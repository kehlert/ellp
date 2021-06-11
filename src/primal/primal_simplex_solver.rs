#![allow(non_snake_case)]

use super::primal_problem::{PrimalPhase1, PrimalPhase2, PrimalProblem};

use crate::error::EllPError;
use crate::problem::{Bound, Problem};
use crate::solver::{EllPResult, OptimalPoint, Solution, SolutionStatus, SolverResult};
use crate::standard_form::{Basic, BasicPoint, Nonbasic, NonbasicBound, StandardForm};
use crate::util::EPS;

use log::{debug, trace};

pub struct PrimalSimplexSolver {
    max_iter: u64,
}

impl std::default::Default for PrimalSimplexSolver {
    fn default() -> Self {
        Self { max_iter: 1000 }
    }
}

impl PrimalSimplexSolver {
    pub fn new(max_iter: Option<u64>) -> Self {
        Self {
            max_iter: max_iter.unwrap_or(u64::MAX),
        }
    }

    pub fn solve(&self, prob: Problem) -> EllPResult {
        let mut phase_1: PrimalPhase1 = match prob.into() {
            Some(phase_1) => phase_1,
            None => return Ok(SolverResult::Infeasible),
        };

        let mut phase_2: PrimalPhase2 = match self.solve_with_initial(&mut phase_1)? {
            SolutionStatus::Optimal => {
                let obj = phase_1.obj();

                assert!(obj > -EPS);

                if obj < EPS {
                    debug!("found feasible point");
                    phase_1.into()
                } else {
                    return Ok(SolverResult::Infeasible);
                }
            }

            SolutionStatus::Infeasible => return Ok(SolverResult::Infeasible),
            SolutionStatus::Unbounded => panic!("phase 1 problem is never unbounded"),
            SolutionStatus::MaxIter => {
                debug!("phase 1 reached maximum iterations");
                return Ok(SolverResult::MaxIter { obj: f64::INFINITY });
            }
        };

        debug!("initial primal feasible point:\n{:#?}", phase_2.pt());

        Ok(match self.solve_with_initial(&mut phase_2)? {
            SolutionStatus::Optimal => {
                let opt_pt = OptimalPoint::new(phase_2.feasible_point.into_pt());
                SolverResult::Optimal(Solution::new(phase_2.std_form, opt_pt))
            }

            SolutionStatus::Infeasible => SolverResult::Infeasible,
            SolutionStatus::Unbounded => SolverResult::Unbounded,
            SolutionStatus::MaxIter => SolverResult::MaxIter { obj: phase_2.obj() },
        })
    }

    pub fn solve_with_initial<P: PrimalProblem>(
        &self,
        prob: &mut P,
    ) -> Result<SolutionStatus, EllPError> {
        let (std_form, pt) = prob.unpack();

        let pt = &mut **pt;
        let x = &mut pt.x;
        let N = &mut pt.N;
        let B = &mut pt.B;

        trace!("c: {}", std_form.c);
        trace!("A: {}", std_form.A);
        trace!("b: {}", std_form.b);

        if std_form.rows() == 0 {
            //trivial problem, and would run into errors if we proceed
            assert_eq!(std_form.c.len(), std_form.bounds.len());
            assert!(B.is_empty());
            N.clear();

            for (i, ((x_i, &c_i), bound)) in x
                .iter_mut()
                .zip(&std_form.c)
                .zip(&std_form.bounds)
                .enumerate()
            {
                if c_i > 0. {
                    match *bound {
                        Bound::Free | Bound::Upper(..) => return Ok(SolutionStatus::Unbounded),
                        Bound::Lower(lb) | Bound::TwoSided(lb, ..) | Bound::Fixed(lb) => {
                            *x_i = lb;
                            N.push(Nonbasic::new(i, NonbasicBound::Lower));
                        }
                    }
                } else if c_i < 0. {
                    match *bound {
                        Bound::Free | Bound::Lower(..) => return Ok(SolutionStatus::Unbounded),
                        Bound::Upper(ub) | Bound::TwoSided(.., ub) | Bound::Fixed(ub) => {
                            *x_i = ub;
                            N.push(Nonbasic::new(i, NonbasicBound::Upper));
                        }
                    }
                } else {
                    *x_i = 0.; //can set to anything, but zero seems reasonable
                }
            }

            return Ok(SolutionStatus::Optimal);
        }

        if B.len() != std_form.rows() {
            return Err(EllPError::new(format!(
                "invalid B, has {} elements but {} expected",
                B.len(),
                std_form.rows(),
            )));
        }

        let expected_N_len = std_form.cols().checked_sub(std_form.rows()).unwrap();

        if N.len() != expected_N_len {
            return Err(EllPError::new(format!(
                "invalid N, has {} elements but {} expected",
                N.len(),
                expected_N_len,
            )));
        }

        let B_cols: Vec<_> = B.iter().map(|i| std_form.A.column(i.index)).collect();
        let mut A_B = nalgebra::DMatrix::from_columns(&B_cols);
        let mut c_B =
            nalgebra::DVector::from_iterator(B.len(), B.iter().map(|i| std_form.c[i.index]));

        let N_cols: Vec<_> = N.iter().map(|i| std_form.A.column(i.index)).collect();

        if N_cols.is_empty() {
            return Ok(SolutionStatus::Optimal);
        }

        let mut A_N = nalgebra::DMatrix::from_columns(&N_cols);
        let mut c_N =
            nalgebra::DVector::from_iterator(N.len(), N.iter().map(|i| std_form.c[i.index]));

        let mut iter = 0u64;

        //TODO set max iterations for the solver
        loop {
            if iter >= self.max_iter {
                debug!("reached max iterations");
                return Ok(SolutionStatus::MaxIter);
            }

            iter += 1;

            //TODO check that objective is nonincreasing
            debug!("obj: {}", std_form.obj(&x));

            //TODO avoid clone, and update LU decomp instead of recomputing it
            let lu_decomp = A_B.clone().lu();

            if lu_decomp.u().diagonal().iter().any(|d| d.abs() < EPS) {
                return Err(EllPError::new(
                    "invalid B, A_B is not invertible".to_string(),
                ));
            }

            //should always have a solution
            //(perhaps add solve_transpose to lu decomp? or .transpose() to lu_decomp?)
            //would be a good contribution to nalgebra
            let u_tilde = lu_decomp.u().tr_solve_upper_triangular(&c_B).unwrap();
            let mut u = lu_decomp.l().tr_solve_lower_triangular(&u_tilde).unwrap();

            lu_decomp.p().inv_permute_rows(&mut u);

            let r = &c_N - A_N.tr_mul(&u);

            let pivot_result = Self::pivot(&lu_decomp, &r, std_form, x, N, B);

            trace!("pivot result: {:?}", pivot_result);

            let pivot = match pivot_result {
                PivotResult::Pivot(pivot) => pivot,

                PivotResult::Optimal => {
                    return Ok(SolutionStatus::Optimal);
                }

                PivotResult::Unbounded => return Ok(SolutionStatus::Unbounded),
            };

            let nonbasic = &mut N[pivot.nonbasic];

            match pivot.basic {
                Some((basic_index, bound)) => {
                    std::mem::swap(&mut B[basic_index].index, &mut nonbasic.index);

                    for (x_i, y_i) in A_N
                        .column_mut(pivot.nonbasic)
                        .iter_mut()
                        .zip(A_B.column_mut(basic_index).iter_mut())
                    {
                        std::mem::swap(x_i, y_i);
                    }

                    std::mem::swap(&mut c_N[pivot.nonbasic], &mut c_B[basic_index]);
                    nonbasic.bound = bound;
                }

                None => {
                    let bound = &mut nonbasic.bound;
                    //nonbasic went from lower to upper bound, or vice versa
                    match bound {
                        NonbasicBound::Lower => *bound = NonbasicBound::Upper,
                        NonbasicBound::Upper => *bound = NonbasicBound::Lower,
                        NonbasicBound::Free => panic!("pivot should have been unbounded"),
                    }
                }
            }

            //basic variable becomes nonbasic, and vice versa
        }
    }

    fn pivot(
        lu_decomp: &nalgebra::LU<f64, nalgebra::Dynamic, nalgebra::Dynamic>,
        r: &nalgebra::DVector<f64>,
        std_form: &StandardForm,
        x: &mut nalgebra::DVector<f64>,
        N: &mut [Nonbasic],
        B: &mut [Basic],
    ) -> PivotResult {
        trace!("\n-------------------\nobj: {}", std_form.obj(x));
        trace!("bounds: {:?}", std_form.bounds);
        trace!("x: {}", x);
        trace!("B: {:?}", B);
        trace!("N: {:?}", N);
        trace!("r: {}", r);

        let pivot = N
            .iter_mut()
            .enumerate()
            .zip(r.iter())
            .filter_map(|((i, nonbasic), &r_i)| {
                if r_i.abs() < EPS {
                    return None;
                }

                trace!("{}, {}, {:?}", r_i, r_i > 0., nonbasic.bound);

                match (r_i > 0., &nonbasic.bound) {
                    (true, NonbasicBound::Upper) => Some((r_i, nonbasic, i)),
                    (false, NonbasicBound::Lower) => Some((-r_i, nonbasic, i)),
                    (_, NonbasicBound::Free) => Some((r_i.abs(), nonbasic, i)),
                    _ => None,
                }
            })
            .max_by(|(r1, N1, _i1), (r2, N2, _i2)| {
                trace!(
                    "{}, {}, {}, {:?}",
                    r1,
                    r2,
                    (r1 - r2).abs(),
                    r1.partial_cmp(r2)
                );

                //this logic breaks cycles (smallest subscript rule)
                if (r1 - r2).abs() >= EPS {
                    r1.partial_cmp(r2).expect("NaN detected")
                } else {
                    N1.index.cmp(&N2.index)
                }
            })
            .map(|(_r_i, pivot, i)| (pivot, i));

        let pivot = match pivot {
            Some(pivot) => pivot,
            None => return PivotResult::Optimal,
        };

        //should always have a solution
        let mut d = lu_decomp.solve(&std_form.A.column(pivot.0.index)).unwrap();
        let nonbasic_at_lower = matches!(pivot.0.bound, NonbasicBound::Lower);

        if nonbasic_at_lower {
            d = -d;
        }

        trace!("d: {}", d);
        trace!("pivot: {:?}", pivot);

        let mut lambda = match std_form.bounds[pivot.0.index] {
            Bound::Free => f64::INFINITY,
            Bound::Lower(..) => f64::INFINITY,
            Bound::Upper(..) => f64::INFINITY,
            Bound::TwoSided(lb, ub) => ub - lb,
            Bound::Fixed(..) => 0.,
        };

        assert!(d.len() == B.len());

        let mut new_basic = None;
        let mut new_basic_index = None;

        assert!(B.len() == d.len());

        for (i, (basic, &d_i)) in B.iter_mut().zip(d.iter()).enumerate() {
            if d_i.abs() < EPS {
                continue;
            }

            let x_i = x[basic.index];

            let lambda_i = match std_form.bounds[basic.index] {
                Bound::Free => f64::INFINITY,

                Bound::Lower(lb) => {
                    if d_i > 0. {
                        f64::INFINITY
                    } else {
                        (lb - x_i) / d_i
                    }
                }

                Bound::Upper(ub) => {
                    if d_i > 0. {
                        (ub - x_i) / d_i
                    } else {
                        f64::INFINITY
                    }
                }

                Bound::TwoSided(lb, ub) => {
                    if d_i > 0. {
                        (ub - x_i) / d_i
                    } else {
                        (lb - x_i) / d_i
                    }
                }

                Bound::Fixed(..) => 0.,
            };

            trace!("i: {}, lambda_i: {}, lambda: {}", i, lambda_i, lambda);

            if lambda_i < lambda - EPS {
                lambda = lambda_i;

                if d_i > 0. {
                    new_basic = Some((i, NonbasicBound::Upper));
                } else {
                    new_basic = Some((i, NonbasicBound::Lower));
                };
            } else if (lambda_i - lambda).abs() < EPS {
                //this logic breaks cycles (smallest subscript rule)
                if new_basic_index.is_none() || basic.index < new_basic_index.unwrap() {
                    new_basic_index = Some(basic.index);
                    lambda = lambda_i;

                    if d_i > 0. {
                        new_basic = Some((i, NonbasicBound::Upper));
                    } else {
                        new_basic = Some((i, NonbasicBound::Lower));
                    };
                }
            }
        }

        assert!(lambda >= 0.);

        if lambda.is_infinite() {
            return PivotResult::Unbounded;
        }

        if lambda > 0. {
            for (B_i, d_i) in B.iter().zip(d.iter()) {
                x[B_i.index] += lambda * d_i;
            }

            if nonbasic_at_lower {
                x[pivot.0.index] += lambda;
            } else {
                x[pivot.0.index] -= lambda;
            }

            //to access the basic varbiable, need to use an index rather than a reference
            //because of Rust's borrowing rules
            PivotResult::Pivot(Pivot {
                basic: new_basic,
                nonbasic: pivot.1,
            })
        } else if lambda == 0. {
            PivotResult::Pivot(Pivot {
                basic: new_basic,
                nonbasic: pivot.1,
            })
        } else {
            //lambda == infinity
            PivotResult::Unbounded
        }
    }
}

#[derive(Debug)]
enum PivotResult {
    Pivot(Pivot),
    Optimal,
    Unbounded,
}

#[derive(Debug)]
struct Pivot {
    basic: Option<(usize, NonbasicBound)>,
    nonbasic: usize,
}
