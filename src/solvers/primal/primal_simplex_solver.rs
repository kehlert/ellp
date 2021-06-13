#![allow(non_snake_case)]

use super::primal_problem::{PrimalFeasiblePoint, PrimalPhase1, PrimalPhase2};
use crate::error::EllPError;
use crate::problem::{Bound, Problem};
use crate::solver::{EllPResult, OptimalPoint, Solution, SolutionStatus, SolverResult};
use crate::solvers::trivial::solve_trivial_problem;
use crate::standard_form::{
    Basic, BasicPoint, Nonbasic, NonbasicBound, StandardForm, StandardizedProblem,
};
use crate::util::{EPS, ITER_WIDTH};

use log::{debug, info, trace};

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

        info!("PRIMAL PHASE 1");

        let mut phase_2: PrimalPhase2 = match self.solve_with_initial(&mut phase_1)? {
            SolutionStatus::Optimal => {
                let obj = phase_1.obj();
                assert!(obj > -EPS);

                if obj < EPS {
                    info!("found feasible point");
                    phase_1.into()
                } else {
                    info!("problem is infeasible");
                    return Ok(SolverResult::Infeasible);
                }
            }

            SolutionStatus::Infeasible => {
                info!("problem is infeasible");
                return Ok(SolverResult::Infeasible);
            }

            SolutionStatus::Unbounded => panic!("primal phase 1 should never be unbounded"),

            SolutionStatus::MaxIter => {
                info!("reached maximum iterations");
                return Ok(SolverResult::MaxIter { obj: f64::INFINITY });
            }
        };

        info!("PRIMAL PHASE 2");

        Ok(match self.solve_with_initial(&mut phase_2)? {
            SolutionStatus::Optimal => {
                let opt_pt = OptimalPoint::new(phase_2.point.into_pt());

                info!(
                    "found optimal point with objective value {}",
                    phase_2.std_form.obj(&opt_pt.x)
                );

                SolverResult::Optimal(Solution::new(phase_2.std_form, opt_pt))
            }

            SolutionStatus::Infeasible => panic!("primal phase 2 should never be infeasible"),

            SolutionStatus::Unbounded => {
                info!("problem is unbounded");
                SolverResult::Unbounded
            }

            SolutionStatus::MaxIter => {
                info!("reached maximum iterations");
                SolverResult::MaxIter { obj: phase_2.obj() }
            }
        })
    }

    pub fn solve_with_initial<P>(&self, prob: &mut P) -> Result<SolutionStatus, EllPError>
    where
        P: StandardizedProblem<FeasiblePoint = PrimalFeasiblePoint>,
    {
        let (std_form, pt) = prob.unpack();

        let pt = &mut **pt;
        let x = &mut pt.x;
        let N = &mut pt.N;
        let B = &mut pt.B;

        info!(
            "solving problem with {} variables and {} constraints",
            std_form.cols(),
            std_form.rows()
        );

        trace!("c: {}", std_form.c);
        trace!("A: {}", std_form.A);
        trace!("b: {}", std_form.b);

        info!("Iteration  |  Objective");

        if std_form.rows() == 0 {
            //trivial problem, and would run into errors if we proceed
            assert!(B.is_empty());
            return Ok(solve_trivial_problem(std_form, x, N, true));
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

        let mut iter = 1u64;

        //TODO set max iterations for the solver
        loop {
            info!("{:it$}  |  {:.8E}", iter, std_form.obj(x), it = ITER_WIDTH,);

            if iter > self.max_iter {
                debug!("reached max iterations");
                return Ok(SolutionStatus::MaxIter);
            }

            iter += 1;

            //TODO check that objective is nonincreasing

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
