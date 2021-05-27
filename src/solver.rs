#![allow(non_snake_case)]

use crate::problem::Bound;
use crate::standard_form::{
    Basic, BasicFeasiblePoint, Nonbasic, NonbasicBound, StandardForm, StandardFormPhase1,
};
use crate::util::EPS;
use crate::{error::EllPError, problem::Problem};

use log::{debug, trace};

pub type EllPResult = Result<SolverResult, EllPError>;

pub struct SimplexSolver {
    max_iter: u64,
    dual: bool,
}

impl std::default::Default for SimplexSolver {
    fn default() -> Self {
        Self {
            max_iter: 1000,
            dual: false,
        }
    }
}

impl SimplexSolver {
    pub fn new(max_iter: Option<u64>, dual: bool) -> Self {
        Self {
            max_iter: max_iter.unwrap_or(u64::MAX),
            dual,
        }
    }

    pub fn solve(&self, prob: &Problem) -> EllPResult {
        if !self.dual {
            self.solve_primal(prob)
        } else {
            self.solve_dual(prob)
        }
    }

    fn solve_primal(&self, prob: &Problem) -> EllPResult {
        let mut std_form: StandardForm = prob.into();

        debug!("finding an initial basic feasible point");

        let feasible_point = match self.find_feasible_point(std_form)? {
            Phase1Result::Feasible(bfp, modified_std_form) => {
                debug!("found feasible point");
                std_form = modified_std_form;
                bfp
            }

            Phase1Result::MaxIter(bfp, phase_1_std_form) => {
                debug!("phase 1 reached maximum iterations");
                let obj = phase_1_std_form.obj(&bfp.x);
                return Ok(SolverResult::MaxIter(Solution { obj, x: bfp.x }));
            }

            Phase1Result::Infeasible => return Ok(SolverResult::Infeasible),
        };

        debug!("initial basic feasible point:\n{:#?}", feasible_point);

        assert!(prob.is_feasible(feasible_point.x().rows(0, prob.vars().len()).as_slice()));

        self.solve_std_form(&std_form, feasible_point)
            .map(|result| match result {
                StandardFormResult::Optimal(bfp) => {
                    let obj = std_form.obj(&bfp.x);
                    let x = std_form.extract_solution(&bfp);
                    SolverResult::Optimal(Solution { obj, x })
                }

                StandardFormResult::MaxIter(bfp) => {
                    let obj = std_form.obj(&bfp.x);
                    let x = std_form.extract_solution(&bfp);
                    SolverResult::MaxIter(Solution { obj, x })
                }

                StandardFormResult::Infeasible => SolverResult::Infeasible,
                StandardFormResult::Unbounded => SolverResult::Unbounded,
            })
    }

    fn solve_dual(&self, prob: &Problem) -> EllPResult {
        todo!()
    }

    fn solve_std_form(
        &self,
        std_form: &StandardForm,
        mut init: BasicFeasiblePoint, //x0 assumed to be feasible
    ) -> Result<StandardFormResult, EllPError> {
        trace!("c: {}", std_form.c);
        trace!("A: {}", std_form.A);
        trace!("b: {}", std_form.b);

        if std_form.rows() == 0 {
            //trivial problem, and would run into errors if we proceed
            assert_eq!(std_form.c.len(), std_form.bounds.len());
            assert!(init.B.is_empty());
            let x = &mut init.x;
            let N = &mut init.N;
            N.clear();

            for (i, ((x_i, &c_i), bound)) in x
                .iter_mut()
                .zip(&std_form.c)
                .zip(&std_form.bounds)
                .enumerate()
            {
                if c_i > 0. {
                    match *bound {
                        Bound::Free | Bound::Upper(..) => return Ok(StandardFormResult::Unbounded),
                        Bound::Lower(lb) | Bound::TwoSided(lb, ..) | Bound::Fixed(lb) => {
                            *x_i = lb;
                            N.push(Nonbasic::new(i, NonbasicBound::Lower));
                        }
                    }
                } else if c_i < 0. {
                    match *bound {
                        Bound::Free | Bound::Lower(..) => return Ok(StandardFormResult::Unbounded),
                        Bound::Upper(ub) | Bound::TwoSided(.., ub) | Bound::Fixed(ub) => {
                            *x_i = ub;
                            N.push(Nonbasic::new(i, NonbasicBound::Upper));
                        }
                    }
                } else {
                    *x_i = 0.; //can set to anything, but zero seems reasonable
                }
            }

            return Ok(StandardFormResult::Optimal(init));
        }

        let mut B = init.B;
        let mut N = init.N;
        let mut x = init.x;

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
            return Ok(StandardFormResult::Optimal(BasicFeasiblePoint { x, B, N }));
        }

        let mut A_N = nalgebra::DMatrix::from_columns(&N_cols);
        let mut c_N =
            nalgebra::DVector::from_iterator(N.len(), N.iter().map(|i| std_form.c[i.index]));

        let mut iter = 0u64;

        //TODO set max iterations for the solver
        loop {
            if iter >= self.max_iter {
                debug!("reached max iterations");
                break Ok(StandardFormResult::MaxIter(BasicFeasiblePoint { x, B, N }));
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

            //should always have a solution\
            //(perhaps add solve_transpose to lu decomp? or .transpose() to lu_decomp?)
            //would be a good contribution to nalgebra
            let u_tilde = lu_decomp
                .u()
                .transpose()
                .solve_lower_triangular(&c_B)
                .unwrap();

            let mut u = lu_decomp
                .l()
                .transpose()
                .solve_upper_triangular(&u_tilde)
                .unwrap();

            lu_decomp.p().inv_permute_rows(&mut u);

            let r = &c_N - A_N.transpose() * u;

            let pivot_result = Self::pivot(&lu_decomp, &mut x, std_form, &r, &mut B, &mut N);

            trace!("pivot result: {:?}", pivot_result);

            let pivot = match pivot_result {
                PivotResult::Pivot(pivot) => pivot,

                PivotResult::Optimal => {
                    return Ok(StandardFormResult::Optimal(BasicFeasiblePoint { x, B, N }));
                }

                PivotResult::Unbounded => return Ok(StandardFormResult::Unbounded),
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

    fn find_feasible_point<'a>(
        &self,
        std_form: StandardForm<'a>,
    ) -> Result<Phase1Result<'a>, EllPError> {
        let (std_form_phase_1, x0) = std_form.phase_1();

        match self.solve_std_form(&std_form_phase_1.std_form, x0)? {
            StandardFormResult::Optimal(mut bfp) => {
                let obj = std_form_phase_1.obj(&bfp.x);

                assert!(obj > -EPS);

                if obj < EPS {
                    let std_form = std_form_phase_1.phase_2(&mut bfp);
                    Ok(Phase1Result::Feasible(bfp, std_form))
                } else {
                    Ok(Phase1Result::Infeasible)
                }
            }

            StandardFormResult::MaxIter(bfp) => Ok(Phase1Result::MaxIter(bfp, std_form_phase_1)),
            StandardFormResult::Infeasible => Ok(Phase1Result::Infeasible),
            StandardFormResult::Unbounded => Err(EllPError::new(
                "phase 1 problem is never unbounded".to_string(),
            )),
        }
    }

    fn pivot<'a>(
        lu_decomp: &nalgebra::LU<f64, nalgebra::Dynamic, nalgebra::Dynamic>,
        x: &mut nalgebra::DVector<f64>,
        std_form: &StandardForm,
        r: &nalgebra::DVector<f64>,
        B: &'a mut [Basic],
        N: &'a mut [Nonbasic],
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

#[derive(Debug)]
pub enum StandardFormResult {
    Optimal(BasicFeasiblePoint),
    Infeasible,
    Unbounded,
    MaxIter(BasicFeasiblePoint),
}

#[derive(Debug)]
pub enum SolverResult {
    Optimal(Solution),
    Infeasible,
    Unbounded,
    MaxIter(Solution),
}

#[derive(Debug)]
pub struct Solution {
    pub obj: f64,
    pub x: nalgebra::DVector<f64>,
}

#[derive(Debug)]
pub enum Phase1Result<'a> {
    Feasible(BasicFeasiblePoint, StandardForm<'a>),
    Infeasible,
    MaxIter(BasicFeasiblePoint, StandardFormPhase1<'a>),
}
