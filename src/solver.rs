#![allow(non_snake_case)]

use crate::problem::Bound;
use crate::standard_form::{Basic, BasicFeasiblePoint, Nonbasic, NonbasicBound, StandardForm};
use crate::util::EPS;
use crate::{error::EllPError, problem::Problem};

use log::{debug, trace};

pub struct Solver {}

impl std::default::Default for Solver {
    fn default() -> Self {
        Self {}
    }
}

impl Solver {
    pub fn new() -> Self {
        Self {}
    }

    pub fn solve(
        &self,
        prob: &Problem,
        x0: Option<nalgebra::DVector<f64>>,
    ) -> Result<SolverResult, EllPError> {
        let n = prob.vars().len();
        let x = x0.unwrap_or_else(|| {
            debug!("initial point not provided, defaulting to the zero vector");
            nalgebra::DVector::zeros(n)
        });

        if x.len() != n {
            return Err(EllPError::new(format!(
                "x0 dimensions invalid: {}, expected: {}",
                x.len(),
                n
            )));
        }

        let mut std_form: StandardForm = prob.into();

        let feasible_point = if !prob.is_feasible(x.as_slice()) {
            debug!("initial point not feasible, will find a feasible point");

            match self.find_feasible_point(std_form)? {
                Phase1Result::Feasible(bfp, modified_std_form) => {
                    debug!("found feasible point");
                    std_form = modified_std_form;
                    bfp
                }

                Phase1Result::Infeasible => return Ok(SolverResult::Infeasible),
            }
        } else {
            todo!()
            //self.solve(std_form, feasible_point)
        };

        debug!("initial basic feasible point:\n{:#?}", feasible_point);
        assert!(prob.is_feasible(feasible_point.x.rows(0, prob.vars().len()).as_slice()));
        self.solve_std_form(&std_form, feasible_point)
    }

    fn solve_std_form(
        &self,
        std_form: &StandardForm,
        init: BasicFeasiblePoint, //x0 assumed to be feasible
    ) -> Result<SolverResult, EllPError> {
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

        if N.len() != std_form.cols().checked_sub(std_form.rows()).unwrap() {
            return Err(EllPError::new(format!(
                "invalid N, has {} elements but {} expected",
                N.len(),
                std_form.rows(),
            )));
        }

        let B_cols: Vec<_> = B.iter().map(|i| std_form.A.column(i.index)).collect();
        let mut A_B = nalgebra::DMatrix::from_columns(&B_cols);

        let N_cols: Vec<_> = N.iter().map(|i| std_form.A.column(i.index)).collect();
        let mut A_N = nalgebra::DMatrix::from_columns(&N_cols);

        let mut c_B =
            nalgebra::DVector::from_iterator(B.len(), B.iter().map(|i| std_form.c[i.index]));

        let mut c_N =
            nalgebra::DVector::from_iterator(N.len(), N.iter().map(|i| std_form.c[i.index]));

        //TODO set max iterations for the solver
        loop {
            debug!("obj: {}", std_form.obj(&x));

            //TODO avoid clone, and update LU decomp instead of recomputing it
            let lu_decomp = A_B.clone().lu();

            if lu_decomp.u().diagonal().iter().any(|d| d.abs() < EPS) {
                return Err(EllPError::new(
                    "invalid B, A_B is not invertible".to_string(),
                ));
            }

            //should always have a solution
            let u = lu_decomp.solve(&c_B).unwrap();
            let r = &c_N - A_N.transpose() * u;
            let pivot_result = Self::pivot(&lu_decomp, &mut x, std_form, &r, &mut B, &mut N);

            trace!("pivot result: {:?}", pivot_result);

            let pivot = match pivot_result? {
                PivotResult::Pivot(pivot) => pivot,

                PivotResult::Optimal => {
                    return Ok(SolverResult::Optimal(
                        std_form.c.dot(&x),
                        BasicFeasiblePoint { x, B, N },
                    ));
                }

                PivotResult::Unbounded => return Ok(SolverResult::Unbounded),
            };

            let nonbasic = &mut N[pivot.nonbasic];

            match pivot.basic {
                Some((basic_index, bound)) => {
                    //basic variable becomes nonbasic, and vice versa
                    std::mem::swap(&mut B[basic_index].index, &mut nonbasic.index);

                    nonbasic.bound = bound;

                    for (x_i, y_i) in A_N
                        .column_mut(pivot.nonbasic)
                        .iter_mut()
                        .zip(A_B.column_mut(basic_index).iter_mut())
                    {
                        std::mem::swap(x_i, y_i);
                    }

                    std::mem::swap(&mut c_N[pivot.nonbasic], &mut c_B[basic_index]);
                }

                None => {
                    //the pivot variable stays nonbasic
                    match nonbasic.bound {
                        NonbasicBound::Lower => nonbasic.bound = NonbasicBound::Upper,
                        NonbasicBound::Upper => nonbasic.bound = NonbasicBound::Lower,
                        NonbasicBound::Free => {
                            return Err(EllPError::new("pivot variable cannot be free".to_string()))
                        }
                    }
                }
            }
        }
    }

    fn find_feasible_point<'a>(
        &self,
        std_form: StandardForm<'a>,
    ) -> Result<Phase1Result<'a>, EllPError> {
        let (std_form_phase_1, x0) = std_form.phase_1();

        match self.solve_std_form(&std_form_phase_1.std_form, x0)? {
            SolverResult::Optimal(obj_val, mut bfp) => {
                assert!(obj_val > -EPS);
                if obj_val < EPS {
                    let std_form = std_form_phase_1.phase_2(&mut bfp);
                    Ok(Phase1Result::Feasible(bfp, std_form))
                } else {
                    Ok(Phase1Result::Infeasible)
                }
            }

            SolverResult::Infeasible => Ok(Phase1Result::Infeasible),
            SolverResult::Unbounded => Err(EllPError::new(
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
    ) -> Result<PivotResult, EllPError> {
        //TODO avoid cycles
        let pivot = N
            .iter_mut()
            .enumerate()
            .zip(r.iter())
            .filter_map(|((i, nonbasic), &r_i)| {
                if r_i.abs() < EPS {
                    return None;
                }

                match (r_i > 0., &nonbasic.bound) {
                    (true, NonbasicBound::Upper) => Some((r_i, nonbasic, i)),
                    (false, NonbasicBound::Lower) => Some((-r_i, nonbasic, i)),
                    (_, NonbasicBound::Free) => Some((r_i.abs(), nonbasic, i)),
                    _ => None,
                }
            })
            .max_by(|(r1, _N1, _i1), (r2, _N2, _i2)| r1.partial_cmp(r2).expect("NaN detected"))
            .map(|(_r_i, pivot, i)| (pivot, i));

        let pivot = match pivot {
            Some(pivot) => pivot,
            None => {
                return Ok(PivotResult::Optimal);
            }
        };

        //should always have a solution
        let mut d = lu_decomp.solve(&std_form.A.column(pivot.0.index)).unwrap();
        let nonbasic_at_lower = matches!(pivot.0.bound, NonbasicBound::Lower);

        if nonbasic_at_lower {
            d = -d;
        }

        let mut lambda = match std_form.bounds[pivot.0.index] {
            Bound::Free => f64::INFINITY,
            Bound::Lower(..) => f64::INFINITY,
            Bound::Upper(..) => f64::INFINITY,
            Bound::TwoSided(lb, ub) => ub - lb,
            Bound::Fixed(..) => 0.,
        };

        assert!(d.len() == B.len());

        let mut new_basic = None;

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

                Bound::Fixed(..) => {
                    return Err(EllPError::new(
                        "fixed variable should not be basic".to_string(),
                    ))
                }
            };

            if lambda_i < lambda {
                lambda = lambda_i;

                if d_i > 0. {
                    new_basic = Some((i, NonbasicBound::Upper));
                } else {
                    new_basic = Some((i, NonbasicBound::Lower));
                };
            }
        }

        assert!(lambda >= 0.);

        if lambda > 0. && lambda.is_finite() {
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
            Ok(PivotResult::Pivot(Pivot {
                basic: new_basic,
                nonbasic: pivot.1,
            }))
        } else if lambda == 0. {
            Ok(PivotResult::Optimal)
        } else {
            //lambda == infinity
            Ok(PivotResult::Unbounded)
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
pub enum SolverResult {
    Optimal(f64, BasicFeasiblePoint),
    Infeasible,
    Unbounded,
}

#[derive(Debug)]
pub enum Phase1Result<'a> {
    Feasible(BasicFeasiblePoint, StandardForm<'a>),
    Infeasible,
}
