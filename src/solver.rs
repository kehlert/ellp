#![allow(non_snake_case)]

use crate::problem::Bound;
use crate::standard_form::{Basic, BasicFeasiblePoint, Nonbasic, NonbasicBound, StandardForm};
use crate::util::EPS;
use crate::{error::EllPError, problem::Problem};

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

    pub fn solve(&self, prob: Problem, x0: Option<nalgebra::DVector<f64>>) -> SolverResult {
        let n = prob.vars().len();
        let mut x = x0.unwrap_or_else(|| nalgebra::DVector::zeros(n));

        if x.len() != n {
            return SolverResult::Error(EllPError::new(format!(
                "x0 dimensions invalid: {}, expected: {}",
                x.len(),
                n
            )));
        }

        println!("initial x: {}", x);

        if !prob.is_feasible(x.as_slice()) {
            let (result, std_form) = self.find_feasible_point(&prob);

            let feasible_point = match result {
                SolverResult::Optimal(obj_val, bfp) => {
                    assert!(obj_val > -EPS);
                    if obj_val < EPS {
                        bfp
                    } else {
                        return SolverResult::Infeasible;
                    }
                }

                SolverResult::Infeasible => return SolverResult::Infeasible,
                SolverResult::Unbounded => panic!("phase 1 problem is never unbounded"),
                SolverResult::Error(err) => return SolverResult::Error(err),
            };

            println!("feasible x: {}", feasible_point.x);
            assert!(prob.is_feasible(feasible_point.x.rows(0, prob.vars().len()).as_slice()));
            // println!("{:#?}", std_form);
            self.solve_std_form(&std_form, feasible_point)
        } else {
            todo!()
            //self.solve(std_form, feasible_point)
        }
    }

    fn solve_std_form(
        &self,
        std_form: &StandardForm,
        init: BasicFeasiblePoint, //x0 assumed to be feasible
    ) -> SolverResult {
        let mut B = init.B;
        let mut N = init.N;
        let mut x = init.x;

        if B.len() != std_form.rows() {
            return SolverResult::Error(EllPError::new(format!(
                "invalid B, has {} elements but {} expected",
                B.len(),
                std_form.rows(),
            )));
        }

        if N.len() != std_form.cols().checked_sub(std_form.rows()).unwrap() {
            return SolverResult::Error(EllPError::new(format!(
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

        for _i in 0..3 {
            //TODO avoid clone, and update LU decomp instead of recomputing it
            let lu_decomp = A_B.clone().lu();

            let L = lu_decomp.l();
            let U = lu_decomp.u();

            if U.diagonal().iter().any(|d| d.abs() < EPS) {
                return SolverResult::Error(EllPError::new(
                    "invalid B, A_B is not invertible".to_string(),
                ));
            }

            //should always have a solution
            let u = lu_decomp.solve(&c_B).unwrap();
            println!("u: {}", u);

            println!(
                "{}, ({}, {}), {}",
                c_N.nrows(),
                A_N.nrows(),
                A_N.ncols(),
                u.nrows()
            );

            let r = &c_N - A_N.transpose() * u;

            //TODO find nonbasic to enter the basis

            println!("r: {}", r);

            println!("old x: {}", x);
            println!("old obj: {}", std_form.c.dot(&x));

            println!("old nonbasic: {:?}", N);
            println!("old basic: {:?}", B);

            let pivot_result = Self::pivot(&lu_decomp, &mut x, std_form, &r, &mut B, &mut N);

            println!("new x: {}", x);
            println!("new obj: {}", std_form.c.dot(&x));

            println!("old A_B: {}", A_B);
            println!("old A_N: {}", A_N);

            println!("pivot result: {:?}", pivot_result);

            let pivot = match pivot_result {
                PivotResult::Pivot(pivot) => pivot,

                PivotResult::Optimal => {
                    return SolverResult::Optimal(
                        std_form.c.dot(&x),
                        BasicFeasiblePoint { x, B, N },
                    );
                }

                PivotResult::Unbounded => return SolverResult::Unbounded,
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

                    println!("c_N: {}", c_N.len());
                    println!("c_B: {}", c_B.len());

                    std::mem::swap(&mut c_N[pivot.nonbasic], &mut c_B[basic_index]);
                }

                None => {
                    //the pivot variable stays nonbasic
                    match nonbasic.bound {
                        NonbasicBound::Lower => nonbasic.bound = NonbasicBound::Upper,
                        NonbasicBound::Upper => nonbasic.bound = NonbasicBound::Lower,
                    }
                }
            }

            println!("A_B: {}", A_B);
            println!("A_N: {}", A_N);

            // println!("d: {}", d);
            // println!("lambda: {}", lambda);
            println!("new_nonbasic: {:?}", N);
            println!("new_basic: {:?}", B);
        }

        todo!()

        // println!("new x: {}", x);
        // println!("new obj: {}", std_form.c.dot(&x));
        // //TODO monitor infeasibility

        // todo!()
    }

    fn find_feasible_point<'a>(&self, prob: &'a Problem) -> (SolverResult, StandardForm<'a>) {
        let (phase_1_std_form, x0) = StandardForm::phase_1(prob);

        println!("{:#?}", x0);

        println!("lhs: {}", &phase_1_std_form.A * &x0.x);
        println!("b: {}", phase_1_std_form.b);

        println!("A: {}", phase_1_std_form.A);
        println!("bounds:\n{:?}", phase_1_std_form.bounds);

        let result = self.solve_std_form(&phase_1_std_form, x0);
        let phase_2_std_form = phase_1_std_form.phase_2();

        (result, phase_2_std_form)
    }

    fn pivot<'a>(
        lu_decomp: &nalgebra::LU<f64, nalgebra::Dynamic, nalgebra::Dynamic>,
        x: &mut nalgebra::DVector<f64>,
        std_form: &StandardForm,
        r: &nalgebra::DVector<f64>,
        B: &'a mut [Basic],
        N: &'a mut [Nonbasic],
    ) -> PivotResult {
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
                    _ => None,
                }
            })
            .max_by(|(r1, _N1, i1), (r2, _N2, i2)| r1.partial_cmp(r2).expect("NaN detected"))
            .map(|(_r_i, pivot, i)| (pivot, i));

        let pivot = match pivot {
            Some(pivot) => pivot,
            None => {
                return PivotResult::Optimal;
            }
        };

        println!("pivot: {:?}", pivot);
        println!("r: {}", r);

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

                Bound::Fixed(..) => panic!("fixed variable should not be basic"),
            };

            if lambda_i < lambda {
                lambda = lambda_i;

                if d_i > 0. {
                    new_basic.replace((i, NonbasicBound::Upper));
                } else {
                    new_basic.replace((i, NonbasicBound::Lower));
                };
            }
        }

        assert!(lambda >= 0.);

        println!("lambda: {}", lambda);

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
            PivotResult::Pivot(Pivot {
                basic: new_basic,
                nonbasic: pivot.1,
            })
        } else if lambda == 0. {
            PivotResult::Optimal
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
pub enum SolverResult {
    Optimal(f64, BasicFeasiblePoint),
    Infeasible,
    Unbounded,
    Error(EllPError),
}
