#![allow(non_snake_case)]

use std::ops::AddAssign;

use super::dual_problem::{DualPhase1, DualPhase2, DualProblem};
use crate::error::EllPError;
use crate::problem::{Bound, Problem};
use crate::solver::{EllPResult, OptimalPoint, Solution, SolutionStatus, SolverResult};
use crate::standard_form::{BasicPoint, Nonbasic, NonbasicBound};
use crate::util::EPS;
use crate::PrimalSimplexSolver;

use log::debug;

pub struct DualSimplexSolver {
    max_iter: u64,
}

impl std::default::Default for DualSimplexSolver {
    fn default() -> Self {
        Self { max_iter: 1000 }
    }
}

impl DualSimplexSolver {
    pub fn new(max_iter: Option<u64>) -> Self {
        Self {
            max_iter: max_iter.unwrap_or(u64::MAX),
        }
    }

    pub fn solve(&self, prob: Problem) -> EllPResult {
        let mut phase_1: DualPhase1 = match prob.into() {
            Some(phase_1) => phase_1,
            None => return Ok(SolverResult::Infeasible),
        };

        println!("\n---------------------------\nPHASE 1\n---------------------------\n");

        let mut phase_2: DualPhase2 = match self.solve_with_initial(&mut phase_1)? {
            SolutionStatus::Optimal => {
                let obj = phase_1.obj();

                assert!(obj < EPS);

                if obj > -EPS {
                    debug!("found feasible point");
                    phase_1.into()
                } else {
                    let primal_solver = PrimalSimplexSolver::default();
                    let result = primal_solver.solve(phase_1.into_orig_prob())?;

                    assert!(matches!(
                        result,
                        SolverResult::Infeasible
                            | SolverResult::Unbounded
                            | SolverResult::MaxIter { .. }
                    ));

                    return Ok(result);
                }
            }

            SolutionStatus::Infeasible => return Ok(SolverResult::Infeasible),

            SolutionStatus::Unbounded => panic!("phase 1 should never be unbounded"),

            SolutionStatus::MaxIter => {
                debug!("phase 1 reached maximum iterations");
                return Ok(SolverResult::MaxIter { obj: f64::INFINITY });
            }
        };

        println!("\n---------------------------\nPHASE 2\n---------------------------\n");

        //debug!("initial dual feasible point:\n{:#?}", phase_2.pt());

        Ok(match self.solve_with_initial(&mut phase_2)? {
            SolutionStatus::Optimal => {
                let opt_pt = OptimalPoint::new(phase_2.point.into_pt());
                println!("optimal point:\n{:?}", opt_pt);
                SolverResult::Optimal(Solution::new(phase_2.std_form, opt_pt))
            }

            SolutionStatus::Infeasible => SolverResult::Infeasible,
            SolutionStatus::Unbounded => SolverResult::Unbounded,
            SolutionStatus::MaxIter => SolverResult::MaxIter { obj: phase_2.obj() },
        })
    }

    pub fn solve_with_initial<P: DualProblem>(
        &self,
        prob: &mut P,
    ) -> Result<SolutionStatus, EllPError> {
        let (std_form, pt) = prob.unpack();

        // let (x, N, B) = pt.unpack();

        let y = &mut pt.y;
        let mut d = &std_form.c - std_form.A.tr_mul(y);

        let pt = &mut pt.point;
        let x = &mut pt.x;
        let N = &mut pt.N;
        let B = &mut pt.B;

        println!("c:{}", std_form.c);
        println!("A:{}", std_form.A);
        println!("y:{}", y);
        println!("d:{}", d);
        println!("B:{:?}", B);
        println!("N:{:?}", N);
        println!("x: {}", x);

        for n in N.as_slice() {
            let d_i = d[n.index];

            let infeasible = match n.bound {
                NonbasicBound::Lower => d_i < -EPS,
                NonbasicBound::Upper => d_i > EPS,
                NonbasicBound::Free => d_i.abs() > EPS,
            };

            if infeasible {
                panic!("initial point of dual phase 2 is dual infeasible");
            }
        }

        println!("std form:\n{}", std_form);
        println!("x:\n{}", x);
        println!("B:\n{:?}", B);
        println!("N:\n{:?}", N);

        if std_form.rows() == 0 {
            //trivial problem, and would run into errors if we proceed
            assert_eq!(std_form.c.len(), std_form.bounds.len());
            assert_eq!(std_form.c.len(), x.len());

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
                        Bound::Free | Bound::Upper(..) => return Ok(SolutionStatus::Infeasible), //dual unbounded
                        Bound::Lower(lb) | Bound::TwoSided(lb, ..) | Bound::Fixed(lb) => {
                            *x_i = lb;
                            N.push(Nonbasic::new(i, NonbasicBound::Lower));
                        }
                    }
                } else if c_i < 0. {
                    match *bound {
                        Bound::Free | Bound::Lower(..) => return Ok(SolutionStatus::Infeasible), //dual unbounded
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

        let mut A_B = std_form.A.select_columns(B.iter().map(|b| &b.index));
        let mut c_B =
            nalgebra::DVector::from_iterator(B.len(), B.iter().map(|i| std_form.c[i.index]));

        if N.is_empty() {
            return Ok(SolutionStatus::Optimal);
        }

        let mut A_N = std_form.A.select_columns(N.iter().map(|n| &n.index));
        let mut c_N =
            nalgebra::DVector::from_iterator(N.len(), N.iter().map(|i| std_form.c[i.index]));

        println!("A:{}", std_form.A);
        println!("b:{}", std_form.b);
        println!("c:{}", std_form.c);
        println!("bounds:\n{:?}", std_form.bounds);
        println!("A_N:{}", A_N);

        let mut iter = 0u64;
        let mut obj = std_form.dual_obj(y, &d);

        let mut zero_vector = nalgebra::DVector::zeros(std_form.rows());

        loop {
            println!("ITER: {}", iter);

            if iter >= self.max_iter {
                debug!("reached max iterations");
                return Ok(SolutionStatus::MaxIter);
            }

            iter += 1;

            println!("y:{}", y);
            println!("d:{}", d);
            println!("B:{:?}", B);
            println!("N:{:?}", N);
            println!("x: {}", x);

            //TODO check that objective is nondecreasing
            //could do an online update of the objective
            debug!("obj: {}", obj);

            let leaving = B.iter().enumerate().find_map(|(i, B_i)| {
                let x_i = x[B_i.index];

                let bound_change = match std_form.bounds[B_i.index] {
                    Bound::Free => None,

                    Bound::Lower(lb) => {
                        if x_i < lb {
                            Some((x_i - lb, NonbasicBound::Lower))
                        } else {
                            None
                        }
                    }

                    Bound::Upper(ub) => {
                        if x_i > ub {
                            Some((x_i - ub, NonbasicBound::Upper))
                        } else {
                            None
                        }
                    }

                    Bound::TwoSided(lb, ub) => {
                        if x_i > ub {
                            Some((x_i - ub, NonbasicBound::Upper))
                        } else if x_i < lb {
                            Some((x_i - lb, NonbasicBound::Lower))
                        } else {
                            None
                        }
                    }

                    Bound::Fixed(_val) => None,
                };

                bound_change.map(|(delta, nonbasic_bound)| (i, delta, nonbasic_bound))
            });

            debug!("leaving variable: {:?}", leaving);

            //TODO avoid the clone and update LU decomp instead
            let lu = A_B.clone().lu();

            let (leaving_index, delta, nonbasic_bound) = match leaving {
                Some(inner) => inner,
                None => return Ok(SolutionStatus::Optimal),
            };

            zero_vector[leaving_index] = 1.;
            let rho_tilde = lu.u().tr_solve_upper_triangular(&zero_vector).unwrap();
            zero_vector[leaving_index] = 0.;

            let mut rho = lu.l().tr_solve_lower_triangular(&rho_tilde).unwrap();
            lu.p().inv_permute_rows(&mut rho);
            debug!("rho:{}", rho);

            let mut alpha = A_N.tr_mul(&rho);

            if delta < 0. {
                alpha.neg_mut();
            }

            debug!("alpha:{}", alpha);

            assert!(alpha.len() == N.len());

            let entering = N
                .iter()
                .enumerate()
                .filter_map(|(i, N_i)| {
                    let keep = match N_i.bound {
                        NonbasicBound::Lower => alpha[i] > EPS,
                        NonbasicBound::Upper => alpha[i] < -EPS,
                        NonbasicBound::Free => true,
                    };

                    if keep {
                        Some((i, d[N_i.index] / alpha[i]))
                    } else {
                        None
                    }
                })
                .min_by(|&(_i1, r1), &(_i2, r2)| r1.partial_cmp(&r2).unwrap());

            let (entering_index, mut theta_dual) = match entering {
                Some(ent) => ent,
                None => return Ok(SolutionStatus::Infeasible), //dual unbounded
            };

            if delta < 0. {
                alpha.neg_mut();
                theta_dual = -theta_dual;
            }

            let leaving = &B[leaving_index];
            let entering = &N[entering_index];

            debug!("leaving: {:?}", leaving);
            debug!("entering: {:?}", entering);

            let alpha_q = lu.solve(&std_form.A.column(entering.index)).unwrap();

            debug!("alpha_q: {}", alpha_q);

            d[leaving_index] = -theta_dual;

            for (i, n) in N.iter().enumerate() {
                d[n.index] -= theta_dual * alpha[i]
            }

            d[entering_index] = 0.;

            y.add_assign(theta_dual * &rho);

            let theta_primal = delta / alpha_q[leaving_index];

            assert!(alpha_q.len() == B.len());

            for (i, b) in B.iter().enumerate() {
                println!(
                    "x[{}]={}, {}, {}",
                    b.index, x[b.index], theta_primal, alpha_q[i]
                );
                x[b.index] -= theta_primal * alpha_q[i];
                println!("after: {}", x[b.index]);
            }

            x[entering.index] += theta_primal;

            obj += theta_dual * delta;

            debug!("delta: {}", delta);
            debug!("theta_dual: {}", theta_dual);
            debug!("theta_primal: {}", theta_primal);

            debug!("obj after: {}", obj);

            debug!("obj after2: {}", std_form.dual_obj(y, &d));

            debug!("x after: {}", x);

            std::mem::swap(&mut B[leaving_index].index, &mut N[entering_index].index);
            N[entering_index].bound = nonbasic_bound;

            for (x, y) in A_B
                .column_mut(leaving_index)
                .iter_mut()
                .zip(A_N.column_mut(entering_index).iter_mut())
            {
                std::mem::swap(x, y);
            }

            std::mem::swap(&mut c_B[leaving_index], &mut c_N[entering_index]);

            debug!("B after: {:?}", B);
            debug!("N after: {:?}", N);
        }
    }
}
