#![allow(non_snake_case)]

use std::ops::AddAssign;

use super::dual_problem::{DualFeasiblePoint, DualPhase1, DualPhase2};
use crate::error::EllPError;
use crate::problem::{Bound, Problem};
use crate::solver::{EllPResult, OptimalPoint, Solution, SolutionStatus, SolverResult};
use crate::solvers::trivial::solve_trivial_problem;
use crate::standard_form::{BasicPoint, NonbasicBound, StandardizedProblem};
use crate::util::{EPS, ITER_WIDTH};
use crate::PrimalSimplexSolver;

use log::{debug, info};

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

        info!("DUAL PHASE 1");

        let mut phase_2: DualPhase2 = match self.solve_with_initial(&mut phase_1)? {
            SolutionStatus::Optimal => {
                let obj = phase_1.obj();

                assert!(obj < EPS);

                if obj > -EPS {
                    info!("found feasible point");
                    phase_1.into()
                } else {
                    info!(
                        "dual problem is infeasible, using primal solver to determine \
                         if problem is infeasible or unbounded"
                    );

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

            SolutionStatus::Infeasible => {
                panic!("dual phase 1 should never be infeasible");
            }

            SolutionStatus::Unbounded => panic!("dual phase 1 should never be unbounded"),

            SolutionStatus::MaxIter => {
                info!("reached maximum iterations");
                return Ok(SolverResult::MaxIter { obj: f64::INFINITY });
            }
        };

        info!("DUAL PHASE 2");

        Ok(match self.solve_with_initial(&mut phase_2)? {
            SolutionStatus::Optimal => {
                let opt_pt = OptimalPoint::new(phase_2.point.into_pt());

                info!(
                    "found optimal point with objective value {}",
                    phase_2.std_form.obj(&opt_pt.x)
                );

                SolverResult::Optimal(Solution::new(phase_2.std_form, opt_pt))
            }

            SolutionStatus::Infeasible => {
                info!("problem is infeasible");
                SolverResult::Infeasible
            }

            SolutionStatus::Unbounded => panic!("dual phase 2 should never return unbounded"),

            SolutionStatus::MaxIter => {
                info!("reached maximum iterations");
                SolverResult::MaxIter { obj: phase_2.obj() }
            }
        })
    }

    pub fn solve_with_initial<P>(&self, prob: &mut P) -> Result<SolutionStatus, EllPError>
    where
        P: StandardizedProblem<FeasiblePoint = DualFeasiblePoint>,
    {
        let (std_form, pt) = prob.unpack();

        info!(
            "solving problem with {} variables and {} constraints",
            std_form.cols(),
            std_form.rows()
        );

        let y = &mut pt.y;
        let d = &mut pt.d;

        let pt = &mut pt.point;
        let x = &mut pt.x;
        let N = &mut pt.N;
        let B = &mut pt.B;

        info!("Iteration  |  Objective");

        if std_form.rows() == 0 {
            //trivial problem, and would run into errors if we proceed
            assert!(B.is_empty());
            return Ok(solve_trivial_problem(std_form, x, N, true));
        }

        //panics if std_form.rows() == 0
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

        let mut iter = 0u64;
        let mut obj = std_form.dual_obj(y, d);

        let mut zero_vector = nalgebra::DVector::zeros(std_form.rows());

        loop {
            debug!("{:it$}  |  {:.8E}", iter, obj, it = ITER_WIDTH,);

            if iter >= self.max_iter {
                debug!("reached max iterations");
                return Ok(SolutionStatus::MaxIter);
            }

            iter += 1;

            //TODO check that objective is nondecreasing

            let leaving = B.iter().enumerate().find_map(|(i, B_i)| {
                let x_i = x[B_i.index];

                let bound_change = match std_form.bounds[B_i.index] {
                    Bound::Free => None,

                    Bound::Lower(lb) => {
                        if x_i < lb - EPS {
                            Some((x_i - lb, NonbasicBound::Lower))
                        } else {
                            None
                        }
                    }

                    Bound::Upper(ub) => {
                        if x_i > ub + EPS {
                            Some((x_i - ub, NonbasicBound::Upper))
                        } else {
                            None
                        }
                    }

                    Bound::TwoSided(lb, ub) => {
                        if x_i > ub + EPS {
                            Some((x_i - ub, NonbasicBound::Upper))
                        } else if x_i < lb - EPS {
                            Some((x_i - lb, NonbasicBound::Lower))
                        } else {
                            None
                        }
                    }

                    Bound::Fixed(_val) => None,
                };

                bound_change.map(|(delta, nonbasic_bound)| (i, delta, nonbasic_bound))
            });

            //debug!("leaving variable: {:?}", leaving);

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

            let mut alpha = A_N.tr_mul(&rho);

            if delta < 0. {
                alpha.neg_mut();
            }

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

            let alpha_q = lu.solve(&std_form.A.column(entering.index)).unwrap();

            d[leaving.index] = -theta_dual;

            for (i, n) in N.iter().enumerate() {
                d[n.index] -= theta_dual * alpha[i]
            }

            d[entering.index] = 0.;

            y.add_assign(theta_dual * &rho);

            let theta_primal = delta / alpha_q[leaving_index];

            assert!(alpha_q.len() == B.len());

            for (i, b) in B.iter().enumerate() {
                x[b.index] -= theta_primal * alpha_q[i];
            }

            x[entering.index] += theta_primal;

            obj += theta_dual * delta;

            // debug!("delta: {}", delta);
            // debug!("theta_dual: {}", theta_dual);
            // debug!("theta_primal: {}", theta_primal);

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
        }
    }
}
