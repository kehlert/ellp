use std::ops::{AddAssign, DerefMut};

use super::dual_problem::{DualFeasiblePoint, DualPhase1, DualPhase2, DualProblem};
use crate::error::EllPError;
use crate::problem::{Bound, Problem};
use crate::solver::{EllPResult, SolutionStatus};
use crate::standard_form::{Nonbasic, NonbasicBound};
use crate::util::EPS;

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
        let mut phase_1: DualPhase1 = prob.into();
        self.solve_with_initial(&mut phase_1)?;
        todo!()
    }

    pub fn solve_with_initial<P: DualProblem>(
        &self,
        prob: &mut P,
    ) -> Result<SolutionStatus, EllPError> {
        let (std_form, pt) = prob.unpack();

        // let (x, N, B) = pt.unpack();

        let y = &mut pt.y;
        let w = &mut pt.w;
        let v = &mut pt.v;

        let pt = &mut pt.point;
        let x = &mut pt.x;
        let N = &mut pt.N;
        let B = &mut pt.B;

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

        let mut d = &std_form.c - std_form.A.tr_mul(y);

        let mut iter = 0u64;

        let mut zero_vector = nalgebra::DVector::zeros(std_form.rows());

        loop {
            if iter >= self.max_iter {
                debug!("reached max iterations");
                return Ok(SolutionStatus::MaxIter);
            }

            iter += 1;

            //TODO check that objective is nondecreasing
            //could do an online update of the objective
            debug!("obj: {}", std_form.dual_obj(y, v, w));

            let leaving = B.iter().enumerate().find_map(|(i, B_i)| {
                let x_i = x[B_i.index];

                let delta = match std_form.bounds[B_i.index] {
                    Bound::Free => None,

                    Bound::Lower(lb) => {
                        if x_i < lb {
                            Some(x_i - lb)
                        } else {
                            None
                        }
                    }
                    Bound::Upper(ub) => {
                        if x_i > ub {
                            Some(x_i - ub)
                        } else {
                            None
                        }
                    }

                    Bound::TwoSided(lb, ub) => {
                        if x_i > ub {
                            Some(x_i - ub)
                        } else if x_i < lb {
                            Some(x_i - lb)
                        } else {
                            None
                        }
                    }

                    Bound::Fixed(val) => {
                        if (x_i - val).abs() > EPS {
                            Some(x_i - val)
                        } else {
                            None
                        }
                    }
                };

                delta.map(|d| (i, d))
            });

            debug!("leaving variable: {:?}", leaving);

            //TODO avoid the clone and update LU decomp instead
            let lu = A_B.clone().lu();

            let (leaving_index, delta) = match leaving {
                Some((i, delta)) => (i, delta),
                None => return Ok(SolutionStatus::Optimal),
            };

            zero_vector[leaving_index] = 1.;
            let rho_tilde = lu.u().tr_solve_upper_triangular(&zero_vector).unwrap();
            zero_vector[leaving_index] = 0.;

            let mut rho = lu.l().tr_solve_lower_triangular(&rho_tilde).unwrap();
            lu.p().inv_permute_rows(&mut rho);
            debug!("rho:{}", rho);

            let mut alpha = A_N.tr_mul(&rho);
            debug!("alpha:{}", alpha);

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

            debug!("entering: {:?}", N[entering_index]);

            let alpha_q = lu.solve(&std_form.A.column(entering.index)).unwrap();

            debug!("alpha_q: {}", alpha_q);

            d[leaving_index] = -theta_dual;

            for (i, n) in N.iter().enumerate() {
                d[n.index] -= theta_dual * alpha[i]
            }

            d[entering_index] = 0.;

            let theta_primal = delta / alpha_q[leaving_index];

            assert!(alpha_q.len() == B.len());

            for (i, b) in B.iter().enumerate() {
                x[b.index] -= theta_primal * alpha_q[i];
            }

            x[entering.index] += theta_primal;

            //TODO update basics and nonbasics, including A_N, c_N, etc.

            debug!("delta: {}, theta_dual: {}", delta, theta_dual);

            debug!(
                "obj after: {}",
                std_form.dual_obj(y, v, w) + theta_dual * delta
            );

            debug!("x after: {}", x);

            todo!()
        }

        todo!()
    }
}
