#![allow(non_snake_case)]

use crate::problem::{Bound, Problem};
use crate::solver::OptimalPoint;
use crate::standard_form::{
    Basic, BasicPoint, Nonbasic, NonbasicBound, Point, StandardForm, StandardizedProblem,
};
use crate::util::EPS;

use log::debug;

#[derive(Debug, Clone)]
pub struct PrimalFeasiblePoint(Point);

impl BasicPoint for PrimalFeasiblePoint {
    #[inline]
    fn into_pt(self) -> Point {
        self.0
    }
}

impl std::ops::Deref for PrimalFeasiblePoint {
    type Target = Point;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl std::ops::DerefMut for PrimalFeasiblePoint {
    #[inline]
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

#[derive(Debug)]
pub struct PrimalPhase1 {
    pub std_form: StandardForm,
    pub point: PrimalFeasiblePoint,
    phase_1_vars: Vec<usize>,
}

impl StandardizedProblem for PrimalPhase1 {
    type FeasiblePoint = PrimalFeasiblePoint;

    #[inline]
    fn obj(&self) -> f64 {
        self.std_form.obj(&self.point.x)
    }

    #[inline]
    fn unpack(&mut self) -> (&StandardForm, &mut Self::FeasiblePoint) {
        (&self.std_form, &mut self.point)
    }
}

#[derive(Debug)]
pub struct PrimalPhase2 {
    pub std_form: StandardForm,
    pub point: PrimalFeasiblePoint,
}

impl StandardizedProblem for PrimalPhase2 {
    type FeasiblePoint = PrimalFeasiblePoint;

    #[inline]
    fn obj(&self) -> f64 {
        self.std_form.obj(&self.point.x)
    }

    #[inline]
    fn unpack(&mut self) -> (&StandardForm, &mut PrimalFeasiblePoint) {
        (&self.std_form, &mut self.point)
    }
}

//returns None if infeasible
impl std::convert::From<Problem> for Option<PrimalPhase1> {
    fn from(prob: Problem) -> Self {
        let mut std_form: StandardForm = match prob.into() {
            Some(std_form) => std_form,
            None => return None,
        };

        debug!("converting standard form to phase 1 standard form");

        let n = std_form.cols();
        let m = std_form.rows();
        let mut N = Vec::with_capacity(n);
        let mut B = Vec::with_capacity(n.checked_sub(m).unwrap_or(10));
        let mut v = nalgebra::DVector::<f64>::zeros(std_form.cols());

        for (i, bound) in std_form.bounds.iter().enumerate() {
            match *bound {
                Bound::Free => (), //will set these values later

                Bound::Lower(lb) => {
                    v[i] = lb;

                    N.push(Nonbasic {
                        index: i,
                        bound: NonbasicBound::Lower,
                    });
                }

                Bound::Upper(ub) => {
                    v[i] = ub;

                    N.push(Nonbasic {
                        index: i,
                        bound: NonbasicBound::Upper,
                    });
                }

                Bound::TwoSided(lb, _ub) => {
                    v[i] = lb;

                    N.push(Nonbasic {
                        index: i,
                        bound: NonbasicBound::Lower,
                    });
                }

                Bound::Fixed(fixed_val) => {
                    v[i] = fixed_val;

                    N.push(Nonbasic {
                        index: i,
                        bound: NonbasicBound::Lower,
                    });
                }
            }
        }

        for c_i in &mut std_form.c {
            *c_i = 0.;
        }

        std_form.c = std_form.c.resize_vertically(n + m, 1.);

        let free_vars: Vec<_> = std_form
            .bounds
            .iter()
            .enumerate()
            .filter_map(|(i, &bound)| {
                if matches!(bound, Bound::Free) {
                    Some(i)
                } else {
                    None
                }
            })
            .collect();

        let mut free_vars = nalgebra::DVector::from_vec(free_vars);

        if !free_vars.is_empty() && !std_form.A.is_empty() {
            let cols: Vec<_> = free_vars.iter().map(|&i| std_form.A.column(i)).collect();
            let A_F = nalgebra::DMatrix::from_columns(&cols);

            let lu_decomp = A_F.full_piv_lu();
            let rank = lu_decomp
                .u()
                .diagonal()
                .iter()
                .enumerate()
                .find(|(_i, d)| d.abs() < EPS)
                .map(|(i, _d)| i)
                .unwrap_or_else(|| free_vars.len());

            let P = lu_decomp.p();
            let Q = lu_decomp.q();

            Q.permute_rows(&mut free_vars);

            for &i in free_vars.iter().take(rank) {
                B.push(Basic::new(i));
            }

            for &i in free_vars.iter().skip(rank) {
                std_form.bounds[i] = Bound::Fixed(0.);

                N.push(Nonbasic {
                    index: i,
                    bound: NonbasicBound::Lower,
                });
            }

            //don't need to use all the rows, but keeping it simpler for now
            let mut b_tilde = &std_form.b - &std_form.A * &v;
            P.permute_rows(&mut b_tilde);
            let len = b_tilde.len();
            let b_tilde = b_tilde.remove_rows(rank, len - rank);

            let L = lu_decomp.l();
            let L = L.slice((0, 0), (rank, rank));

            let U = lu_decomp.u();
            let U = U.slice((0, 0), (rank, rank));

            //TODO might not have a solution?
            let solution = U
                .solve_upper_triangular(&L.solve_lower_triangular(&b_tilde).unwrap())
                .unwrap();

            assert!(solution.len() == rank);

            for (i, v_i) in free_vars.iter().take(rank).zip(solution.iter()) {
                v[*i] = *v_i;
            }

            let mut rows = nalgebra::DVector::from_vec((0..m).collect());
            P.permute_rows(&mut rows);
            let rows = rows.rows(rank, rows.len() - rank);

            let b_tilde = &std_form.b - &std_form.A * &v;
            v = v.resize_vertically(n + m, 0.);

            std_form.A = std_form.A.resize_horizontally(n + rows.len(), 0.);
            let mut cur_col = std_form.A.ncols() - 1;

            for &i in rows.iter() {
                v[cur_col] = b_tilde[i].abs();
                std_form.A[(i, cur_col)] = b_tilde[i].signum();
                B.push(Basic::new(cur_col));
                cur_col -= 1;
            }
        } else {
            //TODO don't need phase 1 variables for equality constraints
            let b_tilde = &std_form.b - &std_form.A * &v;
            v = v.resize_vertically(n + m, 0.);
            std_form.A = std_form.A.resize_horizontally(n + m, 0.);

            for i in 0..m {
                let index = n + i;
                v[index] = b_tilde[i].abs();
                std_form.A[(i, index)] = b_tilde[i].signum();
                B.push(Basic::new(index));
            }
        }

        let mut phase_1_vars = Vec::with_capacity(m);

        for _ in 0..m {
            phase_1_vars.push(std_form.bounds.len());
            std_form.bounds.push(Bound::Lower(0.));
        }

        Some(PrimalPhase1 {
            phase_1_vars,
            std_form,
            point: PrimalFeasiblePoint(Point { x: v, N, B }),
        })
    }
}

impl std::convert::From<PrimalPhase1> for PrimalPhase2 {
    fn from(phase_1: PrimalPhase1) -> Self {
        debug!("converting phase 1 standard form to phase 2 standard form");

        let mut std_form = phase_1.std_form;

        for i in &phase_1.phase_1_vars {
            std_form.c[*i] = 0.;
            std_form.bounds[*i] = Bound::Fixed(0.);
        }

        for (i, var) in std_form.prob.variables.iter().enumerate() {
            std_form.c[i] = var.obj_coeff;

            // some of the free variable bounds may have been set to Fixed(0)
            std_form.bounds[i] = var.bound;
        }

        let mut point = phase_1.point;

        for var in &mut point.N {
            if matches!(std_form.bounds[var.index], Bound::Free) {
                var.bound = NonbasicBound::Free;
            }
        }

        PrimalPhase2 { std_form, point }
    }
}

impl std::convert::From<OptimalPoint> for PrimalFeasiblePoint {
    fn from(opt_pt: OptimalPoint) -> Self {
        Self(opt_pt.into_pt())
    }
}
