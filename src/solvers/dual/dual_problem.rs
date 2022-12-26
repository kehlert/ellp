#![allow(non_snake_case)]

use crate::problem::{Bound, ConstraintOp, Problem};
use crate::standard_form::{
    Basic, BasicPoint, Nonbasic, NonbasicBound, Point, StandardForm, StandardizedProblem,
};
use crate::util::EPS;

use std::collections::HashSet;

#[derive(Debug, Clone)]
pub struct DualFeasiblePoint {
    pub y: nalgebra::DVector<f64>,
    pub d: nalgebra::DVector<f64>,
    pub point: Point,
}

impl BasicPoint for DualFeasiblePoint {
    #[inline]
    fn into_pt(self) -> Point {
        self.point
    }
}

impl std::ops::Deref for DualFeasiblePoint {
    type Target = Point;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.point
    }
}

impl std::ops::DerefMut for DualFeasiblePoint {
    #[inline]
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.point
    }
}

#[derive(Debug)]
pub struct DualPhase1 {
    pub std_form: StandardForm,
    pub point: DualFeasiblePoint,
    orig_std_form: StandardForm,
}

impl DualPhase1 {
    #[inline]
    pub fn into_orig_prob(self) -> Problem {
        self.orig_std_form.prob
    }
}

impl StandardizedProblem for DualPhase1 {
    type FeasiblePoint = DualFeasiblePoint;

    #[inline]
    fn obj(&self) -> f64 {
        self.std_form.dual_obj(&self.point.y, &self.point.d)
    }

    #[inline]
    fn unpack(&mut self) -> (&StandardForm, &mut Self::FeasiblePoint) {
        (&self.std_form, &mut self.point)
    }
}

#[derive(Debug)]
pub struct DualPhase2 {
    pub std_form: StandardForm,
    pub point: DualFeasiblePoint,
}

impl StandardizedProblem for DualPhase2 {
    type FeasiblePoint = DualFeasiblePoint;

    #[inline]
    fn obj(&self) -> f64 {
        self.std_form.dual_obj(&self.point.y, &self.point.d)
    }

    #[inline]
    fn unpack(&mut self) -> (&StandardForm, &mut Self::FeasiblePoint) {
        (&self.std_form, &mut self.point)
    }
}

impl std::convert::From<Problem> for Option<DualPhase1> {
    fn from(prob: Problem) -> Self {
        let orig_std_form: StandardForm = match prob.into() {
            Some(std_form) => std_form,
            None => return None,
        };

        let mut phase_1_prob = Problem::new();
        let mut vars_kept = HashSet::new();

        for i in 0..orig_std_form.cols() {
            if let Some(bound) = match &orig_std_form.bounds[i] {
                Bound::Free => Some(Bound::TwoSided(-1., 1.)),
                Bound::Lower(..) => Some(Bound::TwoSided(0., 1.)),
                Bound::Upper(..) => Some(Bound::TwoSided(-1., 0.)),
                Bound::TwoSided(..) => None,
                Bound::Fixed(..) => None,
            } {
                vars_kept.insert(i);
                phase_1_prob
                    .add_var_with_id(orig_std_form.c[i], bound, i.into(), None)
                    .unwrap();
            }
        }

        for i in 0..orig_std_form.rows() {
            let coeffs: Vec<_> = orig_std_form
                .A
                .row(i)
                .iter()
                .enumerate()
                .filter_map(|(col_index, &coeff)| {
                    if vars_kept.contains(&col_index) {
                        Some((col_index.into(), coeff))
                    } else {
                        None
                    }
                })
                .collect();

            if !coeffs.is_empty() {
                phase_1_prob
                    .add_constraint(coeffs, ConstraintOp::Eq, 0.)
                    .unwrap();
            }
        }

        let std_form: StandardForm = match phase_1_prob.into() {
            Some(std_form) => std_form,
            None => return None,
        };

        let lu = std_form.A.transpose().lu();

        for u_i in lu.u().diagonal().iter() {
            if u_i.abs() < EPS {
                panic!("should always have a basis available");
            }
        }

        //would be good to avoid creating this entire vector
        let n = std_form.A.ncols();
        let mut perm_cols = nalgebra::DVector::<usize>::from_iterator(n, 0..n);
        lu.p().permute_rows(&mut perm_cols);

        let B: Vec<_> = (0..std_form.A.nrows())
            .map(|i| Basic::new(perm_cols[i]))
            .collect();

        let mut N: Vec<_> = (std_form.A.nrows()..perm_cols.len())
            .map(|i| Nonbasic::new(perm_cols[i], NonbasicBound::Lower))
            .collect();

        let c_B = nalgebra::DVector::from_iterator(B.len(), B.iter().map(|i| std_form.c[i.index]));
        let A_B = std_form.A.select_columns(B.iter().map(|b| &b.index));
        let A_B_lu = A_B.lu();

        if !B.is_empty() {
            //would be good to avoid creating this vector
            let y_tilde = A_B_lu.u().tr_solve_upper_triangular(&c_B).unwrap();
            let mut y = A_B_lu.l().tr_solve_lower_triangular(&y_tilde).unwrap();
            A_B_lu.p().inv_permute_rows(&mut y);

            let d = &std_form.c - std_form.A.tr_mul(&y);

            let mut x = nalgebra::DVector::zeros(std_form.bounds.len());

            assert!(d.len() == std_form.bounds.len());

            for n in &mut N {
                let i = n.index;

                match std_form.bounds[i] {
                    Bound::TwoSided(lb, ub) => {
                        if d[i] >= 0. {
                            x[i] = lb;
                            n.bound = NonbasicBound::Lower;
                        } else {
                            x[i] = ub;
                            n.bound = NonbasicBound::Upper;
                        }
                    }

                    Bound::Fixed(val) => {
                        x[i] = val;

                        if d[i] >= 0. {
                            n.bound = NonbasicBound::Lower;
                        } else {
                            n.bound = NonbasicBound::Upper;
                        }
                    }

                    _ => panic!("bounds should always be fixed or two-sided"),
                }
            }

            //assumes that, at this point, the elements of x correspond to x_B are 0
            let b_tilde = &std_form.b - &std_form.A * &x;
            let x_B = A_B_lu.solve(&b_tilde).unwrap();

            assert!(x_B.len() == B.len());

            for (i, val) in B.iter().zip(x_B.iter()) {
                x[i.index] = *val;
            }

            let point = DualFeasiblePoint {
                y,
                d,
                point: Point { x, N, B },
            };

            Some(DualPhase1 {
                std_form,
                point,
                orig_std_form,
            })
        } else {
            let empty_vec = nalgebra::DVector::zeros(0);
            let mut x = nalgebra::DVector::zeros(N.len());
            assert_eq!(N.len(), std_form.bounds.len());

            for n in &mut N {
                let i = n.index;
                n.bound = NonbasicBound::Lower;

                match std_form.bounds[i] {
                    Bound::TwoSided(lb, _ub) => x[i] = lb,
                    Bound::Fixed(val) => x[i] = val,
                    _ => panic!("bounds should always be fixed or two-sided"),
                }
            }

            let point = DualFeasiblePoint {
                y: empty_vec,
                d: std_form.c.clone(),
                point: Point { x, N, B },
            };

            Some(DualPhase1 {
                std_form,
                point,
                orig_std_form,
            })
        }
    }
}

impl std::convert::From<DualPhase1> for DualPhase2 {
    fn from(phase_1: DualPhase1) -> Self {
        let phase_1_prob = phase_1.std_form.prob;
        let std_form = phase_1.orig_std_form;
        let mut is_basic = vec![false; std_form.cols()];

        let B: Vec<_> = phase_1
            .point
            .B
            .iter()
            .map(|b| {
                let index = phase_1_prob.variables[b.index].id.into();
                is_basic[index] = true;
                Basic::new(index)
            })
            .collect();

        if !B.is_empty() {
            let c_B = std_form.c.select_rows(B.iter().map(|b| &b.index));
            let A_B_lu = std_form.A.select_columns(B.iter().map(|b| &b.index)).lu();

            let y_tilde = A_B_lu.u().tr_solve_upper_triangular(&c_B).unwrap();
            let mut y = A_B_lu.l().tr_solve_lower_triangular(&y_tilde).unwrap();
            A_B_lu.p().inv_permute_rows(&mut y);

            let d = &std_form.c - std_form.A.tr_mul(&y);

            let (x_N, N): (Vec<_>, Vec<_>) = is_basic
                .iter()
                .enumerate()
                .filter_map(|(i, b)| {
                    if !b {
                        let d_i = d[i];

                        let (x_i, bound) = match &std_form.bounds[i] {
                            Bound::Free => {
                                assert!(d_i.abs() < EPS);
                                (0., NonbasicBound::Free)
                            }

                            Bound::Lower(lb) => {
                                assert!(d_i > -EPS);
                                (*lb, NonbasicBound::Lower)
                            }

                            Bound::Upper(ub) => {
                                assert!(d_i < EPS);
                                (*ub, NonbasicBound::Upper)
                            }

                            Bound::TwoSided(lb, ub) => {
                                if d_i >= 0. {
                                    (*lb, NonbasicBound::Lower)
                                } else {
                                    (*ub, NonbasicBound::Upper)
                                }
                            }

                            Bound::Fixed(val) => (*val, NonbasicBound::Lower),
                        };

                        Some((x_i, Nonbasic::new(i, bound)))
                    } else {
                        None
                    }
                })
                .unzip();

            let x_N = nalgebra::DVector::from_vec(x_N);
            let A_N = std_form.A.select_columns(N.iter().map(|b| &b.index));
            let x_B = A_B_lu.solve(&(&std_form.b - &A_N * &x_N)).unwrap();

            let mut x = nalgebra::DVector::zeros(std_form.A.ncols());

            assert!(B.len() == x_B.len());

            for (b, val) in B.iter().zip(x_B.iter()) {
                x[b.index] = *val;
            }

            assert!(N.len() == x_N.len());

            for (n, val) in N.iter().zip(x_N.iter()) {
                x[n.index] = *val;
            }

            let point = DualFeasiblePoint {
                y,
                d,
                point: Point { x, N, B },
            };

            DualPhase2 { std_form, point }
        } else {
            let empty_vec = nalgebra::DVector::zeros(0);
            let number_nonbasic = std_form.A.ncols();
            let mut x_N = nalgebra::DVector::zeros(number_nonbasic);

            let N: Vec<_> = is_basic
                .iter()
                .enumerate()
                .filter_map(|(i, b)| {
                    if !b {
                        let bound = match &std_form.bounds[i] {
                            Bound::Free => {
                                x_N[i] = 0.;
                                NonbasicBound::Free
                            }

                            Bound::Lower(lb) => {
                                x_N[i] = *lb;
                                NonbasicBound::Lower
                            }

                            Bound::Upper(ub) => {
                                x_N[i] = *ub;
                                NonbasicBound::Upper
                            }

                            Bound::TwoSided(lb, _ub) => {
                                x_N[i] = *lb;
                                NonbasicBound::Lower
                            }

                            Bound::Fixed(val) => {
                                x_N[i] = *val;
                                NonbasicBound::Lower
                            }
                        };

                        Some(Nonbasic::new(i, bound))
                    } else {
                        None
                    }
                })
                .collect();

            let point = DualFeasiblePoint {
                y: empty_vec,
                d: std_form.c.clone(),
                point: Point { x: x_N, N, B },
            };

            DualPhase2 { std_form, point }
        }
    }
}
