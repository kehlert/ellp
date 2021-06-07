#![allow(non_snake_case)]

use crate::problem::{Bound, ConstraintOp, Problem, VariableId};
use crate::standard_form::{Basic, BasicPoint, Nonbasic, NonbasicBound, Point, StandardForm};
use crate::util::EPS;

use std::collections::HashSet;

pub trait DualProblem {
    fn obj(&self) -> f64;
    fn std_form(&self) -> &StandardForm;
    fn pt(&self) -> &DualFeasiblePoint;
    fn pt_mut(&mut self) -> &mut DualFeasiblePoint;
    fn unpack(&mut self) -> (&StandardForm, &mut DualFeasiblePoint);
}

#[derive(Debug, Clone)]
pub struct DualFeasiblePoint {
    pub y: nalgebra::DVector<f64>,
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

impl DualProblem for DualPhase1 {
    fn obj(&self) -> f64 {
        self.std_form.obj(&self.point.x)
    }

    #[inline]
    fn std_form(&self) -> &StandardForm {
        &self.std_form
    }

    #[inline]
    fn pt(&self) -> &DualFeasiblePoint {
        &self.point
    }

    #[inline]
    fn pt_mut(&mut self) -> &mut DualFeasiblePoint {
        &mut self.point
    }

    #[inline]
    fn unpack(&mut self) -> (&StandardForm, &mut DualFeasiblePoint) {
        (&self.std_form, &mut self.point)
    }
}

#[derive(Debug)]
pub struct DualPhase2 {
    pub std_form: StandardForm,
    pub point: DualFeasiblePoint,
}

impl DualProblem for DualPhase2 {
    fn obj(&self) -> f64 {
        self.std_form.obj(&self.point.x)
    }

    #[inline]
    fn std_form(&self) -> &StandardForm {
        &self.std_form
    }

    #[inline]
    fn pt(&self) -> &DualFeasiblePoint {
        &self.point
    }

    #[inline]
    fn pt_mut(&mut self) -> &mut DualFeasiblePoint {
        &mut self.point
    }

    #[inline]
    fn unpack(&mut self) -> (&StandardForm, &mut DualFeasiblePoint) {
        (&self.std_form, &mut self.point)
    }
}

impl std::convert::From<Problem> for DualPhase1 {
    fn from(prob: Problem) -> Self {
        println!("orig prob:\n{}", prob);

        let orig_std_form: StandardForm = prob.into();
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

            phase_1_prob
                .add_constraint(coeffs, ConstraintOp::Eq, 0.)
                .unwrap();
        }

        println!("{}", phase_1_prob);

        let std_form: StandardForm = phase_1_prob.into();
        let lu = std_form.A.transpose().lu();

        for u_i in lu.u().diagonal().iter() {
            if u_i.abs() < EPS {
                todo!("repair the basis");
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

        //would be good to avoid creating this vector
        let c_B = nalgebra::DVector::from_iterator(B.len(), B.iter().map(|i| std_form.c[i.index]));
        let A_B = std_form.A.select_columns(B.iter().map(|b| &b.index));

        println!("c:{}", std_form.c);
        println!("A:{}", std_form.A);
        println!("A_B:{}", A_B);
        println!("B:{:?}", B);

        let A_B_lu = A_B.lu();
        let y_tilde = A_B_lu.u().tr_solve_upper_triangular(&c_B).unwrap();
        let mut y = A_B_lu.l().tr_solve_lower_triangular(&y_tilde).unwrap();
        A_B_lu.p().inv_permute_rows(&mut y);

        let d = &std_form.c - std_form.A.tr_mul(&y);

        println!("y: {}", y);
        println!("d: {}", d);

        let mut x = nalgebra::DVector::zeros(std_form.bounds.len());

        assert!(d.len() == std_form.bounds.len());

        for n in &mut N {
            let i = n.index;

            if let Bound::TwoSided(lb, ub) = std_form.bounds[i] {
                if d[i] >= 0. {
                    x[i] = lb;
                    n.bound = NonbasicBound::Lower;
                } else {
                    x[i] = ub;
                    n.bound = NonbasicBound::Upper;
                }
            } else {
                panic!("bounds should always be two-sided");
            }
        }

        println!("N:{:?}", N);

        //assumes that, at this point, the elements of x correspond to x_B are 0
        let b_tilde = &std_form.b - &std_form.A * &x;
        let x_B = A_B_lu.solve(&b_tilde).unwrap();

        assert!(x_B.len() == B.len());

        for (i, val) in B.iter().zip(x_B.iter()) {
            x[i.index] = *val;
        }

        println!("{:?}", std_form.bounds);
        println!("y: {}", y);
        println!("x: {}", x);
        println!("d: {}", d);
        println!("Ax:{}", &std_form.A * &x);
        println!("b:{}", std_form.b);

        let A_B_cols: Vec<_> = B.iter().map(|i| std_form.A.column(i.index)).collect();
        let A_B = nalgebra::DMatrix::from_columns(&A_B_cols);
        println!("A_B_T y:{}", A_B.tr_mul(&y));
        println!("c_B: {}", c_B);

        let point = DualFeasiblePoint {
            y,
            point: Point { x, N, B },
        };

        DualPhase1 {
            std_form,
            point,
            orig_std_form,
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

        let d = &std_form.c - std_form.A.tr_mul(&phase_1.point.y);
        println!("d: {}", d);

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

        println!("B: {:?}", B);
        println!("N: {:?}", N);

        let c_B = std_form.c.select_rows(B.iter().map(|b| &b.index));
        let A_B_lu = std_form.A.select_columns(B.iter().map(|b| &b.index)).lu();
        let A_N = std_form.A.select_columns(N.iter().map(|b| &b.index));

        let x_B = A_B_lu.solve(&(&std_form.b - &A_N * &x_N)).unwrap();

        let mut x = nalgebra::DVector::zeros(std_form.A.ncols());

        for (b, val) in B.iter().zip(x_B.iter()) {
            x[b.index] = *val;
        }

        for (n, val) in N.iter().zip(x_N.iter()) {
            x[n.index] = *val;
        }

        println!("x:{}", x);
        println!("obj:{}", std_form.obj(&x));

        let y_tilde = A_B_lu.u().tr_solve_upper_triangular(&c_B).unwrap();
        let mut y = A_B_lu.l().tr_solve_lower_triangular(&y_tilde).unwrap();
        A_B_lu.p().inv_permute_rows(&mut y);

        println!("y:{}", y);

        let point = DualFeasiblePoint {
            y,
            point: Point { x, N, B },
        };

        DualPhase2 { point, std_form }
    }
}
