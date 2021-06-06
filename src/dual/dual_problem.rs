#![allow(non_snake_case)]

use crate::problem::{Bound, ConstraintOp, Problem};
use crate::standard_form::{Basic, BasicPoint, Nonbasic, NonbasicBound, Point, StandardForm};
use crate::util::EPS;

use std::collections::HashMap;

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
        let mut phase_1_prob = Problem::new();
        let mut old_to_new_var_ids = HashMap::new();
        // let mut std_form: StandardForm = prob.into();

        for var in prob.variables {
            if let Some(bound) = match &var.bound {
                Bound::Free => Some(Bound::TwoSided(-1., 1.)),
                Bound::Lower(..) => Some(Bound::TwoSided(0., 1.)),
                Bound::Upper(..) => Some(Bound::TwoSided(-1., 0.)),
                Bound::TwoSided(..) => None,
                Bound::Fixed(..) => None,
            } {
                let new_var_id = phase_1_prob
                    .add_var(var.obj_coeff, bound, var.name)
                    .unwrap();

                old_to_new_var_ids.insert(var.id, new_var_id);
            }
        }

        println!("kept var ids: {:?}", old_to_new_var_ids);

        for constraint in prob.constraints {
            println!("{:?}", constraint);

            let mut coeffs: Vec<_> = constraint
                .coeffs
                .into_iter()
                .filter_map(|(var_id, coeff)| old_to_new_var_ids.get(&var_id).map(|c| (*c, coeff)))
                .collect();

            println!("{:?}\n", coeffs);

            if let Some(slack_coeff) = match constraint.op {
                ConstraintOp::Lte => Some(1.),
                ConstraintOp::Eq => None,
                ConstraintOp::Gte => Some(-1.),
            } {
                //slack vars always have a bound of [0, infty)
                let slack_var = phase_1_prob
                    .add_var(0., Bound::TwoSided(0., 1.), None)
                    .unwrap();

                coeffs.push((slack_var, slack_coeff));
            }

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

        let N_indices: Vec<_> = (std_form.A.nrows()..perm_cols.len())
            .map(|i| perm_cols[i])
            .collect();

        //would be good to avoid creating this vector
        let c_B = nalgebra::DVector::from_iterator(B.len(), B.iter().map(|i| std_form.c[i.index]));
        let A_B_cols: Vec<_> = B.iter().map(|i| std_form.A.column(i.index)).collect();
        let A_B = nalgebra::DMatrix::from_columns(&A_B_cols);

        println!("c:{}", std_form.c);
        println!("A:{}", std_form.A);
        println!("A_B:{}", A_B);
        println!("B:{:?}", B);

        let y = A_B.transpose().lu().solve(&c_B).unwrap();
        let d = &std_form.c - std_form.A.tr_mul(&y);

        println!("y: {}", y);
        println!("d: {}", d);

        let mut x = nalgebra::DVector::zeros(std_form.bounds.len());

        assert!(d.len() == std_form.bounds.len());

        let mut N = Vec::with_capacity(N_indices.len());

        for &i in &N_indices {
            if let Bound::TwoSided(lb, ub) = std_form.bounds[i] {
                if d[i] >= 0. {
                    x[i] = lb;
                    N.push(Nonbasic::new(i, NonbasicBound::Lower))
                } else {
                    x[i] = ub;
                    N.push(Nonbasic::new(i, NonbasicBound::Upper))
                }
            } else {
                panic!("bounds should always be two-sided");
            }
        }

        println!("N:{:?}", N);

        //assumes that, at this point, the elements of x correspond to x_B are 0
        let b_tilde = &std_form.b - &std_form.A * &x;
        let x_B = A_B.lu().solve(&b_tilde).unwrap();

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

        DualPhase1 { std_form, point }
    }
}

impl std::convert::From<DualPhase1> for DualPhase2 {
    fn from(mut _phase_1: DualPhase1) -> Self {
        todo!()
    }
}
