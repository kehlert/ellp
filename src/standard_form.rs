#![allow(non_snake_case)]

use crate::problem::{Bound, ConstraintOp, Problem};
use log::debug;
use std::collections::HashMap;

pub trait BasicPoint: std::ops::Deref<Target = Point> + std::ops::DerefMut<Target = Point> {
    fn into_pt(self) -> Point;
}

#[derive(Debug, Clone)]
pub struct Point {
    pub x: nalgebra::DVector<f64>,
    pub N: Vec<Nonbasic>,
    pub B: Vec<Basic>,
}

#[derive(Debug, Clone)]
pub struct StandardForm {
    pub c: nalgebra::DVector<f64>,
    pub A: nalgebra::DMatrix<f64>,
    pub b: nalgebra::DVector<f64>,
    pub bounds: Vec<Bound>,
    pub prob: Problem,
}

impl StandardForm {
    #[inline]
    pub fn rows(&self) -> usize {
        self.A.nrows()
    }

    #[inline]
    pub fn cols(&self) -> usize {
        self.A.ncols()
    }

    #[inline]
    pub fn obj(&self, x: &nalgebra::DVector<f64>) -> f64 {
        self.c.dot(x)
    }

    pub fn dual_obj(&self, y: &nalgebra::DVector<f64>, d: &nalgebra::DVector<f64>) -> f64 {
        assert!(d.len() == self.bounds.len());

        let mut obj = self.b.dot(y);

        for (i, bound) in self.bounds.iter().enumerate() {
            match bound {
                Bound::Free => (),
                Bound::Lower(lb) => obj += lb * d[i],
                Bound::Upper(ub) => obj += ub * d[i],
                Bound::TwoSided(lb, ub) => obj += if d[i] > 0. { lb * d[i] } else { ub * d[i] },
                Bound::Fixed(val) => obj += val * d[i],
            }
        }

        obj
    }

    #[inline]
    pub fn extract_solution<'a>(&self, point: &'a Point) -> nalgebra::DVectorSlice<'a, f64> {
        let n = self.prob.variables.len();
        point.x.rows(0, n)
    }
}

impl std::convert::From<Problem> for StandardForm {
    fn from(prob: Problem) -> StandardForm {
        debug!("converting problem to standard form");

        let n = prob.variables.len();
        let m = prob.constraints().len();

        let num_slack_vars = prob
            .constraints()
            .iter()
            .map(|constraint| match constraint.op {
                ConstraintOp::Lte | ConstraintOp::Gte => 1,
                ConstraintOp::Eq => 0,
            })
            .sum::<usize>();

        let total_vars = n + num_slack_vars;

        let mut c = nalgebra::DVector::zeros(total_vars);
        let mut A = nalgebra::DMatrix::zeros(m, total_vars);
        let mut b = nalgebra::DVector::zeros(A.nrows());

        //default to Lower(0.), because that's what the slack variable bounds are
        let mut bounds = vec![Bound::Lower(0.); total_vars];
        let mut id_to_index = HashMap::with_capacity(prob.variables.len());

        for (i, var) in prob.variables.iter().enumerate() {
            c[i] = var.obj_coeff;
            bounds[i] = var.bound;
            id_to_index.insert(var.id, i);
        }

        let mut cur_slack_col = A.ncols().saturating_sub(1);

        for (i, constraint) in prob.constraints().iter().enumerate() {
            b[i] = constraint.rhs;

            for (id, coeff) in &constraint.coeffs {
                let var_index = *id_to_index.get(id).unwrap();
                A[(i, var_index)] = *coeff;
            }

            if let Some(slack_coeff) = match constraint.op {
                ConstraintOp::Lte => Some(1.),
                ConstraintOp::Eq => None,
                ConstraintOp::Gte => Some(-1.),
            } {
                A[(i, cur_slack_col)] = slack_coeff;
                cur_slack_col -= 1;
            }
        }

        assert!(cur_slack_col == n.saturating_sub(1));

        StandardForm {
            c,
            A,
            b,
            bounds,
            prob,
        }
    }
}

#[derive(Debug, Clone)]
pub struct Nonbasic {
    pub index: usize,
    pub bound: NonbasicBound,
}

impl Nonbasic {
    pub fn new(index: usize, bound: NonbasicBound) -> Self {
        Self { index, bound }
    }
}

#[derive(Debug, Clone, Copy)]
pub enum NonbasicBound {
    Lower,
    Upper,
    Free,
}

#[derive(Debug, Clone, Copy)]
pub struct Basic {
    pub index: usize,
}

impl Basic {
    pub fn new(index: usize) -> Self {
        Self { index }
    }
}

pub trait Phase2 {}
