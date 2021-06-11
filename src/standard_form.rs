#![allow(non_snake_case)]

use crate::problem::{Bound, ConstraintOp, Problem};
use crate::util::EPS;
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

//returns None if infeasible
impl std::convert::From<Problem> for Option<StandardForm> {
    fn from(prob: Problem) -> Self {
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

            if constraint.coeffs.is_empty() && b[i] != 0. {
                return None;
            }

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

        //remove redundant rows
        let A_qr = A.transpose().col_piv_qr();
        let mut R = A_qr.r();
        A_qr.p().inv_permute_rows(&mut b);

        //sometimes a diagonal element is nonzero when it should actually be zero
        //can lead to issues when detecting infeasibility a couple lines below
        for i in 0..R.nrows() {
            let r_ii = &mut R[(i, i)];

            if r_ii.abs() < EPS {
                *r_ii = 0.
            }
        }

        //need to handle R = [0] as a special case
        //we know that R.ncols() >= R.nrows()

        let is_trivial = !R.is_empty() && R[(0, 0)].abs() < EPS && b[0].abs() < EPS;

        if !is_trivial {
            R.tr_solve_upper_triangular(&b)?;
        }

        A_qr.p().permute_rows(&mut b);

        //now we know that the system is feasible, and we can remove redundant rows of [A, b]

        //R.nrows() = min(A rows, A cols), which is the max rank of A
        let num_indep_rows = (0..R.nrows())
            .map(|i| (i, R[(i, i)]))
            .find(|(_i, r)| r.abs() < EPS)
            .map(|(i, _r)| i)
            .unwrap_or_else(|| R.nrows());

        let mut indep_rows = nalgebra::DVector::from_iterator(A.nrows(), 0..A.nrows());
        A_qr.p().permute_rows(&mut indep_rows);
        let indep_rows = indep_rows.rows(0, num_indep_rows);

        let A = A.select_rows(indep_rows.as_slice());
        let b = b.select_rows(indep_rows.as_slice());

        Some(StandardForm {
            c,
            A,
            b,
            bounds,
            prob,
        })
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

impl std::fmt::Display for StandardForm {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        writeln!(f, "c:{}", self.c)?;
        writeln!(f, "A:{}", self.A)?;
        writeln!(f, "b:{}", self.b)?;

        writeln!(f, "bounds:\n")?;

        for (i, bound) in self.bounds.iter().enumerate() {
            writeln!(f, "x{}: {}", i, bound)?;
        }

        Ok(())
    }
}
