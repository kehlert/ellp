#![allow(non_snake_case)]

use crate::problem::{Bound, Constraint, ConstraintOp, Problem, Variable, VariableId};
use std::borrow::Borrow;

use crate::util::EPS;

#[derive(Debug, Clone)]
pub struct StandardForm {
    pub c: nalgebra::DVector<f64>,
    pub A: nalgebra::DMatrix<f64>,
    pub b: nalgebra::DVector<f64>,
    pub bounds: Vec<Bound>,

    //store the problem so we know that it is not mutated while using its standard form
    prob: Problem,
}

impl std::convert::From<Problem> for StandardForm {
    fn from(prob: Problem) -> StandardForm {
        let n = prob.vars().len();
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

        let mut cur_slack_col = A.ncols() - 1;

        for (i, var) in prob.vars().iter().enumerate() {
            c[i] = var.obj_coeff;
            bounds[i] = var.bound;
        }

        for (i, constraint) in prob.constraints().iter().enumerate() {
            b[i] = constraint.rhs;

            for (id, coeff) in &constraint.coeffs {
                A[(i, id.into())] = *coeff;
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

        assert!(cur_slack_col == n - 1);

        StandardForm {
            c,
            A,
            b,
            bounds,
            prob,
        }
    }
}

impl StandardForm {
    pub fn rows(&self) -> usize {
        self.A.nrows()
    }

    pub fn cols(&self) -> usize {
        self.A.ncols()
    }

    //returns a feasible point
    pub fn phase_1(prob: Problem) -> (StandardForm, InitialPoint) {
        let mut std_form: StandardForm = prob.into();
        let n = std_form.cols();
        let m = std_form.rows();
        let mut N = Vec::with_capacity(m);
        let mut B = Vec::with_capacity(n.checked_sub(m).unwrap_or(10));

        let mut v = nalgebra::DVector::<f64>::zeros(std_form.cols());

        for (i, bound) in std_form.bounds.iter().enumerate() {
            match *bound {
                Bound::Free => (), //will set these values later

                Bound::Lower(lb) => {
                    v[i] = lb;
                    N.push(i);
                }

                Bound::Upper(ub) => {
                    v[i] = ub;
                    N.push(i);
                }

                Bound::TwoSided(lb, _ub) => {
                    v[i] = lb;
                    N.push(i);
                }

                Bound::Fixed(fixed_val) => {
                    v[i] = fixed_val;
                    N.push(i);
                }
            }
        }

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
            B.push(i);
        }

        for &i in free_vars.iter().skip(rank) {
            std_form.bounds[i] = Bound::Fixed(0.);
            N.push(i);
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
        let remaining = b_tilde.len().checked_sub(rank).unwrap();
        let b_tilde = b_tilde.rows(0, remaining);

        std_form.A = std_form.A.resize_horizontally(n + rows.len(), 0.);
        v = v.resize_vertically(n + rows.len(), 0.);
        let mut cur_col = std_form.A.ncols() - 1;

        for &i in rows.iter() {
            v[cur_col] = b_tilde[i].abs();
            std_form.A[(i, cur_col)] = b_tilde[i].signum();
            B.push(cur_col);
            cur_col -= 1;
        }

        for c_i in &mut std_form.c {
            *c_i = 0.;
        }

        std_form.c = std_form.c.resize_vertically(n + rows.len(), 1.);

        for _ in &rows {
            std_form.bounds.push(Bound::Lower(0.));
        }

        (std_form, InitialPoint { x: v, N, B })
    }
}

#[derive(Debug, Clone)]
pub struct InitialPoint {
    pub x: nalgebra::DVector<f64>,
    pub N: Vec<usize>,
    pub B: Vec<usize>,
}
