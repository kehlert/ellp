#![allow(non_snake_case)]

use crate::problem::{Bound, ConstraintOp, Problem};
use crate::util::EPS;

use log::debug;

#[derive(Debug, Clone)]
pub struct StandardForm<'a> {
    pub c: nalgebra::DVector<f64>,
    pub A: nalgebra::DMatrix<f64>,
    pub b: nalgebra::DVector<f64>,
    pub bounds: Vec<Bound>,
    prob: &'a Problem,
}

#[derive(Debug)]
pub struct StandardFormPhase1<'a> {
    pub std_form: StandardForm<'a>,
    phase_1_vars: Vec<usize>,
}

impl<'a> std::convert::From<&'a Problem> for StandardForm<'a> {
    fn from(prob: &'a Problem) -> StandardForm<'a> {
        debug!("converting problem to standard form");

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

        for (i, var) in prob.vars().iter().enumerate() {
            c[i] = var.obj_coeff;
            bounds[i] = var.bound;
        }

        let mut cur_slack_col = A.ncols().saturating_sub(1);

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

impl<'a> StandardForm<'a> {
    pub fn rows(&self) -> usize {
        self.A.nrows()
    }

    pub fn cols(&self) -> usize {
        self.A.ncols()
    }

    pub fn obj(&self, x: &nalgebra::DVector<f64>) -> f64 {
        self.c.dot(x)
    }

    pub fn extract_solution(&self, bfp: &BasicFeasiblePoint) -> nalgebra::DVector<f64> {
        let n = self.prob.vars().len();
        bfp.x.rows(0, n).into()
    }

    pub fn phase_1(mut self) -> (StandardFormPhase1<'a>, BasicFeasiblePoint) {
        debug!("converting standard form to phase 1 standard form");

        let n = self.cols();
        let m = self.rows();
        let mut N = Vec::with_capacity(n);
        let mut B = Vec::with_capacity(n.checked_sub(m).unwrap_or(10));
        let mut v = nalgebra::DVector::<f64>::zeros(self.cols());

        for (i, bound) in self.bounds.iter().enumerate() {
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

        for c_i in &mut self.c {
            *c_i = 0.;
        }

        self.c = self.c.resize_vertically(n + m, 1.);

        let free_vars: Vec<_> = self
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

        if !free_vars.is_empty() && !self.A.is_empty() {
            let cols: Vec<_> = free_vars.iter().map(|&i| self.A.column(i)).collect();
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
                self.bounds[i] = Bound::Fixed(0.);

                N.push(Nonbasic {
                    index: i,
                    bound: NonbasicBound::Lower,
                });
            }

            //don't need to use all the rows, but keeping it simpler for now
            let mut b_tilde = &self.b - &self.A * &v;
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

            let b_tilde = &self.b - &self.A * &v;
            v = v.resize_vertically(n + m, 0.);

            for i in 0..m {
                v[n + i] = b_tilde[i].abs();
            }

            self.A = self.A.resize_horizontally(n + rows.len(), 0.);
            let mut cur_col = self.A.ncols() - 1;

            for &i in rows.iter() {
                self.A[(i, cur_col)] = b_tilde[i].signum();
                B.push(Basic::new(cur_col));
                cur_col -= 1;
            }
        } else {
            //TODO don't need phase 1 variables for equality constraints
            let b_tilde = &self.b - &self.A * &v;
            v = v.resize_vertically(n + m, 0.);
            self.A = self.A.resize_horizontally(n + m, 0.);

            for i in 0..m {
                let index = n + i;
                v[index] = b_tilde[i].abs();
                self.A[(i, index)] = b_tilde[i].signum();
                B.push(Basic::new(index));
            }
        }

        let mut phase_1_vars = Vec::with_capacity(m);

        for _ in 0..m {
            phase_1_vars.push(self.bounds.len());
            self.bounds.push(Bound::Lower(0.));
        }

        (
            StandardFormPhase1 {
                std_form: self,
                phase_1_vars,
            },
            BasicFeasiblePoint { x: v, N, B },
        )
    }
}

impl<'a> StandardFormPhase1<'a> {
    pub fn obj(&self, x: &nalgebra::DVector<f64>) -> f64 {
        self.std_form.obj(x)
    }

    pub fn phase_2(self, bfp: &mut BasicFeasiblePoint) -> StandardForm<'a> {
        debug!("converting phase 1 standard form to phase 2 standard form");

        let mut std_form = self.std_form;

        for i in &self.phase_1_vars {
            std_form.bounds[*i] = Bound::Fixed(0.);
        }

        for (i, var) in std_form.prob.vars().iter().enumerate() {
            std_form.c[i] = var.obj_coeff;

            // some of the free variable bounds may have been set to Fixed(0)
            std_form.bounds[i] = var.bound;
        }

        //needs to be after the above code, because it relies on std_form.bounds
        for var in &mut bfp.N {
            if matches!(std_form.bounds[var.index], Bound::Free) {
                var.bound = NonbasicBound::Free;
            }
        }

        std_form
    }
}

#[derive(Debug, Clone)]
pub struct BasicFeasiblePoint {
    pub x: nalgebra::DVector<f64>,
    pub N: Vec<Nonbasic>,
    pub B: Vec<Basic>,
}

impl BasicFeasiblePoint {
    pub fn x(&self) -> &nalgebra::DVector<f64> {
        &self.x
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
