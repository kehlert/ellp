#![allow(non_snake_case)]

use ndarray::linalg::Dot;

use crate::problem::{Bound, Constraint, ConstraintOp, Problem, Variable, VariableId};
use std::{borrow::Borrow, collections::HashMap};

type VarMap = HashMap<VariableId, VarParts>;

#[derive(Debug, Clone)]
enum VarParts {
    Pos(usize),
    Neg(usize),
    PosNeg(usize, usize),
}

#[derive(Debug, Clone)]
pub struct StandardForm {
    //note: first n variables will always corresponds to the original problem variables
    //and likewise, the first m rows are for the original constraints
    pub A: sprs::CsMat<f64>,
    pub b: sprs::CsVec<f64>,
    pub c: sprs::CsVec<f64>,

    //map from original variable to the columns in A corresponding to the positive
    //and negative parts of the variable
    var_map: VarMap,
}

struct StandardFormBuilder {
    c_ind: Vec<usize>,
    c_data: Vec<f64>,
    A: sprs::TriMat<f64>,
    b_ind: Vec<usize>,
    b_data: Vec<f64>,
    var_map: VarMap,
    cur_row: usize,
    cur_slack_col: usize,
}

impl StandardFormBuilder {
    fn new(prob: &Problem) -> Self {
        let n = prob.vars().len();
        let m = prob.constraints().len();

        let (var_map, num_bounds_constraints, mut total_vars) = Self::split_vars(prob);

        let num_coeffs: usize = prob
            .vars()
            .iter()
            .filter_map(|var| {
                if var.obj_coeff != 0. {
                    Some(match *var_map.get(&var.id).unwrap() {
                        VarParts::Pos(..) => 1,
                        VarParts::Neg(..) => 1,
                        VarParts::PosNeg(..) => 2,
                    })
                } else {
                    None
                }
            })
            .sum();

        total_vars += prob
            .constraints()
            .iter()
            .map(|constraint| match constraint.op {
                ConstraintOp::Lte => 1,

                ConstraintOp::Eq => 0,

                ConstraintOp::Gte => 1,
            })
            .sum::<usize>();

        let A = sprs::TriMat::new((m + num_bounds_constraints, total_vars));

        let b = vec![0.; A.rows()];
        let cur_slack_col = A.cols() - 1;

        println!(
            "n: {}, num_bounds_constraints: {}, total_vars: {}",
            n, num_bounds_constraints, total_vars
        );

        println!("A cols: {}", A.cols());

        StandardFormBuilder {
            c_ind: Vec::with_capacity(num_coeffs),
            c_data: Vec::with_capacity(num_coeffs),
            A,
            b_ind: Vec::with_capacity(num_coeffs),
            b_data: Vec::with_capacity(num_coeffs),
            var_map,
            cur_row: 0,
            cur_slack_col,
        }
    }

    fn build(self) -> StandardForm {
        assert!(self.cur_row == self.A.rows());

        StandardForm {
            A: self.A.to_csc(),
            b: sprs::CsVec::new_from_unsorted(self.A.rows(), self.b_ind, self.b_data).unwrap(),
            c: sprs::CsVec::new_from_unsorted(self.A.cols(), self.c_ind, self.c_data).unwrap(),
            var_map: self.var_map,
        }
    }

    fn add_obj_coeffs(&mut self, var: &Variable) {
        println!("{:?}", self.var_map.get(&var.id).unwrap());
        match *self.var_map.get(&var.id).unwrap() {
            VarParts::Pos(i) => {
                self.c_ind.push(i);
                self.c_data.push(var.obj_coeff);
            }

            VarParts::Neg(i) => {
                self.c_ind.push(i);
                self.c_data.push(-var.obj_coeff);
            }

            VarParts::PosNeg(i1, i2) => {
                self.c_ind.push(i1);
                self.c_data.push(var.obj_coeff);

                self.c_ind.push(i2);
                self.c_data.push(-var.obj_coeff);
            }
        }
    }

    fn add_bound(&mut self, var: &Variable, rhs: f64, lower: bool) {
        if rhs == 0. {
            return;
        }

        let var_parts = self.var_map.get(&var.id).unwrap();
        self.b_ind.push(self.cur_row);
        self.b_data.push(rhs);

        match *var_parts {
            VarParts::Pos(i) => {
                self.A.add_triplet(self.cur_row, i, 1.);
            }

            VarParts::Neg(i) => {
                self.A.add_triplet(self.cur_row, i, -1.);
            }

            VarParts::PosNeg(i1, i2) => {
                self.A.add_triplet(self.cur_row, i1, 1.);
                self.A.add_triplet(self.cur_row, i2, -1.);
            }
        }

        self.A.add_triplet(
            self.cur_row,
            self.cur_slack_col,
            if lower { -1. } else { 1. },
        );

        self.cur_row += 1;
        self.cur_slack_col -= 1;
    }

    fn add_constraint(&mut self, constraint: &Constraint) {
        self.b_ind.push(self.cur_row);
        self.b_data.push(constraint.rhs);

        for (id, coeff) in &constraint.coeffs {
            let var_parts = self.var_map.get(id).unwrap();

            match *var_parts {
                VarParts::Pos(i) => {
                    self.A.add_triplet(self.cur_row, i, *coeff);
                }

                VarParts::Neg(i) => {
                    self.A.add_triplet(self.cur_row, i, -coeff);
                }

                VarParts::PosNeg(i1, i2) => {
                    self.A.add_triplet(self.cur_row, i1, *coeff);
                    self.A.add_triplet(self.cur_row, i2, -coeff);
                }
            }
        }

        if let Some(slack_coeff) = match constraint.op {
            ConstraintOp::Lte => Some(1.),
            ConstraintOp::Eq => None,
            ConstraintOp::Gte => Some(-1.),
        } {
            self.A
                .add_triplet(self.cur_row, self.cur_slack_col, slack_coeff);
            self.cur_slack_col -= 1;
        }

        self.cur_row += 1;
    }

    fn split_vars(prob: &Problem) -> (VarMap, usize, usize) {
        let n = prob.vars().len();
        let mut var_map = HashMap::with_capacity(n);
        let mut next_col = n;
        let mut num_extra_constraints = 0;

        for (i, var) in prob.vars().iter().enumerate() {
            match var.bound {
                Bound::Free => {
                    var_map.insert(var.id, VarParts::PosNeg(i, next_col));
                    next_col += 1;
                }

                Bound::Lower(lb) => {
                    if lb >= 0. {
                        var_map.insert(var.id, VarParts::Pos(i));
                    } else {
                        var_map.insert(var.id, VarParts::PosNeg(i, next_col));
                        next_col += 1;
                    }

                    if lb != 0. {
                        num_extra_constraints += 1;
                    }
                }

                Bound::Upper(ub) => {
                    if ub >= 0. {
                        var_map.insert(var.id, VarParts::PosNeg(i, next_col));
                        next_col += 1;
                    } else {
                        var_map.insert(var.id, VarParts::Neg(i));
                    }

                    if ub != 0. {
                        num_extra_constraints += 1;
                    }
                }

                Bound::TwoSided(lb, ub) => {
                    if lb == 0. {
                        num_extra_constraints += 1;
                    } else {
                        num_extra_constraints += 2;
                    }

                    match (lb >= 0., ub >= 0.) {
                        (true, true) => {
                            var_map.insert(var.id, VarParts::Pos(i));
                        }

                        (true, false) => panic!("invalid bound: ({}, {})", lb, ub),

                        (false, true) => {
                            var_map.insert(var.id, VarParts::PosNeg(i, next_col));
                            next_col += 1;
                        }

                        (false, false) => {
                            var_map.insert(var.id, VarParts::Neg(i));
                        }
                    }
                }
            }
        }

        (
            var_map,
            num_extra_constraints,
            next_col + num_extra_constraints, //add slack variables
        )
    }
}

impl<T> std::convert::From<T> for StandardForm
where
    T: Borrow<Problem>,
{
    fn from(prob: T) -> StandardForm {
        let prob = prob.borrow();
        let mut builder = StandardFormBuilder::new(prob);

        for var in prob.vars() {
            builder.add_obj_coeffs(var);

            match var.bound {
                Bound::Free => (),
                Bound::Lower(lb) => builder.add_bound(var, lb, true),
                Bound::Upper(ub) => builder.add_bound(var, ub, false),
                Bound::TwoSided(lb, ub) => {
                    builder.add_bound(var, lb, true);
                    builder.add_bound(var, ub, false);
                }
            }
        }

        for constraint in prob.constraints() {
            builder.add_constraint(constraint);
        }

        builder.build()
    }
}

impl StandardForm {
    pub fn rows(&self) -> usize {
        self.A.rows()
    }

    pub fn cols(&self) -> usize {
        self.A.cols()
    }

    //returns a feasible point
    pub fn to_phase_1(&self) -> (StandardForm, ndarray::Array1<f64>) {
        let num_vars = self.cols();
        let num_new_vars = self.rows();
        let total_vars = num_vars + num_new_vars;

        let c_ind: Vec<_> = (num_vars..total_vars).collect();
        let c_data = (0..c_ind.len()).map(|_x| 1.).collect();
        let c = sprs::CsVec::new(total_vars, c_ind, c_data);

        //note: appending the columns depends on A having a csc representation
        assert!(self.A.outer_dims() == self.cols());

        let m = self.rows();
        let mut A = self.A.clone();

        println!("phase1 before\n{:?}\n", A.to_dense());

        let mut b_iter = self.b.iter().peekable();

        for i in 0..num_new_vars {
            let elem = if b_iter.peek().map(|&(ind, _b_i)| ind == i).unwrap_or(false) {
                let (_ind, &b_i) = b_iter.next().unwrap();

                if b_i >= 0. {
                    1.
                } else {
                    -1.
                }
            } else {
                1.
            };

            let new_col = sprs::CsVec::new(m, vec![i], vec![elem]);
            A = A.append_outer_csvec(new_col.view());
        }

        println!("phase1\n{:?}", A.to_dense());

        let mut x0 = ndarray::Array1::zeros(total_vars);

        for (i, &b_i) in self.b.iter() {
            x0[num_vars + i] = b_i.abs();
        }

        (
            StandardForm {
                A,
                b: self.b.clone(),
                c,
                var_map: self.var_map.clone(),
            },
            x0,
        )
    }
}
