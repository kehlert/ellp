#![allow(non_snake_case)]

use crate::standard_form::{InitialPoint, StandardForm};
use crate::util::EPS;
use crate::{error::EllPError, problem::Problem};

pub struct Solver {}

impl std::default::Default for Solver {
    fn default() -> Self {
        Self {}
    }
}

impl Solver {
    pub fn new() -> Self {
        Self {}
    }

    pub fn solve(
        &self,
        prob: Problem,
        x0: Option<nalgebra::DVector<f64>>,
    ) -> Result<(), EllPError> {
        let n = prob.vars().len();
        let mut x = x0.unwrap_or_else(|| nalgebra::DVector::zeros(n));

        if x.len() != n {
            return Err(EllPError::new(format!(
                "x0 dimensions invalid: {}, expected: {}",
                x.len(),
                n
            )));
        }

        if !prob.is_feasible(&x) {
            self.find_feasible_point(prob);
        }

        todo!()
    }

    fn solve_std_form(
        &self,
        std_form: &StandardForm,
        init: InitialPoint, //x0 assumed to be feasible
    ) -> Result<(), EllPError> {
        let B = init.B;
        let N = init.N;
        let x0 = init.x;

        if B.len() != std_form.rows() {
            return Err(EllPError::new(format!(
                "invalid B, has {} elements but {} expected",
                B.len(),
                std_form.rows(),
            )));
        }

        if N.len() != std_form.cols().checked_sub(std_form.rows()).unwrap() {
            return Err(EllPError::new(format!(
                "invalid N, has {} elements but {} expected",
                N.len(),
                std_form.rows(),
            )));
        }

        let B_cols: Vec<_> = B.iter().map(|i| std_form.A.column(*i)).collect();
        let A_B = nalgebra::DMatrix::from_columns(&B_cols);

        let N_cols: Vec<_> = N.iter().map(|i| std_form.A.column(*i)).collect();
        let A_N = nalgebra::DMatrix::from_columns(&N_cols);

        let lu_decomp = A_B.lu();

        let L = lu_decomp.l();
        let U = lu_decomp.u();

        if U.diagonal().iter().any(|d| d.abs() < EPS) {
            return Err(EllPError::new(
                "invalid B, A_B is not invertible".to_string(),
            ));
        }

        let c_B = nalgebra::DVector::from_iterator(B.len(), B.iter().map(|i| std_form.c[*i]));
        let c_N = nalgebra::DVector::from_iterator(N.len(), N.iter().map(|i| std_form.c[*i]));

        //should always have a solution
        let u = lu_decomp.solve(&c_B).unwrap();
        println!("u: {}", u);

        println!(
            "{}, ({}, {}), {}",
            c_N.nrows(),
            A_N.nrows(),
            A_N.ncols(),
            u.nrows()
        );

        let r = &c_N - A_N.transpose() * u;

        //TODO find nonbasic to enter the basis

        println!("r: {}", r);

        //TODO monitor infeasibility

        todo!()
    }

    fn find_feasible_point(&self, prob: Problem) -> Result<(), EllPError> {
        let (phase_1_std_form, x0) = StandardForm::phase_1(prob);

        println!("{:#?}", x0);

        println!("lhs: {}", &phase_1_std_form.A * &x0.x);
        println!("b: {}", phase_1_std_form.b);

        println!("A: {}", phase_1_std_form.A);
        println!("bounds:\n{:?}", phase_1_std_form.bounds);

        //should always be feasbible and bounded
        self.solve_std_form(&phase_1_std_form, x0)
    }
}
