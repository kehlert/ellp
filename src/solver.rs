use crate::standard_form::StandardForm;
use crate::{error::EllPError, problem::Problem};

use ndarray::linalg::Dot;

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

    pub fn solve(&self, prob: &Problem, x0: Option<ndarray::Array1<f64>>) -> Result<(), EllPError> {
        let n = prob.vars().len();
        let mut x = x0.unwrap_or_else(|| ndarray::Array1::zeros(n));

        if x.len() != n {
            return Err(EllPError::new(format!(
                "x0 dimensions invalid: {}, expected: {}",
                x.len(),
                n
            )));
        }

        let std_form: StandardForm = prob.into();

        if !prob.is_feasible(&x) {
            self.find_feasible_point(&std_form);
        }

        todo!()
    }

    fn solve_std_form(
        &self,
        std_form: &StandardForm,
        x0: ndarray::Array1<f64>,
    ) -> Result<(), EllPError> {
        todo!()
    }

    fn find_feasible_point(&self, std_form: &StandardForm) -> Result<(), EllPError> {
        let (phase_1_std_form, x0) = std_form.to_phase_1();
        println!("\nx0\n{:?}\n", x0);

        println!("lhs: {:?}", phase_1_std_form.A.dot(&x0));
        println!("b: {:?}", phase_1_std_form.b);

        self.solve_std_form(&phase_1_std_form, x0)
    }
}
