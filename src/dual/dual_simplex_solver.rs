use super::dual_problem::{DualFeasiblePoint, DualPhase1, DualPhase2, DualProblem};
use crate::problem::Problem;
use crate::solver::EllPResult;

pub struct DualSimplexSolver {
    max_iter: u64,
}

impl std::default::Default for DualSimplexSolver {
    fn default() -> Self {
        Self { max_iter: 1000 }
    }
}

impl DualSimplexSolver {
    pub fn new(max_iter: Option<u64>) -> Self {
        Self {
            max_iter: max_iter.unwrap_or(u64::MAX),
        }
    }

    pub fn solve(&self, prob: Problem) -> EllPResult {
        let mut phase_1: DualPhase1 = prob.into();
        todo!()
        //self.solve_with_initial(&mut phase_1)?;
    }
}
