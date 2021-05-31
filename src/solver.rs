use crate::dual::dual_problem::{DualFeasible, DualPhase2};
use crate::error::EllPError;
use crate::primal::primal_problem::PrimalFeasible;
use crate::problem::Constraint;
use crate::standard_form::{BasicPoint, Point, StandardForm};

pub type EllPResult = Result<SolverResult, EllPError>;

#[derive(Debug)]
pub enum SolverResult {
    Optimal(Solution),
    Infeasible,
    Unbounded,
    MaxIter { obj: f64 },
}

#[derive(Debug)]
pub enum SolutionStatus {
    Optimal,
    Infeasible,
    Unbounded,
    MaxIter,
}

#[derive(Debug)]
pub struct Solution {
    std_form: StandardForm,
    point: OptimalPoint,
}

impl Solution {
    pub fn new(std_form: StandardForm, point: OptimalPoint) -> Self {
        Self { std_form, point }
    }

    #[inline]
    pub fn obj(&self) -> f64 {
        self.std_form.obj(&self.point.x)
    }

    pub fn x(&self) -> nalgebra::DVectorSlice<f64> {
        self.std_form.extract_solution(&self.point)
    }

    pub fn add_constraint(self, _constraint: Constraint) -> DualPhase2 {
        todo!()
    }
}

#[derive(Debug)]
pub struct OptimalPoint(Point);

impl OptimalPoint {
    pub fn new(point: Point) -> Self {
        Self(point)
    }
}

impl PrimalFeasible for OptimalPoint {}

impl DualFeasible for OptimalPoint {}

impl BasicPoint for OptimalPoint {
    #[inline]
    fn into_pt(self) -> Point {
        self.0
    }
}

impl std::ops::Deref for OptimalPoint {
    type Target = Point;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl std::ops::DerefMut for OptimalPoint {
    #[inline]
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}
