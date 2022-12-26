use crate::error::EllPError;
use crate::standard_form::{BasicPoint, Point, StandardForm};

pub type EllPResult = Result<SolverResult, EllPError>;

#[derive(Debug)]
pub enum SolverResult {
    Optimal(Solution),
    Infeasible,
    Unbounded,
    MaxIter { obj: f64 },
}

impl std::fmt::Display for SolverResult {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Self::Optimal(sol) => write!(f, "found optimal point with objective {}", sol.obj()),
            Self::Infeasible => write!(f, "problem is infeasible"),
            Self::Unbounded => write!(f, "problem is unbounded"),
            Self::MaxIter { obj } => {
                write!(f, "reached max iterations, current objective = {}", obj)
            }
        }
    }
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
}

#[derive(Debug)]
pub struct OptimalPoint(Point);

impl OptimalPoint {
    pub fn new(point: Point) -> Self {
        Self(point)
    }
}

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
