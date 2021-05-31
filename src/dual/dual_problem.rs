use crate::standard_form::{BasicPoint, Point, StandardForm};

pub trait DualFeasible: BasicPoint {}

pub trait DualProblem<P: DualFeasible> {
    fn obj(&self) -> f64;
    fn pt(&self) -> &P;
}

#[derive(Debug, Clone)]
pub struct DualFeasiblePoint(Point);

impl DualFeasible for DualFeasiblePoint {}

impl BasicPoint for DualFeasiblePoint {
    #[inline]
    fn into_pt(self) -> Point {
        self.0
    }
}

impl std::ops::Deref for DualFeasiblePoint {
    type Target = Point;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl std::ops::DerefMut for DualFeasiblePoint {
    #[inline]
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

#[derive(Debug)]
pub struct DualPhase1 {
    pub std_form: StandardForm,
    feasible_point: DualFeasiblePoint,
}

#[derive(Debug)]
pub struct DualPhase2 {
    std_form: StandardForm,
    feasible_point: DualFeasiblePoint,
}

impl std::convert::From<StandardForm> for DualPhase1 {
    fn from(mut _std_form: StandardForm) -> Self {
        todo!()
    }
}

impl std::convert::From<DualPhase1> for DualPhase2 {
    fn from(mut _phase_1: DualPhase1) -> Self {
        todo!()
    }
}
