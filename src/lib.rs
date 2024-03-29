/*!
Linear programming library that provides primal and dual simplex solvers. Both solvers are currently working for a small set of test problems. This library is an *early work-in-progress*.

## Examples

Here is example code that sets up a linear program, and then solves it with both the primal and dual simplex solvers.

First we setup the problem


```
use ellp::*;

let mut prob = Problem::new();

let x1 = prob
    .add_var(2., Bound::TwoSided(-1., 1.), Some("x1".to_string()))
    .unwrap();

let x2 = prob
    .add_var(10., Bound::Upper(6.), Some("x2".to_string()))
    .unwrap();

let x3 = prob
    .add_var(0., Bound::Lower(0.), Some("x3".to_string()))
    .unwrap();

let x4 = prob
    .add_var(1., Bound::Fixed(0.), Some("x4".to_string()))
    .unwrap();

let x5 = prob
    .add_var(0., Bound::Free, Some("x5".to_string()))
    .unwrap();

prob.add_constraint(vec![(x1, 2.5), (x2, 3.5)], ConstraintOp::Gte, 5.)
    .unwrap();

prob.add_constraint(vec![(x2, 2.5), (x1, 4.5)], ConstraintOp::Lte, 1.)
    .unwrap();

prob.add_constraint(vec![(x3, -1.), (x4, -3.), (x5, -4.)], ConstraintOp::Eq, 2.)
    .unwrap();

println!("{}", prob);

let primal_solver = PrimalSimplexSolver::default();
let dual_solver = DualSimplexSolver::default();

let primal_result = primal_solver.solve(prob.clone()).unwrap();
let dual_result = dual_solver.solve(prob).unwrap();

if let SolverResult::Optimal(sol) = primal_result {
    println!("primal obj: {}", sol.obj());
    println!("primal opt point: {}", sol.x());
} else {
    panic!("should have an optimal point");
}

if let SolverResult::Optimal(sol) = dual_result {
    println!("dual obj: {}", sol.obj());
    println!("dual opt point: {}", sol.x());
} else {
    panic!("should have an optimal point");
}
```

The output is
```console
minimize
+ 2 x1 + 10 x2 + 1 x4

subject to
+ 2.5 x1 + 3.5 x2 ≥ 5
+ 2.5 x2 + 4.5 x1 ≤ 1
- 1 x3 - 3 x4 - 4 x5 = 2

with the bounds
-1 ≤ x1 ≤ 1
x2 ≤ 6
x3 ≥ 0
x4 = 0
x5 free

primal obj: 19.157894736842103
primal opt point:
  ┌                     ┐
  │ -0.9473684210526313 │
  │  2.1052631578947367 │
  │                   0 │
  │                   0 │
  │                -0.5 │
  └                     ┘

dual obj: 19.157894736842103
dual opt point:
  ┌                     ┐
  │ -0.9473684210526313 │
  │  2.1052631578947367 │
  │                   0 │
  │                   0 │
  │                -0.5 │
  └                     ┘
```

If the problem is infeasible or unbounded, then `solve` will return [`SolverResult::Infeasible`] or [`SolverResult::Unbounded`], respectively.
*/

mod error;
pub mod problem;
pub mod solver;
mod solvers;
mod standard_form;
mod util;

#[cfg(feature = "mps")]
mod parse_mps;

#[cfg(feature = "mps")]
pub use parse_mps::parse_mps;

#[cfg(feature = "mps")]
pub use parse_mps::MpsParsingError;

pub use crate::error::EllPError;
pub use crate::problem::{Bound, Constraint, ConstraintOp, Problem, Variable};
pub use crate::solver::{EllPResult, SolverResult};
pub use crate::solvers::dual::dual_simplex_solver::DualSimplexSolver;
pub use crate::solvers::primal::primal_simplex_solver::PrimalSimplexSolver;
