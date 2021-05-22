mod error;
pub mod problem;
pub mod solver;
mod standard_form;
mod util;

pub use crate::problem::{Bound, Constraint, ConstraintOp, Problem, Variable};
pub use crate::solver::Solver;

#[cfg(test)]
mod tests {
    use super::*;

    fn setup_logger(log_level: log::LevelFilter) {
        use fern::colors::{Color, ColoredLevelConfig};
        let colors = ColoredLevelConfig::new()
            .debug(Color::White)
            .info(Color::Green)
            .warn(Color::BrightYellow)
            .error(Color::BrightRed);

        fern::Dispatch::new()
            .format(move |out, message, record| {
                out.finish(format_args!(
                    "{} | {:5} | {}",
                    chrono::Local::now().format("%Y-%m-%d %H:%M:%S%.6f"),
                    colors.color(record.level()),
                    message
                ))
            })
            .level(log_level)
            .chain(std::io::stdout())
            .apply()
            .unwrap();
    }

    #[test]
    fn it_works() {
        setup_logger(log::LevelFilter::Trace);

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

        let solver = Solver::new();
        let result = solver.solve(&prob, None);
        println!("RESULT:\n{:?}", result);

        //TODO test a system where free var constraints are infeasible
    }
}
