mod problems;

use ellp::*;
use paste::paste;

#[allow(dead_code)]
pub fn setup_logger(log_level: log::LevelFilter) {
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

macro_rules! generate_tests {
    ($solver_name:ident, $solver:expr, $($problem:ident,)+) => {
        paste! {
            $(
                #[test]
                fn  [<$solver_name _ $problem>]() {
                    //setup_logger(log::LevelFilter::Trace);
                    let test_prob = problems::$problem();
                    let solver = $solver;
                    let result = solver.solve(test_prob.prob).unwrap();
                    (test_prob.check_result)(&result)
                }
            )+
        }
    };
}

generate_tests! {
    primal,
    PrimalSimplexSolver::default(),
    empty_problem,
    one_variable_no_constraints,
    one_variable_infeasible,
    one_variable_unbounded_upper,
    one_variable_unbounded_free,
    two_variables_unbounded,
    two_variables_infeasible_with_bounds,
    two_variables_infeasible_free,
    infeasible_constraint_without_coeffs,
    feasible_constraint_without_coeffs,
    feasible_constraint_without_coeffs_and_no_vars,
    infeasible_constraint_without_coeffs_and_no_vars,
    linear_system_2d,
    linear_system_3d,
    linear_system_3d_infeasible,
    small_prob_1,
    small_prob_2,
    small_prob_3,
    small_prob_4,
    small_prob_5,
    small_prob_6,
    small_prob_7,
    small_prob_unbounded_1,
    small_prob_unbounded_2,
    beale_cycle,
}

generate_tests! {
    dual,
    DualSimplexSolver::default(),
    empty_problem,
    one_variable_no_constraints,
    one_variable_infeasible,
    one_variable_unbounded_upper,
    one_variable_unbounded_free,
    two_variables_unbounded,
    two_variables_infeasible_with_bounds,
    two_variables_infeasible_free,
    infeasible_constraint_without_coeffs,
    feasible_constraint_without_coeffs,
    feasible_constraint_without_coeffs_and_no_vars,
    infeasible_constraint_without_coeffs_and_no_vars,
    linear_system_2d,
    linear_system_3d,
    linear_system_3d_infeasible,
    //small_prob_1,
    small_prob_2,
    //small_prob_3,
    // small_prob_4,
    // small_prob_5,
    // small_prob_6,
    // small_prob_7,
    small_prob_unbounded_1,
    small_prob_unbounded_2,
    // beale_cycle,
}
