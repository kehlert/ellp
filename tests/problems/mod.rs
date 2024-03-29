#[cfg(feature = "benchmarks")]
use ellp::parse_mps;

use ellp::*;

const EPS: f64 = 0.00000001;
const REL_EPS: f64 = 0.000001;

macro_rules! assert_optimal {
    ($prob:expr, $expected_obj:expr, $expected_x:expr) => {
        TestProblem::new($prob, |result: &SolverResult| match result {
            SolverResult::Optimal(sol) => {
                assert!(
                    (sol.obj() - $expected_obj).abs() < EPS,
                    "obj: {}, expected: {}",
                    sol.obj(),
                    $expected_obj
                );

                let x = sol.x();

                //added this line because compiler couldn't always figure out the type
                let expected_x: &[f64] = $expected_x;

                assert_eq!(x.len(), expected_x.len());

                for (x1, x2) in x.iter().zip(expected_x) {
                    assert!((x1 - x2).abs() < EPS, "x_i: {}, expected: {}", x1, x2);
                }
            }

            _ => panic!("not optimal: {:?}", result),
        })
    };
}

macro_rules! assert_unbounded {
    ($prob:expr) => {
        TestProblem::new($prob, |result: &SolverResult| match result {
            SolverResult::Unbounded => (),
            _ => panic!("not unbounded: {:?}", result),
        })
    };
}

macro_rules! assert_optimal_obj {
    ($prob:expr, $expected_obj:expr) => {
        TestProblem::new($prob, |result: &SolverResult| match result {
            SolverResult::Optimal(sol) => {
                assert!(
                    (sol.obj() - $expected_obj).abs() < EPS
                        || (sol.obj() / $expected_obj - 1.).abs() < REL_EPS,
                    "obj: {}, expected: {}",
                    sol.obj(),
                    $expected_obj
                );
            }

            _ => panic!("not optimal: {:?}", result),
        })
    };
}

macro_rules! assert_infeasible {
    ($prob:expr) => {
        TestProblem::new($prob, |result: &SolverResult| match result {
            SolverResult::Infeasible => (),
            _ => panic!("not infeasible: {:?}", result),
        })
    };
}

#[cfg(feature = "mps")]
fn read_mps_file(problem_name: &str) -> Result<Problem, MpsParsingError> {
    use std::io::Read;

    let mut mps_path = std::env::current_dir().unwrap();
    mps_path.push("tests");
    mps_path.push("benchmark_problems");
    mps_path.push(problem_name);
    mps_path.push(format!("{}.mps", problem_name));

    let f = std::fs::File::open(mps_path).unwrap();
    let mut reader = std::io::BufReader::new(f);
    let mut mps = Vec::new();
    reader.read_to_end(&mut mps).unwrap();
    let mps = std::str::from_utf8(&mps).unwrap();
    parse_mps(mps)
}

#[allow(dead_code)]
fn setup_logger(log_level: log::LevelFilter) {
    use fern::colors::{Color, ColoredLevelConfig};

    let colors = ColoredLevelConfig::new()
        .debug(Color::White)
        .info(Color::Green)
        .warn(Color::BrightYellow)
        .error(Color::BrightRed);

    //ignore the result so setup_logger can be called more than once when running multiple tests
    let _ = fern::Dispatch::new()
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
        .apply();
}

pub struct TestProblem {
    pub prob: Problem,
    pub check_result: Box<dyn FnOnce(&SolverResult)>,
}

impl TestProblem {
    fn new<F: FnOnce(&SolverResult) + 'static>(prob: Problem, check_result: F) -> Self {
        Self {
            prob,
            check_result: Box::new(check_result),
        }
    }
}

pub fn empty_problem() -> TestProblem {
    let prob = Problem::new();
    assert_optimal!(prob, 0., &[])
}

pub fn one_variable_no_constraints() -> TestProblem {
    let mut prob = Problem::new();

    prob.add_var(2., Bound::TwoSided(-1., 1.), Some("x1".to_string()))
        .unwrap();

    assert_optimal!(prob, -2., &[-1.])
}

pub fn one_variable_infeasible() -> TestProblem {
    let mut prob = Problem::new();

    let x1 = prob
        .add_var(2., Bound::Upper(0.), Some("x1".to_string()))
        .unwrap();
    prob.add_constraint(vec![(x1, 1.)], ConstraintOp::Gte, 1.)
        .unwrap();

    assert_infeasible!(prob)
}

pub fn one_variable_unbounded_upper() -> TestProblem {
    let mut prob = Problem::new();

    prob.add_var(2., Bound::Upper(0.), Some("x1".to_string()))
        .unwrap();

    assert_unbounded!(prob)
}

pub fn one_variable_unbounded_free() -> TestProblem {
    let mut prob = Problem::new();

    prob.add_var(2., Bound::Free, Some("x1".to_string()))
        .unwrap();

    assert_unbounded!(prob)
}

pub fn two_variables_unbounded() -> TestProblem {
    let mut prob = Problem::new();

    prob.add_var(2., Bound::Lower(0.), Some("x1".to_string()))
        .unwrap();

    prob.add_var(2., Bound::Upper(1.), Some("x2".to_string()))
        .unwrap();

    assert_unbounded!(prob)
}

pub fn two_variables_infeasible_with_bounds() -> TestProblem {
    let mut prob = Problem::new();

    let x1 = prob
        .add_var(2., Bound::Lower(0.), Some("x1".to_string()))
        .unwrap();

    let x2 = prob
        .add_var(2., Bound::Lower(1.), Some("x2".to_string()))
        .unwrap();

    prob.add_constraint(vec![(x1, 1.), (x2, 1.)], ConstraintOp::Lte, 0.)
        .unwrap();

    assert_infeasible!(prob)
}

pub fn two_variables_infeasible_free() -> TestProblem {
    let mut prob = Problem::new();

    let x1 = prob
        .add_var(2., Bound::Free, Some("x1".to_string()))
        .unwrap();

    let x2 = prob
        .add_var(2., Bound::Free, Some("x2".to_string()))
        .unwrap();

    prob.add_constraint(vec![(x1, 1.), (x2, 1.)], ConstraintOp::Eq, -1.)
        .unwrap();

    prob.add_constraint(vec![(x1, 2.), (x2, 2.)], ConstraintOp::Eq, 1.)
        .unwrap();

    assert_infeasible!(prob)
}

pub fn infeasible_constraint_without_coeffs() -> TestProblem {
    let mut prob = Problem::new();

    prob.add_var(2., Bound::Free, Some("x1".to_string()))
        .unwrap();

    prob.add_constraint(vec![], ConstraintOp::Eq, 1.).unwrap();

    assert_infeasible!(prob)
}

pub fn feasible_constraint_without_coeffs() -> TestProblem {
    let mut prob = Problem::new();

    prob.add_var(2., Bound::Lower(3.), Some("x1".to_string()))
        .unwrap();

    prob.add_constraint(vec![], ConstraintOp::Eq, 0.).unwrap();

    assert_optimal!(prob, 6., &[3.])
}

pub fn feasible_constraint_without_coeffs_and_no_vars() -> TestProblem {
    let mut prob = Problem::new();
    prob.add_constraint(vec![], ConstraintOp::Eq, 0.).unwrap();
    assert_optimal!(prob, 0., &[])
}

pub fn infeasible_constraint_without_coeffs_and_no_vars() -> TestProblem {
    let mut prob = Problem::new();
    prob.add_constraint(vec![], ConstraintOp::Eq, 1.).unwrap();
    assert_infeasible!(prob)
}

pub fn linear_system_2d() -> TestProblem {
    let mut prob = Problem::new();

    let x = prob
        .add_var(0., Bound::Free, Some("x".to_string()))
        .unwrap();

    let y = prob
        .add_var(0., Bound::Free, Some("y".to_string()))
        .unwrap();

    prob.add_constraint(vec![(x, 2.), (y, 1.)], ConstraintOp::Eq, 1.)
        .unwrap();

    prob.add_constraint(vec![(x, 3.), (y, 1.)], ConstraintOp::Eq, 1.)
        .unwrap();

    assert_optimal!(prob, 0., &[0., 1.])
}

pub fn linear_system_3d() -> TestProblem {
    let mut prob = Problem::new();

    let x = prob
        .add_var(0., Bound::Free, Some("x".to_string()))
        .unwrap();

    let y = prob
        .add_var(0., Bound::Free, Some("y".to_string()))
        .unwrap();

    let z = prob
        .add_var(0., Bound::Free, Some("z".to_string()))
        .unwrap();

    prob.add_constraint(vec![(x, 1.), (y, 2.), (z, 4.)], ConstraintOp::Eq, 1.)
        .unwrap();

    prob.add_constraint(vec![(x, 3.), (y, 4.), (z, 8.)], ConstraintOp::Eq, 2.)
        .unwrap();

    prob.add_constraint(vec![(x, 5.), (y, 6.), (z, 13.)], ConstraintOp::Eq, 5.)
        .unwrap();

    assert_optimal!(prob, 0., &[0., -3.5, 2.])
}

pub fn linear_system_3d_infeasible() -> TestProblem {
    let mut prob = Problem::new();

    let x = prob
        .add_var(0., Bound::Free, Some("x".to_string()))
        .unwrap();

    let y = prob
        .add_var(0., Bound::Free, Some("y".to_string()))
        .unwrap();

    let z = prob
        .add_var(0., Bound::Free, Some("z".to_string()))
        .unwrap();

    prob.add_constraint(vec![(x, 1.), (y, 2.), (z, 4.)], ConstraintOp::Eq, 1.)
        .unwrap();

    prob.add_constraint(vec![(x, 3.), (y, 4.), (z, 8.)], ConstraintOp::Eq, 2.)
        .unwrap();

    prob.add_constraint(vec![(x, 5.), (y, 6.), (z, 12.)], ConstraintOp::Eq, 5.)
        .unwrap();

    assert_infeasible!(prob)
}

pub fn small_prob_1() -> TestProblem {
    //setup_logger(log::LevelFilter::Trace);

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

    assert_optimal!(
        prob,
        19.1578947368421,
        &[-0.94736842105, 2.105263157894, 0., 0., -0.5]
    )
}

pub fn small_prob_2() -> TestProblem {
    // setup_logger(log::LevelFilter::Trace);
    let mut prob = Problem::new();

    let x = prob
        .add_var(-5., Bound::Lower(0.), Some("x".to_string()))
        .unwrap();

    let y = prob
        .add_var(-4., Bound::Lower(0.), Some("y".to_string()))
        .unwrap();

    prob.add_constraint(vec![(x, 1.)], ConstraintOp::Lte, 6.)
        .unwrap();

    prob.add_constraint(vec![(x, 0.25), (y, 1.)], ConstraintOp::Lte, 6.)
        .unwrap();

    prob.add_constraint(vec![(x, 3.), (y, 2.)], ConstraintOp::Lte, 22.)
        .unwrap();

    assert_optimal!(prob, -40., &[4., 5.])
}

pub fn small_prob_3() -> TestProblem {
    // setup_logger(log::LevelFilter::Trace);
    let mut prob = Problem::new();

    let x = prob
        .add_var(3., Bound::Lower(0.), Some("x".to_string()))
        .unwrap();

    let y = prob
        .add_var(-6., Bound::Lower(0.), Some("y".to_string()))
        .unwrap();

    prob.add_constraint(vec![(x, 1.), (y, 2.)], ConstraintOp::Gte, -1.)
        .unwrap();

    prob.add_constraint(vec![(x, 2.), (y, 1.)], ConstraintOp::Gte, 0.)
        .unwrap();

    prob.add_constraint(vec![(x, 1.), (y, -1.)], ConstraintOp::Gte, -1.)
        .unwrap();

    prob.add_constraint(vec![(x, 1.), (y, -4.)], ConstraintOp::Gte, -13.)
        .unwrap();

    prob.add_constraint(vec![(x, -4.), (y, 1.)], ConstraintOp::Gte, -23.)
        .unwrap();

    assert_optimal!(prob, -15., &[3., 4.])
}

pub fn small_prob_4() -> TestProblem {
    //NOTE: this problem has multiple optimal points, so we only test the objective value
    let mut prob = Problem::new();

    let x = prob
        .add_var(-1., Bound::Lower(0.), Some("x".to_string()))
        .unwrap();

    let y = prob
        .add_var(-1., Bound::Lower(0.), Some("y".to_string()))
        .unwrap();

    let z = prob
        .add_var(-1., Bound::Lower(0.), Some("z".to_string()))
        .unwrap();

    prob.add_constraint(vec![(x, 1.), (y, -1.), (z, 1.)], ConstraintOp::Gte, -2.)
        .unwrap();

    prob.add_constraint(vec![(x, -1.), (y, 1.), (z, 1.)], ConstraintOp::Gte, -3.)
        .unwrap();

    prob.add_constraint(vec![(x, 1.), (y, 1.), (z, -1.)], ConstraintOp::Gte, -1.)
        .unwrap();

    prob.add_constraint(vec![(x, -1.), (y, -1.), (z, -1.)], ConstraintOp::Gte, -4.)
        .unwrap();

    assert_optimal_obj!(prob, -4.)
}

pub fn small_prob_5() -> TestProblem {
    let mut prob = Problem::new();

    let x = prob
        .add_var(4., Bound::Lower(0.), Some("x".to_string()))
        .unwrap();

    let y = prob
        .add_var(5., Bound::Lower(0.), Some("y".to_string()))
        .unwrap();

    prob.add_constraint(vec![(x, 1.), (y, 1.)], ConstraintOp::Gte, -1.)
        .unwrap();

    prob.add_constraint(vec![(x, 1.), (y, 2.)], ConstraintOp::Gte, 1.)
        .unwrap();

    prob.add_constraint(vec![(x, 4.), (y, 2.)], ConstraintOp::Gte, 8.)
        .unwrap();

    prob.add_constraint(vec![(x, -1.), (y, -1.)], ConstraintOp::Gte, -3.)
        .unwrap();

    prob.add_constraint(vec![(x, -1.), (y, 1.)], ConstraintOp::Gte, 1.)
        .unwrap();

    assert_optimal!(prob, 14., &[1., 2.])
}

pub fn small_prob_6() -> TestProblem {
    let mut prob = Problem::new();

    let x = prob
        .add_var(-2., Bound::Lower(0.), Some("x".to_string()))
        .unwrap();

    let y = prob
        .add_var(-4., Bound::Lower(0.), Some("y".to_string()))
        .unwrap();

    let z = prob
        .add_var(-1., Bound::Lower(0.), Some("z".to_string()))
        .unwrap();

    let w = prob
        .add_var(-1., Bound::Lower(0.), Some("w".to_string()))
        .unwrap();

    prob.add_constraint(vec![(x, -1.), (y, -3.), (w, -1.)], ConstraintOp::Gte, -4.)
        .unwrap();

    prob.add_constraint(vec![(x, -2.), (y, -1.)], ConstraintOp::Gte, -3.)
        .unwrap();

    prob.add_constraint(vec![(y, -1.), (z, -4.), (w, -1.)], ConstraintOp::Gte, -3.)
        .unwrap();

    prob.add_constraint(vec![(x, 1.), (y, 1.), (z, 2.)], ConstraintOp::Gte, 1.)
        .unwrap();

    prob.add_constraint(vec![(x, -1.), (y, 1.), (z, 4.)], ConstraintOp::Gte, 1.)
        .unwrap();

    assert_optimal!(prob, -6.5, &[1., 1., 0.5, 0.])
}

pub fn small_prob_7() -> TestProblem {
    // setup_logger(log::LevelFilter::Trace);

    let mut prob = Problem::new();

    let x = prob
        .add_var(2., Bound::Lower(0.), Some("x".to_string()))
        .unwrap();

    let y = prob
        .add_var(-1., Bound::Lower(0.), Some("y".to_string()))
        .unwrap();

    let z = prob
        .add_var(1., Bound::Free, Some("z".to_string()))
        .unwrap();

    prob.add_constraint(vec![(x, 1.), (y, -1.), (z, 4.)], ConstraintOp::Gte, -1.)
        .unwrap();

    prob.add_constraint(vec![(x, 1.), (y, -1.), (z, -1.)], ConstraintOp::Gte, 2.)
        .unwrap();

    prob.add_constraint(vec![(x, 1.), (y, 3.), (z, 2.)], ConstraintOp::Eq, 3.)
        .unwrap();

    assert_optimal!(prob, 2.9, &[2.1, 0.7, -0.6])
}

pub fn small_prob_unbounded_1() -> TestProblem {
    // setup_logger(log::LevelFilter::Trace);
    let mut prob = Problem::new();

    let x = prob
        .add_var(-2., Bound::Lower(0.), Some("x".to_string()))
        .unwrap();

    let y = prob
        .add_var(-3., Bound::Lower(0.), Some("y".to_string()))
        .unwrap();

    let z = prob
        .add_var(1., Bound::Lower(0.), Some("z".to_string()))
        .unwrap();

    prob.add_constraint(vec![(x, 1.), (y, 1.), (z, 1.)], ConstraintOp::Gte, -3.)
        .unwrap();

    prob.add_constraint(vec![(x, -1.), (y, 1.), (z, -1.)], ConstraintOp::Gte, -4.)
        .unwrap();

    prob.add_constraint(vec![(x, 1.), (y, -1.), (z, -2.)], ConstraintOp::Gte, -1.)
        .unwrap();

    assert_unbounded!(prob)
}

pub fn small_prob_unbounded_2() -> TestProblem {
    // setup_logger(log::LevelFilter::Trace);
    let mut prob = Problem::new();

    let x = prob
        .add_var(-2., Bound::Lower(0.), Some("x".to_string()))
        .unwrap();

    let y = prob
        .add_var(-3., Bound::Lower(0.), Some("y".to_string()))
        .unwrap();

    let z = prob
        .add_var(1., Bound::Lower(0.), Some("z".to_string()))
        .unwrap();

    let w = prob
        .add_var(1., Bound::Lower(0.), Some("w".to_string()))
        .unwrap();

    prob.add_constraint(vec![(y, 1.), (z, -2.), (w, -1.)], ConstraintOp::Gte, -4.)
        .unwrap();

    prob.add_constraint(
        vec![(x, 2.), (y, -1.), (z, -1.), (w, 4.)],
        ConstraintOp::Gte,
        -5.,
    )
    .unwrap();

    prob.add_constraint(vec![(x, -1.), (y, 1.), (w, -2.)], ConstraintOp::Gte, -3.)
        .unwrap();

    assert_unbounded!(prob)
}

pub fn beale_cycle() -> TestProblem {
    // setup_logger(log::LevelFilter::Trace);

    let mut prob = Problem::new();

    let x = prob
        .add_var(-10., Bound::Lower(0.), Some("x".to_string()))
        .unwrap();

    let y = prob
        .add_var(57., Bound::Lower(0.), Some("y".to_string()))
        .unwrap();

    let z = prob
        .add_var(9., Bound::Lower(0.), Some("z".to_string()))
        .unwrap();

    let w = prob
        .add_var(24., Bound::Lower(0.), Some("w".to_string()))
        .unwrap();

    prob.add_constraint(
        vec![(x, -0.5), (y, 5.5), (z, 2.5), (w, -9.)],
        ConstraintOp::Gte,
        0.,
    )
    .unwrap();

    prob.add_constraint(
        vec![(x, -0.5), (y, 1.5), (z, 0.5), (w, -1.)],
        ConstraintOp::Gte,
        0.,
    )
    .unwrap();

    prob.add_constraint(vec![(x, -1.)], ConstraintOp::Gte, -1.)
        .unwrap();

    assert_optimal!(prob, -1., &[1., 0., 1., 0.])
}

#[cfg(feature = "benchmarks")]
pub fn afiro() -> TestProblem {
    //setup_logger(log::LevelFilter::Info);
    let prob = read_mps_file("afiro").unwrap();
    assert_optimal_obj!(prob, -4.6475314286E+02)
}

#[cfg(feature = "benchmarks")]
pub fn adlittle() -> TestProblem {
    let prob = read_mps_file("adlittle").unwrap();
    assert_optimal_obj!(prob, 2.2549496316E+05)
}

#[cfg(feature = "benchmarks")]
pub fn blend() -> TestProblem {
    let prob = read_mps_file("blend").unwrap();
    assert_optimal_obj!(prob, -3.0812149846E+01)
}
