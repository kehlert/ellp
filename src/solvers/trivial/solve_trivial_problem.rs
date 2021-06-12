use crate::problem::Bound;
use crate::solver::SolutionStatus;
use crate::standard_form::{Nonbasic, NonbasicBound, StandardForm};

pub fn solve_trivial_problem(
    std_form: &StandardForm,
    x: &mut nalgebra::DVector<f64>,
    #[allow(non_snake_case)] N: &mut Vec<Nonbasic>,
    minimize: bool,
) -> SolutionStatus {
    N.clear();

    assert_eq!(std_form.c.len(), std_form.bounds.len());

    for (i, ((x_i, &c_i), bound)) in x
        .iter_mut()
        .zip(&std_form.c)
        .zip(&std_form.bounds)
        .enumerate()
    {
        match *bound {
            Bound::Free => {
                N.push(Nonbasic::new(i, NonbasicBound::Free));

                if c_i != 0. {
                    return SolutionStatus::Unbounded;
                } else {
                    *x_i = 0.;
                }
            }

            Bound::Lower(lb) => {
                N.push(Nonbasic::new(i, NonbasicBound::Lower));

                match (c_i > 0., minimize) {
                    (true, true) => *x_i = lb,
                    (true, false) => return SolutionStatus::Unbounded,
                    (false, true) => *x_i = lb,
                    (false, false) => {
                        if c_i != 0. {
                            return SolutionStatus::Unbounded;
                        } else {
                            *x_i = lb
                        }
                    }
                }
            }

            Bound::Upper(ub) => {
                N.push(Nonbasic::new(i, NonbasicBound::Upper));

                match (c_i > 0., minimize) {
                    (true, true) => return SolutionStatus::Unbounded,
                    (true, false) => *x_i = ub,
                    (false, true) => {
                        if c_i != 0. {
                            return SolutionStatus::Unbounded;
                        } else {
                            *x_i = ub
                        }
                    }
                    (false, false) => *x_i = ub,
                }
            }

            Bound::TwoSided(lb, ub) => match (c_i > 0., minimize) {
                (true, true) => {
                    N.push(Nonbasic::new(i, NonbasicBound::Lower));
                    *x_i = lb;
                }

                (true, false) => {
                    N.push(Nonbasic::new(i, NonbasicBound::Upper));
                    *x_i = ub;
                }

                (false, true) => {
                    N.push(Nonbasic::new(i, NonbasicBound::Upper));
                    *x_i = ub;
                }

                (false, false) => {
                    N.push(Nonbasic::new(i, NonbasicBound::Lower));
                    *x_i = lb;
                }
            },

            Bound::Fixed(val) => {
                N.push(Nonbasic::new(i, NonbasicBound::Lower));
                *x_i = val;
            }
        }
    }

    SolutionStatus::Optimal
}
