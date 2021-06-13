use crate::problem::{Bound, Constraint, ConstraintOp, Problem, Variable};

use nom::branch::alt;
use nom::bytes::complete::{tag, take_till1, take_until, take_while1};
use nom::character::complete::{char, multispace0, multispace1, satisfy, space0, space1};
use nom::character::{is_alphanumeric, is_newline};
use nom::multi::{fold_many1, separated_list0};
use nom::number::complete::double;
use nom::Finish;
use nom::IResult;

use log::error;
use thiserror::Error;

use std::collections::HashMap;

type Rows<'a> = HashMap<&'a [u8], Row<'a>>;
type Cols<'a> = HashMap<&'a [u8], Col>;

#[derive(Error, Debug)]
#[error("MPS parsing error. {msg}\n{code:?}\ninput: {input:?}")]
pub struct MpsParsingError<'a> {
    msg: String,
    input: Option<std::borrow::Cow<'a, str>>,
    code: Option<nom::error::ErrorKind>,
}

#[allow(dead_code)]
pub fn parse_mps(mps: &[u8]) -> Result<Problem, MpsParsingError> {
    match parse_mps_impl(mps).finish() {
        Ok((_i, (rows, cols))) => {
            let mut prob = Problem::new();
            let mut var_ids = HashMap::new();

            for (var_name_bytes, col) in cols {
                let var_name = String::from_utf8(var_name_bytes.to_owned()).map_err(|err| {
                    MpsParsingError {
                        msg: format!("invalid variable name: {}", err),
                        input: Some(String::from_utf8_lossy(var_name_bytes)),
                        code: None,
                    }
                })?;

                let obj_coeff = col.obj_coeff.unwrap_or(0.);
                let bound = col.bound.unwrap_or(Bound::Free);

                //should never panic, we know that the var names are unique
                let var_id = prob.add_var(obj_coeff, bound, Some(var_name)).unwrap();
                var_ids.insert(var_name_bytes, var_id);
            }

            for (row_name, row) in rows {
                match row {
                    Row::Objective => continue,
                    Row::Constraint { op, coeffs, rhs } => {
                        let rhs = rhs.ok_or_else(|| MpsParsingError {
                            msg: format!(
                                "did not specify right-hand side for {}",
                                String::from_utf8_lossy(row_name)
                            ),
                            input: None,
                            code: None,
                        })?;

                        let coeffs = coeffs
                            .into_iter()
                            .map(|(var_name, coeff)| match var_ids.get(var_name) {
                                Some(var_id) => Ok((*var_id, coeff)),

                                None => Err(MpsParsingError {
                                    msg: format!(
                                        "for row {}, column {} does not exist",
                                        String::from_utf8_lossy(row_name),
                                        String::from_utf8_lossy(var_name)
                                    ),
                                    input: None,
                                    code: None,
                                }),
                            })
                            .collect::<Result<_, _>>()?;

                        //should never panic, we know the variable ids in coeffs exist
                        prob.add_constraint(coeffs, op, rhs).unwrap();
                    }
                }
            }

            Ok(prob)
        }

        Err(err) => {
            let input = String::from_utf8_lossy(err.input);

            Err(MpsParsingError {
                msg: "parsing error".to_string(),
                input: Some(input),
                code: Some(err.code),
            })
        }
    }
}

fn parse_mps_impl(i: &[u8]) -> IResult<&[u8], (Rows, Cols)> {
    let (i, _whitespace) = multispace0(i)?;
    let (i, _name) = parse_name(i)?;

    let (i, _whitespace) = multispace0(i)?;
    let (i, mut rows) = parse_rows(i)?;

    let (i, _whitespace) = multispace0(i)?;
    let (i, mut cols) = parse_columns(i, &mut rows)?;

    let (i, _whitespace) = multispace0(i)?;
    let (i, ()) = parse_all_rhs(i, &mut rows)?;

    let (i, _whitespace) = multispace0(i)?;
    let (i, ()) = parse_bounds(i, &mut cols)?;

    let (i, _whitespace) = multispace0(i)?;
    let (i, _) = tag("ENDATA")(i)?;

    Ok((i, (rows, cols)))
}

fn parse_name(i: &[u8]) -> IResult<&[u8], &[u8]> {
    let (i, _) = tag("NAME")(i)?;
    let (i, _whitespace) = space1(i)?;
    take_till1(is_newline)(i)
}

fn parse_rows(i: &[u8]) -> IResult<&[u8], Rows> {
    let (i, _) = tag("ROWS")(i)?;
    let (i, _whitespace) = multispace1(i)?;

    fold_many1(parse_row, HashMap::new(), |mut rows, (name, row)| {
        let already_present = rows.insert(name, row).is_some();
        assert!(!already_present);
        rows
    })(i)
}

fn parse_row(i: &[u8]) -> IResult<&[u8], (&[u8], Row)> {
    let (i, _whitespace) = space0(i)?;
    let (i, row_type) = satisfy(|c| c == 'L' || c == 'G' || c == 'E' || c == 'N')(i)?;

    let new_constraint = |op| Row::Constraint {
        op,
        coeffs: HashMap::new(),
        rhs: None,
    };

    let row = match row_type {
        'L' => new_constraint(ConstraintOp::Lte),
        'G' => new_constraint(ConstraintOp::Gte),
        'E' => new_constraint(ConstraintOp::Eq),
        'N' => Row::Objective,
        _ => panic!(
            "internal mps parser error, unexpected row type: {}",
            row_type
        ),
    };

    let (i, _whitespace) = space0(i)?;
    let (i, row_name) = take_till1(is_newline)(i)?;
    let (i, _newline) = char('\n')(i)?;

    Ok((i, (row_name, row)))
}

//returns map from variable names to objective coefficients
fn parse_columns<'a>(i: &'a [u8], rows: &mut Rows<'a>) -> IResult<&'a [u8], Cols<'a>> {
    let (i, _) = tag("COLUMNS")(i)?;
    let (mut i, _whitespace) = multispace1(i)?;

    let mut cols = HashMap::new();

    while let Ok((j, ())) = parse_column_line(i, rows, &mut cols) {
        i = j
    }

    Ok((i, cols))
}

fn parse_column_line<'a>(
    i: &'a [u8],
    rows: &mut Rows<'a>,
    cols: &mut Cols<'a>,
) -> IResult<&'a [u8], ()> {
    let (i, var_name) = take_while1(is_alphanumeric)(i)?;

    let (i, _ws) = space1(i)?;
    let (i, col_info) = separated_list0(space1, parse_column)(i)?;

    for (row_name, coeff) in col_info {
        match rows.get_mut(row_name) {
            Some(row) => match row {
                Row::Objective => {
                    if cols
                        .insert(
                            var_name,
                            Col {
                                obj_coeff: Some(coeff),
                                bound: None,
                            },
                        )
                        .is_some()
                    {
                        error!(
                            "specified objective coefficient for {} more than once",
                            std::str::from_utf8(var_name).unwrap()
                        );

                        return Err(nom::Err::Failure(nom::error::Error::new(
                            i,
                            nom::error::ErrorKind::Tag,
                        )));
                    }
                }

                Row::Constraint { coeffs, .. } => {
                    if coeffs.insert(var_name, coeff).is_some() {
                        error!(
                            "specified constraint coefficient for {} more than once",
                            std::str::from_utf8(var_name).unwrap()
                        );

                        return Err(nom::Err::Failure(nom::error::Error::new(
                            i,
                            nom::error::ErrorKind::Tag,
                        )));
                    }
                }
            },

            None => {
                return Err(nom::Err::Failure(nom::error::Error::new(
                    i,
                    nom::error::ErrorKind::Tag,
                )))
            }
        }
    }

    let (i, _ws) = multispace1(i)?;
    Ok((i, ()))
}

fn parse_column(i: &[u8]) -> IResult<&[u8], (&[u8], f64)> {
    let (i, col_name) = take_while1(is_alphanumeric)(i)?;
    let (i, _whitespace) = space1(i)?;
    let (i, coeff) = double(i)?;
    Ok((i, (col_name, coeff)))
}

fn parse_all_rhs<'a>(i: &'a [u8], rows: &mut Rows) -> IResult<&'a [u8], ()> {
    let (i, _) = tag("RHS")(i)?;
    let (mut i, _whitespace) = multispace1(i)?;

    while let Ok((j, ())) = parse_rhs_line(i, rows) {
        i = j
    }

    Ok((i, ()))
}

fn parse_rhs_line<'a>(i: &'a [u8], rows: &mut Rows) -> IResult<&'a [u8], ()> {
    if i.starts_with(b"BOUNDS") {
        return Err(nom::Err::Failure(nom::error::Error::new(
            i,
            nom::error::ErrorKind::Tag,
        )));
    }

    let (i, _rhs_name) = take_while1(is_alphanumeric)(i)?;
    let (i, _ws) = space0(i)?;
    let (i, rhs_info) = separated_list0(space1, parse_rhs)(i)?;

    for (row_name, rhs_val) in rhs_info {
        match rows.get_mut(row_name) {
            Some(row) => match row {
                Row::Objective => {
                    error!("should not specific rhs value for the objective",);

                    return Err(nom::Err::Failure(nom::error::Error::new(
                        i,
                        nom::error::ErrorKind::Tag,
                    )));
                }

                Row::Constraint { rhs, .. } => {
                    if rhs.is_some() {
                        error!(
                            "specified rhs for {} more than once",
                            std::str::from_utf8(row_name).unwrap()
                        );

                        return Err(nom::Err::Failure(nom::error::Error::new(
                            i,
                            nom::error::ErrorKind::Tag,
                        )));
                    }

                    *rhs = Some(rhs_val);
                }
            },

            None => {
                return Err(nom::Err::Failure(nom::error::Error::new(
                    i,
                    nom::error::ErrorKind::Tag,
                )))
            }
        }
    }

    let (i, _ws) = multispace1(i)?;
    Ok((i, ()))
}

fn parse_rhs(i: &[u8]) -> IResult<&[u8], (&[u8], f64)> {
    let (i, row_name) = take_while1(is_alphanumeric)(i)?;
    let (i, _ws) = space1(i)?;
    let (i, coeff) = double(i)?;
    Ok((i, (row_name, coeff)))
}

fn parse_bounds<'a>(i: &'a [u8], cols: &mut Cols) -> IResult<&'a [u8], ()> {
    let (i, _) = tag("BOUNDS")(i)?;
    let (i, _ws) = multispace1(i)?;
    let (i, bound_info) = separated_list0(multispace1, parse_bound)(i)?;

    for (col_name, bound) in bound_info {
        match cols.get_mut(col_name) {
            Some(col) => match col.bound {
                Some(prev_bound) => match prev_bound {
                    Bound::Lower(lb) => {
                        if let Bound::Upper(ub) = bound {
                            col.bound.replace(Bound::TwoSided(lb, ub));
                        } else {
                            error!(
                                "invalid bounds for {}",
                                std::str::from_utf8(col_name).unwrap()
                            );

                            return Err(nom::Err::Failure(nom::error::Error::new(
                                i,
                                nom::error::ErrorKind::Tag,
                            )));
                        }
                    }

                    Bound::Upper(ub) => {
                        if let Bound::Lower(lb) = bound {
                            col.bound.replace(Bound::TwoSided(lb, ub));
                        } else {
                            error!(
                                "invalid bounds for {}",
                                std::str::from_utf8(col_name).unwrap()
                            );

                            return Err(nom::Err::Failure(nom::error::Error::new(
                                i,
                                nom::error::ErrorKind::Tag,
                            )));
                        }
                    }

                    _ => {
                        error!(
                            "invalid bounds for {}",
                            std::str::from_utf8(col_name).unwrap()
                        );

                        return Err(nom::Err::Failure(nom::error::Error::new(
                            i,
                            nom::error::ErrorKind::Tag,
                        )));
                    }
                },

                None => col.bound = Some(bound),
            },

            None => {
                error!(
                    "found bound for {}, but that column does not exist",
                    std::str::from_utf8(col_name).unwrap()
                );

                return Err(nom::Err::Failure(nom::error::Error::new(
                    i,
                    nom::error::ErrorKind::Tag,
                )));
            }
        }
    }

    Ok((i, ()))
}

fn parse_bound(i: &[u8]) -> IResult<&[u8], (&[u8], Bound)> {
    let (i, bound_type) = alt((tag("UP"), tag("LO"), tag("FR")))(i)?;

    let (i, _ws) = space1(i)?;
    let (i, _bound_name) = take_while1(is_alphanumeric)(i)?;

    let (i, _ws) = space1(i)?;
    let (mut i, col_name) = take_while1(is_alphanumeric)(i)?;

    let bound = match bound_type {
        b"UP" => {
            let (j, _ws) = space1(i)?;
            let res = double(j)?;
            i = res.0;
            Bound::Upper(res.1)
        }

        b"LO" => {
            let (j, _ws) = space1(i)?;
            let res = double(j)?;
            i = res.0;
            Bound::Lower(res.1)
        }

        b"FR" => Bound::Free,

        _ => panic!(
            "internal mps parser error, unexpected bound type: {}",
            std::str::from_utf8(bound_type).unwrap()
        ),
    };

    Ok((i, (col_name, bound)))
}

#[derive(Debug, Clone)]
enum Row<'a> {
    Objective, //keys are variable names

    Constraint {
        op: ConstraintOp,
        coeffs: HashMap<&'a [u8], f64>,
        rhs: Option<f64>,
    },
}

#[derive(Debug, Clone)]
struct Col {
    obj_coeff: Option<f64>,
    bound: Option<Bound>,
}

#[cfg(test)]
mod tests {
    use std::f64::EPSILON;

    use super::*;

    #[test]
    fn parse_example() {
        let mps_str = r#"
            NAME          TESTPROB
            ROWS
            N  COST
            L  LIM1
            G  LIM2
            E  MYEQN
            COLUMNS
                XONE      COST                 1   LIM1                 1
                XONE      LIM2                 1
                YTWO      COST                 4   LIM1                 1
                YTWO      MYEQN               -1
                ZTHREE    COST                 9   LIM2                 1
                ZTHREE    MYEQN                1
            RHS
                RHS1      LIM1                 5   LIM2                10
                RHS1      MYEQN                7
            BOUNDS
            UP BND1      XONE                 4
            LO BND1      YTWO                -1
            UP BND1      YTWO                 1
            ENDATA
        "#;

        let prob: Problem = parse_mps(mps_str.as_bytes()).unwrap();

        for var in &prob.variables {
            match var.name.as_ref().unwrap().as_str() {
                "XONE" => {
                    assert!((var.obj_coeff - 1.) < EPSILON);
                    assert_eq!(var.bound, Bound::Upper(4.));
                }

                "YTWO" => {
                    assert!((var.obj_coeff - 4.) < EPSILON);
                    assert_eq!(var.bound, Bound::TwoSided(-1., 1.));
                }

                "ZTHREE" => {
                    assert!((var.obj_coeff - 9.) < EPSILON);
                    assert_eq!(var.bound, Bound::Free);
                }

                _ => panic!("unexpected variable: {}", var.name.as_ref().unwrap()),
            }
        }

        for constraint in &prob.constraints {
            match constraint.op {
                ConstraintOp::Lte => {
                    assert!((constraint.rhs - 5.) < EPSILON);
                    assert_eq!(constraint.coeffs.len(), 2);
                }

                ConstraintOp::Gte => {
                    assert!((constraint.rhs - 10.) < EPSILON);
                    assert_eq!(constraint.coeffs.len(), 2);
                }

                ConstraintOp::Eq => {
                    assert!((constraint.rhs - 7.) < EPSILON);
                    assert_eq!(constraint.coeffs.len(), 2);
                }
            }
        }
    }
}
