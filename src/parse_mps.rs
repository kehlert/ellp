use crate::problem::{Bound, ConstraintOp, Problem};

use log::{debug, error};
use thiserror::Error;

use std::{collections::HashMap, iter::Peekable};

type Rows<'a> = HashMap<&'a str, Row<'a>>;
type Cols<'a> = HashMap<&'a str, Col>;

#[derive(Error, Debug)]
#[error("MPS parsing error. {msg}")]
pub struct MpsParsingError {
    msg: String,
}

impl MpsParsingError {
    fn new(msg: String) -> Self {
        Self { msg }
    }
}

pub fn parse_mps(mps: &str) -> Result<Problem, MpsParsingError> {
    let (rows, cols) = parse_rows_and_cols(mps)?;

    let mut prob = Problem::new();
    let mut var_ids = HashMap::new();

    for (var_name, col) in cols {
        let var_name = var_name.to_string();
        let bound = col.bound.unwrap_or(Bound::Lower(0.));

        //should never panic, we know that the var names are unique
        let var_id = prob
            .add_var(col.obj_coeff, bound, Some(var_name.clone()))
            .unwrap();

        var_ids.insert(var_name, var_id);
    }

    for (row_name, row) in rows {
        match row {
            Row::Objective => continue,
            Row::Constraint { op, coeffs, rhs } => {
                let rhs = rhs.unwrap_or(0.);

                let coeffs = coeffs
                    .into_iter()
                    .map(|(var_name, coeff)| match var_ids.get(var_name) {
                        Some(var_id) => Ok((*var_id, coeff)),

                        None => Err(MpsParsingError::new(format!(
                            "for row {}, column {} does not exist",
                            row_name, var_name
                        ))),
                    })
                    .collect::<Result<_, _>>()?;

                //should never panic, we know the variable ids in coeffs exist
                prob.add_constraint(coeffs, op, rhs).unwrap();
            }
        }
    }

    Ok(prob)
}

fn parse_rows_and_cols(mps: &str) -> Result<(Rows, Cols), MpsParsingError> {
    let mut lines = mps
        .split_terminator('\n')
        .filter(|s| !s.trim().is_empty())
        .peekable();

    match lines.next() {
        Some(line) => parse_name(line)?,
        None => return Err(MpsParsingError::new("could not find NAME line".to_string())),
    }

    let mut rows = parse_rows(&mut lines)?;
    let mut cols = parse_columns(&mut lines, &mut rows)?;

    parse_rhs(&mut lines, &mut rows)?;
    parse_bounds(&mut lines, &mut cols)?;

    match lines.next() {
        Some(line) => {
            let line = line.trim();

            if line != "ENDATA" {
                return Err(MpsParsingError::new(format!(
                    "expected 'ENDATA', found '{}'",
                    line
                )));
            }
        }

        None => {
            return Err(MpsParsingError::new(
                "could not find ENDATA line".to_string(),
            ))
        }
    }

    if let Some(line) = lines.next() {
        return Err(MpsParsingError::new(format!("unexpected line: {}", line)));
    }

    debug!(
        "parsed {} rows and {} columns from mps file",
        rows.len(),
        cols.len()
    );

    Ok((rows, cols))
}

fn parse_name(line: &str) -> Result<(), MpsParsingError> {
    let mut split_line = line.split_whitespace();

    let name = split_line.next().and_then(|tag| {
        if tag == "NAME" {
            split_line.next()
        } else {
            None
        }
    });

    if name.is_none() {
        return Err(MpsParsingError::new(format!(
            "could not find name in NAME line: {}",
            line
        )));
    }

    //ignore everything after NAME
    Ok(())
}

fn parse_rows<'a, I>(lines: &mut Peekable<I>) -> Result<Rows<'a>, MpsParsingError>
where
    I: Iterator<Item = &'a str>,
{
    match lines.next() {
        Some(line) => {
            let line = line.trim();

            if line != "ROWS" {
                return Err(MpsParsingError::new(format!(
                    "expected 'ROWS', found '{}'",
                    line
                )));
            }
        }

        None => return Err(MpsParsingError::new("could not find ROWS line".to_string())),
    }

    let mut rows = HashMap::new();

    while let Some(&line) = lines.peek() {
        let line = line.trim_start();

        if line.starts_with("COLUMNS") {
            break;
        }

        lines.next();

        let (name, row) = parse_row_line(line)?;

        if rows.insert(name, row).is_some() {
            return Err(MpsParsingError::new(format!("row name repeated: {}", name)));
        }
    }

    Ok(rows)
}

fn parse_row_line(line: &str) -> Result<(&str, Row), MpsParsingError> {
    let mut split_line = line.split_whitespace();

    let constraint_op_char = split_line.next().ok_or_else(|| {
        MpsParsingError::new(format!(
            "expected a row type character in this line: {}",
            line
        ))
    })?;

    let new_constraint = |op| Row::Constraint {
        op,
        coeffs: HashMap::new(),
        rhs: None,
    };

    let row = match constraint_op_char {
        "L" => new_constraint(ConstraintOp::Lte),
        "G" => new_constraint(ConstraintOp::Gte),
        "E" => new_constraint(ConstraintOp::Eq),
        "N" => Row::Objective,
        _ => {
            return Err(MpsParsingError::new(format!(
                "unexpected row type: {}",
                constraint_op_char
            )))
        }
    };

    let row_name = split_line.next().ok_or_else(|| {
        MpsParsingError::new(format!("expected a row name in this line: {}", line))
    })?;

    match split_line.next() {
        Some(s) => Err(MpsParsingError::new(format!(
            "unexpected input in row line: {}",
            s
        ))),

        None => Ok((row_name, row)),
    }
}

fn parse_columns<'a, I>(
    lines: &mut Peekable<I>,
    rows: &mut Rows<'a>,
) -> Result<Cols<'a>, MpsParsingError>
where
    I: Iterator<Item = &'a str>,
{
    match lines.next() {
        Some(line) => {
            let line = line.trim();

            if line != "COLUMNS" {
                return Err(MpsParsingError::new(format!(
                    "expected 'COLUMNS', found '{}'",
                    line
                )));
            }
        }

        None => {
            return Err(MpsParsingError::new(
                "could not find COLUMNS line".to_string(),
            ))
        }
    }

    let mut cols = HashMap::new();

    while let Some(&line) = lines.peek() {
        let line = line.trim_start();

        if line.starts_with("RHS") {
            break;
        }

        lines.next();
        parse_column_line(line, rows, &mut cols)?;
    }

    Ok(cols)
}

fn parse_column_line<'a>(
    line: &'a str,
    rows: &mut Rows<'a>,
    cols: &mut Cols<'a>,
) -> Result<(), MpsParsingError> {
    let mut split_line = line.split_whitespace();

    let var_name = split_line.next().ok_or_else(|| {
        MpsParsingError::new(format!("expected a column name in this line: {}", line))
    })?;

    let row_name = split_line.next().ok_or_else(|| {
        MpsParsingError::new(format!("expected a row name in this line: {}", line))
    })?;

    let coeff = split_line.next().ok_or_else(|| {
        MpsParsingError::new(format!("expected a coefficient in this line: {}", line))
    })?;

    let coeff: f64 = coeff.parse().map_err(|err| {
        MpsParsingError::new(format!(
            "could not parse the coefficient {}\nerror: {}\nline: {}",
            coeff, err, line
        ))
    })?;

    if let Some(s) = split_line.next() {
        return Err(MpsParsingError::new(format!(
            "unexpected input '{}' in column line: {}",
            s, line
        )));
    }

    let col_entry = cols.entry(var_name).or_insert_with(|| Col {
        obj_coeff: 0.,
        bound: None,
    });

    match rows.get_mut(row_name) {
        Some(row) => match row {
            Row::Objective => col_entry.obj_coeff = coeff,
            Row::Constraint { coeffs, .. } => {
                if coeffs.insert(var_name, coeff).is_some() {
                    return Err(MpsParsingError::new(format!(
                        "specified constraint coefficient for the column {} and row {} more than once",
                        var_name, row_name
                    )));
                }
            }
        },

        None => {
            return Err(MpsParsingError::new(format!(
                "could not find the row {}",
                row_name
            )))
        }
    }

    match split_line.next() {
        Some(s) => Err(MpsParsingError::new(format!(
            "unexpected input in column line: {}",
            s
        ))),

        None => Ok(()),
    }
}

fn parse_rhs<'a, I>(lines: &mut Peekable<I>, rows: &mut Rows) -> Result<(), MpsParsingError>
where
    I: Iterator<Item = &'a str>,
{
    match lines.next() {
        Some(line) => {
            let line = line.trim();

            if line != "RHS" {
                return Err(MpsParsingError::new(format!(
                    "expected 'RHS', found '{}'",
                    line
                )));
            }
        }

        None => return Err(MpsParsingError::new("could not find RHS line".to_string())),
    }

    while let Some(&line) = lines.peek() {
        let line = line.trim_start();

        if line.starts_with("BOUNDS") || line.starts_with("ENDATA") {
            break;
        }

        lines.next();
        parse_rhs_line(line, rows)?;
    }

    Ok(())
}

fn parse_rhs_line(line: &str, rows: &mut Rows) -> Result<(), MpsParsingError> {
    let split_line_iter = line.split_whitespace();
    let mut split_line = split_line_iter.clone();

    if split_line_iter.count() == 3 {
        split_line.next(); //skip the name
    }

    let row_name = split_line.next().ok_or_else(|| {
        MpsParsingError::new(format!("expected a row name in this line: {}", line))
    })?;

    let rhs_val = split_line.next().ok_or_else(|| {
        MpsParsingError::new(format!("expected a rhs value in this line: {}", line))
    })?;

    let rhs_val: f64 = rhs_val.parse().map_err(|err| {
        MpsParsingError::new(format!(
            "could not parse the rhs value {}\nerror: {}\nline: {}",
            rhs_val, err, line
        ))
    })?;

    match rows.get_mut(row_name) {
        Some(row) => match row {
            Row::Objective => {
                return Err(MpsParsingError::new(
                    "should not specify rhs value for the objective".to_string(),
                ));
            }

            Row::Constraint { rhs, .. } => {
                if rhs.is_some() {
                    return Err(MpsParsingError::new(format!(
                        "specified rhs for {} more than once",
                        row_name
                    )));
                }

                *rhs = Some(rhs_val);
            }
        },

        None => {
            return Err(MpsParsingError::new(format!(
                "could not find the row {}",
                row_name
            )))
        }
    }

    match split_line.next() {
        Some(s) => Err(MpsParsingError::new(format!(
            "unexpected input in column line: {}",
            s
        ))),

        None => Ok(()),
    }
}

fn parse_bounds<'a, I>(lines: &mut Peekable<I>, cols: &mut Cols) -> Result<(), MpsParsingError>
where
    I: Iterator<Item = &'a str>,
{
    if let Some(line) = lines.peek() {
        if *line == "ENDATA" {
            return Ok(()); //MPS file does not have a BOUNDS section
        }
    }

    match lines.next() {
        Some(line) => {
            let line = line.trim();

            if line == "ENDATA" {
                return Ok(());
            } else if line != "BOUNDS" {
                return Err(MpsParsingError::new(format!(
                    "expected 'BOUNDS', found '{}'",
                    line
                )));
            }
        }

        None => {
            return Err(MpsParsingError::new(
                "could not find BOUNDS line".to_string(),
            ))
        }
    }

    while let Some(&line) = lines.peek() {
        let line = line.trim_start();

        if line.starts_with("ENDATA") {
            break;
        }

        lines.next();
        parse_bound_line(line, cols)?;
    }

    Ok(())
}

fn parse_bound_line(line: &str, cols: &mut Cols) -> Result<(), MpsParsingError> {
    let mut split_line = line.split_whitespace();

    let bound_type = split_line.next().ok_or_else(|| {
        MpsParsingError::new(format!("expected a bound type in this line: {}", line))
    })?;

    split_line.next(); //skip bound name

    let col_name = split_line.next().ok_or_else(|| {
        MpsParsingError::new(format!("expected a column name in this line: {}", line))
    })?;

    let bound_val = split_line
        .next()
        .map(str::parse::<f64>)
        .transpose()
        .map_err(|err| {
            MpsParsingError::new(format!(
                "could not parse the bound value\nerror: {}\nline: {}",
                err, line
            ))
        })?;

    let bound = match (bound_type, bound_val) {
        ("UP", Some(v)) => Bound::Upper(v),
        ("LO", Some(v)) => Bound::Lower(v),
        ("FR", None) => Bound::Free,
        _ => {
            return Err(MpsParsingError::new(format!(
                "invalid bound specification: {}",
                line
            )))
        }
    };

    let col = match cols.get_mut(col_name) {
        Some(col) => col,

        None => {
            return Err(MpsParsingError::new(format!(
                "found bound for the column {}, but it does not exist",
                col_name,
            )))
        }
    };

    match (col.bound, bound) {
        (Some(Bound::Upper(ub)), Bound::Lower(lb)) => {
            col.bound.replace(Bound::TwoSided(lb, ub));
        }

        (Some(Bound::Lower(lb)), Bound::Upper(ub)) => {
            col.bound.replace(Bound::TwoSided(lb, ub));
        }

        (None, _) => col.bound = Some(bound),

        _ => {
            return Err(MpsParsingError::new(format!(
                "invalid bounds for {}",
                col_name,
            )))
        }
    }

    match split_line.next() {
        Some(s) => Err(MpsParsingError::new(format!(
            "unexpected input in column line: {}",
            s
        ))),

        None => Ok(()),
    }
}

#[derive(Debug, Clone)]
enum Row<'a> {
    Objective,

    Constraint {
        op: ConstraintOp,
        coeffs: HashMap<&'a str, f64>, //keys are variable names
        rhs: Option<f64>,
    },
}

#[derive(Debug, Clone)]
struct Col {
    obj_coeff: f64,
    bound: Option<Bound>,
}

#[cfg(test)]
mod tests {
    use std::f64::EPSILON;

    use super::*;

    #[test]
    fn parse_mps_example() {
        let mps_str = r#"
            NAME          TESTPROB
            ROWS
            N  COST
            L  LIM1
            G  LIM2
            E  MYEQN
            COLUMNS
                XONE      COST                 1
                XONE      LIM1                 1
                XONE      LIM2                 1
                YTWO      COST                 4
                YTWO      LIM1                 1
                YTWO      MYEQN               -1
                ZTHREE    COST                 9
                ZTHREE    LIM2                 1
                ZTHREE    MYEQN                1
            RHS
                RHS1      LIM1                 5
                RHS1      LIM2                10
                RHS1      MYEQN                7
            BOUNDS
            UP BND1      XONE                 4
            LO BND1      YTWO                -1
            UP BND1      YTWO                 1
            ENDATA
        "#;

        let prob: Problem = parse_mps(mps_str).unwrap();

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
                    assert_eq!(var.bound, Bound::Lower(0.));
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
