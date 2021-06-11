use crate::error::EllPError;
use crate::util::EPS;

use std::collections::{HashMap, HashSet};

const LTE_STR: &str = "\u{2264}";
const EQ_STR: &str = "\u{003D}";
const GTE_STR: &str = "\u{2265}";
const INF_STR: &str = "\u{221E}";

#[derive(Debug, Clone, Default)]
pub struct Problem {
    pub variables: Vec<Variable>,
    pub constraints: Vec<Constraint>,
    var_names: HashSet<String>, //these strings are duplicated in the variables
    var_ids: HashSet<VariableId>,
}

impl Problem {
    pub fn new() -> Self {
        Default::default()
    }

    pub fn add_var(
        &mut self,
        obj_coeff: f64,
        bound: Bound,
        name: Option<String>,
    ) -> Result<VariableId, EllPError> {
        let id = VariableId(self.variables.len());
        self.add_var_with_id(obj_coeff, bound, id, name)?;
        Ok(id)
    }

    pub fn add_var_with_id(
        &mut self,
        obj_coeff: f64,
        bound: Bound,
        id: VariableId,
        name: Option<String>,
    ) -> Result<VariableId, EllPError> {
        if let Bound::TwoSided(lb, ub) = bound {
            if lb > ub {
                return Err(EllPError::new(format!(
                    "invalid variable bounds: ({}, {})",
                    lb, ub
                )));
            }
        }

        let bound_valid = match bound {
            Bound::Free => true,
            Bound::Lower(lb) => lb.is_finite(),
            Bound::Upper(ub) => ub.is_finite(),
            Bound::TwoSided(lb, ub) => lb.is_finite() && ub.is_finite(),
            Bound::Fixed(fixed_val) => fixed_val.is_finite(),
        };

        if !bound_valid {
            return Err(EllPError::new(format!("invalid bound: {:?}", bound)));
        }

        if let Some(name) = &name {
            if !self.var_names.insert(name.clone()) {
                return Err(EllPError::new(format!(
                    "variable names must be unique, {} was added twice",
                    name
                )));
            }
        }

        let var = Variable::new(id, obj_coeff, bound, name);
        self.variables.push(var);

        if !self.var_ids.insert(id) {
            return Err(EllPError::new(format!(
                "cannot add variable with {:?}, that id is already used",
                id
            )));
        }

        Ok(id)
    }

    pub fn add_constraint(
        &mut self,
        coeffs: Vec<(VariableId, f64)>,
        op: ConstraintOp,
        rhs: f64,
    ) -> Result<(), EllPError> {
        match coeffs
            .iter()
            .find(|(id, _coeff)| !self.var_ids.contains(id))
        {
            Some((invalid_var, _coeff)) => {
                Err(EllPError::new(format!("{:?} is invalid", invalid_var)))
            }

            None => {
                self.constraints.push(Constraint { coeffs, op, rhs });
                Ok(())
            }
        }
    }

    pub fn constraints(&self) -> &[Constraint] {
        self.constraints.as_slice()
    }

    pub fn is_feasible(&self, x: &[f64]) -> bool {
        if x.len() != self.variables.len() {
            return false;
        }

        for (var, &val) in self.variables.iter().zip(x.iter()) {
            match var.bound {
                Bound::Free => (),

                Bound::Lower(lb) => {
                    if val < lb - EPS {
                        return false;
                    }
                }
                Bound::Upper(ub) => {
                    if val > ub + EPS {
                        return false;
                    }
                }

                Bound::TwoSided(lb, ub) => {
                    if val < lb - EPS {
                        return false;
                    }

                    if val > ub + EPS {
                        return false;
                    }
                }

                Bound::Fixed(fixed_val) => {
                    if (val - fixed_val).abs() > EPS {
                        return false;
                    }
                }
            }
        }

        for constraint in &self.constraints {
            if !constraint.is_feasible(x) {
                return false;
            }
        }

        true
    }
}

#[derive(Debug, Clone)]
pub struct Variable {
    pub id: VariableId,
    pub obj_coeff: f64,
    pub bound: Bound,
    pub name: Option<String>,
}

impl Variable {
    fn new(id: VariableId, obj_coeff: f64, bound: Bound, name: Option<String>) -> Self {
        Self {
            id,
            obj_coeff,
            bound,
            name,
        }
    }
}

impl PartialEq for Variable {
    fn eq(&self, other: &Self) -> bool {
        //Problem add_var guarantees that different variables have different ids
        self.id == other.id
    }
}

impl Eq for Variable {}

impl std::hash::Hash for Variable {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.id.hash(state);
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Bound {
    Free,
    Lower(f64),
    Upper(f64),
    TwoSided(f64, f64),
    Fixed(f64),
}

impl Bound {
    fn display(&self, f: &mut std::fmt::Formatter, var: &Variable) -> std::fmt::Result {
        match self {
            Bound::Free => write!(f, "{} free", var),
            Bound::Lower(lb) => write!(f, "{} {gte} {}", var, lb, gte = GTE_STR),
            Bound::Upper(ub) => write!(f, "{} {lte} {}", var, ub, lte = LTE_STR),
            Bound::TwoSided(lb, ub) => {
                write!(f, "{} {lte} {} {lte} {}", lb, var, ub, lte = LTE_STR)
            }
            Bound::Fixed(val) => write!(f, "{} {eq} {}", var, val, eq = EQ_STR),
        }
    }
}

impl std::fmt::Display for Bound {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Bound::Free => write!(f, "(-{inf}, {inf})", inf = INF_STR),
            Bound::Lower(lb) => write!(f, "[{}, {inf})", lb, inf = INF_STR),
            Bound::Upper(ub) => write!(f, "(-{inf}, {}]", ub, inf = INF_STR),
            Bound::TwoSided(lb, ub) => write!(f, "[{}, {}]", lb, ub),
            Bound::Fixed(val) => write!(f, "[{val}, {val}]", val = val),
        }
    }
}

#[derive(Debug, Clone)]
pub struct Constraint {
    pub coeffs: Vec<(VariableId, f64)>,
    pub op: ConstraintOp,
    pub rhs: f64,
}

impl Constraint {
    pub fn add_coeff(&mut self, var: VariableId, coeff: f64) {
        self.coeffs.push((var, coeff));
    }

    fn is_feasible(&self, x: &[f64]) -> bool {
        let mut lhs = 0.;

        for (var, coeff) in &self.coeffs {
            let i: usize = var.into();
            lhs += coeff * x[i];
        }

        match self.op {
            ConstraintOp::Lte => lhs <= self.rhs + EPS,
            ConstraintOp::Eq => (lhs - self.rhs).abs() < EPS,
            ConstraintOp::Gte => lhs >= self.rhs - EPS,
        }
    }

    fn display(
        &self,
        f: &mut std::fmt::Formatter,
        var_names: &HashMap<VariableId, &Variable>,
    ) -> std::fmt::Result {
        for (var_id, coeff) in &self.coeffs {
            if *coeff == 0. {
                continue;
            }

            let var = *var_names.get(&var_id).unwrap();

            write!(
                f,
                "{:+} {} {} ",
                if *coeff >= 0. { "+" } else { "-" },
                coeff.abs(),
                var
            )?;
        }

        write!(f, "{} {}", self.op, self.rhs)
    }
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub struct VariableId(usize);

impl std::convert::From<usize> for VariableId {
    fn from(id: usize) -> Self {
        Self(id)
    }
}

impl std::convert::From<VariableId> for usize {
    fn from(id: VariableId) -> Self {
        id.0
    }
}

impl std::convert::From<&VariableId> for usize {
    fn from(id: &VariableId) -> Self {
        id.0
    }
}

#[derive(Debug, Clone)]
pub enum ConstraintOp {
    Lte,
    Eq,
    Gte,
}

impl std::fmt::Display for Problem {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        writeln!(f, "minimize")?;
        let mut var_id_to_var: HashMap<VariableId, &Variable> = HashMap::new();

        for var in &self.variables {
            let result = var_id_to_var.insert(var.id, &var);
            assert!(result.is_none()); //should have have repeated ids

            if var.obj_coeff == 0. {
                continue;
            }

            write!(
                f,
                "{} {} {} ",
                if var.obj_coeff > 0. { "+" } else { "-" },
                var.obj_coeff.abs(),
                var
            )?;
        }

        writeln!(f, "\n\nsubject to")?;

        for constraint in &self.constraints {
            constraint.display(f, &var_id_to_var)?;
            writeln!(f)?;
        }

        writeln!(f, "\nwith the bounds")?;

        for var in &self.variables {
            var.bound.display(f, var)?;
            writeln!(f)?;
        }

        Ok(())
    }
}

impl std::fmt::Display for Variable {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match &self.name {
            Some(name) => write!(f, "{}", name),
            None => write!(f, "id[{}]", self.id.0),
        }
    }
}

impl std::fmt::Display for ConstraintOp {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            ConstraintOp::Lte => write!(f, "{}", LTE_STR),
            ConstraintOp::Eq => write!(f, "{}", EQ_STR),
            ConstraintOp::Gte => write!(f, "{}", GTE_STR),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::ConstraintOp::*;
    use super::{Bound, Problem, VariableId};

    #[test]
    fn add_var() {
        let mut prob = Problem::new();
        prob.add_var(1., Bound::Free, Some("x".to_string()))
            .unwrap();
        assert_eq!(prob.variables[0].name.as_ref().unwrap(), "x");
    }

    #[test]
    fn add_var_with_bounds() {
        let mut prob = Problem::new();
        prob.add_var(1., Bound::TwoSided(1., 2.), Some("x".to_string()))
            .unwrap();
        let var = &prob.variables[0];
        assert_eq!(var.bound, Bound::TwoSided(1., 2.));
    }

    #[test]
    fn add_var_bad_bounds() {
        let mut prob = Problem::new();
        assert!(prob
            .add_var(1., Bound::TwoSided(1., 0.), Some("x".to_string()))
            .is_err());
    }

    #[test]
    fn add_constraint() {
        let mut prob = Problem::new();
        let var_id = prob.add_var(1., Bound::Free, None).unwrap();
        assert!(prob.add_constraint(vec![(var_id, 1.)], Lte, 0.).is_ok());
    }

    #[test]
    fn add_constraint_invalid_var() {
        let mut prob = Problem::new();
        assert!(prob
            .add_constraint(vec![(VariableId(0), 1.)], Lte, 0.)
            .is_err());
    }

    #[test]
    fn nonunique_var_names() {
        let mut prob = Problem::new();
        prob.add_var(0., Bound::Free, Some("x".to_string()))
            .unwrap();

        let is_err = prob
            .add_var(0., Bound::Free, Some("x".to_string()))
            .is_err();

        assert!(is_err);
    }
}
