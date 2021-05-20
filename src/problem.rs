use crate::error::EllPError;
use crate::util::EPS;

#[derive(Debug, Clone, Default)]
pub struct Problem {
    variables: Vec<Variable>,
    constraints: Vec<Constraint>,
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

        //TODO check that name is unique if it's provided

        let var = Variable::new(VariableId(self.variables.len()), obj_coeff, bound, name);
        self.variables.push(var);
        let var = self.variables.last().unwrap();
        Ok(var.id)
    }

    pub fn var(&self, id: VariableId) -> Option<&Variable> {
        self.variables.get(id.0)
    }

    pub fn vars(&self) -> &[Variable] {
        self.variables.as_slice()
    }

    pub fn add_constraint(
        &mut self,
        coeffs: Vec<(VariableId, f64)>,
        op: ConstraintOp,
        rhs: f64,
    ) -> Result<(), EllPError> {
        match coeffs
            .iter()
            .find(|(id, _coeff)| id.0 >= self.variables.len())
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

    pub fn is_feasible(&self, x: &nalgebra::DVector<f64>) -> bool {
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

#[derive(Debug, Clone)]
pub struct Constraint {
    pub coeffs: Vec<(VariableId, f64)>,
    pub op: ConstraintOp,
    pub rhs: f64,
}

impl Constraint {
    fn is_feasible(&self, x: &nalgebra::DVector<f64>) -> bool {
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
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub struct VariableId(usize);

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

#[cfg(test)]
mod tests {
    use super::ConstraintOp::*;
    use super::{Bound, Problem, VariableId};

    #[test]
    fn add_var() {
        let mut prob = Problem::new();
        let var_id = prob
            .add_var(1., Bound::Free, Some("x".to_string()))
            .unwrap();
        assert_eq!(prob.var(var_id).unwrap().name.as_ref().unwrap(), "x");
    }

    #[test]
    fn add_var_with_bounds() {
        let mut prob = Problem::new();
        let var_id = prob
            .add_var(1., Bound::TwoSided(1., 2.), Some("x".to_string()))
            .unwrap();
        let var = prob.var(var_id).unwrap();
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
}
