use thiserror::Error;

#[derive(Error, Debug)]
#[error("{0}")]
pub struct EllPError(String);

impl EllPError {
    pub fn new(msg: String) -> Self {
        Self(msg)
    }
}
