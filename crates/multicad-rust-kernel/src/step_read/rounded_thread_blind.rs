//! Actual nine-face blind graph correspondence to an unchanged source owner.
use super::*;
use crate::rounded_thread::BlindBore;
pub fn verify_rounded_thread_blind_bore(
    input: &str,
    bore: &BlindBore,
    origin: Point,
    u: Point,
    v: Point,
    tolerance: f64,
) -> Result<f64> {
    crate::step::rounded_thread_blind::dimensions(bore, tolerance).map_err(|e| match e {
        crate::step::StepError::Limit => ReadError::Limit,
        _ => ReadError::Invalid,
    })?;
    let mut members = compound::read_members_by(input, tolerance, true, true, |shell, id, _| {
        rounded_thread::verify_blind_member(shell, id, bore, [origin, u, v])
    })?;
    members.pop().ok_or(ReadError::Invalid)
}

#[cfg(test)]
mod tests;
