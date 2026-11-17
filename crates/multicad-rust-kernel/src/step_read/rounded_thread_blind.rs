//! Actual nine-face blind graph correspondence to an unchanged source owner.
use super::*;
use crate::{rounded_thread::BlindBore, tapered_thread::interval::I};
fn raw_axis_class(shell: &Shell<'_, '_>) -> Result<()> {
    // Narrow source-known contract: normal-scale ratio vectors and nearly
    // orthogonal AXIS2 normal/reference directions. This bounds REAL->cache
    // normalization conditioning; arbitrary projected or tiny-scale STEP
    // frames are not granted an unearned EPS*extent transfer certificate.
    for parts in shell.entities.values() {
        for entity in parts {
            if entity.name == "DIRECTION" {
                let [Value::String(_), Value::List(values)] = entity.args.as_slice() else {
                    return Err(ReadError::Invalid);
                };
                if values.len() != 3 {
                    return Err(ReadError::Unsupported);
                }
                let mut scale = 0_f64;
                for value in values {
                    let Value::Number(value) = value else {
                        return Err(ReadError::Invalid);
                    };
                    if !value.is_finite() {
                        return Err(ReadError::Invalid);
                    }
                    scale = scale.max(value.abs());
                }
                if scale < f64::MIN_POSITIVE {
                    return Err(ReadError::Unsupported);
                }
            }
        }
    }
    for parts in shell.entities.values() {
        for entity in parts {
            if entity.name == "AXIS2_PLACEMENT_3D" {
                let [Value::String(_), _, normal, reference_dir] = entity.args.as_slice() else {
                    return Err(ReadError::Invalid);
                };
                let n = shell.direction(reference(normal)?)?;
                let u = shell.direction(reference(reference_dir)?)?;
                let cosine = n
                    .into_iter()
                    .zip(u)
                    .fold(I::point(0.), |sum, (n, u)| sum + I::point(n) * I::point(u));
                if !cosine.valid() || cosine.magnitude() > 1e-8 {
                    return Err(ReadError::Unsupported);
                }
            }
        }
    }
    Ok(())
}
/// Verify one original-source blind member. Raw analytic frames must have a
/// normal-scale ratio vector and quasiorthogonal normal/reference axes; tiny
/// vector scales and ill-conditioned projected axes are explicitly unsupported.
pub fn verify_rounded_thread_blind_bore(
    input: &str,
    bore: &BlindBore,
    origin: Point,
    u: Point,
    v: Point,
    tolerance: f64,
) -> Result<f64> {
    crate::step::rounded_thread_blind::precision(bore, origin, tolerance).map_err(|e| match e {
        crate::step::StepError::Limit => ReadError::Limit,
        _ => ReadError::Invalid,
    })?;
    let mut members = compound::read_members_by(input, tolerance, true, true, |shell, id, _| {
        raw_axis_class(shell)?;
        rounded_thread::verify_blind_member(shell, id, bore, [origin, u, v])
    })?;
    members.pop().ok_or(ReadError::Invalid)
}

#[cfg(test)]
mod tests;

mod repaired;
pub use repaired::{
    RepairedRoundedBlindCandidate, read_repaired_rounded_blind_bore,
    repaired_rounded_blind_import_matches,
};

pub(super) use repaired::recognize_member;

mod correspondence;
pub use correspondence::{RoundedBlindCorrespondence, verify_rounded_blind_correspondence};
