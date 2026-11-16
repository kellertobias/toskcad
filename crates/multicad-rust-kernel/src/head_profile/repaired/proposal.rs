//! Bounded correspondence to the original proposal, never source identity.
use super::*;
#[cfg(test)]
mod tests;
/// Whole-meridian comparison under a per-arc monotone angular correspondence.
/// Lines use identical parameters. Arc correspondence maps the repaired base
/// circle's angular fraction to the original proposal's uniform-angle fraction.
/// The correction remains part of the repaired geometry. Derivatives are with
/// respect to repaired normalized-ray t, including the correspondence derivative.
#[derive(Clone, Copy, Debug)]
pub struct RepairedProposalTransfer {
    pub position_mm: f64,
    pub derivative_mm: f64,
    pub normal_angle_radians: f64,
    pub work_used: usize,
}
// Proposal angles are cached atan2 constants in[-pi,pi], plus a source
// sweep of magnitude<3. No reduction or transcendental library allowance is
// needed on this complete domain. Taylor degree63 has Lagrange remainder
// |x|^64/64! because every real sine/cosine derivative has magnitude<=1.
pub(crate) fn rotation(angle: I) -> Result<[I; 2], RepairedHeadError> {
    if !angle.valid() || angle.magnitude() > 7. {
        return Err(RepairedHeadError::Unresolved);
    }
    let mut term = [I::point(1.), I::point(0.)];
    let mut value = term;
    for k in 1..=63 {
        term = [
            -term[1] * angle / I::point(k as f64),
            term[0] * angle / I::point(k as f64),
        ];
        for axis in 0..2 {
            value[axis] = value[axis] + term[axis];
        }
    }
    let mut remainder = I::point(1.);
    for k in 1..=64 {
        remainder = remainder * I::point(angle.magnitude()) / I::point(k as f64);
    }
    Ok(value.map(|v| {
        I::new(
            (v - I::point(remainder.hi)).lo.max(-1.),
            (v + I::point(remainder.hi)).hi.min(1.),
        )
    }))
}
impl RepairedHeadMeridian<'_> {
    pub fn proposal_transfer(
        &self,
        maximum_work: usize,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<RepairedProposalTransfer, RepairedHeadError> {
        poll(interrupt)?;
        if maximum_work > WORK {
            return Err(RepairedHeadError::Invalid);
        }
        let required = 1 + self.count * 4096;
        if maximum_work < required {
            return Err(RepairedHeadError::Limit);
        }
        let mut result = RepairedProposalTransfer {
            position_mm: 0.,
            derivative_mm: 0.,
            normal_angle_radians: 0.,
            work_used: required,
        };
        for part in &self.parts[..self.count] {
            poll(interrupt)?;
            if let Part::Arc(arc) = part.expect("constructor retained every qualified part") {
                let bound = arc.proposal_transfer()?;
                result.position_mm = result.position_mm.max(bound[0]);
                result.derivative_mm = result.derivative_mm.max(bound[1]);
                result.normal_angle_radians = result.normal_angle_radians.max(bound[2]);
            }
        }
        poll(interrupt)?;
        Ok(result)
    }
}
