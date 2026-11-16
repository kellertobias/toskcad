//! A genuine closed one-winding literal radial graph, not an inferred repeat
//! of separately fitted long-domain curves. No shell/material certificate.
use super::*;
pub struct PeriodicAngular {
    left: bool,
    curve: AngularCarrier,
    jacobian_lower: f64,
    work: usize,
    bytes: usize,
}
pub(super) fn radial_proof(position_error: f64, derivative_error: f64) -> Result<f64, Error> {
    let tau = I::new(
        std::f64::consts::TAU.next_down(),
        std::f64::consts::TAU.next_up(),
    );
    let ep = p(position_error);
    let ed = p(derivative_error);
    let jac = tau - ep * tau - ed - ep * ed;
    if position_error >= 1. || jac.lo <= 0. {
        return Err(Error::Unresolved);
    }
    Ok(jac.lo)
}
impl PeriodicAngular {
    pub fn prepare(
        left: bool,
        options: CurveOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, Error> {
        let mut w = Work::new(options.work, interrupt)?;
        w.charge(1)?;
        options.validate()?;
        let mut memory = Memory::new(size_of::<Self>(), options.bytes)?;
        let curve = super::super::insert::fit_closed_unit(left, options, &mut memory, &mut w)?;
        w.charge(64)?;
        // det(B,B') differs from signed true tau by at most
        // ep*|C'| + ed*|C| + ep*ed, with |C|=1. Actual decoded-token
        // coefficient uncertainty is already included by the fitter.
        let jacobian_lower = radial_proof(curve.position_error, curve.derivative_error)?;
        // Endpoints have EXACT common literal xy=(1,0). Positive signed
        // angular Jacobian and the nonzero straight homotopy to the true
        // source circle (ep<1) prove degree +/-1, excluding double wraps
        // and nonadjacent contact. Every eventual repeat REUSES these same
        // coefficient bytes; its axial station is a separate source owner.
        w.charge(1)?;
        Ok(Self {
            left,
            curve,
            jacobian_lower,
            work: w.used,
            bytes: memory.used,
        })
    }
    pub fn left_hand(&self) -> bool {
        self.left
    }
    pub fn poles(&self) -> &[[f64; 3]] {
        &self.curve.poles
    }
    pub fn knots(&self) -> &[f64] {
        &self.curve.knots
    }
    pub fn position_error(&self) -> f64 {
        self.curve.position_error
    }
    /// Derivative in physical unit turns, since the original domain is[0,1].
    pub fn derivative_error(&self) -> f64 {
        self.curve.derivative_error
    }
    pub fn signed_angular_jacobian_lower(&self) -> f64 {
        self.jacobian_lower
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_bytes(&self) -> usize {
        self.bytes
    }
}
#[cfg(test)]
mod tests;
