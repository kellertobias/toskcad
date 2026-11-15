//! Verify a cubic rim against its master parameter line on a linear extrusion.
use super::*;

impl Shell<'_, '_> {
    /// `S(u, v) = C(u) + v V`. The master line `(u, v) = (t, v0)` transmits the
    /// swept curve translated by `v0 V` with unchanged poles and knots. Rulings
    /// and reparametrized rims are outside this proof.
    pub(super) fn verify_extrusion_rim_master(
        &self,
        curve_3d: u32,
        surface: u32,
        line: u32,
    ) -> Result<()> {
        let (swept, direction, speed) = self.spline_extrusion_vector(surface)?;
        let ([u, v, _], tangent, scale) = self.line(line, true)?;
        if u != 0. || tangent != [1., 0., 0.] || scale != 1. {
            return Err(ReadError::Unsupported);
        }
        let travel = direction.map(|axis| axis * speed * v);
        if !self.spline_curve(swept)?.translated(
            &self.spline_curve(curve_3d)?,
            travel,
            self.tolerance,
        ) {
            return Err(ReadError::Invalid);
        }
        Ok(())
    }
}
