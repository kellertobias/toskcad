//! Preserve a literal cone angle and bound its finite-radius scalar conversion.
use super::*;
#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct LiteralAngle {
    pub(super) angle: f64,
    pub(super) error: TransportError,
}
impl FiniteBinding {
    /// Attach actual support data. Cache volume and admit a solid only after
    /// attachment, using the resulting spatial/conversion error bounds.
    pub fn with_literal_cone_semi_angle(mut self, angle: f64) -> Result<Self, Error> {
        if !angle.is_finite() || angle <= 0. || angle >= PI / 2. {
            return Err(Error::Invalid);
        }
        let numeric = I::new(angle.next_down(), angle.next_up());
        let sine = numeric.sin();
        if !sine.valid() || sine.lo <= 0. {
            return Err(Error::Unresolved);
        }
        // The trig enclosure subtracts8eps. Its smallest positive sine
        // residual exceeds2^-103, so division of the bounded cosine is finite.
        let numeric_slope = numeric.cos() / sine;
        let k = self.source.cone().depth_per_radius;
        if k < numeric_slope.lo || k > numeric_slope.hi {
            return Err(Error::Invalid);
        }
        let literal = I::point(angle);
        let slope = literal.cos() / literal.sin();
        // Include the retained trim and source tube, not just the nominal rim.
        let radius = I::point(self.source.maximum_radius_bound_mm())
            + I::point(self.nominal_binding_error_mm)
            + I::point(16.) * I::point(self.source.source_error_mm());
        let arithmetic = (radius * I::point((slope - I::point(k)).magnitude())).hi;
        // d(cot a)/da=-csc²a. Bound decimal-to-f64 conversion over its
        // entire adjacent-value interval, independently of geometry policy.
        let input = (numeric - I::point(angle)).magnitude();
        let numeric_conversion = (radius * I::point(input) / sine.square()).hi;
        // The admitted carrier/line limits bound the inflated radius; coupled
        // k<=1e7 keeps this finite. Neither bound uses a caller tolerance.
        self.literal_cone_angle = Some(LiteralAngle {
            angle,
            error: TransportError {
                numeric_conversion_mm: numeric_conversion,
                arithmetic_mm: arithmetic,
            },
        });
        Ok(self)
    }
    pub fn literal_cone_semi_angle_rad(&self) -> Option<f64> {
        self.literal_cone_angle.map(|a| a.angle)
    }
    /// Millimetre floor-height conversion bounds over the complete finite
    /// radius. Coordinate transport remains a separate source quantity.
    pub fn cone_angle_conversion_error(&self) -> TransportError {
        self.literal_cone_angle
            .map_or(TransportError::ZERO, |a| a.error)
    }
}
