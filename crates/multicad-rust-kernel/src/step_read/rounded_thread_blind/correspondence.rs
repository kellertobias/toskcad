//! Total positional allowance for the admitted rounded exchange class.
//! This combines the already checked whole-curve fitting receipts, actual
//! affine amplification, source evaluation/REAL conditioning, exact stock
//! scalar transfer and complete analytic closure checks. No D1/mass follows.
use super::*;
use crate::{rounded_thread::BlindBore, tapered_thread::interval::I};
pub struct RoundedBlindCorrespondence<'input, 'source> {
    input: &'input str,
    source: &'source BlindBore,
    frame: [Point; 3],
    total: f64,
    diagnostic: f64,
}
impl<'input, 'source> RoundedBlindCorrespondence<'input, 'source> {
    pub fn original_input(&self) -> &'input str {
        self.input
    }
    pub fn source(&self) -> &'source BlindBore {
        self.source
    }
    pub fn frame(&self) -> [Point; 3] {
        self.frame
    }
    pub fn maximum_position_error_mm(&self) -> f64 {
        self.total
    }
    /// Numerical decoded spline flux, not sealed raw physical mass.
    pub fn decoded_flux_diagnostic_mm3(&self) -> f64 {
        self.diagnostic
    }
}
pub fn verify_rounded_blind_correspondence<'input, 'source>(
    input: &'input str,
    source: &'source BlindBore,
    origin: Point,
    u: Point,
    v: Point,
    tolerance: f64,
) -> Result<RoundedBlindCorrespondence<'input, 'source>> {
    let diagnostic = verify_rounded_thread_blind_bore(input, source, origin, u, v, tolerance)?;
    let frame = [origin, u, v];
    let total = total_bound(source, frame, tolerance, false)?;
    Ok(RoundedBlindCorrespondence {
        input,
        source,
        frame,
        total,
        diagnostic,
    })
}
pub(super) fn total_bound(
    bore: &BlindBore,
    [origin, u, v]: [Point; 3],
    tolerance: f64,
    repaired: bool,
) -> Result<f64> {
    if !tolerance.is_finite() || tolerance <= 0. {
        return Err(ReadError::Invalid);
    }
    let p = I::point;
    let t = p(tolerance.min(1e-7));
    // prepare() re-earns the complete final controls (including shared ends):
    // radial L1 error <=t/32, unit rotation L1 error <=t/(32*major).
    let splines = crate::step::rounded_thread::prepare(bore.cutter(), origin, t.hi)
        .map_err(helical::step_error)?;
    let mut axial = 0_f64;
    for a in &splines.axial.poles {
        let value = p(a[0].abs()) + p(a[1].abs());
        if !value.valid() {
            return Err(ReadError::Limit);
        }
        axial = axial.max(value.hi);
    }
    let mut amplification = 0_f64;
    for column in [u, v, cross(u, v)] {
        let value = column.into_iter().fold(p(0.), |s, x| s + p(x.abs()));
        if !value.valid() {
            return Err(ReadError::Limit);
        }
        amplification = amplification.max(value.hi);
    }
    let thread = bore.cutter();
    let spec = thread.spec();
    let mut extent = 0_f64;
    let mut center = 0_f64;
    for i in 0..5 {
        if let Some(c) = thread.profile_circle(i) {
            extent = extent.max((p(c.center_radius_mm.abs()) + p(c.radius_mm)).hi);
            center = center.max(c.center_error_mm);
        }
    }
    let o = origin.into_iter().map(f64::abs).fold(0., f64::max);
    // Same finite arithmetic/conditioning guard as prepare(), reconstructed
    // outward. It covers the complete admitted profile/rotation/height and
    // coefficient-expression evaluation, not merely one control coordinate.
    let scale = p(o)
        + p(extent)
        + p(spec.major_radius_mm) * (p(spec.length_mm) / p(thread.helical_rate().abs()))
        + p(spec.major_radius_mm)
        + p(spec.length_mm)
        + p(spec.pitch_mm)
        + p(1.);
    let numeric = p(256. * f64::EPSILON) * scale + p(center);
    let local_fit = t / p(32.) * p(axial)
        + p(bore.cutter_radial_bounds_mm()[1]) * t / (p(32.) * p(spec.major_radius_mm));
    let delta =
        |bounds: [f64; 2], literal: f64| p((I::new(bounds[0], bounds[1]) - p(literal)).magnitude());
    let stock = p(amplification)
        * (p(2.)
            * delta(
                bore.stock_radius_bounds_mm(),
                bore.stock_radius_mm().as_f64(),
            )
            + delta(
                bore.stock_height_bounds_mm(),
                bore.stock_height_mm().as_f64(),
            ));
    // Actual REAL leaves and normalized quasiorthogonal analytic axes were
    // admitted by raw_axis_class/precision. Adding the prepare scale covers
    // thread controls as well as the finite stock/world coordinates.
    let raw = p(256. * f64::EPSILON)
        * (scale + p(bore.stock_radius_bounds_mm()[1]) + p(bore.stock_height_bounds_mm()[1]));
    // Full cap normals and outer-circle/cylinder closure each independently
    // earn t/32 over their actual finite extent. Summing four is conservative
    // even though separate faces would only require their maximum.
    let closure = p(4.) * t / p(32.);
    // Repaired comparisons allow <=t/16 in WORLD control coordinates. Twice
    // that allowance includes subtraction/norm roundoff; it is not remapped.
    let controls = if repaired { t / p(8.) } else { p(0.) };
    let total =
        p(amplification) * (local_fit + numeric) + numeric + stock + raw + closure + controls;
    if !total.valid() || total.hi > t.lo {
        return Err(ReadError::Limit);
    }
    Ok(total.hi)
}
