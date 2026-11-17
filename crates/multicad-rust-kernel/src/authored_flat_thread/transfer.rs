//! Whole-band correspondence; moving crease normals are never compared pointwise.
use super::*;
fn abs(v: I) -> I {
    I::point(v.magnitude())
}
fn positive(v: I) -> Result<I, Error> {
    if v.valid() && v.lo > 0. {
        Ok(v)
    } else {
        Err(Error::Unresolved)
    }
}
pub(super) fn certify(
    spec: Spec,
    phases: [f64; 6],
    radii: [f64; 6],
    legacy: &Thread,
    axial_end: I,
    options: SourceOptions,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<TransferBounds, Error> {
    let true_gain = angular_gain(spec);
    let old_rate = I::point(legacy.helical_rate());
    let true_rate = I::point(1.) / true_gain;
    let gain_error = abs(true_gain - I::point(1.) / old_rate);
    let span = axial_end + I::point(spec.pitch_mm) * I::point(f64::from(spec.starts));
    let angular_error = gain_error * span;
    let mut bands = [BandTransfer {
        position_error_mm: 0.,
        slope_error: 0.,
        normal_angle_error_rad: 0.,
    }; 5];
    for band in 0..5 {
        checkpoint(interrupt)?;
        let [qa, qb] = [phases[band], phases[band + 1]];
        let [ra, rb] = [radii[band], radii[band + 1]];
        let width = positive(I::point(qb) - I::point(qa))?;
        let shared_slope = (I::point(rb) - I::point(ra)) / width;
        let old_slope = I::point(legacy.slope_at_phase(qa, band));
        let old_radius = |q: f64| {
            if band == 0 || band == 4 {
                I::point(spec.major_radius_mm)
            } else if band == 2 {
                I::point(legacy.root_radius_mm())
            } else {
                I::point(spec.major_radius_mm)
                    + old_slope
                        * (I::point(q) - I::point(if band == 1 { phases[1] } else { phases[4] }))
            }
        };
        let a = old_radius(qa);
        let b = old_radius(qb);
        let radius_error = (a - I::point(ra))
            .magnitude()
            .max((b - I::point(rb)).magnitude());
        let slope_error = abs(shared_slope - old_slope);
        let minimum = positive(I::new(
            a.lo.min(b.lo).min(ra.min(rb)),
            a.hi.max(b.hi).max(ra.max(rb)),
        ))?;
        let rmin = I::point(minimum.lo);
        let rmax = I::point(minimum.hi);
        // Coupled algebra avoids subtracting two independent wide beta hulls.
        // beta=true_rate*shared_slope/shared_radius in the radial normal basis.
        let beta_error = abs(true_rate - old_rate) * abs(shared_slope) / rmin
            + slope_error * abs(old_rate) / rmin
            + abs(old_slope * old_rate) * I::point(radius_error) / rmin.square();
        let beta = abs(old_rate * old_slope) / rmin;
        let horizontal = (I::point(1.) + beta.square()).sqrt();
        let normal_difference = horizontal * angular_error + beta_error + slope_error;
        // Both unscaled normals have radial component one, hence norm>=1.
        // Unit chord <=2*D; for D<=.5, angle<=2asin(D)<=4D.
        if !normal_difference.valid() || normal_difference.hi > 0.5 {
            return Err(Error::Unresolved);
        }
        let angle = I::point(4.) * normal_difference;
        let position = I::point(radius_error) + rmax * angular_error;
        if !angle.valid()
            || !position.valid()
            || angle.hi > options.normal_angle_limit_rad
            || position.hi > options.position_limit_mm
        {
            return Err(Error::Unresolved);
        }
        bands[band] = BandTransfer {
            position_error_mm: position.hi,
            slope_error: slope_error.hi,
            normal_angle_error_rad: angle.hi,
        };
    }
    Ok(TransferBounds {
        maximum_position_error_mm: bands
            .into_iter()
            .map(|b| b.position_error_mm)
            .fold(0., f64::max),
        maximum_normal_angle_error_rad: bands
            .into_iter()
            .map(|b| b.normal_angle_error_rad)
            .fold(0., f64::max),
        bands,
        angular_error_rad: angular_error.hi,
    })
}
