//! Bounded approximation to an explicitly supplied original literal Thread.
//! Regular band normals are compared by source correspondence. The old tiny
//! radial closure faces are not normal-equivalent to this new continuous law.
//! This compares real analytic supports with retained literal coefficients,
//! rate and indexed stations. It does NOT certify floating Thread.radius() or
//! contains() evaluation at arbitrary arguments; literal emission and query
//! rounding need independent proofs. Direction: each old regular-band source
//! point maps to the new law at the same normalized band coordinate and z.
//! Same-physical-angle radial/material bounds instead use the NEW continuous
//! profile's global Lipschitz constant sqrt(3), including old one-sided limits.
use super::*;
use crate::helical_thread::Thread;
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TransferOptions {
    pub position_limit_mm: f64,
    pub normal_angle_limit_rad: f64,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BandTransfer {
    pub phase_error_mm: f64,
    pub radius_error_mm: f64,
    pub slope_error: f64,
    pub position_error_mm: f64,
    pub regular_normal_angle_error_rad: f64,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TransferBounds {
    pub bands: [BandTransfer; 5],
    pub maximum_position_error_mm: f64,
    pub maximum_regular_normal_angle_error_rad: f64,
    /// Same physical angle/axial point; covers old one-sided rail limits.
    pub material_radial_error_mm: f64,
    pub material_volume_difference_bound_mm3: f64,
    pub root_radius_error_mm: f64,
}
fn abs(v: I) -> I {
    I::point(v.magnitude())
}
impl AuthoredMetricThread {
    pub fn certify_legacy(
        &self,
        legacy: &Thread,
        options: TransferOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<TransferBounds, Error> {
        checkpoint(interrupt)?;
        if [options.position_limit_mm, options.normal_angle_limit_rad]
            .into_iter()
            .any(|x| !x.is_finite() || x <= 0.)
        {
            return Err(Error::Invalid);
        }
        let a = self.spec;
        let b = legacy.spec();
        if a.major_radius_mm != b.major_radius_mm
            || a.pitch_mm != b.pitch_mm
            || a.length_mm != b.length_mm
            || a.crest_width_mm != b.crest_width_mm
            || b.root_radius_mm.is_none()
            || a.starts != b.starts
            || a.left_hand != b.left_hand
            || b.lower_flank_deg != 30.
            || b.upper_flank_deg != 30.
        {
            return Err(Error::Invalid);
        }
        let old = legacy.profile_breaks();
        let knots = self.knots();
        let root = legacy.root_circle().expect("matched rounded source");
        let rate = self.rate();
        let old_rate = I::point(legacy.helical_rate());
        let span =
            I::point(a.length_mm) + I::point(a.pitch_mm) * I::point(f64::from(a.starts) + 2.);
        let mut indexed = 0_f64;
        let margin = 2 * i64::from(a.starts) + 4;
        let last = (a.length_mm / a.pitch_mm).ceil() as i64 + margin;
        for period in -margin..=last {
            checkpoint(interrupt)?;
            for q in old {
                let exact = I::point(period as f64) * I::point(a.pitch_mm) + I::point(q);
                for literal in [
                    (period as f64).mul_add(a.pitch_mm, q),
                    period as f64 * a.pitch_mm + q,
                ] {
                    indexed = indexed.max((I::point(literal) - exact).magnitude());
                }
            }
        }
        let mut old_min = I::point(root.center_radius) - I::point(root.radius);
        for (band, end, base) in [(1usize, 2usize, 1usize), (3, 3, 4)] {
            let v = I::point(b.major_radius_mm)
                + I::point(legacy.slope_at_phase(old[base], band))
                    * (I::point(old[end]) - I::point(old[base]));
            old_min = I::new(old_min.lo.min(v.lo), old_min.hi.max(v.hi));
        }
        let new_min = interval(self.minimum_radius_mm());
        let radial_support = I::new(new_min.lo.min(old_min.lo), new_min.hi.max(old_min.hi));
        let rmin = positive(radial_support)?;
        let rmax = I::point(
            a.major_radius_mm
                .max(legacy.source_coordinate_bounds_mm()[0][1]),
        );
        let mut bands = [BandTransfer {
            phase_error_mm: 0.,
            radius_error_mm: 0.,
            slope_error: 0.,
            position_error_mm: 0.,
            regular_normal_angle_error_rad: 0.,
        }; 5];
        let angular_convention = abs(I::point(std::f64::consts::TAU) - I::point(2.) * pi());
        let mut radial_error = 0_f64;
        for band in 0..5 {
            checkpoint(interrupt)?;
            let phase = I::point(
                (knots[band][0] - I::point(old[band]))
                    .magnitude()
                    .max((knots[band + 1][0] - I::point(old[band + 1])).magnitude()),
            ) + I::point(indexed);
            let (radius_error, slope_error, old_slope) = if band == 2 {
                let new_x = [-sqrt3() / I::point(2.), sqrt3() / I::point(2.)];
                let old_x = [old[2], old[3]]
                    .map(|q| (I::point(q) - I::point(root.phase)) / I::point(root.radius));
                let dx = I::point(
                    (new_x[0] - old_x[0])
                        .magnitude()
                        .max((new_x[1] - old_x[1]).magnitude()),
                );
                let maximum = new_x
                    .into_iter()
                    .chain(old_x)
                    .map(I::magnitude)
                    .fold(0., f64::max);
                let root_height = (I::point(1.) - I::point(maximum).square()).sqrt();
                let h = positive(root_height)?;
                let slope = I::point(maximum) / h;
                let derivative = I::point(1.) / (h * h * h);
                (
                    abs(self.center_radius() - I::point(root.center_radius))
                        + abs(self.root_radius() - I::point(root.radius))
                        + I::point(root.radius) * slope * dx,
                    derivative * dx,
                    slope,
                )
            } else if band == 1 || band == 3 {
                let base = if band == 1 { old[1] } else { old[4] };
                let m = I::point(legacy.slope_at_phase(base, band));
                let radius = [band, band + 1]
                    .into_iter()
                    .map(|i| {
                        abs(knots[i][1]
                            - (I::point(b.major_radius_mm)
                                + m * (I::point(old[i]) - I::point(base))))
                    })
                    .fold(I::point(0.), |a, b| I::point(a.hi.max(b.hi)));
                let new = if band == 1 { -sqrt3() } else { sqrt3() };
                (radius, abs(new - m), abs(m))
            } else {
                (I::point(0.), I::point(0.), I::point(0.))
            };
            let angular = abs(I::point(1.) / rate - I::point(1.) / old_rate) * span
                + phase / abs(rate)
                + angular_convention;
            let m = if band == 0 || band == 4 {
                I::point(0.)
            } else {
                sqrt3()
            };
            let beta_error = abs(rate - old_rate) * m / I::point(rmin.lo)
                + slope_error * abs(old_rate) / I::point(rmin.lo)
                + old_slope * abs(old_rate) * radius_error / I::point(rmin.lo).square();
            let beta = old_slope * abs(old_rate) / I::point(rmin.lo);
            let normal = (I::point(1.) + beta.square()).sqrt() * angular + beta_error + slope_error;
            let angle = I::point(4.) * normal;
            let position = radius_error + rmax * angular;
            if !normal.valid()
                || normal.hi > 0.5
                || !position.valid()
                || position.hi > options.position_limit_mm
                || angle.hi > options.normal_angle_limit_rad
            {
                return Err(Error::Unresolved);
            }
            bands[band] = BandTransfer {
                phase_error_mm: phase.hi,
                radius_error_mm: radius_error.hi,
                slope_error: slope_error.hi,
                position_error_mm: position.hi,
                regular_normal_angle_error_rad: angle.hi,
            };
            radial_error = radial_error.max(
                (radius_error
                    + sqrt3()
                        * (phase
                            + abs(rate - old_rate) * (I::point(2.) * pi())
                            + abs(old_rate) * angular_convention))
                    .hi,
            );
        }
        let e = I::point(radial_error);
        let volume = pi() * I::point(a.length_mm) * (I::point(2.) * rmax * e + e.square());
        let result = TransferBounds {
            maximum_position_error_mm: bands
                .into_iter()
                .map(|b| b.position_error_mm)
                .fold(0., f64::max),
            maximum_regular_normal_angle_error_rad: bands
                .into_iter()
                .map(|b| b.regular_normal_angle_error_rad)
                .fold(0., f64::max),
            bands,
            root_radius_error_mm: (self.root_radius() - I::point(root.radius)).magnitude(),
            material_radial_error_mm: radial_error,
            material_volume_difference_bound_mm3: volume
                .checked()
                .map_err(|_| Error::Unresolved)?
                .hi,
        };
        checkpoint(interrupt)?;
        Ok(result)
    }
}
