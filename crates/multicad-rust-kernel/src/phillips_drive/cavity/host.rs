//! Actual narrowing envelope containment; no thread source is modified.
use super::*;
use crate::{
    head_profile::{Curve as HeadCurve, Profile, Role},
    helical_thread::Thread,
};

/// Source-envelope containment only. This does not sew the mouth, certify a
/// manufacturing mesh, or validate a foreign file's finite incidence.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HostEnvelopeProof {
    cells: usize,
    minimum_radial_clearance_mm: f64,
    axial_floor_clearance_mm: f64,
}
impl HostEnvelopeProof {
    pub fn cells(&self) -> usize {
        self.cells
    }
    pub fn minimum_radial_clearance_mm(&self) -> f64 {
        self.minimum_radial_clearance_mm
    }
    pub fn axial_floor_clearance_mm(&self) -> f64 {
        self.axial_floor_clearance_mm
    }
}
// Every admitted head curve is a graph on its retained axial interval, with
// positive-circle Crown or negative-circle Neck support. Expand the queried
// axial interval AND radial result by the profile's explicit source reserve.
// This avoids treating rounded trim vertices as exact support intersections.
fn head_radius_lower(profile: &Profile, z: I) -> f64 {
    let reserve = I::point(profile.rounding_error_mm());
    let z = z + I::new(-reserve.hi, reserve.hi);
    let mut lower = f64::INFINITY;
    for s in profile.segments() {
        let start = (I::point(s.start[1]) - reserve).lo;
        let end = (I::point(s.end[1]) + reserve).hi;
        if z.hi < start || z.lo > end {
            continue;
        }
        // The literal endpoint may itself differ from the circle support by
        // the declared source reserve. Preserve that endpoint domain too;
        // clipping to its rounded Z would lose near-tangent support points.
        let range = I::new(z.lo.max(start), z.hi.min(end));
        let radius = match s.curve {
            HeadCurve::Line if s.start[1] == s.end[1] => I::point(s.start[0].min(s.end[0])),
            HeadCurve::Line => {
                let t =
                    (range - I::point(s.start[1])) / (I::point(s.end[1]) - I::point(s.start[1]));
                I::point(s.start[0]) + t * (I::point(s.end[0]) - I::point(s.start[0]))
            }
            HeadCurve::Arc {
                center, radius_mm, ..
            } => {
                let z = range - I::point(center[1]);
                let squared = I::point(radius_mm).square() - z.square();
                // The exact supported circle exists only where this radicand
                // is nonnegative; intersect that analytic domain before sqrt.
                let root = I::new(squared.lo.max(0.), squared.hi.max(0.)).sqrt();
                if s.role == Role::Crown {
                    I::point(center[0]) + root
                } else {
                    I::point(center[0]) - root
                }
            }
        };
        lower = lower.min((radius - reserve).lo);
    }
    if lower.is_finite() { lower } else { 0. }
}
impl Cavity {
    /// Prove this actual depth-dependent cut stays within the retained head
    /// and the original thread's common cylindrical material core. `plain_mm`
    /// is the original thread datum in the profile coordinates; no phase,
    /// root, thread length or visible-prefix trim is changed by this proof.
    /// `wall_mm` specifies radial and axial margins, not a normal-offset wall.
    pub fn certify_host_envelope(
        &self,
        profile: &Profile,
        thread: &Thread,
        plain_mm: f64,
        wall_mm: f64,
        work_limit: usize,
        cancelled: &dyn Fn() -> bool,
    ) -> Result<HostEnvelopeProof, Error> {
        if ![plain_mm, wall_mm].into_iter().all(f64::is_finite)
            || plain_mm < 0.
            || wall_mm <= 0.
            || profile.end()[1] < plain_mm
            || profile.end()[0] < thread.spec().major_radius_mm
        {
            return Err(Error::Invalid);
        }
        if work_limit == 0 || work_limit > 1_000_000 {
            return Err(Error::Limit);
        }
        if cancelled() {
            return Err(Error::Cancelled);
        }
        let reserve = I::point(2.) * I::point(profile.rounding_error_mm());
        let top = I::point(profile.top()[1]);
        let tip = self.source.envelope.tip;
        let end = I::point(plain_mm) + I::point(thread.spec().length_mm);
        let floor = (end - (top + tip) - reserve).lo;
        if floor <= wall_mm {
            return Err(Error::Unresolved);
        }
        let mouth = I::point(self.radial_suffix_bound_mm(0.)?);
        if (I::point(profile.top()[0]) - mouth - reserve).lo <= wall_mm {
            return Err(Error::Unresolved);
        }
        let mut core =
            (I::point(thread.spec().major_radius_mm) - I::point(thread.spec().depth_mm)).lo;
        if let Some(root) = thread.root_circle() {
            core = core.min((I::point(root.center_radius) - I::point(root.radius)).lo);
        }
        let mut pending = vec![[0., tip.hi]];
        let mut cells = 0;
        let mut clearance = f64::INFINITY;
        while let Some([a, b]) = pending.pop() {
            if cancelled() {
                return Err(Error::Cancelled);
            }
            if cells == work_limit {
                return Err(Error::Limit);
            }
            cells += 1;
            let z = top + I::new(a, b);
            let head = if z.lo <= plain_mm {
                head_radius_lower(profile, I::new(z.lo, z.hi.min(plain_mm)))
            } else {
                f64::INFINITY
            };
            let material = if z.hi >= plain_mm {
                head.min(core)
            } else {
                head
            };
            let cut = self.radial_suffix_bound_mm(a)?;
            let gap = (I::point(material) - I::point(cut) - reserve).lo;
            if gap > wall_mm {
                clearance = clearance.min(gap);
            } else {
                let mid = a + (b - a) / 2.;
                if mid <= a || mid >= b {
                    return Err(Error::Unresolved);
                }
                pending.push([mid, b]);
                pending.push([a, mid]);
            }
        }
        Ok(HostEnvelopeProof {
            cells,
            minimum_radial_clearance_mm: clearance,
            axial_floor_clearance_mm: floor,
        })
    }
}

#[cfg(test)]
mod tests;
