//! Conservative whole-prism relations against the actual retained annulus.
use super::*;
use crate::geometry::{
    helical_body::contained_boolean::{chart_bounds, probe_bounds, whole_operand_result},
    internal_thread_treatment_body::interval::Interval as B,
};
use multicad_rust_kernel::helical_thread::BoxClassification;

fn relation(source: &Annulus, bounds: [[f64; 2]; 3]) -> BoxClassification {
    use BoxClassification::*;
    let [x, y, z] = bounds.map(|[lo, hi]| B::new(lo, hi));
    if [x, y, z].iter().any(|v| !v.valid()) {
        return UnresolvedBoundary;
    }
    let radius = x.square() + y.square();
    let inner = B::point(source.inner_radius_mm()).square();
    let outer = B::point(source.outer_radius_mm()).square();
    if z.hi < 0. || z.lo > source.height_mm() {
        return Outside;
    }
    if !radius.valid() {
        UnresolvedBoundary
    } else if radius.hi < inner.lo || radius.lo > outer.hi {
        Outside
    } else if z.lo > 0. && z.hi < source.height_mm() && radius.lo > inner.hi && radius.hi < outer.lo
    {
        Inside
    } else {
        UnresolvedBoundary
    }
}
/// An axial ring of material around its frame origin: radial interval and
/// height along the frame normal. A solid cylinder has inner radius zero.
fn ring(solid: &Solid) -> Option<(PlaneFrame, [f64; 2], f64)> {
    if let Some((body, frame)) = solid.annulus_body() {
        let geometry = &body.geometry;
        return Some((
            frame,
            [geometry.inner_radius_mm(), geometry.outer_radius_mm()],
            geometry.height_mm(),
        ));
    }
    let (cylinder, frame) = solid.axial_cylinder_frame()?;
    Some((
        frame,
        [0., cylinder.circle().radius().as_f64()],
        cylinder.height().as_f64(),
    ))
}

/// Concentric rings whose radial or axial intervals at most touch share no
/// interior. Only rounding of the placement is absorbed, never a tolerance.
fn concentric_rings_disjoint(a: &Solid, b: &Solid) -> Option<bool> {
    let (af, ar, ah) = ring(a)?;
    let (bf, br, bh) = ring(b)?;
    let axis = af.normal();
    let other = bf.normal();
    let magnitude = af
        .origin
        .into_iter()
        .chain(bf.origin)
        .chain([ah, bh, ar[1], br[1]])
        .map(f64::abs)
        .fold(1., f64::max);
    let margin = 64. * f64::EPSILON * magnitude;
    let parallel = dot(axis, other);
    if (parallel.abs() - 1.).abs() > 64. * f64::EPSILON {
        return None;
    }
    let delta: [f64; 3] = std::array::from_fn(|i| bf.origin[i] - af.origin[i]);
    let along = dot(delta, axis);
    let offset: [f64; 3] = std::array::from_fn(|i| delta[i] - axis[i] * along);
    let radial = dot(offset, offset).sqrt();
    if radial > margin {
        return None;
    }
    let end = along + parallel.signum() * bh;
    let [lo, hi] = [along.min(end), along.max(end)];
    Some(hi <= margin || lo >= ah - margin || ar[1] <= br[0] + margin || br[1] <= ar[0] + margin)
}

impl Solid {
    /// Concentric cylinders and annuli that only touch: an empty
    /// intersection, or an unchanged minuend.
    pub(in crate::geometry) fn concentric_ring_boolean(
        &self,
        other: &Self,
        operation: BooleanOperation,
        policy: GeometryPolicy,
    ) -> Option<Result<Option<Self>, GeometryError>> {
        if policy.boolean_fuzzy_mm != 0.
            || self.policy != policy
            || other.policy != policy
            || operation == BooleanOperation::Union
            || !concentric_rings_disjoint(self, other)?
        {
            return None;
        }
        Some(match operation {
            BooleanOperation::Subtract => self.share().map(Some),
            _ => Ok(None),
        })
    }

    pub(in crate::geometry) fn contained_annulus_boolean(
        &self,
        other: &Self,
        operation: BooleanOperation,
        policy: GeometryPolicy,
    ) -> Option<Result<Option<Self>, GeometryError>> {
        if policy.boolean_fuzzy_mm != 0. || self.policy != policy || other.policy != policy {
            return None;
        }
        let (probe, source_is_left, body, frame) = if let Some((body, frame)) = self.annulus_body()
        {
            (other, true, body, frame)
        } else {
            let (body, frame) = other.annulus_body()?;
            (self, false, body, frame)
        };
        let bounds = chart_bounds(
            probe_bounds(probe)?,
            frame.origin,
            [frame.u, frame.v, frame.normal()],
        );
        whole_operand_result(
            self,
            other,
            source_is_left,
            relation(&body.geometry, bounds),
            operation,
        )
    }
}

#[cfg(test)]
mod tests;
