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
impl Solid {
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
