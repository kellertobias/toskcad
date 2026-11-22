//! Conservative whole-prism relations against a retained radial boundary
//! (tapered and custom conical threads, authored pipe threads).
//!
//! No cut surface is computed. An authored pipe thread classifies the probe
//! with its analytic law, after dilating the probe by the recorded
//! law-to-retained error and both retained boundaries' sewing collars, so
//! the relation holds for the retained material itself. Otherwise the
//! retained source's own certified box predicate decides. A proved relation
//! shares one unchanged operand; probes reaching the boundary are left to
//! other providers.
use super::*;
use crate::geometry::helical_body::contained_boolean::{
    chart_bounds, probe_bounds, whole_operand_result,
};
use crate::geometry::internal_thread_treatment_body::interval::Interval as B;
use multicad_rust_kernel::axial_pipe_thread::authored::AuthoredBoxRelation;
use multicad_rust_kernel::helical_thread::BoxClassification;

fn local(world: [B; 3], frame: PlaneFrame) -> [B; 3] {
    chart_bounds(world, frame.origin, [frame.u, frame.v, frame.normal()])
        .map(|[lo, hi]| B::new(lo, hi))
}

impl RadialBoundaryBody {
    fn authored_relation(
        &self,
        frame: PlaneFrame,
        placed: &RadialBoundary,
        world: [B; 3],
    ) -> Option<BoxClassification> {
        let source = &self.authored_pipe.as_ref()?.source;
        let chart = local(local(world, frame), source.construction_frame);
        let magnitude = world
            .iter()
            .chain(&chart)
            .flat_map(|b| [b.lo.abs(), b.hi.abs()])
            .chain(frame.origin.map(f64::abs))
            .chain(source.construction_frame.origin.map(f64::abs))
            .fold(1., f64::max);
        // The law-to-net transfer, both certified material collars, and the
        // rounding of two re-posed pole sets.
        let margin = (source.local_source_error_mm
            + self.source.boundary.boundary_error_mm()
            + placed.boundary_error_mm()
            + 128. * f64::EPSILON * magnitude)
            * (1. + 16. * f64::EPSILON);
        if !margin.is_finite() || chart.iter().any(|b| !b.valid()) {
            return None;
        }
        let dilated = chart.map(|b| [(b.lo - margin).next_down(), (b.hi + margin).next_up()]);
        match source.transport.classify_box(dilated).ok()?? {
            AuthoredBoxRelation::Inside => Some(BoxClassification::Inside),
            AuthoredBoxRelation::Outside => Some(BoxClassification::Outside),
        }
    }
}

impl Solid {
    pub(in crate::geometry) fn contained_radial_boundary_boolean(
        &self,
        other: &Self,
        operation: BooleanOperation,
        policy: GeometryPolicy,
    ) -> Option<Result<Option<Self>, GeometryError>> {
        if policy.boolean_fuzzy_mm != 0. || self.policy != policy || other.policy != policy {
            return None;
        }
        let (probe, source_is_left, (body, frame)) = match self.radial_boundary_body() {
            Some(found) => (other, true, found),
            None => (self, false, other.radial_boundary_body()?),
        };
        let world = probe_bounds(probe)?;
        let placed = match body.placement(frame, None) {
            Ok(placed) => placed,
            Err(error) => return Some(Err(error)),
        };
        let relation = match body.authored_relation(frame, &placed.boundary, world) {
            Some(relation) => relation,
            None => match placed
                .boundary
                .classify_box(world.map(|b| [b.lo, b.hi]), Default::default(), None)
                .ok()?
                .location()
            {
                Location::Inside => BoxClassification::Inside,
                Location::Outside => BoxClassification::Outside,
                Location::Unresolved => return None,
            },
        };
        whole_operand_result(self, other, source_is_left, relation, operation)
    }
}
