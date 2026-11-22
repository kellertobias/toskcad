//! Conservative whole-prism relations against an actual rounded thread.
use super::*;
use crate::geometry::helical_body::contained_boolean::{
    chart_bounds, probe_bounds, whole_operand_result,
};

enum Source<'a> {
    Thread(&'a multicad_rust_kernel::rounded_thread::Thread),
    Bore(&'a multicad_rust_kernel::rounded_thread::ThroughBore),
}

fn source(solid: &Solid) -> Option<(Source<'_>, PlaneFrame)> {
    solid
        .rounded_thread_body()
        .map(|(body, frame)| (Source::Thread(&body.geometry), frame))
        .or_else(|| {
            solid
                .rounded_thread_bore_body()
                .map(|(body, frame)| (Source::Bore(&body.geometry), frame))
        })
}

impl Solid {
    /// A finite prism wholly inside the thread's core or teeth, or wholly
    /// outside them, reuses one unchanged operand. Probes reaching the
    /// thread surface or a cap remain with other providers.
    pub(in crate::geometry) fn contained_rounded_thread_boolean(
        &self,
        other: &Self,
        operation: BooleanOperation,
        policy: GeometryPolicy,
    ) -> Option<Result<Option<Self>, GeometryError>> {
        if policy.boolean_fuzzy_mm != 0. || self.policy != policy || other.policy != policy {
            return None;
        }
        // The thread or its through receiver, whichever operand it is.
        let (probe, source_is_left, (geometry, frame)) = match source(self) {
            Some(found) => (other, true, found),
            None => (self, false, source(other)?),
        };
        let bounds = chart_bounds(
            probe_bounds(probe)?,
            frame.origin,
            [frame.u, frame.v, frame.normal()],
        );
        let relation = match geometry {
            Source::Thread(thread) => thread.classify_box(bounds),
            Source::Bore(bore) => bore.classify_box(bounds),
        };
        whole_operand_result(self, other, source_is_left, relation, operation)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use multicad_rust_kernel::rounded_thread::Spec;

    fn p() -> GeometryPolicy {
        GeometryPolicy::default()
    }

    fn probe(center: [f64; 3]) -> Solid {
        Solid::box_solid([2e-4; 3], p())
            .unwrap()
            .transformed(RigidTransform {
                translation_mm: center.map(|c| c - 1e-4),
                ..Default::default()
            })
            .unwrap()
    }

    #[test]
    fn placed_thread_probes_find_teeth_core_and_clearance() {
        let angle = 27.5_f64.to_radians();
        let spec = Spec {
            major_radius_mm: 5.,
            pitch_mm: 1.,
            length_mm: 4.,
            depth_mm: 1. / (2. * angle.portable_tan()) - 0.274 * (1. / angle.portable_sin() - 1.),
            crest_radius_mm: 0.137,
            root_radius_mm: 0.137,
            flank_angle_deg: 27.5,
            starts: 1,
            left_hand: false,
        };
        // Axis along world +X from (1, 2, 3).
        let frame = PlaneFrame {
            origin: [1., 2., 3.],
            u: [0., 1., 0.],
            v: [0., 0., 1.],
        };
        let thread = Solid::rust_rounded_thread(spec, EntityId::new(), frame, p()).unwrap();
        let normal = frame.normal();
        let world = |[x, y, z]: [f64; 3]| -> [f64; 3] {
            std::array::from_fn(|i| {
                frame.origin[i] + frame.u[i] * x + frame.v[i] * y + normal[i] * z
            })
        };
        let core = probe(world([1., 0., 2.]));
        let hit = thread
            .boolean(&core, BooleanOperation::Intersect, p())
            .unwrap()
            .unwrap();
        assert!((hit.volume_mm3().unwrap() - 8e-12).abs() < 1e-18);
        assert!(
            core.boolean(&thread, BooleanOperation::Intersect, p())
                .unwrap()
                .is_some()
        );
        // Between the teeth at a root phase, and past the far cap.
        let root = probe(world([4.95, 0., 2.5]));
        assert!(
            thread
                .boolean(&root, BooleanOperation::Intersect, p())
                .unwrap()
                .is_none()
        );
        let kept = thread
            .boolean(
                &probe(world([1., 0., 4.5])),
                BooleanOperation::Subtract,
                p(),
            )
            .unwrap()
            .unwrap();
        assert!((kept.volume_mm3().unwrap() - thread.volume_mm3().unwrap()).abs() < 1e-9);
        // On the crest cylinder nothing is claimed.
        assert!(
            thread
                .contained_rounded_thread_boolean(
                    &probe(world([5., 0., 2.])),
                    BooleanOperation::Intersect,
                    p()
                )
                .is_none()
        );
        // The through receiver: tooth material at a root phase, none in
        // the core.
        let bore = Solid::rust_rounded_thread_bore(
            Spec {
                major_radius_mm: 5.05,
                ..spec
            },
            8.,
            EntityId::new(),
            frame,
            p(),
        )
        .unwrap();
        for (point, occupied) in [
            ([5.02, 0., 2.5], true),
            ([1., 0., 2.], false),
            ([6.5, 0., 2.], true),
        ] {
            let small = probe(world(point));
            assert_eq!(
                bore.boolean(&small, BooleanOperation::Intersect, p())
                    .unwrap()
                    .is_some(),
                occupied,
                "{point:?}"
            );
        }
        let fuzzy = GeometryPolicy {
            boolean_fuzzy_mm: 1e-6,
            ..p()
        };
        assert!(
            thread
                .contained_rounded_thread_boolean(&core, BooleanOperation::Intersect, fuzzy)
                .is_none()
        );
        assert!(
            core.contained_rounded_thread_boolean(&core, BooleanOperation::Intersect, p())
                .is_none()
        );
    }
}
