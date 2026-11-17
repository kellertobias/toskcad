//! Conservative whole-prism relations against an actual unpierced sweep.
use super::*;
use crate::geometry::helical_body::contained_boolean::{
    chart_bounds, probe_bounds, whole_operand_result,
};

impl Solid {
    /// A finite rectangular or polygon prism wholly inside or wholly outside
    /// the swept material reuses one unchanged operand. Pierced sweeps and
    /// probes reaching any wall, hole or cap remain with other providers.
    pub(in crate::geometry) fn contained_sweep_boolean(
        &self,
        other: &Self,
        operation: BooleanOperation,
        policy: GeometryPolicy,
    ) -> Option<Result<Option<Self>, GeometryError>> {
        if policy.boolean_fuzzy_mm != 0. || self.policy != policy || other.policy != policy {
            return None;
        }
        let (probe, source_is_left, body, frame) =
            if let Some((body, frame)) = self.normal_sweep_body() {
                (other, true, body, frame)
            } else {
                let (body, frame) = other.normal_sweep_body()?;
                (self, false, body, frame)
            };
        body.openings.is_none().then_some(())?;
        // The sweep's own chart: profile in u/v, material along the signed
        // extrusion direction from z = 0 to its height.
        let normal = frame.normal();
        let axis = if body.reverse {
            normal.map(|n| -n)
        } else {
            normal
        };
        let bounds = chart_bounds(probe_bounds(probe)?, frame.origin, [frame.u, frame.v, axis]);
        whole_operand_result(
            self,
            other,
            source_is_left,
            body.sweep.classify_box(bounds),
            operation,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{CurvedEdge2d, CurvedProfile};

    fn p() -> GeometryPolicy {
        GeometryPolicy::default()
    }

    /// A 10 x 10 plate with a radius-2 hole, extruded from a frame whose
    /// normal is world -Y, so the chart is not world aligned.
    fn plate(height: f64) -> Solid {
        let corners = [[0., 0.], [10., 0.], [10., 10.], [0., 10.]];
        let profile = CurvedProfile {
            frame: PlaneFrame {
                origin: [1., 2., 3.],
                u: [1., 0., 0.],
                v: [0., 0., 1.],
            },
            outer: (0..4)
                .map(|i| CurvedEdge2d::Line {
                    a: corners[i],
                    b: corners[(i + 1) % 4],
                })
                .collect(),
            holes: vec![vec![CurvedEdge2d::Circle {
                center: [5., 5.],
                radius_mm: 2.,
            }]],
        };
        let feature = EntityId::new();
        let origins: Vec<_> = (0..5)
            .map(|i| Origin {
                source: feature,
                role: format!("edge/{i}"),
            })
            .collect();
        Solid::rust_extrude_normal_profile_named(&profile, height, feature, &origins, p()).unwrap()
    }

    fn probe(center: [f64; 3], side: f64) -> Solid {
        Solid::box_solid([side; 3], p())
            .unwrap()
            .transformed(RigidTransform {
                translation_mm: center.map(|c| c - side / 2.),
                ..Default::default()
            })
            .unwrap()
    }

    #[test]
    fn probes_inside_and_outside_both_extrusion_directions() {
        for height in [4., -4.] {
            let sweep = plate(height);
            assert!(sweep.normal_sweep_body().is_some());
            // World y runs along -normal: positive extrusion occupies y in [-2, 2].
            let y = if height > 0. { -0. } else { 4. };
            let inside = probe([2., y, 4.], 0.01);
            let hit = sweep
                .boolean(&inside, BooleanOperation::Intersect, p())
                .unwrap()
                .unwrap();
            assert!((hit.volume_mm3().unwrap() - 1e-6).abs() < 1e-15);
            // Probe ordering does not matter.
            assert!(
                inside
                    .boolean(&sweep, BooleanOperation::Intersect, p())
                    .unwrap()
                    .is_some()
            );
            // Within the hole and past the far cap.
            for point in [[6., y, 8.], [2., y - 3. * height.signum(), 4.]] {
                assert!(
                    sweep
                        .boolean(&probe(point, 0.01), BooleanOperation::Intersect, p())
                        .unwrap()
                        .is_none(),
                    "{height} {point:?}"
                );
            }
        }
    }

    #[test]
    fn probes_on_a_boundary_or_pierced_sweeps_are_not_claimed() {
        let sweep = plate(4.);
        // Straddling the hole wall at local (7, 5).
        let wall = probe([8., -1., 8.], 0.01);
        assert!(
            sweep
                .contained_sweep_boolean(&wall, BooleanOperation::Intersect, p())
                .is_none()
        );
        // Fuzzy and mismatched policies are left alone.
        let fuzzy = GeometryPolicy {
            boolean_fuzzy_mm: 1e-6,
            ..p()
        };
        let inside = probe([2., -1., 4.], 0.01);
        assert!(
            sweep
                .contained_sweep_boolean(&inside, BooleanOperation::Intersect, fuzzy)
                .is_none()
        );
        // Neither operand is a sweep.
        assert!(
            inside
                .contained_sweep_boolean(&inside, BooleanOperation::Intersect, p())
                .is_none()
        );
        // Subtracting an outside probe keeps the sweep; an inside probe in a
        // union keeps the sweep too.
        let outside = probe([6., -1., 8.], 0.01);
        let kept = sweep
            .boolean(&outside, BooleanOperation::Subtract, p())
            .unwrap()
            .unwrap();
        assert!((kept.volume_mm3().unwrap() - sweep.volume_mm3().unwrap()).abs() < 1e-9);
        let union = sweep
            .boolean(&inside, BooleanOperation::Union, p())
            .unwrap()
            .unwrap();
        assert!((union.volume_mm3().unwrap() - sweep.volume_mm3().unwrap()).abs() < 1e-9);
    }
}
