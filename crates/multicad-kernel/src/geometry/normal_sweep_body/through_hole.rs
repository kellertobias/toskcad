//! Cut a straight cylinder through both caps of an unpierced normal sweep.
use super::*;

impl Solid {
    /// A cylinder parallel to the sweep direction that spans both caps and
    /// lies strictly inside the profile material becomes one more circular
    /// hole ring of the same exact sweep. Every existing side keeps its
    /// authored identity; the new wall carries the cylinder's.
    pub(in crate::geometry) fn subtract_cylinder_through_sweep(
        &self,
        tool: &Self,
        operation: BooleanOperation,
        policy: GeometryPolicy,
    ) -> Option<Result<Option<Self>, GeometryError>> {
        if operation != BooleanOperation::Subtract
            || policy.boolean_fuzzy_mm != 0.
            || self.policy != policy
            || tool.policy != policy
        {
            return None;
        }
        let (body, placement) = self.normal_sweep_body()?;
        let (cylinder, frame) = tool.axial_cylinder_frame()?;
        if body.openings.is_some() || body.stock_history.is_some() {
            return None;
        }
        let normal = placement.normal();
        let axis = if body.reverse {
            normal.map(|n| -n)
        } else {
            normal
        };
        let direction = dot(frame.normal(), axis);
        let cross_axis = cross(frame.normal(), axis);
        let tilt = dot(cross_axis, cross_axis).sqrt();
        let delta: [f64; 3] = std::array::from_fn(|i| frame.origin[i] - placement.origin[i]);
        let height = body.sweep.height_mm();
        let length = cylinder.height().as_f64();
        let magnitude = delta
            .into_iter()
            .chain(placement.origin)
            .map(f64::abs)
            .fold(height.max(length).max(1.), f64::max);
        let rounding = 64. * f64::EPSILON * magnitude;
        if direction.is_nan() || direction == 0. || tilt.is_nan() || tilt > 64. * f64::EPSILON {
            return None;
        }
        // The cylinder's axial interval in the sweep chart must cover both
        // caps; its circle is then the same at every height.
        let start = dot(delta, axis);
        let end = start + direction.signum() * length;
        if start.min(end) > rounding || start.max(end) < height - rounding {
            return None;
        }
        let center = [dot(delta, placement.u), dot(delta, placement.v)];
        let radius = cylinder.circle().radius().as_f64();
        let sweep = body
            .sweep
            .with_circular_hole(
                center,
                radius,
                policy.linear_tolerance_mm,
                policy.minimum_feature_mm,
            )
            .ok()?;
        let origin = tool.cylinder_sketch().map_or_else(
            || Origin {
                source: EntityId::from_uuid_bits(cylinder.feature()),
                role: "side.cylinder".into(),
            },
            |sketch| sketch.origin.clone(),
        );
        if body.sides.iter().any(|side| side.origin == origin) {
            return None;
        }
        let mut drilled = body.clone();
        let index = drilled.sides.len();
        drilled.sides.push(Side {
            edge: *sweep.rings().last()?.first()?,
            origin,
            reversed: true,
            previous: index,
            next: index,
            rail: None,
        });
        drilled.rings.push(index..index + 1);
        // Straight sides own rails after the two rims of every side.
        let mut rail = 2 * drilled.sides.len() as u32;
        for side in &mut drilled.sides {
            if side.rail.is_some() {
                side.rail = Some(rail);
                rail += 1;
            }
        }
        drilled.sweep = sweep;
        Some(Ok(Some(Self {
            body: SolidBody::NormalSweep(Rc::new(drilled), placement),
            policy,
            provenance: self.provenance.clone(),
            rust_origins: self.rust_origins.clone(),
            retained_placement: None,
        })))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{CurvedEdge2d, CurvedProfile};
    use multicad_core::references::TopologyKind;

    fn p() -> GeometryPolicy {
        GeometryPolicy::default()
    }

    /// A 100 x 60 plate with a radius-10 hole, 10 high along +Z or -Z.
    fn plate(height: f64) -> Solid {
        let corners = [[0., 0.], [100., 0.], [100., 60.], [0., 60.]];
        let profile = CurvedProfile {
            frame: PlaneFrame::default(),
            outer: (0..4)
                .map(|i| CurvedEdge2d::Line {
                    a: corners[i],
                    b: corners[(i + 1) % 4],
                })
                .collect(),
            holes: vec![vec![CurvedEdge2d::Circle {
                center: [30., 30.],
                radius_mm: 10.,
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

    fn cylinder(center: [f64; 2], z: f64, radius: f64, height: f64) -> Solid {
        Solid::cylinder(
            PlaneFrame {
                origin: [center[0], center[1], z],
                ..Default::default()
            },
            radius,
            height,
            p(),
        )
        .unwrap()
    }

    #[test]
    fn a_through_cylinder_becomes_an_exact_named_hole() {
        for height in [10., -10.] {
            let source = plate(height);
            let z = height.min(0.);
            let tool = cylinder([85., 10.], z, 1.5, 10.);
            let cut = source
                .boolean(&tool, BooleanOperation::Subtract, p())
                .unwrap()
                .unwrap();
            let disc = std::f64::consts::PI * 2.25 * 10.;
            assert!((source.volume_mm3().unwrap() - cut.volume_mm3().unwrap() - disc).abs() < 1e-9);
            assert_eq!(
                cut.faces().unwrap().len(),
                source.faces().unwrap().len() + 1
            );
            // The route itself names the new wall from the cutter; a bare
            // primitive cutter's private identity is then withheld by the
            // public Boolean like every other authored-source retention.
            let direct = source
                .subtract_cylinder_through_sweep(&tool, BooleanOperation::Subtract, p())
                .unwrap()
                .unwrap()
                .unwrap();
            let snapshot = direct.topology_snapshot(EntityId::new(), 0).unwrap();
            let walls = snapshot
                .elements
                .iter()
                .filter(|e| e.key.kind == TopologyKind::Face)
                .filter(|e| {
                    e.origins.iter().any(|o| {
                        o.role.ends_with("primitive.circle") || o.role.ends_with("side.cylinder")
                    })
                })
                .count();
            assert_eq!(walls, 1);
            // A longer cylinder protruding past both caps cuts the same hole.
            let long = cylinder([85., 10.], z - 5., 1.5, 20.);
            let again = source
                .boolean(&long, BooleanOperation::Subtract, p())
                .unwrap()
                .unwrap();
            assert!((again.volume_mm3().unwrap() - cut.volume_mm3().unwrap()).abs() < 1e-9);
        }
    }

    #[test]
    fn blind_crossing_or_misaligned_cylinders_are_not_holes() {
        let source = plate(10.);
        let declined = |tool: &Solid| {
            source
                .subtract_cylinder_through_sweep(tool, BooleanOperation::Subtract, p())
                .is_none()
        };
        // Blind from either cap.
        assert!(declined(&cylinder([85., 10.], 2., 1.5, 20.)));
        assert!(declined(&cylinder([85., 10.], -5., 1.5, 10.)));
        // Touching the outer wall, overlapping the existing hole, outside.
        assert!(declined(&cylinder([98.5, 10.], 0., 1.5, 10.)));
        assert!(declined(&cylinder([30., 41.], 0., 1.5, 10.)));
        assert!(declined(&cylinder([130., 10.], 0., 1.5, 10.)));
        // Not parallel to the sweep.
        let sideways = Solid::cylinder(
            PlaneFrame {
                origin: [50., 0., 5.],
                u: [1., 0., 0.],
                v: [0., 0., 1.],
            },
            1.,
            60.,
            p(),
        )
        .unwrap();
        assert!(declined(&sideways));
        // Only subtraction, with the exact policy, of a cylinder.
        let tool = cylinder([85., 10.], 0., 1.5, 10.);
        assert!(
            source
                .subtract_cylinder_through_sweep(&tool, BooleanOperation::Intersect, p())
                .is_none()
        );
        let fuzzy = GeometryPolicy {
            boolean_fuzzy_mm: 1e-6,
            ..p()
        };
        assert!(
            source
                .subtract_cylinder_through_sweep(&tool, BooleanOperation::Subtract, fuzzy)
                .is_none()
        );
        assert!(
            source
                .subtract_cylinder_through_sweep(&source, BooleanOperation::Subtract, p())
                .is_none()
        );
        // The same cutter twice cannot add a second wall with its identity,
        // while a different cutter can.
        let once = source
            .subtract_cylinder_through_sweep(&tool, BooleanOperation::Subtract, p())
            .unwrap()
            .unwrap()
            .unwrap();
        let moved = tool
            .transformed(RigidTransform {
                translation_mm: [-15., 0., 0.],
                ..Default::default()
            })
            .unwrap();
        if moved.cylinder_sketch().is_some() || moved.axial_cylinder_frame().is_some() {
            assert!(
                once.subtract_cylinder_through_sweep(&moved, BooleanOperation::Subtract, p())
                    .is_none()
            );
        }
        let other = cylinder([70., 10.], 0., 1.5, 10.);
        let twice = once
            .subtract_cylinder_through_sweep(&other, BooleanOperation::Subtract, p())
            .unwrap()
            .unwrap()
            .unwrap();
        assert_eq!(
            twice.faces().unwrap().len(),
            once.faces().unwrap().len() + 1
        );
    }
}
