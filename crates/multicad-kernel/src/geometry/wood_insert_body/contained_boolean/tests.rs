use super::*;
use crate::RigidTransform;

fn probe(frame: PlaneFrame, center: [f64; 3], policy: GeometryPolicy) -> Solid {
    let n = frame.normal();
    let world = std::array::from_fn(|i| {
        frame.origin[i] + frame.u[i] * center[0] + frame.v[i] * center[1] + n[i] * center[2]
    });
    Solid::box_solid([0.02; 3], policy)
        .unwrap()
        .transformed(RigidTransform {
            translation_mm: world.map(|x| x - 0.01),
            ..Default::default()
        })
        .unwrap()
}

#[test]
fn paired_insert_material_bore_and_exterior_use_actual_whole_operand_proofs() {
    let policy = GeometryPolicy::default();
    for left in [false, true] {
        for frame in [
            PlaneFrame::default(),
            PlaneFrame {
                origin: [11., -7., 5.],
                u: [0.8, 0.6, 0.],
                v: [-0.36, 0.48, 0.8],
            },
            PlaneFrame {
                origin: [1e4, -2e4, 3e4],
                u: [0., 1., 0.],
                v: [0., 0., 1.],
            },
        ] {
            let source = crate::geometry::wood_insert_world_tests::source(left);
            let insert =
                Solid::rust_wood_insert(source, EntityId::new(), frame, policy, None).unwrap();
            let before = insert.rust_wood_insert_data();
            for (center, inside) in [
                ([3., 0., 1.], true),
                ([0., 0., 1.], false),
                ([6., 0., 1.], false),
            ] {
                let box_solid = probe(frame, center, policy);
                for (a, b) in [(&insert, &box_solid), (&box_solid, &insert)] {
                    let intersection = a.boolean(b, BooleanOperation::Intersect, policy).unwrap();
                    assert_eq!(intersection.is_some(), inside);
                    if let Some(retained) = intersection {
                        assert_eq!(retained.bounds().unwrap(), box_solid.bounds().unwrap());
                        assert_eq!(
                            retained.volume_mm3().unwrap(),
                            box_solid.volume_mm3().unwrap()
                        );
                    }
                }
                if inside {
                    assert!(
                        box_solid
                            .boolean(&insert, BooleanOperation::Subtract, policy)
                            .unwrap()
                            .is_none()
                    );
                    for (a, b) in [(&insert, &box_solid), (&box_solid, &insert)] {
                        let retained = a
                            .boolean(b, BooleanOperation::Union, policy)
                            .unwrap()
                            .unwrap();
                        assert_eq!(retained.rust_wood_insert_data(), before);
                    }
                } else {
                    for (a, b) in [(&insert, &box_solid), (&box_solid, &insert)] {
                        let retained = a
                            .boolean(b, BooleanOperation::Subtract, policy)
                            .unwrap()
                            .unwrap();
                        assert_eq!(retained.bounds().unwrap(), a.bounds().unwrap());
                        assert_eq!(retained.volume_mm3().unwrap(), a.volume_mm3().unwrap());
                    }
                }
            }
            let crossing = probe(frame, [3., 0., 0.], policy);
            assert!(
                insert
                    .contained_insert_boolean(&crossing, BooleanOperation::Intersect, policy)
                    .is_none()
            );
            assert!(
                insert
                    .boolean(&crossing, BooleanOperation::Intersect, policy)
                    .is_err()
            );
            assert_eq!(insert.rust_wood_insert_data(), before);
        }
    }
}

#[test]
fn insert_relation_refuses_unknown_operands_changed_policy_and_required_new_boundaries() {
    let policy = GeometryPolicy::default();
    let source = crate::geometry::wood_insert_world_tests::source(false);
    let insert =
        Solid::rust_wood_insert(source, EntityId::new(), PlaneFrame::default(), policy, None)
            .unwrap();
    let inside = probe(PlaneFrame::default(), [3., 0., 1.], policy);
    let outside = probe(PlaneFrame::default(), [6., 0., 1.], policy);
    assert!(
        insert
            .contained_insert_boolean(&inside, BooleanOperation::Subtract, policy)
            .is_none()
    );
    assert!(
        insert
            .contained_insert_boolean(&outside, BooleanOperation::Union, policy)
            .is_none()
    );
    assert!(
        insert
            .contained_insert_boolean(&insert, BooleanOperation::Intersect, policy)
            .is_none()
    );
    assert!(
        inside
            .contained_insert_boolean(&outside, BooleanOperation::Intersect, policy)
            .is_none()
    );
    // Every corner is outside, but the box encloses real insert material.
    // Corner sampling must never turn that into an empty intersection.
    let enclosing = Solid::box_solid([12., 12., 3.], policy)
        .unwrap()
        .transformed(RigidTransform {
            translation_mm: [-6., -6., 0.5],
            ..Default::default()
        })
        .unwrap();
    assert!(
        insert
            .contained_insert_boolean(&enclosing, BooleanOperation::Intersect, policy)
            .is_none()
    );
    assert!(
        insert
            .boolean(&enclosing, BooleanOperation::Intersect, policy)
            .is_err()
    );
    let mut changed = policy;
    changed.boolean_fuzzy_mm = 0.001;
    assert!(
        insert
            .contained_insert_boolean(&inside, BooleanOperation::Intersect, changed)
            .is_none()
    );
    changed = policy;
    changed.linear_tolerance_mm *= 2.;
    assert!(
        insert
            .contained_insert_boolean(&inside, BooleanOperation::Intersect, changed)
            .is_none()
    );
}
