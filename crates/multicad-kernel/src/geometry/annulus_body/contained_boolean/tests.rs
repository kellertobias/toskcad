use super::*;
use crate::{GeometryErrorKind, RigidTransform};

fn probe(frame: PlaneFrame, point: [f64; 3], policy: GeometryPolicy) -> Solid {
    Solid::box_solid([0.01; 3], policy)
        .unwrap()
        .transformed(RigidTransform {
            translation_mm: position(frame, point).map(|x| x - 0.005),
            ..Default::default()
        })
        .unwrap()
}
#[test]
fn retained_annular_material_bore_and_exterior_keep_exact_operand_results() {
    let policy = GeometryPolicy::default();
    for posed in [false, true] {
        for reverse in [false, true] {
            let profile = if posed {
                PlaneFrame {
                    origin: [11., -7., 5.],
                    u: [0.8, 0.6, 0.],
                    v: [-0.36, 0.48, 0.8],
                }
            } else {
                PlaneFrame::default()
            };
            let owner = EntityId::new();
            let solid = Solid::rust_annulus_from_sketch(
                [2., 4.],
                if reverse { -3. } else { 3. },
                owner,
                profile,
                [0.7, -0.3],
                ["inner", "outer"].map(|role| Origin {
                    source: owner,
                    role: role.into(),
                }),
                policy,
            )
            .unwrap();
            let (_, frame) = solid.annulus_body().unwrap();
            let volume = solid.volume_mm3().unwrap();
            for (point, inside) in [
                ([3., 0.1, 1.5], true),
                ([0.1, 0.1, 1.5], false),
                ([4.5, 0.1, 1.5], false),
                ([3., 0.1, -0.1], false),
                ([3., 0.1, 3.1], false),
            ] {
                let small = probe(frame, point, policy);
                let small_volume = small.volume_mm3().unwrap();
                for swapped in [false, true] {
                    let (a, b) = if swapped {
                        (&small, &solid)
                    } else {
                        (&solid, &small)
                    };
                    for operation in [
                        BooleanOperation::Intersect,
                        BooleanOperation::Union,
                        BooleanOperation::Subtract,
                    ] {
                        let expected = match (inside, operation, swapped) {
                            (true, BooleanOperation::Intersect, _) => Some(Some(small_volume)),
                            (true, BooleanOperation::Union, _) => Some(Some(volume)),
                            (true, BooleanOperation::Subtract, true)
                            | (false, BooleanOperation::Intersect, _) => Some(None),
                            (false, BooleanOperation::Subtract, _) => {
                                Some(Some(if swapped { small_volume } else { volume }))
                            }
                            _ => None,
                        };
                        if let Some(expected) = expected {
                            let result = a.boolean(b, operation, policy).unwrap();
                            assert_eq!(result.is_some(), expected.is_some());
                            if let Some(expected) = expected {
                                assert_eq!(result.unwrap().volume_mm3().unwrap(), expected);
                            }
                        } else {
                            assert!(a.contained_annulus_boolean(b, operation, policy).is_none());
                        }
                    }
                }
                assert_eq!(solid.volume_mm3().unwrap(), volume);
                assert_eq!(small.volume_mm3().unwrap(), small_volume);
            }
            for point in [[2., 0., 1.5], [4., 0., 1.5], [3., 0., 0.], [3., 0., 3.]] {
                let crossing = probe(frame, point, policy);
                assert!(
                    solid
                        .contained_annulus_boolean(&crossing, BooleanOperation::Intersect, policy)
                        .is_none()
                );
                assert_eq!(
                    solid
                        .boolean(&crossing, BooleanOperation::Intersect, policy)
                        .err()
                        .unwrap()
                        .kind,
                    GeometryErrorKind::Unsupported
                );
            }
            let small = probe(frame, [3., 0., 1.5], policy);
            for changed in [
                GeometryPolicy {
                    boolean_fuzzy_mm: 1e-8,
                    ..policy
                },
                GeometryPolicy {
                    linear_tolerance_mm: 1e-6,
                    ..policy
                },
            ] {
                assert!(
                    solid
                        .contained_annulus_boolean(&small, BooleanOperation::Intersect, changed)
                        .is_none()
                );
                let other = probe(frame, [3., 0., 1.5], changed);
                assert!(
                    solid
                        .contained_annulus_boolean(&other, BooleanOperation::Intersect, policy)
                        .is_none()
                );
            }
            let sphere = Solid::rust_sphere_named(0.01, frame, EntityId::new(), policy).unwrap();
            assert!(
                solid
                    .contained_annulus_boolean(&sphere, BooleanOperation::Intersect, policy)
                    .is_none()
            );
            assert!(
                small
                    .contained_annulus_boolean(&sphere, BooleanOperation::Intersect, policy)
                    .is_none()
            );
        }
    }
}
#[test]
fn invalid_and_overflowing_annular_queries_remain_unresolved() {
    let source = Annulus::new(2., 4., 3.).unwrap();
    for bounds in [
        [[f64::NAN, 1.]; 3],
        [[1., 0.]; 3],
        [[-f64::MAX, f64::MAX]; 3],
    ] {
        assert_eq!(
            relation(&source, bounds),
            BoxClassification::UnresolvedBoundary
        );
    }
}
