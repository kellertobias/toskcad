use super::*;
use multicad_rust_kernel::helical_thread::Thread;
fn source(equal: bool, left: bool) -> PhillipsScrew {
    let source = super::super::tests::fixture(1, 0., left);
    if !equal {
        return source;
    }
    let thread = Thread::new(multicad_rust_kernel::helical_thread::Spec {
        major_radius_mm: source.profile().end()[0],
        starts: 2,
        ..source.thread().spec()
    })
    .unwrap();
    PhillipsScrew::new(
        thread,
        *source.profile(),
        source.spec(),
        source.cavity().clone(),
        100_000,
        &|| false,
    )
    .unwrap()
}
#[test]
fn actual_source_incidence_and_named_faces_survive_both_real_shoulder_kinds() {
    for equal in [false, true] {
        for left in [false, true] {
            let source = source(equal, left);
            let policy = GeometryPolicy::default();
            let feature = EntityId::from_uuid_bits(61);
            let solid = Solid::rust_phillips_source_for_test(
                source.clone(),
                feature,
                PlaneFrame::default(),
                policy,
                &|| false,
            )
            .unwrap();
            let owner = EntityId::from_uuid_bits(62);
            let snapshot = solid.topology_snapshot(owner, 7).unwrap();
            snapshot.validate().unwrap();
            assert_eq!(snapshot.owner, owner);
            assert_eq!(snapshot.revision, 7);
            for (ordinal, role) in source.faces().into_iter().enumerate() {
                let actual = snapshot
                    .elements
                    .iter()
                    .find(|e| {
                        e.key
                            == TopologyKey {
                                kind: TopologyKind::Face,
                                ordinal: ordinal as u32,
                            }
                    })
                    .unwrap();
                assert!(
                    actual
                        .origins
                        .iter()
                        .any(|o| o.source == feature && o.role == face_name(&source, role)),
                    "{role:?} {:?}",
                    actual.origins
                );
            }
            let native = solid.helical_body().unwrap().0;
            let cutter = source.cavity().source().boundary_topology().unwrap();
            let cavity_edges = snapshot
                .elements
                .iter()
                .filter(|e| is_cavity_edge(e))
                .collect::<Vec<_>>();
            assert_eq!(cavity_edges.len(), 76);
            assert_eq!(
                snapshot
                    .elements
                    .iter()
                    .filter(|e| e.key.kind == TopologyKind::Vertex
                        && e.origins.iter().any(|o| o.role.starts_with(PREFIX)))
                    .count(),
                49
            );
            for edge in cutter.edges() {
                let origin = curve_origin(native, &cutter, edge.source_identity());
                let actual = cavity_edges
                    .iter()
                    .find(|e| e.origins.contains(&origin))
                    .unwrap();
                let uses = snapshot
                    .elements
                    .iter()
                    .filter(|e| e.key.kind == TopologyKind::Face && e.boundary.contains(&origin))
                    .count();
                assert_eq!(uses, 2);
                assert!(actual.geometry.is_none());
            }
            let head = source
                .faces()
                .iter()
                .position(|r| *r == PhiFace::HeadTop)
                .unwrap() as u32;
            let land = snapshot
                .elements
                .iter()
                .find(|e| {
                    e.key
                        == TopologyKey {
                            kind: TopologyKind::Face,
                            ordinal: head,
                        }
                })
                .unwrap();
            for edge in cutter
                .edges()
                .iter()
                .filter(|e| matches!(e.source_identity(), CurveKey::Mouth { .. }))
            {
                assert!(land.boundary.contains(&curve_origin(
                    native,
                    &cutter,
                    edge.source_identity()
                )));
            }
            assert!(
                !snapshot
                    .elements
                    .iter()
                    .any(|e| e.origins.iter().any(|o| o.role.contains("face/Mouth")))
            );
        }
    }
}
#[test]
fn directed_source_point_mapping_refuses_invalid_intervals_and_real_precision_floors() {
    let source = source(false, false);
    let cutter = source.cavity().source().boundary_topology().unwrap();
    let point = cutter.vertices()[0].point_mm_bounds();
    let policy = GeometryPolicy::default();
    let mapped = placed_point(PlaneFrame::default(), point, source.head_top_z(), policy).unwrap();
    policy.point(mapped).unwrap();
    let far = PlaneFrame {
        origin: [1_000_000.1; 3],
        ..PlaneFrame::default()
    };
    assert!(
        placed_point(
            far,
            point,
            source.head_top_z(),
            GeometryPolicy {
                linear_tolerance_mm: 1e-12,
                ..policy
            }
        )
        .is_err()
    );
    assert!(
        placed_point(
            PlaneFrame::default(),
            [[f64::NAN, f64::NAN]; 3],
            source.head_top_z(),
            policy
        )
        .is_err()
    );
    assert!(
        placed_point(
            PlaneFrame {
                origin: [policy.maximum_coordinate_mm; 3],
                ..PlaneFrame::default()
            },
            point,
            source.head_top_z(),
            policy
        )
        .is_err()
    );
}

#[test]
fn stale_helical_metadata_is_rejected_before_any_cavity_handoff() {
    let source = source(false, false);
    let policy = GeometryPolicy::default();
    let solid = Solid::rust_phillips_source_for_test(
        source.clone(),
        EntityId::from_uuid_bits(81),
        PlaneFrame::default(),
        policy,
        &|| false,
    )
    .unwrap();
    let (body, frame) = solid.helical_body().unwrap();
    let snapshot = body.topology(frame, body.feature, 0).unwrap();
    let base_count = 2 + body.rails.len();
    let original = snapshot
        .elements
        .iter()
        .filter(|e| e.key.kind != TopologyKind::Face || e.key.ordinal < (base_count as u32))
        .cloned()
        .collect::<Vec<_>>();
    let mut missing = original.clone();
    missing.retain(|e| {
        e.key
            != TopologyKey {
                kind: TopologyKind::Face,
                ordinal: 0,
            }
    });
    let before = missing.clone();
    assert!(extend(body, &source, frame, policy, &mut missing).is_err());
    assert_eq!(missing, before);
    let mut stale = original;
    stale
        .iter_mut()
        .find(|e| e.key.kind == TopologyKind::Face)
        .unwrap()
        .key
        .ordinal = u32::MAX;
    let before = stale.clone();
    assert!(extend(body, &source, frame, policy, &mut stale).is_err());
    assert_eq!(stale, before);
}
