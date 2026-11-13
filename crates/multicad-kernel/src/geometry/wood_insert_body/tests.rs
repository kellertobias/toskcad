use super::*;
#[test]
fn actual_source_solid_publishes_annular_topology_faces_and_closed_manufacturing() {
    for frame in [
        PlaneFrame::default(),
        PlaneFrame {
            origin: [11., -7., 5.],
            u: [0.8, 0.6, 0.],
            v: [-0.36, 0.48, 0.8],
        },
    ] {
        let source = crate::geometry::wood_insert_world_tests::source(false);
        let feature = EntityId::new();
        let timer = std::time::Instant::now();
        let solid =
            Solid::rust_wood_insert(source, feature, frame, GeometryPolicy::default(), None)
                .unwrap();
        assert_eq!(solid.solid_count().unwrap(), 1);
        assert!(solid.volume_mm3().unwrap() > 0.);
        assert_eq!(
            solid.rust_wood_insert_data(),
            Some(WoodInsertData { source, frame })
        );
        let topology = solid
            .topology_with_display(EntityId::new(), 7, Some(0.1))
            .unwrap();
        assert_eq!(
            topology
                .topology
                .elements
                .iter()
                .filter(|e| e.key.kind == TopologyKind::Face)
                .count(),
            source.face_count() as usize
        );
        let faces = solid.faces().unwrap();
        assert_eq!(faces.len(), source.face_count() as usize);
        assert!(faces.iter().all(|f| f.area_mm2 > 0.));
        assert!(faces[0].plane.is_some() && faces[1].plane.is_some());
        assert!(!topology.edges.is_empty() && !topology.vertices.is_empty());
        let proof = solid
            .manufacturing_mesh_with_receipt(MeshOptions::default())
            .unwrap();
        assert!(proof.edge_pairs > 0 && proof.vertex_links > 0);
        assert_eq!(proof.removed_exact_zero_triangles, 0);
        assert_eq!(solid.rust_wood_insert_data().unwrap().source, source);
        eprintln!(
            "actual full insert Solid {}faces {}edges {}triangles {:?}",
            faces.len(),
            topology.edges.len(),
            proof.mesh.indices.len() / 3,
            timer.elapsed()
        );
    }
}

#[test]
fn original_source_ids_and_complete_void_survive_actual_transform_and_failures() {
    let source = crate::geometry::wood_insert_world_tests::source(true);
    let feature = EntityId::new();
    let policy = GeometryPolicy::default();
    assert_eq!(
        Solid::rust_wood_insert(source, EntityId::nil(), PlaneFrame::default(), policy, None)
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Input
    );
    assert_eq!(
        Solid::rust_wood_insert(
            source,
            feature,
            PlaneFrame::default(),
            policy,
            Some(&|| true)
        )
        .err()
        .unwrap()
        .kind,
        GeometryErrorKind::Cancelled
    );
    let shape =
        Solid::rust_wood_insert(source, feature, PlaneFrame::default(), policy, None).unwrap();
    let owner = EntityId::new();
    let original = shape.topology_snapshot(owner, 0).unwrap();
    let moved = shape
        .transformed(RigidTransform {
            translation_mm: [10., -5., 3.],
            rotation: [[0.8, -0.6, 0.], [0.6, 0.8, 0.], [0., 0., 1.]],
        })
        .unwrap();
    let placed = moved.topology_snapshot(owner, 1).unwrap();
    assert_eq!(
        original
            .elements
            .iter()
            .map(|e| &e.origins)
            .collect::<Vec<_>>(),
        placed
            .elements
            .iter()
            .map(|e| &e.origins)
            .collect::<Vec<_>>()
    );
    assert_eq!(moved.rust_wood_insert_data().unwrap().source, source);
    assert!((shape.volume_mm3().unwrap() - moved.volume_mm3().unwrap()).abs() < 1e-6);
    assert!(shape.topology_snapshot(EntityId::nil(), 0).is_err());
    assert_eq!(
        shape.display_edges(0.0001).unwrap_err().kind,
        GeometryErrorKind::Input
    );
    assert!(
        shape
            .manufacturing_mesh_with_receipt(MeshOptions {
                maximum_triangles: 32,
                ..Default::default()
            })
            .is_err()
    );
    assert_eq!(shape.rust_wood_insert_data().unwrap().source, source);
    let bounds = moved.bounds().unwrap();
    assert!(bounds.min[0] < bounds.max[0]);
    assert_eq!(
        source
            .classify_box([[0., 0.], [0., 0.], [0.5, 3.5]], None)
            .unwrap(),
        multicad_rust_kernel::wood_insert::Classification::Outside
    );
}

#[test]
fn source_failures_preserve_error_classes_and_uncertain_vertices_are_refused() {
    for (cause, kind) in [
        (source::Error::Invalid, GeometryErrorKind::Input),
        (source::Error::Limit, GeometryErrorKind::Limit),
        (source::Error::Cancelled, GeometryErrorKind::Cancelled),
        (source::Error::Unsupported, GeometryErrorKind::Unsupported),
        (source::Error::Unresolved, GeometryErrorKind::Unsupported),
    ] {
        assert_eq!(failure(cause).kind, kind);
    }
    let frame = PlaneFrame::default();
    let policy = GeometryPolicy::default();
    for point in [
        [B::new(f64::NAN, f64::NAN), B::ZERO, B::ZERO],
        [B::new(0., 1.), B::ZERO, B::ZERO],
    ] {
        assert_eq!(
            world_point(frame, point, policy).unwrap_err().kind,
            GeometryErrorKind::Unsupported
        );
    }
    assert!(
        world_point(
            PlaneFrame {
                origin: [1e100, 0., 0.],
                ..frame
            },
            [B::ZERO; 3],
            policy
        )
        .is_err()
    );
    let source = crate::geometry::wood_insert_world_tests::source(false);
    let feature = EntityId::new();
    let solid = Solid::rust_wood_insert(source, feature, frame, policy, None).unwrap();
    let (body, _) = solid.wood_insert_body().unwrap();
    assert!(body.topology(frame, EntityId::nil(), 0, policy).is_err());
    for deflection in [0., f64::NAN, f64::INFINITY] {
        assert_eq!(
            body.display_edges(frame, deflection, policy)
                .unwrap_err()
                .kind,
            GeometryErrorKind::Input
        );
    }
    let other = Solid::box_solid([1.; 3], policy).unwrap();
    assert!(other.rust_wood_insert_data().is_none());
}

#[test]
fn cancellation_during_source_graph_and_after_qualification_returns_no_solid() {
    use std::cell::Cell;
    let source = crate::geometry::wood_insert_world_tests::source(false);
    let polls = Cell::new(0usize);
    let count = || {
        polls.set(polls.get() + 1);
        false
    };
    Solid::rust_wood_insert(
        source,
        EntityId::new(),
        PlaneFrame::default(),
        GeometryPolicy::default(),
        Some(&count),
    )
    .unwrap();
    let total = polls.get();
    assert!(total > 3);
    for stop in [3, total] {
        polls.set(0);
        let cancel = || {
            polls.set(polls.get() + 1);
            polls.get() >= stop
        };
        let failed = Solid::rust_wood_insert(
            source,
            EntityId::new(),
            PlaneFrame::default(),
            GeometryPolicy::default(),
            Some(&cancel),
        );
        assert_eq!(failed.err().unwrap().kind, GeometryErrorKind::Cancelled);
    }
}

#[test]
fn semantic_endpoint_roles_do_not_alias_after_constraint_reordering() {
    use multicad_rust_kernel::authored_internal_treatment::BoundaryPointKey;
    let a = BoundaryConstraint::Rail(7);
    let b = BoundaryConstraint::Station(2, 3);
    let key = BoundaryPointKey([a, b]);
    let reversed = BoundaryPointKey([b, a]);
    assert_eq!(
        point_tag(BoundaryPointRole::Source(key)),
        point_tag(BoundaryPointRole::Source(reversed))
    );
    let roles = [
        BoundaryPointRole::Source(key),
        BoundaryPointRole::Tail {
            original: key,
            end: 0,
        },
        BoundaryPointRole::Tail {
            original: key,
            end: 1,
        },
        BoundaryPointRole::StockRim { end: 0 },
        BoundaryPointRole::StockRim { end: 1 },
    ];
    assert_eq!(
        roles
            .into_iter()
            .map(point_tag)
            .collect::<BTreeSet<_>>()
            .len(),
        roles.len()
    );
}

#[test]
fn asymmetric_multistart_outer_law_keeps_distinct_faces_and_reference_ids() {
    use multicad_rust_kernel::authored_flat_thread::{AuthoredThread, SourceOptions};
    let original = crate::geometry::wood_insert_world_tests::source(false);
    let mut outer = original.outer().spec();
    outer.starts = 2;
    outer.lower_flank_deg = 25.;
    outer.upper_flank_deg = 35.;
    let outer = AuthoredThread::new(
        outer,
        SourceOptions {
            position_limit_mm: 1e-7,
            normal_angle_limit_rad: 1e-7,
        },
    )
    .unwrap();
    let source = source::Insert::new(outer, *original.inner(), 0.5, None).unwrap();
    let solid = Solid::rust_wood_insert(
        source,
        EntityId::new(),
        PlaneFrame::default(),
        GeometryPolicy::default(),
        None,
    )
    .unwrap();
    let topology = solid.topology_snapshot(EntityId::new(), 0).unwrap();
    topology.validate().unwrap();
    let faces = solid.faces().unwrap();
    assert_eq!(faces.len(), 14);
    let origins = topology
        .elements
        .iter()
        .flat_map(|e| e.origins.iter())
        .collect::<BTreeSet<_>>();
    assert_eq!(origins.len(), topology.elements.len());
    assert!(
        solid
            .manufacturing_mesh_with_receipt(Default::default())
            .unwrap()
            .vertex_links
            > 0
    );
    assert_eq!(solid.rust_wood_insert_data().unwrap().source, source);
}
