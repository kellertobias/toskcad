use super::*;
use multicad_core::references::{GeometryIntent, Resolution, resolve};

fn profile(starts: u32, taper: f64, left_hand: bool) -> HelicalThread {
    HelicalThread {
        major_radius_mm: 4.,
        pitch_mm: 1.25,
        length_mm: 1.25,
        depth_mm: 0.5,
        crest_width_mm: 0.25,
        lower_flank_deg: 12.,
        upper_flank_deg: 18.,
        radial_taper: taper,
        starts,
        left_hand,
    }
}
fn build(profile: HelicalThread, feature: EntityId, frame: PlaneFrame) -> Solid {
    Solid::rust_retained_helical_thread_with_interrupt(
        profile,
        feature,
        frame,
        GeometryPolicy::default(),
        None,
    )
    .unwrap()
}
#[test]
fn authored_transport_keeps_original_law_and_immutable_source_through_occurrences() {
    let feature = EntityId::new();
    for (starts, taper, left) in [(1, 0., false), (2, 0.03125, true), (8, -0.03125, false)] {
        let law = profile(starts, taper, left);
        let construction = PlaneFrame {
            origin: [3., 4., 5.],
            u: [0., 1., 0.],
            v: [0., 0., 1.],
        };
        let solid = build(law, feature, construction);
        let (source, occurrence) = solid.rust_retained_helical_thread_source().unwrap();
        assert_eq!(source.profile, law);
        assert_eq!(source.transport.radial_taper(), taper);
        assert_eq!(source.construction_frame, construction);
        assert_eq!(source.band_faces.len(), starts as usize);
        assert_eq!(source.gap_faces.len(), starts as usize);
        assert_eq!(occurrence, PlaneFrame::default());
        let copy = solid.share().unwrap();
        assert!(std::ptr::eq(
            source,
            copy.rust_retained_helical_thread_source().unwrap().0
        ));
        let posed = solid
            .transformed(RigidTransform {
                rotation: [[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]],
                translation_mm: [10., 20., 30.],
            })
            .unwrap();
        let (moved, frame) = posed.rust_retained_helical_thread_source().unwrap();
        assert!(std::ptr::eq(source, moved));
        assert_eq!(frame.origin, [10., 20., 30.]);
        let mass = posed.rust_radial_boundary_properties(None).unwrap();
        assert!(mass.volume_mm3[0] > 0.);
        let scalar = posed.volume_mm3().unwrap();
        assert!(mass.volume_mm3[0] <= scalar && scalar <= mass.volume_mm3[1]);
    }
    let blind = Solid::from_rust_radial_boundary(
        build(profile(1, 0.03125, false), feature, PlaneFrame::default())
            .rust_radial_boundary()
            .unwrap()
            .0
            .clone(),
        feature,
        GeometryPolicy::default(),
    )
    .unwrap();
    assert!(blind.rust_retained_helical_thread_source().is_none());
}
#[test]
fn authored_refs_follow_roles_across_law_and_placement_edits() {
    let feature = EntityId::new();
    let owner = EntityId::new();
    let before = build(profile(1, 0.03125, false), feature, PlaneFrame::default());
    let old = before.topology_snapshot(owner, 0).unwrap();
    let roles = [
        (TopologyKind::Face, "cap.end"),
        (TopologyKind::Face, "thread/start/0/crest"),
        (TopologyKind::Edge, "thread/start/0/crest/leading_rail"),
        (TopologyKind::Vertex, "cap.start/corner/start/0/crest"),
    ];
    let intents = roles.map(|(kind, role)| {
        GeometryIntent::capture(
            &old,
            old.elements
                .iter()
                .find(|e| {
                    e.key.kind == kind
                        && e.origins.contains(&Origin {
                            source: feature,
                            role: role.into(),
                        })
                })
                .unwrap()
                .key,
        )
        .unwrap()
    });
    let mut next_law = profile(2, -0.03125, true);
    next_law.length_mm = 1.5;
    let next = build(
        next_law,
        feature,
        PlaneFrame {
            origin: [20., 30., 40.],
            ..Default::default()
        },
    );
    let snapshot = next.topology_snapshot(owner, 1).unwrap();
    for intent in intents {
        assert!(matches!(
            resolve(&intent, &snapshot, 1),
            Resolution::Resolved { .. }
        ));
    }
    let topology = next.topology_with_display(owner, 1, Some(0.05)).unwrap();
    assert_eq!(topology.edges.len(), 28);
    assert_eq!(topology.vertices.len(), 20);
}

#[test]
fn authored_source_manufactures_and_exchanges_without_guessing_law_on_import() {
    let solid = build(
        profile(2, 0.03125, true),
        EntityId::new(),
        PlaneFrame::default(),
    );
    let receipt = solid
        .manufacturing_mesh_with_receipt(MeshOptions {
            linear_deflection_mm: 0.05,
            angular_deflection_rad: 0.4,
            maximum_triangles: 200000,
        })
        .unwrap();
    assert!(receipt.edge_pairs > 0 && receipt.vertex_links > 0);
    let file = tempfile::NamedTempFile::new().unwrap();
    solid.write_step(file.path()).unwrap();
    let text = std::fs::read_to_string(file.path()).unwrap();
    let proof = solid
        .verify_retained_step_source(&text, 1e-7)
        .unwrap()
        .unwrap();
    assert_eq!(proof.retained_members.len(), 1);
    assert!(proof.retained_members[0].decoded_volume_bounds_mm3[0] > 0.);
    let decoded =
        Solid::import_rust_step(&text, EntityId::new(), GeometryPolicy::default()).unwrap();
    assert!(decoded.rust_radial_boundary().is_some());
    assert!(decoded.rust_retained_helical_thread_source().is_none());
}

#[test]
fn authored_factory_validates_identity_policy_law_frame_and_cancellation() {
    use std::cell::Cell;
    let law = profile(1, 0.03125, false);
    let feature = EntityId::new();
    fn run(
        law: HelicalThread,
        feature: EntityId,
        frame: PlaneFrame,
        policy: GeometryPolicy,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> GeometryErrorKind {
        Solid::rust_retained_helical_thread_with_interrupt(law, feature, frame, policy, interrupt)
            .err()
            .unwrap()
            .kind
    }
    assert_eq!(
        run(
            law,
            EntityId::nil(),
            PlaneFrame::default(),
            GeometryPolicy::default(),
            None
        ),
        GeometryErrorKind::Input
    );
    assert_eq!(
        run(
            law,
            feature,
            PlaneFrame::default(),
            GeometryPolicy {
                boolean_fuzzy_mm: 1e-8,
                ..Default::default()
            },
            None
        ),
        GeometryErrorKind::Input
    );
    assert_eq!(
        run(
            law,
            feature,
            PlaneFrame::default(),
            GeometryPolicy {
                linear_tolerance_mm: f64::NAN,
                ..Default::default()
            },
            None
        ),
        GeometryErrorKind::Input
    );
    assert_eq!(
        run(
            law,
            feature,
            PlaneFrame {
                u: [0.; 3],
                ..Default::default()
            },
            GeometryPolicy::default(),
            None
        ),
        GeometryErrorKind::Input
    );
    assert_eq!(
        run(
            HelicalThread {
                pitch_mm: 0.,
                ..law
            },
            feature,
            PlaneFrame::default(),
            GeometryPolicy::default(),
            None
        ),
        GeometryErrorKind::Input
    );
    assert_eq!(
        run(
            HelicalThread {
                radial_taper: 0.5,
                ..law
            },
            feature,
            PlaneFrame::default(),
            GeometryPolicy::default(),
            None
        ),
        GeometryErrorKind::Input
    );
    for after in [0, 2, 80] {
        let count = Cell::new(0);
        let cancel = || {
            let n = count.get();
            count.set(n + 1);
            n >= after
        };
        assert_eq!(
            run(
                law,
                feature,
                PlaneFrame::default(),
                GeometryPolicy::default(),
                Some(&cancel)
            ),
            GeometryErrorKind::Cancelled
        );
        assert_eq!(count.get(), after + 1);
    }
    let solid = build(law, feature, PlaneFrame::default());
    assert!(solid.rust_retained_helical_thread_source().is_some());
}

#[test]
fn authored_ownership_rejects_incomplete_duplicate_incompatible_or_unconnected_receipts() {
    let law = profile(1, 0.03125, false);
    let solid = build(law, EntityId::new(), PlaneFrame::default());
    let (source, _) = solid.rust_retained_helical_thread_source().unwrap();
    let authored = AuthoredBoundary {
        boundary: solid.rust_radial_boundary().unwrap().0.clone(),
        cap_faces: source.cap_faces,
        band_faces: source.band_faces.clone(),
        gap_faces: source.gap_faces.clone(),
    };
    let make =
        |a: &AuthoredBoundary| AuthoredThread::new(law, source.transport, PlaneFrame::default(), a);
    assert!(make(&authored).is_ok());
    for edit in 0..5 {
        let mut raw = authored.clone();
        match edit {
            0 => raw.band_faces.clear(),
            1 => raw.gap_faces.clear(),
            2 => raw.cap_faces[0] = usize::MAX,
            3 => raw.cap_faces[0] = raw.cap_faces[1],
            4 => std::mem::swap(&mut raw.cap_faces[0], &mut raw.band_faces[0][0]),
            _ => unreachable!(),
        }
        assert!(make(&raw).is_err());
    }
    let owned = make(&authored).unwrap();
    let input = authored.boundary.input();
    let mut incidence = vec![BTreeSet::new(); input.edges.len()];
    let vertices = vec![BTreeSet::new(); input.vertices.len()];
    assert!(AuthoredThread::edge_role(input, 0, &incidence, &vertices, &owned.faces, 1).is_err());
    // Raw malformed adjacency is never promoted to a semantic rail or rim.
    for adjacent in [
        [FaceRole::Cap(0), FaceRole::Cap(1)],
        [FaceRole::Cap(0), FaceRole::Band { start: 0, band: 3 }],
        [
            FaceRole::Band { start: 0, band: 0 },
            FaceRole::Band { start: 0, band: 2 },
        ],
        [
            FaceRole::Band { start: 0, band: 1 },
            FaceRole::Band { start: 0, band: 3 },
        ],
        [
            FaceRole::Band { start: 0, band: 2 },
            FaceRole::Band { start: 1, band: 3 },
        ],
    ] {
        incidence[0] = BTreeSet::from([0, 1]);
        assert!(AuthoredThread::edge_role(input, 0, &incidence, &vertices, &adjacent, 2).is_err());
    }
    assert!(
        AuthoredThread::vertex_role(&[], &BTreeSet::new(), &owned.faces, &owned.edges).is_err()
    );
    assert!(
        AuthoredThread::vertex_role(
            &[],
            &BTreeSet::from([authored.cap_faces[0]]),
            &owned.faces,
            &owned.edges
        )
        .is_err()
    );
    assert!(
        AuthoredThread::vertex_role(
            &[0],
            &BTreeSet::from([0, 1]),
            &[FaceRole::Band { start: 0, band: 0 }, FaceRole::Cap(0)],
            &[EdgeRole::Rail { start: 0, band: 0 }]
        )
        .is_ok()
    );
    assert!(
        AuthoredThread::vertex_role(
            &[0],
            &BTreeSet::from([0]),
            &[FaceRole::Cap(0)],
            &[EdgeRole::Rim {
                end: 0,
                start: 0,
                band: 0,
                segment: None
            }]
        )
        .is_err()
    );
    let mut face = TopologyElement {
        key: multicad_core::references::TopologyKey {
            kind: TopologyKind::Face,
            ordinal: 0,
        },
        origins: BTreeSet::new(),
        geometry: None,
        ..solid
            .topology_snapshot(EntityId::new(), 0)
            .unwrap()
            .elements[0]
            .clone()
    };
    owned.apply_origin(EntityId::new(), &mut face);
    assert!(face.origins.is_empty());
}
