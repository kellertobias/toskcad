use multicad_core::{
    document::EntityId,
    references::{GeometryIntent, Resolution, TopologyKind, resolve},
};
use multicad_kernel::{
    BlindThreadCut, BlindThreadSeat, GeometryPolicy, HelicalThread, PlaneFrame, Solid,
    rust_kernel::Rational,
};
fn n(v: i64) -> Rational {
    Rational::new(v, 1).unwrap()
}
fn point(f: PlaneFrame, p: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|i| f.origin[i] + f.u[i] * p[0] + f.v[i] * p[1] + f.normal()[i] * p[2])
}
fn spec(left_hand: bool) -> HelicalThread {
    HelicalThread {
        major_radius_mm: 1.5,
        pitch_mm: 0.5,
        length_mm: 3.5,
        depth_mm: 0.27,
        crest_width_mm: 0.08,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        radial_taper: 0.,
        starts: 1,
        left_hand,
    }
}
fn cut(
    stock: &Solid,
    feature: EntityId,
    entry: PlaneFrame,
    left: bool,
    depth: f64,
) -> Result<Solid, multicad_kernel::GeometryError> {
    stock.rust_cut_counterbored_blind_internal_thread(
        BlindThreadCut {
            profile: spec(left),
            entry,
            receiving_radius_mm: 1.5,
            minimum_wall_mm: 1.,
            minimum_floor_mm: 1.,
        },
        BlindThreadSeat {
            radius_mm: 3.5,
            depth_mm: depth,
        },
        feature,
        GeometryPolicy::default(),
    )
}
#[test]
fn actual_counterbored_blind_source_preserves_stock_and_finite_reference_geometry() {
    for left in [false, true] {
        for top in [false, true] {
            for frame in [
                PlaneFrame::default(),
                PlaneFrame {
                    origin: [100., -30., 50.],
                    u: [0.8, 0.6, 0.],
                    v: [-0.36, 0.48, 0.8],
                },
            ] {
                let stock_id = EntityId::new();
                let feature = EntityId::new();
                let policy = GeometryPolicy::default();
                let stock =
                    Solid::rust_box(stock_id.uuid_bits(), [n(30), n(20), n(20)], frame, policy)
                        .unwrap();
                let entry = PlaneFrame {
                    origin: point(frame, [10., 10., if top { 20. } else { 0. }]),
                    u: frame.u,
                    v: frame.v.map(|v| if top { -v } else { v }),
                };
                let original = stock.topology_snapshot(stock_id, 0).unwrap();
                let seated = cut(&stock, feature, entry, left, 0.8).unwrap();
                let snapshot = seated.topology_snapshot(stock_id, 1).unwrap();
                assert_eq!(seated.faces().unwrap().len(), 13);
                assert_eq!(
                    snapshot
                        .elements
                        .iter()
                        .filter(|e| e.key.kind == TopologyKind::Face)
                        .count(),
                    13
                );
                let shoulder = snapshot
                    .elements
                    .iter()
                    .find(|e| {
                        e.key.kind == TopologyKind::Face
                            && e.origins
                                .iter()
                                .any(|o| o.role == "thread.counterbore/shoulder")
                    })
                    .unwrap();
                let intent = GeometryIntent::capture(&snapshot, shoulder.key).unwrap();
                let edited = cut(&stock, feature, entry, left, 1.1).unwrap();
                assert!(matches!(
                    resolve(&intent, &edited.topology_snapshot(stock_id, 2).unwrap(), 2),
                    Resolution::Resolved { .. }
                ));
                let receipt = seated
                    .manufacturing_mesh_with_receipt(Default::default())
                    .unwrap();
                assert!(receipt.edge_pairs > 0 && receipt.vertex_links > 0);
                assert!(
                    (multicad_kernel::rust_kernel::manufacturing::certify(
                        &multicad_kernel::rust_kernel::display::DisplayMesh {
                            vertices: receipt.mesh.vertices.clone(),
                            indices: receipt.mesh.indices.clone(),
                            face_indices: receipt.mesh.face_indices.clone()
                        }
                    )
                    .unwrap()
                    .signed_volume_mm3
                        - seated.volume_mm3().unwrap())
                    .abs()
                        < 1.
                );
                assert_eq!(stock.topology_snapshot(stock_id, 0).unwrap(), original);
                assert!(cut(&stock, feature, entry, left, 3.5).is_err());
                assert_eq!(stock.topology_snapshot(stock_id, 0).unwrap(), original);
                // No legacy blind STEP support is inferred from a new valid mesh.
                assert!(!seated.supports_step_export());
            }
        }
    }
}

#[test]
fn actual_regular_screw_getter_retains_source_and_reports_thread_start_after_pose() {
    use multicad_kernel::{
        RigidTransform,
        rust_kernel::headed_screw::{ScrewHead, ScrewSpec},
    };
    let policy = GeometryPolicy::default();
    let source_frame = PlaneFrame {
        origin: [5., 6., 7.],
        u: [0.8, 0.6, 0.],
        v: [-0.36, 0.48, 0.8],
    };
    let occurrence = PlaneFrame {
        origin: [100., -30., 50.],
        u: [0., 1., 0.],
        v: [0., 0., 1.],
    };
    let transform = RigidTransform {
        rotation: std::array::from_fn(|i| {
            [occurrence.u[i], occurrence.v[i], occurrence.normal()[i]]
        }),
        translation_mm: occurrence.origin,
    };
    let head = ScrewSpec {
        head: ScrewHead::SocketCap {
            radius_mm: 2.75,
            height_mm: 3.,
            socket_across_flats_mm: 2.5,
            socket_depth_mm: 1.3,
        },
        shank_length_mm: 0.375,
        shank_radius_mm: 1.5,
    };
    for left in [false, true] {
        let mut thread = spec(left);
        thread.starts = 2;
        let screw =
            Solid::rust_machine_screw(thread, None, head, EntityId::new(), source_frame, policy)
                .unwrap();
        let (original, frame) = screw.rust_machine_screw_geometry().unwrap();
        assert_eq!(original.spec(), head);
        assert_eq!(original.thread().spec().starts, 2);
        assert_eq!(original.thread().spec().left_hand, left);
        let original = *original;
        let close = |a: [f64; 3], b: [f64; 3]| {
            for i in 0..3 {
                assert!((a[i] - b[i]).abs() < 1e-12, "{a:?} vs {b:?}");
            }
        };
        close(frame.origin, point(source_frame, [0., 0., 0.375]));
        close(
            point(frame, [0., 0., original.head_seat_z()]),
            source_frame.origin,
        );
        let posed = screw.transformed(transform).unwrap();
        let (actual, placed) = posed.rust_machine_screw_geometry().unwrap();
        assert_eq!(*actual, original);
        close(placed.origin, point(occurrence, frame.origin));
        close(
            point(placed, [0., 0., actual.head_seat_z()]),
            point(occurrence, source_frame.origin),
        );
        assert_eq!(*screw.rust_machine_screw_geometry().unwrap().0, original);
        assert_eq!(screw.rust_machine_screw_geometry().unwrap().1, frame);
        let bare = Solid::rust_helical_thread(thread, None, EntityId::new(), source_frame, policy)
            .unwrap();
        assert!(bare.rust_machine_screw_geometry().is_none());
    }
    let stock = Solid::rust_box(
        EntityId::new().uuid_bits(),
        [n(20); 3],
        PlaneFrame::default(),
        policy,
    )
    .unwrap();
    assert!(stock.rust_machine_screw_geometry().is_none());
}

#[test]
fn counterbore_actual_material_and_mating_proof_keep_the_original_blind_thread() {
    use multicad_kernel::rust_kernel::headed_screw::{ScrewHead, ScrewSpec};
    let policy = GeometryPolicy::default();
    let stock_id = EntityId::new();
    let stock = Solid::rust_box(
        stock_id.uuid_bits(),
        [n(20); 3],
        PlaneFrame::default(),
        policy,
    )
    .unwrap();
    let entry = PlaneFrame {
        origin: [10., 10., 20.],
        u: [1., 0., 0.],
        v: [0., -1., 0.],
    };
    for left in [false, true] {
        let mut profile = spec(left);
        profile.starts = 2;
        let feature = EntityId::new();
        let seated = stock
            .rust_cut_counterbored_blind_internal_thread(
                BlindThreadCut {
                    profile,
                    entry,
                    receiving_radius_mm: 1.5,
                    minimum_wall_mm: 1.,
                    minimum_floor_mm: 1.,
                },
                BlindThreadSeat {
                    radius_mm: 3.5,
                    depth_mm: 0.8,
                },
                feature,
                policy,
            )
            .unwrap();
        let screw = Solid::rust_machine_screw(
            profile,
            None,
            ScrewSpec {
                head: ScrewHead::ExternalHex {
                    across_flats_mm: 5.5,
                    height_mm: 2.,
                },
                shank_length_mm: 0.,
                shank_radius_mm: 1.5,
            },
            EntityId::new(),
            entry,
            policy,
        )
        .unwrap();
        let proof = screw.rust_screw_blind_thread(&seated).unwrap();
        assert_eq!(
            proof.counterbore.unwrap(),
            BlindThreadSeat {
                radius_mm: 3.5,
                depth_mm: 0.8
            }
        );
        assert_eq!(proof.threads[1].0.spec().length_mm, 3.5);
        assert_eq!(proof.threads[1].0.spec().starts, 2);
        assert_eq!(proof.threads[1].1, entry);
        assert_eq!(proof.protected_radius_mm, 4.5);
        assert!(proof.entry_is_support_plane);
        assert!(
            !seated
                .rust_thread_contains(point(entry, [2., 0., 0.4]))
                .unwrap()
        );
        assert!(
            seated
                .rust_thread_contains(point(entry, [2., 0., 1.]))
                .unwrap()
        );
        assert!(
            !seated
                .rust_thread_contains(point(entry, [0., 0., 3.4]))
                .unwrap()
        );
        assert!(
            seated
                .rust_thread_contains(point(entry, [0., 0., 3.6]))
                .unwrap()
        );
        let faces = seated.faces().unwrap();
        let cylinder = &faces[10];
        assert!((cylinder.area_mm2 - std::f64::consts::TAU * 3.5 * 0.8).abs() < 1e-10);
        assert!((cylinder.center_mm[2] - 19.6).abs() < 1e-12);
        assert!(cylinder.plane.is_none());
        let floor = &faces[11];
        assert!((floor.center_mm[2] - 19.2).abs() < 1e-12);
        assert!(floor.plane.is_some());
        assert_eq!(faces.len(), 17);
        let unseated = stock
            .rust_cut_blind_internal_thread(
                BlindThreadCut {
                    profile,
                    entry,
                    receiving_radius_mm: 1.5,
                    minimum_wall_mm: 1.,
                    minimum_floor_mm: 1.,
                },
                EntityId::new(),
                policy,
            )
            .unwrap();
        assert!(
            screw
                .rust_screw_blind_thread(&unseated)
                .unwrap()
                .counterbore
                .is_none()
        );
        assert!(seated.rust_screw_blind_thread(&seated).is_none());
        assert!(screw.rust_screw_blind_thread(&stock).is_none());
        assert!(screw.rust_screw_blind_thread(&screw).is_none());
        assert!(stock.rust_screw_blind_thread(&seated).is_none());
    }
}

#[test]
fn invalid_seated_blind_geometry_cannot_replace_current_stock() {
    let policy = GeometryPolicy::default();
    let owner = EntityId::new();
    let stock =
        Solid::rust_box(owner.uuid_bits(), [n(20); 3], PlaneFrame::default(), policy).unwrap();
    let original = stock.topology_snapshot(owner, 0).unwrap();
    let entry = PlaneFrame {
        origin: [10., 10., 20.],
        u: [1., 0., 0.],
        v: [0., -1., 0.],
    };
    for seat in [
        BlindThreadSeat {
            radius_mm: f64::NAN,
            depth_mm: 0.8,
        },
        BlindThreadSeat {
            radius_mm: 3.5,
            depth_mm: f64::INFINITY,
        },
        BlindThreadSeat {
            radius_mm: 1.5,
            depth_mm: 0.8,
        },
        BlindThreadSeat {
            radius_mm: 3.5,
            depth_mm: 3.5,
        },
        BlindThreadSeat {
            radius_mm: 10.,
            depth_mm: 0.8,
        },
        BlindThreadSeat {
            radius_mm: 1.5 + 1e-8,
            depth_mm: 0.8,
        },
    ] {
        assert!(
            stock
                .rust_cut_counterbored_blind_internal_thread(
                    BlindThreadCut {
                        profile: spec(false),
                        entry,
                        receiving_radius_mm: 1.5,
                        minimum_wall_mm: 1.,
                        minimum_floor_mm: 1.
                    },
                    seat,
                    EntityId::new(),
                    policy
                )
                .is_err()
        );
        assert_eq!(stock.topology_snapshot(owner, 0).unwrap(), original);
    }
}

#[test]
fn seated_blind_rejects_wrong_policy_frame_profile_and_source_atomically() {
    let policy = GeometryPolicy::default();
    let owner = EntityId::new();
    let stock =
        Solid::rust_box(owner.uuid_bits(), [n(20); 3], PlaneFrame::default(), policy).unwrap();
    let original = stock.topology_snapshot(owner, 0).unwrap();
    let entry = PlaneFrame {
        origin: [10., 10., 20.],
        u: [1., 0., 0.],
        v: [0., -1., 0.],
    };
    for case in 0..7 {
        let mut cut = BlindThreadCut {
            profile: spec(false),
            entry,
            receiving_radius_mm: 1.5,
            minimum_wall_mm: 1.,
            minimum_floor_mm: 1.,
        };
        let mut requested = policy;
        let mut feature = EntityId::new();
        match case {
            0 => requested.minimum_feature_mm = -1.,
            1 => cut.entry.u = [0.; 3],
            2 => cut.receiving_radius_mm = -1.,
            3 => requested.minimum_feature_mm *= 2.,
            4 => cut.profile.pitch_mm = 0.,
            5 => feature = owner,
            _ => cut.entry.origin[2] += 0.1,
        }
        assert!(
            stock
                .rust_cut_counterbored_blind_internal_thread(
                    cut,
                    BlindThreadSeat {
                        radius_mm: 3.5,
                        depth_mm: 0.8
                    },
                    feature,
                    requested
                )
                .is_err(),
            "case{case}"
        );
        assert_eq!(stock.topology_snapshot(owner, 0).unwrap(), original);
    }
    let bare =
        Solid::rust_helical_thread(spec(false), None, EntityId::new(), entry, policy).unwrap();
    assert!(cut(&bare, EntityId::new(), entry, false, 0.8).is_err());
}

#[test]
fn intrinsic_counterbore_invalidity_and_stricter_policy_reserves_both_refuse_atomically() {
    let policy = GeometryPolicy::default();
    let owner = EntityId::new();
    let stock =
        Solid::rust_box(owner.uuid_bits(), [n(20); 3], PlaneFrame::default(), policy).unwrap();
    let before = stock.topology_snapshot(owner, 0).unwrap();
    let entry = PlaneFrame {
        origin: [10., 10., 20.],
        u: [1., 0., 0.],
        v: [0., -1., 0.],
    };
    for (radius, depth, message) in [
        (f64::NAN, 0.8, "Invalid finite counterbore"),
        (1.5, 0.8, "Invalid finite counterbore"),
        (3.5, 3.5, "Invalid finite counterbore"),
        (1e8, 0.8, "Invalid finite counterbore"),
        (
            1.5 + policy.minimum_feature_mm / 2.,
            0.8,
            "Length is below minimum feature size",
        ),
        (
            3.5,
            policy.minimum_feature_mm / 2.,
            "Length is below minimum feature size",
        ),
    ] {
        let cause = stock
            .rust_cut_counterbored_blind_internal_thread(
                BlindThreadCut {
                    profile: spec(false),
                    entry,
                    receiving_radius_mm: 1.5,
                    minimum_wall_mm: 1.,
                    minimum_floor_mm: 1.,
                },
                BlindThreadSeat {
                    radius_mm: radius,
                    depth_mm: depth,
                },
                EntityId::new(),
                policy,
            )
            .err()
            .unwrap();
        assert_eq!(cause.kind, multicad_kernel::GeometryErrorKind::Input);
        assert!(cause.message.contains(message), "{cause}");
        assert_eq!(stock.topology_snapshot(owner, 0).unwrap(), before);
    }
    let accepted = cut(&stock, EntityId::new(), entry, false, 0.8).unwrap();
    assert!(
        !accepted
            .rust_thread_contains(point(entry, [2., 0., 0.4]))
            .unwrap()
    );
    assert!(
        accepted
            .rust_thread_contains(point(entry, [2., 0., 1.]))
            .unwrap()
    );
    assert_eq!(stock.topology_snapshot(owner, 0).unwrap(), before);
}
