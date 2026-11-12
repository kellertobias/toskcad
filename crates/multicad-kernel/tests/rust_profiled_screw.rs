use multicad_core::{
    document::EntityId,
    references::{AnalyticGeometry, TopologyKind},
};
use multicad_kernel::{
    GeometryPolicy, HelicalThread, MeshOptions, PlaneFrame, RigidTransform, Solid,
    rust_kernel::{
        head_profile::{Head, Role},
        profiled_screw::Spec,
    },
};

fn thread() -> HelicalThread {
    HelicalThread {
        major_radius_mm: 1.5,
        pitch_mm: 0.5,
        length_mm: 7.5,
        depth_mm: 0.2,
        crest_width_mm: 0.0625,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        radial_taper: 0.,
        starts: 1,
        left_hand: false,
    }
}
fn spec(crown: bool) -> Spec {
    Spec {
        head: if crown {
            Head::Crowned {
                radius_mm: 3.,
                height_mm: 2.5,
                crown_radius_mm: 4.,
                top_land_radius_mm: 1.5,
                neck_radius_mm: 0.1,
            }
        } else {
            Head::Countersunk {
                radius_mm: 3.,
                height_mm: 1.65,
                included_angle_deg: 90.,
                edge_chamfer_mm: 0.1,
                neck_radius_mm: 0.1,
            }
        },
        shank_length_mm: 1.,
        socket_across_flats_mm: 2.,
        socket_depth_mm: 0.6,
        minimum_wall_mm: 0.1,
    }
}

#[test]
fn actual_heads_export_and_import_through_default_rust_facade() {
    let scratch = tempfile::tempdir().unwrap();
    for crown in [false, true] {
        for (starts, left_hand) in [(1, false), (2, true)] {
            let policy = GeometryPolicy::default();
            let owner = EntityId::new();
            let solid = Solid::rust_profiled_machine_screw(
                HelicalThread {
                    starts,
                    left_hand,
                    ..thread()
                },
                None,
                spec(crown),
                owner,
                PlaneFrame {
                    origin: [30., -20., 10.],
                    u: [0.6, 0.8, 0.],
                    v: [0., 0., 1.],
                },
                policy,
            )
            .unwrap();
            assert!(solid.supports_step_export());
            let path = scratch.path().join(format!("head-{crown}-{starts}.step"));
            solid.write_step(&path).unwrap();
            let text = std::fs::read_to_string(&path).unwrap();
            let (read, imported) = multicad_kernel::rust_kernel::step_read::read_import_solid(
                &text,
                policy.linear_tolerance_mm,
            )
            .unwrap();
            let (source, frame) = solid.rust_profiled_screw_geometry().unwrap();
            assert!(
                multicad_kernel::rust_kernel::step_read::profiled_screw_import_matches(
                    &imported,
                    source,
                    [frame.origin, frame.u, frame.v],
                    policy.linear_tolerance_mm
                )
            );
            let reopened = Solid::import_rust_step(&text, EntityId::new(), policy).unwrap();
            let (retained, _) = reopened.rust_profiled_screw_geometry().unwrap();
            assert!(retained.spec().head.is_none());
            assert!((solid.volume_mm3().unwrap() - reopened.volume_mm3().unwrap()).abs() < 1e-6);
            assert!((read.volume_mm3 - reopened.volume_mm3().unwrap()).abs() < 1e-6);
            let mesh = reopened
                .manufacturing_mesh_with_receipt(MeshOptions::default())
                .unwrap();
            assert_eq!(mesh.edge_pairs as usize * 2, mesh.mesh.indices.len());
            reopened
                .topology_snapshot(EntityId::new(), 0)
                .unwrap()
                .validate()
                .unwrap();
        }
    }
}

#[test]
fn fitted_shank_annular_step_exports_and_imports_without_losing_nominal_radius() {
    let scratch = tempfile::tempdir().unwrap();
    for crown in [false, true] {
        for (starts, left_hand) in [(1, false), (2, true)] {
            let policy = GeometryPolicy::default();
            let solid = Solid::rust_profiled_machine_screw_with_shank(
                HelicalThread {
                    major_radius_mm: 1.48,
                    starts,
                    left_hand,
                    ..thread()
                },
                None,
                1.5,
                spec(crown),
                EntityId::new(),
                PlaneFrame {
                    origin: [30., -20., 10.],
                    u: [0.6, 0.8, 0.],
                    v: [0., 0., 1.],
                },
                policy,
            )
            .unwrap();
            assert!(solid.supports_step_export());
            let path = scratch.path().join(format!("fitted-{crown}-{starts}.step"));
            solid.write_step(&path).unwrap();
            let text = std::fs::read_to_string(path).unwrap();
            let imported = Solid::import_rust_step(&text, EntityId::new(), policy).unwrap();
            let (source, _) = imported.rust_profiled_screw_geometry().unwrap();
            assert_eq!(source.profile().end()[0], 1.5);
            assert_eq!(source.thread().spec().major_radius_mm, 1.48);
            assert!(source.spec().head.is_none());
            assert!((solid.volume_mm3().unwrap() - imported.volume_mm3().unwrap()).abs() < 1e-6);
            let mesh = imported
                .manufacturing_mesh_with_receipt(MeshOptions::default())
                .unwrap();
            assert_eq!(mesh.edge_pairs as usize * 2, mesh.mesh.indices.len());
            imported
                .topology_snapshot(EntityId::new(), 0)
                .unwrap()
                .validate()
                .unwrap();
        }
    }
}

#[test]
fn retained_socket_contour_controls_actual_faces_planes_references_and_exchange() {
    use multicad_kernel::rust_kernel::profiled_screw::{Face, Screw};
    let scratch = tempfile::tempdir().unwrap();
    for crown in [false, true] {
        let policy = GeometryPolicy::default();
        let owner = EntityId::new();
        let original = Solid::rust_profiled_machine_screw(
            thread(),
            None,
            spec(crown),
            owner,
            PlaneFrame {
                origin: [30., -20., 10.],
                u: [0.6, 0.8, 0.],
                v: [0., 0., 1.],
            },
            policy,
        )
        .unwrap();
        let (source, at) = original.rust_profiled_screw_geometry().unwrap();
        let (s, c) = 0.37_f64.sin_cos();
        let points = source
            .drive()
            .corners_mm()
            .map(|[x, y]| [c * x - s * y + 1e-14, s * x + c * y]);
        let mut physical = source.spec();
        physical.head = None;
        let source =
            Screw::from_source_with_drive(*source.thread(), *source.profile(), physical, points)
                .unwrap();
        let retained = Solid::rust_profiled_screw_from_source(source, owner, at, policy).unwrap();
        let snapshot = retained.topology_snapshot(owner, 0).unwrap();
        snapshot.validate().unwrap();
        let wall = snapshot
            .elements
            .iter()
            .find(|e| {
                e.key.kind == TopologyKind::Face
                    && e.origins.iter().any(|o| o.role == "screw.drive/edge/0")
            })
            .unwrap();
        let frame = wall.frames.first().unwrap().frame;
        let project = |[x, y, z]: [f64; 3]| {
            std::array::from_fn::<_, 3, _>(|i| {
                at.origin[i] + at.u[i] * x + at.v[i] * y + at.normal()[i] * z
            })
        };
        let expected = project([points[0][0], points[0][1], source.head_top_z()]);
        let d = [points[1][0] - points[0][0], points[1][1] - points[0][1]];
        let length = d[0].hypot(d[1]);
        assert!(
            frame
                .u
                .into_iter()
                .enumerate()
                .all(|(i, x)| (x + (at.u[i] * d[0] + at.v[i] * d[1]) / length).abs() < 1e-12)
        );
        assert!(
            frame
                .origin_mm
                .into_iter()
                .zip(expected)
                .all(|(a, b)| (a - b).abs() < 1e-12)
        );
        let index = source
            .faces()
            .iter()
            .position(|f| *f == Face::SocketWall(0))
            .unwrap();
        let faces = retained.faces().unwrap();
        assert!(
            (faces[index].area_mm2
                - source.drive().side_length_mm(0) * source.spec().socket_depth_mm)
                .abs()
                < 1e-12
        );
        for (i, point) in points.into_iter().enumerate() {
            let vertex = snapshot
                .elements
                .iter()
                .find(|e| {
                    e.key.kind == TopologyKind::Vertex
                        && e.origins
                            .iter()
                            .any(|o| o.role == format!("screw.drive/corner/0/{i}"))
                })
                .unwrap();
            assert!(
                matches!(vertex.geometry,Some(AnalyticGeometry::Point { position_mm }) if position_mm.into_iter().zip(project([point[0],point[1],source.head_top_z()])).all(|(a,b)|(a-b).abs()<1e-12))
            );
        }
        let path = scratch.path().join(format!("rotated-drive-{crown}.step"));
        retained.write_step(&path).unwrap();
        let text = std::fs::read_to_string(path).unwrap();
        let imported = Solid::import_rust_step(&text, EntityId::new(), policy).unwrap();
        let (decoded, _) = imported.rust_profiled_screw_geometry().unwrap();
        assert!(!decoded.drive().is_canonical());
        let mesh = imported
            .manufacturing_mesh_with_receipt(MeshOptions::default())
            .unwrap();
        assert_eq!(mesh.edge_pairs as usize * 2, mesh.mesh.indices.len());
    }
}

#[test]
fn typed_source_and_finite_faces_survive_rotated_placement() {
    for crown in [false, true] {
        let owner = EntityId::new();
        let frame = PlaneFrame {
            origin: [7., -3., 4.],
            u: [0., 0., 1.],
            v: [0., -1., 0.],
        };
        let solid = Solid::rust_profiled_machine_screw(
            thread(),
            None,
            spec(crown),
            owner,
            frame,
            GeometryPolicy::default(),
        )
        .unwrap();
        let (source, at_thread) = solid.rust_profiled_screw_geometry().unwrap();
        assert_eq!(source.spec(), spec(crown).into());
        assert_eq!(at_thread.origin, [8., -3., 4.]);
        assert_eq!(at_thread.u, frame.u);
        assert_eq!(at_thread.v, frame.v);
        assert_eq!(source.thread().spec().length_mm, 7.5);
        assert!(
            source
                .profile()
                .segments()
                .any(|s| s.role == if crown { Role::Crown } else { Role::Cone })
        );
        assert_eq!(solid.volume_mm3().unwrap(), source.volume_mm3());
        assert!(
            !solid
                .rust_thread_contains([7. + source.profile().top()[1] + 0.1, -3., 4.])
                .unwrap()
        );
        let faces = solid.faces().unwrap();
        assert_eq!(faces.len(), source.faces().len());
        assert!(
            faces
                .iter()
                .all(|f| f.area_mm2 > 0. && f.center_mm.into_iter().all(f64::is_finite))
        );
        let topology = solid.topology_snapshot(owner, 3).unwrap();
        topology.validate().unwrap();
        let names = |name: &str| {
            topology
                .elements
                .iter()
                .find(|e| {
                    e.key.kind == TopologyKind::Face && e.origins.iter().any(|o| o.role == name)
                })
                .unwrap()
        };
        assert!(names("screw.head/neck.blend").boundary.len() >= 2);
        assert!(
            names(if crown {
                "screw.head/crown"
            } else {
                "screw.head/cone"
            })
            .geometry
            .is_none()
        );
        assert!(matches!(
            names("screw.drive/end").geometry,
            Some(AnalyticGeometry::Plane { .. })
        ));
        assert_eq!(names("screw.drive/end").boundary.len(), 6);
        let shoulder = names("screw.thread/start");
        assert!(
            shoulder
                .boundary
                .iter()
                .all(|o| !o.role.ends_with("/crest"))
        );
        let mesh = solid
            .manufacturing_mesh_with_receipt(MeshOptions::default())
            .unwrap();
        assert_eq!(mesh.edge_pairs as usize * 2, mesh.mesh.indices.len());
        assert!(solid.display_edges(0.03).unwrap().len() > 30);
        let shifted = solid
            .transformed(RigidTransform {
                translation_mm: [2., 3., 4.],
                ..Default::default()
            })
            .unwrap();
        assert_eq!(
            shifted.rust_profiled_screw_geometry().unwrap().1.origin,
            [10., 0., 8.]
        );
        assert_eq!(shifted.volume_mm3().unwrap(), solid.volume_mm3().unwrap());
        assert!(solid.supports_step_export());
    }
}

#[test]
fn precision_and_invalid_drive_edits_preserve_the_previous_source() {
    let owner = EntityId::new();
    let make = |spec| {
        Solid::rust_profiled_machine_screw(
            thread(),
            None,
            spec,
            owner,
            PlaneFrame::default(),
            GeometryPolicy::default(),
        )
    };
    let valid = make(spec(false)).unwrap();
    let before = valid.volume_mm3().unwrap();
    assert!(
        make(Spec {
            socket_depth_mm: 100.,
            ..spec(false)
        })
        .is_err()
    );
    let thin = Spec {
        head: Head::Countersunk {
            radius_mm: 3.,
            height_mm: 1.65,
            included_angle_deg: 90.,
            edge_chamfer_mm: 0.1,
            neck_radius_mm: 1e-7,
        },
        ..spec(false)
    };
    assert!(make(thin).is_err());
    let unrepresentable = Spec {
        head: Head::Crowned {
            radius_mm: 3.,
            height_mm: 2.5,
            crown_radius_mm: 1e6,
            top_land_radius_mm: 1.5,
            neck_radius_mm: 0.1,
        },
        ..spec(true)
    };
    assert!(make(unrepresentable).is_err());
    assert_eq!(valid.volume_mm3().unwrap(), before);
    assert!(valid.manufacturing_mesh(MeshOptions::default()).is_ok());
    let shaft = Solid::rust_helical_thread(
        thread(),
        None,
        owner,
        PlaneFrame::default(),
        GeometryPolicy::default(),
    )
    .unwrap();
    assert!(shaft.rust_profiled_screw_geometry().is_none());
    let box_solid = Solid::box_solid([1.; 3], GeometryPolicy::default()).unwrap();
    assert!(box_solid.rust_profiled_screw_geometry().is_none());
}

#[test]
fn covered_thread_prefix_keeps_actual_visible_faces_rims_and_phase() {
    let scratch = tempfile::tempdir().unwrap();
    for crown in [false, true] {
        for left_hand in [false, true] {
            let mut thread = thread();
            thread.left_hand = left_hand;
            thread.starts = 2;
            let owner = EntityId::new();
            let frame = PlaneFrame::default();
            let full = Solid::rust_profiled_machine_screw(
                thread,
                None,
                spec(crown),
                owner,
                frame,
                GeometryPolicy::default(),
            )
            .unwrap();
            for plain in [0., 0.03125] {
                let mut spec = spec(crown);
                spec.shank_length_mm = plain;
                let solid = Solid::rust_profiled_machine_screw_with_head_thread_union(
                    thread,
                    None,
                    spec,
                    owner,
                    frame,
                    GeometryPolicy::default(),
                )
                .unwrap();
                let (source, at) = solid.rust_profiled_screw_geometry().unwrap();
                let [start, end] = source.visible_thread_interval_mm();
                assert!(start > 0. && start < end);
                assert_eq!(source.thread().spec().length_mm, thread.length_mm);
                assert_eq!(at.origin[2], plain);
                let faces = solid.faces().unwrap();
                let full_faces = full.faces().unwrap();
                assert!((faces[0].center_mm[2] - plain - start).abs() < 1e-12);
                assert!((faces[1].center_mm[2] - plain - end).abs() < 1e-12);
                for i in 2..2 + 4 * thread.starts as usize {
                    let ratio = (end - start) / end;
                    assert!((faces[i].area_mm2 / full_faces[i].area_mm2 - ratio).abs() < 1e-12);
                    assert!((faces[i].center_mm[2] - plain - start.midpoint(end)).abs() < 1e-12);
                }
                let topology = solid.topology_snapshot(owner, 1).unwrap();
                topology.validate().unwrap();
                let shoulder = topology
                    .elements
                    .iter()
                    .find(|e| {
                        e.key.kind == TopologyKind::Face
                            && e.origins.iter().any(|o| o.role == "screw.thread/start")
                    })
                    .unwrap();
                let Some(AnalyticGeometry::Plane { frame: cap }) = shoulder.geometry else {
                    panic!("actual shoulder plane")
                };
                assert!((cap.origin_mm[2] - plain - start).abs() < 1e-12);
                let count = 4 * thread.starts as usize;
                let edges = solid.display_edges(0.01).unwrap();
                for edge in &edges[..count] {
                    assert!((edge.points_mm.first().unwrap()[2] - plain - start).abs() < 1e-12);
                    assert!((edge.points_mm.last().unwrap()[2] - plain - end).abs() < 1e-12);
                }
                for element in topology.elements.iter().filter(|e| {
                    e.key.kind == TopologyKind::Vertex && e.key.ordinal < 2 * count as u32
                }) {
                    let Some(AnalyticGeometry::Point { position_mm: p }) = element.geometry else {
                        panic!("actual helical rim vertex")
                    };
                    assert!(p[2] >= plain + start - 1e-12);
                    let z = p[2] - plain;
                    let radius = p[0].hypot(p[1]);
                    assert!((radius - source.thread().radius(p[1].atan2(p[0]), z)).abs() < 1e-10);
                }
                assert_eq!(solid.volume_mm3().unwrap(), source.volume_mm3());
                if plain == 0. && !left_hand {
                    let mesh = solid
                        .manufacturing_mesh_with_receipt(MeshOptions::default())
                        .unwrap();
                    assert_eq!(mesh.edge_pairs as usize * 2, mesh.mesh.indices.len());
                }
                assert!(solid.supports_step_export());
                let path = scratch
                    .path()
                    .join(format!("union-{crown}-{left_hand}-{plain}.step"));
                solid.write_step(&path).unwrap();
                let text = std::fs::read_to_string(path).unwrap();
                let imported =
                    Solid::import_rust_step(&text, EntityId::new(), GeometryPolicy::default())
                        .unwrap();
                let (decoded, _) = imported.rust_profiled_screw_geometry().unwrap();
                assert!(decoded.has_head_thread_overlap());
                assert!((decoded.thread().spec().length_mm - thread.length_mm).abs() < 1e-12);
                assert!(
                    (imported.volume_mm3().unwrap() - solid.volume_mm3().unwrap()).abs() < 1e-6
                );
            }
        }
    }
}

#[test]
fn hidden_prefix_cannot_leave_a_subminimum_visible_thread() {
    let owner = EntityId::new();
    let policy = GeometryPolicy::default();
    let mut t = thread();
    t.pitch_mm = 0.05;
    t.depth_mm = 0.02;
    t.crest_width_mm = 0.00625;
    t.length_mm = 0.5;
    let mut spec = spec(true);
    spec.shank_length_mm = 0.;
    if let Head::Crowned {
        ref mut neck_radius_mm,
        ..
    } = spec.head
    {
        *neck_radius_mm = 0.125;
    }
    let make = |t| {
        Solid::rust_profiled_machine_screw_with_head_thread_union(
            t,
            None,
            spec,
            owner,
            PlaneFrame::default(),
            policy,
        )
    };
    let previous = make(t).unwrap();
    let volume = previous.volume_mm3().unwrap();
    let start = previous
        .rust_profiled_screw_geometry()
        .unwrap()
        .0
        .visible_thread_interval_mm()[0];
    t.length_mm = start + policy.minimum_feature_mm / 2.;
    assert!(t.length_mm > policy.minimum_feature_mm);
    assert!(make(t).is_err());
    t.length_mm = start + 2. * policy.minimum_feature_mm;
    assert!(make(t).is_ok());
    assert_eq!(previous.volume_mm3().unwrap(), volume);
    assert_eq!(
        previous
            .rust_profiled_screw_geometry()
            .unwrap()
            .0
            .thread()
            .spec()
            .length_mm,
        0.5
    );
}

#[test]
fn physical_head_source_keeps_literal_supports_transport_budget_and_original_frame() {
    use multicad_kernel::rust_kernel::{
        head_profile::{Profile, TransportError},
        profiled_screw::Screw,
    };
    let policy = GeometryPolicy::default();
    for crown in [false, true] {
        for union in [false, true] {
            let owner = EntityId::new();
            let frame = PlaneFrame {
                origin: [1000., -2000., 4000.],
                u: [0., 0., 1.],
                v: [0., -1., 0.],
            };
            let mut spec = spec(crown);
            if union {
                spec.shank_length_mm = 0.;
            }
            let authored = if union {
                Solid::rust_profiled_machine_screw_with_head_thread_union(
                    thread(),
                    None,
                    spec,
                    owner,
                    frame,
                    policy,
                )
            } else {
                Solid::rust_profiled_machine_screw(thread(), None, spec, owner, frame, policy)
            }
            .unwrap();
            let (original, at) = authored.rust_profiled_screw_geometry().unwrap();
            let error = TransportError {
                numeric_conversion_mm: 1e-11,
                arithmetic_mm: 2e-11,
            };
            let segments: Vec<_> = original.profile().segments().collect();
            let profile = Profile::from_transported_segments(&segments, error).unwrap();
            let mut physical_spec = original.spec();
            physical_spec.head = None;
            let physical = if union {
                Screw::from_source_with_union(*original.thread(), profile, physical_spec)
            } else {
                Screw::from_source(*original.thread(), profile, physical_spec)
            }
            .unwrap();
            let make = |policy| Solid::rust_profiled_screw_from_source(physical, owner, at, policy);
            let solid = make(policy).unwrap();
            let (retained, pose) = solid.rust_profiled_screw_geometry().unwrap();
            assert_eq!(*retained, physical);
            assert_eq!(pose, at);
            assert_eq!(retained.spec().head, None);
            assert_eq!(retained.profile().transport_error(), error);
            assert_eq!(retained.profile().segments().collect::<Vec<_>>(), segments);
            assert_eq!(solid.volume_mm3().unwrap(), authored.volume_mm3().unwrap());
            let snapshot = solid.topology_snapshot(owner, 1).unwrap();
            assert_eq!(snapshot, authored.topology_snapshot(owner, 1).unwrap());
            let mesh = solid
                .manufacturing_mesh_with_receipt(MeshOptions::default())
                .unwrap();
            assert_eq!(mesh.edge_pairs as usize * 2, mesh.mesh.indices.len());
            assert!(
                make(GeometryPolicy {
                    linear_tolerance_mm: 1e-12,
                    ..policy
                })
                .is_err()
            );
            assert!(
                Solid::rust_profiled_screw_from_source(physical, EntityId::nil(), at, policy)
                    .is_err()
            );
            assert_eq!(
                solid
                    .rust_profiled_screw_geometry()
                    .unwrap()
                    .0
                    .profile()
                    .transport_error(),
                error
            );
        }
    }
}

#[test]
fn nominal_shank_and_fitted_rounded_thread_keep_both_actual_shoulder_boundaries() {
    let owner = EntityId::new();
    let policy = GeometryPolicy::default();
    for crown in [false, true] {
        for plain in [0., 1.] {
            let mut t = thread();
            t.major_radius_mm = 1.48;
            let root = 0.05;
            t.depth_mm = (t.pitch_mm - t.crest_width_mm) * 3_f64.sqrt() / 2. - root;
            let mut s = spec(crown);
            s.shank_length_mm = plain;
            let solid = Solid::rust_profiled_machine_screw_with_shank(
                t,
                Some(root),
                1.5,
                s,
                owner,
                PlaneFrame::default(),
                policy,
            )
            .unwrap();
            let (source, _) = solid.rust_profiled_screw_geometry().unwrap();
            assert_eq!(source.profile().end()[0], 1.5);
            assert_eq!(source.thread().spec().major_radius_mm, 1.48);
            assert_eq!(source.thread().spec().root_radius_mm, Some(root));
            let faces = solid.faces().unwrap();
            let removed = source.thread().volume_mm3() / t.length_mm;
            assert!(
                (faces[0].area_mm2 - (std::f64::consts::PI * 1.5_f64.powi(2) - removed)).abs()
                    < 1e-10
            );
            let snapshot = solid.topology_snapshot(owner, 1).unwrap();
            let shoulder = snapshot
                .elements
                .iter()
                .find(|e| {
                    e.key.kind == TopologyKind::Face
                        && e.origins.iter().any(|o| o.role == "screw.thread/start")
                })
                .unwrap();
            let mut outer = 0;
            let mut crest = 0;
            for edge in snapshot.elements.iter().filter(|e| {
                e.key.kind == TopologyKind::Edge && !e.origins.is_disjoint(&shoulder.boundary)
            }) {
                if let Some(AnalyticGeometry::CircularArc { radius_mm, .. }) = edge.geometry {
                    if radius_mm == 1.5 {
                        outer += 1;
                    }
                    if radius_mm == t.major_radius_mm {
                        crest += 1;
                    }
                }
            }
            assert_eq!(outer, 1);
            assert!(crest > 0);
            let mesh = solid
                .manufacturing_mesh_with_receipt(MeshOptions::default())
                .unwrap();
            assert_eq!(mesh.edge_pairs as usize * 2, mesh.mesh.indices.len());
            let original = solid.volume_mm3().unwrap();
            t.major_radius_mm = 1.5 - policy.minimum_feature_mm / 2.;
            assert!(
                Solid::rust_profiled_machine_screw_with_shank(
                    t,
                    Some(root),
                    1.5,
                    s,
                    owner,
                    PlaneFrame::default(),
                    policy
                )
                .is_err()
            );
            assert_eq!(solid.volume_mm3().unwrap(), original);
        }
    }
}
