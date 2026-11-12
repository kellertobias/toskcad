use multicad_core::document::EntityId;
use multicad_kernel::{GeometryPolicy, HelicalThread, MeshOptions, PlaneFrame, Solid};
fn request(left: bool, starts: u32) -> (HelicalThread, f64) {
    let root = 1. / (4. * 3_f64.sqrt());
    (
        HelicalThread {
            major_radius_mm: 3.,
            pitch_mm: 1.,
            length_mm: 6.6,
            depth_mm: 0.875 * 3_f64.sqrt() / 2. - root,
            crest_width_mm: 0.125,
            lower_flank_deg: 30.,
            upper_flank_deg: 30.,
            radial_taper: 0.,
            starts,
            left_hand: left,
        },
        root,
    )
}

#[test]
fn original_multistart_far_rims_and_corners_bind_once_to_complete_source_chains() {
    use multicad_core::references::{
        AnalyticGeometry, GeometryIntent, Resolution, TopologyKind, resolve,
    };
    for left in [false, true] {
        let (spec, root) = request(left, 2);
        let feature = EntityId::new();
        let frame = PlaneFrame::default();
        let original =
            Solid::rust_helical_thread(spec, Some(root), feature, frame, Default::default())
                .unwrap();
        let old = original.topology_snapshot(feature, 0).unwrap();
        let intents = old
            .elements
            .iter()
            .filter(|e| {
                e.origins.iter().any(|o| {
                    o.role.starts_with("cap.end/rim/")
                        || o.role.starts_with("cap.end/corner/")
                        || (e.key.kind == TopologyKind::Face && o.role.starts_with("thread/start/"))
                })
            })
            .map(|e| (GeometryIntent::capture(&old, e.key).unwrap(), e.geometry))
            .collect::<Vec<_>>();
        assert_eq!(intents.len(), 24);
        let body = Solid::rust_helical_thread_entry_chamfer(
            spec,
            Some(root),
            0.3,
            None,
            feature,
            frame,
            Default::default(),
        )
        .unwrap();
        let snapshot = body.topology_snapshot(feature, 1).unwrap();
        for (intent, geometry) in intents {
            let Resolution::Resolved {
                geometry: actual, ..
            } = resolve(&intent, &snapshot, 1)
            else {
                panic!("Lost original source ownership {:?}", intent.origins)
            };
            if let (
                Some(AnalyticGeometry::Point { position_mm: a }),
                Some(AnalyticGeometry::Point { position_mm: b }),
            ) = (geometry, actual)
            {
                assert!((0..3).map(|i| (a[i] - b[i]).abs()).sum::<f64>() < 1e-7);
            }
        }
        for edge in body.rust_thread_entry_chamfer_edges().unwrap().0 {
            if edge.source_edges().len() > 1 {
                let graph = body.rust_thread_entry_chamfer_topology().unwrap().0;
                for pair in edge.source_edges().windows(2) {
                    assert_eq!(
                        graph.edges()[pair[0]].vertices()[1],
                        graph.edges()[pair[1]].vertices()[0]
                    );
                }
            }
        }
    }
}
#[test]
fn retained_entry_material_publishes_actual_faces_mesh_and_source_known_file_proof() {
    let root = 1. / (4. * 3_f64.sqrt());
    let spec = HelicalThread {
        major_radius_mm: 3.,
        pitch_mm: 1.,
        length_mm: 6.6,
        depth_mm: 0.875 * 3_f64.sqrt() / 2. - root,
        crest_width_mm: 0.125,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        radial_taper: 0.,
        starts: 1,
        left_hand: false,
    };
    let body = Solid::rust_helical_thread_entry_chamfer(
        spec,
        Some(root),
        0.3,
        None,
        EntityId::new(),
        PlaneFrame::default(),
        GeometryPolicy::default(),
    )
    .unwrap();
    assert!(body.uses_rust_geometry());
    assert!(body.rust_thread_entry_chamfer_data().is_some());
    let (properties, _) = body.rust_thread_entry_chamfer_properties().unwrap();
    assert!(
        properties.volume_mm3.lower <= body.volume_mm3().unwrap()
            && body.volume_mm3().unwrap() <= properties.volume_mm3.upper
    );
    assert_eq!(body.faces().unwrap().len(), 8);
    let topology = body.topology_snapshot(EntityId::new(), 0).unwrap();
    topology.validate().unwrap();
    assert!(body.rust_thread_contains([0., 0., 0.]).unwrap());
    assert!(!body.rust_thread_contains([2.9, 0., 0.]).unwrap());
    assert_eq!(
        body.planar_face_contours(
            1,
            PlaneFrame {
                origin: [0., 0., 6.6],
                ..Default::default()
            }
        )
        .unwrap_err()
        .kind,
        multicad_kernel::GeometryErrorKind::Unsupported
    );
    assert_eq!(
        multicad_kernel::project_drawing(
            &[multicad_kernel::DrawingSolid {
                solid: &body,
                hatch: None
            }],
            Default::default()
        )
        .unwrap_err()
        .kind,
        multicad_kernel::GeometryErrorKind::Unsupported
    );
    let mesh = body.manufacturing_mesh(MeshOptions::default()).unwrap();
    assert!(!mesh.vertices.is_empty());
    let dir = tempfile::tempdir().unwrap();
    let stock = Solid::box_solid([1.; 3], Default::default()).unwrap();
    let compound =
        Solid::compound_named(&[&stock, &body], EntityId::new(), Default::default()).unwrap();
    let cap = stock.faces().unwrap().len() as u32 + 1;
    assert_eq!(
        compound
            .planar_face_contours_with_sources(
                cap,
                PlaneFrame {
                    origin: [0., 0., 6.6],
                    ..Default::default()
                }
            )
            .unwrap_err()
            .kind,
        multicad_kernel::GeometryErrorKind::Unsupported
    );
    assert!(!compound.supports_step_export());
    let refused = dir.path().join("unqualified-compound.step");
    assert!(compound.write_step(&refused).is_err());
    assert!(!refused.exists());
    let path = dir.path().join("entry.step");
    body.write_step(&path).unwrap();
    let text = std::fs::read_to_string(&path).unwrap();
    let proof = body.verify_entry_chamfer_step(&text).unwrap().unwrap();
    assert_eq!(proof.face_count, 8);
    assert!(body.verify_helical_step(&text).is_none());
    assert!(
        body.verify_entry_chamfer_step(&text.replace("CONICAL_SURFACE", "CYLINDRICAL_SURFACE"))
            .unwrap()
            .is_err()
    );
}

#[test]
fn placed_material_keeps_cached_references_and_embedded_facets_at_far_proper_poses() {
    use multicad_core::references::{GeometryIntent, Resolution, TopologyKind, resolve};
    use multicad_kernel::RigidTransform;
    use multicad_kernel::rust_kernel::thread_end_relief::Intent;
    for (left, starts, terminal, frame) in [
        (
            true,
            1,
            false,
            PlaneFrame {
                origin: [1e6, -2e6, 3e6],
                u: [0.8, 0.6, 0.],
                v: [-0.36, 0.48, 0.8],
            },
        ),
        (
            false,
            2,
            true,
            PlaneFrame {
                origin: [17., -3., 21.],
                u: [1., 0., 0.],
                v: [0., -1., 0.],
            },
        ),
    ] {
        let (spec, root) = request(left, starts);
        let feature = EntityId::new();
        let body = Solid::rust_helical_thread_entry_chamfer(
            spec,
            Some(root),
            0.3,
            terminal.then_some(Intent {
                width_mm: 0.6,
                radial_depth_mm: 0.1,
            }),
            feature,
            frame,
            Default::default(),
        )
        .unwrap();
        let snapshot = body.topology_snapshot(feature, 0).unwrap();
        let cap = snapshot
            .elements
            .iter()
            .find(|e| {
                e.key.kind == TopologyKind::Face && e.origins.iter().any(|o| o.role == "cap.end")
            })
            .unwrap();
        let expected =
            std::array::from_fn::<_, 3, _>(|i| frame.origin[i] + frame.normal()[i] * 6.6);
        for (a, b) in cap.frames[0].frame.origin_mm.into_iter().zip(expected) {
            assert!((a - b).abs() < 1e-7);
        }
        let intent = GeometryIntent::capture(&snapshot, cap.key).unwrap();
        let shared = body.share().unwrap();
        assert!(std::ptr::eq(
            body.rust_thread_entry_chamfer_properties().unwrap().0,
            shared.rust_thread_entry_chamfer_properties().unwrap().0
        ));
        assert!(std::ptr::eq(
            body.rust_thread_entry_chamfer_topology().unwrap().0,
            shared.rust_thread_entry_chamfer_topology().unwrap().0
        ));
        let posed = body
            .transformed(RigidTransform {
                rotation: [[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]],
                translation_mm: [7., -9., 11.],
            })
            .unwrap();
        assert!(std::ptr::eq(
            body.rust_thread_entry_chamfer_properties().unwrap().0,
            posed.rust_thread_entry_chamfer_properties().unwrap().0
        ));
        assert!(matches!(
            resolve(&intent, &posed.topology_snapshot(feature, 1).unwrap(), 1),
            Resolution::Resolved { .. }
        ));
        let display = body.topology_with_display(feature, 0, Some(0.1)).unwrap();
        let tight = body.display_edges(0.001).unwrap();
        assert!(
            tight.iter().map(|e| e.points_mm.len()).sum::<usize>() < 50_000,
            "Chord certificates must not subdivide by curve length"
        );
        assert_eq!(
            display.edges.len(),
            body.rust_thread_entry_chamfer_edges()
                .unwrap()
                .0
                .iter()
                .filter(|edge| !matches!(
                    body.rust_thread_entry_chamfer_topology().unwrap().0.edges()
                        [edge.source_edges()[0]]
                        .source(),
                    multicad_kernel::rust_kernel::thread_entry_chamfer::SourceCurve::NeckSeam { .. }
                ))
                .count()
        );
        assert!(
            !body
                .faces()
                .unwrap()
                .iter()
                .any(|f| !f.area_mm2.is_finite() || f.area_mm2 <= 0.)
        );
        let started = std::time::Instant::now();
        let mesh = body.manufacturing_mesh(MeshOptions::default());
        // Keep the real authored root knot; phase-shift only the artificial
        // sampling grid so its nearly coincident top cut cannot collapse here.
        let mesh = mesh.unwrap();
        eprintln!(
            "placed entry: {} facets in {:?}",
            mesh.indices.len() / 3,
            started.elapsed()
        );
        assert!(mesh.vertices.iter().all(|v| v.is_finite()));
        let before = body.bounds().unwrap();
        assert_eq!(before, body.bounds().unwrap());
    }
}

#[test]
fn tiny_cap_and_cancelled_or_invalid_construction_never_publish_partial_geometry() {
    use multicad_kernel::GeometryErrorKind;
    let (spec, root) = request(false, 1);
    let policy = GeometryPolicy::default();
    let tiny = Solid::rust_helical_thread_entry_chamfer(
        spec,
        Some(root),
        2.998,
        None,
        EntityId::new(),
        PlaneFrame::default(),
        policy,
    )
    .unwrap();
    assert!(tiny.faces().unwrap()[0].area_mm2 > 0.);
    for stop in [1, 40, 200] {
        let calls = std::cell::Cell::new(0);
        let result = Solid::rust_helical_thread_entry_chamfer_with_interrupt(
            spec,
            Some(root),
            0.3,
            None,
            EntityId::new(),
            PlaneFrame::default(),
            policy,
            Some(&|| {
                calls.set(calls.get() + 1);
                calls.get() >= stop
            }),
        );
        assert_eq!(result.err().unwrap().kind, GeometryErrorKind::Cancelled);
    }
    for c in [0., -1., f64::NAN, 3.1] {
        assert!(
            Solid::rust_helical_thread_entry_chamfer(
                spec,
                Some(root),
                c,
                None,
                EntityId::new(),
                PlaneFrame::default(),
                policy
            )
            .is_err()
        );
    }
    let retry = Solid::rust_helical_thread_entry_chamfer(
        spec,
        Some(root),
        0.3,
        None,
        EntityId::new(),
        PlaneFrame::default(),
        policy,
    )
    .unwrap();
    assert!(retry.volume_mm3().unwrap() > 0.);
}
