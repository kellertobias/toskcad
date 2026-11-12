use multicad_core::{
    commands::Session,
    document::{Document, EntityId},
    model::{ThreadFit, ThreadHand, ThreadMode, ThreadParameters, length},
    references::{GeometryIntent, Origin, Resolution, TopologyKind, resolve},
};
use multicad_kernel::{
    GeometryErrorKind, GeometryPolicy, HelicalThread, MeshOptions, PlaneFrame, RigidTransform,
    Solid,
    rust_kernel::{Circle, Rational, cylinder::AxialCylinder},
};
fn tool(left: bool, starts: u32) -> (HelicalThread, Option<f64>) {
    let session = Session::new(Document::new("Metric receiving source").unwrap()).unwrap();
    let t = ThreadParameters {
        major_diameter: length(3.),
        pitch: length(0.5),
        length: length(3.3),
        starts,
        hand: if left {
            ThreadHand::Left
        } else {
            ThreadHand::Right
        },
        fit: if starts == 1 {
            ThreadFit::Metric6gExternal
        } else {
            ThreadFit::NominalBasic
        },
        ..Default::default()
    }
    .resolve(session.parameters(), ThreadMode::CutExternal)
    .unwrap();
    (
        HelicalThread {
            major_radius_mm: t.major_radius_mm,
            pitch_mm: t.pitch_mm,
            length_mm: t.length_mm,
            depth_mm: t.profile.depth_mm,
            crest_width_mm: t.profile.crest_width_mm,
            lower_flank_deg: t.profile.lower_flank_deg,
            upper_flank_deg: t.profile.upper_flank_deg,
            radial_taper: 0.,
            starts,
            left_hand: left,
        },
        t.metric_fit.and_then(|f| f.external_root_radius_mm),
    )
}
fn world(f: PlaneFrame, p: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|i| f.origin[i] + f.u[i] * p[0] + f.v[i] * p[1] + f.normal()[i] * p[2])
}
fn stock(radius: f64, frame: PlaneFrame) -> Solid {
    let n = |v: f64| Rational::parse_decimal(&v.to_string()).unwrap();
    let id = EntityId::new();
    Solid::from_axial_cylinder_sweep(
        AxialCylinder::new(
            id.uuid_bits(),
            Circle::new([n(0.); 2], n(radius)).unwrap(),
            n(0.),
            n(3.3),
        )
        .unwrap(),
        frame,
        [0.; 2],
        false,
        Origin {
            source: id,
            role: "circle".into(),
        },
        GeometryPolicy::default(),
    )
    .unwrap()
}
#[test]
fn terminal_only_receiver_retains_actual_neck_cap_ancestry_and_independent_step() {
    let policy = GeometryPolicy::default();
    for (left, starts, reversed) in [
        (false, 1, false),
        (true, 1, true),
        (false, 2, true),
        (true, 2, false),
    ] {
        let frame = PlaneFrame {
            origin: [100., -30., 50.],
            u: [0.8, 0.6, 0.],
            v: [-0.36, 0.48, 0.8],
        };
        let source = stock(1.7, frame);
        let old = source.topology_snapshot(EntityId::new(), 0).unwrap();
        let cut = PlaneFrame {
            origin: world(frame, [0., 0., if reversed { 3.3 } else { 0. }]),
            u: frame.u,
            v: frame.v.map(|x| if reversed { -x } else { x }),
        };
        let (spec, root) = tool(left, starts);
        let feature = EntityId::new();
        let solid = source
            .rust_cut_external_thread_end_relief_with_interrupt(
                spec, root, 0.3, 0.05, feature, cut, policy, None,
            )
            .unwrap();
        let (relief, placed) = solid.rust_thread_end_relief_data().unwrap();
        assert_eq!(relief.source_thread().spec().length_mm, 3.3);
        assert_eq!(relief.source_thread().spec().starts, starts);
        assert_eq!(relief.source_thread().spec().left_hand, left);
        assert_eq!(relief.full_length_mm(), 3.3);
        assert_eq!(relief.material_length_mm(), 3.3 - 0.3);
        assert_eq!(relief.intent().unwrap().width_mm, 0.3);
        assert_eq!(relief.intent().unwrap().radial_depth_mm, 0.05);
        let intervals = solid.rust_receiving_thread_intervals().unwrap();
        assert_eq!(intervals.frame, placed);
        assert_eq!(intervals.requested_tool_length_mm, 3.3);
        assert_eq!(intervals.stock_length_mm, 3.3);
        assert_eq!(
            intervals.physical_thread_length_mm,
            relief.material_length_mm()
        );
        let current = solid.topology_snapshot(feature, 1).unwrap();
        for original in old
            .elements
            .iter()
            .filter(|e| e.key.kind == TopologyKind::Face && e.key.ordinal < 2)
        {
            let mut intent = GeometryIntent::capture(&old, original.key).unwrap();
            intent.boundary = None;
            let Resolution::Resolved { key, .. } = resolve(&intent, &current, 1) else {
                panic!("actual cap ancestry")
            };
            let actual = current.elements.iter().find(|e| e.key == key).unwrap();
            for anchor in &original.frames {
                let saved = actual
                    .frames
                    .iter()
                    .find(|a| a.origin == anchor.origin)
                    .unwrap();
                for (a, b) in saved
                    .frame
                    .origin_mm
                    .into_iter()
                    .zip(anchor.frame.origin_mm)
                    .chain(saved.frame.u.into_iter().zip(anchor.frame.u))
                    .chain(saved.frame.v.into_iter().zip(anchor.frame.v))
                {
                    assert!((a - b).abs() < policy.linear_tolerance_mm / 16.);
                }
            }
        }
        let pose = RigidTransform {
            rotation: [[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]],
            translation_mm: [7., -11., 3.],
        };
        let moved_stock = source.transformed(pose).unwrap();
        let moved = solid.transformed(pose).unwrap();
        let moved_original = moved_stock.topology_snapshot(old.owner, 2).unwrap();
        let moved_topology = moved.topology_snapshot(feature, 2).unwrap();
        for original in moved_original
            .elements
            .iter()
            .filter(|e| e.key.kind == TopologyKind::Face && e.key.ordinal < 2)
        {
            for anchor in &original.frames {
                let actual = moved_topology
                    .elements
                    .iter()
                    .flat_map(|e| &e.frames)
                    .find(|a| a.origin == anchor.origin)
                    .unwrap();
                for (a, b) in actual
                    .frame
                    .origin_mm
                    .into_iter()
                    .zip(anchor.frame.origin_mm)
                    .chain(actual.frame.u.into_iter().zip(anchor.frame.u))
                    .chain(actual.frame.v.into_iter().zip(anchor.frame.v))
                {
                    assert!((a - b).abs() < policy.linear_tolerance_mm / 16.);
                }
            }
        }
        for removed in old.elements.iter().filter(|e| {
            e.key.kind == TopologyKind::Edge
                || (e.key.kind == TopologyKind::Face && e.key.ordinal == 2)
        }) {
            assert!(
                current
                    .elements
                    .iter()
                    .filter(|e| e.key.kind == removed.key.kind)
                    .all(|e| e.origins.is_disjoint(&removed.origins))
            );
        }
        let fresh =
            Solid::rust_helical_thread_end_relief(spec, root, 0.3, 0.05, feature, placed, policy)
                .unwrap();
        assert_eq!(solid.volume_mm3().unwrap(), fresh.volume_mm3().unwrap());
        assert!(solid.volume_mm3().unwrap() < source.volume_mm3().unwrap());
        let neck = relief.neck_radius_mm();
        assert!(
            solid
                .rust_thread_contains(world(placed, [neck - 0.01, 0., 3.15]))
                .unwrap()
        );
        assert!(
            !solid
                .rust_thread_contains(world(placed, [neck + 0.01, 0., 3.15]))
                .unwrap()
        );
        assert!(
            !solid
                .rust_thread_contains(world(placed, [0., 0., 3.4]))
                .unwrap()
        );
        let receipt = solid
            .manufacturing_mesh_with_receipt(MeshOptions {
                linear_deflection_mm: 0.05,
                angular_deflection_rad: 0.3,
                ..Default::default()
            })
            .unwrap();
        assert!(receipt.edge_pairs > 0 && receipt.vertex_links > 0);
        let dir = tempfile::tempdir().unwrap();
        let file = dir.path().join("terminal.step");
        solid.write_step(&file).unwrap();
        let text = std::fs::read_to_string(&file).unwrap();
        let decoded_volume = solid.verify_helical_step(&text).unwrap().unwrap();
        assert!((decoded_volume - solid.volume_mm3().unwrap()).abs() < 1e-5);
        let imported = Solid::import_rust_step(&text, EntityId::new(), policy).unwrap();
        let (decoded, _) = imported.rust_thread_end_relief_data().unwrap();
        assert!(decoded.intent().is_none());
        assert!((decoded.full_length_mm() - 3.3).abs() < policy.linear_tolerance_mm);
        assert!(
            (decoded.material_length_mm() - relief.material_length_mm()).abs()
                < policy.linear_tolerance_mm
        );
        assert!((imported.volume_mm3().unwrap() - solid.volume_mm3().unwrap()).abs() < 1e-5);
        assert!(
            !imported
                .rust_thread_contains(world(placed, [neck + 0.01, 0., 3.15]))
                .unwrap()
        );
        assert_eq!(source.topology_snapshot(old.owner, 0).unwrap(), old);
    }
}
#[test]
fn malformed_partial_overrun_offaxis_and_cancelled_terminal_cuts_are_atomic() {
    let policy = GeometryPolicy::default();
    let source = stock(1.5, PlaneFrame::default());
    let old = source.topology_snapshot(EntityId::new(), 0).unwrap();
    let (spec, root) = tool(false, 1);
    for (length, width, depth, frame) in [
        (3., 0.3, 0.05, PlaneFrame::default()),
        (3.6, 0.3, 0.05, PlaneFrame::default()),
        (3.3, 0., 0.05, PlaneFrame::default()),
        (3.3, 4., 0.05, PlaneFrame::default()),
        (3.3, 0.3, 0., PlaneFrame::default()),
        (3.3, 0.3, 1.5, PlaneFrame::default()),
        (
            3.3,
            0.3,
            0.05,
            PlaneFrame {
                origin: [0.01, 0., 0.],
                ..Default::default()
            },
        ),
    ] {
        assert!(
            source
                .rust_cut_external_thread_end_relief_with_interrupt(
                    HelicalThread {
                        length_mm: length,
                        ..spec
                    },
                    root,
                    width,
                    depth,
                    EntityId::new(),
                    frame,
                    policy,
                    None
                )
                .is_err()
        );
        assert_eq!(source.topology_snapshot(old.owner, 0).unwrap(), old);
    }
    for after in [1, 2] {
        let calls = std::cell::Cell::new(0);
        let cancel = || {
            calls.set(calls.get() + 1);
            calls.get() >= after
        };
        let error = source
            .rust_cut_external_thread_end_relief_with_interrupt(
                spec,
                root,
                0.3,
                0.05,
                EntityId::new(),
                PlaneFrame::default(),
                policy,
                Some(&cancel),
            )
            .err()
            .unwrap();
        assert_eq!(error.kind, GeometryErrorKind::Cancelled);
        assert_eq!(calls.get(), after);
        assert_eq!(source.topology_snapshot(old.owner, 0).unwrap(), old);
    }
}
