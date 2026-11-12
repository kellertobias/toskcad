use multicad_core::{
    commands::Session,
    document::{Document, EntityId},
    model::{ThreadFit, ThreadHand, ThreadMode, ThreadParameters, length},
    references::{GeometryIntent, Origin, Resolution, TopologyKind, resolve},
};
use multicad_kernel::{
    GeometryErrorKind, GeometryPolicy, HelicalThread, MeshOptions, PlaneFrame, Solid,
    rust_kernel::{Circle, Rational, cylinder::AxialCylinder, thread_end_relief::Intent},
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
fn complete_metric_receiver_retains_actual_material_cap_charts_mesh_and_step() {
    let policy = GeometryPolicy::default();
    for (left, starts, reversed, terminal) in [
        (false, 1, false, false),
        (true, 1, true, false),
        (false, 2, true, true),
        (true, 2, false, true),
    ] {
        let frame = PlaneFrame {
            origin: [100., -30., 50.],
            u: [0.8, 0.6, 0.],
            v: [-0.36, 0.48, 0.8],
        };
        let source = stock(1.7, frame);
        let before = source.topology_snapshot(EntityId::new(), 0).unwrap();
        let cut = PlaneFrame {
            origin: world(frame, [0., 0., if reversed { 3.3 } else { 0. }]),
            u: frame.u,
            v: frame.v.map(|x| if reversed { -x } else { x }),
        };
        let (spec, root) = tool(left, starts);
        let feature = EntityId::new();
        let tail = terminal.then_some(Intent {
            width_mm: 0.3,
            radial_depth_mm: 0.05,
        });
        let solid = source
            .rust_cut_external_thread_entry_chamfer_with_interrupt(
                spec, root, 0.2, tail, feature, cut, policy, None,
            )
            .unwrap();
        let data = solid.rust_thread_entry_chamfer_data().unwrap();
        assert_eq!(data.source.source_thread().spec().length_mm, 3.3);
        assert_eq!(data.source.full_length_mm(), 3.3);
        if let Some(tail) = data.source.terminal_relief() {
            assert_eq!(tail.full_length_mm(), 3.3);
            assert_eq!(tail.intent().unwrap().width_mm, 0.3);
            assert_eq!(data.source.material_length_mm(), 3.3 - 0.3);
        }
        assert_eq!(data.source.source_thread().spec().starts, starts);
        assert_eq!(data.source.source_thread().spec().left_hand, left);
        assert_eq!(data.source.terminal_relief().is_some(), terminal);
        let intervals = solid.rust_receiving_thread_intervals().unwrap();
        assert_eq!(intervals.requested_tool_length_mm, 3.3);
        assert_eq!(intervals.stock_length_mm, 3.3);
        assert_eq!(
            intervals.physical_thread_length_mm,
            data.source.material_length_mm()
        );
        let snapshot = solid.topology_snapshot(feature, 1).unwrap();
        for original in before
            .elements
            .iter()
            .filter(|e| e.key.kind == TopologyKind::Face && e.key.ordinal < 2)
        {
            let mut intent = GeometryIntent::capture(&before, original.key).unwrap();
            intent.boundary = None;
            let Resolution::Resolved { key, .. } = resolve(&intent, &snapshot, 1) else {
                panic!("original cap plane ancestry")
            };
            let current = snapshot.elements.iter().find(|e| e.key == key).unwrap();
            for anchor in &original.frames {
                let actual = current
                    .frames
                    .iter()
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
        // No original full-radius rim or removed cylindrical wall is renamed
        // to a different finite boundary just because its cap plane survives.
        for original in before.elements.iter().filter(|e| {
            e.key.kind == TopologyKind::Edge
                || (e.key.kind == TopologyKind::Face && e.key.ordinal == 2)
        }) {
            assert!(
                snapshot
                    .elements
                    .iter()
                    .filter(|e| e.key.kind == original.key.kind)
                    .all(|e| e.origins.is_disjoint(&original.origins))
            );
        }
        assert!(
            solid
                .rust_thread_contains(world(intervals.frame, [0., 0., 1.]))
                .unwrap()
        );
        assert!(
            !solid
                .rust_thread_contains(world(intervals.frame, [1.69, 0., 1.]))
                .unwrap()
        );
        let fresh = Solid::rust_helical_thread_entry_chamfer(
            spec,
            root,
            0.2,
            tail,
            feature,
            intervals.frame,
            policy,
        )
        .unwrap();
        assert_eq!(solid.volume_mm3().unwrap(), fresh.volume_mm3().unwrap());
        assert!(solid.volume_mm3().unwrap() < source.volume_mm3().unwrap());
        let receipt = solid
            .manufacturing_mesh_with_receipt(MeshOptions {
                linear_deflection_mm: 0.1,
                angular_deflection_rad: 0.3,
                ..Default::default()
            })
            .unwrap();
        assert!(receipt.edge_pairs > 0 && receipt.vertex_links > 0);
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("receiving.step");
        solid.write_step(&path).unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        let proof = solid.verify_entry_chamfer_step(&text).unwrap().unwrap();
        assert!(
            proof.source_volume_bounds_mm3[0] <= solid.volume_mm3().unwrap()
                && solid.volume_mm3().unwrap() <= proof.source_volume_bounds_mm3[1]
        );
        assert_eq!(source.topology_snapshot(before.owner, 0).unwrap(), before);
    }
}
#[test]
fn partial_overrun_offaxis_and_cancelled_cuts_refuse_without_source_changes() {
    let policy = GeometryPolicy::default();
    let source = stock(1.5, PlaneFrame::default());
    let before = source.topology_snapshot(EntityId::new(), 0).unwrap();
    let (spec, root) = tool(false, 1);
    for length in [3., 3.6] {
        let error = source
            .rust_cut_external_thread_entry_chamfer_with_interrupt(
                HelicalThread {
                    length_mm: length,
                    ..spec
                },
                root,
                0.2,
                None,
                EntityId::new(),
                PlaneFrame::default(),
                policy,
                None,
            )
            .err()
            .unwrap();
        assert_eq!(error.kind, GeometryErrorKind::Unsupported);
        assert!(error.message.contains("complete literal stock height"));
    }
    for frame in [
        PlaneFrame {
            origin: [0.01, 0., 0.],
            ..PlaneFrame::default()
        },
        PlaneFrame {
            u: [0.9999995, 0., 0.001],
            v: [0., 1., 0.],
            ..PlaneFrame::default()
        },
    ] {
        assert!(
            source
                .rust_cut_external_thread_entry_chamfer_with_interrupt(
                    spec,
                    root,
                    0.2,
                    None,
                    EntityId::new(),
                    frame,
                    policy,
                    None
                )
                .is_err()
        );
    }
    let error = source
        .rust_cut_external_thread_entry_chamfer_with_interrupt(
            spec,
            root,
            0.2,
            None,
            EntityId::new(),
            PlaneFrame::default(),
            policy,
            Some(&|| true),
        )
        .err()
        .unwrap();
    assert_eq!(error.kind, GeometryErrorKind::Cancelled);
    let calls = std::cell::Cell::new(0);
    let cancel_after_qualification = || {
        calls.set(calls.get() + 1);
        calls.get() >= 2
    };
    let error = source
        .rust_cut_external_thread_entry_chamfer_with_interrupt(
            spec,
            root,
            0.2,
            None,
            EntityId::new(),
            PlaneFrame::default(),
            policy,
            Some(&cancel_after_qualification),
        )
        .err()
        .unwrap();
    assert_eq!(error.kind, GeometryErrorKind::Cancelled);
    assert_eq!(calls.get(), 2);
    assert!(
        stock(1., PlaneFrame::default())
            .rust_cut_external_thread_entry_chamfer_with_interrupt(
                spec,
                root,
                0.2,
                None,
                EntityId::new(),
                PlaneFrame::default(),
                policy,
                None
            )
            .is_err()
    );
    assert_eq!(source.topology_snapshot(before.owner, 0).unwrap(), before);
}
