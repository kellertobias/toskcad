use multicad_core::{
    document::EntityId,
    references::{GeometryIntent, Origin, Resolution, TopologyKind, resolve},
};
use multicad_kernel::{
    CircularPocketTool, GeometryPolicy, PlaneFrame, Solid, rust_kernel::Rational,
};
use multicad_precision::PortableMath as _;
fn point(f: PlaneFrame, p: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|i| f.origin[i] + f.u[i] * p[0] + f.v[i] * p[1] + f.normal()[i] * p[2])
}
fn tool(entry: PlaneFrame, depth: f64) -> CircularPocketTool {
    let n = entry.normal();
    CircularPocketTool {
        frame: PlaneFrame {
            origin: std::array::from_fn(|i| entry.origin[i] - n[i] * 0.01),
            ..entry
        },
        center_uv: [0.; 2],
        radius_mm: 2.,
        displacement_mm: n.map(|x| x * (depth + 0.01)),
        origin: Origin {
            source: EntityId::new(),
            role: "connection.receiver/bore".into(),
        },
    }
}
fn stock(frame: PlaneFrame, owner: EntityId) -> Solid {
    Solid::rust_box(
        owner.uuid_bits(),
        [Rational::new(20, 1).unwrap(); 3],
        frame,
        GeometryPolicy::default(),
    )
    .unwrap()
}
#[test]
fn guarded_blind_bores_keep_actual_caps_floors_material_and_exchange_in_proper_stock_frames() {
    let policy = GeometryPolicy::default();
    let dir = tempfile::tempdir().unwrap();
    for frame in [
        PlaneFrame::default(),
        PlaneFrame {
            origin: [100., -30., 50.],
            u: [0.8, 0.6, 0.],
            v: [-0.36, 0.48, 0.8],
        },
    ] {
        for axis in 0..3 {
            for reverse in [false, true] {
                let owner = EntityId::new();
                let stock = stock(frame, owner);
                let before = stock.topology_snapshot(owner, 0).unwrap();
                let cap = before
                    .elements
                    .iter()
                    .find(|e| e.key.kind == TopologyKind::Face)
                    .unwrap();
                let intent = GeometryIntent::capture(&before, cap.key).unwrap();
                let mut center = [10.; 3];
                center[axis] = if reverse { 20. } else { 0. };
                let axes = [frame.u, frame.v, frame.normal()];
                let n = axes[axis].map(|x| if reverse { -x } else { x });
                let u = axes[(axis + 1) % 3];
                let v = [
                    n[1] * u[2] - n[2] * u[1],
                    n[2] * u[0] - n[0] * u[2],
                    n[0] * u[1] - n[1] * u[0],
                ];
                let entry = PlaneFrame {
                    origin: point(frame, center),
                    u: std::array::from_fn(|i| 0.8 * u[i] + 0.6 * v[i]),
                    v: std::array::from_fn(|i| -0.6 * u[i] + 0.8 * v[i]),
                };
                let (cut, mouth) = stock
                    .rust_connection_blind_bore(&tool(entry, 4.), 1., 1., policy)
                    .unwrap();
                for (actual, expected) in mouth.into_iter().zip(entry.origin) {
                    assert!((actual - expected).abs() < 1e-10);
                }
                assert!(
                    (cut.volume_mm3().unwrap() - (8000. - 16. * std::f64::consts::PI)).abs() < 1e-8
                );
                assert_eq!(cut.solid_count().unwrap(), 1);
                assert_eq!(cut.faces().unwrap().len(), 8);
                let floor = cut
                    .faces()
                    .unwrap()
                    .into_iter()
                    .find(|f| (f.area_mm2 - 4. * std::f64::consts::PI).abs() < 1e-8)
                    .unwrap();
                let expected = point(entry, [0., 0., 4.]);
                for (actual, expected) in floor.center_mm.into_iter().zip(expected) {
                    assert!((actual - expected).abs() < 1e-10);
                }
                let normal = floor.plane.unwrap().1;
                assert!(normal.into_iter().zip(n).map(|(a, b)| a * b).sum::<f64>() < -0.999999999);
                let topology = cut.topology_snapshot(owner, 1).unwrap();
                assert!(matches!(
                    resolve(&intent, &topology, 1),
                    Resolution::Resolved { .. }
                ));
                let receipt = cut
                    .manufacturing_mesh_with_receipt(Default::default())
                    .unwrap();
                assert!(receipt.edge_pairs > 0 && receipt.vertex_links > 0);
                let path = dir.path().join("actual-blind.step");
                cut.write_step(&path).unwrap();
                let back = Solid::read_step(&path, policy).unwrap();
                assert!((back.volume_mm3().unwrap() - cut.volume_mm3().unwrap()).abs() < 1e-7);
                assert_eq!(stock.topology_snapshot(owner, 0).unwrap(), before);
            }
        }
    }
}
#[test]
fn declared_wall_floor_prior_cuts_and_oblique_tools_cannot_publish_a_blind_receiver() {
    let policy = GeometryPolicy::default();
    let owner = EntityId::new();
    let stock = stock(PlaneFrame::default(), owner);
    let before = stock.topology_snapshot(owner, 0).unwrap();
    let entry = PlaneFrame {
        origin: [10., 10., 20.],
        u: [1., 0., 0.],
        v: [0., -1., 0.],
    };
    for case in 0..7 {
        let mut t = tool(entry, 4.);
        let (mut wall, mut floor) = (1., 1.);
        let mut requested = policy;
        match case {
            0 => wall = 20.,
            1 => floor = 20.,
            2 => t.radius_mm = -1.,
            3 => requested.minimum_feature_mm *= 2.,
            4 => {
                let a: f64 = 0.01;
                let tilted = PlaneFrame {
                    origin: entry.origin,
                    u: [a.portable_cos(), 0., a.portable_sin()],
                    v: entry.v,
                };
                t = tool(tilted, 4.);
            }
            5 => t.origin.source = EntityId::nil(),
            _ => t = tool(entry, 19.5),
        };
        assert!(
            stock
                .rust_connection_blind_bore(&t, wall, floor, requested)
                .is_err(),
            "case{case}"
        );
        assert_eq!(stock.topology_snapshot(owner, 0).unwrap(), before);
    }
    let prior = stock
        .rust_drill_circular_pockets(
            &[tool(
                PlaneFrame {
                    origin: [5., 5., 20.],
                    ..entry
                },
                1.,
            )],
            policy,
        )
        .unwrap();
    let prior_before = prior.topology_snapshot(owner, 1).unwrap();
    assert!(
        prior
            .rust_connection_blind_bore(&tool(entry, 4.), 1., 1., policy)
            .is_err()
    );
    assert_eq!(prior.topology_snapshot(owner, 1).unwrap(), prior_before);
    let detached = Solid::rust_box(
        EntityId::new().uuid_bits(),
        [Rational::new(20, 1).unwrap(); 3],
        PlaneFrame {
            origin: [30., 0., 0.],
            ..Default::default()
        },
        GeometryPolicy::default(),
    )
    .unwrap();
    let disconnected = stock
        .rust_union_orthogonal_tools(&[detached], policy)
        .unwrap()
        .unwrap();
    assert_eq!(disconnected.solid_count().unwrap(), 2);
    let disconnected_before = disconnected.topology_snapshot(owner, 0).unwrap();
    assert!(
        disconnected
            .rust_connection_blind_bore(&tool(entry, 4.), 1., 1., policy)
            .is_err()
    );
    assert_eq!(
        disconnected.topology_snapshot(owner, 0).unwrap(),
        disconnected_before
    );
}

#[test]
fn bounding_boxes_cannot_replace_actual_declared_wall_and_floor_material() {
    let policy = GeometryPolicy::default();
    let owner = EntityId::new();
    let solid = stock(PlaneFrame::default(), owner);
    let entry = PlaneFrame {
        origin: [10., 10., 20.],
        u: [1., 0., 0.],
        v: [0., -1., 0.],
    };
    for (origin, size) in [
        (
            [9., 9., 15.],
            [
                Rational::new(2, 1).unwrap(),
                Rational::new(2, 1).unwrap(),
                Rational::new(1, 2).unwrap(),
            ],
        ),
        (
            [12.5, 9., 17.],
            [
                Rational::new(1, 4).unwrap(),
                Rational::new(2, 1).unwrap(),
                Rational::new(2, 1).unwrap(),
            ],
        ),
    ] {
        let cavity = Solid::rust_box(
            EntityId::new().uuid_bits(),
            size,
            PlaneFrame {
                origin,
                ..Default::default()
            },
            policy,
        )
        .unwrap();
        let missing = solid
            .rust_subtract_orthogonal_tools(&[cavity], policy)
            .unwrap()
            .unwrap();
        assert_eq!(missing.bounds().unwrap(), solid.bounds().unwrap());
        assert_eq!(missing.solid_count().unwrap(), 1);
        let before = missing.topology_snapshot(owner, 0).unwrap();
        let drill = tool(entry, 4.);
        // Intrinsic pocket separation permits this cut. Its larger declared wall
        // or floor is missing real material, so the receiving qualification refuses.
        assert!(
            missing
                .rust_drill_circular_pockets(std::slice::from_ref(&drill), policy)
                .is_ok()
        );
        assert!(
            missing
                .rust_connection_blind_bore(&drill, 1., 1., policy)
                .is_err()
        );
        assert_eq!(missing.topology_snapshot(owner, 0).unwrap(), before);
    }
}
