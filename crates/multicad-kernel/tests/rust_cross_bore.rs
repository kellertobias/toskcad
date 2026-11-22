//! Exact perpendicular cross-bores through Rust stepped shafts and plain
//! cylinders, through the public solid API.
use multicad_core::{
    document::EntityId,
    references::{Origin, TopologyKind},
};
use multicad_kernel::{
    BooleanOperation, GeometryErrorKind, GeometryPolicy, MeshOptions, PlaneFrame, RigidTransform,
    Solid,
    rust_kernel::{Circle, Rational, cylinder::AxialCylinder},
};
use multicad_precision::PortableMath as _;
use std::f64::consts::PI;

fn cylinder(radius: Rational, height: i64, frame: PlaneFrame) -> (Solid, EntityId) {
    let feature = EntityId::new();
    let size = |n| Rational::new(n, 1).unwrap();
    let cylinder = AxialCylinder::new(
        feature.uuid_bits(),
        Circle::new([size(0); 2], radius).unwrap(),
        size(0),
        size(height),
    )
    .unwrap();
    let solid = Solid::from_axial_cylinder_sweep(
        cylinder,
        frame,
        [0.; 2],
        false,
        Origin {
            source: feature,
            role: "authored circle".into(),
        },
        GeometryPolicy::default(),
    )
    .unwrap();
    (solid, feature)
}

/// A cylinder along +Y from `origin`.
fn along_y(radius: i64, height: i64, origin: [f64; 3]) -> Solid {
    let frame = PlaneFrame {
        origin,
        u: [1., 0., 0.],
        v: [0., 0., -1.],
    };
    cylinder(Rational::new(radius, 1).unwrap(), height, frame).0
}

/// A Ø5.5 tool along +Z from `origin`.
fn bore(origin: [f64; 3], height: i64) -> (Solid, EntityId) {
    let frame = PlaneFrame {
        origin,
        u: [1., 0., 0.],
        v: [0., 1., 0.],
    };
    cylinder(Rational::new(11, 4).unwrap(), height, frame)
}

/// Hub Ø18 × 30 on y ∈ [0, 30] joined to a Ø10 shank to y = 570.
fn screw() -> Solid {
    let hub = along_y(9, 30, [0.; 3]);
    let shank = along_y(5, 540, [0., 30., 0.]);
    Solid::union_coaxial_cylinders(&[&hub, &shank], EntityId::new(), GeometryPolicy::default())
        .unwrap()
        .unwrap()
}

fn cut(target: &Solid, tool: &Solid) -> Result<Option<Solid>, multicad_kernel::GeometryError> {
    target.boolean(tool, BooleanOperation::Subtract, GeometryPolicy::default())
}

/// Volume of two crossing cylinders of radii `a <= big` with meeting axes.
fn crossing(a: f64, big: f64) -> f64 {
    let n = 4096;
    let step = PI / f64::from(n);
    (0..n)
        .map(|k| {
            let t = -PI / 2. + (f64::from(k) + 0.5) * step;
            t.portable_cos().powi(2) * (big * big - a * a * t.portable_sin().powi(2)).sqrt()
        })
        .sum::<f64>()
        * step
        * 4.
        * a
        * a
}

#[test]
fn radial_through_bore_keeps_exact_volume_bounds_and_topology() {
    let screw = screw();
    let (tool, tool_feature) = bore([0., 15., -50.], 100);
    let result = cut(&screw, &tool).unwrap().unwrap();
    let expected = screw.volume_mm3().unwrap() - crossing(2.75, 9.);
    let volume = result.volume_mm3().unwrap();
    assert!(
        (volume - expected).abs() < 1e-10 * expected,
        "{volume} {expected}"
    );
    assert_eq!(result.solid_count().unwrap(), 1);
    assert_eq!(result.bounds().unwrap(), screw.bounds().unwrap());
    // One retained tool wall joins the shaft faces; the bored wall loses
    // both holes.
    let faces = result.faces().unwrap();
    let before = screw.faces().unwrap();
    assert_eq!(faces.len(), before.len() + 1);
    let lost: f64 = before.iter().map(|f| f.area_mm2).sum::<f64>()
        - faces
            .iter()
            .take(before.len())
            .map(|f| f.area_mm2)
            .sum::<f64>();
    assert!(lost > 2. * PI * 2.75 * 2.75);
    let tool_wall = faces.last().unwrap();
    assert!(tool_wall.plane.is_none());
    assert!((tool_wall.center_mm[1] - 15.).abs() < 1e-9);
    // Stable origins: the new face is the tool's side and both rims carry
    // the tool and the bored wall.
    let feature = EntityId::new();
    let topology = result.topology_snapshot(feature, 0).unwrap();
    let side = topology
        .elements
        .iter()
        .find(|e| e.key.kind == TopologyKind::Face && e.key.ordinal == tool_wall.index)
        .unwrap();
    assert!(side.origins.iter().any(|o| o.source == tool_feature));
    let rims: Vec<_> = topology
        .elements
        .iter()
        .filter(|e| e.key.kind == TopologyKind::Edge && e.geometry.is_none())
        .collect();
    assert_eq!(rims.len(), 2);
    assert!(
        rims.iter()
            .all(|rim| rim.origins.is_superset(&side.origins))
    );
    // A closed manufacturing mesh within the chord bound, and sampled rims.
    let receipt = result
        .manufacturing_mesh_with_receipt(MeshOptions::default())
        .unwrap();
    assert!(receipt.edge_pairs > 0);
    let edges = result.display_edges(0.05).unwrap();
    assert!(edges.len() >= rims.len());
    // The bore is empty, the hub around it is material.
    assert!(!result.rust_revolution_contains([0., 15., 0.]).unwrap());
    assert!(!result.rust_revolution_contains([1., 15., -8.]).unwrap());
    assert!(result.rust_revolution_contains([5., 15., 0.]).unwrap());
    assert!(result.rust_revolution_contains([0., 100., 0.]).unwrap());
    // STEP has no encoding for the quartic rims yet: typed refusal.
    assert!(!result.supports_step_export());
    let dir = std::env::temp_dir().join(format!("cross-bore-{}", EntityId::new()));
    let error = result.write_step(&dir.with_extension("step")).unwrap_err();
    assert_eq!(error.kind, GeometryErrorKind::Unsupported);
    // A rigid placement keeps the exact body.
    let moved = result
        .transformed(RigidTransform {
            rotation: [[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]],
            translation_mm: [10., 20., 30.],
        })
        .unwrap();
    assert!((moved.volume_mm3().unwrap() - volume).abs() < 1e-9);
    assert!(!moved.rust_revolution_contains([-5., 20., 30.]).unwrap());
    assert!(moved.rust_revolution_contains([-5., 25., 30.]).unwrap());
}

#[test]
fn blind_and_offset_bores_and_a_plain_cylinder_target() {
    let screw = screw();
    // Blind: from below, stopping 1 mm past the axis with a planar bottom.
    let blind = cut(&screw, &bore([0., 15., -50.], 51).0).unwrap().unwrap();
    let faces = blind.faces().unwrap();
    assert_eq!(faces.len(), screw.faces().unwrap().len() + 2);
    let bottom = faces.last().unwrap();
    let (point, normal) = bottom.plane.unwrap();
    assert_eq!(normal, [0., 0., -1.]);
    assert!((point[2] - 1.).abs() < 1e-12);
    assert!((bottom.area_mm2 - PI * 2.75 * 2.75).abs() < 1e-9);
    let removed = screw.volume_mm3().unwrap() - blind.volume_mm3().unwrap();
    assert!(removed > crossing(2.75, 9.) / 2. && removed < crossing(2.75, 9.));
    assert!(blind.rust_revolution_contains([0., 15., 2.]).unwrap());
    assert!(!blind.rust_revolution_contains([0., 15., 0.]).unwrap());
    blind
        .manufacturing_mesh_with_receipt(MeshOptions::default())
        .unwrap();
    // Offset through bore removes less than the radial one.
    let offset = cut(&screw, &bore([3., 15., -50.], 100).0).unwrap().unwrap();
    let removed = screw.volume_mm3().unwrap() - offset.volume_mm3().unwrap();
    assert!(removed > 0. && removed < crossing(2.75, 9.));
    // A plain cylinder is bored as its one-step revolution.
    let hub = along_y(9, 30, [0.; 3]);
    let bored = cut(&hub, &bore([0., 15., -50.], 100).0).unwrap().unwrap();
    let expected = PI * 81. * 30. - crossing(2.75, 9.);
    assert!((bored.volume_mm3().unwrap() - expected).abs() < 1e-10 * expected);
}

#[test]
fn misses_contacts_and_tangency() {
    let screw = screw();
    let volume = screw.volume_mm3().unwrap();
    // Apart, and in exact contact with the hub wall: unchanged.
    for origin in [[20., 15., -50.], [11.75, 15., -50.]] {
        let result = cut(&screw, &bore(origin, 100).0).unwrap().unwrap();
        assert_eq!(result.volume_mm3().unwrap(), volume);
    }
    // Internal tangency of the bore to the hub wall is refused explicitly.
    let Err(error) = cut(&screw, &bore([6.25, 15., -50.], 100).0) else {
        panic!("tangency must be refused")
    };
    assert_eq!(error.kind, GeometryErrorKind::Unsupported);
    assert!(error.message.contains("tangency"), "{}", error.message);
    // Crossing the shoulder between hub and shank is refused, too.
    let Err(error) = cut(&screw, &bore([0., 29., -50.], 100).0) else {
        panic!("a shoulder crossing must be refused")
    };
    assert_eq!(error.kind, GeometryErrorKind::Unsupported);
}

#[test]
fn drawings_see_through_the_bore_and_publish_its_rims() {
    use multicad_core::drawings::{VectorCurve, VectorRole};
    use multicad_kernel::{DrawingProjection, DrawingSolid, project_drawing};
    let screw = screw();
    let bored = cut(&screw, &bore([0., 15., -50.], 100).0).unwrap().unwrap();
    // A plate under the hub whose top edges at y = 13 and 17 run beneath
    // the bore.
    let size = |n| Rational::new(n, 1).unwrap();
    let plate = Solid::rust_box(
        7,
        [size(40), size(4), size(2)],
        PlaneFrame {
            origin: [-20., 13., -20.],
            ..PlaneFrame::default()
        },
        GeometryPolicy::default(),
    )
    .unwrap();
    let solids = [
        DrawingSolid {
            solid: &bored,
            hatch: None,
        },
        DrawingSolid {
            solid: &plate,
            hatch: None,
        },
    ];
    let top = DrawingProjection {
        hidden_lines: true,
        ..Default::default()
    };
    let vectors = project_drawing(&solids, top).unwrap();
    let window = (2.75f64 * 2.75 - 4.).sqrt();
    let plate_edge = |role: VectorRole, from: f64, to: f64| {
        vectors.iter().any(|vector| {
            vector.source_index == 1
                && vector.role == role
                && matches!(vector.curve, VectorCurve::Line { a, b }
                    if (a[1] - 13.).abs() < 1e-6 && (b[1] - 13.).abs() < 1e-6
                        && (a[0].min(b[0]) - from).abs() < 1e-6
                        && (a[0].max(b[0]) - to).abs() < 1e-6)
        })
    };
    // Seen through the bore, hidden by the hub around it.
    assert!(
        plate_edge(VectorRole::Visible, -window, window),
        "{vectors:#?}"
    );
    assert!(plate_edge(VectorRole::Hidden, window, 9.));
    assert!(plate_edge(VectorRole::Hidden, -9., -window));
    // The near rim projects onto the tool circle, as visible chords.
    let rim = vectors
        .iter()
        .filter(|vector| vector.source_index == 0 && vector.role == VectorRole::Visible)
        .filter(|vector| match vector.curve {
            VectorCurve::Line { a, b } => [a, b]
                .iter()
                .all(|p| (p[0].portable_hypot(p[1] - 15.) - 2.75).abs() < 1e-9),
            _ => false,
        })
        .count();
    assert!(rim >= 16, "{rim}");
    // An isometric view of the bored shaft alone also evaluates exactly.
    let (s, c) = (0.5f64.sqrt(), (1. / 3f64).sqrt());
    let iso = DrawingProjection {
        frame: PlaneFrame {
            origin: [0.; 3],
            u: [s, -s, 0.],
            v: [-c * s, -c * s, (2. / 3f64).sqrt()],
        },
        hidden_lines: true,
        ..Default::default()
    };
    let vectors = project_drawing(&solids[..1], iso).unwrap();
    assert!(vectors.iter().any(|v| v.role == VectorRole::Hidden));
    assert!(vectors.iter().any(|v| v.role == VectorRole::Visible));
}
