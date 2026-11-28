//! CAD-424: a profile sketched on a flat face of a box whose twelve edges
//! are all filleted joins outward and cuts inward exactly.
//!
//! The independent oracle is closed form: the rounded box is the Minkowski
//! sum of the inner box with a ball, and the tool adds or removes its
//! profile area times the extruded height. Every result is one closed
//! shell with exact cylindrical walls and reads back from STEP as the same
//! solid. Profiles reaching the fillets are refused in plain words.
use multicad_core::{
    document::EntityId,
    references::{Origin, TopologyKind},
};
use multicad_kernel::{
    BooleanOperation, CurvedEdge2d, CurvedProfile, EdgeTreatment, GeometryErrorKind,
    GeometryPolicy, PlaneFrame, Solid, TreatmentEdge,
    rust_kernel::{Circle, Rational, cylinder::AxialCylinder},
};
use std::f64::consts::PI;

const SIZE: [f64; 3] = [40., 30., 20.];
const R: f64 = 3.;

fn policy() -> GeometryPolicy {
    GeometryPolicy::default()
}

/// Closed-form volume of the box with every edge rounded by `R`.
fn rounded_volume() -> f64 {
    let i = SIZE.map(|s| s - 2. * R);
    i[0] * i[1] * i[2]
        + 2. * R * (i[0] * i[1] + i[1] * i[2] + i[0] * i[2])
        + PI * R * R * (i[0] + i[1] + i[2])
        + 4. / 3. * PI * R.powi(3)
}

fn rounded() -> Solid {
    let owner = EntityId::new();
    let feature = EntityId::new();
    let source = Solid::box_named(SIZE, owner, policy()).unwrap();
    let edges: Vec<_> = source
        .topology_snapshot(owner, 0)
        .unwrap()
        .elements
        .iter()
        .filter(|e| e.key.kind == TopologyKind::Edge)
        .map(|e| TreatmentEdge {
            key: e.key,
            origin: Origin {
                source: feature,
                role: format!("edge/{}", e.key.ordinal),
            },
        })
        .collect();
    assert_eq!(edges.len(), 12);
    let solid = source
        .treat_edges_named(
            &edges,
            EdgeTreatment::Fillet { radius_mm: R },
            feature,
            policy(),
        )
        .unwrap();
    near(solid.volume_mm3().unwrap(), rounded_volume());
    assert_eq!(solid.faces().unwrap().len(), 26);
    solid
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-9 * expected.abs().max(1.),
        "{actual} != {expected}"
    );
}

fn top() -> PlaneFrame {
    PlaneFrame {
        origin: [0., 0., SIZE[2]],
        ..Default::default()
    }
}

/// The +x side face, with sketch `u` along y and `v` along z.
fn side() -> PlaneFrame {
    PlaneFrame {
        origin: [SIZE[0], 0., 0.],
        u: [0., 1., 0.],
        v: [0., 0., 1.],
    }
}

/// A profile region on `frame` swept by `height`; negative heights sweep
/// into the body.
fn sweep(frame: PlaneFrame, outer: Vec<CurvedEdge2d>, height: f64) -> Solid {
    let owner = EntityId::new();
    let origins: Vec<_> = (0..outer.len())
        .map(|i| Origin {
            source: owner,
            role: format!("curve/{i}"),
        })
        .collect();
    let profile = CurvedProfile {
        frame,
        outer,
        holes: vec![],
    };
    Solid::rust_extrude_normal_profile_named(&profile, height, owner, &origins, policy()).unwrap()
}

/// The circular cylinder the Extrude command builds from a full circle on
/// the top face.
fn cylinder(center: [i64; 2], radius: i64, bottom: i64, height: i64) -> Solid {
    let int = |n| Rational::new(n, 1).unwrap();
    let circle = Circle::new(center.map(int), int(radius)).unwrap();
    let geometry = AxialCylinder::new(7, circle, int(bottom), int(height)).unwrap();
    let frame = PlaneFrame {
        origin: [center[0] as f64, center[1] as f64, bottom as f64],
        ..Default::default()
    };
    Solid::from_axial_cylinder(geometry, frame, policy()).unwrap()
}

fn circle(center: [f64; 2], radius: f64) -> Vec<CurvedEdge2d> {
    vec![CurvedEdge2d::Circle {
        center,
        radius_mm: radius,
    }]
}

/// A D-shaped region: a half disk of radius 5 about `(20, 15)` closed by its
/// diameter, in the face's own coordinates.
fn half_disk() -> (Vec<CurvedEdge2d>, f64) {
    (
        vec![
            CurvedEdge2d::Line {
                a: [15., 15.],
                b: [25., 15.],
            },
            CurvedEdge2d::Arc {
                center: [20., 15.],
                a: [25., 15.],
                b: [15., 15.],
                clockwise: false,
            },
        ],
        PI * 25. / 2.,
    )
}

fn curved_faces(solid: &Solid) -> usize {
    solid
        .faces()
        .unwrap()
        .iter()
        .filter(|f| f.plane.is_none())
        .count()
}

fn planar_faces(solid: &Solid, axis: usize, level: f64, up: bool) -> Vec<f64> {
    solid
        .faces()
        .unwrap()
        .into_iter()
        .filter_map(|f| {
            let (p, n) = f.plane?;
            ((p[axis] - level).abs() < 1e-7
                && (n[axis].abs() - 1.).abs() < 1e-9
                && (n[axis] > 0.) == up)
                .then_some(f.area_mm2)
        })
        .collect()
}

/// The exchange file reads back as the same solid.
fn round_trips(solid: &Solid) {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("result.step");
    solid.write_step(&path).unwrap();
    let read = Solid::read_step(&path, policy()).unwrap();
    near(read.volume_mm3().unwrap(), solid.volume_mm3().unwrap());
    assert_eq!(read.faces().unwrap().len(), solid.faces().unwrap().len());
    assert_eq!(curved_faces(&read), curved_faces(solid));
}

fn checked(result: Result<Option<Solid>, multicad_kernel::GeometryError>, what: &str) -> Solid {
    let solid = result
        .unwrap_or_else(|e| panic!("{what}: {e}"))
        .unwrap_or_else(|| panic!("{what}: empty result"));
    assert!(solid.uses_rust_geometry(), "{what}");
    assert!(
        !solid
            .manufacturing_mesh_with_receipt(Default::default())
            .unwrap()
            .mesh
            .indices
            .is_empty(),
        "{what}"
    );
    solid
}

/// The flat top face keeps the inner rectangle of the rounded box.
const TOP: f64 = (SIZE[0] - 2. * R) * (SIZE[1] - 2. * R);

#[test]
fn a_circle_inside_the_flat_face_joins_outward_exactly() {
    let body = rounded();
    for (what, tool) in [
        ("sweep", sweep(top(), circle([20., 15.], 6.), 5.)),
        ("cylinder", cylinder([20, 15], 6, 20, 5)),
    ] {
        let union = checked(body.boolean(&tool, BooleanOperation::Union, policy()), what);
        near(
            union.volume_mm3().unwrap(),
            rounded_volume() + 36. * PI * 5.,
        );
        // The flat face keeps the rest of itself around the boss; the boss
        // adds its wall and its own top.
        let face = planar_faces(&union, 2, SIZE[2], true);
        assert_eq!(face.len(), 1, "{what}");
        near(face[0], TOP - 36. * PI);
        near(planar_faces(&union, 2, SIZE[2] + 5., true)[0], 36. * PI);
        assert!(planar_faces(&union, 2, SIZE[2], false).is_empty(), "{what}");
        assert_eq!(curved_faces(&union), 21, "{what}");
        assert_eq!(union.faces().unwrap().len(), 28, "{what}");
        round_trips(&union);
        // The operands commute.
        let swapped = checked(tool.boolean(&body, BooleanOperation::Union, policy()), what);
        near(
            swapped.volume_mm3().unwrap(),
            rounded_volume() + 36. * PI * 5.,
        );
    }
}

#[test]
fn a_circle_inside_the_flat_face_cuts_blind_and_through_exactly() {
    let body = rounded();
    for (what, tool, depth) in [
        ("blind sweep", sweep(top(), circle([20., 15.], 6.), -5.), 5.),
        ("blind cylinder", cylinder([20, 15], 6, 15, 5), 5.),
        // From the top face to 5 mm past the bottom face.
        ("through", sweep(top(), circle([20., 15.], 6.), -25.), 20.),
    ] {
        let cut = checked(
            body.boolean(&tool, BooleanOperation::Subtract, policy()),
            what,
        );
        near(
            cut.volume_mm3().unwrap(),
            rounded_volume() - 36. * PI * depth,
        );
        near(planar_faces(&cut, 2, SIZE[2], true)[0], TOP - 36. * PI);
        if depth < SIZE[2] {
            // The pocket floor faces up.
            let floor = planar_faces(&cut, 2, SIZE[2] - depth, true);
            assert_eq!(floor.len(), 1, "{what}");
            near(floor[0], 36. * PI);
            assert_eq!(cut.faces().unwrap().len(), 28, "{what}");
        } else {
            near(planar_faces(&cut, 2, 0., false)[0], TOP - 36. * PI);
            assert_eq!(cut.faces().unwrap().len(), 27, "{what}");
        }
        assert_eq!(curved_faces(&cut), 21, "{what}");
        round_trips(&cut);
    }
}

#[test]
fn line_and_arc_regions_and_side_faces_extrude_exactly() {
    let body = rounded();
    let (outer, area) = half_disk();
    let union = checked(
        body.boolean(
            &sweep(top(), outer.clone(), 4.),
            BooleanOperation::Union,
            policy(),
        ),
        "half disk join",
    );
    near(union.volume_mm3().unwrap(), rounded_volume() + area * 4.);
    // A planar wall on the diameter and a cylindrical one on the arc.
    assert_eq!(union.faces().unwrap().len(), 29);
    round_trips(&union);
    let cut = checked(
        body.boolean(
            &sweep(top(), outer, -4.),
            BooleanOperation::Subtract,
            policy(),
        ),
        "half disk cut",
    );
    near(cut.volume_mm3().unwrap(), rounded_volume() - area * 4.);
    round_trips(&cut);
    // A circle on the +x face, whose sketch normal points along x.
    let boss = checked(
        body.boolean(
            &sweep(side(), circle([15., 10.], 4.), 6.),
            BooleanOperation::Union,
            policy(),
        ),
        "side join",
    );
    near(boss.volume_mm3().unwrap(), rounded_volume() + 16. * PI * 6.);
    near(
        planar_faces(&boss, 0, SIZE[0], true)[0],
        (SIZE[1] - 2. * R) * (SIZE[2] - 2. * R) - 16. * PI,
    );
    round_trips(&boss);
    let pocket = checked(
        body.boolean(
            &sweep(side(), circle([15., 10.], 4.), -6.),
            BooleanOperation::Subtract,
            policy(),
        ),
        "side cut",
    );
    near(
        pocket.volume_mm3().unwrap(),
        rounded_volume() - 16. * PI * 6.,
    );
    round_trips(&pocket);
}

#[test]
fn a_second_extrusion_onto_the_joined_result_is_exact() {
    let body = rounded();
    let first = checked(
        body.boolean(
            &sweep(top(), circle([12., 15.], 4.), 5.),
            BooleanOperation::Union,
            policy(),
        ),
        "first boss",
    );
    let second = checked(
        first.boolean(
            &sweep(top(), circle([28., 15.], 4.), -6.),
            BooleanOperation::Subtract,
            policy(),
        ),
        "then a pocket",
    );
    near(
        second.volume_mm3().unwrap(),
        rounded_volume() + 16. * PI * 5. - 16. * PI * 6.,
    );
    round_trips(&second);
}

#[test]
fn profiles_reaching_the_fillets_are_refused_in_plain_words() {
    let body = rounded();
    // The flat face ends 3 mm inside each box edge, at y = 27 here.
    for (what, center, radius) in [
        ("into the edge fillet", [20., 24.], 6.),
        ("into the corner blend", [36., 26.], 3.),
        ("past the box edge", [20., 30.], 6.),
        ("touching the fillet seam", [20., 21.], 6.),
    ] {
        for (operation, height) in [
            (BooleanOperation::Union, 5.),
            (BooleanOperation::Subtract, -5.),
        ] {
            let tool = sweep(top(), circle(center, radius), height);
            let error = body
                .boolean(&tool, operation, policy())
                .err()
                .unwrap_or_else(|| panic!("{what} {operation:?} must be refused"));
            assert_eq!(
                error.kind,
                GeometryErrorKind::Unsupported,
                "{what}: {error}"
            );
            assert!(
                error.message.contains("fillets") || error.message.contains("face"),
                "{what}: {error}"
            );
        }
    }
}
