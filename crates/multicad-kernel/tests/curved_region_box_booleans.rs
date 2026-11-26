//! CAD-406 / CAD-407: a box joined with, or cut by, a prism extruded from a
//! sketch region on one of its faces whose boundary has circular arcs (a
//! circle clipped by the face). Results are exact Rust solids: one closed
//! shell with exact cylindrical walls, no face left at the contact, and a
//! volume equal to the closed-form box and region measures. Each result
//! reads back from its STEP file as the same solid.
use multicad_core::{document::EntityId, references::Origin};
use multicad_kernel::{
    BooleanOperation, CurvedEdge2d, CurvedProfile, GeometryErrorKind, GeometryPolicy, PlaneFrame,
    RigidTransform, Solid,
};
use multicad_precision::PortableMath as _;
use std::f64::consts::PI;

fn policy() -> GeometryPolicy {
    GeometryPolicy::default()
}

const SIZE: [f64; 3] = [40., 30., 20.];
const BOX: f64 = SIZE[0] * SIZE[1] * SIZE[2];

fn stock() -> Solid {
    Solid::box_solid(SIZE, policy()).unwrap()
}

/// The sketch frame on the box's top face.
fn top_frame() -> PlaneFrame {
    PlaneFrame {
        origin: [0., 0., SIZE[2]],
        ..Default::default()
    }
}

/// A region on `frame` swept by `height` (negative sweeps into the box).
fn prism(frame: PlaneFrame, outer: Vec<CurvedEdge2d>, height: f64) -> Solid {
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

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-9 * expected.abs().max(1.),
        "{actual} != {expected}"
    );
}

/// Area of the part of a disk of radius `r` beyond a chord at distance `d`
/// from its centre.
fn segment_area(r: f64, d: f64) -> f64 {
    r * r * (d / r).portable_acos() - d * (r * r - d * d).sqrt()
}

/// Area of the part of a disk of radius `r` with x > a and y > b relative to
/// its centre (a, b >= 0 and a^2 + b^2 < r^2).
fn corner_area(r: f64, a: f64, b: f64) -> f64 {
    let f = |x: f64| 0.5 * (x * (r * r - x * x).sqrt() + r * r * (x / r).portable_asin()) - b * x;
    f((r * r - b * b).sqrt()) - f(a)
}

struct Region {
    name: &'static str,
    outer: Vec<CurvedEdge2d>,
    area: f64,
}

/// A circle clipped by the 40 x 30 face, in the face's own coordinates.
fn regions() -> Vec<Region> {
    let arc = |center: [f64; 2], a: [f64; 2], b: [f64; 2]| CurvedEdge2d::Arc {
        center,
        a,
        b,
        clockwise: false,
    };
    let line = |a: [f64; 2], b: [f64; 2]| CurvedEdge2d::Line { a, b };
    // Crossing the edge y = 30 with the centre 3 mm inside it.
    let h = (64f64 - 9.).sqrt();
    // A circle about a point 3 mm inside the corner (40, 30).
    let k = (64f64 - 9.).sqrt();
    vec![
        Region {
            name: "inside",
            outer: vec![CurvedEdge2d::Circle {
                center: [20., 15.],
                radius_mm: 6.,
            }],
            area: PI * 36.,
        },
        Region {
            name: "one edge",
            outer: vec![
                line([20. + h, 30.], [20. - h, 30.]),
                arc([20., 27.], [20. - h, 30.], [20. + h, 30.]),
            ],
            area: PI * 64. - segment_area(8., 3.),
        },
        Region {
            name: "corner sector",
            outer: vec![
                line([40., 30.], [30., 30.]),
                arc([40., 30.], [30., 30.], [40., 20.]),
                line([40., 20.], [40., 30.]),
            ],
            area: PI * 100. / 4.,
        },
        Region {
            name: "around corner",
            outer: vec![
                line([40., 27. - k], [40., 30.]),
                line([40., 30.], [37. - k, 30.]),
                arc([37., 27.], [37. - k, 30.], [40., 27. - k]),
            ],
            area: PI * 64. - 2. * segment_area(8., 3.) + corner_area(8., 3., 3.),
        },
    ]
}

fn planar_faces_at(solid: &Solid, z: f64, up: bool) -> Vec<f64> {
    solid
        .faces()
        .unwrap()
        .into_iter()
        .filter_map(|f| {
            let (p, n) = f.plane?;
            ((p[2] - z).abs() < 1e-7 && n[0].abs() < 1e-9 && n[1].abs() < 1e-9 && (n[2] > 0.) == up)
                .then_some(f.area_mm2)
        })
        .collect()
}

fn curved_faces(solid: &Solid) -> usize {
    solid
        .faces()
        .unwrap()
        .iter()
        .filter(|f| f.plane.is_none())
        .count()
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

#[test]
fn an_outward_join_merges_the_region_prism_into_the_box() {
    for region in regions() {
        let what = format!("join {}", region.name);
        let boss = prism(top_frame(), region.outer.clone(), 5.);
        let union = checked(
            stock().boolean(&boss, BooleanOperation::Union, policy()),
            &what,
        );
        near(union.volume_mm3().unwrap(), BOX + region.area * 5.);
        // No face is left at the contact: the box top keeps the rest of the
        // face, and nothing faces down there.
        let top = planar_faces_at(&union, SIZE[2], true);
        assert_eq!(top.len(), 1, "{what}: {top:?}");
        near(top[0], SIZE[0] * SIZE[1] - region.area);
        assert!(planar_faces_at(&union, SIZE[2], false).is_empty(), "{what}");
        near(planar_faces_at(&union, SIZE[2] + 5., true)[0], region.area);
        // One exact cylinder; coplanar box and region sides are one face.
        assert_eq!(curved_faces(&union), 1, "{what}");
        // Bottom, four sides, the rest of the top, the region top and the
        // cylinder: a region side on a box side adds no face.
        assert_eq!(union.faces().unwrap().len(), 8, "{what}");
        round_trips(&union);
        // The operands commute.
        let swapped = checked(
            boss.boolean(&stock(), BooleanOperation::Union, policy()),
            &what,
        );
        near(swapped.volume_mm3().unwrap(), BOX + region.area * 5.);
    }
}

#[test]
fn an_inward_join_leaves_the_box_unchanged() {
    for region in regions() {
        let what = format!("inward join {}", region.name);
        let inner = prism(top_frame(), region.outer.clone(), -5.);
        let union = checked(
            stock().boolean(&inner, BooleanOperation::Union, policy()),
            &what,
        );
        near(union.volume_mm3().unwrap(), BOX);
        assert_eq!(union.faces().unwrap().len(), 6, "{what}");
        round_trips(&union);
    }
}

#[test]
fn a_blind_cut_removes_the_region_prism_from_the_box() {
    for region in regions() {
        let what = format!("blind cut {}", region.name);
        let pocket = prism(top_frame(), region.outer.clone(), -5.);
        let cut = checked(
            stock().boolean(&pocket, BooleanOperation::Subtract, policy()),
            &what,
        );
        near(cut.volume_mm3().unwrap(), BOX - region.area * 5.);
        near(
            planar_faces_at(&cut, SIZE[2], true)[0],
            SIZE[0] * SIZE[1] - region.area,
        );
        // The pocket floor faces up 5 mm below the face.
        let floor = planar_faces_at(&cut, SIZE[2] - 5., true);
        assert_eq!(floor.len(), 1, "{what}");
        near(floor[0], region.area);
        assert_eq!(curved_faces(&cut), 1, "{what}");
        round_trips(&cut);
    }
}

#[test]
fn a_through_cut_removes_the_region_through_the_box() {
    for region in regions() {
        let what = format!("through cut {}", region.name);
        // From the top face to 5 mm past the bottom.
        let tool = prism(top_frame(), region.outer.clone(), -(SIZE[2] + 5.));
        let cut = checked(
            stock().boolean(&tool, BooleanOperation::Subtract, policy()),
            &what,
        );
        near(cut.volume_mm3().unwrap(), BOX - region.area * SIZE[2]);
        near(
            planar_faces_at(&cut, 0., false)[0],
            SIZE[0] * SIZE[1] - region.area,
        );
        assert_eq!(curved_faces(&cut), 1, "{what}");
        round_trips(&cut);
    }
}

#[test]
fn an_outward_cut_misses_the_box() {
    let region = &regions()[2];
    let tool = prism(top_frame(), region.outer.clone(), 5.);
    let cut = checked(
        stock().boolean(&tool, BooleanOperation::Subtract, policy()),
        "outward cut",
    );
    near(cut.volume_mm3().unwrap(), BOX);
}

#[test]
fn a_join_on_a_rotated_side_face_follows_that_face() {
    // The sketch on the +x side face: u = +y, v = +z, normal +x.
    let frame = PlaneFrame {
        origin: [SIZE[0], 0., 0.],
        u: [0., 1., 0.],
        v: [0., 0., 1.],
    };
    // Quarter disk about the face corner (y, z) = (30, 20).
    let outer = vec![
        CurvedEdge2d::Line {
            a: [30., 20.],
            b: [20., 20.],
        },
        CurvedEdge2d::Arc {
            center: [30., 20.],
            a: [20., 20.],
            b: [30., 10.],
            clockwise: false,
        },
        CurvedEdge2d::Line {
            a: [30., 10.],
            b: [30., 20.],
        },
    ];
    let area = PI * 100. / 4.;
    // The box placed away from the origin keeps the chart exact.
    let moved = RigidTransform {
        translation_mm: [3., -2., 7.],
        ..Default::default()
    };
    let block = stock().transformed(moved).unwrap();
    let frame = PlaneFrame {
        origin: [
            frame.origin[0] + 3.,
            frame.origin[1] - 2.,
            frame.origin[2] + 7.,
        ],
        ..frame
    };
    let union = checked(
        block.boolean(
            &prism(frame, outer.clone(), 4.),
            BooleanOperation::Union,
            policy(),
        ),
        "side join",
    );
    near(union.volume_mm3().unwrap(), BOX + area * 4.);
    assert_eq!(curved_faces(&union), 1);
    round_trips(&union);
    let cut = checked(
        block.boolean(
            &prism(frame, outer, -4.),
            BooleanOperation::Subtract,
            policy(),
        ),
        "side cut",
    );
    near(cut.volume_mm3().unwrap(), BOX - area * 4.);
    round_trips(&cut);
}

#[test]
fn a_tangent_region_is_refused_in_plain_words() {
    // A circle touching the face edge y = 30 from inside.
    let tool = prism(
        top_frame(),
        vec![CurvedEdge2d::Circle {
            center: [20., 22.],
            radius_mm: 8.,
        }],
        -5.,
    );
    let Err(error) = stock().boolean(&tool, BooleanOperation::Subtract, policy()) else {
        panic!("a tangent region is not certified");
    };
    assert_eq!(error.kind, GeometryErrorKind::Unsupported, "{error}");
    assert!(error.message.contains("tangentially"), "{error}");
}

/// A face circle that no edge clips extrudes as an exact placed cylinder.
#[test]
fn a_cylinder_from_the_face_joins_and_cuts_flush_with_it() {
    let disk = |z: f64, height: f64| {
        Solid::cylinder(
            PlaneFrame {
                origin: [20., 15., z],
                ..Default::default()
            },
            6.,
            height,
            policy(),
        )
        .unwrap()
    };
    let area = PI * 36.;
    let join = checked(
        stock().boolean(&disk(SIZE[2], 5.), BooleanOperation::Union, policy()),
        "cylinder join",
    );
    near(join.volume_mm3().unwrap(), BOX + area * 5.);
    assert!(planar_faces_at(&join, SIZE[2], false).is_empty());
    assert_eq!(curved_faces(&join), 1);
    round_trips(&join);
    // Blind from the face, then through the box.
    for (z, height, removed) in [(SIZE[2] - 5., 5., 5.), (-5., SIZE[2] + 5., SIZE[2])] {
        let cut = checked(
            stock().boolean(&disk(z, height), BooleanOperation::Subtract, policy()),
            "cylinder cut",
        );
        near(cut.volume_mm3().unwrap(), BOX - area * removed);
        assert_eq!(curved_faces(&cut), 1);
        round_trips(&cut);
    }
}

/// Cut target discovery measures the exact common volume of the box and
/// the region prism; the planar-only provider is not asked to certify it.
#[test]
fn the_common_volume_of_box_and_region_prism_is_exact() {
    for region in regions() {
        let inward = prism(top_frame(), region.outer.clone(), -5.);
        let volume = stock()
            .rust_occupied_overlap_volume(&inward, policy())
            .unwrap();
        near(volume, region.area * 5.);
        let common = checked(
            stock().boolean(&inward, BooleanOperation::Intersect, policy()),
            region.name,
        );
        near(common.volume_mm3().unwrap(), region.area * 5.);
        // Face contact alone removes nothing.
        let outward = prism(top_frame(), region.outer.clone(), 5.);
        let contact = stock()
            .rust_occupied_overlap_volume(&outward, policy())
            .unwrap();
        assert_eq!(contact, 0., "{}", region.name);
    }
}

#[test]
fn a_residual_floor_below_the_minimum_feature_is_refused_not_cut_through() {
    let region = &regions()[2];
    let tool = prism(top_frame(), region.outer.clone(), -(SIZE[2] - 1e-8));
    let Err(error) = stock().boolean(&tool, BooleanOperation::Subtract, policy()) else {
        panic!("a 1e-8 mm floor is neither kept nor cut through");
    };
    assert_eq!(error.kind, GeometryErrorKind::Unsupported, "{error}");
    assert!(error.message.contains("minimum feature"), "{error}");
}
