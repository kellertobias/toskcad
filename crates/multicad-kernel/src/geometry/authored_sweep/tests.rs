use super::*;
use std::f64::consts::PI;

fn policy() -> GeometryPolicy {
    GeometryPolicy::default()
}

fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-6 * b.abs().max(1.), "{a} != {b}");
}

fn line(a: P2, b: P2) -> CurvedEdge2d {
    CurvedEdge2d::Line { a, b }
}

fn names(profile: &CurvedProfile) -> Vec<Origin> {
    edges(profile)
        .map(|_| Origin {
            source: EntityId::new(),
            role: "curve".into(),
        })
        .collect()
}

/// A half disc of radius 3 bounded by a clockwise arc, then its diameter.
fn half_disc() -> CurvedProfile {
    CurvedProfile {
        frame: PlaneFrame::default(),
        outer: vec![
            CurvedEdge2d::Arc {
                center: [0., 0.],
                a: [3., 0.],
                b: [-3., 0.],
                clockwise: true,
            },
            line([-3., 0.], [3., 0.]),
        ],
        holes: vec![],
    }
}

/// An L-shaped section whose outer corner is rounded.
fn rounded_l() -> CurvedProfile {
    CurvedProfile {
        frame: PlaneFrame::default(),
        outer: vec![
            line([0., 0.], [8., 0.]),
            CurvedEdge2d::Arc {
                center: [8., 2.],
                a: [8., 0.],
                b: [10., 2.],
                clockwise: false,
            },
            line([10., 2.], [10., 4.]),
            line([10., 4.], [4., 4.]),
            line([4., 4.], [4., 10.]),
            line([4., 10.], [0., 10.]),
            line([0., 10.], [0., 0.]),
        ],
        holes: vec![],
    }
}

fn translated(profile: &CurvedProfile, displacement: P) -> Result<Solid, GeometryError> {
    let feature = EntityId::new();
    Solid::rust_authored_translation(profile, displacement, feature, &names(profile), policy())
}

#[test]
fn clockwise_arc_translation_keeps_area_times_normal_height_and_reversed_walls() {
    // The clockwise arc runs below the diameter; the section is a half disc.
    for displacement in [[2., -1., 5.], [-2., 1., -5.]] {
        let solid = translated(&half_disc(), displacement).unwrap();
        near(solid.volume_mm3().unwrap(), 4.5 * PI * 5.);
        let bounds = solid.bounds().unwrap();
        let low = if displacement[2] > 0. { 0. } else { -5. };
        near(bounds.min[2], low);
        near(bounds.min[1], -3. + displacement[1].min(0.));
        assert_eq!(solid.solid_count().unwrap(), 1);
    }
}

#[test]
fn authored_names_follow_rigid_placements() {
    let feature = EntityId::new();
    let profile = half_disc();
    let origins = names(&profile);
    let solid =
        Solid::rust_authored_translation(&profile, [1., 0., 4.], feature, &origins, policy())
            .unwrap();
    let transform = RigidTransform {
        rotation: [[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]],
        translation_mm: [10., 0., 1.],
    };
    let moved = solid.transformed(transform).unwrap();
    let cap = |solid: &Solid| {
        solid
            .topology_snapshot(feature, 0)
            .unwrap()
            .elements
            .into_iter()
            .find(|e| {
                e.origins.contains(&Origin {
                    source: feature,
                    role: "cap.end".into(),
                })
            })
            .unwrap()
    };
    let (before, after) = (cap(&solid), cap(&moved));
    assert_eq!(before.origins, after.origins);
    assert_eq!(before.frames[0].frame.origin_mm, [1., 0., 4.]);
    let origin = after.frames[0].frame.origin_mm;
    for (axis, value) in [10., 1., 5.].into_iter().enumerate() {
        near(origin[axis], value);
    }
    near(after.frames[0].frame.u[1], 1.);
    near(moved.volume_mm3().unwrap(), solid.volume_mm3().unwrap());
}

#[test]
fn plane_sweeps_refuse_invalid_directions_planes_and_reach() {
    let profile = half_disc();
    let feature = EntityId::new();
    let origins = names(&profile);
    let plane = |height: f64| PlaneFrame {
        origin: [0., 0., height],
        ..PlaneFrame::default()
    };
    let sweep = |direction: P, start: Option<PlaneFrame>, end: PlaneFrame, from_sketch: bool| {
        Solid::rust_authored_plane_sweep(
            &profile,
            direction,
            start,
            end,
            from_sketch,
            feature,
            &origins,
            policy(),
        )
    };
    let messages = [
        sweep([0., 0., 2.], None, plane(5.), true),
        sweep([f64::NAN, 0., 1.], None, plane(5.), true),
        sweep([1., 0., 0.], None, plane(5.), true),
        sweep(
            [0., 0., 1.],
            None,
            PlaneFrame {
                origin: [0., 0., 5.],
                u: [1., 0., 0.],
                v: [0., 0., 1.],
            },
            true,
        ),
        sweep([0., 0., 1.], None, plane(-5.), true),
        sweep([0., 0., 1.], Some(plane(1.)), plane(5.), true),
        sweep([0., 0., 1.], Some(plane(6.)), plane(5.), false),
        sweep([0., 0., 1.], None, plane(2e9), true),
    ]
    .map(|result| result.err().unwrap().message);
    assert!(messages[0].contains("unit length"));
    assert!(messages[1].contains("unit length"));
    assert!(messages[2].contains("cross the sketch plane"));
    assert!(messages[3].contains("termination plane"));
    assert!(messages[4].contains("ahead"));
    assert!(messages[5].contains("behind"));
    assert!(messages[6].contains("insufficient material"));
    assert!(!messages[7].is_empty());

    // Between two tilted planes on either side of the section.
    let tilt = |height: f64| {
        let l = 1_f64.hypot(0.25);
        PlaneFrame {
            origin: [0., 0., height],
            u: [1. / l, 0., 0.25 / l],
            v: [0., 1., 0.],
        }
    };
    let between = sweep([0., 0., 1.], Some(tilt(-2.)), tilt(4.), false).unwrap();
    near(between.volume_mm3().unwrap(), 4.5 * PI * 6.);
    let behind = sweep([0., 0., 1.], Some(plane(-2.)), plane(4.), true).unwrap();
    near(behind.volume_mm3().unwrap(), 4.5 * PI * 6.);
}

#[test]
fn tapered_halves_refuse_invalid_controls_and_join_drafted_frusta() {
    let profile = CurvedProfile {
        frame: PlaneFrame::default(),
        outer: vec![CurvedEdge2d::Circle {
            center: [0., 0.],
            radius_mm: 3.,
        }],
        holes: vec![],
    };
    let feature = EntityId::new();
    let origins = names(&profile);
    let halves = |direction: P, lengths: [f64; 2]| {
        Solid::rust_extrude_tapered_halves_named(
            &profile,
            direction,
            lengths,
            0.1_f64.atan(),
            feature,
            &origins,
            policy(),
        )
    };
    assert!(halves([0., 0., 2.], [4., 4.]).is_err());
    assert!(halves([0., 0., f64::INFINITY], [4., 4.]).is_err());
    assert!(halves([0., 0., 1.], [0., 4.]).is_err());
    assert!(halves([0., 0., 1.], [4., -1.]).is_err());
    let joined = halves([0., 0., 1.], [10., 5.]).unwrap();
    near(joined.volume_mm3().unwrap(), 176.25 * PI);
    assert_eq!(joined.solid_count().unwrap(), 1);
}

#[test]
fn rulings_refuse_unsupported_or_mismatched_sections() {
    let arc = half_disc();
    let mut end = arc.clone();
    end.frame.origin = [0., 0., 5.];
    let feature = EntityId::new();
    let origins = names(&arc);
    let refused = Solid::rust_authored_taper(&arc, &end, [0., 0., 5.], feature, &origins, policy());
    assert!(
        refused
            .err()
            .unwrap()
            .message
            .contains("lines and full circles")
    );
    // Every section curve needs its own authored origin.
    let refused =
        Solid::rust_authored_taper(&arc, &end, [0., 0., 5.], feature, &origins[..1], policy());
    assert!(refused.is_err());
}

#[test]
fn crossing_sections_and_bad_names_refuse_before_any_boundary() {
    let crossing = CurvedProfile {
        frame: PlaneFrame::default(),
        outer: vec![
            line([0., 0.], [4., 4.]),
            CurvedEdge2d::Arc {
                center: [2., 4.],
                a: [4., 4.],
                b: [0., 4.],
                clockwise: false,
            },
            line([0., 4.], [4., 0.]),
            line([4., 0.], [0., 0.]),
        ],
        holes: vec![],
    };
    assert!(translated(&crossing, [1., 0., 3.]).is_err());
    let profile = half_disc();
    let feature = EntityId::new();
    let mut origins = names(&profile);
    origins[1] = origins[0].clone();
    assert!(
        Solid::rust_authored_translation(&profile, [1., 0., 3.], feature, &origins, policy())
            .is_err()
    );
    assert!(translated(&profile, [1., 0., 1e-9]).is_err());
    assert!(translated(&profile, [1e12, 0., 1.]).is_err());
}

#[test]
fn glued_unions_refuse_other_operations_and_unshared_or_overlapping_operands() {
    let circle = CurvedProfile {
        frame: PlaneFrame::default(),
        outer: vec![CurvedEdge2d::Circle {
            center: [0., 0.],
            radius_mm: 3.,
        }],
        holes: vec![],
    };
    let up = translated(&circle, [1., 0., 4.]).unwrap();
    let union = BooleanOperation::Union;
    assert!(
        up.glued_brep_union(&up, BooleanOperation::Subtract, policy())
            .is_none()
    );
    let other = Solid::box_solid([1., 1., 1.], policy()).unwrap();
    assert!(up.glued_brep_union(&other, union, policy()).is_none());
    // Both operands start on the sketch plane facing the same way.
    let again = translated(&circle, [-1., 0., 4.]).unwrap();
    let refused = up.glued_brep_union(&again, union, policy()).unwrap();
    assert!(refused.err().unwrap().message.contains("exactly one"));
    let mut loose = policy();
    loose.linear_tolerance_mm = 1e-6;
    let refused = up.glued_brep_union(&again, union, loose).unwrap();
    assert!(refused.err().unwrap().message.contains("geometry policy"));
    // The far cap of a downward sweep and the start of an upward sweep
    // share the sketch section from opposite sides.
    let down = translated(&circle, [1., 0., -4.]).unwrap();
    let joined = up
        .glued_brep_union(&down, union, policy())
        .unwrap()
        .unwrap()
        .unwrap();
    near(joined.volume_mm3().unwrap(), 9. * PI * 8.);
    assert_eq!(joined.faces().unwrap().len(), 4);
}

#[test]
fn half_spaces_reject_a_face_that_does_not_bound_its_solid() {
    let solid = translated(&rounded_l(), [0., 0., 3.]).unwrap();
    let body = solid.imported_brep_body().unwrap();
    let boundary = body.solid.boundary();
    // The cap bounds the solid; the inner wall of the L (x = 4) does not.
    assert!(glue::behind(&body.solid, 0, 1e-7).unwrap());
    let inner = (0..boundary.faces.len())
        .find(|&face| {
            glue::plane(&body.solid, face)
                .is_some_and(|(origin, normal)| normal[0] > 0.5 && (origin[0] - 4.).abs() < 1e-9)
        })
        .unwrap();
    assert!(!glue::behind(&body.solid, inner, 1e-7).unwrap());
    near(solid.volume_mm3().unwrap(), (64. - 4. + PI) * 3.);
}

#[test]
fn boundary_failures_keep_their_kind_and_names_cover_every_face() {
    for (cause, kind) in [
        (brep::Error::Invalid, GeometryErrorKind::InvalidResult),
        (brep::Error::Unsupported, GeometryErrorKind::Unsupported),
        (brep::Error::Limit, GeometryErrorKind::Limit),
        (brep::Error::Cancelled, GeometryErrorKind::Cancelled),
    ] {
        assert_eq!(failure(cause).kind, kind);
    }
    let solid = translated(&half_disc(), [0., 1., 2.]).unwrap();
    let retained = solid.imported_brep_body().unwrap().solid.clone();
    let refused = Solid::from_authored_brep(
        retained,
        EntityId::new(),
        AuthoredNaming { faces: vec![] },
        policy(),
    );
    assert_eq!(
        refused.err().unwrap().kind,
        GeometryErrorKind::InvalidResult
    );
    assert!(unit([0.; 3]).is_err());
}
