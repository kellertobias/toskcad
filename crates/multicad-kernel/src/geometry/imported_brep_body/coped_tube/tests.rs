use super::*;
use crate::{CurvedEdge2d, CurvedProfile};

const R: f64 = 33.7 / 2.;
const WALL: f64 = 3.2;

fn policy() -> GeometryPolicy {
    GeometryPolicy::default()
}

fn near(actual: f64, expected: f64) {
    assert!((actual - expected).abs() <= 1e-6, "{actual} != {expected}");
}

fn names(feature: EntityId) -> CopedTubeNames {
    let origin = |role: &str| Origin {
        source: feature,
        role: role.into(),
    };
    CopedTubeNames {
        outer: origin("member.outer.0"),
        inner: origin("member.inner.0"),
        end: origin("member.end"),
        saddle: origin("connection.cope"),
    }
}

fn branch(angle_deg: f64, grow: f64, length: f64) -> Result<Solid, GeometryError> {
    let feature = EntityId::new();
    let (s, c) = angle_deg.to_radians().portable_sin_cos();
    Solid::rust_coped_round_tube(
        &CopedTube {
            origin_mm: [300., 0., 0.],
            axis: [c, 0., s],
            outer_radius_mm: R,
            inner_radius_mm: R - WALL,
            length_mm: length,
            envelope_axis: [1., 0., 0.],
            envelope_radius_mm: R + grow,
        },
        &names(feature),
        feature,
        policy(),
    )
}

/// The receiving round tube along x through the origin.
fn main_tube(outer: f64) -> Solid {
    let feature = EntityId::new();
    let profile = CurvedProfile {
        frame: PlaneFrame {
            origin: [0., 0., 0.],
            u: [0., 1., 0.],
            v: [0., 0., 1.],
        },
        outer: vec![CurvedEdge2d::Circle {
            center: [0., 0.],
            radius_mm: outer,
        }],
        holes: vec![vec![CurvedEdge2d::Circle {
            center: [0., 0.],
            radius_mm: outer - WALL,
        }]],
    };
    let origins = ["outer", "inner"].map(|role| Origin {
        source: feature,
        role: role.into(),
    });
    Solid::rust_extrude_linear_named(&profile, [600., 0., 0.], 0., feature, &origins, policy())
        .unwrap()
}

#[test]
fn a_square_saddle_with_clearance_reaches_down_the_envelope_sides() {
    let coped = branch(90., 1., 400.).unwrap();
    let bounds = coped.bounds().unwrap();
    near(bounds.min[2], ((R + 1.).powi(2) - R * R).sqrt());
    near(bounds.max[2], 400.);
    near(bounds.min[0], 300. - R);
    near(bounds.max[1], R);
    // Every generator leaves the envelope between the side height and the
    // full envelope radius above its axis.
    let area = PI * (R * R - (R - WALL).powi(2));
    let volume = coped.volume_mm3().unwrap();
    assert!(volume < area * (400. - (2. * R + 1.).sqrt()));
    assert!(volume > area * (400. - (R + 1.)));
    assert_eq!(coped.solid_count().unwrap(), 1);
    assert_eq!(coped.faces().unwrap().len(), 4);
}

#[test]
fn an_equal_angled_branch_meets_the_receiving_wall_in_two_half_ellipses() {
    let coped = branch(60., 0., 300.).unwrap();
    let ellipses = coped
        .topology_snapshot(EntityId::new(), 0)
        .unwrap()
        .elements
        .into_iter()
        .filter(|element| matches!(element.geometry, Some(AnalyticGeometry::EllipticArc { .. })))
        .count();
    assert_eq!(ellipses, 2);
    let plain = PI * (R * R - (R - WALL).powi(2)) * 300.;
    assert!(coped.volume_mm3().unwrap() < plain);
}

#[test]
fn the_coped_branch_only_touches_the_member_it_is_fitted_to() {
    let main = main_tube(R);
    assert!(main.annulus_body().is_some());
    for (angle, grow) in [(90., 1.), (60., 0.)] {
        let coped = branch(angle, grow, 300.).unwrap();
        for (a, b) in [(&coped, &main), (&main, &coped)] {
            let overlap = a.boolean(b, BooleanOperation::Intersect, policy()).unwrap();
            assert!(overlap.is_none());
        }
    }
    // A receiving tube wider than the envelope really overlaps; the
    // contact proof does not claim it.
    let coped = branch(90., 1., 300.).unwrap();
    let wide = main_tube(R + 2.);
    assert!(
        coped
            .cope_contact_intersection(&wide, BooleanOperation::Intersect, policy())
            .is_none()
    );
    assert!(
        coped
            .cope_contact_intersection(&main, BooleanOperation::Subtract, policy())
            .is_none()
    );
}

#[test]
fn wide_short_or_parallel_copes_are_refused() {
    let feature = EntityId::new();
    let base = CopedTube {
        origin_mm: [0.; 3],
        axis: [0., 0., 1.],
        outer_radius_mm: R,
        inner_radius_mm: R - WALL,
        length_mm: 100.,
        envelope_axis: [1., 0., 0.],
        envelope_radius_mm: R,
    };
    let refused = |tube: CopedTube| {
        Solid::rust_coped_round_tube(&tube, &names(feature), feature, policy())
            .err()
            .unwrap()
            .kind
    };
    assert_eq!(
        refused(CopedTube {
            envelope_radius_mm: R - 1.,
            ..base
        }),
        GeometryErrorKind::Unsupported
    );
    assert_eq!(
        refused(CopedTube {
            length_mm: R,
            ..base
        }),
        GeometryErrorKind::Unsupported
    );
    assert_eq!(
        refused(CopedTube {
            envelope_axis: [0.01, 0., 1.],
            ..base
        }),
        GeometryErrorKind::Unsupported
    );
    assert!(Solid::rust_coped_round_tube(&base, &names(feature), feature, policy()).is_ok());
}
