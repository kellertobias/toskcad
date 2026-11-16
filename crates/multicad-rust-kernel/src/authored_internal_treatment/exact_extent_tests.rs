use super::*;
use crate::{
    Circle, Rational, authored_flat_thread::SourceOptions, cylinder::AxialCylinder,
    helical_thread::Spec,
};
use num_bigint::{BigInt, Sign};
use num_rational::BigRational as Q;
fn thread(left: bool, length: f64) -> AuthoredThread {
    AuthoredThread::new(
        Spec {
            major_radius_mm: 1.5,
            pitch_mm: 0.5,
            length_mm: length,
            depth_mm: 0.25,
            crest_width_mm: 0.125,
            lower_flank_deg: 30.,
            upper_flank_deg: 30.,
            starts: 1,
            left_hand: left,
            root_radius_mm: None,
        },
        SourceOptions {
            position_limit_mm: 1e-7,
            normal_angle_limit_rad: 1e-7,
        },
    )
    .unwrap()
}
fn intent() -> Intent {
    Intent {
        entry_depth_mm: Some(0.4),
        terminal: Some(Terminal {
            width_mm: 0.3,
            radial_depth_mm: 0.05,
        }),
    }
}
fn q(r: Rational) -> Q {
    let (n, d) = r.ratio();
    Q::new(n.into(), d.into())
}
fn exact(x: ExactBoundaryRational<'_>) -> Q {
    let sign = if x.signum() < 0 {
        Sign::Minus
    } else {
        Sign::Plus
    };
    Q::new(
        BigInt::from_bytes_le(sign, x.numerator_magnitude_le()),
        BigInt::from_bytes_le(Sign::Plus, x.denominator_magnitude_le()),
    )
}
fn stock(end: Rational) -> AxialCylinder {
    AxialCylinder::new(
        7,
        Circle::new(
            [Rational::new(0, 1).unwrap(); 2],
            Rational::new(7, 2).unwrap(),
        )
        .unwrap(),
        Rational::new(0, 1).unwrap(),
        end,
    )
    .unwrap()
}
#[test]
fn explicit_canonical_receiving_domain_keeps_original_visible_law_and_exact_shared_cap() {
    for left in [false, true] {
        for end in [
            Rational::new(33, 10).unwrap(),
            Rational::new(33000000000000003, 10000000000000000).unwrap(),
        ] {
            let original = Cavity::new(thread(left, 3.3000000000000003), intent()).unwrap();
            let body =
                Cavity::new_with_exact_receiving_extent(*original.source_thread(), intent(), end)
                    .unwrap();
            assert_eq!(body.source_thread(), original.source_thread());
            assert_eq!(body.entry, original.entry);
            assert_eq!(body.tail, original.tail);
            assert_eq!(body.intent(), original.intent());
            assert_eq!(body.exact_receiving_extent(), Some(end));
            for z in [0.01, 0.2, 0.7, 1.23, 2.7, 3.1] {
                for theta in [-2., 0., 0.43, 2.1] {
                    assert_eq!(body.radius_mm(theta, z), original.radius_mm(theta, z));
                }
            }
            let material = Material::new(stock(end), body, 1.).unwrap();
            assert_eq!(*material.stock(), stock(end));
            let boundary = material
                .boundary_with_interrupt(Default::default(), None)
                .unwrap();
            let mut terminal = 0;
            let mut outer = 0;
            for (id, curve) in boundary.curves() {
                if curve.faces().contains(&1) {
                    let data = boundary
                        .exact_curve_coefficients(id, 1_000_000, 0, None)
                        .unwrap();
                    assert_eq!(exact(data.coordinate(0, 1).unwrap()), q(end));
                    assert_eq!(
                        exact(data.coordinate(1, 1).unwrap()),
                        Q::from_integer(0.into())
                    );
                    if curve.faces().contains(&(boundary.face_count() - 1)) {
                        outer += 1
                    } else {
                        terminal += 1
                    }
                }
            }
            assert_eq!(outer, 1);
            assert!(terminal >= 3);
            let face = material
                .local_face_properties_with_interrupt(FaceKey::EndCap, Default::default(), None)
                .unwrap();
            assert!(
                Q::from_float(face.centroid_mm[2].lower).unwrap() <= q(end)
                    && q(end) <= Q::from_float(face.centroid_mm[2].upper).unwrap()
            );
            let cylinder = material
                .local_face_properties_with_interrupt(
                    FaceKey::TerminalCylinder,
                    Default::default(),
                    None,
                )
                .unwrap();
            let height = q(end) - Q::from_float(body.tail.unwrap().station).unwrap();
            let radius = Q::from_float(body.tail.unwrap().radius).unwrap();
            let lower = Q::from_float(std::f64::consts::PI.next_down()).unwrap()
                * Q::from_integer(2.into())
                * &radius
                * &height;
            let upper = Q::from_float(std::f64::consts::PI.next_up()).unwrap()
                * Q::from_integer(2.into())
                * &radius
                * &height;
            assert!(
                Q::from_float(cylinder.area_mm2.lower).unwrap() <= lower
                    && upper <= Q::from_float(cylinder.area_mm2.upper).unwrap()
            );
            let volume = material.volume_bounds(Default::default()).unwrap();
            assert!(volume.lower > 0. && volume.upper - volume.lower <= 1e-7);
            let mesh = material
                .tessellate_candidate(
                    crate::display::Tessellation {
                        deflection: 0.01,
                        angle: 0.3,
                        triangle_limit: 200_000,
                    },
                    None,
                )
                .unwrap();
            assert!(mesh.mesh.face_indices.len() > 100);
            assert!(mesh.source_chord_bound_mm <= 0.01);
        }
    }
}
#[test]
fn exact_end_stock_mismatch_even_with_same_f64_and_strict_legacy_refusal_remain() {
    let end = Rational::new(33, 10).unwrap();
    let thread = thread(false, 3.3);
    let legacy = Cavity::new(thread, intent()).unwrap();
    assert_eq!(
        Material::new(stock(end), legacy, 1.),
        Err(Error::Unresolved)
    );
    let body = Cavity::new_with_exact_receiving_extent(thread, intent(), end).unwrap();
    let other = Rational::new(330000000000000001, 100000000000000000).unwrap();
    assert_eq!(other.as_f64(), end.as_f64());
    assert_eq!(
        Material::new(stock(other), body, 1.),
        Err(Error::Unresolved)
    );
    for end in [
        Rational::new(0, 1).unwrap(),
        Rational::new(-1, 1).unwrap(),
        Rational::new(3, 1).unwrap(),
    ] {
        assert!(Cavity::new_with_exact_receiving_extent(thread, intent(), end).is_err());
    }
    assert!(
        Cavity::new_with_exact_receiving_extent(
            thread,
            Intent {
                terminal: None,
                ..intent()
            },
            end
        )
        .is_err()
    );
    let material = Material::new(stock(end), body, 1.).unwrap();
    assert_eq!(
        material
            .boundary_with_interrupt(Default::default(), Some(&|| true))
            .err(),
        Some(Error::Cancelled)
    );
    assert!(
        material
            .boundary_with_interrupt(Default::default(), None)
            .is_ok()
    );
}

#[test]
fn rounded_cap_representative_does_not_expand_the_exact_void_domain() {
    let end = Rational::new(16, 5).unwrap();
    let c = Cavity::new_with_exact_receiving_extent(thread(false, 3.2), intent(), end).unwrap();
    assert!(Q::from_float(c.full_length_mm()).unwrap() > q(end));
    assert!(!c.contains_void([0., 0., c.full_length_mm()], 0.));
    assert!(c.contains_void([0., 0., c.full_length_mm().next_down()], 0.));
    assert!(c.contains_void([0., 0., c.full_length_mm()], 1e-7));
}
