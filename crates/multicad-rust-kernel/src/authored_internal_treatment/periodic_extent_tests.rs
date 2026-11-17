use super::*;
use crate::{
    Circle, Rational, authored_flat_thread::SourceOptions, cylinder::AxialCylinder,
    helical_thread::Spec,
};
use num_rational::BigRational as Q;
fn exact(q: Rational) -> Q {
    let (n, d) = q.ratio();
    Q::new(n.into(), d.into())
}
fn float(x: f64) -> Q {
    Q::from_float(x).unwrap()
}
fn source(left: bool, end: Rational) -> Cavity {
    let source = AuthoredThread::new_with_exact_axial_domain(
        Spec {
            major_radius_mm: 1.5,
            pitch_mm: 0.5,
            length_mm: 3.,
            depth_mm: 0.25,
            crest_width_mm: 0.125,
            lower_flank_deg: 30.,
            upper_flank_deg: 30.,
            starts: 1,
            left_hand: left,
            root_radius_mm: None,
        },
        end,
        SourceOptions {
            position_limit_mm: 1e-7,
            normal_angle_limit_rad: 1e-7,
        },
        None,
    )
    .unwrap();
    Cavity::new_with_exact_receiving_extent(
        source,
        Intent {
            entry_depth_mm: Some(0.4),
            terminal: None,
        },
        end,
    )
    .unwrap()
}
fn stock(end: Rational) -> AxialCylinder {
    let zero = Rational::new(0, 1).unwrap();
    AxialCylinder::new(
        7,
        Circle::new([zero; 2], Rational::new(7, 2).unwrap()).unwrap(),
        zero,
        end,
    )
    .unwrap()
}
fn enclosed(e: Enclosure, lo: Q, hi: Q) {
    assert!(float(e.lower) <= lo && hi <= float(e.upper), "{:?}", e);
}
#[test]
fn true_exact_entry_only_domain_has_real_end_graph_and_independent_max_column_mass() {
    let pi0 = float(std::f64::consts::PI.next_down());
    let pi1 = float(std::f64::consts::PI.next_up());
    for left in [false, true] {
        for end in [
            Rational::new(33, 10).unwrap(),
            Rational::new(16, 5).unwrap(),
            Rational::new(11, 4).unwrap(),
            Rational::new(33_000_000_000_000_003, 10_000_000_000_000_000).unwrap(),
        ] {
            let cavity = source(left, end);
            let material = Material::new(stock(end), cavity, 0.5).unwrap();
            assert_eq!(cavity.source_thread().spec().length_mm, 3.);
            assert_eq!(cavity.source_thread().exact_axial_extent(), Some(end));
            assert_eq!(cavity.exact_receiving_extent(), Some(end));
            assert_eq!(cavity.intent().terminal, None);
            let depth = float(cavity.intent().entry_depth_mm.unwrap());
            let a = float(cavity.entry_radius_mm().unwrap());
            let k = (&a - float(cavity.source_thread().minimum_radius_mm())) / depth;
            let mut column = Q::from_integer(0.into());
            let mut mean_squared = Q::from_integer(0.into());
            for b in cavity.source_thread().bands() {
                let x = float(b.radius_mm[0]);
                let y = float(b.radius_mm[1]);
                let dq = float(b.phase_mm[1]) - float(b.phase_mm[0]);
                let r2 = (&x * &x + &x * &y + &y * &y) / Q::from_integer(3.into());
                let r3 = (&x * &x * &x + &x * &x * &y + &x * &y * &y + &y * &y * &y)
                    / Q::from_integer(4.into());
                column += &dq
                    * (exact(end) * &r2
                        + (&a * &a * &a - Q::from_integer(3.into()) * &a * &r2
                            + Q::from_integer(2.into()) * r3)
                            / (Q::from_integer(3.into()) * &k));
                mean_squared += dq * r2;
            }
            let pitch = float(cavity.source_thread().spec().pitch_mm);
            column /= pitch.clone();
            mean_squared /= pitch;
            enclosed(
                cavity.removed_volume_mm3(Default::default()).unwrap(),
                &pi0 * &column,
                &pi1 * &column,
            );
            let radius = exact(stock(end).circle().radius());
            let material_factor = &radius * &radius * exact(end) - column;
            enclosed(
                material.volume_bounds(Default::default()).unwrap(),
                &pi0 * &material_factor,
                &pi1 * &material_factor,
            );
            let endcap = material
                .local_face_properties(FaceKey::EndCap, Default::default())
                .unwrap();
            let area = &radius * &radius - mean_squared;
            enclosed(endcap.area_mm2, &pi0 * &area, &pi1 * &area);
            enclosed(endcap.centroid_mm[2], exact(end), exact(end));
            let startcap = material
                .local_face_properties(FaceKey::StartCap, Default::default())
                .unwrap();
            let area = &radius * &radius - &a * &a;
            enclosed(startcap.area_mm2, &pi0 * &area, &pi1 * &area);
            let boundary = material
                .boundary_with_interrupt(Default::default(), None)
                .unwrap();
            let mut actual_ends = 0;
            for (id, curve) in boundary.curves() {
                if curve.faces().contains(&1) {
                    let coefficients = boundary
                        .exact_curve_coefficients(id, 1_000_000, 0, None)
                        .unwrap();
                    let z = coefficients.coordinate(0, 1).unwrap();
                    let zero = coefficients.coordinate(1, 1).unwrap();
                    let convert = |v: ExactBoundaryRational<'_>| {
                        use num_bigint::{BigInt, Sign};
                        Q::new(
                            BigInt::from_bytes_le(
                                if v.signum() < 0 {
                                    Sign::Minus
                                } else {
                                    Sign::Plus
                                },
                                v.numerator_magnitude_le(),
                            ),
                            BigInt::from_bytes_le(Sign::Plus, v.denominator_magnitude_le()),
                        )
                    };
                    assert_eq!(convert(z), exact(end));
                    assert_eq!(convert(zero), Q::from_integer(0.into()));
                    actual_ends += 1;
                }
            }
            assert!(actual_ends >= 6);
            let wrong = Rational::new(33, 10).unwrap();
            if end != wrong {
                assert!(Material::new(stock(wrong), cavity, 0.5).is_err());
                assert!(
                    Cavity::new_with_exact_receiving_extent(
                        *cavity.source_thread(),
                        cavity.intent(),
                        wrong
                    )
                    .is_err()
                );
            }
        }
    }
}
#[test]
fn canonical_entry_only_actual_mesh_retains_every_finite_source_and_closes() {
    for left in [false, true] {
        for end in [
            Rational::new(33, 10).unwrap(),
            Rational::new(16, 5).unwrap(),
        ] {
            let actual = Material::new(stock(end), source(left, end), 0.5).unwrap();
            let mesh = actual
                .tessellate_candidate(
                    crate::display::Tessellation {
                        deflection: 0.02,
                        angle: 0.3,
                        triangle_limit: 500_000,
                    },
                    None,
                )
                .unwrap();
            assert!(
                mesh.source_chord_bound_mm <= 0.02 && mesh.source_normal_angle_bound_rad <= 0.3
            );
            assert!(
                crate::manufacturing::certify(&mesh.mesh)
                    .unwrap()
                    .signed_volume_mm3
                    > 0.
            );
        }
    }
}

#[test]
fn exact_entry_source_refuses_same_float_different_end_and_too_short_actual_domain() {
    let canonical = Rational::new(16, 5).unwrap();
    let binary = Rational::new(3_602_879_701_896_397, 1_125_899_906_842_624).unwrap();
    assert_eq!(canonical.as_f64(), binary.as_f64());
    assert_ne!(canonical, binary);
    let cavity = source(false, canonical);
    assert!(Material::new(stock(binary), cavity, 0.5).is_err());
    assert!(
        Cavity::new_with_exact_receiving_extent(*cavity.source_thread(), cavity.intent(), binary)
            .is_err()
    );
    let original = AuthoredThread::new(
        cavity.source_thread().spec(),
        SourceOptions {
            position_limit_mm: 1e-7,
            normal_angle_limit_rad: 1e-7,
        },
    )
    .unwrap();
    assert!(Cavity::new_with_exact_receiving_extent(original, cavity.intent(), canonical).is_err());
    let short = AuthoredThread::new_with_exact_axial_domain(
        original.spec(),
        Rational::new(3, 5).unwrap(),
        SourceOptions {
            position_limit_mm: 1e-7,
            normal_angle_limit_rad: 1e-7,
        },
        None,
    )
    .unwrap();
    assert!(matches!(
        Cavity::new(short, cavity.intent()),
        Err(Error::Unresolved)
    ));
    assert!(Material::new(stock(canonical), cavity, 0.5).is_ok());
}
