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

#[test]
fn positive_ulp_terminal_band_keeps_original_source_and_bounded_literal_mesh() {
    let length = 3.3000000000000003;
    let binary = Q::from_float(length).unwrap();
    let binary_end = Rational::new(
        binary.numer().to_string().parse().unwrap(),
        binary.denom().to_string().parse().unwrap(),
    )
    .unwrap();
    let decimal_end = Rational::new(33_000_000_000_000_003, 10_000_000_000_000_000).unwrap();
    let options = crate::display::Tessellation {
        deflection: 0.01,
        angle: 0.225,
        triangle_limit: 500_000,
    };
    for left in [false, true] {
        for explicit in [false, true] {
            let original = Cavity::new(thread(left, length), intent()).unwrap();
            assert_eq!(original.terminal_station_mm(), Some(3.0000000000000004));
            let end = if explicit { decimal_end } else { binary_end };
            let cavity = if explicit {
                Cavity::new_with_exact_receiving_extent(*original.source_thread(), intent(), end)
                    .unwrap()
            } else {
                original
            };
            let material = Material::new(stock(end), cavity, 1.).unwrap();
            let before = material;
            let calls = std::cell::Cell::new(0usize);
            let count = || {
                calls.set(calls.get() + 1);
                false
            };
            let candidate = material
                .tessellate_candidate(options, Some(&count))
                .unwrap();
            let complete_calls = calls.get();
            let proof = crate::manufacturing::certify(&candidate.mesh).unwrap();
            assert!(proof.signed_volume_mm3 > 0. && proof.vertex_links > 0);
            assert!(candidate.source_chord_bound_mm <= options.deflection);
            assert!(candidate.source_normal_angle_bound_rad <= options.angle);
            assert!(candidate.work_used <= 16_000_000);
            assert_eq!(material, before);
            assert_eq!(cavity.source_thread(), original.source_thread());
            assert_eq!(cavity.terminal_station_mm(), original.terminal_station_mm());
            assert_eq!(cavity.intent(), original.intent());
            assert_eq!(cavity.exact_receiving_extent().is_some(), explicit);
            assert_eq!(material.stock().height(), end);
            eprintln!(
                "ulpTerminal left={left} explicit={explicit} facets={} work={}",
                candidate.mesh.face_indices.len(),
                candidate.work_used
            );
            if !left && explicit {
                for stop in [0, complete_calls / 2, complete_calls - 1] {
                    calls.set(0);
                    let cancel = || {
                        let n = calls.get();
                        calls.set(n + 1);
                        n >= stop
                    };
                    assert_eq!(
                        material.tessellate_candidate(options, Some(&cancel)).err(),
                        Some(Error::Cancelled)
                    );
                    assert_eq!(material, before);
                }
                assert_eq!(
                    material
                        .tessellate_candidate(
                            crate::display::Tessellation {
                                triangle_limit: 16,
                                ..options
                            },
                            None
                        )
                        .err(),
                    Some(Error::Limit)
                );
                let retry = material.tessellate_candidate(options, None).unwrap();
                assert_eq!(retry.mesh, candidate.mesh);
                assert_eq!(retry.work_used, candidate.work_used);
            }
        }
    }
}

#[test]
fn tail_only_exact_domain_has_original_visible_law_real_cap_owners_and_independent_measures() {
    fn f(value: f64) -> Q {
        Q::from_float(value).unwrap()
    }
    fn pi_oracle(value: Enclosure, coefficient: &Q) {
        let lower = Q::from_float(std::f64::consts::PI.next_down()).unwrap() * coefficient;
        let upper = Q::from_float(std::f64::consts::PI.next_up()).unwrap() * coefficient;
        assert!(
            Q::from_float(value.lower).unwrap() <= lower
                && upper <= Q::from_float(value.upper).unwrap(),
            "{value:?}"
        );
    }
    let only_tail = Intent {
        entry_depth_mm: None,
        ..intent()
    };
    for left in [false, true] {
        for (visible, end) in [
            (4., Rational::new(4, 1).unwrap()),
            (3.3, Rational::new(33, 10).unwrap()),
            (
                3.3000000000000003,
                Rational::new(33_000_000_000_000_003, 10_000_000_000_000_000).unwrap(),
            ),
            (4., Rational::new(41, 10).unwrap()),
            (4., Rational::new(39, 10).unwrap()),
        ] {
            let original = Cavity::new(thread(left, visible), only_tail).unwrap();
            let body =
                Cavity::new_with_exact_receiving_extent(*original.source_thread(), only_tail, end)
                    .unwrap();
            assert_eq!(body.source_thread(), original.source_thread());
            assert_eq!(body.tail, original.tail);
            assert_eq!(body.intent(), original.intent());
            assert_eq!(body.entry_radius_mm(), None);
            let tail = body.tail.unwrap();
            let material = Material::new(stock(end), body, 1.).unwrap();
            // Exact independent integral of the ACTUAL five shared affine
            // bands, integrated over one true angular turn at each station.
            let mut mean = Q::from_integer(0.into());
            for band in body.source_thread().bands() {
                let width = f(band.phase_mm[1]) - f(band.phase_mm[0]);
                let a = f(band.radius_mm[0]);
                let b = f(band.radius_mm[1]);
                mean += width * (&a * &a + &a * &b + &b * &b) / Q::from_integer(3.into());
            }
            mean /= f(body.source_thread().spec().pitch_mm);
            let station = f(tail.station);
            let terminal_radius = f(tail.radius);
            let terminal_squared = &terminal_radius * &terminal_radius;
            let terminal_height = q(end) - &station;
            let cavity_coefficient = &mean * &station + &terminal_squared * &terminal_height;
            pi_oracle(
                body.removed_volume_mm3(Default::default()).unwrap(),
                &cavity_coefficient,
            );
            let radius = q(material.stock().circle().radius());
            let stock_squared = &radius * &radius;
            pi_oracle(
                material.volume_bounds(Default::default()).unwrap(),
                &(&stock_squared * q(end) - &cavity_coefficient),
            );
            pi_oracle(
                material
                    .local_face_properties(FaceKey::StartCap, Default::default())
                    .unwrap()
                    .area_mm2,
                &(&stock_squared - &mean),
            );
            let endcap = material
                .local_face_properties(FaceKey::EndCap, Default::default())
                .unwrap();
            pi_oracle(endcap.area_mm2, &(&stock_squared - &terminal_squared));
            assert!(
                f(endcap.centroid_mm[2].lower) <= q(end)
                    && q(end) <= f(endcap.centroid_mm[2].upper)
            );
            let terminal = material
                .local_face_properties(FaceKey::TerminalCylinder, Default::default())
                .unwrap();
            pi_oracle(
                terminal.area_mm2,
                &(Q::from_integer(2.into()) * &terminal_radius * &terminal_height),
            );
            let middle = (q(end) + &station) / Q::from_integer(2.into());
            assert!(
                f(terminal.centroid_mm[2].lower) <= middle
                    && middle <= f(terminal.centroid_mm[2].upper)
            );
            let boundary = material
                .boundary_with_interrupt(Default::default(), None)
                .unwrap();
            let mut cap_edges = 0;
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
                    cap_edges += 1;
                }
            }
            assert!(cap_edges >= 2);
            let mesh = material
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
            assert!(mesh.work_used <= 16_000_000);
            assert!(
                crate::manufacturing::certify(&mesh.mesh)
                    .unwrap()
                    .signed_volume_mm3
                    > 0.
            );
            if end == Rational::new(41, 10).unwrap() {
                assert!(body.source_thread().point_in_band(0, 0, 0., 4.05).is_err());
                assert!(body.contains_void([1., 0., 4.05], 0.));
            }
        }
    }
}
