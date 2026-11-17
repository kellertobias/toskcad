use super::*;
use crate::Rational;
use num_rational::BigRational as Q;
fn rational(q: Rational) -> Q {
    let (n, d) = q.ratio();
    Q::new(n.into(), d.into())
}
fn options() -> SourceOptions {
    SourceOptions {
        position_limit_mm: 1e-7,
        normal_angle_limit_rad: 1e-7,
    }
}
fn spec(left: bool, starts: u32) -> Spec {
    Spec {
        major_radius_mm: 1.5,
        pitch_mm: 0.5,
        length_mm: 3.,
        depth_mm: 0.25,
        crest_width_mm: 0.125,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts,
        left_hand: left,
        root_radius_mm: None,
    }
}
#[test]
fn exact_finite_domain_preserves_original_law_and_independent_radial_end_witnesses() {
    for left in [false, true] {
        for starts in [1, 3] {
            for end in [
                Rational::new(33, 10).unwrap(),
                Rational::new(16, 5).unwrap(),
                Rational::new(33_000_000_000_000_003, 10_000_000_000_000_000).unwrap(),
                Rational::new(11, 4).unwrap(),
            ] {
                let original = AuthoredThread::new(spec(left, starts), options()).unwrap();
                let source = AuthoredThread::new_with_exact_axial_domain(
                    original.spec(),
                    end,
                    options(),
                    None,
                )
                .unwrap();
                assert_eq!(source.spec(), original.spec());
                assert_eq!(source.bands(), original.bands());
                assert_eq!(source.exact_axial_extent(), Some(end));
                assert_eq!(original.exact_axial_extent(), None);
                let [lo, hi] = source.axial_extent_bounds_mm();
                let q = rational(end);
                assert!(Q::from_float(lo).unwrap() <= q && q <= Q::from_float(hi).unwrap());
                let (approx, transfer) = source.legacy_approximation();
                assert!(Q::from_float(approx.spec().length_mm).unwrap() >= q);
                assert!(
                    transfer.maximum_position_error_mm < 1e-7
                        && transfer.maximum_normal_angle_error_rad < 1e-7
                );
                let last = if Q::from_float(end.as_f64()).unwrap() <= q {
                    end.as_f64()
                } else {
                    end.as_f64().next_down()
                };
                assert!(source.point_in_band(0, 0, 0., last.next_up()).is_err());
                for z in [0., last / 2., last] {
                    for start in 0..starts {
                        for (band, b) in source.bands().into_iter().enumerate() {
                            for u in [0., 0.25, 0.5, 1.] {
                                let sample = source.point_in_band(band, start, u, z).unwrap();
                                assert_eq!(sample.point_mm[2], z);
                                let a = Q::from_float(b.radius_mm[0]).unwrap();
                                let dr = Q::from_float(b.radius_mm[1]).unwrap() - &a;
                                let radius = a + dr * Q::from_float(u).unwrap();
                                let norm = sample.point_mm[0].hypot(sample.point_mm[1]);
                                // Independent exact affine radius, with outward sample/trig
                                // rounding enclosure rather than another source evaluator.
                                let low = (norm - sample.error_mm).next_down();
                                let high = (norm + sample.error_mm).next_up();
                                assert!(
                                    Q::from_float(low).unwrap() <= radius
                                        && radius <= Q::from_float(high).unwrap()
                                );
                                if band < 4 && u == 1. {
                                    assert_eq!(
                                        sample.point_mm,
                                        source
                                            .point_in_band(band + 1, start, 0., z)
                                            .unwrap()
                                            .point_mm
                                    );
                                }
                            }
                        }
                    }
                }
                if q > Q::from_integer(3.into()) {
                    assert!(original.point_in_band(0, 0, 0., last).is_err());
                    assert!(source.point_in_band(0, 0, 0., last).is_ok());
                }
                assert_eq!(original.spec().length_mm, 3.);
            }
        }
    }
}
#[test]
fn exact_domain_rejects_real_invalid_span_turn_budget_and_cancellation_without_changing_legacy() {
    let original = AuthoredThread::new(spec(false, 1), options()).unwrap();
    for end in [
        Rational::new(-1, 1).unwrap(),
        Rational::new(0, 1).unwrap(),
        Rational::new(1, 10_000).unwrap(),
        Rational::new(10_000_001, 1).unwrap(),
    ] {
        assert_eq!(
            AuthoredThread::new_with_exact_axial_domain(original.spec(), end, options(), None),
            Err(Error::Invalid)
        );
    }
    assert_eq!(
        AuthoredThread::new_with_exact_axial_domain(
            original.spec(),
            Rational::new(100, 1).unwrap(),
            options(),
            None
        ),
        Err(Error::Limit)
    );
    let end = Rational::new(33, 10).unwrap();
    let expected =
        AuthoredThread::new_with_exact_axial_domain(original.spec(), end, options(), None).unwrap();
    let calls = std::cell::Cell::new(0usize);
    let count = || {
        calls.set(calls.get() + 1);
        false
    };
    assert_eq!(
        AuthoredThread::new_with_exact_axial_domain(original.spec(), end, options(), Some(&count))
            .unwrap(),
        expected
    );
    let total = calls.get();
    for stop in [1, total / 2, total] {
        calls.set(0);
        let cancel = || {
            calls.set(calls.get() + 1);
            calls.get() >= stop
        };
        assert_eq!(
            AuthoredThread::new_with_exact_axial_domain(
                original.spec(),
                end,
                options(),
                Some(&cancel)
            ),
            Err(Error::Cancelled)
        );
    }
    assert_eq!(
        AuthoredThread::new(original.spec(), options()).unwrap(),
        original
    );
    assert!(original.point_in_band(0, 0, 0., 3.).is_ok());
    assert!(original.point_in_band(0, 0, 0., 3_f64.next_up()).is_err());
    assert_eq!(
        AuthoredThread::new_with_exact_axial_domain(
            Spec {
                length_mm: -1.,
                ..original.spec()
            },
            end,
            options(),
            None
        ),
        Err(Error::Invalid)
    );
    assert_eq!(
        AuthoredThread::new_with_exact_axial_domain(
            Spec {
                root_radius_mm: Some(0.01),
                ..original.spec()
            },
            end,
            options(),
            None
        ),
        Err(Error::Unsupported)
    );
    assert_eq!(
        AuthoredThread::new_with_exact_axial_domain(
            original.spec(),
            end,
            SourceOptions {
                position_limit_mm: 0.,
                ..options()
            },
            None
        ),
        Err(Error::Invalid)
    );
    for z in [f64::NAN, f64::INFINITY, -0.001] {
        assert!(expected.point_in_band(0, 0, 0., z).is_err());
    }
}

#[test]
fn public_bare_exact_source_graph_shell_and_planar_query_keep_the_same_endpoint() {
    use crate::authored_internal_treatment::{ExactBoundaryRational, FaceKey};
    use num_bigint::{BigInt, Sign};
    let q = crate::Rational::new(16, 5).unwrap();
    for left in [false, true] {
        let source =
            AuthoredThread::new_with_exact_axial_domain(spec(left, 1), q, options(), None).unwrap();
        let boundary = source
            .boundary_with_interrupt(Default::default(), None)
            .unwrap();
        assert_eq!(boundary.source(), &source);
        let mut end_edges = 0;
        for (id, curve) in boundary.source_curves() {
            if curve.faces().contains(&1) {
                let coefficients = boundary
                    .exact_curve_coefficients(id, 1_000_000, 0, None)
                    .unwrap();
                let to_q = |v: ExactBoundaryRational<'_>| {
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
                assert_eq!(to_q(coefficients.coordinate(0, 1).unwrap()), rational(q));
                end_edges += 1;
            }
        }
        // z=16/5 cuts the phase period at1/5, strictly inside a flank.
        // That one physical band straddles the angular seam: six rim uses.
        assert_eq!(end_edges, 6);
        let cap = source
            .local_face_properties(FaceKey::EndCap, Default::default(), None)
            .unwrap();
        assert!(
            Q::from_float(cap.centroid_mm[2].lower).unwrap() <= rational(q)
                && rational(q) <= Q::from_float(cap.centroid_mm[2].upper).unwrap()
        );
        let quality = crate::display::Tessellation {
            deflection: 0.02,
            angle: 0.3,
            triangle_limit: 500_000,
        };
        let shell = source.tessellate_shell(quality, 16_000_000, None).unwrap();
        assert_eq!(shell.source(), &source);
        assert!(shell.rim_count(ShellEnd::End) >= 5);
        let [lo, hi] = source.axial_extent_bounds_mm();
        for edge in 0..shell.rim_count(ShellEnd::End) {
            for vertex in shell.rim_edge(ShellEnd::End, edge).unwrap().vertex_indices {
                let z = shell.mesh().vertices[3 * vertex + 2];
                assert!(lo <= z && z <= hi);
            }
        }
        assert!(
            source
                .tessellate_prefix_shell(q.as_f64(), quality, 16_000_000, None)
                .is_err()
        );
        assert!(
            source
                .local_prefix_face_properties(q.as_f64(), FaceKey::EndCap, Default::default(), None)
                .is_err()
        );
        assert!(
            source
                .tessellate_prefix_shell(3., quality, 16_000_000, None)
                .is_ok()
        );
    }
}
