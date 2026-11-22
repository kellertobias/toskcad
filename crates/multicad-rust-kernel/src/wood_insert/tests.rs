use super::*;
use crate::{authored_flat_thread::SourceOptions, helical_thread::Spec};
use std::cell::Cell;
fn source(outer: bool, left: bool, starts: u32, nonbinary: bool) -> AuthoredThread {
    let pitch_mid = (2.675 + 2.775) / 2.;
    let minor_mid = (2.459 + 2.599) / 2.;
    let spec = if outer {
        Spec {
            major_radius_mm: 5.,
            pitch_mm: 2.,
            length_mm: if nonbinary { 4.3 } else { 4. },
            depth_mm: 1.,
            crest_width_mm: 0.3,
            lower_flank_deg: 20.,
            upper_flank_deg: 30.,
            starts,
            left_hand: left,
            root_radius_mm: None,
        }
    } else {
        Spec {
            major_radius_mm: 1.5,
            pitch_mm: 0.5,
            length_mm: if nonbinary { 4.3 } else { 4. },
            depth_mm: (3. - minor_mid) / 2.,
            crest_width_mm: 0.25 - (3. - pitch_mid) / 3_f64.sqrt(),
            lower_flank_deg: 30.,
            upper_flank_deg: 30.,
            starts,
            left_hand: left,
            root_radius_mm: None,
        }
    };
    AuthoredThread::new(
        spec,
        SourceOptions {
            position_limit_mm: 1e-7,
            normal_angle_limit_rad: 1e-7,
        },
    )
    .unwrap()
}
pub(super) fn fixture(left: bool, starts: u32, nonbinary: bool) -> Insert {
    Insert::new(
        source(true, left, starts, nonbinary),
        source(false, !left, starts, nonbinary),
        0.5,
        None,
    )
    .unwrap()
}
fn options() -> QueryOptions {
    QueryOptions {
        volume_error_mm3: 1e-7,
        maximum_cells: 10,
    }
}
#[test]
fn paired_actual_laws_volume_and_box_domains() {
    for left in [false, true] {
        for starts in [1, 2] {
            for nb in [false, true] {
                let body = fixture(left, starts, nb);
                assert_eq!(body.outer().spec(), source(true, left, starts, nb).spec());
                assert_eq!(body.inner().spec(), source(false, !left, starts, nb).spec());
                assert_eq!(body.minimum_wall_mm(), 0.5);
                let v = body.volume(options(), None).unwrap();
                let reference = |t: &AuthoredThread| {
                    let sum: f64 = t
                        .bands()
                        .into_iter()
                        .map(|b| {
                            let [a, c] = b.radius_mm;
                            (b.phase_mm[1] - b.phase_mm[0]) * (a * a + a * c + c * c) / 3.
                        })
                        .sum();
                    std::f64::consts::PI * t.full_length_mm() * sum / t.spec().pitch_mm
                };
                let expected = reference(body.outer()) - reference(body.inner());
                assert!(v.material_mm3[0] < expected && expected < v.material_mm3[1]);
                assert!(v.outer_mm3[0] > v.cavity_mm3[1]);
                let mut sample_sum = 0.;
                let n = 20000;
                for i in 0..n {
                    let u = (i as f64 + 0.5) / n as f64;
                    let square = |t: &AuthoredThread| {
                        let b = t.bands();
                        let q = u * t.spec().pitch_mm;
                        let b = b.iter().find(|b| q <= b.phase_mm[1]).unwrap();
                        let r = b.radius_mm[0]
                            + (b.radius_mm[1] - b.radius_mm[0]) * (q - b.phase_mm[0])
                                / (b.phase_mm[1] - b.phase_mm[0]);
                        r * r
                    };
                    sample_sum += square(body.outer()) - square(body.inner());
                }
                let numerical =
                    std::f64::consts::PI * body.full_length_mm() * sample_sum / n as f64;
                assert!((expected - numerical).abs() < 1e-5);
                for [x, y] in [[2., 0.], [-2., 1.], [-2., -1.], [0., 2.], [0., -2.]] {
                    assert_eq!(
                        body.classify_box([[x, x], [y, y], [1., 1.1]], None),
                        Ok(Classification::Inside)
                    );
                }
                assert_eq!(
                    body.classify_box([[0., 0.], [0., 0.], [1., 2.]], None),
                    Ok(Classification::Outside)
                );
                assert_eq!(
                    body.classify_box([[6., 7.], [0., 0.], [1., 2.]], None),
                    Ok(Classification::Outside)
                );
                assert_eq!(
                    body.classify_box([[2., 2.], [0., 0.], [-2., -1.]], None),
                    Ok(Classification::Outside)
                );
                assert_eq!(
                    body.classify_box([[2., 2.], [0., 0.], [0., 0.]], None),
                    Ok(Classification::UnresolvedBoundary)
                );
                assert_eq!(
                    body.classify_box([[1.4, 1.5], [-0.1, 0.1], [1., 1.1]], None),
                    Ok(Classification::UnresolvedBoundary)
                );
                assert_eq!(
                    body.classify_box([[3., 3.], [0., 0.], [1., 1.]], None),
                    Ok(Classification::Inside)
                );
            }
        }
    }
}
#[test]
fn strict_wall_and_length_domain_are_immutable() {
    let outer = source(true, false, 1, false);
    let inner = source(false, false, 1, false);
    let exact = outer.minimum_radius_mm() - inner.maximum_radius_mm();
    let body = Insert::new(outer, inner, exact, None).unwrap();
    assert_eq!(
        Insert::new(outer, inner, exact.next_up(), None),
        Err(Error::Unsupported)
    );
    for bad in [0., -1., f64::NAN, f64::INFINITY] {
        assert_eq!(Insert::new(outer, inner, bad, None), Err(Error::Invalid));
    }
    assert_eq!(
        Insert::new(outer, source(false, false, 1, true), 0.1, None),
        Err(Error::Unsupported)
    );
    assert_eq!(
        Insert::new(inner, outer, 0.1, None),
        Err(Error::Unsupported)
    );
    assert_eq!(body.outer(), &outer);
    assert_eq!(body.inner(), &inner);
}
#[test]
fn bounds_are_separation_only_and_touching_is_not_disjoint() {
    let body = fixture(false, 1, false);
    let b = body.bounds_mm();
    assert_eq!(b, [[-5., 5.], [-5., 5.], [0., 4.]]);
    assert_eq!(
        body.separated_from_box([[5_f64.next_up(), 6.], [0., 0.], [1., 2.]], None),
        Ok(true)
    );
    assert_eq!(
        body.separated_from_box([[-7., -6.], [0., 0.], [1., 2.]], None),
        Ok(true)
    );
    assert_eq!(
        body.separated_from_box([[5., 6.], [0., 0.], [1., 2.]], None),
        Ok(false)
    );
    assert_eq!(
        body.separated_from_box([[0., 0.], [0., 0.], [1., 2.]], None),
        Ok(false)
    );
    for query in [
        [[1., 0.], [0., 0.], [0., 0.]],
        [[f64::NAN, 1.], [0., 0.], [0., 0.]],
    ] {
        assert_eq!(body.separated_from_box(query, None), Err(Error::Invalid));
        assert_eq!(body.classify_box(query, None), Err(Error::Invalid));
    }
    assert_eq!(
        body.classify_box([[f64::MAX / 2., f64::MAX], [0., 0.], [1., 2.]], None),
        Ok(Classification::UnresolvedBoundary)
    );
}
#[test]
fn real_phase_samples_include_flanks_and_no_proxy_radius() {
    let body = fixture(false, 1, true);
    for source in [body.outer(), body.inner()] {
        for band in 0..5 {
            for u in [0.2, 0.7] {
                let p = source.point_in_band(band, 0, u, 2.).unwrap();
                let r = p.point_mm[0].hypot(p.point_mm[1]);
                let factor = if std::ptr::eq(source, body.outer()) {
                    0.999
                } else {
                    1.001
                };
                let point = [p.point_mm[0] * factor, p.point_mm[1] * factor, 2.];
                assert!(r > 0.);
                assert_eq!(
                    body.classify_box(point.map(|v| [v, v]), None),
                    Ok(Classification::Inside)
                );
            }
        }
    }
}
#[test]
fn query_budgets_cancel_and_retry_never_publish_partial_receipts() {
    let body = fixture(false, 1, false);
    let reference = body.volume(options(), None).unwrap();
    for error in [0., -1., f64::INFINITY, f64::NAN] {
        assert_eq!(
            body.volume(
                QueryOptions {
                    volume_error_mm3: error,
                    ..options()
                },
                None
            ),
            Err(Error::Invalid)
        );
    }
    assert_eq!(
        body.volume(
            QueryOptions {
                maximum_cells: 9,
                ..options()
            },
            None
        ),
        Err(Error::Limit)
    );
    assert_eq!(
        body.volume(
            QueryOptions {
                volume_error_mm3: 1e-30,
                ..options()
            },
            None
        ),
        Err(Error::Unresolved)
    );
    for stop in 1..=12 {
        let calls = Cell::new(0);
        let cancel = || {
            calls.set(calls.get() + 1);
            calls.get() == stop
        };
        assert_eq!(body.volume(options(), Some(&cancel)), Err(Error::Cancelled));
    }
    for stop in [1, 2] {
        let calls = Cell::new(0);
        let cancel = || {
            calls.set(calls.get() + 1);
            calls.get() == stop
        };
        assert_eq!(
            Insert::new(*body.outer(), *body.inner(), 0.1, Some(&cancel)),
            Err(Error::Cancelled)
        );
        let calls = Cell::new(0);
        let cancel = || {
            calls.set(calls.get() + 1);
            calls.get() == stop
        };
        assert_eq!(
            body.separated_from_box([[6., 7.], [0., 0.], [1., 2.]], Some(&cancel)),
            Err(Error::Cancelled)
        );
        let calls = Cell::new(0);
        let cancel = || {
            calls.set(calls.get() + 1);
            calls.get() == stop
        };
        assert_eq!(
            body.classify_box([[2., 2.], [0., 0.], [1., 2.]], Some(&cancel)),
            Err(Error::Cancelled)
        );
    }
    assert_eq!(body.volume(options(), None).unwrap(), reference);
}

#[test]
fn finite_huge_axial_box_cannot_become_a_periodic_membership_receipt() {
    let body = fixture(false, 1, false);
    assert_eq!(
        body.classify_box([[2., 2.], [0., 0.], [-f64::MAX, f64::MAX]], None),
        Ok(Classification::UnresolvedBoundary)
    );
}

#[test]
fn exact_finite_periodic_domain_does_not_borrow_unearned_binary_insert_caps() {
    let outer = source(true, false, 1, false);
    let inner = source(false, true, 1, false);
    assert!(Insert::new(outer, inner, 0.5, None).is_ok());
    let options = SourceOptions {
        position_limit_mm: 1e-7,
        normal_angle_limit_rad: 1e-7,
    };
    let end = crate::Rational::new(43, 10).unwrap();
    let a = AuthoredThread::new_with_exact_axial_domain(outer.spec(), end, options, None).unwrap();
    let b = AuthoredThread::new_with_exact_axial_domain(inner.spec(), end, options, None).unwrap();
    for (a, b) in [(a, inner), (outer, b), (a, b)] {
        assert_eq!(Insert::new(a, b, 0.5, None), Err(Error::Unsupported));
    }
    assert!(Insert::new(outer, inner, 0.5, None).is_ok());
}
