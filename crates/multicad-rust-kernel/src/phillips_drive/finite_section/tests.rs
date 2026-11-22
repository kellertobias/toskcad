use super::*;
use std::cell::Cell;
fn i([lo, hi]: Bounds) -> I {
    I::new(lo, hi)
}
fn contains_zero(x: I) {
    assert!(x.valid() && x.lo <= 0. && x.hi >= 0., "{x:?}")
}
fn midpoint(p: Point2Bounds) -> [f64; 2] {
    p.map(|x| i(x).midpoint())
}
fn rot([x, y]: [f64; 2], q: u8) -> [f64; 2] {
    match q {
        0 => [x, y],
        1 => [-y, x],
        2 => [-x, -y],
        _ => [y, -x],
    }
}
fn support_residual(curve: SectionCurve, p: Point2Bounds) {
    let [x, y] = p.map(i);
    match curve {
        SectionCurve::CircleArc {
            radius_mm_bounds, ..
        } => contains_zero(x.square() + y.square() - i(radius_mm_bounds).square()),
        SectionCurve::Line {
            normal_bounds,
            offset_mm_bounds,
        } => contains_zero(i(normal_bounds[0]) * x + i(normal_bounds[1]) * y - i(offset_mm_bounds)),
        SectionCurve::EllipseArc {
            center_mm_bounds,
            semiaxes_mm_bounds,
            quarter_turn,
            ..
        } => {
            let delta = rotate(
                [x - i(center_mm_bounds[0]), y - i(center_mm_bounds[1])],
                (4 - quarter_turn) % 4,
            );
            contains_zero(
                (delta[0] / i(semiaxes_mm_bounds[0])).square()
                    + (delta[1] / i(semiaxes_mm_bounds[1])).square()
                    - I::point(1.),
            );
        }
    }
}
#[test]
fn metric_sections_have_shared_oriented_real_supports_and_agree_with_independent_source_queries() {
    let mut roles = std::collections::HashSet::new();
    let mut observed_counts = std::collections::HashSet::new();
    for index in 0..=5 {
        let (s, d) = super::super::tests::fixture(index);
        let c = Cutter::new(s, d).unwrap();
        for fraction in [0., 0.00001, 0.03, 0.3, 0.6, 0.85, 0.95] {
            let depth = c.envelope.tip.midpoint() * fraction;
            let section = c.finite_section(depth).unwrap();
            assert_eq!(section.source(), &c);
            assert_eq!(section.depth_mm(), depth);
            if fraction <= 0.03 {
                assert_eq!(section.envelope_role(), EnvelopeRole::EntryRoundover)
            } else if fraction == 0.3 {
                assert_eq!(section.envelope_role(), EnvelopeRole::OuterFlank)
            } else if fraction == 0.95 {
                assert_eq!(section.envelope_role(), EnvelopeRole::Floor)
            }
            assert_eq!(section.edges().len(), section.vertices().len());
            observed_counts.insert(section.edges().len());
            let mut signed_area = 0.;
            for (index, edge) in section.edges().iter().enumerate() {
                assert_eq!(
                    edge.vertices(),
                    [index, (index + 1) % section.vertices().len()]
                );
                let [a, b] = edge.vertices().map(|v| &section.vertices()[v]);
                assert_eq!(a.source_identity().0, edge.source_identity().0);
                let pa = midpoint(a.point_mm_bounds());
                let pb = midpoint(b.point_mm_bounds());
                let curve = edge.curve();
                support_residual(curve, a.point_mm_bounds());
                support_residual(curve, b.point_mm_bounds());
                let (point, tangent, area) = match curve {
                    SectionCurve::Line { .. } => {
                        roles.insert(0);
                        (
                            [(pa[0] + pb[0]) / 2., (pa[1] + pb[1]) / 2.],
                            [pb[0] - pa[0], pb[1] - pa[1]],
                            (pa[0] * pb[1] - pa[1] * pb[0]) / 2.,
                        )
                    }
                    SectionCurve::CircleArc {
                        radius_mm_bounds,
                        start_angle_rad_bounds,
                        sweep_rad_bounds,
                    } => {
                        roles.insert(1);
                        let r = i(radius_mm_bounds).midpoint();
                        let sweep = i(sweep_rad_bounds).midpoint();
                        assert!(sweep > 0.);
                        let angle = i(start_angle_rad_bounds).midpoint() + sweep / 2.;
                        (
                            [r * angle.cos(), r * angle.sin()],
                            [-angle.sin(), angle.cos()],
                            r * r * sweep / 2.,
                        )
                    }
                    SectionCurve::EllipseArc {
                        center_mm_bounds,
                        semiaxes_mm_bounds,
                        quarter_turn,
                        start_angle_rad_bounds,
                        sweep_rad_bounds,
                    } => {
                        roles.insert(2);
                        let center = midpoint(center_mm_bounds);
                        let [rx, ry] = semiaxes_mm_bounds.map(|v| i(v).midpoint());
                        assert!(rx > ry);
                        let sweep = i(sweep_rad_bounds).midpoint();
                        assert!(sweep < 0.);
                        let angle = i(start_angle_rad_bounds).midpoint() + sweep / 2.;
                        let p = rot([rx * angle.cos(), ry * angle.sin()], quarter_turn);
                        let t = rot([rx * angle.sin(), -ry * angle.cos()], quarter_turn);
                        (
                            [center[0] + p[0], center[1] + p[1]],
                            t,
                            (rx * ry * sweep + center[0] * (pb[1] - pa[1])
                                - center[1] * (pb[0] - pa[0]))
                                / 2.,
                        )
                    }
                };
                signed_area += area;
                if depth > 0. {
                    let norm = tangent[0].hypot(tangent[1]);
                    let epsilon = 1e-5_f64.min((pa[0] - pb[0]).hypot(pa[1] - pb[1]) / 100.);
                    assert!(epsilon > 0.);
                    let left = [-tangent[1] / norm * epsilon, tangent[0] / norm * epsilon];
                    let inside = [point[0] + left[0], point[1] + left[1], depth];
                    let outside = [point[0] - left[0], point[1] - left[1], depth];
                    assert!(
                        super::super::tests::oracle(s, d, inside),
                        "{index} {curve:?} {inside:?}"
                    );
                    assert!(
                        !super::super::tests::oracle(s, d, outside),
                        "{index} {curve:?} {outside:?}"
                    );
                    assert_eq!(c.classify_point_mm(inside).unwrap(), Classification::Inside);
                    assert_eq!(
                        c.classify_point_mm(outside).unwrap(),
                        Classification::Outside
                    );
                }
            }
            assert!(signed_area > 0.);
            if depth > 0. {
                let branch = if depth < c.envelope.join.lo {
                    0
                } else if depth < c.envelope.flank.lo {
                    1
                } else {
                    2
                };
                let expected = volume::Sections::new(&c)
                    .unwrap()
                    .area(J::d(I::point(depth)), branch)
                    .unwrap()
                    .v;
                let scale = c.envelope.entry.hi * c.envelope.entry.hi;
                assert!(
                    (signed_area - expected.midpoint()).abs() < 1e-11 * scale,
                    "{signed_area} {expected:?}"
                );
            }
        }
    }
    assert_eq!(roles.len(), 3);
    assert!(observed_counts.contains(&4));
    assert!(observed_counts.contains(&12));
    assert!(observed_counts.contains(&20));
}
#[test]
fn invalid_event_and_finite_domain_sections_refuse_and_cancel_without_mutating_source() {
    let (s, d) = super::super::tests::fixture(1);
    let c = Cutter::new(s, d).unwrap();
    let before = c.clone();
    for depth in [f64::NAN, -1., c.envelope.tip.hi + 1.] {
        assert!(matches!(c.finite_section(depth), Err(Error::Invalid)))
    }
    for depth in [
        c.envelope.join.midpoint(),
        c.envelope.flank.midpoint(),
        c.envelope.tip.midpoint(),
    ] {
        assert!(matches!(c.finite_section(depth), Err(Error::Unresolved)))
    }
    let short = Cutter::new(
        s,
        Domain {
            corner_reach_mm: 0.01,
            ..d
        },
    )
    .unwrap();
    assert!(matches!(short.finite_section(0.), Err(Error::Unresolved)));
    for stop in 1..=7 {
        let calls = Cell::new(0);
        assert!(matches!(
            c.finite_section_with_cancel(0., &|| {
                calls.set(calls.get() + 1);
                calls.get() == stop
            }),
            Err(Error::Cancelled)
        ));
    }
    // No-cut sections use the same four explicit source quarter seams and their
    // cancellation checkpoints rather than accidentally returning a cached wire.
    let depth = c.envelope.tip.lo * 0.95;
    for stop in 1..=6 {
        let calls = Cell::new(0);
        assert!(matches!(
            c.finite_section_with_cancel(depth, &|| {
                calls.set(calls.get() + 1);
                calls.get() == stop
            }),
            Err(Error::Cancelled)
        ));
    }
    assert_eq!(c, before);
    assert!(c.finite_section(0.).is_ok());
}

#[test]
fn real_contact_and_section_join_events_and_unrepresentable_rounds_refuse() {
    let (s, d) = super::super::tests::fixture(1);
    let c = Cutter::new(s, d).unwrap();
    let e = &c.envelope;
    let k = &c.corners[0];
    let cot = e.cos_c.midpoint() / e.sin_c.midpoint();
    let tilt = k.sin_beta.midpoint() / k.cos_beta.midpoint();
    let contact =
        (e.tip.midpoint() * cot - k.b_half.midpoint() - k.b_depth.midpoint() * tilt) / (cot - tilt);
    assert!(matches!(c.finite_section(contact), Err(Error::Unresolved)));
    let (mut lo, mut hi) = (e.flank.hi, contact);
    for _ in 0..80 {
        let z = (lo + hi) / 2.;
        let r = (e.tip.midpoint() - z) * cot;
        let x = k.b_half.midpoint()
            + (k.b_depth.midpoint() - z) * tilt
            + k.join_s.midpoint() / k.cos_beta.midpoint();
        if r * r > x * x + k.join_y.midpoint().powi(2) {
            lo = z
        } else {
            hi = z
        }
    }
    assert!(matches!(
        c.finite_section((lo + hi) / 2.),
        Err(Error::Unresolved)
    ));
    let (h0, domain) = super::super::tests::fixture(0);
    let tiny = Cutter::new(
        Spec {
            inner: Inner::Rounded { radius_mm: 1e-300 },
            ..h0
        },
        domain,
    )
    .unwrap();
    assert!(matches!(tiny.finite_section(0.), Err(Error::Unresolved)));
    // Directed endpoint-chart tests retain a zero-containing contact interval;
    // no nominal midpoint is allowed to turn it into a finite angular edge.
    assert!(intersection_angle(I::point(1.), I::new(0., 1e-15)).is_err());
    assert!(intersection_angle(I::point(1.), I::point(1.)).is_err());
    assert!(intersection_angle(I::point(1.), I::point(1e-300)).is_err());
    let a = intersection_angle(I::point(1.), I::point(0.2)).unwrap();
    assert!(a.lo <= 0.2_f64.atan() && a.hi >= 0.2_f64.atan());
}
