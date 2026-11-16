use super::*;
use crate::authored_hanger_bolt::{BoundaryOptions, Hanger};
use num_rational::BigRational as Q;
use num_traits::ToPrimitive;
#[test]
fn exact_turn_covers_preserve_nonbinary_ulp_boundaries_and_true_station_identity() {
    for left in [false, true] {
        for nb in [false, true] {
            for starts in [1, 3] {
                let mut spec = crate::authored_hanger_bolt::tests::fixture(left, nb);
                spec.metric = AuthoredMetricThread::new(crate::authored_metric_thread::Spec {
                    starts,
                    ..spec.metric.spec()
                })
                .unwrap();
                let s = Hanger::new(spec, None).unwrap();
                let b = s.boundary(BoundaryOptions::default(), None).unwrap();
                let f = SourceField::prepare(&b, DomainOptions::default(), None).unwrap();
                assert!(std::ptr::eq(f.source_boundary(), &b));
                assert!(f.work_used() < 16_000_000 && f.retained_bytes() >= b.retained_bytes());
                for c in [Component::Metric, Component::Wood] {
                    let (pitch, length, n) = match c {
                        Component::Metric => {
                            let x = s.metric().spec();
                            (x.pitch_mm, x.length_mm, x.starts)
                        }
                        Component::Wood => {
                            let x = s.wood().spec();
                            (x.pitch_mm, x.length_mm, x.starts)
                        }
                    };
                    let expected = (Q::from_float(length).unwrap()
                        / (Q::from_float(pitch).unwrap() * Q::from_integer(n.into())))
                    .ceil()
                    .to_integer()
                    .to_i32()
                    .unwrap();
                    assert_eq!(f.physical_turns(c), [-1, expected]);
                    let range = f.cell_periods(c);
                    assert_eq!(range[0], -(n as i32));
                    for cell in range[0]..=range[1] - (n as i32) {
                        for node in [
                            PhaseNode::CrestSplit,
                            PhaseNode::Knot(Knot::CrestEnd),
                            PhaseNode::Knot(Knot::UpperRoot),
                            PhaseNode::Knot(Knot::LowerRoot),
                            PhaseNode::Knot(Knot::CrestStart),
                        ] {
                            let a = f
                                .seam_point(c, cell, true, node, DomainOptions::default(), None)
                                .unwrap();
                            let next = f
                                .seam_point(
                                    c,
                                    cell + n as i32,
                                    false,
                                    node,
                                    DomainOptions::default(),
                                    None,
                                )
                                .unwrap();
                            assert_eq!(
                                a, next,
                                "structurally same source station must have identical data"
                            );
                            assert!(a.work_used <= 16_000_000 && a.source_error_mm < 1e-7);
                            assert_eq!(a.point_mm[1], 0.);
                            let err = (0..3)
                                .map(|k| {
                                    let x = hull(a.source_bounds_mm[k]);
                                    (x.lo - a.point_mm[k])
                                        .abs()
                                        .max((x.hi - a.point_mm[k]).abs())
                                })
                                .sum::<f64>();
                            assert!(err <= a.source_error_mm);
                            if let PhaseNode::Knot(knot) = node {
                                // Independent ORIGINAL finite 3D rail oracle;
                                // catches wrong body-offset/phase/hand even if
                                // both equivalent station calls share a bug.
                                let phase = match c {
                                    Component::Metric => {
                                        use crate::authored_metric_thread::Knot as K;
                                        let kind = match knot {
                                            Knot::CrestEnd => K::CrestEnd,
                                            Knot::UpperRoot => K::UpperRoot,
                                            Knot::LowerRoot => K::LowerRoot,
                                            Knot::CrestStart => K::CrestStart,
                                        };
                                        hull(s.metric().knot(kind).phase_mm)
                                    }
                                    Component::Wood => {
                                        let i = match knot {
                                            Knot::CrestEnd => 0,
                                            Knot::UpperRoot => 1,
                                            Knot::LowerRoot => 2,
                                            Knot::CrestStart => 3,
                                        };
                                        p(s.wood().bands()[i].phase_mm[1])
                                    }
                                };
                                // Planned parameter comes from the immutable
                                // source phase, NEVER from the returned point.
                                let local = phase + p(pitch) * p(f64::from(a.key.total_period));
                                let t = local / p(length);
                                if t.lo >= 0. && t.hi <= 1. {
                                    let original = b
                                        .curve_bounds(
                                            crate::authored_hanger_bolt::CurveKey::Rail {
                                                component: c,
                                                period: a.key.total_period.rem_euclid(n as i32)
                                                    as u32,
                                                knot,
                                            },
                                            [t.lo, t.hi],
                                            None,
                                        )
                                        .unwrap();
                                    for k in 0..3 {
                                        let expected = hull(original.point_mm[k]);
                                        let actual = hull(a.source_bounds_mm[k]);
                                        assert!(
                                            actual.lo <= expected.hi && expected.lo <= actual.hi
                                        );
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    // The floating L/p proposal rounds to exactly24 for this actualbinary
    // length; the exact source ratio is >24 and needs a genuine extra cover.
    assert_eq!(7.2 / 0.3, 24.);
    assert_eq!(cover(0.3, 7.2, 1).last_turn, 25);
    for length in [7.2_f64.next_down(), 7.2, 7.2_f64.next_up()] {
        let p = 0.3;
        let c = cover(p, length, 1);
        let expected = (Q::from_float(length).unwrap() / Q::from_float(p).unwrap())
            .ceil()
            .to_integer()
            .to_i32()
            .unwrap();
        assert_eq!(c.last_turn, expected);
    }
}
#[test]
fn distinct_original_knots_and_periods_never_alias_and_resources_cancel_atomically() {
    use std::cell::Cell;
    let s = Hanger::new(
        crate::authored_hanger_bolt::tests::fixture(false, true),
        None,
    )
    .unwrap();
    let b = s.boundary(BoundaryOptions::default(), None).unwrap();
    let opts = DomainOptions::default();
    let f = SourceField::prepare(&b, opts, None).unwrap();
    let a = f
        .seam_point(
            Component::Metric,
            0,
            false,
            PhaseNode::Knot(Knot::UpperRoot),
            opts,
            None,
        )
        .unwrap();
    let other = f
        .seam_point(
            Component::Metric,
            0,
            false,
            PhaseNode::Knot(Knot::LowerRoot),
            opts,
            None,
        )
        .unwrap();
    assert_ne!(a.key, other.key);
    assert_ne!(a.point_mm[2], other.point_mm[2]);
    let later = f
        .seam_point(
            Component::Metric,
            1,
            false,
            PhaseNode::Knot(Knot::UpperRoot),
            opts,
            None,
        )
        .unwrap();
    assert_ne!(a.key, later.key);
    assert_ne!(a.point_mm[2], later.point_mm[2]);
    for cell in [i32::MIN, i32::MAX] {
        assert_eq!(
            f.seam_point(
                Component::Wood,
                cell,
                false,
                PhaseNode::CrestSplit,
                opts,
                None
            ),
            Err(Error::Invalid)
        );
    }
    assert!(matches!(
        SourceField::prepare(
            &b,
            DomainOptions {
                maximum_work: 1,
                ..opts
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        SourceField::prepare(
            &b,
            DomainOptions {
                maximum_bytes: 1,
                ..opts
            },
            None
        ),
        Err(Error::Limit)
    ));
    let bytes = f.retained_bytes() + SCRATCH;
    assert!(
        f.seam_point(
            Component::Wood,
            0,
            false,
            PhaseNode::CrestSplit,
            DomainOptions {
                maximum_bytes: bytes,
                ..opts
            },
            None
        )
        .is_ok()
    );
    assert_eq!(
        f.seam_point(
            Component::Wood,
            0,
            false,
            PhaseNode::CrestSplit,
            DomainOptions {
                maximum_bytes: bytes - 1,
                ..opts
            },
            None
        ),
        Err(Error::Limit)
    );
    let n = Cell::new(0);
    f.seam_point(
        Component::Wood,
        0,
        false,
        PhaseNode::CrestSplit,
        opts,
        Some(&|| {
            n.set(n.get() + 1);
            false
        }),
    )
    .unwrap();
    for stop in 0..n.get() {
        let at = Cell::new(0);
        assert_eq!(
            f.seam_point(
                Component::Wood,
                0,
                false,
                PhaseNode::CrestSplit,
                opts,
                Some(&|| {
                    let k = at.get();
                    at.set(k + 1);
                    k == stop
                })
            ),
            Err(Error::Cancelled)
        );
    }
    let n = Cell::new(0);
    SourceField::prepare(
        &b,
        opts,
        Some(&|| {
            n.set(n.get() + 1);
            false
        }),
    )
    .unwrap();
    for stop in 0..n.get() {
        let at = Cell::new(0);
        assert!(matches!(
            SourceField::prepare(
                &b,
                opts,
                Some(&|| {
                    let k = at.get();
                    at.set(k + 1);
                    k == stop
                })
            ),
            Err(Error::Cancelled)
        ));
    }
    assert_eq!(*b.source(), s);
}
