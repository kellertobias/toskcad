use super::*;
use crate::authored_hanger_bolt::{BoundaryOptions, Hanger};
fn value(row: &MeridianRow, t: f64) -> [I; 3] {
    let mut p = row
        .controls()
        .iter()
        .map(|v| v.source_bounds.map(hull))
        .collect::<Vec<_>>();
    for n in (1..p.len()).rev() {
        for i in 0..n {
            p[i] = std::array::from_fn(|k| p[i][k] * I::point(1. - t) + p[i + 1][k] * I::point(t));
        }
    }
    std::array::from_fn(|k| p[0][k] / p[0][3])
}
fn overlap(a: I, b: I) {
    assert!(
        a.valid() && b.valid() && a.lo <= b.hi && b.lo <= a.hi,
        "{a:?} {b:?}"
    );
}
#[test]
fn all_actual_profile_rows_and_physical_roles_share_complete_periodic_control_identity() {
    for left in [false, true] {
        for nb in [false, true] {
            let mut spec = crate::authored_hanger_bolt::tests::fixture(left, nb);
            spec.metric = AuthoredMetricThread::new(crate::authored_metric_thread::Spec {
                starts: 3,
                ..spec.metric.spec()
            })
            .unwrap();
            let s = Hanger::new(spec, None).unwrap();
            let b = s.boundary(BoundaryOptions::default(), None).unwrap();
            let f = SourceField::prepare(&b, DomainOptions::default(), None).unwrap();
            for i in 0..s.face_count() {
                let face = s.face_key(i).unwrap();
                let (c, band) = match face {
                    FaceKey::MetricThread { band, .. } => (Component::Metric, band),
                    FaceKey::WoodThread { band, .. } => (Component::Wood, band),
                    _ => continue,
                };
                let owner = b
                    .curves()
                    .iter()
                    .find_map(|v| match v.key {
                        CurveKey::Cap {
                            end: End::Start,
                            period,
                            ..
                        } if v.faces.contains(&face) => Some(period),
                        _ => None,
                    })
                    .unwrap();
                let pieces = if band == ThreadBandRole::Crest {
                    vec![BandPiece::CrestLowPhase, BandPiece::CrestHighPhase]
                } else {
                    vec![BandPiece::Whole]
                };
                let (pitch, length, n) = match c {
                    Component::Metric => {
                        let v = s.metric().spec();
                        (v.pitch_mm, v.length_mm, v.starts)
                    }
                    Component::Wood => {
                        let v = s.wood().spec();
                        (v.pitch_mm, v.length_mm, v.starts)
                    }
                };
                for piece in pieces {
                    let j = owner as i32
                        + if piece == BandPiece::CrestLowPhase {
                            1
                        } else {
                            0
                        };
                    let a = f
                        .meridian_row(
                            face,
                            piece,
                            CellSeam {
                                period: j,
                                upper_turn: true,
                            },
                            DomainOptions::default(),
                            None,
                        )
                        .unwrap();
                    let next = f
                        .meridian_row(
                            face,
                            piece,
                            CellSeam {
                                period: j + n as i32,
                                upper_turn: false,
                            },
                            DomainOptions::default(),
                            None,
                        )
                        .unwrap();
                    assert_eq!(a.key(), next.key());
                    assert_eq!(a.controls(), next.controls());
                    assert_eq!(a.work_used(), next.work_used());
                    assert_eq!(a.key().face, face);
                    let rounded = c == Component::Metric && band == ThreadBandRole::Root;
                    assert_eq!(a.degree(), if rounded { 2 } else { 1 });
                    if rounded {
                        assert_eq!(a.controls()[1].literal[3], 0.5);
                        assert!(
                            a.controls()[1].literal[0] / 0.5 < s.metric().minimum_radius_mm()[0],
                            "middle CONTROL must not be falsely treated as a physical source sample"
                        );
                    }
                    let index = match band {
                        ThreadBandRole::Crest => {
                            if piece == BandPiece::CrestLowPhase {
                                0
                            } else {
                                4
                            }
                        }
                        ThreadBandRole::Upper => 1,
                        ThreadBandRole::Root => 2,
                        ThreadBandRole::Lower => 3,
                    };
                    for t in [0., 0.25, 0.5, 0.75, 1.] {
                        let (q, r) = match c {
                            Component::Metric => {
                                if rounded {
                                    let v =
                                        s.metric().root_meridian().bounds([t, t], None).unwrap();
                                    (hull(v.point_mm[0]), hull(v.point_mm[1]))
                                } else {
                                    let v = s.metric().profile_bounds(index, [t, t]).unwrap();
                                    (hull(v.phase_mm), hull(v.radius_mm))
                                }
                            }
                            Component::Wood => {
                                let v = s.wood().bands()[index];
                                (
                                    p(v.phase_mm[0]) + (p(v.phase_mm[1]) - p(v.phase_mm[0])) * p(t),
                                    p(v.radius_mm[0])
                                        + (p(v.radius_mm[1]) - p(v.radius_mm[0])) * p(t),
                                )
                            }
                        };
                        let z = q + p(pitch) * p(f64::from(a.key().total_period));
                        let z = if c == Component::Metric {
                            z - p(length) - p(s.spec().plain_length_mm)
                        } else {
                            z
                        };
                        let actual = value(&a, t);
                        overlap(actual[0], r);
                        overlap(actual[1], p(0.));
                        overlap(actual[2], z);
                    }
                }
            }
        }
    }
}
#[test]
fn actual_row_domains_role_mismatches_and_resource_cancellation_are_atomic() {
    use std::cell::Cell;
    let mut spec = crate::authored_hanger_bolt::tests::fixture(false, false);
    spec.metric = AuthoredMetricThread::new(crate::authored_metric_thread::Spec {
        starts: 3,
        ..spec.metric.spec()
    })
    .unwrap();
    let s = Hanger::new(spec, None).unwrap();
    let b = s.boundary(BoundaryOptions::default(), None).unwrap();
    let options = DomainOptions::default();
    let f = SourceField::prepare(&b, options, None).unwrap();
    let face = FaceKey::MetricThread {
        start: 0,
        band: ThreadBandRole::Root,
    };
    let seam = CellSeam {
        period: 0,
        upper_turn: false,
    };
    let good = f
        .meridian_row(face, BandPiece::Whole, seam, options, None)
        .unwrap();
    for bad in [
        FaceKey::PlainWall,
        FaceKey::MetricThread {
            start: 99,
            band: ThreadBandRole::Root,
        },
    ] {
        assert!(matches!(
            f.meridian_row(bad, BandPiece::Whole, seam, options, None),
            Err(Error::Invalid)
        ));
    }
    assert!(matches!(
        f.meridian_row(face, BandPiece::CrestLowPhase, seam, options, None),
        Err(Error::Invalid)
    ));
    assert!(matches!(
        f.meridian_row(
            FaceKey::MetricThread {
                start: 0,
                band: ThreadBandRole::Crest
            },
            BandPiece::Whole,
            seam,
            options,
            None
        ),
        Err(Error::Invalid)
    ));
    for period in [i32::MIN, 1, f.cell_periods(Component::Metric)[1]] {
        assert!(matches!(
            f.meridian_row(
                face,
                BandPiece::Whole,
                CellSeam { period, ..seam },
                options,
                None
            ),
            Err(Error::Invalid)
        ));
    }
    let header = b
        .curves()
        .iter()
        .position(|c| {
            c.faces.contains(&face)
                && matches!(
                    c.key,
                    CurveKey::Cap {
                        end: End::Start,
                        ..
                    }
                )
        })
        .unwrap()
        + 2;
    for budget in [1, header, header + 1, good.work_used() - 1] {
        assert!(matches!(
            f.meridian_row(
                face,
                BandPiece::Whole,
                seam,
                DomainOptions {
                    maximum_work: budget,
                    ..options
                },
                None
            ),
            Err(Error::Limit)
        ));
    }
    let bytes = f.retained_bytes() + 2 * SCRATCH;
    assert!(
        f.meridian_row(
            face,
            BandPiece::Whole,
            seam,
            DomainOptions {
                maximum_bytes: bytes,
                ..options
            },
            None
        )
        .is_ok()
    );
    assert!(matches!(
        f.meridian_row(
            face,
            BandPiece::Whole,
            seam,
            DomainOptions {
                maximum_bytes: bytes - 1,
                ..options
            },
            None
        ),
        Err(Error::Limit)
    ));
    let n = Cell::new(0);
    f.meridian_row(
        face,
        BandPiece::Whole,
        seam,
        options,
        Some(&|| {
            n.set(n.get() + 1);
            false
        }),
    )
    .unwrap();
    for stop in [0, 1, header, header + 1, n.get() / 2, n.get() - 1] {
        let at = Cell::new(0);
        assert!(matches!(
            f.meridian_row(
                face,
                BandPiece::Whole,
                seam,
                options,
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
