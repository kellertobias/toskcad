use super::*;
use crate::authored_hanger_bolt::{BoundaryOptions, Hanger};
fn bezier(mut c: Vec<[I; 4]>, t: f64) -> [I; 4] {
    for n in (1..c.len()).rev() {
        for i in 0..n {
            c[i] = std::array::from_fn(|k| c[i][k] * p(1. - t) + c[i + 1][k] * p(t));
        }
    }
    c[0]
}
#[test]
fn complete_actual_fields_preserve_full_seams_and_source_images() {
    for (left, nb, starts) in [(false, false, 1), (true, true, 1), (false, false, 3)] {
        let mut spec = crate::authored_hanger_bolt::tests::fixture(left, nb);
        spec.metric = AuthoredMetricThread::new(crate::authored_metric_thread::Spec {
            starts,
            ..spec.metric.spec()
        })
        .unwrap();
        let s = Hanger::new(spec, None).unwrap();
        let b = s.boundary(BoundaryOptions::default(), None).unwrap();
        let options = CurveOptions::default();
        let field = BodyField::prepare(&b, options, None).unwrap();
        assert!(std::ptr::eq(field.source_field().source_boundary(), &b));
        assert!(field.work_used() < options.work && field.retained_bytes() < options.bytes);
        assert!(field.cells().len() > 100);
        eprintln!(
            "field left={left} nonbinary={nb} cells={} controls={} work={} retained_bytes={}",
            field.cells().len(),
            field
                .cells()
                .iter()
                .map(|c| c.controls().len())
                .sum::<usize>(),
            field.work_used(),
            field.retained_bytes()
        );
        for cell in field.cells() {
            let component = if matches!(cell.key.face, FaceKey::MetricThread { .. }) {
                Component::Metric
            } else {
                Component::Wood
            };
            let starts = match component {
                Component::Metric => s.metric().spec().starts,
                Component::Wood => s.wood().spec().starts,
            };
            if let Some(next) = field.cells().iter().find(|v| {
                v.key
                    == MeridianRowKey {
                        total_period: cell.key.total_period + starts as i32,
                        ..cell.key
                    }
            }) {
                for row in 0..=cell.degree {
                    assert_eq!(
                        cell.controls[row * cell.columns + cell.columns - 1],
                        next.controls[row * next.columns]
                    );
                }
            }
            let knots = field.angular_knots(component);
            let lower = field
                .source
                .meridian_row(
                    cell.key.face,
                    cell.key.piece,
                    CellSeam {
                        period: cell.key.total_period,
                        upper_turn: false,
                    },
                    DomainOptions::default(),
                    None,
                )
                .unwrap();
            let upper = field
                .source
                .meridian_row(
                    cell.key.face,
                    cell.key.piece,
                    CellSeam {
                        period: cell.key.total_period,
                        upper_turn: true,
                    },
                    DomainOptions::default(),
                    None,
                )
                .unwrap();
            for u in [0., 0.37, 1.] {
                let a = bezier(
                    lower
                        .controls()
                        .iter()
                        .map(|v| v.source_bounds.map(hull))
                        .collect(),
                    u,
                );
                let z = bezier(
                    upper
                        .controls()
                        .iter()
                        .map(|v| v.source_bounds.map(hull))
                        .collect(),
                    u,
                );
                for (piece, span) in knots.windows(2).enumerate() {
                    let t = 0.43;
                    let v = p(span[0]) + (p(span[1]) - p(span[0])) * p(t);
                    let angle =
                        v * I::new(
                            std::f64::consts::TAU.next_down(),
                            std::f64::consts::TAU.next_up(),
                        ) * p(if left { -1. } else { 1. });
                    let expected = [
                        a[0] / a[3] * angle.cos(),
                        a[0] / a[3] * angle.sin(),
                        (a[2] + (z[2] - a[2]) * v) / a[3],
                    ];
                    let columns = (0..=10)
                        .map(|j| {
                            bezier(
                                (0..=cell.degree)
                                    .map(|i| {
                                        cell.controls[i * cell.columns + piece * 10 + j]
                                            .literal
                                            .map(p)
                                    })
                                    .collect(),
                                u,
                            )
                        })
                        .collect();
                    let h = bezier(columns, t);
                    for k in 0..3 {
                        let actual = h[k] / h[3];
                        let d = expected[k] - actual;
                        assert!(
                            d.magnitude() <= cell.bounds.position_error_mm + 1e-11,
                            "{:?} {d:?} {:?}",
                            cell.key,
                            cell.bounds
                        );
                    }
                }
            }
        }
    }
}
#[test]
fn aggregate_limits_and_cancellation_do_not_publish_partial_fields() {
    let s = Hanger::new(
        crate::authored_hanger_bolt::tests::fixture(false, false),
        None,
    )
    .unwrap();
    let b = s.boundary(BoundaryOptions::default(), None).unwrap();
    let options = CurveOptions::default();
    assert!(matches!(
        BodyField::prepare(&b, CurveOptions { work: 0, ..options }, None),
        Err(Error::Invalid)
    ));
    assert!(matches!(
        BodyField::prepare(
            &b,
            CurveOptions {
                work: 100,
                ..options
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        BodyField::prepare(
            &b,
            CurveOptions {
                bytes: 1,
                ..options
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        BodyField::prepare(&b, options, Some(&|| true)),
        Err(Error::Cancelled)
    ));
    let n = std::cell::Cell::new(0);
    assert!(matches!(
        BodyField::prepare(
            &b,
            options,
            Some(&|| {
                n.set(n.get() + 1);
                n.get() > 10000
            })
        ),
        Err(Error::Cancelled)
    ));
    assert!(matches!(
        BodyField::prepare(
            &b,
            CurveOptions {
                position_error_mm: 1e-300,
                derivative_error_mm: 1e-300,
                ..options
            },
            None
        ),
        Err(Error::Unresolved) | Err(Error::Limit)
    ));
}
#[test]
fn rational_weight_and_derivative_residuals_bound_actual_axial_amplification() {
    let s = Hanger::new(
        crate::authored_hanger_bolt::tests::fixture(false, true),
        None,
    )
    .unwrap();
    let b = s.boundary(BoundaryOptions::default(), None).unwrap();
    let field = BodyField::prepare(&b, CurveOptions::default(), None).unwrap();
    let cell = field.cells().iter().find(|c| c.degree == 2).unwrap();
    assert_eq!(cell.key(), cell.key);
    assert_eq!(cell.meridian_degree(), cell.degree);
    assert_eq!(cell.angular_columns(), cell.columns);
    assert_eq!(cell.controls(), cell.controls);
    assert_eq!(cell.bounds(), cell.bounds);
    let knots = field.angular_knots(Component::Metric);
    let mut candidate = cell.controls.clone();
    for c in &mut candidate {
        c.literal[3] += 0.001;
    }
    let mut work = Work::new(16_000_000, None).unwrap();
    let domain = || NetDomain {
        degree: cell.degree,
        columns: cell.columns,
        knots,
        angular_error: [0.; 2],
        radial: [3., 3.],
    };
    let loose = CurveOptions {
        position_error_mm: 1.,
        derivative_error_mm: 10.,
        ..CurveOptions::default()
    };
    let result = certify(&candidate, domain(), loose, &mut work).unwrap();
    assert!(result.position_error_mm > 0.005);
    for u in [0., 0.3, 0.7, 1.] {
        for (piece, _) in knots.windows(2).enumerate() {
            let column = |data: &[FieldControl], literal: bool| {
                let rows = (0..=10)
                    .map(|j| {
                        bezier(
                            (0..=cell.degree)
                                .map(|i| {
                                    let c = data[i * cell.columns + piece * 10 + j];
                                    if literal {
                                        c.literal.map(p)
                                    } else {
                                        c.coefficient.map(hull)
                                    }
                                })
                                .collect(),
                            u,
                        )
                    })
                    .collect();
                bezier(rows, 0.47)
            };
            let source = column(&candidate, false);
            let literal = column(&candidate, true);
            let d = l1(std::array::from_fn(|k| {
                source[k] / source[3] - literal[k] / literal[3]
            }));
            assert!(d.hi <= result.position_error_mm);
        }
    }
    assert!(matches!(
        certify(&candidate, domain(), CurveOptions::default(), &mut work),
        Err(Error::Unresolved)
    ));
    candidate[0].literal[3] = -1.;
    assert!(matches!(
        certify(&candidate, domain(), loose, &mut work),
        Err(Error::Unresolved)
    ));
}
#[test]
fn source_precision_floors_and_checkpoint_prefixes_preserve_atomic_refusal() {
    let source = Hanger::new(
        crate::authored_hanger_bolt::tests::fixture(false, false),
        None,
    )
    .unwrap();
    let boundary = source.boundary(BoundaryOptions::default(), None).unwrap();
    let options = CurveOptions::default();
    assert!(matches!(
        BodyField::prepare(
            &boundary,
            CurveOptions {
                position_error_mm: f64::from_bits(1),
                ..options
            },
            None
        ),
        Err(Error::Unresolved)
    ));
    let mut w = Work::new(1, None).unwrap();
    w.charge(1).unwrap();
    assert!(matches!(nested(options, &w), Err(Error::Limit)));
    let huge = BigRational::from_integer(num_bigint::BigInt::from(1) << 1024);
    let mut w = Work::new(2000, None).unwrap();
    assert!(matches!(enclosure(&huge, &mut w), Err(Error::Unresolved)));
    let mut spec = source.spec();
    spec.plain_length_mm = 1e7;
    let far = Hanger::new(spec, None).unwrap();
    let b = far.boundary(BoundaryOptions::default(), None).unwrap();
    assert!(matches!(
        BodyField::prepare(&b, options, None),
        Err(Error::Unresolved)
    ));
    // Exercise cancellations within successive real preparation/row/proof
    // prefixes. No private source is forged and no stage bypass is installed.
    for prefix in (3000..18000).step_by(32) {
        let n = std::cell::Cell::new(0);
        assert!(matches!(
            BodyField::prepare(
                &boundary,
                options,
                Some(&|| {
                    n.set(n.get() + 1);
                    n.get() > prefix
                })
            ),
            Err(Error::Cancelled)
        ));
    }
    assert_eq!(source.spec(), boundary.source().spec());
}
#[test]
fn cancellation_at_measured_first_row_frontier_never_publishes_a_partial_net() {
    let source = Hanger::new(
        crate::authored_hanger_bolt::tests::fixture(false, false),
        None,
    )
    .unwrap();
    let boundary = source.boundary(BoundaryOptions::default(), None).unwrap();
    let options = CurveOptions::default();
    let m = source.metric().spec();
    let scale = (p(m.major_radius_mm) + p(m.pitch_mm) * p(f64::from(m.starts)) + p(1.)).hi;
    let n = std::cell::Cell::new(0);
    let callback = || {
        n.set(n.get() + 1);
        false
    };
    let mut w = Work::new(options.work, Some(&callback)).unwrap();
    let mut memory = Memory::new(0, options.bytes).unwrap();
    let angular = fit_reduced_unit(
        false,
        CurveOptions {
            position_error_mm: options.position_error_mm / (64. * scale),
            derivative_error_mm: options.derivative_error_mm / (64. * scale),
            ..options
        },
        &mut memory,
        &mut w,
    )
    .unwrap();
    let frontier = n.get();
    assert!(angular.knots().len() > 2);
    for prefix in (frontier + 128..frontier + 4096).step_by(8) {
        let count = std::cell::Cell::new(0);
        assert!(matches!(
            BodyField::prepare(
                &boundary,
                options,
                Some(&|| {
                    count.set(count.get() + 1);
                    count.get() > prefix
                })
            ),
            Err(Error::Cancelled)
        ));
    }
}
