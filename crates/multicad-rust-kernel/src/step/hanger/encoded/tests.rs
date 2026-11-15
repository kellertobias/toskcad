use super::*;
use crate::authored_hanger_bolt::{BoundaryOptions, Hanger};
#[test]
fn actual_encoded_scalar_field_is_small_and_all_seams_are_identical() {
    for (left, nb) in [(false, false), (true, true)] {
        let s = Hanger::new(crate::authored_hanger_bolt::tests::fixture(left, nb), None).unwrap();
        let b = s.boundary(BoundaryOptions::default(), None).unwrap();
        let f = EncodedField::prepare(&b, CurveOptions::default(), None).unwrap();
        assert!(std::ptr::eq(f.source_field().source_boundary(), &b));
        assert_eq!(f.meridians().len(), 10);
        eprintln!(
            "encoded left={left} nonbinary={nb} cells={} work={} bytes={} angular_pieces={}",
            f.cells().len(),
            f.work_used(),
            f.retained_bytes(),
            f.angular().coefficients().len()
        );
        for (i, cell) in f.cells().iter().enumerate() {
            let profile = &f.meridians()[cell.meridian];
            let n = if profile.component() == Component::Metric {
                s.metric().spec().starts
            } else {
                s.wood().spec().starts
            };
            let next = f.cells().iter().position(|v| {
                v.source_key
                    == MeridianRowKey {
                        total_period: cell.source_key.total_period + n as i32,
                        ..cell.source_key
                    }
            });
            if let Some(j) = next {
                for row in 0..=profile.degree() {
                    let a = f
                        .control(
                            i,
                            row,
                            f.angular().coefficients().len() - 1,
                            10,
                            DomainOptions::default(),
                            None,
                        )
                        .unwrap();
                    let z = f
                        .control(j, row, 0, 0, DomainOptions::default(), None)
                        .unwrap();
                    assert_eq!(a.homogeneous, z.homogeneous);
                }
            }
        }
    }
}
fn qi(q: &BigRational) -> I {
    use num_traits::ToPrimitive;
    let x = q.to_f64().unwrap();
    let a = x.next_down();
    let b = x.next_up();
    assert!(BigRational::from_float(a).unwrap() <= *q && *q <= BigRational::from_float(b).unwrap());
    I::new(a, b)
}
#[test]
fn actual_original_source_position_and_both_chart_derivatives_bound_encoded_tensors() {
    use super::cap::tests::{bezier, number};
    use num_traits::One;
    type Q = BigRational;
    for left in [false, true] {
        let s = Hanger::new(
            crate::authored_hanger_bolt::tests::fixture(left, true),
            None,
        )
        .unwrap();
        let b = s.boundary(BoundaryOptions::default(), None).unwrap();
        let f = EncodedField::prepare(&b, CurveOptions::default(), None).unwrap();
        assert!(f.angular_jacobian_lower() > 6.);
        for (index, meridian) in f.meridians().iter().enumerate() {
            let cell = f.cells().iter().position(|c| c.meridian == index).unwrap();
            let j = f.cells()[cell].source_key.total_period;
            let band = match meridian.band() {
                ThreadBandRole::Upper => 1,
                ThreadBandRole::Root => 2,
                ThreadBandRole::Lower => 3,
                ThreadBandRole::Crest => {
                    if meridian.piece() == BandPiece::CrestLowPhase {
                        0
                    } else {
                        4
                    }
                }
            };
            assert_eq!(meridian.phase_coefficients().len(), meridian.degree() + 1);
            assert_eq!(
                meridian.radial_coefficients().len(),
                meridian.weights().len()
            );
            let (pitch, starts) = if meridian.component() == Component::Metric {
                let spec = s.metric().spec();
                (spec.pitch_mm, spec.starts)
            } else {
                let spec = s.wood().spec();
                (spec.pitch_mm, spec.starts)
            };
            for u in [0., 0.5, 1.] {
                for piece in 0..f.angular().coefficients().len() {
                    let uexact = Q::from_float(u).unwrap();
                    let half = Q::one() / Q::from_integer(2.into());
                    let va = Q::from_float(f.angular().knots()[piece]).unwrap();
                    let width = Q::from_float(f.angular().knots()[piece + 1]).unwrap() - &va;
                    let vexact = &va + &width * &half;
                    let v = qi(&vexact);
                    let mut h = Vec::new();
                    let mut hu: Vec<[Q; 4]> = Vec::new();
                    for column in 0..=10 {
                        let controls = (0..=meridian.degree())
                            .map(|row| {
                                f.control(cell, row, piece, column, DomainOptions::default(), None)
                                    .unwrap()
                                    .homogeneous
                                    .each_ref()
                                    .map(number)
                            })
                            .collect::<Vec<_>>();
                        let differences: Vec<[Q; 4]> = controls
                            .windows(2)
                            .map(|c| {
                                std::array::from_fn(|k| {
                                    (&c[1][k] - &c[0][k])
                                        * Q::from_integer(meridian.degree().into())
                                })
                            })
                            .collect();
                        hu.push(bezier(differences, &uexact));
                        h.push(bezier(controls, &uexact));
                    }
                    let differences: Vec<[Q; 4]> = h
                        .windows(2)
                        .map(|c| {
                            std::array::from_fn(|k| {
                                (&c[1][k] - &c[0][k]) * Q::from_integer(10.into()) / &width
                            })
                        })
                        .collect();
                    let hv = bezier(differences, &half);
                    let hu = bezier(hu, &half);
                    let h = bezier(h, &half);
                    let position: [Q; 3] = std::array::from_fn(|k| &h[k] / &h[3]);
                    let du: [Q; 3] = std::array::from_fn(|k| {
                        (&hu[k] * &h[3] - &h[k] * &hu[3]) / (&h[3] * &h[3])
                    });
                    let dv: [Q; 3] = std::array::from_fn(|k| {
                        (&hv[k] * &h[3] - &h[k] * &hv[3]) / (&h[3] * &h[3])
                    });
                    let (q, r, dq, dr) = if meridian.component() == Component::Metric && band == 2 {
                        let actual = s.metric().root_meridian().bounds([u, u], None).unwrap();
                        (
                            hull(actual.point_mm[0]),
                            hull(actual.point_mm[1]),
                            hull(actual.first_derivative_mm[0]),
                            hull(actual.first_derivative_mm[1]),
                        )
                    } else if meridian.component() == Component::Metric {
                        use crate::authored_metric_thread::Knot;
                        let knots = [
                            Knot::Begin,
                            Knot::CrestEnd,
                            Knot::UpperRoot,
                            Knot::LowerRoot,
                            Knot::CrestStart,
                            Knot::End,
                        ];
                        let actual = s.metric().profile_bounds(band, [u, u]).unwrap();
                        let dq = hull(s.metric().knot(knots[band + 1]).phase_mm)
                            - hull(s.metric().knot(knots[band]).phase_mm);
                        (
                            hull(actual.phase_mm),
                            hull(actual.radius_mm),
                            dq,
                            hull(actual.slope) * dq,
                        )
                    } else {
                        let b = s.wood().bands()[band];
                        let dq = p(b.phase_mm[1]) - p(b.phase_mm[0]);
                        let dr = p(b.radius_mm[1]) - p(b.radius_mm[0]);
                        (
                            p(b.phase_mm[0]) + dq * p(u),
                            p(b.radius_mm[0]) + dr * p(u),
                            dq,
                            dr,
                        )
                    };
                    let tau = I::new(
                        std::f64::consts::TAU.next_down(),
                        std::f64::consts::TAU.next_up(),
                    );
                    let rate = tau * p(if left { -1. } else { 1. });
                    let angle = rate * v;
                    let lead = p(pitch) * p(f64::from(starts));
                    let mut z = q + p(pitch) * p(f64::from(j)) + lead * v;
                    if meridian.component() == Component::Metric {
                        z = z - p(s.metric().full_length_mm()) - p(s.spec().plain_length_mm);
                    }
                    let expected = [r * angle.cos(), r * angle.sin(), z];
                    let expected_u = [dr * angle.cos(), dr * angle.sin(), dq];
                    let expected_v = [-r * angle.sin() * rate, r * angle.cos() * rate, lead];
                    for (a, b, bound) in [
                        (position, expected, meridian.bounds().position_error_mm),
                        (du, expected_u, meridian.bounds().derivative_error_mm[0]),
                        (dv, expected_v, meridian.bounds().derivative_error_mm[1]),
                    ] {
                        for k in 0..3 {
                            let actual = qi(&a[k]);
                            let gap = (actual.lo - b[k].hi).max(b[k].lo - actual.hi).max(0.);
                            assert!(
                                gap <= bound,
                                "{index} {u} {piece} axis{k} gap{gap} bound{bound}"
                            );
                        }
                    }
                }
            }
        }
    }
}
#[test]
fn scalar_field_resource_precision_and_query_failures_are_atomic() {
    let s = Hanger::new(
        crate::authored_hanger_bolt::tests::fixture(false, false),
        None,
    )
    .unwrap();
    let b = s.boundary(BoundaryOptions::default(), None).unwrap();
    let o = CurveOptions::default();
    for options in [
        CurveOptions { work: 1, ..o },
        CurveOptions { bytes: 1, ..o },
    ] {
        assert!(matches!(
            EncodedField::prepare(&b, options, None),
            Err(Error::Limit)
        ));
    }
    assert!(matches!(
        EncodedField::prepare(&b, CurveOptions { work: 0, ..o }, None),
        Err(Error::Invalid)
    ));
    assert!(matches!(
        EncodedField::prepare(
            &b,
            CurveOptions {
                position_error_mm: f64::from_bits(1),
                ..o
            },
            None
        ),
        Err(Error::Unresolved)
    ));
    assert!(matches!(
        EncodedField::prepare(&b, o, Some(&|| true)),
        Err(Error::Cancelled)
    ));
    let field = EncodedField::prepare(&b, o, None).unwrap();
    assert!(matches!(
        field.control(usize::MAX, 0, 0, 0, DomainOptions::default(), None),
        Err(Error::Invalid)
    ));
    assert!(matches!(
        field.control(0, 10, 0, 0, DomainOptions::default(), None),
        Err(Error::Invalid)
    ));
    assert!(matches!(
        field.control(0, 0, usize::MAX, 0, DomainOptions::default(), None),
        Err(Error::Invalid)
    ));
    assert!(matches!(
        field.control(0, 0, 0, 11, DomainOptions::default(), None),
        Err(Error::Invalid)
    ));
    assert!(matches!(
        field.control(0, 0, 0, 0, DomainOptions::default(), Some(&|| true)),
        Err(Error::Cancelled)
    ));
    assert!(matches!(
        field.cap_trace(usize::MAX, End::Start, DomainOptions::default(), None),
        Err(Error::Invalid)
    ));
    assert!(matches!(
        field.cap_trace(0, End::Start, DomainOptions::default(), Some(&|| true)),
        Err(Error::Cancelled)
    ));
    for prefix in (0..2000).step_by(11) {
        let count = std::cell::Cell::new(0);
        assert!(matches!(
            EncodedField::prepare(
                &b,
                o,
                Some(&|| {
                    count.set(count.get() + 1);
                    count.get() > prefix
                })
            ),
            Err(Error::Cancelled)
        ));
    }
}
#[test]
fn multistart_and_every_actual_publication_checkpoint_preserve_atomicity() {
    let mut spec = crate::authored_hanger_bolt::tests::fixture(true, true);
    spec.metric = AuthoredMetricThread::new(crate::authored_metric_thread::Spec {
        starts: 3,
        ..spec.metric.spec()
    })
    .unwrap();
    let s = Hanger::new(spec, None).unwrap();
    let b = s.boundary(BoundaryOptions::default(), None).unwrap();
    let o = CurveOptions::default();
    let n = std::cell::Cell::new(0);
    let field = EncodedField::prepare(
        &b,
        o,
        Some(&|| {
            n.set(n.get() + 1);
            false
        }),
    )
    .unwrap();
    let checkpoints = n.get();
    assert!(checkpoints > 100);
    for prefix in 0..checkpoints {
        let seen = std::cell::Cell::new(0);
        assert!(matches!(
            EncodedField::prepare(
                &b,
                o,
                Some(&|| {
                    seen.set(seen.get() + 1);
                    seen.get() > prefix
                })
            ),
            Err(Error::Cancelled)
        ));
    }
    let n = std::cell::Cell::new(0);
    field
        .control(
            0,
            0,
            0,
            4,
            DomainOptions::default(),
            Some(&|| {
                n.set(n.get() + 1);
                false
            }),
        )
        .unwrap();
    for prefix in 0..n.get() {
        let seen = std::cell::Cell::new(0);
        assert!(matches!(
            field.control(
                0,
                0,
                0,
                4,
                DomainOptions::default(),
                Some(&|| {
                    seen.set(seen.get() + 1);
                    seen.get() > prefix
                })
            ),
            Err(Error::Cancelled)
        ));
    }
    let target = field
        .cells()
        .iter()
        .enumerate()
        .find_map(|(i, c)| {
            let m = &field.meridians()[c.meridian];
            if m.component() == Component::Metric
                && m.band() == ThreadBandRole::Crest
                && field
                    .cap_trace(i, End::Start, DomainOptions::default(), None)
                    .unwrap()
                    .is_some()
            {
                Some(i)
            } else {
                None
            }
        })
        .unwrap();
    let n = std::cell::Cell::new(0);
    field
        .cap_trace(
            target,
            End::Start,
            DomainOptions::default(),
            Some(&|| {
                n.set(n.get() + 1);
                false
            }),
        )
        .unwrap();
    for prefix in 0..n.get() {
        let seen = std::cell::Cell::new(0);
        assert!(matches!(
            field.cap_trace(
                target,
                End::Start,
                DomainOptions::default(),
                Some(&|| {
                    seen.set(seen.get() + 1);
                    seen.get() > prefix
                })
            ),
            Err(Error::Cancelled)
        ));
    }
    assert_eq!(s.spec(), b.source().spec());
}
#[test]
fn actual_source_precision_floors_and_first_cell_allocation_frontier_refuse() {
    let spec = crate::authored_hanger_bolt::tests::fixture(false, false);
    let source = Hanger::new(spec, None).unwrap();
    let b = source.boundary(BoundaryOptions::default(), None).unwrap();
    let o = CurveOptions::default();
    let field = EncodedField::prepare(&b, o, None).unwrap();
    let wood = source.wood().spec();
    let metric = source.metric().spec();
    let scale = (p(source
        .wood()
        .maximum_radius_mm()
        .max(source.metric().maximum_radius_mm()))
        + p(wood.pitch_mm.max(metric.pitch_mm)) * p(f64::from(wood.starts.max(metric.starts)))
        + p(1.))
    .hi;
    let mut w = Work::new(o.work, None).unwrap();
    w.charge(1 + field.source_field().work_used()).unwrap();
    let mut memory = Memory::new(
        b.retained_bytes() + size_of::<EncodedField<'_>>() + SCRATCH,
        o.bytes,
    )
    .unwrap();
    let _angular = fit_reduced_unit(
        false,
        CurveOptions {
            position_error_mm: o.position_error_mm / (4. * scale),
            derivative_error_mm: o.derivative_error_mm / (4. * scale),
            ..o
        },
        &mut memory,
        &mut w,
    )
    .unwrap();
    w.charge(64).unwrap();
    let frontier = w.used + field.cells().len() + 1;
    assert!(matches!(
        EncodedField::prepare(
            &b,
            CurveOptions {
                work: frontier - 1,
                ..o
            },
            None
        ),
        Err(Error::Limit)
    ));
    // A real accepted source with a tiny positive core cannot use a global
    // normal lower bound at this requested precision; retain that refusal.
    let m = spec.metric.spec();
    let offset = m.major_radius_mm - hull(spec.metric.minimum_radius_mm()).midpoint();
    let metric = AuthoredMetricThread::new(crate::authored_metric_thread::Spec {
        major_radius_mm: offset + 1e-12,
        ..m
    })
    .unwrap();
    let small = Hanger::new(crate::authored_hanger_bolt::Spec { metric, ..spec }, None).unwrap();
    let b = small.boundary(BoundaryOptions::default(), None).unwrap();
    assert!(matches!(
        EncodedField::prepare(&b, o, None),
        Err(Error::Unresolved)
    ));
}

#[test]
fn actual_adjacent_band_rails_share_the_exact_encoded_radius_and_station() {
    for (left, nonbinary) in [(false, false), (true, true)] {
        let source = Hanger::new(
            crate::authored_hanger_bolt::tests::fixture(left, nonbinary),
            None,
        )
        .unwrap();
        let boundary = source.boundary(BoundaryOptions::default(), None).unwrap();
        let field = EncodedField::prepare(&boundary, CurveOptions::default(), None).unwrap();
        for base in [0, 5] {
            for band in 0..5 {
                let a = &field.meridians[base + band];
                let b = &field.meridians[base + (band + 1) % 5];
                assert_eq!(
                    a.radius[a.degree()],
                    b.radius[0],
                    "left={left} nonbinary={nonbinary} band={band}"
                );
                assert_eq!(a.weights[a.degree()], b.weights[0]);
                let i = field
                    .cells
                    .iter()
                    .position(|c| c.meridian == base + band)
                    .unwrap();
                let period = field.cells[i].source_key.total_period + i32::from(band == 4);
                let j = field
                    .cells
                    .iter()
                    .position(|c| {
                        c.meridian == base + (band + 1) % 5 && c.source_key.total_period == period
                    })
                    .unwrap();
                for (piece, column) in [(0, 0), (field.angular.coefficients().len() - 1, 10)] {
                    let x = field
                        .control(i, a.degree(), piece, column, DomainOptions::default(), None)
                        .unwrap();
                    let y = field
                        .control(j, 0, piece, column, DomainOptions::default(), None)
                        .unwrap();
                    assert_eq!(x.homogeneous, y.homogeneous);
                }
            }
        }
    }
}
