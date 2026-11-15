//! Exact affine-row angular chart prototype, never a rounded-root alias.
use super::*;
use crate::authored_hanger_bolt::{BoundaryOptions, Hanger};
use crate::step::hanger::encoded::cap::tests::{bezier, number};
use num_traits::{One, Zero};
type Q = BigRational;
use super::affine::{build, relation};
#[test]
fn exact_angular_chart_retains_actual_affine_tensor_and_reverses_source_parameter() {
    for (left, nb) in [(false, false), (true, true)] {
        let s = Hanger::new(crate::authored_hanger_bolt::tests::fixture(left, nb), None).unwrap();
        let b = s.boundary(BoundaryOptions::default(), None).unwrap();
        let f = EncodedField::prepare(&b, CurveOptions::default(), None).unwrap();
        let mut w = Work::new(16_000_000, None).unwrap();
        w.charge(f.work_used()).unwrap();
        let mut count = 0;
        for index in 0..f.cells().len() {
            for end in [End::Start, End::End] {
                for piece in 0..f.angular.coefficients().len() {
                    let Some(cap) = build(&f, index, end, piece, &mut w).unwrap() else {
                        continue;
                    };
                    count += 1;
                    assert!(cap.reversed);
                    assert!(cap.derivative[0].0 < Q::zero());
                    assert!(cap.derivative[1].0 > Q::zero());
                    let row = &f.meridians[f.cells[index].meridian];
                    assert_eq!(row.degree(), 1);
                    let net = (0..=1)
                        .map(|i| {
                            (0..=10)
                                .map(|j| {
                                    f.control(index, i, piece, j, DomainOptions::default(), None)
                                        .unwrap()
                                        .homogeneous
                                        .each_ref()
                                        .map(number)
                                })
                                .collect::<Vec<_>>()
                        })
                        .collect::<Vec<_>>();
                    for k in 0..=4 {
                        let t = &cap.trim[0].0
                            + (&cap.trim[1].0 - &cap.trim[0].0) * Q::new(k.into(), 4.into());
                        let u = &cap.uv[0][0].0 + &cap.derivative[0].0 * &t;
                        let v = &cap.uv[0][1].0 + &cap.derivative[1].0 * &t;
                        assert!(u >= Q::zero() && u <= Q::one());
                        let h = bezier(
                            cap.controls
                                .iter()
                                .map(|c| c.clone().map(|n| n.0))
                                .collect(),
                            &t,
                        );
                        let tensor =
                            bezier(net.iter().map(|r| bezier(r.clone(), &t)).collect(), &u);
                        assert_eq!(
                            h[..3].iter().map(|x| x / &h[3]).collect::<Vec<_>>(),
                            tensor[..3]
                                .iter()
                                .map(|x| x / &tensor[3])
                                .collect::<Vec<_>>()
                        );
                        let curve = cap
                            .controls
                            .iter()
                            .map(|c| c.clone().map(|n| n.0))
                            .collect::<Vec<_>>();
                        let (_, d) = super::restrict::tests::jet(curve, &t);
                        let hv = bezier(
                            net.iter()
                                .map(|r| {
                                    bezier(
                                        r.windows(2)
                                            .map(|p| {
                                                std::array::from_fn::<_, 4, _>(|axis| {
                                                    (&p[1][axis] - &p[0][axis])
                                                        * Q::from_integer(10.into())
                                                })
                                            })
                                            .collect(),
                                        &t,
                                    )
                                })
                                .collect(),
                            &u,
                        );
                        let hu = bezier(
                            net.windows(2)
                                .map(|r| {
                                    std::array::from_fn::<_, 4, _>(|axis| {
                                        let a = bezier(r[0].clone(), &t);
                                        let b = bezier(r[1].clone(), &t);
                                        &b[axis] - &a[axis]
                                    })
                                })
                                .collect(),
                            &u,
                        );
                        let expected = std::array::from_fn::<_, 3, _>(|axis| {
                            ((&hu[axis] * &tensor[3] - &tensor[axis] * &hu[3])
                                * &cap.derivative[0].0
                                + (&hv[axis] * &tensor[3] - &tensor[axis] * &hv[3]))
                                / (&tensor[3] * &tensor[3])
                        });
                        assert_eq!(d, expected);
                        assert!(v >= cap.uv[0][1].0 && v <= cap.uv[1][1].0);
                    }
                }
            }
        }
        eprintln!(
            "actual affine v-chart left={left} nb={nb} pieces={count} selection+support+field_work={}",
            w.used
        );
        assert!(count >= 36);
    }
}

#[test]
fn only_actual_constant_positive_weights_earn_the_affine_chart() {
    let mut w = Work::new(100_000, None).unwrap();
    let row = |weight| [N::integer(3), N::integer(1), N::integer(weight)];
    assert!(relation(&[row(1), row(2)], &mut w).unwrap().is_none());
    assert!(relation(&[row(0), row(0)], &mut w).unwrap().is_none());
    assert!(relation(&[row(-1), row(-1)], &mut w).unwrap().is_none());
    assert!(
        relation(&[row(1), row(1), row(1)], &mut w)
            .unwrap()
            .is_none()
    );
    assert!(relation(&[row(2), row(2)], &mut w).unwrap().is_some());
}

#[test]
fn finite_mixed_charts_earn_same_tensor_lift_under_one_actual_work_ledger() {
    use super::literal::{Precision, encode_cap_at, quantize_at};
    for (left, nb) in [(false, false), (true, true)] {
        let s = Hanger::new(crate::authored_hanger_bolt::tests::fixture(left, nb), None).unwrap();
        let b = s.boundary(BoundaryOptions::default(), None).unwrap();
        let f = EncodedField::prepare(&b, CurveOptions::default(), None).unwrap();
        let mut w = Work::new(16_000_000, None).unwrap();
        w.charge(f.work_used()).unwrap();
        let mut m = Memory::new(f.retained_bytes(), 64 * 1024 * 1024).unwrap();
        m.add(SCRATCH + size_of::<CapTrace>()).unwrap();
        let mut maximum = [0_f64; 3];
        let mut count = 0;
        let mut root_count = 0;
        let mut root_work = 0;
        for index in 0..f.cells().len() {
            for end in [End::Start, End::End] {
                // The retained original curve owns the source parameter and
                // orientation; the computational v chart never replaces it.
                let trace_memory = m.used;
                let Some(trace) = f
                    .cap_trace_with_ledgers(index, end, &mut m, &mut w)
                    .unwrap()
                else {
                    continue;
                };
                for piece in 0..f.angular.coefficients().len() {
                    let Some(cap) = build(&f, index, end, piece, &mut w).unwrap() else {
                        if f.meridians[f.cells[index].meridian].degree() == 2
                            && root_overlap(&f, index, end, piece, &trace, &mut w).unwrap()
                        {
                            // Exact overlap selection avoids performing an
                            // unreceipted empty standalone root query. The
                            // entire actual positive support/proof increment
                            // is debited before the next attempted piece.
                            let curve = f
                                .cap_curve_piece(
                                    index,
                                    end,
                                    piece,
                                    DomainOptions {
                                        maximum_work: w.remaining(),
                                        ..DomainOptions::default()
                                    },
                                    None,
                                )
                                .unwrap()
                                .unwrap();
                            let proof = super::trim_probe::run_with_limit(
                                &curve,
                                f.work_used() + w.remaining(),
                            )
                            .unwrap();
                            let increment = proof.work - f.work_used();
                            w.charge(increment).unwrap();
                            root_work += increment;
                            root_count += 1;
                            maximum[0] = maximum[0].max(proof.position);
                            maximum[1] = maximum[1].max(proof.derivative);
                        }
                        continue;
                    };
                    count += 1;
                    let before = m.used;
                    let chosen = cap
                        .trim
                        .each_ref()
                        .map(|t| quantize_at(t, Precision::Bits64, &mut w).unwrap());
                    let span = chosen[1].sub(&chosen[0], &mut w).unwrap();
                    let span_i = span.enclosure(&mut w).unwrap();
                    assert!(span_i.lo > 0.);
                    let xyz = encode_cap_at(
                        &cap.controls,
                        None,
                        CurveOptions::default(),
                        &mut m,
                        &mut w,
                        Precision::Bits64,
                    )
                    .unwrap();
                    let uv = chosen.each_ref().map(|t| {
                        let u = cap.uv[0][0]
                            .add(&cap.derivative[0].mul(t, &mut w).unwrap(), &mut w)
                            .unwrap();
                        let v = cap.uv[0][1]
                            .add(&cap.derivative[1].mul(t, &mut w).unwrap(), &mut w)
                            .unwrap();
                        [u, v, N::integer(0), N::integer(1)]
                    });
                    let encoded_uv = encode_cap_at(
                        &uv,
                        None,
                        CurveOptions::default(),
                        &mut m,
                        &mut w,
                        Precision::Bits64,
                    )
                    .unwrap();
                    let mut vrange = I::new(f64::INFINITY, f64::NEG_INFINITY);
                    for p in &encoded_uv.poles {
                        let u = p[0].number(&mut w).unwrap();
                        let v = p[1].number(&mut w).unwrap();
                        // Positive constant weights make the entire finite
                        // literal line stay inside this actual tensor chart.
                        assert!(u.0 >= Q::zero() && u.0 <= Q::one());
                        let ka = binary(f.angular.knots()[piece], &mut w).unwrap();
                        let kb = binary(f.angular.knots()[piece + 1], &mut w).unwrap();
                        assert!(v.0 >= ka.0 && v.0 <= kb.0);
                        let vi = v.enclosure(&mut w).unwrap();
                        vrange = I::new(vrange.lo.min(vi.lo), vrange.hi.max(vi.hi));
                    }
                    for point in &uv {
                        let vi = point[1].enclosure(&mut w).unwrap();
                        vrange = I::new(vrange.lo.min(vi.lo), vrange.hi.max(vi.hi));
                    }
                    let [j, h] =
                        super::trim_probe::field_bounds(&f, index, end, piece, vrange, &mut w)
                            .unwrap();
                    let gp = (p(cap.derivative[0].enclosure(&mut w).unwrap().magnitude())
                        + p(cap.derivative[1].enclosure(&mut w).unwrap().magnitude()))
                        * span_i;
                    let e0 = p(xyz.error[0]) + p(j) * p(encoded_uv.error[0]);
                    // Both terms are now derivatives in finite normalized
                    // lambda. Original-u and lambda errors are never mixed.
                    let e1 = p(xyz.error[1]) * span_i
                        + p(j) * p(encoded_uv.error[1])
                        + p(h) * p(encoded_uv.error[0]) * gp;
                    let phase0 = trace.endpoints[0].phase.number(&mut w).unwrap();
                    let phase1 = trace.endpoints[1].phase.number(&mut w).unwrap();
                    let original0 = trace.endpoints[0]
                        .original_parameter
                        .number(&mut w)
                        .unwrap();
                    let original1 = trace.endpoints[1]
                        .original_parameter
                        .number(&mut w)
                        .unwrap();
                    let row = &f.meridians[f.cells[index].meridian];
                    let q0 = row.phase[0]
                        .number(&mut w)
                        .unwrap()
                        .div(&row.weights[0].number(&mut w).unwrap(), &mut w)
                        .unwrap();
                    let q1 = row.phase[1]
                        .number(&mut w)
                        .unwrap()
                        .div(&row.weights[1].number(&mut w).unwrap(), &mut w)
                        .unwrap();
                    let gamma = q1
                        .sub(&q0, &mut w)
                        .unwrap()
                        .mul(&cap.derivative[0], &mut w)
                        .unwrap()
                        .mul(&original1.sub(&original0, &mut w).unwrap(), &mut w)
                        .unwrap()
                        .div(&phase1.sub(&phase0, &mut w).unwrap(), &mut w)
                        .unwrap();
                    assert!(gamma.0 < Q::zero());
                    let original_span = gamma.mul(&span, &mut w).unwrap();
                    let original_scale = original_span.enclosure(&mut w).unwrap();
                    assert!(original_scale.hi < 0.);
                    let original_e1 = e1 / -original_scale;
                    assert!(e0.valid() && e1.valid() && original_e1.valid());
                    assert!(e0.hi < 1e-7 && e1.hi < 1e-6 && original_e1.hi < 1e-6);
                    maximum[0] = maximum[0].max(e0.hi);
                    maximum[1] = maximum[1].max(e1.hi);
                    maximum[2] = maximum[2].max(original_e1.hi);
                    // Local prototype drops every encoded carrier before
                    // returning its retained allocation credit.
                    drop((xyz, encoded_uv));
                    m.used = before;
                }
                drop(trace);
                m.used = trace_memory;
            }
        }
        eprintln!(
            "finite mixed caps left={left} nb={nb} affine={count} root={root_count} root_increment={root_work} field+selection+trace+support+encoding+lift_work={} [C0all,D1finiteAll,D1originalAffine]={maximum:?}; streaming cap probe, excludes field serialization/closed material",
            w.used
        );
        assert!(count >= 36);
        assert_eq!(count + root_count, if nb { 48 } else { 44 });
    }
}

fn root_overlap(
    f: &EncodedField<'_>,
    index: usize,
    end: End,
    piece: usize,
    trace: &CapTrace,
    w: &mut Work<'_>,
) -> Result<bool, Error> {
    let cell = f.cells[index];
    let source = f.source.source_boundary().source();
    let spec = source.metric().spec();
    let pitch = binary(spec.pitch_mm, w)?;
    let lead = pitch.mul(&N::integer(i64::from(spec.starts)), w)?;
    let station = if end == End::Start {
        N::integer(0)
    } else {
        binary(spec.length_mm, w)?
    };
    let shift = station.sub(
        &pitch.mul(&N::integer(i64::from(cell.source_key.total_period)), w)?,
        w,
    )?;
    let a = binary(f.angular.knots()[piece], w)?;
    let b = binary(f.angular.knots()[piece + 1], w)?;
    let lower = shift.sub(&lead.mul(&b, w)?, w)?;
    let upper = shift.sub(&lead.mul(&a, w)?, w)?;
    let first = trace.endpoints[0].phase.number(w)?;
    let last = trace.endpoints[1].phase.number(w)?;
    let lo = if first.cmp(&lower, w)? == Ordering::Greater {
        first
    } else {
        lower
    };
    let hi = if last.cmp(&upper, w)? == Ordering::Less {
        last
    } else {
        upper
    };
    Ok(lo.cmp(&hi, w)? == Ordering::Less)
}

pub(super) fn bind_vertices(
    f: &EncodedField<'_>,
    index: usize,
    end: End,
    piece: usize,
    targets: [&[EncodedScalar; 3]; 2],
    m: &mut Memory,
    w: &mut Work<'_>,
) -> Result<[f64; 3], Error> {
    let cap = build(f, index, end, piece, w)?.ok_or(Error::Invalid)?;
    let proof = super::bound_affine::bind(f, (index, end, piece), targets, cap, m, w)?;
    let result = proof.bounds;
    m.release(&proof.xyz.poles);
    m.release(&proof.xyz.weights);
    m.release(&proof.uv.poles);
    m.release(&proof.uv.weights);
    drop(proof);
    Ok(result)
}
#[test]
fn actual_affine_and_root_overlap_arithmetic_cancel_at_every_attempted_checkpoint() {
    use std::cell::Cell;
    let s = Hanger::new(
        crate::authored_hanger_bolt::tests::fixture(false, false),
        None,
    )
    .unwrap();
    let b = s.boundary(BoundaryOptions::default(), None).unwrap();
    let f = EncodedField::prepare(&b, CurveOptions::default(), None).unwrap();
    let index = f
        .cells()
        .iter()
        .enumerate()
        .find_map(|(i, c)| {
            (f.meridians[c.meridian].degree() == 1 && c.source_key.total_period == -1).then_some(i)
        })
        .unwrap();
    let calls = Cell::new(0);
    let count = || {
        calls.set(calls.get() + 1);
        false
    };
    let mut w = Work::new(16_000_000, Some(&count)).unwrap();
    // Some angular pieces are empty; pick the first actual finite candidate.
    let piece = (0..f.angular.coefficients().len())
        .find(|&p| {
            build(
                &f,
                index,
                End::Start,
                p,
                &mut Work::new(16_000_000, None).unwrap(),
            )
            .unwrap()
            .is_some()
        })
        .unwrap();
    build(&f, index, End::Start, piece, &mut w)
        .unwrap()
        .unwrap();
    for stop in 1..=calls.get() {
        let n = Cell::new(0);
        let cancel = || {
            n.set(n.get() + 1);
            n.get() >= stop
        };
        let mut w = Work::new(16_000_000, Some(&cancel)).unwrap();
        assert_eq!(
            build(&f, index, End::Start, piece, &mut w).err(),
            Some(Error::Cancelled),
            "affine stop={stop}"
        );
    }
    let root = f
        .cells()
        .iter()
        .enumerate()
        .find_map(|(i, c)| {
            (f.meridians[c.meridian].degree() == 2 && c.source_key.total_period == -1).then_some(i)
        })
        .unwrap();
    let trace = f
        .cap_trace(root, End::Start, DomainOptions::default(), None)
        .unwrap()
        .unwrap();
    calls.set(0);
    let mut w = Work::new(16_000_000, Some(&count)).unwrap();
    root_overlap(&f, root, End::Start, 3, &trace, &mut w).unwrap();
    for stop in 1..=calls.get() {
        let n = Cell::new(0);
        let cancel = || {
            n.set(n.get() + 1);
            n.get() >= stop
        };
        let mut w = Work::new(16_000_000, Some(&cancel)).unwrap();
        assert_eq!(
            root_overlap(&f, root, End::Start, 3, &trace, &mut w).err(),
            Some(Error::Cancelled),
            "root overlap stop={stop}"
        );
    }
    let mut w = Work::new(16_000_000, None).unwrap();
    assert!(
        build(&f, root, End::Start, 3, &mut w).unwrap().is_none(),
        "actual retained rounded root must use its rational path"
    );
}

#[test]
fn actual_nonbinary_pitches_earn_inward_trims_and_every_bound_carrier_checkpoint() {
    use std::cell::Cell;
    let mut spec = crate::authored_hanger_bolt::tests::fixture(true, true);
    let metric = spec.metric.spec();
    spec.metric = AuthoredMetricThread::new(crate::authored_metric_thread::Spec {
        pitch_mm: 0.3,
        crest_width_mm: 0.03,
        ..metric
    })
    .unwrap();
    let wood = spec.wood.spec();
    spec.wood = crate::authored_flat_thread::AuthoredThread::new(
        crate::helical_thread::Spec {
            pitch_mm: 2.3,
            ..wood
        },
        crate::authored_flat_thread::SourceOptions {
            position_limit_mm: 1e-7,
            normal_angle_limit_rad: 1e-6,
        },
    )
    .unwrap();
    let source = Hanger::new(spec, None).unwrap();
    let boundary = source.boundary(BoundaryOptions::default(), None).unwrap();
    let f = EncodedField::prepare(&boundary, CurveOptions::default(), None).unwrap();
    let mut selected = None;
    let mut w = Work::new(16_000_000, None).unwrap();
    w.charge(f.work_used()).unwrap();
    let mut memory = Memory::new(f.retained_bytes(), 64 * 1024 * 1024).unwrap();
    memory.add(SCRATCH).unwrap();
    let mut inward = [false; 2];
    for index in 0..f.cells().len() {
        for end in [End::Start, End::End] {
            for piece in 0..f.angular.coefficients().len() {
                let Some(cap) = build(&f, index, end, piece, &mut w).unwrap() else {
                    continue;
                };
                for (e, changed) in inward.iter_mut().enumerate() {
                    let rounded = super::literal::quantize_at(
                        &cap.trim[e],
                        super::literal::Precision::Bits64,
                        &mut w,
                    )
                    .unwrap();
                    *changed |= rounded.cmp(&cap.trim[e], &mut w).unwrap()
                        == if e == 0 {
                            Ordering::Less
                        } else {
                            Ordering::Greater
                        };
                }
                let mut local = Work::new(16_000_000, None).unwrap();
                let targets = cap.trim.each_ref().map(|t| {
                    let h = super::trim_probe::eval(&cap.controls, t, &mut local).unwrap();
                    std::array::from_fn::<_, 3, _>(|a| {
                        EncodedScalar::retain(
                            &h[a].div(&h[3], &mut local).unwrap(),
                            &mut memory,
                            &mut local,
                        )
                        .unwrap()
                    })
                });
                // The public owner order remains increasing source phase.
                let refs = [&targets[1], &targets[0]];
                bind_vertices(&f, index, end, piece, refs, &mut memory, &mut w).unwrap();
                if selected.is_none() {
                    selected = Some((index, end, piece, targets));
                }
                if inward.into_iter().all(|x| x) {
                    break;
                }
            }
            if inward.into_iter().all(|x| x) {
                break;
            }
        }
        if inward.into_iter().all(|x| x) {
            break;
        }
    }
    assert!(inward.into_iter().all(|x| x));
    let (index, end, piece, targets) = selected.unwrap();
    let refs = [&targets[1], &targets[0]];
    let calls = Cell::new(0);
    let count = || {
        calls.set(calls.get() + 1);
        false
    };
    let mut w = Work::new(16_000_000, Some(&count)).unwrap();
    let mut m = Memory::new(f.retained_bytes(), 64 * 1024 * 1024).unwrap();
    m.add(SCRATCH).unwrap();
    bind_vertices(&f, index, end, piece, refs, &mut m, &mut w).unwrap();
    for stop in 1..=calls.get() {
        let n = Cell::new(0);
        let cancel = || {
            n.set(n.get() + 1);
            n.get() >= stop
        };
        let mut w = Work::new(16_000_000, Some(&cancel)).unwrap();
        let mut m = Memory::new(f.retained_bytes(), 64 * 1024 * 1024).unwrap();
        m.add(SCRATCH).unwrap();
        assert_eq!(
            bind_vertices(&f, index, end, piece, refs, &mut m, &mut w).err(),
            Some(Error::Cancelled),
            "bound carrier stop={stop}"
        );
    }
}
#[test]
fn source_owned_affine_support_api_retains_exact_controls_and_refusals() {
    for (left, nb) in [(false, false), (true, true)] {
        let source =
            Hanger::new(crate::authored_hanger_bolt::tests::fixture(left, nb), None).unwrap();
        let boundary = source.boundary(BoundaryOptions::default(), None).unwrap();
        let field = EncodedField::prepare(&boundary, CurveOptions::default(), None).unwrap();
        assert!(matches!(
            field.affine_cap_piece(
                0,
                End::Start,
                0,
                DomainOptions {
                    maximum_bytes: field.retained_bytes() + 1024 * 1024,
                    ..Default::default()
                },
                None
            ),
            Err(Error::Limit)
        ));
        let scan_work: usize = field
            .meridians
            .iter()
            .flat_map(|row| row.radius.iter().chain(&row.phase).chain(&row.weights))
            .map(|s| 1 + s.n.len() + s.d.len())
            .sum();
        assert!(matches!(
            field.affine_cap_piece(
                0,
                End::Start,
                0,
                DomainOptions {
                    maximum_work: scan_work + 2,
                    ..Default::default()
                },
                None
            ),
            Err(Error::Limit)
        ));
        let mut positive = 0;
        let mut absent = 0;
        for cell in 0..field.cells.len() {
            for end in [End::Start, End::End] {
                for piece in 0..field.angular.coefficients().len() {
                    let options = DomainOptions::default();
                    let result = field
                        .affine_cap_piece(cell, end, piece, options, None)
                        .unwrap();
                    let mut work = Work::new(16_000_000, None).unwrap();
                    let raw = build(&field, cell, end, piece, &mut work).unwrap();
                    match (result, raw) {
                        (Some(literal), Some(raw)) => {
                            positive += 1;
                            assert!(std::ptr::eq(literal.field(), &field));
                            assert_eq!(literal.source_cell(), cell);
                            assert_eq!(literal.end(), end);
                            assert_eq!(literal.angular_piece(), piece);
                            assert!(literal.original_phase_reversed());
                            assert!(literal.retained_bytes() > 0);
                            for (a, b) in literal.controls().iter().zip(&raw.controls) {
                                for (axis, q) in a.iter().zip(b) {
                                    assert_eq!(axis.number(&mut work).unwrap().0, q.0);
                                }
                            }
                            for i in 0..2 {
                                assert_eq!(
                                    literal.trim()[i].number(&mut work).unwrap().0,
                                    raw.trim[i].0
                                );
                                assert_eq!(
                                    literal.chart_derivative()[i].number(&mut work).unwrap().0,
                                    raw.derivative[i].0
                                );
                                for axis in 0..2 {
                                    assert_eq!(
                                        literal.surface_endpoints()[i][axis]
                                            .number(&mut work)
                                            .unwrap()
                                            .0,
                                        raw.uv[i][axis].0
                                    );
                                }
                            }
                            let used = literal.work_used();
                            assert!(
                                field
                                    .affine_cap_piece(
                                        cell,
                                        end,
                                        piece,
                                        DomainOptions {
                                            maximum_work: used,
                                            ..options
                                        },
                                        None
                                    )
                                    .unwrap()
                                    .is_some()
                            );
                            assert!(matches!(
                                field.affine_cap_piece(
                                    cell,
                                    end,
                                    piece,
                                    DomainOptions {
                                        maximum_work: used - 1,
                                        ..options
                                    },
                                    None
                                ),
                                Err(Error::Limit)
                            ));
                        }
                        (None, None) => {
                            absent += 1;
                        }
                        _ => panic!("source carrier differs"),
                    }
                }
            }
        }
        assert!(positive > 0 && absent > 0);
        assert!(matches!(
            field.affine_cap_piece(usize::MAX, End::Start, 0, DomainOptions::default(), None),
            Err(Error::Invalid)
        ));
        assert!(matches!(
            field.affine_cap_piece(0, End::Start, usize::MAX, DomainOptions::default(), None),
            Err(Error::Invalid)
        ));
        assert!(matches!(
            field.affine_cap_piece(
                0,
                End::Start,
                0,
                DomainOptions {
                    maximum_bytes: field.retained_bytes(),
                    ..Default::default()
                },
                None
            ),
            Err(Error::Limit)
        ));
        assert!(matches!(
            field.affine_cap_piece(0, End::Start, 0, DomainOptions::default(), Some(&|| true)),
            Err(Error::Cancelled)
        ));
    }
}
