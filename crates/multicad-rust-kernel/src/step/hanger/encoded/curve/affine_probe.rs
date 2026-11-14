//! Exact affine-row angular chart prototype, never a rounded-root alias.
use super::*;
use crate::authored_hanger_bolt::{BoundaryOptions, Hanger};
use crate::step::hanger::encoded::cap::tests::{bezier, number};
use num_traits::{One, Zero};
type Q = BigRational;
struct AffineCap {
    controls: Vec<[N; 4]>,
    trim: [N; 2],
    uv: [[N; 2]; 2],
    // Original meridian-u and angular-v derivatives in basis t.
    derivative: [N; 2],
    reversed: bool,
}
fn relation(h: &[[N; 3]], w: &mut Work<'_>) -> Result<Option<[[N; 2]; 2]>, Error> {
    if h.len() != 2
        || h[0][2].cmp(&h[1][2], w)? != Ordering::Equal
        || h[0][2].cmp(&N::integer(0), w)? != Ordering::Greater
    {
        return Ok(None);
    }
    Ok(Some([
        [h[0][0].div(&h[0][2], w)?, h[0][1].div(&h[0][2], w)?],
        [h[1][0].div(&h[1][2], w)?, h[1][1].div(&h[1][2], w)?],
    ]))
}
fn build(
    field: &EncodedField<'_>,
    index: usize,
    end: End,
    piece: usize,
    w: &mut Work<'_>,
) -> Result<Option<AffineCap>, Error> {
    w.charge(1)?;
    let cell = field.cells[index];
    let row = &field.meridians[cell.meridian];
    if row.degree() != 1 {
        return Ok(None);
    }
    let input = (0..=1)
        .map(|i| {
            Ok([
                row.radius[i].number(w)?,
                row.phase[i].number(w)?,
                row.weights[i].number(w)?,
            ])
        })
        .collect::<Result<Vec<_>, Error>>()?;
    let Some([[r0, q0], [r1, q1]]) = relation(&input, w)? else {
        return Ok(None);
    };
    let source = field.source.source_boundary().source();
    let (pitch, starts, left, length) = if row.component == Component::Metric {
        let s = source.metric().spec();
        (s.pitch_mm, s.starts, s.left_hand, s.length_mm)
    } else {
        let s = source.wood().spec();
        (s.pitch_mm, s.starts, s.left_hand, s.length_mm)
    };
    let pitch = binary(pitch, w)?;
    let lead = pitch.mul(&N::integer(i64::from(starts)), w)?;
    let station = if end == End::Start {
        N::integer(0)
    } else {
        binary(length, w)?
    };
    let shift = station.sub(
        &pitch.mul(&N::integer(i64::from(cell.source_key.total_period)), w)?,
        w,
    )?;
    let a = binary(field.angular.knots()[piece], w)?;
    let b = binary(field.angular.knots()[piece + 1], w)?;
    let from = shift.sub(&q1, w)?.div(&lead, w)?;
    let to = shift.sub(&q0, w)?.div(&lead, w)?;
    let lo = if from.cmp(&a, w)? == Ordering::Greater {
        from
    } else {
        a.clone()
    };
    let hi = if to.cmp(&b, w)? == Ordering::Less {
        to
    } else {
        b.clone()
    };
    if lo.cmp(&hi, w)? != Ordering::Less {
        return Ok(None);
    }
    let width = b.sub(&a, w)?;
    let qspan = q1.sub(&q0, w)?;
    let ua = shift
        .sub(&lead.mul(&a, w)?, w)?
        .sub(&q0, w)?
        .div(&qspan, w)?;
    let ub = shift
        .sub(&lead.mul(&b, w)?, w)?
        .sub(&q0, w)?
        .div(&qspan, w)?;
    let ra = r0.add(&r1.sub(&r0, w)?.mul(&ua, w)?, w)?;
    let rb = r0.add(&r1.sub(&r0, w)?.mul(&ub, w)?, w)?;
    let ap = if left {
        field.angular.coefficients().len() - 1 - piece
    } else {
        piece
    };
    let coefficients = &field.angular.coefficients()[ap];
    let mut z = station;
    if row.component == Component::Metric {
        z = z
            .sub(&binary(source.metric().full_length_mm(), w)?, w)?
            .sub(&binary(source.spec().plain_length_mm, w)?, w)?;
    }
    let mut controls = Vec::with_capacity(12);
    for i in 0usize..=11 {
        let mut xy = [N::integer(0), N::integer(0)];
        for (j, r, scale) in [(i, &ra, 11 - i), (i.wrapping_sub(1), &rb, i)] {
            if scale == 0 {
                continue;
            }
            let j = if left { 10 - j } else { j };
            for (out, c) in xy.iter_mut().zip(&coefficients[j]) {
                w.charge(512)?;
                *out = out.add(
                    &r.mul(&N(c.clone()), w)?.mul(&N::integer(scale as i64), w)?,
                    w,
                )?;
            }
        }
        controls.push([
            xy[0].clone(),
            xy[1].clone(),
            z.mul(&N::integer(11), w)?,
            N::integer(11),
        ]);
    }
    let trim = [
        lo.sub(&a, w)?.div(&width, w)?,
        hi.sub(&a, w)?.div(&width, w)?,
    ];
    let du = ub.sub(&ua, w)?;
    Ok(Some(AffineCap {
        controls,
        trim,
        uv: [[ua, a], [ub, b]],
        derivative: [du, width],
        reversed: true,
    }))
}
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
