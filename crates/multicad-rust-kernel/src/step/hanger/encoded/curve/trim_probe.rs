//! Bounded actual-source prototype. No file/native admission.
use super::correction::correct;
use super::literal::{Precision, encode_cap_at, quantize_at};
use super::restrict::{chosen, derivatives, restrict};
use super::*;
use crate::authored_hanger_bolt::{BoundaryOptions, Hanger};
fn encode64(
    h: &[[N; 4]],
    endpoints: Option<&[[EncodedScalar; 3]; 2]>,
    options: CurveOptions,
    m: &mut Memory,
    w: &mut Work<'_>,
) -> Result<super::literal::CartesianCap, Error> {
    encode_cap_at(h, endpoints, options, m, w, Precision::Bits64)
}
fn quantize64(value: &N, w: &mut Work<'_>) -> Result<N, Error> {
    quantize_at(value, Precision::Bits64, w)
}
pub(super) fn eval<const D: usize>(
    input: &[[N; D]],
    t: &N,
    w: &mut Work<'_>,
) -> Result<[N; D], Error> {
    w.charge(input.len() * D * 512)?;
    let mut values = input.to_vec();
    let one = N::integer(1).sub(t, w)?;
    for n in (1..values.len()).rev() {
        for i in 0..n {
            let (lo, hi) = values.split_at_mut(i + 1);
            for (a, b) in lo[i].iter_mut().zip(&hi[0]) {
                *a = a.mul(&one, w)?.add(&b.mul(t, w)?, w)?;
            }
        }
    }
    Ok(values.remove(0))
}
pub(super) use super::jets::eval_interval;
fn decode<const D: usize>(
    c: &[[EncodedScalar; D]],
    w: &mut Work<'_>,
) -> Result<Vec<[N; D]>, Error> {
    c.iter()
        .map(|row| {
            let mut result = std::array::from_fn(|_| N::integer(0));
            for (out, x) in result.iter_mut().zip(row) {
                *out = x.number(w)?;
            }
            Ok(result)
        })
        .collect()
}
pub(super) struct Probe {
    pub(super) work: usize,
    pub(super) position: f64,
    pub(super) derivative: f64,
    pub(super) source_position: f64,
    pub(super) source_derivative: f64,
    pub(super) original_source_derivative: f64,
    vertex_error: f64,
    parameter: [N; 2],
    stages: [usize; 7],
    xyz: super::literal::CartesianCap,
    uv: super::literal::CartesianCap,
}
fn run(curve: &CapCurvePiece<'_, '_>) -> Result<Probe, Error> {
    run_with_limit(curve, 16_000_000)
}
pub(super) fn run_with_limit(curve: &CapCurvePiece<'_, '_>, limit: usize) -> Result<Probe, Error> {
    let proof = run_bound(curve, None, limit, None, false)?;
    source_gate(&proof)?;
    Ok(proof)
}
pub(super) fn run_to_vertices(
    curve: &CapCurvePiece<'_, '_>,
    targets: [&[EncodedScalar; 3]; 2],
    limit: usize,
) -> Result<Probe, Error> {
    let proof = run_bound(curve, Some(targets), limit, None, false)?;
    source_gate(&proof)?;
    Ok(proof)
}
fn parameter_envelope(
    bounds: [[f64; 2]; 2],
    chosen: &[N; 2],
    w: &mut Work<'_>,
) -> Result<[N; 2], Error> {
    let low = binary(bounds[0][0], w)?;
    let high = binary(bounds[1][1], w)?;
    Ok([
        if low.cmp(&chosen[0], w)? == Ordering::Greater {
            chosen[0].clone()
        } else {
            low
        },
        if high.cmp(&chosen[1], w)? == Ordering::Less {
            chosen[1].clone()
        } else {
            high
        },
    ])
}
fn source_gate(proof: &Probe) -> Result<(), Error> {
    source_widths([
        proof.source_position,
        proof.source_derivative,
        proof.original_source_derivative,
    ])
}
fn source_widths(widths: [f64; 3]) -> Result<(), Error> {
    if widths
        .into_iter()
        .zip([1e-7, 1e-6, 1e-6])
        .any(|(v, limit)| !v.is_finite() || v < 0. || v > limit)
    {
        Err(Error::Unresolved)
    } else {
        Ok(())
    }
}

fn run_bound(
    curve: &CapCurvePiece<'_, '_>,
    targets: Option<[&[EncodedScalar; 3]; 2]>,
    limit: usize,
    interrupt: Option<&dyn Fn() -> bool>,
    exact_endpoints: bool,
) -> Result<Probe, Error> {
    let f = curve.field;
    let cell = f.cells[curve.source_cell];
    let row = &f.meridians[cell.meridian];
    let s = f.source.source_boundary().source();
    let (pitch, starts, left, length) = if row.component == Component::Metric {
        let v = s.metric().spec();
        (v.pitch_mm, v.starts, v.left_hand, v.length_mm)
    } else {
        let v = s.wood().spec();
        (v.pitch_mm, v.starts, v.left_hand, v.length_mm)
    };
    let mut w = Work::new(limit, interrupt)?;
    w.charge(f.work_used() + curve.work_used)?;
    let mut m = Memory::new(f.retained_bytes() + curve.retained_bytes, 64 * 1024 * 1024)?;
    m.add(SCRATCH)?;
    let pitch = binary(pitch, &mut w)?;
    let lead = pitch.mul(&N::integer(i64::from(starts)), &mut w)?;
    let station = if curve.end == End::Start {
        N::integer(0)
    } else {
        binary(length, &mut w)?
    };
    let shift = station.sub(
        &pitch.mul(&N::integer(i64::from(cell.source_key.total_period)), &mut w)?,
        &mut w,
    )?;
    let plane = curve.surface_trace.plane_z.number(&mut w)?;
    let ka = binary(f.angular.knots()[curve.angular_piece], &mut w)?;
    let kb = binary(f.angular.knots()[curve.angular_piece + 1], &mut w)?;
    let kw = kb.sub(&ka, &mut w)?;
    let angular_piece = if left {
        f.angular.coefficients().len() - 1 - curve.angular_piece
    } else {
        curve.angular_piece
    };
    let angular = (0..=10)
        .map(|i| {
            f.angular.coefficients()[angular_piece][if left { 10 - i } else { i }]
                .clone()
                .map(N)
        })
        .collect::<Vec<_>>();
    let angular_interval = angular
        .iter()
        .map(|c| Ok([c[0].enclosure(&mut w)?, c[1].enclosure(&mut w)?]))
        .collect::<Result<Vec<_>, Error>>()?;
    let meridian = (0..=row.degree())
        .map(|i| {
            Ok([
                row.radius[i].number(&mut w)?,
                row.phase[i].number(&mut w)?,
                row.weights[i].number(&mut w)?,
            ])
        })
        .collect::<Result<Vec<_>, Error>>()?;
    let baseline = w.used;
    let mut stages = [0usize; 7];
    let mut parameter: [N; 2] = std::array::from_fn(|_| N::integer(0));
    let mut vtarget: [N; 2] = std::array::from_fn(|_| N::integer(0));
    let mut points: [[N; 3]; 2] = std::array::from_fn(|_| std::array::from_fn(|_| N::integer(0)));
    let mut delta: [[N; 3]; 2] = std::array::from_fn(|_| std::array::from_fn(|_| N::integer(0)));
    let mut vertex_error = 0_f64;
    let mut uv_delta: [N; 2] = std::array::from_fn(|_| N::integer(0));
    let mut parameter_error = [0_f64; 2];
    for e in 0..2 {
        let phase = curve.phase[e].number(&mut w)?;
        let (selected, error) = chosen(row, &phase, curve.parameter_bounds[e], &pitch, &mut w)?;
        parameter[e] = quantize64(&selected, &mut w)?;
        parameter_error[e] = (p(error)
            + p(parameter[e]
                .sub(&selected, &mut w)?
                .enclosure(&mut w)?
                .magnitude()))
        .hi;
        vtarget[e] = quantize64(&shift.sub(&phase, &mut w)?.div(&lead, &mut w)?, &mut w)?;
        assert!(
            vtarget[e].cmp(&ka, &mut w)? != Ordering::Less
                && vtarget[e].cmp(&kb, &mut w)? != Ordering::Greater
        );
        let h = eval(&meridian, &parameter[e], &mut w)?;
        let radius = h[0].div(&h[2], &mut w)?;
        let q = h[1].div(&h[2], &mut w)?;
        let vold = shift.sub(&q, &mut w)?.div(&lead, &mut w)?;
        let old_t = vold.sub(&ka, &mut w)?.div(&kw, &mut w)?;
        let target_t = vtarget[e].sub(&ka, &mut w)?.div(&kw, &mut w)?;
        let old_i = eval_interval(&angular_interval, old_t.enclosure(&mut w)?, &mut w)?;
        let target_i = eval_interval(&angular_interval, target_t.enclosure(&mut w)?, &mut w)?;
        let radius_i = radius.enclosure(&mut w)?;
        let mut residual = p(0.);
        for axis in 0..2 {
            let point = radius_i * target_i[axis];
            points[e][axis] = if let Some(targets) = targets {
                targets[e][axis].number(&mut w)?
            } else {
                quantize64(&binary(point.midpoint(), &mut w)?, &mut w)?
            };
            let required = points[e][axis].enclosure(&mut w)? - radius_i * old_i[axis];
            delta[e][axis] = quantize64(&binary(required.midpoint(), &mut w)?, &mut w)?;
            residual = residual + p((required - delta[e][axis].enclosure(&mut w)?).magnitude());
        }
        if exact_endpoints || !residual.valid() || residual.hi > 1e-7 {
            // Bounded precision fallback, with all failed proposal work retained.
            let old = eval(&angular, &old_t, &mut w)?;
            let target = eval(&angular, &target_t, &mut w)?;
            residual = p(0.);
            for axis in 0..2 {
                points[e][axis] = if let Some(targets) = targets {
                    targets[e][axis].number(&mut w)?
                } else {
                    quantize64(&radius.mul(&target[axis], &mut w)?, &mut w)?
                };
                let required = points[e][axis].sub(&radius.mul(&old[axis], &mut w)?, &mut w)?;
                delta[e][axis] = quantize64(&required, &mut w)?;
                residual = residual
                    + p(required
                        .sub(&delta[e][axis], &mut w)?
                        .enclosure(&mut w)?
                        .magnitude());
            }
        }
        vertex_error = vertex_error.max(residual.hi);
        if let Some(targets) = targets {
            assert_eq!(targets[e][2].number(&mut w)?.0, plane.0);
        }
        points[e][2] = plane.clone();
        w.charge(512)?;
        uv_delta[e] = vtarget[e].sub(&vold, &mut w)?;
    }
    stages[0] = w.used - baseline;
    let span = parameter[1].sub(&parameter[0], &mut w)?;
    let span_bound = span.enclosure(&mut w)?;
    assert!(span_bound.lo > 0.);
    let h = decode(curve.controls(), &mut w)?;
    let corrected = correct(&h, &delta, &parameter[0], &parameter[1], &mut m, &mut w)?;
    stages[1] = w.used - baseline - stages[..1].iter().sum::<usize>();
    let xyz = encode64(
        &corrected.controls,
        None,
        CurveOptions::default(),
        &mut m,
        &mut w,
    )?;
    stages[2] = w.used - baseline - stages[..2].iter().sum::<usize>();
    vertex_error = (p(vertex_error) + p(xyz.error[0])).hi;
    assert!(vertex_error <= 1e-7);
    let finite = restrict(&meridian, &parameter[0], &parameter[1], &mut m, &mut w)?;
    let degree = row.degree() + 1;
    let mut uv = Vec::new();
    for i in 0..=degree {
        let mut den = N::integer(0);
        let mut num = N::integer(0);
        let mut linear = N::integer(0);
        for (j, scale, e) in [(i, degree - i, 0), (i.wrapping_sub(1), i, 1)] {
            if scale == 0 {
                continue;
            }
            let d = lead.mul(&finite[j][2], &mut w)?;
            let v = shift
                .mul(&finite[j][2], &mut w)?
                .sub(&finite[j][1], &mut w)?;
            den = den.add(&d.mul(&N::integer(scale as i64), &mut w)?, &mut w)?;
            num = num.add(
                &v.add(&uv_delta[e].mul(&d, &mut w)?, &mut w)?
                    .mul(&N::integer(scale as i64), &mut w)?,
                &mut w,
            )?;
            if e == 1 {
                linear = d.mul(&N::integer(i as i64), &mut w)?;
            }
        }
        let u = parameter[0]
            .mul(&den, &mut w)?
            .add(&span.mul(&linear, &mut w)?, &mut w)?;
        uv.push([u, num, N::integer(0), den]);
    }
    let mut encoded_endpoints: [[EncodedScalar; 3]; 2] = std::array::from_fn(|_| {
        std::array::from_fn(|_| EncodedScalar {
            sign: 0,
            n: Box::new([]),
            d: Box::new([1]),
        })
    });
    for e in 0..2 {
        encoded_endpoints[e] = [
            EncodedScalar::retain(&parameter[e], &mut m, &mut w)?,
            EncodedScalar::retain(&vtarget[e], &mut m, &mut w)?,
            EncodedScalar::retain(&N::integer(0), &mut m, &mut w)?,
        ];
    }
    stages[3] = w.used - baseline - stages[..3].iter().sum::<usize>();
    let uv_literal = encode64(
        &uv,
        Some(&encoded_endpoints),
        CurveOptions::default(),
        &mut m,
        &mut w,
    )?;
    stages[4] = w.used - baseline - stages[..4].iter().sum::<usize>();
    for point in &uv_literal.poles {
        let u = point[0].number(&mut w)?;
        let v = point[1].number(&mut w)?;
        assert!(
            u.cmp(&N::integer(0), &mut w)? != Ordering::Less
                && u.cmp(&N::integer(1), &mut w)? != Ordering::Greater
        );
        assert!(v.cmp(&ka, &mut w)? != Ordering::Less && v.cmp(&kb, &mut w)? != Ordering::Greater);
    }
    // Whole-trim lift proof is supplied below by the actual tensor Jacobian
    // and Hessian coefficient bounds, not by the endpoint correspondence.
    let original_uv = decode(curve.surface_trace.controls(), &mut w)?;
    let original_uv = restrict(&original_uv, &parameter[0], &parameter[1], &mut m, &mut w)?;
    let orig4 = original_uv
        .iter()
        .map(|p| [p[0].clone(), p[1].clone(), N::integer(0), p[2].clone()])
        .collect::<Vec<_>>();
    let gp = derivatives(&orig4, &N::integer(1), &mut w)?[0];
    let mut vrange = I::new(f64::INFINITY, f64::NEG_INFINITY);
    for point in &original_uv {
        let v = point[1].div(&point[2], &mut w)?.enclosure(&mut w)?;
        vrange = I::new(vrange.lo.min(v.lo), vrange.hi.max(v.hi));
    }
    for point in &uv_literal.poles {
        let v = point[1].number(&mut w)?.enclosure(&mut w)?;
        vrange = I::new(vrange.lo.min(v.lo), vrange.hi.max(v.hi));
    }
    stages[5] = w.used - baseline - stages[..5].iter().sum::<usize>();
    let [j, hessian] = tensor_bounds(curve, vrange, &mut w)?;
    let duv = p(uv_delta[0]
        .enclosure(&mut w)?
        .magnitude()
        .max(uv_delta[1].enclosure(&mut w)?.magnitude()))
        + p(uv_literal.error[0]);
    let duvp = p(uv_delta[1]
        .sub(&uv_delta[0], &mut w)?
        .enclosure(&mut w)?
        .magnitude())
        + p(uv_literal.error[1]);
    let d3 = corrected.position_error + p(xyz.error[0]);
    let d3p = (corrected.derivative_error + p(xyz.error[1])) * span_bound;
    let lift = d3 + p(j) * duv;
    let lift_prime = d3p + p(j) * duvp + p(hessian) * duv * p(gp);
    assert!(lift.valid() && lift_prime.valid());
    assert!(
        lift.hi <= 1e-7 && lift_prime.hi <= 1e-6,
        "lift={lift:?} derivative={lift_prime:?}"
    );
    // True-event transfer: chosen finite parameters are proposals, not exact
    // algebraic endpoints. Compare the two finite curves under their actual
    // affine-in-root-u lambda correspondence before claiming source accuracy.
    let [low, high] = parameter_envelope(curve.parameter_bounds, &parameter, &mut w)?;
    let jets = super::restrict::factor_derivatives(curve, &low, &high, &mut m, &mut w)?;
    let eps = p(parameter_error[0].max(parameter_error[1]));
    let event_position = p(jets[0]) * eps;
    let event_derivative =
        p(jets[1]) * eps + p(jets[0]) * (p(parameter_error[0]) + p(parameter_error[1]));
    let true_span = I::new(curve.parameter_bounds[1][0], curve.parameter_bounds[1][1])
        - I::new(curve.parameter_bounds[0][0], curve.parameter_bounds[0][1]);
    assert!(true_span.lo > 0.);
    w.charge(512)?;
    let dq = if row.degree() == 2 {
        let b = s.metric().root_meridian().whole_bounds();
        I::new(b.first_derivative_mm[0][0], b.first_derivative_mm[0][1])
    } else {
        meridian[row.degree()][1]
            .div(&meridian[row.degree()][2], &mut w)?
            .sub(&meridian[0][1].div(&meridian[0][2], &mut w)?, &mut w)?
            .enclosure(&mut w)?
    };
    let source_position = lift + event_position + p(row.bounds.position_error_mm);
    let source_derivative = lift_prime
        + event_derivative
        + (p(row.bounds.derivative_error_mm[0])
            + p(row.bounds.derivative_error_mm[1]) * p(dq.magnitude()) / lead.enclosure(&mut w)?)
            * p(true_span.hi.min(1.));
    let original_span = curve.surface_trace.endpoints[1]
        .original_parameter
        .number(&mut w)?
        .sub(
            &curve.surface_trace.endpoints[0]
                .original_parameter
                .number(&mut w)?,
            &mut w,
        )?;
    let phase_span = curve.surface_trace.endpoints[1].phase.number(&mut w)?.sub(
        &curve.surface_trace.endpoints[0].phase.number(&mut w)?,
        &mut w,
    )?;
    let original_scale =
        original_span.div(&phase_span, &mut w)?.enclosure(&mut w)? * dq * true_span;
    assert!(original_scale.lo > 0.);
    let original_source_derivative = source_derivative / original_scale;

    w.charge(1)?;
    stages[6] = w.used - baseline - stages[..6].iter().sum::<usize>();
    // Independent exact-rational point/D1 oracle is added by the caller test.
    Ok(Probe {
        work: w.used,
        position: lift.hi,
        derivative: lift_prime.hi,
        source_position: source_position.hi,
        source_derivative: source_derivative.hi,
        original_source_derivative: original_source_derivative.hi,
        vertex_error,
        parameter,
        stages,
        xyz,
        uv: uv_literal,
    })
}
fn tensor_bounds(
    curve: &CapCurvePiece<'_, '_>,
    vrange: I,
    w: &mut Work<'_>,
) -> Result<[f64; 2], Error> {
    field_bounds(
        curve.field,
        curve.source_cell,
        curve.end,
        curve.angular_piece,
        vrange,
        w,
    )
}
pub(super) use super::jets::field_bounds;
#[test]
fn actual_root_full_basis_and_same_tensor_pcurve_lift_pass_both_hands_nonbinary() {
    for (left, nb, period, end, piece) in [
        (false, false, -1, End::Start, 3),
        (true, true, 25, End::End, 0),
    ] {
        let s = Hanger::new(crate::authored_hanger_bolt::tests::fixture(left, nb), None).unwrap();
        let b = s.boundary(BoundaryOptions::default(), None).unwrap();
        let f = EncodedField::prepare(&b, CurveOptions::default(), None).unwrap();
        let i = f
            .cells()
            .iter()
            .position(|c| {
                c.source_key.total_period == period
                    && matches!(
                        c.source_key.face,
                        FaceKey::MetricThread {
                            band: ThreadBandRole::Root,
                            ..
                        }
                    )
            })
            .unwrap();
        let curve = f
            .cap_curve_piece(i, end, piece, DomainOptions::default(), None)
            .unwrap()
            .unwrap();
        let proof = run(&curve).unwrap();
        independent_oracle(&curve, &proof);
        let targets = std::array::from_fn(|e| {
            let mut memory = Memory::new(0, 64 * 1024 * 1024).unwrap();
            let mut work = Work::new(16_000_000, None).unwrap();
            let u = &proof.parameter[e];
            let xyz = eval(&decode(curve.controls(), &mut work).unwrap(), u, &mut work).unwrap();
            std::array::from_fn(|axis| {
                EncodedScalar::retain(
                    &xyz[axis].div(&xyz[3], &mut work).unwrap(),
                    &mut memory,
                    &mut work,
                )
                .unwrap()
            })
        });
        let bound = run_to_vertices(&curve, targets.each_ref(), 16_000_000).unwrap();
        assert!(bound.vertex_error < 1e-7);
        eprintln!(
            "actual full-basis+UV left={left} nb={nb} cumulative_work={} position={} derivative={}",
            proof.work, proof.position, proof.derivative
        );
    }
}

#[test]
fn old_meridian_u_chart_retains_actual_field_lift_and_real_source_c1_refusal() {
    for (left, nb) in [(false, false), (true, true)] {
        let source =
            Hanger::new(crate::authored_hanger_bolt::tests::fixture(left, nb), None).unwrap();
        let boundary = source.boundary(BoundaryOptions::default(), None).unwrap();
        let field = EncodedField::prepare(&boundary, CurveOptions::default(), None).unwrap();
        let mut count = 0;
        let mut subtotal = field.work_used();
        let mut worst = [0_f64; 2];
        let mut stages = [0usize; 7];
        let mut vertex = 0_f64;
        let mut source_refusals = 0;
        for index in 0..field.cells().len() {
            for end in [End::Start, End::End] {
                for piece in 0..field.angular.coefficients().len() {
                    let Some(curve) = field
                        .cap_curve_piece(index, end, piece, DomainOptions::default(), None)
                        .unwrap()
                    else {
                        continue;
                    };
                    // The historical chart is a paired diagnostic, not the
                    // new Source-qualified carrier route. Preserve its real
                    // tiny-span C1 refusal while checking the field proof.
                    let proof = run_bound(&curve, None, 16_000_000, None, false).unwrap();
                    if source_gate(&proof).is_err() {
                        source_refusals += 1;
                        assert!(nb);
                        assert!(matches!(
                            field.cells[index].source_key.face,
                            FaceKey::WoodThread {
                                band: ThreadBandRole::Crest,
                                ..
                            }
                        ));
                        eprintln!(
                            "old-u source C1 refusal left={left} nb={nb} cell={index} end={end:?} piece={piece} originalD1={}",
                            proof.original_source_derivative
                        );
                    }
                    subtotal += proof.work - field.work_used();
                    for (sum, value) in stages.iter_mut().zip(proof.stages) {
                        *sum += value;
                    }
                    vertex = vertex.max(proof.vertex_error);
                    worst[0] = worst[0].max(proof.position);
                    worst[1] = worst[1].max(proof.derivative);
                    count += 1;
                }
            }
        }
        eprintln!(
            "old-u field-only source_refusals={source_refusals} caps-only positive-query subtotal left={left} nb={nb} count={count} work={subtotal} worst={worst:?} vertex={vertex} stages={stages:?}; excludes selection/field-surfaces/closedfile"
        );
        assert!(count >= 44);
        assert_eq!(source_refusals > 0, nb);
    }
}

fn independent_oracle(curve: &CapCurvePiece<'_, '_>, proof: &Probe) {
    use super::restrict::tests::jet;
    use crate::step::hanger::encoded::cap::tests::{bezier, number};
    use num_rational::BigRational as Q;
    use num_traits::Signed;
    let homogeneous = |data: &super::literal::CartesianCap| {
        data.poles
            .iter()
            .zip(&data.weights)
            .map(|(p, w)| {
                let p = p.each_ref().map(number);
                let w = number(w);
                [&p[0] * &w, &p[1] * &w, &p[2] * &w, w]
            })
            .collect::<Vec<_>>()
    };
    let xyz = homogeneous(&proof.xyz);
    let uv = homogeneous(&proof.uv);
    let f = curve.field;
    let row = &f.meridians[f.cells[curve.source_cell].meridian];
    let controls = (0..=row.degree())
        .map(|i| {
            (0..=10)
                .map(|j| {
                    f.control(
                        curve.source_cell,
                        i,
                        curve.angular_piece,
                        j,
                        DomainOptions::default(),
                        None,
                    )
                    .unwrap()
                    .homogeneous
                    .each_ref()
                    .map(number)
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let ka = Q::from_float(f.angular.knots()[curve.angular_piece]).unwrap();
    let kw = Q::from_float(f.angular.knots()[curve.angular_piece + 1]).unwrap() - &ka;
    let span = &proof.parameter[1].0 - &proof.parameter[0].0;
    for k in 0..=4 {
        let t = Q::new(k.into(), 4.into());
        let u3 = &proof.parameter[0].0 + &span * &t;
        let (point, derivative) = jet(xyz.clone(), &u3);
        let derivative = derivative.map(|d| d * &span);
        let (foot, foot_derivative) = jet(uv.clone(), &t);
        let v = (&foot[1] - &ka) / &kw;
        let hrows = controls
            .iter()
            .map(|r| bezier(r.clone(), &v))
            .collect::<Vec<_>>();
        let h = bezier(hrows.clone(), &foot[0]);
        let hu = bezier(
            hrows
                .windows(2)
                .map(|p| {
                    std::array::from_fn::<_, 4, _>(|axis| {
                        (&p[1][axis] - &p[0][axis]) * Q::from_integer(row.degree().into())
                    })
                })
                .collect(),
            &foot[0],
        );
        let hvrows = controls
            .iter()
            .map(|r| {
                bezier(
                    r.windows(2)
                        .map(|p| {
                            std::array::from_fn::<_, 4, _>(|axis| {
                                (&p[1][axis] - &p[0][axis]) * Q::from_integer(10.into()) / &kw
                            })
                        })
                        .collect(),
                    &v,
                )
            })
            .collect::<Vec<_>>();
        let hv = bezier(hvrows, &foot[0]);
        let surface = std::array::from_fn::<_, 3, _>(|axis| &h[axis] / &h[3]);
        let surface_derivative = std::array::from_fn::<_, 3, _>(|axis| {
            ((&hu[axis] * &h[3] - &h[axis] * &hu[3]) * &foot_derivative[0]
                + (&hv[axis] * &h[3] - &h[axis] * &hv[3]) * &foot_derivative[1])
                / (&h[3] * &h[3])
        });
        let gap = (0..3).map(|i| (&point[i] - &surface[i]).abs()).sum::<Q>();
        let tangent = (0..3)
            .map(|i| (&derivative[i] - &surface_derivative[i]).abs())
            .sum::<Q>();
        assert!(gap <= Q::from_float(proof.position).unwrap());
        assert!(tangent <= Q::from_float(proof.derivative).unwrap());
    }
}

#[test]
fn same_actual_root_exact_endpoint_fallback_and_every_earned_checkpoint_remain_bounded() {
    use std::cell::Cell;
    for x in [
        [f64::NAN, 0., 0.],
        [-1., 0., 0.],
        [0., 2e-6, 0.],
        [0., 0., 2e-6],
    ] {
        assert!(matches!(source_widths(x), Err(Error::Unresolved)));
    }
    source_widths([0., 0., 0.]).unwrap();
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
        .position(|c| {
            c.source_key.total_period == -1
                && matches!(
                    c.source_key.face,
                    FaceKey::MetricThread {
                        band: ThreadBandRole::Root,
                        ..
                    }
                )
        })
        .unwrap();
    let curve = f
        .cap_curve_piece(index, End::Start, 3, DomainOptions::default(), None)
        .unwrap()
        .unwrap();
    let interval = run(&curve).unwrap();
    let exact = run_bound(&curve, None, 16_000_000, None, true).unwrap();
    source_gate(&exact).unwrap();
    independent_oracle(&curve, &exact);
    let mut memory = Memory::new(0, 64 * 1024 * 1024).unwrap();
    let mut work = Work::new(16_000_000, None).unwrap();
    let h = decode(curve.controls(), &mut work).unwrap();
    let targets = interval.parameter.each_ref().map(|u| {
        let p = eval(&h, u, &mut work).unwrap();
        std::array::from_fn::<_, 3, _>(|a| {
            EncodedScalar::retain(&p[a].div(&p[3], &mut work).unwrap(), &mut memory, &mut work)
                .unwrap()
        })
    });
    let exact_target = run_bound(&curve, Some(targets.each_ref()), 16_000_000, None, true).unwrap();
    source_gate(&exact_target).unwrap();
    // A real source-root envelope must also cover an explicitly coarser
    // parameter proposal. This is a numeric hull check, not curve admission.
    let half_root = f
        .cap_curve_piece(index, End::Start, 4, DomainOptions::default(), None)
        .unwrap()
        .unwrap();
    let proposed = [N::integer(0), N::integer(1)];
    let envelope = parameter_envelope(half_root.parameter_bounds, &proposed, &mut work).unwrap();
    assert_eq!(envelope[0].0, proposed[0].0);
    assert_eq!(envelope[1].0, proposed[1].0);
    assert!(
        exact.work > interval.work,
        "attempted interval proposal and exact fallback are both charged"
    );
    let calls = Cell::new(0);
    let count = || {
        calls.set(calls.get() + 1);
        false
    };
    run_bound(&curve, None, 16_000_000, Some(&count), false).unwrap();
    for stop in 1..=calls.get() {
        let n = Cell::new(0);
        let cancel = || {
            n.set(n.get() + 1);
            n.get() >= stop
        };
        assert_eq!(
            run_bound(&curve, None, 16_000_000, Some(&cancel), false).err(),
            Some(Error::Cancelled),
            "stop={stop}"
        );
    }
    eprintln!(
        "actual root carrier checkpoint callbacks={} intervalwork={} exactwork={}",
        calls.get(),
        interval.work,
        exact.work
    );
}
