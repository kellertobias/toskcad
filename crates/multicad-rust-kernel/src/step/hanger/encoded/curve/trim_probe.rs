//! Bounded actual-source prototype. No file/native admission.
use super::correction::correct;
use super::literal::{Precision, encode_cap_at, quantize_at};
use super::restrict::{amplification, chosen, derivatives, restrict};
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
fn eval<const D: usize>(input: &[[N; D]], t: &N, w: &mut Work<'_>) -> Result<[N; D], Error> {
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
fn eval_interval<const D: usize>(
    input: &[[I; D]],
    t: I,
    w: &mut Work<'_>,
) -> Result<[I; D], Error> {
    w.charge(32 * D * input.len() * (input.len() - 1) / 2 + input.len() * D)?;
    let mut values = input.to_vec();
    let one = p(1.) - t;
    for n in (1..values.len()).rev() {
        for i in 0..n {
            let (lo, hi) = values.split_at_mut(i + 1);
            for (a, b) in lo[i].iter_mut().zip(&hi[0]) {
                *a = *a * one + *b * t;
            }
        }
    }
    Ok(values.remove(0))
}
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
fn norm(c: &[N], w: &mut Work<'_>) -> Result<I, Error> {
    let mut s = p(0.);
    for x in c {
        s = s + p(x.enclosure(w)?.magnitude());
    }
    Ok(s)
}
pub(super) struct Probe {
    pub(super) work: usize,
    pub(super) position: f64,
    pub(super) derivative: f64,
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
    let mut w = Work::new(limit, None)?;
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
    for e in 0..2 {
        let phase = curve.phase[e].number(&mut w)?;
        parameter[e] = quantize64(
            &chosen(row, &phase, curve.parameter_bounds[e], &pitch, &mut w)?.0,
            &mut w,
        )?;
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
            points[e][axis] = quantize64(&binary(point.midpoint(), &mut w)?, &mut w)?;
            let required = points[e][axis].enclosure(&mut w)? - radius_i * old_i[axis];
            delta[e][axis] = quantize64(&binary(required.midpoint(), &mut w)?, &mut w)?;
            residual = residual + p((required - delta[e][axis].enclosure(&mut w)?).magnitude());
        }
        if !residual.valid() || residual.hi > 1e-7 {
            // Bounded precision fallback, with all failed proposal work retained.
            let old = eval(&angular, &old_t, &mut w)?;
            let target = eval(&angular, &target_t, &mut w)?;
            residual = p(0.);
            for axis in 0..2 {
                points[e][axis] = quantize64(&radius.mul(&target[axis], &mut w)?, &mut w)?;
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
    stages[6] = w.used - baseline - stages[..6].iter().sum::<usize>();
    // Independent exact-rational point/D1 oracle is added by the caller test.
    Ok(Probe {
        work: w.used,
        position: lift.hi,
        derivative: lift_prime.hi,
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
pub(super) fn field_bounds(
    f: &EncodedField<'_>,
    index: usize,
    end: End,
    angular_piece: usize,
    vrange: I,
    w: &mut Work<'_>,
) -> Result<[f64; 2], Error> {
    let cell = f.cells[index];
    let r = &f.meridians[cell.meridian];
    let source = f.source.source_boundary().source();
    let (pitch, starts, left, length) = if r.component == Component::Metric {
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
    let shift = pitch
        .mul(&N::integer(i64::from(cell.source_key.total_period)), w)?
        .sub(&station, w)?;
    let a = binary(f.angular.knots()[angular_piece], w)?;
    let b = binary(f.angular.knots()[angular_piece + 1], w)?;
    let width = b.sub(&a, w)?;
    let t = (vrange - a.enclosure(w)?) / width.enclosure(w)?;
    let factors = [
        amplification(t, 10, w)?,
        amplification(t, 9, w)?,
        amplification(t, 8, w)?,
    ];
    let ap = if left {
        f.angular.coefficients().len() - 1 - angular_piece
    } else {
        angular_piece
    };
    let mut net = Vec::new();
    let mut weights = Vec::new();
    for i in 0..=r.degree() {
        let radius = r.radius[i].number(w)?;
        let q = r.phase[i].number(w)?;
        let weight = r.weights[i].number(w)?;
        let mut row = Vec::new();
        for j in 0..=10 {
            let c = &f.angular.coefficients()[ap][if left { 10 - j } else { j }];
            let v = a.add(
                &width
                    .mul(&N::integer(j as i64), w)?
                    .div(&N::integer(10), w)?,
                w,
            )?;
            w.charge(512)?;
            row.push([
                radius.mul(&N(c[0].clone()), w)?,
                radius.mul(&N(c[1].clone()), w)?,
                q.add(&weight.mul(&shift.add(&lead.mul(&v, w)?, w)?, w)?, w)?,
            ]);
        }
        net.push(row);
        weights.push(weight);
    }
    let du = r.degree();
    let mut raw = [0_f64; 6];
    for row in &net {
        for v in row {
            raw[0] = raw[0].max(norm(v, w)?.hi);
        }
    }
    for rows in net.windows(2) {
        for (upper, lower) in rows[1].iter().zip(&rows[0]) {
            let d = std::array::from_fn::<_, 3, _>(|_| N::integer(0));
            let mut d = d;
            for (axis, out) in d.iter_mut().enumerate() {
                *out = upper[axis]
                    .sub(&lower[axis], w)?
                    .mul(&N::integer(du as i64), w)?;
            }
            raw[1] = raw[1].max(norm(&d, w)?.hi);
        }
    }
    for row in &net {
        for p in row.windows(2) {
            let mut d: [N; 3] = std::array::from_fn(|_| N::integer(0));
            for axis in 0..3 {
                d[axis] = p[1][axis]
                    .sub(&p[0][axis], w)?
                    .mul(&N::integer(10), w)?
                    .div(&width, w)?;
            }
            raw[2] = raw[2].max(norm(&d, w)?.hi);
        }
    }
    for rows in net.windows(3) {
        for ((upper, middle), lower) in rows[2].iter().zip(&rows[1]).zip(&rows[0]) {
            let mut d: [N; 3] = std::array::from_fn(|_| N::integer(0));
            for axis in 0..3 {
                d[axis] = upper[axis]
                    .sub(&middle[axis].mul(&N::integer(2), w)?, w)?
                    .add(&lower[axis], w)?
                    .mul(&N::integer((du * (du - 1)) as i64), w)?;
            }
            raw[3] = raw[3].max(norm(&d, w)?.hi);
        }
    }
    for rows in net.windows(2) {
        for j in 0..10 {
            let mut d: [N; 3] = std::array::from_fn(|_| N::integer(0));
            for axis in 0..3 {
                d[axis] = rows[1][j + 1][axis]
                    .sub(&rows[1][j][axis], w)?
                    .sub(&rows[0][j + 1][axis].sub(&rows[0][j][axis], w)?, w)?
                    .mul(&N::integer((du * 10) as i64), w)?
                    .div(&width, w)?;
            }
            raw[4] = raw[4].max(norm(&d, w)?.hi);
        }
    }
    let square = width.mul(&width, w)?;
    for row in &net {
        for p in row.windows(3) {
            let mut d: [N; 3] = std::array::from_fn(|_| N::integer(0));
            for axis in 0..3 {
                d[axis] = p[2][axis]
                    .sub(&p[1][axis].mul(&N::integer(2), w)?, w)?
                    .add(&p[0][axis], w)?
                    .mul(&N::integer(90), w)?
                    .div(&square, w)?;
            }
            raw[5] = raw[5].max(norm(&d, w)?.hi);
        }
    }
    let mut minimum = f64::INFINITY;
    for x in &weights {
        minimum = minimum.min(x.enclosure(w)?.lo);
    }
    let mut wu = 0_f64;
    for p in weights.windows(2) {
        wu = wu.max(
            p[1].sub(&p[0], w)?
                .mul(&N::integer(du as i64), w)?
                .enclosure(w)?
                .magnitude(),
        );
    }
    let mut wuu = 0_f64;
    for p in weights.windows(3) {
        wuu = wuu.max(
            p[2].sub(&p[1].mul(&N::integer(2), w)?, w)?
                .add(&p[0], w)?
                .mul(&N::integer((du * (du - 1)) as i64), w)?
                .enclosure(w)?
                .magnitude(),
        );
    }
    w.charge(1024)?;
    let a = p(raw[0]) * factors[0];
    let au = p(raw[1]) * factors[0];
    let av = p(raw[2]) * factors[1];
    let auu = p(raw[3]) * factors[0];
    let auv = p(raw[4]) * factors[1];
    let avv = p(raw[5]) * factors[2];
    let lo = p(minimum);
    let fu = au / lo + a * p(wu) / lo.square();
    let fv = av / lo;
    let fuu = auu / lo
        + p(2.) * au * p(wu) / lo.square()
        + a * p(wuu) / lo.square()
        + p(2.) * a * p(wu).square() / (lo.square() * lo);
    let fuv = auv / lo + av * p(wu) / lo.square();
    let fvv = avv / lo;
    Ok([fu.hi.max(fv.hi), fuu.hi.max(fuv.hi).max(fvv.hi)])
}
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
        eprintln!(
            "actual full-basis+UV left={left} nb={nb} cumulative_work={} position={} derivative={}",
            proof.work, proof.position, proof.derivative
        );
    }
}

#[test]
fn all_actual_cap_pieces_prove_full_basis_uv_lift_under_original_query_limits() {
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
        for index in 0..field.cells().len() {
            for end in [End::Start, End::End] {
                for piece in 0..field.angular.coefficients().len() {
                    let Some(curve) = field
                        .cap_curve_piece(index, end, piece, DomainOptions::default(), None)
                        .unwrap()
                    else {
                        continue;
                    };
                    let proof = run(&curve).unwrap();
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
            "caps-only positive-query subtotal left={left} nb={nb} count={count} work={subtotal} worst={worst:?} vertex={vertex} stages={stages:?}; excludes selection/field-surfaces/closedfile"
        );
        assert!(count >= 44);
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
