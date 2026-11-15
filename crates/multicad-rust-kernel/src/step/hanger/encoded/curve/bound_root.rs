//! Actual finite root-cap numeric candidate, on caller-owned ledgers.
use super::correction::correct;
use super::literal::{Precision, encode_cap_at, quantize_at};
use super::restrict::{chosen, derivatives, restrict};
use super::*;
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
pub(super) fn decode<const D: usize>(
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
pub(super) struct Candidate {
    pub(super) work: usize,
    pub(super) position: f64,
    pub(super) derivative: f64,
    pub(super) source_position: f64,
    pub(super) source_derivative: f64,
    pub(super) original_source_derivative: f64,
    pub(super) vertex_error: f64,
    pub(super) parameter: [N; 2],
    pub(super) stages: [usize; 7],
    pub(super) xyz: super::literal::CartesianCap,
    pub(super) uv: super::literal::CartesianCap,
}
pub(super) fn parameter_envelope(
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
pub(super) fn source_gate(proof: &Candidate) -> Result<(), Error> {
    source_widths([
        proof.source_position,
        proof.source_derivative,
        proof.original_source_derivative,
    ])
}
pub(super) fn source_widths(widths: [f64; 3]) -> Result<(), Error> {
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

pub(super) fn bind(
    curve: &CapCurvePiece<'_, '_>,
    targets: Option<[&[EncodedScalar; 3]; 2]>,
    m: &mut Memory,
    w: &mut Work<'_>,
    arithmetic: EndpointArithmetic,
) -> Result<Candidate, Error> {
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
    let pitch = binary(pitch, w)?;
    let lead = pitch.mul(&N::integer(i64::from(starts)), w)?;
    let station = if curve.end == End::Start {
        N::integer(0)
    } else {
        binary(length, w)?
    };
    let shift = station.sub(
        &pitch.mul(&N::integer(i64::from(cell.source_key.total_period)), w)?,
        w,
    )?;
    let plane = curve.surface_trace.plane_z.number(w)?;
    let ka = binary(f.angular.knots()[curve.angular_piece], w)?;
    let kb = binary(f.angular.knots()[curve.angular_piece + 1], w)?;
    let kw = kb.sub(&ka, w)?;
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
        .map(|c| Ok([c[0].enclosure(w)?, c[1].enclosure(w)?]))
        .collect::<Result<Vec<_>, Error>>()?;
    let meridian = (0..=row.degree())
        .map(|i| {
            Ok([
                row.radius[i].number(w)?,
                row.phase[i].number(w)?,
                row.weights[i].number(w)?,
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
        let phase = curve.phase[e].number(w)?;
        let (selected, error) = chosen(row, &phase, curve.parameter_bounds[e], &pitch, w)?;
        parameter[e] = quantize64(&selected, w)?;
        parameter_error[e] =
            (p(error) + p(parameter[e].sub(&selected, w)?.enclosure(w)?.magnitude())).hi;
        vtarget[e] = quantize64(&shift.sub(&phase, w)?.div(&lead, w)?, w)?;
        if vtarget[e].cmp(&ka, w)? == Ordering::Less || vtarget[e].cmp(&kb, w)? == Ordering::Greater
        {
            return Err(Error::Unresolved);
        }
        let h = eval(&meridian, &parameter[e], w)?;
        let radius = h[0].div(&h[2], w)?;
        let q = h[1].div(&h[2], w)?;
        let vold = shift.sub(&q, w)?.div(&lead, w)?;
        let old_t = vold.sub(&ka, w)?.div(&kw, w)?;
        let target_t = vtarget[e].sub(&ka, w)?.div(&kw, w)?;
        let old_i = eval_interval(&angular_interval, old_t.enclosure(w)?, w)?;
        let target_i = eval_interval(&angular_interval, target_t.enclosure(w)?, w)?;
        let radius_i = radius.enclosure(w)?;
        let mut residual = p(0.);
        for axis in 0..2 {
            let point = radius_i * target_i[axis];
            points[e][axis] = if let Some(targets) = targets {
                targets[e][axis].number(w)?
            } else {
                quantize64(&binary(point.midpoint(), w)?, w)?
            };
            let required = points[e][axis].enclosure(w)? - radius_i * old_i[axis];
            delta[e][axis] = quantize64(&binary(required.midpoint(), w)?, w)?;
            residual = residual + p((required - delta[e][axis].enclosure(w)?).magnitude());
        }
        if arithmetic == EndpointArithmetic::Exact || !residual.valid() || residual.hi > 1e-7 {
            // Bounded precision fallback, with all failed proposal work retained.
            let old = eval(&angular, &old_t, w)?;
            let target = eval(&angular, &target_t, w)?;
            residual = p(0.);
            for axis in 0..2 {
                points[e][axis] = if let Some(targets) = targets {
                    targets[e][axis].number(w)?
                } else {
                    quantize64(&radius.mul(&target[axis], w)?, w)?
                };
                let required = points[e][axis].sub(&radius.mul(&old[axis], w)?, w)?;
                delta[e][axis] = quantize64(&required, w)?;
                residual =
                    residual + p(required.sub(&delta[e][axis], w)?.enclosure(w)?.magnitude());
            }
        }
        vertex_error = vertex_error.max(residual.hi);
        if let Some(targets) = targets
            && targets[e][2].number(w)?.0 != plane.0
        {
            return Err(Error::Unresolved);
        }
        points[e][2] = plane.clone();
        w.charge(512)?;
        uv_delta[e] = vtarget[e].sub(&vold, w)?;
    }
    stages[0] = w.used - baseline;
    let span = parameter[1].sub(&parameter[0], w)?;
    let span_bound = span.enclosure(w)?;
    if span_bound.lo <= 0. {
        return Err(Error::Unresolved);
    }
    let h = decode(curve.controls(), w)?;
    let corrected = correct(&h, &delta, &parameter[0], &parameter[1], m, w)?;
    stages[1] = w.used - baseline - stages[..1].iter().sum::<usize>();
    let xyz = encode64(&corrected.controls, None, CurveOptions::default(), m, w)?;
    stages[2] = w.used - baseline - stages[..2].iter().sum::<usize>();
    vertex_error = (p(vertex_error) + p(xyz.error[0])).hi;
    if vertex_error > 1e-7 || !vertex_error.is_finite() {
        return Err(Error::Unresolved);
    }
    let finite = restrict(&meridian, &parameter[0], &parameter[1], m, w)?;
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
            let d = lead.mul(&finite[j][2], w)?;
            let v = shift.mul(&finite[j][2], w)?.sub(&finite[j][1], w)?;
            den = den.add(&d.mul(&N::integer(scale as i64), w)?, w)?;
            num = num.add(
                &v.add(&uv_delta[e].mul(&d, w)?, w)?
                    .mul(&N::integer(scale as i64), w)?,
                w,
            )?;
            if e == 1 {
                linear = d.mul(&N::integer(i as i64), w)?;
            }
        }
        let u = parameter[0].mul(&den, w)?.add(&span.mul(&linear, w)?, w)?;
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
            EncodedScalar::retain(&parameter[e], m, w)?,
            EncodedScalar::retain(&vtarget[e], m, w)?,
            EncodedScalar::retain(&N::integer(0), m, w)?,
        ];
    }
    stages[3] = w.used - baseline - stages[..3].iter().sum::<usize>();
    let uv_literal = encode64(&uv, Some(&encoded_endpoints), CurveOptions::default(), m, w)?;
    stages[4] = w.used - baseline - stages[..4].iter().sum::<usize>();
    for point in &uv_literal.poles {
        let u = point[0].number(w)?;
        let v = point[1].number(w)?;
        if u.cmp(&N::integer(0), w)? == Ordering::Less
            || u.cmp(&N::integer(1), w)? == Ordering::Greater
        {
            return Err(Error::Unresolved);
        }
        if v.cmp(&ka, w)? == Ordering::Less || v.cmp(&kb, w)? == Ordering::Greater {
            return Err(Error::Unresolved);
        }
    }
    // Whole-trim lift proof is supplied below by the actual tensor Jacobian
    // and Hessian coefficient bounds, not by the endpoint correspondence.
    let original_uv = decode(curve.surface_trace.controls(), w)?;
    let original_uv = restrict(&original_uv, &parameter[0], &parameter[1], m, w)?;
    let orig4 = original_uv
        .iter()
        .map(|p| [p[0].clone(), p[1].clone(), N::integer(0), p[2].clone()])
        .collect::<Vec<_>>();
    let gp = derivatives(&orig4, &N::integer(1), w)?[0];
    let mut vrange = I::new(f64::INFINITY, f64::NEG_INFINITY);
    for point in &original_uv {
        let v = point[1].div(&point[2], w)?.enclosure(w)?;
        vrange = I::new(vrange.lo.min(v.lo), vrange.hi.max(v.hi));
    }
    for point in &uv_literal.poles {
        let v = point[1].number(w)?.enclosure(w)?;
        vrange = I::new(vrange.lo.min(v.lo), vrange.hi.max(v.hi));
    }
    stages[5] = w.used - baseline - stages[..5].iter().sum::<usize>();
    let [j, hessian] = tensor_bounds(curve, vrange, w)?;
    let duv = p(uv_delta[0]
        .enclosure(w)?
        .magnitude()
        .max(uv_delta[1].enclosure(w)?.magnitude()))
        + p(uv_literal.error[0]);
    let duvp =
        p(uv_delta[1].sub(&uv_delta[0], w)?.enclosure(w)?.magnitude()) + p(uv_literal.error[1]);
    let d3 = corrected.position_error + p(xyz.error[0]);
    let d3p = (corrected.derivative_error + p(xyz.error[1])) * span_bound;
    let lift = d3 + p(j) * duv;
    let lift_prime = d3p + p(j) * duvp + p(hessian) * duv * p(gp);
    if !lift.valid() || !lift_prime.valid() || lift.hi > 1e-7 || lift_prime.hi > 1e-6 {
        return Err(Error::Unresolved);
    }
    // True-event transfer: chosen finite parameters are proposals, not exact
    // algebraic endpoints. Compare the two finite curves under their actual
    // affine-in-root-u lambda correspondence before claiming source accuracy.
    let [low, high] = parameter_envelope(curve.parameter_bounds, &parameter, w)?;
    let jets = super::restrict::factor_derivatives(curve, &low, &high, m, w)?;
    let eps = p(parameter_error[0].max(parameter_error[1]));
    let event_position = p(jets[0]) * eps;
    let event_derivative =
        p(jets[1]) * eps + p(jets[0]) * (p(parameter_error[0]) + p(parameter_error[1]));
    let true_span = I::new(curve.parameter_bounds[1][0], curve.parameter_bounds[1][1])
        - I::new(curve.parameter_bounds[0][0], curve.parameter_bounds[0][1]);
    if true_span.lo <= 0. {
        return Err(Error::Unresolved);
    }
    w.charge(512)?;
    let dq = if row.degree() == 2 {
        let b = s.metric().root_meridian().whole_bounds();
        I::new(b.first_derivative_mm[0][0], b.first_derivative_mm[0][1])
    } else {
        meridian[row.degree()][1]
            .div(&meridian[row.degree()][2], w)?
            .sub(&meridian[0][1].div(&meridian[0][2], w)?, w)?
            .enclosure(w)?
    };
    let source_position = lift + event_position + p(row.bounds.position_error_mm);
    let source_derivative = lift_prime
        + event_derivative
        + (p(row.bounds.derivative_error_mm[0])
            + p(row.bounds.derivative_error_mm[1]) * p(dq.magnitude()) / lead.enclosure(w)?)
            * p(true_span.hi.min(1.));
    let original_span = curve.surface_trace.endpoints[1]
        .original_parameter
        .number(w)?
        .sub(
            &curve.surface_trace.endpoints[0]
                .original_parameter
                .number(w)?,
            w,
        )?;
    let phase_span = curve.surface_trace.endpoints[1]
        .phase
        .number(w)?
        .sub(&curve.surface_trace.endpoints[0].phase.number(w)?, w)?;
    let original_scale = original_span.div(&phase_span, w)?.enclosure(w)? * dq * true_span;
    if original_scale.lo <= 0. {
        return Err(Error::Unresolved);
    }
    let original_source_derivative = source_derivative / original_scale;

    w.charge(1)?;
    stages[6] = w.used - baseline - stages[..6].iter().sum::<usize>();
    // Independent exact-rational point/D1 oracle is added by the caller test.
    Ok(Candidate {
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

use super::jets::field_bounds;

/// Endpoint correction arithmetic: interval proposal with an exact fallback,
/// or direct rational evaluation of the same defining angular polynomial.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EndpointArithmetic {
    Adaptive,
    Exact,
}
pub struct BoundRootCap<'r, 'g, 'c, 'd, 'f, 's> {
    source: &'r CapEvents<'g, 'c, 'd, 'f, 's>,
    segment: usize,
    xyz: super::literal::CartesianCap,
    uv: super::literal::CartesianCap,
    parameter: [EncodedScalar; 2],
    field: [f64; 2],
    bounds: [f64; 3],
    vertex_error: f64,
    stages: [usize; 7],
    certificate_work: usize,
    work: usize,
    bytes: usize,
}
impl BoundRootCap<'_, '_, '_, '_, '_, '_> {
    pub fn source(&self) -> &CapEvents<'_, '_, '_, '_, '_> {
        self.source
    }
    pub fn segment_id(&self) -> usize {
        self.segment
    }
    pub fn points(&self) -> &[[EncodedScalar; 3]] {
        &self.xyz.poles
    }
    pub fn weights(&self) -> &[EncodedScalar] {
        &self.xyz.weights
    }
    pub fn surface_points(&self) -> &[[EncodedScalar; 3]] {
        &self.uv.poles
    }
    pub fn surface_weights(&self) -> &[EncodedScalar] {
        &self.uv.weights
    }
    pub fn parameter_range(&self) -> &[EncodedScalar; 2] {
        &self.parameter
    }
    pub fn field_bounds(&self) -> [f64; 2] {
        self.field
    }
    pub fn source_bounds(&self) -> [f64; 3] {
        self.bounds
    }
    pub fn vertex_incidence_error_mm(&self) -> f64 {
        self.vertex_error
    }
    pub fn stage_work(&self) -> [usize; 7] {
        self.stages
    }
    pub fn certificate_work_used(&self) -> usize {
        self.certificate_work
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_bytes(&self) -> usize {
        self.bytes
    }
}
impl<'g, 'c, 'd, 'f, 's> CapEvents<'g, 'c, 'd, 'f, 's> {
    pub fn bound_root_cap<'r>(
        &'r self,
        segment: usize,
        options: CurveOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Option<BoundRootCap<'r, 'g, 'c, 'd, 'f, 's>>, Error> {
        self.bound_root_cap_with_arithmetic(
            segment,
            options,
            EndpointArithmetic::Adaptive,
            interrupt,
        )
    }
    pub fn bound_root_cap_with_arithmetic<'r>(
        &'r self,
        segment: usize,
        options: CurveOptions,
        arithmetic: EndpointArithmetic,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Option<BoundRootCap<'r, 'g, 'c, 'd, 'f, 's>>, Error> {
        options.validate()?;
        let owner = self.segments().get(segment).ok_or(Error::Invalid)?;
        let field = self.source().source().domains().field();
        let (mut w, mut m) = crate::step::hanger::domain::ledgers(
            DomainOptions {
                maximum_work: options.work,
                maximum_bytes: options.bytes,
            },
            self.retained_bytes(),
            interrupt,
        )?;
        w.charge(self.work_used())?;
        m.add(SCRATCH)?;
        let base = m.used;
        m.add(size_of::<BoundRootCap<'_, '_, '_, '_, '_, '_>>())?;
        if field.meridians[field.cells[owner.source_cell()].meridian].degree() != 2 {
            return Ok(None);
        }
        let curve = field
            .cap_curve_piece_with_ledgers(
                owner.source_cell(),
                owner.end(),
                owner.angular_piece(),
                &mut m,
                &mut w,
            )?
            .ok_or(Error::Unresolved)?;
        let ids = owner.event_ids();
        let targets = [
            self.coordinates(ids[0]).ok_or(Error::Unresolved)?,
            self.coordinates(ids[1]).ok_or(Error::Unresolved)?,
        ];
        let candidate = bind(&curve, Some(targets), &mut m, &mut w, arithmetic)?;
        source_gate(&candidate)?;
        let bounds = [
            candidate.source_position,
            candidate.source_derivative,
            candidate.original_source_derivative,
        ];
        if bounds[0] > options.position_error_mm
            || bounds[1] > options.derivative_error_mm
            || bounds[2] > options.derivative_error_mm
        {
            return Err(Error::Unresolved);
        }
        let parameter = [
            EncodedScalar::retain(&candidate.parameter[0], &mut m, &mut w)?,
            EncodedScalar::retain(&candidate.parameter[1], &mut m, &mut w)?,
        ];
        let curve_bytes = curve.retained_bytes();
        drop(curve);
        m.used -= curve_bytes;
        w.charge(1)?;
        Ok(Some(BoundRootCap {
            source: self,
            segment,
            xyz: candidate.xyz,
            uv: candidate.uv,
            parameter,
            field: [candidate.position, candidate.derivative],
            bounds,
            vertex_error: candidate.vertex_error,
            stages: candidate.stages,
            certificate_work: candidate.work,
            work: w.used,
            bytes: m.used - base,
        }))
    }
}
#[cfg(test)]
mod tests;
