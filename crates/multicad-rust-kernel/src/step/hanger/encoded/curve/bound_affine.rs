//! Bound finite affine cap curves and PCurves, retaining actual numeric output.
use super::*;

use num_traits::{One, Zero};
type Q = BigRational;
pub(super) struct Candidate {
    pub(super) xyz: super::literal::CartesianCap,
    pub(super) uv: super::literal::CartesianCap,
    pub(super) parameter: [N; 2],
    pub(super) bounds: [f64; 3],
    pub(super) vertex_error: f64,
}
pub(super) fn bind(
    f: &EncodedField<'_>,
    owner: (usize, End, usize),
    targets: [&[EncodedScalar; 3]; 2],
    cap: super::affine::AffineCap,
    m: &mut Memory,
    w: &mut Work<'_>,
) -> Result<Candidate, Error> {
    let (index, end, piece) = owner;
    use super::literal::{Precision, encode_cap_at, quantize_at};

    // Caller endpoints retain ORIGINAL increasing-phase coedge ownership.
    // This actual affine-v basis has decreasing phase: reverse only its
    // target correspondence, never the original graph edge/coedge IDs.
    // The actual affine constructor always increases v and decreases phase.
    let targets = [targets[1], targets[0]];
    let trace_before = m.used;
    let trace = f
        .cap_trace_with_ledgers(index, end, m, w)?
        .ok_or(Error::Invalid)?;
    let trace_bytes = m.used - trace_before;
    let mut chosen = [
        quantize_at(&cap.trim[0], Precision::Bits64, w)?,
        quantize_at(&cap.trim[1], Precision::Bits64, w)?,
    ];
    // Inward finite proposals keep the complete original UV line in the
    // actual source patch. The omitted true-end strips remain in the error.
    let quantum = N(BigRational::new(
        1.into(),
        num_bigint::BigInt::from(1) << 64,
    ));
    if chosen[0].cmp(&cap.trim[0], w)? == Ordering::Less {
        chosen[0] = chosen[0].add(&quantum, w)?;
    }
    if chosen[1].cmp(&cap.trim[1], w)? == Ordering::Greater {
        chosen[1] = chosen[1].sub(&quantum, w)?;
    }
    let span = chosen[1].sub(&chosen[0], w)?;
    let span_i = span.enclosure(w)?;
    if span_i.lo <= 0. {
        return Err(Error::Unresolved);
    }
    let mut delta =
        std::array::from_fn::<_, 2, _>(|_| std::array::from_fn::<_, 3, _>(|_| N::integer(0)));
    let intervals = cap
        .controls
        .iter()
        .map(|c| Ok([c[0].enclosure(w)?, c[1].enclosure(w)?, c[3].enclosure(w)?]))
        .collect::<Result<Vec<_>, Error>>()?;
    let mut vertex_error = 0_f64;
    let plane = cap.controls[0][2].div(&cap.controls[0][3], w)?;
    for e in 0..2 {
        let actual = super::jets::eval_interval(&intervals, chosen[e].enclosure(w)?, w)?;
        if targets[e][2].number(w)?.0 != plane.0 {
            return Err(Error::Unresolved);
        }
        let mut residual = p(0.);
        for axis in 0..2 {
            let required = targets[e][axis].number(w)?.enclosure(w)? - actual[axis] / actual[2];
            delta[e][axis] = quantize_at(&binary(required.midpoint(), w)?, Precision::Bits64, w)?;
            residual = residual + p((required - delta[e][axis].enclosure(w)?).magnitude());
        }
        vertex_error = vertex_error.max(residual.hi);
    }
    let corrected =
        super::correction::correct(&cap.controls, &delta, &chosen[0], &chosen[1], m, w)?;
    let xyz = encode_cap_at(
        &corrected.controls,
        None,
        CurveOptions::default(),
        m,
        w,
        Precision::Bits64,
    )?;
    // Release defining/corrected rational nets before whole-patch jets.
    // The live exact-object bound of affine::reserve remains512 BigInts:
    // no raw net overlaps the derivative-net construction below.
    let original_jets = super::restrict::derivatives(&cap.controls, &N::integer(1), w)?;
    drop(cap.controls);
    let correction_error = [corrected.position_error, corrected.derivative_error];
    m.release(&corrected.controls);
    drop(corrected);
    let mut original =
        std::array::from_fn::<_, 2, _>(|_| std::array::from_fn::<_, 2, _>(|_| N::integer(0)));
    let mut desired = original.clone();
    for e in 0..2 {
        for axis in 0..2 {
            original[e][axis] =
                cap.uv[0][axis].add(&cap.derivative[axis].mul(&chosen[e], w)?, w)?;
            let truth = cap.uv[0][axis].add(&cap.derivative[axis].mul(&cap.trim[e], w)?, w)?;
            desired[e][axis] = quantize_at(&truth, Precision::Bits64, w)?;
        }
    }
    let uvh = desired
        .each_ref()
        .map(|c| [c[0].clone(), c[1].clone(), N::integer(0), N::integer(1)]);
    let uv = encode_cap_at(&uvh, None, CurveOptions::default(), m, w, Precision::Bits64)?;
    let mut vrange = I::new(f64::INFINITY, f64::NEG_INFINITY);
    for point in &uv.poles {
        let u = point[0].number(w)?;
        let v = point[1].number(w)?;
        if u.0 < Q::zero() || u.0 > Q::one() {
            return Err(Error::Unresolved);
        }
        if v.cmp(&cap.uv[0][1], w)? == Ordering::Less
            || v.cmp(&cap.uv[1][1], w)? == Ordering::Greater
        {
            return Err(Error::Unresolved);
        }
        let i = v.enclosure(w)?;
        vrange = I::new(vrange.lo.min(i.lo), vrange.hi.max(i.hi));
    }
    let mut offset =
        std::array::from_fn::<_, 2, _>(|_| std::array::from_fn::<_, 2, _>(|_| N::integer(0)));
    let mut duv = 0_f64;
    for e in 0..2 {
        let mut norm = p(0.);
        for axis in 0..2 {
            offset[e][axis] = desired[e][axis].sub(&original[e][axis], w)?;
            norm = norm + p(offset[e][axis].enclosure(w)?.magnitude());
        }
        duv = duv.max(norm.hi);
        let v = original[e][1].enclosure(w)?;
        vrange = I::new(vrange.lo.min(v.lo), vrange.hi.max(v.hi));
    }
    let mut duvp = p(0.);
    for (upper, lower) in offset[1].iter().zip(&offset[0]) {
        duvp = duvp + p(upper.sub(lower, w)?.enclosure(w)?.magnitude());
    }
    let [j, h] = super::jets::field_bounds(f, index, end, piece, vrange, w)?;
    let gp = (p(cap.derivative[0].enclosure(w)?.magnitude())
        + p(cap.derivative[1].enclosure(w)?.magnitude()))
        * span_i;
    let d0 = p(duv) + p(uv.error[0]);
    let d1 = duvp + p(uv.error[1]);
    let field_position = correction_error[0] + p(xyz.error[0]) + p(j) * d0;
    let field_derivative =
        (correction_error[1] + p(xyz.error[1])) * span_i + p(j) * d1 + p(h) * d0 * gp;
    if field_position.hi > 1e-7 || field_derivative.hi > 1e-6 || vertex_error + xyz.error[0] > 1e-7
    {
        return Err(Error::Unresolved);
    }
    let eps = chosen
        .each_ref()
        .into_iter()
        .zip(cap.trim.each_ref())
        .map(|(a, b)| Ok(a.sub(b, w)?.enclosure(w)?.magnitude()))
        .collect::<Result<Vec<_>, Error>>()?;
    let jets = original_jets;
    let event_position = p(jets[0]) * p(eps[0].max(eps[1]));
    let event_derivative =
        p(jets[1]) * p(eps[0].max(eps[1])) + p(jets[0]) * (p(eps[0]) + p(eps[1]));
    let row = &f.meridians[f.cells[index].meridian];
    let true_span = cap.trim[1].sub(&cap.trim[0], w)?.enclosure(w)?;
    let source_position = field_position + event_position + p(row.bounds.position_error_mm);
    let source_derivative = field_derivative
        + event_derivative
        + (p(row.bounds.derivative_error_mm[0]) * p(cap.derivative[0].enclosure(w)?.magnitude())
            + p(row.bounds.derivative_error_mm[1])
                * p(cap.derivative[1].enclosure(w)?.magnitude()))
            * true_span;
    let ratio = trace.endpoints[1]
        .original_parameter
        .number(w)?
        .sub(&trace.endpoints[0].original_parameter.number(w)?, w)?
        .div(
            &trace.endpoints[1]
                .phase
                .number(w)?
                .sub(&trace.endpoints[0].phase.number(w)?, w)?,
            w,
        )?;
    let q0 = row.phase[0].number(w)?.div(&row.weights[0].number(w)?, w)?;
    let q1 = row.phase[1].number(w)?.div(&row.weights[1].number(w)?, w)?;
    let gamma = q1.sub(&q0, w)?.mul(&cap.derivative[0], w)?.mul(&ratio, w)?;
    if gamma.0 >= Q::zero() {
        return Err(Error::Unresolved);
    }
    let scale = -(gamma.enclosure(w)? * true_span);
    if scale.lo <= 0. {
        return Err(Error::Unresolved);
    }
    let original_derivative = source_derivative / scale;
    if !source_position.hi.is_finite()
        || !source_derivative.hi.is_finite()
        || !original_derivative.hi.is_finite()
        || source_position.hi > 1e-7
        || source_derivative.hi > 1e-6
        || original_derivative.hi > 1e-6
    {
        return Err(Error::Unresolved);
    }
    drop(trace);
    m.used -= trace_bytes;
    let xyz_error = xyz.error[0];
    Ok(Candidate {
        xyz,
        uv,
        parameter: chosen,
        bounds: [
            source_position.hi,
            source_derivative.hi,
            original_derivative.hi,
        ],
        vertex_error: vertex_error + xyz_error,
    })
}

pub struct BoundAffineCap<'r, 'g, 'c, 'd, 'f, 's> {
    source: &'r CapEvents<'g, 'c, 'd, 'f, 's>,
    segment: usize,
    xyz: super::literal::CartesianCap,
    uv: super::literal::CartesianCap,
    bounds: [f64; 3],
    vertex_error: f64,
    parameter: [EncodedScalar; 2],
    work: usize,
    bytes: usize,
}
impl BoundAffineCap<'_, '_, '_, '_, '_, '_> {
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
    /// Source position, finite-interval D1, ORIGINAL source-coedge D1.
    pub fn source_bounds(&self) -> [f64; 3] {
        self.bounds
    }
    pub fn vertex_incidence_error_mm(&self) -> f64 {
        self.vertex_error
    }
    pub fn original_phase_reversed(&self) -> bool {
        true
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_bytes(&self) -> usize {
        self.bytes
    }
}
impl<'g, 'c, 'd, 'f, 's> CapEvents<'g, 'c, 'd, 'f, 's> {
    pub fn bound_affine_cap<'r>(
        &'r self,
        segment: usize,
        options: CurveOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Option<BoundAffineCap<'r, 'g, 'c, 'd, 'f, 's>>, Error> {
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
        super::affine::reserve(field, options.bytes, &mut m, &mut w)?;
        let base = m.used;
        m.add(size_of::<BoundAffineCap<'_, '_, '_, '_, '_, '_>>())?;
        let Some(cap) = super::affine::build(
            field,
            owner.source_cell(),
            owner.end(),
            owner.angular_piece(),
            &mut w,
        )?
        else {
            return Ok(None);
        };
        let ids = owner.event_ids();
        let targets = [
            self.coordinates(ids[0]).ok_or(Error::Unresolved)?,
            self.coordinates(ids[1]).ok_or(Error::Unresolved)?,
        ];
        let candidate = bind(
            field,
            (owner.source_cell(), owner.end(), owner.angular_piece()),
            targets,
            cap,
            &mut m,
            &mut w,
        )?;
        if candidate.bounds[0] > options.position_error_mm
            || candidate.bounds[1] > options.derivative_error_mm
            || candidate.bounds[2] > options.derivative_error_mm
        {
            return Err(Error::Unresolved);
        }
        let parameter = [
            EncodedScalar::retain(&candidate.parameter[0], &mut m, &mut w)?,
            EncodedScalar::retain(&candidate.parameter[1], &mut m, &mut w)?,
        ];
        w.charge(1)?;
        Ok(Some(BoundAffineCap {
            source: self,
            segment,
            xyz: candidate.xyz,
            uv: candidate.uv,
            bounds: candidate.bounds,
            vertex_error: candidate.vertex_error,
            parameter,
            work: w.used,
            bytes: m.used - base,
        }))
    }
}
#[cfg(test)]
mod tests;
