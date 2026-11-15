//! Finite restriction of the actual rational curve, not a new sampled law.
use super::*;
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RestrictionBounds {
    /// Difference from the true encoded finite interval under linear parameter
    /// correspondence, before Cartesian encoding or endpoint corrections.
    pub position_error_mm: f64,
    pub derivative_error_mm: f64,
    pub meridian_parameter_error: [f64; 2],
}
pub struct RestrictedCapCurve<'c, 'f, 's> {
    source: &'c CapCurvePiece<'f, 's>,
    parameter: [EncodedScalar; 2],
    controls: Vec<[EncodedScalar; 4]>,
    bounds: RestrictionBounds,
    work: usize,
    bytes: usize,
}
impl RestrictedCapCurve<'_, '_, '_> {
    pub fn source(&self) -> &CapCurvePiece<'_, '_> {
        self.source
    }
    pub fn meridian_parameters(&self) -> &[EncodedScalar; 2] {
        &self.parameter
    }
    pub fn controls(&self) -> &[[EncodedScalar; 4]] {
        &self.controls
    }
    pub fn bounds(&self) -> RestrictionBounds {
        self.bounds
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_bytes(&self) -> usize {
        self.bytes
    }
}
fn copy<const D: usize>(
    v: &[[N; D]],
    m: &mut Memory,
    w: &mut Work<'_>,
) -> Result<Vec<[N; D]>, Error> {
    let mut out = m.allocate(v.len(), w)?;
    for row in v {
        for x in row {
            w.charge(1 + ((x.0.numer().bits() + x.0.denom().bits()) as usize).div_ceil(64))?;
        }
        out.push(row.clone());
    }
    Ok(out)
}
fn split<const D: usize>(
    v: &[[N; D]],
    t: &N,
    m: &mut Memory,
    w: &mut Work<'_>,
) -> Result<[Vec<[N; D]>; 2], Error> {
    let mut work = copy(v, m, w)?;
    let mut left = m.allocate(v.len(), w)?;
    let mut right = m.allocate(v.len(), w)?;
    for _ in v {
        left.push(std::array::from_fn(|_| N::integer(0)));
        right.push(std::array::from_fn(|_| N::integer(0)));
    }
    let n = v.len() - 1;
    let one = N::integer(1).sub(t, w)?;
    for level in 0..=n {
        for x in work[0].iter().chain(&work[n - level]) {
            w.charge(1 + ((x.0.numer().bits() + x.0.denom().bits()) as usize).div_ceil(64))?;
        }
        left[level] = work[0].clone();
        right[n - level] = work[n - level].clone();
        for i in 0..n - level {
            let (prefix, suffix) = work.split_at_mut(i + 1);
            for (slot, next) in prefix[i].iter_mut().zip(&suffix[0]) {
                *slot = slot.mul(&one, w)?.add(&next.mul(t, w)?, w)?;
            }
        }
    }
    m.release(&work);
    drop(work);
    Ok([left, right])
}
pub(super) fn restrict<const D: usize>(
    h: &[[N; D]],
    a: &N,
    b: &N,
    m: &mut Memory,
    w: &mut Work<'_>,
) -> Result<Vec<[N; D]>, Error> {
    // All callers have a genuine positive interval inside[0,1]; division by b
    // is therefore exact and nonzero. Keep generic numeric validation here so
    // this scratch primitive cannot accidentally extrapolate a later caller.
    if h.len() < 2
        || h.len() > 23
        || a.cmp(&N::integer(0), w)? == Ordering::Less
        || b.cmp(&N::integer(1), w)? == Ordering::Greater
        || a.cmp(b, w)? != Ordering::Less
    {
        return Err(Error::Invalid);
    }
    if a.cmp(&N::integer(0), w)? == Ordering::Equal && b.cmp(&N::integer(1), w)? == Ordering::Equal
    {
        return copy(h, m, w);
    }
    let low = if b.cmp(&N::integer(1), w)? == Ordering::Equal {
        copy(h, m, w)?
    } else {
        let [low, discard] = split(h, b, m, w)?;
        m.release(&discard);
        drop(discard);
        low
    };
    // Exact De Casteljau endpoints: splitting again at zero only copies the
    // complete left polygon. Retain it directly with the same coefficient law.
    if a.cmp(&N::integer(0), w)? == Ordering::Equal {
        return Ok(low);
    }
    let [discard, out] = split(&low, &a.div(b, w)?, m, w)?;
    m.release(&low);
    m.release(&discard);
    drop(low);
    drop(discard);
    Ok(out)
}
pub(super) fn chosen(
    row: &EncodedMeridian,
    phase: &N,
    enclosure: [f64; 2],
    pitch: &N,
    w: &mut Work<'_>,
) -> Result<(N, f64), Error> {
    let a = row.phase[0].number(w)?;
    let b = row.phase[row.degree()].number(w)?;
    if phase.cmp(&a, w)? == Ordering::Equal {
        return Ok((N::integer(0), 0.));
    }
    if phase.cmp(&b, w)? == Ordering::Equal {
        return Ok((N::integer(1), 0.));
    }
    if row.degree() == 1 {
        return Ok((phase.sub(&a, w)?.div(&b.sub(&a, w)?, w)?, 0.));
    }
    // Actual root Q=[3p/8,p/4,5p/8], W=[1,1/2,1]: q(1/2)=p/2
    // identically. This is a symbolic equality witness, not event snapping.
    let half = N::integer(1).div(&N::integer(2), w)?;
    if phase.cmp(&pitch.mul(&half, w)?, w)? == Ordering::Equal {
        return Ok((half, 0.));
    }
    let proposal = (enclosure[0] + enclosure[1]) / 2.;
    let value = binary(proposal, w)?;
    let error = (value.enclosure(w)? - I::new(enclosure[0], enclosure[1])).magnitude();
    Ok((value, error))
}
/// XY derivative bounds of the actual plane-valued rational curve over an
/// exact restricted domain. Z=constant is structural, so distant axial datums
/// do not create false derivative amplification through a quotient estimate.
#[cfg(test)]
pub(super) fn derivatives(h: &[[N; 4]], span: &N, w: &mut Work<'_>) -> Result<[f64; 2], Error> {
    let n = h.len() - 1;
    let mut a = 0_f64;
    let mut da = 0_f64;
    let mut dda = 0_f64;
    let mut minimum = f64::INFINITY;
    let mut dw = 0_f64;
    let mut ddw = 0_f64;
    for row in h {
        a = a.max((p(row[0].enclosure(w)?.magnitude()) + p(row[1].enclosure(w)?.magnitude())).hi);
        minimum = minimum.min(row[3].enclosure(w)?.lo);
    }
    for pair in h.windows(2) {
        let mut sum = p(0.);
        for axis in [0, 1, 3] {
            let d = pair[1][axis]
                .sub(&pair[0][axis], w)?
                .mul(&N::integer(n as i64), w)?
                .div(span, w)?
                .enclosure(w)?
                .magnitude();
            if axis == 3 {
                dw = dw.max(d);
            } else {
                sum = sum + p(d);
            }
        }
        da = da.max(sum.hi);
    }
    for three in h.windows(3) {
        let mut sum = p(0.);
        for axis in [0, 1, 3] {
            let d = three[2][axis]
                .sub(&three[1][axis].mul(&N::integer(2), w)?, w)?
                .add(&three[0][axis], w)?
                .mul(&N::integer((n * (n - 1)) as i64), w)?
                .div(&span.mul(span, w)?, w)?
                .enclosure(w)?
                .magnitude();
            if axis == 3 {
                ddw = ddw.max(d);
            } else {
                sum = sum + p(d);
            }
        }
        dda = dda.max(sum.hi);
    }
    let lo = p(minimum);
    let first = p(da) / lo + p(a) * p(dw) / lo.square();
    let second = p(dda) / lo
        + p(2.) * p(da) * p(dw) / lo.square()
        + p(a) * p(ddw) / lo.square()
        + p(2.) * p(a) * p(dw).square() / (lo.square() * lo);
    if minimum <= 0. || !first.valid() || !second.valid() {
        return Err(Error::Unresolved);
    }
    Ok([first.hi, second.hi])
}
// Scalar numerator translation is EXACT. It avoids inventing radial/phase
// derivative growth from an arbitrary constant datum in a quotient bound.
struct ScalarJet {
    value: I,
    first: f64,
    second: f64,
}
fn scalar_jet(h: &[[N; 3]], axis: usize, span: &N, w: &mut Work<'_>) -> Result<ScalarJet, Error> {
    let n = h.len() - 1;
    let reference = h[0][axis].div(&h[0][2], w)?;
    let mut a = 0_f64;
    let mut da = 0_f64;
    let mut dda = 0_f64;
    let mut minimum = f64::INFINITY;
    let mut dw = 0_f64;
    let mut ddw = 0_f64;
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    for row in h {
        let value = row[axis].div(&row[2], w)?.enclosure(w)?;
        lo = lo.min(value.lo);
        hi = hi.max(value.hi);
        a = a.max(
            row[axis]
                .sub(&reference.mul(&row[2], w)?, w)?
                .enclosure(w)?
                .magnitude(),
        );
        minimum = minimum.min(row[2].enclosure(w)?.lo);
    }
    for pair in h.windows(2) {
        let wd = pair[1][2].sub(&pair[0][2], w)?;
        let ad = pair[1][axis]
            .sub(&pair[0][axis], w)?
            .sub(&reference.mul(&wd, w)?, w)?;
        da = da.max(
            ad.mul(&N::integer(n as i64), w)?
                .div(span, w)?
                .enclosure(w)?
                .magnitude(),
        );
        dw = dw.max(
            wd.mul(&N::integer(n as i64), w)?
                .div(span, w)?
                .enclosure(w)?
                .magnitude(),
        );
    }
    let square = span.mul(span, w)?;
    for t in h.windows(3) {
        let wd = t[2][2]
            .sub(&t[1][2].mul(&N::integer(2), w)?, w)?
            .add(&t[0][2], w)?;
        let ad = t[2][axis]
            .sub(&t[1][axis].mul(&N::integer(2), w)?, w)?
            .add(&t[0][axis], w)?
            .sub(&reference.mul(&wd, w)?, w)?;
        dda = dda.max(
            ad.mul(&N::integer((n * (n - 1)) as i64), w)?
                .div(&square, w)?
                .enclosure(w)?
                .magnitude(),
        );
        ddw = ddw.max(
            wd.mul(&N::integer((n * (n - 1)) as i64), w)?
                .div(&square, w)?
                .enclosure(w)?
                .magnitude(),
        );
    }
    w.charge(512)?;
    let den = p(minimum);
    let first = p(da) / den + p(a) * p(dw) / den.square();
    let second = p(dda) / den
        + p(2.) * p(da) * p(dw) / den.square()
        + p(a) * p(ddw) / den.square()
        + p(2.) * p(a) * p(dw).square() / (den.square() * den);
    let value = I::new(lo, hi);
    if minimum <= 0. || !value.valid() || !first.valid() || !second.valid() {
        return Err(Error::Unresolved);
    }
    Ok(ScalarJet {
        value,
        first: first.hi,
        second: second.hi,
    })
}
pub(super) fn amplification(t: I, degree: usize, w: &mut Work<'_>) -> Result<I, Error> {
    // Bernstein basis L1 sum is (|t|+|1-t|)^degree, also outside[0,1].
    // That convex piecewise-linear base attains its maximum at a range end.
    w.charge(32 * (degree + 1))?;
    let at = |x| (p(p(x).magnitude()) + p((p(1.) - p(x)).magnitude())).hi;
    let base = at(t.lo).max(at(t.hi));
    let mut out = p(1.);
    for _ in 0..degree {
        out = out * p(base);
    }
    Ok(out)
}
pub(super) fn factor_derivatives(
    curve: &CapCurvePiece<'_, '_>,
    low: &N,
    high: &N,
    m: &mut Memory,
    w: &mut Work<'_>,
) -> Result<[f64; 2], Error> {
    let f = curve.field;
    let cell = f.cells[curve.source_cell];
    let row = &f.meridians[cell.meridian];
    let source = f.source.source_boundary().source();
    let (pitch, starts, left, length) = if row.component == Component::Metric {
        let s = source.metric().spec();
        (s.pitch_mm, s.starts, s.left_hand, s.length_mm)
    } else {
        let s = source.wood().spec();
        (s.pitch_mm, s.starts, s.left_hand, s.length_mm)
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
    let mut meridian = m.allocate(row.degree() + 1, w)?;
    for i in 0..=row.degree() {
        meridian.push([
            row.radius[i].number(w)?,
            row.phase[i].number(w)?,
            row.weights[i].number(w)?,
        ]);
    }
    let restricted = restrict(&meridian, low, high, m, w)?;
    let span = high.sub(low, w)?;
    let r = scalar_jet(&restricted, 0, &span, w)?;
    let q = scalar_jet(&restricted, 1, &span, w)?;
    m.release(&meridian);
    m.release(&restricted);
    drop(meridian);
    drop(restricted);
    let a = binary(f.angular.knots()[curve.angular_piece], w)?;
    let b = binary(f.angular.knots()[curve.angular_piece + 1], w)?;
    let width = b.sub(&a, w)?;
    let lead = lead.enclosure(w)?;
    let angular_range =
        ((shift.enclosure(w)? - q.value) / lead - a.enclosure(w)?) / width.enclosure(w)?;
    let piece = if left {
        f.angular.coefficients().len() - 1 - curve.angular_piece
    } else {
        curve.angular_piece
    };
    let coefficients = &f.angular.coefficients()[piece];
    let mut c = m.allocate(11, w)?;
    for i in 0..=10 {
        let j = if left { 10 - i } else { i };
        w.charge(
            1 + coefficients[j]
                .iter()
                .map(|x| ((x.numer().bits() + x.denom().bits()) as usize).div_ceil(64))
                .sum::<usize>(),
        )?;
        c.push(coefficients[j].clone().map(N));
    }
    let mut norms = [0_f64; 3];
    for pair in &c {
        norms[0] = norms[0]
            .max((p(pair[0].enclosure(w)?.magnitude()) + p(pair[1].enclosure(w)?.magnitude())).hi);
    }
    for pair in c.windows(2) {
        let mut sum = p(0.);
        for (upper, lower) in pair[1].iter().zip(&pair[0]) {
            sum = sum
                + p(upper
                    .sub(lower, w)?
                    .mul(&N::integer(10), w)?
                    .div(&width, w)?
                    .enclosure(w)?
                    .magnitude());
        }
        norms[1] = norms[1].max(sum.hi);
    }
    let width_square = width.mul(&width, w)?;
    for t in c.windows(3) {
        let mut sum = p(0.);
        for ((upper, middle), lower) in t[2].iter().zip(&t[1]).zip(&t[0]) {
            sum = sum
                + p(upper
                    .sub(&middle.mul(&N::integer(2), w)?, w)?
                    .add(lower, w)?
                    .mul(&N::integer(90), w)?
                    .div(&width_square, w)?
                    .enclosure(w)?
                    .magnitude());
        }
        norms[2] = norms[2].max(sum.hi);
    }
    m.release(&c);
    drop(c);
    let b = p(norms[0]) * amplification(angular_range, 10, w)?;
    let bv = p(norms[1]) * amplification(angular_range, 9, w)?;
    let bvv = p(norms[2]) * amplification(angular_range, 8, w)?;
    chain_bound(r, q, lead, [b, bv, bvv], w)
}
fn chain_bound(
    r: ScalarJet,
    q: ScalarJet,
    lead: I,
    bounds: [I; 3],
    w: &mut Work<'_>,
) -> Result<[f64; 2], Error> {
    w.charge(512)?;
    let [b, bv, bvv] = bounds;
    let v1 = p(q.first) / lead;
    let v2 = p(q.second) / lead;
    let first = p(r.first) * b + p(r.value.magnitude()) * bv * v1;
    let second = p(r.second) * b
        + p(2.) * p(r.first) * bv * v1
        + p(r.value.magnitude()) * (bvv * v1.square() + bv * v2);
    if !first.valid() || !second.valid() {
        return Err(Error::Unresolved);
    }
    Ok([first.hi, second.hi])
}
impl<'f, 's> CapCurvePiece<'f, 's> {
    pub fn restrict_finite<'c>(
        &'c self,
        options: CurveOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<RestrictedCapCurve<'c, 'f, 's>, Error> {
        options.validate()?;
        let mut w = Work::new(options.work, interrupt)?;
        w.charge(self.field.work_used() + self.work_used)?;
        let base = self.field.retained_bytes()
            + self.retained_bytes
            + size_of::<RestrictedCapCurve<'_, '_, '_>>();
        let mut memory = Memory::new(base, options.bytes)?;
        memory.add(SCRATCH)?;
        let field = self.field;
        let row = &field.meridians[field.cells[self.source_cell].meridian];
        let s = field.source.source_boundary().source();
        let pitch = binary(
            if row.component == Component::Metric {
                s.metric().spec().pitch_mm
            } else {
                s.wood().spec().pitch_mm
            },
            &mut w,
        )?;
        let phase = self.phase.each_ref().map(|x| x.number(&mut w));
        let [p0, p1] = phase;
        let p0 = p0?;
        let p1 = p1?;
        let (a, ea) = chosen(row, &p0, self.parameter_bounds[0], &pitch, &mut w)?;
        let (b, eb) = chosen(row, &p1, self.parameter_bounds[1], &pitch, &mut w)?;
        let mut input = memory.allocate(self.controls.len(), &mut w)?;
        for x in &self.controls {
            input.push([
                x[0].number(&mut w)?,
                x[1].number(&mut w)?,
                x[2].number(&mut w)?,
                x[3].number(&mut w)?,
            ]);
        }
        let restricted = restrict(&input, &a, &b, &mut memory, &mut w)?;
        let bounds = if ea == 0. && eb == 0. {
            RestrictionBounds {
                position_error_mm: 0.,
                derivative_error_mm: 0.,
                meridian_parameter_error: [0., 0.],
            }
        } else {
            let low = binary(self.parameter_bounds[0][0], &mut w)?;
            let high = binary(self.parameter_bounds[1][1], &mut w)?;
            let jet = factor_derivatives(self, &low, &high, &mut memory, &mut w)?;
            let position = p(jet[0]) * p(ea.max(eb));
            let derivative = p(jet[1]) * p(ea.max(eb)) + p(jet[0]) * (p(ea) + p(eb));
            RestrictionBounds {
                position_error_mm: position.hi,
                derivative_error_mm: derivative.hi,
                meridian_parameter_error: [ea, eb],
            }
        };
        if !bounds.position_error_mm.is_finite()
            || !bounds.derivative_error_mm.is_finite()
            || bounds.position_error_mm > options.position_error_mm
            || bounds.derivative_error_mm > options.derivative_error_mm
        {
            return Err(Error::Unresolved);
        }
        let parameter = [
            EncodedScalar::retain(&a, &mut memory, &mut w)?,
            EncodedScalar::retain(&b, &mut memory, &mut w)?,
        ];
        let mut controls = memory.allocate(restricted.len(), &mut w)?;
        for row in &restricted {
            controls.push([
                EncodedScalar::retain(&row[0], &mut memory, &mut w)?,
                EncodedScalar::retain(&row[1], &mut memory, &mut w)?,
                EncodedScalar::retain(&row[2], &mut memory, &mut w)?,
                EncodedScalar::retain(&row[3], &mut memory, &mut w)?,
            ]);
        }
        memory.release(&input);
        memory.release(&restricted);
        drop(input);
        drop(restricted);
        w.charge(1)?;
        Ok(RestrictedCapCurve {
            source: self,
            parameter,
            controls,
            bounds,
            work: w.used,
            bytes: memory.used - SCRATCH,
        })
    }
}
#[cfg(test)]
pub(super) mod tests;
