//! Actual 3D cap curves composed from the SAME exact scalar and angular owners.
//! These homogeneous supports are not rounded STEP Cartesian poles or a file.
use super::*;
mod affine;
pub use affine::AffineCapPiece;
#[cfg(test)]
mod affine_probe;
mod bound_affine;
mod bound_root;
pub use bound_root::{BoundRootCap, EndpointArithmetic};
mod correction;
mod jets;
pub use bound_affine::BoundAffineCap;
mod literal;
mod side;
pub use side::{SidePole, SideSupport, SideSupportBatch, SideSupportKind};
mod restrict;
mod vertices;
pub use vertices::{CapEvent, CapEvents, CapSegment, LiteralVertex, LiteralVertices};
#[cfg(test)]
mod side_probe;
#[cfg(test)]
mod trim_probe;
#[cfg(test)]
mod unit_probe;
pub use literal::LiteralCapCurve;
pub use restrict::{RestrictedCapCurve, RestrictionBounds};
use std::cmp::Ordering;
const DEGREE: usize = 22;
struct Polynomial {
    c: [N; DEGREE + 1],
    degree: usize,
}
impl Polynomial {
    fn zero(degree: usize) -> Self {
        Self {
            c: std::array::from_fn(|_| N::integer(0)),
            degree,
        }
    }
    fn one() -> Self {
        let mut p = Self::zero(0);
        p.c[0] = N::integer(1);
        p
    }
    fn coefficients(c: &[N], w: &mut Work<'_>) -> Result<Self, Error> {
        if c.is_empty() || c.len() > DEGREE + 1 {
            return Err(Error::Invalid);
        }
        w.charge(1 + c.len() * 512)?;
        let mut p = Self::zero(c.len() - 1);
        p.c[..c.len()].clone_from_slice(c);
        Ok(p)
    }
    fn product(&self, b: &Self, w: &mut Work<'_>) -> Result<Self, Error> {
        let n = self.degree + b.degree;
        if n > DEGREE {
            return Err(Error::Limit);
        }
        let mut out = Self::zero(n);
        for i in 0..=self.degree {
            for j in 0..=b.degree {
                let factor = N::integer(choose(self.degree, i, w)? * choose(b.degree, j, w)?)
                    .div(&N::integer(choose(n, i + j, w)?), w)?;
                let term = self.c[i].mul(&b.c[j], w)?.mul(&factor, w)?;
                out.c[i + j] = out.c[i + j].add(&term, w)?;
            }
        }
        Ok(out)
    }
}
fn choose(n: usize, k: usize, w: &mut Work<'_>) -> Result<i64, Error> {
    w.charge(n + 1)?;
    let mut out = 1_i64;
    for i in 0..k.min(n - k) {
        out = out * (n - i) as i64 / (i + 1) as i64;
    }
    Ok(out)
}
/// Rational Bezier support, parameterized by the ORIGINAL meridian u.
/// Only the stated finite u interval uses the selected angular polynomial.
/// Exact homogeneous controls must be converted to Cartesian poles/weights
/// with a separately earned literal encoding bound before serialization.
pub struct CapCurvePiece<'f, 's> {
    field: &'f EncodedField<'s>,
    surface_trace: CapTrace,
    source_cell: usize,
    angular_piece: usize,
    end: End,
    phase: [EncodedScalar; 2],
    parameter_bounds: [[f64; 2]; 2],
    controls: Vec<[EncodedScalar; 4]>,
    work_used: usize,
    retained_bytes: usize,
}
impl CapCurvePiece<'_, '_> {
    pub fn field(&self) -> &EncodedField<'_> {
        self.field
    }
    pub fn surface_trace(&self) -> &CapTrace {
        &self.surface_trace
    }
    pub fn source_cell(&self) -> usize {
        self.source_cell
    }
    pub fn angular_piece(&self) -> usize {
        self.angular_piece
    }
    pub fn end(&self) -> End {
        self.end
    }
    pub fn phase(&self) -> &[EncodedScalar; 2] {
        &self.phase
    }
    pub fn parameter_bounds(&self) -> [[f64; 2]; 2] {
        self.parameter_bounds
    }
    /// Query-local work; borrowed EncodedField preparation is separate.
    pub fn work_used(&self) -> usize {
        self.work_used
    }
    /// Owned curve storage; borrowed EncodedField storage is separate.
    pub fn retained_bytes(&self) -> usize {
        self.retained_bytes
    }
    pub fn degree(&self) -> usize {
        self.controls.len() - 1
    }
    pub fn controls(&self) -> &[[EncodedScalar; 4]] {
        &self.controls
    }
}
pub(super) fn parameter(
    row: &EncodedMeridian,
    pitch: &N,
    phase: &N,
    w: &mut Work<'_>,
) -> Result<I, Error> {
    let a = row.phase[0].number(w)?.div(&row.weights[0].number(w)?, w)?;
    let b = row.phase[row.degree()]
        .number(w)?
        .div(&row.weights[row.degree()].number(w)?, w)?;
    if phase.cmp(&a, w)? == Ordering::Equal {
        return Ok(p(0.));
    }
    if phase.cmp(&b, w)? == Ordering::Equal {
        return Ok(p(1.));
    }
    let value = if row.component == Component::Metric && row.band == ThreadBandRole::Root {
        let d = phase
            .div(pitch, w)?
            .sub(&N::integer(1).div(&N::integer(2), w)?, w)?
            .mul(&N::integer(8), w)?
            .enclosure(w)?;
        w.charge(64)?;
        (p(1.) + p(3.) * d / (p(2.) + (p(4.) - p(3.) * d.square()).sqrt())) / p(2.)
    } else {
        phase.sub(&a, w)?.div(&b.sub(&a, w)?, w)?.enclosure(w)?
    };
    // Exact phase has already been intersected with this strictly monotone
    // row. Same positive-root inverse as CapTrace, with all N conversion
    // failures propagated before the earned [0,1] intersection.
    Ok(I::new(value.lo.max(0.), value.hi.min(1.)))
}
impl<'s> EncodedField<'s> {
    pub fn cap_curve_piece<'f>(
        &'f self,
        index: usize,
        end: End,
        angular_piece: usize,
        options: DomainOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Option<CapCurvePiece<'f, 's>>, Error> {
        let (mut w, mut memory) = super::super::domain::ledgers(options, self.bytes, interrupt)?;
        memory.add(SCRATCH - super::super::SCRATCH)?;
        let mut output =
            self.cap_curve_piece_with_ledgers(index, end, angular_piece, &mut memory, &mut w)?;
        if let Some(curve) = &mut output {
            curve.work_used = w.used;
        }
        Ok(output)
    }
    /// Aggregate callers retain the same field and reserve SCRATCH once.
    /// Failed/empty selections release only their own memory; spent work
    /// remains debited. Successful receipts report this query's increment.
    pub(super) fn cap_curve_piece_with_ledgers<'f>(
        &'f self,
        index: usize,
        end: End,
        angular_piece: usize,
        memory: &mut Memory,
        w: &mut Work<'_>,
    ) -> Result<Option<CapCurvePiece<'f, 's>>, Error> {
        let base = memory.used;
        let result = self.compute_cap_piece(index, end, angular_piece, memory, w);
        if !matches!(&result, Ok(Some(_))) {
            memory.used = base;
        }
        result
    }
    fn compute_cap_piece<'f>(
        &'f self,
        index: usize,
        end: End,
        angular_piece: usize,
        memory: &mut Memory,
        w: &mut Work<'_>,
    ) -> Result<Option<CapCurvePiece<'f, 's>>, Error> {
        let start_work = w.used;
        let start_memory = memory.used;
        memory.add(size_of::<CapCurvePiece>())?;
        if angular_piece >= self.angular.coefficients().len() {
            return Err(Error::Invalid);
        }
        if w.remaining() == 0 {
            return Err(Error::Limit);
        }
        let Some(trace) = self.cap_trace_with_ledgers(index, end, memory, w)? else {
            w.charge(1)?;
            return Ok(None);
        };
        // CapCurvePiece already contains the moved CapTrace descriptor.
        memory.used -= size_of::<CapTrace>();
        let cell = &self.cells[index];
        let row = &self.meridians[cell.meridian];
        let s = self.source.source_boundary().source();
        let (pitch, starts, left) = if row.component == Component::Metric {
            let t = s.metric().spec();
            (t.pitch_mm, t.starts, t.left_hand)
        } else {
            let t = s.wood().spec();
            (t.pitch_mm, t.starts, t.left_hand)
        };
        let pitch = binary(pitch, w)?;
        let lead = pitch.mul(&N::integer(i64::from(starts)), w)?;
        let station = if end == End::Start {
            N::integer(0)
        } else {
            binary(
                if row.component == Component::Metric {
                    s.metric().full_length_mm()
                } else {
                    s.wood().full_length_mm()
                },
                w,
            )?
        };
        let shifted = station.sub(
            &pitch.mul(&N::integer(i64::from(cell.source_key.total_period)), w)?,
            w,
        )?;
        let a = binary(self.angular.knots()[angular_piece], w)?;
        let b = binary(self.angular.knots()[angular_piece + 1], w)?;
        let lower = shifted.sub(&lead.mul(&b, w)?, w)?;
        let upper = shifted.sub(&lead.mul(&a, w)?, w)?;
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
        if lo.cmp(&hi, w)? != Ordering::Less {
            w.charge(1)?;
            return Ok(None);
        }
        let u0 = parameter(row, &pitch, &lo, w)?;
        let u1 = parameter(row, &pitch, &hi, w)?;
        if u0.hi >= u1.lo {
            return Err(Error::Unresolved);
        }
        let q = row
            .phase
            .iter()
            .map(|x| x.number(w))
            .collect::<Result<Vec<_>, _>>()?;
        let radius = row
            .radius
            .iter()
            .map(|x| x.number(w))
            .collect::<Result<Vec<_>, _>>()?;
        let weights = row
            .weights
            .iter()
            .map(|x| x.number(w))
            .collect::<Result<Vec<_>, _>>()?;
        let scale = lead.mul(&b.sub(&a, w)?, w)?;
        let shift = shifted.sub(&lead.mul(&a, w)?, w)?;
        let mut t = Polynomial::zero(row.degree());
        let mut other = Polynomial::zero(row.degree());
        for i in 0..=row.degree() {
            t.c[i] = shift.mul(&weights[i], w)?.sub(&q[i], w)?.div(&scale, w)?;
            other.c[i] = weights[i].sub(&t.c[i], w)?;
        }
        let mut tp: [Polynomial; 11] = std::array::from_fn(|_| Polynomial::one());
        let mut op: [Polynomial; 11] = std::array::from_fn(|_| Polynomial::one());
        for i in 1..=10 {
            tp[i] = tp[i - 1].product(&t, w)?;
            op[i] = op[i - 1].product(&other, w)?;
        }
        let mut xy = [
            Polynomial::zero(10 * row.degree()),
            Polynomial::zero(10 * row.degree()),
        ];
        let (piece, reverse) = if left {
            (self.angular.coefficients().len() - 1 - angular_piece, true)
        } else {
            (angular_piece, false)
        };
        for i in 0..=10 {
            let term = tp[i].product(&op[10 - i], w)?;
            let coefficient = &self.angular.coefficients()[piece][if reverse { 10 - i } else { i }];
            for axis in 0..2 {
                w.charge(512)?;
                let factor = N(coefficient[axis].clone()).mul(&N::integer(choose(10, i, w)?), w)?;
                for j in 0..=term.degree {
                    xy[axis].c[j] = xy[axis].c[j].add(&term.c[j].mul(&factor, w)?, w)?;
                }
            }
        }
        let r = Polynomial::coefficients(&radius, w)?;
        let x = r.product(&xy[0], w)?;
        let y = r.product(&xy[1], w)?;
        let weight = Polynomial::coefficients(&weights, w)?;
        let mut den = Polynomial::one();
        for _ in 0..11 {
            den = den.product(&weight, w)?;
        }
        let plane = trace.plane_z.number(w)?;
        let mut controls = memory.allocate(den.degree + 1, w)?;
        for i in 0..=den.degree {
            controls.push([
                EncodedScalar::retain(&x.c[i], memory, w)?,
                EncodedScalar::retain(&y.c[i], memory, w)?,
                EncodedScalar::retain(&plane.mul(&den.c[i], w)?, memory, w)?,
                EncodedScalar::retain(&den.c[i], memory, w)?,
            ]);
        }
        let phase = [
            EncodedScalar::retain(&lo, memory, w)?,
            EncodedScalar::retain(&hi, memory, w)?,
        ];
        w.charge(1)?;
        Ok(Some(CapCurvePiece {
            field: self,
            surface_trace: trace,
            source_cell: index,
            angular_piece,
            end,
            phase,
            parameter_bounds: [[u0.lo, u0.hi], [u1.lo, u1.hi]],
            controls,
            work_used: w.used - start_work,
            retained_bytes: memory.used - start_memory,
        }))
    }
}
#[cfg(test)]
mod tests;
