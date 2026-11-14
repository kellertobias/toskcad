//! Finite Cartesian pole encoding of a retained rational cap support.
//! All weights share ONE exact positive scale; no weight-fit error is omitted.
use super::*;
use num_bigint::BigInt;
use num_traits::{One, Signed};
const FRACTION_BITS: usize = 128;
pub struct LiteralCapCurve<'c, 'f, 's> {
    source: &'c CapCurvePiece<'f, 's>,
    poles: Vec<[EncodedScalar; 3]>,
    weights: Vec<EncodedScalar>,
    error: [f64; 2],
    work: usize,
    bytes: usize,
}
impl LiteralCapCurve<'_, '_, '_> {
    pub fn source(&self) -> &CapCurvePiece<'_, '_> {
        self.source
    }
    pub fn poles(&self) -> &[[EncodedScalar; 3]] {
        &self.poles
    }
    pub fn weights(&self) -> &[EncodedScalar] {
        &self.weights
    }
    pub fn position_error_mm(&self) -> f64 {
        self.error[0]
    }
    /// Whole ORIGINAL root/flat meridian-u derivative, before finite restriction.
    pub fn derivative_error_mm(&self) -> f64 {
        self.error[1]
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_bytes(&self) -> usize {
        self.bytes
    }
}
fn quantize(x: &N, w: &mut Work<'_>) -> Result<N, Error> {
    let scale = N(BigRational::from_integer(BigInt::one() << FRACTION_BITS));
    let scaled = x.mul(&scale, w)?;
    let n = scaled.0.numer();
    let d = scaled.0.denom();
    let words = n.bits().max(d.bits()).div_ceil(64) as usize;
    // Two integer divisions, remainder comparison and retained normalization.
    w.charge(8 * (1 + words * words))?;
    let mut integer = n / d;
    let remainder = n % d;
    if remainder.abs() * 2 >= *d {
        integer += if n.sign() == Sign::Minus { -1 } else { 1 };
    }
    Ok(N(BigRational::new(integer, scale.0.to_integer())))
}
fn finite_weights(h: &[[N; 4]], m: &mut Memory, w: &mut Work<'_>) -> Result<Vec<N>, Error> {
    let mut factor = N::integer(1);
    for row in h {
        if row[3].cmp(&N::integer(0), w)? != Ordering::Greater {
            return Err(Error::Unresolved);
        }
        w.charge(512)?;
        factor = factor.mul(&N(BigRational::from_integer(row[3].0.denom().clone())), w)?;
    }
    let mut maximum = 0;
    for row in h {
        maximum = maximum.max(row[3].mul(&factor, w)?.0.numer().bits());
    }
    // Product of ALL denominators makes every scaled weight an integer.
    // Uniform dyadic normalization changes neither the rational curve nor its
    // quotient derivatives, and keeps the largest encoded weight in[1/2,1).
    w.charge(1 + (maximum as usize).div_ceil(64))?;
    let scale = N(BigRational::from_integer(BigInt::one() << maximum));
    let mut out = m.allocate(h.len(), w)?;
    for row in h {
        out.push(row[3].mul(&factor, w)?.div(&scale, w)?);
    }
    Ok(out)
}
impl<'f, 's> CapCurvePiece<'f, 's> {
    /// Every endpoint override is included in the SAME quotient error proof.
    /// It must later be earned by the finite shared-vertex/trim registry; this
    /// support query alone does not prove endpoint or BRep incidence.
    pub fn literal<'c>(
        &'c self,
        endpoints: Option<&[[EncodedScalar; 3]; 2]>,
        options: CurveOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<LiteralCapCurve<'c, 'f, 's>, Error> {
        options.validate()?;
        let mut w = Work::new(options.work, interrupt)?;
        w.charge(self.field.work_used() + self.work_used)?;
        let endpoint_bytes = endpoints.map_or(0, |values| {
            size_of_val(values)
                + values
                    .iter()
                    .flatten()
                    .map(|v| v.n.len() + v.d.len())
                    .sum::<usize>()
        });
        let base = endpoint_bytes
            + self.field.retained_bytes()
            + self.retained_bytes
            + size_of::<LiteralCapCurve<'_, '_, '_>>();
        let mut memory = Memory::new(base, options.bytes)?;
        // <=23 four-coordinate source controls, poles, differences, and <=23
        // denominator products: <512 simultaneously live <=16384-bit values.
        memory.add(SCRATCH)?;
        let mut h = memory.allocate(self.controls.len(), &mut w)?;
        for row in &self.controls {
            h.push([
                row[0].number(&mut w)?,
                row[1].number(&mut w)?,
                row[2].number(&mut w)?,
                row[3].number(&mut w)?,
            ]);
        }
        let plane = h[0][2].div(&h[0][3], &mut w)?;
        // CapCurvePiece structurally owns Z=plane*W. Preserve that exact
        // dyadic station, including subnormal positive plain-length offsets.
        // Endpoint corrections may never turn this plane into a fitted slab.
        if let Some(points) = endpoints {
            for point in points {
                if point[2].number(&mut w)?.cmp(&plane, &mut w)? != Ordering::Equal {
                    return Err(Error::Unresolved);
                }
            }
        }
        let weight = finite_weights(&h, &mut memory, &mut w)?;
        let mut poles = memory.allocate(h.len(), &mut w)?;
        let mut weights = memory.allocate(h.len(), &mut w)?;
        let mut delta = memory.allocate(h.len(), &mut w)?;
        let mut minimum = f64::INFINITY;
        let mut error = 0_f64;
        for (i, row) in h.iter().enumerate() {
            let mut chosen: [N; 3] = std::array::from_fn(|_| N::integer(0));
            let mut diff: [N; 3] = std::array::from_fn(|_| N::integer(0));
            let endpoint = if i == 0 {
                Some(0)
            } else if i + 1 == h.len() {
                Some(1)
            } else {
                None
            };
            let mut sum = p(0.);
            for axis in 0..3 {
                let actual = row[axis].div(&row[3], &mut w)?;
                chosen[axis] = if axis == 2 {
                    w.charge(512)?;
                    plane.clone()
                } else {
                    quantize(
                        &if let (Some(e), Some(v)) = (endpoint, endpoints) {
                            v[e][axis].number(&mut w)?
                        } else {
                            actual
                        },
                        &mut w,
                    )?
                };
                diff[axis] = row[axis].sub(&row[3].mul(&chosen[axis], &mut w)?, &mut w)?;
                sum = sum + p(diff[axis].enclosure(&mut w)?.magnitude());
            }
            minimum = minimum.min(row[3].enclosure(&mut w)?.lo);
            error = error.max(sum.hi);
            poles.push([
                EncodedScalar::retain(&chosen[0], &mut memory, &mut w)?,
                EncodedScalar::retain(&chosen[1], &mut memory, &mut w)?,
                EncodedScalar::retain(&chosen[2], &mut memory, &mut w)?,
            ]);
            weights.push(EncodedScalar::retain(&weight[i], &mut memory, &mut w)?);
            delta.push(diff);
        }
        let mut derivative = 0_f64;
        let mut weight_derivative = 0_f64;
        for (pair, weights) in delta.windows(2).zip(h.windows(2)) {
            let mut sum = p(0.);
            for (upper, lower) in pair[1].iter().zip(&pair[0]) {
                sum = sum
                    + p(upper.sub(lower, &mut w)?.enclosure(&mut w)?.magnitude())
                        * p(self.degree() as f64);
            }
            derivative = derivative.max(sum.hi);
            weight_derivative = weight_derivative.max(
                (weights[1][3]
                    .sub(&weights[0][3], &mut w)?
                    .enclosure(&mut w)?
                    * p(self.degree() as f64))
                .magnitude(),
            );
        }
        let position = p(error) / p(minimum);
        let tangent =
            p(derivative) / p(minimum) + p(error) * p(weight_derivative) / p(minimum).square();
        if minimum <= 0.
            || !position.valid()
            || !tangent.valid()
            || position.hi > options.position_error_mm
            || tangent.hi > options.derivative_error_mm
        {
            return Err(Error::Unresolved);
        }
        memory.release(&h);
        memory.release(&weight);
        memory.release(&delta);
        drop(h);
        drop(weight);
        drop(delta);
        w.charge(1)?;
        Ok(LiteralCapCurve {
            source: self,
            poles,
            weights,
            error: [position.hi, tangent.hi],
            work: w.used,
            bytes: memory.used - SCRATCH,
        })
    }
}
#[cfg(test)]
mod tests;
