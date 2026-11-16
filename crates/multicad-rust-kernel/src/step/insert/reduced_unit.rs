//! Audited unit-circle source -> actual closed degree10 coefficient candidates.
//! No manufacturing, imported-file, occupied-material or winding admission.
use super::*;
use domains::N;
use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{One, Signed, Zero};
const BITS: usize = 96;
type AngularPiece = [[BigRational; 2]; 11];
const SCRATCH: usize = 4 * 1024 * 1024;
const LIMBS: usize = 22 * 512;
/// Immutable actual dyadic coefficients. Bounds refer to the ORIGINAL complete
/// unit-turn coordinate v∈[0,1], not each piece's local normalized coordinate.
pub struct ReducedUnitAngular {
    knots: Vec<f64>,
    coefficients: Vec<AngularPiece>,
    position_error: f64,
    derivative_error: f64,
    left: bool,
}
impl ReducedUnitAngular {
    pub fn prepare(
        left: bool,
        options: CurveOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, Error> {
        let mut work = Work::new(options.work, interrupt)?;
        let mut memory = Memory::new(0, options.bytes)?;
        fit_reduced_unit(left, options, &mut memory, &mut work)
    }
    pub fn knots(&self) -> &[f64] {
        &self.knots
    }
    pub fn coefficients(&self) -> &[AngularPiece] {
        &self.coefficients
    }
    pub fn position_error(&self) -> f64 {
        self.position_error
    }
    pub fn derivative_error(&self) -> f64 {
        self.derivative_error
    }
    pub fn left_hand(&self) -> bool {
        self.left
    }
    pub(crate) fn release_memory(&self, memory: &mut Memory) {
        memory.release(&self.knots);
        memory.release(&self.coefficients);
        memory.used -= self.coefficients.len() * LIMBS + size_of::<Self>();
    }
}
fn round(n: BigInt, d: &BigInt) -> BigInt {
    let negative = n.is_negative();
    let absolute = n.abs();
    let q = &absolute / d;
    let r = &absolute % d;
    let twice = &r << 1;
    let q = if twice > *d || (twice == *d && (&q % BigInt::from(2)) != BigInt::zero()) {
        q + 1
    } else {
        q
    };
    if negative { -q } else { q }
}
fn binary(value: f64, w: &mut Work<'_>) -> Result<BigInt, Error> {
    let value = BigRational::from_float(value).expect("audited finite angular controls");
    // Every finite binary64 numerator/denominator has<=1075bits. This fixed96
    // shift has a derived<2200bit preflight; no unearned arbitrary-rational DTO.
    let words =
        ((value.numer().bits() + BITS as u64).max(value.denom().bits()) as usize).div_ceil(64);
    w.charge(1 + words * words)?;
    Ok(round(value.numer() << BITS, value.denom()))
}
fn upper(value: BigRational, w: &mut Work<'_>) -> Result<f64, Error> {
    Ok(N(value).enclosure(w)?.hi)
}
fn reduce(
    parent: &AngularCarrier,
    left: bool,
    memory: &mut Memory,
    w: &mut Work<'_>,
) -> Result<(Vec<AngularPiece>, [f64; 2]), Error> {
    memory.add(SCRATCH)?;
    let mut pieces = memory.allocate(parent.knots.len() - 1, w)?;
    let mut error = [0_f64; 2];
    for (piece, input) in parent.poles.windows(12).step_by(11).enumerate() {
        let mut q: [[BigInt; 2]; 11] =
            std::array::from_fn(|_| std::array::from_fn(|_| BigInt::zero()));
        for axis in 0..2 {
            q[0][axis] = binary(input[0][axis], w)?;
            for i in 1..10 {
                let bits = binary(input[i][axis], w)?;
                let words = (bits.bits().max(q[i - 1][axis].bits()) as usize + 8).div_ceil(64);
                w.charge(8 * (1 + words * words))?;
                q[i][axis] = round(bits * 11 - &q[i - 1][axis] * i, &BigInt::from(11 - i));
            }
            q[10][axis] = binary(input[11][axis], w)?;
        }
        // Exact authored unit-turn identities, including closure. Actual
        // endpoint changes are INCLUDED in the elevated-net residual below.
        for (index, t) in [(0, parent.knots[piece]), (10, parent.knots[piece + 1])] {
            if (t * 4.).fract() == 0. {
                let mut xy = match (t * 4.) as u8 {
                    0 | 4 => [1., 0.],
                    1 => [0., 1.],
                    2 => [-1., 0.],
                    _ => [0., -1.],
                };
                if left {
                    xy[1] = -xy[1];
                }
                q[index] = [binary(xy[0], w)?, binary(xy[1], w)?];
            }
        }
        let mut residual: [[BigRational; 2]; 12] =
            std::array::from_fn(|_| std::array::from_fn(|_| BigRational::zero()));
        let mut position = BigRational::zero();
        let mut derivative = BigRational::zero();
        for i in 0..12 {
            for axis in 0..2 {
                w.charge(128)?;
                let numerator = if i == 0 {
                    &q[0][axis] * 11
                } else if i == 11 {
                    &q[10][axis] * 11
                } else {
                    &q[i][axis] * (11 - i) + &q[i - 1][axis] * i
                };
                residual[i][axis] = BigRational::new(numerator, BigInt::from(11) << BITS)
                    - BigRational::from_float(input[i][axis]).expect("audited finite control");
            }
            position = position.max(residual[i][0].abs() + residual[i][1].abs());
        }
        let h = BigRational::from_float(parent.knots[piece + 1]).expect("finite knots")
            - BigRational::from_float(parent.knots[piece]).expect("finite knots");
        for pair in residual.windows(2) {
            derivative = derivative.max(
                ((pair[1][0].clone() - &pair[0][0]).abs()
                    + (pair[1][1].clone() - &pair[0][1]).abs())
                    * BigInt::from(11)
                    / &h,
            );
        }
        error[0] = error[0].max(upper(position, w)?);
        error[1] = error[1].max(upper(derivative, w)?);
        memory.add(LIMBS)?;
        w.charge(22 * 128)?;
        pieces.push(q.map(|point| point.map(|n| BigRational::new(n, BigInt::one() << BITS))));
    }
    memory.used -= SCRATCH;
    Ok((
        pieces,
        [
            (I::point(parent.position_error) + I::point(error[0])).hi,
            (I::point(parent.derivative_error) + I::point(error[1])).hi,
        ],
    ))
}
/// Shared Work/Memory remain the authoritative aggregate limits; this function
/// never resets either ledger. Options set the desired source correspondence.
pub(crate) fn fit_reduced_unit(
    left: bool,
    options: CurveOptions,
    memory: &mut Memory,
    w: &mut Work<'_>,
) -> Result<ReducedUnitAngular, Error> {
    options.validate()?;
    w.charge(1)?;
    memory.add(size_of::<ReducedUnitAngular>())?;
    let mut wanted = CurveOptions {
        position_error_mm: options.position_error_mm / 2.,
        derivative_error_mm: options.derivative_error_mm / 2.,
        ..options
    };
    loop {
        if wanted.position_error_mm <= 0. || wanted.derivative_error_mm <= 0. {
            return Err(Error::Unresolved);
        }
        let parent = fit_angular([0., 1.], left, wanted, memory, w)?;
        let (pieces, error) = reduce(&parent, left, memory, w)?;
        memory.release(&parent.poles);
        drop(parent.poles);
        let candidate = ReducedUnitAngular {
            knots: parent.knots,
            coefficients: pieces,
            position_error: error[0],
            derivative_error: error[1],
            left,
        };
        if error[0] <= options.position_error_mm && error[1] <= options.derivative_error_mm {
            return Ok(candidate);
        }
        candidate.release_memory(memory);
        drop(candidate);
        memory.add(size_of::<ReducedUnitAngular>())?;
        wanted.position_error_mm /= 4.;
        wanted.derivative_error_mm /= 4.;
    }
}
#[cfg(test)]
mod tests;
