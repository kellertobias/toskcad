//! Measured reduction of the immutable audited angular carrier. Source-known.
use super::number::{Grid, rounded_ratio, upper};
use super::*;
const BITS: usize = 96;
pub(super) struct Reduced {
    pub poles: Vec<[[Grid; 2]; 11]>,
    pub errors: [f64; 2],
}
fn div_round(value: BigInt, divisor: usize) -> Grid {
    Grid(rounded_ratio(value, &BigInt::from(divisor)))
}
fn scalar(value: f64, w: &mut Work<'_>) -> Result<Grid, Error> {
    Grid::rational(&N(BigRational::from_float(value).unwrap()), BITS, w)
}
pub(super) fn reduce(
    parent: &crate::step::insert::angular::AngularCarrier,
    left: bool,
    memory: &mut Memory,
    w: &mut Work<'_>,
) -> Result<Reduced, Error> {
    memory.add(4 * 1024 * 1024)?;
    let mut poles = memory.allocate(parent.knots.len() - 1, w)?;
    let mut error = [I::point(0.); 2];
    for (piece, input) in parent.poles.windows(12).step_by(11).enumerate() {
        let mut q: [[Grid; 2]; 11] =
            std::array::from_fn(|_| std::array::from_fn(|_| Grid(BigInt::zero())));
        for axis in 0..2 {
            q[0][axis] = scalar(input[0][axis], w)?;
            for i in 1..10 {
                w.charge(16)?;
                q[i][axis] = div_round(
                    scalar(input[i][axis], w)?.0 * 11 - &q[i - 1][axis].0 * i,
                    11 - i,
                );
            }
            q[10][axis] = scalar(input[11][axis], w)?;
        }
        // Actual exact-turn identity. Both incident pieces reference the same
        // cardinal endpoint, and the measured reduction includes this change.
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
                q[index] = [scalar(xy[0], w)?, scalar(xy[1], w)?];
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
                    &q[0][axis].0 * 11
                } else if i == 11 {
                    &q[10][axis].0 * 11
                } else {
                    &q[i][axis].0 * (11 - i) + &q[i - 1][axis].0 * i
                };
                residual[i][axis] = BigRational::new(numerator, BigInt::from(11) << BITS)
                    - BigRational::from_float(input[i][axis]).unwrap();
            }
            position = position.max(residual[i][0].abs() + residual[i][1].abs());
        }
        let h = BigRational::from_float(parent.knots[piece + 1]).unwrap()
            - BigRational::from_float(parent.knots[piece]).unwrap();
        for pair in residual.windows(2) {
            let value = ((pair[1][0].clone() - &pair[0][0]).abs()
                + (pair[1][1].clone() - &pair[0][1]).abs())
                * BigInt::from(11)
                / &h;
            derivative = derivative.max(value);
        }
        error[0].hi = error[0].hi.max(upper(&position, w)?);
        error[1].hi = error[1].hi.max(upper(&derivative, w)?);
        memory.add(22 * 128)?;
        poles.push(q);
    }
    memory.used -= 4 * 1024 * 1024;
    Ok(Reduced {
        poles,
        errors: [
            (I::point(parent.position_error) + error[0]).hi,
            (I::point(parent.derivative_error) + error[1]).hi,
        ],
    })
}
