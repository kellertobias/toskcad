//! Exact Machin bounds and gcd-free reduction for a rigorous rotation.
//! Each pi bound is the alternating32-term rational Machin construction;
//! tests recompute that identity independently. No libm trigonometric call
//! or heuristic ulp margin enters the source correspondence proof.
use super::*;
use num_bigint::BigInt as Z;
use num_traits::{Signed, ToPrimitive};
const PI_LOW_N: &str = "174118290308941112378113826430623356189660704808346636107620488511958189347105442573891134318947595085642963956171759213315555515609276797618972308309360496669357217484010844413919575812751533622854219483882629877787022476";
const PI_LOW_D: &str = "55423573170755267764678847670127618640928415844293113986805896510622486333639072176235601943412170410931973377118988209489449832305322773611681329213437134813548110027581192184201012984612333411860163323581218719482421875";
const PI_HIGH_N: &str = "381029503836025963258070207171231587747196094285332700534943760868988539549066611527487739149816931391630886753632848466806557368935301121273928292086515952217380320479815457699307436447121321805533860996086486663952944";
const PI_HIGH_D: &str = "121285457998711655443442095880078295725145777919445374701961398852047597060711192416614734387117195101039839501056940988186152711580072945177083141956191975835393528710065457940601996167373499702080152928829193115234375";
pub(super) struct Rotation {
    bounds: [(Z, Z); 2],
}
fn fraction(v: f64) -> (Z, u32) {
    let bits = v.to_bits();
    let e = ((bits >> 52) & 2047) as u32;
    let mantissa = (bits & ((1_u64 << 52) - 1)) | if e == 0 { 0 } else { 1_u64 << 52 };
    let m = Z::from(mantissa);
    let m = if bits >> 63 == 0 { m } else { -m };
    // Admitted abs<=1e7 gives a negative binary exponent, including subnormals.
    (m, if e == 0 { 1074 } else { 1075 - e })
}
fn rational_bound(n: &Z, d: &Z) -> I {
    // abs(n/d)<=4. Fixed2^-60 flooring gives a rigorous enclosure without
    // division by floating-point approximations of enormous integer operands.
    let q = ((n.abs() << 60_u32) / d)
        .to_u64()
        .expect("proved bounded fixed-point quotient");
    assert!(q <= 4_u64 << 60);
    let factor = 2_f64.powi(-60);
    let a = ((q as f64) * factor).next_down().max(0.);
    let b = (((q + 1) as f64) * factor).next_up();
    if n.is_negative() {
        I::new(-b, -a)
    } else {
        I::new(a, b)
    }
}
impl Rotation {
    pub(super) fn new(work: &mut Work<'_>) -> Result<Self> {
        work.step(4096)?;
        Ok(Self {
            bounds: [
                (
                    Z::parse_bytes(PI_LOW_N.as_bytes(), 10).expect("fixed exact pi numerator"),
                    Z::parse_bytes(PI_LOW_D.as_bytes(), 10).expect("fixed exact pi denominator"),
                ),
                (
                    Z::parse_bytes(PI_HIGH_N.as_bytes(), 10).expect("fixed exact pi numerator"),
                    Z::parse_bytes(PI_HIGH_D.as_bytes(), 10).expect("fixed exact pi denominator"),
                ),
            ],
        })
    }
    pub(super) fn sin_cos(&self, angle: I, work: &mut Work<'_>) -> Result<[I; 2]> {
        //735-bit constants +binary64 denominator<=2^1074 keep all live
        // numerators/denominators below2048bits. Two bounded integer divisions
        // and fixed degree31 interval arithmetic fit4096 geometry primitives.
        // This is not an instruction/RSS claim; caller reserves fixed scratch.
        work.step(4096)?;
        if !angle.valid() || angle.magnitude() > 1e7 {
            return Err(ReadError::Unsupported);
        }
        if angle.lo == 0. && angle.hi == 0. {
            return Ok([I::point(1.), I::point(0.)]);
        }
        let turns = (angle.midpoint() / std::f64::consts::TAU).round() as i64;
        let reduce = |v: f64, use_upper_pi: bool| {
            let (m, shift) = fraction(v);
            let (n, d) = &self.bounds[usize::from(use_upper_pi)];
            let numerator = m * d - ((Z::from(2 * turns) * n) << shift);
            let denominator = d << shift;
            if numerator.abs() > 4 * &denominator {
                return Err(ReadError::Unsupported);
            }
            Ok(rational_bound(&numerator, &denominator))
        };
        let a = reduce(angle.lo, turns >= 0)?;
        let b = reduce(angle.hi, turns < 0)?;
        let r = I::new(a.lo, b.hi);
        // The exact rational reduction already earned abs<=4. Outward
        // binary conversion may widen by one endpoint ulp;4.1 is a proven
        // finite arithmetic bound, not an acceptance tolerance.
        debug_assert!(r.valid() && r.magnitude() <= 4.1);
        let mut term = [I::point(1.), I::point(0.)];
        let mut value = term;
        for k in 1..=31 {
            term = [
                -term[1] * r / I::point(k as f64),
                term[0] * r / I::point(k as f64),
            ];
            for axis in 0..2 {
                value[axis] = value[axis] + term[axis];
            }
        }
        let mut remainder = I::point(1.);
        for k in 1..=32 {
            remainder = remainder * I::point(r.magnitude()) / I::point(k as f64);
        }
        Ok(value.map(|v| {
            I::new(
                (v - I::point(remainder.hi)).lo.max(-1.),
                (v + I::point(remainder.hi)).hi.min(1.),
            )
        }))
    }
}

#[cfg(test)]
mod tests;
