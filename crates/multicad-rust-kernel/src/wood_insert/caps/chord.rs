//! Whole actual cap source versus each literal affine facet, piece by piece.
use super::*;
type R<T> = Result<T, Error>;
macro_rules! add {
    ($a:expr,$b:expr,$w:expr) => {
        op($a, $b, Op::Add, $w)?
    };
}
macro_rules! sub {
    ($a:expr,$b:expr,$w:expr) => {
        op($a, $b, Op::Sub, $w)?
    };
}
macro_rules! mul {
    ($a:expr,$b:expr,$w:expr) => {
        op($a, $b, Op::Mul, $w)?
    };
}
macro_rules! div {
    ($a:expr,$b:expr,$w:expr) => {
        op($a, $b, Op::Div, $w)?
    };
}
pub(super) fn binary_work(x: f64, w: &mut Work<'_>) -> R<Q> {
    if !x.is_finite() {
        return Err(Error::Unresolved);
    }
    let bits = x.to_bits();
    let exponent = ((bits >> 52) & 0x7ff) as i32;
    let mantissa = bits & ((1_u64 << 52) - 1);
    if exponent == 0 && mantissa == 0 {
        w.charge(1)?;
        return Ok(Q::zero());
    }
    let (mantissa, power) = if exponent == 0 {
        (mantissa, -1074)
    } else {
        (mantissa | (1 << 52), exponent - 1023 - 52)
    };
    let n = u64::from(64 - mantissa.leading_zeros()) + power.max(0) as u64;
    let d = 1 + (-power).max(0) as u64;
    let words = (n + d).div_ceil(64) as usize;
    w.charge(1 + words * words)?;
    Ok(binary(x))
}
fn interval(q: &Q, w: &mut Work<'_>) -> R<I> {
    let words = (q.numer().bits() + q.denom().bits()).div_ceil(64) as usize;
    w.charge(1 + words * words)?;
    let x = q
        .to_f64()
        .filter(|v| v.is_finite())
        .ok_or(Error::Unresolved)?;
    let [lo, hi] = [x.next_down(), x.next_up()];
    if !lo.is_finite() || !hi.is_finite() {
        return Err(Error::Unresolved);
    }
    let a = binary_work(lo, w)?;
    let b = binary_work(hi, w)?;
    if cmp(&a, q, w)? == Ordering::Greater || cmp(&b, q, w)? == Ordering::Less {
        return Err(Error::Unresolved);
    }
    Ok(I::new(lo, hi))
}
fn coefficients(source: &AuthoredThread, z: &Q, band: u8, k: i64, w: &mut Work<'_>) -> R<[I; 2]> {
    // These closed radius formulas use<=16 finite source conversions. Their
    // admitted operands have<=1075bits; reserve the full limb bound first.
    w.charge(16 * 300)?;
    let a = radius(source, z, &Q::zero(), band, k, w)?;
    let b = radius(source, z, &Q::one(), band, k, w)?;
    let b = sub!(&b, &a, w);
    Ok([interval(&a, w)?, interval(&b, w)?])
}
fn affine(uv: &[[Q; 2]; 3], p: [[f64; 3]; 3], w: &mut Work<'_>) -> R<([[I; 3]; 3], bool)> {
    let du = [sub!(&uv[1][0], &uv[0][0], w), sub!(&uv[2][0], &uv[0][0], w)];
    let dt = [sub!(&uv[1][1], &uv[0][1], w), sub!(&uv[2][1], &uv[0][1], w)];
    let determinant = sub!(&mul!(&du[0], &dt[1], w), &mul!(&du[1], &dt[0], w), w);
    if determinant.is_zero() {
        return Err(Error::Unresolved);
    }
    let positive = determinant > Q::zero();
    let mut result = [[I::point(0.); 3]; 3];
    for axis in 0..3 {
        let origin = binary_work(p[0][axis], w)?;
        let dp = [
            sub!(&binary_work(p[1][axis], w)?, &origin, w),
            sub!(&binary_work(p[2][axis], w)?, &origin, w),
        ];
        let a = div!(
            &sub!(&mul!(&dp[0], &dt[1], w), &mul!(&dp[1], &dt[0], w), w),
            &determinant,
            w
        );
        let b = div!(
            &sub!(&mul!(&du[0], &dp[1], w), &mul!(&du[1], &dp[0], w), w),
            &determinant,
            w
        );
        let c = sub!(
            &sub!(&origin, &mul!(&a, &uv[0][0], w), w),
            &mul!(&b, &uv[0][1], w),
            w
        );
        result[axis] = [interval(&a, w)?, interval(&b, w)?, interval(&c, w)?];
    }
    Ok((result, positive))
}
fn cross(a: &[Q; 2], b: &[Q; 2], p: &[Q; 2], w: &mut Work<'_>) -> R<Q> {
    Ok(sub!(
        &mul!(&sub!(&b[0], &a[0], w), &sub!(&p[1], &a[1], w), w),
        &mul!(&sub!(&b[1], &a[1], w), &sub!(&p[0], &a[0], w), w),
        w
    ))
}
fn outside(uv: &[[Q; 2]; 3], bounds: &[[Q; 2]; 2], positive: bool, w: &mut Work<'_>) -> R<bool> {
    for e in 0..3 {
        let mut outside = true;
        for u in &bounds[0] {
            for t in &bounds[1] {
                let c = cross(&uv[e], &uv[(e + 1) % 3], &[u.clone(), t.clone()], w)?;
                if (positive && c >= Q::zero()) || (!positive && c <= Q::zero()) {
                    outside = false;
                }
            }
        }
        if outside {
            return Ok(true);
        }
    }
    Ok(false)
}
fn field(c: [[I; 2]; 2], u: I, t: I, z: I, w: &mut Work<'_>) -> R<[[I; 5]; 3]> {
    // Fixed point/first/Hessian expression count, charged before evaluation.
    w.charge(128)?;
    let tau = I::new(
        std::f64::consts::PI.next_down(),
        std::f64::consts::PI.next_up(),
    ) * I::point(2.);
    let angle = u * tau;
    let co = angle.cos();
    let si = angle.sin();
    let inner = c[0][0] + c[0][1] * u;
    let gap = c[1][0] - c[0][0] + (c[1][1] - c[0][1]) * u;
    let r = inner + t * gap;
    let rw = c[0][1] + t * (c[1][1] - c[0][1]);
    let slope = c[1][1] - c[0][1];
    let out = [
        [
            r * co,
            rw * co - tau * r * si,
            gap * co,
            -I::point(2.) * tau * rw * si - tau.square() * r * co,
            slope * co - tau * gap * si,
        ],
        [
            r * si,
            rw * si + tau * r * co,
            gap * si,
            I::point(2.) * tau * rw * co - tau.square() * r * si,
            slope * si + tau * gap * co,
        ],
        [z, I::point(0.), I::point(0.), I::point(0.), I::point(0.)],
    ];
    if !out.iter().flatten().copied().all(I::valid) {
        return Err(Error::Unresolved);
    }
    Ok(out)
}
fn residual(
    c: [[I; 2]; 2],
    bounds: &[[Q; 2]; 2],
    z: I,
    a: [[I; 3]; 3],
    w: &mut Work<'_>,
) -> R<f64> {
    let mut center = [Q::zero(), Q::zero()];
    let mut ci = [I::point(0.); 2];
    let mut range = ci;
    let mut delta = ci;
    for i in 0..2 {
        center[i] = div!(
            &add!(&bounds[i][0], &bounds[i][1], w),
            &Q::from_integer(2.into()),
            w
        );
        ci[i] = interval(&center[i], w)?;
        let lo = interval(&bounds[i][0], w)?;
        let hi = interval(&bounds[i][1], w)?;
        range[i] = I::new(lo.lo, hi.hi);
        let lo = interval(&sub!(&bounds[i][0], &center[i], w), w)?;
        let hi = interval(&sub!(&bounds[i][1], &center[i], w), w)?;
        delta[i] = I::new(lo.lo, hi.hi);
    }
    let f0 = field(c, ci[0], ci[1], z, w)?;
    let whole = field(c, range[0], range[1], z, w)?;
    w.charge(128)?;
    let mut squares = I::point(0.);
    for axis in 0..3 {
        let l = a[axis][0] * ci[0] + a[axis][1] * ci[1] + a[axis][2];
        let r0 = f0[axis][0] - l;
        let rw = f0[axis][1] - a[axis][0];
        let rt = f0[axis][2] - a[axis][1];
        let b = I::point(r0.magnitude())
            + I::point(rw.magnitude()) * I::point(delta[0].magnitude())
            + I::point(rt.magnitude()) * I::point(delta[1].magnitude())
            + I::point(0.5)
                * I::point(whole[axis][3].magnitude())
                * I::point(delta[0].magnitude()).square()
            + I::point(whole[axis][4].magnitude())
                * I::point(delta[0].magnitude())
                * I::point(delta[1].magnitude());
        squares = squares + b.square();
    }
    let out = squares.sqrt();
    if !out.valid() {
        return Err(Error::Unresolved);
    }
    Ok(out.hi)
}
pub(super) struct CapProof<'a> {
    cap: &'a CapDomains,
    coefficients: Vec<[[I; 2]; 2]>,
    z: I,
}
impl<'a> CapProof<'a> {
    pub fn new(cap: &'a CapDomains, w: &mut Work<'_>) -> R<Self> {
        let z = station(cap.owner(), cap.end());
        let mut cache = reserve(cap.cells().len(), w)?;
        for cell in cap.cells() {
            cache.push([
                coefficients(
                    cap.owner().inner(),
                    &z,
                    cell.inner_band,
                    cell.inner_period,
                    w,
                )?,
                coefficients(
                    cap.owner().outer(),
                    &z,
                    cell.outer_band,
                    cell.outer_period,
                    w,
                )?,
            ]);
        }
        let tau = I::new(
            std::f64::consts::PI.next_down(),
            std::f64::consts::PI.next_up(),
        ) * I::point(2.);
        let jac = tau
            * I::point(cap.owner().inner().minimum_radius_mm())
            * I::point(cap.owner().minimum_wall_mm());
        if !jac.valid() || jac.lo <= 0. {
            return Err(Error::Unresolved);
        }
        Ok(Self {
            cap,
            coefficients: cache,
            z: I::point(if cap.end() == CapEnd::Start {
                0.
            } else {
                cap.owner().full_length_mm()
            }),
        })
    }
    pub fn certify(
        &self,
        uv: &[[Q; 2]; 3],
        p: [[f64; 3]; 3],
        limit: f64,
        w: &mut Work<'_>,
    ) -> R<f64> {
        let (a, positive) = affine(uv, p, w)?;
        let min = uv.iter().map(|p| p[0].clone()).min().unwrap();
        let max = uv.iter().map(|p| p[0].clone()).max().unwrap();
        if min < Q::zero() || max > Q::from_integer(2.into()) {
            return Err(Error::Unresolved);
        }
        let mut maximum = 0_f64;
        let mut checked = 0;
        for shift in 0..=1 {
            let shift_q = Q::from_integer(shift.into());
            for (index, cell) in self.cap.cells().iter().enumerate() {
                let lo =
                    add!(&self.cap.angles[cell.events[0]].decode(w)?, &shift_q, w).max(min.clone());
                let hi =
                    add!(&self.cap.angles[cell.events[1]].decode(w)?, &shift_q, w).min(max.clone());
                if lo >= hi {
                    continue;
                }
                let mut c = self.coefficients[index];
                if shift != 0 {
                    for v in &mut c {
                        v[0] = v[0] - v[1];
                    }
                }
                // True period transfer r(w+1)=r(w), equivalent to k−hand*starts.
                let mut stack = reserve(37, w)?;
                stack.push(([[lo, hi], [Q::zero(), Q::one()]], 0));
                while let Some((bounds, depth)) = stack.pop() {
                    if outside(uv, &bounds, positive, w)? {
                        continue;
                    }
                    checked += 1;
                    if checked > 4096 {
                        return Err(Error::Limit);
                    }
                    let e = residual(c, &bounds, self.z, a, w)?;
                    if e <= limit {
                        maximum = maximum.max(e);
                        continue;
                    }
                    if depth == 12 {
                        return Err(Error::Unresolved);
                    }
                    let mut mid = [Q::zero(), Q::zero()];
                    for i in 0..2 {
                        mid[i] = div!(
                            &add!(&bounds[i][0], &bounds[i][1], w),
                            &Q::from_integer(2.into()),
                            w
                        );
                    }
                    for i in 0..2 {
                        for j in 0..2 {
                            let half = |axis: usize, upper: usize| {
                                if upper == 0 {
                                    [bounds[axis][0].clone(), mid[axis].clone()]
                                } else {
                                    [mid[axis].clone(), bounds[axis][1].clone()]
                                }
                            };
                            stack.push(([half(0, i), half(1, j)], depth + 1));
                        }
                    }
                }
            }
        }
        w.charge(1)?;
        Ok(maximum)
    }
}
