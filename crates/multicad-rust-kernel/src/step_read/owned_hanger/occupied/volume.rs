//! Volume of the healed, decoded support field; never the raw trim occupied set.
use super::*;
#[derive(Clone, Copy)]
pub struct VolumeOptions {
    pub maximum_work: usize,
    pub maximum_bytes: usize,
    pub maximum_width_mm3: f64,
}
impl Default for VolumeOptions {
    fn default() -> Self {
        Self {
            maximum_work: 16_000_000,
            maximum_bytes: 64 * 1024 * 1024,
            maximum_width_mm3: 1e-6,
        }
    }
}
pub struct HealedVolume {
    bounds: I,
    work: usize,
}
impl HealedVolume {
    pub fn bounds_mm3(&self) -> [f64; 2] {
        [self.bounds.lo, self.bounds.hi]
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
}
fn integral(a: &[I], b: &mut Budget<'_>) -> Result<I, Error> {
    let mut sum = p(0.);
    for &v in a {
        b.charge(1)?;
        sum = sum + v;
    }
    Ok(sum / p(a.len() as f64))
}
fn profile(strip: &PolarStrip, algebra: &Algebra, b: &mut Budget<'_>) -> Result<I, Error> {
    let w = &strip.weights;
    b.charge(strip.radius.len() * 2)?;
    let r = strip
        .radius
        .iter()
        .zip(w)
        .map(|(&r, &w)| r * w)
        .collect::<Polynomial>();
    let q = strip
        .phase
        .iter()
        .zip(w)
        .map(|(&q, &w)| q * w)
        .collect::<Polynomial>();
    let dw = poly::derivative(w, b)?;
    let dq = poly::derivative(&q, b)?;
    let qprime = poly::subtract(
        &algebra.multiply(&dq, w, b)?,
        &algebra.multiply(&q, &dw, b)?,
        b,
    )?;
    let r2 = algebra.multiply(&r, &r, b)?;
    let numerator = algebra.multiply(&r2, &qprime, b)?;
    let low = w.iter().map(|v| v.lo).fold(f64::INFINITY, f64::min);
    let high = w.iter().map(|v| v.hi).fold(0., f64::max);
    let center = (low + high) * 0.5;
    if low <= 0. || !center.is_finite() {
        return Err(Error::Unsupported);
    }
    let x = w
        .iter()
        .map(|&w| w / p(center) - p(1.))
        .collect::<Polynomial>();
    let t = p(poly::magnitude(&x));
    if t.hi >= 0.5 {
        return Err(Error::Unsupported);
    }
    // (1+x)^-4 = sum (-1)^k C(k+3,3) x^k. For |x|<=t,
    // tail after K is <= C(K+4,3)t^(K+1)/(1-t)^4: termwise
    // C(K+4+j,3)<=C(K+4,3)C(j+3,3). No quadrature sample claim.
    let mut power = vec![p(1.)];
    let mut total = p(0.);
    for k in 0..=40 {
        let coefficient =
            ((k + 1) * (k + 2) * (k + 3) / 6) as f64 * if k % 2 == 0 { 1. } else { -1. };
        let product = algebra.multiply(&numerator, &power, b)?;
        total = total + p(coefficient) * integral(&product, b)?;
        power = algebra.multiply(&power, &x, b)?;
    }
    let mut tail = p(13244.);
    for _ in 0..41 {
        b.charge(1)?;
        tail = tail * t;
    }
    let one = p(1.) - t;
    tail = tail / (one * one * one * one);
    let error = tail * p(poly::magnitude(&numerator));
    total = (total + I::new(-error.hi, error.hi)) / (p(center) * p(center) * p(center) * p(center));
    if !total.valid() || total.lo <= 0. {
        return Err(Error::Unsupported);
    }
    Ok(total)
}
impl HealedMaterial<'_, '_, '_, '_, '_> {
    /// At any fixed Z, one angular revolution covers exactly one earned phase
    /// cycle. Bounding its actual angular Jacobian uniformly therefore bounds
    /// area by J/(2*lead) times the exact profile integral of r² dq, including
    /// fractional end turns, arbitrary phase origins and either handedness.
    pub fn healed_volume(
        &self,
        options: VolumeOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<HealedVolume, Error> {
        if options.maximum_work == 0
            || options.maximum_work > 64_000_000
            || options.maximum_bytes == 0
            || options.maximum_bytes > 64 * 1024 * 1024
            || !options.maximum_width_mm3.is_finite()
            || options.maximum_width_mm3 <= 0.
        {
            return Err(Error::Invalid);
        }
        let mut b = Budget {
            work: 0,
            remaining: options.maximum_work,
            bytes: self.bytes,
            limit: options.maximum_bytes,
            interrupt,
        };
        b.charge(1)?;
        b.retain(1024 * 1024)?;
        let algebra = Algebra::new(&mut b)?;
        let mut density = I::new(f64::INFINITY, 0.);
        let mut angular_area = p(0.);
        for piece in &self.polar.angular {
            let x = piece.bounds.iter().map(|p| p[0]).collect::<Polynomial>();
            let y = piece.bounds.iter().map(|p| p[1]).collect::<Polynomial>();
            let dx = poly::derivative(&x, &mut b)?;
            let dy = poly::derivative(&y, &mut b)?;
            let jac = poly::subtract(
                &algebra.multiply(&x, &dy, &mut b)?,
                &algebra.multiply(&y, &dx, &mut b)?,
                &mut b,
            )?;
            let width = &piece.knots[1] - &piece.knots[0];
            let width = physical(width.numer(), width.denom(), &mut b)?;
            let normalized = jac.iter().map(|&v| v / width).collect::<Polynomial>();
            let mut pending = vec![(normalized, 0)];
            while let Some((net, depth)) = pending.pop() {
                if depth < 3 {
                    let halves = poly::split(&net, p(0.5), &mut b)?;
                    pending.push((halves[0].clone(), depth + 1));
                    pending.push((halves[1].clone(), depth + 1));
                } else {
                    density.lo = density.lo.min(poly::lower(&net));
                    density.hi = density.hi.max(net.iter().map(|v| v.hi).fold(0., f64::max));
                }
            }
            angular_area = angular_area + integral(&jac, &mut b)?;
        }
        if !density.valid() || density.lo <= 0. {
            return Err(Error::Unsupported);
        }
        let mut profiles = [p(0.); 2];
        for strip in &self.polar.strips {
            let group = usize::from(strip.lower != 0);
            profiles[group] = profiles[group] + profile(strip, &algebra, &mut b)?;
        }
        let unit = integer(self.polar.grid.denominator(), &mut b)?;
        let heights = self
            .polar
            .planes
            .windows(2)
            .map(|z| physical(&(&z[1].1 - &z[0].1), &unit, &mut b))
            .collect::<Result<Vec<_>, _>>()?;
        let leads = [
            physical(&self.polar.leads[0], &unit, &mut b)?,
            physical(&self.polar.leads[1], &unit, &mut b)?,
        ];
        let plain = physical(&self.polar.plain_radius, &unit, &mut b)?;
        let total = density / p(2.)
            * (profiles[0] * heights[0] / leads[0] + profiles[1] * heights[2] / leads[1])
            + angular_area / p(2.) * plain * plain * heights[1];
        if !total.valid()
            || total.lo <= 0.
            || (p(total.hi) - p(total.lo)).hi > options.maximum_width_mm3
        {
            return Err(Error::Unsupported);
        }
        Ok(HealedVolume {
            bounds: total,
            work: b.work,
        })
    }
}
