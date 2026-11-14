use super::*;
use crate::authored_internal_treatment::BoundaryCurve;
type V = [I; 3];
pub(super) type Basis = [[I; 11]; 11];
fn p(x: f64) -> I {
    I::point(x)
}
fn choose(n: usize, k: usize) -> f64 {
    let mut v = 1.;
    for j in 0..k {
        v = v * (n - j) as f64 / (j + 1) as f64;
    }
    v
}
fn power(x: I, n: usize) -> I {
    let mut v = p(1.);
    for _ in 0..n {
        v = v * x;
    }
    v
}
fn l1(v: V) -> I {
    v.into_iter().fold(p(0.), |a, x| a + p(x.magnitude()))
}
fn finite_upper(v: I) -> Result<f64, Error> {
    if v.valid() && v.hi.is_finite() {
        Ok(v.hi)
    } else {
        Err(Error::Unresolved)
    }
}
// The eventual decimal token rounding to x lies between these adjacent floats.
// This transfer is deliberately not an exact raw-file incidence certificate.
fn literal(x: f64) -> I {
    I::new(x.next_down(), x.next_up())
}
// Bernstein-10 controls of (s-1/2)^k, using exact small binomial integers.
pub(super) fn derivative_basis(w: &mut Work<'_>) -> Result<Basis, Error> {
    let mut out = [[p(0.); 11]; 11];
    for (i, row) in out.iter_mut().enumerate() {
        for (k, v) in row.iter_mut().enumerate() {
            for j in 0..=k.min(i) {
                w.charge(1)?;
                *v = *v
                    + p(choose(k, j)) * power(p(-0.5), k - j) * p(choose(i, j)) / p(choose(10, j));
            }
        }
    }
    Ok(out)
}
#[derive(Clone, Copy)]
pub(super) struct Law {
    base: V,
    delta: V,
    omega: I,
    tau: I,
}
impl Law {
    fn new(source: &BoundaryCurve, left: bool) -> Self {
        let [base, delta] = source
            .affine_rzv_bounds()
            .map(|v| v.map(|x| I::new(x.lower, x.upper)));
        Self::from_bounds([base, delta], left)
    }
    pub(super) fn from_bounds([base, delta]: [V; 2], left: bool) -> Self {
        let tau = I::new(
            std::f64::consts::PI.next_down(),
            std::f64::consts::PI.next_up(),
        ) * p(if left { -2. } else { 2. });
        Self {
            base,
            delta,
            omega: tau * delta[2],
            tau,
        }
    }
    pub(super) fn point(&self, t: f64, w: &mut Work<'_>) -> Result<V, Error> {
        Ok(self.jets(t, w)?[0])
    }
    fn jets(&self, t: f64, w: &mut Work<'_>) -> Result<[V; 12], Error> {
        w.charge(48)?;
        let r = self.base[0] + self.delta[0] * p(t);
        let theta = self.tau * (self.base[2] + self.delta[2] * p(t));
        let e = [theta.cos(), theta.sin()];
        let mut out = [[p(0.); 3]; 12];
        out[0] = [r * e[0], r * e[1], self.base[1] + self.delta[1] * p(t)];
        out[1][2] = self.delta[1];
        let mut pe = e;
        for n in 1..12 {
            w.charge(8)?;
            out[n][0] = -self.omega * out[n - 1][1] + self.delta[0] * pe[0];
            out[n][1] = self.omega * out[n - 1][0] + self.delta[0] * pe[1];
            pe = [-self.omega * pe[1], self.omega * pe[0]];
        }
        Ok(out)
    }
    // D^12[(r0+dr*t) exp(i*w*t)] = (iw)^12*r*E + 12*dr*(iw)^11*E.
    // Twice the component bound encloses the vector L1 norm on the whole span.
    fn remainder(&self, a: f64, b: f64) -> I {
        let r = (self.base[0] + self.delta[0] * I::new(a, b)).magnitude();
        let omega = p(self.omega.magnitude());
        p(2.) * (power(omega, 12) * p(r) + p(12.) * power(omega, 11) * p(self.delta[0].magnitude()))
    }
}
fn ideals(jets: &[V; 12], h: f64) -> [V; 6] {
    let mut out = [[p(0.); 3]; 6];
    for (i, v) in out.iter_mut().enumerate() {
        let mut factor = p(1.);
        for (j, jet) in jets.iter().enumerate().take(i + 1) {
            if j > 0 {
                factor = factor * p(h) / p((12 - j) as f64);
            }
            for k in 0..3 {
                v[k] = v[k] + jet[k] * factor * p(choose(i, j));
            }
        }
    }
    out
}
fn certify(
    law: &Law,
    a: f64,
    b: f64,
    ideal: &[V; 12],
    poles: &[[f64; 3]],
    basis: &Basis,
    w: &mut Work<'_>,
) -> Result<(f64, f64), Error> {
    let h = b - a;
    let mut error = p(0.);
    for i in 0..12 {
        w.charge(4)?;
        let e = l1(std::array::from_fn(|k| literal(poles[i][k]) - ideal[i][k]));
        error = I::new(0., error.hi.max(finite_upper(e)?));
    }
    let m = law.remainder(a, b);
    // Six derivatives at each endpoint give the Hermite remainder
    // M12*(t-a)^6*(b-t)^6/12! <= M12*h^12/(4096*12!).
    let pos = finite_upper(error + m * power(p(h), 12) / p(4096. * 479001600.))?;
    // Independently compare the LITERAL polynomial derivative with Taylor-10
    // of C' at the midpoint. Its remainder is M12*(h/2)^11/11!: this does
    // not differentiate a Hermite remainder containing an unknown xi(t).
    let jet = law.jets(a + h * 0.5, w)?;
    let mut derivative = p(0.);
    for (i, row) in basis.iter().enumerate() {
        let mut t = [p(0.); 3];
        let mut factor = p(1.);
        for k in 0..11 {
            w.charge(4)?;
            if k > 0 {
                factor = factor * p(h) / p(k as f64);
            }
            for axis in 0..3 {
                t[axis] = t[axis] + jet[k + 1][axis] * factor * row[k];
            }
        }
        let d = std::array::from_fn(|axis| {
            p(11.) * (literal(poles[i + 1][axis]) - literal(poles[i][axis])) / p(h) - t[axis]
        });
        let e = l1(d);
        derivative = I::new(0., derivative.hi.max(finite_upper(e)?));
    }
    let der = finite_upper(derivative + m * power(p(h * 0.5), 11) / p(39916800.))?;
    Ok((pos, der))
}
#[allow(clippy::too_many_arguments)]
pub(super) fn prepare(
    key: BoundaryCurveKey,
    endpoints: [BoundaryVertexKey; 2],
    source: &BoundaryCurve,
    left: bool,
    ends: [[f64; 3]; 2],
    cap: Option<f64>,
    options: CurveOptions,
    basis: &Basis,
    memory: &mut Memory,
    w: &mut Work<'_>,
) -> Result<CurveCarrier, Error> {
    let law = Law::new(source, left);
    if law.delta[2].lo == 0. && law.delta[2].hi == 0. {
        let j = law.jets(0., w)?;
        let end = law.jets(1., w)?;
        let pos = l1(std::array::from_fn(|k| literal(ends[0][k]) - j[0][k]))
            .hi
            .max(l1(std::array::from_fn(|k| literal(ends[1][k]) - end[0][k])).hi);
        let der = l1(std::array::from_fn(|k| {
            literal(ends[1][k]) - literal(ends[0][k]) - j[1][k]
        }))
        .hi;
        if pos > options.position_error_mm || der > options.derivative_error_mm {
            return Err(Error::Unresolved);
        }
        let mut poles = memory.allocate(2, w)?;
        poles.extend(ends);
        let mut knots = memory.allocate(2, w)?;
        knots.extend([0., 1.]);
        return Ok(CurveCarrier {
            key,
            endpoints,
            degree: 1,
            knots,
            poles,
            position_error_mm: pos,
            derivative_error_mm: der,
        });
    }
    let mut n = 1;
    loop {
        let fit = fixed(&law, n, ends, cap, basis, memory, w)?;
        if fit.position_error_mm <= options.position_error_mm
            && fit.derivative_error_mm <= options.derivative_error_mm
        {
            return Ok(CurveCarrier {
                key,
                endpoints,
                degree: 11,
                knots: fit.knots,
                poles: fit.poles,
                position_error_mm: fit.position_error_mm,
                derivative_error_mm: fit.derivative_error_mm,
            });
        }
        memory.release(&fit.poles);
        memory.release(&fit.knots);
        if n > options.maximum_spans_per_curve / 2 {
            return Err(Error::Unresolved);
        }
        n *= 2;
    }
}

pub(super) struct Fixed {
    pub knots: Vec<f64>,
    pub poles: Vec<[f64; 3]>,
    pub position_error_mm: f64,
    pub derivative_error_mm: f64,
}
/// Fixed dyadic subdivision for coupled rows. Caller supplies n=power-of-two
/// <=1024; accepting a row never accepts a surface or occupied material.
pub(super) fn fixed(
    law: &Law,
    n: usize,
    ends: [[f64; 3]; 2],
    cap: Option<f64>,
    basis: &Basis,
    memory: &mut Memory,
    w: &mut Work<'_>,
) -> Result<Fixed, Error> {
    let mut poles = memory.allocate(11 * n + 1, w)?;
    let mut knots = memory.allocate(n + 1, w)?;
    let mut pos: f64 = 0.;
    let mut der: f64 = 0.;
    for span in 0..n {
        let a = span as f64 / n as f64;
        let b = (span + 1) as f64 / n as f64;
        let h = b - a;
        knots.push(a);
        let ja = law.jets(a, w)?;
        let jb = law.jets(b, w)?;
        let start = ideals(&ja, h);
        let end = ideals(&jb, -h);
        let mut ideal = [[p(0.); 3]; 12];
        for i in 0..6 {
            ideal[i] = start[i];
            ideal[11 - i] = end[i];
        }
        let mut local = ideal.map(|v| v.map(|x| x.midpoint()));
        local[0] = if span == 0 { ends[0] } else { poles[11 * span] };
        if span + 1 == n {
            local[11] = ends[1];
        }
        if let Some(z) = cap {
            for p in &mut local {
                p[2] = z;
            }
        }
        let (ep, ed) = certify(law, a, b, &ideal, &local, basis, w)?;
        pos = pos.max(ep);
        der = der.max(ed);
        if span == 0 {
            poles.push(local[0]);
        }
        poles.extend_from_slice(&local[1..]);
    }
    knots.push(1.);
    Ok(Fixed {
        knots,
        poles,
        position_error_mm: pos,
        derivative_error_mm: der,
    })
}

pub(super) fn recheck(
    source: &BoundaryCurve,
    left: bool,
    curve: &CurveCarrier,
    poles: &[[f64; 3]],
    basis: &Basis,
    w: &mut Work<'_>,
) -> Result<(f64, f64), Error> {
    let law = Law::new(source, left);
    if curve.degree == 1 {
        let ja = law.jets(0., w)?;
        let jb = law.jets(1., w)?;
        let pos = l1(std::array::from_fn(|k| literal(poles[0][k]) - ja[0][k]))
            .hi
            .max(l1(std::array::from_fn(|k| literal(poles[1][k]) - jb[0][k])).hi);
        let der = l1(std::array::from_fn(|k| {
            literal(poles[1][k]) - literal(poles[0][k]) - ja[1][k]
        }))
        .hi;
        return Ok((pos, der));
    }
    let mut pos: f64 = 0.;
    let mut der: f64 = 0.;
    for (span, pair) in curve.knots.windows(2).enumerate() {
        let [a, b] = [pair[0], pair[1]];
        let start = ideals(&law.jets(a, w)?, b - a);
        let end = ideals(&law.jets(b, w)?, a - b);
        let mut ideal = [[p(0.); 3]; 12];
        for i in 0..6 {
            ideal[i] = start[i];
            ideal[11 - i] = end[i];
        }
        let (ep, ed) = certify(
            &law,
            a,
            b,
            &ideal,
            &poles[span * 11..span * 11 + 12],
            basis,
            w,
        )?;
        pos = pos.max(ep);
        der = der.max(ed);
    }
    Ok((pos, der))
}
