//! Sound chart-derivative magnitude bounds over the decoded supports, for
//! transporting local derivative proofs through a decoded placement operator
//! (docs/architecture/hanger-placed-step-transport.md). Convex-hull bounds of
//! derivative control nets; rational supports use (|N'| + |S||W'|) / w_min
//! with positive weights. These are upper bounds, not correspondence errors.
use super::*;
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DerivativeHull {
    /// Upper bound of |dS/du| over every tensor support (rational U span).
    pub tensor_u_mm: f64,
    /// Upper bound of |dS/dv| over every tensor support (degree-10 pieces).
    pub tensor_v_mm: f64,
    /// Upper bound of |dC/dt| over every 3D edge curve.
    pub curve_mm: f64,
    /// Upper bound of either planar chart axis norm.
    pub plane_mm: f64,
    pub work_used: usize,
}
struct Meter<'a> {
    used: usize,
    limit: usize,
    interrupt: Option<&'a dyn Fn() -> bool>,
}
impl Meter<'_> {
    fn step(&mut self, n: usize) -> Result<(), Error> {
        if self.interrupt.is_some_and(|f| f()) {
            return Err(Error::Cancelled);
        }
        self.used = self.used.checked_add(n).ok_or(Error::Limit)?;
        if self.used > self.limit {
            return Err(Error::Limit);
        }
        Ok(())
    }
}
fn norm(v: [I; 3]) -> Result<f64, Error> {
    let s = v.into_iter().fold(I::point(0.), |s, x| s + x.square());
    if !s.valid() {
        return Err(Error::Invalid);
    }
    let r = I::new(s.lo.max(0.), s.hi).sqrt();
    if !r.valid() || !r.hi.is_finite() {
        return Err(Error::Invalid);
    }
    Ok(r.hi)
}
fn point(d: &Decoded<'_>, id: u32) -> Result<[I; 3], Error> {
    let p = d.points().get(&id).ok_or(Error::Invalid)?;
    let c = p.coordinates();
    if c.len() != 3 {
        return Err(Error::Invalid);
    }
    Ok([c[0].interval(), c[1].interval(), c[2].interval()])
}
/// Homogeneous net (w P, w) with positive weights; empty weights mean 1.
fn net(d: &Decoded<'_>, ids: &[u32], weights: &[Scalar<'_>]) -> Result<Vec<([I; 3], I)>, Error> {
    if !weights.is_empty() && weights.len() != ids.len() {
        return Err(Error::Invalid);
    }
    ids.iter()
        .enumerate()
        .map(|(k, &id)| {
            let w = weights.get(k).map_or(I::point(1.), Scalar::interval);
            if !w.valid() || w.lo <= 0. {
                return Err(Error::Invalid);
            }
            let p = point(d, id)?;
            Ok((p.map(|x| x * w), w))
        })
        .collect()
}
/// sup |d(N/W)| from hull maxima of the derivative nets of N and W.
fn rational(dn: f64, dw: f64, p_max: f64, w_min: f64) -> Result<f64, Error> {
    let r = (I::point(dn) + I::point(p_max) * I::point(dw)) / I::point(w_min);
    if !r.valid() || !r.hi.is_finite() {
        return Err(Error::Invalid);
    }
    Ok(r.hi)
}
/// Maximum derivative-net norms for consecutive pairs (a, b) scaled by gain.
fn pair(a: &([I; 3], I), b: &([I; 3], I), gain: I) -> Result<(f64, f64), Error> {
    let dn = norm(std::array::from_fn(|k| (b.0[k] - a.0[k]) * gain))?;
    let dw = ((b.1 - a.1) * gain).magnitude();
    if !dw.is_finite() {
        return Err(Error::Invalid);
    }
    Ok((dn, dw))
}
fn extrema(d: &Decoded<'_>, ids: &[u32], net: &[([I; 3], I)]) -> Result<(f64, f64), Error> {
    let mut p_max = 0f64;
    let mut w_min = f64::INFINITY;
    for (&id, (_, w)) in ids.iter().zip(net) {
        p_max = p_max.max(norm(point(d, id)?)?);
        w_min = w_min.min(w.lo);
    }
    if w_min.is_nan() || w_min <= 0. {
        return Err(Error::Invalid);
    }
    Ok((p_max, w_min))
}
impl Decoded<'_> {
    pub fn derivative_hull(
        &self,
        maximum_work: usize,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<DerivativeHull, Error> {
        let mut m = Meter {
            used: 0,
            limit: maximum_work,
            interrupt,
        };
        let mut hull = DerivativeHull {
            tensor_u_mm: 0.,
            tensor_v_mm: 0.,
            curve_mm: 0.,
            plane_mm: 0.,
            work_used: 0,
        };
        for t in self.tensors() {
            let [du, dv] = t.degrees();
            let [nu, nv] = t.dimensions();
            let ids = t.point_entities();
            m.step(ids.len().checked_mul(4).ok_or(Error::Limit)?)?;
            let h = net(self, ids, t.weights())?;
            let (p_max, w_min) = extrema(self, ids, &h)?;
            let knots = t.knots_v();
            if nu != du + 1 || nv != dv * (knots.len() - 1) + 1 {
                return Err(Error::Invalid);
            }
            // U: one Bezier span on [0, 1], derivative gain du.
            let (mut n_u, mut w_u) = (0f64, 0f64);
            for i in 0..du {
                for j in 0..nv {
                    let (a, b) = pair(&h[i * nv + j], &h[(i + 1) * nv + j], I::point(du as f64))?;
                    n_u = n_u.max(a);
                    w_u = w_u.max(b);
                }
            }
            // V: C0-joined Bezier pieces of degree dv between distinct knots.
            let (mut n_v, mut w_v) = (0f64, 0f64);
            for s in 0..knots.len() - 1 {
                let span = knots[s + 1].interval() - knots[s].interval();
                if !span.valid() || span.lo <= 0. {
                    return Err(Error::Invalid);
                }
                let gain = I::point(dv as f64) / span;
                for i in 0..nu {
                    for k in 0..dv {
                        let j = s * dv + k;
                        let (a, b) = pair(&h[i * nv + j], &h[i * nv + j + 1], gain)?;
                        n_v = n_v.max(a);
                        w_v = w_v.max(b);
                    }
                }
            }
            hull.tensor_u_mm = hull.tensor_u_mm.max(rational(n_u, w_u, p_max, w_min)?);
            hull.tensor_v_mm = hull.tensor_v_mm.max(rational(n_v, w_v, p_max, w_min)?);
        }
        for c in self.curves().values() {
            let ids = c.point_entities();
            // Only 3D edge curves are mapped; 2D pcurves are operator-free.
            if self
                .points()
                .get(ids.first().ok_or(Error::Invalid)?)
                .is_none_or(|p| p.coordinates().len() != 3)
            {
                continue;
            }
            m.step(ids.len().checked_mul(4).ok_or(Error::Limit)?)?;
            let degree = c.degree();
            let knots = c.knots();
            if ids.len() != degree * (knots.len() - 1) + 1 {
                return Err(Error::Invalid);
            }
            let h = net(self, ids, c.weights())?;
            let (p_max, w_min) = extrema(self, ids, &h)?;
            let (mut n, mut w) = (0f64, 0f64);
            for s in 0..knots.len() - 1 {
                let span = knots[s + 1].interval() - knots[s].interval();
                if !span.valid() || span.lo <= 0. {
                    return Err(Error::Invalid);
                }
                let gain = I::point(degree as f64) / span;
                for k in 0..degree {
                    let (a, b) = pair(&h[s * degree + k], &h[s * degree + k + 1], gain)?;
                    n = n.max(a);
                    w = w.max(b);
                }
            }
            hull.curve_mm = hull.curve_mm.max(rational(n, w, p_max, w_min)?);
        }
        // ISO 10303-42 normalizes AXIS2_PLACEMENT_3D directions, so a PLANE's
        // chart axes are unit vectors whatever ratios are written.
        if !self.planes().is_empty() {
            m.step(self.planes().len())?;
            hull.plane_mm = 1.;
        }
        hull.work_used = m.used;
        Ok(hull)
    }
}
