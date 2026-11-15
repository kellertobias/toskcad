//! Paired active-source charts. Input atlases independently own their finite
//! topology proofs. Taylor subtraction compares complete trim expressions;
//! finite cutter domains are never assumed equal or replaced by a recipe.
use super::*;
use crate::tapered_thread::{Scalar, Series};
type S = Series<9>;
const CELLS: usize = 262_144;
impl evaluate::Value for S {
    fn bound(v: I) -> Self {
        Self::interval(v)
    }
    fn parameter(v: I) -> Self {
        Self::variable(v)
    }
    fn value(self) -> I {
        self.c[0]
    }
    fn checked(self) -> Result<Self, Error> {
        checked(self)
    }
    fn sine(self) -> Self {
        self.sin_cos()[0]
    }
}
#[derive(Clone, Copy, Debug)]
pub struct PatchCorrespondence {
    pub position_mm: f64,
    pub du_mm: f64,
    pub dv_mm: f64,
    /// Actual paired leaf attempts, including rejected parent cells.
    pub cells_used: usize,
}
#[derive(Clone, Copy)]
pub struct PatchCorrespondenceOptions {
    pub maximum_position_mm: f64,
    pub maximum_derivative_mm: f64,
    pub maximum_cells: usize,
}
impl Default for PatchCorrespondenceOptions {
    fn default() -> Self {
        Self {
            maximum_position_mm: 1e-7,
            maximum_derivative_mm: 1e-7,
            maximum_cells: CELLS,
        }
    }
}
fn c(x: I) -> S {
    S::interval(x)
}
fn n(x: f64) -> S {
    S::constant(x)
}
fn checked(s: S) -> Result<S, Error> {
    s.checked().map_err(|_| Error::Unresolved)
}
fn atan(x: S) -> Result<S, Error> {
    let dx = S {
        c: std::array::from_fn(|i| {
            if i < 8 {
                x.c[i + 1] * I::point((i + 1) as f64)
            } else {
                I::point(0.)
            }
        }),
    };
    let d = dx / (n(1.) + x * x);
    let mut out = c(x.c[0].atan());
    for i in 1..9 {
        out.c[i] = d.c[i - 1] / I::point(i as f64);
    }
    checked(out)
}
// Canonical scalar functions before exact quarter rotations/reflections.
// Envelope/Floor: radius, phi, z; Gap:radius,theta,0;
// Wall: x,y,z; Trace: top x,top y,0.
fn reduced(atlas: &PatchAtlas<'_>, index: usize, u: I) -> Result<[S; 3], Error> {
    let events = atlas.topology.events();
    let source = events.source();
    let k = &source.corners[0];
    let e = &source.envelope;
    let result = match atlas.patches[index].law {
        Law::Envelope(i, _) => {
            let [x, y, z] = evaluate::trim_state::<S>(events, i, u)?;
            [checked((x * x + y * y).sqrt())?, atan(y / x)?, z]
        }
        Law::Wall(i, ..) => evaluate::trim_state::<S>(events, i, u)?,
        Law::Trace(i, ..) => {
            let [x, y, z] = evaluate::trim_state::<S>(events, i, u)?;
            [x + c(events.beta_tangent) * z, y, n(0.)]
        }
        Law::Gap(_) => {
            let p = events.event(EventRole::MouthEnvelope).point;
            let phi = (p[1] / p[0]).atan();
            [
                c(e.entry),
                c(phi) + S::variable(u) * c(pi() / I::point(2.) - I::point(2.) * phi),
                n(0.),
            ]
        }
        Law::Floor(_) => [
            c(k.b_half) * (n(1.) - S::variable(u)),
            n(0.),
            c(k.b_depth) + S::variable(u) * c(e.tip - k.b_depth),
        ],
    };
    Ok([
        checked(result[0])?,
        checked(result[1])?,
        checked(result[2])?,
    ])
}
fn abs(x: I) -> I {
    I::point(x.magnitude())
}
fn sum(x: [I; 3]) -> I {
    x.into_iter().fold(I::point(0.), |a, b| a + b)
}
fn maximum(x: [I; 3]) -> f64 {
    x.into_iter().map(|x| x.hi).fold(0., f64::max)
}
fn difference(a: S, b: S, ar: S, br: S, h: f64) -> [I; 2] {
    let t = I::new(-h, h);
    let mut value = a.c[7] - b.c[7];
    for i in (0..7).rev() {
        value = value * t + (a.c[i] - b.c[i]);
    }
    let mut derivative = I::point(7.) * (a.c[7] - b.c[7]);
    for i in (1..7).rev() {
        derivative = derivative * t + I::point(i as f64) * (a.c[i] - b.c[i]);
    }
    let remainder = abs(ar.c[8]) + abs(br.c[8]);
    let mut power = I::point(1.);
    for _ in 0..7 {
        power = power * I::point(h);
    }
    [
        abs(value) + remainder * power * I::point(h),
        abs(derivative) + I::point(8.) * remainder * power,
    ]
}
fn lift(
    a: &PatchAtlas<'_>,
    b: &PatchAtlas<'_>,
    index: usize,
    whole: [[S; 3]; 2],
    delta: [[I; 2]; 3],
) -> [I; 3] {
    let [x, y] = whole;
    let d = delta.map(|v| v[0]);
    let du = delta.map(|v| v[1]);
    match a.patches[index].law {
        Law::Wall(..) => {
            let beta = a.topology.events().beta_tangent;
            let db = abs(beta - b.topology.events().beta_tangent);
            let z = abs(y[2].c[0]);
            let dz = abs(y[2].c[1]);
            [
                sum(d) + db * z + abs(beta) * d[2],
                sum(du) + db * dz + abs(beta) * du[2],
                d[2] + db * z + abs(beta) * d[2],
            ]
        }
        Law::Trace(..) => [sum(d), sum(du), sum(d)],
        Law::Envelope(..) | Law::Floor(_) | Law::Gap(_) => {
            // Unit radial/tangent vectors are 1-Lipschitz in angle. Expand
            // the product differences using common source parameters before
            // taking independent interval magnitudes.
            let r = abs(y[0].c[0]);
            let rp = abs(y[0].c[1]);
            let angle = d[1];
            let ap = abs(y[1].c[1]);
            let asp = abs(x[1].c[1]);
            let p = d[0] + r * angle + d[2];
            let u = du[0] + rp * angle + d[0] * asp + r * (du[1] + ap * angle) + du[2];
            let v = if matches!(a.patches[index].law, Law::Gap(_)) {
                d[0] + r * angle
            } else {
                let width = abs(pi() / I::point(2.) - I::point(2.) * x[1].c[0]);
                let other_width = abs(pi() / I::point(2.) - I::point(2.) * y[1].c[0]);
                d[0] * width + r * (I::point(2.) * angle + other_width * angle)
            };
            [p, u, v]
        }
    }
}
impl PatchAtlas<'_> {
    /// Complete paired active charts at identical normalized coordinates.
    /// No retained allocation: depth-first traversal owns a fixed 64-cell
    /// stack; receipt cells count both successful leaves and rejected parents.
    pub fn source_correspondence(
        &self,
        other: &PatchAtlas<'_>,
        options: PatchCorrespondenceOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<PatchCorrespondence, Error> {
        let poll = || {
            if interrupt.is_some_and(|f| f()) {
                Err(Error::Cancelled)
            } else {
                Ok(())
            }
        };
        poll()?;
        if [options.maximum_position_mm, options.maximum_derivative_mm]
            .into_iter()
            .any(|v| !v.is_finite() || v <= 0. || v > 1e-7)
            || options.maximum_cells > CELLS
        {
            return Err(Error::Invalid);
        }
        // Sole private builder consumes BoundaryEvents::trims:[FiniteTrim;4]
        // and emits four quarters of five charts per trim, plus eight
        // gap/floor charts: 4*4*5+8=88. Neither atlas exposes mutation.
        assert_eq!(self.patches.len(), other.patches.len());
        let mut result = PatchCorrespondence {
            position_mm: 0.,
            du_mm: 0.,
            dv_mm: 0.,
            cells_used: 0,
        };
        for index in 0..self.patches.len() {
            if self.patches[index].key != other.patches[index].key {
                return Err(Error::Unresolved);
            }
            let mut stack = [(0., 1., 0usize); 64];
            let mut count = 1;
            while count > 0 {
                poll()?;
                if result.cells_used == options.maximum_cells {
                    return Err(Error::Limit);
                }
                result.cells_used += 1;
                count -= 1;
                let (lo, hi, depth) = stack[count];
                let mid = (lo + hi) * 0.5;
                let h = (hi - lo) * 0.5;
                let left = reduced(self, index, I::point(mid));
                let right = reduced(other, index, I::point(mid));
                let aw = reduced(self, index, I::new(lo, hi));
                let bw = reduced(other, index, I::new(lo, hi));
                let bound = match (left, right, aw, bw) {
                    (Ok(a), Ok(b), Ok(ar), Ok(br)) => Some(lift(
                        self,
                        other,
                        index,
                        [ar, br],
                        std::array::from_fn(|i| difference(a[i], b[i], ar[i], br[i], h)),
                    )),
                    _ => None,
                };
                if let Some(bound) = bound.filter(|b| {
                    b.iter().all(|v| v.valid() && v.hi.is_finite())
                        && b[0].hi <= options.maximum_position_mm
                        && maximum([b[1], b[2], I::point(0.)]) <= options.maximum_derivative_mm
                }) {
                    result.position_mm = result.position_mm.max(bound[0].hi);
                    result.du_mm = result.du_mm.max(bound[1].hi);
                    result.dv_mm = result.dv_mm.max(bound[2].hi);
                } else {
                    if depth == 60 || mid == lo || mid == hi {
                        return Err(Error::Unresolved);
                    }
                    // Depth <=60 implies <=61 pending cells in depth-first order.
                    stack[count] = (mid, hi, depth + 1);
                    stack[count + 1] = (lo, mid, depth + 1);
                    count += 2;
                }
            }
        }
        poll()?;
        Ok(result)
    }
}

#[cfg(test)]
mod tests;
