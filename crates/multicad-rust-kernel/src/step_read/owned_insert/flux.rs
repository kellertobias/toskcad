//! Literal surface flux. Global embedding/material admission remains separate.
//! Each actual polynomial span is integrated in its normalized parameter; the
//! original knot Jacobian cancels the span's integration measure.
use super::*;
use crate::tapered_thread::interval::I;
mod volume;
pub use volume::OccupiedVolume;
const WORK: usize = 64_000_000;
const BYTES: usize = 64 * 1024 * 1024;
// Fixed polynomial buffers, a33x33 integer table, scalar-format temporaries
// and call headers. This is accounted query storage, not a process RSS limit.
const SCRATCH: usize = 128 * 1024;
const STRIDE: usize = 33;
const COEFFICIENTS: usize = 3 * STRIDE;

#[derive(Clone, Copy, Debug)]
pub struct FluxOptions {
    pub maximum_width_mm3: f64,
    /// Includes retained support-binding work, then this independent query.
    pub maximum_work: usize,
    /// Retained packet plus fixed query storage and receipt.
    pub maximum_bytes: usize,
}
impl Default for FluxOptions {
    fn default() -> Self {
        Self {
            maximum_width_mm3: 1e-6,
            maximum_work: WORK,
            maximum_bytes: BYTES,
        }
    }
}
impl FluxOptions {
    fn validate(self) -> Result<()> {
        if !self.maximum_width_mm3.is_finite()
            || self.maximum_width_mm3 <= 0.
            || self.maximum_work == 0
            || self.maximum_work > WORK
            || self.maximum_bytes == 0
            || self.maximum_bytes > BYTES
        {
            Err(Error::Invalid)
        } else {
            Ok(())
        }
    }
}
pub struct OrientedFlux<'s, 'p, 'input> {
    bindings: &'s SupportBindings<'p, 'input>,
    oriented_mm3: [f64; 2],
    work: usize,
    bytes: usize,
}
impl<'s, 'p, 'input> OrientedFlux<'s, 'p, 'input> {
    pub fn bindings(&self) -> &'s SupportBindings<'p, 'input> {
        self.bindings
    }
    /// Signed flux of these literal carriers, not yet occupied material mass.
    pub fn oriented_mm3(&self) -> [f64; 2] {
        self.oriented_mm3
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn query_budget_bytes(&self) -> usize {
        self.bytes
    }
}
struct Weights([[u64; 33]; 33]);
impl Weights {
    fn new(w: &mut Ledger<'_>) -> Result<Self> {
        w.step(33 * 33)?;
        let mut c = [[0; 33]; 33];
        for n in 0..33 {
            c[n][0] = 1;
            c[n][n] = 1;
            for k in 1..n {
                c[n][k] = c[n - 1][k - 1] + c[n - 1][k];
            }
        }
        Ok(Self(c))
    }
    fn product(&self, m: usize, n: usize, i: usize, j: usize) -> I {
        // Through degree32, choose(m,i)*choose(n,j)<=choose(32,16)
        // by Vandermonde. Both integer numerator/denominator are exact f64.
        I::point((self.0[m][i] * self.0[n][j]) as f64) / I::point(self.0[m + n][i + j] as f64)
    }
}
struct Poly {
    u: usize,
    v: usize,
    c: [I; COEFFICIENTS],
}
impl Poly {
    fn zero(u: usize, v: usize) -> Self {
        Self {
            u,
            v,
            c: [I::point(0.); COEFFICIENTS],
        }
    }
    fn get(&self, u: usize, v: usize) -> I {
        self.c[u * STRIDE + v]
    }
    fn du(&self, w: &mut Ledger<'_>) -> Result<Self> {
        let mut p = Self::zero(self.u - 1, self.v);
        for i in 0..=p.u {
            for j in 0..=p.v {
                w.step(4)?;
                p.c[i * STRIDE + j] =
                    I::point(self.u as f64) * (self.get(i + 1, j) - self.get(i, j));
            }
        }
        p.checked()
    }
    fn dv(&self, w: &mut Ledger<'_>) -> Result<Self> {
        let mut p = Self::zero(self.u, self.v - 1);
        for i in 0..=p.u {
            for j in 0..=p.v {
                w.step(4)?;
                p.c[i * STRIDE + j] =
                    I::point(self.v as f64) * (self.get(i, j + 1) - self.get(i, j));
            }
        }
        p.checked()
    }
    fn mul(&self, b: &Self, weights: &Weights, w: &mut Ledger<'_>) -> Result<Self> {
        let mut p = Self::zero(self.u + b.u, self.v + b.v);
        for i in 0..=self.u {
            for j in 0..=self.v {
                for k in 0..=b.u {
                    for l in 0..=b.v {
                        w.step(8)?;
                        let factor =
                            weights.product(self.u, b.u, i, k) * weights.product(self.v, b.v, j, l);
                        let slot = (i + k) * STRIDE + j + l;
                        p.c[slot] = p.c[slot] + factor * self.get(i, j) * b.get(k, l);
                    }
                }
            }
        }
        p.checked()
    }
    fn sub(mut self, b: Self, w: &mut Ledger<'_>) -> Result<Self> {
        // All callers pair the same two derivative-product degrees.
        for i in 0..=self.u {
            for j in 0..=self.v {
                w.step(1)?;
                self.c[i * STRIDE + j] = self.get(i, j) - b.get(i, j);
            }
        }
        self.checked()
    }
    fn checked(self) -> Result<Self> {
        if (0..=self.u).any(|i| (0..=self.v).any(|j| !self.get(i, j).valid())) {
            Err(Error::Unsupported)
        } else {
            Ok(self)
        }
    }
    fn integral(&self, w: &mut Ledger<'_>) -> Result<I> {
        let mut sum = I::point(0.);
        for i in 0..=self.u {
            for j in 0..=self.v {
                w.step(1)?;
                sum = sum + self.get(i, j);
            }
        }
        finite(sum / I::point(((self.u + 1) * (self.v + 1)) as f64))
    }
}
fn finite(i: I) -> Result<I> {
    if i.valid() {
        Ok(i)
    } else {
        Err(Error::Unsupported)
    }
}
fn point(p: &Parsed<'_>, id: u32, w: &mut Ledger<'_>) -> Result<[I; 3]> {
    w.step(256)?;
    let s = p.point(id).ok_or(Error::Invalid)?.coordinates();
    let [x, y, z] = s else {
        return Err(Error::Invalid);
    };
    let mut output = [I::point(0.); 3];
    for (out, scalar) in output.iter_mut().zip([x, y, z]) {
        w.step(scalar.spelling().len())?;
        let [a, b] = scalar.bounds();
        *out = finite(I::new(a, b))?;
    }
    Ok(output)
}
fn wall(p: &Parsed<'_>, t: &TensorSurface<'_>, weights: &Weights, w: &mut Ledger<'_>) -> Result<I> {
    let rows = t.rows();
    let mut total = I::point(0.);
    for span in 0..(rows[0].len() - 1) / 11 {
        let mut xyz = std::array::from_fn::<_, 3, _>(|_| Poly::zero(1, 11));
        for (u, row) in rows.into_iter().enumerate() {
            for v in 0..12 {
                let values = point(p, row[11 * span + v], w)?;
                for (poly, value) in xyz.iter_mut().zip(values) {
                    poly.c[u * STRIDE + v] = value;
                }
            }
        }
        let xu = xyz[0].du(w)?;
        let yv = xyz[1].dv(w)?;
        let xv = xyz[0].dv(w)?;
        let yu = xyz[1].du(w)?;
        let determinant = xu.mul(&yv, weights, w)?.sub(xv.mul(&yu, weights, w)?, w)?;
        total = finite(total + xyz[2].mul(&determinant, weights, w)?.integral(w)?)?;
    }
    Ok(total)
}
fn cap(
    p: &Parsed<'_>,
    face: &Face,
    height: Scalar<'_>,
    weights: &Weights,
    w: &mut Ledger<'_>,
) -> Result<I> {
    w.step(height.spelling().len() + 32)?;
    let [lo, hi] = height.bounds();
    let height = finite(I::new(lo, hi))?;
    let mut area = I::point(0.);
    // FACE_BOUND and EDGE_CURVE orientations are exactly true in this sealed
    // schema. Encoded coedge traversal already carries the face orientation;
    // multiplying same_sense again would incorrectly reverse a bottom cap.
    for boundary in face.loops() {
        for coedge in boundary.coedges() {
            let edge = p.edge(coedge.edge).expect("sealed source edge");
            let c = p.curve(edge.curve_entity()).expect("sealed source curve");
            for span in 0..(c.pole_entities().len() - 1) / 11 {
                let mut xy = std::array::from_fn::<_, 2, _>(|_| Poly::zero(0, 11));
                for v in 0..12 {
                    let values = point(p, c.pole_entities()[11 * span + v], w)?;
                    for (poly, value) in xy.iter_mut().zip(values) {
                        poly.c[v] = value;
                    }
                }
                let dy = xy[1].dv(w)?;
                let dx = xy[0].dv(w)?;
                let integral = xy[0]
                    .mul(&dy, weights, w)?
                    .sub(xy[1].mul(&dx, weights, w)?, w)?
                    .integral(w)?
                    / I::point(2.);
                area = finite(area + if coedge.forward { integral } else { -integral })?;
            }
        }
    }
    finite(height * area)
}
impl<'p, 'input> SupportBindings<'p, 'input> {
    pub fn oriented_flux<'s>(
        &'s self,
        options: FluxOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<OrientedFlux<'s, 'p, 'input>> {
        options.validate()?;
        let mut w = Ledger {
            used: 0,
            remaining: options.maximum_work,
            interrupt,
        };
        w.step(self.work_used())?;
        let bytes = self.retained_budget_bytes() + SCRATCH + size_of::<OrientedFlux<'_, '_, '_>>();
        if bytes > options.maximum_bytes {
            return Err(Error::Limit);
        }
        let weights = Weights::new(&mut w)?;
        let p = self.source();
        let mut value = I::point(0.);
        for binding in self.tensors() {
            let flux = wall(p, &p.tensors()[binding.tensor_index()], &weights, &mut w)?;
            let sign = if p.faces()[binding.face_index()].same_sense() {
                1.
            } else {
                -1.
            };
            value = finite(value + I::point(sign) * flux)?;
        }
        for (face, height) in self.cap_face_indices().into_iter().zip(self.cap_heights()) {
            value = finite(value + cap(p, &p.faces()[face], height, &weights, &mut w)?)?;
        }
        w.step(1)?;
        if (I::point(value.hi) - I::point(value.lo)).hi > options.maximum_width_mm3 {
            return Err(Error::Unsupported);
        }
        Ok(OrientedFlux {
            bindings: self,
            oriented_mm3: [value.lo, value.hi],
            work: w.used,
            bytes,
        })
    }
}

#[cfg(test)]
mod tests;
