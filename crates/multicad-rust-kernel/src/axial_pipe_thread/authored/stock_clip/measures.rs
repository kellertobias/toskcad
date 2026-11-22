//! Integrals of the original clipped source, separate from mesh/file measures.
use super::*;
use crate::tapered_thread::{Scalar, Series};
type S = Series<13>;
fn c(v: I) -> S {
    S::interval(v)
}
#[derive(Clone, Copy, Debug)]
pub struct StockClipMeasureOptions {
    pub area_width_mm2: f64,
    pub moment_width_mm3: f64,
    pub centroid_width_mm: f64,
    pub limits: StockClipOptions,
}
impl Default for StockClipMeasureOptions {
    fn default() -> Self {
        Self {
            area_width_mm2: 1e-7,
            moment_width_mm3: 1e-7,
            centroid_width_mm: 1e-7,
            limits: StockClipOptions::default(),
        }
    }
}
/// Actual authored face area and first moment enclosures. No literal mesh,
/// independently decoded geometry or approximate carrier is integrated.
#[derive(Clone, Copy, Debug)]
pub struct StockClipFaceMeasure {
    face: StockClipFace,
    values: [[f64; 2]; 4],
    centroid: [[f64; 2]; 3],
    work: usize,
    bytes: usize,
}
impl StockClipFaceMeasure {
    pub fn face(&self) -> StockClipFace {
        self.face
    }
    pub fn area_mm2(&self) -> [f64; 2] {
        self.values[0]
    }
    pub fn first_moments_mm3(&self) -> [[f64; 2]; 3] {
        [self.values[1], self.values[2], self.values[3]]
    }
    pub fn centroid_mm(&self) -> [[f64; 2]; 3] {
        self.centroid
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_budget_bytes(&self) -> usize {
        self.bytes
    }
}
struct Density {
    corners: [[I; 3]; 4],
    omega: I,
    rq: I,
    rz: I,
    cap: Option<f64>,
}
impl Density {
    fn field<const N: usize>(&self, u: Series<N>, v: Series<N>) -> [Series<N>; 4] {
        let c = Series::<N>::interval;
        let n = Series::<N>::constant;
        let fields: [Series<N>; 3] = std::array::from_fn(|axis| {
            let a = c(self.corners[0][axis]) + c(self.corners[1][axis] - self.corners[0][axis]) * u;
            let b = c(self.corners[2][axis]) + c(self.corners[3][axis] - self.corners[2][axis]) * u;
            a + (b - a) * v
        });
        let [mut r, q, z] = fields;
        // Every admitted quadrature/derivative domain stays within the
        // normalized unit chart. The exact bilinear radius is a convex
        // combination of the SAME four source corners; retain this bound
        // when interval dependence in a+(b-a)*v crosses the axis.
        let radius_min = self
            .corners
            .iter()
            .map(|p| p[0].lo)
            .fold(f64::INFINITY, f64::min);
        let radius_max = self
            .corners
            .iter()
            .map(|p| p[0].hi)
            .fold(f64::NEG_INFINITY, f64::max);
        r.c[0] = I::new(r.c[0].lo.max(radius_min), r.c[0].hi.min(radius_max));
        let qu = c(self.corners[1][1] - self.corners[0][1]) * (n(1.) - v)
            + c(self.corners[3][1] - self.corners[2][1]) * v;
        let [sin, cos] = ((z - q) * c(self.omega)).sin_cos();
        if self.cap.is_some() {
            let area = qu * c(self.omega) * r * r / n(2.);
            let radial = qu * c(self.omega) * r * r * r / n(3.);
            [area, radial * cos, radial * sin, area * z]
        } else {
            let g = if self.rq.lo == 0. && self.rq.hi == 0. {
                // Exact positive-radius identity on natural crest/root and
                // stock charts, avoiding a dependency-inflated sqrt(r²).
                c(self.omega * (I::point(1.) + self.rz.square()).sqrt()) * r
            } else {
                (c(self.omega.square() * (I::point(1.) + (self.rq + self.rz).square())) * r * r
                    + c(self.rq.square()))
                .sqrt()
            };
            let area = qu * c(self.corners[2][2] - self.corners[0][2]) * g;
            [area, area * r * cos, area * r * sin, area * z]
        }
    }
}
#[derive(Clone, Copy)]
struct Cell {
    uv: [[f64; 2]; 2],
    depth: u32,
    budget: [f64; 2],
}
fn power(x: I, n: usize) -> I {
    (0..n).fold(I::point(1.), |a, _| a * x)
}
/// Immutable six-point Gauss intervals. Exact-rational controls prove each
/// disjoint P6 root bracket and its coupled 2/((1-x²)P6'(x)²) weight range.
/// No caller source/precision is used to construct these constants.
fn gauss_nodes() -> [(I, I); 6] {
    let positive = [
        (
            I::new(0.2386191860831969, 0.23861918608319693),
            I::new(0.46791393457269104, 0.4679139345726911),
        ),
        (
            I::new(0.6612093864662645, 0.6612093864662646),
            I::new(0.3607615730481385, 0.36076157304813866),
        ),
        (
            I::new(0.9324695142031519, 0.932469514203152),
            I::new(0.17132449237917027, 0.17132449237917058),
        ),
    ];
    [
        positive[0],
        positive[1],
        positive[2],
        (-positive[0].0, positive[0].1),
        (-positive[1].0, positive[1].1),
        (-positive[2].0, positive[2].1),
    ]
}
fn evaluate(d: &Density, cell: Cell, meter: &mut Meter<'_>) -> Result<([I; 4], [f64; 2]), Error> {
    // Thirty-six interval point fields and two pure-axis order-twelve fields. Fixed-degree
    // series arithmetic/trigonometry is charged before evaluation.
    meter.charge(32768)?;
    let range = cell.uv.map(|v| I::new(v[0], v[1]));
    let width = range.map(|r| I::point(r.hi) - I::point(r.lo));
    let mid = range.map(|r| (I::point(r.lo) + I::point(r.hi)) / I::point(2.));
    let half = width.map(|w| w / I::point(2.));
    let nodes = gauss_nodes();
    let mut result = [I::point(0.); 4];
    for (a, wa) in nodes {
        for (b, wb) in nodes {
            let f = d.field(
                Series::<1>::interval(mid[0] + half[0] * a),
                Series::<1>::interval(mid[1] + half[1] * b),
            );
            for j in 0..4 {
                result[j] = result[j] + f[j].c[0] * half[0] * half[1] * wa * wb;
            }
        }
    }
    let axes = [
        d.field(S::variable(range[0]), c(range[1])),
        d.field(c(range[0]), S::variable(range[1])),
    ];
    let mut steering = [0.; 2];
    for axis in 0..2 {
        // Six-point Gauss: integral minus quadrature is f^(12)*h^13/(12!*11099088).
        // Series coefficient is f^(12)/12!. Sequential tensor integration
        // gives these TWO pure-axis remainders; no mixed derivative omitted.
        let scale = power(width[axis], 13) * width[1 - axis] / I::point(11099088.);
        for j in 0..4 {
            let error = axes[axis][j].c[12] * scale;
            result[j] = result[j] + error;
            steering[axis] += error.magnitude() / cell.budget[usize::from(j != 0)];
        }
    }
    // Finite sealed-domain handoff (specific to these NPT-flat charts):
    // |r|,|r'|,|qu|,|zv| <2e7, |theta'|<900. For flanks H=A*r²+rq²
    // has degree2, |H1,H2|<4e22 and sqrt(H0).lo>1.7. With M=5e22,
    // the outward sqrt recurrence is majorized by h=M*t+M*t²+h².
    // At R=1/(8M)=2.5e-24, h(R)<1, so every coefficient through12
    // is bounded by R^-n. The affine density prefactor has absolute
    // coefficient sum at R <8e21; sin/cos contributes <2. Thus ALL
    // density products/convolution partial sums are <1.6e22*R^-12
    // (<2.7e305), with ample outward-operation slack below1e307.
    // The separate g0 contribution, rq=0 linear-g branch and degree<=4
    // cap prefactors with theta' bounded by900 stay below1e70.
    // The Gauss weights/nodes and widths<=1 cannot overflow these values.
    // This proof includes intermediate sums, not only the final coefficient.
    // A new source family/chart domain must re-earn this handoff.
    Ok((result, steering))
}
impl StockClipPlan {
    pub fn face_measure(
        &self,
        face: StockClipFace,
        options: StockClipMeasureOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<StockClipFaceMeasure, Error> {
        let tolerances = [
            options.area_width_mm2,
            options.moment_width_mm3,
            options.centroid_width_mm,
        ];
        if tolerances.into_iter().any(|v| !v.is_finite() || v <= 0.) {
            return Err(Error::Invalid);
        }
        // A second fixed workspace covers the bounded 64-cell DFS stack,
        // Density and live child point/series caches in addition to Meter scratch.
        let retained = size_of::<Self>() + size_of::<StockClipFaceMeasure>() + SCRATCH;
        let mut work = Meter::new(options.limits, retained, interrupt)?;
        work.charge(self.work)?;
        let mut sum = [I::point(0.); 4];
        for patch in self.patches() {
            let cap = match face {
                StockClipFace::StartCap => Some(0.),
                StockClipFace::EndCap => Some(1.),
                _ => None,
            };
            if cap.is_none() && patch.face() != face {
                continue;
            }
            if cap == Some(0.) && patch.slab != StockClipSlab::Before
                || cap == Some(1.) && patch.slab != StockClipSlab::After
            {
                continue;
            }
            let mut corners = [[I::point(0.); 3]; 4];
            for (i, (u, v)) in [(0., 0.), (1., 0.), (0., 1.), (1., 1.)]
                .into_iter()
                .enumerate()
            {
                let point = self.point_bounds(
                    patch,
                    [u, u],
                    [v, v],
                    work.child_options(options.limits)?,
                    interrupt,
                )?;
                work.charge(point.work_used())?;
                corners[i] = [
                    interval(point.radius_mm())?,
                    interval(point.phase_mm())?,
                    interval(point.axial_mm())?,
                ];
            }
            let k = I::point(self.source.spec.radial_taper);
            let (rq, rz) = match (patch.slab, patch.band) {
                (_, StockClipBand::Lower) => (I::point(1.) / self.slopes[0], k),
                (_, StockClipBand::Upper) => (-I::point(1.) / self.slopes[1], k),
                (StockClipSlab::After, StockClipBand::Crest) => (I::point(0.), I::point(0.)),
                _ => (I::point(0.), k),
            };
            if let Some(v) = cap {
                let start = if v == 0. { 0 } else { 2 };
                corners = [
                    corners[start],
                    corners[start + 1],
                    corners[start],
                    corners[start + 1],
                ];
            }
            let d = Density {
                corners,
                omega: I::point(2.) * pi() / I::point(self.source.spec.pitch_mm),
                rq,
                rz,
                cap,
            };
            let mut pending = [Cell {
                uv: [[0., 1.], [0., 1.]],
                depth: 0,
                budget: [0.; 2],
            }; 64];
            let parts = match face {
                StockClipFace::StartCap | StockClipFace::EndCap => 4.,
                StockClipFace::LowerFlank | StockClipFace::UpperFlank | StockClipFace::Root => 2.,
                _ => 1.,
            };
            // Half the caller width is distributed across the ACTUAL chart
            // count; the remaining half covers outward accumulation.
            pending[0].budget = [
                options.area_width_mm2 / (2. * parts),
                options.moment_width_mm3 / (2. * parts),
            ];
            let mut len = 1;
            while len > 0 {
                len -= 1;
                let cell = pending[len];
                let (value, steering) = evaluate(&d, cell, &mut work)?;

                if value[0].width() <= cell.budget[0]
                    && value[1..].iter().all(|v| v.width() <= cell.budget[1])
                {
                    for j in 0..4 {
                        sum[j] = sum[j] + value[j];
                    }
                    continue;
                }
                if cell.depth >= 24 {
                    return Err(Error::Unresolved);
                }
                let axis = usize::from(steering[1] > steering[0]);
                let [a, b] = cell.uv[axis];
                let mid = a + (b - a) * 0.5;
                // Unit endpoints are exact dyadics with denominator at most
                // 2^24. This midpoint is exactly representable and strict.
                let mut left = cell;
                left.depth += 1;
                left.budget = left.budget.map(|v| v * 0.5);
                let mut right = left;
                left.uv[axis][1] = mid;
                right.uv[axis][0] = mid;
                pending[len] = left;
                pending[len + 1] = right;
                len += 2;
            }
        }
        if sum[0].lo <= 0. {
            return Err(Error::Unresolved);
        }
        let centroid = [sum[1] / sum[0], sum[2] / sum[0], sum[3] / sum[0]];
        if sum[0].width() > options.area_width_mm2
            || sum[1..]
                .iter()
                .any(|v| v.width() > options.moment_width_mm3)
            || centroid
                .iter()
                .any(|v| !v.valid() || v.width() > options.centroid_width_mm)
        {
            return Err(Error::Unresolved);
        }
        work.charge(1)?;
        Ok(StockClipFaceMeasure {
            face,
            values: sum.map(bounds),
            centroid: centroid.map(bounds),
            work: work.used,
            bytes: retained + SCRATCH,
        })
    }
}
#[cfg(test)]
mod tests;
