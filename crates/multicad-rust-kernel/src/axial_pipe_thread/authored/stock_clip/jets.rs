//! Mixed jets of the original authored atlas, not fitted carrier derivatives.
use super::*;
use crate::tapered_thread::interval::J;

/// A sealed bilinear source chart. The borrowed plan owns the original law and
/// symbolic contact; the four enclosed scalar stations are evaluation caches.
#[derive(Debug)]
pub struct StockClipPatchLaw<'a> {
    plan: &'a StockClipPlan,
    patch: StockClipPatch,
    // (u,v) = (0,0),(1,0),(0,1),(1,1), each [radius,unwrapped phase,z].
    corners: [[I; 3]; 4],
    work: usize,
    bytes: usize,
}
#[derive(Clone, Copy, Debug)]
pub struct StockClipJets {
    sample: StockClipSample,
    du: [[f64; 2]; 3],
    dv: [[f64; 2]; 3],
    duu: [[f64; 2]; 3],
    duv: [[f64; 2]; 3],
    dvv: [[f64; 2]; 3],
    work: usize,
    bytes: usize,
}
impl StockClipJets {
    pub fn sample(&self) -> &StockClipSample {
        &self.sample
    }
    pub fn du(&self) -> [[f64; 2]; 3] {
        self.du
    }
    pub fn dv(&self) -> [[f64; 2]; 3] {
        self.dv
    }
    pub fn duu(&self) -> [[f64; 2]; 3] {
        self.duu
    }
    pub fn duv(&self) -> [[f64; 2]; 3] {
        self.duv
    }
    pub fn dvv(&self) -> [[f64; 2]; 3] {
        self.dvv
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_budget_bytes(&self) -> usize {
        self.bytes
    }
}
impl StockClipPlan {
    pub fn patch_law(
        &self,
        patch: StockClipPatch,
        options: StockClipOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<StockClipPatchLaw<'_>, Error> {
        let retained =
            size_of::<Self>() + size_of::<StockClipPatchLaw<'_>>() + size_of::<StockClipJets>();
        let mut meter = Meter::new(options, retained, interrupt)?;
        let mut corners = [[I::point(0.); 3]; 4];
        for (i, (u, v)) in [(0., 0.), (1., 0.), (0., 1.), (1., 1.)]
            .into_iter()
            .enumerate()
        {
            let point = self.point_bounds(
                patch,
                [u, u],
                [v, v],
                StockClipOptions {
                    maximum_work: meter.remaining,
                    ..options
                },
                interrupt,
            )?;
            meter.charge(point.work_used())?;
            corners[i] = [
                interval(point.radius_mm())?,
                interval(point.phase_mm())?,
                interval(point.axial_mm())?,
            ];
        }
        meter.charge(32)?;
        Ok(StockClipPatchLaw {
            plan: self,
            patch,
            corners,
            work: meter.used,
            bytes: retained + SCRATCH,
        })
    }
}
fn constant(v: I) -> J {
    J {
        v,
        ..J::constant(0.)
    }
}
fn bilinear(c: [I; 4], u: J, v: J) -> J {
    // Radius, unwrapped phase and height are exactly bilinear on each of the
    // eight admitted source charts. Enclosed corner coefficients therefore
    // enclose the SAME law and every first/mixed second derivative.
    let a = constant(c[0]) + constant(c[1] - c[0]) * u;
    let b = constant(c[2]) + constant(c[3] - c[2]) * u;
    a + (b - a) * v
}
impl StockClipPatchLaw<'_> {
    pub fn plan(&self) -> &StockClipPlan {
        self.plan
    }
    pub fn patch(&self) -> StockClipPatch {
        self.patch
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_budget_bytes(&self) -> usize {
        self.bytes
    }
    pub fn jets(
        &self,
        u: [f64; 2],
        v: [f64; 2],
        options: StockClipOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<StockClipJets, Error> {
        let retained = self.bytes - SCRATCH;
        let mut meter = Meter::new(options, retained, interrupt)?;
        meter.charge(self.work)?;
        let sample = self.plan.sample_bounds(
            self.patch,
            u,
            v,
            StockClipOptions {
                maximum_work: meter.remaining,
                ..options
            },
            interrupt,
        )?;
        meter.charge(sample.point().work_used())?;
        let u = J::d(sample::unit_parameter(u)?);
        let v = J::z(sample::unit_parameter(v)?);
        meter.charge(2048)?;
        let fields: [J; 3] = std::array::from_fn(|j| bilinear(self.corners.map(|p| p[j]), u, v));
        let [r, q, z] = fields;
        let theta =
            (z - q) * constant(I::point(2.) * pi() / I::point(self.plan.source.spec.pitch_mm));
        let p = [r * theta.cos(), r * theta.sin(), z];
        // Sealed r,z<=1e7, |q|<3p and L/p<128 bound normalized chart
        // derivatives, including mixed terms, far below f64 overflow.
        meter.charge(1)?;
        Ok(StockClipJets {
            sample,
            du: p.map(|j| bounds(j.d)),
            dv: p.map(|j| bounds(j.z)),
            duu: p.map(|j| bounds(j.dd)),
            duv: p.map(|j| bounds(j.dz)),
            dvv: p.map(|j| bounds(j.zz)),
            work: meter.used,
            bytes: self.bytes,
        })
    }
}
