//! Exact source-edge restrictions and outward Taylor coefficients. These are
//! source evaluations, not literal STEP curves or finite-incidence receipts.
use super::*;

/// The original affine radius, unwrapped phase and axial law along ONE owned
/// chart edge. The borrowed plan retains the symbolic contact expression.
#[derive(Debug)]
pub struct StockClipCurveLaw<'a> {
    plan: &'a StockClipPlan,
    curve: StockClipCurve,
    endpoints: [[I; 3]; 2],
    work: usize,
    bytes: usize,
}

/// Coefficient n is the nth original-whole-edge derivative divided by n!.
/// An interval parameter encloses coefficients at EVERY center in that domain.
#[derive(Clone, Copy, Debug)]
pub struct StockClipCurveTaylor {
    coefficients: [[[f64; 2]; 3]; 13],
    work: usize,
    bytes: usize,
}
impl StockClipCurveTaylor {
    pub fn coefficients_mm(&self) -> [[[f64; 2]; 3]; 13] {
        self.coefficients
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_budget_bytes(&self) -> usize {
        self.bytes
    }
}
impl StockClipPlan {
    pub fn curve_law(
        &self,
        curve: StockClipCurve,
        options: StockClipOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<StockClipCurveLaw<'_>, Error> {
        let retained = size_of::<Self>() + size_of::<StockClipCurveLaw<'_>>();
        let mut meter = Meter::new(options, retained, interrupt)?;
        let mut endpoints = [[I::point(0.); 3]; 2];
        for (index, t) in [0., 1.].into_iter().enumerate() {
            let p = self.curve_bounds(curve, [t, t], meter.child_options(options)?, interrupt)?;
            meter.charge(p.work_used())?;
            endpoints[index] = [
                interval(p.radius_mm())?,
                interval(p.phase_mm())?,
                interval(p.axial_mm())?,
            ];
        }
        // Every retained edge fixes u or v on the bilinear source chart.
        // Hence these SOURCE scalar restrictions are affine. No Cartesian
        // point fit or replacement rounded contact datum is being inferred.
        meter.charge(16)?;
        Ok(StockClipCurveLaw {
            plan: self,
            curve,
            endpoints,
            work: meter.used,
            bytes: retained + SCRATCH,
        })
    }
}
impl StockClipCurveLaw<'_> {
    pub fn plan(&self) -> &StockClipPlan {
        self.plan
    }
    pub fn curve(&self) -> StockClipCurve {
        self.curve
    }
    /// Evaluation caches of actual [radius,unwrapped phase,z] endpoint
    /// expressions. Equality of these enclosures is NEVER an ownership proof.
    pub fn affine_fields_mm(&self) -> [[[f64; 2]; 3]; 2] {
        self.endpoints.map(|p| p.map(bounds))
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_budget_bytes(&self) -> usize {
        self.bytes
    }
    pub fn taylor(
        &self,
        parameter: [f64; 2],
        options: StockClipOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<StockClipCurveTaylor, Error> {
        type S = Series<13>;
        let retained = self.bytes - SCRATCH + size_of::<StockClipCurveTaylor>();
        let mut meter = Meter::new(options, retained, interrupt)?;
        meter.charge(self.work)?;
        let t = S::variable(sample::unit_parameter(parameter)?);
        // Three affine compositions + sin/cos recurrence through order12,
        // coordinate products and coefficient publication, bounded below
        // 8192 scalar interval operations. Fixed arrays fit the16KiB scratch.
        meter.charge(8192)?;
        let fields: [S; 3] = std::array::from_fn(|i| {
            S::interval(self.endpoints[0][i])
                + S::interval(self.endpoints[1][i] - self.endpoints[0][i]) * t
        });
        let [radius, phase, z] = fields;
        let theta = (z - phase)
            * S::interval(I::point(2.) * pi() / I::point(self.plan.source.spec.pitch_mm));
        let [sin, cos] = theta.sin_cos();
        let xyz = [radius * cos, radius * sin, z];
        // Sealed source r,z<=1e7, L/p<128 and |q|<3p imply angular
        // derivative<900. All coefficients through12 and intermediate
        // recurrences are therefore finite (<1e50), even on interval domains.
        meter.charge(1)?;
        Ok(StockClipCurveTaylor {
            coefficients: std::array::from_fn(|n| xyz.map(|x| bounds(x.c[n]))),
            work: meter.used,
            bytes: retained + SCRATCH,
        })
    }
}

#[cfg(test)]
mod tests;
