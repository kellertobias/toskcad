//! Source-known shared tensor candidates. No decoded file or solid admission.
use super::*;
use crate::axial_pipe_thread::authored::{
    StockClipBand, StockClipPatch, StockClipPatchLaw, StockClipSlab,
};
mod refine;
mod tensor;
/// One complete inherited-edge plus surface request. The aggregate point
/// cap can be reduced by a caller, never raised above262144.
#[derive(Clone, Copy, Debug)]
pub struct SurfaceOptions {
    pub query: CurveOptions,
    pub maximum_control_points: usize,
}
impl From<CurveOptions> for SurfaceOptions {
    fn from(query: CurveOptions) -> Self {
        Self {
            query,
            maximum_control_points: 262_144,
        }
    }
}
impl Default for SurfaceOptions {
    fn default() -> Self {
        CurveOptions::default().into()
    }
}
#[derive(Clone, Copy, Debug)]
pub struct SurfaceCellCertificate {
    domain: [[f64; 2]; 2],
    position_error: f64,
    derivative_error: [f64; 2],
}
impl SurfaceCellCertificate {
    pub fn domain(&self) -> [[f64; 2]; 2] {
        self.domain
    }
    pub fn position_error_mm(&self) -> f64 {
        self.position_error
    }
    pub fn original_derivative_error_mm(&self) -> [f64; 2] {
        self.derivative_error
    }
}
#[derive(Debug)]
pub struct SurfacePatch {
    source: StockClipPatch,
    spans: [usize; 2],
    poles: Vec<[f64; 3]>,
    cells: Vec<SurfaceCellCertificate>,
}
impl SurfacePatch {
    pub fn source_patch(&self) -> StockClipPatch {
        self.source
    }
    pub fn spans(&self) -> [usize; 2] {
        self.spans
    }
    pub fn poles_mm(&self) -> &[[f64; 3]] {
        &self.poles
    }
    pub fn cells(&self) -> &[SurfaceCellCertificate] {
        &self.cells
    }
}
#[derive(Debug)]
pub struct SurfacePlan<'a> {
    curves: CurvePlan<'a>,
    patches: Vec<SurfacePatch>,
    work: usize,
    bytes: usize,
}
fn band_index(band: StockClipBand) -> usize {
    match band {
        StockClipBand::Lower => 0,
        StockClipBand::Crest => 1,
        StockClipBand::Upper => 2,
        StockClipBand::Root => 3,
    }
}
fn slab_index(slab: StockClipSlab) -> usize {
    match slab {
        StockClipSlab::Before => 0,
        StockClipSlab::After => 1,
    }
}
impl<'a> SurfacePlan<'a> {
    pub fn new(
        mut curves: CurvePlan<'a>,
        limits: SurfaceOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, Error> {
        let options = limits.query;
        options.validate()?;
        if limits.maximum_control_points == 0 || limits.maximum_control_points > 262_144 {
            return Err(Error::Invalid);
        }
        if curves
            .curves
            .iter()
            .any(|curve| curve.spans.len() > options.maximum_spans_per_curve)
        {
            return Err(Error::Limit);
        }
        let mut work = Work::new(options.maximum_work, interrupt)?;
        // The original edge request is inherited ONCE, not silently renewed.
        work.charge(curves.work)?;
        let mut memory = Memory::new(curves.bytes + size_of::<Self>(), options.maximum_bytes)?;
        let mut u = [1; 4];
        let mut v = [1; 2];
        for curve in &curves.curves {
            match curve.source {
                StockClipCurve::Ring { band, .. } => {
                    u[band_index(band)] = u[band_index(band)].max(curve.spans.len())
                }
                StockClipCurve::Axial { slab, .. } => {
                    v[slab_index(slab)] = v[slab_index(slab)].max(curve.spans.len())
                }
            }
        }
        for curve in &mut curves.curves {
            let count = match curve.source {
                StockClipCurve::Ring { band, .. } => u[band_index(band)],
                StockClipCurve::Axial { slab, .. } => v[slab_index(slab)],
            };
            refine::curve(curve, count, &mut memory, &mut work)?;
        }
        curves.work = work.used;
        curves.bytes = memory.used - size_of::<Self>();
        let basis = curves::basis(&mut work)?;
        let mut patches = memory.allocate(8, &mut work)?;
        let mut control_count = 0usize;
        for source in curves.source.patches() {
            let spans = [u[band_index(source.band)], v[slab_index(source.slab)]];
            let dimensions = spans.map(|n| 11 * n + 1);
            let count = dimensions[0] * dimensions[1];
            control_count += count;
            if control_count > limits.maximum_control_points {
                return Err(Error::Limit);
            }
            let law = curves
                .source
                .patch_law(source, source_options(&work)?, interrupt)
                .map_err(source_error)?;
            work.charge(law.work_used())?;
            let scaled = tensor::scaled_basis(spans, &basis, &mut work)?;
            let mut poles = memory.allocate(count, &mut work)?;
            poles.resize(count, [0.; 3]);
            for (side, edge) in source.coedges().into_iter().enumerate() {
                let curve = curves
                    .curves
                    .iter()
                    .find(|c| c.source == edge.curve)
                    .expect("original source edge");
                let edge_count = if side % 2 == 0 {
                    dimensions[1]
                } else {
                    dimensions[0]
                };
                for i in 0..edge_count {
                    let span = (i / 11).min(curve.spans.len() - 1);
                    let point = curve.spans[span].poles[i - 11 * span];
                    let [x, y] = match side {
                        0 => [0, i],
                        1 => [i, dimensions[1] - 1],
                        2 => [dimensions[0] - 1, i],
                        _ => [i, 0],
                    };
                    poles[x * dimensions[1] + y] = point;
                    work.charge(1)?;
                }
            }
            let mut cells = memory.allocate(spans[0] * spans[1], &mut work)?;
            for i in 0..spans[0] {
                for j in 0..spans[1] {
                    let domain = [
                        [i as f64 / spans[0] as f64, (i + 1) as f64 / spans[0] as f64],
                        [j as f64 / spans[1] as f64, (j + 1) as f64 / spans[1] as f64],
                    ];
                    let candidate = tensor::cell(&law, domain, &scaled, &mut work, interrupt)?;
                    let mut actual = [[[0.; 3]; 12]; 12];
                    for (a, row) in actual.iter_mut().enumerate() {
                        for (b, pole) in row.iter_mut().enumerate() {
                            let x = 11 * i + a;
                            let y = 11 * j + b;
                            let boundary = x == 0
                                || y == 0
                                || x + 1 == dimensions[0]
                                || y + 1 == dimensions[1];
                            let inherited = (a == 0 && i > 0) || (b == 0 && j > 0);
                            if !boundary && !inherited {
                                poles[x * dimensions[1] + y] =
                                    candidate.ideal[a][b].map(I::midpoint);
                            }
                            *pole = poles[x * dimensions[1] + y];
                        }
                    }
                    let certificate = candidate.certify(actual, &mut work)?;
                    if certificate.position_error > options.position_error_mm
                        || certificate
                            .derivative_error
                            .iter()
                            .any(|e| *e > options.derivative_error_mm)
                    {
                        return Err(Error::Unresolved);
                    }
                    cells.push(certificate);
                }
            }
            patches.push(SurfacePatch {
                source,
                spans,
                poles,
                cells,
            });
        }
        work.charge(1)?;
        Ok(Self {
            curves,
            patches,
            work: work.used,
            bytes: memory.used,
        })
    }
    pub fn curves(&self) -> &CurvePlan<'a> {
        &self.curves
    }
    pub fn patches(&self) -> &[SurfacePatch] {
        &self.patches
    }
    /// Whole original-edge and surface construction work, under ONE quota.
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_budget_bytes(&self) -> usize {
        self.bytes
    }
}
