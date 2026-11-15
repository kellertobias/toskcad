//! Candidate curves of the original clipped NPT atlas. These source-known
//! C0/C1 certificates do NOT admit a file, solid, material or decoded mass.
pub use super::insert::Error;
use super::insert::resources::{Memory, Work};
use crate::{
    axial_pipe_thread::{
        Error as SourceError,
        authored::{
            StockClipCurve, StockClipCurveLaw, StockClipOptions, StockClipPlan, StockClipVertex,
        },
    },
    tapered_thread::interval::I,
};
use std::mem::size_of;
mod curves;
mod surfaces;
pub use surfaces::{SurfaceCellCertificate, SurfaceOptions, SurfacePatch, SurfacePlan};
mod writer;
pub use writer::{CandidateFile, EncodingReceipt, write_source_candidate};
#[cfg(test)]
mod tests;

// Source Taylor fixed workspace (16KiB), basis, ideal/literal control arrays
// and endpoint/derivative temporaries together remain below64KiB.
const SCRATCH: usize = 64 * 1024;
#[derive(Clone, Copy, Debug)]
pub struct CurveOptions {
    pub position_error_mm: f64,
    /// Absolute error in the ORIGINAL whole-source-edge parameter.
    pub derivative_error_mm: f64,
    /// Uniform power-of-two subdivision limit, at most256.
    pub maximum_spans_per_curve: usize,
    pub maximum_work: usize,
    /// Borrowed original plan, retained result and construction workspace.
    pub maximum_bytes: usize,
}
impl Default for CurveOptions {
    fn default() -> Self {
        Self {
            position_error_mm: 1e-7,
            derivative_error_mm: 1e-6,
            maximum_spans_per_curve: 256,
            maximum_work: 16_000_000,
            maximum_bytes: 64 * 1024 * 1024,
        }
    }
}
impl CurveOptions {
    fn validate(self) -> Result<(), Error> {
        if !self.position_error_mm.is_finite()
            || self.position_error_mm <= 0.
            || !self.derivative_error_mm.is_finite()
            || self.derivative_error_mm <= 0.
            || !self.maximum_spans_per_curve.is_power_of_two()
            || self.maximum_spans_per_curve > 256
            || self.maximum_work == 0
            || self.maximum_work > 16_000_000
            || self.maximum_bytes == 0
            || self.maximum_bytes > 64 * 1024 * 1024
        {
            Err(Error::Invalid)
        } else {
            Ok(())
        }
    }
}
fn source_error(e: SourceError) -> Error {
    // These private calls use ONLY the sealed plan's twelve vertices/twenty
    // curves/eight patches, finite closed dyadic domains and validated positive
    // remaining work with the fixed64MiB source-query allowance. Admission already
    // handled unsupported/unresolved geometry. The source enforces p>=.001,
    // dimensions<=1e7 and L/p<128, so all retained query arithmetic is finite.
    // A different error here is an internal invariant break, not a new user
    // geometry refusal. Public source/domain validation remains unchanged.
    assert!(
        matches!(e, SourceError::Limit | SourceError::Cancelled),
        "sealed clipped source query invariant: {e:?}"
    );
    if matches!(e, SourceError::Limit) {
        Error::Limit
    } else {
        Error::Cancelled
    }
}
fn source_options(work: &Work<'_>) -> Result<StockClipOptions, Error> {
    if work.remaining() == 0 {
        return Err(Error::Limit);
    }
    Ok(StockClipOptions {
        maximum_work: work.remaining(),
        ..Default::default()
    })
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CurvePointOwner {
    Vertex(StockClipVertex),
    /// Reduced dyadic station on the original source curve, never an inferred
    /// equality of rounded coordinates or interval overlap.
    Interior {
        curve: StockClipCurve,
        numerator: u16,
        denominator: u16,
    },
}
#[derive(Clone, Copy, Debug)]
pub struct CarrierVertex {
    owner: StockClipVertex,
    position: [f64; 3],
}
impl CarrierVertex {
    pub fn owner(&self) -> StockClipVertex {
        self.owner
    }
    pub fn position_mm(&self) -> [f64; 3] {
        self.position
    }
}
#[derive(Clone, Copy, Debug)]
pub struct BezierSpan {
    domain: [f64; 2],
    owners: [CurvePointOwner; 2],
    poles: [[f64; 3]; 12],
    position_error: f64,
    derivative_error: f64,
}
impl BezierSpan {
    pub fn domain(&self) -> [f64; 2] {
        self.domain
    }
    pub fn endpoint_owners(&self) -> [CurvePointOwner; 2] {
        self.owners
    }
    /// Degree11 Cartesian coefficients. This candidate owns them; an eventual
    /// writer must reuse endpoint owner IDs and qualify actual numeric leaves.
    pub fn poles_mm(&self) -> &[[f64; 3]; 12] {
        &self.poles
    }
    pub fn position_error_mm(&self) -> f64 {
        self.position_error
    }
    pub fn original_derivative_error_mm(&self) -> f64 {
        self.derivative_error
    }
}
#[derive(Debug)]
pub struct CurveCarrier {
    source: StockClipCurve,
    spans: Vec<BezierSpan>,
}
impl CurveCarrier {
    pub fn source_curve(&self) -> StockClipCurve {
        self.source
    }
    pub fn spans(&self) -> &[BezierSpan] {
        &self.spans
    }
}
#[derive(Debug)]
pub struct CurvePlan<'a> {
    source: &'a StockClipPlan,
    vertices: [CarrierVertex; 12],
    curves: Vec<CurveCarrier>,
    work: usize,
    bytes: usize,
}
impl<'a> CurvePlan<'a> {
    pub fn new(
        source: &'a StockClipPlan,
        options: CurveOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, Error> {
        options.validate()?;
        let fixed = size_of::<StockClipPlan>() + size_of::<Self>() + SCRATCH;
        let mut memory = Memory::new(fixed, options.maximum_bytes)?;
        let mut work = Work::new(options.maximum_work, interrupt)?;
        let blank = CarrierVertex {
            owner: source.vertices()[0],
            position: [0.; 3],
        };
        let mut vertices = [blank; 12];
        for (index, owner) in source.vertices().into_iter().enumerate() {
            let p = source
                .vertex_bounds(owner, source_options(&work)?, interrupt)
                .map_err(source_error)?;
            work.charge(p.work_used())?;
            vertices[index] = CarrierVertex {
                owner,
                position: p.position_mm().map(|b| b[0].midpoint(b[1])),
            };
        }
        let basis = curves::basis(&mut work)?;
        let mut curves = memory.allocate(20, &mut work)?;
        for key in source.curves() {
            let law = source
                .curve_law(key, source_options(&work)?, interrupt)
                .map_err(source_error)?;
            work.charge(law.work_used())?;
            let spans = curves::fit(
                &law,
                &vertices,
                &basis,
                options,
                &mut memory,
                &mut work,
                interrupt,
            )?;
            curves.push(CurveCarrier { source: key, spans });
        }
        work.charge(1)?;
        Ok(Self {
            source,
            vertices,
            curves,
            work: work.used,
            bytes: memory.used,
        })
    }
    pub fn source(&self) -> &StockClipPlan {
        self.source
    }
    pub fn vertices(&self) -> &[CarrierVertex; 12] {
        &self.vertices
    }
    pub fn curves(&self) -> &[CurveCarrier] {
        &self.curves
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    /// Conservative whole-query peak; includes the original borrowed plan
    /// exactly once and fixed live construction scratch, not allocator RSS.
    pub fn retained_budget_bytes(&self) -> usize {
        self.bytes
    }
    /// The chosen candidate interface plane, not a replacement source datum.
    /// Ring coefficients at that station share this SAME literal height.
    pub fn contact_plane_height_mm(&self) -> f64 {
        let b = self.source.contact_height_bounds_mm();
        b[0].midpoint(b[1])
    }
    /// Outward SOURCE symbolic-contact to candidate-plane displacement.
    /// This is not an independently decoded physical-file uncertainty.
    pub fn contact_plane_error_mm(&self) -> f64 {
        let b = self.source.contact_height_bounds_mm();
        (I::point(self.contact_plane_height_mm()) - I::new(b[0], b[1])).magnitude()
    }
}
