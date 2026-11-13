//! Authored Phillips H cutter with finite true C-C corner prisms.
//! Finite boundary and manufacturing publication require additional proofs.
use crate::tapered_thread::interval::I;
use std::f64::consts::PI;
mod cavity;
mod classify;
pub use cavity::{Cavity, HostEnvelopeProof};
mod events;
mod finite_section;
pub use events::{
    BoundaryEvents, BoundaryMesh, CubicTrimSpan, EventRole, FiniteTrim, PatchAtlas, PatchBounds,
    PatchKey, QuadErrorBound, SourceEvent, SurfacePatch, TrimErrorBound, WallKind,
};
pub use finite_section::{
    EdgeRole, SectionBoundary, SectionCurve, SectionEdge, SectionVertex, VertexRole,
};
mod source;
mod topology;
pub use topology::{
    AnalyticCurve, BoundaryEdge, BoundaryFace, BoundaryTopology, BoundaryVertex, CurveKey, EdgeUse,
    FaceKey, LineSupport, MeridianStation, Surface, TopologyCertificate, VertexKey,
};
#[cfg(test)]
mod tests;
#[cfg(test)]
pub(crate) use tests::fixture as test_fixture;
mod volume;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Spec {
    pub mouth_diameter_mm: f64,
    pub mouth_radius_mm: f64,
    pub b_width_mm: f64,
    pub g_diameter_mm: f64,
    pub tilt_angle_deg: f64,
    pub outer_flank_angle_deg: f64,
    pub floor_angle_deg: f64,
    pub outer_corner_angle_deg: f64,
    pub inner: Inner,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Inner {
    Rounded { radius_mm: f64 },
    ShallowV { angle_deg: f64, e_width_mm: f64 },
}
/// The finite corner prism, measured relative to the cutter mouth plane.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Domain {
    pub corner_start_depth_mm: f64,
    pub corner_end_depth_mm: f64,
    /// Length along each outer section ray, not its horizontal projection.
    pub corner_reach_mm: f64,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Error {
    Invalid,
    Unresolved,
    Limit,
    Cancelled,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Classification {
    Inside,
    Outside,
    UnresolvedBoundary,
}
/// Bounds enclose the exact authored expression, not fitted support geometry.
pub type Bounds = [f64; 2];
pub type Point2Bounds = [Bounds; 2];
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Curve {
    Line,
    CircularArc {
        center_mm_bounds: Point2Bounds,
        radius_mm: f64,
        sweep_rad_bounds: Bounds,
    },
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EnvelopeRole {
    EntryRoundover,
    OuterFlank,
    Floor,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MeridianSegment {
    pub start_mm_bounds: Point2Bounds,
    pub end_mm_bounds: Point2Bounds,
    pub curve: Curve,
    pub role: EnvelopeRole,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SectionRole {
    InnerRound,
    InnerUpper,
    OuterUpper,
    Closure,
    OuterLower,
    InnerLower,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SectionSegment {
    pub start_mm_bounds: Point2Bounds,
    pub end_mm_bounds: Point2Bounds,
    pub curve: Curve,
    pub role: SectionRole,
}
#[derive(Clone, Debug)]
pub struct Envelope {
    segments: [MeridianSegment; 3],
    entry: I,
    join: I,
    flank: I,
    tip: I,
    sin_a: I,
    cos_a: I,
    sin_c: I,
    cos_c: I,
}
impl Envelope {
    pub fn segments(&self) -> &[MeridianSegment; 3] {
        &self.segments
    }
    pub fn entry_radius_mm_bounds(&self) -> Bounds {
        bounds(self.entry)
    }
    pub fn roundover_join_depth_mm_bounds(&self) -> Bounds {
        bounds(self.join)
    }
    pub fn flank_join_depth_mm_bounds(&self) -> Bounds {
        bounds(self.flank)
    }
    pub fn tip_depth_mm_bounds(&self) -> Bounds {
        bounds(self.tip)
    }
}
#[derive(Clone, Debug)]
pub struct Corner {
    quarter_turn: u8,
    domain: Domain,
    section: [Option<SectionSegment>; 5],
    sin_beta: I,
    cos_beta: I,
    sin_theta: I,
    cos_theta: I,
    b_depth: I,
    join_s: I,
    join_y: I,
    far_s: I,
    b_half: I,
    inner: SourceInner,
}
#[derive(Clone, Debug)]
enum SourceInner {
    Rounded { radius_mm: f64 },
    ShallowV { sine: I, cosine: I },
}
impl Corner {
    pub fn quarter_turns(&self) -> u8 {
        self.quarter_turn
    }
    pub fn domain(&self) -> Domain {
        self.domain
    }
    pub fn section_segments(&self) -> impl Iterator<Item = &SectionSegment> {
        self.section.iter().flatten()
    }
    /// Source direction normal to the true C-C section, pointing into depth.
    pub fn extrusion_axis_bounds(&self) -> [Bounds; 3] {
        let [x, y] = rotate([-self.sin_beta, I::point(0.)], self.quarter_turn);
        [bounds(x), bounds(y), bounds(self.cos_beta)]
    }
    /// C-C coordinates (s,y) and the finite extrusion parameter in depth units.
    pub fn section_coordinates_mm_bounds(&self, p: [f64; 3]) -> Result<[Bounds; 3], Error> {
        if !p.into_iter().all(f64::is_finite) {
            return Err(Error::Invalid);
        }
        let values = self.coordinates(p.map(I::point));
        if values.into_iter().all(I::valid) {
            Ok(values.map(bounds))
        } else {
            Err(Error::Unresolved)
        }
    }
}
/// Bounded volume of the exact authored source expressions. This is not a
/// physical-file mass receipt or a finite manufacturing-boundary certificate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VolumeBounds {
    pub volume_mm3: [f64; 2],
    pub section_evaluations: usize,
    /// Outward sum of midpoint value-enclosure widths times their cell lengths.
    pub midpoint_evaluation_width_mm3: f64,
    /// Outward sum of signed curvature-remainder widths and whole-area event
    /// ranges. Diagnostics are not independent additive volume intervals.
    pub curvature_and_event_width_mm3: f64,
}
#[derive(Clone, Debug)]
pub struct Cutter {
    spec: Spec,
    domain: Domain,
    envelope: Envelope,
    corners: [Corner; 4],
}
impl PartialEq for Cutter {
    fn eq(&self, other: &Self) -> bool {
        self.spec == other.spec && self.domain == other.domain
    }
}
impl Cutter {
    pub fn new(spec: Spec, domain: Domain) -> Result<Self, Error> {
        Self::new_with_cancel(spec, domain, &|| false)
    }
    pub fn new_with_cancel(
        spec: Spec,
        domain: Domain,
        cancelled: &dyn Fn() -> bool,
    ) -> Result<Self, Error> {
        source::construct(spec, domain, cancelled)
    }
    /// Real support charts and whole-cell derivative bounds for detached mesh
    /// construction. A patch atlas alone does not authorize manufacturing.
    pub fn surface_patches(&self) -> Result<PatchAtlas<'_>, Error> {
        self.surface_patches_with_cancel(&|| false)
    }
    pub fn surface_patches_with_cancel(
        &self,
        cancelled: &dyn Fn() -> bool,
    ) -> Result<PatchAtlas<'_>, Error> {
        events::patches::construct(self, cancelled)
    }
    /// Closed source incidence with real support ownership. Manufacturing and
    /// literal emitted-file incidence require their separate certificates.
    pub fn boundary_topology(&self) -> Result<BoundaryTopology<'_>, Error> {
        self.boundary_topology_with_cancel(&|| false)
    }
    pub fn boundary_topology_with_cancel(
        &self,
        cancelled: &dyn Fn() -> bool,
    ) -> Result<BoundaryTopology<'_>, Error> {
        topology::construct(self, cancelled)
    }
    /// Shared expression-owned events and finite trim domains of the canonical
    /// three-dimensional boundary. General event orders remain unsupported.
    pub fn boundary_events(&self) -> Result<BoundaryEvents<'_>, Error> {
        self.boundary_events_with_cancel(&|| false)
    }
    pub fn boundary_events_with_cancel(
        &self,
        cancelled: &dyn Fn() -> bool,
    ) -> Result<BoundaryEvents<'_>, Error> {
        events::construct(self, cancelled)
    }
    /// A sealed oriented wire of actual authored section supports. Shared IDs
    /// preserve expression ownership; these bounds do not authorize meshing or
    /// independently certify literal foreign-file support incidence.
    pub fn finite_section(&self, depth_mm: f64) -> Result<SectionBoundary<'_>, Error> {
        self.finite_section_with_cancel(depth_mm, &|| false)
    }
    pub fn finite_section_with_cancel(
        &self,
        depth_mm: f64,
        cancelled: &dyn Fn() -> bool,
    ) -> Result<SectionBoundary<'_>, Error> {
        finite_section::construct(self, depth_mm, cancelled)
    }
    /// Integrates eligible Phillips cutters after proving their finite caps,
    /// closure, and adjacent corner cutouts inactive/disjoint throughout the
    /// cavity. Other authored domains return `Unresolved`.
    pub fn volume_bounds_mm3(
        &self,
        width_mm3: f64,
        work_limit: usize,
    ) -> Result<VolumeBounds, Error> {
        self.volume_bounds_with_cancel(width_mm3, work_limit, &|| false)
    }
    /// Work counts each section-cell evaluation, including initial cells.
    /// Limits above one million evaluations exceed this routine’s resource cap.
    pub fn volume_bounds_with_cancel(
        &self,
        width_mm3: f64,
        work_limit: usize,
        cancelled: &dyn Fn() -> bool,
    ) -> Result<VolumeBounds, Error> {
        volume::compute(self, width_mm3, work_limit, cancelled)
    }
    pub fn spec(&self) -> Spec {
        self.spec
    }
    pub fn domain(&self) -> Domain {
        self.domain
    }
    pub fn envelope(&self) -> &Envelope {
        &self.envelope
    }
    pub fn corners(&self) -> &[Corner; 4] {
        &self.corners
    }
    /// Conservative box for the cavity at and below its mouth, not an exact BBox.
    pub fn bounding_box_mm(&self) -> [[f64; 3]; 2] {
        let r = self.envelope.entry.hi;
        [[-r, -r, 0.], [r, r, self.envelope.tip.hi]]
    }
    pub fn classify_point_mm(&self, p: [f64; 3]) -> Result<Classification, Error> {
        classify::point(self, p)
    }
}
fn bounds(value: I) -> Bounds {
    [value.lo, value.hi]
}
fn point2(p: [I; 2]) -> Point2Bounds {
    p.map(bounds)
}
fn rotate([x, y]: [I; 2], quarter: u8) -> [I; 2] {
    match quarter {
        0 => [x, y],
        1 => [-y, x],
        2 => [-x, -y],
        _ => [y, -x],
    }
}
