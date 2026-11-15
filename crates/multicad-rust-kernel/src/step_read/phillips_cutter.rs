//! Explicitly source-known Phillips finite-graph readback. No foreign authoring
//! Domain, imported Cutter, independent mass or app admission is inferred.
use super::*;
use crate::{
    phillips_drive::{BoundaryTopology, CubicTrimSpan, Cutter},
    step::phillips_cutter::{
        curves::{self, Literal},
        support::{self, Kind, LiteralAxis, LiteralSupport, Pose},
    },
    tapered_thread::interval::I,
};
use helical::require;
use std::cell::Cell;
mod geometry;
mod graph;
mod physical;
pub use physical::{
    AnalyticEdgeLift, AnalyticEdgeLiftKind, AnalyticEdgeLiftOptions, AnalyticEdgeLifts,
    AnalyticSupportBounds as PhillipsAnalyticSupportBounds,
    AnalyticSupportOptions as PhillipsAnalyticSupportOptions, Axis as ParsedPhillipsAxis,
    Body as ParsedPhillipsBody, Boundary as ParsedPhillipsBoundary,
    CubicSupportBounds as PhillipsCubicSupportBounds,
    CubicSupportOptions as PhillipsCubicSupportOptions, Curve as ParsedPhillipsCurve,
    Edge as ParsedPhillipsEdge, EdgeSupportBound as PhillipsEdgeSupportBound,
    EndpointCouplingBounds as PhillipsEndpointCouplingBounds,
    EndpointCouplingOptions as PhillipsEndpointCouplingOptions,
    ExactPointLeaf as PhillipsExactPointLeaf, Face as ParsedPhillipsFace,
    ImplicitSupportEnclosure as PhillipsImplicitSupportEnclosure, IntersectionEdgeTube,
    IntersectionSection, IntersectionTubeOptions, IntersectionTubes,
    MaterialLoop as ParsedPhillipsMaterialLoop, MaterialUse as ParsedPhillipsMaterialUse,
    Mouth as ParsedPhillipsMouth, OriginalNumericLeaf as PhillipsOriginalNumericLeaf,
    OriginalNumericLeaves as PhillipsOriginalNumericLeaves,
    OuterCarrierBounds as PhillipsOuterCarrierBounds,
    OuterCarrierOptions as PhillipsOuterCarrierOptions, SewnEdgeEvaluation, SewnEdgeOptions,
    SewnEdgeParameter, SewnEdgeSource, SewnEdgeSources, Surface as ParsedPhillipsSurface,
    TensorCarrierBound as PhillipsTensorCarrierBound, Vertex as ParsedPhillipsVertex,
    VirtualOuterCap as ParsedPhillipsVirtualOuterCap, read_phillips_body,
    read_phillips_body_with_interrupt, read_phillips_screw_boundary,
    read_phillips_screw_boundary_with_interrupt,
};
#[cfg(test)]
mod tests;
/// Oriented source-decomposition view backed by actual parsed entities.
/// This never changes the file's graph or authorizes an imported recipe.
pub(in crate::step_read) struct SourceFaceView {
    pub surface: u32,
    pub forward: bool,
    pub uses: Vec<(u32, bool)>,
}
pub(in crate::step_read) fn verify_component(
    shell: &mut Shell<'_, '_>,
    views: &[SourceFaceView],
    source: &Cutter,
    pose: &Pose,
    tolerance: f64,
) -> Result<PhillipsSourceKnownProof> {
    let mut work = Work {
        remaining: 1_000_000,
        interrupt: None,
        cancelled: Cell::new(false),
    };
    let topology = source.boundary_topology().map_err(cutter_error)?;
    let graph = graph::read_views(shell, views, &topology, pose, tolerance, &mut work)?;
    geometry::verify(shell, &graph, &topology, pose, tolerance, &mut work)
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PhillipsReadError {
    Invalid,
    Unsupported,
    Limit,
    Cancelled,
}
impl From<ReadError> for PhillipsReadError {
    fn from(e: ReadError) -> Self {
        match e {
            ReadError::Invalid => Self::Invalid,
            ReadError::Unsupported => Self::Unsupported,
            ReadError::Limit => Self::Limit,
        }
    }
}
fn source_error(e: crate::step::PhillipsStepError) -> ReadError {
    match e {
        crate::step::PhillipsStepError::InvalidGeometry
        | crate::step::PhillipsStepError::InvalidPlacement => ReadError::Invalid,
        crate::step::PhillipsStepError::Unsupported => ReadError::Unsupported,
        crate::step::PhillipsStepError::Limit | crate::step::PhillipsStepError::Cancelled => {
            ReadError::Limit
        }
    }
}
fn cutter_error(e: crate::phillips_drive::Error) -> ReadError {
    match e {
        crate::phillips_drive::Error::Invalid => ReadError::Invalid,
        crate::phillips_drive::Error::Unresolved => ReadError::Unsupported,
        crate::phillips_drive::Error::Limit | crate::phillips_drive::Error::Cancelled => {
            ReadError::Limit
        }
    }
}
struct Work<'a> {
    remaining: usize,
    interrupt: Option<&'a dyn Fn() -> bool>,
    cancelled: Cell<bool>,
}
impl Work<'_> {
    fn interrupted(&self) -> bool {
        let yes = self.interrupt.is_some_and(|f| f());
        if yes {
            self.cancelled.set(true);
        }
        yes
    }
    fn step(&mut self, n: usize) -> Result<()> {
        if self.interrupted() {
            return Err(ReadError::Limit);
        }
        self.remaining = self.remaining.checked_sub(n).ok_or(ReadError::Limit)?;
        Ok(())
    }
}
/// Bounds against the supplied authored cutter, retaining actual file graph
/// identity and measured residuals. This is not an independent file-mass or
/// exact edge/support incidence certificate and cannot create an imported body.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhillipsSourceKnownProof {
    vertices: usize,
    edges: usize,
    faces: usize,
    maximum_vertex_error_mm: f64,
    maximum_analytic_support_error_mm: f64,
    maximum_polynomial_error_mm: f64,
    maximum_polynomial_support_error_mm: f64,
    certified_spans: usize,
}
impl PhillipsSourceKnownProof {
    pub fn vertex_count(&self) -> usize {
        self.vertices
    }
    pub fn edge_count(&self) -> usize {
        self.edges
    }
    pub fn face_count(&self) -> usize {
        self.faces
    }
    pub fn maximum_vertex_error_mm(&self) -> f64 {
        self.maximum_vertex_error_mm
    }
    pub fn maximum_analytic_support_error_mm(&self) -> f64 {
        self.maximum_analytic_support_error_mm
    }
    pub fn maximum_polynomial_error_mm(&self) -> f64 {
        self.maximum_polynomial_error_mm
    }
    pub fn maximum_polynomial_support_error_mm(&self) -> f64 {
        self.maximum_polynomial_support_error_mm
    }
    pub fn certified_span_count(&self) -> usize {
        self.certified_spans
    }
}
pub fn verify_phillips_cutter_source_known(
    input: &str,
    source: &Cutter,
    origin: Point,
    u: Point,
    v: Point,
    tolerance: f64,
) -> std::result::Result<PhillipsSourceKnownProof, PhillipsReadError> {
    verify_phillips_cutter_source_known_with_interrupt(input, source, origin, u, v, tolerance, None)
}
pub fn verify_phillips_cutter_source_known_with_interrupt(
    input: &str,
    source: &Cutter,
    origin: Point,
    u: Point,
    v: Point,
    tolerance: f64,
    interrupt: Option<&dyn Fn() -> bool>,
) -> std::result::Result<PhillipsSourceKnownProof, PhillipsReadError> {
    let mut work = Work {
        remaining: 1_000_000,
        interrupt,
        cancelled: Cell::new(false),
    };
    let result = verify(input, source, origin, u, v, tolerance, &mut work);
    result.map_err(|e| {
        if work.cancelled.get() {
            PhillipsReadError::Cancelled
        } else {
            e.into()
        }
    })
}
fn verify(
    input: &str,
    source: &Cutter,
    origin: Point,
    u: Point,
    v: Point,
    tolerance: f64,
    work: &mut Work<'_>,
) -> Result<PhillipsSourceKnownProof> {
    work.step(1)?;
    if !tolerance.is_finite() || tolerance <= 0. {
        return Err(ReadError::Invalid);
    }
    let tolerance = tolerance.min(1e-7);
    let pose = Pose::new(origin, u, v).map_err(source_error)?;
    let topology = source
        .boundary_topology_with_cancel(&|| work.interrupted())
        .map_err(cutter_error)?;
    let entities = part21::parse(input)?;
    let mut shell = Shell::new(&entities, tolerance, true);
    let solids = shell.represented_solids(true)?;
    let faces = shell.member_faces(solids[0])?;
    let graph = graph::read(&mut shell, &faces, &topology, &pose, tolerance, work)?;
    geometry::verify(&shell, &graph, &topology, &pose, tolerance, work)
}

pub use physical::{
    DerivedCavityOptions, DerivedCavityRepresentative, DerivedEdgeCorrespondence,
    DerivedFaceCorrespondence, DerivedVertexCorrespondence,
};

pub use physical::ClosedDerivedCavityCertificate;

pub use physical::{
    CommonMouthChart, CommonMouthHeadEnvelope, CommonMouthOptions, CommonTensorTransfer,
};

pub use physical::{
    ClosedHeadOptions, ClosedRepairedHead, HeadSupportTransfer, RawHeadSupportKind,
};

pub use physical::{
    AnalyticThreadOptions, AnalyticThreadRepair, AnalyticThreadTransfer, AuthoredPeriodicOptions,
    CoherentDegree11ThreadCarrier, CommonChartPeriodicThreadRepair, ContactThreadCarrier,
    Degree11ContactOptions, Degree11DomainEnd, Degree11RimContact, Degree11RimContactId,
    Degree11RimKnot, Degree11SignDomain, Degree11ThreadContacts, Degree11ThreadFace,
    Degree11ThreadOptions, DerivedThreadCarrier, DerivedThreadOptions, DyadicParameter,
    OccupiedDerivedThread, OwnedPeriodicMaterial, PeriodicAuthoredApproximation,
    PeriodicClosedMaterial, PeriodicCurveOwner, PeriodicCurveQuery, PeriodicFaceProperties,
    PeriodicFiniteCurve, PeriodicFiniteCurves, PeriodicFiniteKind, PeriodicFiniteOptions,
    PeriodicGeometryApproximation, PeriodicLiteralCertificate, PeriodicMaterialFace,
    PeriodicMaterialMatch, PeriodicMaterialMatchOptions, PeriodicMaterialOptions,
    PeriodicMaterialProvenance, PeriodicMaterialVolume, PeriodicOwnedMaterialQuery,
    PeriodicOwnedOptions, PeriodicProposalOrigin, PeriodicShoulderOptions, PeriodicShoulderRepair,
    PeriodicThreadBounds, PeriodicThreadCorrespondence, PeriodicThreadRepair, PeriodicThreadSource,
    PeriodicVolumeOptions, RimContact, RimContactId, RimDomainEnd, RimDomainRole, RimKnot,
    RimSignDomain, ThreadCarrierTransfer, ThreadContactOptions,
};

pub use physical::{ClosedMaterialOptions, DerivedClosedMaterial, MaterialFaceSource};

pub use physical::{MaterialVolume, MaterialVolumeOptions};

pub use physical::{
    DelegatedWorkspaceBounds, JoinedMaterialCandidate, JoinedMeshOptions, OwnedDerivedMaterial,
    OwnedMaterialOptions, OwnedMaterialQuery, OwnedMaterialQueryError, OwnedMaterialQueryOptions,
};

pub use physical::{
    PeriodicMaterialRead, PeriodicReadError, PeriodicReadOptions, PeriodicReadStage,
    read_periodic_phillips_material, read_periodic_phillips_material_with_progress,
};
