//! Safe solid operations in millimeters. Display meshes are detached copies;
//! callers rebuild exact solids from source rather than editing triangles.
use crate::{Mesh, Shape, ffi, mesh::valid_parts, rust_prism_topology::PlacedPolygon};
use multicad_core::{
    document::EntityId,
    references::{Origin, ReferenceFrame, TopologySnapshot},
};
use multicad_rust_kernel::{
    Circle, Rational, cylinder::AxialCylinder, polygon::PolygonPrism, prism::Prism,
};
use std::path::Path;
use std::rc::{Rc, Weak};
mod annulus_body;
mod arithmetic;
mod authored_sweep;
mod axial_stack_body;
mod bevel_body;
pub(crate) mod bore_treatment;
mod bored_prism_boolean;
mod coaxial;
mod coaxial_counterbore;
mod contained_void_cut;
pub use bore_treatment::ConicalSeatGeometry;
pub(crate) mod circular_pocket_body;
mod clipped_npt_body;
mod pipe_hardware_body;
pub use multicad_rust_kernel::axial_pipe_thread::authored::hardware::Kind as PipeHardwareKind;
use pipe_hardware_body::PipeHardwareBody;
pub use pipe_hardware_body::{PipeThreadFactory, PipeThreadFactoryStats};
mod actual_clipped_head_body;
use actual_clipped_head_body::ActualClippedHeadBody;
pub use actual_clipped_head_body::{ActualClippedHeadRuntimeEvidence, ActualHeadStepVerification};
mod literal_npt_body;
use literal_npt_body::LiteralNptBody;
pub use literal_npt_body::LiteralNptRuntimeEvidence;
mod external_npt_stock_body;
mod line8_nut_body;
use external_npt_stock_body::ExternalNptStockBody;
pub use external_npt_stock_body::ExternalNptStockTransfer;
use line8_nut_body::Line8NutBody;
pub use line8_nut_body::{Line8NutCut, Line8NutTransfer};
pub(crate) mod construction;
mod cosmetic_fastener_body;
mod cosmetic_profiled_body;
use clipped_npt_body::ClippedNptBody;
use cosmetic_profiled_body::CosmeticProfiledBody;
mod entry_chamfer_body;
mod wood_insert_body;
use wood_insert_body::WoodInsertBody;
pub use wood_insert_body::{
    PartitionedStepVerification, WoodInsertData, WoodInsertStepVerification,
};
mod wood_insert_material;
pub use wood_insert_material::{WoodInsertFaceProperties, WoodInsertMaterial};
mod hanger_material;
pub use hanger_material::{HangerFaceProperties, HangerMaterial};
mod hanger_body;
mod source_face_measure;
use hanger_body::HangerBody;
pub use hanger_body::{HangerData, HealedHangerStepVerification, HealedSupportStepVerification};
mod cross_bore_body;
mod gear_pair_boolean;
mod gear_probe;
mod helical_body;
mod helical_tube_body;
mod herringbone_body;
pub(crate) mod imported_brep_body;
mod internal_ring_body;
mod linear_body;
mod machining_boolean;
mod nonplanar_loft_body;
mod normal_relief_body;
mod normal_sweep_body;
mod oblique_bore_body;
mod orthogonal_body;
mod parallel_cylinder;
mod parallel_extrusion_intersection;
mod planar_boolean_body;
mod radial_boundary_body;
mod repaired_phillips_body;
use repaired_phillips_body::RepairedPhillipsBody;
pub use repaired_phillips_body::{
    DecodedPeriodicPhillipsStepVerification, Degree11RepairedPhillipsMeshReceipt,
    Degree11RepairedPhillipsQuery, Degree11RepairedPhillipsQueryOptions,
    PeriodicPhillipsStepCandidateResources, PeriodicRepairedPhillipsMeshReceipt,
    PeriodicRepairedPhillipsQuery, RepairedPhillipsMeshBudget, RepairedPhillipsMeshReceipt,
    RepairedPhillipsQuery, RepairedPhillipsQueryOptions,
};
mod revolution_body;
mod rounded_blind_bore_body;
mod rounded_cylinder_body;
mod rounded_slot_body;
mod rounded_thread_body;
mod rounded_thread_bore_body;
mod rust_compound;
mod separated_boolean;
mod sphere_body;
mod spline_sweep_body;
mod spur_gear_body;
mod stack_union_boolean;
mod stacked_prism_boolean;
mod step_import;
mod through_profile_cut;
mod torus_body;
mod warped_union_body;
mod wood_contact;
mod worm_body;
use annulus_body::AnnulusBody;
use axial_stack_body::AxialStackBody;
use bevel_body::BevelBody;
pub use circular_pocket_body::CircularPocketTool;
use cosmetic_fastener_body::CosmeticFastenerBody;
use entry_chamfer_body::EntryChamferBody;
pub use entry_chamfer_body::LiteralMeshEmbedding;
pub use entry_chamfer_body::{EntryChamferBoundaryProof, EntryChamferData, EntryChamferEdge};
use helical_body::HelicalBody;
pub use helical_body::PhillipsMeshVolumePair;
pub use helical_body::{
    BlindThreadCut, BlindThreadMating, BlindThreadSeat, InternalThreadData,
    ReceivingThreadIntervals,
};
pub use helical_tube_body::HelicalTube;
use helical_tube_body::HelicalTubeBody;
use herringbone_body::HerringboneBody;
use internal_ring_body::InternalRingBody;
use linear_body::LinearBody;
use nonplanar_loft_body::NonplanarLoftBody;
use normal_relief_body::NormalReliefBody;
mod spline_cap_relief_body;
use cross_bore_body::CrossBoreBody;
pub use normal_relief_body::NormalReliefMotif;
use normal_sweep_body::NormalSweepBody;
pub use normal_sweep_body::SweepOpening;
use oblique_bore_body::ObliqueBoreBody;
use orthogonal_body::OrthogonalBody;
use orthogonal_body::seats::OrthogonalSeatBody;
pub use orthogonal_body::seats::{OrthogonalCountersink, OrthogonalCountersinkMeasures};
pub use orthogonal_body::slots::RoundedSlotTool;
use planar_boolean_body::PlanarBooleanBody;
use radial_boundary_body::RadialBoundaryBody;
pub use radial_boundary_body::{RetainedAxialPipeThreadSource, RetainedHelicalThreadSource};
pub use radial_boundary_body::{RetainedStepMemberMass, RetainedStepVerification};
pub use revolution_body::CylindricalBlindPocketSpec;
use revolution_body::RevolutionBody;
use rounded_cylinder_body::RoundedCylinderBody;
use rounded_slot_body::RoundedSlotBody;
use rounded_thread_body::RoundedThreadBody;
pub use rounded_thread_body::{RoundedThreadData, RoundedThreadMassProperties};
mod internal_thread_treatment_body;
pub use imported_brep_body::AssemblyOccurrence;
use imported_brep_body::ImportedBrepBody;
use internal_thread_treatment_body::InternalThreadTreatmentBody;
pub use internal_thread_treatment_body::{
    InternalThreadFaceProperties, InternalThreadTreatmentCut, InternalThreadTreatmentData,
};
use rounded_blind_bore_body::RoundedBlindBoreBody;
pub use rounded_blind_bore_body::{
    RoundedBlindBoreData, RoundedBlindMeshBudget, RoundedBlindMeshReceipt,
};
use rounded_thread_bore_body::RoundedThreadBoreBody;
pub use rounded_thread_bore_body::RoundedThreadBoreData;
use rust_compound::RustCompound;
use sphere_body::SphereBody;
use spline_cap_relief_body::SplineCapReliefBody;
use spline_sweep_body::SplineSweepBody;
pub use spur_gear_body::GearFlankWitness;
use spur_gear_body::SpurGearBody;
use torus_body::TorusBody;
use warped_union_body::WarpedUnionBody;
use worm_body::WormBody;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GeometryErrorKind {
    Input,
    Native,
    InvalidResult,
    Limit,
    Cancelled,
    Stale,
    Busy,
    Unsupported,
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GeometryError {
    pub operation: &'static str,
    pub kind: GeometryErrorKind,
    pub message: String,
}
impl std::fmt::Display for GeometryError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}: {}", self.operation, self.message)
    }
}
impl std::error::Error for GeometryError {}
fn error(
    operation: &'static str,
    kind: GeometryErrorKind,
    message: impl Into<String>,
) -> GeometryError {
    GeometryError {
        operation,
        kind,
        message: message.into(),
    }
}
fn input(message: impl Into<String>) -> GeometryError {
    error("validate", GeometryErrorKind::Input, message)
}
fn native(operation: &'static str, e: impl std::fmt::Display) -> GeometryError {
    error(operation, GeometryErrorKind::Native, e.to_string())
}

fn rational_coordinate(value: f64, tolerance: f64) -> Option<Rational> {
    // Use the first continued-fraction convergent inside the error budget
    // so floating-point projection noise does not inflate denominators.
    let mut remainder = value.abs();
    let mut previous = [0_i128, 1];
    let mut current = [1_i128, 0];
    for _ in 0..32 {
        if !remainder.is_finite() || remainder > i64::MAX as f64 {
            break;
        }
        let whole = remainder.floor();
        let next = std::array::from_fn(|i| whole as i128 * current[i] + previous[i]);
        let [Ok(n), Ok(d)] = next.map(i64::try_from) else {
            break;
        };
        let Ok(candidate) = Rational::new(if value.is_sign_negative() { -n } else { n }, d) else {
            break;
        };
        if (candidate.as_f64() - value).abs() <= tolerance / 32. {
            return Some(candidate);
        }
        previous = current;
        current = next;
        remainder = 1. / (remainder - whole);
    }
    None
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GeometryPolicy {
    pub linear_tolerance_mm: f64,
    /// Additional Boolean fuzz is explicit. Zero requests no added gap closing.
    pub boolean_fuzzy_mm: f64,
    pub minimum_feature_mm: f64,
    pub maximum_coordinate_mm: f64,
}
impl Default for GeometryPolicy {
    fn default() -> Self {
        Self {
            linear_tolerance_mm: 1e-7,
            boolean_fuzzy_mm: 0.,
            minimum_feature_mm: 1e-5,
            maximum_coordinate_mm: 1e7,
        }
    }
}
impl GeometryPolicy {
    pub fn validate(self) -> Result<(), GeometryError> {
        if ![
            self.linear_tolerance_mm,
            self.boolean_fuzzy_mm,
            self.minimum_feature_mm,
            self.maximum_coordinate_mm,
        ]
        .iter()
        .all(|v| v.is_finite())
            || !(1e-7..=1e-3).contains(&self.linear_tolerance_mm)
            || !(0. ..=self.linear_tolerance_mm).contains(&self.boolean_fuzzy_mm)
            || self.minimum_feature_mm < self.linear_tolerance_mm * 10.
            || self.maximum_coordinate_mm <= self.minimum_feature_mm
            || self.maximum_coordinate_mm > 1e7
        {
            return Err(input(
                "Invalid numerical policy: tolerance, fuzz, feature size or model bounds",
            ));
        }
        Ok(())
    }
    pub(crate) fn point(self, p: [f64; 3]) -> Result<(), GeometryError> {
        if p.iter()
            .any(|v| !v.is_finite() || v.abs() > self.maximum_coordinate_mm)
        {
            Err(input("Coordinate is nonfinite or outside model bounds"))
        } else {
            Ok(())
        }
    }
    pub(crate) fn length(self, v: f64) -> Result<(), GeometryError> {
        if !v.is_finite() || v < self.minimum_feature_mm || v > self.maximum_coordinate_mm {
            Err(input(
                "Length is below minimum feature size or outside model bounds",
            ))
        } else {
            Ok(())
        }
    }
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlaneFrame {
    pub origin: [f64; 3],
    pub u: [f64; 3],
    pub v: [f64; 3],
}
impl Default for PlaneFrame {
    fn default() -> Self {
        Self {
            origin: [0.; 3],
            u: [1., 0., 0.],
            v: [0., 1., 0.],
        }
    }
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.into_iter().zip(b).map(|(a, b)| a * b).sum()
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn point(p: [f64; 3]) -> ffi::Point3 {
    ffi::Point3 {
        x: p[0],
        y: p[1],
        z: p[2],
    }
}
fn xyz(p: ffi::Point3) -> [f64; 3] {
    [p.x, p.y, p.z]
}
impl PlaneFrame {
    fn checked_metric_extent(
        self,
        extent: [f64; 3],
        policy: GeometryPolicy,
    ) -> Result<(), GeometryError> {
        self.checked(policy)?;
        let norm = |p: [f64; 3]| p[0].hypot(p[1]).hypot(p[2]);
        let unit = |p: [f64; 3]| {
            let n = norm(p);
            p.map(|x| x / n)
        };
        let u = unit(self.u);
        let n = unit(self.normal());
        let v = cross(n, u);
        let residual = [self.u, self.v, self.normal()]
            .into_iter()
            .zip([u, v, n])
            .zip(extent)
            .map(|((actual, proper), size)| {
                norm(std::array::from_fn(|i| actual[i] - proper[i])) * size
            })
            .sum::<f64>();
        if !residual.is_finite() || residual > policy.linear_tolerance_mm / 4. {
            return Err(error(
                "placement",
                GeometryErrorKind::Unsupported,
                "Frame metric residual exceeds the geometry precision budget at this size",
            ));
        }
        Ok(())
    }
    // Callers supply the analytic enclosure in the body's retained local
    // chart, not an authored descriptor that may have been recentered.
    fn checked_metric_bounds(
        self,
        bounds: Bounds,
        policy: GeometryPolicy,
    ) -> Result<(), GeometryError> {
        self.checked_metric_extent(
            std::array::from_fn(|i| bounds.min[i].abs().max(bounds.max[i].abs())),
            policy,
        )
    }
    pub fn normal(self) -> [f64; 3] {
        cross(self.u, self.v)
    }
    pub fn point(self, p: [f64; 2]) -> [f64; 3] {
        std::array::from_fn(|i| self.origin[i] + self.u[i] * p[0] + self.v[i] * p[1])
    }
    pub(crate) fn checked(self, policy: GeometryPolicy) -> Result<ffi::Frame3, GeometryError> {
        policy.point(self.origin)?;
        if self.u.iter().chain(&self.v).any(|v| !v.is_finite())
            || (dot(self.u, self.u) - 1.).abs() > 1e-10
            || (dot(self.v, self.v) - 1.).abs() > 1e-10
            || dot(self.u, self.v).abs() > 1e-10
        {
            return Err(input("Plane axes must be finite and orthonormal"));
        }
        Ok(ffi::Frame3 {
            origin: point(self.origin),
            u: point(self.u),
            v: point(self.v),
        })
    }
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidTransform {
    pub rotation: [[f64; 3]; 3],
    pub translation_mm: [f64; 3],
}
impl Default for RigidTransform {
    fn default() -> Self {
        Self {
            rotation: [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
            translation_mm: [0.; 3],
        }
    }
}
impl RigidTransform {
    fn vector(self, vector: [f64; 3]) -> [f64; 3] {
        self.rotation.map(|row| dot(row, vector))
    }
    fn frame(self, frame: PlaneFrame) -> PlaneFrame {
        let origin = self.vector(frame.origin);
        PlaneFrame {
            origin: std::array::from_fn(|i| origin[i] + self.translation_mm[i]),
            u: self.vector(frame.u),
            v: self.vector(frame.v),
        }
    }
    fn validate(self, policy: GeometryPolicy) -> Result<(), GeometryError> {
        policy.point(self.translation_mm)?;
        if self.rotation.iter().flatten().any(|v| !v.is_finite()) {
            return Err(input("Rotation must be finite"));
        }
        for i in 0..3 {
            for j in 0..3 {
                if (dot(self.rotation[i], self.rotation[j]) - if i == j { 1. } else { 0. }).abs()
                    > 1e-10
                {
                    return Err(input(
                        "Rotation must be orthonormal; scaling and shear are not rigid placements",
                    ));
                }
            }
        }
        if (dot(cross(self.rotation[0], self.rotation[1]), self.rotation[2]) - 1.).abs() > 1e-10 {
            return Err(input("Rotation cannot reflect the solid"));
        }
        Ok(())
    }
}
#[derive(Clone, Debug, PartialEq)]
pub enum ProfileEdge {
    Line {
        a: [f64; 2],
        b: [f64; 2],
    },
    Arc {
        center: [f64; 2],
        a: [f64; 2],
        b: [f64; 2],
        clockwise: bool,
    },
    Circle {
        center: [f64; 2],
        radius_mm: f64,
    },
}
fn distance(a: [f64; 2], b: [f64; 2]) -> f64 {
    (a[0] - b[0]).hypot(a[1] - b[1])
}
impl ProfileEdge {
    fn ends(&self) -> ([f64; 2], [f64; 2]) {
        match *self {
            Self::Line { a, b } | Self::Arc { a, b, .. } => (a, b),
            Self::Circle { center, radius_mm } => {
                let p = [center[0] + radius_mm, center[1]];
                (p, p)
            }
        }
    }
    fn arc(&self) -> Option<(f64, f64, f64)> {
        if let Self::Arc {
            center,
            a,
            b,
            clockwise,
        } = *self
        {
            let r = distance(a, center);
            let start = (a[1] - center[1]).atan2(a[0] - center[0]);
            let end = (b[1] - center[1]).atan2(b[0] - center[0]);
            let sweep = if clockwise {
                -((start - end).rem_euclid(std::f64::consts::TAU))
            } else {
                (end - start).rem_euclid(std::f64::consts::TAU)
            };
            Some((r, start, sweep))
        } else {
            None
        }
    }
    fn area(&self) -> f64 {
        match *self {
            Self::Line { a, b } => (a[0] * b[1] - b[0] * a[1]) * 0.5,
            Self::Circle { radius_mm, .. } => std::f64::consts::PI * radius_mm * radius_mm,
            Self::Arc { center, a, b, .. } => {
                let (r, _, sweep) = self.arc().expect("arc");
                (center[0] * (b[1] - a[1]) - center[1] * (b[0] - a[0]) + r * r * sweep) * 0.5
            }
        }
    }
    fn reversed(&self) -> Self {
        match *self {
            Self::Line { a, b } => Self::Line { a: b, b: a },
            Self::Arc {
                center,
                a,
                b,
                clockwise,
            } => Self::Arc {
                center,
                a: b,
                b: a,
                clockwise: !clockwise,
            },
            Self::Circle { .. } => self.clone(),
        }
    }
    fn checked(
        &self,
        frame: PlaneFrame,
        policy: GeometryPolicy,
    ) -> Result<ffi::CurveEdge, GeometryError> {
        let (a, b) = self.ends();
        policy.point(frame.point(a))?;
        policy.point(frame.point(b))?;
        let mut edge = ffi::CurveEdge {
            kind: 0,
            a: point(frame.point(a)),
            mid: point([0.; 3]),
            b: point(frame.point(b)),
            radius: 0.,
            source_index: 0,
        };
        match *self {
            Self::Line { .. } => policy.length(distance(a, b))?,
            Self::Circle { center, radius_mm } => {
                policy.length(radius_mm)?;
                policy.point(frame.point(center))?;
                edge.kind = 2;
                edge.a = point(frame.point(center));
                edge.radius = radius_mm;
            }
            Self::Arc { center, .. } => {
                policy.point(frame.point(center))?;
                let (r, start, sweep) = self.arc().expect("arc");
                policy.length(r)?;
                if (distance(b, center) - r).abs() > policy.linear_tolerance_mm
                    || distance(a, b) <= policy.linear_tolerance_mm
                    || r * sweep.abs() < policy.minimum_feature_mm
                {
                    return Err(input("Arc endpoints/radius are inconsistent or degenerate"));
                }
                let mid = start + sweep * 0.5;
                let mid = [center[0] + r * mid.cos(), center[1] + r * mid.sin()];
                policy.point(frame.point(mid))?;
                edge.kind = 1;
                edge.mid = point(frame.point(mid));
            }
        }
        Ok(edge)
    }
}
#[derive(Clone, Debug, PartialEq)]
pub struct Profile {
    pub frame: PlaneFrame,
    pub outer: Vec<ProfileEdge>,
    pub holes: Vec<Vec<ProfileEdge>>,
}
impl Profile {
    pub fn rectangle(width: f64, depth: f64) -> Self {
        let points = [[0., 0.], [width, 0.], [width, depth], [0., depth]];
        Self {
            frame: PlaneFrame::default(),
            outer: (0..4)
                .map(|i| ProfileEdge::Line {
                    a: points[i],
                    b: points[(i + 1) % 4],
                })
                .collect(),
            holes: vec![],
        }
    }
    fn checked(
        &self,
        policy: GeometryPolicy,
    ) -> Result<(ffi::Frame3, Vec<ffi::CurveEdge>, Vec<u32>), GeometryError> {
        let frame = self.frame.checked(policy)?;
        let count = self.outer.len() + self.holes.iter().map(Vec::len).sum::<usize>();
        if count > 4096 || self.holes.len() > 128 {
            return Err(error(
                "profile",
                GeometryErrorKind::Limit,
                "Profile exceeds edge/hole budget",
            ));
        }
        let mut edges = vec![];
        let mut offsets = vec![0];
        let mut source_offset = 0;
        for wire in std::iter::once(&self.outer).chain(&self.holes) {
            if wire.is_empty() {
                return Err(input("Profile wires cannot be empty"));
            }
            for (i, e) in wire.iter().enumerate() {
                if matches!(e, ProfileEdge::Circle { .. }) && wire.len() != 1 {
                    return Err(input("A full-circle wire must contain only that circle"));
                }
                if distance(e.ends().1, wire[(i + 1) % wire.len()].ends().0)
                    > policy.linear_tolerance_mm
                {
                    return Err(input(
                        "Profile wire is not closed/connected within tolerance",
                    ));
                }
                e.checked(self.frame, policy)?;
            }
            let area = wire.iter().map(ProfileEdge::area).sum::<f64>();
            if !area.is_finite() || area.abs() < policy.minimum_feature_mm.powi(2) {
                return Err(input("Profile has degenerate signed area"));
            }
            if area < 0. {
                for (i, e) in wire.iter().enumerate().rev() {
                    let mut edge = e.reversed().checked(self.frame, policy)?;
                    edge.source_index = (source_offset + i) as u32;
                    edges.push(edge);
                }
            } else {
                for (i, e) in wire.iter().enumerate() {
                    let mut edge = e.checked(self.frame, policy)?;
                    edge.source_index = (source_offset + i) as u32;
                    edges.push(edge);
                }
            }
            source_offset += wire.len();
            offsets.push(edges.len() as u32);
        }
        Ok((frame, edges, offsets))
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BooleanOperation {
    Union,
    Subtract,
    Intersect,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Bounds {
    pub min: [f64; 3],
    pub max: [f64; 3],
}

/// STEP solid count, decoded volume in mm³ and world bounds in mm.
pub type StepCompoundMetrics = (u32, f64, [[f64; 3]; 2]);
/// Earned planar support as origin and outward normal.
pub type FacePlane = ([f64; 3], [f64; 3]);
#[derive(Clone, Debug, PartialEq)]
pub struct Face {
    /// Evaluation-local ordinal, never a persistent topology reference.
    pub index: u32,
    pub plane: Option<([f64; 3], [f64; 3])>,
    pub center_mm: [f64; 3],
    pub area_mm2: f64,
}
#[derive(Clone, Debug, PartialEq)]
pub struct DisplayEdge {
    /// Evaluation-local index into the B-rep edge map, not a mesh edge.
    pub index: u32,
    pub points_mm: Vec<[f64; 3]>,
}
#[derive(Clone, Debug, PartialEq)]
pub struct DisplayVertex {
    pub index: u32,
    pub position_mm: [f64; 3],
}
pub struct TopologyDisplay {
    pub topology: TopologySnapshot,
    pub edges: Vec<DisplayEdge>,
    pub vertices: Vec<DisplayVertex>,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MeshOptions {
    pub linear_deflection_mm: f64,
    pub angular_deflection_rad: f64,
    pub maximum_triangles: u32,
}
impl Default for MeshOptions {
    fn default() -> Self {
        Self {
            linear_deflection_mm: 0.1,
            angular_deflection_rad: 0.3,
            maximum_triangles: 1_000_000,
        }
    }
}
/// Owned solid geometry and its provenance. Each solid stays on its owner thread.
/// ```compile_fail
/// fn requires_send<T: Send>() {}
/// requires_send::<multicad_kernel::Solid>();
/// ```
/// ```compile_fail
/// fn requires_sync<T: Sync>() {}
/// requires_sync::<multicad_kernel::Solid>();
/// ```
pub struct Solid {
    body: SolidBody,
    policy: GeometryPolicy,
    pub(crate) provenance: Vec<multicad_core::references::TopologyElement>,
    pub(crate) rust_origins: Vec<Origin>,
    // The last authored placement retains its original operands. A shared
    // placement can then be factored without inverting rounded world points.
    retained_placement: Option<Rc<RetainedPlacement>>,
}

struct RetainedPlacement {
    source: Solid,
    transform: RigidTransform,
    geometry: RetainedGeometryKey,
}

struct RetainedGeometryKey {
    body: Weak<dyn std::any::Any>,
    frame: PlaneFrame,
}

impl Solid {
    fn retained_geometry_key(&self) -> Option<RetainedGeometryKey> {
        fn key<T: std::any::Any>(body: &Rc<T>, frame: PlaneFrame) -> RetainedGeometryKey {
            let body: Weak<T> = Rc::downgrade(body);
            RetainedGeometryKey { body, frame }
        }
        Some(match &self.body {
            SolidBody::Prism(body, frame) => key(body, *frame),
            SolidBody::PolygonPrism(body) => key(body, body.frame),
            SolidBody::Linear(body, frame) => key(body, *frame),
            SolidBody::Orthogonal(body, frame) => key(body, *frame),
            SolidBody::OrthogonalSeats(body, frame) => key(body, *frame),
            SolidBody::PlanarBoolean(body, frame) => key(body, *frame),
            _ => return None,
        })
    }

    fn current_retained_placement(&self) -> Option<&RetainedPlacement> {
        let retained = self.retained_placement.as_deref()?;
        let actual = self.retained_geometry_key()?;
        // Weak ownership prevents allocation-address reuse. Copy-on-write
        // geometry edits, chart changes and metadata edits invalidate reuse.
        if !actual.body.ptr_eq(&retained.geometry.body)
            || actual.frame != retained.geometry.frame
            || self.rust_origins != retained.source.rust_origins
            || self.provenance
                != retained
                    .source
                    .transformed_author_metadata(retained.transform)
        {
            return None;
        }
        Some(retained)
    }
}

enum SolidBody {
    Occt(Shape),
    Prism(Rc<Prism>, PlaneFrame),
    PolygonPrism(Rc<PlacedPolygon>),
    AxialCylinder(Rc<AxialCylinder>, PlaneFrame, Option<Rc<CylinderSketch>>),
    Annulus(Rc<AnnulusBody>, PlaneFrame),
    CosmeticFastener(Rc<CosmeticFastenerBody>, PlaneFrame),
    CosmeticProfiled(Rc<CosmeticProfiledBody>, PlaneFrame),
    ClippedNpt(Rc<ClippedNptBody>),
    PipeHardware(Rc<PipeHardwareBody>),
    LiteralNpt(Rc<LiteralNptBody>),
    ActualClippedHead(Rc<ActualClippedHeadBody>),
    ExternalNptStock(Rc<ExternalNptStockBody>),
    Line8Nut(Rc<Line8NutBody>),
    Helical(Rc<HelicalBody>, PlaneFrame),
    RoundedThread(Rc<RoundedThreadBody>, PlaneFrame),
    EntryChamfer(Rc<EntryChamferBody>, PlaneFrame),
    RoundedThreadBore(Rc<RoundedThreadBoreBody>, PlaneFrame),
    RoundedBlindBore(Rc<RoundedBlindBoreBody>, PlaneFrame),
    InternalThreadTreatment(Rc<InternalThreadTreatmentBody>, PlaneFrame),
    RepairedPhillips(Rc<RepairedPhillipsBody>, PlaneFrame),
    WoodInsert(Rc<WoodInsertBody>, PlaneFrame),
    Hanger(Rc<HangerBody>, PlaneFrame),
    HelicalTubeCut(Rc<HelicalTubeBody>, PlaneFrame),
    Linear(Rc<LinearBody>, PlaneFrame),
    Revolved(Rc<RevolutionBody>, PlaneFrame),
    RoundedSlot(Rc<RoundedSlotBody>, PlaneFrame),
    NormalSweep(Rc<NormalSweepBody>, PlaneFrame),
    NormalRelief(Rc<NormalReliefBody>, PlaneFrame),
    SplineCapRelief(Rc<SplineCapReliefBody>, PlaneFrame),
    RadialBoundary(Rc<RadialBoundaryBody>, PlaneFrame),
    NonplanarLoft(Rc<NonplanarLoftBody>, PlaneFrame),
    SplineSweep(Rc<SplineSweepBody>, PlaneFrame),
    SpurGear(Rc<SpurGearBody>, PlaneFrame),
    InternalRing(Rc<InternalRingBody>, PlaneFrame),
    Herringbone(Rc<HerringboneBody>, PlaneFrame),
    Bevel(Rc<BevelBody>, PlaneFrame),
    Worm(Rc<WormBody>, PlaneFrame),
    Sphere(Rc<SphereBody>, PlaneFrame),
    Torus(Rc<TorusBody>, PlaneFrame),
    ObliqueBore(Rc<ObliqueBoreBody>, PlaneFrame),
    CrossBore(Rc<CrossBoreBody>, PlaneFrame),
    RoundedCylinder(Rc<RoundedCylinderBody>, PlaneFrame),
    AxialStack(Rc<AxialStackBody>, PlaneFrame),
    RustCompound(Rc<RustCompound>),
    Orthogonal(Rc<OrthogonalBody>, PlaneFrame),
    ImportedBrep(Rc<ImportedBrepBody>),
    OrthogonalSeats(Rc<OrthogonalSeatBody>, PlaneFrame),
    PlanarBoolean(Rc<PlanarBooleanBody>, PlaneFrame),
    WarpedUnion(Rc<WarpedUnionBody>, PlaneFrame),
}

#[derive(Clone)]
pub(crate) struct CylinderSketch {
    pub(crate) origin: Origin,
    pub(crate) center_uv: [f64; 2],
    pub(crate) reverse: bool,
    pub(crate) cap_anchors: Option<[ReferenceFrame; 2]>,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum MeshUse {
    Display,
    Manufacturing,
}

fn polygon_metric_extent(prism: &PolygonPrism) -> [f64; 3] {
    let mut extent = [0., 0., prism.height().as_f64().abs()];
    // Polygon coordinates remain authored UV coordinates, unlike centered
    // cylinders/annuli. Interior holes cannot enlarge the outer envelope.
    for point in &prism.rings()[0] {
        for i in 0..2 {
            extent[i] = f64::max(extent[i], point[i].as_f64().abs());
        }
    }
    extent
}

impl Solid {
    /// Publish a straight-sided sketch extrusion with the original authored
    /// edge anchors, even if winding normalization reversed some edges.
    pub fn from_polygon_prism(
        prism: PolygonPrism,
        frame: PlaneFrame,
        feature: EntityId,
        origins: &[Origin],
        reversed: &[bool],
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        Self::from_polygon_prism_sweep(prism, frame, false, feature, origins, reversed, policy)
    }

    /// Keep the authored start and end cap identities when sweeping against
    /// the frame normal. The local prism always has positive height.
    pub fn from_polygon_prism_sweep(
        prism: PolygonPrism,
        frame: PlaneFrame,
        reverse_sweep: bool,
        feature: EntityId,
        origins: &[Origin],
        reversed: &[bool],
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        policy.validate()?;
        frame.checked(policy)?;
        policy.length(prism.height().as_f64())?;
        frame.checked_metric_extent(polygon_metric_extent(&prism), policy)?;
        let line_count = prism.rings().iter().map(Vec::len).sum::<usize>();
        let edge_count = line_count + prism.circular_holes().len();
        let unique: std::collections::BTreeSet<_> = origins.iter().collect();
        if feature.is_nil()
            || origins.len() != edge_count
            || reversed.len() != edge_count
            || reversed[line_count..].iter().any(|reversed| *reversed)
            || unique.len() != edge_count
            || origins.iter().any(|origin| {
                origin.source.is_nil()
                    || origin.role.trim().is_empty()
                    || format!("side/{}/{}", origin.source, origin.role).len() > 256
            })
        {
            return Err(input(
                "Every polygon edge requires a distinct valid authored origin",
            ));
        }
        for circle in prism.circular_holes() {
            policy.length(circle.radius().as_f64())?;
        }
        let placed = PlacedPolygon {
            prism: Rc::new(prism),
            frame,
            reverse_sweep,
            feature,
            edges: origins
                .iter()
                .cloned()
                .zip(reversed.iter().copied())
                .collect(),
            pockets: Rc::default(),
            pocket_origins: Rc::default(),
            bore_treatments: None,
        };
        for length in placed.side_lengths() {
            policy.length(length)?;
        }
        let bounds = placed.bounds();
        policy.point(bounds.min)?;
        policy.point(bounds.max)?;
        Ok(Self {
            body: SolidBody::PolygonPrism(Rc::new(placed)),
            policy,
            provenance: vec![],
            rust_origins: vec![],
            retained_placement: None,
        })
    }

    pub(crate) fn polygon_prism(&self) -> Option<&PlacedPolygon> {
        match &self.body {
            SolidBody::PolygonPrism(placed) => Some(placed),
            _ => None,
        }
    }

    /// Own a validated analytic prism in the Rust backend.
    pub fn from_rust_prism(prism: Prism, policy: GeometryPolicy) -> Result<Self, GeometryError> {
        Self::from_placed_prism(prism, PlaneFrame::default(), policy)
    }

    pub fn rust_box(
        feature: u128,
        size: [Rational; 3],
        frame: PlaneFrame,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        let prism = Prism::rectangular(feature, size)
            .map_err(|cause| input(format!("Rust box failed: {cause:?}")))?;
        Self::from_placed_prism(prism, frame, policy)
    }

    fn from_placed_prism(
        prism: Prism,
        frame: PlaneFrame,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        policy.validate()?;
        frame.checked(policy)?;
        if prism.feature() == 0 || prism.inner_loops().any(|(id, _)| id == 0) {
            return Err(input("Prism topology requires non-nil source identities"));
        }
        for size in prism.dimensions() {
            policy.length(size.as_f64())?;
        }
        // maximum includes a rounded end beyond the straight profile length.
        frame.checked_metric_extent(prism.maximum().map(Rational::as_f64), policy)?;
        if let Some(end) = prism.rounded_end() {
            policy.length(end.radius().as_f64())?;
        }
        for (_, circle) in prism.inner_loops() {
            policy.length(circle.radius().as_f64())?;
        }
        let solid = Self {
            body: SolidBody::Prism(Rc::new(prism), frame),
            policy,
            provenance: vec![],
            rust_origins: vec![],
            retained_placement: None,
        };
        let bounds = solid.bounds()?;
        policy.point(bounds.min)?;
        policy.point(bounds.max)?;
        Ok(solid)
    }

    pub fn rust_prism(&self) -> Option<&Prism> {
        self.rust_prism_frame().map(|(prism, _)| prism)
    }

    pub(crate) fn rust_prism_frame(&self) -> Option<(&Prism, PlaneFrame)> {
        match &self.body {
            SolidBody::Prism(prism, frame) => Some((prism, *frame)),
            _ => None,
        }
    }

    pub fn from_axial_cylinder(
        cylinder: AxialCylinder,
        frame: PlaneFrame,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        let [x, y] = cylinder.circle().center().map(Rational::as_f64);
        if frame.origin != [x, y, cylinder.bottom().as_f64()]
            || frame.u[2] != 0.
            || frame.v[2] != 0.
            || frame.normal()[2] <= 0.
        {
            return Err(input(
                "Rust cylinder requires a matching frame with positive Z normal",
            ));
        }
        Self::build_axial_cylinder(cylinder, frame, None, policy)
    }

    /// Place a circular sketch sweep along the frame normal. The frame
    /// starts at the first rim center and points along the positive height.
    /// A reverse sweep retains the authored cap attachment basis.
    pub fn from_axial_cylinder_sweep(
        cylinder: AxialCylinder,
        frame: PlaneFrame,
        center_uv: [f64; 2],
        reverse: bool,
        origin: Origin,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        if origin.source.is_nil()
            || origin.role.trim().is_empty()
            || format!("side/{}/{}", origin.source, origin.role).len() > 256
        {
            return Err(input(
                "Circular extrusion requires a valid authored curve origin",
            ));
        }
        policy.point([center_uv[0], center_uv[1], 0.])?;
        let sketch = CylinderSketch {
            origin,
            center_uv,
            reverse,
            cap_anchors: None,
        };
        Self::build_axial_cylinder(cylinder, frame, Some(Rc::new(sketch)), policy)
    }

    fn build_axial_cylinder(
        cylinder: AxialCylinder,
        frame: PlaneFrame,
        sketch: Option<Rc<CylinderSketch>>,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        policy.validate()?;
        if cylinder.feature() == 0 {
            return Err(input(
                "Cylinder topology requires a non-nil source identity",
            ));
        }
        policy.length(cylinder.circle().radius().as_f64())?;
        policy.length(cylinder.height().as_f64())?;
        frame.checked_metric_extent(
            [
                cylinder.circle().radius().as_f64(),
                cylinder.circle().radius().as_f64(),
                cylinder.height().as_f64(),
            ],
            policy,
        )?;
        let solid = Self {
            body: SolidBody::AxialCylinder(Rc::new(cylinder), frame, sketch),
            policy,
            provenance: vec![],
            rust_origins: vec![],
            retained_placement: None,
        };
        let bounds = solid.bounds()?;
        policy.point(bounds.min)?;
        policy.point(bounds.max)?;
        Ok(solid)
    }

    pub fn axial_cylinder(&self) -> Option<&AxialCylinder> {
        let (cylinder, frame) = self.axial_cylinder_frame()?;
        let [x, y] = cylinder.circle().center().map(Rational::as_f64);
        (frame.origin == [x, y, cylinder.bottom().as_f64()]
            && frame.u[2] == 0.
            && frame.v[2] == 0.
            && frame.normal()[2] > 0.)
            .then_some(cylinder)
    }

    /// World base-rim centre, unit axis, radius and height of a placed Rust
    /// cylinder whose frame starts at its first rim and points along its
    /// height. Reverse circular sweeps keep an authored cap basis instead and
    /// are not described here.
    pub fn rust_cylinder_axis(&self) -> Option<([f64; 3], [f64; 3], f64, f64)> {
        let (cylinder, frame) = self.axial_cylinder_frame()?;
        if self.cylinder_sketch().is_some_and(|sketch| sketch.reverse) {
            return None;
        }
        Some((
            frame.origin,
            frame.normal(),
            cylinder.circle().radius().as_f64(),
            cylinder.height().as_f64(),
        ))
    }

    pub(crate) fn axial_cylinder_frame(&self) -> Option<(&AxialCylinder, PlaneFrame)> {
        match &self.body {
            SolidBody::AxialCylinder(cylinder, frame, _) => Some((cylinder, *frame)),
            _ => None,
        }
    }

    pub(crate) fn cylinder_sketch(&self) -> Option<&CylinderSketch> {
        match &self.body {
            SolidBody::AxialCylinder(_, _, sketch) => sketch.as_deref(),
            _ => None,
        }
    }

    /// Subtract an exact axial through-cylinder in the Rust backend.
    pub fn subtract_rust_prism_cylinder(
        &self,
        feature: u128,
        circle: Circle,
        bottom: Rational,
        top: Rational,
    ) -> Result<Self, GeometryError> {
        self.try_subtract_rust_prism_through_cylinder(feature, circle, bottom, top)?
            .ok_or_else(|| error("cylindrical cut", GeometryErrorKind::Input, "NotThrough"))
    }

    fn try_subtract_rust_prism_through_cylinder(
        &self,
        feature: u128,
        circle: Circle,
        bottom: Rational,
        top: Rational,
    ) -> Result<Option<Self>, GeometryError> {
        let Some((prism, frame)) = self.rust_prism_frame() else {
            return Err(error(
                "cylindrical cut",
                GeometryErrorKind::Unsupported,
                "Expected a Rust prism",
            ));
        };
        if frame.u != PlaneFrame::default().u || frame.v != PlaneFrame::default().v {
            return Err(error(
                "cylindrical cut",
                GeometryErrorKind::Unsupported,
                "World-axial through-bores require a prism with an aligned source chart",
            ));
        }
        self.policy.point([0., 0., bottom.as_f64()])?;
        self.policy.point([0., 0., top.as_f64()])?;
        // Translate exact world operands into the retained source chart. Do not
        // use tolerance-based rational fitting here: even a very shallow floor
        // must keep a short cutter from being accepted as a through-bore.
        let offset = |axis: usize| {
            multicad_rust_kernel::parse_length_mm(&format!("{} mm", frame.origin[axis])).map_err(
                |_| {
                    error(
                        "cylindrical cut",
                        GeometryErrorKind::Unsupported,
                        "Prism translation cannot be represented in the exact length domain",
                    )
                },
            )
        };
        let subtract = |a: Rational, b: Rational| {
            a.checked_sub(b).map_err(|_| {
                error(
                    "cylindrical cut",
                    GeometryErrorKind::Limit,
                    "Translated cutter exceeds the exact arithmetic budget",
                )
            })
        };
        let center = circle.center();
        let circle = Circle::new(
            [
                subtract(center[0], offset(0)?)?,
                subtract(center[1], offset(1)?)?,
            ],
            circle.radius(),
        )
        .map_err(|reason| input(format!("Invalid translated cutter: {reason:?}")))?;
        let z = offset(2)?;
        let next =
            match prism.subtract_cylinder(feature, circle, subtract(bottom, z)?, subtract(top, z)?)
            {
                Ok(next) => next,
                // The general blind-pocket provider currently accepts only
                // rectilinear stock. Keep the through-cut refusal on rounded
                // stock instead of losing it in an unrelated slot provider.
                Err(multicad_rust_kernel::prism::PrismError::NotThrough)
                    if prism.rounded_end().is_none() =>
                {
                    return Ok(None);
                }
                Err(reason) => {
                    return Err(error(
                        "cylindrical cut",
                        GeometryErrorKind::Input,
                        format!("{reason:?}"),
                    ));
                }
            };
        Self::from_placed_prism(next, frame, self.policy).map(Some)
    }

    pub(crate) fn occt(&self) -> Result<&ffi::Shape, GeometryError> {
        match &self.body {
            SolidBody::Occt(shape) => Ok(&shape.0),
            _ => Err(error(
                "solid backend",
                GeometryErrorKind::Unsupported,
                "This operation is unavailable for the Rust solid",
            )),
        }
    }

    pub fn geometry_policy(&self) -> GeometryPolicy {
        self.policy
    }
    pub fn uses_rust_geometry(&self) -> bool {
        !matches!(self.body, SolidBody::Occt(_))
    }
    pub(crate) fn rust_compound(&self) -> Option<&RustCompound> {
        match &self.body {
            SolidBody::RustCompound(compound) => Some(compound),
            _ => None,
        }
    }
    pub(crate) fn orthogonal_seat_body(&self) -> Option<(&OrthogonalSeatBody, PlaneFrame)> {
        match &self.body {
            SolidBody::OrthogonalSeats(b, f) => Some((b, *f)),
            _ => None,
        }
    }
    pub(crate) fn orthogonal_body(&self) -> Option<(&OrthogonalBody, PlaneFrame)> {
        match &self.body {
            SolidBody::Orthogonal(body, frame) => Some((body, *frame)),
            _ => None,
        }
    }
    pub(crate) fn make_rust_compound(
        solids: &[&Solid],
        feature: EntityId,
        policy: GeometryPolicy,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, GeometryError> {
        RustCompound::assemble(solids, feature, policy, interrupt)
    }
    /// Another owner of the same immutable exact B-rep (and its provenance),
    /// without copying geometry. Kernel operations never modify their input
    /// shapes, so both owners observe identical geometry. Like every Solid
    /// the share stays on the thread that created it (`!Send`/`!Sync`).
    pub fn share(&self) -> Result<Self, GeometryError> {
        let body = match &self.body {
            SolidBody::ClippedNpt(body) => SolidBody::ClippedNpt(Rc::clone(body)),
            SolidBody::PipeHardware(body) => SolidBody::PipeHardware(Rc::clone(body)),
            SolidBody::LiteralNpt(body) => SolidBody::LiteralNpt(Rc::clone(body)),
            SolidBody::ActualClippedHead(body) => SolidBody::ActualClippedHead(Rc::clone(body)),
            SolidBody::ExternalNptStock(body) => SolidBody::ExternalNptStock(Rc::clone(body)),
            SolidBody::Line8Nut(body) => SolidBody::Line8Nut(Rc::clone(body)),
            SolidBody::Prism(prism, frame) => SolidBody::Prism(Rc::clone(prism), *frame),
            SolidBody::PolygonPrism(prism) => SolidBody::PolygonPrism(Rc::clone(prism)),
            SolidBody::AxialCylinder(cylinder, frame, sketch) => {
                SolidBody::AxialCylinder(Rc::clone(cylinder), *frame, sketch.clone())
            }
            SolidBody::Annulus(body, frame) => SolidBody::Annulus(Rc::clone(body), *frame),
            SolidBody::CosmeticFastener(body, frame) => {
                SolidBody::CosmeticFastener(Rc::clone(body), *frame)
            }
            SolidBody::CosmeticProfiled(body, frame) => {
                SolidBody::CosmeticProfiled(Rc::clone(body), *frame)
            }
            SolidBody::EntryChamfer(body, frame) => {
                SolidBody::EntryChamfer(Rc::clone(body), *frame)
            }
            SolidBody::Helical(body, frame) => SolidBody::Helical(Rc::clone(body), *frame),
            SolidBody::RoundedThreadBore(body, frame) => {
                SolidBody::RoundedThreadBore(Rc::clone(body), *frame)
            }
            SolidBody::RoundedBlindBore(body, frame) => {
                SolidBody::RoundedBlindBore(Rc::clone(body), *frame)
            }
            SolidBody::InternalThreadTreatment(body, frame) => {
                SolidBody::InternalThreadTreatment(Rc::clone(body), *frame)
            }
            SolidBody::RepairedPhillips(body, frame) => {
                SolidBody::RepairedPhillips(Rc::clone(body), *frame)
            }
            SolidBody::WoodInsert(body, frame) => SolidBody::WoodInsert(Rc::clone(body), *frame),
            SolidBody::Hanger(body, frame) => SolidBody::Hanger(Rc::clone(body), *frame),
            SolidBody::RoundedThread(body, frame) => {
                SolidBody::RoundedThread(Rc::clone(body), *frame)
            }
            SolidBody::HelicalTubeCut(body, frame) => {
                SolidBody::HelicalTubeCut(Rc::clone(body), *frame)
            }
            SolidBody::Linear(body, frame) => SolidBody::Linear(Rc::clone(body), *frame),
            SolidBody::Revolved(body, frame) => SolidBody::Revolved(Rc::clone(body), *frame),
            SolidBody::RoundedSlot(body, frame) => SolidBody::RoundedSlot(Rc::clone(body), *frame),
            SolidBody::NormalSweep(body, frame) => SolidBody::NormalSweep(Rc::clone(body), *frame),
            SolidBody::NormalRelief(body, frame) => {
                SolidBody::NormalRelief(Rc::clone(body), *frame)
            }
            SolidBody::SplineCapRelief(body, frame) => {
                SolidBody::SplineCapRelief(Rc::clone(body), *frame)
            }
            SolidBody::RadialBoundary(body, frame) => {
                SolidBody::RadialBoundary(Rc::clone(body), *frame)
            }
            SolidBody::SplineSweep(body, frame) => SolidBody::SplineSweep(Rc::clone(body), *frame),
            SolidBody::SpurGear(body, frame) => SolidBody::SpurGear(Rc::clone(body), *frame),
            SolidBody::InternalRing(body, frame) => {
                SolidBody::InternalRing(Rc::clone(body), *frame)
            }
            SolidBody::Herringbone(body, frame) => SolidBody::Herringbone(Rc::clone(body), *frame),
            SolidBody::Bevel(body, frame) => SolidBody::Bevel(Rc::clone(body), *frame),
            SolidBody::Worm(body, frame) => SolidBody::Worm(Rc::clone(body), *frame),
            SolidBody::Sphere(body, frame) => SolidBody::Sphere(Rc::clone(body), *frame),
            SolidBody::ImportedBrep(body) => SolidBody::ImportedBrep(Rc::clone(body)),
            SolidBody::Torus(body, frame) => SolidBody::Torus(Rc::clone(body), *frame),
            SolidBody::ObliqueBore(body, frame) => SolidBody::ObliqueBore(Rc::clone(body), *frame),
            SolidBody::CrossBore(body, frame) => SolidBody::CrossBore(Rc::clone(body), *frame),
            SolidBody::RoundedCylinder(body, frame) => {
                SolidBody::RoundedCylinder(Rc::clone(body), *frame)
            }
            SolidBody::AxialStack(body, frame) => SolidBody::AxialStack(Rc::clone(body), *frame),
            SolidBody::PlanarBoolean(body, frame) => {
                SolidBody::PlanarBoolean(Rc::clone(body), *frame)
            }
            SolidBody::NonplanarLoft(body, frame) => {
                SolidBody::NonplanarLoft(Rc::clone(body), *frame)
            }
            SolidBody::WarpedUnion(body, frame) => SolidBody::WarpedUnion(Rc::clone(body), *frame),
            SolidBody::RustCompound(compound) => SolidBody::RustCompound(Rc::clone(compound)),
            SolidBody::Orthogonal(body, frame) => SolidBody::Orthogonal(Rc::clone(body), *frame),
            SolidBody::OrthogonalSeats(body, frame) => {
                SolidBody::OrthogonalSeats(Rc::clone(body), *frame)
            }
            SolidBody::Occt(shape) => {
                let handle = ffi::share_shape(&shape.0).map_err(|e| native("share", e))?;
                if handle.is_null() {
                    return Err(error(
                        "share",
                        GeometryErrorKind::InvalidResult,
                        "Kernel returned an empty shape",
                    ));
                }
                SolidBody::Occt(Shape(handle))
            }
        };
        Ok(Self {
            body,
            policy: self.policy,
            provenance: self.provenance.clone(),
            rust_origins: self.rust_origins.clone(),
            retained_placement: self.retained_placement.clone(),
        })
    }
    pub(crate) fn adopt(
        handle: crate::OcctHandle,
        policy: GeometryPolicy,
        operation: &'static str,
    ) -> Result<Self, GeometryError> {
        if handle.is_null() {
            return Err(error(
                operation,
                GeometryErrorKind::InvalidResult,
                "Kernel returned an empty shape",
            ));
        }
        let solid = Self {
            body: SolidBody::Occt(Shape(handle)),
            policy,
            provenance: vec![],
            rust_origins: vec![],
            retained_placement: None,
        };
        if !ffi::valid(solid.occt()?).map_err(|e| native(operation, e))?
            || solid.solid_count()? == 0
        {
            return Err(error(
                operation,
                GeometryErrorKind::InvalidResult,
                "Result is not a valid volumetric solid",
            ));
        }
        let v = solid.volume_mm3()?;
        if !v.is_finite() || v <= 0. {
            return Err(error(
                operation,
                GeometryErrorKind::InvalidResult,
                "Result has no positive finite volume",
            ));
        }
        let bounds = solid.bounds()?;
        policy.point(bounds.min)?;
        policy.point(bounds.max)?;
        Ok(solid)
    }
    pub fn box_solid(size_mm: [f64; 3], policy: GeometryPolicy) -> Result<Self, GeometryError> {
        if !crate::OCCT_LINKED {
            return construction::box_solid(size_mm, EntityId::new(), policy)?
                .without_authored_provenance();
        }
        policy.validate()?;
        for v in size_mm {
            policy.length(v)?;
        }
        Self::adopt(
            ffi::extrude(size_mm[0], size_mm[1], size_mm[2]).map_err(|e| native("box", e))?,
            policy,
            "box",
        )
    }
    pub fn cylinder(
        frame: PlaneFrame,
        radius_mm: f64,
        height_mm: f64,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        if !crate::OCCT_LINKED {
            return construction::cylinder(frame, radius_mm, height_mm, EntityId::new(), policy)?
                .without_authored_provenance();
        }
        policy.validate()?;
        policy.length(radius_mm)?;
        policy.length(height_mm)?;
        Self::adopt(
            ffi::cylinder(&frame.checked(policy)?, radius_mm, height_mm)
                .map_err(|e| native("cylinder", e))?,
            policy,
            "cylinder",
        )
    }
    pub fn extrude_profile(
        profile: &Profile,
        height_mm: f64,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        if !crate::OCCT_LINKED {
            let feature = EntityId::new();
            return construction::profile(
                profile,
                height_mm,
                feature,
                &construction::origins(
                    feature,
                    profile.outer.len() + profile.holes.iter().map(Vec::len).sum::<usize>(),
                ),
                policy,
            )?
            .without_authored_provenance();
        }
        policy.validate()?;
        policy.length(height_mm.abs())?;
        let (frame, edges, offsets) = profile.checked(policy)?;
        Self::adopt(
            ffi::profile_extrude(&frame, &edges, &offsets, height_mm)
                .map_err(|e| native("extrude", e))?,
            policy,
            "extrude",
        )
    }
    pub fn boolean(
        &self,
        other: &Self,
        operation: BooleanOperation,
        policy: GeometryPolicy,
    ) -> Result<Option<Self>, GeometryError> {
        self.boolean_geometry(other, operation, policy)?
            .map(|result| result.retain_authored_boolean_sources(&[self, other]))
            .transpose()
    }

    fn boolean_geometry(
        &self,
        other: &Self,
        operation: BooleanOperation,
        policy: GeometryPolicy,
    ) -> Result<Option<Self>, GeometryError> {
        policy.validate()?;
        let mut cross_refusal = None;
        if operation == BooleanOperation::Subtract
            && policy.boolean_fuzzy_mm == 0.
            && policy == self.policy
            && policy == other.policy
        {
            if let Some(result) = self.subtract_coaxial_cylinder(other, policy) {
                return result;
            }
            // A cylinder oblique to plain box stock in the stock's own chart,
            // whatever the world placement of either operand.
            if let Some(result) = self.oblique_cylinder_cut(other, None) {
                return result;
            }
            // A cylinder perpendicular to a stepped shaft or plain cylinder.
            // Its refusal stands only when no later provider qualifies.
            match self.cross_bore_cut(other, None) {
                Some(Err(error)) if error.kind == GeometryErrorKind::Unsupported => {
                    cross_refusal = Some(error);
                }
                Some(result) => return result,
                None => {}
            }
            if let Some(tool) = other.axial_cylinder() {
                if let Some(target) = self.polygon_prism() {
                    return match Self::subtract_polygon_cylinder(target, tool, policy) {
                        Err(error) if error.kind == GeometryErrorKind::Unsupported => self
                            .subtract_machining_tool(other, operation, policy)
                            .or_else(|| self.coaxial_counterbore_cut(other, operation, policy))
                            .unwrap_or(Err(error)),
                        result => result.map(Some),
                    };
                }
                if self.rust_prism().is_some()
                    && let Some(result) = self.try_subtract_rust_prism_through_cylinder(
                        tool.feature(),
                        tool.circle(),
                        tool.bottom(),
                        tool.top(),
                    )?
                {
                    return Ok(Some(result));
                }
            }
            if let Some(result) = self.subtract_orthogonal_cylinder(other, policy) {
                // A profile with curved holes is outside the rectilinear grid;
                // the exact polygon pocket provider may still own the cut.
                return match result {
                    Err(error) if error.kind == GeometryErrorKind::Unsupported => self
                        .subtract_machining_tool(other, operation, policy)
                        .or_else(|| self.coaxial_counterbore_cut(other, operation, policy))
                        .unwrap_or(Err(error)),
                    result => result.map(Some),
                };
            }
        }
        if self.uses_rust_geometry() && other.uses_rust_geometry() {
            if let Some(result) = self.glued_brep_union(other, operation, policy) {
                return result;
            }
            if let Some(result) = self.separated_treatment_intersection(other, operation, policy) {
                return result;
            }
            if let Some(result) = self.contained_treatment_box_boolean(other, operation, policy) {
                return result;
            }
            if let Some(result) = self.contained_thread_boolean(other, operation, policy) {
                return result;
            }
            if let Some(result) = self.contained_annulus_boolean(other, operation, policy) {
                return result;
            }
            if let Some(result) = self.contained_insert_boolean(other, operation, policy) {
                return result;
            }
            if let Some(result) = self.contained_hanger_boolean(other, operation, policy) {
                return result;
            }
            // Cheap whole-probe predicates of retained thread sources come
            // before the planar provider, which would integrate their faces.
            if let Some(result) = self.contained_radial_boundary_boolean(other, operation, policy) {
                return result;
            }
            if let Some(result) = self.contained_clipped_npt_boolean(other, operation, policy) {
                return result;
            }
            if let Some(result) = self.contained_warped_union_boolean(other, operation, policy) {
                return result;
            }
            if let Some(result) = self.default_certified_boolean(other, operation, policy) {
                return result;
            }
            if let Some(result) = self.contained_sweep_boolean(other, operation, policy) {
                return result;
            }
            if let Some(result) = self.subtract_cylinder_through_sweep(other, operation, policy) {
                return result;
            }
            if let Some(result) = self.concentric_ring_boolean(other, operation, policy) {
                return result;
            }
            if let Some(result) = self.rounded_thread_clearance_boolean(other, operation, policy) {
                return result;
            }
            if let Some(result) = self.contained_rounded_thread_boolean(other, operation, policy) {
                return result;
            }
            if let Some(result) = self.contained_gear_boolean(other, operation, policy) {
                return result;
            }
            if let Some(result) = self.gear_pair_boolean(other, operation, policy) {
                return result;
            }
            if let Some(result) = self.rounded_thread_axial_clip(other, operation, policy) {
                return result;
            }
            return match self.orthogonal_boolean(other, operation, policy) {
                Err(error) if error.kind == GeometryErrorKind::Unsupported => self
                    .separated_boolean(other, operation, policy)
                    .or_else(|| self.subtract_machining_tool(other, operation, policy))
                    .or_else(|| self.coaxial_counterbore_cut(other, operation, policy))
                    .or_else(|| self.through_profile_cut(other, operation, policy))
                    .or_else(|| self.bored_prism_subtract(other, operation, policy))
                    .or_else(|| self.stacked_prism_boolean(other, operation, policy))
                    .or_else(|| self.parallel_extrusion_intersection(other, operation, policy))
                    .or_else(|| self.stacked_cylinder_union(other, operation, policy))
                    .or_else(|| self.pairwise_warped_boolean(other, operation, policy))
                    .or_else(|| self.warped_shell_probe_boolean(other, operation, policy))
                    .or_else(|| self.contained_void_cut(other, operation, policy))
                    .unwrap_or_else(|| {
                        Err(cross_refusal
                            .unwrap_or_else(|| self.no_exact_boolean(other, operation, error)))
                    }),
                result => result,
            };
        }
        let op = match operation {
            BooleanOperation::Union => 0,
            BooleanOperation::Subtract => 1,
            BooleanOperation::Intersect => 2,
        };
        let handle = ffi::boolean_shapes(self.occt()?, other.occt()?, op, policy.boolean_fuzzy_mm)
            .map_err(|e| native("boolean", e))?;
        if handle.is_null() {
            return Ok(None);
        }
        let mut result = Self::adopt(handle, policy, "boolean")?;
        result.inherit_provenance(&[self, other], None)?;
        Ok(Some(result))
    }

    fn subtract_polygon_cylinder(
        target: &PlacedPolygon,
        tool: &AxialCylinder,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        let operation = "cylindrical cut";
        if target.frame != PlaneFrame::default() {
            return Err(error(
                operation,
                GeometryErrorKind::Unsupported,
                "Rust polygon through-cuts currently require an unplaced XY profile",
            ));
        }
        let zero = Rational::new(0, 1).expect("zero is an exact length");
        let height = target.prism.height();
        let (bottom, top) = if target.reverse_sweep {
            (
                zero.checked_sub(height)
                    .expect("negation fits the rational domain"),
                zero,
            )
        } else {
            (zero, height)
        };
        // These endpoints are exact rationals. Cap contact opens the bore,
        // but even a sub-tolerance shortfall leaves an uncut floor.
        if tool.bottom() > bottom || tool.top() < top {
            return Err(error(
                operation,
                GeometryErrorKind::Input,
                "The tool must reach or pass both polygon caps",
            ));
        }
        let result = Self::with_polygon_through_holes(
            target,
            [(
                tool.circle(),
                Origin {
                    source: EntityId::from_uuid_bits(tool.feature()),
                    role: "side.cylinder".into(),
                },
            )],
            policy,
        )?;
        if target.pockets.is_empty() {
            Ok(result)
        } else {
            result.with_polygon_pockets(target.pockets.to_vec(), target.pocket_origins.to_vec())
        }
    }
    /// Native COMMON query using OCCT's internal worker pool. Rust handles stay
    /// on their owner thread; the non-destructive operation owns its workers.
    /// Native validity, fuzzy tolerance and provenance rules match `boolean`.
    pub fn intersect_parallel(
        &self,
        other: &Self,
        policy: GeometryPolicy,
    ) -> Result<Option<Self>, GeometryError> {
        policy.validate()?;
        if self.uses_rust_geometry() || other.uses_rust_geometry() {
            // The exact Rust intersection has no parallel variant; its serial
            // result is the same certified material.
            return self.boolean(other, BooleanOperation::Intersect, policy);
        }
        let handle = ffi::intersect_parallel(self.occt()?, other.occt()?, policy.boolean_fuzzy_mm)
            .map_err(|e| native("parallel intersection", e))?;
        if handle.is_null() {
            return Ok(None);
        }
        let mut result = Self::adopt(handle, policy, "parallel intersection")?;
        result.inherit_provenance(&[self, other], None)?;
        Ok(Some(result))
    }
    pub fn transformed(&self, transform: RigidTransform) -> Result<Self, GeometryError> {
        transform.validate(self.policy)?;
        if self.uses_rust_geometry() {
            if transform == RigidTransform::default() {
                return self.share();
            }
            let body = match &self.body {
                SolidBody::Line8Nut(body) => {
                    SolidBody::Line8Nut(Rc::new(body.transformed(transform, self.policy)?))
                }
                SolidBody::ExternalNptStock(body) => {
                    SolidBody::ExternalNptStock(Rc::new(body.transformed(transform, self.policy)?))
                }
                SolidBody::LiteralNpt(body) => {
                    SolidBody::LiteralNpt(Rc::new(body.transformed(transform, self.policy)?))
                }
                SolidBody::ActualClippedHead(body) => {
                    SolidBody::ActualClippedHead(Rc::new(body.transformed(transform, self.policy)?))
                }
                SolidBody::ClippedNpt(body) => {
                    SolidBody::ClippedNpt(Rc::new(body.transformed(transform, self.policy)?))
                }
                SolidBody::PipeHardware(body) => {
                    SolidBody::PipeHardware(Rc::new(body.transformed(transform, self.policy)?))
                }
                SolidBody::Prism(prism, frame) => {
                    let frame = transform.frame(*frame);
                    frame.checked_metric_extent(
                        prism.maximum().map(Rational::as_f64),
                        self.policy,
                    )?;
                    SolidBody::Prism(Rc::clone(prism), frame)
                }
                SolidBody::PolygonPrism(placed) => {
                    let mut placed = placed.as_ref().clone();
                    placed.frame = transform.frame(placed.frame);
                    placed
                        .frame
                        .checked_metric_extent(polygon_metric_extent(&placed.prism), self.policy)?;
                    SolidBody::PolygonPrism(Rc::new(placed))
                }
                SolidBody::AxialCylinder(cylinder, frame, sketch) => {
                    let frame = transform.frame(*frame);
                    frame.checked_metric_extent(
                        [
                            cylinder.circle().radius().as_f64(),
                            cylinder.circle().radius().as_f64(),
                            cylinder.height().as_f64(),
                        ],
                        self.policy,
                    )?;
                    SolidBody::AxialCylinder(Rc::clone(cylinder), frame, sketch.clone())
                }
                SolidBody::Annulus(body, frame) => {
                    let frame = transform.frame(*frame);
                    frame.checked_metric_extent(
                        [
                            body.geometry.outer_radius_mm(),
                            body.geometry.outer_radius_mm(),
                            body.geometry.height_mm(),
                        ],
                        self.policy,
                    )?;
                    SolidBody::Annulus(Rc::clone(body), frame)
                }
                SolidBody::Hanger(body, frame) => {
                    let frame = transform.frame(*frame);
                    body.material(frame, self.policy)?;
                    SolidBody::Hanger(Rc::clone(body), frame)
                }
                SolidBody::WoodInsert(body, frame) => {
                    let frame = transform.frame(*frame);
                    body.material(frame, self.policy)?;
                    SolidBody::WoodInsert(Rc::clone(body), frame)
                }
                SolidBody::RepairedPhillips(body, frame) => {
                    let frame = transform.frame(*frame);
                    body.validate_frame(frame)?;
                    SolidBody::RepairedPhillips(Rc::clone(body), frame)
                }
                SolidBody::InternalThreadTreatment(body, frame) => {
                    let frame = transform.frame(*frame);
                    body.validate_frame(frame, self.policy)?;
                    SolidBody::InternalThreadTreatment(Rc::clone(body), frame)
                }
                SolidBody::RoundedBlindBore(body, frame) => {
                    let frame = transform.frame(*frame);
                    body.validate_frame(frame, self.policy)?;
                    SolidBody::RoundedBlindBore(Rc::clone(body), frame)
                }
                SolidBody::RoundedThreadBore(body, frame) => {
                    let frame = transform.frame(*frame);
                    frame.checked_metric_extent(
                        [
                            body.geometry.radius_mm(),
                            body.geometry.radius_mm(),
                            body.geometry.length_mm(),
                        ],
                        self.policy,
                    )?;
                    SolidBody::RoundedThreadBore(Rc::clone(body), frame)
                }
                SolidBody::RoundedThread(body, frame) => {
                    let frame = transform.frame(*frame);
                    let spec = body.geometry.spec();
                    frame.checked_metric_extent(
                        [spec.major_radius_mm, spec.major_radius_mm, spec.length_mm],
                        self.policy,
                    )?;
                    SolidBody::RoundedThread(Rc::clone(body), frame)
                }
                SolidBody::EntryChamfer(body, frame) => {
                    let frame = transform.frame(*frame);
                    let spec = body.geometry.source_thread().spec();
                    frame.checked_metric_extent(
                        [spec.major_radius_mm, spec.major_radius_mm, spec.length_mm],
                        self.policy,
                    )?;
                    SolidBody::EntryChamfer(Rc::clone(body), frame)
                }
                SolidBody::Helical(body, frame) => {
                    let frame = transform.frame(*frame);
                    frame.checked_metric_bounds(body.bounds(PlaneFrame::default()), self.policy)?;
                    SolidBody::Helical(Rc::clone(body), frame)
                }
                SolidBody::HelicalTubeCut(body, frame) => {
                    SolidBody::HelicalTubeCut(Rc::clone(body), transform.frame(*frame))
                }
                SolidBody::CosmeticFastener(body, frame) => {
                    SolidBody::CosmeticFastener(Rc::clone(body), transform.frame(*frame))
                }
                SolidBody::CosmeticProfiled(body, frame) => {
                    let frame = transform.frame(*frame);
                    frame.checked_metric_bounds(body.bounds(PlaneFrame::default()), self.policy)?;
                    SolidBody::CosmeticProfiled(Rc::clone(body), frame)
                }
                SolidBody::Linear(body, frame) => {
                    let frame = transform.frame(*frame);
                    frame.checked_metric_bounds(body.bounds(PlaneFrame::default()), self.policy)?;
                    SolidBody::Linear(Rc::clone(body), frame)
                }
                SolidBody::Revolved(body, frame) => {
                    SolidBody::Revolved(Rc::clone(body), transform.frame(*frame))
                }
                SolidBody::RoundedSlot(body, frame) => {
                    SolidBody::RoundedSlot(Rc::clone(body), transform.frame(*frame))
                }
                SolidBody::SplineSweep(body, frame) => {
                    SolidBody::SplineSweep(Rc::clone(body), transform.frame(*frame))
                }
                SolidBody::NormalSweep(body, frame) => {
                    let frame = transform.frame(*frame);
                    frame.checked_metric_bounds(body.bounds(PlaneFrame::default()), self.policy)?;
                    SolidBody::NormalSweep(Rc::clone(body), frame)
                }
                SolidBody::NormalRelief(body, frame) => {
                    let frame = transform.frame(*frame);
                    frame.checked_metric_bounds(body.bounds(PlaneFrame::default()), self.policy)?;
                    SolidBody::NormalRelief(Rc::clone(body), frame)
                }
                SolidBody::SplineCapRelief(body, frame) => {
                    let frame = transform.frame(*frame);
                    frame.checked_metric_bounds(body.bounds(PlaneFrame::default()), self.policy)?;
                    SolidBody::SplineCapRelief(Rc::clone(body), frame)
                }

                SolidBody::RadialBoundary(body, frame) => {
                    let frame = transform.frame(*frame);
                    frame.checked(self.policy)?;
                    let body = Rc::new(body.new_occurrence());
                    let bounds = body.bounds(frame, None)?;
                    self.policy.point(bounds.min)?;
                    self.policy.point(bounds.max)?;
                    SolidBody::RadialBoundary(body, frame)
                }
                SolidBody::SpurGear(body, frame) => {
                    SolidBody::SpurGear(Rc::clone(body), transform.frame(*frame))
                }
                SolidBody::InternalRing(body, frame) => {
                    SolidBody::InternalRing(Rc::clone(body), transform.frame(*frame))
                }
                SolidBody::Herringbone(body, frame) => {
                    SolidBody::Herringbone(Rc::clone(body), transform.frame(*frame))
                }
                SolidBody::Bevel(body, frame) => {
                    SolidBody::Bevel(Rc::clone(body), transform.frame(*frame))
                }
                SolidBody::Worm(body, frame) => {
                    SolidBody::Worm(Rc::clone(body), transform.frame(*frame))
                }
                SolidBody::Sphere(body, frame) => {
                    let frame = transform.frame(*frame);
                    frame.checked_metric_extent([body.geometry.radius_mm(); 3], self.policy)?;
                    SolidBody::Sphere(Rc::clone(body), frame)
                }
                SolidBody::ImportedBrep(body) => {
                    SolidBody::ImportedBrep(Rc::new(body.transformed(transform)?))
                }
                SolidBody::Torus(body, frame) => {
                    let frame = transform.frame(*frame);
                    let major = body.geometry.major_radius_mm();
                    let minor = body.geometry.minor_radius_mm();
                    frame.checked_metric_extent(
                        [major + minor, major + minor, minor],
                        self.policy,
                    )?;
                    SolidBody::Torus(Rc::clone(body), frame)
                }
                SolidBody::ObliqueBore(body, frame) => {
                    let frame = transform.frame(*frame);
                    frame.checked_metric_extent(body.geometry.size_mm(), self.policy)?;
                    SolidBody::ObliqueBore(Rc::clone(body), frame)
                }
                SolidBody::CrossBore(body, frame) => {
                    let frame = transform.frame(*frame);
                    frame.checked_metric_bounds(body.bounds(PlaneFrame::default()), self.policy)?;
                    SolidBody::CrossBore(Rc::clone(body), frame)
                }
                SolidBody::RoundedCylinder(body, frame) => {
                    let frame = transform.frame(*frame);
                    frame.checked_metric_bounds(body.bounds(PlaneFrame::default()), self.policy)?;
                    SolidBody::RoundedCylinder(Rc::clone(body), frame)
                }
                SolidBody::AxialStack(body, frame) => {
                    let frame = transform.frame(*frame);
                    frame.checked_metric_bounds(body.bounds(PlaneFrame::default()), self.policy)?;
                    SolidBody::AxialStack(Rc::clone(body), frame)
                }
                SolidBody::PlanarBoolean(body, frame) => {
                    SolidBody::PlanarBoolean(Rc::clone(body), transform.frame(*frame))
                }
                SolidBody::NonplanarLoft(body, frame) => {
                    SolidBody::NonplanarLoft(Rc::clone(body), transform.frame(*frame))
                }
                SolidBody::WarpedUnion(body, frame) => {
                    SolidBody::WarpedUnion(Rc::clone(body), transform.frame(*frame))
                }
                SolidBody::RustCompound(compound) => {
                    SolidBody::RustCompound(compound.transformed(transform, self.policy)?)
                }
                SolidBody::Orthogonal(body, frame) => {
                    SolidBody::Orthogonal(Rc::clone(body), transform.frame(*frame))
                }
                SolidBody::OrthogonalSeats(body, frame) => {
                    let frame = transform.frame(*frame);
                    body.validate_frame(frame, self.policy)?;
                    SolidBody::OrthogonalSeats(Rc::clone(body), frame)
                }
                SolidBody::Occt(_) => unreachable!("checked Rust body"),
            };
            let mut result = Self {
                body,
                policy: self.policy,
                provenance: self.transformed_author_metadata(transform),
                rust_origins: self.rust_origins.clone(),
                retained_placement: None,
            };
            let bounds = result.bounds()?;
            self.policy.point(bounds.min)?;
            self.policy.point(bounds.max)?;
            if let Some(geometry) = result.retained_geometry_key() {
                let mut source = self.share()?;
                // Keep only the most recent placement, not a recursive history.
                source.retained_placement = None;
                result.retained_placement = Some(Rc::new(RetainedPlacement {
                    source,
                    transform,
                    geometry,
                }));
            }
            return Ok(result);
        }
        let r: Vec<_> = transform.rotation.into_iter().flatten().collect();
        let mut result = Self::adopt(
            ffi::transform_shape(self.occt()?, &r, &point(transform.translation_mm))
                .map_err(|e| native("transform", e))?,
            self.policy,
            "transform",
        )?;
        result.inherit_provenance(&[self], Some(transform))?;
        Ok(result)
    }
    pub fn volume_mm3(&self) -> Result<f64, GeometryError> {
        if let Some(body) = self.pipe_hardware_body() {
            return body.volume_mm3();
        }
        if let Some((body, frame)) = self.repaired_phillips_body() {
            let [lo, hi] = body.volume(frame)?;
            return Ok(lo.midpoint(hi));
        }
        if let Some((body, frame)) = self.hanger_body() {
            return body.volume(frame, self.policy);
        }
        if let Some((body, frame)) = self.wood_insert_body() {
            let material = body.material(frame, self.policy)?;
            let b = material.source();
            let scale = b
                .outer()
                .maximum_radius_mm()
                .max(b.full_length_mm())
                .max(1.);
            let v = material.volume_bounds(
                multicad_rust_kernel::wood_insert::QueryOptions {
                    volume_error_mm3: (self.policy.linear_tolerance_mm * scale * scale).max(1e-7),
                    maximum_cells: 10,
                },
                None,
            )?;
            return Ok(v[0].midpoint(v[1]));
        }
        if let Some((body, frame)) = self.helical_body()
            && let Some(interval) = body.phillips_volume(frame)
        {
            let [lo, hi] = interval?;
            return Ok(lo.midpoint(hi));
        }
        if let Some((body, frame)) = self.internal_thread_treatment_body() {
            let interval = body.volume(frame)?;
            return Ok(interval.lower.midpoint(interval.upper));
        }
        if let Some((body, _)) = self.rounded_cylinder_body() {
            return Ok(body.geometry.volume_mm3());
        }
        if let Some((body, _)) = self.helical_tube_cut_body() {
            return Ok(body.cut.volume_mm3());
        }
        if let Some((body, _)) = self.torus_body() {
            return Ok(body.geometry.volume_mm3());
        }
        if let Some((body, _)) = self.oblique_bore_body() {
            return Ok(body.geometry.volume_mm3());
        }
        if let Some((body, _)) = self.cross_bore_body() {
            return Ok(body.geometry.volume_mm3());
        }
        if let Some((body, _)) = self.sphere_body() {
            return Ok(body.geometry.volume_mm3());
        }
        if let Some(body) = self.imported_brep_body() {
            return Ok(body.solid.volume_mm3());
        }
        if let Some((body, _)) = self.spur_gear_body() {
            return Ok(body.profile.volume_mm3());
        }
        if let Some((body, _)) = self.internal_ring_body() {
            return Ok(body.ring.volume_mm3());
        }
        if let Some((body, _)) = self.herringbone_body() {
            return Ok(body.geometry.volume_mm3());
        }
        if let Some((body, _)) = self.bevel_body() {
            return Ok(body.geometry.volume_mm3());
        }
        if let Some((body, _)) = self.worm_body() {
            return Ok(body.volume_mm3());
        }
        if let Some((body, _)) = self.warped_union_body() {
            return Ok(body.volume_mm3());
        }
        if let Some((body, _)) = self.planar_boolean_body() {
            return Ok(body.volume_mm3());
        }
        if let Some((body, _)) = self.nonplanar_loft_body() {
            return Ok(body.geometry.volume_mm3());
        }
        if let Some((body, _)) = self.axial_stack_body() {
            return Ok(body.volume_mm3());
        }
        if let Some((body, _)) = self.spline_cap_relief_body() {
            return Ok(body.volume_mm3());
        }
        if let Some((body, _)) = self.normal_relief_body() {
            return Ok(body.volume_mm3());
        }
        if let Some((body, frame)) = self.radial_boundary_body() {
            return body.volume_mm3(frame);
        }
        if let Some((body, _)) = self.normal_sweep_body() {
            return Ok(body.volume_mm3());
        }
        if let Some((body, _)) = self.spline_sweep_body() {
            return Ok(body.sweep.volume_mm3());
        }
        if let Some(body) = self.line8_nut_body() {
            return Ok(body.volume_mm3());
        }
        if let Some((body, frame)) = self.orthogonal_seat_body() {
            return body.volume_mm3(frame);
        }
        if let Some(body) = self.external_npt_stock_body() {
            return Ok(body.volume_mm3());
        }
        if let Some(body) = self.literal_npt_body() {
            return Ok(body.volume_mm3());
        }
        if let Some(body) = self.actual_clipped_head_body() {
            return Ok(body.volume_mm3());
        }
        if let Some(body) = self.clipped_npt_body() {
            return Ok(body.volume_mm3());
        }
        if let Some((body, _)) = self.cosmetic_profiled_body_ref() {
            return Ok(body.volume_mm3());
        }
        if let Some((body, _)) = self.cosmetic_fastener_body() {
            return Ok(body.geometry.volume_mm3());
        }
        if let Some((body, _)) = self.rounded_slot_body() {
            return Ok(body.volume_mm3());
        }
        if let Some((body, _)) = self.linear_body() {
            return Ok(body.volume_mm3());
        }
        if let Some((body, _)) = self.revolution_body() {
            return Ok(body.geometry.volume_mm3());
        }
        if let Some((body, _)) = self.annulus_body() {
            return Ok(body.geometry.volume_mm3());
        }
        if let Some((body, frame)) = self.rounded_blind_bore_body() {
            return Ok(body.volume(frame));
        }
        if let Some((body, _)) = self.rounded_thread_bore_body() {
            return Ok(body.geometry.volume_mm3());
        }
        if let Some((body, _)) = self.rounded_thread_body() {
            return Ok(body.geometry.volume_mm3());
        }
        if let Some((body, _)) = self.entry_chamfer_body() {
            return Ok(body.volume());
        }
        if let Some((body, _)) = self.helical_body() {
            return Ok(body.volume_mm3());
        }
        if let Some((body, _)) = self.orthogonal_body() {
            return Ok(body.volume_mm3());
        }
        if let Some(compound) = self.rust_compound() {
            return Ok(compound.volume_mm3());
        }
        if let Some(placed) = self.polygon_prism() {
            if let Some(treatments) = &placed.bore_treatments {
                return Ok(treatments.geometry.volume_mm3());
            }
            return Ok(multicad_rust_kernel::circular_pockets::volume_mm3(
                &placed.prism,
                &placed.pockets,
            ));
        }
        if let Some((cylinder, _)) = self.axial_cylinder_frame() {
            return Ok(cylinder.volume_mm3());
        }
        if let Some(prism) = self.rust_prism() {
            return Ok(prism.volume().as_f64());
        }
        ffi::volume(self.occt()?).map_err(|e| native("volume", e))
    }
    pub fn solid_count(&self) -> Result<u32, GeometryError> {
        if let Some((body, _)) = self.warped_union_body() {
            return Ok(body.solid_count());
        }
        if let Some((body, _)) = self.planar_boolean_body() {
            return Ok(body.solid_count());
        }
        if let Some((body, _)) = self.orthogonal_body() {
            return Ok(body.geometry.solid_count());
        }
        if let Some(compound) = self.rust_compound() {
            return Ok(compound.solid_count());
        }
        if self.uses_rust_geometry() {
            return Ok(1);
        }
        ffi::solid_count(self.occt()?).map_err(|e| native("solid count", e))
    }
    pub fn bounds(&self) -> Result<Bounds, GeometryError> {
        if let Some(body) = self.pipe_hardware_body() {
            return body.bounds();
        }
        if let Some((body, frame)) = self.repaired_phillips_body() {
            return body.bounds(frame);
        }
        if let Some((body, frame)) = self.hanger_body() {
            return Ok(body.material(frame, self.policy)?.bounds());
        }
        if let Some((body, frame)) = self.wood_insert_body() {
            return body.material(frame, self.policy)?.bounds();
        }
        if let Some((body, frame)) = self.internal_thread_treatment_body() {
            return Ok(body.bounds(frame));
        }
        if let Some((body, frame)) = self.rounded_cylinder_body() {
            return Ok(body.bounds(frame));
        }
        if let Some((body, frame)) = self.helical_tube_cut_body() {
            return body.bounds(frame);
        }
        if let Some((body, frame)) = self.torus_body() {
            return Ok(body.bounds(frame));
        }
        if let Some((body, frame)) = self.oblique_bore_body() {
            return Ok(body.bounds(frame));
        }
        if let Some((body, frame)) = self.cross_bore_body() {
            return Ok(body.bounds(frame));
        }
        if let Some((body, frame)) = self.sphere_body() {
            return Ok(body.bounds(frame));
        }
        if let Some(body) = self.imported_brep_body() {
            return Ok(body.bounds());
        }
        if let Some((body, frame)) = self.spur_gear_body() {
            return Ok(body.bounds(frame));
        }
        if let Some((body, frame)) = self.internal_ring_body() {
            return Ok(body.bounds(frame));
        }
        if let Some((body, frame)) = self.herringbone_body() {
            return Ok(body.bounds(frame));
        }
        if let Some((body, frame)) = self.bevel_body() {
            return Ok(body.bounds(frame));
        }
        if let Some((body, frame)) = self.worm_body() {
            return Ok(body.bounds(frame));
        }
        if let Some((body, frame)) = self.warped_union_body() {
            return Ok(body.bounds(frame));
        }
        if let Some((body, frame)) = self.planar_boolean_body() {
            return Ok(body.bounds(frame));
        }
        if let Some((body, frame)) = self.nonplanar_loft_body() {
            return Ok(body.bounds(frame));
        }
        if let Some((body, frame)) = self.axial_stack_body() {
            return Ok(body.bounds(frame));
        }
        if let Some((body, frame)) = self.spline_cap_relief_body() {
            return Ok(body.bounds(frame));
        }
        if let Some((body, frame)) = self.normal_relief_body() {
            return Ok(body.bounds(frame));
        }
        if let Some((body, frame)) = self.radial_boundary_body() {
            return body.bounds(frame, None);
        }
        if let Some((body, frame)) = self.normal_sweep_body() {
            return Ok(body.bounds(frame));
        }
        if let Some((body, frame)) = self.spline_sweep_body() {
            return Ok(body.bounds(frame));
        }
        if let Some(body) = self.line8_nut_body() {
            return Ok(body.bounds());
        }
        if let Some((body, frame)) = self.orthogonal_seat_body() {
            return body.bounds(frame);
        }
        if let Some(body) = self.external_npt_stock_body() {
            return Ok(body.bounds());
        }
        if let Some(body) = self.literal_npt_body() {
            return Ok(body.bounds());
        }
        if let Some(body) = self.actual_clipped_head_body() {
            return body.bounds();
        }
        if let Some(body) = self.clipped_npt_body() {
            return Ok(body.bounds());
        }
        if let Some((body, frame)) = self.cosmetic_profiled_body_ref() {
            return Ok(body.bounds(frame));
        }
        if let Some((body, frame)) = self.cosmetic_fastener_body() {
            return Ok(body.bounds(frame));
        }
        if let Some((body, frame)) = self.rounded_slot_body() {
            return Ok(body.bounds(frame));
        }
        if let Some((body, frame)) = self.linear_body() {
            return Ok(body.bounds(frame));
        }
        if let Some((body, frame)) = self.revolution_body() {
            return Ok(body.bounds(frame));
        }
        if let Some((body, frame)) = self.annulus_body() {
            return Ok(body.bounds(frame));
        }
        if let Some((body, frame)) = self.rounded_blind_bore_body() {
            return body.bounds(frame);
        }
        if let Some((body, frame)) = self.rounded_thread_bore_body() {
            return body.bounds(frame);
        }
        if let Some((body, frame)) = self.rounded_thread_body() {
            return body.bounds(frame);
        }
        if let Some((body, frame)) = self.entry_chamfer_body() {
            return body.bounds(frame);
        }
        if let Some((body, frame)) = self.helical_body() {
            return Ok(body.bounds(frame));
        }
        if let Some((body, frame)) = self.orthogonal_body() {
            return Ok(body.bounds(frame));
        }
        if let Some(compound) = self.rust_compound() {
            return Ok(compound.bounds());
        }
        if let Some(placed) = self.polygon_prism() {
            return Ok(placed.bounds());
        }
        if let Some(cylinder) = self.axial_cylinder() {
            let [min, max] = cylinder.bounds().map(|corner| corner.map(Rational::as_f64));
            return Ok(Bounds { min, max });
        }
        if let Some((cylinder, frame)) = self.axial_cylinder_frame() {
            let radius = cylinder.circle().radius().as_f64();
            let height = cylinder.height().as_f64();
            let axis = frame.normal();
            let extent: [f64; 3] = std::array::from_fn(|i| radius * frame.u[i].hypot(frame.v[i]));
            return Ok(Bounds {
                min: std::array::from_fn(|i| {
                    frame.origin[i] + (height * axis[i]).min(0.) - extent[i]
                }),
                max: std::array::from_fn(|i| {
                    frame.origin[i] + (height * axis[i]).max(0.) + extent[i]
                }),
            });
        }
        if let Some((prism, frame)) = self.rust_prism_frame() {
            if frame == PlaneFrame::default() {
                return Ok(Bounds {
                    min: [0.; 3],
                    max: prism.maximum().map(Rational::as_f64),
                });
            }
            let [length, width, height] = prism.dimensions().map(Rational::as_f64);
            let axis = frame.normal();
            let support = |u: f64, v: f64, n: f64| {
                let rectangle = (u * length).max(0.) + (v * width).max(0.);
                let cap = prism.rounded_end().map_or(rectangle, |circle| {
                    let radius = circle.radius().as_f64();
                    let arc = u * length + v * radius + radius * u.max(0.).hypot(v);
                    rectangle.max(arc)
                });
                cap + (n * height).max(0.)
            };
            return Ok(Bounds {
                min: std::array::from_fn(|i| {
                    frame.origin[i] - support(-frame.u[i], -frame.v[i], -axis[i])
                }),
                max: std::array::from_fn(|i| {
                    frame.origin[i] + support(frame.u[i], frame.v[i], axis[i])
                }),
            });
        }
        let b = ffi::bounds(self.occt()?).map_err(|e| native("bounds", e))?;
        Ok(Bounds {
            min: xyz(b.min),
            max: xyz(b.max),
        })
    }
    /// Measures the actual selected finite subshapes. Native handles remain on
    /// their owner thread; no display tessellation participates.
    pub fn minimum_distance(
        &self,
        a: crate::BRepSelection,
        other: &Self,
        b: crate::BRepSelection,
    ) -> Result<crate::MinimumDistance, GeometryError> {
        crate::distance::measure(
            self,
            a,
            other,
            b,
            self.policy
                .linear_tolerance_mm
                .max(other.policy.linear_tolerance_mm),
        )
    }
    /// Extract ordered exact wires on the explicitly chosen finite planar face.
    /// The frame must be coplanar; unsupported boundaries never become polylines.
    pub fn planar_face_contours(
        &self,
        face: u32,
        frame: PlaneFrame,
    ) -> Result<multicad_core::manufacturing::PlanarContours, GeometryError> {
        crate::planar_contours::extract(self, face, frame, self.policy)
    }
    pub fn faces(&self) -> Result<Vec<Face>, GeometryError> {
        if let Some(body) = self.pipe_hardware_body() {
            return body.faces(self.policy);
        }
        if let Some((body, frame)) = self.repaired_phillips_body() {
            return if body.has_degree11_source() {
                body.degree11_faces(frame)
            } else {
                body.periodic_faces(frame)
            };
        }
        if let Some((body, frame)) = self.hanger_body() {
            return body.material(frame, self.policy)?.faces(self.policy, None);
        }
        if let Some((body, frame)) = self.wood_insert_body() {
            return body.material(frame, self.policy)?.faces(self.policy, None);
        }
        if let Some((body, frame)) = self.helical_body()
            && let Some(faces) = body.phillips_faces(frame, self.policy)
        {
            return faces;
        }
        if let Some((body, frame)) = self.internal_thread_treatment_body() {
            return body.faces(frame, self.policy);
        }
        if let Some((body, frame)) = self.rounded_cylinder_body() {
            return Ok(body.faces(frame));
        }
        if let Some((body, frame)) = self.helical_tube_cut_body() {
            return Ok(body.faces(frame));
        }
        if let Some((body, frame)) = self.torus_body() {
            return Ok(vec![body.face(frame)]);
        }
        if let Some((body, frame)) = self.oblique_bore_body() {
            return Ok(body.faces(frame));
        }
        if let Some((body, frame)) = self.cross_bore_body() {
            return Ok(body.faces(frame));
        }
        if let Some((body, frame)) = self.sphere_body() {
            return Ok(vec![body.face(frame)]);
        }
        if let Some(body) = self.imported_brep_body() {
            return Ok(body.faces());
        }
        if let Some((body, frame)) = self.spur_gear_body() {
            return Ok(body.faces(frame));
        }
        if let Some((body, frame)) = self.internal_ring_body() {
            return Ok(body.faces(frame));
        }
        if let Some((body, frame)) = self.herringbone_body() {
            return Ok(body.faces(frame));
        }
        if let Some((body, frame)) = self.bevel_body() {
            return Ok(body.faces(frame));
        }
        if let Some((body, frame)) = self.worm_body() {
            return Ok(body.faces(frame));
        }
        if let Some((body, frame)) = self.warped_union_body() {
            return Ok(body.faces(frame));
        }
        if let Some((body, frame)) = self.planar_boolean_body() {
            return Ok(body.faces(frame));
        }
        if let Some((body, frame)) = self.nonplanar_loft_body() {
            return body.faces(frame);
        }
        if let Some((body, frame)) = self.axial_stack_body() {
            return Ok(body.faces(frame));
        }
        if let Some((body, frame)) = self.spline_cap_relief_body() {
            return Ok(body.faces(frame));
        }
        if let Some((body, frame)) = self.normal_relief_body() {
            return Ok(body.faces(frame));
        }
        if let Some((body, frame)) = self.radial_boundary_body() {
            return body.faces(frame);
        }
        if let Some((body, frame)) = self.normal_sweep_body() {
            return Ok(body.faces(frame));
        }
        if let Some((body, frame)) = self.spline_sweep_body() {
            return Ok(body.faces(frame));
        }
        if let Some(body) = self.line8_nut_body() {
            return Ok(body.faces());
        }
        if let Some((body, frame)) = self.orthogonal_seat_body() {
            return body.faces(frame);
        }
        if let Some(body) = self.external_npt_stock_body() {
            return Ok(body.faces());
        }
        if let Some(body) = self.literal_npt_body() {
            return body.faces();
        }
        if let Some(body) = self.actual_clipped_head_body() {
            return body.faces();
        }
        if let Some(body) = self.clipped_npt_body() {
            return Ok(body.faces());
        }
        if let Some((body, frame)) = self.cosmetic_profiled_body_ref() {
            return Ok(body.faces(frame));
        }
        if let Some((body, frame)) = self.cosmetic_fastener_body() {
            return Ok(body.faces(frame));
        }
        if let Some((body, frame)) = self.rounded_slot_body() {
            return Ok(body.faces(frame));
        }
        if let Some((body, frame)) = self.linear_body() {
            return Ok(body.faces(frame));
        }
        if let Some((body, frame)) = self.revolution_body() {
            return Ok(body.faces(frame));
        }
        if let Some((body, frame)) = self.annulus_body() {
            return Ok(body.faces(frame));
        }
        if let Some((body, frame)) = self.rounded_blind_bore_body() {
            return body.faces(frame);
        }
        if let Some((body, frame)) = self.rounded_thread_bore_body() {
            return body.faces(frame);
        }
        if let Some((body, frame)) = self.rounded_thread_body() {
            return body.faces(frame);
        }
        if let Some((body, frame)) = self.entry_chamfer_body() {
            return Ok(body.faces(frame));
        }
        if let Some((body, frame)) = self.helical_body() {
            return Ok(body.faces(frame));
        }
        if let Some((body, frame)) = self.orthogonal_body() {
            return Ok(body.faces(frame));
        }
        if let Some(compound) = self.rust_compound() {
            return compound.faces();
        }
        if let Some(placed) = self.polygon_prism() {
            return Ok(placed.faces());
        }
        if let Some((cylinder, frame)) = self.axial_cylinder_frame() {
            return Ok(crate::rust_prism_topology::cylinder_faces(cylinder, frame));
        }
        if let Some((prism, frame)) = self.rust_prism_frame() {
            return Ok(crate::rust_prism_topology::faces(prism, frame));
        }
        // An exact affine image of a proven plane stays planar even when OCCT
        // stores it as a NURBS surface. Preserve that typed proof for picking;
        // no surface fitting or display triangles establish this descriptor.
        let planes: std::collections::BTreeMap<_, _> = self
            .provenance
            .iter()
            .filter_map(|e| match e.geometry {
                Some(multicad_core::references::AnalyticGeometry::Plane { frame })
                    if e.key.kind == multicad_core::references::TopologyKind::Face =>
                {
                    Some((e.key.ordinal, (frame.origin_mm, frame.normal())))
                }
                _ => None,
            })
            .collect();
        ffi::faces(self.occt()?)
            .map_err(|e| native("faces", e))
            .map(|faces| {
                faces
                    .into_iter()
                    .map(|f| Face {
                        index: f.index,
                        plane: f
                            .planar
                            .then_some((xyz(f.origin), xyz(f.normal)))
                            .or_else(|| planes.get(&f.index).copied()),
                        center_mm: xyz(f.center),
                        area_mm2: f.area,
                    })
                    .collect()
            })
    }
    fn mesh_options_valid(&self, options: MeshOptions) -> bool {
        (self.policy.linear_tolerance_mm.max(0.001)..=self.policy.maximum_coordinate_mm)
            .contains(&options.linear_deflection_mm)
            && (0.01..=1.).contains(&options.angular_deflection_rad)
            && (1..=5_000_000).contains(&options.maximum_triangles)
    }

    fn rust_mesh_options(options: MeshOptions) -> multicad_rust_kernel::display::Tessellation {
        multicad_rust_kernel::display::Tessellation {
            deflection: options.linear_deflection_mm,
            angle: options.angular_deflection_rad,
            triangle_limit: options.maximum_triangles as usize,
        }
    }

    fn rust_mesh(&self, options: MeshOptions, usage: MeshUse) -> Result<Mesh, GeometryError> {
        let rust_options = Self::rust_mesh_options(options);
        match &self.body {
            SolidBody::Prism(prism, frame) => crate::rust_display::placed_mesh(
                multicad_rust_kernel::display::rust_prism(prism, rust_options),
                *frame,
                false,
            ),
            SolidBody::PolygonPrism(placed) => {
                let mesh = if let Some(treatments) = &placed.bore_treatments {
                    multicad_rust_kernel::display::chamfered_polygon_prism(
                        &placed.prism,
                        &treatments.geometry,
                        rust_options,
                    )
                } else if placed.pockets.is_empty() {
                    multicad_rust_kernel::display::polygon_prism(&placed.prism, rust_options)
                } else {
                    multicad_rust_kernel::display::pocketed_polygon_prism(
                        &placed.prism,
                        &placed.pockets,
                        rust_options,
                    )
                };
                crate::rust_display::placed_mesh(mesh, placed.frame, placed.reverse_sweep)
            }
            SolidBody::AxialCylinder(cylinder, frame, _) => {
                let mut mesh = crate::rust_display::mesh(
                    multicad_rust_kernel::display::axial_cylinder(cylinder, rust_options),
                )?;
                let [x, y] = cylinder.circle().center().map(Rational::as_f64);
                let origin = [x, y, cylinder.bottom().as_f64()];
                if self.axial_cylinder().is_none() {
                    for point in mesh.vertices.as_chunks_mut::<3>().0 {
                        let local = std::array::from_fn(|i| point[i] - origin[i]);
                        *point = crate::rust_prism_topology::position(*frame, local);
                    }
                }
                Ok(mesh)
            }
            SolidBody::Annulus(body, frame) => crate::rust_display::placed_mesh(
                multicad_rust_kernel::annulus::tessellate(&body.geometry, rust_options),
                *frame,
                false,
            ),
            SolidBody::EntryChamfer(body, frame) => body.mesh(*frame, rust_options),
            SolidBody::InternalThreadTreatment(body, frame) => body.mesh(*frame, rust_options),
            SolidBody::RepairedPhillips(body, frame) => body.mesh(*frame, rust_options, usage),
            SolidBody::Hanger(body, frame) => {
                body.material(*frame, self.policy)?.mesh(options, None)
            }
            SolidBody::WoodInsert(body, frame) => body.material(*frame, self.policy)?.mesh(
                MeshOptions {
                    linear_deflection_mm: rust_options.deflection,
                    angular_deflection_rad: rust_options.angle,
                    maximum_triangles: rust_options.triangle_limit as u32,
                },
                None,
            ),
            SolidBody::Helical(body, frame) => body.mesh(*frame, rust_options),
            SolidBody::RoundedThreadBore(body, frame) => body.mesh(*frame, rust_options),
            SolidBody::RoundedBlindBore(body, frame) => body.mesh(*frame, rust_options),
            SolidBody::RoundedThread(body, frame) => body.mesh(*frame, rust_options),
            SolidBody::HelicalTubeCut(body, frame) => {
                crate::rust_display::placed_mesh(body.cut.mesh(rust_options), *frame, false)
            }
            SolidBody::ClippedNpt(body) => body.mesh(rust_options),
            SolidBody::PipeHardware(body) => body.mesh(rust_options),
            SolidBody::LiteralNpt(body) => body.mesh(rust_options),
            SolidBody::ActualClippedHead(body) => body.mesh(rust_options),
            SolidBody::ExternalNptStock(body) => body.mesh(rust_options),
            SolidBody::Line8Nut(body) => body.mesh(rust_options),
            SolidBody::OrthogonalSeats(body, frame) => body.mesh(*frame, rust_options),
            SolidBody::CosmeticProfiled(body, frame) => match usage {
                MeshUse::Display => body.mesh(*frame, rust_options),
                MeshUse::Manufacturing => Err(error(
                    "manufacturing mesh",
                    GeometryErrorKind::Unsupported,
                    "Cosmetic thread envelopes require modeled thread geometry for manufacturing",
                )),
            },
            SolidBody::CosmeticFastener(body, frame) => match usage {
                MeshUse::Display => body.mesh(*frame, rust_options),
                MeshUse::Manufacturing => Err(error(
                    "manufacturing mesh",
                    GeometryErrorKind::Unsupported,
                    "Cosmetic thread envelopes require modeled thread geometry for manufacturing",
                )),
            },
            SolidBody::Linear(body, frame) => {
                crate::rust_display::placed_mesh(body.mesh(rust_options), *frame, false)
            }
            SolidBody::Revolved(body, frame) => body.mesh(*frame, rust_options),
            SolidBody::RoundedCylinder(body, frame) => body.mesh(*frame, rust_options),
            SolidBody::RoundedSlot(body, frame) => body.mesh(*frame, rust_options),
            SolidBody::NormalSweep(body, frame) => body.mesh(*frame, rust_options),
            SolidBody::NormalRelief(body, frame) => body.mesh(*frame, rust_options),
            SolidBody::SplineCapRelief(body, frame) => body.mesh(*frame, rust_options),
            SolidBody::RadialBoundary(body, frame) => body.mesh(*frame, rust_options, usage),
            SolidBody::SplineSweep(body, frame) => body.mesh(*frame, rust_options),
            SolidBody::SpurGear(body, frame) => {
                crate::rust_display::placed_mesh(body.profile.mesh(rust_options), *frame, false)
            }
            SolidBody::InternalRing(body, frame) => {
                crate::rust_display::placed_mesh(body.ring.mesh(rust_options), *frame, false)
            }
            SolidBody::Herringbone(body, frame) => {
                crate::rust_display::placed_mesh(body.geometry.mesh(rust_options), *frame, false)
            }
            SolidBody::Bevel(body, frame) => body.mesh(*frame, rust_options),
            SolidBody::Worm(body, frame) => body.mesh(*frame, rust_options),
            SolidBody::ImportedBrep(body) => body.mesh(rust_options),
            SolidBody::Sphere(body, frame) => {
                crate::rust_display::placed_mesh(body.geometry.mesh(rust_options), *frame, false)
            }
            SolidBody::Torus(body, frame) => {
                crate::rust_display::placed_mesh(body.geometry.mesh(rust_options), *frame, false)
            }
            SolidBody::ObliqueBore(body, frame) => {
                crate::rust_display::placed_mesh(body.geometry.mesh(rust_options), *frame, false)
            }
            SolidBody::CrossBore(body, frame) => body.mesh(*frame, rust_options),
            SolidBody::AxialStack(body, frame) => body.mesh(*frame, rust_options),
            SolidBody::PlanarBoolean(body, frame) => body.mesh(*frame, rust_options),
            SolidBody::NonplanarLoft(body, frame) => {
                crate::rust_display::placed_mesh(body.geometry.mesh(rust_options), *frame, false)
            }
            SolidBody::WarpedUnion(body, frame) => body.mesh(*frame, rust_options),
            SolidBody::RustCompound(compound) => compound.mesh(options, usage),
            SolidBody::Orthogonal(body, frame) => body.mesh(*frame, rust_options, usage),
            SolidBody::Occt(_) => unreachable!("Rust-only tessellation"),
        }
    }

    pub fn mesh(&self, options: MeshOptions) -> Result<Mesh, GeometryError> {
        if !self.mesh_options_valid(options) {
            return Err(input("Invalid display tessellation settings/budget"));
        }
        let validate = |vertices: &[f64], indices: &[u32], faces: &[u32]| {
            valid_parts(vertices, indices, faces)
                .then_some(())
                .ok_or_else(|| {
                    error(
                        "mesh",
                        GeometryErrorKind::InvalidResult,
                        "Invalid display mesh",
                    )
                })
        };
        let mesh = match &self.body {
            SolidBody::Occt(shape) => {
                let mesh = ffi::tessellate(
                    &shape.0,
                    options.linear_deflection_mm,
                    options.angular_deflection_rad,
                    options.maximum_triangles,
                )
                .map_err(|e| native("mesh", e))?;
                validate(&mesh.vertices, &mesh.indices, &mesh.face_indices)?;
                // Retain only after validating all detached arrays. Optional
                // receipt allocation failure does not invalidate the mesh.
                let _ = ffi::remember_validated_mesh(
                    &shape.0,
                    options.linear_deflection_mm,
                    options.angular_deflection_rad,
                    options.maximum_triangles,
                    &mesh,
                );
                return Ok(mesh.into());
            }
            _ => self.rust_mesh(options, MeshUse::Display)?,
        };
        validate(&mesh.vertices, &mesh.indices, &mesh.face_indices)?;
        Ok(mesh)
    }
    /// Fresh manufacturing triangulation with exact boundary identity.
    pub fn manufacturing_mesh(&self, options: MeshOptions) -> Result<Mesh, GeometryError> {
        Ok(self.manufacturing_mesh_with_receipt(options)?.mesh)
    }

    pub fn has_retained_radial_source(&self) -> bool {
        self.radial_boundary_body().is_some()
            || self
                .rust_compound()
                .is_some_and(|compound| compound.members().any(Self::has_retained_radial_source))
    }

    /// These sources prove actual facets and literal embedding independently of
    /// an analytic scalar. Keep analytic children separately checked in mixed
    /// compounds; a source interval is not a mesh-volume equality receipt.
    fn has_certified_source_mesh(&self) -> bool {
        self.helical_body()
            .is_some_and(|(body, _)| body.has_phillips_source())
            || self.has_retained_radial_source()
            || self.internal_thread_treatment_body().is_some()
            || self.repaired_phillips_body().is_some()
            || self.clipped_npt_body().is_some()
            || self.pipe_hardware_body().is_some()
            || self.literal_npt_body().is_some()
            || self.actual_clipped_head_body().is_some()
            || self.external_npt_stock_body().is_some()
            || self.line8_nut_body().is_some()
            || self.wood_insert_body().is_some()
            || self.hanger_body().is_some()
            || self
                .rust_compound()
                .is_some_and(|compound| compound.members().any(Self::has_certified_source_mesh))
    }
    fn rust_manufacturing_mesh(
        &self,
        options: MeshOptions,
    ) -> Result<crate::ManufacturingMesh, GeometryError> {
        use multicad_rust_kernel::{
            display::{DisplayMesh, MeshError},
            manufacturing,
        };
        const OPERATION: &str = "manufacturing mesh";
        if let Some((body, frame)) = self.rounded_blind_bore_body() {
            return body.manufacturing(frame, options);
        }
        if let Some(body) = self.pipe_hardware_body() {
            return body.manufacturing(Self::rust_mesh_options(options));
        }
        if let Some(body) = self.line8_nut_body() {
            return body.manufacturing(Self::rust_mesh_options(options));
        }
        if let Some(body) = self.external_npt_stock_body() {
            return body.manufacturing(Self::rust_mesh_options(options));
        }
        if let Some(body) = self.literal_npt_body() {
            return body.manufacturing(Self::rust_mesh_options(options));
        }
        if let Some(body) = self.actual_clipped_head_body() {
            return body.manufacturing(Self::rust_mesh_options(options));
        }
        if let Some(body) = self.clipped_npt_body() {
            return body.manufacturing(Self::rust_mesh_options(options));
        }
        if let Some((body, frame)) = self.repaired_phillips_body() {
            // This provider carries its fresh actual-world manifold proof and
            // consumes one source/placement/embedding work ledger. Do not run a
            // second unmetered certification or an analytic-volume heuristic.
            return body.manufacturing(frame, options);
        }
        let Mesh {
            vertices,
            indices,
            face_indices,
        } = self.rust_mesh(options, MeshUse::Manufacturing)?;
        let tessellation = DisplayMesh {
            vertices,
            indices,
            face_indices,
        };
        let proof = manufacturing::certify(&tessellation).map_err(|cause| {
            let kind = match cause {
                MeshError::InvalidOptions => GeometryErrorKind::Input,
                MeshError::BudgetExceeded => GeometryErrorKind::Limit,
                MeshError::UnresolvedBoundary => GeometryErrorKind::InvalidResult,
            };
            error(
                OPERATION,
                kind,
                format!("Rust indexed mesh certification failed: {cause:?}"),
            )
        })?;
        // Retained radial providers have already proved every actual facet
        // against the complete finite source and certified literal embedding.
        // Their independent mass enclosure is not an exact analytic scalar.
        // Mixed compounds verify analytic children before merging instead.
        if !self.has_certified_source_mesh() {
            let exact_volume = self.volume_mm3()?;
            let volume_tolerance = self.rust_mesh_volume_tolerance(options, exact_volume)?;
            if !volume_tolerance.is_finite()
                || (proof.signed_volume_mm3 - exact_volume).abs() > volume_tolerance
            {
                return Err(error(
                    OPERATION,
                    GeometryErrorKind::InvalidResult,
                    "Rust manufacturing mesh volume differs from the exact solid beyond the chord error bound",
                ));
            }
        }
        let DisplayMesh {
            vertices,
            indices,
            face_indices,
        } = tessellation;
        Ok(crate::ManufacturingMesh {
            mesh: Mesh {
                vertices,
                indices,
                face_indices,
            },
            maximum_displacement_mm: 0.,
            edge_pairs: proof.edge_pairs,
            vertex_links: proof.vertex_links,
            removed_native_pole_triangles: 0,
            removed_exact_zero_triangles: 0,
            exact_vertex_pairs: 0,
            internal_constraint_edges: 0,
            removed_meshing_internal_edges: 0,
            removed_meshing_constraint_wire_containers: 0,
            removed_meshing_internal_vertices: 0,
            remeshed_without_internal_constraints: false,
        })
    }

    fn rust_mesh_volume_tolerance(
        &self,
        options: MeshOptions,
        exact_volume: f64,
    ) -> Result<f64, GeometryError> {
        if let Some(compound) = self.rust_compound() {
            return compound.members().try_fold(0., |total, member| {
                Ok(total + member.rust_mesh_volume_tolerance(options, member.volume_mm3()?)?)
            });
        }
        let mut surface_area = 0.;
        let mut curved_area = 0.;
        for face in self.faces()? {
            surface_area += face.area_mm2;
            if face.plane.is_none() {
                curved_area += face.area_mm2;
            }
        }
        let chord_error = if let Some((body, _)) = self.rounded_thread_bore_body() {
            body.volume_error_bound(Self::rust_mesh_options(options))?
        } else if let Some((body, _)) = self.rounded_thread_body() {
            body.volume_error_bound(Self::rust_mesh_options(options))?
        } else if let Some((body, _)) = self.helical_body() {
            body.volume_error_bound(Self::rust_mesh_options(options))
        } else if let Some((body, _)) = self.cross_bore_body() {
            body.volume_error_bound(Self::rust_mesh_options(options))
        } else if let Some((body, _)) = self.revolution_body() {
            multicad_rust_kernel::revolution_mesh::volume_error_bound(
                &body.geometry,
                Self::rust_mesh_options(options),
            )
        } else if let Some((body, _)) = self.annulus_body() {
            multicad_rust_kernel::annulus::volume_error_bound(
                &body.geometry,
                Self::rust_mesh_options(options),
            )
        } else {
            // Arc/chord area loss is bounded by sagitta times arc length.
            // Extruding this loss bounds volume by the curved wall area.
            options.linear_deflection_mm * curved_area
        };
        Ok(
            chord_error
                + self.policy.linear_tolerance_mm * surface_area * 2.
                + exact_volume * 1e-10,
        )
    }

    pub fn manufacturing_mesh_with_receipt(
        &self,
        options: MeshOptions,
    ) -> Result<crate::ManufacturingMesh, GeometryError> {
        if !self.mesh_options_valid(options) {
            return Err(input("Invalid manufacturing tessellation settings/budget"));
        }
        if self.uses_rust_geometry() {
            return self.rust_manufacturing_mesh(options);
        }
        let tolerance = self
            .policy
            .linear_tolerance_mm
            .min(options.linear_deflection_mm / 16.);
        let receipt = ffi::tessellate_manufacturing(
            self.occt()?,
            options.linear_deflection_mm,
            options.angular_deflection_rad,
            options.maximum_triangles,
            tolerance,
        )
        .map_err(|e| native("manufacturing mesh", e))?;
        let mesh = &receipt.mesh;
        if !valid_parts(&mesh.vertices, &mesh.indices, &mesh.face_indices)
            || !receipt.maximum_displacement_mm.is_finite()
            || receipt.maximum_displacement_mm < 0.
            || receipt.maximum_displacement_mm > tolerance
        {
            return Err(error(
                "manufacturing mesh",
                GeometryErrorKind::InvalidResult,
                "Invalid native manufacturing mesh or seam receipt",
            ));
        }
        Ok(receipt.into())
    }
    /// Publish exact topology and its viewport overlays together. Pass no
    /// deflection for topology alone, without sampling edges or vertices.
    pub fn topology_with_display(
        &self,
        owner: EntityId,
        revision: u64,
        deflection_mm: Option<f64>,
    ) -> Result<TopologyDisplay, GeometryError> {
        let topology = self.topology_snapshot(owner, revision)?;
        let (edges, vertices) = match deflection_mm {
            None => (Vec::new(), Vec::new()),
            Some(deflection) if self.uses_rust_geometry() => {
                self.validate_display_deflection(deflection)?;
                let edges = if let Some(body) = self.pipe_hardware_body() {
                    body.display_edges(deflection, self.policy)?
                } else if let Some((body, frame)) = self.repaired_phillips_body() {
                    if body.has_degree11_source() {
                        body.degree11_display_edges(frame, deflection)?
                    } else {
                        body.periodic_display_edges(frame, deflection)?
                    }
                } else if let Some((body, frame)) = self.hanger_body() {
                    body.display_edges(frame, deflection, self.policy)?
                } else if let Some((body, frame)) = self.wood_insert_body() {
                    body.display_edges(frame, deflection, self.policy)?
                } else if let Some((body, frame)) = self.internal_thread_treatment_body() {
                    body.display_edges(frame, deflection, self.policy)?
                } else if let Some((body, frame)) = self.rounded_blind_bore_body() {
                    body.display_edges(frame, deflection)?
                } else if let Some((body, frame)) = self.rounded_thread_bore_body() {
                    body.display_edges(frame, deflection)?
                } else if let Some((body, frame)) = self.rounded_thread_body() {
                    body.display_edges(frame, deflection)?
                } else if let Some((body, frame)) = self.entry_chamfer_body() {
                    body.display_edges(frame, deflection)?
                } else if let Some((body, frame)) = self.helical_body() {
                    body.display_edges(frame, deflection)?
                } else if let Some(compound) = self.rust_compound() {
                    compound.display_edges(deflection)?
                } else {
                    crate::rust_display::edges(Some(self), &topology, deflection)?
                };
                (edges, crate::rust_display::vertices(&topology))
            }
            Some(deflection) => (self.display_edges(deflection)?, self.display_vertices()?),
        };
        Ok(TopologyDisplay {
            topology,
            edges,
            vertices,
        })
    }

    fn validate_display_deflection(&self, deflection_mm: f64) -> Result<(), GeometryError> {
        if !deflection_mm.is_finite()
            || deflection_mm < 0.001
            || deflection_mm > self.policy.maximum_coordinate_mm
        {
            return Err(input("Invalid edge display deflection"));
        }
        Ok(())
    }

    pub fn display_edges(&self, deflection_mm: f64) -> Result<Vec<DisplayEdge>, GeometryError> {
        self.validate_display_deflection(deflection_mm)?;
        if let Some(body) = self.pipe_hardware_body() {
            return body.display_edges(deflection_mm, self.policy);
        }
        if let Some((body, frame)) = self.repaired_phillips_body() {
            return if body.has_degree11_source() {
                body.degree11_display_edges(frame, deflection_mm)
            } else {
                body.periodic_display_edges(frame, deflection_mm)
            };
        }
        if let Some((body, frame)) = self.hanger_body() {
            return body.display_edges(frame, deflection_mm, self.policy);
        }
        if let Some((body, frame)) = self.wood_insert_body() {
            return body.display_edges(frame, deflection_mm, self.policy);
        }
        if let Some((body, frame)) = self.internal_thread_treatment_body() {
            return body.display_edges(frame, deflection_mm, self.policy);
        }
        if let Some((body, frame)) = self.rounded_blind_bore_body() {
            return body.display_edges(frame, deflection_mm);
        }
        if let Some((body, frame)) = self.rounded_thread_bore_body() {
            return body.display_edges(frame, deflection_mm);
        }
        if let Some((body, frame)) = self.rounded_thread_body() {
            return body.display_edges(frame, deflection_mm);
        }
        if let Some((body, frame)) = self.entry_chamfer_body() {
            return body.display_edges(frame, deflection_mm);
        }
        if let Some((body, frame)) = self.helical_body() {
            return body.display_edges(frame, deflection_mm);
        }
        if let Some(compound) = self.rust_compound() {
            return compound.display_edges(deflection_mm);
        }
        if let Some(topology) = self.rust_topology() {
            return crate::rust_display::edges(Some(self), &topology?, deflection_mm);
        }
        let edges = ffi::display_edges(self.occt()?, deflection_mm)
            .map_err(|e| native("display edges", e))?;
        edges
            .into_iter()
            .map(|edge| {
                if edge.points.len() < 6
                    || !edge.points.len().is_multiple_of(3)
                    || edge.points.iter().any(|v| !v.is_finite())
                {
                    return Err(error(
                        "display edges",
                        GeometryErrorKind::InvalidResult,
                        "Invalid native edge display data",
                    ));
                }
                Ok(DisplayEdge {
                    index: edge.index,
                    points_mm: edge.points.as_chunks::<3>().0.to_vec(),
                })
            })
            .collect()
    }
    pub(crate) fn sample_non_analytic_edge(
        &self,
        ordinal: u32,
        options: multicad_rust_kernel::display::Tessellation,
    ) -> Option<Result<Vec<[f64; 3]>, multicad_rust_kernel::display::MeshError>> {
        if let Some(body) = self.pipe_hardware_body() {
            return body.sample_edge(ordinal, options);
        }
        if let Some(body) = self.line8_nut_body() {
            return body.sample_edge(ordinal, options);
        }
        if let Some(body) = self.external_npt_stock_body() {
            return body.sample_edge(ordinal, options);
        }
        if let Some(body) = self.literal_npt_body() {
            return body.sample_edge(ordinal, options);
        }
        if let Some(body) = self.actual_clipped_head_body() {
            return body.sample_edge(ordinal, options);
        }
        if let Some(body) = self.clipped_npt_body() {
            return body.sample_edge(ordinal, options);
        }
        if let Some((body, frame)) = self.cosmetic_profiled_body_ref() {
            return body.sample_edge(frame, ordinal, options);
        }
        if let Some((body, frame)) = self.spline_cap_relief_body() {
            return body.sample_rim(frame, ordinal, options);
        }
        if let Some((body, frame)) = self.spline_sweep_body() {
            return body.sample_rim(frame, ordinal, options);
        }
        if let Some((body, frame)) = self.herringbone_body() {
            return body.sample_rail(frame, ordinal, options);
        }
        if let Some((body, frame)) = self.bevel_body() {
            return body.sample_edge(frame, ordinal, options);
        }
        if let Some((body, frame)) = self.worm_body() {
            return body.sample_edge(frame, ordinal, options);
        }
        if let Some((body, frame)) = self.helical_tube_cut_body() {
            return body.sample_edge(frame, ordinal, options);
        }
        if let Some((body, frame)) = self.radial_boundary_body() {
            return body.sample_edge(frame, ordinal, options);
        }
        if let Some((body, frame)) = self.cross_bore_body() {
            return body.sample_edge(frame, ordinal, options);
        }
        if let Some(body) = self.imported_brep_body() {
            return body.sample_edge(ordinal, options);
        }
        let (member, local) = self.rust_compound()?.edge_source(ordinal)?;
        member.sample_non_analytic_edge(local, options)
    }

    /// Query one physical face without asking an explicit provider to measure
    /// every other face. Other providers retain their existing bulk fallback.
    pub fn face(&self, index: u32) -> Result<Face, GeometryError> {
        if let Some(body) = self.actual_clipped_head_body() {
            return body.face(index);
        }
        if let Some((body, frame)) = self.repaired_phillips_body()
            && body.has_degree11_source()
        {
            return body.degree11_face(frame, index);
        }
        self.faces()?
            .into_iter()
            .find(|f| f.index == index)
            .ok_or_else(|| input("Face index is out of range"))
    }
    /// Cheap earned planar support; None denotes an actual curved face.
    pub fn face_plane(&self, index: u32) -> Result<Option<FacePlane>, GeometryError> {
        if let Some(body) = self.literal_npt_body() {
            return body.face_plane(index);
        }
        if let Some(body) = self.actual_clipped_head_body() {
            return Ok(body.plane(index)?.map(|f| (f.origin_mm, f.normal())));
        }
        if let Some((body, frame)) = self.repaired_phillips_body()
            && body.has_degree11_source()
        {
            return body.degree11_face_plane(frame, index);
        }
        Ok(self.face(index)?.plane)
    }

    pub fn display_vertices(&self) -> Result<Vec<DisplayVertex>, GeometryError> {
        if let Some(topology) = self.rust_topology() {
            return Ok(crate::rust_display::vertices(&topology?));
        }
        ffi::vertices(self.occt()?)
            .map_err(|e| native("display vertices", e))?
            .into_iter()
            .map(|v| {
                let position_mm = xyz(v.point);
                self.policy.point(position_mm)?;
                Ok(DisplayVertex {
                    index: v.index,
                    position_mm,
                })
            })
            .collect()
    }
    pub(crate) fn rust_topology(
        &self,
    ) -> Option<Result<multicad_core::references::TopologySnapshot, GeometryError>> {
        let feature = self.rust_feature()?;
        Some(self.topology_snapshot(feature, 0))
    }

    /// The source feature a Rust body names its topology after.
    pub(crate) fn rust_feature(&self) -> Option<EntityId> {
        Some(match &self.body {
            SolidBody::PolygonPrism(placed) => placed.feature,
            SolidBody::Prism(prism, _) => EntityId::from_uuid_bits(prism.feature()),
            SolidBody::AxialCylinder(cylinder, ..) => EntityId::from_uuid_bits(cylinder.feature()),
            SolidBody::Annulus(body, _) => body.feature,
            SolidBody::EntryChamfer(body, _) => body.feature,
            SolidBody::InternalThreadTreatment(body, _) => body.feature,
            SolidBody::RepairedPhillips(body, _) => body.feature,
            SolidBody::WoodInsert(body, _) => body.feature,
            SolidBody::Hanger(body, _) => body.feature,
            SolidBody::Helical(body, _) => body.feature,
            SolidBody::RoundedThreadBore(body, _) => body.feature,
            SolidBody::RoundedBlindBore(body, _) => body.feature,
            SolidBody::RoundedThread(body, _) => body.feature,
            SolidBody::HelicalTubeCut(body, _) => body.feature,
            SolidBody::CosmeticFastener(body, _) => body.feature,
            SolidBody::CosmeticProfiled(body, _) => body.feature,
            SolidBody::ClippedNpt(body) => body.feature,
            SolidBody::PipeHardware(body) => body.feature,
            SolidBody::LiteralNpt(body) => body.feature,
            SolidBody::ActualClippedHead(body) => body.feature,
            SolidBody::ExternalNptStock(body) => body.feature,
            SolidBody::Line8Nut(body) => body.feature,
            SolidBody::Linear(body, _) => body.feature,
            SolidBody::Revolved(body, _) => body.feature,
            SolidBody::RoundedSlot(body, _) => body.feature,
            SolidBody::NormalSweep(body, _) => body.feature,
            SolidBody::NormalRelief(body, _) => body.feature,
            SolidBody::SplineCapRelief(body, _) => body.feature,
            SolidBody::RadialBoundary(body, _) => body.feature,
            SolidBody::SplineSweep(body, _) => body.feature,
            SolidBody::SpurGear(body, _) => body.feature,
            SolidBody::InternalRing(body, _) => body.feature,
            SolidBody::Herringbone(body, _) => body.feature,
            SolidBody::Bevel(body, _) => body.feature,
            SolidBody::Worm(body, _) => body.feature,
            SolidBody::Sphere(body, _) => body.feature,
            SolidBody::ImportedBrep(body) => body.feature,
            SolidBody::Torus(body, _) => body.feature,
            SolidBody::ObliqueBore(body, _) => body.feature,
            SolidBody::CrossBore(body, _) => body.feature,
            SolidBody::RoundedCylinder(body, _) => body.feature,
            SolidBody::AxialStack(body, _) => body.feature,
            SolidBody::PlanarBoolean(body, _) => body.feature,
            SolidBody::NonplanarLoft(body, _) => body.feature,
            SolidBody::WarpedUnion(body, _) => body.feature,
            SolidBody::RustCompound(compound) => compound.feature,
            SolidBody::Orthogonal(body, _) => body.feature,
            SolidBody::OrthogonalSeats(body, _) => body.feature,
            SolidBody::Occt(_) => return None,
        })
    }
    /// Whether the active backend can encode this exact solid in STEP.
    pub fn supports_step_export(&self) -> bool {
        self.supports_step_shell(true)
    }

    fn supports_step_shell(&self, standalone: bool) -> bool {
        match &self.body {
            SolidBody::NormalRelief(..)
            | SolidBody::SplineCapRelief(..)
            | SolidBody::RadialBoundary(..)
            | SolidBody::RoundedThread(..)
            | SolidBody::RoundedThreadBore(..)
            | SolidBody::RoundedBlindBore(..) => true,
            SolidBody::Occt(_)
            | SolidBody::PolygonPrism(_)
            | SolidBody::AxialCylinder(..)
            | SolidBody::Annulus(..)
            | SolidBody::Linear(..)
            | SolidBody::NormalSweep(..)
            | SolidBody::NonplanarLoft(..)
            | SolidBody::SplineSweep(..)
            | SolidBody::SpurGear(..)
            | SolidBody::InternalRing(..)
            | SolidBody::AxialStack(..)
            | SolidBody::Sphere(..)
            | SolidBody::Torus(..)
            | SolidBody::ImportedBrep(..)
            | SolidBody::RoundedSlot(..) => true,
            SolidBody::Prism(prism, _) => prism.rounded_end().is_none(),
            SolidBody::RustCompound(compound) => compound
                .members()
                .all(|member| member.supports_step_shell(false)),
            SolidBody::EntryChamfer(..) => true,
            SolidBody::Helical(body, frame) => {
                (standalone || body.has_importable_step_closure())
                    && body.supports_step(*frame, self.policy)
            }
            SolidBody::Herringbone(body, frame) => body.supports_step(*frame, self.policy),
            SolidBody::Bevel(body, frame) => body.supports_step(*frame, self.policy),
            SolidBody::Worm(body, frame) => body.supports_step(*frame, self.policy),
            SolidBody::HelicalTubeCut(body, frame) => {
                multicad_rust_kernel::step::supports_helical_tube_cut(
                    body.cut,
                    frame.origin,
                    self.policy.linear_tolerance_mm,
                )
            }
            SolidBody::WoodInsert(body, frame) => body.supports_step(*frame, self.policy),
            SolidBody::ClippedNpt(body) => body.supports_step(self.policy),
            SolidBody::Hanger(body, frame) => body.supports_step(*frame, self.policy),
            SolidBody::OrthogonalSeats(body, frame) => body.supports_step(*frame, self.policy),
            SolidBody::LiteralNpt(body) => body.supports_step(),
            SolidBody::ActualClippedHead(_) => standalone,
            SolidBody::RepairedPhillips(..) => {
                self.rust_periodic_repaired_phillips_source().is_some()
            }
            SolidBody::PipeHardware(..)
            | SolidBody::Line8Nut(..)
            | SolidBody::ExternalNptStock(..)
            | SolidBody::CosmeticFastener(..)
            | SolidBody::CosmeticProfiled(..)
            | SolidBody::InternalThreadTreatment(..)
            | SolidBody::ObliqueBore(..)
            | SolidBody::CrossBore(..) => false,
            SolidBody::WarpedUnion(body, _) => {
                multicad_rust_kernel::step::supports_warped_union(&body.geometry)
            }
            SolidBody::PlanarBoolean(body, _) => {
                multicad_rust_kernel::step::supports_planar_union(&body.geometry)
            }
            SolidBody::RoundedCylinder(body, frame) => {
                multicad_rust_kernel::step::write_rounded_cylinder(
                    &body.geometry,
                    frame.origin,
                    frame.u,
                    frame.v,
                    self.policy.linear_tolerance_mm,
                )
                .is_ok()
            }
            SolidBody::Revolved(body, _) => {
                multicad_rust_kernel::step::supports_cylindrical_revolution(&body.geometry)
                    || multicad_rust_kernel::step::supports_conical_revolution(
                        &body.geometry,
                        self.policy.linear_tolerance_mm,
                    )
            }
            SolidBody::Orthogonal(body, _) => body.supports_step(),
        }
    }

    pub fn write_step(&self, path: &Path) -> Result<(), GeometryError> {
        let path_text = crate::step_path(path).map_err(|e| input(e.to_string()))?;
        if let SolidBody::Occt(shape) = &self.body {
            return ffi::write_step(&shape.0, path_text).map_err(|e| native("write STEP", e));
        }
        if let Some(result) =
            self.write_actual_clipped_head_step_with_verification(path, Default::default(), None)
        {
            return result.map(|_| ());
        }
        if let Some(result) = self.write_periodic_phillips_step_with_verification(
            path,
            self.policy.linear_tolerance_mm,
            None,
        ) {
            return result.map(|_| ());
        }
        let operation = "write STEP";
        let encoded = self.rust_step_text(operation)?;
        std::fs::write(path, encoded)
            .map_err(|e| error(operation, GeometryErrorKind::InvalidResult, e.to_string()))
    }

    /// Exact, deterministic text of this shape for lossless equality checks:
    /// OCCT BRep text for an OCCT shape, otherwise the Rust kernel's STEP
    /// encoding (round-trip f64 coordinates, fixed header). Shapes the Rust
    /// STEP writer cannot encode are refused rather than approximated.
    pub fn exact_text(&self) -> Result<String, GeometryError> {
        if matches!(self.body, SolidBody::Occt(_)) {
            return self.to_brep();
        }
        self.rust_step_text("exact text")
    }

    /// Rebuild an independent shape from [`Solid::exact_text`]: Rust STEP
    /// text goes through the strict Rust reader with fresh identity, and
    /// anything else is read as OCCT BRep text.
    pub fn from_exact_text(text: &str, policy: GeometryPolicy) -> Result<Self, GeometryError> {
        if text.starts_with("ISO-10303-21;") {
            policy.validate()?;
            return Self::import_rust_step(text, EntityId::new(), policy)?
                .without_authored_provenance();
        }
        Self::from_brep(text, policy)
    }

    fn rust_step_text(&self, operation: &'static str) -> Result<String, GeometryError> {
        use multicad_rust_kernel::step;

        if !self.supports_step_export() {
            return Err(error(
                operation,
                GeometryErrorKind::Unsupported,
                "Rust STEP export requires supported geometry with one closed exterior shell per solid",
            ));
        }
        self.imported_assembly_step()
            .unwrap_or_else(|| {
                let mut writer = step::Writer::new();
                self.append_step(&mut writer).and_then(|()| writer.finish())
            })
            .map_err(|e| {
                let kind = match e {
                    step::StepError::Unsupported => GeometryErrorKind::Unsupported,
                    step::StepError::Limit => GeometryErrorKind::Limit,
                    _ => GeometryErrorKind::InvalidResult,
                };
                error(operation, kind, format!("{e:?}"))
            })
    }

    pub(crate) fn append_step(
        &self,
        writer: &mut multicad_rust_kernel::step::Writer,
    ) -> Result<(), multicad_rust_kernel::step::StepError> {
        use multicad_rust_kernel::step;
        match &self.body {
            SolidBody::NormalRelief(body, frame) => body.append_step(writer, *frame),
            SolidBody::SplineCapRelief(body, frame) => body.append_step(writer, *frame),
            SolidBody::RadialBoundary(body, frame) => {
                body.append_step(writer, *frame, self.policy.linear_tolerance_mm)
            }
            SolidBody::PolygonPrism(placed) => {
                if let Some(treatments) = &placed.bore_treatments {
                    writer.append_treated_polygon_prism(
                        &placed.prism,
                        &treatments.geometry,
                        placed.frame.origin,
                        placed.frame.u,
                        placed.frame.v,
                        placed.reverse_sweep,
                    )
                } else if placed.pockets.is_empty() {
                    writer.append_polygon_prism(
                        &placed.prism,
                        placed.frame.origin,
                        placed.frame.u,
                        placed.frame.v,
                        placed.reverse_sweep,
                    )
                } else {
                    writer.append_pocketed_polygon_prism(
                        &placed.prism,
                        &placed.pockets,
                        placed.frame.origin,
                        placed.frame.u,
                        placed.frame.v,
                        placed.reverse_sweep,
                    )
                }
            }
            SolidBody::Prism(prism, frame) => {
                writer.append_rectangular_prism(prism, frame.origin, frame.u, frame.v)
            }
            SolidBody::AxialCylinder(cylinder, frame, _) => {
                if self.axial_cylinder().is_some() {
                    writer.append_axial_cylinder(cylinder, frame.u)
                } else {
                    writer.append_placed_cylinder(cylinder, frame.origin, frame.u, frame.v)
                }
            }
            SolidBody::Annulus(body, frame) => {
                writer.append_placed_annulus(&body.geometry, frame.origin, frame.u, frame.v)
            }
            SolidBody::Linear(body, frame) => body.append_step(writer, *frame),
            SolidBody::RoundedSlot(body, frame) => body.append_step(writer, *frame),
            SolidBody::NormalSweep(body, frame) => body.append_step(writer, *frame),
            SolidBody::SplineSweep(body, frame) => writer.append_spline_sweep(
                &body.sweep,
                frame.origin,
                frame.u,
                frame.v,
                self.policy.linear_tolerance_mm,
            ),
            SolidBody::Herringbone(body, frame) => body.append_step(writer, *frame, self.policy),
            SolidBody::Bevel(body, frame) => body.append_step(writer, *frame, self.policy),
            SolidBody::Worm(body, frame) => body.append_step(writer, *frame, self.policy),
            SolidBody::SpurGear(body, frame) => writer.append_spur_gear(
                &body.profile,
                frame.origin,
                frame.u,
                frame.v,
                self.policy.linear_tolerance_mm,
            ),
            SolidBody::InternalRing(body, frame) => writer.append_internal_ring(
                &body.ring,
                frame.origin,
                frame.u,
                frame.v,
                self.policy.linear_tolerance_mm,
            ),
            SolidBody::AxialStack(body, frame) => body.append_step(writer, *frame),
            SolidBody::PlanarBoolean(body, frame) => {
                writer.append_planar_union(&body.geometry, frame.origin, frame.u, frame.v)
            }
            SolidBody::NonplanarLoft(body, frame) => {
                writer.append_nonplanar_loft(&body.geometry, frame.origin, frame.u, frame.v)
            }
            SolidBody::WarpedUnion(body, frame) => {
                writer.append_warped_union(&body.geometry, frame.origin, frame.u, frame.v)
            }
            SolidBody::Revolved(body, frame) => {
                if let Some((base, radius, height, bore)) = body.circular_sweep(*frame) {
                    writer.append_placed_circular_sweep(
                        radius,
                        bore,
                        height,
                        base.origin,
                        base.u,
                        base.v,
                    )
                } else if step::supports_conical_revolution(
                    &body.geometry,
                    self.policy.linear_tolerance_mm,
                ) {
                    writer.append_conical_revolution(
                        &body.geometry,
                        frame.origin,
                        frame.u,
                        frame.v,
                        self.policy.linear_tolerance_mm,
                    )
                } else {
                    writer.append_cylindrical_revolution(
                        &body.geometry,
                        frame.origin,
                        frame.u,
                        frame.v,
                    )
                }
            }
            SolidBody::RustCompound(compound) => {
                for member in compound.members() {
                    member.append_step(writer)?;
                }
                Ok(())
            }
            SolidBody::EntryChamfer(body, frame) => body.append_step(writer, *frame, self.policy),
            SolidBody::Helical(body, frame) => body.append_step(writer, *frame, self.policy),
            SolidBody::RoundedThreadBore(body, frame) => {
                body.append_step(writer, *frame, self.policy)
            }
            SolidBody::RoundedBlindBore(body, frame) => {
                body.append_step(writer, *frame, self.policy)
            }
            SolidBody::RoundedThread(body, frame) => body.append_step(writer, *frame, self.policy),
            SolidBody::Sphere(body, frame) => {
                writer.append_sphere(body.geometry, frame.origin, frame.u, frame.v)
            }
            SolidBody::ImportedBrep(body) => writer.append_brep(&body.solid),
            SolidBody::RoundedCylinder(body, frame) => writer.append_rounded_cylinder(
                &body.geometry,
                frame.origin,
                frame.u,
                frame.v,
                self.policy.linear_tolerance_mm,
            ),
            SolidBody::Torus(body, frame) => {
                writer.append_torus(body.geometry, frame.origin, frame.u, frame.v)
            }
            SolidBody::HelicalTubeCut(body, frame) => writer.append_helical_tube_cut(
                body.cut,
                frame.origin,
                frame.u,
                frame.v,
                self.policy.linear_tolerance_mm,
            ),
            SolidBody::WoodInsert(body, frame) => body.append_step(writer, *frame, self.policy),
            SolidBody::ClippedNpt(body) => body.append_step(writer, self.policy),
            SolidBody::Hanger(body, frame) => body.append_step(writer, *frame, self.policy),
            SolidBody::OrthogonalSeats(body, frame) => {
                body.append_step(writer, *frame, self.policy)
            }
            SolidBody::LiteralNpt(body) => body.append_step(writer),
            SolidBody::ActualClippedHead(_) => Err(step::StepError::Unsupported),
            SolidBody::RepairedPhillips(..) => self
                .append_periodic_phillips_step(writer, None)
                .ok_or(step::StepError::Unsupported)?
                .map(|_| ())
                .map_err(|e| match e.kind {
                    GeometryErrorKind::Limit => step::StepError::Limit,
                    GeometryErrorKind::Unsupported => step::StepError::Unsupported,
                    _ => step::StepError::InvalidGeometry,
                }),
            SolidBody::PipeHardware(..)
            | SolidBody::Line8Nut(..)
            | SolidBody::ExternalNptStock(..)
            | SolidBody::CosmeticFastener(..)
            | SolidBody::CosmeticProfiled(..)
            | SolidBody::InternalThreadTreatment(..)
            | SolidBody::ObliqueBore(..)
            | SolidBody::CrossBore(..) => Err(step::StepError::Unsupported),
            SolidBody::Orthogonal(body, frame) => body.append_step(writer, *frame),
            SolidBody::Occt(_) => unreachable!("native STEP dispatched above"),
        }
    }

    /// Every solid in a STEP file as (representation name, exact BRep text).
    pub fn step_solids(path: &Path) -> Result<Vec<(String, String)>, GeometryError> {
        let path = crate::step_path(path).map_err(|e| input(e.to_string()))?;
        Ok(ffi::step_solids(path)
            .map_err(|e| native("read STEP", e))?
            .into_iter()
            .map(|s| (s.name, s.brep))
            .collect())
    }
    /// The exact shape as OCCT BRep text.
    pub fn to_brep(&self) -> Result<String, GeometryError> {
        ffi::shape_to_brep(self.occt()?).map_err(|e| native("write BRep", e))
    }
    /// Lossless native BRep representation digest for exact proof-cache keys.
    ///
    /// Hashes bounded (16 MiB) OCCT binary v4 output without triangulations or
    /// normals. Includes native geometry, tolerances, topology orientations and
    /// locations; excludes Rust provenance and display data. Equal physical
    /// geometry need not have an equal digest if its native representation or
    /// flags differ. This is not a canonical or cross-version geometry identity.
    /// Serialization/size failures return an error: callers must run the native
    /// query rather than substituting a lossy key or accepting a cached proof.
    pub fn exact_geometry_digest(
        &self,
    ) -> Result<multicad_core::project_file::ContentHash, GeometryError> {
        let bytes = ffi::shape_binary_geometry(self.occt()?, 16 * 1024 * 1024)
            .map_err(|e| native("binary geometry digest", e))?;
        Ok(multicad_core::project_file::ContentHash::of(&bytes))
    }
    /// Rebuild an exact shape from OCCT BRep text.
    pub fn from_brep(text: &str, policy: GeometryPolicy) -> Result<Self, GeometryError> {
        policy.validate()?;
        Self::adopt(
            ffi::shape_from_brep(text).map_err(|e| native("read BRep", e))?,
            policy,
            "read BRep",
        )
    }
    pub fn read_step(path: &Path, policy: GeometryPolicy) -> Result<Self, GeometryError> {
        policy.validate()?;
        if !crate::OCCT_LINKED {
            return Self::read_rust_step_file(path, policy);
        }
        let path = crate::step_path(path).map_err(|e| input(e.to_string()))?;
        Self::adopt(
            ffi::read_step(path).map_err(|e| native("read STEP", e))?,
            policy,
            "read STEP",
        )
    }
}

#[cfg(test)]
mod wood_insert_world_tests;

// Typed refusal for operand pairs outside every exact Rust provider.
impl Solid {
    /// The last Rust provider's operand refusal names both operand kinds,
    /// so an unsupported pair is never reported as a rectilinear limitation.
    fn no_exact_boolean(
        &self,
        other: &Self,
        operation: BooleanOperation,
        refusal: GeometryError,
    ) -> GeometryError {
        if refusal.operation != "orthogonal Boolean"
            || !refusal.message.contains("requires a rectilinear prism")
        {
            return refusal;
        }
        error(
            "Boolean",
            GeometryErrorKind::Unsupported,
            format!(
                "Rust Boolean has no exact provider for {} {operation:?} {}",
                self.body_kind_name(),
                other.body_kind_name()
            ),
        )
    }

    fn body_kind_name(&self) -> &'static str {
        match &self.body {
            SolidBody::Occt(..) => "Occt",
            SolidBody::Prism(..) => "Prism",
            SolidBody::PolygonPrism(..) => "PolygonPrism",
            SolidBody::AxialCylinder(..) => "AxialCylinder",
            SolidBody::Annulus(..) => "Annulus",
            SolidBody::CosmeticFastener(..) => "CosmeticFastener",
            SolidBody::CosmeticProfiled(..) => "CosmeticProfiled",
            SolidBody::ClippedNpt(..) => "ClippedNpt",
            SolidBody::PipeHardware(..) => "PipeHardware",
            SolidBody::LiteralNpt(..) => "LiteralNpt",
            SolidBody::ActualClippedHead(..) => "ActualClippedHead",
            SolidBody::ExternalNptStock(..) => "ExternalNptStock",
            SolidBody::Line8Nut(..) => "Line8Nut",
            SolidBody::Helical(..) => "Helical",
            SolidBody::RoundedThread(..) => "RoundedThread",
            SolidBody::EntryChamfer(..) => "EntryChamfer",
            SolidBody::RoundedThreadBore(..) => "RoundedThreadBore",
            SolidBody::RoundedBlindBore(..) => "RoundedBlindBore",
            SolidBody::InternalThreadTreatment(..) => "InternalThreadTreatment",
            SolidBody::RepairedPhillips(..) => "RepairedPhillips",
            SolidBody::WoodInsert(..) => "WoodInsert",
            SolidBody::Hanger(..) => "Hanger",
            SolidBody::HelicalTubeCut(..) => "HelicalTubeCut",
            SolidBody::Linear(..) => "Linear",
            SolidBody::Revolved(..) => "Revolved",
            SolidBody::RoundedSlot(..) => "RoundedSlot",
            SolidBody::NormalSweep(..) => "NormalSweep",
            SolidBody::NormalRelief(..) => "NormalRelief",
            SolidBody::SplineCapRelief(..) => "SplineCapRelief",
            SolidBody::RadialBoundary(..) => "RadialBoundary",
            SolidBody::NonplanarLoft(..) => "NonplanarLoft",
            SolidBody::SplineSweep(..) => "SplineSweep",
            SolidBody::SpurGear(..) => "SpurGear",
            SolidBody::InternalRing(..) => "InternalRing",
            SolidBody::Herringbone(..) => "Herringbone",
            SolidBody::Bevel(..) => "Bevel",
            SolidBody::Worm(..) => "Worm",
            SolidBody::Sphere(..) => "Sphere",
            SolidBody::Torus(..) => "Torus",
            SolidBody::ObliqueBore(..) => "ObliqueBore",
            SolidBody::CrossBore(..) => "CrossBore",
            SolidBody::RoundedCylinder(..) => "RoundedCylinder",
            SolidBody::AxialStack(..) => "AxialStack",
            SolidBody::RustCompound(..) => "RustCompound",
            SolidBody::Orthogonal(..) => "Orthogonal",
            SolidBody::ImportedBrep(..) => "ImportedBrep",
            SolidBody::OrthogonalSeats(..) => "OrthogonalSeats",
            SolidBody::PlanarBoolean(..) => "PlanarBoolean",
            SolidBody::WarpedUnion(..) => "WarpedUnion",
        }
    }
}
