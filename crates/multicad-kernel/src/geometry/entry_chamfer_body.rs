//! Owned conical-entry thread material; no serialized recipe or mesh is geometry.
use super::*;
use crate::HelicalThread;
type Point = [f64; 3];
use crate::rust_prism_topology::position;
use multicad_core::references::{
    AnalyticGeometry, AnchoredFrame, ReferenceFrame, TopologyElement, TopologyKey, TopologyKind,
    TopologySnapshot,
};
use multicad_rust_kernel::thread_end_relief::Intent as TerminalIntent;
use multicad_rust_kernel::thread_entry_chamfer::{
    self as source, BoundaryTopology, Chamfer, Enclosure, MassProperties, PropertiesOptions,
    QueryOptions,
};
use std::cell::RefCell;
use std::collections::BTreeSet;
mod display;
mod embedding;
mod placement;
mod references;
mod topology;
pub use references::EntryChamferEdge;

/// Complete retained material, including the original tool and optional neck.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EntryChamferData {
    pub source: Chamfer,
    pub frame: PlaneFrame,
}
/// A full source-known finite-file proof, with source-integral mass bounds.
/// This deliberately contains no decoded-file volume or blind-import claim.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EntryChamferBoundaryProof {
    pub source_volume_bounds_mm3: [f64; 2],
    pub verified_tolerance_mm: f64,
    pub maximum_polynomial_trim_error_mm: f64,
    pub analytic_closure_tolerance_mm: f64,
    pub face_count: usize,
}
pub(crate) struct EntryChamferBody {
    pub(crate) geometry: Chamfer,
    pub(super) feature: EntityId,
    topology: BoundaryTopology,
    references: references::ReferenceView,
    properties: MassProperties,
    tolerance: f64,
    bounds_cache: RefCell<Option<(PlaneFrame, Bounds)>>,
}
fn failure(cause: source::Error) -> GeometryError {
    error(
        "Rust thread entry chamfer",
        match cause {
            source::Error::Invalid => GeometryErrorKind::Input,
            source::Error::Limit => GeometryErrorKind::Limit,
            source::Error::Unresolved => GeometryErrorKind::Unsupported,
            source::Error::Cancelled => GeometryErrorKind::Cancelled,
        },
        format!("Entry chamfer source proof failed: {cause:?}"),
    )
}
fn middle(value: Enclosure) -> f64 {
    value.lower.midpoint(value.upper)
}
fn failure_at(stage: &'static str, cause: source::Error) -> GeometryError {
    let mut error = failure(cause);
    error.operation = stage;
    error
}
impl Solid {
    /// Actual retained metric entry material, constructed without exchange I/O.
    #[allow(clippy::too_many_arguments)]
    pub fn rust_helical_thread_entry_chamfer(
        spec: HelicalThread,
        root_radius_mm: Option<f64>,
        radial_depth_mm: f64,
        terminal: Option<TerminalIntent>,
        feature: EntityId,
        frame: PlaneFrame,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        Self::rust_helical_thread_entry_chamfer_with_interrupt(
            spec,
            root_radius_mm,
            radial_depth_mm,
            terminal,
            feature,
            frame,
            policy,
            None,
        )
    }
    #[allow(clippy::too_many_arguments)]
    pub fn rust_helical_thread_entry_chamfer_with_interrupt(
        spec: HelicalThread,
        root_radius_mm: Option<f64>,
        radial_depth_mm: f64,
        terminal: Option<TerminalIntent>,
        feature: EntityId,
        frame: PlaneFrame,
        policy: GeometryPolicy,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, GeometryError> {
        if interrupt.is_some_and(|f| f()) {
            return Err(failure(source::Error::Cancelled));
        }
        policy.validate()?;
        frame.checked(policy)?;
        if feature.is_nil() {
            return Err(input("Entry chamfer identity must not be nil"));
        }
        policy.length(radial_depth_mm)?;
        if let Some(tail) = terminal {
            policy.length(tail.width_mm)?;
            policy.length(tail.radial_depth_mm)?;
        }
        let tool = helical_body::checked_thread(spec, root_radius_mm, policy)?;
        let geometry = Chamfer::new(tool, radial_depth_mm, terminal)
            .map_err(|e| failure_at("Entry chamfer material", e))?;
        frame.checked_metric_extent(
            [spec.major_radius_mm, spec.major_radius_mm, spec.length_mm],
            policy,
        )?;
        policy.length(geometry.entry_radius_mm())?;
        policy.length(geometry.material_length_mm())?;
        if let Some(relief) = geometry.terminal_relief() {
            policy.length(relief.neck_radius_mm())?;
        }
        let tolerance = policy.linear_tolerance_mm;
        let topology = geometry
            .boundary_topology_with_interrupt(
                QueryOptions {
                    tolerance: tolerance / 16.,
                    ..Default::default()
                },
                interrupt,
            )
            .map_err(|e| failure_at("Entry chamfer source topology", e))?;
        // Physical integrals have their own documented dimensional accuracy;
        // their midpoint face centers are not attachment/geometry anchors.
        // Keep the complete immutable enclosures, rather than manufacture a
        // tighter claim from rounded face centroids or a viewport mesh.
        let properties = geometry
            .properties_with_interrupt(PropertiesOptions::default(), interrupt)
            .map_err(|e| failure_at("Entry chamfer physical integrals", e))?;
        let references = references::ReferenceView::new(&geometry, &topology)?;
        let body = EntryChamferBody {
            geometry,
            feature,
            topology,
            references,
            properties,
            tolerance,
            bounds_cache: RefCell::new(None),
        };
        let bounds = body.bounds_with_interrupt(frame, interrupt)?;
        policy.point(bounds.min)?;
        policy.point(bounds.max)?;
        Ok(Self {
            body: SolidBody::EntryChamfer(Rc::new(body), frame),
            policy,
            provenance: vec![],
            rust_origins: vec![],
            retained_placement: None,
        })
    }
    pub(crate) fn entry_chamfer_body(&self) -> Option<(&EntryChamferBody, PlaneFrame)> {
        match &self.body {
            SolidBody::EntryChamfer(body, frame) => Some((body, *frame)),
            _ => None,
        }
    }
    /// Typed geometry handoff; a consumer must include the cone and whole neck.
    pub fn rust_thread_entry_chamfer_data(&self) -> Option<EntryChamferData> {
        self.entry_chamfer_body()
            .map(|(body, frame)| EntryChamferData {
                source: body.geometry,
                frame,
            })
    }
    /// Cached unit-density integral enclosures in the source chart, with its
    /// current placement. Midpoint UI properties do not replace these bounds.
    pub fn rust_thread_entry_chamfer_properties(&self) -> Option<(&MassProperties, PlaneFrame)> {
        self.entry_chamfer_body()
            .map(|(body, frame)| (&body.properties, frame))
    }
    pub fn rust_thread_entry_chamfer_topology(&self) -> Option<(&BoundaryTopology, PlaneFrame)> {
        self.entry_chamfer_body()
            .map(|(body, frame)| (&body.topology, frame))
    }
    /// Each published reference edge owns one complete source curve or an
    /// exact ordered chain of source coedges. Raw and reference IDs differ.
    pub fn rust_thread_entry_chamfer_edges(&self) -> Option<(&[EntryChamferEdge], PlaneFrame)> {
        self.entry_chamfer_body()
            .map(|(body, frame)| (body.references.edges.as_slice(), frame))
    }
    pub fn verify_entry_chamfer_step(
        &self,
        input: &str,
    ) -> Option<Result<EntryChamferBoundaryProof, GeometryError>> {
        self.verify_entry_chamfer_step_with_tolerance(input, self.policy.linear_tolerance_mm)
    }
    pub fn verify_entry_chamfer_step_with_tolerance(
        &self,
        input: &str,
        tolerance: f64,
    ) -> Option<Result<EntryChamferBoundaryProof, GeometryError>> {
        let (body, frame) = self.entry_chamfer_body()?;
        if !tolerance.is_finite() || tolerance <= 0. {
            return Some(Err(failure(source::Error::Invalid)));
        }
        let tolerance = tolerance.min(self.policy.linear_tolerance_mm);
        Some(multicad_rust_kernel::step_read::verify_thread_entry_chamfer(input,&body.geometry,frame.origin,frame.u,frame.v,tolerance).map(|proof|EntryChamferBoundaryProof {
            source_volume_bounds_mm3:[body.properties.volume_mm3.lower,body.properties.volume_mm3.upper],
            verified_tolerance_mm:tolerance,
            maximum_polynomial_trim_error_mm:proof.maximum_polynomial_trim_error_mm,
            analytic_closure_tolerance_mm:tolerance/16.,
            face_count:proof.face_count,
        }).map_err(|cause|error("Entry chamfer STEP proof", match cause {
            multicad_rust_kernel::step_read::EntryChamferReadError::Limit=>GeometryErrorKind::Limit,
            multicad_rust_kernel::step_read::EntryChamferReadError::Cancelled=>GeometryErrorKind::Cancelled,
            _=>GeometryErrorKind::InvalidResult,
        },format!("Actual finite STEP boundary does not prove the current entry material: {cause:?}"))))
    }
}
impl EntryChamferBody {
    fn origin(&self, role: impl Into<String>) -> Origin {
        Origin {
            source: self.feature,
            role: role.into(),
        }
    }
    pub(super) fn bounds(&self, frame: PlaneFrame) -> Result<Bounds, GeometryError> {
        self.bounds_with_interrupt(frame, None)
    }
    fn bounds_with_interrupt(
        &self,
        frame: PlaneFrame,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Bounds, GeometryError> {
        if interrupt.is_some_and(|f| f()) {
            return Err(failure(source::Error::Cancelled));
        }
        if let Some((key, bounds)) = *self.bounds_cache.borrow()
            && key == frame
        {
            return Ok(bounds);
        }
        let axes = [frame.u, frame.v, frame.normal()];
        let mut min = [0.; 3];
        let mut max = [0.; 3];
        for (sign, target) in [(-1., &mut min), (1., &mut max)] {
            for i in 0..3 {
                let objective = axes.map(|a| sign * a[i]);
                let bound = self
                    .geometry
                    .support_with_interrupt(
                        objective,
                        QueryOptions {
                            tolerance: self.tolerance / 16.,
                            ..Default::default()
                        },
                        interrupt,
                    )
                    .map_err(failure)?;
                let value = frame.origin[i] + sign * bound.upper;
                target[i] = if sign < 0. {
                    value.next_down()
                } else {
                    value.next_up()
                };
            }
        }
        let result = Bounds { min, max };
        *self.bounds_cache.borrow_mut() = Some((frame, result));
        Ok(result)
    }
    fn cap_frame(&self, frame: PlaneFrame, role: source::Face) -> Option<ReferenceFrame> {
        let (z, upward) = match role {
            source::Face::EntryCap => (0., false),
            source::Face::EndCap => (self.geometry.material_length_mm(), true),
            source::Face::FarCap => (self.geometry.full_length_mm(), true),
            _ => return None,
        };
        Some(ReferenceFrame {
            origin_mm: position(frame, [0., 0., z]),
            u: frame.u,
            v: frame.v.map(|v| if upward { v } else { -v }),
        })
    }
    pub(super) fn faces(&self, frame: PlaneFrame) -> Vec<Face> {
        self.properties
            .faces
            .iter()
            .enumerate()
            .map(|(i, p)| Face {
                index: i as u32,
                plane: self
                    .cap_frame(frame, p.face)
                    .map(|f| (f.origin_mm, f.normal())),
                center_mm: position(frame, p.centroid_mm.map(middle)),
                area_mm2: middle(p.area_mm2),
            })
            .collect()
    }
    pub(super) fn volume(&self) -> f64 {
        middle(self.properties.volume_mm3)
    }
    pub(super) fn contains(&self, frame: PlaneFrame, point: Point) -> bool {
        let n = frame.normal();
        let rows = [cross(frame.v, n), cross(n, frame.u), n];
        let det = dot(frame.u, rows[0]);
        let delta = std::array::from_fn(|i| point[i] - frame.origin[i]);
        let p = rows.map(|r| dot(r, delta) / det);
        p[2] >= 0.
            && p[2] <= self.geometry.full_length_mm()
            && p[0].hypot(p[1]) <= self.geometry.radius(p[1].atan2(p[0]), p[2])
    }
    pub(super) fn mesh(
        &self,
        frame: PlaneFrame,
        options: multicad_rust_kernel::display::Tessellation,
    ) -> Result<Mesh, GeometryError> {
        placement::mesh(&self.geometry, frame, options)
    }
    pub(super) fn append_step(
        &self,
        writer: &mut multicad_rust_kernel::step::Writer,
        frame: PlaneFrame,
        policy: GeometryPolicy,
    ) -> Result<(), multicad_rust_kernel::step::StepError> {
        use multicad_rust_kernel::step::{EntryChamferStepError as E, StepError as S};
        writer
            .append_thread_entry_chamfer(
                &self.geometry,
                frame.origin,
                frame.u,
                frame.v,
                policy.linear_tolerance_mm,
            )
            .map_err(|e| match e {
                E::InvalidGeometry => S::InvalidGeometry,
                E::InvalidPlacement => S::InvalidPlacement,
                E::Limit => S::Limit,
                // This ordinary append installs no interrupt callback. Controlled
                // source clients use the pure typed append API directly.
                E::Unsupported | E::Cancelled => S::Unsupported,
            })
    }
}
