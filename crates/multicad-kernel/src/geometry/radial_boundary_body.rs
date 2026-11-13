//! Facade ownership for an actual certified finite tensor/cone boundary.
use super::*;
use crate::rust_prism_topology::Builder;
use multicad_core::references::{AnalyticGeometry, TopologySnapshot};
use multicad_rust_kernel::radial_tensor::{
    Axis, Curve, Location, MassProperties, RadialBoundary, Surface, VolumeProperties,
};
use std::cell::RefCell;
use std::collections::BTreeSet;
mod display_cache;
mod identity;
mod pipe;
mod thread;
pub use pipe::RetainedAxialPipeThreadSource;
pub use thread::RetainedHelicalThreadSource;

const OPERATION: &str = "retained radial boundary";
const MASS_ACCURACY: f64 = 1e-4;
#[cfg(test)]
mod tests;

struct Placement {
    boundary: Rc<RadialBoundary>,
    properties: RefCell<Option<MassProperties>>,
    volume: RefCell<Option<VolumeProperties>>,
    bounds: RefCell<Option<Bounds>>,
}

pub(crate) struct RadialBoundaryBody {
    pub(super) feature: EntityId,
    source: Rc<Placement>,
    // One successful occurrence pose, rather than an unbounded pose history.
    placed: RefCell<Option<(PlaneFrame, Rc<Placement>)>>,
    policy: GeometryPolicy,
    face_keys: Vec<multicad_core::project_file::ContentHash>,
    edge_keys: Vec<multicad_core::project_file::ContentHash>,
    vertex_keys: Vec<multicad_core::project_file::ContentHash>,
    authored_thread: Option<Rc<thread::AuthoredThread>>,
    authored_pipe: Option<Rc<pipe::AuthoredPipe>>,
}

pub(super) fn source_signature(
    boundary: &RadialBoundary,
) -> multicad_core::project_file::ContentHash {
    identity::keys(boundary.input()).member
}

#[derive(Clone, Debug, PartialEq)]
pub struct RetainedStepMemberMass {
    pub member_index: u32,
    pub decoded_volume_bounds_mm3: [f64; 2],
    pub source_volume_bounds_mm3: [f64; 2],
}

#[derive(Clone, Debug, PartialEq)]
pub struct RetainedStepVerification {
    pub solid_count: u32,
    pub bounds_mm: [[f64; 3]; 2],
    pub retained_members: Vec<RetainedStepMemberMass>,
    /// Existing decoded scalar contracts of non-retained members are kept
    /// separate; they are not promoted to rigorous zero-width enclosures.
    pub decoded_scalar_members: Vec<(u32, f64)>,
}

fn failure(cause: multicad_rust_kernel::radial_tensor::Error) -> GeometryError {
    use multicad_rust_kernel::radial_tensor::Error;
    let kind = match cause {
        Error::Invalid => GeometryErrorKind::Input,
        Error::Unresolved | Error::Unsupported => GeometryErrorKind::Unsupported,
        Error::Limit => GeometryErrorKind::Limit,
        Error::Cancelled => GeometryErrorKind::Cancelled,
    };
    error(
        OPERATION,
        kind,
        format!("Retained source proof failed: {cause:?}"),
    )
}

fn checkpoint(interrupt: Option<&dyn Fn() -> bool>) -> Result<(), GeometryError> {
    if interrupt.is_some_and(|probe| probe()) {
        return Err(error(
            OPERATION,
            GeometryErrorKind::Cancelled,
            "Retained boundary work was superseded",
        ));
    }
    Ok(())
}

fn unit(v: [f64; 3]) -> [f64; 3] {
    let length = v.into_iter().map(|x| x * x).sum::<f64>().sqrt();
    v.map(|x| x / length)
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.into_iter().zip(b).map(|(a, b)| a * b).sum()
}

fn source_frame(axis: Axis, forward: bool) -> ReferenceFrame {
    let normal = unit(axis.normal).map(|v| if forward { v } else { -v });
    let projection = dot(axis.reference, normal);
    let u = unit(std::array::from_fn(|i| {
        axis.reference[i] - projection * normal[i]
    }));
    let v = [
        normal[1] * u[2] - normal[2] * u[1],
        normal[2] * u[0] - normal[0] * u[2],
        normal[0] * u[1] - normal[1] * u[0],
    ];
    ReferenceFrame {
        origin_mm: axis.origin,
        u,
        v,
    }
}

impl Solid {
    pub(super) fn qualify_imported_radial(
        owned: multicad_rust_kernel::radial_tensor::BoundaryOwnership,
        policy: GeometryPolicy,
    ) -> Result<RadialBoundary, GeometryError> {
        let tolerance = policy.linear_tolerance_mm;
        owned
            .bind_surfaces(tolerance, None)
            .map_err(failure)?
            .certify_axial(tolerance, policy.minimum_feature_mm, None)
            .map_err(failure)?
            .certify_caps(None)
            .map_err(failure)?
            .certify_ordering(None)
            .map_err(failure)?
            .certify_material(tolerance, None)
            .map_err(failure)
    }
    /// Admit already certified occupied material. This does not reconstruct a
    /// thread recipe from finite spline poles or replace the retained source.
    pub fn from_rust_radial_boundary(
        boundary: RadialBoundary,
        feature: EntityId,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        Self::from_rust_radial_boundary_with_interrupt(boundary, feature, policy, None)
    }

    pub fn from_rust_radial_boundary_with_interrupt(
        boundary: RadialBoundary,
        feature: EntityId,
        policy: GeometryPolicy,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, GeometryError> {
        Self::from_certified_radial_boundary(boundary, feature, policy, interrupt, None, None)
    }

    fn from_certified_radial_boundary(
        boundary: RadialBoundary,
        feature: EntityId,
        policy: GeometryPolicy,
        interrupt: Option<&dyn Fn() -> bool>,
        authored_thread: Option<Rc<thread::AuthoredThread>>,
        authored_pipe: Option<Rc<pipe::AuthoredPipe>>,
    ) -> Result<Self, GeometryError> {
        policy.validate()?;
        checkpoint(interrupt)?;
        if feature.is_nil() || policy.boolean_fuzzy_mm != 0. {
            return Err(error(
                OPERATION,
                GeometryErrorKind::Input,
                "Retained boundaries require a source identity and exact geometry policy",
            ));
        }
        if boundary.collars().radial_sewing_mm > policy.linear_tolerance_mm {
            return Err(error(
                OPERATION,
                GeometryErrorKind::Unsupported,
                "The retained source precision is coarser than the geometry policy",
            ));
        }
        let input = boundary.input();
        let points = input
            .vertices
            .iter()
            .chain(input.curves.iter().flat_map(|curve| match curve {
                Curve::Spline { poles, .. } => poles.as_slice(),
                Curve::Circle { axis, .. } => std::slice::from_ref(&axis.origin),
            }))
            .chain(input.surfaces.iter().flat_map(|surface| match surface {
                Surface::Tensor(t) => t.poles(),
                Surface::Plane(axis) | Surface::Cone { axis, .. } => {
                    std::slice::from_ref(&axis.origin)
                }
            }));
        for (i, point) in points.enumerate() {
            if i % 128 == 0 {
                checkpoint(interrupt)?;
            }
            policy.point(*point)?;
        }
        let keys = identity::keys(boundary.input());
        let body = RadialBoundaryBody {
            feature,
            source: Rc::new(Placement {
                boundary: Rc::new(boundary),
                properties: RefCell::new(None),
                volume: RefCell::new(None),
                bounds: RefCell::new(None),
            }),
            placed: RefCell::new(None),
            policy,
            face_keys: keys.faces,
            edge_keys: keys.edges,
            vertex_keys: keys.vertices,
            authored_thread,
            authored_pipe,
        };
        let bounds = body.bounds(PlaneFrame::default(), interrupt)?;
        policy.point(bounds.min)?;
        policy.point(bounds.max)?;
        checkpoint(interrupt)?;
        Ok(Self {
            body: SolidBody::RadialBoundary(Rc::new(body), PlaneFrame::default()),
            policy,
            provenance: vec![],
            rust_origins: vec![],
            retained_placement: None,
        })
    }

    /// The actual immutable source and its occurrence pose. Finite source table
    /// ordinals are evaluation data, rather than a reconstructed feature recipe.
    pub fn rust_radial_boundary(&self) -> Option<(&RadialBoundary, PlaneFrame)> {
        self.radial_boundary_body()
            .map(|(body, frame)| (&*body.source.boundary, frame))
    }

    pub fn rust_radial_boundary_classify(
        &self,
        point: [f64; 3],
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Location, GeometryError> {
        let (body, frame) = self.radial_boundary_body().ok_or_else(|| {
            error(
                OPERATION,
                GeometryErrorKind::Unsupported,
                "This solid has no retained radial boundary",
            )
        })?;
        body.placement(frame, interrupt)?
            .boundary
            .classify(point, interrupt)
            .map_err(failure)
    }

    /// Source-integral enclosures in the active occurrence's world coordinates.
    pub fn rust_radial_boundary_properties(
        &self,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<MassProperties, GeometryError> {
        let (body, frame) = self.radial_boundary_body().ok_or_else(|| {
            error(
                OPERATION,
                GeometryErrorKind::Unsupported,
                "This solid has no retained radial boundary",
            )
        })?;
        body.properties(frame, interrupt)
    }

    pub(crate) fn radial_boundary_body(&self) -> Option<(&RadialBoundaryBody, PlaneFrame)> {
        match &self.body {
            SolidBody::RadialBoundary(body, frame) => Some((body, *frame)),
            _ => None,
        }
    }
}

impl RadialBoundaryBody {
    pub(super) fn placed_boundary(
        &self,
        frame: PlaneFrame,
    ) -> Result<Rc<RadialBoundary>, GeometryError> {
        Ok(Rc::clone(&self.placement(frame, None)?.boundary))
    }
    pub(super) fn new_occurrence(&self) -> Self {
        Self {
            feature: self.feature,
            source: Rc::clone(&self.source),
            placed: RefCell::new(None),
            policy: self.policy,
            face_keys: self.face_keys.clone(),
            edge_keys: self.edge_keys.clone(),
            vertex_keys: self.vertex_keys.clone(),
            authored_thread: self.authored_thread.clone(),
            authored_pipe: self.authored_pipe.clone(),
        }
    }
    fn placement(
        &self,
        frame: PlaneFrame,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Rc<Placement>, GeometryError> {
        checkpoint(interrupt)?;
        if frame == PlaneFrame::default() {
            return Ok(Rc::clone(&self.source));
        }
        if let Some((cached, placement)) = self.placed.borrow().as_ref()
            && *cached == frame
        {
            return Ok(Rc::clone(placement));
        }
        frame.checked(self.policy)?;
        let boundary = self
            .source
            .boundary
            .transformed_with_interrupt(
                [frame.origin, frame.u, frame.v],
                self.policy.linear_tolerance_mm,
                self.policy.minimum_feature_mm,
                interrupt,
            )
            .map_err(failure)?;
        let placement = Rc::new(Placement {
            boundary: Rc::new(boundary),
            properties: RefCell::new(None),
            volume: RefCell::new(None),
            bounds: RefCell::new(None),
        });
        checkpoint(interrupt)?;
        *self.placed.borrow_mut() = Some((frame, Rc::clone(&placement)));
        Ok(placement)
    }

    pub(super) fn properties(
        &self,
        frame: PlaneFrame,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<MassProperties, GeometryError> {
        let placement = self.placement(frame, interrupt)?;
        if let Some(properties) = placement.properties.borrow().as_ref() {
            return Ok(properties.clone());
        }
        let mut properties = placement
            .boundary
            .properties_with_interrupt(MASS_ACCURACY, interrupt)
            .map_err(failure)?;
        checkpoint(interrupt)?;
        // Keep the scalar report independent of whether full properties were
        // requested first. Each property retains a valid independent enclosure.
        properties.volume_mm3 = Self::volume_at(&placement, interrupt)?.volume_mm3;
        checkpoint(interrupt)?;
        *placement.properties.borrow_mut() = Some(properties.clone());
        Ok(properties)
    }

    fn volume_at(
        placement: &Placement,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<VolumeProperties, GeometryError> {
        checkpoint(interrupt)?;
        if let Some(volume) = *placement.volume.borrow() {
            return Ok(volume);
        }
        let volume = placement
            .boundary
            .volume_with_interrupt(MASS_ACCURACY, interrupt)
            .map_err(failure)?;
        checkpoint(interrupt)?;
        *placement.volume.borrow_mut() = Some(volume);
        Ok(volume)
    }

    pub(super) fn volume_properties(
        &self,
        frame: PlaneFrame,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<VolumeProperties, GeometryError> {
        let placement = self.placement(frame, interrupt)?;
        Self::volume_at(&placement, interrupt)
    }
    pub(super) fn volume_mm3(&self, frame: PlaneFrame) -> Result<f64, GeometryError> {
        let placement = self.placement(frame, None)?;
        let [a, b] = Self::volume_at(&placement, None)?.volume_mm3;
        Ok(a + (b - a) * 0.5)
    }

    pub(super) fn bounds(
        &self,
        frame: PlaneFrame,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Bounds, GeometryError> {
        let placement = self.placement(frame, interrupt)?;
        if let Some(bounds) = *placement.bounds.borrow() {
            return Ok(bounds);
        }
        let accuracy = 64. * self.policy.linear_tolerance_mm;
        let mut min = [0.; 3];
        let mut max = [0.; 3];
        for i in 0..3 {
            let direction = std::array::from_fn(|k| if k == i { 1. } else { 0. });
            max[i] = placement
                .boundary
                .support_bounds(direction, accuracy, interrupt)
                .map_err(failure)?[1];
            min[i] = -placement
                .boundary
                .support_bounds(direction.map(|v| -v), accuracy, interrupt)
                .map_err(failure)?[1];
        }
        let bounds = Bounds { min, max };
        checkpoint(interrupt)?;
        *placement.bounds.borrow_mut() = Some(bounds);
        Ok(bounds)
    }

    pub(super) fn faces(&self, frame: PlaneFrame) -> Result<Vec<Face>, GeometryError> {
        let placement = self.placement(frame, None)?;
        Ok(self
            .properties(frame, None)?
            .faces
            .into_iter()
            .map(|face| Face {
                index: face.face_index as u32,
                plane: face.planar_support.map(|axis| {
                    (
                        axis.origin,
                        source_frame(
                            axis,
                            placement.boundary.input().faces[face.face_index].forward,
                        )
                        .normal(),
                    )
                }),
                center_mm: face.centroid_mm.map(|[a, b]| a + (b - a) * 0.5),
                area_mm2: face.area_mm2[0] + (face.area_mm2[1] - face.area_mm2[0]) * 0.5,
            })
            .collect())
    }

    pub(crate) fn planar_loops(&self, face: u32) -> Vec<Vec<(u32, bool)>> {
        let input = self.source.boundary.input();
        let Some(face) = input.faces.get(face as usize) else {
            return vec![];
        };
        if !matches!(input.surfaces[face.surface], Surface::Plane(_)) {
            return vec![];
        }
        vec![
            face.outer
                .iter()
                .map(|&(edge, forward)| (edge as u32, !forward))
                .collect(),
        ]
    }

    pub(crate) fn topology(
        &self,
        frame: PlaneFrame,
        owner: EntityId,
        revision: u64,
    ) -> Result<TopologySnapshot, GeometryError> {
        let placement = self.placement(frame, None)?;
        let input = placement.boundary.input();
        let mut builder = Builder::new(
            (0..input.faces.len())
                .map(|i| {
                    if let Some(authored) = &self.authored_thread {
                        authored.face_origins(self.feature, i)
                    } else if let Some(pipe) = &self.authored_pipe {
                        pipe.face_origins(self.feature, i)
                    } else {
                        vec![Origin {
                            source: self.feature,
                            role: format!("retained.face/{}", self.face_keys[i]),
                        }]
                    }
                })
                .collect(),
        );
        let mut edge_faces = vec![BTreeSet::new(); input.edges.len()];
        let mut vertex_faces = vec![BTreeSet::new(); input.vertices.len()];
        for (i, face) in input.faces.iter().enumerate() {
            for &(edge, _) in &face.outer {
                edge_faces[edge].insert(i);
                for vertex in input.edges[edge].vertices {
                    vertex_faces[vertex].insert(i);
                }
            }
        }
        for (i, edge) in input.edges.iter().enumerate() {
            let incident = edge_faces[i].iter().copied().collect::<Vec<_>>();
            match input.curves[edge.curve] {
                Curve::Spline { .. } => builder.non_analytic_edge(&incident),
                Curve::Circle { axis, radius_mm } => {
                    let frame = source_frame(axis, true);
                    let phase = |vertex: usize| {
                        let p =
                            std::array::from_fn(|k| input.vertices[vertex][k] - frame.origin_mm[k]);
                        dot(p, frame.v).atan2(dot(p, frame.u))
                    };
                    let start = phase(edge.vertices[0]);
                    let end = phase(edge.vertices[1]);
                    let sweep = if edge.forward {
                        (end - start).rem_euclid(std::f64::consts::TAU)
                    } else {
                        -(start - end).rem_euclid(std::f64::consts::TAU)
                    };
                    builder.non_face(
                        AnalyticGeometry::CircularArc {
                            frame,
                            radius_mm,
                            start_rad: start,
                            sweep_rad: sweep,
                        },
                        &incident,
                    );
                }
            }
        }
        for (i, &point) in input.vertices.iter().enumerate() {
            builder.non_face(
                AnalyticGeometry::Point { position_mm: point },
                &vertex_faces[i].iter().copied().collect::<Vec<_>>(),
            );
        }
        for (i, face) in input.faces.iter().enumerate() {
            let plane = match input.surfaces[face.surface] {
                Surface::Plane(axis) => Some(source_frame(axis, face.forward)),
                _ => None,
            };
            builder.face(
                i,
                plane,
                plane,
                face.outer
                    .iter()
                    .flat_map(|&(e, _)| edge_faces[e].iter().copied())
                    .filter(|&neighbor| neighbor != i),
            );
        }
        let mut elements = builder.into_elements();
        for element in &mut elements {
            if let Some(authored) = &self.authored_thread {
                authored.apply_origin(self.feature, element);
                continue;
            }
            if let Some(pipe) = &self.authored_pipe {
                pipe.apply_origin(self.feature, element);
                continue;
            }
            let keys = match element.key.kind {
                multicad_core::references::TopologyKind::Face => None,
                multicad_core::references::TopologyKind::Edge => Some((&self.edge_keys, "edge")),
                multicad_core::references::TopologyKind::Vertex => {
                    Some((&self.vertex_keys, "vertex"))
                }
            };
            if let Some((keys, kind)) = keys {
                element.origins.insert(Origin {
                    source: self.feature,
                    role: format!("retained.{kind}/{}", keys[element.key.ordinal as usize]),
                });
            }
        }
        Ok(TopologySnapshot {
            owner,
            revision,
            elements,
        })
    }

    pub(super) fn mesh(
        &self,
        frame: PlaneFrame,
        options: multicad_rust_kernel::display::Tessellation,
        usage: MeshUse,
    ) -> Result<Mesh, GeometryError> {
        let options = multicad_rust_kernel::display::Tessellation {
            triangle_limit: options.triangle_limit.min(4_000_000),
            ..options
        };
        let placement = self.placement(frame, None)?;
        display_cache::mesh(&placement, options, usage, None, || {
            let mesh = placement
                .boundary
                .tessellate_with_interrupt(options, None)
                .map_err(failure)?;
            Ok(Mesh {
                vertices: mesh.vertices,
                indices: mesh.indices,
                face_indices: mesh.face_indices,
            })
        })
    }

    pub(super) fn sample_edge(
        &self,
        frame: PlaneFrame,
        ordinal: u32,
        options: multicad_rust_kernel::display::Tessellation,
    ) -> Option<Result<Vec<[f64; 3]>, multicad_rust_kernel::display::MeshError>> {
        let placement = match self.placement(frame, None) {
            Ok(placement) => placement,
            Err(error) => {
                return Some(Err(match error.kind {
                    GeometryErrorKind::Limit => {
                        multicad_rust_kernel::display::MeshError::BudgetExceeded
                    }
                    GeometryErrorKind::Input => {
                        multicad_rust_kernel::display::MeshError::InvalidOptions
                    }
                    _ => multicad_rust_kernel::display::MeshError::UnresolvedBoundary,
                }));
            }
        };
        if ordinal as usize >= placement.boundary.input().edges.len() {
            return None;
        }
        Some(
            placement
                .boundary
                .sample_edge_with_interrupt(
                    ordinal as usize,
                    options.deflection,
                    options.triangle_limit.min(1_000_000),
                    None,
                )
                .map_err(|cause| match cause {
                    multicad_rust_kernel::radial_tensor::Error::Invalid => {
                        multicad_rust_kernel::display::MeshError::InvalidOptions
                    }
                    multicad_rust_kernel::radial_tensor::Error::Limit => {
                        multicad_rust_kernel::display::MeshError::BudgetExceeded
                    }
                    _ => multicad_rust_kernel::display::MeshError::UnresolvedBoundary,
                }),
        )
    }

    pub(super) fn append_step(
        &self,
        writer: &mut multicad_rust_kernel::step::Writer,
        frame: PlaneFrame,
        tolerance: f64,
    ) -> Result<(), multicad_rust_kernel::step::StepError> {
        let placement = self.placement(frame, None).map_err(|cause| {
            if cause.kind == GeometryErrorKind::Limit {
                multicad_rust_kernel::step::StepError::Limit
            } else {
                multicad_rust_kernel::step::StepError::Unsupported
            }
        })?;
        writer.append_radial_boundary(&placement.boundary, tolerance)
    }
}
