//! True finite tapered pipe thread material joined to hex stock. Original
//! authored intent remains distinct from the actual retained thread boundary.
use super::*;
use crate::geometry::internal_thread_treatment_body::interval::Interval as I;
use multicad_core::references::{
    AnalyticGeometry, AnchoredFrame, TopologyElement, TopologyKey, TopologyKind,
};
use multicad_rust_kernel::{
    axial_pipe_thread::{
        Error as Cause,
        authored::{
            AuthoredSpec, AuthoredThread,
            hardware::{Material, Options as MaterialOptions},
        },
    },
    display::Tessellation,
};
use std::{
    cell::RefCell,
    collections::{BTreeMap, BTreeSet},
};
mod factory;
mod mesh;
pub use factory::{PipeThreadFactory, PipeThreadFactoryStats};
mod placement;
mod query;
mod source;
#[cfg(test)]
mod tests;
mod topology;
use placement::Placement;
#[derive(Clone)]
enum FaceMeasure {
    Ready { area: I, center: [I; 3] },
    Retained(u32),
}
#[derive(Clone)]
struct LocalFace {
    measure: FaceMeasure,
    frame: Option<ReferenceFrame>,
    normal: Option<[I; 3]>,
}
type FaceMeasures = BTreeMap<u32, (I, [I; 3])>;
struct ThreadSource {
    feature: EntityId,
    spec: AuthoredSpec,
    policy: GeometryPolicy,
    thread: Solid,
    topology: TopologySnapshot,
    measured: FaceMeasures,
    lazy_measured: RefCell<Option<FaceMeasures>>,
    #[cfg(test)]
    descriptor_queries: std::cell::Cell<usize>,
    caps: [u32; 2],
    volume: I,
    radius: I,
}
impl ThreadSource {
    fn measure(
        &self,
        index: u32,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<(I, [I; 3]), GeometryError> {
        checkpoint(interrupt)?;
        if let Some(value) = self.measured.get(&index) {
            return Ok(*value);
        }
        if self.lazy_measured.borrow().is_none() {
            let (boundary, pose) = self
                .thread
                .rust_radial_boundary()
                .ok_or_else(|| failure(Cause::Unsupported))?;
            if pose != PlaneFrame::default() {
                return Err(failure(Cause::Unsupported));
            }
            #[cfg(test)]
            self.descriptor_queries
                .set(self.descriptor_queries.get() + 1);
            let properties = boundary
                .face_properties_with_interrupt(source::MEASURE_ACCURACY, interrupt)
                .map_err(radial_failure)?;
            let measured = properties
                .into_iter()
                .map(|face| {
                    (
                        face.face_index as u32,
                        (interval(face.area_mm2), face.centroid_mm.map(interval)),
                    )
                })
                .collect();
            checkpoint(interrupt)?;
            *self.lazy_measured.borrow_mut() = Some(measured);
        }
        self.lazy_measured
            .borrow()
            .as_ref()
            .and_then(|m| m.get(&index).copied())
            .ok_or_else(|| failure(Cause::Unsupported))
    }
    fn retained_bytes(&self) -> usize {
        // Actual complete retained source capacities/unique Arc poles, plus a
        // conservative 1 MiB for <=34 face/112 edge/80 vertex Native descriptors,
        // topology/role strings and the small measured face map. Clipped
        // plans and their <=34 finite face measures fit this same fixed bound.
        1024 * 1024
            + self
                .thread
                .rust_radial_boundary()
                .map_or(0, |(b, _)| b.retained_budget_bytes())
    }
}
struct Source {
    feature: EntityId,
    intent: Material,
    thread: Solid,
    _prepared: Rc<ThreadSource>,
    faces: Vec<LocalFace>,
    face_map: BTreeMap<u32, u32>,
    topology: TopologySnapshot,
    loops: BTreeMap<u32, Vec<Vec<(u32, bool)>>>,
    hex_vertices: [[I; 3]; 12],
    edge_offset: u32,
    vertex_offset: u32,
    bounds: [[f64; 2]; 3],
    volume: I,
    clearance: I,
    minimum_wall: f64,
}
pub(crate) struct PipeHardwareBody {
    pub(super) feature: EntityId,
    source: Rc<Source>,
    placement: Placement,
    display: RefCell<Option<mesh::Cached>>,
    mesh_builds: std::cell::Cell<usize>,
}
fn failure(cause: Cause) -> GeometryError {
    error(
        "authored pipe hardware",
        match cause {
            Cause::Invalid => GeometryErrorKind::Input,
            Cause::Limit => GeometryErrorKind::Limit,
            Cause::Cancelled => GeometryErrorKind::Cancelled,
            Cause::Unsupported | Cause::Unresolved => GeometryErrorKind::Unsupported,
        },
        format!("True pipe hardware qualification failed: {cause:?}"),
    )
}
fn checkpoint(interrupt: Option<&dyn Fn() -> bool>) -> Result<(), GeometryError> {
    if interrupt.is_some_and(|f| f()) {
        Err(failure(Cause::Cancelled))
    } else {
        Ok(())
    }
}
fn radial_failure(cause: multicad_rust_kernel::radial_tensor::Error) -> GeometryError {
    error(
        "pipe hardware retained boundary",
        match cause {
            multicad_rust_kernel::radial_tensor::Error::Invalid => GeometryErrorKind::Input,
            multicad_rust_kernel::radial_tensor::Error::Limit => GeometryErrorKind::Limit,
            multicad_rust_kernel::radial_tensor::Error::Cancelled => GeometryErrorKind::Cancelled,
            _ => GeometryErrorKind::Unsupported,
        },
        format!("Actual retained thread proof failed: {cause:?}"),
    )
}
fn interval(b: [f64; 2]) -> I {
    I::new(b[0], b[1])
}
fn midpoint(p: [I; 3]) -> [f64; 3] {
    p.map(|v| v.lo.midpoint(v.hi))
}
fn norm(p: [I; 3]) -> I {
    let sum = p.into_iter().fold(I::ZERO, |s, v| s + v.square());
    // Directed addition can put the lower endpoint one ulp below zero.
    // The exact sum of these three squares is nonnegative on the whole box.
    I::new(sum.lo.max(0.), sum.hi).sqrt()
}
fn quality(q: Tessellation) -> Result<(), GeometryError> {
    if !q.deflection.is_finite()
        || q.deflection <= 0.
        || !q.angle.is_finite()
        || q.angle <= 0.
        || q.angle > std::f64::consts::PI
        || q.triangle_limit == 0
        || q.triangle_limit > 4_000_000
    {
        Err(failure(Cause::Invalid))
    } else {
        Ok(())
    }
}
impl Solid {
    /// Head construction uses the ORIGINAL neck frame; the exact internal
    /// [X,-Y,-Z],z=L factor couples its tip and neck without a rounded datum.
    /// A through nut uses the original mouth frame and unchanged negative taper.
    pub fn rust_authored_axial_pipe_hardware_with_interrupt(
        spec: AuthoredSpec,
        kind: PipeHardwareKind,
        minimum_wall_mm: f64,
        feature: EntityId,
        frame: PlaneFrame,
        policy: GeometryPolicy,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, GeometryError> {
        construct(
            spec,
            kind,
            minimum_wall_mm,
            feature,
            frame,
            policy,
            interrupt,
            None,
        )
        .map(|(solid, _)| solid)
    }
    pub(crate) fn pipe_hardware_body(&self) -> Option<&PipeHardwareBody> {
        match &self.body {
            SolidBody::PipeHardware(body) => Some(body),
            _ => None,
        }
    }
    /// Original source intent only; actual retained-source measures are used
    /// by this solid's volume/faces and manufacturing paths.
    pub fn rust_pipe_hardware_intent(&self) -> Option<&Material> {
        self.pipe_hardware_body().map(|b| &b.source.intent)
    }
    pub fn rust_pipe_hardware_mesh_build_count(&self) -> Option<usize> {
        self.pipe_hardware_body().map(|b| b.mesh_builds.get())
    }
}
#[allow(clippy::too_many_arguments)]
fn construct(
    spec: AuthoredSpec,
    kind: PipeHardwareKind,
    minimum_wall_mm: f64,
    feature: EntityId,
    frame: PlaneFrame,
    policy: GeometryPolicy,
    interrupt: Option<&dyn Fn() -> bool>,
    prepared: Option<Rc<ThreadSource>>,
) -> Result<(Solid, Rc<ThreadSource>), GeometryError> {
    policy.validate()?;
    frame.checked(policy)?;
    checkpoint(interrupt)?;
    if feature.is_nil()
        || !minimum_wall_mm.is_finite()
        || minimum_wall_mm < 0.
        || policy.boolean_fuzzy_mm != 0.
    {
        return Err(failure(Cause::Invalid));
    }
    let authored = AuthoredThread::new(spec, policy.linear_tolerance_mm, policy.minimum_feature_mm)
        .map_err(failure)?;
    let intent =
        Material::new(authored, kind, MaterialOptions::default(), interrupt).map_err(failure)?;
    let prepared = if let Some(prepared) = prepared {
        prepared
    } else {
        let thread = Solid::rust_authored_axial_pipe_thread_with_interrupt(
            spec,
            feature,
            PlaneFrame::default(),
            policy,
            interrupt,
        )?;
        let across = match kind {
            PipeHardwareKind::HexHead {
                across_flats_mm, ..
            }
            | PipeHardwareKind::HexNut { across_flats_mm } => across_flats_mm,
        };
        source::prepare(
            thread,
            intent.thread(),
            across,
            minimum_wall_mm.max(policy.minimum_feature_mm),
            feature,
            policy,
            interrupt,
        )?
    };
    let source = source::build(
        Rc::clone(&prepared),
        intent,
        minimum_wall_mm.max(policy.minimum_feature_mm),
        feature,
        policy,
        interrupt,
    )?;
    let placement = Placement::new(frame, kind, spec.length_mm, policy)?;
    let body = PipeHardwareBody {
        feature,
        source: Rc::new(source),
        placement,
        display: RefCell::new(None),
        mesh_builds: std::cell::Cell::new(0),
    };
    body.validate(policy)?;
    checkpoint(interrupt)?;
    let solid = Solid {
        body: SolidBody::PipeHardware(Rc::new(body)),
        policy,
        provenance: vec![],
        rust_origins: vec![],
        retained_placement: None,
    };
    Ok((solid, prepared))
}
impl PipeHardwareBody {
    pub(super) fn transformed(
        &self,
        t: RigidTransform,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        let result = Self {
            feature: self.feature,
            source: Rc::clone(&self.source),
            placement: self.placement.transformed(t, policy)?,
            display: RefCell::new(None),
            mesh_builds: std::cell::Cell::new(0),
        };
        result.validate(policy)?;
        Ok(result)
    }
    fn validate(&self, policy: GeometryPolicy) -> Result<(), GeometryError> {
        if let PipeHardwareKind::HexHead { height_mm, .. } = self.source.intent.kind()
            && (I::point(height_mm) * I::point(self.placement.sigma[0])).lo
                <= policy.minimum_feature_mm
        {
            return Err(failure(Cause::Unsupported));
        }
        if (self.source.clearance * I::point(self.placement.sigma[0])).lo
            <= self.source.minimum_wall
        {
            return Err(failure(Cause::Unsupported));
        }
        self.curved_centroid_shift(policy)?;
        let bounds = self.bounds()?;
        policy.point(bounds.min)?;
        policy.point(bounds.max)?;
        // Complete physical bounds enclose every source point and therefore
        // its area centroid. Descriptor integration is deferred until queried.
        Ok(())
    }
    pub(super) fn bounds(&self) -> Result<Bounds, GeometryError> {
        // A linear map reaches each hex-prism support at an actual vertex.
        // Keep the long thread's box separate from the head stock, avoiding
        // a fictitious full-length extrusion of the larger head radius.
        let mut b = [I::new(f64::INFINITY, f64::NEG_INFINITY); 3];
        for p in self.source.hex_vertices {
            let world = self.placement.point_bounds(p)?;
            for k in 0..3 {
                b[k] = I::new(b[k].lo.min(world[k].lo), b[k].hi.max(world[k].hi));
            }
        }
        // Retained cap sewing collars also enclose their SAME literal rim
        // vertices. Keep that actual source enclosure for both roles rather
        // than assuming every declared bore endpoint has exactly zero lift.
        {
            let original = self.source.thread.bounds()?;
            let world = self.placement.point_bounds(std::array::from_fn(|k| {
                I::new(original.min[k], original.max[k])
            }))?;
            for k in 0..3 {
                b[k] = I::new(b[k].lo.min(world[k].lo), b[k].hi.max(world[k].hi));
            }
        }
        Ok(Bounds {
            min: b.map(|p| p.lo),
            max: b.map(|p| p.hi),
        })
    }
    pub(super) fn volume_mm3(&self) -> Result<f64, GeometryError> {
        let v = self.source.volume * self.placement.determinant;
        if !v.valid() || v.lo <= 0. {
            return Err(failure(Cause::Unresolved));
        }
        Ok(v.lo.midpoint(v.hi))
    }
    fn curved_centroid_shift(&self, policy: GeometryPolicy) -> Result<I, GeometryError> {
        let span = self
            .source
            .bounds
            .map(|b| I::point(b[1]) - I::point(b[0]))
            .into_iter()
            .fold(I::ZERO, |s, v| s + v);
        let scale = self.placement.area_scale();
        let delta = scale / I::point(scale.lo) - I::ONE;
        let shift = I::point(self.placement.sigma[1]) * span * delta / (I::ONE - delta);
        if self.source.faces.iter().any(|f| f.frame.is_none())
            && (!shift.valid() || shift.hi > policy.linear_tolerance_mm / 8.)
        {
            return Err(failure(Cause::Unresolved));
        }
        Ok(shift)
    }
    pub(super) fn faces(&self, policy: GeometryPolicy) -> Result<Vec<Face>, GeometryError> {
        let shift = self.curved_centroid_shift(policy)?;
        let scale = self.placement.area_scale();
        self.source
            .faces
            .iter()
            .enumerate()
            .map(|(index, f)| {
                let (local_area, local_center) = match f.measure {
                    FaceMeasure::Ready { area, center } => (area, center),
                    FaceMeasure::Retained(old) => self.source._prepared.measure(old, None)?,
                };
                let (area, center, plane) = if let Some(frame) = f.frame {
                    let n = f.normal.expect("sealed actual planar source normal");
                    let factor = self.placement.normal_factor(n)?;
                    let normal = self.placement.normal(n)?;
                    let origin = self.placement.point_bounds(frame.origin_mm.map(I::point))?;
                    (
                        local_area * factor,
                        self.placement.point_bounds(local_center)?,
                        Some((midpoint(origin), midpoint(normal))),
                    )
                } else {
                    let center = self
                        .placement
                        .point_bounds(local_center)?
                        .map(|v| v + I::new(-shift.hi, shift.hi));
                    (local_area * scale, center, None)
                };
                if !area.valid() || area.lo <= 0. || center.iter().any(|v| !v.valid()) {
                    return Err(failure(Cause::Unresolved));
                }
                policy.point(midpoint(center))?;
                Ok(Face {
                    index: index as u32,
                    plane,
                    center_mm: midpoint(center),
                    area_mm2: area.lo.midpoint(area.hi),
                })
            })
            .collect()
    }
}
