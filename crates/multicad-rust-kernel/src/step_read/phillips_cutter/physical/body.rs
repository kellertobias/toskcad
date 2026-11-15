//! Same-file outer/thread headers and literal raw arena. This packet is NOT a
//! Solid, decoded mass, verified tensor law, cavity containment or virtual-body
//! publication receipt. No expected authoring source or gauge is accepted.
use super::*;
use crate::step_read::helical_import::{law::Law, seed::TensorHeader};
mod outer;
mod retained;
pub use outer::{OuterCarrierBounds, OuterCarrierOptions, TensorCarrierBound};
pub use retained::{ExactPointLeaf, OriginalNumericLeaf, OriginalNumericLeaves};
#[cfg(test)]
mod tests;

/// A virtual filled outer cap, expressly distinguished from the unchanged
/// actual material face and its retained inner hole. No raw entity is changed.
#[derive(Clone, Debug)]
pub struct VirtualOuterCap {
    pub material_face: u32,
    pub surface: u32,
    pub retained_outer_bound: u32,
    pub removed_inner_bound: u32,
    pub removed_cavity_faces: Vec<u32>,
}
#[derive(Debug)]
pub struct Body<'input> {
    input: &'input str,
    entities: BTreeMap<u32, Vec<Entity<'input>>>,
    material_faces: Vec<u32>,
    cavity: Boundary,
    thread: crate::helical_thread::Thread,
    frame: [Point; 3],
    visible: [f64; 2],
    profile: crate::head_profile::Profile,
    profile_faces: Vec<u32>,
    cap: VirtualOuterCap,
    tensor_face_count: usize,
    tolerance: f64,
    retained_budget_bytes: usize,
}
impl<'input> Body<'input> {
    pub(super) fn source_entities(&self) -> &BTreeMap<u32, Vec<Entity<'input>>> {
        &self.entities
    }
    pub fn cavity(&self) -> &Boundary {
        &self.cavity
    }
    pub fn literal_input(&self) -> &'input str {
        self.input
    }
    pub fn literal_entity_count(&self) -> usize {
        self.entities.len()
    }
    pub fn material_face_entities(&self) -> &[u32] {
        &self.material_faces
    }
    pub fn tensor_face_entities(&self) -> &[u32] {
        &self.material_faces[..self.tensor_face_count]
    }
    /// Inferred law/header only; actual literal nets remain authoritative and
    /// require their independent whole-support and finite-trim proof.
    pub fn inferred_thread_spec(&self) -> crate::helical_thread::Spec {
        self.thread.spec()
    }
    pub fn inferred_thread_frame(&self) -> [Point; 3] {
        self.frame
    }
    pub fn inferred_visible_interval_mm(&self) -> [f64; 2] {
        self.visible
    }
    pub fn retained_meridian_header(&self) -> &crate::head_profile::Profile {
        &self.profile
    }
    pub fn profile_face_entities(&self) -> &[u32] {
        &self.profile_faces
    }
    pub fn virtual_outer_cap(&self) -> &VirtualOuterCap {
        &self.cap
    }
    pub fn literal_tensor_control_points(
        &self,
        index: usize,
    ) -> std::result::Result<Vec<Vec<Point>>, PhillipsReadError> {
        self.literal_tensor_control_points_raw(index)
            .map_err(PhillipsReadError::from)
    }
    pub(super) fn literal_tensor_control_points_raw(
        &self,
        index: usize,
    ) -> Result<Vec<Vec<Point>>> {
        let id = *self
            .tensor_face_entities()
            .get(index)
            .ok_or(ReadError::Invalid)?;
        let shell = Shell::new(&self.entities, self.tolerance, true);
        let (surface, _, _) = helical::face(&shell, id)?;
        let args = shell.args(surface, "B_SPLINE_SURFACE_WITH_KNOTS")?;
        let Value::List(rows) = args.get(3).ok_or(ReadError::Invalid)? else {
            return Err(ReadError::Invalid);
        };
        rows.iter()
            .map(|row| refs(row)?.into_iter().map(|id| shell.point(id)).collect())
            .collect()
    }
}
fn plane(shell: &Shell<'_, '_>, id: u32) -> Result<Plane> {
    let (surface, forward, _) = helical::face(shell, id)?;
    require(forward)?;
    let [Value::String(_), axis] = shell.args(surface, "PLANE")? else {
        return Err(ReadError::Invalid);
    };
    shell.placement(reference(axis)?)
}
fn inferred_thread(
    shell: &Shell<'_, '_>,
    faces: &[u32],
) -> Result<(crate::helical_thread::Thread, [Point; 3], [f64; 2], usize)> {
    let mut count: usize = 0;
    for &id in faces {
        let (surface, _, _) = helical::face(shell, id)?;
        if shell.simple(surface)?.name != "B_SPLINE_SURFACE_WITH_KNOTS" {
            break;
        }
        count += 1;
    }
    if count < 4 || !count.is_multiple_of(4) {
        return Err(ReadError::Unsupported);
    }
    if count > 32 {
        return Err(ReadError::Limit);
    }
    let first = helical::face(shell, faces[0])?;
    let bottom = plane(shell, *faces.get(count).ok_or(ReadError::Invalid)?)?;
    let top = plane(shell, *faces.get(count + 1).ok_or(ReadError::Invalid)?)?;
    let tensor = TensorHeader::read(shell, faces)?;
    let (u, v, n) = (bottom.u, bottom.v, bottom.normal);
    let origin = tensor.retained_origin(shell, faces[0], bottom.origin, n)?;
    let height = dot(sub(top.origin, origin), n);
    require(height.is_finite() && height > shell.tolerance)?;
    require(
        length(sub(
            top.origin,
            std::array::from_fn(|i| origin[i] + n[i] * height),
        )) <= shell.tolerance / 16.,
    )?;
    let [start, spin] = tensor.axial;
    require(start == 0. && spin > 0.)?;
    let rate = height / spin;
    let starts = (count / 4) as u32;
    let pitch = std::f64::consts::TAU * rate / f64::from(starts);
    let major = TensorHeader::profiled_major(shell, faces[3], first.1)?;
    let upper = tensor.radius(shell, origin, 0, true)?;
    let lower = tensor.radius(shell, origin, 2, false)?;
    let thread = Law::infer(major, upper, lower, &tensor.radial, rate, pitch)?
        .thread(height, starts, first.1)?;
    Ok((
        thread,
        [origin, u, v],
        [dot(sub(bottom.origin, origin), n), height],
        count,
    ))
}
pub fn read_phillips_body(
    input: &str,
    tolerance: f64,
    work_limit: usize,
) -> std::result::Result<Body<'_>, PhillipsReadError> {
    read_phillips_body_with_interrupt(input, tolerance, work_limit, None)
}
pub fn read_phillips_body_with_interrupt<'input>(
    input: &'input str,
    tolerance: f64,
    work_limit: usize,
    interrupt: Option<&dyn Fn() -> bool>,
) -> std::result::Result<Body<'input>, PhillipsReadError> {
    let mut work = Work {
        remaining: work_limit,
        interrupt,
        cancelled: Cell::new(false),
    };
    let result =
        compound::read_members_retaining_arena(input, tolerance, true, true, |shell, id, _| {
            work.step(1)?;
            let faces = shell.member_faces(id)?;
            // The cavity extractor checks the ORIGINAL whole material incidence
            // before its cut, while this same Shell/arena retains all outer nets.
            let cavity = read_member(shell, &faces, &mut work)?;
            work.step(1)?;
            let (thread, frame, visible, count) = inferred_thread(shell, &faces)?;
            let outer = helical::profiled::read_phillips_outer_header(
                shell, &faces, &thread, frame, visible[0], &cavity,
            )?;
            work.step(1)?;
            Ok((faces, cavity, thread, frame, visible, count, outer))
        });
    if work.cancelled.get() {
        return Err(PhillipsReadError::Cancelled);
    }
    let (entities, mut members) = result.map_err(PhillipsReadError::from)?;
    let (faces, cavity, thread, frame, visible, count, outer) =
        members.pop().ok_or(PhillipsReadError::Invalid)?;
    let cap = VirtualOuterCap {
        material_face: outer.head_top_face,
        surface: outer.head_top_surface,
        retained_outer_bound: outer.outer_bound,
        removed_inner_bound: outer.removed_hole,
        removed_cavity_faces: cavity
            .faces()
            .iter()
            .filter_map(|f| f.material_entity)
            .collect(),
    };
    let mut body = Body {
        input,
        entities,
        material_faces: faces,
        cavity,
        thread,
        frame,
        visible,
        profile: outer.profile,
        profile_faces: outer.profile_faces,
        cap,
        tensor_face_count: count,
        tolerance,
        retained_budget_bytes: 0,
    };
    let budget = retained::budget(&body, &mut work);
    if work.cancelled.get() {
        return Err(PhillipsReadError::Cancelled);
    }
    body.retained_budget_bytes = budget.map_err(PhillipsReadError::from)?;
    Ok(body)
}

mod common_mouth;
pub use common_mouth::{
    CommonMouthChart, CommonMouthHeadEnvelope, CommonMouthOptions, CommonTensorTransfer,
};

mod closed_head;
pub use closed_head::{
    ClosedHeadOptions, ClosedRepairedHead, HeadSupportTransfer, RawHeadSupportKind,
};

mod derived_thread;
pub use derived_thread::{
    CoherentDegree11ThreadCarrier, ContactThreadCarrier, Degree11ContactOptions, Degree11DomainEnd,
    Degree11RimContact, Degree11RimContactId, Degree11RimKnot, Degree11SignDomain,
    Degree11ThreadContacts, Degree11ThreadFace, Degree11ThreadOptions, DerivedThreadCarrier,
    DerivedThreadOptions, DyadicParameter, OccupiedDerivedThread, RimContact, RimContactId,
    RimDomainEnd, RimDomainRole, RimKnot, RimSignDomain, ThreadCarrierTransfer,
    ThreadContactOptions,
};

mod material;

pub use material::{ClosedMaterialOptions, DerivedClosedMaterial, MaterialFaceSource};

pub use material::{MaterialVolume, MaterialVolumeOptions};

pub use material::{
    DelegatedWorkspaceBounds, JoinedMaterialCandidate, JoinedMeshOptions, OwnedDerivedMaterial,
    OwnedMaterialOptions, OwnedMaterialQuery, OwnedMaterialQueryError, OwnedMaterialQueryOptions,
};
