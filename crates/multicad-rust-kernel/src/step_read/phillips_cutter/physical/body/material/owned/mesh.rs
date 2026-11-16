//! Source-certified joined facet candidate. Global literal embedding is a
//! separate consumer gate, not inferred from indexed manifold connectivity.
use super::*;
use crate::display::{DisplayMesh, Tessellation};
mod delegated;
mod geometry;
pub use delegated::DelegatedWorkspaceBounds;
#[cfg(test)]
pub(super) mod tests;
use geometry::*;
const MESH_WORK: usize = 512_000_000;
const MESH_BYTES: usize = 64 * 1024 * 1024;
// All fixed local headers and bounded degree-nine interval scratch. Delegated
// provider workspaces are separate; this accounts NEW producer-owned storage.
const MESH_SCRATCH: usize = 65536;
#[derive(Clone, Copy)]
pub struct JoinedMeshOptions {
    pub tessellation: Tessellation,
    pub maximum_work: usize,
    pub maximum_bytes: usize,
}
impl Default for JoinedMeshOptions {
    fn default() -> Self {
        Self {
            tessellation: Tessellation {
                deflection: 0.1,
                angle: 0.3,
                triangle_limit: 1_000_000,
            },
            maximum_work: MESH_WORK,
            maximum_bytes: MESH_BYTES,
        }
    }
}
/// Bound to the actual owned repaired source. Source distance, actual normals
/// and closed indexed manifold are certified; no global embedding or Native
/// manufacturing publication is earned by this value alone.
#[derive(Clone)]
pub struct JoinedMaterialCandidate<'s> {
    source: &'s OwnedDerivedMaterial,
    mesh: DisplayMesh,
    connectivity: crate::manufacturing::MeshProof,
    distance: f64,
    quality: Tessellation,
    work: usize,
    remaining: usize,
    peak: usize,
    delegated: DelegatedWorkspaceBounds,
    combined_peak: usize,
}
impl<'s> JoinedMaterialCandidate<'s> {
    pub fn delegated_workspace_bounds(&self) -> DelegatedWorkspaceBounds {
        self.delegated
    }
    /// Conservative sequential sum of the accounted peak plus the largest
    /// delegated bound. Components do not run concurrently across stages.
    pub fn combined_live_upper_bound_bytes(&self) -> usize {
        self.combined_peak
    }
    pub fn source(&self) -> &'s OwnedDerivedMaterial {
        self.source
    }
    pub fn mesh(&self) -> &DisplayMesh {
        &self.mesh
    }
    pub fn connectivity(&self) -> crate::manufacturing::MeshProof {
        self.connectivity
    }
    /// Maximum computed transfer bound for NEW head, joins and translated
    /// cavity. The unchanged retained-thread provider separately earns the
    /// requested deflection through its mandatory whole-source certificates.
    pub fn new_surface_error_mm(&self) -> f64 {
        self.distance
    }
    pub fn certified_tessellation(&self) -> Tessellation {
        self.quality
    }
    /// Balance that a later whole-candidate embedding gate must consume;
    /// this producer never renews the caller's query budget.
    pub fn remaining_work(&self) -> usize {
        self.remaining
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    /// Accounted owner/live-component/new-buffer peak, including fixed scratch.
    /// Existing delegated mesher/CDT temporaries retain separate bounds. This
    /// is neither whole-process RSS nor allocator-fragmentation accounting.
    pub fn accounted_peak_bytes(&self) -> usize {
        self.peak
    }
}
#[derive(Clone)]
struct Work<'a> {
    remaining: usize,
    held: usize,
    peak: usize,
    maximum_bytes: usize,
    maximum_triangles: usize,
    delegated: DelegatedWorkspaceBounds,
    interrupt: Option<&'a dyn Fn() -> bool>,
}
impl Work<'_> {
    fn charge(&mut self, n: usize) -> R<()> {
        if self.interrupt.is_some_and(|f| f()) {
            return Err(PhillipsReadError::Cancelled);
        }
        self.remaining = self
            .remaining
            .checked_sub(n)
            .ok_or(PhillipsReadError::Limit)?;
        Ok(())
    }
    fn storage(&mut self, additional: usize) -> R<()> {
        self.charge(1)?;
        let bytes = self
            .held
            .saturating_add(additional)
            .saturating_add(MESH_SCRATCH);
        if bytes > self.maximum_bytes {
            return Err(PhillipsReadError::Limit);
        }
        self.peak = self.peak.max(bytes);
        Ok(())
    }
}
fn mesh_bytes(mesh: &DisplayMesh) -> usize {
    mesh.vertices.capacity() * 8 + mesh.indices.capacity() * 4 + mesh.face_indices.capacity() * 4
}
fn radial_error(e: crate::tapered_thread::Error) -> PhillipsReadError {
    match e {
        crate::tapered_thread::Error::Cancelled => PhillipsReadError::Cancelled,
        crate::tapered_thread::Error::Limit => PhillipsReadError::Limit,
        _ => PhillipsReadError::Unsupported,
    }
}

fn unresolved_mesh(_: crate::display::MeshError) -> PhillipsReadError {
    PhillipsReadError::Unsupported
}
fn append_mesh(
    target: &mut crate::display::DisplayMesh,
    mesh: &crate::display::DisplayMesh,
    face: &mut dyn FnMut(u32) -> Option<u32>,
    z: f64,
    work: &mut Work<'_>,
) -> R<u32> {
    let old = mesh_bytes(target);
    let v = target.vertices.len() + mesh.vertices.len();
    let i = target.indices.len() + mesh.indices.len();
    let f = target.face_indices.len() + mesh.face_indices.len();
    work.storage(old + v * 8 + (i + f) * 4)?;
    target.vertices.reserve_exact(v - target.vertices.len());
    target.indices.reserve_exact(i - target.indices.len());
    target
        .face_indices
        .reserve_exact(f - target.face_indices.len());
    work.storage(old + mesh_bytes(target))?;
    let offset = (target.vertices.len() / 3) as u32;
    for p in mesh.vertices.as_chunks::<3>().0 {
        work.charge(1)?;
        target.vertices.extend([p[0], p[1], p[2] + z]);
    }
    for (t, &f) in mesh
        .indices
        .as_chunks::<3>()
        .0
        .iter()
        .zip(&mesh.face_indices)
    {
        work.charge(1)?;
        if let Some(f) = face(f) {
            if target.face_indices.len() == work.maximum_triangles {
                return Err(PhillipsReadError::Limit);
            }
            target.indices.extend(t.map(|i| i + offset));
            target.face_indices.push(f);
        }
    }
    Ok(offset)
}
struct PlaneProof {
    face: u32,
    positive: bool,
    source_z: f64,
    boundary_error: f64,
    options: Tessellation,
}
fn plane(
    target: &mut crate::display::DisplayMesh,
    outer: &[u32],
    inner: &[u32],
    inner_forward: bool,
    proof: PlaneProof,
    work: &mut Work<'_>,
) -> R<f64> {
    let PlaneProof {
        face,
        positive,
        source_z,
        boundary_error,
        options,
    } = proof;
    let n = outer.len() + inner.len();
    work.delegated.planar_bytes = work.delegated.planar_bytes.max(delegated::planar(n)?);
    work.charge(n)?;
    if target
        .face_indices
        .len()
        .checked_add(n)
        .ok_or(PhillipsReadError::Limit)?
        > work.maximum_triangles
    {
        return Err(PhillipsReadError::Limit);
    }
    // New vectors plus the existing output, before allocation. CDT internal
    // workspace is a separately delegated provider bound.
    work.storage(mesh_bytes(target) + n * 64)?;
    let vertices = outer.iter().chain(inner).copied().collect::<Vec<_>>();
    let points = vertices
        .iter()
        .map(|&id| {
            spade::Point2::new(
                target.vertices[3 * id as usize],
                target.vertices[3 * id as usize + 1],
            )
        })
        .collect::<Vec<_>>();
    let mut edges = (0..outer.len())
        .map(|i| [i, (i + 1) % outer.len()])
        .collect::<Vec<_>>();
    for i in 0..inner.len() {
        let a = outer.len() + i;
        let b = outer.len() + (i + 1) % inner.len();
        edges.push(if inner_forward { [a, b] } else { [b, a] });
    }
    let region = crate::display::triangulate_region(points, &edges, 1).map_err(unresolved_mesh)?;
    let old = mesh_bytes(target);
    work.storage(
        old + n * 64
            + (target.indices.len() + 3 * region.triangles.len()) * 4
            + (target.face_indices.len() + region.triangles.len()) * 4,
    )?;
    target.indices.reserve_exact(3 * region.triangles.len());
    target.face_indices.reserve_exact(region.triangles.len());
    work.storage(old + mesh_bytes(target) + n * 64)?;
    let mut distance = 0f64;
    for t in region.triangles {
        work.charge(16)?;
        let mut ids = t.map(|i| vertices[i]);
        if !positive {
            ids.swap(1, 2);
        }
        let points: [P; 3] =
            ids.map(|id| std::array::from_fn(|k| target.vertices[3 * id as usize + k]));
        let height = points
            .iter()
            .map(|p| (I::point(p[2]) - I::point(source_z)).magnitude())
            .fold(0., f64::max);
        // Outside the actual straight paired-loop homotopy, winding is
        // unchanged. Every point in that homotopy is within boundary_error
        // of BOTH corresponding loop endpoints. Thus the two planar regions
        // are two-sided close; literal triangle lift adds at most height.
        let error = (I::point(boundary_error) + I::point(height)).hi;
        distance = distance.max(error);
        if error > options.deflection {
            return Err(PhillipsReadError::Unsupported);
        }
        let source_normal = [
            I::point(0.),
            I::point(0.),
            I::point(if positive { 1. } else { -1. }),
        ];
        if !normal_match(points, source_normal, options.angle) {
            return Err(PhillipsReadError::Unsupported);
        }
        target.indices.extend(ids);
        target.face_indices.push(face);
    }
    Ok(distance)
}
fn compact(mesh: &DisplayMesh, work: &mut Work<'_>) -> R<DisplayMesh> {
    work.storage(2 * mesh_bytes(mesh) + mesh.vertices.len() / 3 * 4)?;
    let mut remap = vec![u32::MAX; mesh.vertices.len() / 3];
    let mut out = crate::display::DisplayMesh {
        vertices: Vec::with_capacity(mesh.vertices.len()),
        indices: Vec::with_capacity(mesh.indices.len()),
        face_indices: mesh.face_indices.clone(),
    };
    work.storage(mesh_bytes(mesh) + mesh_bytes(&out) + remap.capacity() * 4)?;
    for &id in &mesh.indices {
        work.charge(1)?;
        let mapped = &mut remap[id as usize];
        if *mapped == u32::MAX {
            *mapped = (out.vertices.len() / 3) as u32;
            out.vertices
                .extend(&mesh.vertices[3 * id as usize..3 * id as usize + 3]);
        }
        out.indices.push(*mapped);
    }
    Ok(out)
}

// RadialBoundary::certify_caps runs caps::positive_spline/positive_curve:
// every degree17 Bernstein coefficient of x*y'-y*x' is strictly positive,
// including each refined source span, and the complete rim has one wrap.
// graph::axis gives the owned derived source literal common XY and±Z axes;
// the checks below also keep this helper honest for other retained sources.
// The original edge/coedge reversal is copied exactly from the same trace.
// Therefore every returned/restricted Angular follows increasing true angle;
// inside its certified endpoint-bisector halfplane it follows the MINOR arc
// between those endpoints, not an arbitrary excursion within the halfplane.
fn angular_from_trace<'s>(
    source: &'s RadialBoundary,
    trace: &crate::radial_tensor::TracedMesh<'s>,
    work: &mut Work<'_>,
) -> R<Vec<Angular<'s>>> {
    if !std::ptr::eq(source, trace.source()) {
        return Err(PhillipsReadError::Unsupported);
    }
    let axis = source.cap_axes()[0];
    if axis.origin[0] != 0.
        || axis.origin[1] != 0.
        || axis.reference != [1., 0., 0.]
        || axis.normal[0] != 0.
        || axis.normal[1] != 0.
        || axis.normal[2].abs() != 1.
    {
        return Err(PhillipsReadError::Unsupported);
    }
    let ring = trace.cap_rims()[0];
    work.storage(ring.len().max(4) * 2 * std::mem::size_of::<Angular<'_>>())?;
    let mut angular = Vec::with_capacity(ring.len());
    work.storage(angular.capacity() * std::mem::size_of::<Angular<'_>>())?;
    for (i, p) in ring.iter().enumerate() {
        work.charge(1)?;
        let edge = source.input().edges[p.edge];
        let crate::radial_tensor::Curve::Spline { spans, poles } =
            &source.input().curves[edge.curve]
        else {
            return Err(PhillipsReadError::Unsupported);
        };
        angular.push(Angular {
            poles,
            spans: *spans,
            reverse: edge.forward != p.edge_forward,
            range: [
                p.uv[0],
                if ring[(i + 1) % ring.len()].side == p.side {
                    ring[(i + 1) % ring.len()].uv[0]
                } else {
                    1.
                },
            ],
        });
    }
    Ok(angular)
}
// A failed angular proposal doubles only this new head grid. It never
// changes the retained thread/cavity grids or accepts any uncertified facet.
fn refine_angular(angular: &mut Vec<Angular<'_>>, work: &mut Work<'_>) -> R<()> {
    let old = angular.capacity() * std::mem::size_of::<Angular<'_>>();
    work.storage(2 * angular.len() * std::mem::size_of::<Angular<'_>>())?;
    let mut children = Vec::with_capacity(2 * angular.len());
    work.storage(children.capacity() * std::mem::size_of::<Angular<'_>>())?;
    for &a in angular.iter() {
        work.charge(1)?;
        let m = a.range[0] + (a.range[1] - a.range[0]) / 2.;
        if m <= a.range[0] || m >= a.range[1] {
            return Err(PhillipsReadError::Unsupported);
        }
        children.push(Angular {
            range: [a.range[0], m],
            ..a
        });
        children.push(Angular {
            range: [m, a.range[1]],
            ..a
        });
    }
    work.held = work.held - old + children.capacity() * std::mem::size_of::<Angular<'_>>();
    work.storage(old)?;
    *angular = children;
    Ok(())
}
// The real bounded state machine is shared by normal assembly and coarse
// source proposals; neither path bypasses any whole-source facet proof.
fn make_head(
    source: &OwnedDerivedMaterial,
    angular: &mut Vec<Angular<'_>>,
    options: Tessellation,
    work: &mut Work<'_>,
) -> R<Head> {
    let mut counts = vec![1; source.head().view().segment_count()];
    work.held += counts.capacity() * std::mem::size_of::<usize>();
    work.storage(0)?;
    let mut head = None;
    for _ in 0..64 {
        match head_attempt(source, angular, &counts, options, work)? {
            Ok(h) => {
                head = Some(h);
                break;
            }
            Err(split) => match split {
                Split::T(p) => counts[p] *= 2,
                Split::U => refine_angular(angular, work)?,
            },
        }
    }
    head.ok_or(PhillipsReadError::Limit)
}
#[derive(Clone, Copy)]
struct SourceQuality {
    requested: Tessellation,
    source_angle: f64,
    angular_reserve: f64,
}
fn source_quality(options: Tessellation) -> R<SourceQuality> {
    // Derived lower choices keep source-angle + transfer-angle <= the
    // exact caller angle, including product/subtraction rounding. The
    // actual measured distance is checked independently for both joins.
    let source_angle = (I::point(options.angle) * I::point(0.99)).lo;
    let angular_reserve = (I::point(options.angle) - I::point(source_angle)).lo;
    if source_angle < 0.01 {
        return Err(PhillipsReadError::Unsupported);
    }
    Ok(SourceQuality {
        requested: options,
        source_angle,
        angular_reserve,
    })
}
// Thread receipt may be reused only for the same immutable source and exact
// earned quality. Cavity/atlas are freshly earned by this private constructor.
// All references point into the external OwnedDerivedMaterial, never into this
// value; no self-referential owner or guessed geometry is introduced.
struct Components<'s, 't> {
    source: &'s OwnedDerivedMaterial,
    thread: &'t crate::radial_tensor::TracedMesh<'s>,
    atlas: crate::phillips_drive::PatchAtlas<'s>,
    cavity: crate::phillips_drive::TracedMouthMesh<'s>,
    quality: SourceQuality,
}
fn components<'s, 't>(
    source: &'s OwnedDerivedMaterial,
    thread: &'t crate::radial_tensor::TracedMesh<'s>,
    quality: SourceQuality,
    work: &mut Work<'_>,
) -> R<Components<'s, 't>> {
    let options = quality.requested;
    let actual = thread.tessellation();
    if !std::ptr::eq(source.thread(), thread.source())
        || actual.deflection.to_bits() != options.deflection.to_bits()
        || actual.angle.to_bits() != options.angle.to_bits()
        || actual.triangle_limit != options.triangle_limit
    {
        return Err(PhillipsReadError::Unsupported);
    }
    let interrupt = work.interrupt;
    // The sealed Cutter schema has at most88 patch descriptors,49 vertices,
    //76 edges and29 faces. One MiB conservatively preflights all their Vec
    // growth plus construction temporaries; the actual retained payload is
    // inventoried below, including the nested mouth SectionBoundary.
    work.storage(1024 * 1024)?;
    work.charge(4096)?;
    let atlas = source.cavity().surface_patches().map_err(cavity_error)?;
    work.held += atlas.retained_budget_bytes();
    work.storage(0)?;
    work.delegated.cutter_bytes = delegated::cutter(options.triangle_limit)?;
    let cavity = atlas
        .boundary_mesh_with_mouth_trace(
            Tessellation {
                deflection: options.deflection * 0.99,
                angle: quality.source_angle,
                ..options
            },
            crate::phillips_drive::MouthTraceOptions {
                maximum_work: work.remaining,
                ..Default::default()
            },
            &|| interrupt.is_some_and(|f| f()),
        )
        .map_err(cavity_error)?;
    work.charge(cavity.work_used())?;
    work.held += mesh_bytes(cavity.boundary_mesh().mesh()) + cavity.retained_metadata_bytes();
    work.storage(0)?;
    Ok(Components {
        source,
        thread,
        atlas,
        cavity,
        quality,
    })
}
fn assemble<'s>(
    parts: &Components<'s, '_>,
    maximum_work: usize,
    mut work: Work<'_>,
) -> R<JoinedMaterialCandidate<'s>> {
    let source = parts.source;
    let options = parts.quality.requested;
    let thread = parts.thread;
    let atlas = &parts.atlas;
    let cavity = &parts.cavity;
    // Actual immutable start-cap ordering is owned by the trace; this
    // derived graph's start face is0, as earned by graph::build/material.
    let ring = thread.cap_rims()[0];
    let mut angular = angular_from_trace(source.thread(), thread, &mut work)?;
    work.held += angular.capacity() * std::mem::size_of::<Angular<'_>>();
    work.storage(0)?;
    let head = make_head(source, &mut angular, options, &mut work)?;
    work.held += mesh_bytes(&head.mesh) + (head.top.capacity() + head.bottom.capacity()) * 4;
    work.storage(0)?;
    let h = source.head().view().segment_count();
    let sides = source.thread().side_faces();
    let mut mesh = crate::display::DisplayMesh {
        vertices: Vec::new(),
        indices: Vec::new(),
        face_indices: Vec::new(),
    };
    let thread_offset = append_mesh(
        &mut mesh,
        thread.mesh(),
        &mut |f| {
            if f == 0 {
                None
            } else if f == 1 {
                Some((h + sides.len()) as u32)
            } else {
                Some((h + sides.iter().position(|&i| i == f as usize).unwrap()) as u32)
            }
        },
        0.,
        &mut work,
    )?;
    let head_offset = append_mesh(&mut mesh, &head.mesh, &mut Some, 0., &mut work)?;
    let cavity_faces = atlas.topology().faces();
    let cavity_offset = append_mesh(
        &mut mesh,
        cavity.boundary_mesh().mesh(),
        &mut |f| {
            let key = cavity_faces[f as usize].source_identity();
            if key == CavityFace::Mouth {
                None
            } else {
                Some(
                    (h + sides.len()
                        + 3
                        + source.cavity_faces[..source.cavity_count]
                            .iter()
                            .position(|&k| k == key)
                            .unwrap()) as u32,
                )
            }
        },
        source.chart().outer_top_z,
        &mut work,
    )?;
    work.charge(
        cavity.boundary_mesh().mesh().vertices.len() + cavity.boundary_mesh().mesh().indices.len(),
    )?;
    let translated_error = translation_transfer(
        cavity.boundary_mesh().mesh(),
        &mesh,
        cavity_offset,
        source.chart().outer_top_z,
        options.deflection * 0.01,
        parts.quality.angular_reserve,
        &mut work,
    )?;
    // Reverse the actual Cutter side winding, retaining original source indices.
    for (tri, &f) in mesh
        .indices
        .as_chunks_mut::<3>()
        .0
        .iter_mut()
        .zip(&mesh.face_indices)
    {
        work.charge(1)?;
        if f as usize >= h + sides.len() + 3 {
            tri.swap(1, 2);
        }
    }
    work.storage(
        mesh_bytes(&mesh)
            + (head.bottom.len() + head.top.len() + ring.len() + cavity.mouth().len()) * 4,
    )?;
    let bottom = head
        .bottom
        .iter()
        .map(|&i| i + head_offset)
        .collect::<Vec<_>>();
    let top = head
        .top
        .iter()
        .map(|&i| i + head_offset)
        .collect::<Vec<_>>();
    let inner = ring
        .iter()
        .map(|p| p.vertex + thread_offset)
        .collect::<Vec<_>>();
    work.held += (bottom.capacity() + top.capacity() + inner.capacity()) * 4;
    let join_distance = plane(
        &mut mesh,
        &bottom,
        &inner,
        false,
        PlaneProof {
            face: (h + sides.len() + 1) as u32,
            positive: true,
            source_z: source.head().view().station(h).unwrap()[1],
            boundary_error: ring
                .iter()
                .map(|p| p.outgoing_chord_error_mm)
                .fold(head.distance, f64::max),
            options,
        },
        &mut work,
    )?;
    let mouth = cavity
        .mouth()
        .iter()
        .map(|p| p.vertex + cavity_offset)
        .collect::<Vec<_>>();
    work.held += mouth.capacity() * 4;
    let top_distance = plane(
        &mut mesh,
        &top,
        &mouth,
        true,
        PlaneProof {
            face: (h + sides.len() + 2) as u32,
            positive: false,
            source_z: source.chart().outer_top_z,
            boundary_error: (I::point(cavity.boundary_mesh().source_distance_mm())
                + I::point(translated_error))
            .hi
            .max(head.distance),
            options,
        },
        &mut work,
    )?;
    let mesh = compact(&mesh, &mut work)?;
    work.charge(mesh.vertices.len() + mesh.indices.len())?;
    work.delegated.connectivity_bytes =
        delegated::connectivity(mesh.vertices.len() / 3, mesh.face_indices.len())?;
    let proof = crate::manufacturing::certify(&mesh).map_err(unresolved_mesh)?;
    let result = JoinedMaterialCandidate {
        source,
        mesh,
        connectivity: proof,
        quality: options,
        distance: head
            .distance
            .max(join_distance)
            .max(top_distance)
            .max(
                (I::point(cavity.boundary_mesh().source_distance_mm())
                    + I::point(translated_error))
                .hi,
            )
            .max(
                ring.iter()
                    .map(|p| p.outgoing_chord_error_mm)
                    .fold(0., f64::max),
            ),
        work: maximum_work - work.remaining,
        remaining: work.remaining,
        peak: work.peak,
        delegated: work.delegated,
        combined_peak: work
            .peak
            .checked_add(work.delegated.maximum())
            .ok_or(PhillipsReadError::Limit)?,
    };
    publish(result, maximum_work, work)
}
impl OwnedDerivedMaterial {
    pub fn joined_mesh_candidate(
        &self,
        query: JoinedMeshOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> R<JoinedMaterialCandidate<'_>> {
        let options = query.tessellation;
        options.check().map_err(|_| PhillipsReadError::Invalid)?;
        if query.maximum_work > MESH_WORK || query.maximum_bytes > MESH_BYTES {
            return Err(PhillipsReadError::Invalid);
        }
        let quality = source_quality(options)?;
        let source = self;
        let mut work = Work {
            remaining: query.maximum_work,
            held: self.retained_budget_bytes(),
            peak: 0,
            maximum_bytes: query.maximum_bytes,
            maximum_triangles: options.triangle_limit,
            delegated: Default::default(),
            interrupt,
        };
        work.storage(0)?;
        work.charge(source.thread().input().surfaces.len())?;
        work.delegated.retained_thread_bytes =
            delegated::retained(source.thread(), options.triangle_limit)?;
        let thread = source
            .thread()
            .tessellate_with_cap_trace_budget(
                options,
                Default::default(),
                work.remaining,
                work.interrupt,
            )
            .map_err(radial_error)?;
        work.charge(thread.work_used())?;
        work.held += mesh_bytes(thread.mesh()) + thread.retained_metadata_bytes();
        work.storage(0)?;
        let parts = components(source, &thread, quality, &mut work)?;
        assemble(&parts, query.maximum_work, work)
    }
}

// Actual source-certified geometry has already been assembled and checked.
// This final phase independently checks its live capacities and cancellation
// before any candidate escapes; failed publication cannot cache a receipt.
fn publish<'s>(
    mut result: JoinedMaterialCandidate<'s>,
    maximum_work: usize,
    mut work: Work<'_>,
) -> R<JoinedMaterialCandidate<'s>> {
    work.storage(mesh_bytes(&result.mesh))?;
    if result.mesh.face_indices.len() > work.maximum_triangles {
        return Err(PhillipsReadError::Limit);
    }
    work.charge(1)?;
    result.work = maximum_work - work.remaining;
    result.remaining = work.remaining;
    result.peak = work.peak;
    result.delegated = work.delegated;
    result.combined_peak = work
        .peak
        .checked_add(work.delegated.maximum())
        .ok_or(PhillipsReadError::Limit)?;
    Ok(result)
}
