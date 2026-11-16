//! ACTUAL retained thread plus its ACTUAL finite head polygon. Global axial
//! sewing, literal rounding, cap chord and contact lift have separate bounds.
//! Indexed closure is earned here; whole literal embedding remains a separate
//! required consumer gate before any manufacturing capability is granted.
use super::*;
use crate::{
    display::{DisplayMesh, Tessellation},
    manufacturing::{self, MeshProof},
    tapered_thread::interval::I,
};
mod caps;
mod normal;
mod world;
pub use world::{MappedHeadMeshCandidate, MappedHeadMeshReceipt};
#[cfg(test)]
mod tests;
const HEAD_SCRATCH: usize = 64 * 1024;
const MAX_MESH_WORK: usize = 1_600_000_000;
#[derive(Clone, Copy)]
pub struct HeadMeshOptions {
    pub tessellation: Tessellation,
    pub maximum_work: usize,
    pub maximum_bytes: usize,
    pub additional_resident_bytes: usize,
}
impl Default for HeadMeshOptions {
    fn default() -> Self {
        Self {
            tessellation: Tessellation {
                deflection: 0.1,
                angle: 0.3,
                triangle_limit: 1_000_000,
            },
            maximum_work: MAX_MESH_WORK,
            maximum_bytes: MAX_BYTES,
            additional_resident_bytes: 0,
        }
    }
}
struct Work<'a> {
    used: usize,
    remaining: usize,
    base: usize,
    held: usize,
    peak: usize,
    maximum_bytes: usize,
    interrupt: Option<&'a dyn Fn() -> bool>,
}
impl Work<'_> {
    fn step(&mut self, n: usize) -> Result<()> {
        if self.interrupt.is_some_and(|f| f()) {
            return Err(Error::Cancelled);
        }
        self.remaining = self.remaining.checked_sub(n).ok_or(Error::Limit)?;
        self.used = self.used.checked_add(n).ok_or(Error::Limit)?;
        Ok(())
    }
    fn storage(&mut self, temporary: usize) -> Result<()> {
        let n = self
            .base
            .checked_add(self.held)
            .and_then(|n| n.checked_add(temporary))
            .ok_or(Error::Limit)?;
        if n > self.maximum_bytes {
            return Err(Error::Limit);
        }
        self.peak = self.peak.max(n);
        Ok(())
    }
    fn allocate<T>(&mut self, count: usize) -> Result<Vec<T>> {
        self.step(count + 1)?;
        let bytes = count.checked_mul(size_of::<T>()).ok_or(Error::Limit)?;
        self.storage(bytes)?;
        let mut v = Vec::new();
        v.try_reserve_exact(count).map_err(|_| Error::Limit)?;
        self.held = self
            .held
            .checked_add(v.capacity() * size_of::<T>())
            .ok_or(Error::Limit)?;
        self.storage(0)?;
        Ok(v)
    }
    fn release<T>(&mut self, v: &Vec<T>) {
        self.held -= v.capacity() * size_of::<T>();
    }
}
#[derive(Clone, Copy, Debug)]
pub struct HeadMeshReceipt {
    pub inherited_head_work: usize,
    pub shared_auxiliary_work: usize,
    pub auxiliary_mesh_work: usize,
    pub joining_work: usize,
    pub aggregate_work: usize,
    pub accounted_peak_bytes: usize,
    pub unchanged_surface_error_mm: f64,
    pub literal_transfer_rounding_mm: f64,
    pub cap_curve_chord_error_mm: f64,
    pub annulus_plane_lift_mm: f64,
    pub source_position_error_mm: f64,
    /// Positive-dot tan(theta) bound, hence also an angle upper bound.
    pub normal_tangent_cone: f64,
}
pub struct ActualHeadMeshCandidate<'owner> {
    source: &'owner OwnedActualHead,
    mesh: DisplayMesh,
    connectivity: MeshProof,
    receipt: HeadMeshReceipt,
}
impl<'owner> ActualHeadMeshCandidate<'owner> {
    pub fn source(&self) -> &'owner OwnedActualHead {
        self.source
    }
    pub fn mesh(&self) -> &DisplayMesh {
        &self.mesh
    }
    pub fn into_mesh(self) -> DisplayMesh {
        self.mesh
    }
    pub fn connectivity(&self) -> MeshProof {
        self.connectivity
    }
    pub fn receipt(&self) -> HeadMeshReceipt {
        self.receipt
    }
}
fn vertex(mesh: &DisplayMesh, id: u32) -> Result<[f64; 3]> {
    mesh.vertices
        .get(3 * id as usize..3 * id as usize + 3)
        .and_then(|p| p.try_into().ok())
        .ok_or(Error::Invalid)
}
impl OwnedActualHead {
    pub fn mesh_candidate(
        &self,
        options: HeadMeshOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<ActualHeadMeshCandidate<'_>> {
        options.tessellation.check().map_err(|_| Error::Invalid)?;
        if options.maximum_work == 0
            || options.maximum_work > MAX_MESH_WORK
            || options.maximum_bytes == 0
            || options.maximum_bytes > MAX_BYTES
        {
            return Err(Error::Invalid);
        }
        let auxiliary = self.auxiliary_thread();
        let shared = auxiliary.work_used();
        let fresh_head = self
            .receipt()
            .aggregate_work
            .checked_sub(shared)
            .ok_or(Error::Limit)?;
        let base = self
            .receipt()
            .retained_bytes
            .checked_add(options.additional_resident_bytes)
            .and_then(|n| n.checked_add(size_of::<ActualHeadMeshCandidate<'_>>()))
            .and_then(|n| n.checked_add(HEAD_SCRATCH))
            .ok_or(Error::Limit)?;
        let inherited_peak = self
            .receipt()
            .accounted_peak_bytes
            .checked_sub(self.receipt().additional_resident_bytes)
            .and_then(|n| n.checked_add(options.additional_resident_bytes))
            .and_then(|n| n.checked_add(size_of::<ActualHeadMeshCandidate<'_>>()))
            .and_then(|n| n.checked_add(HEAD_SCRATCH))
            .ok_or(Error::Limit)?;
        if base > options.maximum_bytes || inherited_peak > options.maximum_bytes {
            return Err(Error::Limit);
        }
        let mut w = Work {
            used: self.receipt().aggregate_work,
            remaining: options
                .maximum_work
                .checked_sub(self.receipt().aggregate_work)
                .ok_or(Error::Limit)?,
            base,
            held: 0,
            peak: base.max(inherited_peak),
            maximum_bytes: options.maximum_bytes,
            interrupt,
        };
        let sewing = self.axial_sewing_transfer();
        let factor = sewing.positive_axial_scale_bounds();
        let factor = I::new(factor[0], factor[1]);
        let child_angle = normal::child_angle(factor, options.tessellation.angle, &mut w)?;
        let sigma = factor.hi.max(1.);
        let local_deflection = (I::point(options.tessellation.deflection) / I::point(sigma)).lo;
        let outside = base
            .checked_sub(auxiliary.retained_budget_bytes())
            .ok_or(Error::Limit)?;
        let pre_child = w.used - self.receipt().aggregate_work;
        let literal = auxiliary.mesh_candidate(
            super::super::LiteralMeshOptions {
                tessellation: Tessellation {
                    deflection: local_deflection,
                    angle: child_angle,
                    ..options.tessellation
                },
                maximum_work: options
                    .maximum_work
                    .checked_sub(fresh_head)
                    .and_then(|n| n.checked_sub(pre_child))
                    .ok_or(Error::Limit)?,
                maximum_bytes: options
                    .maximum_bytes
                    .checked_sub(outside)
                    .filter(|n| *n > 0)
                    .ok_or(Error::Limit)?,
            },
            interrupt,
        )?;
        let child_work = literal.work_used();
        w.step(child_work.checked_sub(shared).ok_or(Error::Limit)?)?;
        w.peak = w.peak.max(
            literal
                .accounted_peak_bytes()
                .checked_add(outside)
                .ok_or(Error::Limit)?,
        );
        w.held = literal.retained_facet_bytes();
        w.storage(0)?;
        let old = literal.mesh();
        let traces = literal.cap_traces();
        let upper = traces
            .iter()
            .find(|t| auxiliary.faces()[t.face].entity() == self.head().auxiliary_cap_entity)
            .ok_or(Error::Unsupported)?;
        if upper.rim.len() != upper.edges.len() || upper.rim.len() < 3 {
            return Err(Error::Unsupported);
        }
        let vertices = old.vertices.len() / 3;
        let expected = old.face_indices.len().checked_add(22).ok_or(Error::Limit)?;
        if expected > options.tessellation.triangle_limit {
            return Err(Error::Limit);
        }
        let count = vertices
            .checked_add(self.head().vertices.len())
            .ok_or(Error::Limit)?;
        let mut mesh = DisplayMesh {
            vertices: w.allocate(count.checked_mul(3).ok_or(Error::Limit)?)?,
            indices: w.allocate(expected.checked_mul(3).ok_or(Error::Limit)?)?,
            face_indices: w.allocate(expected)?,
        };
        let mut rounding = 0f64;
        for &p in old.vertices.as_chunks::<3>().0 {
            w.step(128)?;
            let enclosure = sewing
                .map_raw_box_mm(p.map(|x| [x, x]))?
                .map(|x| I::new(x[0], x[1]));
            let literal = enclosure.map(I::midpoint);
            if literal.iter().any(|v| !v.is_finite()) {
                return Err(Error::Unsupported);
            }
            let error = (0..3).fold(I::point(0.), |s, k| {
                s + I::point((enclosure[k] - I::point(literal[k])).magnitude())
            });
            if !error.valid() {
                return Err(Error::Unsupported);
            }
            rounding = rounding.max(error.hi);
            mesh.vertices.extend(literal);
        }
        let transported = normal::transported(factor, literal.normal_error_radians(), &mut w)?;
        let roles = self.physical_roles();
        let mut normal_bound = 0f64;
        let mut removed = 0usize;
        for (&ids, &face) in old.indices.as_chunks::<3>().0.iter().zip(&old.face_indices) {
            w.step(64)?;
            if face as usize == upper.face {
                if !ids.contains(&upper.apex) {
                    return Err(Error::Unsupported);
                }
                removed += 1;
                continue;
            }
            let entity = auxiliary
                .faces()
                .get(face as usize)
                .ok_or(Error::Invalid)?
                .entity();
            let mapped = roles
                .thread_face_entities
                .iter()
                .position(|id| *id == entity)
                .ok_or(Error::Unsupported)? as u32;
            let raw = [
                vertex(old, ids[0])?,
                vertex(old, ids[1])?,
                vertex(old, ids[2])?,
            ];
            let actual = [
                vertex(&mesh, ids[0])?,
                vertex(&mesh, ids[1])?,
                vertex(&mesh, ids[2])?,
            ];
            let delta = normal::rounding(raw, actual, factor, &mut w)?;
            let cone = normal::combine(transported, delta, &mut w)?;
            if cone > options.tessellation.angle {
                return Err(Error::Unsupported);
            }
            normal_bound = normal_bound.max(cone);
            mesh.indices.extend(ids);
            mesh.face_indices.push(mapped);
        }
        if removed != upper.rim.len() {
            return Err(Error::Unsupported);
        }
        let head_base = u32::try_from(mesh.vertices.len() / 3).map_err(|_| Error::Limit)?;
        for p in &self.head().vertices {
            w.step(4)?;
            mesh.vertices.extend(p);
        }
        let stations = self.head_stations_mm();
        let mut outer = [0u32; 6];
        for (i, p) in self.head_polygon_mm().iter().enumerate() {
            w.step(64)?;
            let index = self
                .head()
                .vertices
                .iter()
                .position(|v| *v == [p[0], p[1], stations[0]])
                .ok_or(Error::Unsupported)?;
            outer[i] = head_base
                .checked_add(u32::try_from(index).map_err(|_| Error::Limit)?)
                .ok_or(Error::Limit)?;
        }
        let mut chord = 0f64;
        let mut lift = 0f64;
        for (i, e) in upper.edges.iter().enumerate() {
            w.step(64)?;
            if e.vertices != [upper.rim[i], upper.rim[(i + 1) % upper.rim.len()]] {
                return Err(Error::Unsupported);
            }
            let a = vertex(&mesh, e.vertices[0])?;
            let b = vertex(&mesh, e.vertices[1])?;
            chord = chord.max(caps::chord(
                auxiliary,
                sewing,
                e.edge,
                e.span,
                e.parameters,
                [a, b],
                &mut w,
            )?);
            lift = lift.max((I::point(a[2]) - I::point(stations[0])).magnitude());
        }
        let mut annuli = 0;
        for (i, face) in self.head().faces.iter().enumerate() {
            let outward = face.plane.normal.map(|n| if face.forward { n } else { -n });
            if !face.inner.is_empty() {
                annuli += 1;
                normal_bound = normal_bound.max(caps::annulus(
                    &mut mesh,
                    &outer,
                    upper.rim,
                    (9 + i) as u32,
                    outward,
                    options.tessellation.angle,
                    &mut w,
                )?);
            } else {
                let mut ring = [0u32; 6];
                for (j, &(edge, direction)) in face.outer.iter().enumerate() {
                    w.step(16)?;
                    let endpoints = self.head().edges.get(edge).ok_or(Error::Invalid)?;
                    ring[j] = head_base
                        .checked_add(
                            u32::try_from(endpoints[usize::from(!direction)])
                                .map_err(|_| Error::Limit)?,
                        )
                        .ok_or(Error::Limit)?;
                }
                for j in 1..face.outer.len() - 1 {
                    let ids = [ring[0], ring[j], ring[j + 1]];
                    let triangle =
                        ids.map(|id| std::array::from_fn(|k| mesh.vertices[3 * id as usize + k]));
                    let cone = normal::planar(outward, triangle, &mut w)?;
                    if cone > options.tessellation.angle {
                        return Err(Error::Unsupported);
                    }
                    normal_bound = normal_bound.max(cone);
                    mesh.indices.extend(ids);
                    mesh.face_indices.push((9 + i) as u32);
                }
            }
        }
        if annuli != 1 || mesh.face_indices.len() != expected {
            return Err(Error::Unsupported);
        }
        let unchanged =
            (I::point(literal.position_error_mm()) * I::point(sigma) + I::point(rounding)).hi;
        let annulus = (I::point(chord) + I::point(lift) + I::point(self.planar_sewing_mm())).hi;
        let distance = unchanged.max(annulus).max(self.planar_sewing_mm());
        if !distance.is_finite() || distance > options.tessellation.deflection {
            return Err(Error::Unsupported);
        }
        let literal_bytes = literal.retained_facet_bytes();
        drop(literal);
        w.held -= literal_bytes;
        // Remove only unused actual indices (the DERIVED cap apex). No weld,
        // coordinate fitting, new cap or altered physical vertex is invented.
        let mut remap = w.allocate(mesh.vertices.len() / 3)?;
        remap.resize(mesh.vertices.len() / 3, u32::MAX);
        let mut compact = DisplayMesh {
            vertices: w.allocate(mesh.vertices.len())?,
            indices: w.allocate(mesh.indices.len())?,
            face_indices: mesh.face_indices,
        };
        for &index in &mesh.indices {
            w.step(12)?;
            let id = &mut remap[index as usize];
            if *id == u32::MAX {
                *id = u32::try_from(compact.vertices.len() / 3).map_err(|_| Error::Limit)?;
                compact
                    .vertices
                    .extend_from_slice(&mesh.vertices[3 * index as usize..3 * index as usize + 3]);
            }
            compact.indices.push(*id);
        }
        w.release(&mesh.vertices);
        w.release(&mesh.indices);
        w.release(&remap);
        drop(mesh.vertices);
        drop(mesh.indices);
        drop(remap);
        let proof = manufacturing::certify_counted(&compact, &mut |n, scratch| {
            w.step(n)?;
            if let Some(bytes) = scratch {
                w.storage(bytes)?;
            }
            Ok(())
        })
        .map_err(|e| match e {
            manufacturing::CountedError::Mesh(_) => Error::Unsupported,
            manufacturing::CountedError::Interrupted(e) => e,
        })?;
        w.step(0)?;
        let joining = w
            .used
            .checked_sub(fresh_head)
            .and_then(|n| n.checked_sub(child_work))
            .ok_or(Error::Limit)?;
        Ok(ActualHeadMeshCandidate {
            source: self,
            mesh: compact,
            connectivity: proof,
            receipt: HeadMeshReceipt {
                inherited_head_work: self.receipt().aggregate_work,
                shared_auxiliary_work: shared,
                auxiliary_mesh_work: child_work,
                joining_work: joining,
                aggregate_work: w.used,
                accounted_peak_bytes: w.peak,
                unchanged_surface_error_mm: unchanged,
                literal_transfer_rounding_mm: rounding,
                cap_curve_chord_error_mm: chord,
                annulus_plane_lift_mm: lift,
                source_position_error_mm: distance,
                normal_tangent_cone: normal_bound,
            },
        })
    }
}
