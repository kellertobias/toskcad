//! Conforming literal facets of the actual retained degree-eleven annulus.
use super::*;
use crate::tapered_thread::{
    NormalPolynomial11, ReducedNormal11, ReducedNormalCell11, TensorCell11,
};
use std::sync::Arc;
const C0_WORK: usize = 8192;
const NORMAL_WORK: usize = ReducedNormal11::CELL_BOUND_WORK;
const SPLIT_WORK: usize = 16_384;
// At depth24 recursion retains <=49 position16KiB +normal4KiB cells.
// Eight16KiB position split/transpose temporaries, fixed normal workspaces,
// original root/cached-field constructors and metadata fit within2MiB.
pub(super) const SCRATCH: usize = 2 * 1024 * 1024;
pub(super) struct Thread {
    pub mesh: DisplayMesh,
    pub entry: Vec<u32>,
    pub tip: Vec<u32>,
    pub curves: Vec<[P; 12]>,
    pub divisions: Vec<usize>,
    pub error: f64,
}
fn tensor_error(e: crate::tapered_thread::Error) -> PhillipsReadError {
    match e {
        crate::tapered_thread::Error::Invalid => PhillipsReadError::Invalid,
        crate::tapered_thread::Error::Limit => PhillipsReadError::Limit,
        crate::tapered_thread::Error::Cancelled => PhillipsReadError::Cancelled,
        crate::tapered_thread::Error::Unresolved | crate::tapered_thread::Error::Unsupported => {
            PhillipsReadError::Unsupported
        }
    }
}
// Two actual position splits (2*16384), four N21 point evaluations
// (each <=8192: two Bernstein21 bases plus484 vector contractions),
// corners/triangle scores/validation and callback bookkeeping: <=69632.
const AXIS_WORK: usize = 69_632;
pub(super) fn preferred_axis(
    cell: &TensorCell11,
    field: &NormalPolynomial11,
    work: &mut Work<'_>,
) -> R<usize> {
    work.charge(AXIS_WORK)?;
    let mut scores = [f64::NEG_INFINITY; 2];
    for (axis, score) in scores.iter_mut().enumerate() {
        let children = cell.subdivide(axis).map_err(tensor_error)?;
        *score = f64::INFINITY;
        for child in children {
            let mid = child.domain().map(|[a, b]| [a + (b - a) * 0.5; 2]);
            let normal = field
                .at(mid)
                .map_err(tensor_error)?
                .map(|[a, b]| I::new(a, b).midpoint());
            let points = child
                .corner_bounds()
                .map(|p| p.map(|[a, b]| I::new(a, b).midpoint()));
            for ids in [[0, 1, 2], [0, 2, 3]] {
                let p = ids.map(|i| points[i]);
                let a: P = std::array::from_fn(|k| p[1][k] - p[0][k]);
                let b: P = std::array::from_fn(|k| p[2][k] - p[0][k]);
                let n: P = std::array::from_fn(|k| {
                    a[(k + 1) % 3] * b[(k + 2) % 3] - a[(k + 2) % 3] * b[(k + 1) % 3]
                });
                let projection = (0..3).map(|k| n[k] * normal[k]).sum::<f64>();
                let length = (n.into_iter().map(|x| x * x).sum::<f64>()
                    * normal.into_iter().map(|x| x * x).sum::<f64>())
                .sqrt();
                let cosine = projection / length;
                let cosine = if cosine.is_finite() {
                    cosine
                } else {
                    f64::NEG_INFINITY
                };
                *score = score.min(cosine);
            }
        }
    }
    Ok(usize::from(scores[1] > scores[0]))
}
struct Normals<'a> {
    source: &'a OwnedDegree11Material,
    values: Vec<Option<Arc<ReducedNormal11>>>,
    face_offsets: Vec<usize>,
    vspans: usize,
}
impl Normals<'_> {
    fn get(&mut self, patch: Patch, work: &mut Work<'_>) -> R<&Arc<ReducedNormal11>> {
        work.charge(1)?;
        let index =
            (self.face_offsets[patch.face as usize] + patch.root[0]) * self.vspans + patch.root[1];
        if self.values[index].is_none() {
            work.storage(SCRATCH + ReducedNormal11::WORKSPACE_BYTES + 128)?;
            work.charge(8 + 4 * self.source.thread_faces().len() + ReducedNormal11::WORK_BOUND)?;
            // This immutable cache was constructed from this SAME original
            // face/root by the earned volume and detached with shared Arcs.
            let field = self
                .source
                .thread_normal(patch.face as usize, patch.root)
                .ok_or(PhillipsReadError::Unsupported)?;
            let normal = Arc::new(ReducedNormal11::new(Arc::clone(field)).map_err(tensor_error)?);
            // The N21 payload and its Arc metadata are already counted by
            // source.retained_budget_bytes(); only this new B7 owner is fresh.
            let bytes = std::mem::size_of::<ReducedNormal11>() + 64;
            work.storage(SCRATCH + bytes)?;
            work.held += bytes;
            self.values[index] = Some(normal);
        }
        Ok(self.values[index]
            .as_ref()
            .expect("same-source normal initialized"))
    }
}
#[derive(Clone, Copy)]
struct Patch {
    offset: [usize; 2],
    face: u32,
    forward: bool,
    root: [usize; 2],
}
struct Grid {
    vertices: Vec<Option<P>>,
    indices: Vec<u32>,
    faces: Vec<u32>,
    shape: [usize; 2],
    divisions: [usize; 2],
    error: f64,
}
impl Grid {
    fn bytes(&self) -> usize {
        self.vertices.capacity() * std::mem::size_of::<Option<P>>()
            + (self.indices.capacity() + self.faces.capacity()) * 4
    }
    fn leaf(
        &mut self,
        cell: &TensorCell11,
        patch: Patch,
        normal: &ReducedNormalCell11,
        quality: Tessellation,
        work: &mut Work<'_>,
    ) -> R<std::result::Result<(), usize>> {
        let Patch {
            offset,
            face,
            forward,
            ..
        } = patch;
        // C0: <=8192 credits cover actual bilinear residual coefficients,
        // corner/literal/warp bounds, secants and loop/callback bookkeeping.
        // Refuse this proposal before any unnecessary normal calculation.
        work.charge(C0_WORK)?;
        let corners = cell
            .corner_bounds()
            .map(|p| p.map(|[lo, hi]| I::new(lo, hi)));
        let uv = cell.domain();
        let u = offset[0] + (uv[0][0] * self.divisions[0] as f64) as usize;
        let v = offset[1] + (uv[1][0] * self.divisions[1] as f64) as usize;
        let ids = [[u, v], [u + 1, v], [u + 1, v + 1], [u, v + 1]]
            .map(|[u, v]| ((u % self.shape[0]) * (self.shape[1] + 1) + v) as u32);
        let points: [P; 4] = std::array::from_fn(|j| {
            *self.vertices[ids[j] as usize].get_or_insert_with(|| corners[j].map(I::midpoint))
        });
        let literal = (0..4)
            .map(|j| {
                sum_norm(std::array::from_fn(|k| {
                    corners[j][k] - I::point(points[j][k])
                }))
                .hi
            })
            .fold(0., f64::max);
        // Actual coefficient residual bounds the whole source-to-bilinear
        // map. The diagonal triangles differ from that bilinear interpolant
        // by <=||P00-P10+P11-P01||/4; literal corner transfer is convex.
        let residual = cell.bilinear_error_bound_mm().map_err(tensor_error)?;
        let warp = sum_norm(std::array::from_fn(|k| {
            corners[0][k] - corners[1][k] + corners[2][k] - corners[3][k]
        })) / I::point(4.);
        let error = (I::point(residual) + warp + I::point(literal)).hi;
        // Secants choose a finite conditioning scalar only. For ANY such k
        // the actual coefficient identity preserves the complete normal.
        let du: P = std::array::from_fn(|k| {
            (points[1][k] + points[2][k] - points[0][k] - points[3][k]) * 0.5
        });
        let dv: P = std::array::from_fn(|k| {
            (points[3][k] + points[2][k] - points[0][k] - points[1][k]) * 0.5
        });
        let axis = || {
            let extent = |p: P| p.into_iter().map(f64::abs).sum::<f64>();
            usize::from(extent(dv) > extent(du))
        };
        if error > quality.deflection {
            return Ok(Err(axis()));
        }
        // SAME original N21 Arc owner and exact global coefficient residual
        // accompany the paired B7 restriction. The residual never shrinks.
        // This whole-leaf field enclosure, not samples, earns the normal gate.
        work.charge(NORMAL_WORK)?;
        debug_assert_eq!(normal.domain(), cell.domain());
        let bounds = normal.bounds().map_err(tensor_error)?;
        // Root-unit normal differs from this leaf's normal by the positive
        // factor h_u*h_v, which leaves the oriented normal angle unchanged.
        let mut normal_box = bounds.map(|[a, b]| I::new(a, b));
        if !forward {
            normal_box = normal_box.map(|n| -n);
        }
        let triangles = [[0, 1, 2], [0, 2, 3]].map(|mut t| {
            if !forward {
                t.swap(1, 2);
            }
            t
        });
        let mut refused = false;
        for tri in triangles {
            if normal_match(tri.map(|j| points[j]), normal_box, quality.angle) {
                continue;
            }
            work.charge(ReducedNormal11::CELL_CONE_WORK)?;
            let p = tri.map(|j| points[j].map(I::point));
            let n = super::super::geometry::cross(
                std::array::from_fn(|k| p[1][k] - p[0][k]),
                std::array::from_fn(|k| p[2][k] - p[0][k]),
            );
            let direction = n.map(|n| {
                let n = if forward { n } else { -n };
                [n.lo, n.hi]
            });
            if !normal
                .certifies_normal_cone(direction, quality.angle)
                .map_err(tensor_error)?
            {
                refused = true;
                break;
            }
        }
        if refused {
            // Actual child facets and point normals choose refinement ONLY.
            // Full C0 and coefficient-cone predicates still authorize leaves.
            let axis = preferred_axis(cell, normal.source().source(), work)?;
            return Ok(Err(axis));
        }
        for triangle in triangles {
            self.indices.extend(triangle.map(|j| ids[j]));
            self.faces.push(face);
        }
        self.error = self.error.max(error);
        Ok(Ok(()))
    }
    fn visit(
        &mut self,
        cell: TensorCell11,
        levels: [usize; 2],
        patch: Patch,
        normal: ReducedNormalCell11,
        quality: Tessellation,
        work: &mut Work<'_>,
    ) -> R<std::result::Result<(), usize>> {
        let axis = if levels[0] > 0 {
            0
        } else if levels[1] > 0 {
            1
        } else {
            return self.leaf(&cell, patch, &normal, quality, work);
        };
        work.charge(SPLIT_WORK + ReducedNormal11::CELL_SPLIT_WORK)?;
        let children = cell.subdivide(axis).map_err(tensor_error)?;
        let normal_children = normal.subdivide(axis).map_err(tensor_error)?;
        let mut next = levels;
        next[axis] -= 1;
        for (child, normal) in children.into_iter().zip(normal_children) {
            if let Err(axis) = self.visit(child, next, patch, normal, quality, work)? {
                // This original cell already requires a shared span split.
                // Remaining siblings cannot rescue this failed proposal.
                // The eventual successful pass still visits every leaf.
                return Ok(Err(axis));
            }
        }
        Ok(Ok(()))
    }
}
pub(super) fn build(
    source: &OwnedDegree11Material,
    quality: Tessellation,
    work: &mut Work<'_>,
) -> R<Thread> {
    let faces = source.thread_faces();
    let vspans = faces[0].source().spans()[1];
    let uspans: usize = faces.iter().map(|f| f.source().spans()[0]).sum();
    // The sealed occupied source has complete shared rails and one closed
    // seam. All adjacent nets share these exact original V coefficients.
    if faces.iter().any(|f| f.source().spans()[1] != vspans) {
        return Err(PhillipsReadError::Unsupported);
    }
    // One immutable global knot sequence per axis. Refinement marks only
    // original spans that fail, while every source row/face uses the SAME
    // resulting axis sequence. No hanging nodes or equal-v graph assumption.
    let metadata = (uspans + vspans + faces.len() + 4) * 128;
    work.storage(metadata + SCRATCH)?;
    work.held += metadata;
    let mut levels_u = vec![0usize; uspans];
    let mut levels_v = vec![0usize; vspans];
    let mut face_offsets = Vec::with_capacity(faces.len());
    let mut offset = 0;
    for face in faces {
        face_offsets.push(offset);
        offset += face.source().spans()[0];
    }
    let normal_count = uspans * vspans;
    let slot_bytes = normal_count * std::mem::size_of::<Option<Arc<ReducedNormal11>>>()
        + face_offsets.len() * std::mem::size_of::<usize>();
    work.storage(slot_bytes + SCRATCH)?;
    work.charge(normal_count)?;
    let mut normals = Normals {
        source,
        values: (0..normal_count).map(|_| None).collect(),
        face_offsets: face_offsets.clone(),
        vspans,
    };
    let actual_slots = normals.values.capacity()
        * std::mem::size_of::<Option<Arc<ReducedNormal11>>>()
        + normals.face_offsets.capacity() * std::mem::size_of::<usize>();
    work.storage(actual_slots + SCRATCH)?;
    work.held += actual_slots;
    for _ in 0..25 {
        if levels_u.iter().max().unwrap() + levels_v.iter().max().unwrap() > 24 {
            return Err(PhillipsReadError::Limit);
        }
        let divisions_u: Vec<_> = levels_u.iter().map(|&n| 1usize << n).collect();
        let divisions_v: Vec<_> = levels_v.iter().map(|&n| 1usize << n).collect();
        let prefix = |d: &[usize]| {
            let mut p = Vec::with_capacity(d.len() + 1);
            p.push(0usize);
            for &n in d {
                p.push(p.last().unwrap() + n);
            }
            p
        };
        let up = prefix(&divisions_u);
        let vp = prefix(&divisions_v);
        let shape = [*up.last().unwrap(), *vp.last().unwrap()];
        let vertices = shape[0]
            .checked_mul(shape[1] + 1)
            .ok_or(PhillipsReadError::Limit)?;
        let triangles = shape[0]
            .checked_mul(shape[1])
            .and_then(|n| n.checked_mul(2))
            .ok_or(PhillipsReadError::Limit)?;
        if triangles > quality.triangle_limit || vertices > u32::MAX as usize {
            return Err(PhillipsReadError::Limit);
        }
        let bytes = vertices * std::mem::size_of::<Option<P>>() + triangles * 16;
        work.storage(bytes + SCRATCH)?;
        work.charge(vertices + uspans + vspans)?;
        let mut grid = Grid {
            vertices: vec![None; vertices],
            indices: Vec::with_capacity(3 * triangles),
            faces: Vec::with_capacity(triangles),
            shape,
            divisions: [1, 1],
            error: 0.,
        };
        work.storage(grid.bytes() + SCRATCH)?;
        let grid_bytes = grid.bytes();
        work.held += grid_bytes;
        let mut mark_u = vec![false; uspans];
        let mut mark_v = vec![false; vspans];
        for (index, face) in faces.iter().enumerate() {
            let [us, vs] = face.source().spans();
            for u in 0..us {
                let global_u = face_offsets[index] + u;
                for v in 0..vs {
                    work.charge(1024)?;
                    let cell = face.source().cell([u, v]).map_err(tensor_error)?;
                    let levels = [levels_u[global_u], levels_v[v]];
                    grid.divisions = [divisions_u[global_u], divisions_v[v]];
                    let patch = Patch {
                        offset: [up[global_u], vp[v]],
                        face: index as u32,
                        forward: face.forward(),
                        root: [u, v],
                    };
                    work.charge(ReducedNormal11::CELL_ROOT_WORK)?;
                    let normal = normals.get(patch, work)?.root_cell();
                    if let Err(axis) = grid.visit(cell, levels, patch, normal, quality, work)? {
                        if axis == 0 {
                            mark_u[global_u] = true;
                        } else {
                            mark_v[v] = true;
                        }
                    }
                }
            }
        }
        if mark_u.iter().chain(&mark_v).any(|&b| b) {
            work.held -= grid_bytes;
            for (level, mark) in levels_u.iter_mut().zip(mark_u) {
                if mark {
                    *level += 1;
                }
            }
            for (level, mark) in levels_v.iter_mut().zip(mark_v) {
                if mark {
                    *level += 1;
                }
            }
            continue;
        }
        let final_bytes = vertices * 24 + grid.indices.capacity() * 4 + grid.faces.capacity() * 4;
        work.storage(final_bytes + 2 * shape[0] * 4 + uspans * 12 * 24 + SCRATCH)?;
        let mut mesh = DisplayMesh {
            vertices: Vec::with_capacity(3 * vertices),
            indices: grid.indices,
            face_indices: grid.faces,
        };
        for p in grid.vertices {
            work.charge(1)?;
            mesh.vertices
                .extend(p.ok_or(PhillipsReadError::Unsupported)?);
        }
        let mut entry: Vec<_> = (0..shape[0]).map(|u| (u * (shape[1] + 1)) as u32).collect();
        let mut tip: Vec<_> = (0..shape[0])
            .map(|u| (u * (shape[1] + 1) + shape[1]) as u32)
            .collect();
        if !faces[0].forward() {
            entry[1..].reverse();
            tip[1..].reverse();
        }
        let mut curves = Vec::with_capacity(uspans);
        let mut curve_divisions = Vec::with_capacity(uspans);
        for use_ in source.join_inner_uses().iter().rev().copied() {
            let net = faces[use_.face_index].source();
            let [_, height] = net.dimensions();
            for span in 0..net.spans()[0] {
                let span = if !use_.reversed {
                    net.spans()[0] - 1 - span
                } else {
                    span
                };
                curves.push(std::array::from_fn(|j| {
                    let j = if !use_.reversed { 11 - j } else { j };
                    net.poles()[(11 * span + j) * height]
                }));
                curve_divisions.push(divisions_u[face_offsets[use_.face_index] + span]);
            }
        }
        work.charge(uspans * 12)?;
        work.storage(
            mesh_bytes(&mesh)
                + (entry.capacity() + tip.capacity()) * 4
                + curves.capacity() * 12 * 24
                + SCRATCH,
        )?;
        work.held -= grid_bytes;
        return Ok(Thread {
            mesh,
            entry,
            tip,
            curves,
            divisions: curve_divisions,
            error: grid.error,
        });
    }
    Err(PhillipsReadError::Limit)
}
