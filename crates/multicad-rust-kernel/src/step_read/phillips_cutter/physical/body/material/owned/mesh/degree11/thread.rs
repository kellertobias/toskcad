//! Conforming literal facets of the actual retained degree-eleven annulus.
use super::*;
use crate::tapered_thread::{
    NormalPolynomial11, ReducedNormal11, ReducedNormalCell11, TensorCell11,
};
use std::sync::Arc;
mod accepted;
mod fast;
mod position;
const C0_WORK: usize = 8192;
const NORMAL_WORK: usize = ReducedNormal11::CELL_BOUND_WORK;
const SPLIT_WORK: usize = 16_384;
// At depth24 recursion retains <=49 position16KiB +normal4KiB cells.
// Eight16KiB position split/transpose temporaries, fixed normal workspaces,
// original root/cached-field constructors, the26-slot U-only MVT stack,
// prepared-cone buffer and row workspace fit within2MiB. Dynamic lattices,
// successful-root receipts and hard-leaf prefixes are separately resident.
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
fn chord_score(points: [P; 4], normal: P) -> f64 {
    [[0, 1, 2], [0, 2, 3]]
        .into_iter()
        .map(|ids| {
            let p = ids.map(|i| points[i]);
            let a: P = std::array::from_fn(|k| p[1][k] - p[0][k]);
            let b: P = std::array::from_fn(|k| p[2][k] - p[0][k]);
            let n: P = std::array::from_fn(|k| {
                a[(k + 1) % 3] * b[(k + 2) % 3] - a[(k + 2) % 3] * b[(k + 1) % 3]
            });
            let projection = (0..3).map(|k| n[k] * normal[k]).sum::<f64>();
            let cosine = projection
                / (n.into_iter().map(|x| x * x).sum::<f64>()
                    * normal.into_iter().map(|x| x * x).sum::<f64>())
                .sqrt();
            if cosine.is_finite() {
                cosine
            } else {
                f64::NEG_INFINITY
            }
        })
        .fold(f64::INFINITY, f64::min)
}
fn planned_normal_fit(
    normal: &ReducedNormalCell11,
    points: [P; 4],
    angle: f64,
    work: &mut Work<'_>,
) -> R<bool> {
    work.charge(NORMAL_WORK + 512)?;
    let bounds = match normal.bounds() {
        Ok(bounds) => bounds.map(|[a, b]| I::new(a, b)),
        Err(crate::tapered_thread::Error::Unresolved) => return Ok(false),
        Err(e) => return Err(tensor_error(e)),
    };
    for tri in [[0, 1, 2], [0, 2, 3]] {
        let points = tri.map(|j| points[j]);
        if normal_match(points, bounds, angle) {
            continue;
        }
        work.charge(ReducedNormal11::CELL_CONE_WORK)?;
        let p = points.map(|p| p.map(I::point));
        let direction = super::super::geometry::cross(
            std::array::from_fn(|k| p[1][k] - p[0][k]),
            std::array::from_fn(|k| p[2][k] - p[0][k]),
        )
        .map(|v| [v.lo, v.hi]);
        match normal.certifies_normal_cone(direction, angle) {
            Ok(true) => (),
            Ok(false) | Err(crate::tapered_thread::Error::Unresolved) => return Ok(false),
            Err(e) => return Err(tensor_error(e)),
        }
    }
    Ok(true)
}
fn preferred_root_axis(
    root: &position::Root,
    domain: [[f64; 2]; 2],
    normal: &ReducedNormalCell11,
    angle: f64,
    work: &mut Work<'_>,
) -> R<usize> {
    let points = root.priority_points(domain, work)?;
    // Planning only: prefer the first axis whose ACTUAL two child fields
    // already prove BOTH proposed facet normals. Final canonical-vertex C0,
    // full normal and closure gates still run on every resulting leaf.
    for axis in 0..2 {
        work.charge(ReducedNormal11::CELL_SPLIT_WORK + 512)?;
        let children = normal.subdivide(axis).map_err(tensor_error)?;
        let mut fits = true;
        for (child, field) in children.iter().enumerate() {
            let corners = if axis == 0 {
                [[child, 0], [child + 1, 0], [child + 1, 2], [child, 2]]
            } else {
                [[0, child], [2, child], [2, child + 1], [0, child + 1]]
            };
            if !planned_normal_fit(field, corners.map(|[u, v]| points[u][v]), angle, work)? {
                fits = false;
                break;
            }
        }
        if fits {
            return Ok(axis);
        }
    }
    // If neither complete child-angle proof decided, retain the prior actual
    // N21 center heuristic and its full independent work debit.
    work.charge(4 * 8192 + 2048)?;
    let field = normal.source().source();
    let mut scores = [f64::INFINITY; 2];
    for (axis, score) in scores.iter_mut().enumerate() {
        for child in 0..2 {
            let corners = if axis == 0 {
                [[child, 0], [child + 1, 0], [child + 1, 2], [child, 2]]
            } else {
                [[0, child], [2, child], [2, child + 1], [0, child + 1]]
            };
            let mut child_domain = domain;
            let mid = domain[axis][0] + (domain[axis][1] - domain[axis][0]) * 0.5;
            child_domain[axis][usize::from(child == 0)] = mid;
            let mid = child_domain.map(|[a, b]| [a + (b - a) * 0.5; 2]);
            let normal = field
                .at(mid)
                .map_err(tensor_error)?
                .map(|[a, b]| I::new(a, b).midpoint());
            *score = score.min(chord_score(corners.map(|[u, v]| points[u][v]), normal));
        }
    }
    Ok(usize::from(scores[1] > scores[0]))
}
#[derive(Clone, Copy)]
struct PositionContext<'a> {
    root: &'a position::Root,
    lattice: Option<&'a position::Lattice>,
}
struct Normals<'a> {
    source: &'a OwnedDegree11Material,
    values: Vec<Option<Arc<ReducedNormal11>>>,
    positions: Vec<Option<Arc<position::Root>>>,
    face_offsets: Vec<usize>,
    vspans: usize,
}
impl Normals<'_> {
    fn get(
        &mut self,
        patch: Patch,
        work: &mut Work<'_>,
    ) -> R<(&Arc<ReducedNormal11>, &Arc<position::Root>)> {
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
        if self.positions[index].is_none() {
            work.storage(SCRATCH + std::mem::size_of::<position::Root>() + 128)?;
            work.charge(position::ROOT_WORK)?;
            let root = Arc::new(position::Root::new(
                self.source.thread_faces()[patch.face as usize].source(),
                patch.root,
            )?);
            work.storage(SCRATCH + root.bytes())?;
            work.held += root.bytes();
            self.positions[index] = Some(root);
        }
        Ok((
            self.values[index]
                .as_ref()
                .expect("same-source normal initialized"),
            self.positions[index]
                .as_ref()
                .expect("same-source position initialized"),
        ))
    }
}
#[derive(Clone, Copy)]
struct Patch {
    offset: [usize; 2],
    face: u32,
    forward: bool,
    root: [usize; 2],
}
struct VisitNode {
    cell: Option<TensorCell11>,
    normal: ReducedNormalCell11,
    levels: [usize; 2],
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
        cell: Option<&TensorCell11>,
        positions: PositionContext<'_>,
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
        let PositionContext { root, lattice } = positions;
        // Complete-root Hessian guard includes both tensor interpolation and
        // its diagonal triangle warp. Actual shared-literal transfer is added.
        work.charge(position::C0_WORK)?;
        let uv = normal.domain();
        let corners = if let Some(lattice) = lattice {
            lattice.corners(uv)?
        } else {
            cell.expect("actual P11 traversal supplies corners")
                .corner_bounds()
                .map(|p| p.map(|[a, b]| I::new(a, b)))
        };
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
        let mut materialized = None;
        let mut error = (I::point(root.interpolation_error(uv)) + I::point(literal)).hi;
        if error > quality.deflection {
            if cell.is_none() {
                materialized = Some(root.materialize(uv, work)?);
            }
            let exact = cell
                .or(materialized.as_ref())
                .expect("actual source fallback initialized");
            work.charge(C0_WORK)?;
            let residual = exact.bilinear_error_bound_mm().map_err(tensor_error)?;
            let c = exact.corner_bounds().map(|p| p.map(|[a, b]| I::new(a, b)));
            let warp = sum_norm(std::array::from_fn(|k| {
                c[0][k] - c[1][k] + c[2][k] - c[3][k]
            })) / I::point(4.);
            error = (I::point(residual) + warp + I::point(literal)).hi;
        }
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
            let p = tri.map(|j| points[j].map(I::point));
            let n = super::super::geometry::cross(
                std::array::from_fn(|k| p[1][k] - p[0][k]),
                std::array::from_fn(|k| p[2][k] - p[0][k]),
            );
            let direction = n.map(|n| {
                let n = if forward { n } else { -n };
                [n.lo, n.hi]
            });
            // The source+residual projections are computed ONCE. Global and
            // grouped hulls only short-circuit on a sufficient cone proof;
            // undecided cases retain the SAME complete coefficient predicate.
            work.charge(ReducedNormal11::CELL_PREPARE_WORK)?;
            let prepared = normal
                .prepare_normal_cones(direction)
                .map_err(tensor_error)?;
            work.charge(ReducedNormal11::CELL_HULL_WORK)?;
            if prepared.hull_cone(quality.angle).map_err(tensor_error)? {
                continue;
            }
            work.charge(ReducedNormal11::CELL_GROUP_WORK)?;
            if prepared
                .grouped_cones(quality.angle)
                .map_err(tensor_error)?
            {
                continue;
            }
            work.charge(ReducedNormal11::CELL_REFINE_WORK)?;
            if !prepared
                .coefficient_cones(quality.angle)
                .map_err(tensor_error)?
            {
                refused = true;
                break;
            }
        }
        if refused {
            // Actual child facets and point normals choose refinement ONLY.
            // Full C0 and coefficient-cone predicates still authorize leaves.
            let axis = if let Some(exact) = cell.or(materialized.as_ref()) {
                preferred_axis(exact, normal.source().source(), work)?
            } else {
                preferred_root_axis(root, uv, normal, quality.angle, work)?
            };
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
        node: VisitNode,
        needs: Option<&fast::Needs>,
        positions: PositionContext<'_>,
        patch: Patch,
        quality: Tessellation,
        work: &mut Work<'_>,
    ) -> R<std::result::Result<(), usize>> {
        let VisitNode {
            cell,
            normal,
            levels,
        } = node;
        if let Some(needs) = needs {
            work.charge(64)?;
            if !needs.contains(normal.domain()) {
                return Ok(Ok(()));
            }
        }
        let axis = if levels[0] > 0 {
            0
        } else if levels[1] > 0 {
            1
        } else {
            return self.leaf(cell.as_ref(), positions, patch, &normal, quality, work);
        };
        work.charge(ReducedNormal11::CELL_SPLIT_WORK)?;
        let children = if let Some(cell) = cell {
            work.charge(SPLIT_WORK)?;
            cell.subdivide(axis).map_err(tensor_error)?.map(Some)
        } else {
            [None, None]
        };
        let normal_children = normal.subdivide(axis).map_err(tensor_error)?;
        let mut next = levels;
        next[axis] -= 1;
        for (child, normal) in children.into_iter().zip(normal_children) {
            if let Err(axis) = self.visit(
                VisitNode {
                    cell: child,
                    normal,
                    levels: next,
                },
                needs,
                positions,
                patch,
                quality,
                work,
            )? {
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
    work.storage(normal_count * std::mem::size_of::<Option<position::Lattice>>() + SCRATCH)?;
    work.charge(normal_count)?;
    let mut lattices: Vec<Option<position::Lattice>> = (0..normal_count).map(|_| None).collect();
    work.held += lattices.capacity() * std::mem::size_of::<Option<position::Lattice>>();
    work.storage(SCRATCH)?;
    work.storage(normal_count * std::mem::size_of::<Option<accepted::Accepted>>() + SCRATCH)?;
    work.charge(normal_count)?;
    let mut accepted: Vec<Option<accepted::Accepted>> = (0..normal_count).map(|_| None).collect();
    work.held += accepted.capacity() * std::mem::size_of::<Option<accepted::Accepted>>();
    work.storage(SCRATCH)?;
    let slot_bytes = normal_count
        * (std::mem::size_of::<Option<Arc<ReducedNormal11>>>()
            + std::mem::size_of::<Option<Arc<position::Root>>>())
        + face_offsets.len() * std::mem::size_of::<usize>();
    work.storage(slot_bytes + SCRATCH)?;
    work.charge(normal_count)?;
    let mut normals = Normals {
        source,
        values: (0..normal_count).map(|_| None).collect(),
        positions: (0..normal_count).map(|_| None).collect(),
        face_offsets: face_offsets.clone(),
        vspans,
    };
    let actual_slots = normals.values.capacity()
        * std::mem::size_of::<Option<Arc<ReducedNormal11>>>()
        + normals.positions.capacity() * std::mem::size_of::<Option<Arc<position::Root>>>()
        + normals.face_offsets.capacity() * std::mem::size_of::<usize>();
    work.storage(actual_slots + SCRATCH)?;
    work.held += actual_slots;
    // A density proposal is NOT a material/mesh certificate. Probe only the
    // actual entry-adjacent root of each U band; all other roots remain
    // untrusted until the complete pass below. No equal-v law is assumed.
    let mut planning = true;
    // Preserve the original25 complete-proof proposals IN ADDITION to a
    // separately bounded planning phase. Both phases share the SAME ledger.
    let mut phase_passes = 25;
    while phase_passes > 0 {
        phase_passes -= 1;
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
                for v in 0..if planning { 1 } else { vs } {
                    // A marked source span already makes this entire proposed
                    // grid insufficient. Defer its other roots until the next
                    // conforming proposal; no successful final pass skips one.
                    work.charge(1)?;
                    if !planning && (mark_u[global_u] || mark_v[v]) {
                        continue;
                    }
                    let levels = [levels_u[global_u], levels_v[v]];
                    grid.divisions = [divisions_u[global_u], divisions_v[v]];
                    let patch = Patch {
                        offset: [up[global_u], vp[v]],
                        face: index as u32,
                        forward: face.forward(),
                        root: [u, v],
                    };
                    let receipt_index = global_u * vspans + v;
                    if !planning
                        && let Some(receipt) = &accepted[receipt_index]
                        && receipt.reuse(&mut grid, patch, levels, work)?
                    {
                        continue;
                    }
                    work.charge(ReducedNormal11::CELL_ROOT_WORK)?;
                    let (field, position) = normals.get(patch, work)?;
                    let normal = field.root_cell();
                    let leaf_domain = [
                        [0., 1. / grid.divisions[0] as f64],
                        [0., 1. / grid.divisions[1] as f64],
                    ];
                    let fast = position.interpolation_error(leaf_domain) < quality.deflection;
                    let lattice = if fast {
                        Some(position.lattice_cached(
                            grid.divisions,
                            lattices[receipt_index].as_ref(),
                            work,
                        )?)
                    } else {
                        None
                    };
                    let lattice_bytes = lattice.as_ref().map_or(0, position::Lattice::bytes);
                    work.held += lattice_bytes;
                    let cell = if fast {
                        None
                    } else {
                        Some(position.materialize([[0., 1.]; 2], work)?)
                    };
                    let needs = if let Some(lattice) = &lattice {
                        Some(fast::Needs::discover(
                            &mut grid, position, lattice, &normal, patch, quality, work,
                        )?)
                    } else {
                        None
                    };
                    let needs_bytes = needs.as_ref().map_or(0, fast::Needs::bytes);
                    let positions = PositionContext {
                        root: position,
                        lattice: lattice.as_ref(),
                    };
                    let result = grid.visit(
                        VisitNode {
                            cell,
                            normal,
                            levels,
                        },
                        needs.as_ref(),
                        positions,
                        patch,
                        quality,
                        work,
                    )?;
                    drop(needs);
                    work.held -= needs_bytes;
                    if let Some(lattice) = lattice {
                        // Move the fully earned numeric enclosure into this
                        // SAME root's cache even if its normal gate refused.
                        // The new payload was already charged while local.
                        let old = lattices[receipt_index].replace(lattice);
                        let old_bytes = old.as_ref().map_or(0, position::Lattice::bytes);
                        drop(old);
                        work.held -= old_bytes;
                    }
                    if let Err(axis) = result {
                        if axis == 0 {
                            mark_u[global_u] = true;
                        } else {
                            mark_v[v] = true;
                        }
                    } else if !planning {
                        // Old+new payload overlap is preflighted by capture;
                        // replacing only this fixed original root preserves all
                        // inherited proof work and immutable source association.
                        let receipt = accepted::Accepted::capture(&grid, patch, levels, work)?;
                        let old = accepted[receipt_index]
                            .as_ref()
                            .map_or(0, accepted::Accepted::bytes);
                        work.held = work.held - old + receipt.bytes();
                        accepted[receipt_index] = Some(receipt);
                    }
                }
            }
        }
        if planning {
            let refine = mark_u.iter().chain(&mark_v).any(|&b| b);
            // Even successful probe roots publish no reusable geometry receipt.
            // Their immutable numerical caches remain fully charged above.
            drop(grid);
            work.held -= grid_bytes;
            work.charge(2 * (uspans + vspans) + 8)?;
            for (level, mark) in levels_u.iter_mut().zip(mark_u) {
                if mark {
                    *level += 1;
                }
            }
            if mark_v[0] {
                let proposal = levels_v[0] + 1;
                levels_v.fill(proposal);
            }
            if !refine || phase_passes == 0 {
                // Planning is optional: even when its own iteration ceiling
                // stops it, only the complete phase may publish a candidate.
                planning = false;
                phase_passes = 25;
            }
            continue;
        }
        if mark_u.iter().chain(&mark_v).any(|&b| b) {
            drop(grid);
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
        let lattice_cache_bytes = lattices.capacity()
            * std::mem::size_of::<Option<position::Lattice>>()
            + lattices
                .iter()
                .flatten()
                .map(position::Lattice::bytes)
                .sum::<usize>();
        drop(lattices);
        work.held -= lattice_cache_bytes;
        let accepted_bytes = accepted.capacity()
            * std::mem::size_of::<Option<accepted::Accepted>>()
            + accepted
                .iter()
                .flatten()
                .map(accepted::Accepted::bytes)
                .sum::<usize>();
        drop(accepted);
        work.held -= accepted_bytes;
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
