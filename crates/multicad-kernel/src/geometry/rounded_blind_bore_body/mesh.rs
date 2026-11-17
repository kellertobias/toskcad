//! One fresh producer/placement/closure/embedding query. No display proof is
//! reused for manufacture. Producer arenas retain their separate old scope.
use super::*;
use crate::geometry::internal_thread_treatment_body::interval::Interval as I;
use multicad_rust_kernel::{
    display::{DisplayMesh, Tessellation},
    manufacturing,
};
const WORK: usize = 512_000_000;
const BYTES: usize = 64 * 1024 * 1024;
const PRODUCER_WORK: usize = 24_000_000;
const SCRATCH: usize = 65536;
#[derive(Clone, Copy, Debug)]
pub struct RoundedBlindMeshBudget {
    pub maximum_work: usize,
    pub maximum_native_buffer_bytes: usize,
}
impl Default for RoundedBlindMeshBudget {
    fn default() -> Self {
        Self {
            maximum_work: WORK,
            maximum_native_buffer_bytes: BYTES,
        }
    }
}
pub struct RoundedBlindMeshReceipt<'s> {
    source: &'s BlindBore,
    remaining: usize,
    native: usize,
    peak: usize,
    pairs: usize,
    transfer: f64,
}
impl RoundedBlindMeshReceipt<'_> {
    pub fn source(&self) -> &BlindBore {
        self.source
    }
    pub fn producer_reserved_work(&self) -> usize {
        PRODUCER_WORK
    }
    /// Native balance debits, including the delegated connectivity reservation;
    /// not a count of executed CPU instructions.
    pub fn native_work_used(&self) -> usize {
        self.native
    }
    pub fn remaining_work(&self) -> usize {
        self.remaining
    }
    /// Native local/result buffers only. Producer and delegated validators
    /// keep their independently bounded arenas; this is not process RSS.
    pub fn native_buffer_peak_bytes(&self) -> usize {
        self.peak
    }
    pub fn embedding_pairs(&self) -> usize {
        self.pairs
    }
    pub fn maximum_literal_placement_error_mm(&self) -> f64 {
        self.transfer
    }
}
struct Work<'a> {
    remaining: usize,
    maximum_bytes: usize,
    peak: usize,
    interrupt: Option<&'a dyn Fn() -> bool>,
}
impl Work<'_> {
    fn charge(&mut self, n: usize) -> Result<(), GeometryError> {
        if self.interrupt.is_some_and(|f| f()) {
            return Err(failure(
                GeometryErrorKind::Cancelled,
                "Blind mesh cancelled",
            ));
        }
        self.remaining = self.remaining.checked_sub(n).ok_or_else(|| {
            failure(
                GeometryErrorKind::Limit,
                "Blind mesh aggregate work exhausted",
            )
        })?;
        Ok(())
    }
    fn storage(&mut self, n: usize) -> Result<(), GeometryError> {
        self.charge(1)?;
        let n = n.saturating_add(SCRATCH);
        if n > self.maximum_bytes {
            return Err(failure(
                GeometryErrorKind::Limit,
                "Blind native buffers exhausted",
            ));
        }
        self.peak = self.peak.max(n);
        Ok(())
    }
}
fn box_point(p: Point) -> [I; 3] {
    p.map(I::point)
}
fn delta(a: Point, b: Point) -> [I; 3] {
    std::array::from_fn(|k| I::point(a[k]) - I::point(b[k]))
}
fn cross_box(a: [I; 3], b: [I; 3]) -> [I; 3] {
    std::array::from_fn(|k| a[(k + 1) % 3] * b[(k + 2) % 3] - a[(k + 2) % 3] * b[(k + 1) % 3])
}
fn dot_box(a: [I; 3], b: [I; 3]) -> I {
    (0..3).fold(I::ZERO, |sum, k| sum + a[k] * b[k])
}
struct Placement {
    columns: [Point; 3],
    origin: Point,
    cofactor: [[I; 3]; 3],
    sigma: f64,
    condition: f64,
}
impl Placement {
    fn new(frame: PlaneFrame) -> Result<Self, GeometryError> {
        let columns = [frame.u, frame.v, frame.normal()];
        let a = columns.map(box_point);
        let cofactor = [
            cross_box(a[1], a[2]),
            cross_box(a[2], a[0]),
            cross_box(a[0], a[1]),
        ];
        if dot_box(a[0], cofactor[0]).lo <= 0. {
            return Err(failure(
                GeometryErrorKind::Unsupported,
                "Blind placement orientation is unresolved",
            ));
        }
        let (mut lower, mut upper) = (f64::INFINITY, 0_f64);
        for i in 0..3 {
            let diagonal = dot_box(a[i], a[i]);
            let mut off = I::ZERO;
            for j in 0..3 {
                if i != j {
                    off = off + I::point(dot_box(a[i], a[j]).magnitude());
                }
            }
            lower = lower.min((diagonal - off).lo);
            upper = upper.max((diagonal + off).hi);
        }
        if lower <= 0. || !upper.is_finite() {
            return Err(failure(
                GeometryErrorKind::Unsupported,
                "Blind placement metric is unresolved",
            ));
        }
        let sigma = I::point(upper).sqrt().hi;
        let condition = (I::point(upper) / I::point(lower)).sqrt().hi;
        Ok(Self {
            columns,
            origin: frame.origin,
            cofactor,
            sigma,
            condition,
        })
    }
    fn point(&self, p: Point, reserve: f64) -> Result<(Point, f64), GeometryError> {
        let q: [I; 3] = std::array::from_fn(|i| {
            (0..3).fold(I::point(self.origin[i]), |sum, k| {
                sum + I::point(self.columns[k][i]) * I::point(p[k])
            })
        });
        if q.iter().any(|v| !v.valid()) {
            return Err(failure(
                GeometryErrorKind::Unsupported,
                "Blind point transfer is unresolved",
            ));
        }
        let literal = q.map(|v| v.lo + (v.hi - v.lo) * 0.5);
        let error = (0..3)
            .fold(I::ZERO, |sum, k| {
                sum + I::point((q[k] - I::point(literal[k])).magnitude())
            })
            .hi;
        if error > reserve {
            return Err(failure(
                GeometryErrorKind::Unsupported,
                "Blind literal point exceeds transfer reserve",
            ));
        }
        Ok((literal, error))
    }
    fn normal(
        &self,
        local: [Point; 3],
        world: [Point; 3],
        reserve: f64,
    ) -> Result<(), GeometryError> {
        let n = cross_box(delta(local[1], local[0]), delta(local[2], local[0]));
        let exact = std::array::from_fn(|i| {
            (0..3).fold(I::ZERO, |sum, j| sum + self.cofactor[j][i] * n[j])
        });
        let actual = cross_box(delta(world[1], world[0]), delta(world[2], world[0]));
        let projection = dot_box(exact, actual);
        let perpendicular = cross_box(exact, actual)
            .into_iter()
            .fold(I::ZERO, |sum, v| sum + I::point(v.magnitude()));
        let ratio = perpendicular / projection;
        if projection.valid() && projection.lo > 0. && ratio.valid() && ratio.hi <= reserve {
            Ok(())
        } else {
            Err(failure(
                GeometryErrorKind::Unsupported,
                "Blind literal normals exceed transfer reserve",
            ))
        }
    }
}
fn buffers(mesh: &DisplayMesh) -> usize {
    mesh.vertices.capacity() * 8 + (mesh.indices.capacity() + mesh.face_indices.capacity()) * 4
}
impl RoundedBlindBoreBody {
    fn fresh<'s>(
        &'s self,
        frame: PlaneFrame,
        options: Tessellation,
        budget: RoundedBlindMeshBudget,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<(Mesh, manufacturing::MeshProof, RoundedBlindMeshReceipt<'s>), GeometryError> {
        if !options.deflection.is_finite()
            || options.deflection <= 0.
            || !options.angle.is_finite()
            || options.angle < 0.01
            || options.angle > 1.
            || !(1..=5_000_000).contains(&options.triangle_limit)
        {
            return Err(failure(
                GeometryErrorKind::Input,
                "Invalid blind mesh quality",
            ));
        }
        if budget.maximum_work > WORK
            || budget.maximum_work == 0
            || budget.maximum_native_buffer_bytes > BYTES
            || budget.maximum_native_buffer_bytes == 0
        {
            return Err(failure(
                GeometryErrorKind::Input,
                "Invalid blind mesh budget",
            ));
        }
        let mut work = Work {
            remaining: budget.maximum_work,
            maximum_bytes: budget.maximum_native_buffer_bytes,
            peak: 0,
            interrupt,
        };
        work.charge(1024)?;
        let placement = Placement::new(frame)?;
        let quality = Tessellation {
            deflection: (I::point(options.deflection) / (I::point(2.) * I::point(placement.sigma)))
                .lo,
            angle: (I::point(0.99) * I::point(options.angle) / I::point(placement.condition)).lo,
            ..options
        };
        if quality.angle < 0.01 {
            return Err(failure(
                GeometryErrorKind::Unsupported,
                "Blind world quality leaves no source angular reserve",
            ));
        }
        // The unchanged producer owns its24M hard work cap. Reserve that whole
        // delegated scope once; do not relabel it as measured native work.
        work.charge(PRODUCER_WORK)?;
        let local =
            multicad_rust_kernel::rounded_thread_mesh::tessellate_blind_bore_with_interrupt(
                &self.geometry,
                quality,
                interrupt,
            )
            .map_err(|e| {
                failure(
                    match e {
                        multicad_rust_kernel::rounded_thread_mesh::Error::Cancelled => {
                            GeometryErrorKind::Cancelled
                        }
                        _ => GeometryErrorKind::Limit,
                    },
                    format!("Actual rounded blind source mesh refused: {e:?}"),
                )
            })?;
        let local_bytes = buffers(&local);
        work.storage(local_bytes + local.vertices.len() * 8)?;
        let mut vertices = Vec::with_capacity(local.vertices.len());
        work.storage(local_bytes + vertices.capacity() * 8)?;
        let mut transfer = 0_f64;
        for &p in local.vertices.as_chunks::<3>().0 {
            work.charge(256)?;
            let (p, error) = placement.point(p, options.deflection / 4.)?;
            transfer = transfer.max(error);
            vertices.extend(p);
        }
        for ids in local.indices.as_chunks::<3>().0 {
            work.charge(256)?;
            let a = ids.map(|id| std::array::from_fn(|k| local.vertices[3 * id as usize + k]));
            let b = ids.map(|id| std::array::from_fn(|k| vertices[3 * id as usize + k]));
            placement.normal(a, b, options.angle / 200.)?;
        }
        let world = DisplayMesh {
            vertices,
            indices: local.indices,
            face_indices: local.face_indices,
        };
        drop(local.vertices);
        work.storage(buffers(&world))?;
        // Same explicit delegated connectivity reservation as the existing
        // Native repaired-source adapter; validator arenas remain separate.
        work.charge(
            world
                .face_indices
                .len()
                .saturating_mul(52)
                .saturating_add(world.vertices.len() / 3),
        )?;
        let proof = manufacturing::certify(&world).map_err(|e| {
            failure(
                GeometryErrorKind::InvalidResult,
                format!("Blind world manifold refused: {e:?}"),
            )
        })?;
        let mesh = Mesh {
            vertices: world.vertices,
            indices: world.indices,
            face_indices: world.face_indices,
        };
        work.charge(1)?;
        let (pairs, used) = {
            let proof = mesh.prove_literal_embedding(work.remaining, interrupt)?;
            (proof.pairs(), proof.work_used())
        };
        work.charge(used)?;
        let receipt = RoundedBlindMeshReceipt {
            source: &self.geometry,
            remaining: work.remaining,
            native: budget.maximum_work - work.remaining - PRODUCER_WORK,
            peak: work.peak,
            pairs,
            transfer,
        };
        Ok((mesh, proof, receipt))
    }
    pub(in crate::geometry) fn mesh(
        &self,
        frame: PlaneFrame,
        options: Tessellation,
    ) -> Result<Mesh, GeometryError> {
        self.fresh(frame, options, Default::default(), None)
            .map(|x| x.0)
    }
    pub(in crate::geometry) fn manufacturing(
        &self,
        frame: PlaneFrame,
        options: MeshOptions,
    ) -> Result<crate::ManufacturingMesh, GeometryError> {
        let (mesh, proof, _) = self.fresh(
            frame,
            Solid::rust_mesh_options(options),
            Default::default(),
            None,
        )?;
        Ok(crate::ManufacturingMesh {
            mesh,
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
}
impl Solid {
    pub fn rust_rounded_blind_bore_mesh_with_receipt(
        &self,
        options: MeshOptions,
        budget: RoundedBlindMeshBudget,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<(Mesh, RoundedBlindMeshReceipt<'_>), GeometryError> {
        let (body, frame) = self.rounded_blind_bore_body().ok_or_else(|| {
            failure(
                GeometryErrorKind::Unsupported,
                "Solid has no authored rounded blind source",
            )
        })?;
        body.fresh(frame, Self::rust_mesh_options(options), budget, interrupt)
            .map(|(mesh, _, receipt)| (mesh, receipt))
    }
}
