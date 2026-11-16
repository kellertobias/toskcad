//! Source-free exact Native occurrence transport of a certified literal mesh.
//! Original decimal surfaces and their original STEP operator remain retained.
use super::*;
pub struct PlacedLiteralMesh {
    source: OwnedMaterial,
    placement: LiteralOccurrencePlacement,
    mesh: DisplayMesh,
    embedding: LiteralEmbeddingProof,
    distance: f64,
    normal_chord: f64,
    work: usize,
    preparation: usize,
    remaining: usize,
    peak: usize,
}
impl PlacedLiteralMesh {
    pub fn source(&self) -> &OwnedMaterial {
        &self.source
    }
    pub fn placement(&self) -> &LiteralOccurrencePlacement {
        &self.placement
    }
    pub fn mesh(&self) -> &DisplayMesh {
        &self.mesh
    }
    pub fn into_mesh(self) -> DisplayMesh {
        self.mesh
    }
    pub fn embedding(&self) -> LiteralEmbeddingProof {
        self.embedding
    }
    /// Euclidean correspondence to the original physical surface under this
    /// exact chain. The original candidate's stronger L1 contract is unchanged.
    pub fn position_error_mm(&self) -> f64 {
        self.distance
    }
    pub fn normal_chord(&self) -> f64 {
        self.normal_chord
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn preparation_work_used(&self) -> usize {
        self.preparation
    }
    pub fn query_work_used(&self) -> usize {
        self.work - self.preparation
    }
    pub fn remaining_work(&self) -> usize {
        self.remaining
    }
    pub fn accounted_peak_bytes(&self) -> usize {
        self.peak
    }
}
fn norm(p: [I; 3]) -> I {
    p.into_iter()
        .fold(I::point(0.), |s, v| s + v.square())
        .sqrt()
}
fn sub(a: [I; 3], b: [I; 3]) -> [I; 3] {
    std::array::from_fn(|i| a[i] - b[i])
}
fn cross(a: [I; 3], b: [I; 3]) -> [I; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn facet(mesh: &DisplayMesh, ids: [u32; 3]) -> Result<[I; 3]> {
    let [a, b, c] =
        ids.map(|id| std::array::from_fn(|i| I::point(mesh.vertices[id as usize * 3 + i])));
    let normal = cross(sub(b, a), sub(c, a));
    let length = norm(normal);
    if !length.valid() || length.lo <= 0. {
        return Err(Error::Unsupported);
    }
    Ok(normal.map(|v| v / length))
}
impl LiteralMeshCandidate {
    /// Continue this SAME caller ledger. Neither the original preparation nor
    /// the actual map/facet/embedding operations can reset its total quota.
    pub fn certify_placed(
        self,
        placement: LiteralOccurrencePlacement,
        quality: Tessellation,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<PlacedLiteralMesh> {
        quality.check().map_err(|_| Error::Invalid)?;
        if self.mesh.face_indices.len() > quality.triangle_limit {
            return Err(Error::Limit);
        }
        let mut budget = Budget {
            used: self.work,
            remaining: self.remaining,
            bytes: self.retained_phase_bytes,
            peak: self.peak,
            limit: self.byte_limit,
            interrupt,
        };
        hold(&mut budget, placement.retained_budget_bytes())?;
        budget.step(placement.work_used())?;
        let preparation = self
            .preparation
            .checked_add(placement.work_used())
            .ok_or(Error::Limit)?;
        let mut vertices = budget.allocate(self.mesh.vertices.len())?;
        let mut rounding = 0f64;
        for point in self.mesh.vertices.as_chunks::<3>().0 {
            budget.step(256)?;
            let exact = placement.map_box_mm(point.map(|v| [v, v]))?.map(interval);
            let literal = exact.map(I::midpoint);
            rounding = rounding.max(cells::l1(sub(exact, literal.map(I::point))));
            vertices.extend(literal);
        }
        let distance = (I::point(placement.linear_l2_norm_upper()) * I::point(self.position)
            + I::point(rounding))
        .hi;
        if distance > quality.deflection {
            return Err(Error::Unsupported);
        }
        let mut world = DisplayMesh {
            vertices,
            indices: vec![],
            face_indices: vec![],
        };
        // Every original normal certificate is a positive-dot cone with
        // tan(theta)<=r. Convert that SAME cone to a unit-normal chord without
        // cancellation in 1-1/sqrt(1+r²).
        let r = I::point(self.angle);
        let s = (I::point(1.) + r.square()).sqrt();
        let original_chord = (I::point(2.).sqrt() * r / (s * (s + I::point(1.))).sqrt()).hi;
        let allowed = (I::point(2.) * (I::point(quality.angle) / I::point(2.)).sin()).lo;
        let mut normal_chord = 0f64;
        for &ids in self.mesh.indices.as_chunks::<3>().0 {
            budget.step(1024)?;
            let original = facet(&self.mesh, ids)?;
            let exact = placement.normal_bounds(original.map(bounds))?.map(interval);
            let literal = facet(&world, ids)?;
            // Exact cofactor maps original source and original facet normals.
            // Its condition bounds normalization of both. Measure actual
            // literal facet rounding separately before retaining connectivity.
            let error = (I::point(placement.normal_condition_upper()) * I::point(original_chord)
                + norm(sub(literal, exact)))
            .hi;
            if error > allowed {
                return Err(Error::Unsupported);
            }
            normal_chord = normal_chord.max(error);
        }
        budget.release(&self.mesh.vertices);
        world.indices = self.mesh.indices;
        world.face_indices = self.mesh.face_indices;
        let embedding = embedding::certify(&world, &self.caps, &mut budget)?;
        Ok(PlacedLiteralMesh {
            source: self.source,
            placement,
            mesh: world,
            embedding,
            distance,
            normal_chord,
            work: budget.used,
            preparation,
            remaining: budget.remaining,
            peak: budget.peak,
        })
    }
}
impl LiteralGeometryPlan {
    /// The original complete-field L1 error also bounds its Euclidean error.
    /// Transport by the proved spectral norm; the placed gate measures every
    /// additional literal vertex rounding and retains no unearned L1 claim.
    pub fn placed_mesh(
        &self,
        placement: LiteralOccurrencePlacement,
        options: LiteralMeshOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<PlacedLiteralMesh> {
        options.check()?;
        self.mesh_candidate(options, interrupt)?.certify_placed(
            placement,
            options.tessellation,
            interrupt,
        )
    }
}
#[cfg(test)]
mod tests;
