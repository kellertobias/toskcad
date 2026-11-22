//! Exact affine placement of an already source-certified literal candidate.
//! Connectivity is retained; whole literal embedding remains a separate gate.
use super::*;
#[derive(Clone, Copy, Debug)]
pub struct StockClipPlacement {
    matrix: [[f64; 3]; 3],
    translation: [f64; 3],
    cofactor: [[I; 3]; 3],
    norm: I,
    condition: I,
    determinant: I,
    metric_error: f64,
    area_scale: [f64; 2],
}
fn dot(a: [I; 3], b: [I; 3]) -> I {
    a.into_iter()
        .zip(b)
        .fold(I::point(0.), |s, (a, b)| s + a * b)
}
impl StockClipPlacement {
    /// The literal matrix is authoritative: it is neither normalized nor fitted.
    /// Only bounded, proper, near-isometric placements enter this provider.
    pub fn new(matrix: [[f64; 3]; 3], translation: [f64; 3]) -> Result<Self, Error> {
        if matrix
            .iter()
            .flatten()
            .chain(&translation)
            .any(|v| !v.is_finite())
        {
            return Err(Error::Invalid);
        }
        if translation.into_iter().any(|v| v.abs() > 1e7)
            || matrix.iter().flatten().any(|v| v.abs() > 2.)
        {
            return Err(Error::Limit);
        }
        let m = matrix.map(|r| r.map(I::point));
        let columns: [[I; 3]; 3] = std::array::from_fn(|j| std::array::from_fn(|i| m[i][j]));
        let mut delta = I::point(0.);
        for i in 0..3 {
            for j in 0..3 {
                delta = delta
                    + absolute(
                        dot(columns[i], columns[j]) - I::point(if i == j { 1. } else { 0. }),
                    );
            }
        }
        if delta.hi > 1e-6 {
            return Err(Error::Unsupported);
        }
        let pairs = [(1, 2), (2, 0), (0, 1)];
        let cofactor = std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                let (a, b) = pairs[i];
                let (c, d) = pairs[j];
                m[a][c] * m[b][d] - m[a][d] * m[b][c]
            })
        });
        let determinant = dot(m[0], cofactor[0]);
        if determinant.lo <= 0. {
            return Err(Error::Unsupported);
        }
        // ||RᵀR-I||₂ <= the sum of absolute entries. Singular values of R
        // lie in sqrt(1±delta); cofactor singular values are pairwise products.
        let lower = I::point(1.) - I::point(delta.hi);
        let upper = I::point(1.) + I::point(delta.hi);
        Ok(Self {
            matrix,
            translation,
            cofactor,
            norm: upper.sqrt(),
            condition: upper / lower,
            determinant,
            metric_error: delta.hi,
            area_scale: [lower.lo, upper.hi],
        })
    }
    pub fn matrix(&self) -> [[f64; 3]; 3] {
        self.matrix
    }
    pub fn translation_mm(&self) -> [f64; 3] {
        self.translation
    }
    pub fn determinant_bounds(&self) -> [f64; 2] {
        bounds(self.determinant)
    }
    pub fn metric_error(&self) -> f64 {
        self.metric_error
    }
    /// Bounds every local surface-area multiplier of this affine map.
    pub fn area_scale_bounds(&self) -> [f64; 2] {
        self.area_scale
    }
    /// Fixed-size outward image of an actual local box under the literal map.
    pub fn point_bounds(&self, local: [[f64; 2]; 3]) -> Result<[[f64; 2]; 3], Error> {
        let local = coordinate_box(local)?;
        bounded(std::array::from_fn(|i| {
            dot(self.matrix[i].map(I::point), local) + I::point(self.translation[i])
        }))
    }
    /// The true affine inverse, including cofactor/determinant rounding bounds.
    /// Transposing a near-orthogonal matrix is not substituted for its inverse.
    pub fn inverse_point_bounds(&self, world: [[f64; 2]; 3]) -> Result<[[f64; 2]; 3], Error> {
        let world = coordinate_box(world)?;
        let shifted = std::array::from_fn(|i| world[i] - I::point(self.translation[i]));
        bounded(std::array::from_fn(|i| {
            dot(std::array::from_fn(|j| self.cofactor[j][i]), shifted) / self.determinant
        }))
    }
    fn point(&self, local: [f64; 3]) -> [I; 3] {
        std::array::from_fn(|i| {
            dot(self.matrix[i].map(I::point), local.map(I::point)) + I::point(self.translation[i])
        })
    }
    fn normal(&self, n: [I; 3]) -> Result<[I; 3], Error> {
        let out = self.cofactor.map(|r| dot(r, n));
        let length = norm(out);
        if length.lo <= 0. || !length.hi.is_finite() {
            return Err(Error::Unresolved);
        }
        Ok(out.map(|v| v / length))
    }
}
fn coordinate_box(box_: [[f64; 2]; 3]) -> Result<[I; 3], Error> {
    if box_
        .into_iter()
        .any(|b| !b[0].is_finite() || !b[1].is_finite() || b[0] > b[1])
    {
        return Err(Error::Invalid);
    }
    if box_.iter().flatten().any(|x| x.abs() > 1e7) {
        return Err(Error::Limit);
    }
    Ok(box_.map(|b| I::new(b[0], b[1])))
}
fn bounded(point: [I; 3]) -> Result<[[f64; 2]; 3], Error> {
    if point
        .into_iter()
        .any(|p| !p.lo.is_finite() || !p.hi.is_finite() || p.magnitude() > 1e7)
    {
        return Err(Error::Limit);
    }
    Ok(point.map(bounds))
}
#[derive(Debug)]
pub struct PlacedStockClipMeshCandidate {
    mesh: DisplayMesh,
    proof: MeshProof,
    distance: f64,
    normal: f64,
    source_volume: [f64; 2],
    literal_volume: [f64; 2],
    work: usize,
    bytes: usize,
}
impl PlacedStockClipMeshCandidate {
    pub fn mesh(&self) -> &DisplayMesh {
        &self.mesh
    }
    pub fn into_mesh(self) -> DisplayMesh {
        self.mesh
    }
    pub fn connectivity(&self) -> MeshProof {
        self.proof
    }
    pub fn source_distance_mm(&self) -> f64 {
        self.distance
    }
    pub fn normal_chord(&self) -> f64 {
        self.normal
    }
    pub fn source_volume_mm3(&self) -> [f64; 2] {
        self.source_volume
    }
    pub fn literal_volume_mm3(&self) -> [f64; 2] {
        self.literal_volume
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_budget_bytes(&self) -> usize {
        self.bytes
    }
}
impl StockClipMeshCandidate {
    /// The same caller total includes this candidate's inherited work once.
    /// Quality is certified in the actual placed space, including rounding.
    pub fn placed(
        mut self,
        placement: StockClipPlacement,
        quality: Tessellation,
        maximum_work: usize,
        maximum_bytes: usize,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<PlacedStockClipMeshCandidate, Error> {
        quality.check().map_err(|_| Error::Invalid)?;
        if maximum_work == 0
            || maximum_work > 512_000_000
            || maximum_bytes == 0
            || maximum_bytes > 64 * 1024 * 1024
        {
            return Err(Error::Invalid);
        }
        let mut work = Work {
            used: 0,
            limit: maximum_work,
            bytes: maximum_bytes,
            peak: 0,
            interrupt,
        };
        work.charge(self.work)?;
        work.memory(self.bytes + self.mesh.vertices.len() * 8 + 16 * 1024)?;
        let mut world = Vec::new();
        world
            .try_reserve_exact(self.mesh.vertices.len())
            .map_err(|_| Error::Limit)?;
        work.observe(self.bytes + world.capacity() * 8 + 16 * 1024);
        let mut rounding = 0_f64;
        for point in self.mesh.vertices.as_chunks::<3>().0 {
            work.charge(128)?;
            let exact = placement.point(*point);
            let literal = exact.map(I::midpoint);
            if literal.into_iter().any(|v| !v.is_finite() || v.abs() > 1e7) {
                return Err(Error::Limit);
            }
            rounding =
                rounding.max(norm(std::array::from_fn(|i| I::point(literal[i]) - exact[i])).hi);
            world.extend(literal);
        }
        let distance = (placement.norm * I::point(self.chord) + I::point(rounding)).hi;
        if distance > quality.deflection {
            return Err(Error::Unresolved);
        }
        let mut placed = DisplayMesh {
            vertices: world,
            indices: vec![],
            face_indices: vec![],
        };
        let allowed = (I::point(2.) * (I::point(quality.angle) / I::point(2.)).sin()).lo;
        let mut maximum_normal = 0_f64;
        let mut volume = I::point(0.);
        for ids in self.mesh.indices.as_chunks::<3>().0 {
            work.charge(640)?;
            let local_normal = facet(&self.mesh, *ids)?;
            let exact_normal = placement.normal(local_normal)?;
            let literal_normal = facet(&placed, *ids)?;
            // For equal-unit u,v, ||Au/|Au|-Av/|Av||| <= kappa(A)||u-v||:
            // |a-b|²=(|a|-|b|)²+|a||b|*|a/|a|-b/|b||². Apply this to
            // the actual cofactor, then add measured literal facet rounding.
            let error = (placement.condition * I::point(self.normal)
                + norm(difference(literal_normal, exact_normal)))
            .hi;
            if error > allowed {
                return Err(Error::Unresolved);
            }
            maximum_normal = maximum_normal.max(error);
            let [a, b, c] = ids.map(|id| at(&placed, id));
            volume = volume + dot(a, cross(b, c)) / I::point(6.);
        }
        // SAME index arrays and vertex count retain the earned connectivity.
        // Every new literal facet has independently proved nonzero area above.
        placed.indices = std::mem::take(&mut self.mesh.indices);
        placed.face_indices = std::mem::take(&mut self.mesh.face_indices);
        let source = I::new(self.source_volume[0], self.source_volume[1]) * placement.determinant;
        work.charge(1)?;
        Ok(PlacedStockClipMeshCandidate {
            mesh: placed,
            proof: self.proof,
            distance,
            normal: maximum_normal,
            source_volume: bounds(source),
            literal_volume: bounds(volume),
            work: work.used,
            bytes: work.peak,
        })
    }
}
#[cfg(test)]
mod tests;
