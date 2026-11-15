//! Actual literal facets checked against their original source cells.
use super::*;
use crate::{
    display::{DisplayMesh, Tessellation},
    manufacturing::{self, MeshProof},
};
mod placement;
pub use placement::{PlacedStockClipMeshCandidate, StockClipPlacement};
#[derive(Debug)]
pub struct StockClipMeshCandidate {
    mesh: DisplayMesh,
    proof: MeshProof,
    stock_radius: f64,
    chord: f64,
    normal: f64,
    source_volume: [f64; 2],
    literal_volume: [f64; 2],
    work: usize,
    bytes: usize,
}
impl StockClipMeshCandidate {
    pub fn mesh(&self) -> &DisplayMesh {
        &self.mesh
    }
    /// Move these SAME certified literal arrays without duplicating storage.
    pub fn into_mesh(self) -> DisplayMesh {
        self.mesh
    }
    pub fn connectivity(&self) -> MeshProof {
        self.proof
    }
    /// Directed maximum radius of actual literal vertices. Convexity then
    /// confines their triangles to the unchanged source stock cylinder.
    pub fn maximum_literal_radius_mm(&self) -> f64 {
        self.stock_radius
    }
    pub fn source_distance_mm(&self) -> f64 {
        self.chord
    }
    pub fn normal_chord(&self) -> f64 {
        self.normal
    }
    pub fn source_volume_mm3(&self) -> [f64; 2] {
        self.source_volume
    }
    /// Outward signed polyhedral flux; occupied volume additionally requires
    /// the separate whole literal embedding proof.
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
fn norm(p: [I; 3]) -> I {
    (p[0].square() + p[1].square() + p[2].square()).sqrt()
}
fn cross(a: [I; 3], b: [I; 3]) -> [I; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn at(mesh: &DisplayMesh, id: u32) -> [I; 3] {
    std::array::from_fn(|i| I::point(mesh.vertices[id as usize * 3 + i]))
}
fn difference(a: [I; 3], b: [I; 3]) -> [I; 3] {
    std::array::from_fn(|i| a[i] - b[i])
}
fn normal(s: I, r: I, rs: I, rz: I, omega: I) -> [I; 3] {
    let angle = -omega * s;
    let t = rs / (omega * r);
    let n = (I::point(1.) + t.square() + rz.square()).sqrt();
    [
        (angle.cos() - t * angle.sin()) / n,
        (angle.sin() + t * angle.cos()) / n,
        -rz / n,
    ]
}
fn slopes(p: StockClipPatch, k: I) -> [I; 2] {
    let m = I::point(3.).sqrt();
    match p.band {
        StockClipBand::Lower => [m - k, m],
        StockClipBand::Upper => [-m - k, -m],
        StockClipBand::Crest if p.slab == StockClipSlab::After => [I::point(0.); 2],
        _ => [I::point(0.), k],
    }
}
fn absolute(p: I) -> I {
    I::point(p.lo.abs().max(p.hi.abs()))
}
// Raw literal numeric checks: these do not manufacture source authority.
fn stock_radius(point: [f64; 3], limit: f64) -> Result<f64, Error> {
    let bound = (I::point(point[0]).square() + I::point(point[1]).square())
        .sqrt()
        .hi;
    if !bound.is_finite() || bound > limit {
        return Err(Error::Unresolved);
    }
    Ok(bound)
}
fn positive_area(points: [[I; 2]; 3]) -> Result<(), Error> {
    let a = [points[1][0] - points[0][0], points[1][1] - points[0][1]];
    let b = [points[2][0] - points[0][0], points[2][1] - points[0][1]];
    if (a[0] * b[1] - a[1] * b[0]).lo <= 0. {
        return Err(Error::Unresolved);
    }
    Ok(())
}
fn facet(mesh: &DisplayMesh, ids: [u32; 3]) -> Result<[I; 3], Error> {
    let [a, b, c] = ids.map(|id| at(mesh, id));
    let raw = cross(difference(b, a), difference(c, a));
    let length = norm(raw);
    if !length.lo.is_finite() || !length.hi.is_finite() || length.lo <= 0. {
        return Err(Error::Unresolved);
    }
    Ok(raw.map(|v| v / length))
}
fn push(mesh: &mut DisplayMesh, ids: [u32; 3], face: u32) {
    mesh.indices.extend(ids);
    mesh.face_indices.push(face)
}
fn face(p: StockClipPatch) -> u32 {
    match (p.band, p.slab) {
        (StockClipBand::Lower, _) => 0,
        (StockClipBand::Crest, StockClipSlab::Before) => 1,
        (StockClipBand::Crest, StockClipSlab::After) => 2,
        (StockClipBand::Upper, _) => 3,
        (StockClipBand::Root, _) => 4,
    }
}

#[derive(Clone, Copy)]
struct CertifiedTriangle {
    ids: [u32; 3],
    distance: f64,
    normal: f64,
}
struct TriangleContext<'a> {
    plan: &'a StockClipCellPlan,
    mesh: &'a DisplayMesh,
    errors: &'a [f64],
    omega: I,
    k: I,
    quality: Tessellation,
    allowed: f64,
}
impl TriangleContext<'_> {
    fn certify(
        &self,
        cell: &SourceCell,
        local: [usize; 3],
        work: &mut Work<'_>,
    ) -> Result<CertifiedTriangle, Error> {
        let [rs, rz] = slopes(cell.patch, self.k);
        work.charge(1024)?;
        let ids = local.map(|j| cell.ids[j]);
        let points = local.map(|j| cell.parameters[j]);
        positive_area(points)?;
        let box_: [I; 2] = std::array::from_fn(|j| {
            I::new(
                points.iter().map(|p| p[j].lo).fold(f64::INFINITY, f64::min),
                points
                    .iter()
                    .map(|p| p[j].hi)
                    .fold(f64::NEG_INFINITY, f64::max),
            )
        });
        let radius = I::new(
            ids.iter()
                .map(|id| self.plan.vertices[*id as usize].radius.lo)
                .fold(f64::INFINITY, f64::min),
            ids.iter()
                .map(|id| self.plan.vertices[*id as usize].radius.hi)
                .fold(f64::NEG_INFINITY, f64::max),
        );
        let ds = I::point(box_[0].hi) - I::point(box_[0].lo);
        let dz = I::point(box_[1].hi) - I::point(box_[1].lo);
        let hss = self.omega.square() * radius + I::point(2.) * self.omega * absolute(rs);
        let hsz = self.omega * absolute(rz);
        let error = ids
            .iter()
            .map(|id| self.errors[*id as usize])
            .fold(0., f64::max);
        let distance = ((hss * ds.square() + I::point(2.) * hsz * ds * dz) / I::point(2.)
            + I::point(error))
        .hi;
        let center = normal(
            I::point(box_[0].lo + (box_[0].hi - box_[0].lo) * 0.5),
            I::point(radius.lo + (radius.hi - radius.lo) * 0.5),
            rs,
            rz,
            self.omega,
        );
        let whole = normal(box_[0], radius, rs, rz, self.omega);
        let spread = norm(difference(whole, center));
        let directed = [ids[0], ids[2], ids[1]];
        let n = facet(self.mesh, directed)?;
        let angle = (norm(difference(n, center)) + spread).hi;
        if distance > self.quality.deflection || angle > self.allowed {
            return Err(Error::Unresolved);
        }

        Ok(CertifiedTriangle {
            ids: directed,
            distance,
            normal: angle,
        })
    }
}
impl StockClipCellPlan {
    /// Checks one fixed source-owned cell plan. Connectivity/source quality are
    /// earned here; whole literal embedding remains a separate unchanged gate.
    pub fn mesh_candidate(
        &self,
        quality: Tessellation,
        maximum_work: usize,
        maximum_bytes: usize,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<StockClipMeshCandidate, Error> {
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
        let count = self.vertices.len() + 2;
        let triangles = self.cells.iter().map(|c| c.count - 2).sum::<usize>() + self.caps.len();
        if count > 131_072 || triangles > quality.triangle_limit.min(4_000_000) {
            return Err(Error::Limit);
        }
        // Exact indexed certifier peak: 3T*[u32;3] edges +3T corner
        // union words +T shell words, conservatively plusV link words.
        work.memory(
            self.bytes + 24 * count + 16 * triangles + 52 * triangles + 4 * count + SCRATCH,
        )?;
        work.charge(count * 32 + triangles * 32)?;
        let mut mesh = DisplayMesh {
            vertices: Vec::new(),
            indices: Vec::new(),
            face_indices: Vec::new(),
        };
        mesh.vertices
            .try_reserve_exact(3 * count)
            .map_err(|_| Error::Limit)?;
        mesh.indices
            .try_reserve_exact(3 * triangles)
            .map_err(|_| Error::Limit)?;
        mesh.face_indices
            .try_reserve_exact(triangles)
            .map_err(|_| Error::Limit)?;
        // Pinned RawVec exact non-ZST capacities equal the requested
        // layouts. The earlier full peak preflight dominates these actual
        // observations, including the subsequently allocated error buffer.
        work.observe(
            self.bytes
                + mesh.vertices.capacity() * 8
                + mesh.indices.capacity() * 4
                + mesh.face_indices.capacity() * 4
                + 52 * triangles
                + 4 * count
                + SCRATCH,
        );
        let mut errors = Vec::new();
        let error_base = work.peak;
        work.memory(error_base + 8 * count)?;
        errors.try_reserve_exact(count).map_err(|_| Error::Limit)?;
        work.observe(error_base + 8 * errors.capacity());
        let mut maximum_stock_radius: f64 = 0.;
        for source in &self.vertices {
            work.charge(256)?;
            let mut literal = source.position.map(|p| p.lo + (p.hi - p.lo) * 0.5);
            // This is a literal approximation choice inside the earned
            // finite axial source range, not a changed source event/datum.
            // Every displacement remains in the full corner error below.
            literal[2] = literal[2].clamp(0., self.source.source.spec.length_mm);
            let radial = (I::point(literal[0]).square() + I::point(literal[1]).square()).sqrt();
            if radial.hi > self.source.stock_radius_limit_mm() {
                let factor = (I::point(self.source.stock_radius_limit_mm()) / radial).lo
                    * (1. - 16. * f64::EPSILON);
                literal[0] *= factor;
                literal[1] *= factor;
            }
            maximum_stock_radius = maximum_stock_radius
                .max(stock_radius(literal, self.source.stock_radius_limit_mm())?);
            errors.push(
                norm(std::array::from_fn(|i| {
                    I::point(literal[i]) - source.position[i]
                }))
                .hi,
            );
            mesh.vertices.extend(literal);
        }
        mesh.vertices
            .extend([0., 0., 0., 0., 0., self.source.source.spec.length_mm]);
        let omega = I::point(2.) * pi() / I::point(self.source.source.spec.pitch_mm);
        let k = I::point(self.source.source.spec.radial_taper);
        let allowed = (I::point(2.) * (I::point(quality.angle) / I::point(2.)).sin()).lo;
        let mut chord: f64 = 0.;
        let mut normal_error: f64 = 0.;
        for cell in &self.cells {
            let mut root = 0;
            let mut best = -1.;
            for candidate in 0..cell.count {
                let mut worst = f64::INFINITY;
                for i in 1..cell.count - 1 {
                    let ids = [
                        candidate,
                        (candidate + i) % cell.count,
                        (candidate + i + 1) % cell.count,
                    ]
                    .map(|j| cell.ids[j]);
                    let [a, b, c] = ids.map(|id| {
                        std::array::from_fn::<_, 3, _>(|j| mesh.vertices[id as usize * 3 + j])
                    });
                    let ab = std::array::from_fn::<_, 3, _>(|j| b[j] - a[j]);
                    let ac = std::array::from_fn::<_, 3, _>(|j| c[j] - a[j]);
                    let area = [
                        ab[1] * ac[2] - ab[2] * ac[1],
                        ab[2] * ac[0] - ab[0] * ac[2],
                        ab[0] * ac[1] - ab[1] * ac[0],
                    ]
                    .into_iter()
                    .map(|x| x * x)
                    .sum::<f64>()
                    .sqrt();
                    let denominator = (0..3)
                        .map(|j| ab[j] * ab[j] + ac[j] * ac[j] + (b[j] - c[j]) * (b[j] - c[j]))
                        .sum::<f64>();
                    worst = worst.min(area / denominator);
                }
                if worst > best {
                    best = worst;
                    root = candidate;
                }
            }
            work.charge(cell.count * cell.count * 64)?;
            // The aspect winner remains the first proposal, preserving every
            // old successful literal mesh. Skinny clipped cylinder cells can
            // have a good aspect score but a genuinely tilted chord facet.
            // Certify an ENTIRE fan before publishing any of its triangles.
            let context = TriangleContext {
                plan: self,
                mesh: &mesh,
                errors: &errors,
                omega,
                k,
                quality,
                allowed,
            };
            let mut chosen = None;
            for offset in 0..cell.count {
                let root = (root + offset) % cell.count;
                let mut fan = [CertifiedTriangle {
                    ids: [0; 3],
                    distance: 0.,
                    normal: 0.,
                }; 6];
                let mut valid = true;
                for i in 1..cell.count - 1 {
                    let local = [root, (root + i) % cell.count, (root + i + 1) % cell.count];
                    match context.certify(cell, local, &mut work) {
                        Ok(triangle) => fan[i - 1] = triangle,
                        Err(Error::Unresolved) => {
                            valid = false;
                            break;
                        }
                        Err(cause) => return Err(cause),
                    }
                }
                if valid {
                    chosen = Some(fan);
                    break;
                }
            }
            let Some(fan) = chosen else {
                return Err(Error::Unresolved);
            };
            for triangle in fan.into_iter().take(cell.count - 2) {
                chord = chord.max(triangle.distance);
                normal_error = normal_error.max(triangle.normal);
                push(&mut mesh, triangle.ids, face(cell.patch));
            }
        }
        for edge in &self.caps {
            work.charge(512)?;
            let start = edge.station == StockClipStation::Start;
            let distance = cap_edge(self, edge, &mesh, quality)?;
            chord = chord.max(distance);
            let center = if start { count - 2 } else { count - 1 };
            let ids = [center as u32, edge.ids[0], edge.ids[1]];
            let angle = cap_normal(&mesh, ids, start, allowed)?;
            normal_error = normal_error.max(angle);
            push(&mut mesh, ids, if start { 5 } else { 6 });
        }
        work.charge(triangles * (128 + 3 * (triangles.max(1).ilog2() as usize + 1)))?;
        let proof = manufacturing::certify(&mesh).map_err(|_| Error::Unresolved)?;
        let mut volume = I::point(0.);
        for ids in mesh.indices.as_chunks::<3>().0 {
            work.charge(64)?;
            let [a, b, c] = std::array::from_fn(|i| at(&mesh, ids[i]));
            let area = cross(b, c);
            volume = volume + (a[0] * area[0] + a[1] * area[1] + a[2] * area[2]) / I::point(6.);
        }
        let source_volume = self.source.volume_bounds(
            (quality.deflection * self.source.radius_range.hi * self.source.radius_range.hi)
                .min(f64::MAX),
            work.query()?,
            interrupt,
        )?;
        work.charge(source_volume.work_used())?;
        work.charge(1)?;
        Ok(StockClipMeshCandidate {
            mesh,
            proof,
            stock_radius: maximum_stock_radius,
            chord,
            normal: normal_error,
            source_volume: source_volume.source_volume_mm3(),
            literal_volume: bounds(volume),
            work: work.used,
            bytes: work.peak,
        })
    }
}

// This numerical stage also accepts unqualified literal coordinates. It
// returns only its local source-transfer bound, not a valid-body receipt.
fn cap_edge(
    plan: &StockClipCellPlan,
    edge: &CapEdge,
    mesh: &DisplayMesh,
    quality: Tessellation,
) -> Result<f64, Error> {
    let omega = I::point(2.) * pi() / I::point(plan.source.source.spec.pitch_mm);
    let k = I::point(plan.source.source.spec.radial_taper);
    let [a, b] = edge.ids.map(|id| at(mesh, id));
    if a.into_iter()
        .chain(b)
        .any(|p| !p.lo.is_finite() || !p.hi.is_finite())
    {
        return Err(Error::Unresolved);
    }
    let start = edge.station == StockClipStation::Start;
    let cross = a[0] * b[1] - a[1] * b[0];
    if (start && cross.hi >= 0.) || (!start && cross.lo <= 0.) {
        return Err(Error::Unresolved);
    }
    let radius = I::new(
        edge.ids
            .iter()
            .map(|id| plan.vertices[*id as usize].radius.lo)
            .fold(f64::INFINITY, f64::min),
        edge.ids
            .iter()
            .map(|id| plan.vertices[*id as usize].radius.hi)
            .fold(f64::NEG_INFINITY, f64::max),
    );
    let sr = I::new(
        edge.parameters[0][0].lo.min(edge.parameters[1][0].lo),
        edge.parameters[0][0].hi.max(edge.parameters[1][0].hi),
    );
    let theta = -omega * sr;
    // Strict common halfplane supplies a nonzero boundary-to-chord homotopy.
    let source = [radius * theta.cos(), radius * theta.sin()];
    let halfplane = source[0] * (a[0] + b[0]) + source[1] * (a[1] + b[1]);
    if halfplane.lo <= 0.
        || (a[0] * (a[0] + b[0]) + a[1] * (a[1] + b[1])).lo <= 0.
        || (b[0] * (a[0] + b[0]) + b[1] * (a[1] + b[1])).lo <= 0.
    {
        return Err(Error::Unresolved);
    }
    let rs = slopes(edge.patch, k)[0];
    let hss = omega.square() * radius + I::point(2.) * omega * absolute(rs);
    let ds = I::point(sr.hi) - I::point(sr.lo);
    let error = edge
        .ids
        .iter()
        .map(|id| {
            norm(difference(
                at(mesh, *id),
                plan.vertices[*id as usize].position,
            ))
            .hi
        })
        .fold(0., f64::max);
    let distance = (I::point(3.) * (hss * ds.square() / I::point(8.) + I::point(error))).hi;
    if distance > quality.deflection {
        return Err(Error::Unresolved);
    }
    Ok(distance)
}
fn cap_normal(mesh: &DisplayMesh, ids: [u32; 3], start: bool, allowed: f64) -> Result<f64, Error> {
    let n = facet(mesh, ids)?;
    let expected = [
        I::point(0.),
        I::point(0.),
        I::point(if start { -1. } else { 1. }),
    ];
    let angle = norm(difference(n, expected)).hi;
    if angle > allowed {
        return Err(Error::Unresolved);
    }
    Ok(angle)
}

#[cfg(test)]
mod tests;
