//! Independent literal-f64 embedding, after source correspondence and rounding.
use super::*;
use std::collections::BTreeSet;

#[derive(Clone)]
pub(super) struct Cap {
    pub face: usize,
    pub apex: u32,
    /// Actual outward fan winding (bottom clockwise, top counter-clockwise).
    pub rim: Vec<u32>,
    pub edges: Vec<LiteralCapEdge>,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LiteralEmbeddingProof {
    connectivity: crate::manufacturing::MeshProof,
    radial_lower: f64,
    height_clearance: f64,
}
impl LiteralEmbeddingProof {
    pub fn connectivity(self) -> crate::manufacturing::MeshProof {
        self.connectivity
    }
    pub fn minimum_radial_determinant(self) -> f64 {
        self.radial_lower
    }
    pub fn minimum_axis_height_clearance_mm(self) -> f64 {
        self.height_clearance
    }
}
fn root(parents: &mut [u32], mut node: u32, budget: &mut Budget<'_>) -> Result<u32> {
    let mut result = node;
    while parents[result as usize] != result {
        budget.step(16)?;
        result = parents[result as usize];
    }
    while node != result {
        budget.step(16)?;
        let next = parents[node as usize];
        parents[node as usize] = result;
        node = next;
    }
    Ok(result)
}
fn connected_sphere(
    mesh: &DisplayMesh,
    proof: crate::manufacturing::MeshProof,
    budget: &mut Budget<'_>,
) -> Result<()> {
    let count = mesh.vertices.len() / 3;
    if proof.vertex_links as usize != count
        || count + mesh.face_indices.len() != proof.edge_pairs as usize + 2
    {
        return Err(Error::Unsupported);
    }
    let mut parents = budget.allocate(count)?;
    parents.extend(0..count as u32);
    let mut ranks = budget.allocate(count)?;
    ranks.resize(count, 0u8);
    for triangle in mesh.indices.as_chunks::<3>().0 {
        for index in [1, 2] {
            budget.step(32)?;
            let mut a = root(&mut parents, triangle[0], budget)?;
            let mut b = root(&mut parents, triangle[index], budget)?;
            if a != b {
                if ranks[a as usize] < ranks[b as usize] {
                    std::mem::swap(&mut a, &mut b);
                }
                parents[b as usize] = a;
                if ranks[a as usize] == ranks[b as usize] {
                    ranks[a as usize] += 1;
                }
            }
        }
    }
    let first = root(&mut parents, 0, budget)?;
    for index in 1..count {
        budget.step(16)?;
        if root(&mut parents, index as u32, budget)? != first {
            return Err(Error::Unsupported);
        }
    }
    budget.release(&parents);
    budget.release(&ranks);
    Ok(())
}
fn cross(a: [I; 3], b: [I; 3]) -> [I; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn dot(a: [I; 3], b: [I; 3]) -> I {
    a.into_iter()
        .zip(b)
        .fold(I::point(0.), |sum, (a, b)| sum + a * b)
}
fn sub(a: [I; 3], b: [I; 3]) -> [I; 3] {
    std::array::from_fn(|k| a[k] - b[k])
}
fn vertex(mesh: &DisplayMesh, index: u32) -> Result<[I; 3]> {
    let i = index as usize;
    let p = mesh
        .vertices
        .get(i * 3..i * 3 + 3)
        .ok_or(Error::Unsupported)?;
    if p.iter().any(|v| !v.is_finite()) {
        return Err(Error::Unsupported);
    }
    Ok(std::array::from_fn(|k| I::point(p[k])))
}
fn unit(vector: [I; 3]) -> Result<([I; 3], I)> {
    let length = dot(vector, vector).sqrt();
    if !length.valid() || length.lo <= 0. {
        return Err(Error::Unsupported);
    }
    Ok((vector.map(|v| v / length), length))
}

pub(super) fn certify(
    mesh: &DisplayMesh,
    caps: &[Cap; 2],
    budget: &mut Budget<'_>,
) -> Result<LiteralEmbeddingProof> {
    if caps[0].apex == caps[1].apex
        || caps[0].face == caps[1].face
        || caps.iter().any(|cap| cap.rim.len() < 3)
    {
        return Err(Error::Unsupported);
    }
    let triangles = mesh.face_indices.len();
    let scratch = triangles
        .checked_mul(52)
        .and_then(|n| n.checked_add(mesh.vertices.len() / 3 * 4))
        .ok_or(Error::Limit)?;
    hold(budget, scratch)?;
    let comparisons = (mesh.indices.len().max(1).ilog2() as usize + 1) * 32;
    budget.step(mesh.indices.len() * (128 + comparisons))?;
    let connectivity = crate::manufacturing::certify(mesh).map_err(|_| Error::Unsupported)?;
    budget.bytes -= scratch;
    connected_sphere(mesh, connectivity, budget)?;
    let origins = [vertex(mesh, caps[0].apex)?, vertex(mesh, caps[1].apex)?];
    let (axis, length) = unit(sub(origins[1], origins[0]))?;
    let reference = (0..3)
        .min_by(|&a, &b| axis[a].magnitude().total_cmp(&axis[b].magnitude()))
        .expect("three literal coordinate axes");
    let reference = std::array::from_fn(|k| I::point(if k == reference { 1. } else { 0. }));
    let (v, _) = unit(cross(axis, reference))?;
    let u = cross(v, axis);
    let mut cap_edges = [BTreeSet::new(), BTreeSet::new()];
    let entries = caps.iter().map(|cap| cap.rim.len()).sum::<usize>();
    let set_bytes = entries
        .checked_mul(3 * size_of::<(u32, u32)>() + 64)
        .and_then(|n| n.checked_add(8192))
        .ok_or(Error::Limit)?;
    hold(budget, set_bytes)?;
    let mut bottom_max = f64::NEG_INFINITY;
    let mut top_min = f64::INFINITY;
    let mut clearance = f64::INFINITY;
    for (end, cap) in caps.iter().enumerate() {
        let mut angle = I::point(0.);
        for index in 0..cap.rim.len() {
            budget.step(512)?;
            let a = cap.rim[index];
            let b = cap.rim[(index + 1) % cap.rim.len()];
            if a == cap.apex || b == cap.apex || !cap_edges[end].insert((a, b)) {
                return Err(Error::Unsupported);
            }
            let pa = sub(vertex(mesh, a)?, origins[0]);
            let pb = sub(vertex(mesh, b)?, origins[0]);
            let height = dot(pa, axis);
            if !height.valid() || height.lo <= 0. || height.hi >= length.lo {
                return Err(Error::Unsupported);
            }
            clearance = clearance.min(height.lo).min((length - height).lo);
            if end == 0 {
                bottom_max = bottom_max.max(height.hi);
            } else {
                top_min = top_min.min(height.lo);
            }
            let a = [dot(pa, u), dot(pa, v)];
            let b = [dot(pb, u), dot(pb, v)];
            let delta = crate::polar_loop11::positive_chord_angle(a, b, end == 1)
                .ok_or(Error::Unsupported)?;
            angle = angle + delta;
        }
        // Positive literal chord increments and exact indexed closure make
        // the angular sum an integer multiple of2pi. This interval permits1.
        if !angle.valid() || angle.lo <= 4. || angle.hi >= 8. {
            return Err(Error::Unsupported);
        }
    }
    if bottom_max >= top_min {
        return Err(Error::Unsupported);
    }
    let mut radial_lower = f64::INFINITY;
    for (triangle, face) in mesh
        .indices
        .as_chunks::<3>()
        .0
        .iter()
        .zip(&mesh.face_indices)
    {
        budget.step(512)?;
        let cap = caps.iter().position(|cap| cap.face == *face as usize);
        if let Some(end) = cap {
            let apex = caps[end].apex;
            let start = triangle
                .iter()
                .position(|v| *v == apex)
                .ok_or(Error::Unsupported)?;
            let edge = (triangle[(start + 1) % 3], triangle[(start + 2) % 3]);
            if !cap_edges[end].remove(&edge) {
                return Err(Error::Unsupported);
            }
        } else if triangle
            .iter()
            .any(|v| caps.iter().any(|cap| cap.apex == *v))
        {
            return Err(Error::Unsupported);
        }
        let points = [
            vertex(mesh, triangle[0])?,
            vertex(mesh, triangle[1])?,
            vertex(mesh, triangle[2])?,
        ];
        let normal = cross(sub(points[1], points[0]), sub(points[2], points[0]));
        for (&index, point) in triangle.iter().zip(points) {
            if caps.iter().any(|cap| cap.apex == index) {
                continue;
            }
            let d = sub(point, origins[0]);
            let z = dot(d, axis);
            if !z.valid() || z.lo <= 0. || z.hi >= length.lo {
                return Err(Error::Unsupported);
            }
            clearance = clearance.min(z.lo).min((length - z).lo);
            let radial = std::array::from_fn(|k| d[k] - axis[k] * z);
            let determinant = dot(normal, radial);
            if !determinant.valid() || determinant.lo <= 0. {
                return Err(Error::Unsupported);
            }
            radial_lower = radial_lower.min(determinant.lo);
        }
    }
    if cap_edges.iter().any(|edges| !edges.is_empty()) {
        return Err(Error::Unsupported);
    }
    budget.bytes -= set_bytes;
    budget.step(1)?;
    // Actual indexed vertex links, connectedness and Euler2 give an oriented
    // sphere. The two disjoint, completely matched cap fans are disks; their
    // removal is an annulus. Removing its
    // two verified fan apices yields a proper map to the open(theta,z)
    // cylinder: every other vertex has0<z<L, and every COMPLETE triangle has
    // positive radial projection (an affine determinant; at a fan apex it is
    // exactly0 and positive everywhere else). Each fan has exact indexed
    // closure and angular winding1, so the global degree is1. At seams,
    // opposite indexed coedges put adjacent positive charts on opposite
    // sides. A vertex fan can only have positive integral local degree;
    // an extra turn or a second preimage would contradict global degree1.
    // Consequently the entire punctured surface is a one-sheet radial graph.
    // The fan apices close its two ends. This excludes self-intersection of
    // the literal rounded coordinates, without source-normal inference,
    // tolerance welding, a sampled mesh test or a nominal body recipe.
    Ok(LiteralEmbeddingProof {
        connectivity,
        radial_lower,
        height_clearance: clearance,
    })
}

#[cfg(test)]
mod tests;
