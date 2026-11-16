//! Certify the exact indexed connectivity of a fresh analytic tessellation.
//!
//! This check preserves vertex identities and the input winding. It proves
//! closed, oriented manifold connectivity and measures signed material
//! volume. The analytic tessellator establishes the geometric embedding and
//! chord tolerance before certification.

use crate::display::{DisplayMesh, MeshError};
mod counting_sort;
mod fan_cycles;

const MAX_TRIANGLES: usize = 5_000_000;
const MAX_VERTICES: usize = MAX_TRIANGLES * 3;
const INVALID: MeshError = MeshError::UnresolvedBoundary;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MeshProof {
    pub edge_pairs: u32,
    pub vertex_links: u32,
    pub signed_volume_mm3: f64,
}

// Path compression keeps the corner and shell unions compact. Numbering
// roots by the lowest input index makes local volume origins deterministic.
struct Sets(Vec<u32>);

impl Sets {
    fn new(count: usize) -> Self {
        Self((0..count as u32).collect())
    }

    fn root<E, const COUNTED: bool>(
        &mut self,
        node: u32,
        budget: &mut Budget<E, COUNTED>,
    ) -> Result<u32, CountedError<E>> {
        budget.work(1)?;
        let mut root = node;
        while self.0[root as usize] != root {
            budget.work(1)?;
            root = self.0[root as usize];
        }
        let mut node = node;
        while node != root {
            budget.work(1)?;
            let parent = self.0[node as usize];
            self.0[node as usize] = root;
            node = parent;
        }
        Ok(root)
    }

    fn join<E, const COUNTED: bool>(
        &mut self,
        a: u32,
        b: u32,
        budget: &mut Budget<E, COUNTED>,
    ) -> Result<(), CountedError<E>> {
        let (a, b) = (self.root(a, budget)?, self.root(b, budget)?);
        budget.work(1)?;
        self.0[a.max(b) as usize] = a.min(b);
        Ok(())
    }
}

fn shape_budget(vertices: usize, triangles: usize) -> Result<(), MeshError> {
    if vertices > MAX_VERTICES || triangles > MAX_TRIANGLES {
        Err(MeshError::BudgetExceeded)
    } else {
        Ok(())
    }
}

fn next(corner: u32) -> u32 {
    corner / 3 * 3 + (corner % 3 + 1) % 3
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|i| a[i] - b[i])
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Check indexed triangles without changing their coordinates or topology.
/// Distinct closed shells may have separate index identities at coincident
/// points. Inward cavity shells contribute negative signed volume.
pub fn certify(mesh: &DisplayMesh) -> Result<MeshProof, MeshError> {
    let mut noop = |_: usize, _: Option<usize>| Ok::<_, std::convert::Infallible>(());
    match certify_impl::<_, false>(mesh, &mut Budget(&mut noop)) {
        Ok(proof) => Ok(proof),
        Err(CountedError::Mesh(error)) => Err(error),
        Err(CountedError::Interrupted(never)) => match never {},
    }
}

/// Private budgeted form of the same exact indexed-connectivity algorithm.
/// The callback charges work, or checks total live certification scratch bytes,
/// before the corresponding operation/allocation. It never authorizes geometry.
pub(crate) fn certify_counted<E>(
    mesh: &DisplayMesh,
    callback: &mut dyn FnMut(usize, Option<usize>) -> Result<(), E>,
) -> Result<MeshProof, CountedError<E>> {
    certify_impl::<_, true>(mesh, &mut Budget(callback))
}
#[derive(Debug, PartialEq, Eq)]
pub(crate) enum CountedError<E> {
    Mesh(MeshError),
    Interrupted(E),
}
impl<E> From<MeshError> for CountedError<E> {
    fn from(error: MeshError) -> Self {
        Self::Mesh(error)
    }
}
struct Budget<'a, E, const COUNTED: bool>(&'a mut dyn FnMut(usize, Option<usize>) -> Result<(), E>);
impl<E, const COUNTED: bool> Budget<'_, E, COUNTED> {
    #[inline(always)]
    fn work(&mut self, n: usize) -> Result<(), CountedError<E>> {
        if !COUNTED {
            return Ok(());
        }
        (self.0)(n, None).map_err(CountedError::Interrupted)
    }
    #[inline(always)]
    fn storage(&mut self, bytes: usize) -> Result<(), CountedError<E>> {
        if !COUNTED {
            return Ok(());
        }
        let bytes = bytes
            .checked_add(4 * std::mem::size_of::<Vec<u32>>())
            .ok_or(MeshError::BudgetExceeded)?;
        (self.0)(0, Some(bytes)).map_err(CountedError::Interrupted)
    }
}
fn sort_counted<E, const COUNTED: bool>(
    edges: &mut Vec<[u32; 3]>,
    vertices: usize,
    budget: &mut Budget<E, COUNTED>,
) -> Result<(), CountedError<E>> {
    // certify_impl has validated every endpoint and generated strictly
    // increasing coedge ordinals. Stable b/a passes therefore also retain
    // the full third-key ordering. Avoid vertex-sized storage for sparse data.
    if vertices <= edges.len() {
        return counting_sort::sort(edges, vertices, budget);
    }
    // Both buffers have exactly the original edge count. Every comparison and
    // literal edge-record copy is charged; cancellation returns normally.
    budget.storage(
        edges
            .capacity()
            .checked_add(edges.len())
            .and_then(|n| n.checked_mul(12))
            .ok_or(MeshError::BudgetExceeded)?,
    )?;
    budget.work(edges.len())?;
    let mut scratch = edges.to_vec();
    budget.storage(
        (edges.capacity() + scratch.capacity())
            .checked_mul(12)
            .ok_or(MeshError::BudgetExceeded)?,
    )?;
    let mut source = edges.as_mut_slice();
    let mut target = scratch.as_mut_slice();
    let mut flipped = false;
    let mut width = 1usize;
    while width < source.len() {
        for start in (0..source.len()).step_by(width * 2) {
            let middle = (start + width).min(source.len());
            let end = (middle + width).min(source.len());
            let (mut left, mut right) = (start, middle);
            for destination in &mut target[start..end] {
                let take_left = if left == middle {
                    false
                } else if right == end {
                    true
                } else {
                    budget.work(1)?;
                    source[left] <= source[right]
                };
                budget.work(1)?;
                *destination = if take_left {
                    let value = source[left];
                    left += 1;
                    value
                } else {
                    let value = source[right];
                    right += 1;
                    value
                };
            }
        }
        std::mem::swap(&mut source, &mut target);
        flipped = !flipped;
        width = width.checked_mul(2).ok_or(MeshError::BudgetExceeded)?;
    }
    if flipped {
        budget.work(source.len())?;
        target.copy_from_slice(source);
    }
    Ok(())
}
fn certify_impl<E, const COUNTED: bool>(
    mesh: &DisplayMesh,
    budget: &mut Budget<E, COUNTED>,
) -> Result<MeshProof, CountedError<E>> {
    let (vertices, vertex_tail) = mesh.vertices.as_chunks::<3>();
    let (triangles, index_tail) = mesh.indices.as_chunks::<3>();
    if !vertex_tail.is_empty()
        || !index_tail.is_empty()
        || triangles.is_empty()
        || mesh.face_indices.len() != triangles.len()
    {
        return Err(INVALID.into());
    }
    shape_budget(vertices.len(), triangles.len())?;
    for coordinate in &mesh.vertices {
        budget.work(1)?;
        if !coordinate.is_finite() {
            return Err(INVALID.into());
        }
    }

    budget.storage(
        mesh.indices
            .len()
            .checked_mul(12)
            .ok_or(MeshError::BudgetExceeded)?,
    )?;
    let mut edges = Vec::with_capacity(mesh.indices.len());
    budget.storage(
        edges
            .capacity()
            .checked_mul(12)
            .ok_or(MeshError::BudgetExceeded)?,
    )?;
    for (face, &triangle) in triangles.iter().enumerate() {
        // Logical vector/record units: three vertex fetches, two differences,
        // one cross product, finite-area and nonzero-area validations (8).
        budget.work(8)?;
        let [Some(a), Some(b), Some(c)] = triangle.map(|i| vertices.get(i as usize)) else {
            return Err(INVALID.into());
        };
        let area = cross(sub(*b, *a), sub(*c, *a));
        if area.iter().any(|c| !c.is_finite()) || area == [0.; 3] {
            return Err(INVALID.into());
        }
        for (i, a) in triangle.into_iter().enumerate() {
            budget.work(1)?;
            let b = triangle[(i + 1) % 3];
            edges.push([a.min(b), a.max(b), (face * 3 + i) as u32]);
        }
    }
    if COUNTED {
        sort_counted(&mut edges, vertices.len(), budget)?;
    } else {
        edges.sort_unstable();
    }
    let graph_bytes = edges
        .capacity()
        .checked_mul(12)
        .and_then(|b| b.checked_add(mesh.indices.len() * 4))
        .and_then(|b| b.checked_add(triangles.len() * 4))
        .ok_or(MeshError::BudgetExceeded)?;
    budget.storage(graph_bytes)?;
    budget.work(mesh.indices.len() + triangles.len())?;

    let mut fans = Sets::new(mesh.indices.len());
    let mut shells = Sets::new(triangles.len());
    budget.storage(edges.capacity() * 12 + (fans.0.capacity() + shells.0.capacity()) * 4)?;
    let mut uses = edges.iter().peekable();
    while let Some([a, b, start]) = uses.next().copied() {
        // Two edge-use fetches, multiplicity and winding validations (4).
        budget.work(4)?;
        let [c, d, opposite] = uses.next().copied().ok_or(INVALID)?;
        if [a, b] != [c, d] || uses.peek().is_some_and(|e| e[..2] == [a, b]) {
            return Err(INVALID.into());
        }
        if mesh.indices[start as usize] == mesh.indices[opposite as usize]
            || mesh.indices[next(next(start)) as usize]
                == mesh.indices[next(next(opposite)) as usize]
        {
            return Err(INVALID.into());
        }
        if COUNTED {
            // Pairing is an involution and next() permutes triangle corners.
            // Their composition permutes ALL corners, preserving the vertex
            // by the opposite-winding test above. Its cycles are precisely
            // the components of the original fan-union graph.
            budget.work(2)?;
            fans.0[start as usize] = next(opposite);
            fans.0[opposite as usize] = next(start);
        } else {
            fans.join(start, next(opposite), budget)?;
            fans.join(next(start), opposite, budget)?;
        }
        shells.join(start / 3, opposite / 3, budget)?;
    }
    let edge_pairs = (edges.len() / 2) as u32;
    drop(edges);

    let vertex_links = if COUNTED {
        fan_cycles::validate(
            &mut fans.0,
            &mesh.indices,
            vertices.len(),
            shells.0.capacity(),
            budget,
        )?
    } else {
        budget.storage(
            (mesh.indices.len() + triangles.len() + vertices.len())
                .checked_mul(4)
                .ok_or(MeshError::BudgetExceeded)?,
        )?;
        budget.work(vertices.len())?;
        let mut links = vec![u32::MAX; vertices.len()];
        budget.storage((fans.0.capacity() + shells.0.capacity() + links.capacity()) * 4)?;
        let mut vertex_links = 0;
        for (corner, &vertex) in mesh.indices.iter().enumerate() {
            budget.work(1)?;
            let root = fans.root(corner as u32, budget)?;
            let link = &mut links[vertex as usize];
            if *link == u32::MAX {
                *link = root;
                vertex_links += 1;
            } else if *link != root {
                return Err(INVALID.into());
            }
        }
        vertex_links
    };
    drop(fans);

    // Recenter each shell before integration so distant placements do not
    // turn volume into a subtraction of large world-coordinate products.
    let mut volume = 0.;
    let mut correction = 0.;
    for (face, &triangle) in triangles.iter().enumerate() {
        // Three vector differences, one cross/dot volume expression and
        // the four compensated accumulation actions (8 logical units).
        budget.work(8)?;
        let shell = shells.root(face as u32, budget)? as usize;
        let origin = vertices[triangles[shell][0] as usize];
        let [a, b, c] = triangle.map(|i| sub(vertices[i as usize], origin));
        let area = cross(b, c);
        let term = (a[0] * area[0] + a[1] * area[1] + a[2] * area[2]) / 6.;
        let adjusted = term - correction;
        let next = volume + adjusted;
        correction = (next - volume) - adjusted;
        volume = next;
    }
    if !volume.is_finite() || volume <= 0. {
        return Err(INVALID.into());
    }
    Ok(MeshProof {
        edge_pairs,
        vertex_links,
        signed_volume_mm3: volume,
    })
}

#[cfg(test)]
mod tests;
