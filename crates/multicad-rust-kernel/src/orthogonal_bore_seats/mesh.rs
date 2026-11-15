use super::*;
use crate::display::{DisplayMesh, MeshError, Tessellation, bored_orthogonal};
fn failure(e: MeshError) -> Error {
    match e {
        MeshError::InvalidOptions => Error::Invalid,
        MeshError::BudgetExceeded => Error::Limit,
        MeshError::UnresolvedBoundary => Error::Unsupported,
    }
}
use std::collections::{BTreeMap, BTreeSet};
fn vertex(mesh: &DisplayMesh, i: u32) -> [f64; 3] {
    let i = i as usize * 3;
    [mesh.vertices[i], mesh.vertices[i + 1], mesh.vertices[i + 2]]
}
fn memory(mesh: &DisplayMesh) -> usize {
    mesh.vertices.capacity() * std::mem::size_of::<f64>()
        + mesh.indices.capacity() * std::mem::size_of::<u32>()
        + mesh.face_indices.capacity() * std::mem::size_of::<u32>()
}
impl Seats {
    /// Mesh only this retained material. Original source face ordinals remain
    /// unchanged; one conical face follows all original cylindrical walls.
    pub fn mesh(
        &self,
        options: Tessellation,
        interrupt: &dyn Fn() -> bool,
    ) -> Result<DisplayMesh, Error> {
        options.check().map_err(failure)?;
        let cancelled = || checkpoint(interrupt);
        cancelled()?;
        let maximum_radius = self.seats.iter().map(|s| s.radius_mm).fold(0_f64, f64::max);
        let minimum_radius = self
            .bores
            .iter()
            .map(|b| b.radius_mm)
            .fold(f64::INFINITY, f64::min);
        let tightened = Tessellation {
            deflection: options.deflection * (minimum_radius / maximum_radius).min(1.) * 0.25,
            triangle_limit: options.triangle_limit.min(65_536),
            ..options
        };
        let mut mesh =
            bored_orthogonal(&self.boundary, &[], &self.bores, tightened).map_err(failure)?;
        let source_triangles = mesh.face_indices.len();
        for (seat_index, seat) in self.seats.iter().enumerate() {
            cancelled()?;
            let bore = self.bores[seat.bore];
            let end = usize::from(seat.end);
            let center = bore.centers_mm[end];
            let face = &self.boundary.faces[bore.faces[end]];
            let axis = (0..3)
                .find(|&i| face.normal[i] != 0.)
                .ok_or(Error::Unsupported)?;
            let wall = (self.boundary.faces.len() + seat.bore) as u32;
            let mut edges = Vec::new();
            let mut rim = BTreeSet::new();
            for (triangle, owner) in mesh.indices[..3 * source_triangles]
                .as_chunks::<3>()
                .0
                .iter()
                .zip(&mesh.face_indices[..source_triangles])
            {
                if *owner != wall {
                    continue;
                }
                for j in 0..3 {
                    let a = triangle[j];
                    let b = triangle[(j + 1) % 3];
                    if vertex(&mesh, a)[axis] == center[axis]
                        && vertex(&mesh, b)[axis] == center[axis]
                    {
                        edges.push([a, b]);
                        rim.insert(a);
                        rim.insert(b);
                    }
                }
            }
            // One complete original rim, with a unique directed edge at each
            // vertex. No coordinate welding, sorting-by-angle or new seam.
            if edges.len() < 7 || rim.len() != edges.len() {
                return Err(Error::Unsupported);
            }
            let count = edges.len();
            let mut outgoing = BTreeMap::new();
            let mut incoming = BTreeSet::new();
            for &[a, b] in &edges {
                if outgoing.insert(a, b).is_some() || !incoming.insert(b) {
                    return Err(Error::Unsupported);
                }
            }
            if outgoing.keys().copied().collect::<BTreeSet<_>>() != incoming {
                return Err(Error::Unsupported);
            }
            let first = edges[0][0];
            let mut at = first;
            let mut visited = BTreeSet::new();
            for _ in 0..count {
                if !visited.insert(at) {
                    return Err(Error::Unsupported);
                }
                at = *outgoing.get(&at).ok_or(Error::Unsupported)?;
            }
            if at != first {
                return Err(Error::Unsupported);
            }
            let pi = I {
                lo: PI.next_down(),
                hi: PI.next_up(),
            };
            let n = I::point(count as f64);
            let coordinate = self
                .stock
                .bounds()
                .min
                .into_iter()
                .chain(self.stock.bounds().max)
                .fold(maximum_radius, |a, b| a.max(b.abs()));
            // Scaling an already rounded original circle offset amplifies its
            // subtraction error by R/r. Retain this factor even for tiny bores
            // far from the origin; a fixed epsilon reserve would be unsound.
            let ratio = I::point(seat.radius_mm) / I::point(bore.radius_mm);
            let roundoff = (I::point(4096. * f64::EPSILON)
                * (I::point(coordinate)
                    + I::point(maximum_radius)
                    + I::point(seat.depth_mm)
                    + I::point(1.))
                * (I::point(1.) + ratio))
                .hi;
            // Every conical trapezoid is a homothetic radial sector. At each
            // axial station its radial sag is <=R*delta_angle²/8. The same
            // bound covers the new cap chord. Include input trig/affine reserve.
            let error =
                I::point(seat.radius_mm) * pi * pi / (I::point(2.) * n * n) + I::point(roundoff);
            if !error.hi.is_finite() || error.hi > options.deflection {
                return Err(Error::Unsupported);
            }
            let extra = count * 2;
            if mesh
                .face_indices
                .len()
                .checked_add(extra)
                .is_none_or(|n| n > options.triangle_limit.min(65_536))
            {
                return Err(Error::Limit);
            }
            if self.accounted_peak_bytes + memory(&mesh) + count * 512 > 64 * 1024 * 1024 {
                return Err(Error::Limit);
            }
            mesh.vertices
                .try_reserve_exact(count * 3)
                .map_err(|_| Error::Limit)?;
            mesh.indices
                .try_reserve_exact(extra * 3)
                .map_err(|_| Error::Limit)?;
            mesh.face_indices
                .try_reserve_exact(extra)
                .map_err(|_| Error::Limit)?;
            let mut inner = BTreeMap::new();
            for old in rim {
                let original = vertex(&mesh, old);
                let new = mesh.vertices.len() / 3;
                let mut p = original;
                p[axis] -= face.normal[axis] * seat.depth_mm;
                mesh.vertices.extend(p);
                inner.insert(old, new as u32);
                for i in 0..3 {
                    if i != axis {
                        mesh.vertices[old as usize * 3 + i] = center[i]
                            + (original[i] - center[i]) * (seat.radius_mm / bore.radius_mm);
                    }
                }
            }
            for (triangle, owner) in mesh.indices[..3 * source_triangles]
                .as_chunks_mut::<3>()
                .0
                .iter_mut()
                .zip(&mesh.face_indices[..source_triangles])
            {
                if *owner == wall {
                    for p in triangle {
                        if let Some(new) = inner.get(p) {
                            *p = *new;
                        }
                    }
                }
            }
            let cone = (self.boundary.faces.len() + self.bores.len() + seat_index) as u32;
            for [a, b] in edges {
                let [c, d] = [inner[&a], inner[&b]];
                mesh.indices.extend([a, b, d, a, d, c]);
                mesh.face_indices.extend([cone, cone]);
            }
        }
        // Relocating a mouth may conservatively invalidate an old constrained
        // planar triangulation. Prove every resulting stock triangle retains
        // positive orientation against its ORIGINAL support; never drop it.
        for (triangle, &owner) in mesh
            .indices
            .as_chunks::<3>()
            .0
            .iter()
            .zip(&mesh.face_indices)
        {
            cancelled()?;
            if let Some(face) = self.boundary.faces.get(owner as usize) {
                let [a, b, c] = [
                    vertex(&mesh, triangle[0]),
                    vertex(&mesh, triangle[1]),
                    vertex(&mesh, triangle[2]),
                ];
                let u = std::array::from_fn::<_, 3, _>(|i| I::point(b[i]) - I::point(a[i]));
                let v = std::array::from_fn::<_, 3, _>(|i| I::point(c[i]) - I::point(a[i]));
                let cross = [
                    u[1] * v[2] - u[2] * v[1],
                    u[2] * v[0] - u[0] * v[2],
                    u[0] * v[1] - u[1] * v[0],
                ];
                let orientation =
                    (0..3).fold(I::point(0.), |s, i| s + cross[i] * I::point(face.normal[i]));
                if orientation.lo <= 0. {
                    return Err(Error::Unsupported);
                }
            }
        }
        if self.accounted_peak_bytes + memory(&mesh) > 64 * 1024 * 1024 {
            return Err(Error::Limit);
        }
        cancelled()?;
        Ok(mesh)
    }
}
