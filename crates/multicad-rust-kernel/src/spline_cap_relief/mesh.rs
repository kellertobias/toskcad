//! Chord-bounded detached display over one shared analytic face graph.
use super::*;
use crate::display::{DisplayMesh, MeshError, Tessellation, triangulate_region};
use spade::{Point2, Triangulation};
use std::collections::BTreeMap;
fn key(point: [f64; 2]) -> [u64; 2] {
    point.map(|v| if v == 0. { 0 } else { v.to_bits() })
}
impl Relief {
    pub fn mesh(&self, mut options: Tessellation) -> Result<DisplayMesh, MeshError> {
        options.check()?;
        options.deflection = remaining_budget(options.deflection, self.boundary_error);
        options.angle = remaining_budget(options.angle, 1e-6);
        if options.deflection <= 0. || options.angle <= 0. {
            return Err(MeshError::UnresolvedBoundary);
        }
        for _ in 0..6 {
            match once(self, options) {
                Err(MeshError::UnresolvedBoundary) => options.deflection *= 0.25,
                result => return result,
            }
        }
        Err(MeshError::UnresolvedBoundary)
    }
}
fn insert_cap_point(
    ids: &mut BTreeMap<[u64; 2], usize>,
    point: [f64; 2],
    id: usize,
) -> Result<(), MeshError> {
    if ids.insert(key(point), id).is_some() {
        return Err(MeshError::UnresolvedBoundary);
    }
    Ok(())
}

fn once(relief: &Relief, options: Tessellation) -> Result<DisplayMesh, MeshError> {
    let graph = relief.boundary();
    if graph.vertices.len() > 3 * options.triangle_limit {
        return Err(MeshError::BudgetExceeded);
    }
    let mut mesh = DisplayMesh {
        vertices: graph.vertices.iter().flatten().copied().collect(),
        indices: vec![],
        face_indices: vec![],
    };
    let mut samples = Vec::with_capacity(graph.edges.len());
    for edge in &graph.edges {
        let points = match &edge.curve {
            Curve::Profile { edge, z_mm } => edge
                .sample(options)?
                .into_iter()
                .map(|p| [p[0], p[1], *z_mm])
                .collect::<Vec<_>>(),
            Curve::Line { a_mm, b_mm } => vec![*a_mm, *b_mm],
        };
        let count = points.len() - 1;
        if mesh.vertices.len() / 3 + count > 3 * options.triangle_limit {
            return Err(MeshError::BudgetExceeded);
        }
        let mut ring = vec![edge.vertices[0]];
        for point in points.iter().take(count).skip(1) {
            let id = mesh.vertices.len() / 3;
            mesh.vertices.extend(*point);
            ring.push(id);
        }
        ring.push(edge.vertices[1]);
        samples.push(ring);
    }
    for (ordinal, face) in graph.faces.iter().enumerate() {
        match &face.surface {
            Surface::Plane { z_mm, upward } => {
                let origin_mm = [0., 0., *z_mm];
                let u = [1., 0., 0.];
                let v = [0., if *upward { 1. } else { -1. }, 0.];
                let mut points = vec![];
                let mut edges = vec![];
                let mut ids = BTreeMap::new();
                for loop_ in &face.loops {
                    let start = points.len();
                    for use_ in loop_ {
                        let mut directed = samples[use_.edge].clone();
                        if !use_.forward {
                            directed.reverse();
                        }
                        directed.pop();
                        for id in directed {
                            let p = graph_point(&mesh, id);
                            let uv = [u, v].map(|axis| {
                                (0..3).map(|i| (p[i] - origin_mm[i]) * axis[i]).sum::<f64>()
                            });
                            insert_cap_point(&mut ids, uv, id)?;
                            points.push(Point2::new(uv[0], uv[1]));
                        }
                    }
                    let end = points.len();
                    for i in start..end {
                        edges.push([i, if i + 1 == end { start } else { i + 1 }]);
                    }
                }
                let result = triangulate_region(points, &edges, face.loops.len() - 1)?;
                let mapping = result
                    .cdt
                    .vertices()
                    .map(|p| {
                        let p = p.position();
                        ids.get(&key([p.x, p.y]))
                            .copied()
                            .ok_or(MeshError::UnresolvedBoundary)
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                for triangle in result.triangles {
                    push(&mut mesh, triangle.map(|i| mapping[i]), ordinal, options)?;
                }
            }
            Surface::Extrusion { inward, .. } => {
                // Constructor-owned ruled walls always contain the same source
                // edge at both levels. Sharing its sampling options proves equal
                // station counts without a second geometry approximation.
                let lower = face.loops[0][if *inward { 3 } else { 0 }];
                let upper = face.loops[0][if *inward { 1 } else { 2 }];
                let low = &samples[lower.edge];
                let high = &samples[upper.edge];
                for i in 0..low.len() - 1 {
                    let mut a = [low[i], low[i + 1], high[i + 1]];
                    let mut b = [low[i], high[i + 1], high[i]];
                    if !lower.forward {
                        a.reverse();
                        b.reverse();
                    }
                    push(&mut mesh, a, ordinal, options)?;
                    push(&mut mesh, b, ordinal, options)?;
                }
            }
        }
    }
    // Every graph seam and cap trim uses the same exact sample indices.
    Ok(mesh)
}
fn graph_point(mesh: &DisplayMesh, id: usize) -> Point {
    [
        mesh.vertices[3 * id],
        mesh.vertices[3 * id + 1],
        mesh.vertices[3 * id + 2],
    ]
}

fn push(
    mesh: &mut DisplayMesh,
    triangle: [usize; 3],
    face: usize,
    options: Tessellation,
) -> Result<(), MeshError> {
    if mesh.face_indices.len() == options.triangle_limit {
        return Err(MeshError::BudgetExceeded);
    }
    mesh.triangle(triangle, face as u32);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn cap_point_ownership_rejects_duplicate_vertices_and_signed_zero_aliases() {
        let mut ids = BTreeMap::new();
        insert_cap_point(&mut ids, [0., 1.], 3).unwrap();
        assert_eq!(
            insert_cap_point(&mut ids, [-0., 1.], 4),
            Err(MeshError::UnresolvedBoundary)
        );
        let mut ids = BTreeMap::new();
        insert_cap_point(&mut ids, [1., 2.], 3).unwrap();
        assert_eq!(
            insert_cap_point(&mut ids, [1., 2.], 3),
            Err(MeshError::UnresolvedBoundary)
        );
        insert_cap_point(&mut ids, [2., 2.], 4).unwrap();
    }
}
