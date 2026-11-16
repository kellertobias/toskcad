//! Vertices re-derived from the supports that meet there.
//!
//! A file may state its vertices only to its declared uncertainty. Each
//! vertex is then moved to the point where its edge curves and face
//! surfaces actually meet, so that the retained boundary agrees with itself
//! exactly instead of within the wider tolerance.
use super::{Boundary, Edge, Error, Point, Result, add, distance, ensure, scale};
use std::collections::BTreeSet;

/// Projection sweeps before a vertex must have settled.
const SWEEPS: usize = 64;

/// The supports of every vertex: incident edges and faces.
fn supports(boundary: &Boundary) -> Vec<(BTreeSet<usize>, BTreeSet<usize>)> {
    let mut incident = vec![(BTreeSet::new(), BTreeSet::new()); boundary.vertices.len()];
    for (index, edge) in boundary.edges.iter().enumerate() {
        for vertex in edge.vertices {
            incident[vertex].0.insert(index);
        }
    }
    for (index, face) in boundary.faces.iter().enumerate() {
        for trim in &face.loops {
            if let Some(vertex) = trim.vertex {
                incident[vertex].1.insert(index);
            }
            for &(edge, _) in &trim.uses {
                for vertex in boundary.edges[edge].vertices {
                    incident[vertex].1.insert(index);
                }
            }
        }
    }
    incident
}

/// Move every vertex onto the common point of its supports, found by
/// averaging the closest points on them until it settles. Returns the
/// boundary and the largest remaining distance from a vertex to one of its
/// supports. A vertex whose supports do not settle within `tolerance` of
/// each of them and of the stated vertex keeps its stated point if that
/// lies within `tolerance` of every support, and is refused as invalid
/// otherwise.
pub fn rederive(boundary: &Boundary, tolerance: f64) -> Result<(Boundary, f64)> {
    let incident = supports(boundary);
    let mut vertices = boundary.vertices.clone();
    let mut residual: f64 = 0.;
    for (vertex, (edges, faces)) in incident.iter().enumerate() {
        let stated = vertices[vertex];
        let mut point = stated;
        let mut worst = f64::INFINITY;
        let mut at_stated = f64::INFINITY;
        for sweep in 0..SWEEPS {
            let mut projections: Vec<Point> = Vec::with_capacity(edges.len() + faces.len());
            for &edge in edges {
                let curve = &boundary.edges[edge].curve;
                projections.push(curve.evaluate(curve.invert(point, None)).0);
            }
            for &face in faces {
                let surface = &boundary.faces[face].surface;
                projections.push(surface.evaluate(surface.invert(point, None)).0);
            }
            worst = projections
                .iter()
                .map(|&p| distance(p, point))
                .fold(0., f64::max);
            if sweep == 0 {
                at_stated = worst;
            }
            let count = projections.len() as f64;
            let next = scale(projections.into_iter().fold([0.; 3], add), 1. / count);
            let settled = distance(next, point) <= 1e-15 * (1. + super::norm(point));
            // Every vertex the reader keeps bounds an edge or a vertex loop.
            point = next;
            if settled {
                break;
            }
        }
        if !(worst <= tolerance && distance(point, stated) <= tolerance) {
            // Supports that disagree within the tolerance have no common
            // point, and their balance may lie farther from one of them
            // than the stated vertex does, which then stands.
            ensure(at_stated <= tolerance, Error::Invalid)?;
            (point, worst) = (stated, at_stated);
        }
        residual = residual.max(worst);
        vertices[vertex] = point;
    }
    let edges = boundary
        .edges
        .iter()
        .map(|edge| {
            let points = edge.vertices.map(|vertex| vertices[vertex]);
            Edge::new(edge.curve.clone(), edge.vertices, points, tolerance)
        })
        .collect::<Result<Vec<_>>>()?;
    Ok((
        Boundary {
            vertices,
            edges,
            faces: boundary.faces.clone(),
            shells: boundary.shells.clone(),
        },
        residual,
    ))
}
