//! Source-known independent proof of every finite curve, support and wire.
use super::normal_relief::Frame;
use super::spline_sweep::{Edges, FaceWire, Geometry};
use super::*;
use crate::{
    spline_cap_relief::{Curve as ExpectedCurve, Relief, Surface},
    spline_sweep::Edge,
};
const MAX_WORK: usize = 24_000_000;
fn require(value: bool) -> Result<()> {
    value.then_some(()).ok_or(ReadError::Invalid)
}
fn charge(work: &mut usize, n: usize) -> Result<()> {
    *work = work.checked_sub(n).ok_or(ReadError::Limit)?;
    Ok(())
}
fn near(a: Point, b: Point, tol: f64) -> bool {
    length(sub(a, b)) <= tol
}
fn aligned(a: Point, b: Point, extent: f64, tol: f64) -> bool {
    dot(a, b) > 0. && length(cross(a, b)) * extent.max(1.) <= tol
}
fn metric_frame(frame: Frame, bounds: [Point; 2], tolerance: f64) -> Result<()> {
    use crate::normal_relief::precision::{add_up, dot_error, mul_up};
    let mut gram = 0.;
    let mut extent = 0.;
    for (i, a) in frame.axes.iter().enumerate() {
        for (j, b) in frame.axes.iter().enumerate() {
            let e = dot_error(*a, *b, if i == j { 1. } else { 0. });
            gram = add_up(gram, mul_up(e, e));
        }
        let e = bounds[0][i].abs().max(bounds[1][i].abs());
        extent = add_up(extent, mul_up(e, e));
    }
    let sqrt = |v: f64| if v == 0. { 0. } else { v.sqrt().next_up() };
    require(mul_up(sqrt(gram), sqrt(extent)) <= tolerance)
}
fn cyclic(mut wire: Vec<(u32, bool)>) -> Vec<(u32, bool)> {
    if let Some((i, _)) = wire.iter().enumerate().min_by_key(|(_, x)| **x) {
        wire.rotate_left(i);
    }
    wire
}
fn curve_matches(
    actual: &Geometry,
    expected: &ExpectedCurve,
    forward: bool,
    frame: Frame,
    tol: f64,
) -> bool {
    match (actual, expected) {
        (
            Geometry::Analytic(Curve::Line),
            ExpectedCurve::Line { .. }
            | ExpectedCurve::Profile {
                edge: Edge::Analytic(crate::normal_sweep::Edge::Line { .. }),
                ..
            },
        ) => true,
        (
            Geometry::Cubic { curve, .. },
            ExpectedCurve::Profile {
                edge: Edge::Cubic(source),
                z_mm,
            },
        ) => {
            if forward == source.reversed()
                || curve.knots() != source.knots()
                || curve.poles_mm().len() != source.poles().len()
            {
                return false;
            }
            curve
                .poles_mm()
                .iter()
                .zip(source.poles())
                .all(|(&actual, &p)| near(actual, frame.world([p[0], p[1], *z_mm]), tol))
        }
        (
            Geometry::Analytic(actual),
            ExpectedCurve::Profile {
                edge: Edge::Analytic(source),
                z_mm,
            },
        ) => {
            let Some((center, axis, radius, sweep)) = super::normal_relief::proof::circle(*actual)
            else {
                return false;
            };
            let (xy, clockwise) = match source {
                crate::normal_sweep::Edge::Arc {
                    center, clockwise, ..
                }
                | crate::normal_sweep::Edge::Circle {
                    center, clockwise, ..
                } => (*center, *clockwise),
                _ => return false,
            };
            let class = matches!(
                (actual, source),
                (Curve::Circle(_), crate::normal_sweep::Edge::Circle { .. })
                    | (
                        Curve::Arc(_) | Curve::Semicircle(_),
                        crate::normal_sweep::Edge::Arc { .. }
                    )
            );
            let expected_radius = source.radius_mm().expect("circular source");
            class
                && near(center, frame.world([xy[0], xy[1], *z_mm]), tol)
                && (radius - expected_radius).abs() <= tol
                && aligned(
                    axis.map(|x| if forward { x } else { -x }),
                    frame.axes[2].map(|x| if clockwise { -x } else { x }),
                    radius,
                    tol,
                )
                && (sweep - source.signed_sweep_rad().expect("circular source").abs()).abs()
                    * radius
                    <= tol
        }
        _ => false,
    }
}
pub fn verify_spline_cap_relief(
    input: &str,
    relief: &Relief,
    origin: Point,
    u: Point,
    v: Point,
    tolerance: f64,
) -> Result<f64> {
    compound::read_members_by(input, tolerance, true, true, |shell, solid, _| {
        verify_member(shell, solid, relief, [origin, u, v])
    })?
    .pop()
    .ok_or(ReadError::Invalid)
}
pub(super) fn verify_member(
    shell: &mut Shell<'_, '_>,
    solid: u32,
    relief: &Relief,
    [origin, u, v]: [Point; 3],
) -> Result<f64> {
    let frame = Frame::new(origin, u, v)?;
    metric_frame(frame, relief.bounds(), shell.tolerance / 32.)?;
    let ids = shell.member_faces(solid)?;
    let graph = relief.boundary();
    require(ids.len() == graph.faces.len())?;
    let mut edges = Edges::new();
    let faces = ids
        .iter()
        .map(|&id| spline_sweep::read_face(shell, id, &mut edges))
        .collect::<Result<Vec<_>>>()?;
    require(edges.len() == graph.edges.len() && edges.values().all(|edge| edge.uses == [1, 1]))?;
    let tol = shell.tolerance / 16.;
    let mut work = MAX_WORK;
    // Match unique physical vertices, retaining actual vertex-ID ownership.
    require(shell.points.len() == graph.vertices.len())?;
    let mut points = shell
        .points
        .iter()
        .map(|(&id, &p)| (p, id))
        .collect::<Vec<_>>();
    points.sort_unstable_by(|a, b| a.0[0].total_cmp(&b.0[0]));
    let mut used = BTreeSet::new();
    let mut vertices = Vec::new();
    for &p in &graph.vertices {
        charge(&mut work, 1)?;
        let p = frame.world(p);
        let a = points.partition_point(|(q, _)| q[0] < p[0] - tol);
        let b = points.partition_point(|(q, _)| q[0] <= p[0] + tol);
        charge(&mut work, b - a)?;
        let matches = points[a..b]
            .iter()
            .filter(|(q, _)| near(*q, p, tol))
            .collect::<Vec<_>>();
        let &[&(_, id)] = matches.as_slice() else {
            return Err(ReadError::Invalid);
        };
        require(used.insert(id))?;
        vertices.push(id);
    }
    let mut by_ends: BTreeMap<[u32; 2], Vec<u32>> = BTreeMap::new();
    for (&id, edge) in &edges {
        let mut key = edge.vertices;
        key.sort_unstable();
        by_ends.entry(key).or_default().push(id);
        charge(&mut work, 1)?;
    }
    let mut mapped = Vec::new();
    used.clear();
    for edge in &graph.edges {
        let ends = edge.vertices.map(|v| vertices[v]);
        let mut key = ends;
        key.sort_unstable();
        let candidates = by_ends.get(&key).ok_or(ReadError::Invalid)?;
        charge(&mut work, 2 * candidates.len())?;
        let mut matches = vec![];
        for &id in candidates {
            for forward in [true, false] {
                let actual = &edges[&id];
                if (if forward {
                    actual.vertices
                } else {
                    [actual.vertices[1], actual.vertices[0]]
                }) == ends
                    && curve_matches(&actual.geometry, &edge.curve, forward, frame, tol)
                {
                    matches.push((id, forward));
                }
            }
        }
        let [(id, forward)] = matches[..] else {
            return Err(ReadError::Invalid);
        };
        require(used.insert(id))?;
        mapped.push((id, forward));
    }
    let mut by_wire = BTreeMap::new();
    for (i, face) in faces.iter().enumerate() {
        let wire = cyclic(face.outer.iter().map(|e| (e.edge, e.forward)).collect());
        require(by_wire.insert(wire, i).is_none())?;
    }
    let extent = relief
        .bounds()
        .into_iter()
        .flatten()
        .map(f64::abs)
        .fold(1., f64::max);
    for expected in &graph.faces {
        let mut wires = expected
            .loops
            .iter()
            .map(|ring| {
                cyclic(
                    ring.iter()
                        .map(|e| {
                            let (id, forward) = mapped[e.edge];
                            (id, e.forward == forward)
                        })
                        .collect(),
                )
            })
            .collect::<Vec<_>>();
        charge(&mut work, wires.iter().map(Vec::len).sum())?;
        let index = by_wire.remove(&wires[0]).ok_or(ReadError::Invalid)?;
        let actual = &faces[index];
        let mut holes = actual
            .holes
            .iter()
            .map(|r| cyclic(r.iter().map(|e| (e.edge, e.forward)).collect()))
            .collect::<Vec<_>>();
        wires.remove(0);
        wires.sort();
        holes.sort();
        require(wires == holes)?;
        support(shell, actual, &expected.surface, frame, tol, extent, &edges)?;
    }
    require(by_wire.is_empty())?;
    Ok(relief.volume_mm3())
}
fn support(
    shell: &Shell<'_, '_>,
    face: &FaceWire,
    expected: &Surface,
    frame: Frame,
    tol: f64,
    extent: f64,
    edges: &Edges,
) -> Result<()> {
    match expected {
        Surface::Plane { z_mm, upward } => {
            let plane = spline_sweep::face_plane(shell, face)?;
            require(
                aligned(
                    plane.normal,
                    frame.axes[2].map(|n| if *upward { n } else { -n }),
                    extent,
                    tol,
                ) && dot(
                    sub(frame.world([0., 0., *z_mm]), plane.origin),
                    plane.normal,
                )
                .abs()
                    <= tol,
            )
        }
        Surface::Extrusion {
            edge,
            span_mm,
            inward,
        } => {
            match edge {
                Edge::Analytic(crate::normal_sweep::Edge::Line { a, b }) => {
                    let plane = spline_sweep::face_plane(shell, face)?;
                    let tangent = frame.vector([b[0] - a[0], b[1] - a[1], 0.]);
                    let normal =
                        unit(cross(tangent, frame.axes[2]).map(|n| if *inward { -n } else { n }))?;
                    require(
                        aligned(plane.normal, normal, extent, tol)
                            && dot(
                                sub(frame.world([a[0], a[1], span_mm[0]]), plane.origin),
                                plane.normal,
                            )
                            .abs()
                                <= tol,
                    )
                }
                Edge::Analytic(
                    edge @ (crate::normal_sweep::Edge::Arc {
                        center, clockwise, ..
                    }
                    | crate::normal_sweep::Edge::Circle {
                        center, clockwise, ..
                    }),
                ) => {
                    let (plane, radius) = shell.cylinder(face.surface)?;
                    let expected_radius = edge.radius_mm().expect("circular source");
                    require(
                        face.forward == (!*clockwise != *inward)
                            && (radius - expected_radius).abs() <= tol
                            && length(cross(plane.normal, frame.axes[2])) * extent <= tol
                            && length(cross(
                                sub(
                                    frame.world([center[0], center[1], span_mm[0]]),
                                    plane.origin,
                                ),
                                plane.normal,
                            )) <= tol,
                    )
                }
                Edge::Cubic(curve) => {
                    let (support, axis) = shell.spline_extrusion(face.surface)?;
                    require(
                        face.forward == (!curve.reversed() != *inward)
                            && aligned(axis, frame.axes[2], extent, tol),
                    )?;
                    let low =
                        face.outer
                            .iter()
                            .filter_map(|use_| match &edges[&use_.edge].geometry {
                                Geometry::Cubic { id, curve } => Some((*id, curve)),
                                _ => None,
                            })
                            .find(|(_, actual)| {
                                actual.poles_mm().iter().zip(curve.poles()).all(|(&a, &p)| {
                                    near(a, frame.world([p[0], p[1], span_mm[0]]), tol)
                                })
                            })
                            .ok_or(ReadError::Invalid)?;
                    require(support == low.0)
                }
            }
        }
    }
}

mod recognize;
pub use recognize::{
    ImportedSplineCapRelief, read_spline_cap_relief, read_spline_cap_relief_compound,
};

pub(super) fn recognize_member(shell: &Shell<'_, '_>, id: u32) -> Result<(Readback, ImportSolid)> {
    let mut clean = Shell::new(shell.entities, shell.tolerance, true);
    let (proof, body) = recognize::read_member(&mut clean, id)?;
    Ok((
        proof,
        ImportSolid {
            origin_mm: body.origin_mm,
            u: body.u,
            v: body.v,
            height_mm: body.relief.base().height_mm(),
            profile: ImportProfile::SplineCapRelief {
                relief: Box::new(body.relief),
            },
        },
    ))
}
mod matcher;
pub use matcher::spline_cap_relief_import_matches;

#[cfg(test)]
mod tests;
