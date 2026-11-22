use super::*;
use std::f64::consts::{PI, TAU};

use crate::spline_sweep::{ConstructionError, validation::*};
fn certificate_error(error: ConstructionError) -> ReadError {
    match error {
        ConstructionError::Invalid => ReadError::Invalid,
        ConstructionError::Limit => ReadError::Limit,
        ConstructionError::Unresolved | ConstructionError::Cancelled => ReadError::Unsupported,
    }
}

pub(super) fn verify(
    reader: &Shell<'_, '_>,
    edges: &Edges,
    plane: &Plane,
    face: &FaceWire,
) -> Result<f64> {
    if let ([rim], [hole]) = (face.outer.as_slice(), face.holes.as_slice())
        && let Geometry::Analytic(Curve::Circle(circle)) = edges[&rim.edge].geometry
    {
        return verify_ring(reader, edges, plane, (circle, rim.forward), hole);
    }
    let circular_hole = face.holes.first().and_then(|ring| match ring.as_slice() {
        [edge] => match edges[&edge.edge].geometry {
            Geometry::Analytic(Curve::Circle(circle)) => Some((circle, edge)),
            _ => None,
        },
        _ => None,
    });
    let legacy_hole = face.holes.is_empty() || (face.holes.len() == 1 && circular_hole.is_some());
    if !legacy_hole {
        return verify_general(reader, edges, plane, face);
    }
    if face.outer.len() < 2
        || face
            .outer
            .iter()
            .map(|edge| edge.start)
            .collect::<BTreeSet<_>>()
            .len()
            != face.outer.len()
        || !face
            .outer
            .iter()
            .any(|edge| matches!(edges[&edge.edge].geometry, Geometry::Cubic { .. }))
    {
        return Err(ReadError::Unsupported);
    }
    let tolerance = reader.tolerance;
    let section = face
        .outer
        .iter()
        .map(|edge| boundary(reader, &edges[&edge.edge].geometry, plane, edge))
        .collect::<Result<Vec<_>>>()?;
    let mut budget = Budget {
        remaining: MAX_PAIR_WORK,
        interrupt: None,
    };
    if section.len() == 2 {
        convex_pair(&section, tolerance)?;
    } else {
        simple(&section, tolerance, &mut budget).map_err(certificate_error)?;
    }
    let mut area = face
        .outer
        .iter()
        .map(|edge| section_area(reader, &edges[&edge.edge].geometry, plane, edge))
        .sum::<Result<f64>>()?;
    if let Some((circle, edge)) = circular_hole {
        if !circle.within_plane(plane, tolerance)
            || (dot(circle.normal, plane.normal) > 0.) == edge.forward
        {
            return Err(ReadError::Invalid);
        }
        let center = project(plane, circle.center_mm);
        let winding: f64 = section
            .iter()
            .flat_map(|edge| &edge.pieces)
            .map(|&piece| clear_disk(piece, center, circle.radius_mm + tolerance, 0, &mut budget))
            .sum::<std::result::Result<f64, ConstructionError>>()
            .map_err(certificate_error)?;
        if (winding - TAU).abs() > 1e-7 {
            return Err(ReadError::Invalid);
        }
        area -= PI * circle.radius_mm.powi(2);
    }
    if !area.is_finite() || area <= tolerance.powi(2) {
        return Err(ReadError::Invalid);
    }
    Ok(area)
}

/// Two edges meet at both of their vertices, so no direction of common
/// travel separates them. Their closed control polygon instead turns
/// strictly one way and exactly once: every line then meets the polygon, and
/// by variation diminution the curve, at most twice, so the ring is convex.
fn convex_pair(section: &[Boundary], tolerance: f64) -> Result<()> {
    let mut polygon = Vec::new();
    for piece in section.iter().flat_map(|edge| &edge.pieces) {
        match *piece {
            Piece::Line([start, _]) => polygon.push(start),
            Piece::Bezier([start, first, second, _]) => polygon.extend([start, first, second]),
            Piece::Rational(span) => {
                let points = span.points();
                polygon.extend(&points[..points.len() - 1]);
            }
            Piece::Arc(_) => return Err(ReadError::Unsupported),
        }
    }
    let corner = |index: usize| {
        let [a, b, c] = [0, 1, 2].map(|step| polygon[(index + step) % polygon.len()]);
        [sub2(b, a), sub2(c, b)]
    };
    let side = cross2(corner(0)[0], corner(0)[1]).signum();
    let mut turning = 0.;
    for index in 0..polygon.len() {
        let [before, after] = corner(index);
        let turn = side * cross2(before, after);
        if turn.is_nan() || turn <= tolerance * norm(before).max(norm(after)) {
            return Err(ReadError::Unsupported);
        }
        turning += turn.portable_atan2(dot2(before, after));
    }
    if (turning - TAU).abs() > 1e-9 {
        return Err(ReadError::Unsupported);
    }
    Ok(())
}

/// An internal ring cap: one counterclockwise outside circle around a
/// simple clockwise cavity with at least one cubic flank. Every cavity piece
/// stays strictly inside the circle by its analytic radial reach.
fn verify_ring(
    reader: &Shell<'_, '_>,
    edges: &Edges,
    plane: &Plane,
    (circle, forward): (CircleEdge, bool),
    hole: &[Traversal],
) -> Result<f64> {
    let tolerance = reader.tolerance;
    let section = hole
        .iter()
        .map(|edge| boundary(reader, &edges[&edge.edge].geometry, plane, edge))
        .collect::<Result<Vec<_>>>()?;
    let mut budget = Budget {
        remaining: MAX_PAIR_WORK,
        interrupt: None,
    };
    simple(&section, tolerance, &mut budget).map_err(certificate_error)?;
    let center = project(plane, circle.center_mm);
    let distance = |p: Point2| norm(sub2(p, center));
    let reach = |piece: &Piece| match *piece {
        Piece::Line(points) => points.map(distance).into_iter().fold(0., f64::max),
        Piece::Arc(arc) => distance(arc.center) + arc.radius,
        Piece::Bezier(points) => points.map(distance).into_iter().fold(0., f64::max),
        Piece::Rational(span) => span.points().into_iter().map(distance).fold(0., f64::max),
    };
    let cavity = hole
        .iter()
        .map(|edge| section_area(reader, &edges[&edge.edge].geometry, plane, edge))
        .sum::<Result<f64>>()?;
    let area = PI * circle.radius_mm.powi(2) + cavity;
    if !circle.within_plane(plane, tolerance)
        || (dot(circle.normal, plane.normal) > 0.) != forward
        || hole
            .iter()
            .map(|edge| edge.start)
            .collect::<BTreeSet<_>>()
            .len()
            != hole.len()
        || !hole
            .iter()
            .any(|edge| matches!(edges[&edge.edge].geometry, Geometry::Cubic { .. }))
        || section
            .iter()
            .flat_map(|edge| &edge.pieces)
            .any(|piece| reach(piece) >= circle.radius_mm - tolerance)
        || cavity >= 0.
        || !area.is_finite()
        || area <= tolerance.powi(2)
    {
        return Err(ReadError::Invalid);
    }
    Ok(area)
}

fn section_area(
    reader: &Shell<'_, '_>,
    geometry: &Geometry,
    plane: &Plane,
    edge: &Traversal,
) -> Result<f64> {
    let [a, b] = [edge.start, edge.end].map(|v| project(plane, reader.points[&v]));
    match geometry {
        Geometry::Analytic(Curve::Line) => Ok(cross2(a, b) / 2.),
        Geometry::Analytic(curve) => {
            let round = round_geometry(*curve).ok_or(ReadError::Invalid)?;
            let center = project(plane, round.center);
            let theta = directed_sweep(round, plane, edge.forward);
            Ok((cross2(center, sub2(b, a)) + round.radius.powi(2) * theta) / 2.)
        }
        Geometry::Cubic { curve, .. } => {
            let signed = curve.signed_area(plane.origin, plane.normal);
            Ok(if edge.forward { signed } else { -signed })
        }
    }
}

fn boundary(
    reader: &Shell<'_, '_>,
    geometry: &Geometry,
    plane: &Plane,
    edge: &Traversal,
) -> Result<Boundary> {
    let points = [edge.start, edge.end].map(|id| reader.points[&id]);
    let tolerance = reader.tolerance;
    let on_plane = |p: Point| dot(sub(p, plane.origin), plane.normal).abs() <= tolerance;
    if !points.into_iter().all(on_plane) {
        return Err(ReadError::Invalid);
    }
    let [a, b] = points.map(|p| project(plane, p));
    let pieces = match geometry {
        Geometry::Analytic(Curve::Line) => vec![Piece::Line([a, b])],
        Geometry::Analytic(curve) => {
            let round = round_geometry(*curve).ok_or(ReadError::Invalid)?;
            if round.full
                || !on_plane(round.center)
                || length(cross(round.normal, plane.normal)) * round.radius > tolerance
            {
                return Err(ReadError::Invalid);
            }
            let center = project(plane, round.center);
            let [x, y] = sub2(a, center);
            let arc = Arc {
                center,
                radius: round.radius,
                start: y.portable_atan2(x),
                sweep: directed_sweep(round, plane, edge.forward),
            };
            vec![Piece::Arc(arc)]
        }
        Geometry::Cubic { curve, .. } if curve.weights().is_some() => {
            if !curve.poles_mm().iter().copied().all(on_plane) {
                return Err(ReadError::Invalid);
            }
            let poles = curve
                .poles_mm()
                .iter()
                .map(|&p| project(plane, p))
                .collect();
            let weights = curve.weights().unwrap_or_default().to_vec();
            let rational = crate::spline_sweep::Rational::new(poles, weights, !edge.forward)
                .map_err(|_| ReadError::Invalid)?;
            let piece = Piece::Rational(rational.span());
            let chord = unit2(sub2(b, a)).ok_or(ReadError::Invalid)?;
            if !piece.forward_along(chord) {
                return Err(ReadError::Unsupported);
            }
            vec![piece]
        }
        Geometry::Cubic { curve, .. } => {
            if !curve.poles_mm().iter().copied().all(on_plane) {
                return Err(ReadError::Invalid);
            }
            let mut spans: Vec<_> = curve
                .bezier_spans()
                .map(|span| span.map(|p| project(plane, p)))
                .collect();
            if !edge.forward {
                spans.reverse();
                for span in &mut spans {
                    span.reverse();
                }
            }
            let chord = unit2(sub2(b, a)).ok_or(ReadError::Invalid)?;
            if spans
                .iter()
                .any(|span| !Piece::Bezier(*span).forward_along(chord))
            {
                return Err(ReadError::Unsupported);
            }
            spans.into_iter().map(Piece::Bezier).collect()
        }
    };
    let bounds = pieces
        .iter()
        .map(|p| p.bounds())
        .reduce(merge)
        .ok_or(ReadError::Invalid)?;
    Ok(Boundary { pieces, bounds })
}

fn directed_sweep(round: Round, plane: &Plane, forward: bool) -> f64 {
    if (dot(round.normal, plane.normal) > 0.) == forward {
        round.sweep
    } else {
        -round.sweep
    }
}

fn project(plane: &Plane, p: Point) -> Point2 {
    let d = sub(p, plane.origin);
    [dot(d, plane.u), dot(d, plane.v)]
}

fn verify_general(
    reader: &Shell<'_, '_>,
    edges: &Edges,
    plane: &Plane,
    face: &FaceWire,
) -> Result<f64> {
    if face.outer.len() < 3 {
        return Err(ReadError::Unsupported);
    }

    for ring in std::iter::once(&face.outer).chain(&face.holes) {
        if ring.is_empty()
            || ring
                .iter()
                .map(|edge| edge.start)
                .collect::<BTreeSet<_>>()
                .len()
                != ring.len()
        {
            return Err(ReadError::Invalid);
        }
        for edge in ring {
            if edge.start == edge.end
                && matches!(edges[&edge.edge].geometry, Geometry::Analytic(Curve::Line))
            {
                return Err(ReadError::Invalid);
            }
            match edges[&edge.edge].geometry {
                Geometry::Analytic(Curve::Circle(circle)) => {
                    if ring.len() != 1 || !circle.within_plane(plane, reader.tolerance) {
                        return Err(ReadError::Invalid);
                    }
                }
                _ => {
                    boundary(reader, &edges[&edge.edge].geometry, plane, edge)?;
                }
            }
        }
    }
    let outer = reconstruction::section_edges(reader, edges, plane, &face.outer, false)?;
    let holes = face
        .holes
        .iter()
        .map(|ring| reconstruction::section_edges(reader, edges, plane, ring, false))
        .collect::<Result<Vec<_>>>()?;
    let sweep = crate::spline_sweep::Sweep::new_rings(outer, holes, 1., reader.tolerance)
        .map_err(certificate_error)?;
    Ok(sweep.area_mm2())
}

#[cfg(test)]
mod error_tests {
    use super::*;
    #[test]
    fn section_certificate_keeps_typed_decoder_failure_contract() {
        for (source, expected) in [
            (ConstructionError::Invalid, ReadError::Invalid),
            (ConstructionError::Limit, ReadError::Limit),
            (ConstructionError::Unresolved, ReadError::Unsupported),
            (ConstructionError::Cancelled, ReadError::Unsupported),
        ] {
            assert_eq!(certificate_error(source), expected);
        }
    }
}

#[cfg(test)]
mod convex_pair_tests {
    use super::*;

    fn edge(pieces: Vec<Piece>) -> Boundary {
        let bounds = [[0.; 2]; 2];
        Boundary { pieces, bounds }
    }

    fn ring(points: &[Point2]) -> Vec<Boundary> {
        let next = |index: usize| points[(index + 1) % points.len()];
        let lines = (0..points.len()).map(|index| Piece::Line([points[index], next(index)]));
        vec![edge(lines.collect()), edge(Vec::new())]
    }

    #[test]
    fn two_edges_need_a_convex_control_polygon_that_turns_once() {
        let cubic = [[20., 0.], [20., 12.], [0., 12.], [0., 0.]];
        let chord = Piece::Line([[0., 0.], [20., 0.]]);
        let section = vec![edge(vec![Piece::Bezier(cubic)]), edge(vec![chord])];
        assert_eq!(convex_pair(&section, 1e-7), Ok(()));
        // The reversed ring is convex too; its winding is judged by its area.
        let mut reversed = cubic;
        reversed.reverse();
        let back = Piece::Line([[20., 0.], [0., 0.]]);
        let section = vec![edge(vec![back]), edge(vec![Piece::Bezier(reversed)])];
        assert_eq!(convex_pair(&section, 1e-7), Ok(()));
        let arc = Piece::Arc(Arc {
            center: [10., 0.],
            radius: 10.,
            start: 0.,
            sweep: PI,
        });
        let section = vec![edge(vec![arc]), edge(vec![back])];
        assert_eq!(convex_pair(&section, 1e-7), Err(ReadError::Unsupported));
        // A reflex corner, a straight corner and a ring that turns twice.
        let reflex = [[0., 0.], [4., 0.], [4., 4.], [2., 1.], [0., 4.]];
        let straight = [[0., 0.], [2., 0.], [4., 0.], [4., 4.]];
        let star = [0., 2., 4., 1., 3.].map(|step: f64| {
            let (sin, cos) = (step * TAU / 5.).portable_sin_cos();
            [cos, sin]
        });
        for points in [&reflex[..], &straight, &star] {
            assert_eq!(
                convex_pair(&ring(points), 1e-7),
                Err(ReadError::Unsupported)
            );
        }
    }
}
