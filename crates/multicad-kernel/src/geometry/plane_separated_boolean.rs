//! Booleans of straight-edged planar operands split by one face plane.
//!
//! Each operand lies in the convex hull of its own vertices when all of its
//! faces are planar and all of its edges are straight. A face plane of either
//! operand with every vertex of one operand on one side and every vertex of
//! the other on the opposite side, up to rounding, therefore proves that the
//! operands share no interior: their intersection is empty and a difference
//! leaves the first operand. A union across such a contact would glue the
//! operands and is left to the other providers. No surface is intersected.
use super::*;
use multicad_core::references::AnalyticGeometry;

type Point = [f64; 3];

/// Vertices and face planes of a straight-edged planar body.
fn planar_hull(solid: &Solid) -> Option<(Vec<Point>, Vec<(Point, Point)>)> {
    let topology = solid.rust_topology()?.ok()?;
    let mut vertices = Vec::new();
    for element in &topology.elements {
        match element.geometry {
            Some(AnalyticGeometry::Point { position_mm }) => vertices.push(position_mm),
            Some(AnalyticGeometry::Line { .. } | AnalyticGeometry::Plane { .. }) | None => {}
            Some(_) => return None,
        }
    }
    let planes = solid
        .faces()
        .ok()?
        .into_iter()
        .map(|face| face.plane)
        .collect::<Option<Vec<_>>>()?;
    (!vertices.is_empty() && !planes.is_empty()).then_some((vertices, planes))
}

fn rounding(points: &[Point]) -> f64 {
    64. * f64::EPSILON * points.iter().flatten().map(|x| x.abs()).fold(1., f64::max)
}

/// Whether `plane` puts `a` on one closed side and `b` on the other.
fn splits((origin, normal): (Point, Point), a: &[Point], b: &[Point], margin: f64) -> bool {
    let side = |p: &Point| dot(std::array::from_fn(|i| p[i] - origin[i]), normal);
    let range = |points: &[Point]| {
        points
            .iter()
            .map(side)
            .fold([f64::INFINITY, f64::NEG_INFINITY], |[lo, hi], s| {
                [lo.min(s), hi.max(s)]
            })
    };
    let ([a_lo, a_hi], [b_lo, b_hi]) = (range(a), range(b));
    (a_hi <= margin && b_lo >= -margin) || (b_hi <= margin && a_lo >= -margin)
}

impl Solid {
    pub(super) fn plane_separated_boolean(
        &self,
        other: &Self,
        operation: BooleanOperation,
        policy: GeometryPolicy,
    ) -> Option<Result<Option<Self>, GeometryError>> {
        if operation == BooleanOperation::Union
            || policy.boolean_fuzzy_mm != 0.
            || policy != self.policy
            || policy != other.policy
        {
            return None;
        }
        let (a, a_planes) = planar_hull(self)?;
        let (b, b_planes) = planar_hull(other)?;
        let margin = rounding(&a).max(rounding(&b));
        if !a_planes
            .into_iter()
            .chain(b_planes)
            .any(|plane| splits(plane, &a, &b, margin))
        {
            return None;
        }
        Some(match operation {
            BooleanOperation::Intersect => Ok(None),
            _ => self.share().map(Some),
        })
    }
}

#[cfg(test)]
mod tests;
