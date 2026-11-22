//! Recognize quarter-cylinder fillet sweeps and preserve their exact volume in
//! the common height coordinate used by the Boolean.
use super::*;
use multicad_rust_kernel::normal_sweep::Edge;
use std::f64::consts::FRAC_PI_2;

#[derive(Clone, Copy)]
pub(super) struct QuarterFillet {
    pub(super) center: Point,
    pub(super) radius: f64,
    axis: usize,
    span: [f64; 2],
    direction: f64,
    back: f64,
    bottom: f64,
}

impl QuarterFillet {
    pub(super) fn recognize(solid: &Solid, tolerance: f64) -> Option<Self> {
        Self::read(solid, tolerance, false)
    }

    /// Also accept a bare fillet whose base or back wall has zero thickness:
    /// the arc ends then meet the remaining straight sides directly.
    pub(super) fn recognize_open(solid: &Solid, tolerance: f64) -> Option<Self> {
        Self::read(solid, tolerance, true)
    }

    fn read(solid: &Solid, tolerance: f64, open: bool) -> Option<Self> {
        let (body, frame) = solid.normal_sweep_body()?;
        let (rings, height) = body.unpierced_profile()?;
        let [ring] = rings else { return None };
        if ring.len() != 5 && !(open && (3..5).contains(&ring.len())) {
            return None;
        }
        let normal = frame.normal();
        let axis = (0..2).find(|&i| normal[i].abs() == 1.)?;
        if normal[(axis + 1) % 2] != 0. || normal[2] != 0. {
            return None;
        }
        let transverse = 1 - axis;
        let mut arcs = ring.iter().filter(|edge| matches!(edge, Edge::Arc { .. }));
        let &arc = arcs.next()?;
        if arcs.next().is_some() {
            return None;
        }
        let Edge::Arc { center, a, b, .. } = arc else {
            return None;
        };
        let center = frame.point(center);
        let radius = arc.radius_mm()?;
        let [mut low, mut high] = [frame.point(a), frame.point(b)];
        if low[2] > high[2] {
            std::mem::swap(&mut low, &mut high);
        }
        let direction = (high[transverse] - center[transverse]).signum();
        let mut expected_low = center;
        expected_low[2] -= radius;
        let mut expected_high = center;
        expected_high[transverse] += radius * direction;
        let near = |a: Point, b: Point| norm(subtract(a, b)) <= tolerance;
        if radius <= tolerance
            || !near(low, expected_low)
            || !near(high, expected_high)
            || (arc.signed_sweep_rad()?.abs() - FRAC_PI_2).abs() * radius > tolerance
        {
            return None;
        }
        let lines: Vec<_> = ring
            .iter()
            .filter_map(|edge| match *edge {
                Edge::Line { a, b } => Some([frame.point(a), frame.point(b)]),
                _ => None,
            })
            .collect();
        if lines.len() + 1 != ring.len() {
            return None;
        }
        let mut bottom = lines.iter().flatten().map(|p| p[2]).reduce(f64::min)?;
        let mut back = if direction > 0. {
            lines
                .iter()
                .flatten()
                .map(|p| p[transverse])
                .reduce(f64::max)?
        } else {
            lines
                .iter()
                .flatten()
                .map(|p| p[transverse])
                .reduce(f64::min)?
        };
        let base = expected_low[2] - bottom;
        let wall = (back - expected_high[transverse]) * direction;
        if !open && (base <= tolerance || wall <= tolerance)
            || base < -tolerance
            || wall < -tolerance
        {
            return None;
        }
        // A side within the tolerance of the arc end has no thickness.
        if base <= tolerance {
            bottom = expected_low[2];
        }
        if wall <= tolerance {
            back = expected_high[transverse];
        }
        let mut at_bottom = center;
        at_bottom[2] = bottom;
        let mut outer_bottom = at_bottom;
        outer_bottom[transverse] = back;
        let mut outer_top = center;
        outer_top[transverse] = back;
        let expected = [
            [expected_low, at_bottom],
            [at_bottom, outer_bottom],
            [outer_bottom, outer_top],
            [outer_top, expected_high],
        ];
        let expected: Vec<_> = expected.into_iter().filter(|&[a, b]| a != b).collect();
        if expected.len() != lines.len()
            || !expected.iter().all(|&[a, b]| {
                lines
                    .iter()
                    .any(|&[c, d]| near(a, c) && near(b, d) || near(a, d) && near(b, c))
            })
        {
            return None;
        }
        let a = center[axis];
        let b = a + normal[axis] * height;
        Some(Self {
            center,
            radius,
            axis,
            span: [a.min(b), a.max(b)],
            direction,
            back,
            bottom,
        })
    }

    pub(super) fn pieces(
        self,
        center_z: f64,
        radius: f64,
        tolerance: f64,
    ) -> Result<Vec<ConvexPiece>, planar_boolean::Error> {
        let axis = self.axis;
        let transverse = 1 - axis;
        let center = self.center[transverse];
        let floor = center_z - radius;
        let mut min = [0.; 3];
        let mut max = [0.; 3];
        [min[axis], max[axis]] = self.span;
        min[transverse] = center.min(self.back);
        max[transverse] = center.max(self.back);
        min[2] = self.bottom - floor;
        let mut pieces = Vec::with_capacity(2);
        if min[2] < 0. {
            pieces.push(ConvexPiece::from_box(min, max, tolerance)?);
        }
        let mut section = vec![[center, 0.], [self.back, 0.], [self.back, radius]];
        let rim = center + self.direction * radius;
        if rim != self.back {
            section.push([rim, radius]);
        }
        let unit = std::array::from_fn(|i| f64::from(i == axis));
        let at = |axis_value, [other, z]: [f64; 2]| {
            let mut point = [0., 0., z];
            point[axis] = axis_value;
            point[transverse] = other;
            point
        };
        let place = |section: &[[f64; 2]], axis_value| -> Vec<Point> {
            section.iter().map(|&p| at(axis_value, p)).collect()
        };
        if winding(&place(&section, self.span[0]), unit) < 0. {
            section.reverse();
        }
        let lower = place(&section, self.span[0]);
        let upper = place(&section, self.span[1]);
        let mut faces = vec![lower.iter().rev().copied().collect(), upper.clone()];
        for i in 0..section.len() {
            let next = (i + 1) % section.len();
            faces.push(vec![lower[i], lower[next], upper[next], upper[i]]);
        }
        pieces.push(ConvexPiece::new(faces, tolerance)?);
        Ok(pieces)
    }
}
