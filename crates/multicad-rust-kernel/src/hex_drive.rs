//! Literal finite flat-hex drive contours. A size label never replaces corners.
use crate::{headed_screw, helical_thread::ThreadError, tapered_thread::interval::I};
type Point = [f64; 2];
fn cross(a: Point, b: Point) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}
fn difference(a: Point, b: Point) -> Point {
    [a[0] - b[0], a[1] - b[1]]
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HexDrive {
    corners: [Point; 6],
    across_flats: f64,
    canonical: bool,
}
impl HexDrive {
    pub fn canonical(across_flats: f64, source_error_mm: f64) -> Result<Self, ThreadError> {
        Self::from_retained_corners(
            headed_screw::hex_vertices(across_flats),
            across_flats,
            source_error_mm,
        )
    }
    /// Retain the caller's ordered physical corners. Size, centering and
    /// regularity are qualified within the separately exposed source budget;
    /// corners are never rotated, fitted or replaced with ideal vertices.
    pub fn from_retained_corners(
        corners: [Point; 6],
        across_flats: f64,
        source_error_mm: f64,
    ) -> Result<Self, ThreadError> {
        if !across_flats.is_finite()
            || across_flats <= 0.
            || !source_error_mm.is_finite()
            || source_error_mm <= 0.
            || !corners.into_iter().flatten().all(f64::is_finite)
        {
            return Err(ThreadError::Invalid);
        }
        if across_flats > 1e7 || corners.into_iter().flatten().any(|x| x.abs() > 1e7) {
            return Err(ThreadError::Limit);
        }
        let radius = across_flats / 3_f64.sqrt();
        let budget = 8. * source_error_mm;
        let minimum = 64. * source_error_mm * source_error_mm;
        let mut center = [0.; 2];
        for i in 0..6 {
            let a = corners[i];
            let b = corners[(i + 1) % 6];
            let c = corners[(i + 2) % 6];
            let edge = difference(b, a);
            let next = difference(c, b);
            if (a[0].hypot(a[1]) - radius).abs() > budget
                || (edge[0].hypot(edge[1]) - radius).abs() > budget
                || cross(a, b) <= minimum
                || cross(edge, next) <= minimum
            {
                return Err(ThreadError::Invalid);
            }
            for axis in 0..2 {
                center[axis] += a[axis] / 6.;
            }
        }
        if center[0].hypot(center[1]) > budget {
            return Err(ThreadError::Invalid);
        }
        Ok(Self {
            corners,
            across_flats,
            canonical: corners == headed_screw::hex_vertices(across_flats),
        })
    }
    pub fn corners_mm(&self) -> [Point; 6] {
        self.corners
    }
    pub fn is_canonical(&self) -> bool {
        self.canonical
    }
    pub fn edge_mm(&self, side: usize) -> [Point; 2] {
        let side = side % 6;
        [self.corners[side], self.corners[(side + 1) % 6]]
    }
    pub fn side_length_mm(&self, side: usize) -> f64 {
        if self.canonical {
            self.across_flats / 3_f64.sqrt()
        } else {
            let [a, b] = self.edge_mm(side);
            let d = difference(b, a);
            d[0].hypot(d[1])
        }
    }
    pub fn area_mm2(&self) -> f64 {
        if self.canonical {
            headed_screw::hex_area(self.across_flats)
        } else {
            (0..6)
                .map(|i| cross(self.corners[i], self.corners[(i + 1) % 6]) * 0.5)
                .sum()
        }
    }
    pub fn center_mm(&self) -> Point {
        if self.canonical {
            return [0.; 2];
        }
        let mut moment = [0.; 2];
        for i in 0..6 {
            let a = self.corners[i];
            let b = self.corners[(i + 1) % 6];
            let weight = cross(a, b);
            for axis in 0..2 {
                moment[axis] += (a[axis] + b[axis]) * weight;
            }
        }
        moment.map(|x| x / (6. * self.area_mm2()))
    }
    pub fn maximum_radius_mm(&self) -> f64 {
        self.corners
            .into_iter()
            .map(|[x, y]| (I::point(x).square() + I::point(y).square()).sqrt().hi)
            .fold(0., f64::max)
    }
    pub fn area_bounds_mm2(&self, error_mm: f64) -> Option<[f64; 2]> {
        if !error_mm.is_finite() || error_mm < 0. {
            return None;
        }
        let point = |p: Point| p.map(|x| I::point(x) + I::new(-error_mm, error_mm));
        let mut area = I::point(0.);
        for i in 0..6 {
            let [ax, ay] = point(self.corners[i]);
            let [bx, by] = point(self.corners[(i + 1) % 6]);
            area = area + (ax * by - ay * bx) * I::point(0.5);
        }
        let exact = self.area_mm2();
        area.valid()
            .then_some([area.lo.min(exact).max(0.), area.hi.max(exact)])
    }
    pub fn side_at_angle(&self, angle: f64) -> Option<usize> {
        if !angle.is_finite() {
            return None;
        }
        let d = [angle.cos(), angle.sin()];
        let mut best = (f64::INFINITY, 0);
        for i in 0..6 {
            let [a, b] = self.edge_mm(i);
            let denominator = cross(difference(b, a), d);
            if denominator < 0. {
                let r = -cross(a, b) / denominator;
                if r < best.0 {
                    best = (r, i);
                }
            }
        }
        Some(best.1)
    }
    pub fn radius_at_angle_mm(&self, angle: f64) -> f64 {
        if !angle.is_finite() {
            return f64::NAN;
        }
        if self.canonical {
            return headed_screw::hex_radius(self.across_flats, angle);
        }
        let [a, b] = self.edge_mm(
            self.side_at_angle(angle)
                .expect("validated drive contains the axis"),
        );
        -cross(a, b) / cross(difference(b, a), [angle.cos(), angle.sin()])
    }
}

#[cfg(test)]
mod tests;
