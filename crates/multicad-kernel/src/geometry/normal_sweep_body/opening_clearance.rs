//! A cylinder inside a certified through opening meets no sheet material.
//!
//! A certified opening removes its whole outline prism from the first side
//! of a wall to the second. A cylinder parallel to the opening's through
//! direction, whose disk lies inside the outline and whose axial interval
//! lies between the two sides, is therefore entirely void: its intersection
//! with the sweep is empty and cutting it removes nothing. A cylinder that
//! also reaches beyond a side is void within the wall, and each part beyond
//! it must miss the unpierced profile by the rectangle test below. The interval may
//! overshoot a side only by coordinate rounding, far below any feature of
//! the sheet, and the disk keeps a rounding margin from the outline.
//!
//! A cylinder whose axis lies in the profile plane projects along the sweep
//! direction onto one rectangle of that plane. When no profile edge of the
//! unpierced sweep enters that rectangle's interior and its center lies
//! outside the profile, the rectangle's interior is entirely outside the
//! profile, so the cylinder shares no interior with the sweep or with any
//! pierced version of it. The rectangle is shrunk by the rounding margin the
//! enclosing cylinders were widened by, so contact within rounding remains
//! contact.
use super::*;

fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|i| a[i] - b[i])
}

/// Whether the disk of `radius` about `p` lies in the outline with `margin`.
fn disk_inside(outline: &[Edge], p: [f64; 2], radius: f64, margin: f64) -> bool {
    if let [
        Edge::Circle {
            center, radius_mm, ..
        },
    ] = outline
    {
        return (p[0] - center[0]).hypot(p[1] - center[1]) + radius + margin <= *radius_mm;
    }
    // Otherwise the certificate holds four straight edges; a convex
    // counter-clockwise polygon contains the disk when the center keeps the
    // disk's radius from every edge line on its inner side.
    if outline.len() != 4 {
        return false;
    }
    let mut lines = Vec::with_capacity(4);
    for edge in outline {
        let Edge::Line { a, b } = *edge else {
            return false;
        };
        lines.push((a, b));
    }
    let n = lines.len();
    let convex = (0..n).all(|i| {
        let (a, b) = lines[i];
        let (_, c) = lines[(i + 1) % n];
        lines[(i + 1) % n].0 == b
            && (b[0] - a[0]) * (c[1] - b[1]) - (b[1] - a[1]) * (c[0] - b[0]) > 0.
    });
    convex
        && lines.iter().all(|&(a, b)| {
            let d = [b[0] - a[0], b[1] - a[1]];
            let length = d[0].hypot(d[1]);
            length > 0. && (d[0] * (p[1] - a[1]) - d[1] * (p[0] - a[0])) / length >= radius + margin
        })
}

/// A profile point in the frame of a rectangle: distance along its axis
/// from its base, and signed distance from its center line.
#[derive(Clone, Copy)]
struct Rectangle {
    base: [f64; 2],
    axis: [f64; 2],
    length: f64,
    half_width: f64,
}

impl Rectangle {
    fn local(&self, p: [f64; 2]) -> [f64; 2] {
        let d = [p[0] - self.base[0], p[1] - self.base[1]];
        [
            d[0] * self.axis[0] + d[1] * self.axis[1],
            d[1] * self.axis[0] - d[0] * self.axis[1],
        ]
    }

    fn center(&self) -> [f64; 2] {
        std::array::from_fn(|i| self.base[i] + self.length / 2. * self.axis[i])
    }

    /// Whether the closed segment `a`-`b` meets the open rectangle.
    fn meets_segment(&self, a: [f64; 2], b: [f64; 2]) -> bool {
        let [p, q] = [self.local(a), self.local(b)];
        let lo = [0., -self.half_width];
        let hi = [self.length, self.half_width];
        let (mut t0, mut t1) = (0_f64, 1_f64);
        for i in 0..2 {
            let d = q[i] - p[i];
            if d == 0. {
                if p[i] <= lo[i] || p[i] >= hi[i] {
                    return false;
                }
                continue;
            }
            let [s, e] = [(lo[i] - p[i]) / d, (hi[i] - p[i]) / d];
            t0 = t0.max(s.min(e));
            t1 = t1.min(s.max(e));
        }
        // A positive clipped length, or a point strictly inside.
        t0 < t1 || {
            let m: [f64; 2] = std::array::from_fn(|i| p[i] + t0 * (q[i] - p[i]));
            (0..2).all(|i: usize| lo[i] < m[i] && m[i] < hi[i])
        }
    }

    /// Whether a circle of `radius` about `center` (or any arc of it) could
    /// meet the open rectangle: false only when the whole circle is clear.
    fn may_meet_circle(&self, center: [f64; 2], radius: f64) -> bool {
        let c = self.local(center);
        let nearest = [
            c[0].clamp(0., self.length),
            c[1].clamp(-self.half_width, self.half_width),
        ];
        let near = (c[0] - nearest[0]).hypot(c[1] - nearest[1]);
        let far = [0., self.length]
            .into_iter()
            .flat_map(|s| [-self.half_width, self.half_width].map(|t| (c[0] - s).hypot(c[1] - t)))
            .fold(0., f64::max);
        !(near >= radius || far <= radius)
    }

    /// Whether an arc could meet the open rectangle. Clear only when its
    /// circle is clear, or when neither end lies in the closed rectangle and
    /// no point where the circle crosses the closed boundary lies on the arc:
    /// an arc that reaches the interior must start there or cross the
    /// boundary on its way.
    fn may_meet_arc(&self, center: [f64; 2], a: [f64; 2], b: [f64; 2], clockwise: bool) -> bool {
        let radius = (a[0] - center[0]).hypot(a[1] - center[1]);
        if !self.may_meet_circle(center, radius) {
            return false;
        }
        let [c, p, q] = [center, a, b].map(|point| self.local(point));
        let closed = |x: [f64; 2]| {
            (0. ..=self.length).contains(&x[0])
                && (-self.half_width..=self.half_width).contains(&x[1])
        };
        if closed(p) || closed(q) {
            return true;
        }
        // Counter-clockwise angle from the start, in [0, τ).
        let turn = |x: [f64; 2]| {
            let start = (p[1] - c[1]).atan2(p[0] - c[0]);
            let angle = (x[1] - c[1]).atan2(x[0] - c[0]) - start;
            let angle = if clockwise { -angle } else { angle };
            angle.rem_euclid(std::f64::consts::TAU)
        };
        let span = turn(q);
        let mut crossings = Vec::new();
        // Each boundary line, fixed coordinate `axis` at `level`, with the
        // other coordinate within its range.
        for (axis, level, range) in [
            (0, 0., [-self.half_width, self.half_width]),
            (0, self.length, [-self.half_width, self.half_width]),
            (1, -self.half_width, [0., self.length]),
            (1, self.half_width, [0., self.length]),
        ] {
            let other = 1 - axis;
            let offset = level - c[axis];
            let reach = radius * radius - offset * offset;
            if reach < 0. {
                continue;
            }
            for sign in [-1., 1.] {
                let value = c[other] + sign * reach.sqrt();
                if value >= range[0] && value <= range[1] {
                    let mut point = [0.; 2];
                    point[axis] = level;
                    point[other] = value;
                    crossings.push(point);
                }
            }
        }
        // Rounding near an end is treated as reaching the rectangle.
        let slack = 1e-9;
        crossings
            .into_iter()
            .any(|x| turn(x) <= span + slack || turn(x) >= std::f64::consts::TAU - slack)
    }

    fn meets(&self, edge: &Edge) -> bool {
        match *edge {
            Edge::Line { a, b } => self.meets_segment(a, b),
            Edge::Arc {
                center,
                a,
                b,
                clockwise,
            } => self.may_meet_arc(center, a, b, clockwise),
            Edge::Circle {
                center, radius_mm, ..
            } => self.may_meet_circle(center, radius_mm),
        }
    }
}

/// The eight corners of the box `min`-`max`, mapped by `place`.
fn box_corners(
    min: [f64; 3],
    max: [f64; 3],
    place: impl Fn([f64; 3]) -> [f64; 3],
) -> Vec<[f64; 3]> {
    (0..8)
        .map(|n: usize| {
            place(std::array::from_fn(|i| {
                if n >> i & 1 == 0 { min[i] } else { max[i] }
            }))
        })
        .collect()
}

/// Whether `p` lies in the outline with `margin` to spare.
fn point_inside(outline: &[Edge], p: [f64; 2], margin: f64) -> bool {
    disk_inside(outline, p, 0., margin)
}

/// Whether the cylinder of `radius` from `base` along unit `axis` for
/// `height` misses the unpierced profile of the sweep: its axis lies in the
/// profile plane, no profile edge enters its projected rectangle, shrunk by
/// a rounding margin, and the rectangle's center lies outside the profile.
fn cylinder_misses_profile(
    body: &NormalSweepBody,
    frame: PlaneFrame,
    base: [f64; 3],
    axis: [f64; 3],
    height: f64,
    radius: f64,
) -> bool {
    if dot(axis, frame.normal()).abs() > 64. * f64::EPSILON {
        return false;
    }
    let in_plane = |w: [f64; 3]| [dot(w, frame.u), dot(w, frame.v)];
    let direction = in_plane(axis);
    let norm = direction[0].hypot(direction[1]);
    let scale = base
        .into_iter()
        .chain(frame.origin)
        .chain([height, radius])
        .map(f64::abs)
        .fold(1., f64::max);
    let margin = 128. * f64::EPSILON * scale;
    if (norm - 1.).abs() > 64. * f64::EPSILON || height <= 2. * margin || radius <= margin {
        return false;
    }
    let start = in_plane(sub3(base, frame.origin));
    let rectangle = Rectangle {
        base: std::array::from_fn(|i| start[i] + margin * direction[i] / norm),
        axis: direction.map(|value| value / norm),
        length: height - 2. * margin,
        half_width: radius - margin,
    };
    let sweep = &body.sweep;
    if sweep
        .rings()
        .iter()
        .flatten()
        .any(|edge| rectangle.meets(edge))
    {
        return false;
    }
    let [x, y] = rectangle.center();
    !sweep.contains([x, y, sweep.height_mm() / 2.], 0.)
}

impl Solid {
    /// The unpierced profile of a sweep shares no interior with a cylinder
    /// lying across the profile plane.
    pub(in crate::geometry) fn profile_clearance(
        &self,
        other: &Self,
        operation: BooleanOperation,
        policy: GeometryPolicy,
    ) -> Option<Result<Option<Self>, GeometryError>> {
        if operation == BooleanOperation::Union
            || policy.boolean_fuzzy_mm != 0.
            || self.policy != policy
            || other.policy != policy
        {
            return None;
        }
        let (sheet, tool) = match (self.normal_sweep_body(), other.axial_cylinder_frame()) {
            (Some(_), Some(_)) => (self, other),
            _ if operation == BooleanOperation::Intersect => (other, self),
            _ => return None,
        };
        let (body, frame) = sheet.normal_sweep_body()?;
        let (cylinder, tool_frame) = tool.axial_cylinder_frame()?;
        if !cylinder_misses_profile(
            body,
            frame,
            tool_frame.origin,
            tool_frame.normal(),
            cylinder.height().as_f64(),
            cylinder.circle().radius().as_f64(),
        ) {
            return None;
        }
        Some(match operation {
            BooleanOperation::Intersect => Ok(None),
            _ => self.share().map(Some),
        })
    }

    /// Each exact rectilinear cell of the other operand is clear of the
    /// unpierced profile or lies inside one certified opening, so the
    /// operands share no interior. Both tests are convex: a profile
    /// rectangle encloses the projection of the cell's corners, and a convex
    /// opening prism holding all eight corners holds the whole cell.
    pub(in crate::geometry) fn cell_clearance(
        &self,
        other: &Self,
        operation: BooleanOperation,
        policy: GeometryPolicy,
    ) -> Option<Result<Option<Self>, GeometryError>> {
        if operation == BooleanOperation::Union
            || policy.boolean_fuzzy_mm != 0.
            || self.policy != policy
            || other.policy != policy
        {
            return None;
        }
        let (sheet, probe) = match self.normal_sweep_body() {
            Some(_) => (self, other),
            None if operation == BooleanOperation::Intersect => (other, self),
            None => return None,
        };
        let (body, frame) = sheet.normal_sweep_body()?;
        // Corners of boxes that together enclose the probe. Cuts only remove
        // material, so the stock cells of orthogonal stock enclose it in any
        // placement; other rectilinear bodies give their exact world cells.
        let boxes: Vec<Vec<[f64; 3]>> = if let Some((stock, placement)) = probe.orthogonal_body() {
            stock
                .geometry
                .cells()
                .iter()
                .map(|cell| box_corners(cell.min, cell.max, |p| position(placement, p)))
                .collect()
        } else {
            probe
                .orthogonal_cells()
                .ok()??
                .into_iter()
                .map(|cell| box_corners(cell.min, cell.max, |p| p))
                .collect()
        };
        if boxes.is_empty() {
            return None;
        }
        let sign = if body.reverse { -1. } else { 1. };
        let world = |[x, y, z]: [f64; 3]| position(frame, [x, y, sign * z]);
        let vector = |d: [f64; 3]| sub3(world(d), world([0.; 3]));
        let in_plane = |w: [f64; 3]| [dot(w, frame.u), dot(w, frame.v)];
        for corners in boxes {
            let scale = corners
                .iter()
                .flatten()
                .chain(&frame.origin)
                .map(|v| v.abs())
                .fold(1., f64::max);
            let margin = 128. * f64::EPSILON * scale;
            let projected: Vec<[f64; 2]> = corners
                .iter()
                .map(|&p| in_plane(sub3(p, frame.origin)))
                .collect();
            let [lo, hi] = [f64::min, f64::max].map(|pick| {
                [0, 1].map(|i| projected.iter().map(|p| p[i]).fold(projected[0][i], pick))
            });
            let rectangle = Rectangle {
                base: [lo[0] + margin, (lo[1] + hi[1]) / 2.],
                axis: [1., 0.],
                length: hi[0] - lo[0] - 2. * margin,
                half_width: (hi[1] - lo[1]) / 2. - margin,
            };
            let sweep = &body.sweep;
            // Beyond either cap along the sweep, touching allowed.
            let heights = corners
                .iter()
                .map(|&p| sign * dot(sub3(p, frame.origin), frame.normal()));
            let [low, high] = heights.fold([f64::INFINITY, f64::NEG_INFINITY], |[lo, hi], z| {
                [lo.min(z), hi.max(z)]
            });
            let beyond_caps = high <= margin || low >= sweep.height_mm() - margin;
            let profile_clear = beyond_caps
                || rectangle.length > 0.
                    && rectangle.half_width > 0.
                    && !sweep
                        .rings()
                        .iter()
                        .flatten()
                        .any(|edge| rectangle.meets(edge))
                    && {
                        let [x, y] = rectangle.center();
                        !sweep.contains([x, y, sweep.height_mm() / 2.], 0.)
                    };
            let in_opening = || {
                body.openings.as_ref().is_some_and(|openings| {
                    openings.geometry.openings().iter().any(|opening| {
                        let spec = opening.spec();
                        let origin = world(spec.origin);
                        let [u, v, through] = [spec.u, spec.v, spec.through].map(vector);
                        let depth2 = dot(through, through);
                        depth2 > 0.
                            && corners.iter().all(|&p| {
                                let d = sub3(p, origin);
                                let along = dot(d, through) / depth2;
                                let slack = margin / depth2.sqrt();
                                along >= -slack
                                    && along <= 1. + slack
                                    && point_inside(&spec.outline, [dot(d, u), dot(d, v)], margin)
                            })
                    })
                })
            };
            if !profile_clear && !in_opening() {
                return None;
            }
        }
        Some(match operation {
            BooleanOperation::Intersect => Ok(None),
            _ => self.share().map(Some),
        })
    }

    pub(in crate::geometry) fn opening_clearance(
        &self,
        other: &Self,
        operation: BooleanOperation,
        policy: GeometryPolicy,
    ) -> Option<Result<Option<Self>, GeometryError>> {
        if operation == BooleanOperation::Union
            || policy.boolean_fuzzy_mm != 0.
            || self.policy != policy
            || other.policy != policy
        {
            return None;
        }
        let (sheet, tool) = match (self.normal_sweep_body(), other.axial_cylinder_frame()) {
            (Some(_), Some(_)) => (self, other),
            _ if operation == BooleanOperation::Intersect => (other, self),
            _ => return None,
        };
        let (body, frame) = sheet.normal_sweep_body()?;
        let openings = body.openings.as_ref()?;
        let (cylinder, tool_frame) = tool.axial_cylinder_frame()?;
        let sign = if body.reverse { -1. } else { 1. };
        let world = |[x, y, z]: [f64; 3]| position(frame, [x, y, sign * z]);
        let vector = |d: [f64; 3]| sub3(world(d), world([0.; 3]));
        let axis = tool_frame.normal();
        let height = cylinder.height().as_f64();
        let radius = cylinder.circle().radius().as_f64();
        let base = tool_frame.origin;
        let top: [f64; 3] = std::array::from_fn(|i| base[i] + height * axis[i]);
        let scale = base
            .into_iter()
            .chain(top)
            .chain(frame.origin)
            .map(f64::abs)
            .fold(1., f64::max);
        // Wide enough to absorb the rounding margin of enclosing cylinders,
        // whose own coordinates may be larger than this pair's.
        let margin = 256. * f64::EPSILON * scale;
        let clear = openings.geometry.openings().iter().any(|opening| {
            let spec = opening.spec();
            let origin = world(spec.origin);
            let [u, v, through] = [spec.u, spec.v, spec.through].map(vector);
            let depth = dot(through, through).sqrt();
            let unit = |w: [f64; 3]| (dot(w, w).sqrt() - 1.).abs() <= 64. * f64::EPSILON;
            let square = |a: [f64; 3], b: [f64; 3]| {
                dot(a, b).abs()
                    <= 64. * f64::EPSILON * dot(a, a).sqrt().max(1.) * dot(b, b).sqrt().max(1.)
            };
            let c = cross(axis, through);
            if !(depth > 0.
                && unit(u)
                && unit(v)
                && square(u, v)
                && square(u, through)
                && square(v, through)
                && dot(c, c).sqrt() <= 64. * f64::EPSILON * depth)
            {
                return false;
            }
            let along = |p: [f64; 3]| dot(sub3(p, origin), through) / (depth * depth);
            let slack = margin / depth;
            let [low, high] = [along(base), along(top)];
            let local = sub3(base, origin);
            if !disk_inside(
                &spec.outline,
                [dot(local, u), dot(local, v)],
                radius,
                margin,
            ) {
                return false;
            }
            // Within the wall the opening is void; each part of the cylinder
            // beyond either side must miss the unpierced profile.
            let unit = through.map(|value| value / depth);
            let point =
                |t: f64| -> [f64; 3] { std::array::from_fn(|i| base[i] + (t - low) * through[i]) };
            let [lo, hi] = [low.min(high), low.max(high)];
            let mut parts = Vec::new();
            if lo < -slack {
                parts.push([lo, hi.min(0.)]);
            }
            if hi > 1. + slack {
                parts.push([lo.max(1.), hi]);
            }
            parts.into_iter().all(|[a, b]| {
                b > a
                    && cylinder_misses_profile(body, frame, point(a), unit, (b - a) * depth, radius)
            })
        });
        clear.then(|| match operation {
            BooleanOperation::Intersect => Ok(None),
            _ => self.share().map(Some),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rectangles_detect_entering_segments_and_circles_but_not_contact() {
        let rectangle = Rectangle {
            base: [0., 0.],
            axis: [1., 0.],
            length: 4.,
            half_width: 1.,
        };
        // Along an edge, outside, crossing and strictly inside.
        assert!(!rectangle.meets_segment([0., 1.], [4., 1.]));
        assert!(!rectangle.meets_segment([5., -3.], [5., 3.]));
        assert!(rectangle.meets_segment([2., -3.], [2., 3.]));
        assert!(rectangle.meets_segment([1., 0.], [1., 0.5]));
        // A circle around the whole rectangle, a far one and a crossing one.
        assert!(!rectangle.may_meet_circle([2., 0.], 3.));
        assert!(!rectangle.may_meet_circle([2., 5.], 3.));
        assert!(rectangle.may_meet_circle([2., 3.], 2.5));
        assert_eq!(rectangle.center(), [2., 0.]);
        // A quarter arc whose circle crosses the rectangle but which turns
        // away from it, and the opposite quarter that enters it.
        assert!(!rectangle.may_meet_arc([-0.5, 0.], [-0.5, -1.], [-1.5, 0.], true));
        assert!(rectangle.may_meet_arc([-0.5, 0.], [-0.5, -1.], [0.5, 0.], false));
    }

    #[test]
    fn a_disk_must_keep_its_radius_inside_a_circle_or_convex_quadrilateral() {
        let circle = [Edge::Circle {
            center: [1., 1.],
            radius_mm: 3.,
            clockwise: false,
        }];
        assert!(disk_inside(&circle, [1., 1.], 2.9, 0.));
        assert!(!disk_inside(&circle, [1.5, 1.], 2.9, 0.));
        let corners = [[0., 0.], [4., 0.], [4., 2.], [0., 2.]];
        let square: Vec<_> = (0..4)
            .map(|i| Edge::Line {
                a: corners[i],
                b: corners[(i + 1) % 4],
            })
            .collect();
        assert!(disk_inside(&square, [2., 1.], 1., 0.));
        assert!(!disk_inside(&square, [2., 1.], 1., 1e-12));
        assert!(!disk_inside(&square, [3.5, 1.], 1., 0.));
        let clockwise: Vec<_> = square
            .iter()
            .rev()
            .map(|edge| match *edge {
                Edge::Line { a, b } => Edge::Line { a: b, b: a },
                _ => unreachable!(),
            })
            .collect();
        assert!(!disk_inside(&clockwise, [2., 1.], 0.5, 0.));
    }
}
