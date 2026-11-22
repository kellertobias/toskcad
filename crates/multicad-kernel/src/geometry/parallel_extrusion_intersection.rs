//! Intersection of a circular cylinder with a prism swept along its axis.
//!
//! When a straight prism's sweep direction is parallel to a cylinder's axis,
//! every section normal to that axis inside both height intervals is the
//! same planar disk-polygon intersection. The common solid is that region
//! swept over the overlap of the two intervals:
//!
//! * a disk clear of every profile boundary is either inside the material,
//!   giving the cylinder trimmed to the overlap, or outside (or inside a
//!   round hole), giving nothing;
//! * a disk crossing one convex outline gives a line-and-arc profile.
//!
//! Tangent or near-tangent contacts, ends that touch within tolerance and
//! disks that cross non-convex or holed outlines are refused, so a later
//! provider or the typed Boolean refusal decides them.
use super::*;
use crate::{CurvedEdge2d, CurvedProfile};

type Point2 = [f64; 2];

/// The prism section in the cylinder chart and its axial interval.
struct Section {
    rings: Vec<Vec<Point2>>,
    holes: Vec<(Point2, f64)>,
    interval: [f64; 2],
}

/// The relative sense of two unit axes, or `None` when they are not
/// parallel to rounding.
fn parallel(a: [f64; 3], b: [f64; 3]) -> Option<f64> {
    let c = cross(a, b);
    (dot(c, c).sqrt() <= 64. * f64::EPSILON).then(|| dot(a, b).signum())
}

fn distance(a: Point2, b: Point2) -> f64 {
    (a[0] - b[0]).portable_hypot(a[1] - b[1])
}

/// Distance from the chart origin to the closed segment `a`-`b`.
fn segment_distance(a: Point2, b: Point2) -> f64 {
    let d = [b[0] - a[0], b[1] - a[1]];
    let length = d[0] * d[0] + d[1] * d[1];
    let t = (-(a[0] * d[0] + a[1] * d[1]) / length).clamp(0., 1.);
    distance([a[0] + t * d[0], a[1] + t * d[1]], [0.; 2])
}

fn signed_area(ring: &[Point2]) -> f64 {
    ring.iter()
        .zip(ring.iter().cycle().skip(1))
        .map(|(a, b)| a[0] * b[1] - a[1] * b[0])
        .sum::<f64>()
        / 2.
}

fn contains(ring: &[Point2], p: Point2) -> bool {
    let mut inside = false;
    for (a, b) in ring.iter().zip(ring.iter().cycle().skip(1)) {
        if (a[1] > p[1]) != (b[1] > p[1])
            && p[0] < a[0] + (p[1] - a[1]) * (b[0] - a[0]) / (b[1] - a[1])
        {
            inside = !inside;
        }
    }
    inside
}

fn strictly_convex(ring: &[Point2]) -> bool {
    let n = ring.len();
    (0..n).all(|i| {
        let [a, b, c] = [ring[i], ring[(i + 1) % n], ring[(i + 2) % n]];
        (b[0] - a[0]) * (c[1] - b[1]) - (b[1] - a[1]) * (c[0] - b[0]) > 0.
    })
}

/// The boundary of a counter-clockwise convex polygon intersected with the
/// disk of `radius` about the origin, or `None` at a tangency within
/// `tolerance`. Arcs run counter-clockwise from an exit to the next entry.
fn clip(ring: &[Point2], radius: f64, tolerance: f64) -> Option<Vec<CurvedEdge2d>> {
    enum Mark {
        Corner(Point2),
        Enter(Point2),
        Exit(Point2),
    }
    let r2 = radius * radius;
    if ring
        .iter()
        .any(|&p| (distance(p, [0.; 2]) - radius).abs() <= tolerance)
    {
        return None;
    }
    let mut marks = Vec::new();
    let mut inside = distance(ring[0], [0.; 2]) < radius;
    for (&a, &b) in ring.iter().zip(ring.iter().cycle().skip(1)) {
        if distance(a, [0.; 2]) < radius {
            marks.push(Mark::Corner(a));
        }
        let d = [b[0] - a[0], b[1] - a[1]];
        let qa = d[0] * d[0] + d[1] * d[1];
        let qb = 2. * (a[0] * d[0] + a[1] * d[1]);
        let qc = a[0] * a[0] + a[1] * a[1] - r2;
        let foot = (-qb / (2. * qa)).clamp(0., 1.);
        let closest = distance([a[0] + foot * d[0], a[1] + foot * d[1]], [0.; 2]);
        if (closest - radius).abs() <= tolerance {
            return None;
        }
        let discriminant = qb * qb - 4. * qa * qc;
        if discriminant <= 0. {
            continue;
        }
        let root = discriminant.sqrt();
        for s in [(-qb - root) / (2. * qa), (-qb + root) / (2. * qa)] {
            if s <= 0. || s >= 1. {
                continue;
            }
            let p = [a[0] + s * d[0], a[1] + s * d[1]];
            marks.push(if inside {
                Mark::Exit(p)
            } else {
                Mark::Enter(p)
            });
            inside = !inside;
        }
    }
    let point = |mark: &Mark| match *mark {
        Mark::Corner(p) | Mark::Enter(p) | Mark::Exit(p) => p,
    };
    let edges = (0..marks.len())
        .map(|i| {
            let (from, to) = (&marks[i], &marks[(i + 1) % marks.len()]);
            let (a, b) = (point(from), point(to));
            match (from, to) {
                (Mark::Exit(_), Mark::Enter(_)) => Some(CurvedEdge2d::Arc {
                    center: [0.; 2],
                    a,
                    b,
                    clockwise: false,
                }),
                (Mark::Exit(_), _) | (_, Mark::Enter(_)) => None,
                _ => Some(CurvedEdge2d::Line { a, b }),
            }
        })
        .collect::<Option<Vec<_>>>()?;
    (edges.len() >= 2).then_some(edges)
}

/// The axial interval, in the cylinder frame, where world-aligned cells hold
/// the whole disk; `Some(None)` when no cell meets the cylinder. Cells that
/// only partly cover the disk, several separate intervals and ends within
/// tolerance of each other are not resolved.
fn cell_interval(
    cells: &[multicad_rust_kernel::orthogonal::Box3],
    frame: PlaneFrame,
    radius: f64,
    height: f64,
    tolerance: f64,
) -> Option<Option<[f64; 2]>> {
    let normal = frame.normal();
    let k = (0..3).find(|&k| normal[k].abs() == 1.)?;
    if (0..3).any(|i| i != k && normal[i] != 0.) {
        return None;
    }
    let sense = normal[k];
    let axial = |world: f64| (world - frame.origin[k]) * sense;
    let transverse = [(k + 1) % 3, (k + 2) % 3];
    let mut spans = Vec::new();
    for cell in cells {
        let apart = transverse.iter().any(|&i| {
            cell.max[i] + tolerance < frame.origin[i] - radius
                || cell.min[i] - tolerance > frame.origin[i] + radius
        });
        let [a, b] = [axial(cell.min[k]), axial(cell.max[k])];
        let [low, high] = [a.min(b), a.max(b)];
        if apart || high < -tolerance || low > height + tolerance {
            continue;
        }
        let covers = transverse.iter().all(|&i| {
            cell.min[i] < frame.origin[i] - radius - tolerance
                && cell.max[i] > frame.origin[i] + radius + tolerance
        });
        let near_end = [low, high].into_iter().any(|e| {
            [0., height]
                .into_iter()
                .any(|end| e != end && (e - end).abs() <= tolerance)
        });
        if !covers || near_end {
            return None;
        }
        spans.push([low.max(0.), high.min(height)]);
    }
    spans.sort_by(|a, b| a[0].total_cmp(&b[0]));
    let mut merged: Option<[f64; 2]> = None;
    for [low, high] in spans {
        merged = match merged {
            None => Some([low, high]),
            Some([a, b]) if low <= b => Some([a, b.max(high)]),
            Some(_) => return None,
        };
    }
    match merged {
        Some([low, high]) if high - low <= tolerance => None,
        merged => Some(merged),
    }
}

/// The part of `span` outside every void span, when it is empty or one
/// interval; ends within tolerance of each other are not resolved.
fn remove_voids(span: [f64; 2], voids: &[[f64; 2]], tolerance: f64) -> Option<Option<[f64; 2]>> {
    let [mut low, mut high] = span;
    let mut voids = voids.to_vec();
    voids.sort_by(|a, b| a[0].total_cmp(&b[0]));
    for [a, b] in voids {
        if b <= low || a >= high {
            continue;
        }
        if [a, b].into_iter().any(|e| {
            [low, high]
                .into_iter()
                .any(|end| e != end && (e - end).abs() <= tolerance)
        }) {
            return None;
        }
        match (a <= low, b >= high) {
            (true, true) => return Some(None),
            (true, false) => low = b,
            (false, true) => high = a,
            (false, false) => return None,
        }
    }
    Some(Some([low, high]))
}

impl Solid {
    /// The section of a straight prism swept parallel to `axis`, in the
    /// chart of `frame` (whose normal is `axis`).
    fn parallel_section(&self, frame: PlaneFrame) -> Option<Section> {
        let axis = frame.normal();
        let chart = |p: [f64; 3]| {
            let d = std::array::from_fn(|i| p[i] - frame.origin[i]);
            [dot(d, frame.u), dot(d, frame.v)]
        };
        let along = |p: [f64; 3]| {
            let d: [f64; 3] = std::array::from_fn(|i| p[i] - frame.origin[i]);
            dot(d, axis)
        };
        let sorted = |a: f64, b: f64| [a.min(b), a.max(b)];
        if let Some((prism, placement)) = self.rust_prism_frame() {
            if prism.rounded_end().is_some() || prism.inner_loops().len() != 0 {
                return None;
            }
            let size = prism.maximum().map(Rational::as_f64);
            let axes = [placement.u, placement.v, placement.normal()];
            let k = (0..3).find(|&k| parallel(axes[k], axis).is_some())?;
            let [i, j] = [(k + 1) % 3, (k + 2) % 3];
            let corner = |a: f64, b: f64| {
                chart(std::array::from_fn(|n| {
                    placement.origin[n] + a * axes[i][n] + b * axes[j][n]
                }))
            };
            let mut ring = vec![
                corner(0., 0.),
                corner(size[i], 0.),
                corner(size[i], size[j]),
                corner(0., size[j]),
            ];
            if signed_area(&ring) < 0. {
                ring.reverse();
            }
            let start = along(placement.origin);
            let end: [f64; 3] = std::array::from_fn(|n| placement.origin[n] + size[k] * axes[k][n]);
            return Some(Section {
                rings: vec![ring],
                holes: vec![],
                interval: sorted(start, along(end)),
            });
        }
        let placed = self.polygon_prism()?;
        if !placed.pockets.is_empty() || placed.bore_treatments.is_some() {
            return None;
        }
        parallel(placed.frame.normal(), axis)?;
        let world = |[x, y]: [Rational; 2]| {
            let [x, y] = [x.as_f64(), y.as_f64()];
            std::array::from_fn(|n| {
                placed.frame.origin[n] + x * placed.frame.u[n] + y * placed.frame.v[n]
            })
        };
        let rings = placed
            .prism
            .rings()
            .iter()
            .map(|ring| {
                let mut ring: Vec<_> = ring.iter().map(|&p| chart(world(p))).collect();
                // The chart may mirror the profile; keep the outer ring
                // counter-clockwise.
                if signed_area(&ring) < 0. {
                    ring.reverse();
                }
                ring
            })
            .collect();
        let holes = placed
            .prism
            .circular_holes()
            .iter()
            .map(|circle| (chart(world(circle.center())), circle.radius().as_f64()))
            .collect();
        let height = placed.prism.height().as_f64();
        let sign = if placed.reverse_sweep { -1. } else { 1. };
        let end: [f64; 3] = std::array::from_fn(|n| {
            placed.frame.origin[n] + sign * height * placed.frame.normal()[n]
        });
        Some(Section {
            rings,
            holes,
            interval: sorted(along(placed.frame.origin), along(end)),
        })
    }

    /// This cylinder between the axial heights `low` and `high` of its frame.
    fn trimmed_cylinder(&self, low: f64, high: f64) -> Result<Self, GeometryError> {
        let (cylinder, frame) = self.axial_cylinder_frame().expect("checked cylinder");
        if low == 0. && high == cylinder.height().as_f64() {
            return self.share();
        }
        let normal = frame.normal();
        let frame = PlaneFrame {
            origin: std::array::from_fn(|i| frame.origin[i] + low * normal[i]),
            ..frame
        };
        let [x, y, z] = construction::cylinder_coordinates(frame);
        let circle = Circle::new([x, y], cylinder.circle().radius())
            .map_err(|cause| input(format!("Invalid trimmed cylinder: {cause:?}")))?;
        let trimmed = AxialCylinder::new(
            cylinder.feature(),
            circle,
            z,
            construction::exact(high - low)?,
        )
        .map_err(|cause| {
            error(
                "Boolean",
                GeometryErrorKind::Limit,
                format!("Trimmed cylinder exceeds exact arithmetic admission: {cause:?}"),
            )
        })?;
        match self.cylinder_sketch() {
            Some(sketch) => Self::from_axial_cylinder_sweep(
                trimmed,
                frame,
                sketch.center_uv,
                sketch.reverse,
                sketch.origin.clone(),
                self.policy,
            ),
            None => Self::build_axial_cylinder(trimmed, frame, None, self.policy),
        }
    }

    pub(super) fn parallel_extrusion_intersection(
        &self,
        other: &Self,
        operation: BooleanOperation,
        policy: GeometryPolicy,
    ) -> Option<Result<Option<Self>, GeometryError>> {
        if operation != BooleanOperation::Intersect
            || policy.boolean_fuzzy_mm != 0.
            || policy != self.policy
            || policy != other.policy
        {
            return None;
        }
        let (cylinder, prism) = if self.axial_cylinder_frame().is_some() {
            (self, other)
        } else {
            (other, self)
        };
        let (geometry, frame) = cylinder.axial_cylinder_frame()?;
        let tolerance = policy.linear_tolerance_mm;
        let radius = geometry.circle().radius().as_f64();
        let height = geometry.height().as_f64();
        let Some(section) = prism.parallel_section(frame) else {
            let (cells, voids) = prism
                .orthogonal_cells_near(frame, radius, height, tolerance)
                .ok()??;
            let normal = frame.normal();
            let axial = |p: [f64; 3]| {
                let d: [f64; 3] = std::array::from_fn(|i| p[i] - frame.origin[i]);
                dot(d, normal)
            };
            let voids: Vec<_> = voids
                .into_iter()
                .map(|[a, b]| {
                    let [a, b] = [axial(a), axial(b)];
                    [a.min(b), a.max(b)]
                })
                .collect();
            let material = match cell_interval(&cells, frame, radius, height, tolerance)? {
                Some(span) => remove_voids(span, &voids, tolerance)?,
                None => None,
            };
            return Some(match material {
                Some([low, high]) => cylinder.trimmed_cylinder(low, high).map(Some),
                None => Ok(None),
            });
        };
        let [low, high] = [section.interval[0].max(0.), section.interval[1].min(height)];
        if (high - low).abs() <= tolerance
            || [0., height].into_iter().any(|end| {
                section
                    .interval
                    .iter()
                    .any(|&e| e != end && (e - end).abs() <= tolerance)
            })
        {
            return None;
        }
        if high < low {
            return Some(Ok(None));
        }
        // A disk inside a round hole meets no material.
        if section
            .holes
            .iter()
            .any(|&(center, r)| distance(center, [0.; 2]) + radius < r - tolerance)
        {
            return Some(Ok(None));
        }
        let clear = section.rings.iter().all(|ring| {
            ring.iter()
                .zip(ring.iter().cycle().skip(1))
                .all(|(&a, &b)| segment_distance(a, b) > radius + tolerance)
        }) && section
            .holes
            .iter()
            .all(|&(center, r)| distance(center, [0.; 2]) > radius + r + tolerance);
        if clear {
            let inside = contains(&section.rings[0], [0.; 2])
                && section.rings[1..]
                    .iter()
                    .all(|ring| !contains(ring, [0.; 2]));
            return Some(if inside {
                cylinder.trimmed_cylinder(low, high).map(Some)
            } else {
                Ok(None)
            });
        }
        let [ring] = &section.rings[..] else {
            return None;
        };
        if !section.holes.is_empty() || !strictly_convex(ring) {
            return None;
        }
        let outer = clip(ring, radius, tolerance)?;
        let normal = frame.normal();
        let profile = CurvedProfile {
            frame: PlaneFrame {
                origin: std::array::from_fn(|i| frame.origin[i] + low * normal[i]),
                ..frame
            },
            outer,
            holes: vec![],
        };
        Some(
            Self::extrude_curved_profile_vector(&profile, normal.map(|n| n * (high - low)), policy)
                .map(Some),
        )
    }
}

#[cfg(test)]
mod tests;
