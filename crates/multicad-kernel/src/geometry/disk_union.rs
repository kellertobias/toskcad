//! Pairwise union of parallel cylinders over one common axial span.
//!
//! Parallel cylinders that occupy the same span along their common axis
//! unite to the extrusion of the planar union of their disks. That union's
//! boundary consists of circle arcs between pairwise crossing points; every
//! pair is classified as clearly apart, clearly nested or clearly crossing,
//! and every boundary vertex clears each third circle, so tangencies and
//! triple points are refused rather than guessed. A swept profile made only
//! of such arcs is accepted as an operand when its own disks rebuild exactly
//! its boundary, so a hub takes several lobes one Boolean at a time.
use super::*;
use crate::{CurvedEdge2d, CurvedProfile};
use multicad_rust_kernel::normal_sweep::Edge;
use std::collections::BTreeMap;
use std::f64::consts::TAU;

type Point = [f64; 3];

#[derive(Clone)]
struct Disk {
    center: [f64; 2],
    radius: f64,
    source: EntityId,
}

/// The common chart of all operands: base plane frame, signed height.
struct Chart {
    frame: PlaneFrame,
    height: f64,
}

/// One boundary arc, counter-clockwise about its disk, between two
/// crossing vertices, or a whole uncovered circle.
#[derive(Clone, Copy, Debug)]
struct BoundaryArc {
    disk: usize,
    ends: Option<[usize; 2]>,
}

struct Union {
    vertices: Vec<[f64; 2]>,
    loops: Vec<Vec<BoundaryArc>>,
}

fn sub2(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] - b[0], a[1] - b[1]]
}

fn length2(a: [f64; 2]) -> f64 {
    a[0].portable_hypot(a[1])
}

fn cross2(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}

fn angle(disk: &Disk, point: [f64; 2]) -> f64 {
    let d = sub2(point, disk.center);
    d[1].portable_atan2(d[0]).rem_euclid(TAU)
}

/// The planar union boundary of `disks`, or `None` for any tangency,
/// coincidence, triple point or several separate regions.
fn union(disks: &[Disk], tolerance: f64) -> Option<Union> {
    let mut vertices = Vec::new();
    let mut crossings = BTreeMap::new();
    // Covered angular intervals per disk: (start, end vertex) along CCW.
    let mut covered: Vec<Vec<[usize; 2]>> = vec![Vec::new(); disks.len()];
    let mut inside = vec![false; disks.len()];
    for i in 0..disks.len() {
        for j in i + 1..disks.len() {
            let [a, b] = [&disks[i], &disks[j]];
            let delta = sub2(b.center, a.center);
            let d = length2(delta);
            let [sum, gap] = [a.radius + b.radius, (a.radius - b.radius).abs()];
            if d >= sum + tolerance {
                continue;
            }
            if d + tolerance <= gap {
                if a.radius < b.radius {
                    inside[i] = true;
                } else {
                    inside[j] = true;
                }
                continue;
            }
            if d <= gap + tolerance || d >= sum - tolerance {
                return None;
            }
            // Crossing points: `left` lies counter-clockwise of a -> b.
            let along = (a.radius * a.radius - b.radius * b.radius + d * d) / (2. * d);
            let half = (a.radius * a.radius - along * along).max(0.).sqrt();
            let unit = [delta[0] / d, delta[1] / d];
            let base = [a.center[0] + along * unit[0], a.center[1] + along * unit[1]];
            let left = [base[0] - half * unit[1], base[1] + half * unit[0]];
            let right = [base[0] + half * unit[1], base[1] - half * unit[0]];
            let [l, r] = [vertices.len(), vertices.len() + 1];
            vertices.extend([left, right]);
            crossings.insert((i, j), [l, r]);
            // On disk a, the arc inside b runs CCW from `right` to `left`;
            // on disk b, the arc inside a runs CCW from `left` to `right`.
            covered[i].push([r, l]);
            covered[j].push([l, r]);
        }
    }
    // Every vertex that bounds the union must clear every third circle.
    let mut boundary_vertex = vec![true; vertices.len()];
    for (&(i, j), ends) in &crossings {
        for &v in ends {
            for (k, disk) in disks.iter().enumerate() {
                if k == i || k == j {
                    continue;
                }
                let reach = length2(sub2(vertices[v], disk.center));
                if (reach - disk.radius).abs() <= tolerance {
                    return None;
                }
                if reach < disk.radius {
                    boundary_vertex[v] = false;
                }
            }
        }
    }
    let mut arcs = Vec::new();
    for (index, disk) in disks.iter().enumerate() {
        if inside[index] {
            continue;
        }
        if covered[index].is_empty() {
            arcs.push(BoundaryArc {
                disk: index,
                ends: None,
            });
            continue;
        }
        // An interval end that lies in no other disk starts an uncovered
        // arc, which runs to the next interval start along the circle.
        for &[_, from] in &covered[index] {
            if !boundary_vertex[from] {
                continue;
            }
            let at = angle(disk, vertices[from]);
            let mut next = None;
            for &[start, _] in &covered[index] {
                let delta = (angle(disk, vertices[start]) - at).rem_euclid(TAU);
                if delta <= 0. {
                    return None;
                }
                if next.is_none_or(|(best, _)| delta < best) {
                    next = Some((delta, start));
                }
            }
            let (_, to) = next?;
            if !boundary_vertex[to] {
                return None;
            }
            arcs.push(BoundaryArc {
                disk: index,
                ends: Some([from, to]),
            });
        }
    }
    // Chain the arcs at their shared vertices.
    let mut loops = Vec::new();
    let mut outgoing = BTreeMap::new();
    for (index, arc) in arcs.iter().enumerate() {
        match arc.ends {
            None => loops.push(vec![*arc]),
            Some([from, _]) => {
                if outgoing.insert(from, index).is_some() {
                    return None;
                }
            }
        }
    }
    let mut used = vec![false; arcs.len()];
    for start in 0..arcs.len() {
        if used[start] || arcs[start].ends.is_none() {
            continue;
        }
        let mut ring = Vec::new();
        let mut at = start;
        while !used[at] {
            used[at] = true;
            ring.push(arcs[at]);
            let [_, to] = arcs[at].ends?;
            at = *outgoing.get(&to)?;
        }
        if at != start {
            return None;
        }
        loops.push(ring);
    }
    Some(Union { vertices, loops })
}

/// Signed area of one loop of counter-clockwise disk arcs.
fn loop_area(union: &Union, disks: &[Disk], ring: &[BoundaryArc]) -> f64 {
    ring.iter()
        .map(|arc| {
            let disk = &disks[arc.disk];
            let r = disk.radius;
            let Some([from, to]) = arc.ends else {
                return std::f64::consts::PI * r * r;
            };
            let [a, b] = [union.vertices[from], union.vertices[to]];
            let mut sweep = angle(disk, b) - angle(disk, a);
            if sweep <= 0. {
                sweep += TAU;
            }
            // Chord triangle from the origin plus the circular segment.
            0.5 * cross2(a, b) + 0.5 * r * r * (sweep - sweep.portable_sin())
        })
        .sum()
}

impl Solid {
    /// The disks of one operand in `chart`, or `None` when it is not a
    /// parallel cylinder or arc-only sweep over the same span.
    fn chart_disks(&self, chart: &Chart, tolerance: f64) -> Option<Vec<Disk>> {
        let normal = chart.frame.normal();
        let local = |p: Point| {
            let d: Point = std::array::from_fn(|i| p[i] - chart.frame.origin[i]);
            [dot(d, chart.frame.u), dot(d, chart.frame.v), dot(d, normal)]
        };
        let same_span = |start: Point, axis: Point, height: f64| {
            let base = local(start)[2];
            let top = base + dot(axis, normal) * height;
            let [lo, hi] = [base.min(top), base.max(top)];
            let [want_lo, want_hi] = [chart.height.min(0.), chart.height.max(0.)];
            (1. - dot(axis, normal).abs()) <= 1e-12
                && (lo - want_lo).abs() <= tolerance
                && (hi - want_hi).abs() <= tolerance
        };
        if let Some((cylinder, frame)) = self.axial_cylinder_frame() {
            // The placement origin is the start rim centre; the cylinder's
            // own bounds must agree with that reading.
            let bottom = frame.origin;
            let [radius, height] =
                [cylinder.circle().radius(), cylinder.height()].map(Rational::as_f64);
            let axis = frame.normal();
            let bounds = self.bounds().ok()?;
            for i in 0..3 {
                let end = bottom[i] + axis[i] * height;
                let reach = radius * (1. - axis[i] * axis[i]).max(0.).sqrt();
                if (bounds.min[i] - (bottom[i].min(end) - reach)).abs() > tolerance
                    || (bounds.max[i] - (bottom[i].max(end) + reach)).abs() > tolerance
                {
                    return None;
                }
            }
            if !same_span(bottom, axis, height) {
                return None;
            }
            let [x, y, _] = local(bottom);
            return Some(vec![Disk {
                center: [x, y],
                radius: cylinder.circle().radius().as_f64(),
                source: EntityId::from_uuid_bits(cylinder.feature()),
            }]);
        }
        let (body, frame) = self.normal_sweep_body()?;
        let (rings, height) = body.unpierced_profile()?;
        if !same_span(frame.origin, frame.normal(), height) {
            return None;
        }
        let source = self.rust_feature()?;
        let mut disks: Vec<Disk> = Vec::new();
        for edge in rings.iter().flatten() {
            let (center, radius) = match *edge {
                Edge::Arc {
                    center, clockwise, ..
                } if !clockwise => (center, edge.radius_mm()?),
                Edge::Circle {
                    center,
                    radius_mm,
                    clockwise: false,
                } => (center, radius_mm),
                _ => return None,
            };
            let [x, y, _] = local(frame.point(center));
            let found = disks.iter().any(|disk| {
                length2(sub2(disk.center, [x, y])) <= tolerance
                    && (disk.radius - radius).abs() <= tolerance
            });
            if !found {
                disks.push(Disk {
                    center: [x, y],
                    radius,
                    source,
                });
            }
        }
        // The disks must rebuild the sweep's own boundary edge for edge.
        let rebuilt = union(&disks, tolerance)?;
        let edges: Vec<&Edge> = rings.iter().flatten().collect();
        let arcs: Vec<&BoundaryArc> = rebuilt.loops.iter().flatten().collect();
        if rebuilt.loops.len() != rings.len() || arcs.len() != edges.len() {
            return None;
        }
        let near = |a: [f64; 2], b: [f64; 2]| length2(sub2(a, b)) <= tolerance;
        let planar = |p: [f64; 2]| {
            let [x, y, _] = local(frame.point(p));
            [x, y]
        };
        let mut matched = vec![false; edges.len()];
        for arc in arcs {
            let disk = &disks[arc.disk];
            let found = edges.iter().enumerate().position(|(i, edge)| {
                !matched[i]
                    && match (**edge, arc.ends) {
                        (
                            Edge::Circle {
                                center, radius_mm, ..
                            },
                            None,
                        ) => {
                            near(planar(center), disk.center)
                                && (radius_mm - disk.radius).abs() <= tolerance
                        }
                        (Edge::Arc { center, a, b, .. }, Some([from, to])) => {
                            near(planar(center), disk.center)
                                && near(planar(a), rebuilt.vertices[from])
                                && near(planar(b), rebuilt.vertices[to])
                        }
                        _ => false,
                    }
            })?;
            matched[found] = true;
        }
        Some(disks)
    }

    pub(super) fn parallel_disk_union(
        &self,
        other: &Self,
        operation: BooleanOperation,
        policy: GeometryPolicy,
    ) -> Option<Result<Option<Self>, GeometryError>> {
        if operation != BooleanOperation::Union
            || policy.boolean_fuzzy_mm != 0.
            || policy != self.policy
            || policy != other.policy
            || (self.axial_cylinder_frame().is_none() && other.axial_cylinder_frame().is_none())
        {
            return None;
        }
        let tolerance = policy.linear_tolerance_mm;
        let (cylinder, frame) = self
            .axial_cylinder_frame()
            .or_else(|| other.axial_cylinder_frame())?;
        let chart = Chart {
            frame,
            height: cylinder.height().as_f64(),
        };
        let mut disks = self.chart_disks(&chart, tolerance)?;
        disks.extend(other.chart_disks(&chart, tolerance)?);
        let union = union(&disks, tolerance)?;
        let areas: Vec<f64> = union
            .loops
            .iter()
            .map(|ring| loop_area(&union, &disks, ring))
            .collect();
        // One region: a single counter-clockwise outer loop.
        if areas.iter().filter(|&&a| a > 0.).count() != 1 {
            return None;
        }
        let mut rings: Vec<(f64, Vec<CurvedEdge2d>, Vec<Origin>)> = Vec::new();
        let mut roles = BTreeMap::<EntityId, usize>::new();
        for (ring, area) in union.loops.iter().zip(areas) {
            let mut edges = Vec::new();
            let mut origins = Vec::new();
            for arc in ring {
                let disk = &disks[arc.disk];
                edges.push(match arc.ends {
                    None => CurvedEdge2d::Circle {
                        center: disk.center,
                        radius_mm: disk.radius,
                    },
                    Some([from, to]) => CurvedEdge2d::Arc {
                        center: disk.center,
                        a: union.vertices[from],
                        b: union.vertices[to],
                        clockwise: false,
                    },
                });
                let count = roles.entry(disk.source).or_default();
                origins.push(Origin {
                    source: disk.source,
                    role: format!("side.cylinder.arc.{count}"),
                });
                *count += 1;
            }
            rings.push((area, edges, origins));
        }
        rings.sort_by(|a, b| b.0.total_cmp(&a.0));
        let origins: Vec<Origin> = rings.iter().flat_map(|r| r.2.clone()).collect();
        let mut edges = rings.into_iter().map(|r| r.1);
        let profile = CurvedProfile {
            frame: chart.frame,
            outer: edges.next()?,
            holes: edges.collect(),
        };
        let feature = EntityId::from_uuid_bits(cylinder.feature());
        Some(
            Self::rust_extrude_normal_profile_with_path(
                &profile,
                chart.height,
                feature,
                &origins,
                None,
                policy,
            )
            .map(Some),
        )
    }
}

#[cfg(test)]
mod tests;
