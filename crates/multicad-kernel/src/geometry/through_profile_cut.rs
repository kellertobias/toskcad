//! Through cuts that change only the section of a straight prism.
//!
//! A straight prism (a box with through bores, a plain polygon prism or a
//! plain normal sweep) has one section along its whole height. A tool that
//! passes both caps and is itself a straight prism along the same axis, such
//! as a cylinder or a box, removes its own section at every height. Two such
//! section changes are certified here, and the result is the exact normal
//! sweep of the new section:
//!
//! * a tool section strictly inside the outer ring and clear of every hole
//!   becomes a new hole;
//! * a cylinder centered on a straight hole edge whose length is its
//!   diameter, between two perpendicular edges that run into the hole, rounds
//!   that end of the hole exactly: the edge becomes the outward half circle,
//!   which meets both neighbours tangentially at their shared corners.
//!
//! Every clearance is proved with a margin of the linear tolerance. Other
//! overlaps, tangencies and tools that stop inside the stock are left to the
//! other providers.
use super::*;
use crate::{CurvedEdge2d, CurvedProfile};

type Point2 = [f64; 2];
type Ring = Vec<(CurvedEdge2d, Origin)>;

/// The section of a straight prism stock in its own placement chart.
struct Stock {
    frame: PlaneFrame,
    /// Signed sweep height along the frame normal.
    height: f64,
    feature: EntityId,
    /// The counter-clockwise outer ring, then clockwise or circular holes.
    rings: Vec<Ring>,
}

enum Tool {
    Circle(Point2, f64, Origin),
    /// A clockwise rectangle.
    Ring(Ring),
}

fn sub(a: Point2, b: Point2) -> Point2 {
    [a[0] - b[0], a[1] - b[1]]
}

fn length(a: Point2) -> f64 {
    a[0].portable_hypot(a[1])
}

fn segment_distance(p: Point2, a: Point2, b: Point2) -> f64 {
    let d = sub(b, a);
    let t = (dot2(sub(p, a), d) / dot2(d, d)).clamp(0., 1.);
    length(sub(p, [a[0] + t * d[0], a[1] + t * d[1]]))
}

fn dot2(a: Point2, b: Point2) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}

/// A lower bound of the distance from `p` to the edge.
fn edge_distance(p: Point2, edge: &CurvedEdge2d) -> f64 {
    match *edge {
        CurvedEdge2d::Line { a, b } => segment_distance(p, a, b),
        CurvedEdge2d::Arc { center, a, .. } => {
            (length(sub(p, center)) - length(sub(a, center))).abs()
        }
        CurvedEdge2d::Circle { center, radius_mm } => (length(sub(p, center)) - radius_mm).abs(),
        CurvedEdge2d::BSpline(_) => 0.,
    }
}

/// Conservative chart bounds of an edge.
fn edge_bounds(edge: &CurvedEdge2d) -> Option<[Point2; 2]> {
    Some(match *edge {
        CurvedEdge2d::Line { a, b } => [
            [a[0].min(b[0]), a[1].min(b[1])],
            [a[0].max(b[0]), a[1].max(b[1])],
        ],
        CurvedEdge2d::Arc { center, a, .. } => {
            let r = length(sub(a, center));
            [
                [center[0] - r, center[1] - r],
                [center[0] + r, center[1] + r],
            ]
        }
        CurvedEdge2d::Circle { center, radius_mm } => [
            [center[0] - radius_mm, center[1] - radius_mm],
            [center[0] + radius_mm, center[1] + radius_mm],
        ],
        CurvedEdge2d::BSpline(_) => return None,
    })
}

fn ring_bounds(ring: &Ring) -> Option<[Point2; 2]> {
    ring.iter().try_fold(
        [[f64::INFINITY; 2], [f64::NEG_INFINITY; 2]],
        |[lo, hi], (edge, _)| {
            let [a, b] = edge_bounds(edge)?;
            Some([
                [lo[0].min(a[0]), lo[1].min(a[1])],
                [hi[0].max(b[0]), hi[1].max(b[1])],
            ])
        },
    )
}

fn apart(a: [Point2; 2], b: [Point2; 2], tolerance: f64) -> bool {
    (0..2).any(|i| a[1][i] + tolerance < b[0][i] || b[1][i] + tolerance < a[0][i])
}

/// Whether the box lies inside a convex straight outer ring with margin.
fn inside_outer(outer: &Ring, [lo, hi]: [Point2; 2], tolerance: f64) -> bool {
    let corners = [lo, [hi[0], lo[1]], hi, [lo[0], hi[1]]];
    outer.iter().all(|(edge, _)| {
        let CurvedEdge2d::Line { a, b } = *edge else {
            return false;
        };
        let d = sub(b, a);
        let n = length(d);
        // The counter-clockwise ring keeps its material on the left.
        corners
            .iter()
            .all(|&p| (d[0] * (p[1] - a[1]) - d[1] * (p[0] - a[0])) / n > tolerance)
    }) && {
        let points: Vec<Point2> = outer
            .iter()
            .map(|(edge, _)| match *edge {
                CurvedEdge2d::Line { a, .. } => a,
                _ => [0.; 2],
            })
            .collect();
        let n = points.len();
        (0..n).all(|i| {
            let [a, b, c] = [points[i], points[(i + 1) % n], points[(i + 2) % n]];
            (b[0] - a[0]) * (c[1] - b[1]) - (b[1] - a[1]) * (c[0] - b[0]) > 0.
        })
    }
}

impl Solid {
    fn through_stock(&self) -> Option<Stock> {
        if let Some((frame, height, feature, rings)) = self.normal_sweep_section() {
            return Some(Stock {
                frame,
                height,
                feature,
                rings,
            });
        }
        let line = |a: Point2, b: Point2| CurvedEdge2d::Line { a, b };
        if let Some((prism, frame)) = self.rust_prism_frame() {
            if prism.rounded_end().is_some() {
                return None;
            }
            let feature = EntityId::from_uuid_bits(prism.feature());
            let [l, w, h] = prism.dimensions().map(Rational::as_f64);
            let corners = [[0., 0.], [l, 0.], [l, w], [0., w]];
            let outer = (0..4)
                .map(|i| {
                    (
                        line(corners[i], corners[(i + 1) % 4]),
                        Origin {
                            source: feature,
                            role: format!("profile.edge.{i}"),
                        },
                    )
                })
                .collect();
            let holes = prism.inner_loops().map(|(source, circle)| {
                vec![(
                    CurvedEdge2d::Circle {
                        center: circle.center().map(Rational::as_f64),
                        radius_mm: circle.radius().as_f64(),
                    },
                    Origin {
                        source: EntityId::from_uuid_bits(source),
                        role: "side.cylinder".into(),
                    },
                )]
            });
            return Some(Stock {
                frame,
                height: h,
                feature,
                rings: std::iter::once(outer).chain(holes).collect(),
            });
        }
        let placed = self.polygon_prism()?;
        if !placed.pockets.is_empty() || placed.bore_treatments.is_some() {
            return None;
        }
        let mut origins = placed.edges.iter().map(|(origin, _)| origin.clone());
        let mut rings: Vec<Ring> = placed
            .prism
            .rings()
            .iter()
            .map(|ring| {
                let points: Vec<Point2> = ring.iter().map(|p| p.map(Rational::as_f64)).collect();
                (0..points.len())
                    .map(|i| {
                        let edge = line(points[i], points[(i + 1) % points.len()]);
                        Some((edge, origins.next()?))
                    })
                    .collect::<Option<_>>()
            })
            .collect::<Option<_>>()?;
        for circle in placed.prism.circular_holes() {
            rings.push(vec![(
                CurvedEdge2d::Circle {
                    center: circle.center().map(Rational::as_f64),
                    radius_mm: circle.radius().as_f64(),
                },
                origins.next()?,
            )]);
        }
        let height = placed.prism.height().as_f64();
        Some(Stock {
            frame: placed.frame,
            height: if placed.reverse_sweep {
                -height
            } else {
                height
            },
            feature: placed.feature,
            rings,
        })
    }

    /// This solid as a tool through both caps of `stock`, in its chart.
    fn through_tool(&self, stock: &Stock, tolerance: f64) -> Option<Tool> {
        let frame = stock.frame;
        let normal = frame.normal();
        let chart = |p: [f64; 3]| -> Point2 {
            let d = std::array::from_fn(|i| p[i] - frame.origin[i]);
            [dot(d, frame.u), dot(d, frame.v)]
        };
        let along = |p: [f64; 3]| {
            let d: [f64; 3] = std::array::from_fn(|i| p[i] - frame.origin[i]);
            dot(d, normal)
        };
        let parallel = |axis: [f64; 3]| {
            let c = cross(axis, normal);
            dot(c, c).sqrt() <= 64. * f64::EPSILON
        };
        let through = |a: f64, b: f64| {
            let [lo, hi] = [stock.height.min(0.), stock.height.max(0.)];
            a.min(b) < lo - tolerance && a.max(b) > hi + tolerance
        };
        if let Some((cylinder, placement)) = self.axial_cylinder_frame() {
            let axis = placement.normal();
            let end: [f64; 3] =
                std::array::from_fn(|i| placement.origin[i] + cylinder.height().as_f64() * axis[i]);
            if !parallel(axis) || !through(along(placement.origin), along(end)) {
                return None;
            }
            let origin = self.cylinder_sketch().map_or_else(
                || Origin {
                    source: EntityId::from_uuid_bits(cylinder.feature()),
                    role: "side.cylinder".into(),
                },
                |sketch| sketch.origin.clone(),
            );
            return Some(Tool::Circle(
                chart(placement.origin),
                cylinder.circle().radius().as_f64(),
                origin,
            ));
        }
        let (prism, placement) = self.rust_prism_frame()?;
        if prism.rounded_end().is_some() || prism.inner_loops().len() != 0 {
            return None;
        }
        let size = prism.dimensions().map(Rational::as_f64);
        let axes = [placement.u, placement.v, placement.normal()];
        let k = (0..3).find(|&k| parallel(axes[k]))?;
        let [i, j] = [(k + 1) % 3, (k + 2) % 3];
        let at = |a: f64, b: f64, c: f64| -> [f64; 3] {
            std::array::from_fn(|n| {
                placement.origin[n] + a * axes[i][n] + b * axes[j][n] + c * axes[k][n]
            })
        };
        if !through(along(at(0., 0., 0.)), along(at(0., 0., size[k]))) {
            return None;
        }
        let mut corners = [
            chart(at(0., 0., 0.)),
            chart(at(size[i], 0., 0.)),
            chart(at(size[i], size[j], 0.)),
            chart(at(0., size[j], 0.)),
        ];
        let area: f64 = (0..4)
            .map(|n| {
                let [a, b] = [corners[n], corners[(n + 1) % 4]];
                a[0] * b[1] - a[1] * b[0]
            })
            .sum();
        // Holes are clockwise.
        if area > 0. {
            corners.reverse();
        }
        let feature = EntityId::from_uuid_bits(prism.feature());
        Some(Tool::Ring(
            (0..4)
                .map(|n| {
                    (
                        CurvedEdge2d::Line {
                            a: corners[n],
                            b: corners[(n + 1) % 4],
                        },
                        Origin {
                            source: feature,
                            role: format!("profile.edge.{n}"),
                        },
                    )
                })
                .collect(),
        ))
    }

    /// The hole and edge a circle rounds off, when it does so exactly.
    fn rounded_end(
        stock: &Stock,
        center: Point2,
        radius: f64,
        tolerance: f64,
    ) -> Option<[usize; 2]> {
        let mut found = None;
        for (index, ring) in stock.rings.iter().enumerate().skip(1) {
            let n = ring.len();
            for e in 0..n {
                let CurvedEdge2d::Line { a, b } = ring[e].0 else {
                    continue;
                };
                let d = sub(b, a);
                if [a[0] + b[0], a[1] + b[1]] != [2. * center[0], 2. * center[1]]
                    || dot2(d, d) != 4. * radius * radius
                {
                    continue;
                }
                let (CurvedEdge2d::Line { a: p, b: pa }, CurvedEdge2d::Line { a: nb, b: q }) =
                    (ring[(e + n - 1) % n].0.clone(), ring[(e + 1) % n].0.clone())
                else {
                    return None;
                };
                // Clockwise holes lie right of `d`; the neighbours must run
                // perpendicular into the hole from both corners.
                let left = [-d[1], d[0]];
                if pa != a
                    || nb != b
                    || dot2(sub(a, p), d) != 0.
                    || dot2(sub(q, b), d) != 0.
                    || dot2(sub(p, a), left) >= 0.
                    || dot2(sub(q, b), left) >= 0.
                    || found.is_some()
                {
                    return None;
                }
                found = Some([index, e]);
            }
        }
        let [hole, e] = found?;
        let n = stock.rings[hole].len();
        let neighbours = [(e + n - 1) % n, e, (e + 1) % n];
        let clear = stock.rings.iter().enumerate().all(|(index, ring)| {
            ring.iter().enumerate().all(|(i, (edge, _))| {
                (index == hole && neighbours.contains(&i))
                    || edge_distance(center, edge) > radius + tolerance
            })
        });
        clear.then_some([hole, e])
    }

    pub(super) fn through_profile_cut(
        &self,
        tool: &Self,
        operation: BooleanOperation,
        policy: GeometryPolicy,
    ) -> Option<Result<Option<Self>, GeometryError>> {
        if operation != BooleanOperation::Subtract
            || policy.boolean_fuzzy_mm != 0.
            || policy != self.policy
            || policy != tool.policy
        {
            return None;
        }
        let tolerance = policy.linear_tolerance_mm;
        let mut stock = self.through_stock()?;
        let cut = tool.through_tool(&stock, tolerance)?;
        let ring = match cut {
            Tool::Circle(center, radius, origin) => {
                if let Some([hole, e]) = Self::rounded_end(&stock, center, radius, tolerance) {
                    let CurvedEdge2d::Line { a, b } = stock.rings[hole][e].0 else {
                        unreachable!("a rounded end replaces a line");
                    };
                    stock.rings[hole][e] = (
                        CurvedEdge2d::Arc {
                            center,
                            a,
                            b,
                            clockwise: true,
                        },
                        origin,
                    );
                    None
                } else {
                    Some(vec![(
                        CurvedEdge2d::Circle {
                            center,
                            radius_mm: radius,
                        },
                        origin,
                    )])
                }
            }
            Tool::Ring(ring) => Some(ring),
        };
        if let Some(ring) = ring {
            let bounds = ring_bounds(&ring)?;
            if !inside_outer(&stock.rings[0], bounds, tolerance)
                || !stock.rings[1..]
                    .iter()
                    .all(|hole| ring_bounds(hole).is_some_and(|b| apart(bounds, b, tolerance)))
            {
                return None;
            }
            stock.rings.push(ring);
        }
        let origins: Vec<_> = stock
            .rings
            .iter()
            .flatten()
            .map(|(_, origin)| origin.clone())
            .collect();
        if origins
            .iter()
            .collect::<std::collections::BTreeSet<_>>()
            .len()
            != origins.len()
        {
            return None;
        }
        let mut rings = stock
            .rings
            .into_iter()
            .map(|ring| ring.into_iter().map(|(edge, _)| edge).collect());
        let profile = CurvedProfile {
            frame: stock.frame,
            outer: rings.next()?,
            holes: rings.collect(),
        };
        Some(
            Self::rust_extrude_normal_profile_with_path(
                &profile,
                stock.height,
                stock.feature,
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
