//! Exact polygon-prism and rounded-slot readings of plain normal sweeps.
//!
//! A sweep whose rings are straight polygons, optionally with full circular
//! holes, is the same exact solid as a polygon prism. Reading it that way
//! lets the polygon machining providers act on it. A sweep whose single ring
//! is a closed stadium is the exact envelope of a rounded slot tool. Both
//! readings require every coordinate to be an exact rational and refuse
//! openings, stock edits and arcs that are not part of a stadium.
use super::*;

/// The exact rational for `value`, or `None` when it has no exact
/// small-denominator representation.
fn exact(value: f64) -> Option<Rational> {
    let rational = rational_coordinate(value, 0.)?;
    (rational.as_f64() == value).then_some(rational)
}

impl NormalSweepBody {
    fn plain(&self) -> bool {
        self.openings.is_none() && self.stock_history.is_none()
    }

    /// The prism, edge origins and authored reversal flags of a straight
    /// sweep with optional circular holes, in polygon face order.
    fn polygon_reading(&self) -> Option<(PolygonPrism, Vec<Origin>, Vec<bool>)> {
        if !self.plain() {
            return None;
        }
        let mut rings = Vec::new();
        let mut circles = Vec::new();
        let mut straight = Vec::new();
        let mut round = Vec::new();
        for (index, range) in self.rings.iter().enumerate() {
            let sides = &self.sides[range.clone()];
            if let [side] = sides
                && index > 0
                && let Edge::Circle {
                    center, radius_mm, ..
                } = side.edge
            {
                let [x, y] = center.map(exact);
                circles.push(Circle::new([x?, y?], exact(radius_mm)?).ok()?);
                round.push((side.origin.clone(), false));
                continue;
            }
            let mut ring = Vec::with_capacity(sides.len());
            for side in sides {
                let Edge::Line { a, .. } = side.edge else {
                    return None;
                };
                let [x, y] = a.map(exact);
                ring.push([x?, y?]);
                straight.push((side.origin.clone(), side.reversed));
            }
            rings.push(ring);
        }
        let prism =
            PolygonPrism::with_circular_holes(rings, circles, exact(self.sweep.height_mm())?)
                .ok()?;
        let (origins, reversed) = straight.into_iter().chain(round).unzip();
        Some((prism, origins, reversed))
    }

    /// Centers, radius and signed height of a closed stadium sweep.
    fn stadium_reading(&self) -> Option<([[f64; 2]; 2], f64, f64)> {
        if !self.plain() {
            return None;
        }
        let [ring] = self.sweep.rings() else {
            return None;
        };
        let arcs: Vec<_> = ring
            .iter()
            .filter_map(|edge| match *edge {
                Edge::Arc { center, a, b, .. } => Some((center, a, b)),
                _ => None,
            })
            .collect();
        let lines = ring
            .iter()
            .filter(|edge| matches!(edge, Edge::Line { .. }))
            .count();
        let [(c0, a0, b0), (c1, a1, b1)] = arcs[..] else {
            return None;
        };
        if ring.len() != 4 || lines != 2 {
            return None;
        }
        let distance = |p: [f64; 2], q: [f64; 2]| (p[0] - q[0]).hypot(p[1] - q[1]);
        let radius = distance(c0, a0);
        // Each arc is an exact half turn: its ends are antipodal about the
        // center, and both arcs share one radius, so the two straight edges
        // are the parallel tangents of the slot.
        let half_turn = |c: [f64; 2], a: [f64; 2], b: [f64; 2]| {
            a[0] + b[0] == 2. * c[0] && a[1] + b[1] == 2. * c[1]
        };
        if !half_turn(c0, a0, b0)
            || !half_turn(c1, a1, b1)
            || distance(c0, b0) != radius
            || distance(c1, a1) != radius
            || distance(c1, b1) != radius
            || c0 == c1
        {
            return None;
        }
        let height = self.sweep.height_mm();
        Some((
            [c0, c1],
            radius,
            if self.reverse { -height } else { height },
        ))
    }
}

impl Solid {
    /// The polygon-prism reading of a plain straight normal sweep.
    pub(in crate::geometry) fn normal_sweep_polygon_stock(&self) -> Option<Self> {
        let (body, frame) = self.normal_sweep_body()?;
        let (prism, origins, reversed) = body.polygon_reading()?;
        Self::from_polygon_prism_sweep(
            prism,
            frame,
            body.reverse,
            body.feature,
            &origins,
            &reversed,
            self.policy,
        )
        .ok()
    }

    /// The rounded slot tool whose swept envelope is this stadium sweep.
    pub(in crate::geometry) fn normal_sweep_slot_tool(&self) -> Option<RoundedSlotTool> {
        let (body, frame) = self.normal_sweep_body()?;
        let (centers_uv, radius_mm, height) = body.stadium_reading()?;
        let normal = frame.normal();
        Some(RoundedSlotTool {
            frame,
            centers_uv,
            radius_mm,
            displacement_mm: normal.map(|n| n * height),
            origin: Origin {
                source: body.feature,
                role: "slot".into(),
            },
        })
    }
}

/// The placement frame, signed height, feature and section rings (each edge
/// with its authored origin) of a plain normal sweep.
pub(in crate::geometry) type SweepSection =
    (PlaneFrame, f64, EntityId, Vec<Vec<(CurvedEdge2d, Origin)>>);

impl Solid {
    /// The section rings, with each edge's authored origin, of a plain
    /// normal sweep: the outer ring first and then the holes, in the sweep's
    /// normalized orientation. Also the frame, signed height and feature.
    pub(in crate::geometry) fn normal_sweep_section(&self) -> Option<SweepSection> {
        let (body, frame) = self.normal_sweep_body()?;
        if !body.plain() {
            return None;
        }
        let rings = body
            .rings
            .iter()
            .map(|range| {
                body.sides[range.clone()]
                    .iter()
                    .map(|side| {
                        let edge = match side.edge {
                            Edge::Line { a, b } => CurvedEdge2d::Line { a, b },
                            Edge::Arc {
                                center,
                                a,
                                b,
                                clockwise,
                            } => CurvedEdge2d::Arc {
                                center,
                                a,
                                b,
                                clockwise,
                            },
                            Edge::Circle {
                                center, radius_mm, ..
                            } => CurvedEdge2d::Circle { center, radius_mm },
                        };
                        (edge, side.origin.clone())
                    })
                    .collect()
            })
            .collect();
        Some((frame, body.z(true), body.feature, rings))
    }
}
