//! Source-owned finite curve display and selection, independent of mesh edges.
use super::*;
use multicad_core::references::{AnalyticGeometry, TopologyKind, TopologySnapshot};
const POINTS: usize = 1_000_000;
fn origin_ordinal(snapshot: &TopologySnapshot, role: &str) -> Result<u32, GeometryError> {
    snapshot
        .elements
        .iter()
        .find(|e| e.key.kind == TopologyKind::Edge && e.origins.iter().any(|o| o.role == role))
        .map(|e| e.key.ordinal)
        .ok_or_else(|| unsupported("Source curve topology handoff is incomplete"))
}
fn source_vertex(snapshot: &TopologySnapshot, role: &str) -> Result<Point, GeometryError> {
    snapshot
        .elements
        .iter()
        .find_map(|e| {
            if e.key.kind == TopologyKind::Vertex
                && e.origins.iter().any(|o| o.role == role)
                && let Some(AnalyticGeometry::Point { position_mm }) = e.geometry
            {
                return Some(position_mm);
            }
            None
        })
        .ok_or_else(|| unsupported("Source curve endpoint handoff is incomplete"))
}
fn placed(
    placement: Placement,
    p: [[f64; 2]; 3],
    existing: Option<Point>,
    reserve: f64,
) -> Result<(Point, f64), GeometryError> {
    let world = placement.point_bounds(p.map(|v| I::new(v[0], v[1])));
    if world.iter().any(|v| !v.valid()) {
        return Err(unsupported("Source display curve point is unresolved"));
    }
    let point = existing.unwrap_or_else(|| world.map(|v| v.lo.midpoint(v.hi)));
    let error = (0..3).fold(I::ZERO, |s, i| {
        s + I::point((world[i] - I::point(point[i])).magnitude())
    });
    if !error.valid() || error.hi > reserve {
        return Err(unsupported(
            "Source display curve reaches the precision floor",
        ));
    }
    Ok((point, error.hi))
}
struct SpanSampler {
    placement: Placement,
    deflection: f64,
    reserve: f64,
    remaining: usize,
}
impl SpanSampler {
    fn append(
        &mut self,
        second: f64,
        output: &mut Vec<Point>,
        start: Option<Point>,
        end: Option<Point>,
        point: impl FnMut(f64) -> Result<[[f64; 2]; 3], GeometryError>,
    ) -> Result<(), GeometryError> {
        let mut point = point;
        let count = (I::point(second) * I::point(self.placement.sigma_max)
            / (I::point(4.) * I::point(self.deflection)))
        .hi
        .sqrt()
        .next_up()
        .ceil()
        .max(1.);
        if !second.is_finite()
            || second < 0.
            || !count.is_finite()
            || count + 1. > self.remaining as f64
        {
            return Err(failure(
                GeometryErrorKind::Limit,
                "Actual source curves exceed display point budget",
            ));
        }
        let count = count as usize;
        let mut previous: Option<(f64, f64)> = None;
        for i in 0..=count {
            let t = i as f64 / count as f64;
            let existing = if i == 0 {
                output.last().copied().or(start)
            } else if i == count {
                end
            } else {
                None
            };
            let (p, error) = placed(self.placement, point(t)?, existing, self.reserve)?;
            if let Some((before, before_error)) = previous {
                let gap = I::point(t) - I::point(before);
                let bound = I::point(self.placement.sigma_max) * I::point(second) * gap.square()
                    / I::point(8.)
                    + I::point(error.max(before_error));
                if !bound.valid() || bound.hi > self.deflection {
                    return Err(unsupported(
                        "Actual source chord exceeds display deflection",
                    ));
                }
            }
            if i != 0 || output.is_empty() {
                output.push(p);
                self.remaining -= 1;
            }
            previous = Some((t, error));
        }
        Ok(())
    }
}
impl RepairedPhillipsBody {
    pub(crate) fn periodic_display_edges(
        &self,
        frame: PlaneFrame,
        deflection: f64,
    ) -> Result<Vec<DisplayEdge>, GeometryError> {
        let source = self.source.periodic().ok_or_else(|| {
            unsupported("Retained polynomial source curves have not been admitted")
        })?;
        let placement = self.placement(frame)?;
        let topology = self.periodic_topology(frame, self.feature, 0)?;
        let mut query = source
            .source_curve_query(16_000_000, None)
            .map_err(source_failure)?;
        let mut sampler = SpanSampler {
            placement,
            deflection,
            reserve: (self.policy.linear_tolerance_mm / 8.).min(deflection / 4.),
            remaining: POINTS,
        };
        let mut output = Vec::new();
        for band in 0..4 {
            let role = format!("thread/rail/{band}");
            let start = source_vertex(&topology, &format!("thread/cap/0/vertex/{band}"))?;
            let end = source_vertex(&topology, &format!("thread/cap/1/vertex/{band}"))?;
            let mut points = Vec::new();
            let second = query
                .thread_rail_second_bound_mm(band)
                .map_err(source_failure)?;
            sampler.append(second, &mut points, Some(start), Some(end), |t| {
                query.thread_rail_point_mm(band, t).map_err(source_failure)
            })?;
            output.push(DisplayEdge {
                index: origin_ordinal(&topology, &role)?,
                points_mm: points,
            });
            for cap in [1, 0] {
                let role = format!("thread/cap/{cap}/rim/{band}");
                let start = source_vertex(&topology, &format!("thread/cap/{cap}/vertex/{band}"))?;
                let end = source_vertex(
                    &topology,
                    &format!("thread/cap/{cap}/vertex/{}", (band + 1) % 4),
                )?;
                let mut points = Vec::new();
                let second = query
                    .thread_rim_second_bound_mm(band)
                    .map_err(source_failure)?;
                sampler.append(second, &mut points, Some(start), Some(end), |t| {
                    query
                        .thread_rim_point_mm(cap, band, t)
                        .map_err(source_failure)
                })?;
                output.push(DisplayEdge {
                    index: origin_ordinal(&topology, &role)?,
                    points_mm: points,
                });
            }
            if band != 3 {
                let role = format!("screw.head/thread.junction/{band}");
                let vertex = |band| {
                    if band == 0 || band == 3 {
                        format!("thread/cap/0/vertex/{band}")
                    } else {
                        format!("screw.head/thread.junction/vertex/{band}")
                    }
                };
                let start = source_vertex(&topology, &vertex(band))?;
                let end = source_vertex(&topology, &vertex((band + 1) % 4))?;
                let mut points = Vec::new();
                let second = query
                    .join_rim_second_bound_mm(band)
                    .map_err(source_failure)?;
                sampler.append(second, &mut points, Some(start), Some(end), |t| {
                    query.join_rim_point_mm(band, t).map_err(source_failure)
                })?;
                output.push(DisplayEdge {
                    index: origin_ordinal(&topology, &role)?,
                    points_mm: points,
                });
            }
        }
        for (station, segment) in source
            .head()
            .view()
            .profile_proposal()
            .segments()
            .enumerate()
        {
            let mut points = Vec::new();
            for band in 0..4 {
                let second = query
                    .head_rim_second_bound_mm(station, band)
                    .map_err(source_failure)?;
                let end = if band == 3 {
                    points.first().copied()
                } else {
                    None
                };
                sampler.append(second, &mut points, None, end, |t| {
                    query
                        .head_rim_point_mm(station, band, t)
                        .map_err(source_failure)
                })?;
            }
            let role = format!("screw.head/{}/start", segment.role.source_role());
            output.push(DisplayEdge {
                index: origin_ordinal(&topology, &role)?,
                points_mm: points,
            });
        }
        let cavity = source
            .cavity()
            .boundary_topology()
            .map_err(|e| unsupported(format!("Actual cavity display topology: {e:?}")))?;
        let mut ordinals = Vec::with_capacity(cavity.edges().len());
        let mut endpoints = Vec::with_capacity(cavity.edges().len());
        for edge in cavity.edges() {
            ordinals.push(origin_ordinal(
                &topology,
                &super::topology::cavity_curve_name(&cavity, edge.source_identity()),
            )?);
            let vertices = edge.vertices();
            endpoints.push([
                source_vertex(
                    &topology,
                    &format!(
                        "screw.drive/phillips/vertex/{:?}",
                        cavity.vertices()[vertices[0]].source_identity()
                    ),
                )?,
                source_vertex(
                    &topology,
                    &format!(
                        "screw.drive/phillips/vertex/{:?}",
                        cavity.vertices()[vertices[1]].source_identity()
                    ),
                )?,
            ]);
        }
        let mouth = source.cavity_mouth_expression_mm();
        let mouth = I::point(mouth[0]) + I::point(mouth[1]);
        let place = |mut p: [[f64; 2]; 3], existing: Option<Point>| {
            let z = I::new(p[2][0], p[2][1]) + mouth;
            p[2] = [z.lo, z.hi];
            placed(placement, p, existing, sampler.reserve)
        };
        output.extend(super::cavity_edges::edges(
            &cavity,
            &ordinals,
            &endpoints,
            deflection,
            placement.sigma_max,
            &mut sampler.remaining,
            &place,
        )?);
        output.sort_unstable_by_key(|e| e.index);
        Ok(output)
    }
}
