//! Display/selection curves from original D11 boundary coefficients and the
//! actual repaired-head circles, independent of any triangle mesh or recipe.
use super::*;
use multicad_core::references::{AnalyticGeometry, TopologyKind, TopologySnapshot};
const POINTS: usize = 1_000_000;
const METADATA_BYTES: usize = 8 * 1024 * 1024;
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

fn curve_point(poles: [Point; 12], t: f64) -> [[f64; 2]; 3] {
    let mut controls = poles.map(|p| p.map(I::point));
    let parameter = I::point(t);
    let complement = I::ONE - parameter;
    for degree in (1..12).rev() {
        for i in 0..degree {
            controls[i] = std::array::from_fn(|k| {
                complement * controls[i][k] + parameter * controls[i + 1][k]
            });
        }
    }
    controls[0].map(|v| [v.lo, v.hi])
}
fn curve_second(poles: [Point; 12]) -> f64 {
    (0..10)
        .map(|i| {
            (0..3)
                .fold(I::ZERO, |sum, k| {
                    let second = I::point(110.)
                        * (I::point(poles[i + 2][k]) - I::point(2.) * I::point(poles[i + 1][k])
                            + I::point(poles[i][k]));
                    sum + I::point(second.magnitude())
                })
                .hi
        })
        .fold(0., f64::max)
}
fn append_curve(
    spans: usize,
    pole: impl Fn(usize) -> Point,
    sampler: &mut SpanSampler,
    start: Point,
    end: Point,
) -> Result<Vec<Point>, GeometryError> {
    if spans == 0 || spans >= sampler.remaining {
        return Err(failure(
            GeometryErrorKind::Limit,
            "Original D11 source spans exceed display budget",
        ));
    }
    let mut output = Vec::new();
    for span in 0..spans {
        // Each original Bezier span has its OWN normalized unit interval.
        // Its coefficient second differences already include its chart scale;
        // no extra original-span count factor belongs in this chord proof.
        let controls = std::array::from_fn(|i| pole(11 * span + i));
        sampler.append(
            curve_second(controls),
            &mut output,
            if span == 0 { Some(start) } else { None },
            if span + 1 == spans { Some(end) } else { None },
            |t| Ok(curve_point(controls, t)),
        )?;
    }
    Ok(output)
}
fn circle_point(radius: f64, z: f64, quarter: usize, t: f64) -> [[f64; 2]; 3] {
    // Exact rational quadratic circle: w=sqrt(1/2), weighted controls
    // (1,0),(w,w),(0,1). The numerator identity Nx²+Ny²=D² follows
    // from w²=1/2. Directed weights enclose the SAME exact circle.
    let parameter = I::point(t);
    let complement = I::ONE - parameter;
    let middle = I::point(2.) * I::point(0.5).sqrt() * complement * parameter;
    let denominator = complement.square() + middle + parameter.square();
    let x = (complement.square() + middle) / denominator;
    let y = (parameter.square() + middle) / denominator;
    let [x, y] = match quarter {
        0 => [x, y],
        1 => [-y, x],
        2 => [-x, -y],
        _ => [y, -x],
    };
    [I::point(radius) * x, I::point(radius) * y, I::point(z)].map(|v| [v.lo, v.hi])
}
fn circle_second(radius: f64) -> f64 {
    let weight = I::point(0.5).sqrt();
    // Positive Bernstein denominator D>=w. For each numerator:
    // (N/D)''=N''/D-ND''/D²-2N'D'/D²+2N(D')²/D³.
    // Whole coefficient hulls bound every term on the complete quadrant.
    let lower = I::point(weight.lo);
    let denominator_first = I::point((I::point(2.) * (I::ONE - weight)).magnitude());
    let denominator_second = I::point((I::point(4.) * (I::ONE - weight)).magnitude());
    let numerator_first = I::point((I::point(2.) * weight).magnitude());
    let numerator_second = I::point((I::point(2.) * (I::ONE - I::point(2.) * weight)).magnitude());
    let per = numerator_second / lower
        + denominator_second / lower.square()
        + I::point(2.) * numerator_first * denominator_first / lower.square()
        + I::point(2.) * denominator_first.square() / (lower.square() * lower);
    (I::point(2.) * I::point(radius) * per).hi
}
fn append_circle(
    radius: f64,
    z: f64,
    sampler: &mut SpanSampler,
) -> Result<Vec<Point>, GeometryError> {
    let mut output = Vec::new();
    let second = circle_second(radius);
    for quarter in 0..4 {
        let end = if quarter == 3 {
            output.first().copied()
        } else {
            None
        };
        sampler.append(second, &mut output, None, end, |t| {
            Ok(circle_point(radius, z, quarter, t))
        })?;
    }
    Ok(output)
}
impl RepairedPhillipsBody {
    pub(crate) fn degree11_display_edges(
        &self,
        frame: PlaneFrame,
        deflection: f64,
    ) -> Result<Vec<DisplayEdge>, GeometryError> {
        if !deflection.is_finite() || deflection <= 0. {
            return Err(input("Invalid D11 edge display deflection"));
        }
        let source = self
            .source
            .degree11()
            .ok_or_else(|| unsupported("No actual D11 material curves"))?;
        self.validate_frame(frame)?;
        let placement = self.placement(frame)?;
        let available = BYTES
            .checked_sub(self.resident_bytes())
            .and_then(|n| n.checked_sub(METADATA_BYTES))
            .ok_or_else(query::exhausted)?;
        // Four point payloads cover geometric Vec growth and old/new allocation
        // overlap. This is logical-capacity admission, not an allocator/RSS claim.
        let points = (available / (4 * std::mem::size_of::<Point>())).min(POINTS);
        if points == 0 {
            return Err(query::exhausted());
        }
        let cavity = source
            .cavity()
            .boundary_topology()
            .map_err(|e| unsupported(format!("Actual D11 Cutter curves: {e:?}")))?;
        let h = source.head().view().segment_count();
        let n = source.thread_faces().len();
        // Actual accepted head count, never an assumed five-segment recipe.
        // At most two bounded Cutter atlases overlap while topology is built.
        // 8KiB per actual graph item covers strings/sets/frames/descriptor vectors.
        let items = source
            .face_count()
            .checked_add(5 * n + h + 1)
            .and_then(|m| m.checked_add(cavity.edges().len() + cavity.vertices().len()))
            .ok_or_else(query::exhausted)?;
        let metadata = items
            .checked_mul(8192)
            .and_then(|m| m.checked_add(2 * 1024 * 1024 + 64 * 1024))
            .ok_or_else(query::exhausted)?;
        if metadata > METADATA_BYTES {
            return Err(query::exhausted());
        }
        let topology = self.degree11_topology(frame, self.feature, 0)?;
        let mut sampler = SpanSampler {
            placement,
            deflection,
            reserve: (self.policy.linear_tolerance_mm / 8.).min(deflection / 4.),
            remaining: points,
        };
        let mut output = Vec::new();
        for (index, face) in source.thread_faces().iter().enumerate() {
            let net = face.source();
            let [width, height] = net.dimensions();
            let start = source_vertex(&topology, &format!("thread/cap/0/vertex/{index}"))?;
            let end = source_vertex(&topology, &format!("thread/cap/1/vertex/{index}"))?;
            output.push(DisplayEdge {
                index: origin_ordinal(&topology, &format!("thread/rail/{index}"))?,
                points_mm: append_curve(
                    net.spans()[1],
                    |v| net.poles()[v],
                    &mut sampler,
                    start,
                    end,
                )?,
            });
            for cap in [1, 0] {
                let start = source_vertex(&topology, &format!("thread/cap/{cap}/vertex/{index}"))?;
                let end = source_vertex(
                    &topology,
                    &format!("thread/cap/{cap}/vertex/{}", (index + 1) % n),
                )?;
                let station = if cap == 0 { 0 } else { height - 1 };
                output.push(DisplayEdge {
                    index: origin_ordinal(&topology, &format!("thread/cap/{cap}/rim/{index}"))?,
                    points_mm: append_curve(
                        (width - 1) / 11,
                        |u| net.poles()[u * height + station],
                        &mut sampler,
                        start,
                        end,
                    )?,
                });
            }
        }
        let head = source.head().view();
        for (station, segment) in head.profile_proposal().segments().enumerate() {
            let [radius, z] = head
                .station(station)
                .ok_or_else(|| unsupported("Actual head station missing"))?;
            output.push(DisplayEdge {
                index: origin_ordinal(
                    &topology,
                    &format!("screw.head/{}/start", segment.role.source_role()),
                )?,
                points_mm: append_circle(radius, z, &mut sampler)?,
            });
        }
        let [radius, z] = head
            .station(h)
            .ok_or_else(|| unsupported("Actual head bottom missing"))?;
        output.push(DisplayEdge {
            index: origin_ordinal(&topology, "screw.head/thread.junction")?,
            points_mm: append_circle(radius, z, &mut sampler)?,
        });
        let mut ordinals = Vec::with_capacity(cavity.edges().len());
        let mut endpoints = Vec::with_capacity(cavity.edges().len());
        for edge in cavity.edges() {
            ordinals.push(origin_ordinal(
                &topology,
                &super::topology::cavity_curve_name(&cavity, edge.source_identity()),
            )?);
            let vertices = edge.vertices();
            let endpoint = |v: usize| {
                source_vertex(
                    &topology,
                    &format!(
                        "screw.drive/phillips/vertex/{:?}",
                        cavity.vertices()[v].source_identity()
                    ),
                )
            };
            endpoints.push([endpoint(vertices[0])?, endpoint(vertices[1])?]);
        }
        let mouth = I::point(source.chart().outer_top_z);
        let place = |mut p: [[f64; 2]; 3], existing: Option<Point>| {
            let z = I::new(p[2][0], p[2][1]) + mouth;
            p[2] = [z.lo, z.hi];
            placed(placement, p, existing, sampler.reserve)
        };
        output.extend(super::super::periodic::cavity_edges::edges(
            &cavity,
            &ordinals,
            &endpoints,
            deflection,
            placement.sigma_max,
            &mut sampler.remaining,
            &place,
        )?);
        output.sort_unstable_by_key(|edge| edge.index);
        if output.len()
            != topology
                .elements
                .iter()
                .filter(|e| e.key.kind == TopologyKind::Edge)
                .count()
        {
            return Err(unsupported("Actual D11 display edge inventory incomplete"));
        }
        Ok(output)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn contains(bounds: [[f64; 2]; 3], point: Point) {
        for k in 0..3 {
            assert!(
                bounds[k][0] <= point[k] && point[k] <= bounds[k][1],
                "coordinate{k}: {:?} vs{}",
                bounds[k],
                point[k]
            );
        }
    }
    #[test]
    fn original_n12_coefficients_enclose_independent_integer_polynomial_and_second_derivative() {
        // Bernstein11 controls i(i-1), i, i(i-1)(i-2) define EXACT
        // P(t)=(110t²,11t,990t³); every control is an exact small integer.
        let poles = std::array::from_fn(|i| {
            let k = i as i32;
            [
                (k * (k - 1)) as f64,
                k as f64,
                (k * (k - 1) * (k - 2)) as f64,
            ]
        });
        for t in [0., 0.125, 0.5, 0.75, 1.] {
            let bounds = curve_point(poles, t);
            contains(bounds, [110. * t * t, 11. * t, 990. * t * t * t]);
            assert!(bounds.iter().all(|v| v[1] - v[0] < 1e-8));
        }
        // ||P''||_1=220+5940t attains6160 on the whole interval.
        let second = curve_second(poles);
        assert!((6160. ..6160. * (1. + 1e-12)).contains(&second));
        let line = std::array::from_fn(|i| [i as f64, 2. * i as f64, 3. * i as f64]);
        assert!(curve_second(line) < 1e-9);
    }
    #[test]
    fn rational_head_quadrants_preserve_actual_circle_endpoints_and_positive_denominator() {
        let rotate = |p: [f64; 2], quarter| match quarter {
            0 => p,
            1 => [-p[1], p[0]],
            2 => [-p[0], -p[1]],
            _ => [p[1], -p[0]],
        };
        for quarter in 0..4 {
            for (t, p) in [(0., [2., 0.]), (1., [0., 2.])] {
                let p = rotate(p, quarter);
                contains(circle_point(2., 3., quarter, t), [p[0], p[1], 3.]);
            }
            for t in [0., 0.125, 0.5, 0.875, 1.] {
                let point = circle_point(2., 3., quarter, t).map(|v| I::new(v[0], v[1]));
                assert!(point.iter().all(|v| v.valid()));
                let radius = point[0].square() + point[1].square();
                assert!(radius.lo <= 4. && 4. <= radius.hi);
                assert!(radius.width() < 1e-10);
            }
        }
        let bound = circle_second(2.);
        // At the quadrant start, independent exact differentiation gives
        // ||P''||_1=2*(2+4sqrt(1/2)-2)=8sqrt(1/2).
        assert!(bound >= 8. * 0.5_f64.sqrt() && bound < 100.);
        assert!((bound - 4. * circle_second(0.5)).abs() < 1e-12);
    }
}
