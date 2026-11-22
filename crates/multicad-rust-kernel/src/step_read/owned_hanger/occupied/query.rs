//! Conservative point queries against the retained decoded, healed field.
use super::*;
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PointLocation {
    Inside,
    Outside,
    /// Includes the sewing collar and any unresolved interval predicate.
    Unresolved,
}
#[derive(Clone, Copy)]
pub struct QueryOptions {
    pub maximum_work: usize,
    pub maximum_bytes: usize,
}
impl Default for QueryOptions {
    fn default() -> Self {
        Self {
            maximum_work: 2_000_000,
            maximum_bytes: 64 * 1024 * 1024,
        }
    }
}
pub struct PointQuery {
    location: PointLocation,
    work: usize,
}
impl PointQuery {
    pub fn location(&self) -> PointLocation {
        self.location
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
}
struct Node {
    xy: [Polynomial; 2],
    a: I,
    z: I,
    depth: usize,
}
fn hull(a: &[I]) -> I {
    I::new(
        a.iter().map(|v| v.lo).fold(f64::INFINITY, f64::min),
        a.iter().map(|v| v.hi).fold(f64::NEG_INFINITY, f64::max),
    )
}
fn rational(values: &[I], weights: &[I], u: I, b: &mut Budget<'_>) -> Result<I, Error> {
    b.charge(values.len())?;
    let weighted = values
        .iter()
        .zip(weights)
        .map(|(&v, &w)| v * w)
        .collect::<Polynomial>();
    let n = poly::evaluate(&weighted, u, b)?;
    let w = poly::evaluate(weights, u, b)?;
    if w.lo <= 0. {
        return Err(Error::Unsupported);
    }
    Ok(n / w)
}
fn row_radius(strip: &PolarStrip, target: I, b: &mut Budget<'_>) -> Result<I, Error> {
    let mut lo = 0.;
    let mut hi = 1.;
    for _ in 0..48 {
        b.charge(1)?;
        let mid = (lo + hi) * 0.5;
        let q = rational(&strip.phase, &strip.weights, p(mid), b)?;
        if q.hi < target.lo {
            lo = mid
        } else if q.lo > target.hi {
            hi = mid
        } else {
            break;
        }
    }
    let mut n = Vec::with_capacity(strip.radius.len());
    for (&r, &w) in strip.radius.iter().zip(&strip.weights) {
        b.charge(1)?;
        n.push(r * w)
    }
    let nr = poly::restrict(&n, p(lo), p(hi), b)?;
    let wr = poly::restrict(&strip.weights, p(lo), p(hi), b)?;
    let w = hull(&wr);
    if w.lo <= 0. {
        return Err(Error::Unsupported);
    }
    Ok(hull(&nr) / w)
}
fn compare(lambda: I, radius: I, margin: I) -> PointLocation {
    if (lambda + margin).hi < radius.lo {
        PointLocation::Inside
    } else if (lambda - margin).lo > radius.hi {
        PointLocation::Outside
    } else {
        PointLocation::Unresolved
    }
}
impl HealedMaterial<'_, '_, '_, '_, '_> {
    /// Uses actual decoded angular controls and rational meridian rows. Nominal
    /// source dimensions and STEP headers are never used to reconstruct a part.
    pub fn classify_point(
        &self,
        point: [f64; 3],
        options: QueryOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<PointQuery, Error> {
        if point.iter().any(|v| !v.is_finite())
            || options.maximum_work == 0
            || options.maximum_work > 16_000_000
            || options.maximum_bytes == 0
            || options.maximum_bytes > 64 * 1024 * 1024
        {
            return Err(Error::Invalid);
        }
        let mut b = Budget {
            work: 0,
            remaining: options.maximum_work,
            bytes: self.bytes,
            limit: options.maximum_bytes,
            interrupt,
        };
        b.charge(1)?;
        b.retain(1024 * 1024)?;
        let location = self.classify_core(point, &mut b)?;
        Ok(PointQuery {
            location,
            work: b.work,
        })
    }
    fn classify_core(&self, point: [f64; 3], b: &mut Budget<'_>) -> Result<PointLocation, Error> {
        let unit = integer(self.polar.grid.denominator(), b)?;
        let heights = self
            .polar
            .planes
            .iter()
            .map(|(_, z)| physical(z, &unit, b))
            .collect::<Result<Vec<_>, _>>()?;
        let z = p(point[2]);
        let collar = p(self.sewing);
        if z.hi < (heights[0] - collar).lo || z.lo > (heights[3] + collar).hi {
            return Ok(PointLocation::Outside);
        }
        if heights
            .iter()
            .any(|&h| z.lo <= (h + collar).hi && z.hi >= (h - collar).lo)
        {
            return Ok(PointLocation::Unresolved);
        }
        let group = (0..3)
            .find(|&i| z.lo > heights[i].hi && z.hi < heights[i + 1].lo)
            .ok_or(Error::Unsupported)?;
        let distance = (p(point[0]) * p(point[0]) + p(point[1]) * p(point[1])).sqrt();
        let coarse = compare(distance, self.radii[group], collar);
        if coarse != PointLocation::Unresolved {
            return Ok(coarse);
        }
        let thread = if group == 0 {
            Some(0)
        } else if group == 2 {
            Some(1)
        } else {
            None
        };
        let lead = thread
            .map(|g| physical(&self.polar.leads[g], &unit, b))
            .transpose()?;
        let origin = thread
            .map(|g| physical(&self.polar.origins[g], &unit, b))
            .transpose()?;
        let plain = physical(&self.polar.plain_radius, &unit, b)?;
        let margin = collar / p(self.norm.lo);
        let mut result = None;
        let mut nodes = 0;
        for piece in &self.polar.angular {
            let mut pending = vec![Node {
                xy: [
                    piece.bounds.iter().map(|v| v[0]).collect(),
                    piece.bounds.iter().map(|v| v[1]).collect(),
                ],
                a: physical(piece.knots[0].numer(), piece.knots[0].denom(), b)?,
                z: physical(piece.knots[1].numer(), piece.knots[1].denom(), b)?,
                depth: 0,
            }];
            while let Some(node) = pending.pop() {
                b.charge(node.xy[0].len() * 8 + 1)?;
                nodes += 1;
                if nodes > 4096 {
                    return Err(Error::Limit);
                }
                let cross = node.xy[0]
                    .iter()
                    .zip(&node.xy[1])
                    .map(|(&x, &y)| p(point[0]) * y - p(point[1]) * x)
                    .collect::<Polynomial>();
                let dot = node.xy[0]
                    .iter()
                    .zip(&node.xy[1])
                    .map(|(&x, &y)| p(point[0]) * x + p(point[1]) * y)
                    .collect::<Polynomial>();
                let c = hull(&cross);
                if c.lo > 0. || c.hi < 0. || hull(&dot).hi <= 0. {
                    continue;
                }
                if node.depth < 40 {
                    let xs = poly::split(&node.xy[0], p(0.5), b)?;
                    let ys = poly::split(&node.xy[1], p(0.5), b)?;
                    let mid = (node.a + node.z) / p(2.);
                    pending.push(Node {
                        xy: [xs[1].clone(), ys[1].clone()],
                        a: mid,
                        z: node.z,
                        depth: node.depth + 1,
                    });
                    pending.push(Node {
                        xy: [xs[0].clone(), ys[0].clone()],
                        a: node.a,
                        z: mid,
                        depth: node.depth + 1,
                    });
                    continue;
                }
                let x = hull(&node.xy[0]);
                let y = hull(&node.xy[1]);
                let squared = x * x + y * y;
                if squared.lo <= 0. {
                    return Ok(PointLocation::Unresolved);
                }
                let lambda = (p(point[0]) * x + p(point[1]) * y) / squared;
                let mut local = None;
                if let (Some(g), Some(lead), Some(origin)) = (thread, lead, origin) {
                    let left = self
                        .polar
                        .strips
                        .iter()
                        .find(|s| s.lower == group)
                        .ok_or(Error::Invalid)?
                        .left;
                    let v = I::new(node.a.lo, node.z.hi);
                    let v = if left { p(1.) - v } else { v };
                    let raw = z - origin - lead * v;
                    let turns = raw / lead;
                    if !turns.valid() || turns.lo.abs() > 1e9 || turns.hi.abs() > 1e9 {
                        return Ok(PointLocation::Unresolved);
                    }
                    let lo = turns.lo.floor() as i64;
                    let hi = turns.hi.floor() as i64;
                    if hi - lo > 2 {
                        return Ok(PointLocation::Unresolved);
                    }
                    for turn in lo..=hi {
                        for strip in self
                            .polar
                            .strips
                            .iter()
                            .filter(|s| s.lower == if g == 0 { 0 } else { 2 })
                        {
                            b.charge(1)?;
                            let target = raw - lead * p(turn as f64) - strip.phase_start;
                            if target.hi < 0. || target.lo > strip.phase_width.hi {
                                continue;
                            }
                            let target =
                                I::new(target.lo.max(0.), target.hi.min(strip.phase_width.hi));
                            let radius = row_radius(strip, target, b)?;
                            let class = compare(lambda, radius, margin);
                            if local.is_some_and(|old| old != class) {
                                return Ok(PointLocation::Unresolved);
                            }
                            local = Some(class);
                        }
                    }
                } else {
                    local = Some(compare(lambda, plain, margin));
                }
                let Some(class) = local else {
                    return Ok(PointLocation::Unresolved);
                };
                if result.is_some_and(|old| old != class) {
                    return Ok(PointLocation::Unresolved);
                }
                result = Some(class);
            }
        }
        Ok(result.unwrap_or(PointLocation::Unresolved))
    }
}
