//! Chords of an original source edge; no manufacturing surface is substituted.
use super::*;
#[derive(Debug)]
pub struct StockClipPolyline {
    curve: StockClipCurve,
    points: Vec<[f64; 3]>,
    error: f64,
    work: usize,
    bytes: usize,
}
impl StockClipPolyline {
    pub fn curve(&self) -> StockClipCurve {
        self.curve
    }
    pub fn points_mm(&self) -> &[[f64; 3]] {
        &self.points
    }
    pub fn source_chord_error_mm(&self) -> f64 {
        self.error
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_budget_bytes(&self) -> usize {
        self.bytes
    }
}
#[derive(Clone, Copy)]
struct Point {
    p: [f64; 3],
    error: f64,
}
#[derive(Clone, Copy)]
struct Edge {
    t: [f64; 2],
    ends: [Point; 2],
    depth: u32,
}
fn chosen(p: StockClipPoint) -> Point {
    let xyz = p.position_mm().map(|v| I::new(v[0], v[1]));
    let literal = xyz.map(I::midpoint);
    let error = (0..3)
        .fold(I::point(0.), |sum, i| {
            sum + I::point((xyz[i] - I::point(literal[i])).magnitude())
        })
        .hi;
    Point { p: literal, error }
}
impl StockClipPlan {
    pub fn curve_polyline(
        &self,
        curve: StockClipCurve,
        chord_error_mm: f64,
        point_limit: usize,
        options: StockClipOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<StockClipPolyline, Error> {
        if !chord_error_mm.is_finite()
            || chord_error_mm <= 0.
            || !(2..=131072).contains(&point_limit)
        {
            return Err(Error::Invalid);
        }
        let fixed = size_of::<Self>() + size_of::<StockClipPolyline>() + 2 * SCRATCH;
        let mut work = Meter::new(options, fixed - SCRATCH, interrupt)?;
        work.charge(self.work)?;
        let (patch, axis, vfixed) = match curve {
            StockClipCurve::Axial { slab, rail } => (
                StockClipPatch {
                    slab,
                    band: BANDS[rail_index(rail)],
                },
                1,
                0.,
            ),
            StockClipCurve::Ring { station, band } => {
                let (slab, v) = match station {
                    StockClipStation::Start => (StockClipSlab::Before, 0.),
                    StockClipStation::Contact => (StockClipSlab::Before, 1.),
                    StockClipStation::End => (StockClipSlab::After, 1.),
                };
                (StockClipPatch { slab, band }, 0, v)
            }
        };
        let law = self.patch_law(patch, work.child_options(options)?, interrupt)?;
        work.charge(law.work_used())?;
        let mut ends = [Point {
            p: [0.; 3],
            error: 0.,
        }; 2];
        for (i, key) in curve.endpoints().into_iter().enumerate() {
            let p = self.vertex_bounds(key, work.child_options(options)?, interrupt)?;
            work.charge(p.work_used())?;
            ends[i] = chosen(p);
        }
        let blank = Edge {
            t: [0., 1.],
            ends,
            depth: 0,
        };
        let mut pending = [blank; 64];
        let mut len = 1;
        let mut out = Vec::new();
        let mut max_error: f64 = 0.;
        let mut peak_bytes = fixed;
        reserve(
            &mut out,
            point_limit,
            fixed,
            options,
            &mut work,
            &mut peak_bytes,
        )?;
        out.push(ends[0].p);
        while len > 0 {
            len -= 1;
            let edge = pending[len];
            let (u, v) = if axis == 0 {
                (edge.t, [vfixed, vfixed])
            } else {
                ([0., 0.], edge.t)
            };
            let jets = law.jets(u, v, work.child_options(options)?, interrupt)?;
            work.charge(jets.work_used())?;
            let second = if axis == 0 { jets.duu() } else { jets.dvv() };
            let curvature = second.into_iter().fold(I::point(0.), |a, v| {
                a + I::point(I::new(v[0], v[1]).magnitude())
            });
            let h = I::point(edge.t[1]) - I::point(edge.t[0]);
            let error = (curvature * h.square() / I::point(8.)
                + I::point(edge.ends[0].error.max(edge.ends[1].error)))
            .hi;
            if error <= chord_error_mm {
                reserve(
                    &mut out,
                    point_limit,
                    fixed,
                    options,
                    &mut work,
                    &mut peak_bytes,
                )?;
                out.push(edge.ends[1].p);
                max_error = max_error.max(error);
                continue;
            }
            if edge.depth >= 48 {
                return Err(Error::Unresolved);
            }
            let mid = edge.t[0] + (edge.t[1] - edge.t[0]) * 0.5;
            // Source t starts at [0,1] and has fewer than48 dyadic
            // bisections here: midpoint representation is exact and strict.
            let p =
                self.curve_bounds(curve, [mid, mid], work.child_options(options)?, interrupt)?;
            work.charge(p.work_used())?;
            let point = chosen(p);
            // Reverse stack order emits increasing exact dyadic source t.
            pending[len] = Edge {
                t: [mid, edge.t[1]],
                ends: [point, edge.ends[1]],
                depth: edge.depth + 1,
            };
            pending[len + 1] = Edge {
                t: [edge.t[0], mid],
                ends: [edge.ends[0], point],
                depth: edge.depth + 1,
            };
            len += 2;
        }
        work.charge(1)?;
        let bytes = peak_bytes.max(fixed + out.capacity() * size_of::<[f64; 3]>());
        Ok(StockClipPolyline {
            curve,
            points: out,
            error: max_error,
            work: work.used,
            bytes,
        })
    }
}
fn reserve(
    out: &mut Vec<[f64; 3]>,
    limit: usize,
    fixed: usize,
    options: StockClipOptions,
    work: &mut Meter<'_>,
    peak_bytes: &mut usize,
) -> Result<(), Error> {
    work.charge(32)?;
    if out.len() == limit {
        return Err(Error::Limit);
    }
    if out.len() < out.capacity() {
        return Ok(());
    }
    let next = (out.capacity().max(16) * 2).min(limit);
    let peak = fixed.saturating_add(
        out.capacity()
            .saturating_add(next)
            .saturating_mul(size_of::<[f64; 3]>()),
    );
    *peak_bytes = (*peak_bytes).max(peak);
    if peak > options.maximum_bytes {
        return Err(Error::Limit);
    }
    out.try_reserve_exact(next - out.len())
        .map_err(|_| Error::Limit)?;
    // Pinned non-ZST RawVec exact-growth invariant (rustc58cf8f9c35b3):
    // requested logical capacity is stored verbatim. Allocation/RSS slack
    // is outside this retained-budget contract; observe actual capacity.
    debug_assert_eq!(out.capacity(), next);
    Ok(())
}
#[cfg(test)]
mod tests;
