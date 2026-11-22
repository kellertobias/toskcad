//! Actual repaired source integral; no legacy Profile volume or angle proxy.
use super::*;
use crate::tapered_thread::Series;
const CELLS: usize = 65536;
#[derive(Clone, Copy)]
pub struct RepairedHeadVolumeOptions {
    pub maximum_width_mm3: f64,
    pub maximum_cells: usize,
    pub maximum_work: usize,
    pub maximum_bytes: usize,
}
impl Default for RepairedHeadVolumeOptions {
    fn default() -> Self {
        Self {
            maximum_width_mm3: 1e-6,
            maximum_cells: CELLS,
            maximum_work: WORK,
            maximum_bytes: BYTES,
        }
    }
}
/// Local-chart volume of this distinct embedded revolution source. Affine world
/// transport, other material components and raw-file matching remain separate.
pub struct RepairedHeadVolume<'h, 'p> {
    source: &'h RepairedHeadMeridian<'p>,
    volume: [f64; 2],
    cells: usize,
    work: usize,
}
impl<'h, 'p> RepairedHeadVolume<'h, 'p> {
    pub fn source(&self) -> &'h RepairedHeadMeridian<'p> {
        self.source
    }
    pub fn volume_mm3(&self) -> [f64; 2] {
        self.volume
    }
    pub fn cells(&self) -> usize {
        self.cells
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    /// Internal aggregate callers validate1<=cells<=CELLS before querying.
    pub(crate) fn scratch_budget_bytes(cells: usize) -> usize {
        4 * (cells + 1) * std::mem::size_of::<Cell>() + SCRATCH
    }
}
#[derive(Clone, Copy)]
struct Cell {
    part: usize,
    range: [f64; 2],
}
fn pi() -> I {
    I::new(
        std::f64::consts::PI.next_down(),
        std::f64::consts::PI.next_up(),
    )
}
fn density(part: ArcSource, t: I) -> Series<5> {
    let [r, z] = part.position_jet(Series::<6>::variable(t));
    let r = Series::<5> {
        c: std::array::from_fn(|i| r.c[i]),
    };
    // Integrand fourth derivative requires the FIFTH z derivative.
    let dz = Series::<5> {
        c: std::array::from_fn(|i| z.c[i + 1] * I::point((i + 1) as f64)),
    };
    r * r * dz * Series::interval(pi())
}
fn line(s: Segment) -> I {
    if s.start[1] == s.end[1] {
        return I::point(0.);
    }
    let [a, b] = [s.start[0], s.end[0]].map(I::point);
    pi() * (I::point(s.end[1]) - I::point(s.start[1])) * (a.square() + a * b + b.square())
        / I::point(3.)
}
fn enclosure(part: ArcSource, [a, b]: [f64; 2]) -> I {
    let width = I::point(b) - I::point(a);
    let mid = (I::point(a) + I::point(b)) / I::point(2.);
    let offset = width / (I::point(2.) * I::point(3.).sqrt());
    let quadrature = width * (density(part, mid - offset).c[0] + density(part, mid + offset).c[0])
        / I::point(2.);
    let fourth = I::point(density(part, I::new(a, b)).c[4].magnitude());
    // Gauss2 error h^5*f''''/4320 = h^5*Taylor_c4/180.
    let h2 = width.square();
    let error = (h2 * h2 * width * fourth / I::point(180.)).hi;
    quadrature + I::new(-error, error)
}
impl<'p> RepairedHeadMeridian<'p> {
    pub fn volume_bounds<'h>(
        &'h self,
        options: RepairedHeadVolumeOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<RepairedHeadVolume<'h, 'p>, RepairedHeadError> {
        poll(interrupt)?;
        if !options.maximum_width_mm3.is_finite()
            || options.maximum_width_mm3 <= 0.
            || options.maximum_cells == 0
            || options.maximum_cells > CELLS
            || options.maximum_work > WORK
            || options.maximum_bytes > BYTES
        {
            return Err(RepairedHeadError::Invalid);
        }
        let peak = std::mem::size_of::<RepairedHeadVolume<'_, '_>>()
            + self.retained_budget_bytes()
            + RepairedHeadVolume::scratch_budget_bytes(options.maximum_cells);
        if peak > options.maximum_bytes {
            return Err(RepairedHeadError::Limit);
        }
        let mut remaining = options
            .maximum_work
            .checked_sub(self.work)
            .ok_or(RepairedHeadError::Limit)?;
        let mut volume = I::point(0.);
        let mut pending = Vec::new();
        // At most five immutable arc descriptors fit the fixed scratch reserve.
        // Only arcs enter numerical integration; lines were integrated exactly.
        let mut arcs = Vec::with_capacity(self.count);
        let mut cells = 0usize;
        for i in 0..self.count {
            poll(interrupt)?;
            match self.parts[i].unwrap() {
                Part::Line(s) => {
                    if cells == options.maximum_cells {
                        return Err(RepairedHeadError::Limit);
                    }
                    remaining = remaining.checked_sub(128).ok_or(RepairedHeadError::Limit)?;
                    cells += 1;
                    volume = volume + line(s);
                }
                Part::Arc(arc) => {
                    arcs.push(arc);
                    pending.push(Cell {
                        part: arcs.len() - 1,
                        range: [0., 1.],
                    });
                }
            }
        }
        // These positive-source contributions cannot lose their current width
        // under subsequent outward addition. Refuse an impossible line floor.
        poll(interrupt)?;
        if !finite(volume)
            || (I::point(volume.hi) - I::point(volume.lo)).hi > options.maximum_width_mm3
        {
            return Err(RepairedHeadError::Unresolved);
        }
        let available =
            I::point(options.maximum_width_mm3) - (I::point(volume.hi) - I::point(volume.lo));
        let quota = (available / I::point(4. * self.count as f64)).lo;
        while let Some(cell) = pending.pop() {
            poll(interrupt)?;
            if cells == options.maximum_cells {
                return Err(RepairedHeadError::Limit);
            }
            // Three Series6/5 evaluations, fixed order<=5 and at most two
            // normalized coordinates; bound all coefficient convolutions.
            remaining = remaining
                .checked_sub(4096)
                .ok_or(RepairedHeadError::Limit)?;
            cells += 1;
            let value = enclosure(arcs[cell.part], cell.range);
            poll(interrupt)?;
            ensure(finite(value))?;
            let width = (I::point(value.hi) - I::point(value.lo)).hi;
            let target = (I::point(quota) * (I::point(cell.range[1]) - I::point(cell.range[0]))).lo;
            if width <= target && value.lo >= 0. {
                volume = volume + value;
            } else {
                let [a, b] = cell.range;
                let mid = a + (b - a) / 2.;
                if mid <= a || mid >= b {
                    return Err(RepairedHeadError::Unresolved);
                }
                pending.push(Cell {
                    part: cell.part,
                    range: [mid, b],
                });
                pending.push(Cell {
                    part: cell.part,
                    range: [a, mid],
                });
            }
        }
        let finish = || {
            ensure(
                finite(volume)
                    && volume.lo > 0.
                    && (I::point(volume.hi) - I::point(volume.lo)).hi <= options.maximum_width_mm3,
            )?;
            Ok([volume.lo, volume.hi])
        };
        let result = finish();
        poll(interrupt)?;
        let remaining = remaining.checked_sub(1).ok_or(RepairedHeadError::Limit)?;
        Ok(RepairedHeadVolume {
            source: self,
            volume: result?,
            cells,
            work: options.maximum_work - remaining,
        })
    }
}
#[cfg(test)]
mod tests;
