//! Actual repaired meridian surface integrals, independent of mesh triangles.
use super::*;
use crate::tapered_thread::{Scalar, Series};
const MAX_CELLS: usize = 65_536;
const DEPTH: usize = 48;
const CELL_WORK: usize = 8192;
const SURFACE_SCRATCH: usize = 32 * 1024;
#[derive(Clone, Copy, Debug)]
pub struct RepairedHeadSurfaceOptions {
    pub area_tolerance_mm2: f64,
    pub first_moment_tolerance_mm3: f64,
    pub centroid_tolerance_mm: f64,
    pub maximum_cells: usize,
    pub maximum_work: usize,
    pub maximum_bytes: usize,
}
impl Default for RepairedHeadSurfaceOptions {
    fn default() -> Self {
        Self {
            area_tolerance_mm2: 1e-7,
            first_moment_tolerance_mm3: 1e-7,
            centroid_tolerance_mm: 1e-7,
            maximum_cells: MAX_CELLS,
            maximum_work: WORK,
            maximum_bytes: BYTES,
        }
    }
}
#[derive(Clone, Copy, Debug)]
pub struct RepairedHeadSurface {
    pub segment: usize,
    pub area_mm2: [f64; 2],
    pub first_moment_mm3: [[f64; 2]; 3],
    pub centroid_mm: [[f64; 2]; 3],
}
pub struct RepairedHeadSurfaces<'h, 'p> {
    source: &'h RepairedHeadMeridian<'p>,
    surfaces: Vec<RepairedHeadSurface>,
    cells: usize,
    work: usize,
    peak: usize,
}
impl<'h, 'p> RepairedHeadSurfaces<'h, 'p> {
    pub fn source(&self) -> &'h RepairedHeadMeridian<'p> {
        self.source
    }
    pub fn surfaces(&self) -> &[RepairedHeadSurface] {
        &self.surfaces
    }
    pub fn cells(&self) -> usize {
        self.cells
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn accounted_peak_bytes(&self) -> usize {
        self.peak
    }
    pub fn into_surfaces(self) -> Vec<RepairedHeadSurface> {
        self.surfaces
    }
}
#[derive(Clone, Copy)]
struct Cell {
    range: [f64; 2],
    depth: usize,
}
fn pair(x: I) -> [f64; 2] {
    [x.lo, x.hi]
}
fn pi() -> I {
    I::new(
        std::f64::consts::PI.next_down(),
        std::f64::consts::PI.next_up(),
    )
}
fn line(s: Segment) -> [I; 2] {
    let [r, z] = s.start.map(I::point);
    let [dr, dz] = std::array::from_fn::<_, 2, _>(|i| I::point(s.end[i]) - I::point(s.start[i]));
    let scale = I::point(2.) * pi() * (dr.square() + dz.square()).sqrt();
    [
        scale * (r + dr / I::point(2.)),
        scale * (r * z + (r * dz + z * dr) / I::point(2.) + dr * dz / I::point(3.)),
    ]
}
fn density(arc: ArcSource, t: I) -> [Series<5>; 2] {
    let [r, z] = arc.position_jet(Series::<6>::variable(t));
    let dr = Series::<5> {
        c: std::array::from_fn(|i| r.c[i + 1] * I::point((i + 1) as f64)),
    };
    let dz = Series::<5> {
        c: std::array::from_fn(|i| z.c[i + 1] * I::point((i + 1) as f64)),
    };
    let r = Series::<5> {
        c: std::array::from_fn(|i| r.c[i]),
    };
    let z = Series::<5> {
        c: std::array::from_fn(|i| z.c[i]),
    };
    let area = Series::interval(I::point(2.) * pi()) * r * (dr * dr + dz * dz).sqrt();
    [area, z * area]
}
fn enclosure(arc: ArcSource, [a, b]: [f64; 2]) -> [I; 2] {
    let h = I::point(b) - I::point(a);
    let mid = (I::point(a) + I::point(b)) / I::point(2.);
    let d = h / (I::point(2.) * I::point(3.).sqrt());
    let left = density(arc, mid - d);
    let right = density(arc, mid + d);
    let all = density(arc, I::new(a, b));
    std::array::from_fn(|i| {
        let q = h * (left[i].c[0] + right[i].c[0]) / I::point(2.);
        // Gauss2: h^5 sup|f''''|/4320 = h^5 sup|Taylor_c4|/180.
        let e = (h.square().square() * h * I::point(all[i].c[4].magnitude()) / I::point(180.)).hi;
        q + I::new(-e, e)
    })
}
fn width(x: I) -> f64 {
    (I::point(x.hi) - I::point(x.lo)).hi
}
impl<'p> RepairedHeadMeridian<'p> {
    /// Each meridian segment is revolved completely. Cap/annulus subtraction
    /// and affine transport belong to the caller's actual material geometry.
    pub fn surface_properties<'h>(
        &'h self,
        o: RepairedHeadSurfaceOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<RepairedHeadSurfaces<'h, 'p>, RepairedHeadError> {
        poll(interrupt)?;
        if [
            o.area_tolerance_mm2,
            o.first_moment_tolerance_mm3,
            o.centroid_tolerance_mm,
        ]
        .into_iter()
        .any(|x| !x.is_finite() || x <= 0.)
            || o.maximum_cells == 0
            || o.maximum_cells > MAX_CELLS
            || o.maximum_work > WORK
            || o.maximum_bytes > BYTES
        {
            return Err(RepairedHeadError::Invalid);
        }
        let mut remaining = o
            .maximum_work
            .checked_sub(512)
            .ok_or(RepairedHeadError::Limit)?;
        let mut peak = self
            .retained_budget_bytes()
            .checked_add(SURFACE_SCRATCH)
            .and_then(|x| {
                x.checked_add(2 * self.count * std::mem::size_of::<RepairedHeadSurface>())
            })
            .ok_or(RepairedHeadError::Limit)?;
        if peak > o.maximum_bytes {
            return Err(RepairedHeadError::Limit);
        }
        let mut surfaces = Vec::new();
        surfaces
            .try_reserve_exact(self.count)
            .map_err(|_| RepairedHeadError::Limit)?;
        peak = peak.max(
            self.retained_budget_bytes()
                + SURFACE_SCRATCH
                + surfaces.capacity() * std::mem::size_of::<RepairedHeadSurface>(),
        );
        if peak > o.maximum_bytes {
            return Err(RepairedHeadError::Limit);
        }
        let mut cells = 0usize;
        for segment in 0..self.count {
            let certificate = self.certificates[segment].expect("sealed segment certificate");
            // Certified strict dz/dt (or radial bearing), integrated on [0,1].
            // A ZeroAtStart certificate bounds dz/dt divided by t.
            let progress = I::point(certificate.minimum_progress_mm)
                / I::point(if certificate.progress == AxialProgress::ZeroAtStart {
                    2.
                } else {
                    1.
                });
            let lower =
                (I::point(2.) * pi() * I::point(certificate.minimum_radius_mm) * progress).lo;
            ensure(lower > 0. && lower.is_finite())?;
            let zmag = self.station(segment).unwrap()[1]
                .abs()
                .max(self.station(segment + 1).unwrap()[1].abs());
            let a_target = o.area_tolerance_mm2.min(
                (I::point(o.centroid_tolerance_mm) * I::point(lower)
                    / (I::point(8.) * (I::point(zmag) + I::point(1.))))
                .lo,
            );
            let m_target = o
                .first_moment_tolerance_mm3
                .min((I::point(o.centroid_tolerance_mm) * I::point(lower) / I::point(8.)).lo);
            ensure(a_target > 0. && m_target > 0.)?;
            let part = self.parts[segment].expect("sealed segment source");
            let mut total = [I::point(0.); 2];
            let mut pending = [None; DEPTH + 2];
            pending[0] = Some(Cell {
                range: [0., 1.],
                depth: 0,
            });
            let mut count = 1;
            while count > 0 {
                poll(interrupt)?;
                if cells == o.maximum_cells {
                    return Err(RepairedHeadError::Limit);
                }
                remaining = remaining
                    .checked_sub(CELL_WORK)
                    .ok_or(RepairedHeadError::Limit)?;
                cells += 1;
                count -= 1;
                let cell = pending[count]
                    .take()
                    .expect("bounded pending segment range");
                let values = match part {
                    Part::Line(s) => line(s),
                    Part::Arc(a) => enclosure(a, cell.range),
                };
                let h = (I::point(cell.range[1]) - I::point(cell.range[0])).lo;
                if values.iter().all(|&x| finite(x))
                    && values[0].lo > 0.
                    && width(values[0]) <= (I::point(a_target) * I::point(h) / I::point(4.)).lo
                    && width(values[1]) <= (I::point(m_target) * I::point(h) / I::point(4.)).lo
                {
                    total = std::array::from_fn(|i| total[i] + values[i]);
                } else {
                    // Exact line expressions have no refinable parameter error.
                    if matches!(part, Part::Line(_)) || cell.depth == DEPTH {
                        return Err(RepairedHeadError::Unresolved);
                    }
                    let [a, b] = cell.range;
                    let mid = a + (b - a) / 2.;
                    ensure(mid > a && mid < b)?;
                    pending[count] = Some(Cell {
                        range: [mid, b],
                        depth: cell.depth + 1,
                    });
                    pending[count + 1] = Some(Cell {
                        range: [a, mid],
                        depth: cell.depth + 1,
                    });
                    count += 2;
                }
            }
            let center = total[1] / total[0];
            ensure(
                total.iter().all(|&x| finite(x))
                    && total[0].lo > 0.
                    && finite(center)
                    && width(total[0]) <= o.area_tolerance_mm2
                    && width(total[1]) <= o.first_moment_tolerance_mm3
                    && width(center) <= o.centroid_tolerance_mm,
            )?;
            surfaces.push(RepairedHeadSurface {
                segment,
                area_mm2: pair(total[0]),
                first_moment_mm3: [pair(I::point(0.)), pair(I::point(0.)), pair(total[1])],
                centroid_mm: [pair(I::point(0.)), pair(I::point(0.)), pair(center)],
            });
        }
        poll(interrupt)?;
        remaining = remaining.checked_sub(1).ok_or(RepairedHeadError::Limit)?;
        Ok(RepairedHeadSurfaces {
            source: self,
            surfaces,
            cells,
            work: o.maximum_work - remaining,
            peak,
        })
    }
}

#[cfg(test)]
mod tests;
