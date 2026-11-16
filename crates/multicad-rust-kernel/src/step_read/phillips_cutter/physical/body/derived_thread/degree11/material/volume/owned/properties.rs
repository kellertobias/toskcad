//! Actual D11 planar end-disk measures. No nominal circle or mesh supplies them.
use super::*;
mod planar;
const END_WORK: usize = 512_000_000;
#[derive(Clone, Copy, Debug)]
pub struct Degree11EndDiskOptions {
    pub maximum_area_width_mm2: f64,
    pub maximum_first_moment_width_mm3: f64,
    pub maximum_centroid_width_mm: f64,
    /// Fresh measurement work; immutable material construction is not repeated.
    pub maximum_work: usize,
    pub maximum_bytes: usize,
}
impl Default for Degree11EndDiskOptions {
    fn default() -> Self {
        Self {
            maximum_area_width_mm2: 1e-7,
            maximum_first_moment_width_mm3: 1e-7,
            maximum_centroid_width_mm: 1e-7,
            maximum_work: END_WORK,
            maximum_bytes: MAX_BYTES,
        }
    }
}
#[derive(Clone, Copy, Debug)]
pub struct Degree11EndDiskProperties {
    pub station_mm: f64,
    pub area_mm2: [f64; 2],
    pub first_moment_mm3: [[f64; 2]; 3],
    pub centroid_mm: [[f64; 2]; 3],
}
pub struct Degree11EndDisks<'s> {
    source: &'s OwnedDegree11Material,
    ends: [Degree11EndDiskProperties; 2],
    work: usize,
    peak: usize,
}
impl Degree11EndDisks<'_> {
    pub fn source(&self) -> &OwnedDegree11Material {
        self.source
    }
    /// Actual entry disk followed by actual tip disk. The entry disk is a
    /// numerical operand for the joined head annulus, not an exposed face.
    pub fn ends(&self) -> &[Degree11EndDiskProperties; 2] {
        &self.ends
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn accounted_peak_bytes(&self) -> usize {
        self.peak
    }
}
fn width(v: I) -> f64 {
    (I::point(v.hi) - I::point(v.lo)).hi
}
fn pair(v: I) -> [f64; 2] {
    [v.lo, v.hi]
}
struct EndWork<'a> {
    remaining: usize,
    interrupt: Option<&'a dyn Fn() -> bool>,
}
impl EndWork<'_> {
    fn step(&mut self, n: usize) -> std::result::Result<(), PhillipsReadError> {
        if self.interrupt.is_some_and(|f| f()) {
            return Err(PhillipsReadError::Cancelled);
        }
        self.remaining = self
            .remaining
            .checked_sub(n)
            .ok_or(PhillipsReadError::Limit)?;
        Ok(())
    }
}
impl OwnedDegree11Material {
    pub fn end_disk_properties(
        &self,
        o: Degree11EndDiskOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> std::result::Result<Degree11EndDisks<'_>, PhillipsReadError> {
        let valid = [
            o.maximum_area_width_mm2,
            o.maximum_first_moment_width_mm3,
            o.maximum_centroid_width_mm,
        ]
        .into_iter()
        .all(|v| v.is_finite() && v > 0.);
        if !valid || o.maximum_work > END_WORK || o.maximum_bytes > MAX_BYTES {
            return Err(PhillipsReadError::Invalid);
        }
        let mut work = EndWork {
            remaining: o.maximum_work,
            interrupt,
        };
        work.step(planar::SETUP_WORK)?;
        let peak = self
            .retained_budget_bytes()
            .checked_add(planar::SCRATCH_BYTES)
            .and_then(|n| n.checked_add(std::mem::size_of::<Degree11EndDisks<'_>>()))
            .ok_or(PhillipsReadError::Limit)?;
        if peak > o.maximum_bytes {
            return Err(PhillipsReadError::Limit);
        }
        let weights = planar::Weights::new();
        let mut totals = [[I::point(0.); 3]; 2];
        for face in self.thread_faces() {
            let tensor = face.source();
            let [_, height] = tensor.dimensions();
            let sign = I::point(if face.forward() { 1. } else { -1. });
            for (end, total) in totals.iter_mut().enumerate() {
                let row = if end == 0 { 0 } else { height - 1 };
                for span in 0..tensor.spans()[0] {
                    work.step(planar::CURVE_WORK)?;
                    let p = std::array::from_fn(|i| {
                        let p = tensor.poles()[(11 * span + i) * height + row];
                        [I::point(p[0]), I::point(p[1])]
                    });
                    let value = weights.curve(p);
                    *total = std::array::from_fn(|i| total[i] + sign * value[i]);
                }
            }
        }
        let mut ends = [Degree11EndDiskProperties {
            station_mm: 0.,
            area_mm2: [0.; 2],
            first_moment_mm3: [[0.; 2]; 3],
            centroid_mm: [[0.; 2]; 3],
        }; 2];
        for (end, v) in totals.into_iter().enumerate() {
            let station = self.cap_stations_mm()[end];
            let moments = [v[1], v[2], I::point(station) * v[0]];
            let centroid = [v[1] / v[0], v[2] / v[0], I::point(station)];
            if !v[0].valid()
                || v[0].lo <= 0.
                || width(v[0]) > o.maximum_area_width_mm2
                || moments
                    .iter()
                    .any(|&x| !x.valid() || width(x) > o.maximum_first_moment_width_mm3)
                || centroid
                    .iter()
                    .any(|&x| !x.valid() || width(x) > o.maximum_centroid_width_mm)
            {
                return Err(PhillipsReadError::Unsupported);
            }
            ends[end] = Degree11EndDiskProperties {
                station_mm: station,
                area_mm2: pair(v[0]),
                first_moment_mm3: moments.map(pair),
                centroid_mm: centroid.map(pair),
            };
        }
        work.step(1)?;
        Ok(Degree11EndDisks {
            source: self,
            ends,
            work: o.maximum_work - work.remaining,
            peak,
        })
    }
}
#[cfg(test)]
mod tests;
