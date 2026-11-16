//! Actual retained face integrals. Planar head measures use the earned polygon
//! extrusion and the physical annulus' identical projected XY chart. Its
//! contact/plane sewing collar remains separate, as in the owned reference.
use super::*;
use crate::tapered_thread::interval::I;
use std::mem::size_of;
mod curved;
pub use curved::{HeadGeometryOptions, HeadGeometryPlan};
const MEASURE_WORK: usize = 512_000_000;
const MEASURE_SCRATCH: usize = 16 * 1024;
#[derive(Clone, Copy, Default)]
pub struct HeadMeasureOptions {
    pub measure: FaceMeasureOptions,
    pub additional_resident_bytes: usize,
}
pub struct HeadFaceMeasure<'source> {
    source: &'source OwnedActualHead,
    ordinal: usize,
    entity: u32,
    area: [f64; 2],
    centroid: [[f64; 2]; 3],
    work: usize,
    inherited_plan_work: usize,
    peak: usize,
}
impl<'source> HeadFaceMeasure<'source> {
    pub fn source(&self) -> &'source OwnedActualHead {
        self.source
    }
    /// SAME physical ordinal convention as the actual-head mesh: nine
    /// original thread faces, then eight actual head faces; no derived cap.
    pub fn physical_face_ordinal(&self) -> usize {
        self.ordinal
    }
    pub fn physical_face_entity(&self) -> u32 {
        self.entity
    }
    pub fn area_mm2(&self) -> [f64; 2] {
        self.area
    }
    pub fn centroid_mm(&self) -> [[f64; 2]; 3] {
        self.centroid
    }
    /// Fresh integration work, not immutable material construction.
    pub fn work_used(&self) -> usize {
        self.work
    }
    /// Preparation of the immutable geometry plan, counted once in the
    /// composed property allowance; zero for the direct planar helper.
    pub fn inherited_plan_work_used(&self) -> usize {
        self.inherited_plan_work
    }
    pub fn aggregate_work_used(&self) -> usize {
        self.inherited_plan_work + self.work
    }
    pub fn accounted_peak_bytes(&self) -> usize {
        self.peak
    }
}
fn checked(options: HeadMeasureOptions) -> Result<()> {
    let o = options.measure;
    if !o.area_width_mm2.is_finite()
        || o.area_width_mm2 <= 0.
        || !o.centroid_width_mm.is_finite()
        || o.centroid_width_mm <= 0.
        || o.maximum_work == 0
        || o.maximum_work > MEASURE_WORK
        || o.maximum_bytes == 0
        || o.maximum_bytes > MAX_BYTES
    {
        return Err(Error::Invalid);
    }
    Ok(())
}
fn polygon(source: &OwnedActualHead, b: &mut Budget<'_>) -> Result<(I, [I; 2])> {
    let p = source.head_polygon_mm();
    let mut twice = I::point(0.);
    let mut moment = [I::point(0.); 2];
    for i in 0..6 {
        b.step(128)?;
        let a = p[i].map(I::point);
        let z = p[(i + 1) % 6].map(I::point);
        let cross = a[0] * z[1] - a[1] * z[0];
        twice = twice + cross;
        for k in 0..2 {
            moment[k] = moment[k] + (a[k] + z[k]) * cross / I::point(6.);
        }
    }
    let area = twice / I::point(2.);
    if !area.valid() || area.lo <= 0. || moment.iter().any(|v| !v.valid()) {
        return Err(Error::Unsupported);
    }
    Ok((area, moment))
}
impl OwnedActualHead {
    /// Select one actual planar physical face without preparing or integrating
    /// any curved descriptor. Curved selections return Unsupported.
    pub fn planar_face_measure(
        &self,
        ordinal: usize,
        options: HeadMeasureOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<HeadFaceMeasure<'_>> {
        checked(options)?;
        if ordinal >= 17 {
            return Err(Error::Invalid);
        }
        let extra = options
            .additional_resident_bytes
            .checked_add(size_of::<HeadFaceMeasure<'_>>())
            .and_then(|n| n.checked_add(MEASURE_SCRATCH))
            .ok_or(Error::Limit)?;
        let bytes = self
            .receipt()
            .retained_bytes
            .checked_add(extra)
            .ok_or(Error::Limit)?;
        if bytes > options.measure.maximum_bytes {
            return Err(Error::Limit);
        }
        let mut b = Budget {
            used: 0,
            remaining: options.measure.maximum_work,
            bytes,
            peak: bytes,
            limit: options.measure.maximum_bytes,
            interrupt,
        };
        b.step(256)?;
        let roles = self.physical_roles();
        let entity = if ordinal < 9 {
            roles.thread_face_entities[ordinal]
        } else {
            roles.head_face_entities[ordinal - 9]
        };
        let cap = |index: usize, b: &mut Budget<'_>| -> Result<FaceMeasure> {
            if b.remaining == 0 {
                return Err(Error::Limit);
            }
            if options.measure.area_width_mm2 / 4. == 0.
                || options.measure.centroid_width_mm / 4. == 0.
            {
                return Err(Error::Unsupported);
            }
            let outside = bytes
                .checked_sub(self.auxiliary_thread().retained_budget_bytes())
                .ok_or(Error::Limit)?;
            let value = self.auxiliary_thread().cap_face_measure(
                index,
                FaceMeasureOptions {
                    maximum_work: b.remaining,
                    maximum_bytes: options
                        .measure
                        .maximum_bytes
                        .checked_sub(outside)
                        .ok_or(Error::Limit)?,
                    area_width_mm2: options.measure.area_width_mm2 / 4.,
                    centroid_width_mm: options.measure.centroid_width_mm / 4.,
                },
                interrupt,
            )?;
            b.step(value.work_used())?;
            b.peak = b.peak.max(
                value
                    .accounted_peak_bytes()
                    .checked_add(outside)
                    .ok_or(Error::Limit)?,
            );
            Ok(value)
        };
        let (area, centroid) = if ordinal < 9 {
            let index = self
                .auxiliary_thread()
                .faces()
                .iter()
                .position(|f| f.entity() == entity)
                .ok_or(Error::Invalid)?;
            if !matches!(
                self.auxiliary_thread().faces()[index].support(),
                Support::Plane { .. }
            ) {
                return Err(Error::Unsupported);
            }
            let value = cap(index, &mut b)?;
            b.step(128)?;
            let centroid = self
                .axial_sewing_transfer()
                .map_raw_box_mm(value.centroid_mm())?
                .map(|v| I::new(v[0], v[1]));
            let area = value.area_mm2();
            (I::new(area[0], area[1]), centroid)
        } else {
            let face = &self.head().faces[ordinal - 9];
            let [low, high] = self.head_stations_mm();
            if face.outer.len() == 6 {
                let (outer, mut moments) = polygon(self, &mut b)?;
                let (area, z) = if !face.inner.is_empty() {
                    let index = self
                        .auxiliary_thread()
                        .faces()
                        .iter()
                        .position(|f| f.entity() == self.head().auxiliary_cap_entity)
                        .ok_or(Error::Invalid)?;
                    let value = cap(index, &mut b)?;
                    b.step(256)?;
                    let a = value.area_mm2();
                    let a = I::new(a[0], a[1]);
                    let anchor = value.moment_anchor_mm();
                    let m = value.first_moments_about_anchor_mm3();
                    for k in 0..2 {
                        moments[k] =
                            moments[k] - (I::new(m[k][0], m[k][1]) + a * I::point(anchor[k]));
                    }
                    (outer - a, low)
                } else {
                    (outer, high)
                };
                if !area.valid() || area.lo <= 0. {
                    return Err(Error::Unsupported);
                }
                (area, [moments[0] / area, moments[1] / area, I::point(z)])
            } else {
                b.step(256)?;
                // Earned prism topology gives exactly two distinct XY
                // stations on each side, repeated at low/high Z.
                let edge = self.head().edges[face.outer[0].0];
                let a = self.head().vertices[edge[0]];
                let other = face
                    .outer
                    .iter()
                    .flat_map(|&(i, _)| self.head().edges[i])
                    .map(|i| self.head().vertices[i])
                    .find(|p| p[..2] != a[..2])
                    .ok_or(Error::Unsupported)?;
                let dx = I::point(other[0]) - I::point(a[0]);
                let dy = I::point(other[1]) - I::point(a[1]);
                let area = (dx.square() + dy.square()).sqrt() * (I::point(high) - I::point(low));
                (
                    area,
                    [
                        (I::point(a[0]) + I::point(other[0])) / I::point(2.),
                        (I::point(a[1]) + I::point(other[1])) / I::point(2.),
                        (I::point(low) + I::point(high)) / I::point(2.),
                    ],
                )
            }
        };
        b.step(128)?;
        if !area.valid()
            || area.lo <= 0.
            || area.width() > options.measure.area_width_mm2
            || centroid
                .iter()
                .any(|v| !v.valid() || v.width() > options.measure.centroid_width_mm)
        {
            return Err(Error::Unsupported);
        }
        b.step(0)?;
        Ok(HeadFaceMeasure {
            source: self,
            ordinal,
            entity,
            area: [area.lo, area.hi],
            centroid: centroid.map(|v| [v.lo, v.hi]),
            work: b.used,
            inherited_plan_work: 0,
            peak: b.peak,
        })
    }
}
#[cfg(test)]
mod tests;
