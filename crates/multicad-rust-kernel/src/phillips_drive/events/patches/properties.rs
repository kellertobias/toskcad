//! Actual source surface integrals with analytic integration in the second coordinate.
use super::*;
use crate::tapered_thread::{Scalar, Series};
use std::collections::BTreeMap;
mod density;
mod integrate;
#[cfg(test)]
mod tests;
type Values = [I; 4];
const MAXIMUM_CELLS: usize = 65536;
#[derive(Clone, Copy, Debug)]
pub struct FacePropertyOptions {
    pub area_tolerance_mm2: f64,
    pub first_moment_tolerance_mm3: f64,
    pub centroid_tolerance_mm: f64,
    /// Cumulative cell proofs, including unsuccessful derivative enclosures.
    pub maximum_cells: usize,
}
impl Default for FacePropertyOptions {
    fn default() -> Self {
        Self {
            area_tolerance_mm2: 1e-7,
            first_moment_tolerance_mm3: 1e-7,
            centroid_tolerance_mm: 1e-7,
            maximum_cells: 16384,
        }
    }
}
#[derive(Clone, Copy, Debug)]
pub struct FaceProperties {
    pub face: FaceKey,
    pub area_mm2: Bounds,
    pub first_moment_mm3: [Bounds; 3],
    pub centroid_mm: [Bounds; 3],
}
struct Work<'a> {
    remaining: usize,
    cancelled: &'a dyn Fn() -> bool,
}
impl Work<'_> {
    fn check(&self) -> Result<(), Error> {
        if (self.cancelled)() {
            Err(Error::Cancelled)
        } else {
            Ok(())
        }
    }
    fn step(&mut self) -> Result<(), Error> {
        self.check()?;
        if self.remaining == 0 {
            return Err(Error::Limit);
        }
        self.remaining -= 1;
        Ok(())
    }
}
impl FacePropertyOptions {
    fn validate(self) -> Result<(), Error> {
        if [
            self.area_tolerance_mm2,
            self.first_moment_tolerance_mm3,
            self.centroid_tolerance_mm,
        ]
        .into_iter()
        .any(|v| !v.is_finite() || v <= 0.)
            || self.maximum_cells == 0
            || self.maximum_cells > MAXIMUM_CELLS
        {
            return Err(Error::Invalid);
        }
        Ok(())
    }
}
fn receipt(face: FaceKey, raw: Values, options: FacePropertyOptions) -> Option<FaceProperties> {
    if !raw.into_iter().all(I::valid)
        || raw[0].lo <= 0.
        || raw[0].width() > options.area_tolerance_mm2
        || raw[1..]
            .iter()
            .any(|v| v.width() > options.first_moment_tolerance_mm3)
    {
        return None;
    }
    let centroid = [raw[1], raw[2], raw[3]].map(|v| v / raw[0]);
    if centroid
        .into_iter()
        .any(|v| !v.valid() || v.width() > options.centroid_tolerance_mm)
    {
        return None;
    }
    Some(FaceProperties {
        face,
        area_mm2: bounds(raw[0]),
        first_moment_mm3: [raw[1], raw[2], raw[3]].map(bounds),
        centroid_mm: centroid.map(bounds),
    })
}
impl Cavity {
    /// Includes the cutter's actual Mouth closure. Material screw publication
    /// must omit that face and retain the head's remaining annular land.
    pub fn face_properties_with_cancel(
        &self,
        options: FacePropertyOptions,
        cancelled: &dyn Fn() -> bool,
    ) -> Result<Vec<FaceProperties>, Error> {
        options.validate()?;
        let mut work = Work {
            remaining: options.maximum_cells,
            cancelled,
        };
        work.check()?;
        let atlas = self.source().surface_patches_with_cancel(cancelled)?;
        let mut result = Vec::with_capacity(atlas.topology.faces().len());
        for (face, source) in atlas.topology.faces().iter().enumerate() {
            let value = integrate::face(&atlas, face, source.source_identity(), options, &mut work);
            result.push(value?);
        }
        work.check()?;
        Ok(result)
    }
    pub fn local_face_properties_with_cancel(
        &self,
        face: FaceKey,
        options: FacePropertyOptions,
        cancelled: &dyn Fn() -> bool,
    ) -> Result<FaceProperties, Error> {
        options.validate()?;
        let mut work = Work {
            remaining: options.maximum_cells,
            cancelled,
        };
        work.check()?;
        let atlas = self.source().surface_patches_with_cancel(cancelled)?;
        let index = atlas
            .topology
            .faces()
            .iter()
            .position(|f| f.source_identity() == face)
            .ok_or(Error::Invalid)?;
        let result = integrate::face(&atlas, index, face, options, &mut work)?;
        work.check()?;
        Ok(result)
    }
}
