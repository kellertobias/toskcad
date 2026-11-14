//! Actual local source integrals, independent of tessellation and placement.
use super::*;
use crate::authored_internal_treatment::{self as measure, properties::bare_properties};
pub use measure::FacePropertyOptions;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FaceProperties {
    pub face: FaceKey,
    pub area_mm2: measure::Enclosure,
    pub first_moment_mm3: [measure::Enclosure; 3],
    pub centroid_mm: [measure::Enclosure; 3],
}
impl Insert {
    /// Bounded area and area-weighted moments of the actual finite source.
    /// These local measures do not imply area invariance after placement.
    pub fn local_face_properties(
        &self,
        face: FaceKey,
        options: FacePropertyOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<FaceProperties, Error> {
        let (thread, subtract, key) = match face {
            FaceKey::StartCap => (self.outer(), Some(self.inner()), measure::FaceKey::StartCap),
            FaceKey::EndCap => (self.outer(), Some(self.inner()), measure::FaceKey::EndCap),
            FaceKey::OuterThread { start, band } => {
                (self.outer(), None, measure::FaceKey::Thread { start, band })
            }
            FaceKey::InnerThread { start, band } => {
                (self.inner(), None, measure::FaceKey::Thread { start, band })
            }
        };
        let result = bare_properties(thread, subtract, key, options, interrupt).map_err(failure)?;
        Ok(FaceProperties {
            face,
            area_mm2: result.area_mm2,
            first_moment_mm3: result.first_moment_mm3,
            centroid_mm: result.centroid_mm,
        })
    }
}
fn failure(error: measure::Error) -> Error {
    match error {
        measure::Error::Invalid => Error::Invalid,
        measure::Error::Unresolved => Error::Unresolved,
        measure::Error::Limit => Error::Limit,
        measure::Error::Cancelled => Error::Cancelled,
    }
}
#[cfg(test)]
mod tests;
