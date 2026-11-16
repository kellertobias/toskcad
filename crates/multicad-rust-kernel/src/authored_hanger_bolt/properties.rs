//! Actual local area measures; a later affine placement needs a separate
//! normal-dependent area and centroid transport certificate.
use super::FaceKey;
use crate::authored_internal_treatment::Enclosure;
pub use crate::authored_internal_treatment::FacePropertyOptions;
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FaceProperties {
    pub face: FaceKey,
    pub area_mm2: Enclosure,
    pub first_moment_mm3: [Enclosure; 3],
    pub centroid_mm: [Enclosure; 3],
}

#[cfg(test)]
mod tests;
