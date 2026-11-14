//! Semantic faces of the actual paired authored material source.
use crate::authored_internal_treatment::ThreadBandRole;

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum FaceKey {
    StartCap,
    EndCap,
    OuterThread { start: u32, band: ThreadBandRole },
    InnerThread { start: u32, band: ThreadBandRole },
}

use super::Insert;
impl Insert {
    /// Semantic roles only; the finite boundary separately owns incidence.
    pub fn face_count(&self) -> u32 {
        2 + 4 * (self.outer().spec().starts + self.inner().spec().starts)
    }
    pub fn face_key(&self, index: u32) -> Option<FaceKey> {
        if index == 0 {
            return Some(FaceKey::StartCap);
        }
        if index == 1 {
            return Some(FaceKey::EndCap);
        }
        if index >= self.face_count() {
            return None;
        }
        let outer_count = 4 * self.outer().spec().starts;
        let index = index - 2;
        let outer = index < outer_count;
        let index = if outer { index } else { index - outer_count };
        let band = match index % 4 {
            0 => ThreadBandRole::Lower,
            1 => ThreadBandRole::Crest,
            2 => ThreadBandRole::Upper,
            _ => ThreadBandRole::Root,
        };
        let start = index / 4;
        Some(if outer {
            FaceKey::OuterThread { start, band }
        } else {
            FaceKey::InnerThread { start, band }
        })
    }
    pub fn face_index(&self, key: FaceKey) -> Option<u32> {
        (0..self.face_count()).find(|&index| self.face_key(index) == Some(key))
    }
}
