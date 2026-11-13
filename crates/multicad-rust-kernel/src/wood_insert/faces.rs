//! Semantic faces of the actual paired authored material source.
use crate::authored_internal_treatment::ThreadBandRole;

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum FaceKey {
    StartCap,
    EndCap,
    OuterThread { start: u32, band: ThreadBandRole },
    InnerThread { start: u32, band: ThreadBandRole },
}
