//! Analytic normal descriptors and actual source second-order fields.
use super::{jet::Jet, *};
#[derive(Clone, Copy, Debug)]
pub(super) enum Kind {
    Cone,
    Wall,
}
#[derive(Clone, Copy)]
pub(super) enum Normal {
    Radial { angle: I, beta: I, z: I },
    Axial(I),
}
impl Normal {
    pub fn components(self, yaw: f64) -> [I; 3] {
        match self {
            Self::Axial(z) => [I::point(0.), I::point(0.), z],
            Self::Radial { angle, beta, z } => {
                let shifted = angle - I::point(yaw);
                let c = shifted.cos();
                let s = shifted.sin();
                [c - beta * s, s + beta * c, z]
            }
        }
    }
}
pub(super) struct Field {
    pub point: [Jet; 3],
    pub normal: Normal,
}
