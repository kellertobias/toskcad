//! Source proof of a radial mouth sweep, distinct from literal ring checks.
use super::*;
/// The authored mouth bounds a star-shaped disk, with one ordered angular
/// sweep. This witness concerns the exact source, not sampled or emitted
/// vertices: a literal zipper must separately retain and qualify their order.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MouthSweepProof {
    minimum_front_tangent_intercept_mm: f64,
    maximum_radius_mm: f64,
    source_edges: usize,
}
impl MouthSweepProof {
    /// Lower bound on x-y*dx/dy for every inner/outer front, including H0's
    /// inclined cylinder section, throughout the whole cavity depth domain.
    pub fn minimum_front_tangent_intercept_mm(&self) -> f64 {
        self.minimum_front_tangent_intercept_mm
    }
    pub fn maximum_radius_mm(&self) -> f64 {
        self.maximum_radius_mm
    }
    pub fn source_edge_count(&self) -> usize {
        self.source_edges
    }
}
impl Cavity {
    pub fn mouth_sweep(&self) -> MouthSweepProof {
        // Volume admission has already proved these exact shared intercepts
        // strictly positive, finite caps/closure inactive, and adjacent corner
        // cutouts disjoint. rb decreases with depth, so the tip is worst.
        //
        // Each V front has x-y*x' equal to its positive intercept. For the
        // H0 ellipse it decreases monotonically from rb to the outer tangent
        // intercept at the shared join. Every circle piece has radial cross
        // product R²>0. The sealed section wire joins these supports with
        // exact quarter rotations. Thus each ray meets its boundary once,
        // in source wire order, rather than merely in sampled-angle order.
        let [inner, _, outer] = volume::radial_intercepts(&self.source);
        MouthSweepProof {
            minimum_front_tangent_intercept_mm: inner.lo.min(outer.lo),
            maximum_radius_mm: self.source.envelope.entry.hi,
            source_edges: self.mouth_edge_count,
        }
    }
}

#[cfg(test)]
mod tests;
