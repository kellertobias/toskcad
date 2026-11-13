//! Authored cavity cache and monotone whole-suffix envelope, without admission.
use super::*;
use std::sync::Arc;
mod host;
pub use host::HostEnvelopeProof;
mod screw;
pub use screw::{PhillipsScrew, ScrewFace, SourceVolumeBounds};
mod mouth;
pub use mouth::MouthSweepProof;

/// Qualified authored cavity. This retains the finite cutter, not an imported
/// body. Its volume is a source-expression enclosure; native publication,
/// head-mouth sewing and physical-file mass require separate proofs.
#[derive(Clone, Debug)]
pub struct Cavity {
    source: Arc<Cutter>,
    volume: VolumeBounds,
    side_faces: Arc<[FaceKey]>,
    mouth_edge_count: usize,
}
impl Cavity {
    pub fn new(source: Cutter, width_mm3: f64, work_limit: usize) -> Result<Self, Error> {
        Self::with_cancel(source, width_mm3, work_limit, &|| false)
    }
    pub fn with_cancel(
        source: Cutter,
        width_mm3: f64,
        work_limit: usize,
        cancelled: &dyn Fn() -> bool,
    ) -> Result<Self, Error> {
        // These independent source gates prove finite-cap eligibility and the
        // actual closed incidence, retaining both H0 and shallow-V supports.
        let topology = source.boundary_topology_with_cancel(cancelled)?;
        let mouth_edge_count = topology
            .edges()
            .iter()
            .filter(|e| matches!(e.source_identity(), CurveKey::Mouth { .. }))
            .count();
        let side_faces = topology
            .faces()
            .iter()
            .map(|f| f.source_identity())
            .filter(|key| *key != FaceKey::Mouth)
            .collect::<Vec<_>>()
            .into();
        let volume = source.volume_bounds_with_cancel(width_mm3, work_limit, cancelled)?;
        Ok(Self {
            source: Arc::new(source),
            volume,
            side_faces,
            mouth_edge_count,
        })
    }
    pub fn source(&self) -> &Cutter {
        &self.source
    }
    pub fn volume(&self) -> VolumeBounds {
        self.volume
    }
    pub fn side_face_keys(&self) -> &[FaceKey] {
        &self.side_faces
    }
    /// Upper radius of the ENTIRE cavity at and below `depth_mm`. All three
    /// actual envelope meridians decrease strictly with depth. Removing the
    /// finite corner prisms can only shrink their radial envelope. This is a
    /// conservative section bound, never a claim that the cavity is circular.
    /// At an uncertain source join both owning analytic expressions contribute.
    pub fn radial_suffix_bound_mm(&self, depth_mm: f64) -> Result<f64, Error> {
        if !depth_mm.is_finite() || depth_mm < 0. {
            return Err(Error::Invalid);
        }
        let e = &self.source.envelope;
        if depth_mm >= e.tip.hi {
            return Ok(0.);
        }
        let d = I::point(depth_mm);
        let mut upper: f64 = 0.;
        if depth_mm <= e.join.hi {
            let r = I::point(self.source.spec.mouth_radius_mm);
            let square = r.square() - (d - r).square();
            // The active roundover domain is 0<=d<=join<r. Its exact
            // radicand is nonnegative, including the mouth endpoint zero.
            let root = I::new(square.lo.max(0.), square.hi.max(0.)).sqrt();
            upper = upper.max((e.entry - root).hi);
        }
        if depth_mm >= e.join.lo && depth_mm <= e.flank.hi {
            upper = upper.max(
                (I::point(self.source.spec.mouth_diameter_mm) / I::point(2.)
                    - d * e.sin_a / e.cos_a)
                    .hi,
            );
        }
        if depth_mm >= e.flank.lo {
            upper = upper.max(((e.tip - d) * e.cos_c / e.sin_c).hi);
        }
        // The validated source has finite positive lengths and acute angles;
        // these branches cover every depth below tip, with finite denominators.
        Ok(upper)
    }
}

#[cfg(test)]
mod tests;
