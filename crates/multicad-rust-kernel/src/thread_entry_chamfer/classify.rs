//! Whole-box predicates on actual conical entry and constant terminal material.
use super::*;
use crate::helical_thread::BoxClassification;
use crate::tapered_thread::interval::I;
impl Chamfer {
    /// Classify EVERY point of a local closed box when an actual treatment
    /// dominates the original source over its whole axial interval. This
    /// bounded predicate does not construct an intersected boundary. Ordinary
    /// helical regions, touching boxes and uncertain arithmetic stay unresolved.
    pub fn classify_treatment_box(
        &self,
        bounds: [[f64; 2]; 3],
    ) -> Result<BoxClassification, Error> {
        use BoxClassification::*;
        let p = bounds.map(|[lo, hi]| I::new(lo, hi));
        if p.iter().any(|p| !p.valid()) {
            return Err(Error::Invalid);
        }
        let [x, y, z] = p;
        let rho = (x.square() + y.square()).sqrt();
        if !rho.valid() {
            return Ok(UnresolvedBoundary);
        }
        if z.hi < 0. || z.lo > self.full_length_mm() {
            return Ok(Outside);
        }
        if z.lo <= 0. || z.hi >= self.full_length_mm() {
            return Ok(UnresolvedBoundary);
        }
        let radius = if let Some(tail) = self.terminal_relief()
            && z.lo > tail.material_length_mm()
        {
            I::point(tail.neck_radius_mm())
        } else if z.hi < self.entry_length_mm() {
            let cone = I::point(self.entry_radius_mm())
                + (I::point(self.tool.spec().major_radius_mm) - I::point(self.entry_radius_mm()))
                    / I::point(self.entry_length_mm())
                    * z;
            // Include BOTH actual one-sided laws at all original profile
            // joins. A nominal root-radius or sampled corner is insufficient.
            let base = self
                .tool
                .radius_at_phase_interval(I::new(0., self.tool.spec().pitch_mm))
                .ok_or(Error::Unresolved)?;
            if !cone.valid() || cone.hi >= base.lo {
                return Ok(UnresolvedBoundary);
            }
            cone
        } else {
            return Ok(UnresolvedBoundary);
        };
        Ok(if rho.hi < radius.lo {
            Inside
        } else if rho.lo > radius.hi {
            Outside
        } else {
            UnresolvedBoundary
        })
    }
}
