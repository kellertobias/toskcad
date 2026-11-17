//! Actual stock-minus-MAX-cavity membership of complete local boxes.
use super::*;
use crate::helical_thread::BoxClassification;
impl Material {
    /// Classify every point of a closed stock-relative box using its actual
    /// conical entry or terminal cylinder. Inside requires the ENTIRE box in
    /// the actual finite stock and outside the void; no point sampling or
    /// display geometry is involved. Other regions remain unresolved.
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
        let end = self.cavity().end_interval();
        let stock = self.radius_interval();
        if z.hi < 0. || z.lo > end.hi || rho.lo > stock.hi {
            return Ok(Outside);
        }
        if z.lo <= 0. || z.hi >= end.lo {
            return Ok(UnresolvedBoundary);
        }
        let inside_stock = rho.hi < stock.lo;
        let body = self.cavity();
        let radius = if let Some(tail) = body.tail
            && z.lo > tail.station
        {
            I::point(tail.radius)
        } else if let Some(entry) = body.entry
            && z.hi < entry.depth
        {
            let cone = I::point(entry.radius)
                - (I::point(entry.radius) - I::point(body.thread.minimum_radius_mm()))
                    / I::point(entry.depth)
                    * z;
            // All shared-band radius endpoints are <= this actual major
            // radius, so strict domination earns the complete MAX branch.
            if !cone.valid() || cone.lo <= body.thread.maximum_radius_mm() {
                return Ok(UnresolvedBoundary);
            }
            cone
        } else {
            return Ok(UnresolvedBoundary);
        };
        Ok(if rho.hi < radius.lo {
            Outside
        } else if inside_stock && rho.lo > radius.hi {
            Inside
        } else {
            UnresolvedBoundary
        })
    }
}

#[cfg(test)]
mod tests;
