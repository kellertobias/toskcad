//! Regular closed difference of finite stock and an unchanged shorter cutter.
//! The floor is bounded by the actual noncircular cutter-tip section. No
//! nominal major-radius circle defines that planar trimmed face.
//! The finite entry/floor/stock caps belong to the regular boundary; whole-box
//! predicates deliberately leave all touching queries unresolved.
use super::{Error, Thread};
use crate::{Rational, helical_thread::BoxClassification, tapered_thread::interval::I};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BlindBore {
    cutter: Thread,
    radius: Rational,
    height: Rational,
    radial_support: [f64; 2],
}

fn meridian_bounds(thread: &Thread) -> Result<I, Error> {
    let mut result = I::new(f64::INFINITY, f64::NEG_INFINITY);
    for band in 0..5 {
        let phases = I::new(thread.breaks[band], thread.breaks[band + 1]);
        let r = if let Some(circle) = thread.profile_circle(band) {
            // Retain the true symbolic-circle coefficient uncertainty, not
            // just its rounded cached centre or nominal major radius.
            let center = I::point(circle.center_radius_mm)
                + I::new(-circle.center_error_mm, circle.center_error_mm);
            let x = phases - I::point(circle.center_phase_mm);
            let radius = I::point(circle.radius_mm);
            // Squaring the SAME phase interval preserves its zero crossing;
            // two independent (r-x)/(r+x) hulls would needlessly inflate the
            // crest maximum and reject otherwise sound thin stock walls.
            let h2 = radius.square() - x.square();
            if !h2.valid() || h2.lo < 0. {
                return Err(Error::Unresolved);
            }
            center + I::point(circle.radial_sign) * h2.sqrt()
        } else {
            let (anchor, radius, slope) = if band == 1 {
                (thread.breaks[1], thread.tangent_radius[0], thread.slope)
            } else {
                (thread.breaks[3], thread.tangent_radius[1], -thread.slope)
            };
            // The entire affine interval includes BOTH one-sided endpoint
            // laws, independently of radius_at_phase's branch tie-breaking.
            I::point(radius) + I::point(slope) * (phases - I::point(anchor))
        };
        if !r.valid() || r.lo <= 0. {
            return Err(Error::Unresolved);
        }
        result.lo = result.lo.min(r.lo);
        result.hi = result.hi.max(r.hi);
    }
    Ok(result)
}

fn exact(q: Rational) -> I {
    let (n, d) = q.ratio();
    let bits = n.unsigned_abs();
    if (d as u64).is_power_of_two()
        && (bits == 0 || 64 - bits.leading_zeros() - bits.trailing_zeros() <= 53)
    {
        I::point(q.as_f64())
    } else {
        // Rational::as_f64 rounds the actual ratio once; its neighbours
        // enclose it even when either integer operand exceeds binary64.
        let x = q.as_f64();
        I::new(x.next_down(), x.next_up())
    }
}

impl BlindBore {
    /// Preserve the original full cutter, independently of the exact stock
    /// extent. A strictly positive wall and floor must be proved before this
    /// source exists. Through and overrun cases remain separate contracts.
    pub fn new(cutter: Thread, radius: Rational, height: Rational) -> Result<Self, Error> {
        let r = exact(radius);
        let h = exact(height);
        if !r.valid() || !h.valid() || r.lo <= 0. || h.lo <= 0. || r.hi > 1e7 || h.hi > 1e7 {
            return Err(Error::Invalid);
        }
        let radial = meridian_bounds(&cutter)?;
        if r.lo <= radial.hi || h.lo <= cutter.spec().length_mm {
            return Err(Error::Unresolved);
        }
        Ok(Self {
            cutter,
            radius,
            height,
            radial_support: [radial.lo, radial.hi],
        })
    }
    pub fn cutter(&self) -> &Thread {
        &self.cutter
    }
    pub fn stock_radius_mm(&self) -> Rational {
        self.radius
    }
    pub fn stock_height_mm(&self) -> Rational {
        self.height
    }
    /// Complete actual circle/affine meridian support, with symbolic circle
    /// centre uncertainty and both one-sided laws at every finite join.
    pub fn cutter_radial_bounds_mm(&self) -> [f64; 2] {
        self.radial_support
    }
    pub fn stock_radius_bounds_mm(&self) -> [f64; 2] {
        let r = exact(self.radius);
        [r.lo, r.hi]
    }
    pub fn stock_height_bounds_mm(&self) -> [f64; 2] {
        let h = exact(self.height);
        [h.lo, h.hi]
    }
    pub fn floor_thickness_bounds_mm(&self) -> [f64; 2] {
        let h = exact(self.height) - I::point(self.cutter.spec().length_mm);
        [h.lo, h.hi]
    }
    /// Directed support of the actual finite stock. Every outer cylindrical
    /// support point remains occupied: the cutter's complete radial support
    /// is strictly smaller, and the whole stock top lies above its finite tip.
    /// Thus this is a material support proof, not just a containing box.
    pub fn support_bounds(&self, direction: [f64; 3], tolerance: f64) -> Result<[f64; 2], Error> {
        if direction.iter().any(|v| !v.is_finite()) || !tolerance.is_finite() || tolerance <= 0. {
            return Err(Error::Invalid);
        }
        if direction == [0.; 3] {
            return Ok([0.; 2]);
        }
        let [x, y, z] = direction.map(I::point);
        let radial = (x.square() + y.square()).sqrt() * exact(self.radius);
        let axial = if direction[2] > 0. {
            z * exact(self.height)
        } else {
            I::point(0.)
        };
        let support = radial + axial;
        if !support.valid() || support.width() > tolerance {
            return Err(Error::Unresolved);
        }
        Ok([support.lo, support.hi])
    }
    /// Whole-box classification only where stock or the solid floor decides
    /// membership. The original helical region and every touching boundary
    /// stay unresolved; these shortcuts never replace its finite rim proof.
    pub fn classify_floor_box(&self, bounds: [[f64; 2]; 3]) -> Result<BoxClassification, Error> {
        use BoxClassification::*;
        let [x, y, z] = bounds.map(|[lo, hi]| I::new(lo, hi));
        if [x, y, z].iter().any(|v| !v.valid()) {
            return Err(Error::Invalid);
        }
        let rho2 = x.square() + y.square();
        let radius2 = exact(self.radius).square();
        let height = exact(self.height);
        if !rho2.valid() || !radius2.valid() {
            return Ok(UnresolvedBoundary);
        }
        Ok(if z.hi < 0. || z.lo > height.hi || rho2.lo > radius2.hi {
            Outside
        } else if z.lo > self.cutter.spec().length_mm && z.hi < height.lo && rho2.hi < radius2.lo {
            Inside
        } else {
            UnresolvedBoundary
        })
    }
}

#[cfg(test)]
mod tests;
