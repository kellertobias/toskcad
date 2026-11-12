//! Single right-hand axial pipe sections, distinct from cone-spine transport.
//! Finite line/circle supports and their radial junctions are retained literally.
//! This source is not an admission certificate for a manufacturing boundary.
pub use crate::tapered_thread::Error;
use crate::tapered_thread::interval::I;
use std::f64::consts::PI;
#[cfg(test)]
mod boundary;
mod properties;
#[cfg(test)]
mod tests;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Section {
    NptFlat {
        depth_mm: f64,
        flat_width_mm: f64,
    },
    BsptRounded {
        depth_mm: f64,
        rounding_radius_mm: f64,
    },
}
/// Physical resolved dimensions. Size, gauge and standards metadata belong to
/// the authoring layer and do not redefine these retained supports.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Spec {
    pub major_radius_mm: f64,
    pub pitch_mm: f64,
    pub length_mm: f64,
    /// Signed radius change per axial millimetre, exactly +/-1/32.
    pub radial_taper: f64,
    pub section: Section,
    /// Actual cylindrical pipe-stock ceiling, supported for NPT only.
    pub crest_radius_limit_mm: Option<f64>,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Curve {
    Line {
        phase_mm: f64,
        radius_mm: f64,
        slope: f64,
    },
    Circle {
        phase_mm: f64,
        radius_mm: f64,
        circle_radius_mm: f64,
        sign: f64,
    },
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Band {
    pub phase_mm: [f64; 2],
    pub curve: Curve,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Location {
    Inside,
    Boundary,
    Outside,
    Unresolved,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Junction {
    pub phase_mm: f64,
    /// Both incident literal supports, before the common axial cone shift.
    pub radial_bounds_mm: [f64; 2],
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Thread {
    spec: Spec,
    bands: [Band; 5],
    reserve: f64,
    tolerance: f64,
    minimum: f64,
}
fn interval(v: [f64; 2]) -> Result<I, Error> {
    if v.into_iter().all(f64::is_finite) && v[0] <= v[1] {
        Ok(I::new(v[0], v[1]))
    } else {
        Err(Error::Invalid)
    }
}
fn bounds(v: I) -> [f64; 2] {
    [v.lo, v.hi]
}
fn hull(a: I, b: I) -> I {
    I::new(a.lo.min(b.lo), a.hi.max(b.hi))
}
fn ceiling(a: I, limit: Option<f64>) -> I {
    limit.map_or(a, |r| I::new(a.lo.min(r), a.hi.min(r)))
}
impl Curve {
    fn radius(self, phase: I) -> I {
        match self {
            Self::Line {
                phase_mm,
                radius_mm,
                slope,
            } => I::point(radius_mm) + I::point(slope) * (phase - I::point(phase_mm)),
            Self::Circle {
                phase_mm,
                radius_mm,
                circle_radius_mm: r,
                sign,
            } => {
                I::point(radius_mm)
                    + I::point(sign)
                        * (I::point(r).square() - (phase - I::point(phase_mm)).square()).sqrt()
            }
        }
    }
}
impl Thread {
    pub fn new(spec: Spec, tolerance: f64, minimum: f64) -> Result<Self, Error> {
        let Spec {
            major_radius_mm: r,
            pitch_mm: p,
            length_mm: l,
            radial_taper: k,
            ..
        } = spec;
        let (depth, extra) = match spec.section {
            Section::NptFlat {
                depth_mm,
                flat_width_mm,
            } => (depth_mm, flat_width_mm),
            Section::BsptRounded {
                depth_mm,
                rounding_radius_mm,
            } => (depth_mm, rounding_radius_mm),
        };
        if [r, p, l, depth, extra, tolerance, minimum]
            .into_iter()
            .any(|x| !x.is_finite() || x <= 0.)
            || !k.is_finite()
            || k.abs() != 1. / 32.
            || r <= depth
        {
            return Err(Error::Invalid);
        }
        if [r, p, l, depth, extra].into_iter().any(|x| x > 1e7)
            || p < 0.001
            || l < p
            || l / p + 2. > 128.
        {
            return Err(Error::Limit);
        }
        let reserve = 256. * f64::EPSILON * (r + p + l + depth + extra);
        if tolerance <= 4. * reserve || r - depth + (k * l).min(0.) <= minimum + 4. * reserve {
            return Err(Error::Unresolved);
        }
        if let Some(limit) = spec.crest_radius_limit_mm {
            if !limit.is_finite() || limit <= 0. {
                return Err(Error::Invalid);
            }
            if limit > 1e7 {
                return Err(Error::Limit);
            }
            if !matches!(spec.section, Section::NptFlat { .. }) {
                return Err(Error::Unsupported);
            }
            if limit <= r - depth + (k * l).max(0.) + minimum + 4. * reserve {
                return Err(Error::Unresolved);
            }
        }
        let line = |a, b, q, rad, slope| Band {
            phase_mm: [a, b],
            curve: Curve::Line {
                phase_mm: q,
                radius_mm: rad,
                slope,
            },
        };
        let circle = |a, b, q, rad, cr, sign| Band {
            phase_mm: [a, b],
            curve: Curve::Circle {
                phase_mm: q,
                radius_mm: rad,
                circle_radius_mm: cr,
                sign,
            },
        };
        let bands = match spec.section {
            Section::NptFlat {
                flat_width_mm: w, ..
            } => {
                let m = 3_f64.sqrt();
                let h = w / 2.;
                let a = h + depth / (m + k);
                let b = p - h - depth / (m - k);
                [
                    line(0., h, 0., r, k),
                    line(h, a, h, r + k * h, -m),
                    line(a, b, 0., r - depth, k),
                    line(b, p - h, p - h, r + k * (p - h), m),
                    line(p - h, p, 0., r, k),
                ]
            }
            Section::BsptRounded {
                rounding_radius_mm: cr,
                ..
            } => {
                let angle = 27.5_f64.to_radians();
                let (sin, cos) = angle.sin_cos();
                let m = 1. / angle.tan();
                let q = (1. + k * k).sqrt();
                let fundamental = p * (m * m - k * k) / (2. * m);
                let expected = fundamental - 2. * cr * (1. / sin - q);
                if (expected - depth).abs() > reserve {
                    return Err(Error::Invalid);
                }
                let h = cr * cos;
                let v = fundamental / (m + k);
                let center = r - cr * q;
                let tangent = center + cr * sin;
                [
                    circle(0., h, 0., center, cr, 1.),
                    line(h, v - h, h, tangent, -m),
                    circle(v - h, v + h, v, r - depth + k * v + cr * q, cr, -1.),
                    line(v + h, p - h, p - h, tangent + k * p, m),
                    circle(p - h, p, p, center + k * p, cr, 1.),
                ]
            }
        };
        if bands
            .iter()
            .any(|b| b.phase_mm[1] - b.phase_mm[0] <= reserve)
        {
            return Err(Error::Unresolved);
        }
        let thread = Self {
            spec,
            bands,
            reserve,
            tolerance,
            minimum,
        };
        // Keep both literal endpoint values. No average or tolerance snap changes
        // a circle support; a finite boundary must retain these radial closures.
        if thread
            .junctions()
            .iter()
            .any(|j| j.radial_bounds_mm[1] - j.radial_bounds_mm[0] > reserve)
        {
            return Err(Error::Unresolved);
        }
        Ok(thread)
    }
    pub fn spec(&self) -> Spec {
        self.spec
    }
    pub fn bands(&self) -> &[Band; 5] {
        &self.bands
    }
    /// Intrinsic construction budget bounding the actual junction enclosures.
    /// It never moves a support or identifies distinct positive material.
    pub fn coefficient_error_mm(&self) -> f64 {
        self.reserve
    }
    pub fn tolerance_mm(&self) -> f64 {
        self.tolerance
    }
    pub fn minimum_feature_mm(&self) -> f64 {
        self.minimum
    }
    pub fn junctions(&self) -> [Junction; 5] {
        std::array::from_fn(|i| {
            let q = self.bands[i].phase_mm[1];
            let a = self.unclipped(i, I::point(q), I::point(0.));
            let next = (i + 1) % 5;
            let b = self.unclipped(next, I::point(if next == 0 { 0. } else { q }), I::point(0.));
            Junction {
                phase_mm: q,
                radial_bounds_mm: bounds(hull(a, b)),
            }
        })
    }
    fn unclipped(&self, band: usize, q: I, z: I) -> I {
        self.bands[band].curve.radius(q) + I::point(self.spec.radial_taper) * (z - q)
    }
    fn radius_interval(&self, q: I, z: I) -> Result<I, Error> {
        let mut result = None;
        for (i, b) in self.bands.iter().enumerate() {
            let lo = q.lo.max(b.phase_mm[0]);
            let hi = q.hi.min(b.phase_mm[1]);
            if lo <= hi {
                let r = self.unclipped(i, I::new(lo, hi), z);
                result = Some(result.map_or(r, |old| hull(old, r)));
            }
        }
        if q.lo == 0. || q.hi == self.spec.pitch_mm {
            let ends = hull(
                self.unclipped(0, I::point(0.), z),
                self.unclipped(4, I::point(self.spec.pitch_mm), z),
            );
            result = Some(result.map_or(ends, |r| hull(r, ends)));
        }
        result
            .map(|r| ceiling(r, self.spec.crest_radius_limit_mm))
            .ok_or(Error::Invalid)?
            .checked()
    }
    pub fn radial_bounds(&self, phase_mm: [f64; 2], axial_mm: [f64; 2]) -> Result<[f64; 2], Error> {
        let q = interval(phase_mm)?;
        let z = interval(axial_mm)?;
        if q.lo < 0. || q.hi > self.spec.pitch_mm || z.lo < 0. || z.hi > self.spec.length_mm {
            return Err(Error::Invalid);
        }
        Ok(bounds(self.radius_interval(q, z)?))
    }
    /// Material query in the retained helical chart q=z-pitch*theta/(2*pi).
    /// Ambiguous numerical contact remains Unresolved, never fabricated contact.
    pub fn locate_radial(&self, phase_mm: f64, z: f64, radius_mm: f64) -> Result<Location, Error> {
        if !radius_mm.is_finite() || radius_mm < 0. || !z.is_finite() || !phase_mm.is_finite() {
            return Err(Error::Invalid);
        }
        if z < 0. || z > self.spec.length_mm {
            return Ok(Location::Outside);
        }
        let r = self.radial_bounds([phase_mm; 2], [z; 2])?;
        Ok(if radius_mm > r[1] {
            Location::Outside
        } else if radius_mm < r[0] {
            if z == 0. || z == self.spec.length_mm {
                Location::Boundary
            } else {
                Location::Inside
            }
        } else if r[0] == r[1] && radius_mm == r[0] {
            Location::Boundary
        } else {
            Location::Unresolved
        })
    }
    /// Actual wall parameterization; the radius ceiling is a physical pipe OD.
    pub fn surface_point_bounds(
        &self,
        band: usize,
        phase_mm: [f64; 2],
        axial_mm: [f64; 2],
    ) -> Result<[[f64; 2]; 3], Error> {
        let b = self.bands.get(band).ok_or(Error::Invalid)?;
        let q = interval(phase_mm)?;
        let z = interval(axial_mm)?;
        if q.lo < b.phase_mm[0] || q.hi > b.phase_mm[1] || z.lo < 0. || z.hi > self.spec.length_mm {
            return Err(Error::Invalid);
        }
        let r = ceiling(self.unclipped(band, q, z), self.spec.crest_radius_limit_mm).checked()?;
        let pi = I::new(PI.next_down(), PI.next_up());
        let angle = (z - q) * (I::point(2.) * pi / I::point(self.spec.pitch_mm));
        Ok([bounds(r * angle.cos()), bounds(r * angle.sin()), bounds(z)])
    }
}
