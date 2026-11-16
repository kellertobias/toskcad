//! Authored pipe laws retain shared expressions, not rounded derived supports.
//! Literal imported/retained sections remain a different source type.
use super::{Error, I, bounds, ceiling, interval};
use crate::tapered_thread::Series;
use std::f64::consts::PI;
mod boundary;
pub mod hardware;
mod stock_clip;
pub use boundary::RawAuthoredBoundary;
pub use stock_clip::{
    PlacedStockClipMeshCandidate, StockClipAffineCurveTaylor, StockClipBand, StockClipCapTraceEdge,
    StockClipCellPlan, StockClipClassification, StockClipCoedge, StockClipCurve, StockClipCurveLaw,
    StockClipCurveTaylor, StockClipFace, StockClipFaceLoop, StockClipFaceMeasure,
    StockClipGridLine, StockClipJets, StockClipMeasureOptions, StockClipMeshCandidate,
    StockClipOptions, StockClipPatch, StockClipPatchLaw, StockClipPatchTaylor, StockClipPlacement,
    StockClipPlan, StockClipPoint, StockClipPointOwner, StockClipPolyline, StockClipQuery,
    StockClipRail, StockClipSample, StockClipSlab, StockClipStation, StockClipTracedMesh,
    StockClipVertex, StockClipVolume,
};
#[cfg(test)]
mod tests;
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum AuthoredSection {
    NptFlat {
        depth_mm: f64,
        flat_width_mm: f64,
    },
    /// Radius is derived from the exact closure equation; it is not an
    /// independently rounded physical input or an inferred imported circle.
    BsptRounded {
        depth_mm: f64,
    },
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AuthoredSpec {
    pub major_radius_mm: f64,
    pub pitch_mm: f64,
    pub length_mm: f64,
    pub radial_taper: f64,
    pub section: AuthoredSection,
    pub crest_radius_limit_mm: Option<f64>,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AuthoredThread {
    spec: AuthoredSpec,
    tolerance: f64,
    minimum: f64,
}
#[derive(Clone, Copy)]
struct Rounded {
    radius: I,
    sin: I,
    cos: I,
    m: I,
    normal: I,
    beta: I,
    valley: I,
}
fn pi() -> I {
    I::new(PI.next_down(), PI.next_up())
}
fn min(a: I, b: I) -> I {
    I::new(a.lo.min(b.lo), a.hi.min(b.hi))
}
fn max(a: I, b: I) -> I {
    I::new(a.lo.max(b.lo), a.hi.max(b.hi))
}
impl AuthoredThread {
    pub fn new(spec: AuthoredSpec, tolerance: f64, minimum: f64) -> Result<Self, Error> {
        let AuthoredSpec {
            major_radius_mm: r,
            pitch_mm: p,
            length_mm: l,
            radial_taper: k,
            ..
        } = spec;
        let d = match spec.section {
            AuthoredSection::NptFlat { depth_mm, .. }
            | AuthoredSection::BsptRounded { depth_mm } => depth_mm,
        };
        if [r, p, l, d, tolerance, minimum]
            .into_iter()
            .any(|x| !x.is_finite() || x <= 0.)
            || !k.is_finite()
            || k.abs() != 1. / 32.
            || r <= d
        {
            return Err(Error::Invalid);
        }
        if [r, p, l, d].into_iter().any(|x| x > 1e7) || p < 0.001 || l < p || l / p + 2. > 128. {
            return Err(Error::Limit);
        }
        let reserve = 256. * f64::EPSILON * (r + p + l + d);
        if tolerance <= 4. * reserve || r - d + (k * l).min(0.) <= minimum + 4. * reserve {
            return Err(Error::Unresolved);
        }
        if let Some(limit) = spec.crest_radius_limit_mm {
            if !limit.is_finite() || limit <= 0. {
                return Err(Error::Invalid);
            }
            if limit > 1e7 {
                return Err(Error::Limit);
            }
            if limit <= r - d + (k * l).max(0.) + minimum {
                return Err(Error::Unresolved);
            }
        }
        let result = Self {
            spec,
            tolerance,
            minimum,
        };
        match spec.section {
            AuthoredSection::NptFlat {
                flat_width_mm: w, ..
            } => {
                if !w.is_finite() || w <= 0. {
                    return Err(Error::Invalid);
                }
                let (a, _) = result.endpoint(0, false)?;
                let (b, _) = result.endpoint(2, true)?;
                if (I::point(p) - (b - a)).lo <= minimum {
                    return Err(Error::Unresolved);
                }
            }
            AuthoredSection::BsptRounded { .. } => {
                let b = result.rounded()?;
                if b.radius.lo <= minimum
                    || (min(b.valley, I::point(p) - b.valley) - I::point(2.) * b.radius * b.cos).lo
                        <= minimum
                {
                    return Err(Error::Unresolved);
                }
            }
        }
        Ok(result)
    }
    pub fn spec(&self) -> AuthoredSpec {
        self.spec
    }
    pub fn tolerance_mm(&self) -> f64 {
        self.tolerance
    }
    pub fn minimum_feature_mm(&self) -> f64 {
        self.minimum
    }
    pub fn band_count(&self) -> usize {
        match self.spec.section {
            AuthoredSection::NptFlat { .. } => 3,
            AuthoredSection::BsptRounded { .. } => 5,
        }
    }
    fn depth(&self) -> f64 {
        match self.spec.section {
            AuthoredSection::NptFlat { depth_mm, .. }
            | AuthoredSection::BsptRounded { depth_mm } => depth_mm,
        }
    }
    fn rounded(&self) -> Result<Rounded, Error> {
        let a = I::point(27.5) * pi() / I::point(180.);
        let sin = a.sin();
        let cos = a.cos();
        let m = cos / sin;
        let k = I::point(self.spec.radial_taper);
        let normal = (I::point(1.) + k.square()).sqrt();
        let h = I::point(self.spec.pitch_mm) * (m.square() - k.square()) / (I::point(2.) * m);
        let radius = ((h - I::point(self.depth()))
            / (I::point(2.) * (I::point(1.) / sin - normal)))
            .positive()?;
        Ok(Rounded {
            radius,
            sin,
            cos,
            m,
            normal,
            beta: pi() / I::point(2.) - a,
            valley: h / (m + k),
        })
    }
    pub fn derived_rounding_radius_bounds_mm(&self) -> Result<Option<[f64; 2]>, Error> {
        match self.spec.section {
            AuthoredSection::NptFlat { .. } => Ok(None),
            AuthoredSection::BsptRounded { .. } => Ok(Some(bounds(self.rounded()?.radius))),
        }
    }
    fn endpoint(&self, band: usize, end: bool) -> Result<(I, I), Error> {
        if band >= self.band_count() {
            return Err(Error::Invalid);
        }
        let i = band + usize::from(end);
        let r = I::point(self.spec.major_radius_mm);
        let d = I::point(self.depth());
        let p = I::point(self.spec.pitch_mm);
        let k = I::point(self.spec.radial_taper);
        Ok(match self.spec.section {
            AuthoredSection::NptFlat {
                flat_width_mm: w, ..
            } => {
                let h = I::point(w) / I::point(2.);
                let m = I::point(3.).sqrt();
                match i {
                    0 => (-h - d / (m - k), r - d),
                    1 => (-h, r),
                    2 => (h, r),
                    _ => (h + d / (m + k), r - d),
                }
            }
            AuthoredSection::BsptRounded { .. } => {
                let b = self.rounded()?;
                let h = b.radius * b.cos;
                let crest = r - b.radius * b.normal;
                let tangent = crest + b.radius * b.sin;
                match i {
                    0 => (I::point(0.), crest + b.radius),
                    1 => (h, tangent - k * h),
                    2 => (
                        b.valley - h,
                        r - d + b.radius * b.normal - b.radius * b.sin + k * h,
                    ),
                    3 => (
                        b.valley + h,
                        r - d + b.radius * b.normal - b.radius * b.sin - k * h,
                    ),
                    4 => (p - h, tangent + k * h),
                    _ => (p, crest + b.radius),
                }
            }
        })
    }
    fn profile<const N: usize>(
        &self,
        band: usize,
        t: Series<N>,
    ) -> Result<(Series<N>, Series<N>), Error> {
        if band >= self.band_count() {
            return Err(Error::Invalid);
        }
        let c = Series::interval;
        let r = c(I::point(self.spec.major_radius_mm));
        let d = c(I::point(self.depth()));
        let k = c(I::point(self.spec.radial_taper));
        let p = c(I::point(self.spec.pitch_mm));
        let (mut q, mut a) = match self.spec.section {
            AuthoredSection::NptFlat {
                flat_width_mm: w, ..
            } => {
                let h = c(I::point(w) / I::point(2.));
                let m = c(I::point(3.).sqrt());
                match band {
                    0 => (-h - d / (m - k) + d / (m - k) * t, r - d + d * t),
                    1 => (-h + c(I::point(w)) * t, r),
                    _ => (h + d / (m + k) * t, r - d * t),
                }
            }
            AuthoredSection::BsptRounded { .. } => {
                let b = self.rounded()?;
                let radius = c(b.radius);
                let h = radius * c(b.cos);
                let v = c(b.valley);
                let center = r - radius * c(b.normal);
                let tangent = center + radius * c(b.sin);
                match band {
                    0 | 2 | 4 => {
                        let angle = match band {
                            0 => c(b.beta) * t,
                            2 => -c(b.beta) + c(I::point(2.) * b.beta) * t,
                            _ => -c(b.beta) + c(b.beta) * t,
                        };
                        let [sin, cos] = angle.sin_cos();
                        let x = radius * sin;
                        match band {
                            0 => (x, center + radius * cos - k * x),
                            2 => (v + x, r - d + radius * c(b.normal) - radius * cos - k * x),
                            _ => (p + x, center + radius * cos - k * x),
                        }
                    }
                    1 => {
                        let q = h + (v - c(I::point(2.)) * h) * t;
                        (q, tangent - k * h - (c(b.m) + k) * (q - h))
                    }
                    _ => {
                        let q = v + h + (p - v - c(I::point(2.)) * h) * t;
                        (q, tangent + k * h + (c(b.m) - k) * (q - (p - h)))
                    }
                }
            }
        };
        // These are algebraic identities of the authored law, not endpoint
        // fitting of the separate literal-support Thread.
        if t.c[0].lo == t.c[0].hi && (t.c[0].lo == 0. || t.c[0].lo == 1.) {
            let (eq, ea) = self.endpoint(band, t.c[0].lo == 1.)?;
            q.c[0] = eq;
            a.c[0] = ea;
        }
        Ok((q.checked()?, a.checked()?))
    }
    pub fn surface_point_bounds(
        &self,
        band: usize,
        parameter: [f64; 2],
        axial_mm: [f64; 2],
    ) -> Result<[[f64; 2]; 3], Error> {
        let t = interval(parameter)?;
        let z = interval(axial_mm)?;
        if t.lo < 0. || t.hi > 1. || z.lo < 0. || z.hi > self.spec.length_mm {
            return Err(Error::Invalid);
        }
        let (q, a) = self.profile(band, Series::<2>::interval(t))?;
        let radius = ceiling(
            a.c[0] + I::point(self.spec.radial_taper) * z,
            self.spec.crest_radius_limit_mm,
        );
        let theta = (z - q.c[0]) * (I::point(2.) * pi() / I::point(self.spec.pitch_mm));
        Ok([
            bounds(radius * theta.cos()),
            bounds(radius * theta.sin()),
            bounds(z),
        ])
    }
    pub fn radial_bounds(&self, phase: [f64; 2], axial: [f64; 2]) -> Result<[f64; 2], Error> {
        let q = interval(phase)?;
        let z = interval(axial)?;
        if q.lo < 0. || q.hi > self.spec.pitch_mm || z.lo < 0. || z.hi > self.spec.length_mm {
            return Err(Error::Invalid);
        }
        let r = I::point(self.spec.major_radius_mm);
        let d = I::point(self.depth());
        let k = I::point(self.spec.radial_taper);
        let major = r + k * z;
        let value = match self.spec.section {
            AuthoredSection::NptFlat {
                flat_width_mm: w, ..
            } => {
                let h = I::point(w) / I::point(2.);
                let m = I::point(3.).sqrt();
                max(
                    major - d,
                    min(
                        major,
                        max(
                            major - (m + k) * (q - h),
                            major + (m - k) * (q - (I::point(self.spec.pitch_mm) - h)),
                        ),
                    ),
                )
            }
            AuthoredSection::BsptRounded { .. } => {
                let b = self.rounded()?;
                let h = b.radius * b.cos;
                let center = r - b.radius * b.normal;
                let tangent = center + b.radius * b.sin;
                let stations = [
                    I::point(0.),
                    h,
                    b.valley - h,
                    b.valley + h,
                    I::point(self.spec.pitch_mm) - h,
                    I::point(self.spec.pitch_mm),
                ];
                let mut all = None;
                for i in 0..5 {
                    let lo = q.lo.max(stations[i].lo);
                    let hi = q.hi.min(stations[i + 1].hi);
                    if lo > hi {
                        continue;
                    }
                    let phase = I::new(lo, hi);
                    let a = match i {
                        0 | 4 => {
                            let x = phase
                                - if i == 0 {
                                    I::point(0.)
                                } else {
                                    I::point(self.spec.pitch_mm)
                                };
                            center + (b.radius.square() - x.square()).sqrt() - k * x
                        }
                        1 => tangent - k * h - (b.m + k) * (phase - h),
                        2 => {
                            let x = phase - b.valley;
                            r - d + b.radius * b.normal
                                - (b.radius.square() - x.square()).sqrt()
                                - k * x
                        }
                        _ => {
                            tangent
                                + k * h
                                + (b.m - k) * (phase - (I::point(self.spec.pitch_mm) - h))
                        }
                    } + k * z;
                    all = Some(all.map_or(a, |old: I| I::new(old.lo.min(a.lo), old.hi.max(a.hi))));
                }
                all.ok_or(Error::Unresolved)?
            }
        };
        Ok(bounds(
            ceiling(value, self.spec.crest_radius_limit_mm).checked()?,
        ))
    }
}
