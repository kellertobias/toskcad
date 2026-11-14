//! Finite root-face trimming of the SAME owned Hanger Boundary. No literal
//! curves, sewing, closed shell, file or occupied-material certificate here.
use super::*;
use crate::{
    authored_hanger_bolt::{Boundary, Component, CurveKey, End, FaceKey, Knot},
    authored_internal_treatment::ThreadBandRole,
    authored_metric_thread::UnitInterval,
};
#[derive(Clone, Copy)]
pub struct DomainOptions {
    pub maximum_work: usize,
    pub maximum_bytes: usize,
}
impl Default for DomainOptions {
    fn default() -> Self {
        Self {
            maximum_work: 16_000_000,
            maximum_bytes: 64 * 1024 * 1024,
        }
    }
}
pub struct RootFaceDomain<'a> {
    boundary: &'a Boundary,
    face: FaceKey,
    period: u32,
    turns: [f64; 2],
    offset: I,
    work: usize,
    bytes: usize,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TrimBounds {
    pub point_uv: [[f64; 2]; 2],
    pub first_derivative: [[f64; 2]; 2],
    pub second_derivative: [[f64; 2]; 2],
    /// Original Boundary curve parameter. Root cap curves use original
    /// LINEAR-PHASE u, while their rational surface PCurve uses root t.
    pub original_parameter: [f64; 2],
    pub original_parameter_first: [f64; 2],
    pub original_parameter_second: [f64; 2],
}
fn pair(v: I) -> [f64; 2] {
    [v.lo, v.hi]
}
enum TrimKind {
    Cap(End),
    Rail(bool),
}
fn work<'a>(
    options: DomainOptions,
    bytes: usize,
    interrupt: Option<&'a dyn Fn() -> bool>,
) -> Result<Work<'a>, Error> {
    let mut w = Work::new(options.maximum_work, interrupt)?;
    w.charge(1)?;
    if options.maximum_bytes == 0 || options.maximum_bytes > 64 * 1024 * 1024 {
        return Err(Error::Invalid);
    }
    let mut memory = Memory::new(bytes, options.maximum_bytes)?;
    memory.add(SCRATCH)?;
    Ok(w)
}
impl<'a> RootFaceDomain<'a> {
    pub fn prepare(
        boundary: &'a Boundary,
        face: FaceKey,
        options: DomainOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, Error> {
        let bytes = boundary.retained_bytes() + size_of::<Self>();
        let mut w = work(options, bytes, interrupt)?;
        if !matches!(
            face,
            FaceKey::MetricThread {
                band: ThreadBandRole::Root,
                ..
            }
        ) {
            return Err(Error::Invalid);
        }
        let mut period = None;
        for curve in boundary.curves() {
            w.charge(1)?;
            if curve.faces.contains(&face)
                && let CurveKey::Cap {
                    component: Component::Metric,
                    end: End::Start,
                    period: actual,
                    band: ThreadBandRole::Root,
                } = curve.key
            {
                period = Some(actual);
                break;
            }
        }
        let period = period.ok_or(Error::Invalid)?;
        w.charge(512)?;
        let s = boundary.source().metric();
        let spec = s.spec();
        let lead = p(spec.pitch_mm) * p(f64::from(spec.starts));
        let index = p(spec.pitch_mm) * p(f64::from(period));
        let root = s.root_meridian().homogeneous_controls();
        let q0 = hull(root[0][0]);
        let q1 = hull(root[2][0]);
        // Outward cover, not a new rounded authored station. The cap maps
        // below retain the original symbolic p,L and actual source period.
        let turns = [
            ((-q1 - index) / lead).lo,
            ((p(spec.length_mm) - q0 - index) / lead).hi,
        ];
        let offset = index - p(spec.length_mm) - p(boundary.source().spec().plain_length_mm);
        w.charge(1)?;
        Ok(Self {
            boundary,
            face,
            period,
            turns,
            offset,
            work: w.used,
            bytes,
        })
    }
    pub fn source_boundary(&self) -> &'a Boundary {
        self.boundary
    }
    pub fn face_key(&self) -> FaceKey {
        self.face
    }
    pub fn phase_period(&self) -> u32 {
        self.period
    }
    pub fn physical_turns(&self) -> [f64; 2] {
        self.turns
    }
    /// Add this exact-source enclosure to principal RootSurface's raw z.
    pub fn axial_offset_mm(&self) -> [f64; 2] {
        pair(self.offset)
    }
    /// q'>0 and radial(Pu×Pv)=−q'*(signed tau)*r prove this orientation.
    pub fn same_sense(&self) -> bool {
        self.boundary.source().metric().spec().left_hand
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_bytes(&self) -> usize {
        self.bytes
    }
    /// Parameter bounds for an ACTUAL edge of this source face. Cap t is the
    /// rational meridian parameter; rail t retains original axial parameter.
    /// This does not certify independently rounded/encoded PCurve controls.
    pub fn trim_bounds(
        &self,
        key: CurveKey,
        parameter: [f64; 2],
        options: DomainOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<TrimBounds, Error> {
        let mut w = work(options, self.bytes, interrupt)?;
        let t = UnitInterval::checked(parameter).ok_or(Error::Invalid)?;
        let kind = match key {
            CurveKey::Cap { end, .. } => TrimKind::Cap(end),
            CurveKey::Rail {
                knot: Knot::UpperRoot,
                ..
            } => TrimKind::Rail(false),
            CurveKey::Rail {
                knot: Knot::LowerRoot,
                ..
            } => TrimKind::Rail(true),
            _ => return Err(Error::Invalid),
        };
        let mut found = false;
        for c in self.boundary.curves() {
            w.charge(1)?;
            if c.key == key && c.faces.contains(&self.face) {
                found = true;
                break;
            }
        }
        if !found {
            return Err(Error::Invalid);
        }
        w.charge(1024)?;
        let s = self.boundary.source().metric();
        let spec = s.spec();
        let lead = p(spec.pitch_mm) * p(f64::from(spec.starts));
        let index = p(spec.pitch_mm) * p(f64::from(self.period));
        let span = p(self.turns[1]) - p(self.turns[0]);
        let divisor = lead * span;
        let result = match kind {
            TrimKind::Cap(end) => {
                let root = s.root_meridian().in_domain(t);
                let q = hull(root.point_mm[0]);
                let dq = hull(root.first_derivative_mm[0]);
                let ddq = hull(root.second_derivative_mm[0]);
                let z = p(if end == End::Start {
                    0.
                } else {
                    spec.length_mm
                });
                let v = (z - q - index - lead * p(self.turns[0])) / divisor;
                let original = (q - p(spec.pitch_mm) * p(3.) / p(8.)) / (p(spec.pitch_mm) / p(4.));
                TrimBounds {
                    point_uv: [pair(t.interval()), pair(v)],
                    first_derivative: [pair(p(1.)), pair(-dq / divisor)],
                    second_derivative: [pair(p(0.)), pair(-ddq / divisor)],
                    // Exact monotone endpoint/range theorem earns this
                    // intersection; no near-event geometry is merged.
                    original_parameter: [original.lo.max(0.), original.hi.min(1.)],
                    original_parameter_first: pair(p(4.) * dq / p(spec.pitch_mm)),
                    original_parameter_second: pair(p(4.) * ddq / p(spec.pitch_mm)),
                }
            }
            TrimKind::Rail(end) => {
                // Membership in the sealed Root face proves these are its
                // two original root rails, never a guessed crest endpoint.
                let u = if end { 1. } else { 0. };
                let root = s.root_meridian().in_domain(UnitInterval::endpoint(end));
                let q = hull(root.point_mm[0]);
                let z = p(spec.length_mm) * t.interval();
                let v = (z - q - index - lead * p(self.turns[0])) / divisor;
                TrimBounds {
                    point_uv: [pair(p(u)), pair(v)],
                    first_derivative: [pair(p(0.)), pair(p(spec.length_mm) / divisor)],
                    second_derivative: [pair(p(0.)); 2],
                    original_parameter: pair(t.interval()),
                    original_parameter_first: pair(p(1.)),
                    original_parameter_second: pair(p(0.)),
                }
            }
        };
        // Keep UV and original-parameter enclosures finite; validated source
        // dimensions and positive lead/span bound all fixed operations above.
        // The true UV image lies in the closed rectangle by construction;
        // outward arithmetic may extend slightly beyond its boundary.
        w.charge(1)?;
        Ok(result)
    }
}
#[cfg(test)]
mod tests;
