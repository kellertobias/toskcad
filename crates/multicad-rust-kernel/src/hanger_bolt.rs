//! Connected flat-wood/rounded-metric source with literal finite end laws.
//! No compound, mesh, receiving cavity or exchange representation is implied.
use crate::{
    authored_flat_thread::AuthoredThread, helical_thread::Thread, tapered_thread::interval::I,
};
use num_rational::BigRational;
mod classify;
#[cfg(test)]
mod tests;
mod volume;
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Error {
    Invalid,
    Unsupported,
    Unresolved,
    Limit,
    Cancelled,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Spec {
    pub wood: AuthoredThread,
    pub metric: Thread,
    pub plain_length_mm: f64,
    pub plain_radius_mm: f64,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SectionKind {
    Metric,
    Plain,
    Wood,
}
/// Endpoints enclose the exact authored sums, not substituted rounded stations.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AxialSection {
    pub kind: SectionKind,
    pub start_mm: [f64; 2],
    pub end_mm: [f64; 2],
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Hanger {
    spec: Spec,
    metric_bottom: [f64; 2],
    core_radius_mm: f64,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Classification {
    Inside,
    Outside,
    UnresolvedBoundary,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct QueryOptions {
    pub volume_error_mm3: f64,
    pub maximum_cells: usize,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Volume {
    pub material_mm3: [f64; 2],
    pub wood_mm3: [f64; 2],
    pub plain_mm3: [f64; 2],
    pub metric_mm3: [f64; 2],
    /// Already included in metric_mm3, not a caller-selected source tolerance.
    pub metric_phase_error_mm3: f64,
}
impl Hanger {
    pub fn new(spec: Spec, interrupt: Option<&dyn Fn() -> bool>) -> Result<Self, Error> {
        checkpoint(interrupt)?;
        // This connected graph still owns binary finite end stations. An
        // explicit canonical domain needs its own shoulder/cap qualification.
        if spec.wood.exact_axial_extent().is_some() {
            return Err(Error::Unsupported);
        }
        if [spec.plain_length_mm, spec.plain_radius_mm]
            .into_iter()
            .any(|v| !v.is_finite() || v <= 0. || v > 1e7)
        {
            return Err(Error::Invalid);
        }
        if spec.metric.spec().root_radius_mm.is_none() {
            return Err(Error::Unsupported);
        }
        let radial = metric_radius(spec.metric)?
            .positive()
            .map_err(|_| Error::Unresolved)?;
        let q = |x| BigRational::from_float(x).expect("sealed finite source support");
        // Strict complete-radius dominance gives two genuine annular shoulders.
        // Crossing/equal radii require their own finite trim graph, not fitting.
        if q(spec.plain_radius_mm) <= q(radial.hi)
            || q(spec.plain_radius_mm) >= q(spec.wood.minimum_radius_mm())
        {
            return Err(Error::Unsupported);
        }
        let bottom = (-I::point(spec.plain_length_mm) - I::point(spec.metric.spec().length_mm))
            .checked()
            .map_err(|_| Error::Unresolved)?;
        let core_radius_mm = radial
            .lo
            .min(spec.plain_radius_mm)
            .min(spec.wood.minimum_radius_mm());
        let result = Self {
            spec,
            metric_bottom: pair(bottom),
            core_radius_mm,
        };
        checkpoint(interrupt)?;
        Ok(result)
    }
    pub fn spec(&self) -> Spec {
        self.spec
    }
    pub fn wood(&self) -> &AuthoredThread {
        &self.spec.wood
    }
    pub fn metric(&self) -> &Thread {
        &self.spec.metric
    }
    pub fn connection_core_radius_mm(&self) -> f64 {
        self.core_radius_mm
    }
    pub fn sections(&self) -> [AxialSection; 3] {
        let p = -self.spec.plain_length_mm;
        [
            AxialSection {
                kind: SectionKind::Metric,
                start_mm: self.metric_bottom,
                end_mm: [p, p],
            },
            AxialSection {
                kind: SectionKind::Plain,
                start_mm: [p, p],
                end_mm: [0., 0.],
            },
            AxialSection {
                kind: SectionKind::Wood,
                start_mm: [0., 0.],
                end_mm: [self.wood().full_length_mm(); 2],
            },
        ]
    }
    /// Complete source envelope for separation only, never material membership.
    pub fn bounds_mm(&self) -> [[f64; 2]; 3] {
        let r = self.wood().maximum_radius_mm();
        [
            [-r, r],
            [-r, r],
            [self.metric_bottom[0], self.wood().full_length_mm()],
        ]
    }
}
fn metric_radius(thread: Thread) -> Result<I, Error> {
    let b = thread.profile_breaks();
    let mut result = I::new(f64::INFINITY, f64::NEG_INFINITY);
    for band in 0..5 {
        let phase = I::new(b[band], b[band + 1]);
        let r = metric_profile(thread, band, phase)
            .checked()
            .map_err(|_| Error::Unresolved)?;
        result.lo = result.lo.min(r.lo);
        result.hi = result.hi.max(r.hi);
    }
    Ok(result)
}
fn metric_profile(thread: Thread, band: usize, phase: I) -> I {
    let b = thread.profile_breaks();
    let major = I::point(thread.spec().major_radius_mm);
    match band {
        1 | 3 => {
            let base = if band == 1 { b[1] } else { b[4] };
            major + I::point(thread.slope_at_phase(base, band)) * (phase - I::point(base))
        }
        2 => {
            let root = thread.root_circle().expect("sealed rounded metric source");
            I::point(root.center_radius)
                - (I::point(root.radius).square() - (phase - I::point(root.phase)).square()).sqrt()
        }
        _ => major,
    }
}
fn checkpoint(interrupt: Option<&dyn Fn() -> bool>) -> Result<(), Error> {
    if interrupt.is_some_and(|f| f()) {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
fn pair(i: I) -> [f64; 2] {
    [i.lo, i.hi]
}
fn pi() -> I {
    I::new(
        std::f64::consts::PI.next_down(),
        std::f64::consts::PI.next_up(),
    )
}
