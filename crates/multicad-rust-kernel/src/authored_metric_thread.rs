//! Continuous authored60° metric profile with canonical6g rounded-root law, distinct from literal legacy Thread.
//! One shared Q(sqrt(3)) expression owns each tangent point. Outward evaluated
//! descriptors never become independent physical line/circle inputs.
use crate::tapered_thread::interval::I;
mod classify;
#[cfg(test)]
mod tests;
mod transfer;
pub use transfer::{BandTransfer, TransferBounds, TransferOptions};
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Error {
    Invalid,
    Limit,
    Unresolved,
    Cancelled,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Spec {
    pub major_radius_mm: f64,
    pub pitch_mm: f64,
    pub length_mm: f64,
    pub crest_width_mm: f64,
    pub starts: u32,
    pub left_hand: bool,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Knot {
    Begin,
    CrestEnd,
    UpperRoot,
    LowerRoot,
    CrestStart,
    End,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PointBounds {
    pub phase_mm: [f64; 2],
    pub radius_mm: [f64; 2],
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ProfileBounds {
    pub phase_mm: [f64; 2],
    pub radius_mm: [f64; 2],
    pub slope: [f64; 2],
    pub second_derivative_per_mm: [f64; 2],
    /// Positive homogeneous radial coefficient of the meridian normal.
    /// On the circle this is H=sqrt(1-X²), not an independent support input.
    pub normal_radial: [f64; 2],
    pub normal_axial: [f64; 2],
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AuthoredMetricThread {
    spec: Spec,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Classification {
    Inside,
    Outside,
    UnresolvedBoundary,
}
impl AuthoredMetricThread {
    pub fn new(spec: Spec) -> Result<Self, Error> {
        if [
            spec.major_radius_mm,
            spec.pitch_mm,
            spec.length_mm,
            spec.crest_width_mm,
        ]
        .into_iter()
        .any(|v| !v.is_finite() || v <= 0.)
            || spec.crest_width_mm >= spec.pitch_mm
        {
            return Err(Error::Invalid);
        }
        if [spec.major_radius_mm, spec.pitch_mm, spec.length_mm]
            .into_iter()
            .any(|v| v > 1e7)
            || !(1..=8).contains(&spec.starts)
            || spec.pitch_mm < 0.001
            || spec.length_mm < spec.pitch_mm
            || spec.length_mm / spec.pitch_mm + 2. * f64::from(spec.starts) > 128.
        {
            return Err(Error::Limit);
        }
        let value = Self { spec };
        // Strict evaluated separation certifies the exact shared expressions. The
        // caller's accuracy never moves a tangent, radius or crest.
        let q = value.knots();
        for pair in q.windows(2) {
            positive(pair[1][0] - pair[0][0])?;
        }
        positive(value.center_radius() - value.root_radius())?;
        Ok(value)
    }
    pub fn spec(&self) -> Spec {
        self.spec
    }
    pub fn full_length_mm(&self) -> f64 {
        self.spec.length_mm
    }
    pub fn minimum_radius_mm(&self) -> [f64; 2] {
        pair(self.center_radius() - self.root_radius())
    }
    /// Exact canonical engineering relation R=p/(4sqrt(3)), enclosed rather
    /// than retained as an independently rounded physical radius.
    pub fn root_radius_mm(&self) -> [f64; 2] {
        pair(self.root_radius())
    }
    fn root_radius(&self) -> I {
        I::point(self.spec.pitch_mm) * sqrt3() / I::point(12.)
    }
    pub fn maximum_radius_mm(&self) -> f64 {
        self.spec.major_radius_mm
    }
    pub fn knot(&self, knot: Knot) -> PointBounds {
        let index = match knot {
            Knot::Begin => 0,
            Knot::CrestEnd => 1,
            Knot::UpperRoot => 2,
            Knot::LowerRoot => 3,
            Knot::CrestStart => 4,
            Knot::End => 5,
        };
        let [q, r] = self.knots()[index];
        PointBounds {
            phase_mm: pair(q),
            radius_mm: pair(r),
        }
    }
    pub fn profile_bounds(&self, band: usize, u: [f64; 2]) -> Result<ProfileBounds, Error> {
        let u = I::new(u[0], u[1]);
        if band >= 5 || !u.valid() || u.lo < 0. || u.hi > 1. {
            return Err(Error::Invalid);
        }
        let k = self.knots();
        let phase = k[band][0] + (k[band + 1][0] - k[band][0]) * u;
        let root = self.root_radius();
        let r3 = sqrt3();
        let (radius, slope, second, radial, axial) = match band {
            1 => (
                k[1][1] + (k[2][1] - k[1][1]) * u,
                -r3,
                I::point(0.),
                I::point(1.),
                r3,
            ),
            3 => (
                k[3][1] + (k[4][1] - k[3][1]) * u,
                r3,
                I::point(0.),
                I::point(1.),
                -r3,
            ),
            2 => {
                // x/R=sqrt(3)*(u−1/2); h/R>=1/2 over the entire true trim.
                // Retain this normalized form to avoid subtracting large phase centers.
                let x = r3 * (u - I::point(0.5));
                let h = (I::point(1.) - x.square()).sqrt();
                let h = positive(h)?;
                (
                    self.center_radius() - root * h,
                    x / h,
                    I::point(1.) / (root * h * h * h),
                    h,
                    -x,
                )
            }
            _ => (
                I::point(self.spec.major_radius_mm),
                I::point(0.),
                I::point(0.),
                I::point(1.),
                I::point(0.),
            ),
        };
        let radius = if u.lo == 0. && u.hi == 0. {
            k[band][1]
        } else if u.lo == 1. && u.hi == 1. {
            k[band + 1][1]
        } else {
            radius
        };
        for value in [phase, radius, slope, second, radial, axial] {
            value.checked().map_err(|_| Error::Unresolved)?;
        }
        Ok(ProfileBounds {
            phase_mm: pair(phase),
            radius_mm: pair(radius),
            slope: pair(slope),
            second_derivative_per_mm: pair(second),
            normal_radial: pair(radial),
            normal_axial: pair(axial),
        })
    }
    /// Exact source envelope. It is not a replacement for material membership.
    pub fn bounds_mm(&self) -> [[f64; 2]; 3] {
        let r = self.spec.major_radius_mm;
        [[-r, r], [-r, r], [0., self.spec.length_mm]]
    }
    fn center_radius(&self) -> I {
        I::point(self.spec.major_radius_mm)
            - (I::point(self.spec.pitch_mm) - I::point(self.spec.crest_width_mm)) * sqrt3()
                / I::point(2.)
            + I::point(2.) * self.root_radius()
    }
    fn knots(&self) -> [[I; 2]; 6] {
        let root = self.root_radius();
        let m = I::point(self.spec.major_radius_mm);
        let tangent_r = self.center_radius() - root / I::point(2.);
        let points = [
            Knot::Begin,
            Knot::CrestEnd,
            Knot::UpperRoot,
            Knot::LowerRoot,
            Knot::CrestStart,
            Knot::End,
        ];
        points.map(|knot| {
            let [(pn, pd), (cn, cd)] = phase_coefficients(knot);
            let q = I::point(self.spec.pitch_mm) * I::point(pn as f64) / I::point(pd as f64)
                + I::point(self.spec.crest_width_mm) * I::point(cn as f64) / I::point(cd as f64);
            [
                q,
                if matches!(knot, Knot::UpperRoot | Knot::LowerRoot) {
                    tangent_r
                } else {
                    m
                },
            ]
        })
    }
    fn rate(&self) -> I {
        let lead = I::point(self.spec.pitch_mm) * I::point(f64::from(self.spec.starts));
        let rate = lead / (I::point(2.) * pi());
        if self.spec.left_hand { -rate } else { rate }
    }
}
fn sqrt3() -> I {
    I::point(3.).sqrt()
}
fn pi() -> I {
    I::new(
        std::f64::consts::PI.next_down(),
        std::f64::consts::PI.next_up(),
    )
}
fn positive(v: I) -> Result<I, Error> {
    v.positive().map_err(|_| Error::Unresolved)
}
fn pair(v: I) -> [f64; 2] {
    [v.lo, v.hi]
}
fn interval(v: [f64; 2]) -> I {
    I::new(v[0], v[1])
}
fn checkpoint(f: Option<&dyn Fn() -> bool>) -> Result<(), Error> {
    if f.is_some_and(|f| f()) {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}

/// The actual owned phase expression q=(pn/pd)*pitch+(cn/cd)*crest.
/// Exact clipping evaluates these same rational expressions; it must never
/// infer source vertices from the midpoint of a point enclosure.
pub(crate) fn phase_coefficients(knot: Knot) -> [(i64, i64); 2] {
    match knot {
        Knot::Begin => [(0, 1), (0, 1)],
        Knot::CrestEnd => [(0, 1), (1, 2)],
        Knot::UpperRoot => [(3, 8), (0, 1)],
        Knot::LowerRoot => [(5, 8), (0, 1)],
        Knot::CrestStart => [(1, 1), (-1, 2)],
        Knot::End => [(1, 1), (0, 1)],
    }
}

// The geometry remains this actual source. The detached open shell is an
// approximation receipt, not a closed manufacturing or exchange authority.
pub use crate::authored_internal_treatment::{MetricRimEdge, MetricShell};
impl AuthoredMetricThread {
    pub fn tessellate_shell(
        &self,
        options: crate::display::Tessellation,
        maximum_work: usize,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<MetricShell, Error> {
        crate::authored_internal_treatment::metric_shell(self, options, maximum_work, interrupt)
            .map_err(shell_error)
    }
}

fn shell_error(e: crate::authored_internal_treatment::Error) -> Error {
    match e {
        crate::authored_internal_treatment::Error::Invalid => Error::Invalid,
        crate::authored_internal_treatment::Error::Unresolved => Error::Unresolved,
        crate::authored_internal_treatment::Error::Limit => Error::Limit,
        crate::authored_internal_treatment::Error::Cancelled => Error::Cancelled,
    }
}
