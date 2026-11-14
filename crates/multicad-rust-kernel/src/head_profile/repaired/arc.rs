use super::*;
#[derive(Clone, Copy)]
pub(super) struct ArcSource {
    center: [f64; 2],
    radius: f64,
    a: [I; 2],
    b: [I; 2],
    delta: [[I; 2]; 2],
    dot: I,
    det: I,
}
fn ray(p: [f64; 2], center: [f64; 2]) -> Result<[I; 2], RepairedHeadError> {
    ensure(p != center)?;
    if p[0] == center[0] {
        return Ok([I::point(0.), I::point((p[1] - center[1]).signum())]);
    }
    if p[1] == center[1] {
        return Ok([I::point((p[0] - center[0]).signum()), I::point(0.)]);
    }
    let v = std::array::from_fn::<_, 2, _>(|i| I::point(p[i]) - I::point(center[i]));
    let n = (v[0].square() + v[1].square()).sqrt();
    ensure(finite(n) && n.lo > 0.)?;
    let r = v.map(|v| v / n);
    ensure(r.into_iter().all(finite))?;
    Ok(r)
}
fn norm(v: [I; 2]) -> Result<f64, RepairedHeadError> {
    let n = (v[0].square() + v[1].square()).sqrt();
    ensure(finite(n))?;
    Ok(n.hi)
}
impl ArcSource {
    pub(super) fn new(
        s: Segment,
        center: [f64; 2],
        radius: f64,
        sweep: f64,
    ) -> Result<(Self, RepairedSegmentCertificate), RepairedHeadError> {
        let a = ray(s.start, center)?;
        let b = ray(s.end, center)?;
        let det = a[0] * b[1] - a[1] * b[0];
        let sign = if sweep > 0. { 1. } else { -1. };
        ensure(finite(det) && (det * I::point(sign)).lo > 0.)?;
        // The actual endpoints must lie on the appropriate x-side for strict
        // axial order of the whole normalized-ray minor arc, not only samples.
        let ordered_side = if sign > 0. {
            s.start[0] >= center[0] && s.end[0] >= center[0]
        } else {
            s.start[0] <= center[0] && s.end[0] <= center[0]
        };
        ensure(ordered_side)?;
        let dot = a[0] * b[0] + a[1] * b[1];
        ensure(finite(dot) && dot.lo > -1.)?;
        let dot = I::new(dot.lo.max(-1.), dot.hi.min(1.));
        let delta = [(s.start, a), (s.end, b)].map(|(p, u)| {
            std::array::from_fn(|i| I::point(p[i]) - I::point(center[i]) - I::point(radius) * u[i])
        });
        let dz = (delta[1][1] - delta[0][1]).magnitude();
        let det_lower = if sign > 0. { det.lo } else { -det.hi };
        let (progress, coefficient) = if s.start[0] == center[0] {
            (
                AxialProgress::ZeroAtStart,
                (I::point(radius) * det * b[0] - I::point(6.) * I::point(dz)).lo,
            )
        } else {
            let x = if sign > 0. {
                a[0].lo.min(b[0].lo)
            } else {
                (-a[0].hi).min(-b[0].hi)
            };
            (
                AxialProgress::Strict,
                (I::point(radius) * I::point(det_lower) * I::point(x)
                    - I::point(1.5) * I::point(dz))
                .lo,
            )
        };
        ensure(coefficient > 0. && coefficient.is_finite())?;
        // |w|<=1 for convex combinations of exact unit rays. On the positive
        // x-side w_x/|w|>=min(endpoint x); on the negative side >=-1.
        let raw_r = if sign > 0. {
            I::point(center[0]) + I::point(radius) * I::point(a[0].lo.min(b[0].lo))
        } else {
            I::point(center[0]) - I::point(radius)
        };
        let minimum_radius = (raw_r + I::point(delta[0][0].lo.min(delta[1][0].lo))).lo;
        ensure(minimum_radius > 0. && minimum_radius.is_finite())?;
        let transfer = norm(delta[0])?.max(norm(delta[1])?);
        // Same circle, same orientation and minor sector. A monotone angular
        // correspondence bounds whole finite-arc set transfer by the maximum
        // endpoint angular discrepancy; normalized t need not equal old t.
        let old_start = super::super::source_volume::angle(s, center);
        let old_end = super::super::source_volume::angle(Segment { start: s.end, ..s }, center);
        let start_width = (I::point(old_start.hi) - I::point(old_start.lo)).hi;
        let end_error = (old_start + I::point(sweep) - old_end).magnitude();
        let trim = (I::point(radius) * I::point(start_width.max(end_error))).hi;
        ensure(trim.is_finite())?;
        Ok((
            Self {
                center,
                radius,
                a,
                b,
                delta,
                dot,
                det,
            },
            RepairedSegmentCertificate {
                role: s.role,
                stations: [0, 0],
                progress,
                minimum_progress_mm: coefficient,
                minimum_radius_mm: minimum_radius,
                support_transfer_mm: transfer,
                finite_trim_transfer_mm: trim,
            },
        ))
    }
    pub(super) fn evaluate(self, t: I) -> [[I; 2]; 2] {
        let w = std::array::from_fn::<_, 2, _>(|i| self.a[i] + t * (self.b[i] - self.a[i]));
        // Exact unit-ray identity preserves the positive denominator when an
        // AABB of w contains zero. t(1-t) lies in[0,1/4] on the unit domain.
        let q = I::point(0.25) - (t - I::point(0.5)).square();
        let q = I::new(q.lo.max(0.), q.hi.min(0.25));
        let s2 = I::point(1.) - I::point(2.) * (I::point(1.) - self.dot) * q;
        let lower = ((I::point(1.) + I::point(self.dot.lo)) / I::point(2.)).lo;
        let n = I::new(s2.lo.max(lower), s2.hi.min(1.)).sqrt();
        let h = t.square() * (I::point(3.) - I::point(2.) * t);
        let h = I::new(h.lo.max(0.), h.hi.min(1.));
        let change = std::array::from_fn::<_, 2, _>(|i| self.delta[1][i] - self.delta[0][i]);
        let p = std::array::from_fn(|i| {
            I::point(self.center[i])
                + I::point(self.radius) * w[i] / n
                + self.delta[0][i]
                + change[i] * h
        });
        let perpendicular = [-w[1], w[0]];
        let derivative = std::array::from_fn(|i| {
            I::point(self.radius) * self.det * perpendicular[i] / (n * n * n)
                + I::point(6.) * q * change[i]
        });
        [p, derivative]
    }
}

impl ArcSource {
    /// Fifth geometry jet of the SAME normalized-ray+cubic source. Only c0
    /// intersections use the already proved unit-domain identities; derivative
    /// coefficients are retained, never clamped or replaced by samples.
    pub(super) fn position_jet(
        self,
        t: crate::tapered_thread::Series<6>,
    ) -> [crate::tapered_thread::Series<6>; 2] {
        use crate::tapered_thread::{Scalar, Series};
        let c = Series::<6>::interval;
        let w = std::array::from_fn::<_, 2, _>(|i| c(self.a[i]) + t * c(self.b[i] - self.a[i]));
        let u = t - c(I::point(0.5));
        let mut q = c(I::point(0.25)) - u * u;
        q.c[0] = I::new(q.c[0].lo.max(0.), q.c[0].hi.min(0.25));
        let mut s2 = c(I::point(1.)) - c(I::point(2.) * (I::point(1.) - self.dot)) * q;
        let lower = ((I::point(1.) + I::point(self.dot.lo)) / I::point(2.)).lo;
        s2.c[0] = I::new(s2.c[0].lo.max(lower), s2.c[0].hi.min(1.));
        let n = s2.sqrt();
        let mut h = t * t * (c(I::point(3.)) - c(I::point(2.)) * t);
        h.c[0] = I::new(h.c[0].lo.max(0.), h.c[0].hi.min(1.));
        std::array::from_fn(|i| {
            c(I::point(self.center[i]))
                + c(I::point(self.radius)) * w[i] / n
                + c(self.delta[0][i])
                + c(self.delta[1][i] - self.delta[0][i]) * h
        })
    }
}
