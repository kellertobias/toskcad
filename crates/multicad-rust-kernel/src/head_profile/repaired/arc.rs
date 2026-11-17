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
    // Original proposal constants for a later explicit correspondence. They
    // do not define or replace the repaired source expression.
    proposal_angles: [f64; 2],
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
    /// Whole normalized-ray parameter speed and smoothstep repair bounds.
    /// The base ray is q(t)=((1-t)a+tb)/|...|, not uniform angle.
    pub(super) fn raw_support_envelope(self) -> [f64; 3] {
        let n2 = (I::point(1.) + I::point(self.dot.lo)) / I::point(2.);
        let speed = I::point(self.det.magnitude()) / I::point(n2.lo);
        let correction = norm(self.delta[0])
            .unwrap()
            .max(norm(self.delta[1]).unwrap());
        let change = norm(std::array::from_fn(|i| self.delta[1][i] - self.delta[0][i])).unwrap();
        [speed.hi, correction, (I::point(1.5) * I::point(change)).hi]
    }
    pub(super) fn proposal_transfer(self) -> Result<[f64; 3], RepairedHeadError> {
        // The proposal's cached atan2 is a defining binary64 constant, not an
        // estimate of the repaired source angle. Enclose its actual circle rays.
        let [start, sweep] = self.proposal_angles;
        let angles = [I::point(start), I::point(start) + I::point(sweep)];
        let old = [
            super::proposal::rotation(angles[0])?,
            super::proposal::rotation(angles[1])?,
        ];
        let errors =
            [(self.a, old[0]), (self.b, old[1])].map(|(a, b)| norm([a[0] - b[0], a[1] - b[1]]));
        let [e0, e1] = [errors[0]?, errors[1]?];
        // For chord e<1/16, minor angular error 2 asin(e/2)<2e. Both
        // source sectors have the same strict orientation and sweep<3<pi;
        // endpoint rays this close therefore select the same unwrapped sector.
        ensure(e0 < 1. / 16. && e1 < 1. / 16. && sweep.abs() < 3.)?;
        let angular = I::point(2.) * I::point(e0.max(e1));
        let det_min = self.det.lo.abs().min(self.det.hi.abs());
        let delta_angle = I::point(2.) * (I::point(e0) + I::point(e1));
        let ratio = delta_angle / I::point(det_min);
        ensure(finite(ratio) && ratio.hi < 0.25)?;
        // Actual minor span>=|sin(span)|=|det|. The angular correspondence
        // derivative differs from identity by <=sum(endpoint errors)/|det|.
        let n2_min = ((I::point(1.) + I::point(self.dot.lo)) / I::point(2.)).lo;
        let speed = I::point(self.det.magnitude()) / I::point(n2_min);
        let correction = norm(self.delta[0])?.max(norm(self.delta[1])?);
        let change = norm(std::array::from_fn(|i| self.delta[1][i] - self.delta[0][i]))?;
        let position = I::point(self.radius) * angular + I::point(correction);
        let derivative =
            I::point(self.radius) * speed * (angular + ratio) + I::point(1.5) * I::point(change);
        let old_min = I::point(self.radius) * I::point(det_min) * (I::point(1.) - ratio);
        let normal_ratio = derivative / old_min;
        ensure(finite(normal_ratio) && normal_ratio.hi < 0.25)?;
        // |a-b|/|b| bounds sin(angle) after perturbation; the conservative
        // 2q/(1-q) bound also includes normalization and is valid for q<1/4.
        let normal = I::point(2.) * normal_ratio / (I::point(1.) - normal_ratio);
        ensure([position, derivative, normal].into_iter().all(finite))?;
        Ok([position.hi, derivative.hi, normal.hi])
    }
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
                proposal_angles: [
                    (s.start[1] - center[1]).atan2(s.start[0] - center[0]),
                    sweep,
                ],
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

impl ArcSource {
    pub(super) fn correspondence(
        self,
        other: Self,
        z_shift: I,
    ) -> Result<[I; 2], RepairedHeadError> {
        let mag = |v: [I; 2]| I::point(v[0].magnitude()) + I::point(v[1].magnitude());
        let da = std::array::from_fn(|i| self.a[i] - other.a[i]);
        let db = std::array::from_fn(|i| self.b[i] - other.b[i]);
        let dw = I::point(mag(da).hi.max(mag(db).hi));
        let dv = mag(std::array::from_fn(|i| db[i] - da[i]));
        let lower =
            ((I::point(1.) + I::point(self.dot.lo.min(other.dot.lo))) / I::point(2.)).sqrt();
        ensure(finite(lower) && lower.lo > 0.)?;
        let lower = I::point(lower.lo);
        // Normalizing either chord obeys ||qa-qb||<=2||wa-wb||/nmin.
        let direction = I::point(2.) * dw / lower;
        let va = mag(std::array::from_fn(|i| self.b[i] - self.a[i]));
        let vb = mag(std::array::from_fn(|i| other.b[i] - other.a[i]));
        // q'=(I-qq^T)v/|w|. Unit-vector projectors have norm<=1;
        // their difference has norm<=2||qa-qb||, and inverse norms
        // differ by at most ||wa-wb||/nmin².
        let derivative =
            dv / lower + vb * dw / lower.square() + I::point(2.) * direction * vb / lower;
        let radius_delta = I::point((I::point(self.radius) - I::point(other.radius)).magnitude());
        let center = mag([
            I::point(self.center[0]) - I::point(other.center[0]),
            I::point(self.center[1]) - I::point(other.center[1]) - z_shift,
        ]);
        let corrections = [0, 1].map(|end| {
            mag(std::array::from_fn(|i| {
                self.delta[end][i] - other.delta[end][i]
            }))
        });
        let correction_change = mag(std::array::from_fn(|i| {
            (self.delta[1][i] - other.delta[1][i]) - (self.delta[0][i] - other.delta[0][i])
        }));
        let position = center
            + radius_delta
            + I::point(other.radius) * direction
            + I::point(corrections[0].hi.max(corrections[1].hi));
        let derivative = radius_delta * va / lower
            + I::point(other.radius) * derivative
            + I::point(1.5) * correction_change;
        ensure([position, derivative].into_iter().all(finite))?;
        Ok([position, derivative])
    }
}
