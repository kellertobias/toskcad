//! Whole derived-face distance to actual raw analytic SUPPORTS, not face fills.
use super::*;
use crate::head_profile::{Curve, Segment};
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RawHeadSupportKind {
    Plane,
    Cylinder,
    Cone,
    Sphere,
    Torus,
}
#[derive(Clone, Copy)]
pub(super) enum Kind {
    Plane,
    Cylinder(I),
    Cone(I, I),
    Sphere(I),
    Torus(I, I),
}
pub(super) struct RawSupport {
    pub kind: Kind,
    pub frame: support::Frame,
}
fn valid(x: I) -> bool {
    x.valid() && x.lo.is_finite() && x.hi.is_finite()
}
fn norm<const N: usize>(v: [I; N]) -> Result<f64> {
    let n = v
        .into_iter()
        .fold(I::point(0.), |a, v| a + v.square())
        .sqrt();
    require(valid(n))?;
    Ok(n.hi)
}
fn dot(a: [I; 3], b: [I; 3]) -> I {
    (0..3).fold(I::point(0.), |v, i| v + a[i] * b[i])
}
fn vector(f: &support::Frame, p: Point) -> [I; 3] {
    std::array::from_fn(|i| dot(f.axes[i], p.map(I::point)))
}
fn origin(f: &support::Frame, c: CommonMouthChart, z: I) -> [I; 3] {
    let p = std::array::from_fn(|i| {
        I::point(c.columns[0][i]) + I::point(c.columns[3][i]) * (z - I::point(c.outer_top_z))
            - f.origin[i]
    });
    std::array::from_fn(|i| dot(f.axes[i], p))
}
pub(super) fn read(shell: &Shell<'_, '_>, id: u32) -> Result<RawSupport> {
    let e = shell.simple(id)?;
    let (axis, kind) = match (e.name, e.args.as_slice()) {
        ("PLANE", [Value::String(_), a]) => (reference(a)?, Kind::Plane),
        ("CYLINDRICAL_SURFACE", [Value::String(_), a, Value::Number(r)]) => {
            (reference(a)?, Kind::Cylinder(support::adjacent(*r)))
        }
        ("CONICAL_SURFACE", [Value::String(_), a, Value::Number(r), Value::Number(angle)]) => {
            let a0 = support::adjacent(*angle);
            let c = a0.cos();
            require(c.lo > 0.)?;
            (
                reference(a)?,
                Kind::Cone(support::adjacent(*r), a0.sin() / c),
            )
        }
        ("SPHERICAL_SURFACE", [Value::String(_), a, Value::Number(r)]) => {
            (reference(a)?, Kind::Sphere(support::adjacent(*r)))
        }
        ("TOROIDAL_SURFACE", [Value::String(_), a, Value::Number(r), Value::Number(t)]) => (
            reference(a)?,
            Kind::Torus(support::adjacent(*r), support::adjacent(*t)),
        ),
        _ => return Err(ReadError::Unsupported),
    };
    let radii = match kind {
        Kind::Plane => [I::point(1.); 2],
        Kind::Cylinder(r) | Kind::Cone(r, _) | Kind::Sphere(r) => [r; 2],
        Kind::Torus(r, t) => [r, t],
    };
    require(radii.into_iter().all(|r| valid(r) && r.lo > 0.))?;
    let frame = support::Frame::literal(&geometry::axis(shell, axis)?).map_err(source_error)?;
    Ok(RawSupport { kind, frame })
}
/// Chosen support-parameter basis only. It does not normalize or replace the
/// actual derived affine frame. Every finite phase and either sine orientation
/// parameterize the exact same raw circle/cylinder/cone support.
fn phase(u: [I; 3]) -> [I; 2] {
    let mut seed = [u[0].lo / 2. + u[0].hi / 2., u[1].lo / 2. + u[1].hi / 2.];
    if seed == [0.; 2] {
        seed = [1., 0.]
    }
    let n = (I::point(seed[0]).square() + I::point(seed[1]).square()).sqrt();
    seed.map(|v| I::point(v) / n)
}
fn polar_errors(u: [I; 3], v: [I; 3], e: [I; 2], sign: f64) -> Result<f64> {
    let a = norm([u[0] - e[0], u[1] - e[1]])?;
    let b = norm([v[0] + I::point(sign) * e[1], v[1] - I::point(sign) * e[0]])?;
    Ok((I::point(a) + I::point(b)).hi)
}
fn stretch(chart: CommonMouthChart) -> Result<[f64; 2]> {
    let mut error = 0_f64;
    for i in 0..3 {
        let mut row = I::point(0.);
        for j in 0..3 {
            let g = dot(
                chart.columns[i + 1].map(I::point),
                chart.columns[j + 1].map(I::point),
            ) - I::point(if i == j { 1. } else { 0. });
            row = row + I::point(g.magnitude());
        }
        error = error.max(row.hi);
    }
    require(error.is_finite() && error < 1.)?;
    let lower = (I::point(1.) - I::point(error)).sqrt();
    let upper = (I::point(1.) + I::point(error)).sqrt();
    require(valid(lower) && valid(upper) && lower.lo > 0.)?;
    Ok([lower.lo, upper.hi])
}
fn repair_world(sigma: f64, correction: f64) -> I {
    I::point(sigma) * I::point(correction)
}
fn plane(raw: &RawSupport, c: CommonMouthChart, bounds: [[f64; 2]; 2]) -> Result<f64> {
    let center = origin(&raw.frame, c, I::new(bounds[1][0], bounds[1][1]));
    let u = vector(&raw.frame, c.columns[1]);
    let v = vector(&raw.frame, c.columns[2]);
    let b = I::point(center[2].magnitude())
        + I::point(bounds[0][0].abs().max(bounds[0][1].abs()))
            * (I::point(u[2].magnitude()) + I::point(v[2].magnitude()));
    require(valid(b))?;
    Ok(b.hi)
}
fn tube(raw: &RawSupport, c: CommonMouthChart, s: Segment, r: I, k: I) -> Result<f64> {
    let u = vector(&raw.frame, c.columns[1]);
    let v = vector(&raw.frame, c.columns[2]);
    let e = phase(u);
    let a = origin(&raw.frame, c, I::point(s.start[1]));
    let b = origin(&raw.frame, c, I::point(s.end[1]));
    let center =
        std::array::from_fn::<_, 3, _>(|i| I::new(a[i].lo.min(b[i].lo), a[i].hi.max(b[i].hi)));
    let q0 = I::point(s.start[0]) - r - k * a[2];
    let q1 = I::point(s.end[0]) - r - k * b[2];
    let q = q0.magnitude().max(q1.magnitude());
    let radius = I::point(s.start[0].max(s.end[0]));
    let tilt = radius * (I::point(u[2].magnitude()) + I::point(v[2].magnitude()));
    let z = center[2] + I::new(-tilt.hi, tilt.hi);
    require((r + k * z).lo > 0.)?;
    let common =
        I::point(norm([center[0], center[1]])?) + I::point(q) + I::point(k.magnitude()) * tilt;
    let mut best = f64::INFINITY;
    for sign in [-1., 1.] {
        let b = common + radius * I::point(polar_errors(u, v, e, sign)?);
        require(valid(b))?;
        best = best.min(b.hi);
    }
    Ok(best)
}
fn sphere(
    raw: &RawSupport,
    c: CommonMouthChart,
    center: [f64; 2],
    radius: f64,
    correction: f64,
    r: I,
    sigma: [f64; 2],
) -> Result<f64> {
    let local = std::array::from_fn::<_, 3, _>(|i| {
        I::point(c.columns[0][i])
            + I::point(c.columns[3][i]) * (I::point(center[1]) - I::point(c.outer_top_z))
            - raw.frame.origin[i]
    });
    // r_local = center.r*e_theta + radius*u, |u|=1. This includes a
    // nonzero proposed major radius conservatively; the usual sphere has zero.
    let span = I::point(radius) + I::new(-center[0].abs(), center[0].abs());
    let radial = I::new(span.lo.max(0.), span.hi) * I::new(sigma[0], sigma[1]) - r;
    let b =
        I::point(norm(local)?) + I::point(radial.magnitude()) + repair_world(sigma[1], correction);
    require(valid(b))?;
    Ok(b.hi)
}
fn torus(
    raw: &RawSupport,
    c: CommonMouthChart,
    center: [f64; 2],
    radius: f64,
    correction: f64,
    radii: [I; 2],
    sigma: [f64; 2],
) -> Result<f64> {
    let [major, minor] = radii;
    let u = vector(&raw.frame, c.columns[1]);
    let e = phase(u);
    let mut best = f64::INFINITY;
    for sign in [-1., 1.] {
        for axial in [-1., 1.] {
            let x: [I; 3] =
                std::array::from_fn(|i| raw.frame.axes[0][i] * e[0] + raw.frame.axes[1][i] * e[1]);
            let y: [I; 3] = std::array::from_fn(|i| {
                I::point(sign) * (-raw.frame.axes[0][i] * e[1] + raw.frame.axes[1][i] * e[0])
            });
            let n = raw.frame.axes[2].map(|v| v * I::point(axial));
            // Each ur/uz is a true unit-minor-circle component. Independent
            // [-1,1] coefficients overbound the complete finite source face.
            let delta = std::array::from_fn::<_, 3, _>(|i| {
                let d = I::point(c.columns[0][i])
                    + I::point(c.columns[3][i]) * (I::point(center[1]) - I::point(c.outer_top_z))
                    - raw.frame.origin[i];
                let a = I::point(c.columns[1][i]) * I::point(center[0]) - x[i] * major;
                let b = I::point(c.columns[2][i]) * I::point(center[0]) - y[i] * major;
                let ar = I::point(c.columns[1][i]) * I::point(radius) - x[i] * minor;
                let br = I::point(c.columns[2][i]) * I::point(radius) - y[i] * minor;
                let z = I::point(c.columns[3][i]) * I::point(radius) - n[i] * minor;
                d + I::new(-1., 1.) * a
                    + I::new(-1., 1.) * b
                    + I::new(-1., 1.) * ar
                    + I::new(-1., 1.) * br
                    + I::new(-1., 1.) * z
            });
            let b = I::point(norm(delta)?) + repair_world(sigma[1], correction);
            require(valid(b))?;
            best = best.min(b.hi);
        }
    }
    Ok(best)
}
pub(super) fn distance(
    raw: &RawSupport,
    chart: CommonMouthChart,
    s: Segment,
    correction: f64,
) -> Result<(RawHeadSupportKind, f64)> {
    let sigma = stretch(chart)?;
    let out = match (raw.kind, s.curve) {
        (Kind::Plane, Curve::Line) => {
            let bounds =
                std::array::from_fn(|i| [s.start[i].min(s.end[i]), s.start[i].max(s.end[i])]);
            (RawHeadSupportKind::Plane, plane(raw, chart, bounds)?)
        }
        (Kind::Cylinder(r), Curve::Line) => (
            RawHeadSupportKind::Cylinder,
            tube(raw, chart, s, r, I::point(0.))?,
        ),
        (Kind::Cone(r, k), Curve::Line) => (RawHeadSupportKind::Cone, tube(raw, chart, s, r, k)?),
        (
            Kind::Sphere(r),
            Curve::Arc {
                center, radius_mm, ..
            },
        ) => (
            RawHeadSupportKind::Sphere,
            sphere(raw, chart, center, radius_mm, correction, r, sigma)?,
        ),
        (
            Kind::Torus(r, t),
            Curve::Arc {
                center, radius_mm, ..
            },
        ) => (
            RawHeadSupportKind::Torus,
            torus(raw, chart, center, radius_mm, correction, [r, t], sigma)?,
        ),
        _ => return Err(ReadError::Unsupported),
    };
    require(out.1.is_finite() && out.1 >= 0.)?;
    Ok(out)
}

pub(super) fn top_cap_distance(
    raw: &RawSupport,
    chart: CommonMouthChart,
    radius: f64,
    z: f64,
) -> Result<f64> {
    require(matches!(raw.kind, Kind::Plane))?;
    plane(raw, chart, [[0., radius], [z, z]])
}
