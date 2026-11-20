//! Closed filled cylinders and annuli with parallel axes.
//!
//! Each body is a radial section (a disk or an annulus) times an axial
//! interval in one shared axial direction, so the squared distance splits
//! into a radial and an axial part. A rotationally symmetric annulus `A`
//! meets a connected section `B` exactly when the radii `|q|` of points
//! `q ∈ B` about A's center, an interval `[m, M]`, meet A's radii; the
//! radial distance is the gap between those two intervals.
use super::*;

type Point = [f64; 3];

struct Ring {
    origin: Point,
    u: Point,
    n: Point,
    radii: [f64; 2],
    height: f64,
    /// Chart drift bound: unit axes differ from the stored ones.
    drift: f64,
}

fn dot(a: Point, b: Point) -> f64 {
    (0..3).map(|i| a[i] * b[i]).sum()
}
fn norm(a: Point) -> f64 {
    a[0].hypot(a[1]).hypot(a[2])
}
fn sub(a: Point, b: Point) -> Point {
    std::array::from_fn(|i| a[i] - b[i])
}
fn cross(a: Point, b: Point) -> Point {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn unit(a: Point) -> Option<Point> {
    let n = norm(a);
    (n.is_finite() && n > 0.).then(|| a.map(|x| x / n))
}

fn ring(solid: &Solid) -> Option<Ring> {
    let (radii, height, frame) = if let Some((cylinder, frame)) = solid.axial_cylinder_frame() {
        (
            [0., cylinder.circle().radius().as_f64()],
            cylinder.height().as_f64(),
            frame,
        )
    } else {
        let (annulus, frame) = solid.rust_annulus_data()?;
        (
            [annulus.inner_radius_mm(), annulus.outer_radius_mm()],
            annulus.height_mm(),
            frame,
        )
    };
    let u = unit(frame.u)?;
    let n = unit(frame.normal())?;
    let v = cross(n, u);
    let drift = (norm(sub(frame.u, u)) + norm(sub(frame.v, v))) * radii[1]
        + norm(sub(frame.normal(), n)) * height;
    Some(Ring {
        origin: frame.origin,
        u,
        n,
        radii,
        height,
        drift,
    })
}

pub(super) fn measure(
    a: &Solid,
    sa: BRepSelection,
    b: &Solid,
    sb: BRepSelection,
    tolerance: f64,
) -> Option<Result<MinimumDistance, GeometryError>> {
    if (sa, sb) != (BRepSelection::Whole, BRepSelection::Whole)
        || (a.rust_annulus_data().is_none() && b.rust_annulus_data().is_none())
    {
        return None;
    }
    let rings = (ring(a)?, ring(b)?);
    Some(solve(rings.0, rings.1, tolerance))
}

fn solve(a: Ring, b: Ring, tolerance: f64) -> Result<MinimumDistance, GeometryError> {
    let unproved = |message: &str| error(GeometryErrorKind::Unsupported, message);
    let n = a.n;
    let delta = sub(b.origin, a.origin);
    let reach = norm(delta) + a.radii[1] + b.radii[1] + a.height + b.height;
    let magnitude = a
        .origin
        .into_iter()
        .chain(b.origin)
        .map(f64::abs)
        .fold(reach.max(1.), f64::max);
    let uncertainty =
        a.drift + b.drift + norm(cross(a.n, b.n)) * reach + 128. * f64::EPSILON * magnitude;
    if !uncertainty.is_finite() || uncertainty > tolerance / 32. {
        return Err(unproved(
            "Parallel ring distance requires parallel axes within the numerical budget",
        ));
    }
    let sign = if dot(a.n, b.n) < 0. { -1. } else { 1. };
    let z0 = dot(delta, n);
    let z1 = z0 + sign * b.height;
    let (low, high) = (z0.min(z1), z0.max(z1));
    let lateral: Point = std::array::from_fn(|i| delta[i] - z0 * n[i]);
    let e = norm(lateral);
    let centered = e <= uncertainty;
    let w = if centered {
        a.u
    } else {
        lateral.map(|x| x / e)
    };
    let p = cross(n, w);
    let [ia, oa] = a.radii;
    let [ib, ob] = b.radii;
    // Radii about A's center reached by B's section.
    let reached = [if e < ib { ib - e } else { (e - ob).max(0.) }, e + ob];
    let radial = (reached[0] - oa).max(ia - reached[1]).max(0.);
    let axial = (low - a.height).max(-high).max(0.);
    // Radial witnesses as coefficients of (w, p) about A's center.
    let (qa, qb, radial_unique) = if ia - reached[1] > 0. {
        ([ia, 0.], [e + ob, 0.], !centered)
    } else if reached[0] - oa > 0. {
        if e < ib {
            ([-oa, 0.], [e - ib, 0.], !centered)
        } else {
            ([oa, 0.], [e - ob, 0.], true)
        }
    } else {
        let r = (reached[0].max(ia) + reached[1].min(oa)) / 2.;
        let q = if centered {
            [r, 0.]
        } else {
            let lo = ib.max((r - e).abs());
            let hi = ob.min(r + e);
            let rho = (lo + hi) / 2.;
            let x = (r * r - rho * rho + e * e) / (2. * e);
            [x, (r * r - x * x).max(0.).sqrt()]
        };
        (q, q, false)
    };
    let (za, zb) = if low - a.height > 0. {
        (a.height, low)
    } else if high < 0. {
        (0., high)
    } else {
        let z = (low.max(0.) + high.min(a.height)) / 2.;
        (z, z)
    };
    let place = |q: [f64; 2], z: f64| -> Point {
        std::array::from_fn(|i| a.origin[i] + q[0] * w[i] + q[1] * p[i] + z * n[i])
    };
    let pa = place(qa, za);
    let pb = place(qb, zb);
    let distance = radial.hypot(axial);
    let actual = norm(sub(pa, pb));
    if !actual.is_finite() || (actual - distance).abs() > tolerance / 4. {
        return Err(unproved(
            "Parallel ring distance witness exceeds the numerical budget",
        ));
    }
    let unique = radial > 0. && axial > 0. && radial_unique;
    let overlap_radial = reached[1].min(oa) - reached[0].max(ia);
    let overlap_axial = high.min(a.height) - low.max(0.);
    Ok(MinimumDistance {
        distance_mm: if distance == 0. { 0. } else { actual },
        witnesses: vec![DistanceWitness { a_mm: pa, b_mm: pb }],
        witnesses_complete: unique,
        total_solutions: unique.then_some(1),
        inner_solution: distance == 0.
            && overlap_radial > uncertainty
            && overlap_axial > uncertainty,
    })
}
