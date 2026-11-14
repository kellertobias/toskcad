//! Untrusted numerical proposal from literal supports, never an acceptance proof.
use super::*;
use crate::phillips_drive::{Domain, Inner, Spec};
fn dot(a: Point, b: Point) -> f64 {
    a.into_iter().zip(b).map(|(a, b)| a * b).sum()
}
fn sub(a: Point, b: Point) -> Point {
    std::array::from_fn(|i| a[i] - b[i])
}
fn cross(a: Point, b: Point) -> Point {
    std::array::from_fn(|i| a[(i + 1) % 3] * b[(i + 2) % 3] - a[(i + 2) % 3] * b[(i + 1) % 3])
}
fn unit(a: Point) -> Point {
    let n = dot(a, a).sqrt();
    a.map(|a| a / n)
}
pub(super) struct Candidate {
    pub source: Cutter,
    pub pose: Pose,
}
// Scalar fitting is an untrusted proposal, not exact decoding of an authoring
// recipe. The original raw Body remains unchanged. Only compare() earns the
// separate finite-graph and whole-source measured correspondence. Selection
// thresholds below are hints/refusals; they never accept geometric incidence.
pub(super) fn candidate(b: &Boundary, work: &mut Work<'_>) -> Result<Candidate> {
    work.step(4096)?;
    let fail = ReadError::Unsupported;
    // Sealed physical::incidence::read already requires29views,49vertices
    // and76edges before publishing this immutable Boundary. Clones preserve
    // those sizes; the final graph separately proves candidate bijection.
    let (axis, minor) = b
        .faces()
        .iter()
        .find_map(|f| match f.support {
            Surface::Torus {
                axis,
                minor_radius_mm,
                ..
            } => Some((axis, minor_radius_mm)),
            _ => None,
        })
        .ok_or(fail)?;
    let n = unit(axis.normal);
    let v = unit(cross(n, axis.reference));
    let u = cross(v, n);
    let origin = std::array::from_fn(|i| axis.origin[i] - minor * n[i]);
    let local = |p: Point| {
        [
            dot(sub(p, origin), u),
            dot(sub(p, origin), v),
            dot(sub(p, origin), n),
        ]
    };
    let vector = |p: Point| [dot(p, u), dot(p, v), dot(p, n)];
    let mut cones = b
        .faces()
        .iter()
        .filter_map(|f| match f.support {
            Surface::Cone {
                axis,
                radius_mm,
                semi_angle_rad,
            } => Some((axis, radius_mm, semi_angle_rad)),
            _ => None,
        })
        .collect::<Vec<_>>();
    cones.sort_by(|a, b| a.2.total_cmp(&b.2));
    let (_, mhalf, alpha) = *cones.first().ok_or(fail)?;
    let (_, floor0, phi) = *cones.last().ok_or(fail)?;
    let k = alpha.tan();
    let kc = phi.tan();
    // Intersect the two actual meridian cone laws r=a-k*z.
    let floor = (mhalf - floor0) / (k - kc);
    let g = 2. * (mhalf - k * floor);
    let tip = floor0 / kc;
    let (caxis, inner) = b
        .faces()
        .iter()
        .find_map(|f| match f.support {
            Surface::Cylinder { axis, radius_mm }
                if local(axis.origin)[0] > 0. && local(axis.origin)[1].abs() < 1e-8 =>
            {
                Some((axis, radius_mm))
            }
            _ => None,
        })
        .ok_or(fail)?;
    let cn = vector(unit(caxis.normal));
    let beta = (-cn[0]).atan2(cn[2]);
    let theta = b
        .faces()
        .iter()
        .find_map(|f| match f.support {
            Surface::Plane { axis } => {
                let x = vector(unit(axis.normal));
                (x[0] > 0. && x[1] < 0. && x[2] > 0.).then(|| x[0].hypot(x[2]).atan2(-x[1]))
            }
            _ => None,
        })
        .ok_or(fail)?;
    // Cylinder center x = b/2 + tan(beta)*(tip-b/(2*kc))
    //                    + inner_radius/cos(beta).
    let bwidth = 2. * (local(caxis.origin)[0] - tip * beta.tan() - inner / beta.cos())
        / (1. - beta.tan() / kc);
    let spec = Spec {
        mouth_diameter_mm: 2. * mhalf,
        mouth_radius_mm: minor,
        b_width_mm: bwidth,
        g_diameter_mm: g,
        tilt_angle_deg: beta.to_degrees(),
        outer_flank_angle_deg: alpha.to_degrees(),
        floor_angle_deg: (1. / kc).atan().to_degrees(),
        outer_corner_angle_deg: 2. * theta.to_degrees(),
        inner: Inner::Rounded { radius_mm: inner },
    };
    // The independent actual torus major radius is deliberately NOT inferred
    // into this proposal. The final whole-support proof must bound its actual
    // disagreement with this new coherent source; no heuristic can accept it.
    let reach = b
        .vertices()
        .iter()
        .flat_map(|v| local(v.point))
        .map(f64::abs)
        .fold(0., f64::max)
        .max(1.);
    // These are deliberately new inactive padding bounds, not recovered
    // original prism limits. Full topology/trim matching below must reject
    // any padding that changes the actual visible finite component.
    let domain = Domain {
        corner_start_depth_mm: -4. * reach,
        corner_end_depth_mm: 4. * reach,
        corner_reach_mm: 16. * reach,
    };
    // Constant authored construction and its fixed four-corner algebra are
    // covered by the caller's fixed source-stage work reservation.
    let source =
        Cutter::new_with_cancel(spec, domain, &|| work.interrupted()).map_err(cutter_error)?;
    let pose = Pose::new(origin, u, v).map_err(|_| fail)?;
    Ok(Candidate { source, pose })
}
