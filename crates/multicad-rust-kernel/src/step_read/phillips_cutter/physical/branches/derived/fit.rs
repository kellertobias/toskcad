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
    let k = alpha.portable_tan();
    let kc = phi.portable_tan();
    // Intersect the two actual meridian cone laws r=a-k*z.
    let floor = (mhalf - floor0) / (k - kc);
    let g = 2. * (mhalf - k * floor);
    let tip = floor0 / kc;
    let rounded = b.faces().iter().find_map(|f| match f.support {
        Surface::Cylinder { axis, radius_mm }
            if local(axis.origin)[0] > 0. && local(axis.origin)[1].abs() < 1e-8 =>
        {
            Some((axis, radius_mm))
        }
        _ => None,
    });
    let (beta, theta, bwidth, inner) = if let Some((caxis, inner)) = rounded {
        let cn = vector(unit(caxis.normal));
        let beta = (-cn[0]).portable_atan2(cn[2]);
        let theta = b
            .faces()
            .iter()
            .find_map(|f| match f.support {
                Surface::Plane { axis } => {
                    let x = vector(unit(axis.normal));
                    (x[0] > 0. && x[1] < 0. && x[2] > 0.)
                        .then(|| x[0].portable_hypot(x[2]).portable_atan2(-x[1]))
                }
                _ => None,
            })
            .ok_or(fail)?;
        // Cylinder center x = b/2 + tan(beta)*(tip-b/(2*kc))
        //                    + inner_radius/cos(beta).
        let bwidth = 2.
            * (local(caxis.origin)[0] - tip * beta.portable_tan() - inner / beta.portable_cos())
            / (1. - beta.portable_tan() / kc);

        (beta, theta, bwidth, Inner::Rounded { radius_mm: inner })
    } else {
        // These hints propose a shallow-V from its actual inner/outer plane
        // supports. No fixture index, head family or standard table is used.
        // Full graph bijection and whole-support comparisons remain mandatory.
        work.step(4096)?;
        let mut planes = b
            .faces()
            .iter()
            .filter_map(|f| match f.support {
                Surface::Plane { axis } => {
                    let n = vector(unit(axis.normal));
                    (n[0] > 0. && n[1] < 0. && n[2] > 0. && n[0] > -n[1])
                        .then(|| (n[0].portable_hypot(n[2]).portable_atan2(-n[1]), axis, n))
                }
                _ => None,
            })
            .collect::<Vec<_>>();
        planes.sort_by(|a, b| a.0.total_cmp(&b.0));
        let (theta, outer_axis, outer_n) = *planes.first().ok_or(fail)?;
        let (half, inner_axis, inner_n) = *planes.last().ok_or(fail)?;
        // The offset of an inner V plane divided by sin(half) is
        // S=cos(beta)*b/2+sin(beta)*(tip-b/(2kc)). The outer plane's
        // additional offset is e/2*(cot(half)*sin(theta)-cos(theta)).
        let beta = inner_n[2].portable_atan2(inner_n[0]);
        let inner_offset = dot(inner_n, local(inner_axis.origin));
        let outer_offset = dot(outer_n, local(outer_axis.origin));
        let center = inner_offset / half.portable_sin();
        let bwidth = 2. * (center - beta.portable_sin() * tip)
            / (beta.portable_cos() - beta.portable_sin() / kc);
        let ewidth = 2. * (outer_offset - theta.portable_sin() * center)
            / (half.portable_cos() / half.portable_sin() * theta.portable_sin()
                - theta.portable_cos());
        (
            beta,
            theta,
            bwidth,
            Inner::ShallowV {
                angle_deg: 2. * half.to_degrees(),
                e_width_mm: ewidth,
            },
        )
    };
    let spec = Spec {
        mouth_diameter_mm: 2. * mhalf,
        mouth_radius_mm: minor,
        b_width_mm: bwidth,
        g_diameter_mm: g,
        tilt_angle_deg: beta.to_degrees(),
        outer_flank_angle_deg: alpha.to_degrees(),
        floor_angle_deg: (1. / kc).portable_atan().to_degrees(),
        outer_corner_angle_deg: 2. * theta.to_degrees(),
        inner,
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
