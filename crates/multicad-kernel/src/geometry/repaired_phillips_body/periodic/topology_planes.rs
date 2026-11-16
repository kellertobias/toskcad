//! Actual source plane support and bounded attachment-frame publication.
use super::*;
use multicad_core::references::ReferenceFrame;
use multicad_rust_kernel::phillips_drive::Surface;
use multicad_rust_kernel::step_read::PhillipsMaterialFaceSource as FaceSource;
fn sum(p: [f64; 2]) -> I {
    I::point(p[0]) + I::point(p[1])
}
fn dot(a: [I; 3], b: [I; 3]) -> I {
    (0..3).fold(I::ZERO, |s, i| s + a[i] * b[i])
}
fn cross(a: [I; 3], b: [I; 3]) -> [I; 3] {
    std::array::from_fn(|i| a[(i + 1) % 3] * b[(i + 2) % 3] - a[(i + 2) % 3] * b[(i + 1) % 3])
}
fn norm(a: [I; 3]) -> Result<[I; 3], GeometryError> {
    let d = dot(a, a);
    let d = I::new(d.lo.max(0.), d.hi).sqrt();
    if !d.valid() || d.lo <= 0. {
        return Err(unsupported("Source plane normal is unresolved"));
    }
    Ok(a.map(|v| v / d))
}
fn midpoint(a: [I; 3]) -> Point {
    a.map(|v| v.lo.midpoint(v.hi))
}
pub(super) fn plane(
    source: &OwnedPeriodicMaterial,
    face: FaceSource,
    placement: Placement,
    policy: GeometryPolicy,
) -> Result<Option<ReferenceFrame>, GeometryError> {
    let axial = match face {
        FaceSource::ThreadTip => Some((sum(source.cap_expressions_mm()[1]), 1.)),
        FaceSource::JoinAnnulus => Some((sum(source.cap_expressions_mm()[0]), 1.)),
        FaceSource::TopAnnulus => Some((sum(source.cavity_mouth_expression_mm()), -1.)),
        FaceSource::HeadMeridian { segment } => {
            let view = source.head().view();
            let s = view
                .profile_proposal()
                .segments()
                .nth(segment)
                .expect("retained segment");
            // A repaired line retains its original endpoints and is exactly
            // linear. A curved meridian cannot acquire a planar descriptor.
            if matches!(s.curve, multicad_rust_kernel::head_profile::Curve::Line)
                && s.start[1] == s.end[1]
            {
                Some((
                    I::point(s.start[1]) - I::point(source.head_z_shift_mm()),
                    -(s.end[0] - s.start[0]).signum(),
                ))
            } else {
                None
            }
        }
        _ => None,
    };
    let (anchor, normal) = if let Some((z, sign)) = axial {
        ([I::ZERO, I::ZERO, z], [I::ZERO, I::ZERO, I::point(sign)])
    } else if let FaceSource::ReversedCavity { face } = face {
        let topology = source
            .cavity()
            .boundary_topology()
            .map_err(|e| unsupported(format!("Actual cavity plane: {e:?}")))?;
        let index = topology
            .faces()
            .iter()
            .position(|f| f.source_identity() == face)
            .expect("retained cavity face");
        let Surface::Plane {
            normal_bounds,
            offset_mm_bounds,
        } = topology
            .face_surface(index)
            .map_err(|e| unsupported(format!("Actual cavity plane: {e:?}")))?
        else {
            return Ok(None);
        };
        let n = normal_bounds.map(|v| I::new(v[0], v[1]));
        let d = dot(n, n);
        if !d.valid() || d.lo <= 0. {
            return Err(unsupported("Source plane denominator is unresolved"));
        }
        let mut p = n.map(|v| v * I::new(offset_mm_bounds[0], offset_mm_bounds[1]) / d);
        p[2] = p[2] + sum(source.cavity_mouth_expression_mm());
        (p, n.map(|v| -v))
    } else {
        return Ok(None);
    };
    let axes: [[I; 3]; 3] = std::array::from_fn(|i| {
        let mut e = [I::ZERO; 3];
        e[i] = I::ONE;
        placement.linear_bounds(e)
    });
    let cofactor = [
        cross(axes[1], axes[2]),
        cross(axes[2], axes[0]),
        cross(axes[0], axes[1]),
    ];
    let actual = norm(std::array::from_fn(|i| {
        (0..3).fold(I::ZERO, |s, j| s + cofactor[j][i] * normal[j])
    }))?;
    let n = midpoint(actual);
    let local_n = midpoint(norm(normal)?);
    let index = (0..3)
        .min_by(|&i, &j| local_n[i].abs().total_cmp(&local_n[j].abs()))
        .expect("three axes");
    let mut axis = [I::ZERO; 3];
    axis[index] = I::ONE;
    let local_u = midpoint(norm(cross(local_n.map(I::point), axis))?);
    let u = midpoint(norm(placement.linear_bounds(local_u.map(I::point)))?);
    let v = midpoint(norm(cross(n.map(I::point), u.map(I::point)))?);
    let published = cross(u.map(I::point), v.map(I::point));
    let error = dot(
        std::array::from_fn(|i| I::point((published[i] - actual[i]).magnitude())),
        [I::ONE; 3],
    );
    let origin = placement.bounded_point(anchor, policy.linear_tolerance_mm / 32.)?;
    let world = placement.point_bounds(anchor);
    let ep = (0..3).fold(I::ZERO, |s, i| {
        s + I::point((world[i] - I::point(origin[i])).magnitude())
    });
    let domain = source.local_enclosure_bounds().map(|v| I::new(v[0], v[1]));
    let relative = placement.linear_bounds(std::array::from_fn(|i| domain[i] - anchor[i]));
    let radius = relative
        .into_iter()
        .fold(I::ZERO, |s, v| s + I::point(v.magnitude()));
    let deviation = (ep * (I::ONE + error) + error * radius) / dot(published, published).sqrt();
    if !deviation.valid() || deviation.hi > policy.linear_tolerance_mm / 8. {
        return Err(unsupported(
            "Source plane attachment reaches the precision floor",
        ));
    }
    let frame = ReferenceFrame {
        origin_mm: origin,
        u,
        v,
    };
    frame.validate().map_err(unsupported)?;
    Ok(Some(frame))
}
