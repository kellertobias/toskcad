//! Plane supports from the actual immutable degree-eleven material inventory.
use super::*;
use multicad_core::references::ReferenceFrame;
use multicad_rust_kernel::{
    phillips_drive::Surface, step_read::PhillipsDegree11MaterialFace as FaceSource,
};
pub(super) fn plane(
    source: &Owned,
    index: usize,
    placement: Placement,
    policy: GeometryPolicy,
) -> Result<Option<ReferenceFrame>, GeometryError> {
    let face = source
        .face_source(index)
        .ok_or_else(|| input("Unknown degree-eleven face"))?;
    let axial = match face {
        FaceSource::ThreadTip { .. } => Some((source.cap_stations_mm()[1], 1.)),
        FaceSource::JoinAnnulus { .. } => Some((source.cap_stations_mm()[0], 1.)),
        FaceSource::TopAnnulus { .. } => Some((source.chart().outer_top_z, -1.)),
        FaceSource::HeadMeridian { segment, .. } => {
            let head = source.head().view();
            let s = head
                .profile_proposal()
                .segments()
                .nth(segment)
                .ok_or_else(|| unsupported("Missing retained head segment"))?;
            if matches!(s.curve, multicad_rust_kernel::head_profile::Curve::Line)
                && s.start[1] == s.end[1]
            {
                Some((s.start[1], -(s.end[0] - s.start[0]).signum()))
            } else {
                None
            }
        }
        _ => None,
    };
    let (anchor, normal) = if let Some((z, sign)) = axial {
        (
            [I::ZERO, I::ZERO, I::point(z)],
            [I::ZERO, I::ZERO, I::point(sign)],
        )
    } else if let FaceSource::ReversedCavity { face, .. } = face {
        let topology = source
            .cavity()
            .boundary_topology()
            .map_err(|e| unsupported(format!("Actual cavity plane: {e:?}")))?;
        let index = topology
            .faces()
            .iter()
            .position(|f| f.source_identity() == face)
            .ok_or_else(|| unsupported("Missing retained cavity face"))?;
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
        let d = n.into_iter().fold(I::ZERO, |sum, n| sum + n * n);
        if !d.valid() || d.lo <= 0. {
            return Err(unsupported("Actual plane denominator unresolved"));
        }
        let mut p = n.map(|n| n * I::new(offset_mm_bounds[0], offset_mm_bounds[1]) / d);
        p[2] = p[2] + I::point(source.chart().outer_top_z);
        (p, n.map(|v| -v))
    } else {
        return Ok(None);
    };
    super::super::periodic::topology_planes::publish(
        anchor,
        normal,
        source.local_enclosure_bounds_mm(),
        placement,
        policy,
    )
    .map(Some)
}
