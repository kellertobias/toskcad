//! Directed raw-frame transport; area factors are never assumed invariant.
use super::*;
use crate::geometry::internal_thread_treatment_body::interval::Interval as I;
fn p(v: f64) -> I {
    I::point(v)
}
fn b(v: [f64; 2]) -> I {
    I::new(v[0], v[1])
}
fn refused() -> GeometryError {
    error(
        OPERATION,
        GeometryErrorKind::Unsupported,
        "Actual Phillips world face precision is unresolved",
    )
}
fn dot(a: [I; 3], b: [I; 3]) -> I {
    (0..3).fold(I::ZERO, |s, i| s + a[i] * b[i])
}
fn cross(a: [I; 3], b: [I; 3]) -> [I; 3] {
    [
        (a[1] * b[2] - a[2] * b[1]),
        (a[2] * b[0] - a[0] * b[2]),
        (a[0] * b[1] - a[1] * b[0]),
    ]
}
fn positive(v: I) -> Result<I, GeometryError> {
    if v.valid() && v.lo > 0. {
        Ok(v)
    } else {
        Err(refused())
    }
}
fn normalized(v: [I; 3]) -> Result<[I; 3], GeometryError> {
    let n = positive(dot(v, v).sqrt())?;
    Ok(v.map(|v| v / n))
}
fn mid(v: [I; 3]) -> Point {
    v.map(|b| b.lo.midpoint(b.hi))
}
fn axes(frame: PlaneFrame) -> [[I; 3]; 3] {
    [frame.u, frame.v, frame.normal()].map(|v| v.map(p))
}
fn world(frame: PlaneFrame, v: [I; 3]) -> [I; 3] {
    let a = axes(frame);
    std::array::from_fn(|i| p(frame.origin[i]) + (0..3).fold(I::ZERO, |s, k| s + a[k][i] * v[k]))
}
fn extent(source: &PhillipsScrew) -> [I; 3] {
    source.source_coordinate_bounds_mm().map(b)
}
pub(super) fn bounds(source: &PhillipsScrew, frame: PlaneFrame) -> Bounds {
    let domain = world(frame, extent(source));
    Bounds {
        min: domain.map(|v| v.lo),
        max: domain.map(|v| v.hi),
    }
}

fn metric(frame: PlaneFrame, planar: bool) -> Result<(I, bool), GeometryError> {
    let a = axes(frame);
    let mut e = I::ZERO;
    for i in 0..3 {
        for j in 0..3 {
            e = e + p((dot(a[i], a[j]) - p(if i == j { 1. } else { 0. })).magnitude());
        }
    }
    let determinant = dot(cross(a[0], a[1]), a[2]);
    if !e.valid() || e.hi > 1e-6 || determinant.lo <= 0. {
        return Err(refused());
    }
    if planar {
        return Ok((
            positive(dot(cross(a[0], a[1]), cross(a[0], a[1])).sqrt())?,
            true,
        ));
    }
    Ok((
        I::new((1. - e.hi).next_down(), (1. + e.hi).next_up()),
        false,
    ))
}
fn determinant(frame: PlaneFrame) -> Result<I, GeometryError> {
    let a = axes(frame);
    positive(dot(cross(a[0], a[1]), a[2]))
}
pub(super) fn volume_bounds(
    source: source::SourceVolumeBounds,
    frame: PlaneFrame,
) -> Result<[f64; 2], GeometryError> {
    metric(frame, false)?;
    let value = positive(b(source.volume_mm3()) * determinant(frame)?)?;
    Ok([value.lo, value.hi])
}
pub(super) fn mesh_pair(
    volume: source::SourceVolumeBounds,
    frame: PlaneFrame,
    mesh: &Mesh,
    local_literal: [f64; 2],
    local_difference: f64,
    options: multicad_rust_kernel::display::Tessellation,
) -> Result<PhillipsMeshVolumePair, GeometryError> {
    let source_world = b(volume_bounds(volume, frame)?);
    let origin = mesh
        .vertices
        .as_chunks::<3>()
        .0
        .first()
        .copied()
        .ok_or_else(refused)?;
    let mut literal = I::ZERO;
    for tri in mesh.indices.as_chunks::<3>().0 {
        let v = tri.map(|id| {
            let q = mesh.vertices.as_chunks::<3>().0[id as usize];
            std::array::from_fn(|i| p(q[i]) - p(origin[i]))
        });
        literal = literal + dot(v[0], cross(v[1], v[2])) / p(6.);
    }
    let literal = positive(literal)?;
    let det = determinant(frame)?;
    let transfer = (literal - b(local_literal) * det).magnitude();
    let bound = (p(local_difference) * det + p(transfer)).hi;
    // Triangle/source pairing is independently retained, not inferred from
    // this volume diagnostic. These intervals describe the actual emitted
    // world vertices, including every literal transfer-rounding term.
    if !bound.is_finite()
        || !local_difference.is_finite()
        || local_difference < 0.
        || (literal - source_world).magnitude() > bound
    {
        return Err(refused());
    }
    Ok(PhillipsMeshVolumePair {
        source_world_mm3: [source_world.lo, source_world.hi],
        frame,
        linear_deflection_mm: options.deflection,
        angular_deflection_rad: options.angle,
        triangle_limit: options.triangle_limit,
        literal_world_mm3: [literal.lo, literal.hi],
        source_to_literal_bound_mm3: bound,
        literal_vertex_transfer_bound_mm3: transfer,
    })
}
fn axial_plane(source: &PhillipsScrew, role: ScrewFace) -> Option<(f64, f64)> {
    match role {
        ScrewFace::ThreadShoulder => Some((source.visible_thread_interval_mm()[0], 1.)),
        ScrewFace::ThreadTip => Some((source.visible_thread_interval_mm()[1], 1.)),
        ScrewFace::HeadTop => Some((source.head_top_z(), -1.)),
        ScrewFace::ProfileWall(i) => {
            let s = source.profile().segments().nth(i as usize)?;
            (s.start[1] == s.end[1]).then_some((
                s.start[1] - source.spec().shank_length_mm,
                -(s.end[0] - s.start[0]).signum(),
            ))
        }
        _ => None,
    }
}
pub(in crate::geometry::helical_body) fn plane(
    source: &PhillipsScrew,
    role: ScrewFace,
    frame: PlaneFrame,
    policy: GeometryPolicy,
) -> Result<Option<ReferenceFrame>, GeometryError> {
    frame.checked(policy)?;
    let (anchor, normal) = if let Some((z, sign)) = axial_plane(source, role) {
        ([I::ZERO, I::ZERO, p(z)], [I::ZERO, I::ZERO, p(sign)])
    } else if let ScrewFace::Cavity(face) = role {
        let topology = source
            .cavity()
            .source()
            .boundary_topology()
            .map_err(failure)?;
        let index = topology
            .faces()
            .iter()
            .position(|f| f.source_identity() == face)
            .ok_or_else(refused)?;
        let source::Surface::Plane {
            normal_bounds,
            offset_mm_bounds,
        } = topology.face_surface(index).map_err(failure)?
        else {
            return Ok(None);
        };
        let n = normal_bounds.map(b);
        let denominator = positive(dot(n, n))?;
        let mut anchor = n.map(|n| n * b(offset_mm_bounds) / denominator);
        anchor[2] = anchor[2] + p(source.head_top_z());
        (anchor, n.map(|v| -v))
    } else {
        return Ok(None);
    };
    let a = axes(frame);
    let cof = [cross(a[1], a[2]), cross(a[2], a[0]), cross(a[0], a[1])];
    positive(dot(cof[2], a[2]))?;
    let actual = normalized(std::array::from_fn(|i| {
        (0..3).fold(I::ZERO, |s, k| s + cof[k][i] * normal[k])
    }))?;
    let n = mid(actual);
    let index = (0..3)
        .min_by(|&i, &j| n[i].abs().total_cmp(&n[j].abs()))
        .expect("three coordinates");
    let mut axis = [I::ZERO; 3];
    axis[index] = I::ONE;
    let u = mid(normalized(cross(n.map(p), axis))?);
    let v = mid(normalized(cross(n.map(p), u.map(p)))?);
    let published = cross(u.map(p), v.map(p));
    let en = (0..3).fold(I::ZERO, |sum, i| {
        sum + p((published[i] - actual[i]).magnitude())
    });
    let anchor_world = world(frame, anchor);
    let origin = mid(anchor_world);
    let ep = (0..3).fold(I::ZERO, |sum, i| {
        sum + p((anchor_world[i] - p(origin[i])).magnitude())
    });
    let domain = extent(source);
    let relative: [I; 3] = std::array::from_fn(|i| domain[i] - anchor[i]);
    let delta = std::array::from_fn::<_, 3, _>(|i| {
        (0..3).fold(I::ZERO, |sum, k| sum + a[k][i] * relative[k])
    });
    let radius = delta
        .into_iter()
        .fold(I::ZERO, |sum, v| sum + p(v.magnitude()));
    let error = (ep * (I::ONE + en) + en * radius) / positive(dot(published, published).sqrt())?;
    if !error.valid() || error.hi > policy.linear_tolerance_mm / 8. {
        return Err(refused());
    }
    policy.point(origin)?;
    let result = ReferenceFrame {
        origin_mm: origin,
        u,
        v,
    };
    verified_reference(result).map(Some)
}
fn verified_reference(frame: ReferenceFrame) -> Result<ReferenceFrame, GeometryError> {
    frame.validate().map_err(|_| refused())?;
    Ok(frame)
}

impl Body {
    pub(super) fn faces(
        &self,
        frame: PlaneFrame,
        policy: GeometryPolicy,
    ) -> Result<Vec<Face>, GeometryError> {
        frame.checked(policy)?;
        let domain = extent(&self.source);
        let world_domain = world(frame, domain);
        policy.point(world_domain.map(|v| v.lo))?;
        policy.point(world_domain.map(|v| v.hi))?;
        let scale = domain.into_iter().map(I::magnitude).fold(1., f64::max);
        let mut result = Vec::with_capacity(self.properties.len());
        for (index, local) in self.properties.iter().enumerate() {
            let (jac, constant) = metric(frame, axial_plane(&self.source, local.face).is_some())?;
            let area = b(local.area_mm2) * jac;
            let c = local.centroid_mm.map(b);
            let center = jac.lo.midpoint(jac.hi);
            let delta = if constant {
                I::ZERO
            } else {
                p((p(center) - jac).magnitude())
            };
            let mut weighted = c;
            for i in 0..3 {
                let error = delta / p(jac.lo) * p((domain[i] - c[i]).magnitude());
                weighted[i] = (c[i] + I::new(-error.hi, error.hi))
                    .intersect(domain[i])
                    .ok_or_else(refused)?;
            }
            let centroid = world(frame, weighted);
            if !area.valid()
                || area.lo <= 0.
                || area.width() > policy.linear_tolerance_mm * scale / 8.
                || centroid
                    .iter()
                    .any(|v| !v.valid() || v.width() > policy.linear_tolerance_mm / 8.)
            {
                return Err(refused());
            }
            result.push(Face {
                index: index as u32,
                plane: plane(&self.source, local.face, frame, policy)?
                    .map(|p| (p.origin_mm, p.normal())),
                center_mm: mid(centroid),
                area_mm2: area.lo.midpoint(area.hi),
            });
        }
        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn standalone_actual_frame_transport_refuses_ambiguous_metric_plane_and_volume_bounds() {
        let s = super::super::tests::fixture(1, 1., false);
        assert!(positive(I::new(0., 1.)).is_err());
        assert!(
            verified_reference(ReferenceFrame {
                origin_mm: [0.; 3],
                u: [0.; 3],
                v: [1., 0., 0.]
            })
            .is_err()
        );
        assert!(normalized([I::ZERO; 3]).is_err());
        let collapse = PlaneFrame {
            u: [0.; 3],
            ..Default::default()
        };
        assert!(metric(collapse, false).is_err());
        let strain = PlaneFrame {
            u: [1.01, 0., 0.],
            ..Default::default()
        };
        assert!(metric(strain, false).is_err());
        assert!(volume_bounds(s.source_volume_bounds().unwrap(), collapse).is_err());
        assert!(
            plane(
                &s,
                ScrewFace::HeadTop,
                PlaneFrame {
                    origin: [1e6, -2e6, 3e6],
                    ..Default::default()
                },
                GeometryPolicy {
                    linear_tolerance_mm: 1e-14,
                    ..Default::default()
                }
            )
            .is_err()
        );
        let bounds = bounds(&s, PlaneFrame::default());
        assert!(bounds.min[2] <= s.head_top_z() && bounds.max[2] >= s.thread().spec().length_mm);
    }
    #[test]
    fn literal_world_pair_refuses_reversed_actual_mesh_and_nonfinite_pair_reserves() {
        let source = super::super::tests::fixture(1, 0., false);
        let body = Body::new(source, GeometryPolicy::default(), &|| false).unwrap();
        let options = Solid::rust_mesh_options(Default::default());
        let frame = PlaneFrame::default();
        let mesh = body.mesh(frame, options).unwrap();
        let pair = (*body.last_mesh_pair.borrow()).unwrap();
        let mut changed = mesh;
        for tri in changed.indices.as_chunks_mut::<3>().0 {
            tri.swap(1, 2);
        }
        assert!(
            mesh_pair(
                body.volume,
                frame,
                &changed,
                pair.literal_world_mm3,
                pair.source_to_literal_bound_mm3,
                options
            )
            .is_err()
        );
        for tri in changed.indices.as_chunks_mut::<3>().0 {
            tri.swap(1, 2);
        }
        assert!(
            mesh_pair(
                body.volume,
                frame,
                &changed,
                pair.literal_world_mm3,
                f64::NAN,
                options
            )
            .is_err()
        );
    }
}
