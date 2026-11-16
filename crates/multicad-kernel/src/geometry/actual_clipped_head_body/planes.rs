//! Cheap earned support planes, without any scalar surface integration.
use super::*;
fn dot(a: [I; 3], b: [I; 3]) -> I {
    (0..3).fold(I::ZERO, |s, k| s + a[k] * b[k])
}
fn cross(a: [I; 3], b: [I; 3]) -> [I; 3] {
    std::array::from_fn(|k| a[(k + 1) % 3] * b[(k + 2) % 3] - a[(k + 2) % 3] * b[(k + 1) % 3])
}
fn unit(a: [I; 3]) -> Result<[I; 3], GeometryError> {
    let length = dot(a, a).sqrt();
    if !length.valid() || length.lo <= 0. {
        return Err(failure(read::Error::Unsupported));
    }
    let v = a.map(|v| v / length);
    if v.iter().any(|v| !v.valid()) {
        return Err(failure(read::Error::Unsupported));
    }
    Ok(v)
}
fn proposed_unit(a: [f64; 3]) -> Result<[f64; 3], GeometryError> {
    let length = a.into_iter().map(|v| v * v).sum::<f64>().sqrt();
    if !length.is_finite() || length <= 0. {
        return Err(failure(read::Error::Unsupported));
    }
    Ok(a.map(|v| v / length))
}
impl ActualClippedHeadBody {
    fn vector(&self, v: [I; 3], normal: bool) -> Result<[I; 3], GeometryError> {
        let frame = self.source.actual_frame();
        let axes = frame.axes();
        let scale = frame.scale_bounds();
        let value = std::array::from_fn(|k| {
            (0..3).fold(I::ZERO, |s, j| {
                s + I::new(axes[j][k][0], axes[j][k][1]) * v[j]
            })
        });
        let value = if normal {
            value
        } else {
            value.map(|v| v * I::new(scale[0], scale[1]))
        };
        let value = if let Some(a) = &self.occurrence {
            ranges(
                if normal {
                    a.normal_bounds(output(value))
                } else {
                    a.vector_bounds(output(value))
                }
                .map_err(failure)?,
            )
        } else {
            value
        };
        if value.iter().any(|v| !v.valid()) {
            return Err(failure(read::Error::Unsupported));
        }
        Ok(value)
    }
    pub(in crate::geometry) fn plane(
        &self,
        index: u32,
    ) -> Result<Option<ReferenceFrame>, GeometryError> {
        let index = index as usize;
        if index >= 17 {
            return Err(failure(read::Error::Invalid));
        }
        if index < 8 {
            return Ok(None);
        }
        let source = self.source.local_geometry();
        let (origin, normal, reference) = if index < 9 {
            let id = source.physical_roles().thread_face_entities[index];
            let aux = source.auxiliary_thread();
            let face = aux
                .faces()
                .iter()
                .find(|f| f.entity() == id)
                .ok_or_else(|| failure(read::Error::Invalid))?;
            let read::Support::Plane { height } = face.support() else {
                return Err(failure(read::Error::Invalid));
            };
            let p = aux
                .point(height)
                .ok_or_else(|| failure(read::Error::Invalid))?
                .coordinates()
                .map(|v| v.bounds());
            (
                source
                    .axial_sewing_transfer()
                    .map_raw_box_mm(p)
                    .map_err(failure)?,
                [0., 0., if face.forward() { 1. } else { -1. }],
                [1., 0., 0.],
            )
        } else {
            let face = &source.head().faces[index - 9];
            (
                face.plane.origin.map(|v| [v, v]),
                face.plane.normal.map(|v| if face.forward { v } else { -v }),
                face.plane.reference,
            )
        };
        let n = unit(normal.map(I::point))?;
        let r = reference.map(I::point);
        let along = dot(n, r);
        let u = unit(std::array::from_fn(|k| r[k] - n[k] * along))?;
        let world_origin = self.map_local_box(origin)?;
        let world_u = self.vector(u, false)?;
        let world_n = unit(self.vector(n, true)?)?;
        let u = proposed_unit(world_u.map(|v| v.lo.midpoint(v.hi)))?;
        let n = proposed_unit(world_n.map(|v| v.lo.midpoint(v.hi)))?;
        let v = proposed_unit(super::super::cross(n, u))?;
        let origin_mm = self.literal_point(world_origin)?;
        let frame = ReferenceFrame { origin_mm, u, v };
        frame
            .validate()
            .map_err(|_| failure(read::Error::Unsupported))?;
        let proposed = unit(cross(u.map(I::point), v.map(I::point)))?;
        let delta = (0..3).fold(I::ZERO, |s, k| {
            s + I::point((proposed[k] - world_n[k]).magnitude())
        });
        let radius = self
            .world_box()?
            .into_iter()
            .zip(ranges(world_origin))
            .fold(I::ZERO, |s, (b, o)| {
                s + I::point((I::new(b[0], b[1]) - o).magnitude())
            });
        let rounding = ranges(world_origin)
            .into_iter()
            .zip(origin_mm)
            .fold(I::ZERO, |s, (o, p)| {
                s + I::point((o - I::point(p)).magnitude())
            });
        let scale = I::point(self.source.actual_frame().scale_bounds()[1])
            * I::point(
                self.occurrence
                    .as_ref()
                    .map_or(1., |a| a.linear_l2_norm_upper()),
            );
        let collar = I::point(source.planar_sewing_mm()) * scale;
        let error = delta * radius + rounding + collar;
        if !error.valid() || error.hi > self.policy.linear_tolerance_mm / 8. {
            return Err(failure(read::Error::Unsupported));
        }
        Ok(Some(frame))
    }
}
