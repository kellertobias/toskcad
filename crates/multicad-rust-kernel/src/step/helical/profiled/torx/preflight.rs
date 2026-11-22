//! Qualify the actual literal supports emitted by the canonical STEP family.
//! A positional trim tube alone never proves raw support/vertex incidence.
use super::*;
use crate::{hexalobular_drive::BoundarySource, tapered_thread::interval::I};

fn point(p: Point) -> [I; 3] {
    p.map(|v| I::new(v.next_down(), v.next_up()))
}
fn dot(a: [I; 3], b: [I; 3]) -> I {
    a.into_iter()
        .zip(b)
        .fold(I::point(0.), |s, (a, b)| s + a * b)
}
fn cross(a: [I; 3], b: [I; 3]) -> [I; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn zero(values: impl IntoIterator<Item = I>) -> Result<(), StepError> {
    if values
        .into_iter()
        .all(|v| v.valid() && v.lo <= 0. && v.hi >= 0.)
    {
        Ok(())
    } else {
        Err(StepError::Unsupported)
    }
}
fn unit(vector: [I; 3]) -> Result<[I; 3], StepError> {
    let scale = vector.into_iter().map(I::magnitude).fold(0_f64, f64::max);
    if !scale.is_finite() || scale == 0. {
        return Err(StepError::InvalidGeometry);
    }
    let scaled = vector.map(|v| v / I::point(scale));
    let norm = dot(scaled, scaled).sqrt();
    if !norm.valid() || norm.lo <= 0. {
        return Err(StepError::Unsupported);
    }
    Ok(scaled.map(|v| v / norm))
}
fn delta(a: [I; 3], b: [I; 3]) -> [I; 3] {
    std::array::from_fn(|i| a[i] - b[i])
}

pub(in super::super) fn admission(
    screw: &Screw,
    recess: &ConicalHexalobular,
    tolerance: f64,
) -> Result<(), StepError> {
    let source = recess.source();
    if recess.finite_binding().is_none()
        && matches!(source.boundary_source(), BoundarySource::RetainedFamily)
    {
        return Err(StepError::Unsupported);
    }
    let angle = recess
        .finite_binding()
        .and_then(|b| b.literal_cone_semi_angle_rad())
        .unwrap_or_else(|| source.cone().depth_per_radius.recip().atan());
    if angle <= 0. || angle >= std::f64::consts::FRAC_PI_2 {
        return Err(StepError::Unsupported);
    }
    let error = screw.recess().source_error_mm();
    let radius = (screw.recess().maximum_radius_bound_mm() + error).next_up();
    let depth = (screw.recess().maximum_depth_mm() + error).next_up();
    if error > tolerance / 8.
        || !screw.admits_recess_envelope(radius, depth, 10. * tolerance)
        || source.arcs().iter().any(|arc| {
            arc.radius_mm <= 10. * tolerance
                || arc.radius_mm * arc.sweep_rad.abs() <= 10. * tolerance
        })
    {
        return Err(StepError::Limit);
    }
    Ok(())
}

pub(in super::super) fn placement(screw: &Screw, p: &Placement) -> Result<(), StepError> {
    let crate::screw_recess::Recess::ConicalHexalobular(recess) = screw.recess() else {
        return Ok(());
    };
    let source = recess.source();
    let binding = recess.finite_binding();
    let mouth = source
        .arcs()
        .map(|arc| [arc.start_mm[0], arc.start_mm[1], 0.]);
    let floor = mouth.map(|[x, y, z]| {
        [
            x,
            y,
            z + source.cone().tip_depth_mm - source.cone().depth_per_radius * x.hypot(y),
        ]
    });
    let vertices = binding.map_or(
        FiniteVertices {
            mouth_mm: mouth,
            floor_mm: floor,
            apex_mm: [0., 0., source.cone().tip_depth_mm],
        },
        |b| *b.vertices(),
    );
    let top = screw.head_top_z();
    let world = |[x, y, z]: Point| p.point_f64([x, y], top + z);
    let n = unit(point(p.normal))?;
    let origin = point(p.origin);
    let top_origin = point(world([0.; 3]));
    let apex = point(world(vertices.apex_mm));
    zero(cross(delta(apex, origin), n))?;
    for (i, arc) in source.arcs().iter().enumerate() {
        let center = point(p.point_f64(arc.center_mm, top));
        zero([dot(delta(center, top_origin), n)])?;
        let radius = I::new(arc.radius_mm.next_down(), arc.radius_mm.next_up());
        for j in [i, (i + 1) % 12] {
            let mouth = point(world(vertices.mouth_mm[j]));
            let d = delta(mouth, center);
            zero([dot(d, n), dot(d, d) - radius.square()])?;
        }
        let post = binding.map_or(
            FinitePost {
                origin_mm: vertices.mouth_mm[i],
                direction_mm_per_parameter: [0., 0., 1.],
                parameters: [0., vertices.floor_mm[i][2] - vertices.mouth_mm[i][2]],
            },
            |b| b.post_lines()[i],
        );
        let raw_vector = p.vector(post.direction_mm_per_parameter);
        let magnitude = raw_vector[0].hypot(raw_vector[1]).hypot(raw_vector[2]);
        let direction = unit(point(raw_vector))?;
        let vector = direction.map(|v| v * I::new(magnitude.next_down(), magnitude.next_up()));
        zero(cross(vector, n))?;
        let index = (0..3)
            .max_by(|&a, &b| raw_vector[a].abs().total_cmp(&raw_vector[b].abs()))
            .ok_or(StepError::InvalidGeometry)?;
        let origin = point(world(post.origin_mm));
        let mut parameters = [I::point(0.); 2];
        for (j, vertex) in [vertices.mouth_mm[i], vertices.floor_mm[i]]
            .into_iter()
            .enumerate()
        {
            let vertex = point(world(vertex));
            let parameter = (vertex[index] - origin[index]) / vector[index];
            zero((0..3).map(|k| origin[k] + vector[k] * parameter - vertex[k]))?;
            parameters[j] = parameter;
        }
        let advance = parameters[1] - parameters[0];
        if !(advance.valid()
            && if post.parameters[1] > post.parameters[0] {
                advance.lo > 0.
            } else {
                advance.hi < 0.
            })
        {
            return Err(StepError::Unsupported);
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn directed_normalization_encloses_finite_supports_and_refuses_a_zero_or_unresolved_normal() {
        assert!(matches!(
            unit([I::point(0.); 3]),
            Err(StepError::InvalidGeometry)
        ));
        assert!(matches!(
            unit([I::new(-1., 1.); 3]),
            Err(StepError::Unsupported)
        ));
        let vector = [1e200, -2e200, 3e200];
        let normalized = unit(point(vector)).unwrap();
        for (v, bounds) in [1. / 14_f64.sqrt(), -2. / 14_f64.sqrt(), 3. / 14_f64.sqrt()]
            .into_iter()
            .zip(normalized)
        {
            assert!(bounds.lo <= v && bounds.hi >= v);
        }
    }
}
