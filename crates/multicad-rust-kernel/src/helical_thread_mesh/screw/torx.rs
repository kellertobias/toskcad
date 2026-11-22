//! Exact-source arc stations and the actual cone floor, never a flat cap.
use super::*;
use crate::hexalobular_drive::Drive;
use std::f64::consts::PI;
#[cfg(test)]
mod tests;

pub(super) fn chord_request(source: &Drive, options: Tessellation) -> f64 {
    use crate::tapered_thread::interval::I;
    (I::point(options.deflection)
        / (I::point(4.) * (I::point(1.) + I::point(source.cone().depth_per_radius))))
    .lo
}
pub(super) fn angular_step(source: &Drive, options: Tessellation, previous: f64) -> f64 {
    previous.min(
        (options.deflection / (source.cone().depth_per_radius * source.maximum_radius_bound_mm()))
            .sqrt(),
    )
}
pub(super) fn stations(source: &Drive, options: Tessellation) -> Result<Vec<f64>, MeshError> {
    let delta = chord_request(source, options);
    let mut result = Vec::new();
    for arc in source.arcs() {
        let step = (options.angle * 0.5).min((4. * delta / arc.radius_mm).sqrt());
        let count = (arc.sweep_rad.abs() / step).ceil();
        if !count.is_finite() || (result.len() as f64 + count) * 2. > options.triangle_limit as f64
        {
            return Err(MeshError::BudgetExceeded);
        }
        for i in 0..count as usize {
            let point = if i == 0 {
                arc.start_mm
            } else {
                let t = arc.start_angle_rad + arc.sweep_rad * i as f64 / count;
                [
                    arc.center_mm[0] + arc.radius_mm * t.portable_cos(),
                    arc.center_mm[1] + arc.radius_mm * t.portable_sin(),
                ]
            };
            result.push(point[1].portable_atan2(point[0]).rem_euclid(TAU));
        }
    }
    Ok(result)
}

pub(super) fn volume_error(
    source: &Drive,
    options: Tessellation,
    reserve: f64,
    width: f64,
) -> Result<f64, MeshError> {
    use crate::tapered_thread::interval::I;
    let xy = (I::point(chord_request(source, options)) + I::point(reserve)).hi;
    let floor = (I::point(options.deflection) / I::point(4.)
        + I::point(source.cone().depth_per_radius) * I::point(xy)
        + I::point(reserve))
    .hi;
    let geometry = source
        .boundary_volume_error_bound_mm3(xy, floor)
        .ok_or(MeshError::UnresolvedBoundary)?;
    Ok((I::point(geometry) + I::point(width)).hi)
}

fn validate_ring(
    points: &[[f64; 3]],
    source: &Drive,
    options: Tessellation,
) -> Result<(), MeshError> {
    use crate::tapered_thread::interval::I;
    let reserve = 4. * source.source_error_mm();
    for (a, b) in points
        .iter()
        .zip(points.iter().cycle().skip(1))
        .take(points.len())
    {
        let theta = a[1].portable_atan2(a[0]).rem_euclid(TAU);
        let delta = (b[1].portable_atan2(b[0]) - theta).rem_euclid(TAU);
        let index = source
            .arc_index_at_angle(theta + delta * 0.5)
            .ok_or(MeshError::UnresolvedBoundary)?;
        let arc = &source.arcs()[index];
        let pa = [a[0] - arc.center_mm[0], a[1] - arc.center_mm[1]];
        let pb = [b[0] - arc.center_mm[0], b[1] - arc.center_mm[1]];
        let span = ((pb[1].portable_atan2(pb[0]) - pa[1].portable_atan2(pa[0]))
            * arc.sweep_rad.signum())
        .rem_euclid(TAU);
        let parameter_reserve = reserve / arc.radius_mm;
        let finite_endpoint = |p: [f64; 2]| {
            let offset = ((p[1].portable_atan2(p[0]) - arc.start_angle_rad)
                * arc.sweep_rad.signum())
            .rem_euclid(TAU);
            offset <= arc.sweep_rad.abs() + parameter_reserve || TAU - offset <= parameter_reserve
        };
        let chord = (I::point(arc.radius_mm)
            * (I::point(span) + I::point(parameter_reserve)).square()
            / I::point(8.))
        .hi;
        let floor = (I::point(source.cone().depth_per_radius)
            * I::point(source.maximum_radius_bound_mm())
            * I::point(delta).square()
            / I::point(8.))
        .hi;
        if !(delta > 0.
            && delta < PI
            && delta <= options.angle
            && span > 0.
            && span < PI
            && span + parameter_reserve <= options.angle
            && span <= arc.sweep_rad.abs() + parameter_reserve
            && finite_endpoint(pa)
            && finite_endpoint(pb)
            && (pa[0].portable_hypot(pa[1]) - arc.radius_mm).abs() <= reserve
            && (pb[0].portable_hypot(pb[1]) - arc.radius_mm).abs() <= reserve
            && chord <= chord_request(source, options)
            && floor <= options.deflection / 4.)
        {
            return Err(MeshError::UnresolvedBoundary);
        }
    }
    Ok(())
}

impl HeadMesh<'_, '_> {
    pub(super) fn torx_ring(
        &mut self,
        top: f64,
        source: &Drive,
        coupling: f64,
        floor: bool,
        options: Tessellation,
    ) -> Result<Vec<u32>, MeshError> {
        let points = source.arcs().map(|arc| arc.start_mm);
        let corners = self.corner_stations(&points, coupling)?;
        let mut positions = Vec::with_capacity(self.angles.len());
        for (i, &angle) in self.angles.iter().enumerate() {
            let point = if let Some(p) = corners[i] {
                p
            } else {
                let r = source
                    .radius_at_angle_mm(angle)
                    .ok_or(MeshError::UnresolvedBoundary)?;
                [r * angle.portable_cos(), r * angle.portable_sin()]
            };
            let z = if floor {
                top + source.cone().tip_depth_mm
                    - point[0].portable_hypot(point[1]) * source.cone().depth_per_radius
            } else {
                top
            };
            positions.push([point[0], point[1], z]);
        }
        validate_ring(&positions, source, options)?;
        Ok(positions.into_iter().map(|p| self.push(p)).collect())
    }
}
