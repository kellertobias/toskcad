//! Display/selectable descriptors from the actual finite source, never mesh areas.
use super::interval::Interval as B;
use super::*;
use source::{FaceKey, FacePropertyOptions};
#[cfg(test)]
mod tests;
fn refused() -> GeometryError {
    unsupported("Actual internal-treatment face descriptor exceeds its precision bound")
}
fn midpoint(values: [B; 3]) -> Result<([f64; 3], f64), GeometryError> {
    if values.into_iter().any(|v| !v.valid()) {
        return Err(refused());
    }
    let points = values.map(|b| b.lo.midpoint(b.hi));
    let error = values
        .into_iter()
        .zip(points)
        .fold(0., |sum, (b, p)| up(sum, (B::point(p) - b).magnitude()));
    Ok((points, error))
}
fn centroid(values: [B; 3], policy: GeometryPolicy) -> Result<[f64; 3], GeometryError> {
    let (point, error) = midpoint(values)?;
    if error > policy.linear_tolerance_mm / 8. {
        return Err(refused());
    }
    Ok(point)
}
type PlaneSupport = ([f64; 3], [f64; 3]);
fn plane(
    material: &Material,
    frame: PlaneFrame,
    face: FaceKey,
    policy: GeometryPolicy,
) -> Result<Option<PlaneSupport>, GeometryError> {
    let z = match face {
        FaceKey::StartCap => B::ZERO,
        FaceKey::EndCap => B::from(material.cavity().full_length_bounds_mm()),
        FaceKey::TerminalShoulder => material
            .cavity()
            .terminal_station_mm()
            .map(B::point)
            .ok_or_else(refused)?,
        _ => return Ok(None),
    };
    let u = frame.u.map(B::point);
    let v = frame.v.map(B::point);
    let cross = [(1, 2), (2, 0), (0, 1)].map(|(j, k)| u[j] * v[k] - u[k] * v[j]);
    let length = cross
        .into_iter()
        .fold(B::ZERO, |s, v| s + v.square())
        .sqrt();
    if !length.valid() || length.lo <= 0. {
        return Err(refused());
    }
    let orientation = if face == FaceKey::StartCap { -1. } else { 1. };
    let unit = cross.map(|v| B::point(orientation) * v / length);
    let (normal, normal_error) = midpoint(unit)?;
    let point =
        std::array::from_fn(|i| B::point(frame.origin[i]) + B::point(frame.normal()[i]) * z);
    let (point, point_error) = midpoint(point)?;
    // The support plane spans the literal U,V columns. The rounded axial
    // column only locates its station; it is never substituted for its normal.
    let extent = B::point(radius_upper(material))
        * (frame
            .u
            .into_iter()
            .chain(frame.v)
            .fold(B::ZERO, |s, v| s + B::point(v.abs())));
    // The published normal is itself rounded: its component norm is bounded
    // by 1+normal_error when charging the support-point error.
    let error =
        B::point(point_error) * (B::ONE + B::point(normal_error)) + B::point(normal_error) * extent;
    if normal_error > 1e-12 || !error.valid() || error.hi > policy.linear_tolerance_mm / 8. {
        return Err(refused());
    }
    policy.point(point)?;
    Ok(Some((point, normal)))
}
impl InternalThreadTreatmentBody {
    pub(in crate::geometry) fn faces(
        &self,
        frame: PlaneFrame,
        policy: GeometryPolicy,
    ) -> Result<Vec<Face>, GeometryError> {
        self.validate_frame(frame, policy)?;
        let boundary = self.qualified_boundary()?;
        let scale = radius_upper(&self.material)
            .max(self.material.cavity().full_length_bounds_mm().upper)
            .max(1.);
        let extent = frame.origin.into_iter().map(f64::abs).fold(1., f64::max)
            + 2. * (2. * radius_upper(&self.material)
                + self.material.cavity().full_length_bounds_mm().upper);
        // Units: linear error times a source length is an area budget; the
        // first-moment budget also includes translation times area uncertainty.
        let area_tolerance_mm2 = policy.linear_tolerance_mm * scale / 8.;
        let options = FacePropertyOptions {
            area_tolerance_mm2,
            centroid_tolerance_mm: policy.linear_tolerance_mm / 8.,
            first_moment_tolerance_mm3: area_tolerance_mm2 * (extent + scale),
            ..Default::default()
        };
        // Count comes from the resource-qualified immutable source graph.
        let mut result = Vec::with_capacity(boundary.face_count() as usize);
        for index in 0..boundary.face_count() {
            let key = self.material.face_key(index).ok_or_else(refused)?;
            let receipt = self.world_face_properties(frame, key, options, policy, None)?;
            let center_mm = centroid(receipt.world.centroid_mm.map(B::from), policy)?;
            let descriptor = Face {
                index,
                plane: plane(&self.material, frame, key, policy)?,
                center_mm,
                area_mm2: receipt
                    .world
                    .area_mm2
                    .lower
                    .midpoint(receipt.world.area_mm2.upper),
            };
            policy.point(descriptor.center_mm)?;
            result.push(descriptor);
        }
        Ok(result)
    }
}
