//! Measures and descriptors of the actual raw affine source image.
//! These APIs do not establish finite boundary incidence or Solid admission.
use super::*;
use crate::geometry::internal_thread_treatment_body::interval::Interval as B;
use crate::geometry::source_face_measure::positive;
use multicad_rust_kernel::authored_internal_treatment::Enclosure;
use multicad_rust_kernel::wood_insert::{FaceKey, FaceProperties, FacePropertyOptions};
#[cfg(test)]
mod tests;
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WoodInsertFaceProperties {
    pub local: FaceProperties,
    pub world: FaceProperties,
    pub area_jacobian: Enclosure,
    pub local_centroid_reweighting_error_mm: [f64; 3],
}
fn refused() -> GeometryError {
    error(
        "Rust wood insert face measure",
        GeometryErrorKind::Unsupported,
        "Actual world face measure exceeds its precision bound",
    )
}
impl WoodInsertMaterial {
    /// Bounded world measures; the literal frame is never normalized.
    pub fn world_face_properties(
        &self,
        face: FaceKey,
        options: FacePropertyOptions,
        policy: GeometryPolicy,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<super::WoodInsertFaceProperties, GeometryError> {
        checkpoint(interrupt)?;
        policy.validate()?;
        policy.length(self.source.minimum_wall_mm())?;
        let r = self.source.outer().maximum_radius_mm();
        self.frame
            .checked_metric_extent([r, r, self.source.full_length_mm()], policy)?;
        transport(&self.source, self.frame, face, options, policy, interrupt)
    }
}
fn planar(face: FaceKey) -> bool {
    matches!(face, FaceKey::StartCap | FaceKey::EndCap)
}
fn coordinates(material: &Insert, face: FaceKey) -> [B; 3] {
    let r = material.outer().maximum_radius_mm();
    let z = match face {
        FaceKey::StartCap => B::ZERO,
        FaceKey::EndCap => B::point(material.full_length_mm()),
        _ => B::new(0., material.full_length_mm()),
    };
    [B::new(-r, r), B::new(-r, r), z]
}
fn area_upper(material: &Insert, face: FaceKey) -> Result<f64, GeometryError> {
    let pi = B::new(
        std::f64::consts::PI.next_down(),
        std::f64::consts::PI.next_up(),
    );
    let thread = match face {
        FaceKey::StartCap | FaceKey::EndCap => {
            return Ok(positive(pi * B::point(material.outer().maximum_radius_mm()).square())?.hi);
        }
        FaceKey::OuterThread { .. } => material.outer(),
        FaceKey::InnerThread { .. } => material.inner(),
    };
    let spec = thread.spec();
    let gain = B::point(2.) * pi / (B::point(spec.pitch_mm) * B::point(f64::from(spec.starts)));
    // Complete actual wall is a conservative upper bound for any one role.
    // This bound supplies numerical budgets, never face ownership or geometry.
    let mut sum = B::ZERO;
    for b in thread.bands() {
        let dq = B::point(b.phase_mm[1]) - B::point(b.phase_mm[0]);
        let m = (B::point(b.radius_mm[1]) - B::point(b.radius_mm[0])) / dq;
        let r = B::point(b.radius_mm[0].max(b.radius_mm[1]));
        sum = sum
            + B::point(spec.length_mm)
                * dq
                * ((r * gain).square() * (B::ONE + m.square()) + m.square()).sqrt();
    }
    Ok(positive(sum * B::point(f64::from(spec.starts)))?.hi)
}
fn transport(
    material: &Insert,
    placement: PlaneFrame,
    face: FaceKey,
    options: FacePropertyOptions,
    policy: GeometryPolicy,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<WoodInsertFaceProperties, GeometryError> {
    use crate::geometry::source_face_measure::{self as shared, Measures, Request};
    checkpoint(interrupt)?;
    if material.face_index(face).is_none() {
        return Err(input("Face does not belong to the actual material source"));
    }
    let r = material.outer().maximum_radius_mm();
    let receipt = shared::transport(
        Request {
            frame: placement,
            options,
            planar: planar(face),
            domain: coordinates(material, face),
            whole_domain: [
                B::new(-r, r),
                B::new(-r, r),
                B::new(0., material.full_length_mm()),
            ],
            area_upper_mm2: area_upper(material, face)?,
        },
        policy,
        |options| {
            let local = material
                .local_face_properties(face, options, interrupt)
                .map_err(failure)?;
            let values = Measures {
                area_mm2: B::from(local.area_mm2),
                first_moment_mm3: local.first_moment_mm3.map(B::from),
                centroid_mm: local.centroid_mm.map(B::from),
            };
            Ok((local, values))
        },
    )?;
    let result = WoodInsertFaceProperties {
        local: receipt.local,
        world: FaceProperties {
            face,
            area_mm2: receipt.world.area_mm2.source(),
            first_moment_mm3: receipt.world.first_moment_mm3.map(B::source),
            centroid_mm: receipt.world.centroid_mm.map(B::source),
        },
        area_jacobian: receipt.area_jacobian.source(),
        local_centroid_reweighting_error_mm: receipt.local_centroid_reweighting_error_mm,
    };
    checkpoint(interrupt)?;
    Ok(result)
}

fn up(a: f64, b: f64) -> f64 {
    (a + b).next_up()
}
fn midpoint(values: [B; 3]) -> Result<([f64; 3], f64), GeometryError> {
    if values.into_iter().any(|v| !v.valid()) {
        return Err(refused());
    }
    let point = values.map(|v| v.lo.midpoint(v.hi));
    let e = values
        .into_iter()
        .zip(point)
        .fold(0., |s, (b, p)| up(s, (B::point(p) - b).magnitude()));
    Ok((point, e))
}
fn center(values: [B; 3], policy: GeometryPolicy) -> Result<[f64; 3], GeometryError> {
    let (p, e) = midpoint(values)?;
    if e > policy.linear_tolerance_mm / 8. {
        return Err(refused());
    }
    Ok(p)
}
type PlaneSupport = ([f64; 3], [f64; 3]);
fn plane(
    material: &Insert,
    frame: PlaneFrame,
    face: FaceKey,
    policy: GeometryPolicy,
) -> Result<Option<PlaneSupport>, GeometryError> {
    let z = match face {
        FaceKey::StartCap => 0.,
        FaceKey::EndCap => material.full_length_mm(),
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
    let point = std::array::from_fn(|i| {
        B::point(frame.origin[i]) + B::point(frame.normal()[i]) * B::point(z)
    });
    let (point, point_error) = midpoint(point)?;
    let extent = B::point(material.outer().maximum_radius_mm())
        * frame
            .u
            .into_iter()
            .chain(frame.v)
            .fold(B::ZERO, |s, v| s + B::point(v.abs()));
    let error =
        B::point(point_error) * (B::ONE + B::point(normal_error)) + B::point(normal_error) * extent;
    if normal_error > 1e-12 || !error.valid() || error.hi > policy.linear_tolerance_mm / 8. {
        return Err(refused());
    }
    policy.point(point)?;
    Ok(Some((point, normal)))
}
impl WoodInsertMaterial {
    /// Source-role descriptors only; no finite boundary incidence is implied.
    pub fn faces(
        &self,
        policy: GeometryPolicy,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Vec<Face>, GeometryError> {
        checkpoint(interrupt)?;
        policy.validate()?;
        let r = self.source.outer().maximum_radius_mm();
        let l = self.source.full_length_mm();
        policy.length(self.source.minimum_wall_mm())?;
        self.frame.checked_metric_extent([r, r, l], policy)?;
        let scale = r.max(l).max(1.);
        let extent = self
            .frame
            .origin
            .into_iter()
            .map(f64::abs)
            .fold(1., f64::max)
            + 2. * (2. * r + l);
        let area_tolerance_mm2 = policy.linear_tolerance_mm * scale / 8.;
        let options = FacePropertyOptions {
            area_tolerance_mm2,
            centroid_tolerance_mm: policy.linear_tolerance_mm / 8.,
            first_moment_tolerance_mm3: area_tolerance_mm2 * (extent + scale),
            ..Default::default()
        };
        let mut result = Vec::with_capacity(self.source.face_count() as usize);
        for index in 0..self.source.face_count() {
            let key = self
                .source
                .face_key(index)
                .expect("sealed source face ordinal");
            let receipt = self.world_face_properties(key, options, policy, interrupt)?;
            let center_mm = center(receipt.world.centroid_mm.map(B::from), policy)?;
            let face = Face {
                index,
                plane: plane(&self.source, self.frame, key, policy)?,
                center_mm,
                area_mm2: receipt
                    .world
                    .area_mm2
                    .lower
                    .midpoint(receipt.world.area_mm2.upper),
            };
            policy.point(face.center_mm)?;
            result.push(face);
        }
        checkpoint(interrupt)?;
        Ok(result)
    }
}
