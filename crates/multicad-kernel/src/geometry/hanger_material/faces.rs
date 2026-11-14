//! Measured face descriptors of the actual placed continuous hanger source.
use super::*;
use crate::geometry::source_face_measure::{self as shared, Measures, Request};
use multicad_rust_kernel::authored_hanger_bolt::{FaceKey, FaceProperties, FacePropertyOptions};
use multicad_rust_kernel::authored_internal_treatment::Enclosure;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HangerFaceProperties {
    pub local: FaceProperties,
    pub world: FaceProperties,
    pub area_jacobian: Enclosure,
    pub local_centroid_reweighting_error_mm: [f64; 3],
}
fn pair([lo, hi]: [f64; 2]) -> I {
    I::new(lo, hi)
}
fn domain(source: Hanger, face: FaceKey) -> [I; 3] {
    let mut bounds = source.bounds_mm().map(pair);
    let [metric, plain, wood] = source.sections();
    bounds[2] = match face {
        FaceKey::MetricTip => pair(metric.start_mm),
        FaceKey::MetricShoulder => pair(metric.end_mm),
        FaceKey::WoodShoulder => pair(wood.start_mm),
        FaceKey::WoodTip => pair(wood.end_mm),
        FaceKey::MetricThread { .. } => I::new(metric.start_mm[0], metric.end_mm[1]),
        FaceKey::WoodThread { .. } => I::new(wood.start_mm[0], wood.end_mm[1]),
        FaceKey::PlainWall => I::new(plain.start_mm[0], plain.end_mm[1]),
    };
    bounds
}
fn planar(face: FaceKey) -> bool {
    matches!(
        face,
        FaceKey::MetricTip | FaceKey::MetricShoulder | FaceKey::WoodShoulder | FaceKey::WoodTip
    )
}
fn area_upper(source: Hanger, face: FaceKey) -> f64 {
    let tau = I::new(
        std::f64::consts::TAU.next_down(),
        std::f64::consts::TAU.next_up(),
    );
    let wall = |radius: f64, length: f64, lead: I, slope: I| {
        // For r(q), |r'|<=m and k=tau/lead, the complete wall is bounded
        // by L*tau*sqrt(R²(1+m²)+(m/k)²) <= L*(tau*R*(1+m)+m*lead).
        // This is only a conservative numerical budget, never face geometry.
        (I::point(length) * (tau * I::point(radius) * (I::ONE + slope) + slope * lead)).hi
    };
    match face {
        FaceKey::PlainWall => {
            (tau * I::point(source.spec().plain_radius_mm)
                * I::point(source.spec().plain_length_mm))
            .hi
        }
        FaceKey::MetricThread { .. } => {
            let s = source.metric().spec();
            wall(
                s.major_radius_mm,
                s.length_mm,
                I::point(s.pitch_mm) * I::point(f64::from(s.starts)),
                I::point(3.).sqrt(),
            )
        }
        FaceKey::WoodThread { .. } => {
            let s = source.wood().spec();
            let mut slope = 0_f64;
            for b in source.wood().bands() {
                let dr = I::point(b.radius_mm[1]) - I::point(b.radius_mm[0]);
                let dq = I::point(b.phase_mm[1]) - I::point(b.phase_mm[0]);
                slope = slope.max((dr / dq).magnitude());
            }
            wall(
                source.wood().maximum_radius_mm(),
                s.length_mm,
                I::point(s.pitch_mm) * I::point(f64::from(s.starts)),
                I::point(slope),
            )
        }
        _ => (tau * I::point(source.wood().maximum_radius_mm()).square() / I::point(2.)).hi,
    }
}
impl HangerMaterial {
    /// Actual affine surface measures, with explicit normal-dependent area
    /// scaling and centroid reweighting. The frame is never normalized.
    pub fn world_face_properties(
        &self,
        face: FaceKey,
        options: FacePropertyOptions,
        policy: GeometryPolicy,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<HangerFaceProperties, GeometryError> {
        Self::new(self.source, self.frame, policy, interrupt)?;
        if self.source.face_index(face).is_none() {
            return Err(input("Face does not belong to the actual hanger source"));
        }
        let receipt = shared::transport(
            Request {
                frame: self.frame,
                options,
                planar: planar(face),
                domain: domain(self.source, face),
                whole_domain: self.source.bounds_mm().map(pair),
                area_upper_mm2: area_upper(self.source, face),
            },
            policy,
            |options| {
                let local = self
                    .source
                    .local_face_properties(face, options, interrupt)
                    .map_err(failure)?;
                let values = Measures {
                    area_mm2: I::from(local.area_mm2),
                    first_moment_mm3: local.first_moment_mm3.map(I::from),
                    centroid_mm: local.centroid_mm.map(I::from),
                };
                Ok((local, values))
            },
        )?;
        checkpoint(interrupt)?;
        Ok(HangerFaceProperties {
            local: receipt.local,
            world: FaceProperties {
                face,
                area_mm2: receipt.world.area_mm2.source(),
                first_moment_mm3: receipt.world.first_moment_mm3.map(I::source),
                centroid_mm: receipt.world.centroid_mm.map(I::source),
            },
            area_jacobian: receipt.area_jacobian.source(),
            local_centroid_reweighting_error_mm: receipt.local_centroid_reweighting_error_mm,
        })
    }
    pub fn faces(
        &self,
        policy: GeometryPolicy,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Vec<Face>, GeometryError> {
        Self::new(self.source, self.frame, policy, interrupt)?;
        let bounds = self.source.bounds_mm();
        let scale = bounds
            .into_iter()
            .flatten()
            .map(f64::abs)
            .fold(1., f64::max);
        let world = self.bounds();
        let extent = world
            .min
            .into_iter()
            .chain(world.max)
            .map(f64::abs)
            .fold(1., f64::max);
        let tolerance = policy.linear_tolerance_mm * scale / 8.;
        let options = FacePropertyOptions {
            area_tolerance_mm2: tolerance,
            centroid_tolerance_mm: policy.linear_tolerance_mm / 8.,
            first_moment_tolerance_mm3: tolerance * (extent + scale),
            ..Default::default()
        };
        let mut result = Vec::with_capacity(self.source.face_count());
        for index in 0..self.source.face_count() {
            let key = self
                .source
                .face_key(index)
                .expect("actual source face ordinal");
            let receipt = self.world_face_properties(key, options, policy, interrupt)?;
            let plane = if planar(key) {
                Some(shared::plane(
                    self.frame,
                    domain(self.source, key)[2],
                    self.source.wood().maximum_radius_mm(),
                    key == FaceKey::WoodTip,
                    policy,
                )?)
            } else {
                None
            };
            let center_mm = shared::center(receipt.world.centroid_mm.map(I::from), policy)?;
            policy.point(center_mm)?;
            result.push(Face {
                index: index as u32,
                plane,
                center_mm,
                area_mm2: receipt
                    .world
                    .area_mm2
                    .lower
                    .midpoint(receipt.world.area_mm2.upper),
            });
        }
        checkpoint(interrupt)?;
        Ok(result)
    }
}
