//! Authored exact cylindrical stock minus the unchanged finite rounded cutter.
//! This is a separate closure; a through receiver never changes its domain.
use super::rounded_thread_body::{RoundedThreadBody, RoundedThreadMassProperties};
use super::*;
use crate::rust_prism_topology::position;
use multicad_core::references::{
    AnalyticGeometry, Origin, TopologyElement, TopologyKey, TopologyKind,
};
use multicad_rust_kernel::{
    Rational,
    rounded_thread::{BlindBore, Spec},
};
use std::collections::BTreeSet;
mod mesh;
mod repaired;
mod stock;
mod topology;
pub use mesh::{RoundedBlindMeshBudget, RoundedBlindMeshReceipt};
type Point = [f64; 3];
const OPERATION: &str = "Rust rounded blind receiver";

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RoundedBlindBoreData {
    pub bore: BlindBore,
    pub frame: PlaneFrame,
}
pub(crate) struct RoundedBlindBoreBody {
    pub(crate) geometry: BlindBore,
    pub(super) feature: EntityId,
    inner: Rc<RoundedThreadBody>,
    inherited: Option<(TopologySnapshot, PlaneFrame, bool)>,
    import_identity: Option<Origin>,
    import_source: Option<multicad_rust_kernel::step_read::RepairedRoundedBlindCandidate>,
}
fn failure(kind: GeometryErrorKind, message: impl Into<String>) -> GeometryError {
    error(OPERATION, kind, message)
}
impl Solid {
    pub fn rust_rounded_blind_bore(
        spec: Spec,
        radius: Rational,
        height: Rational,
        feature: EntityId,
        frame: PlaneFrame,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        let tool = Self::rust_rounded_thread(spec, feature, frame, policy)?;
        let SolidBody::RoundedThread(inner, _) = tool.body else {
            unreachable!("rounded source factory")
        };
        let geometry = BlindBore::new(inner.geometry, radius, height).map_err(|e| {
            failure(
                GeometryErrorKind::Unsupported,
                format!("Actual blind source refused: {e:?}"),
            )
        })?;
        policy.length(
            (geometry.stock_radius_bounds_mm()[0] - geometry.cutter_radial_bounds_mm()[1])
                .next_down(),
        )?;
        policy.length(geometry.floor_thickness_bounds_mm()[0])?;
        let body = RoundedBlindBoreBody {
            geometry,
            feature,
            inner,
            inherited: None,
            import_identity: None,
            import_source: None,
        };
        body.validate_frame(frame, policy)?;
        Ok(Self {
            body: SolidBody::RoundedBlindBore(Rc::new(body), frame),
            policy,
            provenance: vec![],
            rust_origins: vec![],
            retained_placement: None,
        })
    }
    pub(crate) fn rounded_blind_bore_body(&self) -> Option<(&RoundedBlindBoreBody, PlaneFrame)> {
        match &self.body {
            SolidBody::RoundedBlindBore(body, frame) => Some((body, *frame)),
            _ => None,
        }
    }
    pub fn rust_rounded_blind_bore_data(&self) -> Option<RoundedBlindBoreData> {
        self.rounded_blind_bore_body()
            .map(|(b, frame)| RoundedBlindBoreData {
                bore: b.geometry,
                frame,
            })
    }
    /// Point membership in the authored finite cutter difference, at this
    /// solid's numerical policy. This is not a whole-box Boolean certificate.
    pub fn rust_rounded_blind_bore_contains(&self, point: Point) -> Option<bool> {
        self.rounded_blind_bore_body()
            .map(|(b, f)| b.contains(f, point))
    }
    /// Source-known correspondence check against this retained authored bore.
    /// This does not import arbitrary files or certify exact raw REAL mass.
    pub fn verify_source_known_rounded_blind_step(
        &self,
        input: &str,
    ) -> Result<f64, GeometryError> {
        let (body, frame) = self.rounded_blind_bore_body().ok_or_else(|| {
            failure(
                GeometryErrorKind::Unsupported,
                "No authored rounded blind source",
            )
        })?;
        multicad_rust_kernel::step_read::verify_rounded_thread_blind_bore(
            input,
            &body.geometry,
            frame.origin,
            frame.u,
            frame.v,
            self.policy.linear_tolerance_mm,
        )
        .map_err(|e| {
            failure(
                GeometryErrorKind::Unsupported,
                format!("Source-known blind STEP refused: {e:?}"),
            )
        })
    }
    /// Earn the TOTAL source-known positional allowance at the requested
    /// tolerance. Diagnostic spline flux supplies no physical mass authority.
    pub fn verify_source_known_rounded_blind_position_bound(
        &self,
        input: &str,
        tolerance: f64,
    ) -> Result<f64, GeometryError> {
        if !tolerance.is_finite() || tolerance <= 0. {
            return Err(failure(
                GeometryErrorKind::Input,
                "Invalid blind correspondence tolerance",
            ));
        }
        let (body, frame) = self.rounded_blind_bore_body().ok_or_else(|| {
            failure(
                GeometryErrorKind::Unsupported,
                "No canonical rounded blind body",
            )
        })?;
        multicad_rust_kernel::step_read::verify_rounded_blind_correspondence(
            input,
            &body.geometry,
            frame.origin,
            frame.u,
            frame.v,
            tolerance.min(self.policy.linear_tolerance_mm),
        )
        .map(|p| p.maximum_position_error_mm())
        .map_err(|e| {
            failure(
                if e == multicad_rust_kernel::step_read::ReadError::Limit {
                    GeometryErrorKind::Limit
                } else {
                    GeometryErrorKind::Unsupported
                },
                format!("Complete blind correspondence refused: {e:?}"),
            )
        })
    }
    pub fn rust_rounded_blind_bore_mass_properties(
        &self,
    ) -> Result<Option<RoundedThreadMassProperties>, GeometryError> {
        self.rounded_blind_bore_body()
            .map(|(b, f)| b.mass(f))
            .transpose()
    }
}
impl RoundedBlindBoreBody {
    fn extents(&self) -> [f64; 3] {
        let r = self.geometry.stock_radius_bounds_mm()[1];
        [r, r, self.geometry.stock_height_bounds_mm()[1]]
    }
    pub(super) fn validate_frame(
        &self,
        frame: PlaneFrame,
        policy: GeometryPolicy,
    ) -> Result<(), GeometryError> {
        frame.checked_metric_extent(self.extents(), policy)?;
        let bounds = self.bounds(frame)?;
        policy.point(bounds.min)?;
        policy.point(bounds.max)?;
        Ok(())
    }
    pub(super) fn bounds(&self, frame: PlaneFrame) -> Result<Bounds, GeometryError> {
        let n = frame.normal();
        let mut out = Bounds {
            min: [0.; 3],
            max: [0.; 3],
        };
        for axis in 0..3 {
            for (sign, target) in [(-1., &mut out.min), (1., &mut out.max)] {
                let high = self
                    .geometry
                    .support_bounds(
                        [sign * frame.u[axis], sign * frame.v[axis], sign * n[axis]],
                        self.inner.tolerance() / 8.,
                    )
                    .map_err(|e| {
                        failure(
                            GeometryErrorKind::Limit,
                            format!("Actual stock support refused: {e:?}"),
                        )
                    })?[1];
                let value = frame.origin[axis] + sign * high;
                target[axis] = if sign < 0. {
                    value.next_down()
                } else {
                    value.next_up()
                };
            }
        }
        Ok(out)
    }
    pub(super) fn contains(&self, frame: PlaneFrame, point: Point) -> bool {
        let n = frame.normal();
        let rows = [cross(frame.v, n), cross(n, frame.u), n];
        let det = dot(frame.u, rows[0]);
        let delta = std::array::from_fn(|i| point[i] - frame.origin[i]);
        let [x, y, z] = rows.map(|r| dot(r, delta) / det);
        let tolerance = self.inner.tolerance();
        let rho = x.hypot(y);
        let height = self.geometry.stock_height_mm().as_f64();
        if !rho.is_finite()
            || !z.is_finite()
            || z < -tolerance
            || z > height + tolerance
            || rho > self.geometry.stock_radius_mm().as_f64() + tolerance
        {
            return false;
        }
        z >= self.geometry.cutter().spec().length_mm
            || rho >= self.geometry.cutter().radius(y.atan2(x), z.max(0.)) - tolerance
    }
    fn local_volume(&self) -> f64 {
        std::f64::consts::PI
            * self.geometry.stock_radius_mm().as_f64().powi(2)
            * self.geometry.stock_height_mm().as_f64()
            - self.geometry.cutter().volume_mm3()
    }
    pub(super) fn volume(&self, frame: PlaneFrame) -> f64 {
        self.local_volume() * dot(frame.u, cross(frame.v, frame.normal()))
    }
    fn mass(&self, frame: PlaneFrame) -> Result<RoundedThreadMassProperties, GeometryError> {
        let m = self.inner.properties()?;
        let r = self.geometry.stock_radius_mm().as_f64();
        let h = self.geometry.stock_height_mm().as_f64();
        let stock = std::f64::consts::PI * r * r * h;
        let hole = self.geometry.cutter().volume_mm3();
        let volume = self.local_volume();
        if !volume.is_finite() || volume <= 0. {
            return Err(failure(
                GeometryErrorKind::InvalidResult,
                "Blind material mass is unresolved",
            ));
        }
        let center: Point = std::array::from_fn(|i| {
            ((if i == 2 { stock * h / 2. } else { 0. }) - hole * m.centroid_mm[i]) / volume
        });
        let second: [[f64; 3]; 3] = std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                let outer = if i != j {
                    0.
                } else if i == 2 {
                    stock * h * h / 3.
                } else {
                    stock * r * r / 4.
                };
                outer
                    - m.central_second_moment_mm5[i][j]
                    - hole * m.centroid_mm[i] * m.centroid_mm[j]
                    - volume * center[i] * center[j]
            })
        });
        let axes = [frame.u, frame.v, frame.normal()];
        let det = dot(axes[0], cross(axes[1], axes[2]));
        let covariance: [[f64; 3]; 3] = std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                det * (0..3)
                    .flat_map(|a| (0..3).map(move |b| axes[a][i] * second[a][b] * axes[b][j]))
                    .sum::<f64>()
            })
        });
        let trace = (0..3).map(|i| covariance[i][i]).sum::<f64>();
        Ok(RoundedThreadMassProperties {
            volume_mm3: volume * det,
            centroid_mm: position(frame, center),
            inertia_mm5: std::array::from_fn(|i| {
                std::array::from_fn(|j| {
                    if i == j {
                        trace - covariance[i][j]
                    } else {
                        -covariance[i][j]
                    }
                })
            }),
        })
    }
    pub(super) fn faces(&self, frame: PlaneFrame) -> Result<Vec<Face>, GeometryError> {
        let mut faces = self.inner.faces(frame)?;
        let m = self.inner.properties()?;
        let r = self.geometry.stock_radius_mm().as_f64();
        let h = self.geometry.stock_height_mm().as_f64();
        let area = std::f64::consts::PI * r * r;
        let (center, hole) = m.faces[0];
        let remaining = area - hole;
        faces[0].center_mm = position(
            frame,
            [
                -hole * center[0] / remaining,
                -hole * center[1] / remaining,
                0.,
            ],
        );
        faces[0].area_mm2 = remaining;
        // Original finite tip section, with the opposite occupied-side normal.
        if let Some((origin, normal)) = &mut faces[1].plane {
            *normal = normal.map(|v| -v);
            let _ = origin;
        }
        faces.push(Face {
            index: 7,
            plane: None,
            center_mm: position(frame, [0., 0., h / 2.]),
            area_mm2: std::f64::consts::TAU * r * h,
        });
        faces.push(Face {
            index: 8,
            plane: Some((position(frame, [0., 0., h]), frame.normal())),
            center_mm: position(frame, [0., 0., h]),
            area_mm2: area,
        });
        Ok(faces)
    }
    pub(super) fn append_step(
        &self,
        writer: &mut multicad_rust_kernel::step::Writer,
        frame: PlaneFrame,
        policy: GeometryPolicy,
    ) -> Result<(), multicad_rust_kernel::step::StepError> {
        writer.append_rounded_thread_blind_bore(
            &self.geometry,
            frame.origin,
            frame.u,
            frame.v,
            policy.linear_tolerance_mm,
        )
    }
}
