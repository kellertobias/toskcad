//! Typed retained Phillips source with bounded native material receipts.
//! Physical STEP exchange is not admitted by this source constructor.
use super::*;
use multicad_rust_kernel::phillips_drive::{self as source, PhillipsScrew, ScrewFace};
pub(super) mod display;
mod properties;
pub(super) mod queries;
pub(super) mod topology;
pub(super) use properties::plane;

#[derive(Clone, Copy, Debug)]
pub struct PhillipsMeshVolumePair {
    pub source_world_mm3: [f64; 2],
    pub frame: PlaneFrame,
    pub linear_deflection_mm: f64,
    pub angular_deflection_rad: f64,
    pub triangle_limit: usize,
    pub literal_world_mm3: [f64; 2],
    pub source_to_literal_bound_mm3: f64,
    pub literal_vertex_transfer_bound_mm3: f64,
}
pub(super) struct Body {
    source: PhillipsScrew,
    policy: GeometryPolicy,
    properties: Vec<source::ScrewFaceProperties>,
    volume: source::SourceVolumeBounds,
    last_mesh_pair: std::cell::RefCell<Option<PhillipsMeshVolumePair>>,
}
fn failure(cause: source::Error) -> GeometryError {
    let kind = match cause {
        source::Error::Cancelled => GeometryErrorKind::Cancelled,
        source::Error::Limit => GeometryErrorKind::Limit,
        source::Error::Invalid => GeometryErrorKind::Input,
        source::Error::Unresolved => GeometryErrorKind::Unsupported,
    };
    error(
        OPERATION,
        kind,
        format!("Actual retained Phillips source proof failed: {cause:?}"),
    )
}
impl Body {
    fn new(
        source: PhillipsScrew,
        policy: GeometryPolicy,
        interrupt: &dyn Fn() -> bool,
    ) -> Result<Self, GeometryError> {
        policy.validate()?;
        policy.length(source.spec().minimum_wall_mm)?;
        if source.profile().transport_error_mm() != 0. {
            return Err(error(
                OPERATION,
                GeometryErrorKind::Unsupported,
                "Transported Phillips face measures require an actual carrier metric certificate",
            ));
        }
        if source.profile().rounding_error_mm() > policy.linear_tolerance_mm / 8. {
            return Err(input(
                "Phillips head source precision exceeds the requested policy",
            ));
        }
        let extent = source
            .profile()
            .support([1., 0., 0.])
            .max(source.thread().spec().length_mm)
            .max(1.);
        let properties = source
            .face_properties_with_cancel(
                source::FacePropertyOptions {
                    area_tolerance_mm2: policy.linear_tolerance_mm * extent / 16.,
                    first_moment_tolerance_mm3: policy.linear_tolerance_mm * extent * extent / 16.,
                    centroid_tolerance_mm: policy.linear_tolerance_mm / 16.,
                    maximum_cells: 65536,
                },
                interrupt,
            )
            .map_err(failure)?;
        let volume = source.source_volume_bounds().map_err(failure)?;
        if interrupt() {
            return Err(failure(source::Error::Cancelled));
        }
        Ok(Self {
            source,
            policy,
            properties,
            volume,
            last_mesh_pair: Default::default(),
        })
    }
    pub(super) fn bounds(&self, frame: PlaneFrame) -> Bounds {
        properties::bounds(&self.source, frame)
    }
    pub(super) fn source(&self) -> &PhillipsScrew {
        &self.source
    }
    pub(super) fn policy(&self) -> GeometryPolicy {
        self.policy
    }
    pub(super) fn extend_properties(
        &self,
        faces: &mut Vec<(Point, f64)>,
    ) -> Result<(), GeometryError> {
        *faces = self
            .properties
            .iter()
            .map(|p| {
                (
                    p.centroid_mm.map(|b| b[0].midpoint(b[1])),
                    p.area_mm2[0].midpoint(p.area_mm2[1]),
                )
            })
            .collect();
        Ok(())
    }
    pub(super) fn volume_mm3(&self) -> f64 {
        let [a, b] = self.volume.volume_mm3();
        a.midpoint(b)
    }
    pub(super) fn mesh(
        &self,
        frame: PlaneFrame,
        options: multicad_rust_kernel::display::Tessellation,
    ) -> Result<Mesh, GeometryError> {
        self.mesh_with_interrupt(frame, options, &|| false)
    }
    fn mesh_with_interrupt(
        &self,
        frame: PlaneFrame,
        options: multicad_rust_kernel::display::Tessellation,
        interrupt: &dyn Fn() -> bool,
    ) -> Result<Mesh, GeometryError> {
        if interrupt() {
            return Err(failure(source::Error::Cancelled));
        }
        let local_pair = std::cell::Cell::new(None);
        let mesh = super::super::entry_chamfer_body::placement::qualified_source_mesh_with_small_normal_reserve(
            frame,
            options,
            |quality| {
                let candidate =
                    helical_thread_mesh::tessellate_phillips_screw(&self.source, quality, interrupt)
                    .map_err(|cause| match cause {
                        helical_thread_mesh::PhillipsMeshError::Cancelled => {
                            failure(source::Error::Cancelled)
                        }
                        helical_thread_mesh::PhillipsMeshError::Geometry(cause) => {
                            crate::rust_display::mesh(Err(cause)).unwrap_err()
                        }
                    })?;
                // Keep source interval and literal mesh volume distinct. This
                // candidate has its own complete facet/ownership/closure proof;
                // the placement helper additionally proves literal embedding.
                // The producer borrows this exact immutable source. Retain
                // its owned pairing directly; volume equality is not a source
                // identity test and no detached candidate is admitted here.
                local_pair.set(Some((candidate.source_volume(),candidate.literal_volume_mm3(),candidate.volume_difference_bound_mm3())));
                let m = candidate.mesh();
                Ok(Mesh {
                    vertices: m.vertices.clone(),
                    indices: m.indices.clone(),
                    face_indices: m.face_indices.clone(),
                })
            },
        )?;
        if interrupt() {
            return Err(failure(source::Error::Cancelled));
        }
        // Both successful transfer branches invoke the local producer exactly
        // once. A future mesh cache must retain this source receipt as well.
        let (volume, literal, difference) = local_pair
            .get()
            .expect("successful source transfer invoked its owned producer");
        let pair = properties::mesh_pair(volume, frame, &mesh, literal, difference, options)?;
        if interrupt() {
            return Err(failure(source::Error::Cancelled));
        }
        *self.last_mesh_pair.borrow_mut() = Some(pair);
        Ok(mesh)
    }
    pub(super) fn volume_error_bound(
        &self,
        _options: multicad_rust_kernel::display::Tessellation,
    ) -> f64 {
        // This source has directed source/literal pairing, not an exact scalar
        // tube formula. Manufacturing uses its independently earned source gate.
        f64::INFINITY
    }
}
impl Solid {
    /// Source/literal mesh intervals from the last successful native proof.
    /// No file was independently decoded by this query.
    pub fn rust_phillips_mesh_volume_pair(&self) -> Option<PhillipsMeshVolumePair> {
        let (body, frame) = self.helical_body()?;
        let HelicalForm::PhillipsScrew(body) = &body.form else {
            return None;
        };
        let pair = (*body.last_mesh_pair.borrow())?;
        if [frame.origin, frame.u, frame.v]
            .into_iter()
            .flatten()
            .map(f64::to_bits)
            .eq([pair.frame.origin, pair.frame.u, pair.frame.v]
                .into_iter()
                .flatten()
                .map(f64::to_bits))
        {
            Some(pair)
        } else {
            None
        }
    }
    pub fn rust_phillips_screw_geometry(&self) -> Option<(&PhillipsScrew, PlaneFrame)> {
        let (body, frame) = self.helical_body()?;
        let HelicalForm::PhillipsScrew(source) = &body.form else {
            return None;
        };
        Some((source.source(), frame))
    }
    /// Retain the exact supplied thread, head meridian and finite cavity in
    /// their original thread frame. Physical STEP remains unsupported.
    pub fn rust_phillips_screw_from_source(
        source: PhillipsScrew,
        feature: EntityId,
        frame: PlaneFrame,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        Self::make_phillips_source(source, feature, frame, policy, &|| false)
    }
    pub fn rust_phillips_screw_from_source_with_interrupt(
        source: PhillipsScrew,
        feature: EntityId,
        frame: PlaneFrame,
        policy: GeometryPolicy,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, GeometryError> {
        Self::make_phillips_source(source, feature, frame, policy, &|| {
            interrupt.is_some_and(|f| f())
        })
    }
    #[cfg(test)]
    fn rust_phillips_source_for_test(
        source: PhillipsScrew,
        feature: EntityId,
        frame: PlaneFrame,
        policy: GeometryPolicy,
        interrupt: &dyn Fn() -> bool,
    ) -> Result<Self, GeometryError> {
        Self::make_phillips_source(source, feature, frame, policy, interrupt)
    }
    fn make_phillips_source(
        source: PhillipsScrew,
        feature: EntityId,
        frame: PlaneFrame,
        policy: GeometryPolicy,
        interrupt: &dyn Fn() -> bool,
    ) -> Result<Self, GeometryError> {
        policy.validate()?;
        frame.checked(policy)?;
        let geometry = *source.thread();
        checked_thread(
            imported_request(geometry.spec()),
            geometry.spec().root_radius_mm,
            policy,
        )?;
        let body = Body::new(source, policy, interrupt)?;
        Self::publish_helical_geometry(geometry, feature, frame, policy, |_| {
            Ok(HelicalForm::PhillipsScrew(Box::new(body)))
        })
    }
}
impl HelicalBody {
    pub(in crate::geometry) fn phillips_volume(
        &self,
        frame: PlaneFrame,
    ) -> Option<Result<[f64; 2], GeometryError>> {
        let HelicalForm::PhillipsScrew(body) = &self.form else {
            return None;
        };
        Some(properties::volume_bounds(body.volume, frame))
    }
    pub(in crate::geometry) fn has_phillips_source(&self) -> bool {
        matches!(&self.form, HelicalForm::PhillipsScrew(_))
    }
    pub(in crate::geometry) fn phillips_faces(
        &self,
        frame: PlaneFrame,
        policy: GeometryPolicy,
    ) -> Option<Result<Vec<Face>, GeometryError>> {
        let HelicalForm::PhillipsScrew(body) = &self.form else {
            return None;
        };
        Some(body.faces(frame, policy))
    }
}

#[cfg(test)]
mod tests;
