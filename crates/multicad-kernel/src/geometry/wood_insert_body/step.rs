//! Normal retained-insert export. Actual file material/volume and measured
//! source correspondence are independent prerequisites, never header claims.
use super::*;
use multicad_rust_kernel::{
    step::{self, insert},
    step_read::owned_insert as read,
};
const BYTES: usize = 64 * 1024 * 1024;
const WORK: usize = 16_000_000;
const SCRATCH: usize = 4096;
#[derive(Clone, Debug, PartialEq)]
pub struct WoodInsertStepVerification {
    pub decoded_bounds_mm: [[f64; 3]; 2],
    pub decoded_volume_bounds_mm3: [f64; 2],
    pub source_volume_bounds_mm3: [f64; 2],
    pub maximum_position_error_mm: f64,
    pub maximum_row_derivative_error_mm: f64,
    pub maximum_ruling_derivative_error_mm: f64,
    pub verified_tolerance_mm: f64,
}
pub(super) fn fail(kind: GeometryErrorKind, message: impl Into<String>) -> GeometryError {
    error("Rust insert STEP", kind, message)
}
fn plan_error(e: insert::Error) -> GeometryError {
    fail(
        match e {
            insert::Error::Invalid => GeometryErrorKind::Input,
            insert::Error::Limit => GeometryErrorKind::Limit,
            insert::Error::Cancelled => GeometryErrorKind::Cancelled,
            insert::Error::Unresolved => GeometryErrorKind::Unsupported,
        },
        format!("Actual insert encoder proof failed:{e:?}"),
    )
}
pub(super) fn read_error(e: read::Error) -> GeometryError {
    fail(
        match e {
            read::Error::Invalid => GeometryErrorKind::InvalidResult,
            read::Error::Limit => GeometryErrorKind::Limit,
            read::Error::Cancelled => GeometryErrorKind::Cancelled,
            read::Error::Unsupported => GeometryErrorKind::Unsupported,
        },
        format!("Actual insert file proof failed:{e:?}"),
    )
}
fn step_error(e: GeometryError) -> step::StepError {
    match e.kind {
        GeometryErrorKind::Limit => step::StepError::Limit,
        GeometryErrorKind::Unsupported | GeometryErrorKind::Cancelled => {
            step::StepError::Unsupported
        }
        _ => step::StepError::InvalidGeometry,
    }
}
impl WoodInsertBody {
    pub(in crate::geometry) fn supports_step(
        &self,
        frame: PlaneFrame,
        policy: GeometryPolicy,
    ) -> bool {
        self.material(frame, policy).is_ok()
    }
    pub(super) fn with_plan<T>(
        &self,
        policy: GeometryPolicy,
        interrupt: Option<&dyn Fn() -> bool>,
        f: impl FnOnce(&insert::TopologyPlan<'_, '_, '_>) -> Result<T, GeometryError>,
    ) -> Result<T, GeometryError> {
        policy.validate()?;
        let domain = insert::DomainPlan::prepare(
            self.boundary(),
            insert::DomainOptions::default(),
            interrupt,
        )
        .map_err(plan_error)?;
        let surface = insert::SurfacePlan::prepare(
            &domain,
            insert::CurveOptions {
                position_error_mm: policy.linear_tolerance_mm * 0.25,
                ..Default::default()
            },
            interrupt,
        )
        .map_err(plan_error)?;
        let topology =
            insert::TopologyPlan::prepare(&surface, WORK, BYTES, interrupt).map_err(plan_error)?;
        f(&topology)
    }
    fn verify_candidate(
        &self,
        topology: &insert::TopologyPlan<'_, '_, '_>,
        text: &str,
        capacity_slack: usize,
        frame: PlaneFrame,
        policy: GeometryPolicy,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<WoodInsertStepVerification, GeometryError> {
        let fixed = SCRATCH + size_of::<WoodInsertStepVerification>() + capacity_slack;
        let read_bytes = BYTES
            .checked_sub(topology.retained_bytes().saturating_add(fixed))
            .filter(|&n| n > 0)
            .ok_or_else(|| {
                fail(
                    GeometryErrorKind::Limit,
                    "Source plan leaves no room for the actual file arena",
                )
            })?;
        let options = read::Options {
            maximum_work: WORK,
            retained_budget_bytes: read_bytes,
        };
        // Charge the original document once in its explicit global scope,
        // exactly as mixed verification does. The selected proof below keeps
        // its unchanged 16M/64MiB limits, including the live source topology.
        let file = read::CompoundDocument::parse(text, read::DocumentOptions::default(), interrupt)
            .map_err(read_error)?;
        let [member] = file.represented_members() else {
            return Err(fail(
                GeometryErrorKind::InvalidResult,
                "Insert candidate must have one actual occurrence",
            ));
        };
        let packet = file
            .member_packet(member.item_entity(), options, interrupt)
            .map_err(read_error)?;
        let actual = packet.parsed();
        let proof = self.verify_parsed(topology, actual, read_bytes, fixed, policy, interrupt)?;
        if let Some(link) = member.mapping_link() {
            let mapping = actual
                .physical_mapping(link.operator_entity(), options, interrupt)
                .map_err(read_error)?;
            super::step_pose::transport(&self.source, frame, topology, &mapping, proof, interrupt)
        } else if frame == PlaneFrame::default() {
            Ok(proof)
        } else {
            Err(fail(
                GeometryErrorKind::InvalidResult,
                "Placed source requires an explicit STEP map",
            ))
        }
    }
    pub(super) fn verify_parsed(
        &self,
        topology: &insert::TopologyPlan<'_, '_, '_>,
        actual: &read::Parsed<'_>,
        read_bytes: usize,
        fixed: usize,
        policy: GeometryPolicy,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<WoodInsertStepVerification, GeometryError> {
        let options = read::Options {
            maximum_work: WORK,
            retained_budget_bytes: read_bytes,
        };
        let binding = actual
            .support_bindings(options, interrupt)
            .map_err(read_error)?;
        let correspondence = binding
            .source_correspondence(
                topology,
                read::CorrespondenceOptions {
                    position_limit_mm: policy.linear_tolerance_mm,
                    retained_budget_bytes: BYTES - fixed,
                    ..Default::default()
                },
                interrupt,
            )
            .map_err(read_error)?;
        let errors = [
            correspondence.position_error_mm(),
            correspondence.row_derivative_error_mm(),
            correspondence.ruling_derivative_error_mm(),
        ];
        let material = binding
            .material_bounds(options, interrupt)
            .map_err(read_error)?;
        let decoded_bounds = material.bounds_mm();
        let volume = material
            .occupied_volume(
                read::FluxOptions {
                    maximum_bytes: read_bytes,
                    ..Default::default()
                },
                interrupt,
            )
            .map_err(read_error)?
            .volume_mm3();
        // This scalar source query is a separate authority from decoded mass.
        let source = self
            .source
            .volume(
                source::QueryOptions {
                    volume_error_mm3: 1e-6,
                    maximum_cells: 65536,
                },
                interrupt,
            )
            .map_err(failure)?
            .material_mm3;
        Ok(WoodInsertStepVerification {
            decoded_bounds_mm: [decoded_bounds.map(|b| b[0]), decoded_bounds.map(|b| b[1])],
            decoded_volume_bounds_mm3: volume,
            source_volume_bounds_mm3: source,
            maximum_position_error_mm: errors[0],
            maximum_row_derivative_error_mm: errors[1],
            maximum_ruling_derivative_error_mm: errors[2],
            verified_tolerance_mm: policy.linear_tolerance_mm,
        })
    }
    pub(in crate::geometry) fn append_step(
        &self,
        writer: &mut step::Writer,
        frame: PlaneFrame,
        policy: GeometryPolicy,
    ) -> Result<(), step::StepError> {
        if !self.supports_step(frame, policy) {
            return Err(step::StepError::Unsupported);
        }
        self.with_plan(policy, None, |topology| {
            let local = frame == PlaneFrame::default();
            let text = if local {
                insert::write_local_candidate(topology, WORK, BYTES, None)
            } else {
                insert::write_mapped_candidate(
                    topology,
                    frame.origin,
                    frame.u,
                    frame.v,
                    WORK,
                    BYTES,
                    None,
                )
            }
            .map_err(plan_error)?;
            self.verify_candidate(
                topology,
                &text,
                text.capacity() - text.len(),
                frame,
                policy,
                None,
            )?;
            drop(text);
            // The same immutable topology and encoder emit identical numeric
            // leaves under entity renumbering; existing members retain their
            // own identities through the writer's typed transaction.
            if local {
                writer.append_owned_insert_local(topology, WORK, BYTES, None)
            } else {
                writer.append_owned_insert_mapped_candidate(
                    topology,
                    frame.origin,
                    frame.u,
                    frame.v,
                    WORK,
                    BYTES,
                    None,
                )
            }
            .map_err(plan_error)
        })
        .map_err(step_error)
    }
}
impl Solid {
    /// Independent actual file mass plus a distinct measured correspondence
    /// receipt. Exactly one actual local or explicitly mapped occurrence.
    pub fn verify_wood_insert_step_source(
        &self,
        text: &str,
        tolerance: f64,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Option<Result<WoodInsertStepVerification, GeometryError>> {
        let (body, frame) = self.wood_insert_body()?;
        if !tolerance.is_finite() || tolerance <= 0. {
            return Some(Err(fail(
                GeometryErrorKind::Input,
                "STEP source verification needs positive finite tolerance",
            )));
        }
        let policy = GeometryPolicy {
            linear_tolerance_mm: self.policy.linear_tolerance_mm.min(tolerance),
            ..self.policy
        };
        Some(body.with_plan(policy, interrupt, |topology| {
            body.verify_candidate(topology, text, 0, frame, policy, interrupt)
        }))
    }
}
