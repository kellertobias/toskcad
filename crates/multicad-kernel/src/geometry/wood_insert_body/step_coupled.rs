//! Additive exact-unit insert path. Callers must earn the independent actual
//! material and whole-source correspondence before publishing any STEP bytes.
use super::step::{fail, read_error};
use super::*;
use multicad_rust_kernel::{
    step::{
        self,
        insert::{
            self,
            coupled::{CoupledPlan, topology::CoupledTopology},
        },
    },
    step_read::owned_insert as read,
};
const BYTES: usize = 64 * 1024 * 1024;
const WORK: usize = 16_000_000;
const SCRATCH: usize = 4096;
fn plan_error(e: insert::Error) -> GeometryError {
    fail(
        match e {
            insert::Error::Invalid => GeometryErrorKind::Input,
            insert::Error::Limit => GeometryErrorKind::Limit,
            insert::Error::Cancelled => GeometryErrorKind::Cancelled,
            insert::Error::Unresolved => GeometryErrorKind::Unsupported,
        },
        format!("Coupled insert source proof failed:{e:?}"),
    )
}
/// A representation alternative must preserve an earlier specific refusal.
/// Cancellation and a measured new mismatch after an unsupported old format
/// remain authoritative; a generic alternative refusal adds no information.
pub(super) fn preserve_failure(
    previous: GeometryError,
    alternative: GeometryError,
) -> GeometryError {
    match (previous.kind, alternative.kind) {
        (_, GeometryErrorKind::Cancelled)
        | (GeometryErrorKind::Unsupported, GeometryErrorKind::InvalidResult) => alternative,
        _ => previous,
    }
}
pub(super) struct ParsedProofScope {
    pub available: usize,
    pub fixed: usize,
    pub frame: PlaneFrame,
    pub policy: GeometryPolicy,
}
impl WoodInsertBody {
    pub(super) fn with_coupled_plan<T>(
        &self,
        policy: GeometryPolicy,
        interrupt: Option<&dyn Fn() -> bool>,
        f: impl FnOnce(&CoupledTopology<'_, '_>) -> Result<T, GeometryError>,
    ) -> Result<T, GeometryError> {
        policy.validate()?;
        let plan = CoupledPlan::prepare(
            self.boundary(),
            insert::CurveOptions {
                position_error_mm: policy.linear_tolerance_mm * 0.25,
                ..Default::default()
            },
            interrupt,
        )
        .map_err(plan_error)?;
        let topology =
            CoupledTopology::prepare(&plan, WORK, BYTES, interrupt).map_err(plan_error)?;
        f(&topology)
    }
    fn verify_coupled_candidate(
        &self,
        topology: &CoupledTopology<'_, '_>,
        text: &str,
        slack: usize,
        frame: PlaneFrame,
        policy: GeometryPolicy,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<WoodInsertStepVerification, GeometryError> {
        let fixed = SCRATCH + size_of::<WoodInsertStepVerification>() + slack;
        let available = BYTES
            .checked_sub(topology.retained_bytes().saturating_add(fixed))
            .filter(|&n| n > 0)
            .ok_or_else(|| {
                fail(
                    GeometryErrorKind::Limit,
                    "Coupled source leaves no selected-file arena",
                )
            })?;
        let options = read::Options {
            maximum_work: WORK,
            retained_budget_bytes: available,
        };
        let document =
            read::CompoundDocument::parse(text, read::DocumentOptions::default(), interrupt)
                .map_err(read_error)?;
        let [member] = document.represented_members() else {
            return Err(fail(
                GeometryErrorKind::InvalidResult,
                "Coupled candidate requires one actual member",
            ));
        };
        let link = member.mapping_link().ok_or_else(|| {
            fail(
                GeometryErrorKind::InvalidResult,
                "Coupled candidate requires its actual coefficient-unit mapping",
            )
        })?;
        let packet = document
            .member_packet(member.item_entity(), options, interrupt)
            .map_err(read_error)?;
        let actual = packet.parsed();
        let mapping = actual
            .physical_mapping(link.operator_entity(), options, interrupt)
            .map_err(read_error)?;
        self.verify_coupled_parsed(
            topology,
            actual,
            &mapping,
            ParsedProofScope {
                available,
                fixed,
                frame,
                policy,
            },
            interrupt,
        )
    }
    pub(super) fn verify_coupled_parsed(
        &self,
        topology: &CoupledTopology<'_, '_>,
        actual: &read::Parsed<'_>,
        mapping: &read::PhysicalMapping<'_, '_>,
        scope: ParsedProofScope,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<WoodInsertStepVerification, GeometryError> {
        let ParsedProofScope {
            available,
            fixed,
            frame,
            policy,
        } = scope;
        let options = read::Options {
            maximum_work: WORK,
            retained_budget_bytes: available,
        };
        let binding = actual
            .support_bindings(options, interrupt)
            .map_err(read_error)?;
        let material = binding
            .material_bounds(options, interrupt)
            .map_err(read_error)?;
        // This opaque proof also owns actual cap Jordan/nesting prerequisites.
        // A valid differently scaled body is not the expected source body.
        let correspondence = material
            .coupled_source_correspondence(
                topology,
                mapping,
                read::CorrespondenceOptions {
                    position_limit_mm: policy.linear_tolerance_mm,
                    retained_budget_bytes: BYTES - fixed,
                    ..Default::default()
                },
                interrupt,
            )
            .map_err(read_error)?;
        let scale = mapping.volume_scale_bounds().map_err(read_error)?;
        let raw_width = (B::point(1e-6) / B::point(scale[1])).lo;
        let raw_volume = material
            .occupied_volume(
                read::FluxOptions {
                    maximum_width_mm3: raw_width,
                    maximum_bytes: available,
                    ..Default::default()
                },
                interrupt,
            )
            .map_err(read_error)?
            .volume_mm3();
        let source_volume = self
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
        let bounds = material.bounds_mm();
        let proof = WoodInsertStepVerification {
            decoded_bounds_mm: [bounds.map(|b| b[0]), bounds.map(|b| b[1])],
            decoded_volume_bounds_mm3: raw_volume,
            source_volume_bounds_mm3: source_volume,
            maximum_position_error_mm: correspondence.position_error_mm(),
            maximum_row_derivative_error_mm: correspondence.row_derivative_error_mm(),
            maximum_ruling_derivative_error_mm: correspondence.ruling_derivative_error_mm(),
            verified_tolerance_mm: policy.linear_tolerance_mm,
        };
        super::step_coupled_pose::transport(
            &self.source,
            frame,
            self.boundary(),
            mapping,
            proof,
            interrupt,
        )
    }
    pub(super) fn append_coupled_step(
        &self,
        writer: &mut step::Writer,
        frame: PlaneFrame,
        policy: GeometryPolicy,
    ) -> Result<(), GeometryError> {
        self.material(frame, policy)?;
        self.with_coupled_plan(policy, None, |topology| {
            let text = topology
                .write_mapped_candidate(frame.origin, frame.u, frame.v, WORK, BYTES, None)
                .map_err(plan_error)?;
            self.verify_coupled_candidate(
                topology,
                &text,
                text.capacity() - text.len(),
                frame,
                policy,
                None,
            )?;
            drop(text);
            writer
                .append_coupled_insert_candidate(
                    topology,
                    frame.origin,
                    frame.u,
                    frame.v,
                    WORK,
                    BYTES,
                    None,
                )
                .map_err(plan_error)
        })
    }
    pub(super) fn verify_coupled_step(
        &self,
        text: &str,
        frame: PlaneFrame,
        policy: GeometryPolicy,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<WoodInsertStepVerification, GeometryError> {
        self.material(frame, policy)?;
        self.with_coupled_plan(policy, interrupt, |topology| {
            self.verify_coupled_candidate(topology, text, 0, frame, policy, interrupt)
        })
    }
}
