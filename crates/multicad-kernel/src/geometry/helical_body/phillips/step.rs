//! Physical exchange of an authored Phillips source through an independently
//! decoded, distinctly repaired degree-11 material matched back to the SAME
//! author in its occurrence frame. Written in world coordinates; a mapped
//! member is refused. Raw-file exact occupied mass is not asserted.
use super::*;
use multicad_rust_kernel::{step::degree11_phillips as exchange, step_read::owned_insert as read};
const STEP_OPERATION: &str = "authored Phillips STEP";
#[derive(Clone, Debug, PartialEq)]
pub struct DecodedAuthoredPhillipsStepVerification {
    pub represented_member_count: usize,
    pub member_index: usize,
    pub decoded_bounds_mm: [[f64; 3]; 2],
    /// Independently decoded repaired material, not raw occupied file mass.
    pub decoded_volume_bounds_mm3: [f64; 2],
    pub source_volume_bounds_mm3: [f64; 2],
    pub maximum_world_position_error_mm: f64,
    /// Whole-source paired normalized chart derivatives, not unit normals.
    pub maximum_world_chart_derivative_error_mm: f64,
    /// Author minus decoded repaired material.
    pub independent_volume_difference_mm3: [f64; 2],
    pub verified_tolerance_mm: f64,
    pub reader_reserved_work: usize,
    pub decoded_material_proof_work: usize,
    pub correspondence_work: usize,
    pub correspondence_cavity_cells: usize,
    pub document_work: usize,
    pub selection_work: usize,
}
fn fail(kind: GeometryErrorKind, message: impl Into<String>) -> GeometryError {
    error(STEP_OPERATION, kind, message)
}
fn cancelled(interrupt: Option<&dyn Fn() -> bool>) -> Result<(), GeometryError> {
    if interrupt.is_some_and(|f| f()) {
        Err(fail(GeometryErrorKind::Cancelled, "Exchange superseded"))
    } else {
        Ok(())
    }
}
fn exchange_error(e: exchange::Error) -> GeometryError {
    use multicad_rust_kernel::step_read::PhillipsReadError as R;
    let kind = match e {
        exchange::Error::Invalid => GeometryErrorKind::Input,
        exchange::Error::Limit => GeometryErrorKind::Limit,
        exchange::Error::Cancelled => GeometryErrorKind::Cancelled,
        exchange::Error::Unsupported => GeometryErrorKind::Unsupported,
        exchange::Error::Reader(e) => match e.cause {
            R::Invalid => GeometryErrorKind::InvalidResult,
            R::Limit => GeometryErrorKind::Limit,
            R::Cancelled => GeometryErrorKind::Cancelled,
            R::Unsupported => GeometryErrorKind::Unsupported,
        },
    };
    fail(
        kind,
        format!("Independent decoded authored Phillips qualification failed:{e:?}"),
    )
}
fn document_error(e: read::Error) -> GeometryError {
    fail(
        match e {
            read::Error::Invalid => GeometryErrorKind::InvalidResult,
            read::Error::Limit => GeometryErrorKind::Limit,
            read::Error::Cancelled => GeometryErrorKind::Cancelled,
            read::Error::Unsupported => GeometryErrorKind::Unsupported,
        },
        format!("Actual represented STEP document failed:{e:?}"),
    )
}
fn options(
    policy: GeometryPolicy,
    tolerance: f64,
) -> Result<(exchange::AuthoredOptions, f64), GeometryError> {
    policy.validate()?;
    if !tolerance.is_finite() || tolerance <= 0. {
        return Err(fail(
            GeometryErrorKind::Input,
            "Positive finite exchange tolerance required",
        ));
    }
    let mut options = exchange::AuthoredOptions::default();
    let tol = tolerance.min(policy.linear_tolerance_mm).min(1e-7);
    options.matching.maximum_distance_mm = tol;
    options.matching.maximum_derivative_mm = tol;
    Ok((options, tol))
}
fn frame3(frame: PlaneFrame) -> [[f64; 3]; 3] {
    [frame.origin, frame.u, frame.v]
}
/// Fresh proof of one ordinary singleton text against the author.
fn verify_singleton(
    body: &Body,
    frame: PlaneFrame,
    policy: GeometryPolicy,
    text: &str,
    tolerance: f64,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<DecodedAuthoredPhillipsStepVerification, GeometryError> {
    cancelled(interrupt)?;
    let (options, tol) = options(policy, tolerance)?;
    let document = read::CompoundDocument::parse(text, Default::default(), interrupt)
        .map_err(document_error)?;
    if document.represented_members().len() != 1
        || document.represented_members()[0].mapping_link().is_some()
    {
        return Err(fail(
            GeometryErrorKind::InvalidResult,
            "Exactly one ordinary represented member required",
        ));
    }
    let qualified =
        exchange::qualify_authored_file(&body.source, frame3(frame), text, options, interrupt)
            .map_err(exchange_error)?;
    let proof = qualified.verification();
    let bounds =
        crate::geometry::repaired_phillips_body::degree11_world_bounds(qualified.material())?;
    cancelled(interrupt)?;
    Ok(DecodedAuthoredPhillipsStepVerification {
        represented_member_count: 1,
        member_index: 0,
        decoded_bounds_mm: [bounds.min, bounds.max],
        decoded_volume_bounds_mm3: proof.decoded_repaired_volume_bounds_mm3,
        source_volume_bounds_mm3: properties::volume_bounds(body.volume, frame)?,
        maximum_world_position_error_mm: proof.maximum_world_distance_mm,
        maximum_world_chart_derivative_error_mm: proof.maximum_world_chart_derivative_mm,
        independent_volume_difference_mm3: proof.independent_volume_difference_mm3,
        verified_tolerance_mm: tol,
        reader_reserved_work: proof.reader_reserved_work,
        decoded_material_proof_work: proof.actual_decoded_material_proof_work,
        correspondence_work: proof.fresh_correspondence_work,
        correspondence_cavity_cells: proof.correspondence_cavity_cells,
        document_work: document.work_used(),
        selection_work: 0,
    })
}
impl Body {
    pub(in crate::geometry) fn supports_step(&self) -> bool {
        exchange::supports_author(&self.source)
    }
    pub(in crate::geometry) fn append_step(
        &self,
        writer: &mut multicad_rust_kernel::step::Writer,
        frame: PlaneFrame,
        policy: GeometryPolicy,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<exchange::AuthoredAppendReceipt, GeometryError> {
        cancelled(interrupt)?;
        if !self.supports_step() {
            return Err(fail(
                GeometryErrorKind::Unsupported,
                "Authored Phillips source is outside the exchange family",
            ));
        }
        let (options, tol) = options(policy, policy.linear_tolerance_mm)?;
        let candidate =
            exchange::export_authored(&self.source, frame3(frame), tol, options, interrupt)
                .map_err(exchange_error)?;
        writer
            .append_qualified_authored_degree11_phillips(&candidate, Default::default(), interrupt)
            .map_err(exchange_error)
    }
    /// Verify exactly this represented member of the actual parsed document.
    pub(in crate::geometry) fn verify_step_member(
        &self,
        actual: &read::CompoundDocument<'_>,
        member: read::Member,
        frame: PlaneFrame,
        policy: GeometryPolicy,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<DecodedAuthoredPhillipsStepVerification, GeometryError> {
        cancelled(interrupt)?;
        let index = actual
            .represented_members()
            .iter()
            .position(|m| m.item_entity() == member.item_entity())
            .ok_or_else(|| {
                fail(
                    GeometryErrorKind::Input,
                    "Selected member is absent from actual document",
                )
            })?;
        let member = actual.represented_members()[index];
        if member.mapping_link().is_some() {
            // The authored proposal is written in world coordinates.
            return Err(fail(
                GeometryErrorKind::InvalidResult,
                "Authored Phillips member must not carry an occurrence mapping",
            ));
        }
        let selected = actual
            .local_member_document(member.item_entity(), Default::default(), interrupt)
            .map_err(document_error)?;
        let mut proof = verify_singleton(
            self,
            frame,
            policy,
            selected.as_str(),
            policy.linear_tolerance_mm,
            interrupt,
        )?;
        proof.represented_member_count = actual.represented_members().len();
        proof.member_index = index;
        proof.selection_work = selected.work_used();
        cancelled(interrupt)?;
        Ok(proof)
    }
}
impl Solid {
    /// The authored Phillips source when its exchange family admits STEP.
    pub fn rust_authored_phillips_step_source(&self) -> Option<(&PhillipsScrew, PlaneFrame)> {
        let (body, frame) = self.helical_body()?;
        let HelicalForm::PhillipsScrew(body) = &body.form else {
            return None;
        };
        body.supports_step().then_some((&body.source, frame))
    }
    pub(in crate::geometry) fn authored_phillips_body(&self) -> Option<(&Body, PlaneFrame)> {
        let (body, frame) = self.helical_body()?;
        let HelicalForm::PhillipsScrew(body) = &body.form else {
            return None;
        };
        body.supports_step().then_some((body.as_ref(), frame))
    }
    /// Fresh proof of a whole singleton file against the current author.
    pub fn verify_authored_phillips_step_source(
        &self,
        text: &str,
        tolerance: f64,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Option<Result<DecodedAuthoredPhillipsStepVerification, GeometryError>> {
        let (body, frame) = self.authored_phillips_body()?;
        Some(verify_singleton(
            body,
            frame,
            self.policy,
            text,
            tolerance,
            interrupt,
        ))
    }
    /// Qualify the authored proposal, then append it by exact replay.
    pub fn append_authored_phillips_step(
        &self,
        writer: &mut multicad_rust_kernel::step::Writer,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Option<Result<exchange::AuthoredAppendReceipt, GeometryError>> {
        let (body, frame) = self.authored_phillips_body()?;
        Some(body.append_step(writer, frame, self.policy, interrupt))
    }
}
