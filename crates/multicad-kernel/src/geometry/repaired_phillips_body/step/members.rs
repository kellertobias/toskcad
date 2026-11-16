//! Explicit actual represented-member selection. The existing local view
//! preserves original dependency IDs/components; no geometry recipe is rebuilt.
use super::*;
impl RepairedPhillipsBody {
    pub(in crate::geometry) fn verify_step_member(
        &self,
        actual: &read::CompoundDocument<'_>,
        member: read::Member,
        frame: PlaneFrame,
        policy: GeometryPolicy,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<DecodedPeriodicPhillipsStepVerification, GeometryError> {
        checkpoint(interrupt)?;
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
            return Err(fail(
                GeometryErrorKind::Unsupported,
                "Mapped periodic Phillips member needs an independent occurrence proof",
            ));
        }
        let selected = actual
            .local_member_document(member.item_entity(), Default::default(), interrupt)
            .map_err(document_error)?;
        let mut proof = self.verify_singleton_step(
            selected.as_str(),
            frame,
            policy,
            policy.linear_tolerance_mm,
            interrupt,
        )?;
        proof.represented_member_count = actual.represented_members().len();
        proof.member_index = index;
        proof.selection_work = selected.work_used();
        proof.selection_query_peak_bytes = selected.query_budget_bytes();
        // Full input/AST is an already completed caller scope; its work is
        // charged by the public parse wrapper or by the enclosing batch once.
        // Both arenas are resident here, so report their retained upper sum.
        proof.document_retained_bytes += actual.retained_budget_bytes();
        checkpoint(interrupt)?;
        Ok(proof)
    }
}
impl Solid {
    /// Verify exactly the requested represented-member index in these actual
    /// bytes. Global count/order is reported, never inferred from a candidate.
    pub fn verify_periodic_phillips_step_member(
        &self,
        text: &str,
        member_index: usize,
        tolerance: f64,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Option<Result<super::super::DecodedPeriodicPhillipsStepVerification, GeometryError>> {
        let (body, frame) = self.repaired_phillips_body()?;
        body.source.periodic()?;
        Some((|| {
            checkpoint(interrupt)?;
            let options = options(self.policy, tolerance)?;
            let actual = read::CompoundDocument::parse(text, Default::default(), interrupt)
                .map_err(document_error)?;
            let member = *actual
                .represented_members()
                .get(member_index)
                .ok_or_else(|| {
                    fail(
                        GeometryErrorKind::Input,
                        "Represented member index is out of range",
                    )
                })?;
            let policy = GeometryPolicy {
                linear_tolerance_mm: options.correspondence.maximum_distance_mm,
                ..self.policy
            };
            let mut proof = body.verify_step_member(&actual, member, frame, policy, interrupt)?;
            proof.document_work += actual.work_used();
            Ok(proof)
        })())
    }
    /// Candidate material is independently decoded and geometrically matched
    /// before any append. The actual rebased payload is then parsed and bound;
    /// the enclosing writer transaction restores its previous state on failure.
    pub fn append_periodic_phillips_step(
        &self,
        writer: &mut multicad_rust_kernel::step::Writer,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Option<Result<exchange::AppendReceipt, GeometryError>> {
        let (body, frame) = self.repaired_phillips_body()?;
        let source = body.source.periodic()?;
        Some((|| {
            checkpoint(interrupt)?;
            let candidate = exchange::export(
                source,
                occurrence(frame),
                options(self.policy, self.policy.linear_tolerance_mm)?,
                interrupt,
            )
            .map_err(exchange_error)?;
            writer
                .append_qualified_periodic_phillips(&candidate, Default::default(), interrupt)
                .map_err(exchange_error)
        })())
    }
}
