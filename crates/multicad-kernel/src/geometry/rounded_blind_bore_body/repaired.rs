//! Closed canonical interpretation and immutable raw-member provenance.
use super::*;
use multicad_rust_kernel::step_read::{ImportProfile, ImportSolid, RepairedRoundedBlindCandidate};
impl Solid {
    pub(in crate::geometry) fn import_repaired_rounded_blind(
        candidate: RepairedRoundedBlindCandidate,
        frame: PlaneFrame,
        feature: EntityId,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        if candidate.frame() != [frame.origin, frame.u, frame.v] {
            return Err(failure(
                GeometryErrorKind::Unsupported,
                "Repaired blind operator must retain its earned source frame",
            ));
        }
        let source = *candidate.canonical_bore();
        let identity = super::super::step_import::feature_origin(
            feature,
            "repaired-rounded-blind.v1",
            super::super::rounded_thread_body::identity_values(source.cutter().spec(), frame)
                .chain([
                    source.stock_radius_mm().as_f64(),
                    source.stock_height_mm().as_f64(),
                ]),
        );
        let mut solid = Self::rust_rounded_blind_bore(
            source.cutter().spec(),
            source.stock_radius_mm(),
            source.stock_height_mm(),
            feature,
            frame,
            policy,
        )?;
        let SolidBody::RoundedBlindBore(body, _) = &mut solid.body else {
            unreachable!("blind source constructor")
        };
        let body = Rc::get_mut(body).expect("unshared repaired blind");
        body.import_source = Some(candidate);
        body.import_identity = Some(identity);
        Ok(solid)
    }
    /// Immutable original member records and repair contract. Physical queries
    /// use the closed canonical body, not unsealed raw STEP occupied material.
    /// The candidate retains its own original frame; moving this Solid does
    /// not rewrite that input frame or claim a new raw-file correspondence.
    pub fn rust_repaired_rounded_blind_source(&self) -> Option<&RepairedRoundedBlindCandidate> {
        self.rounded_blind_bore_body()?.0.import_source.as_ref()
    }
    pub(in crate::geometry) fn matches_repaired_rounded_blind(
        &self,
        imported: &ImportSolid,
        tolerance: f64,
    ) -> bool {
        let Some((body, frame)) = self.rounded_blind_bore_body() else {
            return false;
        };
        let ImportProfile::RepairedRoundedBlind { candidate } = &imported.profile else {
            return false;
        };
        candidate.frame() == [imported.origin_mm, imported.u, imported.v]
            && multicad_rust_kernel::step_read::repaired_rounded_blind_import_matches(
                candidate,
                &body.geometry,
                [frame.origin, frame.u, frame.v],
                tolerance,
            )
    }
}
