//! Explicit authored pipe transport; receiving tools retain separate admission.
use super::*;
use multicad_kernel::rust_kernel::axial_pipe_thread::authored::{AuthoredSection, AuthoredSpec};
pub(super) fn shaft(
    feature: EntityId,
    family: ThreadFamily,
    thread: &ResolvedThread,
    policy: GeometryPolicy,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<Solid, ModelError> {
    Solid::rust_authored_axial_pipe_thread_with_interrupt(
        authored_spec(feature, family, thread)?,
        feature,
        PlaneFrame {
            origin: thread.frame.origin,
            u: thread.frame.u,
            v: thread.frame.v,
        },
        policy,
        interrupt,
    )
    .map_err(|e| pipe_error(feature, e))
}

pub(super) fn cut_external_stock(
    feature: EntityId,
    thread: &ResolvedThread,
    stock: &Solid,
    policy: GeometryPolicy,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<Solid, ModelError> {
    stock
        .rust_cut_external_npt_stock(
            authored_spec(feature, ThreadFamily::Npt, thread)?,
            feature,
            PlaneFrame {
                origin: thread.frame.origin,
                u: thread.frame.u,
                v: thread.frame.v,
            },
            policy,
            interrupt,
        )
        .map_err(|e| pipe_error(feature, e))
}

fn pipe_error(feature: EntityId, cause: multicad_kernel::GeometryError) -> ModelError {
    if cause.kind == multicad_kernel::GeometryErrorKind::Cancelled {
        error(Some(feature), "Model evaluation superseded")
    } else {
        error(Some(feature), cause)
    }
}

fn authored_spec(
    feature: EntityId,
    family: ThreadFamily,
    thread: &ResolvedThread,
) -> Result<AuthoredSpec, ModelError> {
    let fail = |message| error(Some(feature), message);
    if thread.hand != ThreadHand::Right
        || thread.starts != 1
        || thread
            .end_treatments
            .as_ref()
            .is_some_and(|e| e.start_chamfer_radial_depth_mm != 0. || e.end_relief_length_mm != 0.)
    {
        return Err(fail(
            "Authored axial pipe threads require a right-hand single-start source without end treatments",
        ));
    }
    let (section, limit) = match family {
        ThreadFamily::Npt => {
            let npt = thread
                .npt
                .as_ref()
                .ok_or_else(|| fail("Canonical NPT source is missing"))?;
            (
                AuthoredSection::NptFlat {
                    depth_mm: npt.depth_mm,
                    flat_width_mm: npt.flat_width_mm,
                },
                npt.stock_crest_radius_limit_mm,
            )
        }
        ThreadFamily::Bspt if thread.pipe.is_some() => (
            AuthoredSection::BsptRounded {
                depth_mm: thread.profile.depth_mm,
            },
            None,
        ),
        _ => return Err(fail("Canonical BSPT source is missing")),
    };
    Ok(AuthoredSpec {
        major_radius_mm: thread.major_radius_mm,
        pitch_mm: thread.pitch_mm,
        length_mm: thread.length_mm,
        radial_taper: thread.profile.radial_taper,
        section,
        crest_radius_limit_mm: limit,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use multicad_core::{document::Document, model::ThreadParameters, parameters};
    #[test]
    fn pipe_transport_refuses_missing_or_unsupported_physical_source_without_mutation() {
        let document = Document::new("pipe guards").unwrap();
        let evaluated = parameters::evaluate(&document).unwrap();
        let resolved = ThreadParameters::default()
            .resolve(&evaluated, ThreadMode::NewExternal)
            .unwrap();
        let feature = EntityId::new();
        let run = |r: &ResolvedThread, family| {
            shaft(feature, family, r, GeometryPolicy::default(), None)
                .err()
                .unwrap()
        };
        assert!(
            run(&resolved, ThreadFamily::Npt)
                .to_string()
                .contains("NPT source is missing")
        );
        assert!(
            run(&resolved, ThreadFamily::Bspt)
                .to_string()
                .contains("BSPT source is missing")
        );
        let mut unsupported = resolved.clone();
        unsupported.hand = ThreadHand::Left;
        assert!(
            run(&unsupported, ThreadFamily::Bspt)
                .to_string()
                .contains("right-hand")
        );
        unsupported = resolved.clone();
        unsupported.starts = 2;
        assert!(
            run(&unsupported, ThreadFamily::Bspt)
                .to_string()
                .contains("single-start")
        );
        assert_eq!(document, evaluated.document);
    }
    #[test]
    fn resolved_pipe_cancel_and_end_geometry_refuse_without_repairing_source() {
        use multicad_core::model::{
            self, ResolvedThreadEndTreatments, ThreadEndTreatment, ThreadEndTreatments,
            ThreadStartTreatment,
        };
        let document = Document::new("canonical pipe guards").unwrap();
        let evaluated = parameters::evaluate(&document).unwrap();
        let intent = model::PipeThreadIntent {
            form: model::PipeThreadForm::R,
            size: model::PipeThreadSize::OneEighth,
            gauge: model::PipeGaugeConvention::Iso7Nominal1994,
        };
        let parameters = ThreadParameters {
            family: ThreadFamily::Bspt,
            pipe: Some(intent),
            fit: model::ThreadFit::BsptNominal,
            major_diameter: model::length(intent.base_major_mm()),
            pitch: model::length(intent.preset().pitch_mm),
            length: model::length(intent.default_length_mm()),
            ..Default::default()
        };
        let mut resolved = parameters
            .resolve(&evaluated, ThreadMode::NewExternal)
            .unwrap();
        let feature = EntityId::new();
        let cancel = shaft(
            feature,
            ThreadFamily::Bspt,
            &resolved,
            GeometryPolicy::default(),
            Some(&|| true),
        )
        .err()
        .unwrap();
        assert!(cancel.to_string().contains("superseded"));
        let npt = model::NptThreadIntent {
            form: model::NptThreadForm::External,
            size: model::NptThreadSize::OneEighth,
            gauge: model::NptGaugeConvention::Nbs1959PlainFace,
            truncation: model::NptTruncationTarget::MidpointConeParallel,
        };
        let npt_parameters = ThreadParameters {
            family: ThreadFamily::Npt,
            npt: Some(npt),
            fit: model::ThreadFit::NptHistoricalMidpoint,
            major_diameter: model::length(npt.base_major_mm()),
            pitch: model::length(npt.preset().pitch_mm),
            length: model::length(npt.default_length_mm()),
            ..Default::default()
        };
        let npt_source = npt_parameters
            .resolve(&evaluated, ThreadMode::NewExternal)
            .unwrap();
        let clipped = shaft(
            feature,
            ThreadFamily::Npt,
            &npt_source,
            GeometryPolicy::default(),
            None,
        )
        .unwrap();
        assert!(clipped.rust_clipped_npt_source().is_some());
        // Local STEP admission independently proves literal material and its
        // whole-domain correspondence to this unchanged source.
        assert!(clipped.supports_step_export());
        resolved.end_treatments = Some(ResolvedThreadEndTreatments {
            intent: ThreadEndTreatments {
                start: ThreadStartTreatment::FlatCut,
                end: ThreadEndTreatment::FlatCut,
            },
            start_chamfer_radial_depth_mm: 0.,
            end_relief_length_mm: 1.,
            end_relief_radial_depth_mm: 0.1,
            full_profile_start_mm: 0.,
            full_profile_end_mm: 5.5,
            full_profile_length_mm: 5.5,
            maximum_cutter_radius_mm: 5.,
            scope: "test actual unsupported relief",
        });
        assert!(
            shaft(
                feature,
                ThreadFamily::Bspt,
                &resolved,
                GeometryPolicy::default(),
                None
            )
            .err()
            .unwrap()
            .to_string()
            .contains("without end treatments")
        );
    }
}
