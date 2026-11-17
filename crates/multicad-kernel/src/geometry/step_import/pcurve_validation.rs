//! Preserve actual trim disagreement evidence at the public import boundary.
use super::*;

pub(in crate::geometry) fn validate<'i>(
    input: &'i str,
    policy: GeometryPolicy,
) -> Result<step_read::brep::PcurveValidatedInput<'i>, GeometryError> {
    step_read::brep::validate_pcurve_input(input, policy.linear_tolerance_mm, None).map_err(
        |cause| {
            if !cause.cancelled
                && cause.kind == ReadError::Unsupported
                && matches!(
                    step_read::brep::retained_geometry_refusal(
                        input,
                        policy.linear_tolerance_mm,
                        None
                    ),
                    Some(multicad_rust_kernel::brep::Error::Invalid)
                )
            {
                // A missing trim certificate cannot hide a separately proved
                // invalid solid. This diagnostic path never returns a body.
                return crate::geometry::imported_brep_body::failure(
                    multicad_rust_kernel::brep::Error::Invalid,
                );
            }
            validation_error(cause, policy.linear_tolerance_mm)
        },
    )
}

fn validation_error(
    cause: step_read::brep::PcurveValidationError,
    tolerance: f64,
) -> GeometryError {
    if cause.cancelled {
        return error(
            IMPORT,
            GeometryErrorKind::Cancelled,
            "STEP trim verification cancelled",
        );
    }
    let kind = match cause.kind {
        ReadError::Invalid => GeometryErrorKind::Input,
        ReadError::Unsupported => GeometryErrorKind::Unsupported,
        ReadError::Limit => GeometryErrorKind::Limit,
    };
    let message = if cause.disagreements.is_empty() {
        match cause.kind {
            ReadError::Invalid => "Invalid STEP curve/surface trim data".to_owned(),
            ReadError::Unsupported => {
                "STEP curve/surface trim agreement is outside the certified Rust families"
                    .to_owned()
            }
            ReadError::Limit => {
                "STEP curve/surface trim verification exceeds the Rust budget".to_owned()
            }
        }
    } else {
        let details = cause
            .disagreements
            .iter()
            .map(|d| {
                format!(
                    "{} #{}; edge #{}; pcurve #{}; deviation_mm>={:.17E}; tolerance_mm={:.17E}",
                    d.surface_kind.to_ascii_lowercase(),
                    d.surface_entity,
                    d.edge_entity,
                    d.pcurve_entity,
                    d.deviation_mm,
                    tolerance,
                )
            })
            .collect::<Vec<_>>()
            .join(" | ");
        format!("STEP trim disagreement: {details}")
    };
    error(IMPORT, kind, message)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unqualified_trim_does_not_admit_a_body_or_mask_independent_invalid_geometry() {
        let original = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../multicad-app/tests/fixtures/exchange/corpus/occt-placed-box.step"
        ));
        let unqualified = original.replace(
            "#38 = LINE('',#39,#40);",
            "#38 = OFFSET_CURVE_2D('',#900001,0.,.F.);\n#900001 = LINE('',#39,#40);",
        );
        assert_ne!(original, unqualified);
        let refusal = Solid::import_rust_step(&unqualified, EntityId::new(), Default::default())
            .err()
            .expect("uncertified trims cannot import");
        assert_eq!(refusal.kind, GeometryErrorKind::Unsupported);
        let invalid = unqualified.replace(
            "#17 = ADVANCED_FACE('',(#18),#32,.F.);",
            "#17 = ADVANCED_FACE('',(#18),#32,.T.);",
        );
        assert_ne!(invalid, unqualified);
        let refusal = Solid::import_rust_step(&invalid, EntityId::new(), Default::default())
            .err()
            .expect("an inconsistent shell cannot import");
        assert_eq!(refusal.kind, GeometryErrorKind::InvalidResult);
        assert!(refusal.message.contains("not a valid closed solid"));
    }
}
