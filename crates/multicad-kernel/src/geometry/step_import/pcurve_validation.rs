//! Preserve actual trim disagreement evidence at the public import boundary.
use super::*;

pub(in crate::geometry) fn validate<'i>(
    input: &'i str,
    policy: GeometryPolicy,
) -> Result<step_read::brep::PcurveValidatedInput<'i>, GeometryError> {
    step_read::brep::validate_pcurve_input(input, policy.linear_tolerance_mm, None)
        .map_err(|cause| validation_error(cause, policy.linear_tolerance_mm))
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
