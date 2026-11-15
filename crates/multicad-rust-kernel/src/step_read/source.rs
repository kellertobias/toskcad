//! Flat exchange members retain their actual kind of qualification evidence.
use super::*;
mod literal_npt;
#[cfg(test)]
mod tests;

/// Finite radial ownership still needs cold surface, ordering and occupied
/// material qualification. It deliberately carries no analytic volume receipt.
pub enum SourceImport {
    Analytic(Box<(Readback, ImportSolid)>),
    Radial(Box<crate::radial_tensor::BoundaryOwnership>),
    LiteralNpt(Box<super::npt_stock::OwnedMaterial>),
}

/// Decode bounded flat mixed representations without inferring a construction
/// recipe from retained tensor controls. Global topology ownership is shared
/// across every member, including analytic and retained members together.
pub fn read_source_compound(input: &str, tolerance: f64) -> Result<Vec<SourceImport>> {
    match ordinary(input, tolerance) {
        Err(ReadError::Unsupported) => literal_npt::read(input, tolerance),
        other => other,
    }
}
fn ordinary(input: &str, tolerance: f64) -> Result<Vec<SourceImport>> {
    compound::read_members_by(input, tolerance, true, false, ordinary_member)
}
pub(super) fn ordinary_member(
    shell: &mut Shell<'_, '_>,
    id: u32,
    may_be_warped: bool,
) -> Result<SourceImport> {
    let faces = if shell.simple(id)?.name == "MANIFOLD_SOLID_BREP" {
        shell.member_faces(id)?
    } else {
        Vec::new()
    };
    let mut literal = false;
    for face in faces {
        let surface = helical::face(shell, face)?.0;
        if shell.simple(surface)?.name == "B_SPLINE_SURFACE_WITH_KNOTS" {
            let args = shell.args(surface, "B_SPLINE_SURFACE_WITH_KNOTS")?;
            // A degree-eleven support may have an independent literal reader.
            // This dispatch cue grants no geometric admission.
            if matches!(args.get(1),Some(Value::Number(d)) if *d==11.)
                && matches!(args.get(2),Some(Value::Number(d)) if *d==11.)
            {
                literal = true;
            }
            if matches!(args.get(1),Some(Value::Number(d)) if *d==9.)
                && matches!(args.get(2),Some(Value::Number(d)) if *d==9.)
            {
                return radial_boundary::read_member(shell, id)
                    .map(|body| SourceImport::Radial(Box::new(body)));
            }
        }
    }
    match compound::read_import_member(shell, id, may_be_warped) {
        Ok(body) => Ok(SourceImport::Analytic(Box::new(body))),
        // Preserve existing successful analytic routes (including a finite
        // curved screw cavity). A degree-eleven mismatch is handed to the
        // independent literal verifier, which still must earn every premise.
        Err(ReadError::Invalid) if literal => Err(ReadError::Unsupported),
        Err(error) => Err(error),
    }
}
