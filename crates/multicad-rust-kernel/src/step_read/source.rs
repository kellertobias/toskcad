//! Flat exchange members retain their actual kind of qualification evidence.
use super::*;
mod actual_head;
pub use actual_head::read_actual_mapped_head_source;
mod literal_npt;
mod periodic;
#[cfg(test)]
mod tests;

/// Finite radial ownership still needs cold surface, ordering and occupied
/// material qualification. It deliberately carries no analytic volume receipt.
pub enum SourceImport {
    PeriodicPhillips(Box<super::PhillipsOwnedPeriodicMaterial>),
    Analytic(Box<(Readback, ImportSolid)>),
    Radial(Box<crate::radial_tensor::BoundaryOwnership>),
    /// Cylindrical stock around one through radial cavity, both uncertified.
    PipeSleeve(Box<super::ImportedPipeSleeve>),
    LiteralNpt(Box<super::npt_stock::OwnedMaterial>),
    ActualMappedHead(std::sync::Arc<super::npt_stock::full_head::transport::MappedActualHead>),
}

/// Decode bounded flat mixed representations without inferring a construction
/// recipe from retained tensor controls. Global topology ownership is shared
/// across every member, including analytic and retained members together.
pub fn read_source_compound(input: &str, tolerance: f64) -> Result<Vec<SourceImport>> {
    if let Some(periodic) = periodic::read(input, tolerance)? {
        return Ok(periodic);
    }
    match ordinary(input, tolerance) {
        // The literal NPT reader recognises only members with a degree-11
        // B-spline face. This is a candidate/performance cue only; legal
        // representation ownership is independently resolved by its document
        // reader, including non-transform shape-representation relationships.
        Err(ReadError::Unsupported) if has_degree_eleven_surface(input) => {
            literal_npt::read(input, tolerance)
        }
        other => other,
    }
}
/// True when the file states a degree-11 by degree-11 B-spline surface, the
/// necessary premise of every literal NPT candidate member.
fn has_degree_eleven_surface(input: &str) -> bool {
    super::part21::parse(input).is_ok_and(|entities| {
        entities.values().flatten().any(|entity| {
            entity.name == "B_SPLINE_SURFACE_WITH_KNOTS"
                && matches!(entity.args.get(1), Some(Value::Number(d)) if *d == 11.)
                && matches!(entity.args.get(2), Some(Value::Number(d)) if *d == 11.)
        })
    })
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
                // Degree-9 retained tensors only: a cylindrical stock with two
                // annular caps is a pipe sleeve, never a helical receiver
                // (those keep their degree-11 analytic routes below).
                if let Some(sleeve) = radial_boundary::sleeve::recognize(shell, id)? {
                    return Ok(SourceImport::PipeSleeve(Box::new(sleeve)));
                }
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
