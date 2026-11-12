//! Flat exchange members retain their actual kind of qualification evidence.
use super::*;
#[cfg(test)]
mod tests;

/// Finite radial ownership still needs cold surface, ordering and occupied
/// material qualification. It deliberately carries no analytic volume receipt.
pub enum SourceImport {
    Analytic(Box<(Readback, ImportSolid)>),
    Radial(Box<crate::radial_tensor::BoundaryOwnership>),
}

/// Decode bounded flat mixed representations without inferring a construction
/// recipe from retained tensor controls. Global topology ownership is shared
/// across every member, including analytic and retained members together.
pub fn read_source_compound(input: &str, tolerance: f64) -> Result<Vec<SourceImport>> {
    compound::read_members_by(input, tolerance, true, false, |shell, id, may_be_warped| {
        let faces = if shell.simple(id)?.name == "MANIFOLD_SOLID_BREP" {
            shell.member_faces(id)?
        } else {
            Vec::new()
        };
        for face in faces {
            let surface = helical::face(shell, face)?.0;
            if shell.simple(surface)?.name == "B_SPLINE_SURFACE_WITH_KNOTS" {
                let args = shell.args(surface, "B_SPLINE_SURFACE_WITH_KNOTS")?;
                if matches!(args.get(1),Some(Value::Number(d)) if *d==9.)
                    && matches!(args.get(2),Some(Value::Number(d)) if *d==9.)
                {
                    return radial_boundary::read_member(shell, id)
                        .map(|body| SourceImport::Radial(Box::new(body)));
                }
            }
        }
        compound::read_import_member(shell, id, may_be_warped)
            .map(|body| SourceImport::Analytic(Box::new(body)))
    })
}
