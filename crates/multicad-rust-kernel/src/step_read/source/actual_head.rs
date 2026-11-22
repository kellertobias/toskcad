//! Narrow public descriptor of an independently qualified WHOLE original file.
//! Cues below only select this reader; they never authorize source geometry.
use super::*;
use crate::step_read::{
    npt_stock::full_head::{OwnedHeadOptions, transport::MappedActualHead},
    owned_insert::CompoundDocument,
};
use std::sync::Arc;
pub(super) const DESCRIPTOR_WORK: usize = 128;
pub(super) const DESCRIPTOR_BYTES: usize =
    64 + std::mem::size_of::<SourceImport>() + std::mem::size_of::<Vec<SourceImport>>();
/// Explicit lowerable per-kind import boundary. The original entire input is
/// retained and independently decoded. Imported geometry remains STEP-export
/// unsupported until its separate actual-output/source proof is earned.
pub fn read_actual_mapped_head_source(
    input: Arc<str>,
    mut options: OwnedHeadOptions,
    interrupt: Option<&dyn Fn() -> bool>,
) -> std::result::Result<SourceImport, crate::step_read::npt_stock::Error> {
    use crate::step_read::npt_stock::Error as E;
    if interrupt.is_some_and(|f| f()) {
        return Err(E::Cancelled);
    }
    if options.resources.maximum_work == 0
        || options.resources.maximum_work > 272_000_000
        || options.resources.maximum_bytes == 0
        || options.resources.maximum_bytes > 64 * 1024 * 1024
    {
        return Err(E::Invalid);
    }
    options.resources.maximum_work = options
        .resources
        .maximum_work
        .checked_sub(DESCRIPTOR_WORK)
        .filter(|n| *n > 0)
        .ok_or(E::Limit)?;
    options.resources.additional_resident_bytes = options
        .resources
        .additional_resident_bytes
        .checked_add(DESCRIPTOR_BYTES)
        .ok_or(E::Limit)?;
    let actual = MappedActualHead::read(input, options, interrupt)?;
    if interrupt.is_some_and(|f| f()) {
        return Err(E::Cancelled);
    }
    Ok(SourceImport::ActualMappedHead(Arc::new(actual)))
}
fn step(work: &mut usize, amount: usize, maximum: usize) -> Result<()> {
    *work = work
        .checked_add(amount)
        .filter(|n| *n <= maximum)
        .ok_or(ReadError::Limit)?;
    Ok(())
}
fn simple<'a, 'input>(
    document: &'a CompoundDocument<'input>,
    id: u32,
    work: &mut usize,
    maximum: usize,
) -> Result<&'a Entity<'input>> {
    // Bounded document cardinality makes128 cover actual map comparisons,
    // component/argument reads and dispatch bookkeeping before each lookup.
    step(work, 128, maximum)?;
    let parts = document.records().get(&id).ok_or(ReadError::Invalid)?;
    let [entity] = parts.as_slice() else {
        return Err(ReadError::Unsupported);
    };
    Ok(entity)
}
/// Fresh discovery work only. No AST/Parsed/source object is fabricated, and
/// the caller drops this entire arena before the actual whole-file reader.
pub(super) fn recognized(document: &CompoundDocument<'_>, maximum: usize) -> Result<(bool, usize)> {
    let mut work = 0;
    for member in document.represented_members() {
        let solid = simple(document, member.local_solid_entity(), &mut work, maximum)?;
        if solid.name != "MANIFOLD_SOLID_BREP" {
            continue;
        }
        let [Value::String(_), Value::Reference(shell)] = solid.args.as_slice() else {
            return Err(ReadError::Invalid);
        };
        let shell = simple(document, *shell, &mut work, maximum)?;
        if shell.name != "CLOSED_SHELL" {
            continue;
        }
        let [Value::String(_), Value::List(faces)] = shell.args.as_slice() else {
            return Err(ReadError::Invalid);
        };
        if faces.len() != 17 {
            continue;
        }
        for face in faces {
            let Value::Reference(face) = face else {
                return Err(ReadError::Invalid);
            };
            let f = simple(document, *face, &mut work, maximum)?;
            if f.name != "ADVANCED_FACE" {
                continue;
            }
            let [Value::String(_), _, Value::Reference(surface), _] = f.args.as_slice() else {
                return Err(ReadError::Invalid);
            };
            let surface = simple(document, *surface, &mut work, maximum)?;
            if surface.name == "B_SPLINE_SURFACE_WITH_KNOTS"
                && matches!(surface.args.get(1), Some(Value::Number(11.)))
                && matches!(surface.args.get(2), Some(Value::Number(11.)))
            {
                return Ok((true, work));
            }
        }
    }
    Ok((false, work))
}

#[cfg(test)]
mod tests;
