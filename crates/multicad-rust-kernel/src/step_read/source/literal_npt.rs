//! Strict original-source fallback; no degree or header cue admits geometry.
use super::*;
use crate::step_read::{
    npt_stock as npt,
    owned_insert::{CompoundDocument, DocumentOptions},
};
fn failure(cause: npt::Error) -> ReadError {
    match cause {
        npt::Error::Invalid => ReadError::Invalid,
        npt::Error::Limit => ReadError::Limit,
        npt::Error::Unsupported | npt::Error::Cancelled => ReadError::Unsupported,
    }
}
fn document_failure(cause: crate::step_read::owned_insert::Error) -> ReadError {
    match cause {
        crate::step_read::owned_insert::Error::Invalid => ReadError::Invalid,
        crate::step_read::owned_insert::Error::Limit => ReadError::Limit,
        _ => ReadError::Unsupported,
    }
}
fn candidate(shell: &Shell<'_, '_>, solid: u32) -> Result<bool> {
    if shell.simple(solid)?.name != "MANIFOLD_SOLID_BREP" {
        return Ok(false);
    }
    for face in shell.member_faces(solid)? {
        let surface = helical::face(shell, face)?.0;
        if shell.simple(surface)?.name == "B_SPLINE_SURFACE_WITH_KNOTS" {
            let args = shell.args(surface, "B_SPLINE_SURFACE_WITH_KNOTS")?;
            if matches!(args.get(1), Some(Value::Number(11.)))
                && matches!(args.get(2), Some(Value::Number(11.)))
            {
                return Ok(true);
            }
        }
    }
    Ok(false)
}
pub(super) fn read(input: &str, tolerance: f64) -> Result<Vec<SourceImport>> {
    // The previous ordinary discovery probe has already ended and dropped its
    // arena. This one actual neutral arena is shared by all selected members.
    // Its historical work/bytes remain separately available from each retained
    // literal source; packet construction never resets an inherited ledger.
    let document = CompoundDocument::parse(input, DocumentOptions::default(), None)
        .map_err(document_failure)?;
    let (topology_count, may_be_warped) =
        compound::inspect_literal_entities(document.records().values().flatten())?;
    let mut owners = BTreeMap::new();
    let mut result = Vec::with_capacity(document.represented_members().len());
    let mut found = false;
    // New fallback's retained NPT members share a bounded aggregate512MiB and
    //512M source-work envelope. Ordinary member recognizers retain their
    // existing contracts, independently of this new source-owned phase.
    let mut retained = 0usize;
    let mut work = document.work_used();
    for member in document.represented_members() {
        let id = member.local_solid_entity();
        let mut shell = Shell::new(document.records(), tolerance, true);
        // Preserve every successful existing recognizer, including other
        // degree-eleven families, before selecting the strict literal reader.
        let ordinary = if member.mapping_link().is_none() {
            ordinary_member(&mut shell, id, may_be_warped)
        } else {
            Err(ReadError::Unsupported)
        };
        let imported = match ordinary {
            Ok(imported) => imported,
            Err(ReadError::Unsupported) if candidate(&shell, id)? => {
                let remaining = 512_000_000usize
                    .checked_sub(work)
                    .filter(|v| *v > 0)
                    .ok_or(ReadError::Limit)?;
                let selected = npt::Selected::parse(
                    &document,
                    member.item_entity(),
                    npt::Options {
                        maximum_work: remaining.min(16_000_000),
                        ..Default::default()
                    },
                    None,
                )
                .map_err(failure)?;
                let source = selected
                    .retain(
                        npt::MaterialOptions {
                            maximum_work: remaining.min(64_000_000),
                            ..Default::default()
                        },
                        None,
                    )
                    .map_err(failure)?;
                work = work
                    .checked_add(source.work_used())
                    .filter(|v| *v <= 512_000_000)
                    .ok_or(ReadError::Limit)?;
                retained = retained
                    .checked_add(source.retained_budget_bytes())
                    .filter(|v| *v <= 512 * 1024 * 1024)
                    .ok_or(ReadError::Limit)?;
                found = true;
                SourceImport::LiteralNpt(Box::new(source))
            }
            Err(error) => return Err(error),
        };
        compound::claim_topology(&shell, id, &mut owners)?;
        result.push(imported);
    }
    if !found {
        return Err(ReadError::Unsupported);
    }
    if owners.len() != topology_count {
        return Err(ReadError::Invalid);
    }
    Ok(result)
}
