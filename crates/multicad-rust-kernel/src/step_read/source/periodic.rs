//! New bounded flat mixed-periodic import. Non-periodic-only files retain the
//! old whole-file routes and their original limits. Borrowed discovery is never
//! admission: every recognized Phi receives the existing independent owner.
use super::*;
use crate::step_read::{
    DocumentOwnershipError,
    owned_insert::{CompoundDocument, DocumentOptions},
    validate_document_topology_ownership,
};
const MEMBERS: usize = 16;
const RESERVED_WORK: usize = 4_000_000_000;
const RETAINED: usize = 512 * 1024 * 1024;
fn document_error(e: crate::step_read::owned_insert::Error) -> ReadError {
    match e {
        crate::step_read::owned_insert::Error::Invalid => ReadError::Invalid,
        crate::step_read::owned_insert::Error::Limit => ReadError::Limit,
        _ => ReadError::Unsupported,
    }
}
fn material_error(e: super::super::PhillipsPeriodicReadError) -> ReadError {
    match e.cause {
        super::super::PhillipsReadError::Invalid => ReadError::Invalid,
        super::super::PhillipsReadError::Limit => ReadError::Limit,
        _ => ReadError::Unsupported,
    }
}
fn cue(shell: &Shell<'_, '_>, id: u32) -> Result<bool> {
    if shell.simple(id)?.name != "MANIFOLD_SOLID_BREP" {
        return Ok(false);
    }
    for face in shell.member_faces(id)? {
        let surface = helical::face(shell, face)?.0;
        // Only a cheap non-admitting dispatch cue. Degree is not a family
        // identity: current Phi and NPT writers both use degree eleven.
        if shell.simple(surface)?.name == "B_SPLINE_SURFACE_WITH_KNOTS" {
            return Ok(true);
        }
    }
    Ok(false)
}
pub(super) fn read(input: &str, tolerance: f64) -> Result<Option<Vec<SourceImport>>> {
    let document = match CompoundDocument::parse(input, DocumentOptions::default(), None) {
        Ok(d) => d,
        Err(_) => return Ok(None),
    };
    let count = document.represented_members().len();
    // This new route is deliberately flat/world-space. Existing mapped NPT
    // and all other previous routes remain authoritative outside this scope.
    if !(2..=MEMBERS).contains(&count)
        || document
            .represented_members()
            .iter()
            .any(|m| m.mapping_link().is_some())
    {
        return Ok(None);
    }
    let mut phi = Vec::new();
    phi.try_reserve_exact(count).map_err(|_| ReadError::Limit)?;
    let mut found = false;
    for member in document.represented_members() {
        let shell = Shell::new(document.records(), tolerance, true);
        let hinted = match cue(&shell, member.local_solid_entity()) {
            Ok(v) => v,
            Err(_) => return Ok(None),
        };
        let recognized = if hinted {
            let selected = match document.local_member_document(
                member.item_entity(),
                Default::default(),
                None,
            ) {
                Ok(v) => v,
                Err(_) => return Ok(None),
            };
            super::super::read_phillips_body(selected.as_str(), tolerance, 1_000_000).is_ok()
        } else {
            false
        };
        found |= recognized;
        phi.push(recognized);
    }
    if !found {
        return Ok(None);
    }
    // Reserve the SUM before entering any material reader. Default neutral
    // document<=512M, ownership<=16M, two packet views+discovery<=33M/member,
    // each Phi17-stage reader224M. Old singleton calls are independently
    // bounded in count by MEMBERS and retain ALL their original contracts;
    // this reservation does not invent a scalar conversion for their work.
    // This new4B staged contract never claims the old16M whole-file succeeded.
    let reader_options = super::super::PhillipsPeriodicReadOptions::default();
    let phi_reservation = reader_options.stage_work_limits.iter().sum::<usize>();
    let view_options = crate::step_read::owned_insert::Options::default();
    let mut reserved = DocumentOptions::default().maximum_work
        + crate::step_read::DocumentOwnershipOptions::default().maximum_work
        + count * (2 * view_options.maximum_work + 1_000_000);
    for &p in &phi {
        reserved = reserved
            .checked_add(if p { phi_reservation } else { 0 })
            .filter(|n| *n <= RESERVED_WORK)
            .ok_or(ReadError::Limit)?;
    }
    let ownership =
        validate_document_topology_ownership(&document, tolerance, Default::default(), None)
            .map_err(|e| match e {
                DocumentOwnershipError::Invalid => ReadError::Invalid,
                DocumentOwnershipError::Limit => ReadError::Limit,
                _ => ReadError::Unsupported,
            })?;
    assert_eq!(ownership.represented_member_count(), count);
    let mut output = Vec::new();
    output
        .try_reserve_exact(count)
        .map_err(|_| ReadError::Limit)?;
    let mut retained = 0usize;
    for (member, periodic) in document.represented_members().iter().zip(phi) {
        let selected = document
            .local_member_document(member.item_entity(), Default::default(), None)
            .map_err(document_error)?;
        let result = if periodic {
            SourceImport::PeriodicPhillips(Box::new(
                super::super::read_periodic_phillips_material(
                    selected.as_str(),
                    reader_options,
                    None,
                )
                .map_err(material_error)?
                .into_material(),
            ))
        } else {
            let mut previous = super::read_source_compound(selected.as_str(), tolerance)?;
            if previous.len() != 1 {
                return Err(ReadError::Invalid);
            }
            previous.pop().expect("one original selected member")
        };
        let bytes = match &result {
            SourceImport::PeriodicPhillips(s) => s.retained_budget_bytes(),
            SourceImport::LiteralNpt(s) => s.retained_budget_bytes(),
            _ => 0,
        };
        retained = retained
            .checked_add(bytes)
            .filter(|n| *n <= RETAINED)
            .ok_or(ReadError::Limit)?;
        output.push(result);
    }
    Ok(Some(output))
}
