//! Verify membership in a flat STEP advanced BREP representation.

use super::*;
pub(super) mod ownership;

pub(super) fn read_import_members(
    input: &str,
    tolerance: f64,
    single: bool,
) -> Result<Vec<(Readback, ImportSolid)>> {
    read_members_by(input, tolerance, true, single, read_import_member)
}

pub(super) fn read_import_member(
    shell: &mut Shell<'_, '_>,
    id: u32,
    may_be_warped: bool,
) -> Result<(Readback, ImportSolid)> {
    if let Some(member) = box_void::read_member(shell, id)? {
        return Ok(member);
    }
    let face_ids = shell.member_faces(id)?;
    if let Some(sphere) = sphere::read_member(shell, &face_ids)? {
        return Ok(sphere);
    }
    if let Some(torus) = torus::read_member(shell, &face_ids)? {
        return Ok(torus);
    }
    if let Some(thread) = helical_import::read_member(shell, id, &face_ids)? {
        return Ok(thread);
    }
    if let Some(rounded) = rounded_cylinder::recognize_member(shell, id, &face_ids)? {
        return Ok(rounded);
    }
    if let Some(thread) = rounded_thread::recognize_member(shell, id, &face_ids)? {
        return Ok(thread);
    }
    if let Some(bore) = rounded_thread_blind::recognize_member(shell, id, &face_ids)? {
        return Ok(bore);
    }
    if let Some(bore) = rounded_thread_bore::recognize_member(shell, id, &face_ids)? {
        return Ok(bore);
    }
    if let Some(cut) = import::read_helical_cut_member(shell, &face_ids)? {
        return Ok(cut);
    }
    if let Some((proof, loft)) = nonplanar_loft::read_member(shell, &face_ids)? {
        return Ok(import::recognize_nonplanar_loft(proof, loft));
    }
    if face_ids
        .iter()
        .any(|&face| spline_sweep::has_extrusion_support(shell, face))
    {
        return match spline_sweep::import_member(shell, face_ids) {
            Ok((proof, body)) => Ok(import::recognize_spline(proof, body)),
            Err(error @ (ReadError::Invalid | ReadError::Unsupported)) => {
                spline_cap_relief::recognize_member(shell, id)
                    .map_err(|cause| spline_error(cause, error))
            }
            Err(error) => Err(error),
        };
    }
    let regular = shell
        .validate_faces(face_ids)
        .and_then(|faces| import::recognize_member(shell, &faces));
    let regular = match regular {
        Err(error @ (ReadError::Invalid | ReadError::Unsupported)) if may_be_warped => {
            fallback_warped(shell, id, error)
                .map(|(readback, import)| import::recognize_warped(readback, import))
        }
        _ => regular,
    };
    match regular {
        Err(ReadError::Unsupported) => {
            normal_relief::read_member(shell, id)?.ok_or(ReadError::Unsupported)
        }
        other => other,
    }
}

pub(super) fn read_analytic(input: &str, tolerance: f64, single: bool) -> Result<Vec<Readback>> {
    read_analytic_with_helical_cuts(input, tolerance, single)
        .map(|members| members.into_iter().map(|(proof, _)| proof).collect())
}

pub(super) fn read_analytic_with_helical_cuts(
    input: &str,
    tolerance: f64,
    single: bool,
) -> Result<Vec<(Readback, Option<ImportSolid>)>> {
    read_members_by(input, tolerance, true, single, read_analytic_member)
}

pub(super) fn read_analytic_member(
    shell: &mut Shell<'_, '_>,
    id: u32,
    may_be_warped: bool,
) -> Result<(Readback, Option<ImportSolid>)> {
    if let Some((readback, _)) = box_void::read_member(shell, id)? {
        return Ok((readback, None));
    }
    let faces = shell.member_faces(id)?;
    if let Some((proof, _)) = sphere::read_member(shell, &faces)? {
        return Ok((proof, None));
    }
    if let Some((proof, _)) = torus::read_member(shell, &faces)? {
        return Ok((proof, None));
    }
    if let Some((proof, thread)) = helical_import::read_member(shell, id, &faces)? {
        return Ok((proof, Some(thread)));
    }
    if let Some((proof, _)) = rounded_cylinder::recognize_member(shell, id, &faces)? {
        return Ok((proof, None));
    }
    if let Some((proof, thread)) = rounded_thread::recognize_member(shell, id, &faces)? {
        return Ok((proof, Some(thread)));
    }
    if let Some((proof, bore)) = rounded_thread_blind::recognize_member(shell, id, &faces)? {
        return Ok((proof, Some(bore)));
    }
    if let Some((proof, bore)) = rounded_thread_bore::recognize_member(shell, id, &faces)? {
        return Ok((proof, Some(bore)));
    }
    if let Some((proof, cut)) = import::read_helical_cut_member(shell, &faces)? {
        return Ok((proof, Some(cut)));
    }
    if let Some((proof, _)) = nonplanar_loft::read_member(shell, &faces)? {
        return Ok((proof, None));
    }
    if faces
        .iter()
        .any(|&face| spline_sweep::has_extrusion_support(shell, face))
    {
        return match spline_sweep::read_member(shell, faces) {
            Ok(proof) => Ok((proof, None)),
            Err(error @ (ReadError::Invalid | ReadError::Unsupported)) => {
                spline_cap_relief::recognize_member(shell, id)
                    .map(|(proof, body)| (proof, Some(body)))
                    .map_err(|cause| spline_error(cause, error))
            }
            Err(error) => Err(error),
        };
    }
    let regular = shell
        .validate_faces(faces)
        .and_then(|faces| shell.readback(&faces));
    match regular {
        Err(error @ (ReadError::Invalid | ReadError::Unsupported)) if may_be_warped => {
            fallback_warped(shell, id, error).map(|(proof, _)| (proof, None))
        }
        _ => regular.map(|proof| (proof, None)),
    }
}

// A later provider cannot mask either an earlier malformed-geometry diagnosis
// or its own resource refusal. Both analytic and import entry points share this.
pub(super) fn spline_error(cause: ReadError, original: ReadError) -> ReadError {
    match cause {
        ReadError::Limit => ReadError::Limit,
        _ => original,
    }
}

fn fallback_warped(
    shell: &Shell<'_, '_>,
    id: u32,
    error: ReadError,
) -> Result<(Readback, ImportedWarped)> {
    // The independent verifier does not reuse a partially parsed shell.
    let mut clean = Shell::new(shell.entities, shell.tolerance, true);
    match warped::read_member(&mut clean, id) {
        Ok(result) => Ok(result),
        Err(ReadError::Limit) => Err(ReadError::Limit),
        Err(_) => Err(error),
    }
}

pub(super) fn read_members_by<T>(
    input: &str,
    tolerance: f64,
    allow_circles: bool,
    single: bool,
    member: impl FnMut(&mut Shell<'_, '_>, u32, bool) -> Result<T>,
) -> Result<Vec<T>> {
    read_members_retaining_arena(input, tolerance, allow_circles, single, member)
        .map(|(_, members)| members)
}

/// Retain the actual immutable parsed file when later boundary proofs need
/// its literal tensors/curves. Ownership and orphan-topology checks are exactly
/// the same as the ordinary reader; no virtual geometry is inserted here.
type RetainedMembers<'input, T> = (BTreeMap<u32, Vec<Entity<'input>>>, Vec<T>);
pub(super) fn read_members_retaining_arena<'input, T>(
    input: &'input str,
    tolerance: f64,
    allow_circles: bool,
    single: bool,
    mut member: impl FnMut(&mut Shell<'_, 'input>, u32, bool) -> Result<T>,
) -> Result<RetainedMembers<'input, T>> {
    let mut entities = parse_entities(input, tolerance)?;
    let transmitted = read_arena(&entities, tolerance, allow_circles, single, &mut member);
    let error = match transmitted {
        Ok(members) => return Ok((entities, members)),
        Err(ReadError::Limit) => return Err(ReadError::Limit),
        Err(error) => error,
    };
    // A refused shell may only differ from a recognized family by faces that
    // partition one plane. Only a recognized merged boundary replaces the
    // verdict on the transmitted one.
    let solids = Shell::new(&entities, tolerance, allow_circles).represented_solids(single);
    let merged = match solids {
        Ok(solids) => same_domain::unify(&mut entities, &solids, tolerance),
        Err(_) => false,
    };
    if merged
        && let Ok(members) = read_arena(&entities, tolerance, allow_circles, single, &mut member)
    {
        return Ok((entities, members));
    }
    Err(error)
}

fn read_arena<'input, T>(
    entities: &BTreeMap<u32, Vec<Entity<'input>>>,
    tolerance: f64,
    allow_circles: bool,
    single: bool,
    member: &mut impl FnMut(&mut Shell<'_, 'input>, u32, bool) -> Result<T>,
) -> Result<Vec<T>> {
    let solids = Shell::new(entities, tolerance, allow_circles).represented_solids(single)?;
    let (topology_count, analytic_candidate) = inspect_flat_entities(entities.values().flatten())?;
    reject_reference_cycles(entities, &solids)?;
    let may_be_warped = allow_circles && analytic_candidate;
    let mut owners = BTreeMap::new();
    let members = solids
        .into_iter()
        .map(|id| {
            let mut shell = Shell::new(entities, tolerance, allow_circles);
            let result = member(&mut shell, id, may_be_warped)?;
            claim_topology(&shell, id, &mut owners)?;
            Ok(result)
        })
        .collect::<Result<Vec<_>>>()?;
    if owners.len() != topology_count {
        return Err(ReadError::Invalid);
    }
    Ok(members)
}

/// The topology and geometry below a solid form a hierarchy. A reference
/// back to an entity that is still being expanded is malformed, whatever
/// entity kind the referring attribute expects.
fn reject_reference_cycles(
    entities: &BTreeMap<u32, Vec<Entity<'_>>>,
    solids: &[u32],
) -> Result<()> {
    fn collect(value: &Value<'_>, pending: &mut Vec<(u32, bool)>) {
        match value {
            Value::Reference(id) => pending.push((*id, false)),
            Value::List(items) => items.iter().for_each(|item| collect(item, pending)),
            Value::Typed(entity) => entity.args.iter().for_each(|arg| collect(arg, pending)),
            _ => {}
        }
    }
    let mut closed = BTreeMap::new();
    let mut pending = solids.iter().map(|&id| (id, false)).collect::<Vec<_>>();
    while let Some((id, leaving)) = pending.pop() {
        if leaving {
            closed.insert(id, true);
            continue;
        }
        match closed.get(&id) {
            Some(true) => continue,
            Some(false) => return Err(ReadError::Invalid),
            None => {}
        }
        closed.insert(id, false);
        pending.push((id, true));
        // A dangling reference is reported by the member reader.
        for entity in entities.get(&id).into_iter().flatten() {
            entity
                .args
                .iter()
                .for_each(|arg| collect(arg, &mut pending));
        }
    }
    Ok(())
}

pub(super) fn inspect_flat_entities<'a>(
    entities: impl Iterator<Item = &'a Entity<'a>>,
) -> Result<(usize, bool)> {
    let mut count = 0;
    let mut ellipse = false;
    let mut cylinder = false;
    for entity in entities {
        if matches!(
            entity.name,
            "REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION"
                | "CONTEXT_DEPENDENT_SHAPE_REPRESENTATION"
                | "NEXT_ASSEMBLY_USAGE_OCCURRENCE"
                | "MAPPED_ITEM"
        ) {
            return Err(ReadError::Unsupported);
        }
        count += usize::from(topology(entity.name));
        ellipse |= entity.name == "ELLIPSE";
        cylinder |= entity.name == "CYLINDRICAL_SURFACE";
    }
    Ok((count, ellipse || cylinder))
}

pub(super) fn inspect_literal_entities<'a>(
    entities: impl Iterator<Item = &'a Entity<'a>>,
) -> Result<(usize, bool)> {
    let mut count = 0;
    let mut candidate = false;
    for entity in entities {
        if matches!(
            entity.name,
            "REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION"
                | "CONTEXT_DEPENDENT_SHAPE_REPRESENTATION"
                | "NEXT_ASSEMBLY_USAGE_OCCURRENCE"
        ) {
            return Err(ReadError::Unsupported);
        }
        count += usize::from(topology(entity.name));
        candidate |= matches!(entity.name, "ELLIPSE" | "CYLINDRICAL_SURFACE");
    }
    Ok((count, candidate))
}
pub(super) fn topology(name: &str) -> bool {
    matches!(
        name,
        "MANIFOLD_SOLID_BREP"
            | "CLOSED_SHELL"
            | "ADVANCED_FACE"
            | "FACE_OUTER_BOUND"
            | "FACE_BOUND"
            | "EDGE_LOOP"
            | "ORIENTED_EDGE"
            | "SEAM_EDGE"
            | "EDGE_CURVE"
            | "VERTEX_POINT"
            | "BREP_WITH_VOIDS"
            | "OPEN_SHELL"
            | "ORIENTED_CLOSED_SHELL"
            | "FACE_SURFACE"
            | "VERTEX_LOOP"
            | "SHELL_BASED_SURFACE_MODEL"
    )
}

pub(super) fn claim_topology(
    shell: &Shell<'_, '_>,
    solid: u32,
    owners: &mut BTreeMap<u32, u32>,
) -> Result<()> {
    let mut pending = Vec::new();
    claim_topology_by(shell, solid, owners, &mut pending, &mut || {
        Ok::<_, ReadError>(())
    })
}
fn claim_topology_by<E: From<ReadError>>(
    shell: &Shell<'_, '_>,
    solid: u32,
    owners: &mut BTreeMap<u32, u32>,
    pending: &mut Vec<u32>,
    step: &mut impl FnMut() -> std::result::Result<(), E>,
) -> std::result::Result<(), E> {
    use std::collections::btree_map::Entry;
    pending.clear();
    pending.push(solid);
    while let Some(id) = pending.pop() {
        step()?;
        let parts = shell.components(id).map_err(E::from)?;
        let [entity] = parts else {
            // A verified rational support is a complex geometric entity.
            // Topology remains restricted to single, exclusively owned entities.
            if parts.iter().any(|part| topology(part.name)) {
                return Err(E::from(ReadError::Unsupported));
            }
            continue;
        };
        if !topology(entity.name) {
            continue;
        }
        match owners.entry(id) {
            Entry::Occupied(owner) if *owner.get() == solid => continue,
            Entry::Occupied(_) => return Err(E::from(ReadError::Invalid)),
            Entry::Vacant(owner) => owner.insert(solid),
        };
        // The member shell has already checked every topological reference.
        // Only topology has an exclusive owner; geometric supports can be shared.
        for arg in &entity.args {
            match arg {
                Value::Reference(id) => pending.push(*id),
                Value::List(items) => {
                    pending.extend(items.iter().filter_map(|item| reference(item).ok()));
                }
                _ => {}
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests;
