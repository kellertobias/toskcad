//! Snap-time projection (CAD-383) composes several independently prepared
//! reference captures of one sketch into one transaction. Every capture was
//! checked against the same base document; composing keeps exactly their
//! additions and refuses anything else, so no capture can silently change
//! existing source.
use super::source::invalid;
use crate::{
    constraint_sketch::source::Source,
    document::{Document, DocumentError, EntityData, EntityId, Mutation},
};
use std::collections::BTreeSet;

/// Mutations that apply every capture preview of `sketch` to `base` at once:
/// one insert per new reference entity and one replacement of the sketch
/// whose source and inputs gain the union of the captured additions.
///
/// A preview may only insert new entities and append to the sketch. Changed
/// or removed entities, edited existing sketch elements, conflicting
/// identities and a different sketch frame are refused.
pub fn merge_capture_previews(
    base: &Document,
    sketch: EntityId,
    previews: &[&Document],
) -> Result<Vec<Mutation>, DocumentError> {
    if previews.is_empty() {
        return Err(invalid("No reference capture to compose"));
    }
    let entity = base.get(sketch).ok_or(DocumentError::Missing(sketch))?;
    let EntityData::Sketch {
        definition: old,
        inputs: old_inputs,
        ..
    } = &entity.data
    else {
        return Err(invalid("Choose a sketch"));
    };
    let original = Source::from_definition(old)?;
    let mut merged = original.clone();
    let mut entity = entity.clone();
    let mut schema = old.schema_version;
    let mut inserts = Vec::new();
    let mut inserted = BTreeSet::new();
    let mut ids: BTreeSet<EntityId> = original
        .sketch
        .points
        .iter()
        .map(|p| p.id)
        .chain(original.sketch.curves.iter().map(|c| c.id))
        .chain(original.sketch.constraints.iter().map(|c| c.id))
        .collect();
    for preview in previews {
        if preview.id() != base.id() {
            return Err(invalid("A reference capture belongs to another document"));
        }
        for old in base.entities() {
            if old.id != sketch && preview.get(old.id) != Some(old) {
                return Err(invalid(
                    "A reference capture changed existing project entities; nothing was composed",
                ));
            }
        }
        for new in preview.entities() {
            if base.get(new.id).is_some() {
                continue;
            }
            if !inserted.insert(new.id) {
                return Err(invalid("Two reference captures share one identity"));
            }
            inserts.push(Mutation::Insert(new.clone()));
        }
        let Some(EntityData::Sketch {
            definition,
            inputs,
            plane,
            ..
        }) = preview.get(sketch).map(|e| &e.data)
        else {
            return Err(invalid("A reference capture removed its sketch"));
        };
        if !matches!(&entity.data, EntityData::Sketch { plane: base_plane, .. } if base_plane == plane)
        {
            return Err(invalid("A reference capture changed the sketch plane"));
        }
        let captured = Source::from_definition(definition)?;
        if captured.u != original.u
            || captured.sketch.texts != original.sketch.texts
            || !captured.sketch.points.starts_with(&original.sketch.points)
            || !captured.sketch.curves.starts_with(&original.sketch.curves)
            || !captured
                .sketch
                .constraints
                .starts_with(&original.sketch.constraints)
            || !inputs.starts_with(old_inputs)
        {
            return Err(invalid(
                "A reference capture edited existing sketch source; nothing was composed",
            ));
        }
        let mut claim = |id: EntityId| {
            if ids.insert(id) {
                Ok(())
            } else {
                Err(invalid("Two reference captures share one sketch identity"))
            }
        };
        for point in &captured.sketch.points[original.sketch.points.len()..] {
            claim(point.id)?;
            merged.sketch.points.push(point.clone());
        }
        for curve in &captured.sketch.curves[original.sketch.curves.len()..] {
            claim(curve.id)?;
            merged.sketch.curves.push(curve.clone());
        }
        for constraint in &captured.sketch.constraints[original.sketch.constraints.len()..] {
            claim(constraint.id)?;
            merged.sketch.constraints.push(constraint.clone());
        }
        let EntityData::Sketch { inputs: target, .. } = &mut entity.data else {
            unreachable!("checked sketch")
        };
        for link in &inputs[old_inputs.len()..] {
            if target.contains(link) {
                return Err(invalid("Two reference captures share one sketch input"));
            }
            target.push(link.clone());
        }
        schema = schema.max(definition.schema_version);
    }
    merged.validate()?;
    let EntityData::Sketch { definition, .. } = &mut entity.data else {
        unreachable!("checked sketch")
    };
    *definition = merged.definition()?;
    definition.schema_version = definition.schema_version.max(schema);
    inserts.push(Mutation::Replace(entity));
    Ok(inserts)
}
