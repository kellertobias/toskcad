//! Author-defined fixed plan symbols of generic library parts and assemblies.
//! A plan projection draws the authored linework at each instance's current
//! occurrence pose, composed with the definition's explicit plan alignment,
//! instead of projecting the detailed 3D content. Elevations, sections and 3D
//! views keep the physical geometry.
use crate::{
    Model, ModelError,
    compiler::error,
    drawings::DrawingBody,
    plant_plans::{PlantPlan, project_sketch},
};
use multicad_core::{
    document::{Document, EntityId},
    drawings::*,
    sketch, symbols,
};
use std::collections::{BTreeMap, BTreeSet};

pub(crate) fn plans(
    document: &Document,
    model: &Model,
    revision: u64,
    view: &DrawingView,
    bodies: &[DrawingBody<'_>],
    handled: &BTreeSet<EntityId>,
) -> Result<BTreeMap<EntityId, PlantPlan>, ModelError> {
    if view.kind != ViewKind::Projection
        || view.frame.normal()[2] < 1. - 1e-8
        || !matches!(view.source, ViewSource::Model { .. })
    {
        return Ok(BTreeMap::new());
    }
    let selected: BTreeSet<_> = bodies.iter().map(|b| b.identity).collect();
    let records: Vec<_> = symbols::authored_plan_instances(document)
        .map_err(|e| error(None, e))?
        .into_iter()
        .filter(|plan| selected.contains(&plan.instance) && !handled.contains(&plan.instance))
        .collect();
    if records.is_empty() {
        return Ok(BTreeMap::new());
    }
    let parameters = multicad_core::parameters::evaluate(document).map_err(|e| error(None, e))?;
    let mut output = BTreeMap::new();
    for plan in records {
        let occurrence = model
            .snapshot()
            .occurrences
            .get(&plan.instance)
            .ok_or_else(|| {
                error(
                    Some(plan.instance),
                    "Symbol has no current occurrence pose for its plan",
                )
            })?
            .pose;
        let pose = plan
            .alignment
            .pose()
            .then(occurrence)
            .map_err(|e| error(Some(plan.instance), e))?;
        let evaluated = multicad_solver::evaluate_document_sketch_with_parameters(
            document,
            revision,
            plan.sketch,
            Some(&model.snapshot().references),
            &parameters,
        )
        .map_err(|e| error(Some(plan.instance), e))?;
        // A tilted instance remains a physical projection; the fixed symbol
        // is only substituted when its plane faces the drawing.
        if sketch::dot(pose.vector(evaluated.frame.n), view.frame.normal()) < 1. - 1e-8 {
            continue;
        }
        output.insert(
            plan.instance,
            PlantPlan {
                vectors: project_sketch(view, plan.instance, pose, &evaluated, 1.)?,
                labels: vec![],
            },
        );
    }
    Ok(output)
}
