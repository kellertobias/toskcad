//! World linework of the fixed plan symbol a form is authoring or editing,
//! composed exactly as plan drawings compose it: sketch frame, then the
//! authored alignment, then each placed instance's pose.
use super::forms::{Form, Operation};
use multicad_core::{
    commands::Session,
    document::{EntityData, EntityId},
    model::RigidPose,
    sketch::Geometry,
    symbols::{self, PlanAlignment},
};
use multicad_precision::PortableMath as _;
use std::f64::consts::TAU;

/// Each painted polyline in world millimeters; empty when nothing applies.
pub(super) fn world_lines(form: &Form, session: &Session) -> Result<Vec<Vec<[f64; 3]>>, String> {
    let document = session.document();
    let placed = |definition: EntityId,
                  sketch: EntityId,
                  alignment: PlanAlignment|
     -> Result<Vec<(EntityId, RigidPose)>, String> {
        let mut out = vec![];
        for (_, instance) in symbols::instances(document, definition).map_err(|e| e.to_string())? {
            let Some(EntityData::Instance { placement, .. }) =
                document.get(instance.instance).map(|e| &e.data)
            else {
                continue;
            };
            let pose = RigidPose::from_placement(*placement).map_err(|e| e.to_string())?;
            let mapped = *instance
                .ids
                .get(&sketch)
                .ok_or("A placed instance has no copy of the plan sketch")?;
            out.push((
                mapped,
                alignment.pose().then(pose).map_err(|e| e.to_string())?,
            ));
        }
        Ok(out)
    };
    let sketches = match &form.operation {
        // The authored source lives in the symbol's own coordinates.
        Operation::Author(request) => match &request.plan {
            Some(plan) => vec![(plan.sketch, plan.alignment.pose())],
            None => vec![],
        },
        Operation::PlanEdit {
            definition,
            sketch: Some(sketch),
            alignment,
            ..
        } => placed(*definition, *sketch, *alignment)?,
        _ => vec![],
    };
    let mut lines = vec![];
    for (sketch, pose) in sketches {
        let evaluated = multicad_solver::evaluate_document_sketch_with_parameters(
            document,
            session.revision(),
            sketch,
            None,
            session.parameters(),
        )
        .map_err(|e| e.to_string())?;
        let world = |p: [f64; 2]| pose.point(evaluated.frame.point(p));
        for geometry in evaluated.geometry().map_err(|e| e.to_string())? {
            lines.push(match geometry {
                Geometry::Line { a, b, .. } => vec![world(a), world(b)],
                Geometry::Rectangle { a, b, .. } => [a, [b[0], a[1]], b, [a[0], b[1]], a]
                    .into_iter()
                    .map(world)
                    .collect(),
                Geometry::Circle {
                    center, radius_mm, ..
                } => arc(center, radius_mm, 0., TAU)
                    .into_iter()
                    .map(world)
                    .collect(),
                Geometry::Arc {
                    center,
                    a,
                    b,
                    clockwise,
                    ..
                } => {
                    let start = (a[1] - center[1]).portable_atan2(a[0] - center[0]);
                    let mut end = (b[1] - center[1]).portable_atan2(b[0] - center[0]);
                    if clockwise {
                        while end >= start {
                            end -= TAU;
                        }
                    } else {
                        while end <= start {
                            end += TAU;
                        }
                    }
                    let radius = (a[0] - center[0]).portable_hypot(a[1] - center[1]);
                    arc(center, radius, start, end)
                        .into_iter()
                        .map(world)
                        .collect()
                }
                Geometry::BSpline { poles, basis, .. } => basis
                    .display_segments(&poles)
                    .map_err(|e| e.to_string())?
                    .into_iter()
                    .flat_map(|segment| [segment[0], segment[1]])
                    .map(world)
                    .collect(),
            });
        }
    }
    Ok(lines)
}
fn arc(center: [f64; 2], radius: f64, start: f64, end: f64) -> Vec<[f64; 2]> {
    const STEPS: usize = 48;
    (0..=STEPS)
        .map(|i| {
            let t = start + (end - start) * i as f64 / STEPS as f64;
            let (s, c) = t.portable_sin_cos();
            [center[0] + radius * c, center[1] + radius * s]
        })
        .collect()
}
