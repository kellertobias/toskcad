//! Exact, source-neutral cut target discovery on a current owned native model.
//! Scope and explicit overrides select candidate Parts; cached bounds only
//! prune work, and positive-volume native common geometry decides every hit.
//! Join discovery (`contact`) also counts certified face-to-face contact.
use crate::{Model, ModelError, compiler::error};
pub use multicad_core::model::{CutTargetQuery, CutTargetSelection, MAX_CUT_TARGET_POOL};
use multicad_core::{
    document::{Document, EntityData, EntityId},
    hierarchy,
    model::{self, BodyInput, MAX_BOOLEAN_TARGETS, MAX_BOOLEAN_TOOLS},
    project_file::ContentHash,
    references::{ReferenceAddress, ReferenceReport, Resolution},
    sketch_regions::supporting_face,
};
use multicad_kernel::{BooleanOperation, Bounds, RigidTransform, Solid};
use std::collections::{BTreeMap, BTreeSet};

/// How an extrusion tool's start cap meets the face its profile lies on.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum FaceRelation {
    /// The tool leaves the face: it shares the profile area as a contact face.
    Touches,
    /// The tool sweeps into the material behind the face.
    Enters,
}
fn checkpoint(interrupt: Option<&dyn Fn() -> bool>) -> Result<(), ModelError> {
    if interrupt.is_some_and(|f| f()) {
        Err(error(None, "Cut target query superseded"))
    } else {
        Ok(())
    }
}
fn overlap(a: Bounds, b: Bounds) -> bool {
    (0..3).all(|i| a.max[i] > b.min[i] && b.max[i] > a.min[i])
}
fn touch(a: Bounds, b: Bounds, tolerance: f64) -> bool {
    (0..3).all(|i| a.max[i] + tolerance >= b.min[i] && b.max[i] + tolerance >= a.min[i])
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    (0..3).map(|i| a[i] * b[i]).sum()
}
/// The face a current extrusion tool's profile lies on, and the producer
/// that owns that face. A face-sketch region is a cell of its supporting face
/// (outside areas and holes are never regions) and a face profile is the face
/// itself, so a zero-offset start cap on that plane shares the whole profile
/// area with the face. Anything else returns None and keeps the general test.
fn profile_face(
    document: &Document,
    revision: u64,
    report: &ReferenceReport,
    tool: EntityId,
    tolerance: f64,
) -> Option<(EntityId, ReferenceAddress)> {
    let definition = model::read_extrude(document, tool).ok()?;
    let context = multicad_solver::SketchEvaluationContext::new(document, revision).ok()?;
    let controls = definition
        .parameters
        .resolve_controls(context.parameters())
        .ok()?;
    if controls.start_offset_mm.abs() > tolerance {
        return None;
    }
    if let Some(sketch) = definition.profile.sketch_id() {
        report.sketch_supports.get(&sketch)?.as_ref().ok()?;
        let (plane, link) = supporting_face(document, sketch)?;
        Some((
            link.target,
            ReferenceAddress {
                owner: plane,
                reference: link.geometry.as_ref()?.id,
            },
        ))
    } else {
        let face = definition.profile.face()?;
        let Some(EntityData::Feature { inputs, .. }) = document.get(tool).map(|e| &e.data) else {
            return None;
        };
        let link = inputs.iter().find(|link| {
            link.role == model::PROFILE_FACE_ROLE
                && link.geometry.as_ref().is_some_and(|g| g.id == face)
        })?;
        Some((
            link.target,
            ReferenceAddress {
                owner: tool,
                reference: face,
            },
        ))
    }
}
fn unique(ids: &[EntityId], budget: usize, what: &str) -> Result<BTreeSet<EntityId>, ModelError> {
    if ids.len() > budget {
        return Err(error(
            None,
            format!("Cut target {what} exceeds the {budget}-item budget"),
        ));
    }
    let result: BTreeSet<_> = ids.iter().copied().collect();
    if result.len() != ids.len() || result.iter().any(EntityId::is_nil) {
        return Err(error(
            None,
            format!("Cut target {what} require unique non-nil identities"),
        ));
    }
    Ok(result)
}
impl Model {
    fn tolerance(&self) -> f64 {
        self.policy.linear_tolerance_mm.max(1e-7) * 10.
    }
    /// Relation of a tool to the current face its profile was taken from.
    /// The sketch or face frame must lie on that face plane; the side follows
    /// from the tool's own coplanar planar cap, never from display geometry.
    fn face_relation(
        &self,
        document: &Document,
        revision: u64,
        tool: EntityId,
    ) -> Option<(EntityId, FaceRelation)> {
        let snapshot = self.snapshot();
        let tolerance = self.tolerance();
        let (producer, address) =
            profile_face(document, revision, &snapshot.references, tool, tolerance)?;
        let Resolution::Resolved { key, .. } = &snapshot.references.get(address)?.resolution else {
            return None;
        };
        let (origin, normal) = self.solid(producer)?.face_plane(key.ordinal).ok()??;
        if let Some(sketch) = document
            .get(tool)
            .and_then(|_| model::read_extrude(document, tool).ok())
            .and_then(|definition| definition.profile.sketch_id())
        {
            let (_, frame, _) =
                multicad_core::constraint_sketch::source::Source::read_for_evaluation(
                    document,
                    sketch,
                    revision,
                    Some(&snapshot.references),
                )
                .ok()?;
            let offset = dot(std::array::from_fn(|i| frame.origin[i] - origin[i]), normal);
            if offset.abs() > tolerance || dot(frame.n, normal).abs() < 1. - 1e-8 {
                return None;
            }
        }
        let cap = self
            .solid(tool)?
            .faces()
            .ok()?
            .into_iter()
            .find_map(|face| {
                let (point, n) = face.plane?;
                let offset = dot(std::array::from_fn(|i| point[i] - origin[i]), normal);
                (offset.abs() <= tolerance && dot(n, normal).abs() >= 1. - 1e-8).then_some(n)
            })?;
        Some((
            producer,
            if dot(cap, normal) < 0. {
                FaceRelation::Touches
            } else {
                FaceRelation::Enters
            },
        ))
    }
    /// Certified face-to-face contact of planar tool faces with a target that
    /// has zero common volume. The tool moves a small distance across each of
    /// its planar faces: a contact face gains volume linearly with the move,
    /// while edge or vertex contact gains none or only higher-order volume.
    fn touches_face_to_face(&self, target: &Solid, tool: &Solid) -> Result<bool, String> {
        let step = (self.policy.minimum_feature_mm * 10.).max(self.tolerance() * 100.);
        let mut normals: Vec<[f64; 3]> = vec![];
        for face in tool.faces().map_err(|e| e.to_string())? {
            let Some((_, n)) = face.plane else { continue };
            if normals.iter().any(|m| dot(*m, n) >= 1. - 1e-8) {
                continue;
            }
            normals.push(n);
            let moved = |distance: f64| -> Result<f64, String> {
                let shifted = tool
                    .transformed(RigidTransform {
                        translation_mm: n.map(|c| c * distance),
                        ..Default::default()
                    })
                    .map_err(|e| e.to_string())?;
                if target.uses_rust_geometry() && shifted.uses_rust_geometry() {
                    target.rust_occupied_overlap_volume(&shifted, self.policy)
                } else {
                    target
                        .boolean(&shifted, BooleanOperation::Intersect, self.policy)
                        .and_then(|solid| solid.map(|s| s.volume_mm3()).transpose())
                        .map(|volume| volume.unwrap_or(0.))
                }
                .map_err(|e| e.to_string())
            };
            let near = moved(step)?;
            if near > 0. && moved(step * 2.)? < near * 3. {
                return Ok(true);
            }
        }
        Ok(false)
    }
    pub fn cut_targets(
        &self,
        document: &Document,
        revision: u64,
        query: &CutTargetQuery,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<CutTargetSelection, ModelError> {
        checkpoint(interrupt)?;
        let snapshot = self.snapshot();
        if snapshot.scope != crate::EvaluationScope::Full {
            return Err(error(
                None,
                "A source-only inspection model cannot discover cut targets",
            ));
        }
        if snapshot.stamp.document != document.id()
            || snapshot.stamp.revision != revision
            || snapshot.stamp.source_hash
                != ContentHash::document(document).map_err(|e| error(None, e))?
        {
            return Err(error(
                None,
                "Cut target query requires the current document/revision/source model; rebuild before selecting targets",
            ));
        }
        let tools = unique(&query.tools, MAX_BOOLEAN_TOOLS, "tools")?;
        if tools.is_empty() {
            return Err(error(None, "Choose at least one current tool producer"));
        }
        let includes = unique(&query.include, MAX_CUT_TARGET_POOL, "includes")?;
        let excludes = unique(&query.exclude, MAX_CUT_TARGET_POOL, "excludes")?;
        let mut tool_heads = BTreeMap::new();
        for (i, (body, producer)) in snapshot.bodies.iter().enumerate() {
            if i % 256 == 0 {
                checkpoint(interrupt)?;
            }
            if tools.contains(producer) && tool_heads.insert(*producer, *body).is_some() {
                return Err(error(
                    Some(*producer),
                    "Tool must uniquely identify one current unconsumed body's final producer",
                ));
            }
        }
        let mut tool_bodies = BTreeSet::new();
        for tool in &tools {
            let Some(body) = tool_heads.get(tool) else {
                return Err(error(
                    Some(*tool),
                    "Tool must uniquely identify one current unconsumed body's final producer",
                ));
            };
            if self.solid(*tool).is_none() || !snapshot.shapes.contains_key(tool) {
                return Err(error(
                    Some(*tool),
                    "Tool producer has no current exact native shape",
                ));
            }
            tool_bodies.insert(*body);
        }
        let index = hierarchy::Index::new(document).map_err(|e| error(None, e))?;
        let defaults: BTreeSet<_> = index
            .default_cut_bodies(query.component)
            .map_err(|e| error(query.component, e))?
            .into_iter()
            .filter(|id| snapshot.bodies.contains_key(id) && !tool_bodies.contains(id))
            .collect();
        let mut explicit = BTreeSet::new();
        let mut empty_components = BTreeSet::new();
        for id in includes {
            checkpoint(interrupt)?;
            match document.get(id).map(|e| &e.data) {
                Some(EntityData::Part { .. }) => {
                    if tool_bodies.contains(&id) {
                        return Err(error(
                            Some(id),
                            "A tool Part cannot explicitly include itself as a cut target",
                        ));
                    }
                    if !snapshot.bodies.contains_key(&id) {
                        return Err(error(
                            Some(id),
                            "Included Part has no current unconsumed final body output",
                        ));
                    }
                    explicit.insert(id);
                }
                Some(EntityData::Assembly { .. }) => {
                    let members: Vec<_> = index
                        .default_cut_bodies(Some(id))
                        .map_err(|e| error(Some(id), e))?
                        .into_iter()
                        .filter(|id| snapshot.bodies.contains_key(id) && !tool_bodies.contains(id))
                        .collect();
                    if members.is_empty() {
                        empty_components.insert(id);
                    }
                    explicit.extend(members);
                }
                _ => {
                    return Err(error(
                        Some(id),
                        "Cut target includes must identify an existing Part or component",
                    ));
                }
            }
        }
        for id in &excludes {
            if !matches!(
                document.get(*id).map(|e| &e.data),
                Some(EntityData::Part { .. })
            ) {
                return Err(error(
                    Some(*id),
                    "Cut target excludes must identify existing Parts",
                ));
            }
        }
        let pool: BTreeSet<_> = defaults.union(&explicit).copied().collect();
        if pool.len() > MAX_CUT_TARGET_POOL {
            return Err(error(
                None,
                format!(
                    "Cut target pool exceeds {MAX_CUT_TARGET_POOL} bodies; narrow its component scope or explicit includes"
                ),
            ));
        }
        let relations: BTreeMap<_, _> = tools
            .iter()
            .filter_map(|tool| {
                self.face_relation(document, revision, *tool)
                    .map(|relation| (*tool, relation))
            })
            .collect();
        let mut hits = BTreeMap::new();
        let mut producer_hits = BTreeMap::new();
        for body in pool {
            checkpoint(interrupt)?;
            let producer = snapshot.bodies[&body];
            if let Some(hit) = producer_hits.get(&producer) {
                hits.insert(body, *hit);
                continue;
            }
            let shape = snapshot.shapes.get(&producer).ok_or_else(|| {
                error(
                    Some(producer),
                    "Current target producer has no exact shape snapshot",
                )
            })?;
            let target = self.solid(producer).ok_or_else(|| {
                error(
                    Some(producer),
                    "Current target producer has no exact native solid",
                )
            })?;
            let mut hit = false;
            for tool in &tools {
                checkpoint(interrupt)?;
                match relations.get(tool) {
                    Some((owner, FaceRelation::Enters)) if *owner == producer => {
                        hit = true;
                        break;
                    }
                    Some((owner, FaceRelation::Touches)) if *owner == producer && query.contact => {
                        hit = true;
                        break;
                    }
                    _ => {}
                }
                let tool_shape = &snapshot.shapes[tool];
                if !overlap(shape.bounds, tool_shape.bounds) {
                    if query.contact
                        && touch(shape.bounds, tool_shape.bounds, self.tolerance())
                        && self
                            .touches_face_to_face(target, self.solid(*tool).expect("validated tool"))
                            .map_err(|e| {
                                error(
                                    Some(body),
                                    format!(
                                        "Cannot certify face contact with this body ({e}); choose explicit targets"
                                    ),
                                )
                            })?
                    {
                        hit = true;
                        break;
                    }
                    continue;
                }
                let tool = self.solid(*tool).expect("validated tool");
                let overlaps = if target.uses_rust_geometry() && tool.uses_rust_geometry() {
                    target
                        .rust_occupied_overlap_volume(tool, self.policy)
                        .map(|volume| volume > 0.)
                } else {
                    target
                        .boolean(tool, BooleanOperation::Intersect, self.policy)
                        .map(|solid| solid.is_some())
                }
                .map_err(|e| error(Some(body), e))?;
                checkpoint(interrupt)?;
                if overlaps {
                    hit = true;
                    break;
                }
            }
            producer_hits.insert(producer, hit);
            hits.insert(body, hit);
        }
        let body_input = |id: &EntityId| BodyInput {
            body: *id,
            producer: snapshot.bodies[id],
        };
        let default_hits = defaults
            .iter()
            .filter(|id| hits[id])
            .map(body_input)
            .collect();
        let targets: Vec<_> = hits
            .iter()
            .filter(|(id, hit)| **hit && !excludes.contains(id))
            .map(|(id, _)| body_input(id))
            .collect();
        if targets.len() > MAX_BOOLEAN_TARGETS {
            return Err(error(
                None,
                format!(
                    "More than {MAX_BOOLEAN_TARGETS} exact cut targets remain; narrow the scope or exclude bodies before applying"
                ),
            ));
        }
        let mut not_intersecting = empty_components;
        not_intersecting.extend(explicit.into_iter().filter(|id| !hits[id]));
        checkpoint(interrupt)?;
        Ok(CutTargetSelection {
            defaults: default_hits,
            targets,
            not_intersecting: not_intersecting.into_iter().collect(),
        })
    }
}
