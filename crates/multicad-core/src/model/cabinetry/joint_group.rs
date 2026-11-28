//! Joined panels form one component (CAD-432). Creating a joint puts its
//! panels and the joint into one canonical Assembly component in the model
//! hierarchy: the component one of its panels already belongs to, or a new
//! "Joined panels N" component placed where the first panel was. Panels
//! that belong to another component stay there. Cabinet-generated joints
//! and joints between carcass or Wood Box panels are organised by their
//! owning recipe and are not grouped.
//!
//! Deleting a joint keeps its panels where they are. With
//! [`JointGroupOnDelete::Keep`] the component stays; with
//! [`JointGroupOnDelete::Ungroup`] a component left without joints is
//! dissolved and its members move to its parent in its place.
use super::*;
use crate::{
    commands::CommandError,
    document::Mutation,
    hierarchy::{Hierarchy, Index},
    model::common::error,
};
use serde::{Deserialize, Serialize};

pub const JOINED_PANELS: &str = "Joined panels";

/// What deleting a joint does with the component holding its panels.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum JointGroupOnDelete {
    /// The component and its panels stay as they are.
    #[default]
    Keep,
    /// A component left without joints is dissolved into its parent.
    Ungroup,
}

fn is_assembly(document: &Document, id: EntityId) -> bool {
    matches!(
        document.get(id).map(|e| &e.data),
        Some(EntityData::Assembly { .. })
    )
}
/// The container listing `id` explicitly, if any.
fn explicit_parent(state: &Hierarchy, id: EntityId) -> Option<EntityId> {
    state
        .members
        .iter()
        .find(|(_, members)| members.contains(&id))
        .map(|(parent, _)| *parent)
}
fn unique_name(document: &Document) -> String {
    let taken: BTreeSet<&str> = document.entities().map(|e| e.name.as_str()).collect();
    (1..)
        .map(|n| format!("{JOINED_PANELS} {n}"))
        .find(|name| !taken.contains(name.as_str()))
        .expect("an unused name")
}
/// Whether the joint's panels are organised by a carcass or Wood Box.
fn recipe_owned(document: &Document, joint: &PanelJointDefinition) -> bool {
    joint.generated_by.is_some()
        || joint.pairs.iter().flat_map(|p| [p.a, p.b]).any(|part| {
            panel_feature(document, part)
                .and_then(|f| read_panel(document, f).ok())
                .is_some_and(|p| {
                    matches!(
                        p.source,
                        PanelSource::Carcass { .. } | PanelSource::WoodBox { .. }
                    )
                })
        })
}

/// Put a newly inserted joint and its panels into one component. Appends
/// the hierarchy edits to `edits` and applies them to `document`.
pub(super) fn group_joined_panels(
    document: &mut Document,
    edits: &mut Vec<Mutation>,
    joint_id: EntityId,
    joint: &PanelJointDefinition,
) -> Result<(), CommandError> {
    if recipe_owned(document, joint) {
        return Ok(());
    }
    let mut panels: Vec<EntityId> = vec![];
    for part in joint.pairs.iter().flat_map(|p| [p.a, p.b]) {
        if !panels.contains(&part) {
            panels.push(part);
        }
    }
    let index = Index::new(document).map_err(|e| error(e.to_string()))?;
    let parents = panels
        .iter()
        .map(|p| index.parent(*p).map_err(|e| error(e.to_string())))
        .collect::<Result<Vec<_>, _>>()?;
    let mut state = document.hierarchy().clone();
    let root = document.id();
    let mut steps = vec![];
    let component = match parents
        .iter()
        .flatten()
        .copied()
        .find(|p| is_assembly(document, *p))
    {
        Some(existing) => existing,
        None => {
            // A new component where the first panel is now.
            let entity = Entity::new(
                unique_name(document),
                EntityData::Assembly { instances: vec![] },
            );
            let id = entity.id;
            let container = parents[0].unwrap_or(root);
            let list = state.members.entry(container).or_default();
            match list.iter().position(|m| *m == panels[0]) {
                Some(i) => list.insert(i, id),
                None => list.push(id),
            }
            state.members.entry(id).or_default();
            steps.push(Mutation::Insert(entity));
            id
        }
    };
    // Panels outside any component join it; members of another one stay.
    for (panel, parent) in panels.iter().zip(&parents) {
        if *parent == Some(component) || parent.is_some_and(|p| is_assembly(document, p)) {
            continue;
        }
        if let Some(from) = explicit_parent(&state, *panel) {
            state
                .members
                .get_mut(&from)
                .expect("listed parent")
                .retain(|m| m != panel);
        }
        state.members.entry(component).or_default().push(*panel);
    }
    state.members.entry(component).or_default().push(joint_id);
    steps.push(Mutation::SetHierarchy(state));
    *document = document.edited(&steps)?;
    edits.extend(steps);
    Ok(())
}

/// Detach a joint that is being deleted from the hierarchy and apply the
/// chosen group behaviour. Appends edits and applies them to `document`.
pub(super) fn ungroup_deleted_joint(
    document: &mut Document,
    edits: &mut Vec<Mutation>,
    joint_id: EntityId,
    group: JointGroupOnDelete,
) -> Result<(), CommandError> {
    let mut state = document.hierarchy().clone();
    let component = explicit_parent(&state, joint_id);
    for members in state.members.values_mut() {
        members.retain(|m| *m != joint_id);
    }
    state.hidden.remove(&joint_id);
    let mut steps = vec![];
    if group == JointGroupOnDelete::Ungroup
        && let Some(component) = component.filter(|c| is_assembly(document, *c))
    {
        let members = state.members.get(&component).cloned().unwrap_or_default();
        let name = document.get(component).map_or("", |e| e.name.as_str());
        if members.iter().any(|m| {
            document
                .get(*m)
                .is_some_and(|e| is_kind(e, PANEL_JOINT_KIND))
        }) {
            return Err(error(format!(
                "Other joints still join panels in {name}; keep the group or delete those joints first"
            )));
        }
        if matches!(document.get(component).map(|e| &e.data), Some(EntityData::Assembly { instances }) if !instances.is_empty())
        {
            return Err(error(format!(
                "{name} holds placed instances; keep the group and move them first"
            )));
        }
        let index = Index::new(document).map_err(|e| error(e.to_string()))?;
        let parent = index
            .parent(component)
            .map_err(|e| error(e.to_string()))?
            .unwrap_or(document.id());
        state.members.remove(&component);
        state.hidden.remove(&component);
        if state.active_component == Some(component) {
            state.active_component = None;
        }
        let list = state.members.entry(parent).or_default();
        let at = list.iter().position(|m| *m == component);
        list.retain(|m| *m != component);
        let at = at.unwrap_or(list.len());
        for (k, m) in members.into_iter().enumerate() {
            list.insert(at + k, m);
        }
        steps.push(Mutation::SetHierarchy(state));
        steps.push(Mutation::Remove(component));
    } else if &state != document.hierarchy() {
        steps.push(Mutation::SetHierarchy(state));
    }
    if !steps.is_empty() {
        *document = document.edited(&steps)?;
        edits.extend(steps);
    }
    Ok(())
}
