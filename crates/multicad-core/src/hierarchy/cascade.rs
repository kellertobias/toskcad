//! CAD-426: the complete set of sources one Delete removes. The UI and MCP
//! share this closure, so a confirmed or explicitly cascading delete removes
//! exactly the entries the confirmation named, in one transaction.
//!
//! An entry depends on a removed one when it consumes it through a plain
//! (non-geometry) source link: later features of a part history, Booleans and
//! patterns that take it as an input, instances of a removed part, sketches on
//! a removed plane, and module helpers that only describe removed entries.
//! Links that only capture geometry (a sketch on a face, a projected edge) and
//! drawings or export sets are not deleted: they keep their link and report it
//! through the existing missing-reference diagnostics.
use super::*;
use crate::document::{Entity, Link};

/// The removal one Delete performs.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DeletionPlan {
    pub target: EntityId,
    pub target_name: String,
    /// Every source removed, including the target and the helpers it owns.
    pub removed: BTreeSet<EntityId>,
    /// Model entries removed besides the target, sources first. Helpers that
    /// belong to the target or to a listed entry (its creating operation,
    /// manufacturing properties, placement records) are not listed.
    pub dependents: Vec<Dependent>,
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Dependent {
    pub id: EntityId,
    pub name: String,
    pub feature: bool,
}
impl DeletionPlan {
    fn noun(&self) -> &'static str {
        let features = self.dependents.iter().all(|d| d.feature);
        match (features, self.dependents.len()) {
            (true, 1) => "dependent feature",
            (true, _) => "dependent features",
            (false, 1) => "dependent item",
            (false, _) => "dependent items",
        }
    }
    fn names(&self) -> String {
        self.dependents
            .iter()
            .map(|d| d.name.as_str())
            .collect::<Vec<_>>()
            .join(", ")
    }
    /// The confirmation question, e.g. "Delete Extrude 1 and 2 dependent
    /// features (Fillet 1, Chamfer 1)?".
    pub fn question(&self) -> String {
        if self.dependents.is_empty() {
            return format!("Delete {}?", self.target_name);
        }
        format!(
            "Delete {} and {} {} ({})?",
            self.target_name,
            self.dependents.len(),
            self.noun(),
            self.names()
        )
    }
    /// Why a delete without the cascade option is refused.
    pub fn refusal(&self) -> String {
        let quoted = self
            .dependents
            .iter()
            .map(|d| format!("“{}”", d.name))
            .collect::<Vec<_>>()
            .join(", ");
        let verb = if self.dependents.len() == 1 {
            "depends"
        } else {
            "depend"
        };
        format!(
            "{quoted} {verb} on this entry. Delete it together with its {} {}, or remove or update {} first",
            self.dependents.len(),
            self.noun(),
            if self.dependents.len() == 1 {
                "that dependency"
            } else {
                "those dependencies"
            }
        )
    }
}

fn hard(entity: &Entity) -> impl Iterator<Item = Link> {
    entity
        .references()
        .into_iter()
        .filter(|link| link.geometry.is_none())
}
fn is_extrude_group(entity: &Entity) -> bool {
    matches!(&entity.data, EntityData::Extension { definition, .. }
        if definition.type_id == crate::model::EXTRUDE_GROUP_TYPE)
}

/// Compute the removal closure of `target`. The caller validates that the
/// target itself is a deletable model entry.
pub(super) fn closure(document: &Document, target: EntityId) -> Result<DeletionPlan, String> {
    let index = Index::new(document).map_err(|error| error.to_string())?;
    let mut removed = BTreeSet::from([target]);
    loop {
        let before = removed.len();
        // Ownership: what a removed entry carries with it.
        for id in removed.clone() {
            let Some(entity) = document.get(id) else {
                continue;
            };
            match &entity.data {
                EntityData::Part { features, .. } => removed.extend(features.iter().copied()),
                EntityData::Assembly { .. } | EntityData::Group {} => removed.extend(
                    index
                        .children(Some(id))
                        .map_err(|error| error.to_string())?,
                ),
                EntityData::Instance {
                    placement_source, ..
                } => {
                    let mut driver = *placement_source;
                    while let Some(source) = driver {
                        let Ok(placement) = crate::model::read_instance_placement(document, source)
                        else {
                            break;
                        };
                        if placement.owner != id || !removed.insert(source) {
                            break;
                        }
                        driver = match placement.base {
                            crate::model::PlacementBase::Previous { driver } => Some(driver),
                            crate::model::PlacementBase::Literal { .. } => None,
                        };
                    }
                }
                // An extrusion operation is atomic: its tools and results go
                // together with the operation record.
                _ if is_extrude_group(entity) => removed.extend(
                    hard(entity)
                        .filter(|link| {
                            link.role.ends_with("/producer") || link.role.ends_with("/feature")
                        })
                        .map(|link| link.target),
                ),
                _ => {}
            }
        }
        // Dependency: what consumes a removed entry.
        for entity in document.entities() {
            if removed.contains(&entity.id) {
                continue;
            }
            let consumes = |link: &Link| removed.contains(&link.target);
            let dependent = match &entity.data {
                // A part without its creating operation has no geometry.
                EntityData::Part { features, .. } => {
                    features.first().is_some_and(|id| removed.contains(id))
                }
                EntityData::Feature { .. } => hard(entity).any(|link| consumes(&link)),
                EntityData::Instance { component, .. } => removed.contains(component),
                EntityData::Sketch { plane, .. } => removed.contains(plane),
                EntityData::Extension { .. } if is_extrude_group(entity) => {
                    hard(entity).any(|link| consumes(&link))
                }
                // A module record is removed only when every model entry it
                // describes is removed; otherwise it keeps a lost reference.
                EntityData::Extension { .. } => {
                    let links: Vec<_> = hard(entity).collect();
                    links.iter().any(consumes)
                        && links.iter().all(|link| {
                            consumes(link)
                                || !document
                                    .get(link.target)
                                    .is_some_and(|target| is_model_entry(&target.data))
                        })
                }
                _ => false,
            };
            if dependent {
                removed.insert(entity.id);
            }
        }
        if removed.len() == before {
            break;
        }
    }
    removed.retain(|id| document.get(*id).is_some());

    // List model entries only; creating operations stay with their part.
    let mut hidden = BTreeSet::from([target]);
    for id in &removed {
        if let Some(EntityData::Part { features, .. }) = document.get(*id).map(|e| &e.data)
            && let Some(first) = features.first()
            && *first != target
        {
            hidden.insert(*first);
        }
    }
    let order = document
        .topological_order()
        .unwrap_or_else(|_| document.entities().map(|e| e.id).collect());
    let dependents = order
        .into_iter()
        .filter(|id| removed.contains(id) && !hidden.contains(id))
        .filter_map(|id| document.get(id))
        .filter(|entity| is_model_entry(&entity.data))
        .map(|entity| Dependent {
            id: entity.id,
            name: entity.name.clone(),
            feature: matches!(entity.data, EntityData::Feature { .. }),
        })
        .collect();
    Ok(DeletionPlan {
        target,
        target_name: document
            .get(target)
            .map(|e| e.name.clone())
            .unwrap_or_default(),
        removed,
        dependents,
    })
}
