//! What a drawing view shows (CAD-438). A view shows every part of its
//! sources unless the page holds a filter for it: include tags keep only
//! parts carrying one of them, exclude tags drop parts carrying one, and
//! per-part toggles hide or show single parts over the tag rule. The filter
//! uses the canonical Part tags, so tag edits in Model update the view.
use super::{DrawingView, Page, ViewSource, source::invalid};
use crate::{
    document::{Document, DocumentError, EntityData, EntityId},
    materials::checked_tags,
};
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;

/// Files with view part filters need a reader that applies them.
pub const PART_FILTER_CAPABILITY: &str = "drawings.part_filters.1";

#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ViewPartFilter {
    pub view: EntityId,
    /// Show only parts with at least one of these tags; empty shows all.
    #[serde(default, skip_serializing_if = "BTreeSet::is_empty")]
    pub include_tags: BTreeSet<String>,
    /// Hide parts with any of these tags.
    #[serde(default, skip_serializing_if = "BTreeSet::is_empty")]
    pub exclude_tags: BTreeSet<String>,
    /// Parts toggled off in this view.
    #[serde(default, skip_serializing_if = "BTreeSet::is_empty")]
    pub hidden: BTreeSet<EntityId>,
    /// Parts toggled on in this view although the tag rule hides them.
    #[serde(default, skip_serializing_if = "BTreeSet::is_empty")]
    pub shown: BTreeSet<EntityId>,
}
impl ViewPartFilter {
    pub fn is_empty(&self) -> bool {
        self.include_tags.is_empty()
            && self.exclude_tags.is_empty()
            && self.hidden.is_empty()
            && self.shown.is_empty()
    }
    /// Whether the tag rule alone shows this part.
    pub fn tags_show(&self, document: &Document, part: EntityId) -> bool {
        let empty = BTreeSet::new();
        let tags = match document.get(part).map(|e| &e.data) {
            Some(EntityData::Part { tags, .. }) => tags,
            _ => &empty,
        };
        (self.include_tags.is_empty() || !self.include_tags.is_disjoint(tags))
            && self.exclude_tags.is_disjoint(tags)
    }
    pub fn shows(&self, document: &Document, part: EntityId) -> bool {
        if self.hidden.contains(&part) {
            return false;
        }
        self.shown.contains(&part) || self.tags_show(document, part)
    }
}
impl Page {
    pub fn part_filter(&self, view: EntityId) -> Option<&ViewPartFilter> {
        self.part_filters.iter().find(|f| f.view == view)
    }
    fn part_filter_mut(&mut self, view: EntityId) -> &mut ViewPartFilter {
        if let Some(index) = self.part_filters.iter().position(|f| f.view == view) {
            return &mut self.part_filters[index];
        }
        self.part_filters.push(ViewPartFilter {
            view,
            ..Default::default()
        });
        self.part_filters.last_mut().expect("just pushed")
    }
    fn prune_part_filters(&mut self) {
        self.part_filters.retain(|f| !f.is_empty());
    }
    /// Toggle one part of a view on or off.
    pub fn set_part_shown(
        &mut self,
        document: &Document,
        view: EntityId,
        part: EntityId,
        shown: bool,
    ) -> Result<(), DocumentError> {
        if !self.views.iter().any(|v| v.id == view) {
            return Err(DocumentError::Missing(view));
        }
        let filter = self.part_filter_mut(view);
        filter.hidden.remove(&part);
        filter.shown.remove(&part);
        if filter.tags_show(document, part) != shown {
            if shown {
                filter.shown.insert(part);
            } else {
                filter.hidden.insert(part);
            }
        }
        self.prune_part_filters();
        Ok(())
    }
    /// Replace a view's include and exclude tags.
    pub fn set_tag_filter(
        &mut self,
        view: EntityId,
        include: BTreeSet<String>,
        exclude: BTreeSet<String>,
    ) -> Result<(), DocumentError> {
        if !self.views.iter().any(|v| v.id == view) {
            return Err(DocumentError::Missing(view));
        }
        checked_tags(&include)?;
        checked_tags(&exclude)?;
        let filter = self.part_filter_mut(view);
        filter.include_tags = include;
        filter.exclude_tags = exclude;
        self.prune_part_filters();
        Ok(())
    }
}

pub fn has_part_filters(book: &super::DrawingBook) -> bool {
    book.pages.iter().any(|p| !p.part_filters.is_empty())
}

pub(super) fn validate_part_filters(page: &Page) -> Result<(), DocumentError> {
    let mut seen = BTreeSet::new();
    for filter in &page.part_filters {
        if !page.views.iter().any(|v| v.id == filter.view) {
            return Err(invalid("A part filter must belong to a view on its page"));
        }
        if !seen.insert(filter.view) {
            return Err(invalid("A view has one part filter"));
        }
        checked_tags(&filter.include_tags)?;
        checked_tags(&filter.exclude_tags)?;
        if filter.is_empty()
            || filter.hidden.len() + filter.shown.len() > 4096
            || filter
                .hidden
                .iter()
                .chain(&filter.shown)
                .any(|id| id.is_nil())
            || !filter.hidden.is_disjoint(&filter.shown)
        {
            return Err(invalid("A part filter needs tags or distinct part toggles"));
        }
    }
    Ok(())
}

/// The model Parts a view's sources contain, in source order without
/// repeats: a Part itself, a placed Part's definition, and the Parts of
/// assemblies and groups.
pub fn view_parts(document: &Document, view: &DrawingView) -> Vec<EntityId> {
    let ViewSource::Model { entities } = &view.source else {
        return vec![];
    };
    let hierarchy = crate::hierarchy::Index::new(document).ok();
    let mut parts = vec![];
    let mut seen = BTreeSet::new();
    let mut visited = BTreeSet::new();
    let mut pending: Vec<EntityId> = entities.iter().rev().copied().collect();
    while let Some(id) = pending.pop() {
        if !visited.insert(id) {
            continue;
        }
        match document.get(id).map(|e| &e.data) {
            Some(EntityData::Part { .. }) => {
                if seen.insert(id) {
                    parts.push(id);
                }
            }
            Some(EntityData::Instance { component, .. }) => pending.push(*component),
            Some(EntityData::Assembly { instances }) => {
                pending.extend(instances.iter().rev().copied());
            }
            Some(EntityData::Group { .. }) => {
                if let Some(children) = hierarchy.as_ref().and_then(|h| h.children(Some(id)).ok()) {
                    pending.extend(children.into_iter().rev());
                }
            }
            _ => {}
        }
    }
    parts
}
