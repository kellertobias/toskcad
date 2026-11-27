//! Projected views (CAD-437): a view projected from a base view shares one
//! paper axis with it. A top view below or above a front view keeps the same
//! paper column; a side view keeps the same row. Moving the base moves its
//! projected views with it, and a projected view only slides along its free
//! axis, so model points stay aligned across the views.
use super::{DrawingView, Page, ViewKind, source::invalid};
use crate::document::{DocumentError, EntityId};
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;

/// Files with projected-view links need a reader that keeps them aligned.
pub const PROJECTION_CAPABILITY: &str = "drawings.projections.1";

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ProjectedView {
    /// The projected view.
    pub view: EntityId,
    /// The view it is projected from.
    pub base: EntityId,
}

/// The paper axis a projected view shares with its base.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ProjectionAxis {
    /// Same right axis: the views share paper X (one above the other).
    Column,
    /// Same up axis: the views share paper Y (side by side).
    Row,
}

fn same(a: [f64; 3], b: [f64; 3]) -> bool {
    a.iter().zip(b).all(|(x, y)| (x - y).abs() < 1e-9)
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// The axis `view` shares with `base`, if both are orthographic projections
/// with a common right or up direction (but not the same view direction).
pub fn projection_axis(base: &DrawingView, view: &DrawingView) -> Option<ProjectionAxis> {
    if base.kind != ViewKind::Projection || view.kind != ViewKind::Projection {
        return None;
    }
    let (u, v) = (
        same(base.frame.u, view.frame.u),
        same(base.frame.v, view.frame.v),
    );
    match (u, v) {
        (true, false) => Some(ProjectionAxis::Column),
        (false, true) => Some(ProjectionAxis::Row),
        _ => None,
    }
}

pub fn has_projections(book: &super::DrawingBook) -> bool {
    book.pages.iter().any(|p| !p.projections.is_empty())
}

pub(super) fn validate_projections(page: &Page) -> Result<(), DocumentError> {
    let mut projected = BTreeSet::new();
    let find = |id| page.views.iter().find(|v| v.id == id);
    for link in &page.projections {
        let (Some(view), Some(base)) = (find(link.view), find(link.base)) else {
            return Err(invalid(
                "A projected view and its base must be on the same page",
            ));
        };
        if !projected.insert(link.view) {
            return Err(invalid("A view is projected from one base only"));
        }
        if page.projections.iter().any(|l| l.view == link.base) {
            return Err(invalid(
                "Project views from a base view, not from another projection",
            ));
        }
        if projection_axis(base, view).is_none() {
            return Err(invalid(
                "A projected view must share the right or up direction of its base",
            ));
        }
    }
    Ok(())
}

/// Put every projected view back on its base's scale and shared paper axis.
/// Model points then land on the same paper column or row in both views.
pub fn align_projections(page: &mut Page) {
    for link in page.projections.clone() {
        let Some(base) = page.views.iter().find(|v| v.id == link.base).cloned() else {
            continue;
        };
        let Some(view) = page.views.iter_mut().find(|v| v.id == link.view) else {
            continue;
        };
        let Some(axis) = projection_axis(&base, view) else {
            continue;
        };
        view.scale = base.scale;
        let s = super::scale_factor(base.scale);
        let d = [
            view.frame.origin_mm[0] - base.frame.origin_mm[0],
            view.frame.origin_mm[1] - base.frame.origin_mm[1],
            view.frame.origin_mm[2] - base.frame.origin_mm[2],
        ];
        match axis {
            ProjectionAxis::Column => {
                view.position_mm[0] = base.position_mm[0] + dot(d, base.frame.u) * s;
            }
            ProjectionAxis::Row => {
                view.position_mm[1] = base.position_mm[1] - dot(d, base.frame.v) * s;
            }
        }
    }
}

/// Move a view on its page by `delta_mm`, keeping projections aligned. A base
/// carries its projected views; a projected view moves only along its free
/// axis. The move is limited so every moved anchor stays printable. Returns
/// the applied delta.
pub fn move_view(
    page: &mut Page,
    id: EntityId,
    delta_mm: [f64; 2],
) -> Result<[f64; 2], DocumentError> {
    if delta_mm.iter().any(|v| !v.is_finite()) {
        return Err(invalid("Drawing view move must be finite"));
    }
    let view = page
        .views
        .iter()
        .find(|v| v.id == id)
        .ok_or(DocumentError::Missing(id))?;
    let mut delta = delta_mm;
    let moved: Vec<EntityId> = match page.projections.iter().find(|l| l.view == id) {
        Some(link) => {
            let base = page
                .views
                .iter()
                .find(|v| v.id == link.base)
                .ok_or(DocumentError::Missing(link.base))?;
            match projection_axis(base, view) {
                Some(ProjectionAxis::Column) => delta[0] = 0.,
                Some(ProjectionAxis::Row) => delta[1] = 0.,
                None => {}
            }
            vec![id]
        }
        None => std::iter::once(id)
            .chain(
                page.projections
                    .iter()
                    .filter(|l| l.base == id)
                    .map(|l| l.view),
            )
            .collect(),
    };
    let [w, h] = page.settings.dimensions();
    let m = page.settings.margin_mm;
    // Keep anchors strictly inside, as page validation requires.
    let (lo, hi) = ([m, m], [w - m - 0.01, h - m - 0.01]);
    for (i, d) in delta.iter_mut().enumerate() {
        let (mut min, mut max) = (f64::NEG_INFINITY, f64::INFINITY);
        for view in page.views.iter().filter(|v| moved.contains(&v.id)) {
            min = min.max(lo[i] - view.position_mm[i]);
            max = max.min(hi[i] - view.position_mm[i]);
        }
        *d = if min <= max { d.clamp(min, max) } else { 0. };
    }
    for view in page.views.iter_mut().filter(|v| moved.contains(&v.id)) {
        view.position_mm[0] += delta[0];
        view.position_mm[1] += delta[1];
    }
    align_projections(page);
    Ok(delta)
}

/// Remove a view and any projection link that names it. Views projected from
/// it become free-standing views at their current positions.
pub fn remove_view(page: &mut Page, id: EntityId) {
    page.views.retain(|v| v.id != id);
    page.projections.retain(|l| l.view != id && l.base != id);
}
