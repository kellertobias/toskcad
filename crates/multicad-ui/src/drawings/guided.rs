//! Guided view placement (CAD-436). Add view starts a front view of the
//! selection, or of the whole model, at the largest standard scale that fits
//! the sheet; it follows the pointer until a click places it. Standard views
//! arranges front, top, right and iso views in one click, with top and right
//! projected from the front. On a saved page the result is applied as one
//! undoable edit as soon as its native preview is ready.
use super::*;
use multicad_core::{
    document::{Document, EntityData},
    drawings::{DrawingView, ProjectedView, Scale, ViewEvaluation, ViewKind, ViewSource},
    references::ReferenceFrame,
};

/// One-click orientations offered while adding a view.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StandardView {
    Front,
    Top,
    Right,
    Iso,
}
impl StandardView {
    pub const ALL: [Self; 4] = [Self::Front, Self::Top, Self::Right, Self::Iso];
    pub fn label(self) -> &'static str {
        match self {
            Self::Front => "Front",
            Self::Top => "Top",
            Self::Right => "Right",
            Self::Iso => "Iso",
        }
    }
    /// The view's right (u) and up (v) directions in world coordinates.
    pub fn frame(self) -> ([f64; 3], [f64; 3]) {
        match self {
            Self::Front => ([1., 0., 0.], [0., 0., 1.]),
            Self::Top => ([1., 0., 0.], [0., 1., 0.]),
            Self::Right => ([0., 1., 0.], [0., 0., 1.]),
            Self::Iso => (
                [2f64.sqrt().recip(), -2f64.sqrt().recip(), 0.],
                [6f64.sqrt().recip(), 6f64.sqrt().recip(), 2. / 6f64.sqrt()],
            ),
        }
    }
}

#[derive(Clone, Debug)]
pub(super) struct ViewPlacement {
    /// The pending views in the draft: one, or front, top, right and iso.
    views: Vec<EntityId>,
    standard: bool,
    orientation: StandardView,
    /// The scale the user chose; otherwise the fitted standard scale.
    scale: Option<Scale>,
    /// Scale and frame are fitted to the current preview of these views.
    measured: bool,
    /// The clicked anchor of a single view.
    position: Option<[f64; 2]>,
    /// The page before placement, restored on cancel.
    before: Page,
    /// The page was saved when placement began: apply when ready.
    apply: bool,
    /// What the views show, for the hint.
    target: String,
    /// Half the fitted single view's paper size, to keep it on the paper.
    half_mm: Option<[f64; 2]>,
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// The model's top-level sources: parts and assemblies not used as a
/// component of an instance, and instances not owned by an assembly.
pub fn whole_model(document: &Document) -> Vec<EntityId> {
    let mut used = std::collections::BTreeSet::new();
    for entity in document.entities() {
        match &entity.data {
            EntityData::Instance { component, .. } => {
                used.insert(*component);
            }
            EntityData::Assembly { instances } => used.extend(instances.iter().copied()),
            _ => {}
        }
    }
    document
        .entities()
        .filter(|e| {
            matches!(
                e.data,
                EntityData::Part { .. } | EntityData::Instance { .. } | EntityData::Assembly { .. }
            ) && !used.contains(&e.id)
        })
        .map(|e| e.id)
        .collect()
}

impl DrawingEditor {
    /// The model selection a new view shows instead of the whole model.
    pub fn set_model_selection(&mut self, selection: Vec<EntityId>) {
        self.model_selection = selection;
    }
    /// Where a paper point was drawn in the last frame, for automation.
    pub fn sheet_point(&self, mm: [f64; 2]) -> Option<egui::Pos2> {
        let (paper, scale) = self.last_paper?;
        Some(paper.min + egui::vec2(mm[0] as f32, mm[1] as f32) * scale)
    }
    pub fn view_placement_active(&self) -> bool {
        self.view_placement.is_some()
    }
    fn view_sources(&self, session: &Session) -> Result<(Vec<EntityId>, String), String> {
        let document = session.document();
        let selected: Vec<_> = self
            .model_selection
            .iter()
            .copied()
            .filter(|id| {
                document.get(*id).is_some_and(|e| {
                    matches!(
                        e.data,
                        EntityData::Part { .. }
                            | EntityData::Instance { .. }
                            | EntityData::Assembly { .. }
                    )
                })
            })
            .collect();
        let count = |n: usize| {
            if n == 1 {
                "1 part".into()
            } else {
                format!("{n} parts")
            }
        };
        if !selected.is_empty() {
            let n = selected.len();
            return Ok((selected, format!("the selection ({})", count(n))));
        }
        let all = whole_model(document);
        if all.is_empty() {
            return Err("Add a part in Model before adding a drawing view".into());
        }
        let n = all.len();
        Ok((all, format!("the whole model ({})", count(n))))
    }
    fn new_view(
        &self,
        sources: &[EntityId],
        orientation: StandardView,
        at: [f64; 2],
    ) -> DrawingView {
        let (u, v) = orientation.frame();
        DrawingView {
            id: EntityId::new(),
            name: format!("{} view", orientation.label()),
            source: ViewSource::Model {
                entities: sources.to_vec(),
            },
            frame: ReferenceFrame {
                origin_mm: [0.; 3],
                u,
                v,
            },
            kind: ViewKind::Projection,
            position_mm: at,
            scale: self.draft.as_ref().map_or(
                Scale {
                    numerator: 1,
                    denominator: 1,
                },
                |p| p.settings.scale,
            ),
            hidden_lines: false,
            tangent_edges: false,
            paper_chord_tolerance_mm: 0.01,
            guide: None,
            equipment_categories: None,
        }
    }
    /// Start Add view (one view under the pointer) or Standard views.
    pub(super) fn begin_view_placement(
        &mut self,
        session: &Session,
        standard: bool,
    ) -> Result<(), String> {
        if self.view_placement.is_some() {
            return Err("Place or cancel the current view first".into());
        }
        if self.placing() || self.section_pick.is_some() {
            return Err("Finish or cancel the current sheet placement first".into());
        }
        let (sources, target) = self.view_sources(session)?;
        let apply = !self.dirty() && self.revision == Some(session.revision());
        let before = self.draft.clone().ok_or("Add a page first")?;
        let [w, h] = before.settings.dimensions();
        let center = [w * 0.5, h * 0.45];
        let orientations: &[StandardView] = if standard {
            &StandardView::ALL
        } else {
            &[StandardView::Front]
        };
        let views: Vec<_> = orientations
            .iter()
            .map(|o| self.new_view(&sources, *o, center))
            .collect();
        let ids: Vec<_> = views.iter().map(|v| v.id).collect();
        let page = self.draft.as_mut().expect("draft checked");
        page.views.extend(views);
        if standard {
            for projected in [ids[1], ids[2]] {
                page.projections.push(ProjectedView {
                    view: projected,
                    base: ids[0],
                });
            }
        }
        self.selected_view = None;
        self.selected_annotation = None;
        self.view_placement = Some(ViewPlacement {
            views: ids,
            standard,
            orientation: StandardView::Front,
            scale: None,
            measured: false,
            position: None,
            before,
            apply,
            target,
            half_mm: None,
        });
        self.error = None;
        Ok(())
    }
    pub fn cancel_view_placement(&mut self) -> bool {
        let Some(placement) = self.view_placement.take() else {
            return false;
        };
        self.draft = Some(placement.before);
        true
    }
    /// Re-orient the pending single view; it is fitted again.
    pub fn set_placement_orientation(&mut self, orientation: StandardView) {
        let Some(placement) = self.view_placement.as_mut() else {
            return;
        };
        if placement.standard || placement.orientation == orientation {
            return;
        }
        placement.orientation = orientation;
        placement.measured = false;
        let id = placement.views[0];
        if let Some(view) = self
            .draft
            .as_mut()
            .and_then(|p| p.views.iter_mut().find(|v| v.id == id))
        {
            let (u, v) = orientation.frame();
            view.frame = ReferenceFrame {
                origin_mm: [0.; 3],
                u,
                v,
            };
            view.name = format!("{} view", orientation.label());
        }
    }
    /// Override the fitted scale of the pending views.
    pub fn set_placement_scale(&mut self, scale: Scale) {
        let Some(placement) = self.view_placement.as_mut() else {
            return;
        };
        placement.scale = Some(scale);
        placement.measured = false;
    }
    pub fn placement_scale(&self) -> Option<Scale> {
        let placement = self.view_placement.as_ref()?;
        let page = self.draft.as_ref()?;
        page.views
            .iter()
            .find(|v| v.id == placement.views[0])
            .map(|v| v.scale)
    }
    /// Place the pending single view at this paper point.
    pub fn place_view(&mut self, at: [f64; 2]) -> Result<(), String> {
        let at = self.clamp_to_paper(at);
        let placement = self
            .view_placement
            .as_mut()
            .ok_or("No view is being placed")?;
        if placement.standard {
            return Err("Standard views place themselves".into());
        }
        placement.position = Some(at);
        let id = placement.views[0];
        if let Some(view) = self
            .draft
            .as_mut()
            .and_then(|p| p.views.iter_mut().find(|v| v.id == id))
        {
            view.position_mm = at;
        }
        Ok(())
    }
    /// Keep the anchor printable and, while it fits, the whole pending view.
    fn clamp_to_paper(&self, p: [f64; 2]) -> [f64; 2] {
        let Some(page) = &self.draft else {
            return p;
        };
        let [w, h] = page.settings.dimensions();
        let m = page.settings.margin_mm;
        let half = self
            .view_placement
            .as_ref()
            .and_then(|v| v.half_mm)
            .unwrap_or([0.; 2]);
        let (lo, hi) = ([m, m], [w - m - 0.01, h - m - 0.01]);
        [0, 1].map(|i| {
            let (a, b) = (lo[i] + half[i], hi[i] - half[i]);
            if a <= b {
                p[i].clamp(a, b)
            } else {
                ((lo[i] + hi[i]) * 0.5).clamp(lo[i], hi[i])
            }
        })
    }
    /// The draft as shown while a single view follows the pointer.
    pub(super) fn following(&self, pointer: Option<[f64; 2]>) -> Option<Page> {
        let placement = self.view_placement.as_ref()?;
        if placement.standard || placement.position.is_some() {
            return None;
        }
        let at = self.clamp_to_paper(pointer?);
        let mut page = self.draft.clone()?;
        let view = page.views.iter_mut().find(|v| v.id == placement.views[0])?;
        view.position_mm = at;
        Some(page)
    }
    /// Fit the pending views once their preview is current, then apply a
    /// saved page when the placed result is ready.
    pub(super) fn tick_view_placement(&mut self, session: &mut Session) {
        let Some(placement) = self.view_placement.clone() else {
            return;
        };
        if !placement.measured {
            if !self.preview_matches(session) || self.feed.pending() {
                return;
            }
            let Some(result) = self.feed.result.clone() else {
                return;
            };
            let mut bounds = vec![];
            for id in &placement.views {
                match result.views.get(id) {
                    Some(ViewEvaluation::Ready { vectors, .. }) => {
                        let mut min = [f64::INFINITY; 2];
                        let mut max = [f64::NEG_INFINITY; 2];
                        for p in vectors.iter().flat_map(|v| &v.preview_mm) {
                            for i in 0..2 {
                                min[i] = min[i].min(p[i]);
                                max[i] = max[i].max(p[i]);
                            }
                        }
                        if !min[0].is_finite() {
                            self.error = Some("The view shows no geometry".into());
                            return;
                        }
                        bounds.push((min, max));
                    }
                    Some(ViewEvaluation::Unresolved { message }) => {
                        self.error = Some(message.clone());
                        return;
                    }
                    None => return,
                }
            }
            if let Err(error) = self.fit_placement(&placement, &bounds) {
                self.error = Some(error);
                self.cancel_view_placement();
                return;
            }
            let half = (!placement.standard).then(|| {
                let (min, max) = bounds[0];
                let s = self
                    .placement_scale_of(placement.views[0])
                    .map_or(1., drawings::scale_factor);
                [(max[0] - min[0]) * s * 0.5, (max[1] - min[1]) * s * 0.5]
            });
            if let Some(p) = self.view_placement.as_mut() {
                p.measured = true;
                p.half_mm = half;
                if p.standard {
                    p.position = Some([0.; 2]);
                }
            }
            return;
        }
        if placement.position.is_none() {
            return;
        }
        if !placement.apply {
            // Joins the pending page draft for Apply page.
            self.view_placement = None;
            self.selected_view = placement.views.first().copied();
            return;
        }
        if !self.enabled(DrawingTool::Apply, session) {
            return;
        }
        self.view_placement = None;
        let scale = self.placement_scale_of(placement.views[0]);
        match self.invoke(session, DrawingTool::Apply, Channel::Ui) {
            Ok(()) => {
                self.selected_view = placement.views.first().copied();
                self.notice = Some(if placement.standard {
                    format!(
                        "Front, top, right and iso views placed at {} · drag a view to move it",
                        scale.map_or_else(String::new, drawings::scale_text)
                    )
                } else {
                    format!(
                        "{} view placed at {} · drag it to move it",
                        placement.orientation.label(),
                        scale.map_or_else(String::new, drawings::scale_text)
                    )
                });
            }
            Err(error) => self.error = Some(error),
        }
    }
    fn placement_scale_of(&self, id: EntityId) -> Option<Scale> {
        self.draft
            .as_ref()?
            .views
            .iter()
            .find(|v| v.id == id)
            .map(|v| v.scale)
    }
    /// Choose scale and frame origins from the measured model extents.
    fn fit_placement(
        &mut self,
        placement: &ViewPlacement,
        bounds: &[([f64; 2], [f64; 2])],
    ) -> Result<(), String> {
        let page = self.draft.as_mut().ok_or("Select a page")?;
        let layout = drawings::layout_page(page, 1, 1).map_err(|e| e.to_string())?;
        let find = |page: &mut Page, id: EntityId| -> Option<usize> {
            page.views.iter().position(|v| v.id == id)
        };
        if !placement.standard {
            let index = find(page, placement.views[0]).ok_or("View is missing")?;
            let (min, max) = bounds[0];
            let view = &mut page.views[index];
            let fit = drawings::fit_view(view, &layout, min, max, drawings::FitMode::Standard)?;
            let center = fit.center_mm;
            for i in 0..3 {
                view.frame.origin_mm[i] +=
                    view.frame.u[i] * center[0] + view.frame.v[i] * center[1];
            }
            view.scale = placement.scale.unwrap_or(fit.scale);
            return Ok(());
        }
        // Front (x, z), top (x, y): the model center all four views share.
        let mid = |i: usize, axis: usize| (bounds[i].0[axis] + bounds[i].1[axis]) * 0.5;
        let center = [mid(0, 0), mid(1, 1), mid(0, 1)];
        let extents: Vec<[f64; 2]> = bounds
            .iter()
            .map(|(min, max)| [max[0] - min[0], max[1] - min[1]])
            .collect();
        let arranged = drawings::standard_layout(
            &layout,
            page.settings.projection,
            [extents[0], extents[1], extents[2], extents[3]],
        )?;
        let scale = placement.scale.unwrap_or(arranged.scale);
        for (slot, id) in placement.views.iter().enumerate() {
            let index = find(page, *id).ok_or("View is missing")?;
            let view = &mut page.views[index];
            let (u, v) = (view.frame.u, view.frame.v);
            let mut origin = center;
            if slot == 3 {
                // The iso view centers on its own extent.
                let (du, dv) = (mid(3, 0) - dot(center, u), mid(3, 1) - dot(center, v));
                for i in 0..3 {
                    origin[i] += u[i] * du + v[i] * dv;
                }
            }
            view.frame.origin_mm = origin;
            view.scale = scale;
            view.position_mm = arranged.positions[slot];
        }
        drawings::align_projections(page);
        Ok(())
    }
    /// What to do next, shown on the sheet.
    pub fn next_step_hint(&self) -> Option<String> {
        if let Some(placement) = &self.view_placement {
            return Some(if !placement.measured {
                format!("Preparing a view of {}…", placement.target)
            } else if placement.standard {
                "Placing front, top, right and iso views…".into()
            } else if placement.position.is_none() {
                format!(
                    "Click the sheet to place the {} view of {} · Front, Top, Right or Iso changes it · Escape cancels",
                    placement.orientation.label().to_lowercase(),
                    placement.target
                )
            } else {
                "Placing the view…".into()
            });
        }
        let page = self.draft.as_ref()?;
        if self.placing() || self.section_pick.is_some() {
            return None;
        }
        Some(if page.views.is_empty() {
            "Next: Add view places a front view of the model — then click the sheet. Standard views adds front, top, right and iso at once.".into()
        } else {
            "Drag a view to move it · Add view adds another view · click a view to edit it".into()
        })
    }
    /// Inspector controls while a view is being added.
    pub(super) fn view_placement_inspector(&mut self, ui: &mut egui::Ui) {
        let Some(placement) = self.view_placement.clone() else {
            return;
        };
        ui.heading(if placement.standard {
            "Standard views"
        } else {
            "Add view"
        });
        ui.label(format!("Shows {}", placement.target));
        if !placement.standard {
            ui.label("Orientation");
            ui.horizontal_wrapped(|ui| {
                for orientation in StandardView::ALL {
                    if ui
                        .selectable_label(placement.orientation == orientation, orientation.label())
                        .on_hover_text("Show the model from this standard direction")
                        .clicked()
                    {
                        self.set_placement_orientation(orientation);
                    }
                }
            });
        }
        let current = self.placement_scale();
        ui.label(if placement.scale.is_some() {
            "Scale (chosen)"
        } else {
            "Scale (largest standard scale that fits the sheet)"
        });
        crate::controls::inspector_combo(ui, "placement-scale")
            .selected_text(current.map_or("Fitting…".into(), drawings::scale_text))
            .show_ui(ui, |ui| {
                for scale in drawings::standard_scales() {
                    if ui
                        .selectable_label(current == Some(scale), drawings::scale_text(scale))
                        .clicked()
                    {
                        self.set_placement_scale(scale);
                        ui.close();
                    }
                }
            });
        if ui
            .button("Cancel")
            .on_hover_text("Stop adding the view; the page stays as it was")
            .clicked()
        {
            self.cancel_view_placement();
        }
        ui.separator();
    }
}
