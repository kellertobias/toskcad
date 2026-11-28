//! Sketch-mode left sidebar (CAD-415): Model and Sketch tabs. The Sketch tab
//! is an icon tree of the open sketch's elements, in the model tree's row
//! style, built on the canonical `Sketch::shapes` record (CAD-398). Elements
//! that modify a shape (fillets, chamfers) are listed as its children;
//! constraints and dimensions sit in a collapsible group.
use super::shape_list::{ElementKey, ListedElement, canonical, elements};
use super::*;
use multicad_core::constraint_sketch::{CurveKind, Shape, ShapeKind};
use std::collections::BTreeSet;

/// Which sidebar tab is shown while a sketch is open.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub(crate) enum SidebarTab {
    Model,
    #[default]
    Sketch,
}
impl SidebarTab {
    pub(crate) fn label(self) -> &'static str {
        match self {
            Self::Model => "Model",
            Self::Sketch => "Sketch",
        }
    }
}
pub(crate) fn tab_id(tab: SidebarTab) -> egui::Id {
    egui::Id::new(("sketch-sidebar-tab", tab))
}

/// Element type shown by the Sketch tab icon.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) enum ElementType {
    Rectangle,
    Circle,
    Arc,
    Line,
    Polyline,
    Spline,
    #[allow(
        dead_code,
        reason = "Shown once the core shape record has slots (CAD-412)"
    )]
    Slot,
    Text,
    Point,
    Projected,
    Fillet,
    Chamfer,
}
impl ElementType {
    #[cfg(test)]
    pub(crate) const ALL: [Self; 12] = [
        Self::Rectangle,
        Self::Circle,
        Self::Arc,
        Self::Line,
        Self::Polyline,
        Self::Spline,
        Self::Slot,
        Self::Text,
        Self::Point,
        Self::Projected,
        Self::Fillet,
        Self::Chamfer,
    ];
    pub(crate) fn of_shape(kind: ShapeKind) -> Self {
        match kind {
            ShapeKind::Line => Self::Line,
            ShapeKind::Rectangle => Self::Rectangle,
            ShapeKind::Circle => Self::Circle,
            ShapeKind::Arc => Self::Arc,
            ShapeKind::Polyline => Self::Polyline,
            ShapeKind::Spline => Self::Spline,
            ShapeKind::Fillet => Self::Fillet,
            ShapeKind::Chamfer => Self::Chamfer,
        }
    }
    pub(crate) fn of_curve(kind: &CurveKind) -> Self {
        match kind {
            CurveKind::Line { .. } => Self::Line,
            CurveKind::Arc { .. } => Self::Arc,
            CurveKind::Circle { .. } => Self::Circle,
            CurveKind::BSpline { .. } => Self::Spline,
        }
    }
    pub(crate) fn icon(self) -> Icon {
        match self {
            Self::Rectangle => Icon::Rectangle,
            Self::Circle => Icon::Circle,
            Self::Arc => Icon::Arc,
            Self::Line => Icon::Line,
            Self::Polyline => Icon::Polyline,
            Self::Spline => Icon::Spline,
            Self::Slot => Icon::Slot,
            Self::Text => Icon::SketchText,
            Self::Point => Icon::SketchPoint,
            Self::Projected => Icon::ProjectedReference,
            Self::Fillet => Icon::SketchFillet,
            Self::Chamfer => Icon::SketchChamfer,
        }
    }
    pub(crate) fn label(self) -> &'static str {
        match self {
            Self::Rectangle => "Rectangle",
            Self::Circle => "Circle",
            Self::Arc => "Arc",
            Self::Line => "Line",
            Self::Polyline => "Polyline",
            Self::Spline => "Spline",
            Self::Slot => "Slot",
            Self::Text => "Text",
            Self::Point => "Point",
            Self::Projected => "Projected reference",
            Self::Fillet => "Fillet",
            Self::Chamfer => "Chamfer",
        }
    }
}
/// The shape a recorded fillet or chamfer modifies, read from the canonical
/// record (CAD-414); `None` keeps the shape at the top level.
pub(crate) fn modified_shape(shape: &Shape) -> Option<EntityId> {
    shape.modifies
}

/// One painted tree row.
#[derive(Clone, Debug, PartialEq)]
pub(crate) enum TreeEntry {
    Element {
        key: ElementKey,
        kind: ElementType,
        label: String,
        depth: usize,
        children: bool,
    },
    ConstraintGroup {
        count: usize,
    },
    Constraint {
        id: EntityId,
        dimension: bool,
        label: String,
    },
}
/// Flatten `listed` into rows: each element is followed by the elements
/// that modify it, one level deeper, unless it is collapsed. An element
/// whose parent is missing stays top-level.
pub(crate) fn element_rows(
    listed: &[ListedElement],
    collapsed: &BTreeSet<ElementKey>,
) -> Vec<TreeEntry> {
    let present: BTreeSet<_> = listed.iter().map(|e| e.key).collect();
    let parent = |e: &ListedElement| e.parent.filter(|p| present.contains(p) && *p != e.key);
    let mut out = vec![];
    fn push(
        out: &mut Vec<TreeEntry>,
        listed: &[ListedElement],
        parent: &dyn Fn(&ListedElement) -> Option<ElementKey>,
        collapsed: &BTreeSet<ElementKey>,
        element: &ListedElement,
        depth: usize,
        seen: &mut BTreeSet<ElementKey>,
    ) {
        if !seen.insert(element.key) {
            return;
        }
        let children: Vec<_> = listed
            .iter()
            .filter(|c| parent(c) == Some(element.key))
            .collect();
        out.push(TreeEntry::Element {
            key: element.key,
            kind: element.kind,
            label: element.label.clone(),
            depth,
            children: !children.is_empty(),
        });
        if !collapsed.contains(&element.key) {
            for child in children {
                push(out, listed, parent, collapsed, child, depth + 1, seen);
            }
        }
    }
    let mut seen = BTreeSet::new();
    for element in listed.iter().filter(|e| parent(e).is_none()) {
        push(&mut out, listed, &parent, collapsed, element, 0, &mut seen);
    }
    out
}

/// Sidebar state kept by the editor.
#[derive(Default)]
pub(crate) struct SketchSidebar {
    pub tab: SidebarTab,
    /// Sketch the tab choice belongs to; opening another resets to Sketch.
    open: Option<EntityId>,
    collapsed: BTreeSet<ElementKey>,
    constraints_open: bool,
    /// Last selection scrolled into view, so canvas picks reveal their row.
    revealed: Option<ElementKey>,
    /// Rows painted in the last frame (tests read icons and depth here).
    pub rows: Vec<TreeEntry>,
}
impl SketchSidebar {
    /// Follow the open sketch: entering one selects the Sketch tab.
    pub(crate) fn follow(&mut self, open: Option<EntityId>) {
        if self.open != open {
            self.open = open;
            self.tab = SidebarTab::Sketch;
            self.collapsed.clear();
            self.revealed = None;
            self.rows.clear();
        }
    }
}

fn tab_button(ui: &mut egui::Ui, tab: SidebarTab, selected: bool) -> egui::Response {
    let label = tab.label();
    let galley = egui::WidgetText::from(label).into_galley(
        ui,
        Some(egui::TextWrapMode::Extend),
        f32::INFINITY,
        egui::TextStyle::Button,
    );
    let (rect, _) =
        ui.allocate_exact_size(egui::vec2(galley.size().x + 20., 24.), egui::Sense::hover());
    let response = ui.interact(rect, tab_id(tab), egui::Sense::click());
    response.widget_info(|| {
        egui::WidgetInfo::selected(egui::WidgetType::SelectableLabel, true, selected, label)
    });
    let visuals = ui.visuals();
    if selected {
        ui.painter().rect_filled(rect, 3, visuals.selection.bg_fill);
    } else if response.hovered() {
        ui.painter()
            .rect_filled(rect, 3, visuals.widgets.hovered.weak_bg_fill);
    }
    let color = if selected {
        visuals.strong_text_color()
    } else {
        visuals.text_color()
    };
    ui.painter()
        .galley(rect.center() - galley.size() * 0.5, galley, color);
    response.on_hover_text(match tab {
        SidebarTab::Model => "Show the model tree",
        SidebarTab::Sketch => "Show the elements of the open sketch",
    })
}

impl SketchEditor {
    /// The sketch open for editing, if any.
    fn open_sketch(&self) -> Option<EntityId> {
        self.selected.filter(|_| self.editing)
    }
    /// Model / Sketch tab bar. Returns the tab to show; without an open
    /// sketch there is no bar and the model tree shows.
    pub(crate) fn sidebar_tabs(&mut self, ui: &mut egui::Ui) -> SidebarTab {
        let open = self.open_sketch();
        self.sketch_sidebar.follow(open);
        if open.is_none() {
            return SidebarTab::Model;
        }
        ui.horizontal(|ui| {
            for tab in [SidebarTab::Model, SidebarTab::Sketch] {
                let selected = self.sketch_sidebar.tab == tab;
                if tab_button(ui, tab, selected).clicked() {
                    self.sketch_sidebar.tab = tab;
                }
            }
        });
        ui.separator();
        self.sketch_sidebar.tab
    }
    /// Sketch tab: the open sketch's elements as an icon tree.
    pub(crate) fn sketch_tree(&mut self, ui: &mut egui::Ui, session: &Session) {
        let Some(sketch) = self.open_sketch() else {
            return;
        };
        let Some(source) = canonical(session.document(), sketch) else {
            ui.small("This sketch has no editable elements");
            return;
        };
        let listed = elements(session.document(), sketch, &source);
        if self
            .shape
            .is_some_and(|key| !listed.iter().any(|e| e.key == key))
        {
            self.shape = None;
        }
        // A canvas pick of a collapsed child reveals it.
        if let Some(parent) = self
            .shape
            .and_then(|key| listed.iter().find(|e| e.key == key)?.parent)
        {
            self.sketch_sidebar.collapsed.remove(&parent);
        }
        let mut rows = element_rows(&listed, &self.sketch_sidebar.collapsed);
        let constraints: Vec<_> = source
            .sketch
            .constraints
            .iter()
            .map(|c| (c.id, c.kind.dimension().is_some(), c.kind.clone()))
            .collect();
        if !constraints.is_empty() {
            rows.push(TreeEntry::ConstraintGroup {
                count: constraints.len(),
            });
            if self.sketch_sidebar.constraints_open {
                let mut counts = std::collections::BTreeMap::<&str, usize>::new();
                for (id, dimension, kind) in &constraints {
                    let name = crate::sketch_constraints::label(kind);
                    let n = counts.entry(name).or_default();
                    *n += 1;
                    rows.push(TreeEntry::Constraint {
                        id: *id,
                        dimension: *dimension,
                        label: format!("{name} {n}"),
                    });
                }
            }
        }
        ui.strong("Sketch elements");
        if listed.is_empty() {
            ui.small("Nothing drawn yet");
        }
        let reveal = (self.shape != self.sketch_sidebar.revealed)
            .then_some(self.shape)
            .flatten();
        self.sketch_sidebar.revealed = self.shape;
        let mut clicked = None;
        let mut toggle = None;
        let mut toggle_constraints = false;
        egui::ScrollArea::vertical()
            .id_salt("sketch-element-tree")
            .auto_shrink([false, true])
            .show(ui, |ui| {
                for row in &rows {
                    let (row_rect, _) = ui.allocate_exact_size(
                        egui::vec2(ui.available_width().max(0.), 26.),
                        egui::Sense::hover(),
                    );
                    let mut row_ui = ui.new_child(
                        egui::UiBuilder::new()
                            .max_rect(row_rect)
                            .layout(egui::Layout::left_to_right(egui::Align::Center)),
                    );
                    row_ui.set_clip_rect(row_rect.intersect(ui.clip_rect()));
                    let ui = &mut row_ui;
                    let trailing = 4.;
                    match row {
                        TreeEntry::Element {
                            key,
                            kind,
                            label,
                            depth,
                            children,
                        } => {
                            ui.push_id(("sketch-tree-row", *key), |ui| {
                                ui.add_space((*depth as f32 * 10.).min(80.));
                                if *children {
                                    let expanded =
                                        !self.sketch_sidebar.collapsed.contains(key);
                                    if crate::model_tree::disclosure_button(ui, expanded, label)
                                        .clicked()
                                    {
                                        toggle = Some(*key);
                                    }
                                } else {
                                    ui.add_space(16.);
                                }
                                let icon = icons::model_tree_icon(
                                    ui,
                                    kind.icon(),
                                    &format!("Select {} {label}", kind.label()),
                                    false,
                                );
                                let width = (ui.available_width() - trailing).max(0.);
                                let selected = self.shape == Some(*key);
                                let color = ui.visuals().text_color();
                                let response = crate::model_tree::flat_row_label(
                                    ui,
                                    width,
                                    label,
                                    true,
                                    selected,
                                    color,
                                    false,
                                );
                                if reveal == Some(*key) {
                                    response.scroll_to_me(None);
                                }
                                let response = response.on_hover_text(format!(
                                    "{label}\n{}\nSelect it in the sketch and edit it in the Inspector",
                                    kind.label()
                                ));
                                if response.clicked() || icon.clicked() {
                                    clicked = Some((!selected).then_some(*key));
                                }
                            });
                        }
                        TreeEntry::ConstraintGroup { count } => {
                            let label = format!("Constraints and dimensions ({count})");
                            if crate::model_tree::disclosure_button(
                                ui,
                                self.sketch_sidebar.constraints_open,
                                &label,
                            )
                            .clicked()
                            {
                                toggle_constraints = true;
                            }
                            let icon = icons::model_tree_icon(
                                ui,
                                Icon::SketchConstraint,
                                &label,
                                false,
                            );
                            let width = (ui.available_width() - trailing).max(0.);
                                let color = ui.visuals().text_color();
                            let response = crate::model_tree::flat_row_label(
                                ui,
                                width,
                                &label,
                                true,
                                false,
                                color,
                                false,
                            );
                            if response.clicked() || icon.clicked() {
                                toggle_constraints = true;
                            }
                        }
                        TreeEntry::Constraint {
                            id,
                            dimension,
                            label,
                        } => {
                            ui.push_id(("sketch-tree-constraint", *id), |ui| {
                                ui.add_space(26.);
                                icons::model_tree_icon(
                                    ui,
                                    if *dimension {
                                        Icon::SketchDimension
                                    } else {
                                        Icon::SketchConstraint
                                    },
                                    label,
                                    false,
                                );
                                let width = (ui.available_width() - trailing).max(0.);
                                let color = ui.visuals().text_color();
                                crate::model_tree::flat_row_label(
                                    ui,
                                    width,
                                    label,
                                    false,
                                    false,
                                    color,
                                    false,
                                )
                                .on_hover_text(
                                    "Select its geometry in the canvas to edit or remove it",
                                );
                            });
                        }
                    }
                }
            });
        self.sketch_sidebar.rows = rows;
        if let Some(key) = toggle
            && !self.sketch_sidebar.collapsed.remove(&key)
        {
            self.sketch_sidebar.collapsed.insert(key);
        }
        if toggle_constraints {
            self.sketch_sidebar.constraints_open = !self.sketch_sidebar.constraints_open;
        }
        if let Some(key) = clicked {
            self.select_element(key);
            self.tool = None;
            self.clear_drawing_progress();
        }
    }
}
