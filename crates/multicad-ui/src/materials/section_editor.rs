//! Modal section-sketch editor for "Based on sketch" profile stock (CAD-385).
//!
//! The product owner chose a modal containing the sketch elements, opened
//! from Project → Materials (an explicit exception to the rule that sketch
//! editing is contextual to Model). The modal edits a local copy of the
//! stock's owned constraint sketch with the sketch editor's command
//! pipeline: the same typed draw command (`sketch_draw`), the same
//! constraint solver for dimensions, and the same exact region extraction
//! that validates extrusion profiles. Nothing reaches the document from the
//! modal: Finish returns the solved sketch to the stock draft, and the
//! Materials Apply commits it in one transaction. Cancel discards the copy.
use crate::{egui, modeling::dimension::DimensionField};
use eframe::egui::{Color32, Pos2, Rect, Sense, Stroke, Vec2};
use multicad_core::{
    constraint_sketch::{Constraint, ConstraintKind, CurveKind, Dimension, Sketch},
    document::EntityId,
    materials::{ProfileShape, StockKind},
    model::{SectionEdge, SectionLoop, SketchSection, resolve_sketch_section},
    sketch_draw::{DrawInput, DrawRequest, PointInput},
};
use multicad_precision::PortableMath as _;
use multicad_solver::profile_section::{SolvedSection, draw_section, solve_section_sketch};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Tool {
    Select,
    Line,
    Rectangle,
    Circle,
    Arc,
}
impl Tool {
    pub(crate) const ALL: [Self; 5] = [
        Self::Select,
        Self::Line,
        Self::Rectangle,
        Self::Circle,
        Self::Arc,
    ];
    pub(crate) fn label(self) -> &'static str {
        match self {
            Self::Select => "Select",
            Self::Line => "Line",
            Self::Rectangle => "Rectangle",
            Self::Circle => "Circle",
            Self::Arc => "Arc",
        }
    }
    fn hint(self) -> &'static str {
        match self {
            Self::Select => {
                "Click an element to select it; Delete removes it, Dimension adds a dimension."
            }
            Self::Line => {
                "Click points to draw connected lines; click the first point to close. Escape ends the chain."
            }
            Self::Rectangle => "Click two opposite corners.",
            Self::Circle => "Click the center, then a point on the circle.",
            Self::Arc => "Click the center, the start and the end (counterclockwise).",
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) enum Outcome {
    Finished(Box<SolvedSection>),
    Cancelled,
}

/// One pending tool click: a snapped position, or an existing sketch point.
#[derive(Clone, Copy, Debug, PartialEq)]
struct Click {
    at: [f64; 2],
    point: Option<EntityId>,
}
impl Click {
    fn input(self) -> PointInput {
        self.point
            .map_or(PointInput::Position(self.at), PointInput::Existing)
    }
}

pub(crate) struct SectionEditor {
    sketch: Sketch,
    undo: Vec<Sketch>,
    redo: Vec<Sketch>,
    tool: Tool,
    clicks: Vec<Click>,
    selected: Option<EntityId>,
    dimensions: Vec<(EntityId, DimensionField)>,
    /// Section-plane center and pixels per millimeter of the canvas.
    view: Option<([f64; 2], f64)>,
    message: Option<String>,
    touched: bool,
}

const CANVAS: Vec2 = Vec2::new(480., 340.);
const SNAP_PX: f32 = 8.;

impl SectionEditor {
    pub(crate) fn new(sketch: Sketch) -> Self {
        let mut editor = Self {
            sketch,
            undo: vec![],
            redo: vec![],
            tool: Tool::Select,
            clicks: vec![],
            selected: None,
            dimensions: vec![],
            view: None,
            message: None,
            touched: false,
        };
        editor.sync_dimensions();
        if editor.sketch.curves.is_empty() {
            editor.tool = Tool::Rectangle;
        }
        editor
    }
    /// A first sketch for stock switching to "Based on sketch": the current
    /// preset outline, dimensioned, so the author edits rather than redraws.
    pub(crate) fn seeded(specification: &StockKind) -> Self {
        let StockKind::Profile {
            width_mm,
            height_mm,
            wall_mm,
            shape,
            ..
        } = specification
        else {
            return Self::new(Sketch::default());
        };
        let (w, h) = (*width_mm, *height_mm);
        let mut editor = Self::new(Sketch::default());
        let ok = |w: f64| w.is_finite() && w > 0.;
        if !(ok(w) && ok(h)) {
            return editor;
        }
        let mut seed = |input| {
            if let Ok(next) = draw_section(
                &editor.sketch,
                &DrawRequest {
                    input,
                    construction: false,
                },
            ) {
                editor.sketch = next;
            }
        };
        let p = PointInput::Position;
        match shape {
            ProfileShape::Rectangular => {
                seed(DrawInput::Rectangle {
                    a: p([0., 0.]),
                    b: p([w, h]),
                });
                if let Some(t) = wall_mm.filter(|t| 2. * t < w.min(h)) {
                    seed(DrawInput::Rectangle {
                        a: p([t, t]),
                        b: p([w - t, h - t]),
                    });
                }
            }
            ProfileShape::Round => {
                seed(DrawInput::Circle {
                    center: p([w / 2., w / 2.]),
                    radius_mm: w / 2.,
                });
                if let Some(t) = wall_mm.filter(|t| 2. * t < w) {
                    seed(DrawInput::Circle {
                        center: p([w / 2., w / 2.]),
                        radius_mm: w / 2. - t,
                    });
                }
            }
            ProfileShape::Angle => {
                let t = wall_mm.unwrap_or(2.).min(w.min(h) / 2.);
                seed(DrawInput::Polyline {
                    vertices: [[0., 0.], [w, 0.], [w, t], [t, t], [t, h], [0., h]]
                        .into_iter()
                        .map(p)
                        .collect(),
                    closed: true,
                });
            }
            ProfileShape::ItemLine8_40x40 | ProfileShape::Sketch => {}
        }
        // Dimension each outside edge length and circle diameter once.
        let curves: Vec<_> = editor.sketch.curves.iter().map(|c| c.id).collect();
        for id in curves.into_iter().take(2) {
            editor.selected = Some(id);
            editor.add_dimension();
        }
        editor.selected = None;
        editor.undo.clear();
        editor.touched = false;
        editor.message = None;
        editor.sync_dimensions();
        editor
    }
    #[cfg(test)]
    pub(crate) fn sketch(&self) -> &Sketch {
        &self.sketch
    }
    #[cfg(test)]
    pub(crate) fn touched(&self) -> bool {
        self.touched
    }
    #[cfg(test)]
    pub(crate) fn tool(&self) -> Tool {
        self.tool
    }
    #[cfg(test)]
    pub(crate) fn message(&self) -> Option<&str> {
        self.message.as_deref()
    }
    #[cfg(test)]
    pub(crate) fn dimension_fields(&self) -> &[(EntityId, DimensionField)] {
        &self.dimensions
    }
    pub(crate) fn status(&self) -> Result<SketchSection, String> {
        resolve_sketch_section(&self.sketch).map_err(|e| e.to_string())
    }
    fn commit(&mut self, next: Sketch) {
        self.undo.push(std::mem::replace(&mut self.sketch, next));
        self.redo.clear();
        self.touched = true;
        self.message = None;
        self.sync_dimensions();
    }
    pub(crate) fn undo(&mut self) {
        if let Some(previous) = self.undo.pop() {
            self.redo
                .push(std::mem::replace(&mut self.sketch, previous));
            self.selected = None;
            self.sync_dimensions();
        }
    }
    pub(crate) fn redo(&mut self) {
        if let Some(next) = self.redo.pop() {
            self.undo.push(std::mem::replace(&mut self.sketch, next));
            self.selected = None;
            self.sync_dimensions();
        }
    }
    fn sync_dimensions(&mut self) {
        let previous = std::mem::take(&mut self.dimensions);
        self.dimensions = self
            .sketch
            .constraints
            .iter()
            .filter_map(|c| {
                let dimension = c.kind.dimension()?;
                let unit = if matches!(c.kind, ConstraintKind::Angle { .. }) {
                    "deg"
                } else {
                    "mm"
                };
                let field = previous
                    .iter()
                    .find(|(id, f)| *id == c.id && f.resolve_literal() == literal(dimension))
                    .map(|(_, f)| f.clone())
                    .or_else(|| DimensionField::from_dimension(dimension, unit).ok())?;
                Some((c.id, field))
            })
            .collect();
    }
    /// Run a drawing command; refused input leaves the sketch unchanged.
    pub(crate) fn draw(&mut self, input: DrawInput) {
        match draw_section(
            &self.sketch,
            &DrawRequest {
                input,
                construction: false,
            },
        ) {
            Ok(next) => self.commit(next),
            Err(error) => self.message = Some(error),
        }
    }
    #[cfg(test)]
    pub(crate) fn select(&mut self, curve: Option<EntityId>) {
        self.selected = curve;
    }
    /// Remove the selected element, its constraints and orphaned points.
    pub(crate) fn delete_selected(&mut self) {
        let Some(id) = self.selected.take() else {
            return;
        };
        let mut next = self.sketch.clone();
        next.curves.retain(|c| c.id != id);
        next.constraints
            .retain(|c| !c.kind.references().contains(&id));
        let used: Vec<EntityId> = next.curves.iter().flat_map(curve_points).collect();
        let removed: Vec<EntityId> = next
            .points
            .iter()
            .map(|p| p.id)
            .filter(|p| !used.contains(p))
            .collect();
        next.points.retain(|p| !removed.contains(&p.id));
        next.constraints
            .retain(|c| !c.kind.references().iter().any(|r| removed.contains(r)));
        match solve_section_sketch(&next) {
            Ok(solved) => self.commit(solved),
            Err(error) => self.message = Some(error),
        }
    }
    /// Dimension the selected element at its current size: a line's length,
    /// a circle's diameter or an arc's radius.
    pub(crate) fn add_dimension(&mut self) {
        let Some(id) = self.selected else {
            self.message = Some("Select a line, circle or arc to dimension".into());
            return;
        };
        let Some(curve) = self.sketch.curves.iter().find(|c| c.id == id) else {
            return;
        };
        let mm = |value: f64| Dimension::Literal {
            value: round(value),
            unit: "mm".into(),
        };
        let kind = match curve.kind {
            CurveKind::Line { start, end } => ConstraintKind::Length {
                line: id,
                dimension: mm(distance(self.at(start), self.at(end))),
            },
            CurveKind::Circle { radius_mm, .. } => ConstraintKind::Diameter {
                curve: id,
                dimension: mm(2. * radius_mm),
            },
            CurveKind::Arc { center, start, .. } => ConstraintKind::Radius {
                curve: id,
                dimension: mm(distance(self.at(center), self.at(start))),
            },
            CurveKind::BSpline { .. } => {
                self.message = Some("Splines cannot be dimensioned in a section".into());
                return;
            }
        };
        if self
            .sketch
            .constraints
            .iter()
            .any(|c| c.kind.dimension().is_some() && c.kind.references() == kind.references())
        {
            self.message = Some("This element is already dimensioned".into());
            return;
        }
        let mut next = self.sketch.clone();
        next.constraints.push(Constraint {
            id: EntityId::new(),
            reference: false,
            kind,
        });
        match solve_section_sketch(&next) {
            Ok(solved) => self.commit(solved),
            Err(error) => self.message = Some(error),
        }
    }
    /// Change one dimension to a number and re-solve. Expressions are
    /// refused: shared stock dimensions are literal.
    pub(crate) fn set_dimension(&mut self, constraint: EntityId, text: &str) {
        let Ok(value) = text.trim().parse::<f64>() else {
            self.message = Some(
                "Section dimensions must be numbers; project-parameter expressions are not supported in shared stock"
                    .into(),
            );
            self.sync_dimensions();
            return;
        };
        if !value.is_finite() || value <= 0. {
            self.message = Some("Enter a dimension above zero".into());
            self.sync_dimensions();
            return;
        }
        let mut next = self.sketch.clone();
        let Some(dimension) = next
            .constraints
            .iter_mut()
            .find(|c| c.id == constraint)
            .and_then(|c| c.kind.dimension_mut())
        else {
            return;
        };
        if literal(dimension) == Some(value) {
            return;
        }
        let unit = dimension.unit().to_string();
        *dimension = Dimension::Literal { value, unit };
        match solve_section_sketch(&next) {
            Ok(solved) => self.commit(solved),
            Err(error) => {
                self.message = Some(error);
                self.sync_dimensions();
            }
        }
    }
    pub(crate) fn remove_dimension(&mut self, constraint: EntityId) {
        let mut next = self.sketch.clone();
        next.constraints.retain(|c| c.id != constraint);
        match solve_section_sketch(&next) {
            Ok(solved) => self.commit(solved),
            Err(error) => self.message = Some(error),
        }
    }
    pub(crate) fn set_tool(&mut self, tool: Tool) {
        self.tool = tool;
        self.clicks.clear();
    }
    fn at(&self, point: EntityId) -> [f64; 2] {
        self.sketch
            .points
            .iter()
            .find(|p| p.id == point)
            .map_or([0.; 2], |p| p.position_mm)
    }
    /// Apply one canvas click in section millimeters.
    pub(crate) fn click(&mut self, click: [f64; 2], snap_mm: f64) {
        let point = self
            .sketch
            .points
            .iter()
            .filter(|p| {
                self.sketch
                    .curves
                    .iter()
                    .any(|c| curve_points(c).contains(&p.id))
            })
            .map(|p| (p.id, p.position_mm, distance(p.position_mm, click)))
            .filter(|(_, _, d)| *d <= snap_mm)
            .min_by(|a, b| a.2.total_cmp(&b.2));
        let click = match point {
            Some((id, at, _)) => Click {
                at,
                point: Some(id),
            },
            None => Click {
                at: click,
                point: None,
            },
        };
        match self.tool {
            Tool::Select => {
                self.selected = self.nearest_curve(click.at, snap_mm);
            }
            Tool::Line => {
                if let Some(start) = self.clicks.last().copied() {
                    if start.at == click.at {
                        return;
                    }
                    self.draw(DrawInput::Line {
                        start: start.input(),
                        end: click.input(),
                    });
                    // Continue the chain from the end point unless it closed.
                    let closed = click.point.is_some()
                        && self.clicks.first().is_some_and(|f| f.at == click.at);
                    self.clicks.clear();
                    if !closed && self.message.is_none() {
                        let end = self
                            .sketch
                            .points
                            .iter()
                            .find(|p| p.position_mm == click.at)
                            .map(|p| p.id);
                        self.clicks.push(Click {
                            at: click.at,
                            point: end,
                        });
                    }
                } else {
                    self.clicks.push(click);
                }
            }
            Tool::Rectangle => {
                if let Some(a) = self.clicks.pop() {
                    self.draw(DrawInput::Rectangle {
                        a: a.input(),
                        b: click.input(),
                    });
                } else {
                    self.clicks.push(click);
                }
            }
            Tool::Circle => {
                if let Some(center) = self.clicks.pop() {
                    let radius_mm = round(distance(center.at, click.at));
                    if radius_mm > 0. {
                        self.draw(DrawInput::Circle {
                            center: center.input(),
                            radius_mm,
                        });
                    }
                } else {
                    self.clicks.push(click);
                }
            }
            Tool::Arc => {
                self.clicks.push(click);
                if let [center, start, end] = self.clicks[..] {
                    self.clicks.clear();
                    // Place the end on the start radius so the arc is exact.
                    let r = distance(center.at, start.at);
                    let d = distance(center.at, end.at);
                    let end = if end.point.is_some() || d == 0. {
                        end
                    } else {
                        Click {
                            at: [
                                center.at[0] + (end.at[0] - center.at[0]) * r / d,
                                center.at[1] + (end.at[1] - center.at[1]) * r / d,
                            ],
                            point: None,
                        }
                    };
                    self.draw(DrawInput::Arc {
                        center: center.input(),
                        start: start.input(),
                        end: end.input(),
                        clockwise: false,
                    });
                }
            }
        }
    }
    fn nearest_curve(&self, at: [f64; 2], within: f64) -> Option<EntityId> {
        self.sketch
            .curves
            .iter()
            .filter_map(|c| {
                let d = match c.kind {
                    CurveKind::Line { start, end } => {
                        segment_distance(at, self.at(start), self.at(end))
                    }
                    CurveKind::Circle { center, radius_mm } => {
                        (distance(at, self.at(center)) - radius_mm).abs()
                    }
                    CurveKind::Arc { center, start, .. } => (distance(at, self.at(center))
                        - distance(self.at(start), self.at(center)))
                    .abs(),
                    CurveKind::BSpline { .. } => return None,
                };
                (d <= within).then_some((c.id, d))
            })
            .min_by(|a, b| a.1.total_cmp(&b.1))
            .map(|(id, _)| id)
    }

    /// Paint the modal. Returns the outcome once the author finishes or
    /// cancels; `None` while editing continues.
    pub(crate) fn show(&mut self, ctx: &egui::Context) -> Option<Outcome> {
        let mut outcome = None;
        let width = (ctx.content_rect().width() - 48.).clamp(300., 820.);
        let modal = egui::Modal::new(egui::Id::new("profile-section-sketch"))
            .backdrop_color(Color32::from_black_alpha(80))
            .show(ctx, |ui| {
                ui.set_width(width);
                ui.heading("Profile section sketch");
                ui.small("Draw one closed outline and optional holes. Dimensions are numbers; Finish returns the section to this stock.");
                ui.add_space(6.);
                ui.horizontal_wrapped(|ui| {
                    for tool in Tool::ALL {
                        if ui
                            .selectable_label(self.tool == tool, tool.label())
                            .on_hover_text(tool.hint())
                            .clicked()
                        {
                            self.set_tool(tool);
                        }
                    }
                    ui.separator();
                    if ui.add_enabled(!self.undo.is_empty(), egui::Button::new("Undo")).clicked() {
                        self.undo();
                    }
                    if ui.add_enabled(!self.redo.is_empty(), egui::Button::new("Redo")).clicked() {
                        self.redo();
                    }
                    if ui
                        .add_enabled(self.selected.is_some(), egui::Button::new("Delete"))
                        .clicked()
                    {
                        self.delete_selected();
                    }
                    if ui
                        .add_enabled(self.selected.is_some(), egui::Button::new("Dimension"))
                        .on_hover_text("Dimension the selected line, circle or arc at its current size")
                        .clicked()
                    {
                        self.add_dimension();
                    }
                    if ui.button("Fit").clicked() {
                        self.view = None;
                    }
                });
                ui.small(self.tool.hint());
                let narrow = width < 620.;
                if narrow {
                    self.canvas(ui, Vec2::new(width, CANVAS.y * width / CANVAS.x));
                    self.dimension_list(ui);
                } else {
                    ui.horizontal_top(|ui| {
                        self.canvas(ui, Vec2::new(width - 240., CANVAS.y));
                        ui.vertical(|ui| {
                            ui.set_width(228.);
                            self.dimension_list(ui);
                        });
                    });
                }
                let status = self.status();
                match &status {
                    Ok(section) => {
                        ui.colored_label(
                            ui.visuals().text_color(),
                            format!(
                                "Valid section · {} × {} mm · area {:.1} mm²{}",
                                trim(section.size_mm[0]),
                                trim(section.size_mm[1]),
                                section.area_mm2,
                                if section.holes.is_empty() {
                                    String::new()
                                } else {
                                    format!(" · {} hole(s)", section.holes.len())
                                }
                            ),
                        );
                    }
                    Err(error) => {
                        ui.colored_label(ui.visuals().error_fg_color, error);
                    }
                }
                if let Some(message) = &self.message {
                    ui.colored_label(ui.visuals().warn_fg_color, message);
                }
                ui.add_space(6.);
                ui.horizontal(|ui| {
                    if ui.button("Cancel").clicked() {
                        outcome = Some(Outcome::Cancelled);
                    }
                    let finish = ui
                        .add_enabled(status.is_ok(), egui::Button::new("Finish"))
                        .on_disabled_hover_text("Finish needs one closed outline with optional holes");
                    if finish.clicked()
                        && let Ok(section) = status.clone()
                    {
                        outcome = Some(Outcome::Finished(Box::new(SolvedSection {
                            sketch: self.sketch.clone(),
                            section,
                        })));
                    }
                });
            });
        let _ = modal;
        if outcome.is_none()
            && ctx.input_mut(|i| i.consume_key(egui::Modifiers::NONE, egui::Key::Escape))
        {
            if self.clicks.is_empty() {
                outcome = Some(Outcome::Cancelled);
            } else {
                self.clicks.clear();
            }
        }
        if outcome.is_none()
            && self.selected.is_some()
            && ctx.memory(|m| m.focused().is_none())
            && ctx.input_mut(|i| {
                i.consume_key(egui::Modifiers::NONE, egui::Key::Delete)
                    || i.consume_key(egui::Modifiers::NONE, egui::Key::Backspace)
            })
        {
            self.delete_selected();
        }
        outcome
    }

    fn dimension_list(&mut self, ui: &mut egui::Ui) {
        ui.strong("Dimensions");
        if self.dimensions.is_empty() {
            ui.small("Select an element and press Dimension to give it an exact size.");
        }
        let mut apply = None;
        let mut remove = None;
        let selected = self.selected;
        let constraints = self.sketch.constraints.clone();
        for (index, (id, field)) in self.dimensions.iter_mut().enumerate() {
            let Some(constraint) = constraints.iter().find(|c| c.id == *id) else {
                continue;
            };
            let label = match constraint.kind {
                ConstraintKind::Length { .. } => "Length",
                ConstraintKind::Diameter { .. } => "Diameter",
                ConstraintKind::Radius { .. } => "Radius",
                ConstraintKind::Angle { .. } => "Angle",
                _ => "Distance",
            };
            let highlighted = selected.is_some_and(|s| constraint.kind.references().contains(&s));
            ui.vertical(|ui| {
                let response = field.show_linked(
                    ui,
                    ("section-dimension", *id),
                    &format!("{label} {}", index + 1),
                    "A number; Enter re-solves the sketch",
                    highlighted,
                    false,
                );
                if response.lost_focus() {
                    apply = Some((*id, field.text.clone()));
                }
                if ui
                    .small_button("✕")
                    .on_hover_text("Remove this dimension")
                    .clicked()
                {
                    remove = Some(*id);
                }
            });
        }
        if let Some((id, text)) = apply {
            self.set_dimension(id, &text);
        }
        if let Some(id) = remove {
            self.remove_dimension(id);
        }
    }

    fn fit(&self, size: Vec2) -> ([f64; 2], f64) {
        let mut min = [f64::INFINITY; 2];
        let mut max = [f64::NEG_INFINITY; 2];
        for c in &self.sketch.curves {
            let mut add = |p: [f64; 2], r: f64| {
                for i in 0..2 {
                    min[i] = min[i].min(p[i] - r);
                    max[i] = max[i].max(p[i] + r);
                }
            };
            match c.kind {
                CurveKind::Circle { center, radius_mm } => add(self.at(center), radius_mm),
                CurveKind::Arc { center, start, .. } => {
                    add(self.at(center), distance(self.at(center), self.at(start)))
                }
                _ => curve_points(c)
                    .into_iter()
                    .for_each(|p| add(self.at(p), 0.)),
            }
        }
        if !min[0].is_finite() {
            min = [-10., -10.];
            max = [90., 70.];
        }
        let span = [(max[0] - min[0]).max(1.), (max[1] - min[1]).max(1.)];
        let scale = (f64::from(size.x) / (span[0] * 1.3)).min(f64::from(size.y) / (span[1] * 1.3));
        ([(min[0] + max[0]) / 2., (min[1] + max[1]) / 2.], scale)
    }

    fn canvas(&mut self, ui: &mut egui::Ui, size: Vec2) {
        let (rect, response) = ui.allocate_exact_size(size, Sense::click());
        if self.view.is_none() {
            self.view = Some(self.fit(size));
        }
        let (center, scale) = self.view.expect("fitted view");
        let to_screen = |p: [f64; 2]| {
            Pos2::new(
                rect.center().x + ((p[0] - center[0]) * scale) as f32,
                rect.center().y - ((p[1] - center[1]) * scale) as f32,
            )
        };
        let to_mm = |p: Pos2| {
            [
                center[0] + f64::from(p.x - rect.center().x) / scale,
                center[1] - f64::from(p.y - rect.center().y) / scale,
            ]
        };
        let painter = ui.painter_at(rect);
        let v = ui.visuals();
        painter.rect_filled(rect, 4., v.extreme_bg_color);
        // Grid: 10 mm major lines with the origin axes emphasized.
        let step = grid_step(scale);
        let [x0, y1] = to_mm(rect.min);
        let [x1, y0] = to_mm(rect.max);
        let mut x = (x0 / step).floor() * step;
        while x <= x1 {
            let p = to_screen([x, 0.]).x;
            painter.line_segment(
                [Pos2::new(p, rect.top()), Pos2::new(p, rect.bottom())],
                Stroke::new(
                    if x == 0. { 1.2 } else { 0.4 },
                    v.weak_text_color().gamma_multiply(0.5),
                ),
            );
            x += step;
        }
        let mut y = (y0 / step).floor() * step;
        while y <= y1 {
            let p = to_screen([0., y]).y;
            painter.line_segment(
                [Pos2::new(rect.left(), p), Pos2::new(rect.right(), p)],
                Stroke::new(
                    if y == 0. { 1.2 } else { 0.4 },
                    v.weak_text_color().gamma_multiply(0.5),
                ),
            );
            y += step;
        }
        // A valid section's boundaries are traced in the accent color so the
        // author sees exactly which loops form the material and its holes.
        if let Ok(section) = self.status() {
            let accent = Stroke::new(4., v.selection.bg_fill.gamma_multiply(0.6));
            let offset = section.offset_mm;
            for boundary in std::iter::once(&section.outer).chain(&section.holes) {
                let points: Vec<Pos2> = loop_points(boundary)
                    .into_iter()
                    .map(|p| to_screen([p[0] + offset[0], p[1] + offset[1]]))
                    .collect();
                painter.add(egui::Shape::closed_line(points, accent));
            }
        }
        for curve in self.sketch.curves.iter() {
            let selected = self.selected == Some(curve.id);
            let stroke = Stroke::new(
                if selected { 2.5 } else { 1.6 },
                if selected {
                    v.selection.stroke.color
                } else if curve.construction {
                    v.weak_text_color()
                } else {
                    v.strong_text_color()
                },
            );
            match curve.kind {
                CurveKind::Line { start, end } => {
                    painter
                        .line_segment([to_screen(self.at(start)), to_screen(self.at(end))], stroke);
                }
                CurveKind::Circle { center, radius_mm } => {
                    painter.circle_stroke(
                        to_screen(self.at(center)),
                        (radius_mm * scale) as f32,
                        stroke,
                    );
                }
                CurveKind::Arc {
                    center,
                    start,
                    end,
                    clockwise,
                } => {
                    let edge = SectionEdge::Arc {
                        center: self.at(center),
                        a: self.at(start),
                        b: self.at(end),
                        clockwise,
                    };
                    let points = loop_points(&SectionLoop::Curved(vec![edge]));
                    let mut screen: Vec<Pos2> = points.into_iter().map(to_screen).collect();
                    screen.push(to_screen(self.at(end)));
                    painter.add(egui::Shape::line(screen, stroke));
                }
                CurveKind::BSpline { .. } => {}
            }
        }
        for point in &self.sketch.points {
            painter.circle_filled(to_screen(point.position_mm), 2.5, v.strong_text_color());
        }
        // Rubber band for the pending tool.
        let hover = response.hover_pos().map(|p| snap(to_mm(p), step / 10.));
        if let (Some(first), Some(hover)) = (self.clicks.first(), hover) {
            let guide = Stroke::new(1., v.selection.stroke.color);
            let a = to_screen(first.at);
            let b = to_screen(hover);
            match self.tool {
                Tool::Line => {
                    painter.line_segment([to_screen(self.clicks.last().unwrap().at), b], guide);
                }
                Tool::Rectangle => {
                    painter.rect_stroke(
                        Rect::from_two_pos(a, b),
                        0.,
                        guide,
                        egui::StrokeKind::Middle,
                    );
                }
                Tool::Circle | Tool::Arc => {
                    painter.circle_stroke(a, a.distance(b), guide);
                }
                Tool::Select => {}
            }
        }
        if response.clicked()
            && let Some(p) = response.interact_pointer_pos()
        {
            self.click(snap(to_mm(p), step / 10.), f64::from(SNAP_PX) / scale);
        }
        if response.hovered() {
            let zoom = ui.input(|i| i.zoom_delta());
            if zoom != 1.
                && let Some(view) = &mut self.view
            {
                view.1 *= f64::from(zoom);
            }
        }
    }
}

fn literal(d: &Dimension) -> Option<f64> {
    match d {
        Dimension::Literal { value, .. } => Some(*value),
        Dimension::Expression { .. } => None,
    }
}
trait LiteralField {
    fn resolve_literal(&self) -> Option<f64>;
}
impl LiteralField for DimensionField {
    fn resolve_literal(&self) -> Option<f64> {
        self.text.trim().parse().ok()
    }
}
fn curve_points(curve: &multicad_core::constraint_sketch::Curve) -> Vec<EntityId> {
    match &curve.kind {
        CurveKind::Line { start, end } => vec![*start, *end],
        CurveKind::Arc {
            center, start, end, ..
        } => vec![*center, *start, *end],
        CurveKind::Circle { center, .. } => vec![*center],
        CurveKind::BSpline { poles, .. } => poles.clone(),
    }
}
fn distance(a: [f64; 2], b: [f64; 2]) -> f64 {
    (a[0] - b[0]).portable_hypot(a[1] - b[1])
}
fn segment_distance(p: [f64; 2], a: [f64; 2], b: [f64; 2]) -> f64 {
    let d = [b[0] - a[0], b[1] - a[1]];
    let len = d[0] * d[0] + d[1] * d[1];
    let t = if len == 0. {
        0.
    } else {
        (((p[0] - a[0]) * d[0] + (p[1] - a[1]) * d[1]) / len).clamp(0., 1.)
    };
    distance(p, [a[0] + d[0] * t, a[1] + d[1] * t])
}
/// Six significant decimals: typed and clicked sizes stay exact.
fn round(value: f64) -> f64 {
    (value * 1e6).round() / 1e6
}
fn trim(value: f64) -> String {
    let text = format!("{value:.3}");
    text.trim_end_matches('0').trim_end_matches('.').to_string()
}
fn grid_step(scale: f64) -> f64 {
    // About 30 px between grid lines: 1, 2, 5 × 10ⁿ millimeters.
    let target = 30. / scale;
    let mut step = 1e-3;
    loop {
        for m in [1., 2., 5.] {
            if step * m >= target {
                return step * m;
            }
        }
        step *= 10.;
    }
}
fn snap(p: [f64; 2], step: f64) -> [f64; 2] {
    p.map(|v| round((v / step).round() * step))
}
/// Boundary points for display: vertices and points along arcs.
fn loop_points(section_loop: &SectionLoop) -> Vec<[f64; 2]> {
    let arc = |center: [f64; 2], a: [f64; 2], b: [f64; 2], clockwise: bool| {
        let r = distance(center, a);
        let start = (a[1] - center[1]).portable_atan2(a[0] - center[0]);
        let end = (b[1] - center[1]).portable_atan2(b[0] - center[0]);
        let mut sweep = (end - start).rem_euclid(std::f64::consts::TAU);
        if sweep == 0. {
            sweep = std::f64::consts::TAU;
        }
        if clockwise {
            sweep -= std::f64::consts::TAU;
        }
        (0..32)
            .map(|i| {
                let t = start + sweep * f64::from(i) / 32.;
                [
                    center[0] + r * t.portable_cos(),
                    center[1] + r * t.portable_sin(),
                ]
            })
            .collect::<Vec<_>>()
    };
    match section_loop {
        SectionLoop::Polygon(points) => points.clone(),
        SectionLoop::Circle { center, radius_mm } => {
            let a = [center[0] + radius_mm, center[1]];
            arc(*center, a, a, false)
        }
        SectionLoop::Curved(edges) => edges
            .iter()
            .flat_map(|e| match *e {
                SectionEdge::Line { a, .. } => vec![a],
                SectionEdge::Arc {
                    center,
                    a,
                    b,
                    clockwise,
                } => arc(center, a, b, clockwise),
            })
            .collect(),
    }
}

#[cfg(test)]
mod tests;
