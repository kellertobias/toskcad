//! Sketch elements as drawn shapes (CAD-398): the left-panel list, the
//! two-way selection with the canvas, and the Inspector's size and position
//! fields. Shape identity is the canonical `Sketch::shapes` record written by
//! the draw request; curves without a record (older sources, operation
//! results) are listed one by one, never grouped by geometry.
use super::*;
use multicad_core::constraint_sketch::{
    CurveKind, Dimension, Shape, ShapeKind, Sketch as SourceSketch,
    binding::DimensionInput,
    shape_edit::{self, ShapeSize},
    source::Source,
};
use multicad_solver::SolveState;
use std::collections::{BTreeMap, BTreeSet};

/// What one list entry stands for.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) enum ElementKey {
    Shape(EntityId),
    /// A curve with no recorded shape.
    Curve(EntityId),
    Text(EntityId),
    /// A standalone point used by no curve (CAD-415).
    Point(EntityId),
}
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct ListedElement {
    pub key: ElementKey,
    pub label: String,
    pub curves: Vec<EntityId>,
    /// Element type, for the Sketch tab icon (CAD-415).
    pub kind: super::sketch_tree::ElementType,
    /// The shape this element modifies, e.g. the rectangle a fillet rounds.
    /// The Sketch tab lists it as that shape's child.
    pub parent: Option<ElementKey>,
}
/// Canonical source straight from the stored definition; positions in it are
/// authored guesses, so callers use solved geometry for any displayed value.
pub(super) fn canonical(document: &Document, sketch: EntityId) -> Option<Source> {
    match &document.get(sketch)?.data {
        EntityData::Sketch { definition, .. }
            if definition.type_id == multicad_core::constraint_sketch::source::TYPE =>
        {
            Source::from_definition(definition).ok()
        }
        _ => None,
    }
}
fn curve_label(kind: &CurveKind) -> &'static str {
    match kind {
        CurveKind::Line { .. } => "Line",
        CurveKind::Arc { .. } => "Arc",
        CurveKind::Circle { .. } => "Circle",
        CurveKind::BSpline { .. } => "Spline",
    }
}
/// Every element in list order: recorded shapes, then unrecorded curves
/// (projected references named as such), then texts, then standalone points.
pub(crate) fn elements(
    document: &Document,
    sketch: EntityId,
    source: &Source,
) -> Vec<ListedElement> {
    let mut out: Vec<_> = source
        .sketch
        .shapes
        .iter()
        .map(|shape| ListedElement {
            key: ElementKey::Shape(shape.id),
            // A fillet or chamfer shows its size ("Fillet R10"); every
            // other shape keeps its name.
            label: multicad_core::constraint_sketch::shape_edit::corner_label(
                &source.sketch,
                shape,
            ),
            curves: shape.curves.clone(),
            kind: super::sketch_tree::ElementType::of_shape(shape.kind),
            parent: super::sketch_tree::modified_shape(shape).map(ElementKey::Shape),
        })
        .collect();
    let owned: BTreeSet<_> = source
        .sketch
        .shapes
        .iter()
        .flat_map(|s| s.curves.iter().copied())
        .collect();
    let mut counts = BTreeMap::<&str, usize>::new();
    for curve in source
        .sketch
        .curves
        .iter()
        // A fillet's or chamfer's construction helpers stay in the sketch
        // but are not listed (CAD-415).
        .filter(|c| !owned.contains(&c.id) && !source.sketch.is_corner_helper(c.id))
    {
        let projected =
            crate::external_references::Controller::protected(document, sketch, curve.id);
        let kind = if projected {
            "Projected reference"
        } else {
            curve_label(&curve.kind)
        };
        let n = counts.entry(kind).or_default();
        *n += 1;
        out.push(ListedElement {
            key: ElementKey::Curve(curve.id),
            label: if projected {
                format!("{kind} {n}")
            } else {
                format!("{kind} {n} (single curve)")
            },
            curves: vec![curve.id],
            kind: if projected {
                super::sketch_tree::ElementType::Projected
            } else {
                super::sketch_tree::ElementType::of_curve(&curve.kind)
            },
            parent: None,
        });
    }
    for (index, text) in source.sketch.texts.iter().enumerate() {
        let mut content: String = text.style.text.lines().next().unwrap_or("").into();
        if content.chars().count() > 24 {
            content = content.chars().take(23).collect::<String>() + "…";
        }
        out.push(ListedElement {
            key: ElementKey::Text(text.id),
            label: format!("Text {} \u{201c}{content}\u{201d}", index + 1),
            curves: text
                .contours()
                .map(|c| c.into_iter().map(|c| c.id).collect())
                .unwrap_or_default(),
            kind: super::sketch_tree::ElementType::Text,
            parent: None,
        });
    }
    let used: BTreeSet<_> = source
        .sketch
        .curves
        .iter()
        .flat_map(|curve| match &curve.kind {
            CurveKind::Line { start, end } => vec![*start, *end],
            CurveKind::Arc {
                center, start, end, ..
            } => vec![*center, *start, *end],
            CurveKind::Circle { center, .. } => vec![*center],
            CurveKind::BSpline { poles, .. } => poles.clone(),
        })
        .collect();
    for (index, point) in source
        .sketch
        .points
        .iter()
        .filter(|p| !used.contains(&p.id))
        .enumerate()
    {
        out.push(ListedElement {
            key: ElementKey::Point(point.id),
            label: format!("Point {}", index + 1),
            curves: vec![],
            kind: super::sketch_tree::ElementType::Point,
            parent: None,
        });
    }
    out
}
fn dimension_input(text: &str) -> Result<DimensionInput, String> {
    let text = text.trim();
    if text.starts_with('=') {
        return Ok(DimensionInput::Expression {
            text: text.into(),
            unit: "mm".into(),
        });
    }
    let value: f64 = text
        .parse()
        .map_err(|_| "Enter a number in mm or an expression beginning with =")?;
    if !value.is_finite() {
        return Err("Values must be finite".into());
    }
    Ok(DimensionInput::Literal {
        value,
        unit: "mm".into(),
    })
}
fn length(a: [f64; 2], b: [f64; 2]) -> f64 {
    (b[0] - a[0]).portable_hypot(b[1] - a[1])
}
/// Current solved value of a size, from the solved sketch.
pub(super) fn measured(solved: &SourceSketch, shape: &Shape, size: ShapeSize) -> Option<f64> {
    let curve = shape_edit::size_curve(shape, size)?;
    let at = |id| {
        solved
            .points
            .iter()
            .find(|p| p.id == id)
            .map(|p| p.position_mm)
    };
    match &solved.curves.iter().find(|c| c.id == curve)?.kind {
        CurveKind::Line { start, end } => Some(length(at(*start)?, at(*end)?)),
        CurveKind::Circle { radius_mm, .. } if size == ShapeSize::Diameter => Some(2. * radius_mm),
        CurveKind::Arc { center, start, .. } => Some(length(at(*center)?, at(*start)?)),
        _ => None,
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) enum ShapeField {
    Size(u8),
    X,
    Y,
}
/// Field state for the selected shape, rebuilt whenever the selection or the
/// document revision changes.
#[derive(Default)]
pub(crate) struct ShapeFields {
    scope: Option<(ElementKey, u64)>,
    text: BTreeMap<ShapeField, String>,
    shown: BTreeMap<ShapeField, String>,
    /// Sizes that already have a driving dimension.
    driven: BTreeSet<ShapeField>,
    pub error: Option<String>,
}
pub(crate) fn field_id(field: ShapeField) -> egui::Id {
    egui::Id::new(("sketch-shape-field", field))
}
fn format_mm(value: f64) -> String {
    let rounded = (value * 1e4).round() / 1e4;
    if rounded == rounded.trunc() {
        format!("{rounded:.0}")
    } else {
        format!("{rounded}")
    }
}

impl SketchEditor {
    /// Select a list entry: its curves highlight in the canvas.
    pub(crate) fn select_element(&mut self, key: Option<ElementKey>) {
        self.shape = key;
        self.shape_fields.error = None;
        if key.is_some() {
            self.element = None;
            self.constraints.clear_geometry_selection();
            self.profiles.clear();
        }
        if let Some(ElementKey::Point(point)) = key {
            self.constraints.select_geometry(point, None, false);
        }
    }
    /// Canvas → list: the element that owns `curve`, if any.
    pub(crate) fn select_element_of_curve(
        &mut self,
        document: &Document,
        sketch: EntityId,
        curve: Option<EntityId>,
    ) {
        self.select_element_of_hit(document, sketch, curve, None);
    }
    /// Canvas → tree: the element that owns the clicked curve, or the
    /// standalone point that was clicked.
    pub(crate) fn select_element_of_hit(
        &mut self,
        document: &Document,
        sketch: EntityId,
        curve: Option<EntityId>,
        point: Option<EntityId>,
    ) {
        if let Some(point) = point {
            self.shape = canonical(document, sketch).and_then(|source| {
                elements(document, sketch, &source)
                    .into_iter()
                    .find(|e| e.key == ElementKey::Point(point))
                    .map(|e| e.key)
            });
            self.shape_fields.error = None;
            return;
        }
        self.shape = curve.and_then(|curve| {
            let source = canonical(document, sketch)?;
            // A hidden fillet or chamfer helper selects the treatment that
            // owns it (CAD-415).
            if source.sketch.is_corner_helper(curve) {
                return source
                    .sketch
                    .shapes
                    .iter()
                    .find(|s| s.helpers.contains(&curve))
                    .map(|s| ElementKey::Shape(s.id));
            }
            elements(document, sketch, &source)
                .into_iter()
                .find(|e| e.curves.contains(&curve))
                .map(|e| e.key)
        });
        self.shape_fields.error = None;
    }
    /// Curves of the selected element, for canvas highlighting.
    pub(crate) fn element_highlight(&self, document: &Document) -> BTreeSet<EntityId> {
        let (Some(key), Some(sketch)) = (self.shape, self.selected) else {
            return BTreeSet::new();
        };
        let Some(source) = canonical(document, sketch) else {
            return BTreeSet::new();
        };
        elements(document, sketch, &source)
            .into_iter()
            .find(|e| e.key == key)
            .map(|e| e.curves.into_iter().collect())
            .unwrap_or_default()
    }
    fn solved_sketch(&self, sketch: EntityId, revision: u64) -> Option<SourceSketch> {
        let presentation = self.sketch_evaluation.current(sketch)?;
        (presentation.evaluated.stamp.revision == revision)
            .then(|| presentation.evaluated.solved().cloned())
            .flatten()
    }
    /// Inspector: sizes and position of the selected shape. Each committed
    /// field is one canonical source transaction.
    pub(crate) fn shape_inspector(
        &mut self,
        ui: &mut egui::Ui,
        session: &mut Session,
        sketch: EntityId,
    ) {
        let Some(key) = self.shape else {
            return;
        };
        let Some(source) = canonical(session.document(), sketch) else {
            return;
        };
        let Some(element) = elements(session.document(), sketch, &source)
            .into_iter()
            .find(|e| e.key == key)
        else {
            return;
        };
        ui.strong(format!("Selected: {}", element.label));
        let ElementKey::Shape(shape_id) = key else {
            ui.small(match key {
                ElementKey::Text(_) => "Edit text in the Text tool.",
                ElementKey::Point(_) => "Drag the point in the canvas, or constrain or fix it below.",
                _ => "This curve was not drawn as a shape; select it in the canvas to dimension or constrain it.",
            });
            ui.separator();
            return;
        };
        let shape = source
            .sketch
            .shapes
            .iter()
            .find(|s| s.id == shape_id)
            .expect("listed shape")
            .clone();
        let revision = session.revision();
        let solved = self.solved_sketch(sketch, revision);
        if self.shape_fields.scope != Some((key, revision)) {
            let error = self.shape_fields.error.take();
            self.shape_fields = ShapeFields {
                scope: Some((key, revision)),
                error,
                ..Default::default()
            };
        }
        let mut shown = BTreeMap::new();
        self.shape_fields.driven.clear();
        for (index, size) in shape_edit::sizes(&shape).iter().enumerate() {
            let constraint = shape_edit::size_constraint(&source.sketch, &shape, *size);
            if constraint.is_some() {
                self.shape_fields
                    .driven
                    .insert(ShapeField::Size(index as u8));
            }
            let text = match constraint.and_then(|c| c.kind.dimension()) {
                Some(Dimension::Expression { expression }) => expression.text.clone(),
                _ => solved
                    .as_ref()
                    .and_then(|s| measured(s, &shape, *size))
                    .map(format_mm)
                    .unwrap_or_default(),
            };
            shown.insert(ShapeField::Size(index as u8), text);
        }
        let anchor = shape_edit::anchor(&source.sketch, &shape);
        let anchor_at = anchor.and_then(|a| {
            solved
                .as_ref()?
                .points
                .iter()
                .find(|p| p.id == a)
                .map(|p| p.position_mm)
        });
        if anchor.is_some_and(|a| {
            source
                .sketch
                .constraints
                .iter()
                .any(|c| matches!(c.kind, multicad_core::constraint_sketch::ConstraintKind::Fix { point } if point == a))
        }) {
            self.shape_fields.driven.extend([ShapeField::X, ShapeField::Y]);
        }
        if let Some(at) = anchor_at {
            shown.insert(ShapeField::X, format_mm(at[0]));
            shown.insert(ShapeField::Y, format_mm(at[1]));
        }
        for (field, value) in &shown {
            if self.shape_fields.shown.get(field) != Some(value) {
                self.shape_fields.text.insert(*field, value.clone());
            }
        }
        self.shape_fields.shown = shown.clone();
        let ready = solved.is_some() && !self.has_pending_edits();
        let mut commit = None;
        let mut row = |ui: &mut egui::Ui,
                       fields: &mut ShapeFields,
                       field: ShapeField,
                       label: &str,
                       hover: &str| {
            ui.horizontal(|ui| {
                ui.label(label);
                let text = fields.text.entry(field).or_default();
                let response = ui.add_enabled(
                    ready,
                    egui::TextEdit::singleline(text)
                        .id(field_id(field))
                        .desired_width(90.),
                );
                ui.label("mm");
                let response = response.on_hover_text(hover);
                if response.lost_focus() {
                    let enter = ui.input(|i| i.key_pressed(egui::Key::Enter));
                    if ui.input(|i| i.key_pressed(egui::Key::Escape)) {
                        *text = fields.shown.get(&field).cloned().unwrap_or_default();
                    } else if fields.shown.get(&field) != Some(text)
                        || (enter && !fields.driven.contains(&field))
                    {
                        // Enter on an unchanged free size dimensions it at its
                        // current value.
                        commit = Some((field, text.clone()));
                    }
                }
            });
        };
        let sizes = shape_edit::sizes(&shape);
        if !sizes.is_empty() {
            ui.label(egui::RichText::new("Size").small().weak());
        }
        for (index, size) in sizes.iter().enumerate() {
            row(
                ui,
                &mut self.shape_fields,
                ShapeField::Size(index as u8),
                size.label(),
                "Sets a dimension on this shape only; Enter applies as one undoable edit. =name uses a parameter",
            );
        }
        if anchor_at.is_some() {
            ui.label(
                egui::RichText::new(match shape.kind {
                    ShapeKind::Circle | ShapeKind::Arc => "Center",
                    _ => "First corner",
                })
                .small()
                .weak(),
            );
            row(
                ui,
                &mut self.shape_fields,
                ShapeField::X,
                "X",
                "Moves the whole shape and fixes this point",
            );
            row(
                ui,
                &mut self.shape_fields,
                ShapeField::Y,
                "Y",
                "Moves the whole shape and fixes this point",
            );
        }
        if sizes.is_empty() && anchor_at.is_none() {
            ui.small("This shape has no editable sizes.");
        }
        if let Some(presentation) = self.sketch_evaluation.current(sketch) {
            let diagnostics = &presentation.evaluated.evaluation.diagnostics;
            let mine: BTreeSet<_> = shape_edit::points(&source.sketch, &shape)
                .into_iter()
                .chain(shape.curves.iter().copied())
                .collect();
            let free = diagnostics
                .free_entities
                .iter()
                .filter(|e| mine.contains(e))
                .count();
            if diagnostics.state == SolveState::FullyConstrained || free == 0 {
                ui.small("Fully defined");
            } else {
                ui.colored_label(
                    Color32::LIGHT_YELLOW,
                    format!(
                        "Under-defined · {free} free element(s); add sizes or a position to fix it"
                    ),
                );
            }
        } else if !ready {
            ui.small("Solving the sketch…");
        }
        if let Some(error) = &self.shape_fields.error {
            ui.colored_label(Color32::LIGHT_RED, error);
        }
        ui.separator();
        if let Some((field, text)) = commit {
            let result = self.edit_shape(session, sketch, &shape, field, &text, solved.as_ref());
            match result {
                Ok(message) => {
                    self.shape_fields.error = None;
                    self.message = message;
                }
                Err(error) => {
                    if let Some(value) = self.shape_fields.shown.get(&field).cloned() {
                        self.shape_fields.text.insert(field, value);
                    }
                    self.message = error.clone();
                    self.shape_fields.error = Some(error);
                }
            }
        }
    }
    /// Apply one field as a canonical edit, or explain why nothing changed.
    pub(crate) fn edit_shape(
        &mut self,
        session: &mut Session,
        sketch: EntityId,
        shape: &Shape,
        field: ShapeField,
        text: &str,
        solved: Option<&SourceSketch>,
    ) -> Result<String, String> {
        if self.has_pending_edits() {
            return Err(self.pending_edits_refusal());
        }
        let source = canonical(session.document(), sketch).ok_or("Select a constrained sketch")?;
        let (draft, what) = match field {
            ShapeField::Size(index) => {
                let size = *shape_edit::sizes(shape)
                    .get(index as usize)
                    .ok_or("This shape has no such size")?;
                let input = dimension_input(text)?;
                if let DimensionInput::Literal { value, .. } = input
                    && value <= 0.
                {
                    return Err(format!("{} must be greater than 0 mm", size.label()));
                }
                (
                    shape_edit::set_size(session.document(), &source, shape.id, size, &input)
                        .map_err(|e| e.to_string())?,
                    format!("{} {}", shape.name, size.label().to_lowercase()),
                )
            }
            ShapeField::X | ShapeField::Y => {
                let solved = solved.ok_or("Wait for the current sketch solution")?;
                let value: f64 = text
                    .trim()
                    .parse()
                    .map_err(|_| "Enter the position as a number in mm")?;
                let anchor = shape_edit::anchor(&source.sketch, shape)
                    .ok_or("This shape has no position")?;
                let mut at = solved
                    .points
                    .iter()
                    .find(|p| p.id == anchor)
                    .map(|p| p.position_mm)
                    .ok_or("Wait for the current sketch solution")?;
                at[usize::from(field == ShapeField::Y)] = value;
                (
                    shape_edit::set_position(&source, solved, shape.id, at)
                        .map_err(|e| e.to_string())?,
                    format!("{} position", shape.name),
                )
            }
        };
        let edit = multicad_solver::prepare_canonical_edit(
            session,
            session.revision(),
            Actor::local(Channel::Ui),
            sketch,
            draft,
            self.evaluated.as_ref(),
        )
        .map_err(|e| e.to_string())?;
        let diagnostics = edit.evaluation.diagnostics.clone();
        // A move must not resize the shape: an attached shape can pull the
        // corners of an undimensioned one while the solver balances them.
        if matches!(field, ShapeField::X | ShapeField::Y)
            && let (Some(before), Some(after)) =
                (solved, edit.semantic.as_ref().and_then(|s| s.solved()))
            && let Some(size) = shape_edit::sizes(shape).iter().find(|size| {
                match (
                    measured(before, shape, **size),
                    measured(after, shape, **size),
                ) {
                    (Some(old), Some(new)) => (old - new).abs() > 1e-6,
                    _ => false,
                }
            })
        {
            return Err(format!(
                "Moving {} would change its {} because attached geometry pulls on it · set its sizes first · nothing changed",
                shape.name,
                size.label().to_lowercase()
            ));
        }
        match (edit.preview, diagnostics.state) {
            (Some(preview), SolveState::FullyConstrained | SolveState::Underconstrained) => {
                session.commit(preview).map_err(|e| e.to_string())?;
                Ok(if diagnostics.state == SolveState::FullyConstrained {
                    format!("Set {what} · sketch fully defined · one undoable edit")
                } else {
                    format!(
                        "Set {what} · sketch under-defined, {} degree(s) of freedom remain · one undoable edit",
                        diagnostics.degrees_of_freedom
                    )
                })
            }
            (_, SolveState::Conflict) => Err(format!(
                "Over-constrained: {what} conflicts with {} existing constraint(s) · nothing changed",
                diagnostics.conflicts.len().max(1)
            )),
            (_, state) => Err(format!(
                "{what} could not be solved ({}) · nothing changed",
                match state {
                    SolveState::NotConverged => "the solver did not converge",
                    SolveState::Flipped => "the shape would flip",
                    _ => "no valid solution",
                }
            )),
        }
    }
}
