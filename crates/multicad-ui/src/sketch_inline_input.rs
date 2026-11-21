//! Creation-time numeric input beside the geometry being drawn (CAD-384).
//!
//! While Line, Polyline, Rectangle or Circle is active, editable fields float
//! in the viewport: the next point's X/Y in sketch coordinates before the
//! first click, then the size (rectangle width and height, circle diameter,
//! line length) beside the preview. Entered values fix those quantities while
//! the pointer still places the rest, and the same preview is committed by a
//! click or Enter. Typing a number, or Tab, enters the fields; unfocused fields
//! never intercept the pointer. Sizes become ordinary driving dimensions in the same single
//! draw transaction, so Undo, editing and save/reopen keep them; positions only
//! place points. Values accept a number in millimeters or the parameter
//! expression syntax used elsewhere, such as `=$width / 2`.
use super::*;
use multicad_core::constraint_sketch::{ConstraintKind, Dimension, source::Source};

/// One field of the inline input.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Field {
    X,
    Y,
    Width,
    Height,
    Diameter,
    Length,
}
impl Field {
    pub(crate) fn label(self) -> &'static str {
        match self {
            Self::X => "X",
            Self::Y => "Y",
            Self::Width => "W",
            Self::Height => "H",
            Self::Diameter => "⌀",
            Self::Length => "L",
        }
    }
    fn size(self) -> bool {
        !matches!(self, Self::X | Self::Y)
    }
}

/// The fields offered for a tool at its current drawing stage.
pub(crate) fn fields(tool: DrawTool, started: bool) -> &'static [Field] {
    match (tool, started) {
        (DrawTool::Line | DrawTool::Polyline | DrawTool::Rectangle | DrawTool::Circle, false) => {
            &[Field::X, Field::Y]
        }
        (DrawTool::Rectangle, true) => &[Field::Width, Field::Height],
        (DrawTool::Circle, true) => &[Field::Diameter],
        (DrawTool::Line | DrawTool::Polyline, true) => &[Field::Length],
        _ => &[],
    }
}

/// The pointer-placed point with every entered quantity held fixed. Sizes keep
/// the side of the first point the pointer is on.
pub(crate) fn constrain(
    tool: DrawTool,
    first: Option<[f64; 2]>,
    pointer: [f64; 2],
    values: &[Option<f64>],
) -> [f64; 2] {
    let value = |i: usize| values.get(i).copied().flatten();
    let Some(a) = first else {
        return [
            value(0).unwrap_or(pointer[0]),
            value(1).unwrap_or(pointer[1]),
        ];
    };
    let along = |length: f64| {
        let d = [pointer[0] - a[0], pointer[1] - a[1]];
        let n = d[0].hypot(d[1]);
        let d = if n > 1e-12 {
            [d[0] / n, d[1] / n]
        } else {
            [1., 0.]
        };
        [a[0] + d[0] * length, a[1] + d[1] * length]
    };
    let side = |axis: usize| if pointer[axis] < a[axis] { -1. } else { 1. };
    match tool {
        DrawTool::Rectangle => [
            value(0).map_or(pointer[0], |w| a[0] + side(0) * w),
            value(1).map_or(pointer[1], |h| a[1] + side(1) * h),
        ],
        DrawTool::Circle => value(0).map_or(pointer, |d| along(d / 2.)),
        DrawTool::Line | DrawTool::Polyline => value(0).map_or(pointer, along),
        _ => pointer,
    }
}

/// The current measured value of each field for the pointer-placed preview.
fn measured(tool: DrawTool, first: Option<[f64; 2]>, p: [f64; 2]) -> Vec<f64> {
    let Some(a) = first else {
        return vec![p[0], p[1]];
    };
    let distance = (p[0] - a[0]).hypot(p[1] - a[1]);
    match tool {
        DrawTool::Rectangle => vec![(p[0] - a[0]).abs(), (p[1] - a[1]).abs()],
        DrawTool::Circle => vec![2. * distance],
        _ => vec![distance],
    }
}

/// Parse one entry in millimeters. Plain numbers stay literal; anything else
/// is bound as a parameter expression. Empty text means "not entered".
pub(crate) fn parse(
    session: &Session,
    field: Field,
    text: &str,
) -> Result<Option<(f64, Dimension)>, String> {
    let text = text.trim();
    if text.is_empty() {
        return Ok(None);
    }
    let dimension = match text.parse::<f64>() {
        Ok(value) => Dimension::Literal {
            value,
            unit: "mm".into(),
        },
        Err(_) => Dimension::Expression {
            expression: multicad_core::parameters::BoundExpression::bind(
                session.document(),
                text,
                "mm",
                None,
            )
            .map_err(|e| e.to_string())?,
        },
    };
    let mm = multicad_core::parameters::Unit::parse("mm").expect("millimeters");
    let value = dimension.value_in(mm, Some(session.parameters()))?;
    if field.size() && value <= 1e-7 {
        return Err("Enter a positive size".into());
    }
    Ok(Some((value, dimension)))
}

/// The canonical authored source of a sketch, without resolving its plane.
pub(crate) fn authored(document: &Document, id: EntityId) -> Result<Source, String> {
    match document.get(id).map(|e| &e.data) {
        Some(EntityData::Sketch { definition, .. }) => {
            Source::from_definition(definition).map_err(|e| e.to_string())
        }
        _ => Err("The sketch is missing".into()),
    }
}

#[derive(Clone, Debug, Default)]
pub(crate) struct InlineInput {
    stage: Option<(DrawTool, bool)>,
    pub(crate) texts: Vec<String>,
    errors: Vec<Option<String>>,
    focus_request: Option<usize>,
    /// The last placement shown, kept so the fields and preview stay while
    /// the pointer is over a field instead of the canvas.
    pub(crate) last_point: Option<[f64; 2]>,
    pub(crate) last_pointer: Option<Pos2>,
    pub(crate) commit_requested: bool,
    had_focus: bool,
    /// The fields were painted in the last frame.
    shown: bool,
}

fn field_id(index: usize) -> egui::Id {
    egui::Id::new(("sketch-inline-field", index))
}

impl SketchEditor {
    pub(super) fn inline_stage(&self) -> Option<(DrawTool, bool)> {
        let tool = self.tool.filter(|_| self.editing)?;
        let started = self.first.is_some();
        (!fields(tool, started).is_empty()).then_some((tool, started))
    }
    /// Starts a frame: a new tool or drawing stage starts with empty fields.
    pub(super) fn sync_inline_input(&mut self) {
        let stage = self.inline_stage();
        if stage != self.inline.stage {
            let count = stage.map_or(0, |(tool, started)| fields(tool, started).len());
            let keep_focus = self.inline.had_focus && stage.is_some();
            self.inline = InlineInput {
                stage,
                texts: vec![String::new(); count],
                errors: vec![None; count],
                focus_request: keep_focus.then_some(0),
                last_pointer: self.inline.last_pointer,
                ..Default::default()
            };
        }
    }
    /// Current entered values, recording feedback for invalid entries.
    pub(super) fn inline_values(&mut self, session: &Session) -> Vec<Option<(f64, Dimension)>> {
        let Some((tool, started)) = self.inline.stage else {
            return vec![];
        };
        let mut values = vec![];
        for (i, field) in fields(tool, started).iter().enumerate() {
            match parse(session, *field, &self.inline.texts[i]) {
                Ok(value) => {
                    self.inline.errors[i] = None;
                    values.push(value);
                }
                Err(error) => {
                    self.inline.errors[i] = Some(error);
                    values.push(None);
                }
            }
        }
        values
    }
    /// Apply entered quantities to a pointer placement.
    pub(super) fn inline_constrain(
        &mut self,
        session: &Session,
        tool: DrawTool,
        pointer: [f64; 2],
    ) -> [f64; 2] {
        let values: Vec<_> = self
            .inline_values(session)
            .into_iter()
            .map(|v| v.map(|v| v.0))
            .collect();
        constrain(tool, self.first, pointer, &values)
    }
    /// Whether an entry is invalid, which blocks committing the creation.
    pub(super) fn inline_invalid(&self) -> Option<&str> {
        self.inline
            .errors
            .iter()
            .flatten()
            .next()
            .map(String::as_str)
    }
    /// True while one of the inline fields owns keyboard focus.
    pub(crate) fn inline_input_active(&self) -> bool {
        self.inline.had_focus && self.inline.stage.is_some()
    }
    /// Whether the fields are shown, so Tab moves between them.
    pub(crate) fn inline_fields_visible(&self) -> bool {
        self.inline.shown && self.inline.stage.is_some()
    }
    /// Tab from the canvas enters the first (or, backwards, the last) field.
    pub(crate) fn focus_inline_field(&mut self, forward: bool) {
        let count = self.inline.texts.len();
        if count > 0 {
            self.inline.focus_request = Some(if forward { 0 } else { count - 1 });
        }
    }
    #[cfg(test)]
    pub(crate) fn pending_first(&self) -> Option<[f64; 2]> {
        self.first
    }
    /// Escape from a field drops the pending creation; the tool stays active.
    pub(crate) fn cancel_inline_creation(&mut self) {
        self.clear_drawing_progress();
        self.inline = InlineInput::default();
        self.snap_assist.clear();
        self.message = "Creation cancelled · document unchanged".into();
    }
    /// Draws the creation with the entered sizes as driving dimensions, in
    /// one transaction. Without entered sizes it is the ordinary draw.
    pub(super) fn draw_dimensioned(
        &mut self,
        session: &mut Session,
        request: &multicad_core::sketch_draw::DrawRequest,
        tool: DrawTool,
        channel: Channel,
    ) -> Result<(), String> {
        if let Some(error) = self.inline_invalid() {
            return Err(format!("{error}; the sketch is unchanged"));
        }
        let sizes: Vec<(Field, Dimension)> = match self.inline.stage {
            Some((stage_tool, true)) if stage_tool == tool => {
                let values = self.inline_values(session);
                fields(tool, true)
                    .iter()
                    .zip(values)
                    .filter_map(|(field, value)| value.map(|v| (*field, v.1)))
                    .collect()
            }
            _ => vec![],
        };
        if sizes.is_empty() {
            return self.draw(session, request, channel);
        }
        if !self.editing {
            return Err("Create or open a sketch first".into());
        }
        if self.has_pending_edits() {
            return Err("Apply or discard pending edits before drawing".into());
        }
        let id = self.selected.ok_or("Select a sketch")?;
        let drawn = multicad_solver::draw::prepare_draw(
            session,
            session.revision(),
            Actor::local(channel),
            id,
            self.evaluated.as_ref(),
            request,
        )
        .map_err(|e| e.to_string())?;
        let preview = drawn.edit.preview.ok_or_else(|| {
            format!(
                "Sketch rejected: {:?}",
                drawn.edit.evaluation.diagnostics.state
            )
        })?;
        let mut source = authored(preview.document(), id)?;
        let position = |source: &Source, point| {
            source
                .sketch
                .points
                .iter()
                .find(|p| p.id == point)
                .map(|p| p.position_mm)
        };
        // (curve, (start id, position), (end id, position)) of each new line.
        type End = (EntityId, [f64; 2]);
        let lines: Vec<(EntityId, End, End)> = drawn
            .created_curves
            .iter()
            .filter_map(|curve| {
                let curve = source.sketch.curves.iter().find(|c| c.id == *curve)?;
                let multicad_core::constraint_sketch::CurveKind::Line { start, end } = curve.kind
                else {
                    return None;
                };
                Some((
                    curve.id,
                    (start, position(&source, start)?),
                    (end, position(&source, end)?),
                ))
            })
            .collect();
        let horizontal = |vertical: bool| {
            lines
                .iter()
                .find(|(_, (_, a), (_, b))| {
                    let (along, across) = if vertical { (1, 0) } else { (0, 1) };
                    (a[across] - b[across]).abs() < 1e-9 && (a[along] - b[along]).abs() > 1e-9
                })
                .map(|l| l.0)
        };
        for (field, dimension) in sizes {
            let kind = match field {
                Field::Width => ConstraintKind::Length {
                    line: horizontal(false).ok_or("The rectangle has no horizontal side")?,
                    dimension,
                },
                Field::Height => ConstraintKind::Length {
                    line: horizontal(true).ok_or("The rectangle has no vertical side")?,
                    dimension,
                },
                Field::Diameter => ConstraintKind::Diameter {
                    curve: *drawn.created_curves.first().ok_or("No circle was drawn")?,
                    dimension,
                },
                Field::Length => ConstraintKind::Length {
                    line: lines.last().ok_or("No line was drawn")?.0,
                    dimension,
                },
                Field::X | Field::Y => continue,
            };
            source.add_constraint(kind);
        }
        let edit = multicad_solver::prepare_edit_with_report(
            session,
            session.revision(),
            Actor::local(channel),
            id,
            source,
            self.evaluated.as_ref(),
        )
        .map_err(|e| e.to_string())?;
        session
            .commit(
                edit.preview
                    .ok_or("The entered sizes conflict with the sketch; nothing was drawn")?,
            )
            .map_err(|e| e.to_string())?;
        self.last_draw_line = lines.last().map(|(_, start, end)| [*start, *end]);
        self.clear_drawing_progress();
        self.element = drawn.created_curves.first().copied();
        self.revision = session.revision();
        self.inline.texts.iter_mut().for_each(String::clear);
        self.message = "Geometry and entered dimensions committed · Undo removes both".into();
        Ok(())
    }
    /// No fields this frame: forget keyboard ownership.
    pub(super) fn inline_fields_ui_hidden(&mut self) {
        self.inline.had_focus = false;
        self.inline.shown = false;
    }
    /// Paints the floating fields and handles Tab, Enter and Escape in them.
    /// `point` is the current constrained placement.
    pub(super) fn inline_fields_ui(
        &mut self,
        ui: &egui::Ui,
        frame: Frame,
        rect: Rect,
        point: [f64; 2],
    ) {
        // A click in this frame may have advanced the drawing stage.
        self.sync_inline_input();
        let Some((tool, started)) = self.inline.stage else {
            self.inline_fields_ui_hidden();
            return;
        };
        let kinds = fields(tool, started);
        let ctx = ui.ctx().clone();
        let focused = (0..kinds.len()).find(|i| ctx.memory(|m| m.has_focus(field_id(*i))));
        // Typing a number while drawing starts entry in the first field, so
        // the fields never need a pointer click (Fusion-style).
        if focused.is_none() && !ctx.text_edit_focused() && self.inline.focus_request.is_none() {
            let typed = ctx.input_mut(|input| {
                let mut typed = String::new();
                input.events.retain(|event| match event {
                    egui::Event::Text(text)
                        if typed.is_empty()
                            && text
                                .chars()
                                .next()
                                .is_some_and(|c| c.is_ascii_digit() || "-.=".contains(c)) =>
                    {
                        typed.push_str(text);
                        false
                    }
                    egui::Event::Text(text) if !typed.is_empty() => {
                        typed.push_str(text);
                        false
                    }
                    _ => true,
                });
                typed
            });
            if !typed.is_empty() {
                self.inline.texts[0] = typed;
                self.inline.focus_request = Some(0);
            }
        }
        if let Some(i) = focused {
            let (back, tab, enter) = ctx.input_mut(|input| {
                (
                    input.consume_key(egui::Modifiers::SHIFT, egui::Key::Tab),
                    input.consume_key(egui::Modifiers::NONE, egui::Key::Tab),
                    input.consume_key(egui::Modifiers::NONE, egui::Key::Enter),
                )
            });
            let n = kinds.len();
            if tab {
                self.inline.focus_request = Some((i + 1) % n);
            } else if back {
                self.inline.focus_request = Some((i + n - 1) % n);
            }
            if enter {
                self.inline.commit_requested = true;
                ctx.request_repaint();
            }
        }
        let screen = |p: [f64; 2]| {
            let s = self.camera.project(frame.point(p).map(|v| v as f32), rect);
            Pos2::new(s[0], s[1])
        };
        let current = measured(tool, self.first, point);
        let anchors: Vec<Pos2> = match (self.first, tool) {
            (None, _) => {
                // Position fields follow the pointer; until it has been over
                // the canvas there is no place for them.
                let Some(at) = self.inline.last_pointer else {
                    self.inline_fields_ui_hidden();
                    return;
                };
                vec![at + Vec2::new(70., 36.), at + Vec2::new(70., 66.)]
            }
            (Some(a), DrawTool::Rectangle) => {
                let center = screen([(a[0] + point[0]) / 2., (a[1] + point[1]) / 2.]);
                let outward = |m: Pos2| {
                    let d = m - center;
                    m + if d.length() > 1e-3 {
                        d.normalized() * 26.
                    } else {
                        Vec2::new(0., 26.)
                    }
                };
                vec![
                    outward(screen([(a[0] + point[0]) / 2., a[1]])),
                    outward(screen([point[0], (a[1] + point[1]) / 2.])),
                ]
            }
            (Some(a), _) => {
                let (s, e) = (screen(a), screen(point));
                let d = e - s;
                let normal = if d.length() > 1e-3 {
                    d.normalized().rot90()
                } else {
                    Vec2::new(0., -1.)
                };
                vec![s + d * 0.5 + normal * 22.]
            }
        };
        let mut had_focus = false;
        // A field never sits under the pointer, which must keep placing.
        let size = Vec2::new(132., 34.);
        let pointer = self.inline.last_pointer;
        let clear = |at: Pos2| match pointer {
            Some(p) if Rect::from_min_size(at, size).expand(6.).contains(p) => {
                Pos2::new(at.x, p.y + 26.)
            }
            _ => at,
        };
        for (i, field) in kinds.iter().enumerate() {
            let at = clear(anchors[i.min(anchors.len() - 1)] - Vec2::new(46., 11.));
            let request = self.inline.focus_request == Some(i);
            let error = self.inline.errors.get(i).cloned().flatten();
            let text = &mut self.inline.texts[i];
            // Unfocused fields let the pointer through: it always places.
            let interactive = request || focused == Some(i);
            egui::Area::new(egui::Id::new(("sketch-inline-area", i)))
                .order(egui::Order::Foreground)
                .interactable(interactive)
                .fixed_pos(at)
                .show(&ctx, |ui| {
                    egui::Frame::popup(ui.style()).show(ui, |ui| {
                        ui.horizontal(|ui| {
                            ui.label(egui::RichText::new(field.label()).strong());
                            let response = ui.add(
                                egui::TextEdit::singleline(text)
                                    .id(field_id(i))
                                    .lock_focus(true)
                                    .desired_width(64.)
                                    .hint_text(format!("{:.2}", current[i])),
                            );
                            if request {
                                response.request_focus();
                            }
                            had_focus |= response.has_focus() || request;
                            let tip = if field.size() {
                                "Fix this size; the pointer still places the rest. Numbers are millimeters; =$name expressions bind parameters."
                            } else {
                                "Fix this sketch coordinate in millimeters; Enter places the point."
                            };
                            response.on_hover_text(tip);
                            ui.label("mm");
                        });
                        if let Some(error) = &error {
                            ui.colored_label(Color32::LIGHT_RED, error);
                        }
                    });
                });
        }
        self.inline.focus_request = None;
        self.inline.had_focus = had_focus;
        self.inline.shown = true;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn entered_quantities_hold_while_the_pointer_places_the_rest() {
        let a = Some([10., 5.]);
        // Rectangle: width and height keep the pointer's side.
        assert_eq!(
            constrain(DrawTool::Rectangle, a, [3., 40.], &[Some(20.), None]),
            [-10., 40.]
        );
        assert_eq!(
            constrain(DrawTool::Rectangle, a, [30., 0.], &[Some(20.), Some(8.)]),
            [30., -3.]
        );
        // Circle: the diameter fixes the radius along the pointer direction.
        let p = constrain(DrawTool::Circle, a, [10., 50.], &[Some(12.)]);
        assert!((p[0] - 10.).abs() < 1e-12 && (p[1] - 11.).abs() < 1e-12);
        // Line: the length fixes the distance; a pointer on the start uses +X.
        assert_eq!(
            constrain(DrawTool::Line, a, [10., 5.], &[Some(7.)]),
            [17., 5.]
        );
        // Before the first point the coordinates are fixed independently.
        assert_eq!(
            constrain(DrawTool::Circle, None, [1., 2.], &[None, Some(-4.)]),
            [1., -4.]
        );
        // Nothing entered leaves the pointer placement.
        assert_eq!(
            constrain(DrawTool::Rectangle, a, [1., 2.], &[None, None]),
            [1., 2.]
        );
        assert!(fields(DrawTool::Arc, true).is_empty());
        assert_eq!(
            fields(DrawTool::Rectangle, true),
            &[Field::Width, Field::Height]
        );
    }
}
