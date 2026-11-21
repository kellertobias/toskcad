//! Shared explanatory diagrams for authoring forms (CAD-386, CAD-389,
//! CAD-394, CAD-385).
//!
//! A diagram is a labeled schematic of the selected family or variant. It is
//! display-only: every callout names one canonical form field, and clicking it
//! focuses that same field, so the diagram and the field list always edit one
//! staged value. Proportions follow the current values where they resolve, but
//! no geometry, preview or manufacturing output is derived from a diagram.
//!
//! The module provides four pieces that forms combine:
//! - [`Diagram`]: marks plus dimension callouts, painted to fit the Inspector.
//! - [`DiagramLink`]: links callouts and fields (click → focus, focus →
//!   highlight, validation → highlight and reveal).
//! - [`choice_grid`]: an illustrated, keyboard-accessible family selector.
//! - [`section`]: a named progressive-disclosure group that opens itself when
//!   it holds the field a callout or validation points to.
use super::dimension::DimensionField;
use crate::egui::{self, Align2, Color32, FontId, Pos2, Rect, Sense, Stroke, Vec2};
use multicad_core::{document::Document, parameters::EvaluatedParameters, parameters::Unit};

/// Stable identity of one canonical form field, such as `"washer.bore"`.
pub(crate) type Key = &'static str;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Line {
    /// Visible edge.
    Outline,
    /// Hidden edge or construction reference, dashed.
    Hidden,
    /// Axis or symmetry line, long dashes.
    Center,
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) enum Mark {
    /// Filled shapes must be convex; compose concave sections from pieces.
    Path {
        points: Vec<[f32; 2]>,
        closed: bool,
        line: Line,
        fill: bool,
    },
    Circle {
        center: [f32; 2],
        radius: f32,
        line: Line,
        fill: bool,
    },
    /// Short explanatory text placed in diagram coordinates.
    Note { at: [f32; 2], text: String },
}

/// How a callout's value is read: the symbol makes diameter versus radius
/// explicit on the drawing itself.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Notation {
    Length,
    Diameter,
    Radius,
    Angle,
}
impl Notation {
    pub(crate) fn prefix(self) -> &'static str {
        match self {
            Self::Diameter => "Ø ",
            Self::Radius => "R ",
            Self::Angle => "∠ ",
            Self::Length => "",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum Placement {
    /// A dimension line between `from` and `to`, offset perpendicular by this
    /// many screen points, with extension lines back to the endpoints.
    Linear(f32),
    /// A leader from `from` to a label at `to`, for angles and features.
    Leader,
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct Callout {
    pub key: Key,
    pub from: [f32; 2],
    pub to: [f32; 2],
    pub placement: Placement,
    pub notation: Notation,
    pub name: String,
}
impl Callout {
    pub(crate) fn label(&self, value: &str) -> String {
        format!("{}{} {value}", self.notation.prefix(), self.name)
    }
}

/// Diagram coordinates are millimetre-like with +Y up; painting fits them to
/// the available width.
#[derive(Clone, Debug, Default, PartialEq)]
pub(crate) struct Diagram {
    pub marks: Vec<Mark>,
    pub callouts: Vec<Callout>,
}
impl Diagram {
    pub(crate) fn new() -> Self {
        Self::default()
    }
    pub(crate) fn path(mut self, points: &[[f32; 2]], closed: bool, line: Line) -> Self {
        self.marks.push(Mark::Path {
            points: points.to_vec(),
            closed,
            line,
            fill: false,
        });
        self
    }
    /// A filled convex section outlined as a visible edge.
    pub(crate) fn solid(mut self, points: &[[f32; 2]]) -> Self {
        self.marks.push(Mark::Path {
            points: points.to_vec(),
            closed: true,
            line: Line::Outline,
            fill: true,
        });
        self
    }
    pub(crate) fn rect(self, x0: f32, y0: f32, x1: f32, y1: f32) -> Self {
        self.solid(&[[x0, y0], [x1, y0], [x1, y1], [x0, y1]])
    }
    pub(crate) fn outline_rect(self, x0: f32, y0: f32, x1: f32, y1: f32, line: Line) -> Self {
        self.path(&[[x0, y0], [x1, y0], [x1, y1], [x0, y1]], true, line)
    }
    pub(crate) fn line(self, a: [f32; 2], b: [f32; 2], line: Line) -> Self {
        self.path(&[a, b], false, line)
    }
    pub(crate) fn circle(mut self, center: [f32; 2], radius: f32, line: Line, fill: bool) -> Self {
        self.marks.push(Mark::Circle {
            center,
            radius,
            line,
            fill,
        });
        self
    }
    pub(crate) fn note(mut self, at: [f32; 2], text: impl Into<String>) -> Self {
        self.marks.push(Mark::Note {
            at,
            text: text.into(),
        });
        self
    }
    pub(crate) fn dim(
        mut self,
        key: Key,
        from: [f32; 2],
        to: [f32; 2],
        offset: f32,
        notation: Notation,
        name: impl Into<String>,
    ) -> Self {
        self.callouts.push(Callout {
            key,
            from,
            to,
            placement: Placement::Linear(offset),
            notation,
            name: name.into(),
        });
        self
    }
    pub(crate) fn leader(
        mut self,
        key: Key,
        from: [f32; 2],
        to: [f32; 2],
        notation: Notation,
        name: impl Into<String>,
    ) -> Self {
        self.callouts.push(Callout {
            key,
            from,
            to,
            placement: Placement::Leader,
            notation,
            name: name.into(),
        });
        self
    }
    pub(crate) fn keys(&self) -> Vec<Key> {
        self.callouts.iter().map(|c| c.key).collect()
    }
    fn bounds(&self) -> Option<Rect> {
        let mut points = vec![];
        for mark in &self.marks {
            match mark {
                Mark::Path { points: p, .. } => points.extend(p.iter().copied()),
                Mark::Circle { center, radius, .. } => {
                    points.push([center[0] - radius, center[1] - radius]);
                    points.push([center[0] + radius, center[1] + radius]);
                }
                Mark::Note { at, .. } => points.push(*at),
            }
        }
        for callout in &self.callouts {
            points.extend([callout.from, callout.to]);
        }
        let mut iter = points
            .into_iter()
            .filter(|p| p[0].is_finite() && p[1].is_finite());
        let first = iter.next()?;
        let mut rect = Rect::from_min_max(first.into(), first.into());
        for p in iter {
            rect.extend_with(p.into());
        }
        Some(rect)
    }
}

/// Model-to-screen mapping with +Y up.
#[derive(Clone, Copy, Debug)]
struct Fit {
    origin: Pos2,
    min: Pos2,
    scale: f32,
    height: f32,
}
impl Fit {
    fn new(bounds: Rect, area: Rect, pad: Vec2) -> Self {
        let size = bounds.size().max(Vec2::splat(1e-3));
        let inner = (area.size() - 2. * pad).max(Vec2::splat(1.));
        let scale = (inner.x / size.x).min(inner.y / size.y);
        let used = size * scale;
        let origin = area.min + pad + (inner - used) / 2.;
        Self {
            origin,
            min: bounds.min,
            scale,
            height: size.y,
        }
    }
    fn at(&self, p: [f32; 2]) -> Pos2 {
        Pos2::new(
            self.origin.x + (p[0] - self.min.x) * self.scale,
            self.origin.y + (self.height - (p[1] - self.min.y)) * self.scale,
        )
    }
}

/// Palette taken from the active theme so diagrams read in light and dark.
#[derive(Clone, Copy)]
struct Ink {
    edge: Color32,
    fill: Color32,
    dim: Color32,
    active: Color32,
    invalid: Color32,
    paper: Color32,
}
impl Ink {
    fn of(ui: &egui::Ui) -> Self {
        let v = ui.visuals();
        Self {
            edge: v.strong_text_color(),
            fill: v.faint_bg_color.gamma_multiply(1.6),
            dim: v.weak_text_color(),
            active: v.selection.stroke.color,
            invalid: v.error_fg_color,
            paper: v.extreme_bg_color,
        }
    }
}

fn paint_marks(painter: &egui::Painter, fit: &Fit, marks: &[Mark], ink: Ink, width: f32) {
    for mark in marks {
        match mark {
            Mark::Path {
                points,
                closed,
                line,
                fill,
            } => {
                let mut screen: Vec<Pos2> = points.iter().map(|p| fit.at(*p)).collect();
                if *fill && screen.len() > 2 {
                    painter.add(egui::Shape::convex_polygon(
                        screen.clone(),
                        ink.fill,
                        Stroke::NONE,
                    ));
                }
                if *closed && let Some(first) = screen.first().copied() {
                    screen.push(first);
                }
                stroke_path(painter, &screen, *line, ink, width);
            }
            Mark::Circle {
                center,
                radius,
                line,
                fill,
            } => {
                let c = fit.at(*center);
                let r = radius * fit.scale;
                if *fill {
                    painter.circle_filled(c, r, ink.fill);
                }
                if *line == Line::Outline {
                    painter.circle_stroke(c, r, Stroke::new(width, ink.edge));
                } else {
                    let steps = 48;
                    let points: Vec<Pos2> = (0..=steps)
                        .map(|i| {
                            let a = i as f32 / steps as f32 * std::f32::consts::TAU;
                            c + Vec2::angled(a) * r
                        })
                        .collect();
                    stroke_path(painter, &points, *line, ink, width);
                }
            }
            Mark::Note { at, text } => {
                painter.text(
                    fit.at(*at),
                    Align2::CENTER_CENTER,
                    text,
                    FontId::proportional(10.),
                    ink.dim,
                );
            }
        }
    }
}

fn stroke_path(painter: &egui::Painter, points: &[Pos2], line: Line, ink: Ink, width: f32) {
    match line {
        Line::Outline => {
            painter.add(egui::Shape::line(
                points.to_vec(),
                Stroke::new(width, ink.edge),
            ));
        }
        Line::Hidden => {
            painter.extend(egui::Shape::dashed_line(
                points,
                Stroke::new(width * 0.8, ink.dim),
                4.,
                3.,
            ));
        }
        Line::Center => {
            painter.extend(egui::Shape::dashed_line(
                points,
                Stroke::new(width * 0.7, ink.dim),
                9.,
                3.,
            ));
        }
    }
}

fn arrow(painter: &egui::Painter, tip: Pos2, toward: Vec2, color: Color32) {
    let dir = toward.normalized();
    if !dir.x.is_finite() {
        return;
    }
    let side = dir.rot90() * 2.6;
    let base = tip - dir * 6.;
    painter.add(egui::Shape::convex_polygon(
        vec![tip, base + side, base - side],
        color,
        Stroke::NONE,
    ));
}

/// What happened to a painted diagram this frame.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct Painted {
    pub clicked: Option<Key>,
    pub hovered: Option<Key>,
    /// Painted label rectangle of every callout, for tests and layout checks.
    pub labels: Vec<(Key, Rect)>,
    pub rect: Rect,
}

impl Default for Painted {
    fn default() -> Self {
        Self {
            clicked: None,
            hovered: None,
            labels: vec![],
            rect: Rect::NOTHING,
        }
    }
}

#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct Emphasis {
    pub active: Option<Key>,
    pub invalid: Option<Key>,
}

/// Paint `diagram` across the available Inspector width. `value` supplies the
/// staged text of each callout's field so labels show exactly what is typed.
pub(crate) fn show(
    ui: &mut egui::Ui,
    id_salt: impl egui::AsIdSalt,
    diagram: &Diagram,
    max_height: f32,
    emphasis: Emphasis,
    value: &dyn Fn(Key) -> String,
) -> Painted {
    let id = ui.make_persistent_id(("diagram", id_salt));
    let width = crate::controls::inspector_width(ui).max(40.);
    let font = FontId::proportional(11.);
    let pad = Vec2::new(10., 18.);
    let Some(bounds) = diagram.bounds() else {
        return Painted::default();
    };
    let aspect = bounds.height().max(1e-3) / bounds.width().max(1e-3);
    let height = ((width - 2. * pad.x) * aspect + 2. * pad.y).clamp(60., max_height.max(60.));
    let (rect, _) = ui.allocate_exact_size(Vec2::new(width, height), Sense::hover());
    let mut painted = Painted {
        rect,
        ..Default::default()
    };
    if !ui.is_rect_visible(rect) {
        return painted;
    }
    let ink = Ink::of(ui);
    let painter = ui.painter_at(rect.expand(1.));
    painter.rect_filled(rect, 4., ink.paper);
    let fit = Fit::new(bounds, rect, pad);
    paint_marks(&painter, &fit, &diagram.marks, ink, 1.4);
    let mut placed: Vec<Rect> = vec![];
    for (index, callout) in diagram.callouts.iter().enumerate() {
        let a = fit.at(callout.from);
        let b = fit.at(callout.to);
        let color = if emphasis.invalid == Some(callout.key) {
            ink.invalid
        } else if emphasis.active == Some(callout.key) {
            ink.active
        } else {
            ink.dim
        };
        let width = if color == ink.dim { 1. } else { 2. };
        let stroke = Stroke::new(width, color);
        let (anchor, normal, hit_line) = match callout.placement {
            Placement::Linear(offset) => {
                let along = b - a;
                let normal = if along.length() > 1e-3 {
                    along.normalized().rot90()
                } else {
                    Vec2::Y
                };
                let (pa, pb) = (a + normal * offset, b + normal * offset);
                for (end, start) in [(pa, a), (pb, b)] {
                    let gap = (end - start).normalized() * 2.;
                    if (end - start).length() > 3. {
                        painter.line_segment(
                            [start + gap, end + normal * offset.signum() * 3.],
                            Stroke::new(0.8, color),
                        );
                    }
                }
                painter.line_segment([pa, pb], stroke);
                if callout.notation != Notation::Radius {
                    arrow(&painter, pa, pa - pb, color);
                }
                arrow(&painter, pb, pb - pa, color);
                (pa + (pb - pa) / 2., normal * offset.signum(), [pa, pb])
            }
            Placement::Leader => {
                painter.line_segment([a, b], stroke);
                painter.circle_filled(a, 2., color);
                (b, (b - a).normalized(), [a, b])
            }
        };
        let text = callout.label(&value(callout.key));
        let galley = painter.layout_no_wrap(text, font.clone(), color);
        let size = galley.size() + Vec2::new(6., 2.);
        let mut center = anchor;
        if callout.placement != Placement::Leader {
            // Labels sit on their dimension line; shift outward on overlap.
            let normal = if normal.x.is_finite() {
                normal
            } else {
                Vec2::Y
            };
            for _ in 0..6 {
                let candidate = Rect::from_center_size(center, size);
                if !placed.iter().any(|r| r.intersects(candidate)) {
                    break;
                }
                center += normal * (size.y + 1.) + Vec2::new(0., 0.);
            }
        } else {
            center += normal * (size.y / 2.);
        }
        let mut label = Rect::from_center_size(center, size);
        // Keep every label readable inside the diagram, never clipped.
        let dx =
            (rect.left() + 1. - label.left()).max(0.) + (rect.right() - 1. - label.right()).min(0.);
        let dy =
            (rect.top() + 1. - label.top()).max(0.) + (rect.bottom() - 1. - label.bottom()).min(0.);
        label = label.translate(Vec2::new(dx, dy));
        placed.push(label);
        let response = ui
            .interact(
                label.union(Rect::from_two_pos(hit_line[0], hit_line[1]).expand(3.)),
                id.with(index),
                Sense::click(),
            )
            .on_hover_text(format!("Edit {}", callout.name))
            .on_hover_cursor(egui::CursorIcon::PointingHand);
        response.widget_info(|| {
            egui::WidgetInfo::labeled(
                egui::WidgetType::Button,
                true,
                callout.label(&value(callout.key)),
            )
        });
        let hot = response.hovered() || response.has_focus();
        painter.rect(
            label,
            3.,
            ink.paper,
            Stroke::new(
                if hot { 1.5 } else { 0.6 },
                if hot { ink.active } else { color },
            ),
            egui::StrokeKind::Inside,
        );
        painter.galley(label.min + Vec2::new(3., 1.), galley, color);
        if response.clicked() {
            painted.clicked = Some(callout.key);
        }
        if hot {
            painted.hovered = Some(callout.key);
        }
        painted.labels.push((callout.key, label));
    }
    painted
}

/// Paint only the marks of `diagram` inside `rect` (choice thumbnails).
pub(crate) fn thumbnail(ui: &egui::Ui, rect: Rect, diagram: &Diagram) {
    let Some(bounds) = diagram.bounds() else {
        return;
    };
    let fit = Fit::new(bounds, rect, Vec2::splat(4.));
    paint_marks(&ui.painter_at(rect), &fit, &diagram.marks, Ink::of(ui), 1.2);
}

/// Callout ↔ field linking state, retained by its form between frames.
#[derive(Clone, Debug, Default)]
pub(crate) struct DiagramLink {
    active: Option<Key>,
    focus: Option<Key>,
    /// The first field the current validation refers to.
    pub invalid: Option<Key>,
}
impl DiagramLink {
    #[cfg(test)]
    pub(crate) fn active(&self) -> Option<Key> {
        self.active
    }
    #[cfg(test)]
    pub(crate) fn pending_focus(&self) -> Option<Key> {
        self.focus
    }
    /// Focus `key`'s field on the next time it is painted (callout clicks).
    pub(crate) fn focus(&mut self, key: Key) {
        self.focus = Some(key);
        self.active = Some(key);
    }
    /// Whether a group holding `keys` must be revealed: it contains the
    /// requested focus or the field that currently fails validation.
    pub(crate) fn reveals(&self, keys: &[Key]) -> bool {
        [self.focus, self.invalid]
            .into_iter()
            .flatten()
            .any(|key| keys.contains(&key))
    }
    pub(crate) fn emphasis(&self) -> Emphasis {
        Emphasis {
            active: self.active,
            invalid: self.invalid,
        }
    }
    pub(crate) fn diagram(
        &mut self,
        ui: &mut egui::Ui,
        id_salt: impl egui::AsIdSalt,
        diagram: &Diagram,
        max_height: f32,
        value: &dyn Fn(Key) -> String,
    ) -> Painted {
        let painted = show(ui, id_salt, diagram, max_height, self.emphasis(), value);
        if let Some(key) = painted.clicked {
            self.focus(key);
        }
        painted
    }
    /// Paint a linked dimension field. Returns whether its text changed.
    pub(crate) fn field(
        &mut self,
        ui: &mut egui::Ui,
        key: Key,
        field: &mut DimensionField,
        id_salt: impl egui::AsIdSalt,
        label: &str,
        help: &str,
    ) -> bool {
        let response = field.show_linked(
            ui,
            id_salt,
            label,
            help,
            self.active == Some(key),
            self.invalid == Some(key),
        );
        self.track(key, &response);
        response.changed()
    }
    /// Link any focusable control (for example a drag value) to `key`.
    pub(crate) fn track(&mut self, key: Key, response: &egui::Response) {
        if self.focus == Some(key) {
            response.request_focus();
            response.scroll_to_me(None);
            self.focus = None;
        }
        if response.has_focus() {
            self.active = Some(key);
        } else if response.lost_focus() && self.active == Some(key) {
            self.active = None;
        }
    }
}

/// Text a callout shows for a field: the staged literal with its unit, or the
/// authored expression itself.
pub(crate) fn field_value(field: &DimensionField) -> String {
    let text = field.text.trim();
    if text.is_empty() {
        "?".into()
    } else if text.starts_with('=') || text.chars().any(char::is_alphabetic) {
        text.into()
    } else {
        match field.unit() {
            "count" | "ratio" => text.into(),
            "deg" => format!("{text}°"),
            "teeth_per_in" => format!("{text} teeth/in"),
            unit => format!("{text} {unit}"),
        }
    }
}

/// Current numeric value of a staged field — millimetres for lengths,
/// degrees for angles, the field's own unit otherwise — used only to
/// proportion schematics.
pub(crate) fn value_of(
    field: &DimensionField,
    document: &Document,
    parameters: &EvaluatedParameters,
) -> Option<f64> {
    let dimension = field.resolve(document).ok()?;
    let unit = Unit::parse(field.unit())?;
    // Lengths in millimetres, angles in degrees, counts and ratios as typed.
    let target = if unit.dimension() == Unit::Degree.dimension() {
        Unit::Degree
    } else if unit.dimension() == Unit::Millimeter.dimension() {
        Unit::Millimeter
    } else {
        unit
    };
    dimension
        .value_in(target, Some(parameters))
        .ok()
        .filter(|v| v.is_finite())
}

/// First key in display order whose field text does not resolve.
pub(crate) fn first_invalid<'a>(
    document: &Document,
    fields: impl IntoIterator<Item = (Key, &'a DimensionField)>,
) -> Option<Key> {
    fields
        .into_iter()
        .find(|(_, field)| field.resolve(document).is_err())
        .map(|(key, _)| key)
}

/// Clamp a proportion so a schematic stays legible at any authored value.
pub(crate) fn legible(value: Option<f64>, fallback: f32, min: f32, max: f32) -> f32 {
    value
        .map(|v| v as f32)
        .filter(|v| v.is_finite() && *v > 0.)
        .unwrap_or(fallback)
        .clamp(min, max)
}

pub(crate) struct Choice<'a, T> {
    pub value: T,
    pub label: &'a str,
    pub hint: &'a str,
    pub art: Diagram,
}

/// Illustrated single-choice selector. Tiles reflow into as many columns as
/// fit; Tab focuses tiles, Space/Enter selects, arrows move the selection.
pub(crate) fn choice_grid<T: Copy + PartialEq>(
    ui: &mut egui::Ui,
    id_salt: impl egui::AsIdSalt,
    selected: &mut T,
    choices: &[Choice<'_, T>],
    enabled: bool,
) -> bool {
    if choices.is_empty() {
        return false;
    }
    let id = ui.make_persistent_id(("choice_grid", id_salt));
    let width = crate::controls::inspector_width(ui).max(40.);
    let gap = 4.;
    let min_tile = 74.;
    let columns = (((width + gap) / (min_tile + gap)).floor() as usize).clamp(1, choices.len());
    let tile_w = (width - gap * (columns - 1) as f32) / columns as f32;
    let art_h = (tile_w * 0.55).clamp(30., 56.);
    let font = FontId::proportional(11.5);
    let galleys: Vec<_> = choices
        .iter()
        .map(|c| {
            ui.painter().layout(
                c.label.to_owned(),
                font.clone(),
                ui.visuals().text_color(),
                tile_w - 8.,
            )
        })
        .collect();
    let rows: Vec<f32> = galleys
        .chunks(columns)
        .map(|row| art_h + 10. + row.iter().map(|g| g.size().y).fold(0., f32::max))
        .collect();
    let total = rows.iter().sum::<f32>() + gap * (rows.len().saturating_sub(1)) as f32;
    let (rect, _) = ui.allocate_exact_size(Vec2::new(width, total), Sense::hover());
    let mut changed = false;
    let mut move_to = None;
    let mut y = rect.top();
    for (row, height) in rows.iter().enumerate() {
        for column in 0..columns {
            let index = row * columns + column;
            let Some(choice) = choices.get(index) else {
                break;
            };
            let tile = Rect::from_min_size(
                Pos2::new(rect.left() + column as f32 * (tile_w + gap), y),
                Vec2::new(tile_w, *height),
            );
            let is_selected = *selected == choice.value;
            let response = ui
                .interact(
                    tile,
                    id.with(index),
                    if enabled {
                        Sense::click()
                    } else {
                        Sense::hover()
                    },
                )
                .on_hover_text(choice.hint);
            response.widget_info(|| {
                egui::WidgetInfo::selected(
                    egui::WidgetType::SelectableLabel,
                    enabled,
                    is_selected,
                    choice.label,
                )
            });
            if enabled && response.clicked() {
                // Keep keyboard focus on the chosen tile so arrows continue.
                response.request_focus();
                if !is_selected {
                    *selected = choice.value;
                    changed = true;
                }
            }
            if enabled && response.has_focus() {
                let (prev, next) = ui.input(|i| {
                    (
                        i.key_pressed(egui::Key::ArrowLeft) || i.key_pressed(egui::Key::ArrowUp),
                        i.key_pressed(egui::Key::ArrowRight) || i.key_pressed(egui::Key::ArrowDown),
                    )
                });
                if prev && index > 0 {
                    move_to = Some(index - 1);
                } else if next && index + 1 < choices.len() {
                    move_to = Some(index + 1);
                }
            }
            if ui.is_rect_visible(tile) {
                let v = ui.visuals();
                let (fill, stroke) = if is_selected {
                    (
                        v.selection.bg_fill.gamma_multiply(0.35),
                        Stroke::new(2., v.selection.stroke.color),
                    )
                } else if response.hovered() {
                    (v.widgets.hovered.weak_bg_fill, v.widgets.hovered.bg_stroke)
                } else {
                    (v.faint_bg_color, v.widgets.inactive.bg_stroke)
                };
                let stroke = if response.has_focus() {
                    Stroke::new(2., v.selection.stroke.color)
                } else {
                    stroke
                };
                let painter = ui.painter();
                painter.rect(tile, 4., fill, stroke, egui::StrokeKind::Inside);
                let art = Rect::from_min_size(
                    tile.min + Vec2::new(3., 3.),
                    Vec2::new(tile_w - 6., art_h),
                );
                thumbnail(ui, art, &choice.art);
                let galley = galleys[index].clone();
                let at = Pos2::new(tile.center().x - galley.size().x / 2., art.bottom() + 3.);
                let color = if enabled {
                    ui.visuals().text_color()
                } else {
                    ui.visuals().weak_text_color()
                };
                ui.painter().galley(at, galley, color);
            }
        }
        y += height + gap;
    }
    if let Some(index) = move_to {
        *selected = choices[index].value;
        changed = true;
        ui.memory_mut(|m| m.request_focus(id.with(index)));
    }
    changed
}

/// A named progressive-disclosure group. Closed by default; `reveal` opens it
/// (and it stays open) when it holds the field a callout or validation needs.
/// Collapsing never resets the fields inside; they remain in their form.
pub(crate) fn section<R>(
    ui: &mut egui::Ui,
    id_salt: impl egui::AsIdSalt,
    title: &str,
    reveal: bool,
    body: impl FnOnce(&mut egui::Ui) -> R,
) -> Option<R> {
    let width =
        (ui.available_width() - ui.spacing().indent - 2. * ui.spacing().button_padding.x).max(0.);
    let text = egui::WidgetText::from(egui::RichText::new(title).strong()).into_galley(
        ui,
        Some(egui::TextWrapMode::Truncate),
        width,
        egui::TextStyle::Button,
    );
    let mut header = egui::CollapsingHeader::new(text)
        .id_salt(id_salt)
        .default_open(false);
    if reveal {
        header = header.open(Some(true));
    }
    let response = header.show(ui, body);
    response.header_response.on_hover_text(title);
    response.body_returned
}

#[cfg(test)]
mod tests;
