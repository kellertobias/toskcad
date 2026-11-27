//! CAD-416: the one numeric Inspector field layout.
//!
//! The label sits above a full-width input. The unit is a suffix atom painted
//! inside the input's own frame, right-aligned; it is never part of the edited
//! text, so it cannot be typed over, selected or copied. Field validation
//! messages show directly below the input.
use crate::egui;
use egui::emath::Numeric;
use std::ops::RangeInclusive;

/// How the label is drawn. Diagram-linked fields highlight the measurement
/// the diagram identifies and mark the one validation refers to.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) enum Emphasis {
    #[default]
    Normal,
    Highlighted,
    Invalid,
}

/// The responses of one painted field.
pub(crate) struct Shown {
    #[cfg_attr(not(test), allow(dead_code))]
    pub label: egui::Response,
    pub input: egui::Response,
}

/// Text shown inside the field for a canonical unit name. Counts and ratios
/// have no unit.
pub(crate) fn unit_symbol(unit: &str) -> &str {
    match unit {
        "deg" => "°",
        "count" | "ratio" | "" => "",
        "teeth_per_in" => "teeth/in",
        unit => unit,
    }
}

/// A numeric Inspector field description; paint it with [`Self::show`] for
/// free text (values, `=` expressions, `$parameter` references) or
/// [`Self::value`] for a plain bounded number.
pub(crate) struct NumericField<'a> {
    id: Option<egui::Id>,
    label: &'a str,
    unit: &'a str,
    help: &'a str,
    hint: &'a str,
    emphasis: Emphasis,
    error: Option<&'a str>,
}

impl<'a> NumericField<'a> {
    /// `unit` is a canonical unit name (`mm`, `deg`, `count`, …) or a symbol.
    pub fn new(label: &'a str, unit: &'a str) -> Self {
        Self {
            id: None,
            label,
            unit,
            help: "",
            hint: "",
            emphasis: Emphasis::Normal,
            error: None,
        }
    }
    /// Stable identity; without one the input takes egui's next auto ID.
    pub fn id(mut self, id: egui::Id) -> Self {
        self.id = Some(id);
        self
    }
    pub fn help(mut self, help: &'a str) -> Self {
        self.help = help;
        self
    }
    pub fn hint(mut self, hint: &'a str) -> Self {
        self.hint = hint;
        self
    }
    pub fn emphasis(mut self, emphasis: Emphasis) -> Self {
        self.emphasis = emphasis;
        self
    }
    pub fn error(mut self, error: Option<&'a str>) -> Self {
        self.error = error;
        self
    }

    /// Paint the field around `text`, which the caller owns and validates.
    pub fn show(self, ui: &mut egui::Ui, text: &mut String) -> Shown {
        ui.vertical(|ui| {
            let width = crate::controls::inspector_width(ui);
            let label = self.paint_label(ui);
            let mut edit = egui::TextEdit::singleline(text)
                .desired_width(width)
                .hint_text(self.hint);
            if let Some(id) = self.id {
                edit = edit.id(id);
            }
            let symbol = unit_symbol(self.unit);
            if !symbol.is_empty() {
                edit = edit.suffix(
                    egui::RichText::new(symbol.to_owned()).color(ui.visuals().weak_text_color()),
                );
            }
            let mut input = ui.add(edit).labelled_by(label.id);
            if !self.help.is_empty() {
                input = input.on_hover_text(self.help);
            }
            if input.changed() {
                // Validation of the new text runs after this paint; show its
                // message on the next pass.
                ui.ctx().request_repaint();
            }
            if let Some(error) = self.error.filter(|error| !error.is_empty()) {
                ui.add(
                    egui::Label::new(
                        egui::RichText::new(error)
                            .small()
                            .color(ui.visuals().error_fg_color),
                    )
                    .wrap(),
                );
            }
            Shown { label, input }
        })
        .inner
    }

    /// A plain bounded number (no expressions). The value updates while
    /// typing, as egui's drag value did: parseable text is clamped to `range`
    /// and applied on every change; unparseable text is ignored and reported
    /// below the field until focus leaves, which restores the value's text.
    pub fn value<N: Numeric>(
        self,
        ui: &mut egui::Ui,
        value: &mut N,
        range: RangeInclusive<N>,
    ) -> Shown {
        #[derive(Clone, Default)]
        struct Editing {
            text: String,
            value: f64,
        }
        let id = self
            .id
            .unwrap_or_else(|| ui.next_auto_id().with("numeric-value"));
        let current = value.to_f64();
        let focused = ui.memory(|memory| memory.has_focus(id));
        let mut text = ui
            .data(|data| data.get_temp::<Editing>(id))
            // The value may round-trip through another unit (degrees stored
            // as radians); tolerate that instead of discarding the typing.
            .filter(|editing| {
                focused && (editing.value - current).abs() <= 1e-9 * current.abs().max(1.)
            })
            .map_or_else(|| format_value::<N>(current), |editing| editing.text);
        let parse_error = focused && parse::<N>(&text).is_none();
        let error = if parse_error {
            Some("Enter a number")
        } else {
            self.error
        };
        let shown = NumericField {
            error,
            id: Some(id),
            ..self
        }
        .show(ui, &mut text);
        if shown.input.gained_focus() {
            select_all(ui, id, &text);
        }
        if shown.input.changed()
            && let Some(parsed) = parse::<N>(&text)
        {
            let (min, max) = (range.start().to_f64(), range.end().to_f64());
            *value = N::from_f64(parsed.clamp(min.min(max), max.max(min)));
        }
        if shown.input.has_focus() {
            let value = value.to_f64();
            ui.data_mut(|data| data.insert_temp(id, Editing { text, value }));
        } else {
            ui.data_mut(|data| data.remove_temp::<Editing>(id));
        }
        shown
    }

    fn paint_label(&self, ui: &mut egui::Ui) -> egui::Response {
        let text = egui::RichText::new(self.label);
        let text = match self.emphasis {
            Emphasis::Normal => text,
            Emphasis::Highlighted => text.color(ui.visuals().selection.stroke.color).strong(),
            Emphasis::Invalid => text.color(ui.visuals().error_fg_color).strong(),
        };
        ui.add(egui::Label::new(text).wrap())
    }
}

/// Plain numbers only; integer fields round to the nearest whole value.
fn parse<N: Numeric>(text: &str) -> Option<f64> {
    let text: String = text
        .trim()
        .chars()
        .map(|c| if c == '−' { '-' } else { c })
        .filter(|c| !c.is_whitespace())
        .collect();
    let value = text.parse::<f64>().ok().filter(|value| value.is_finite())?;
    Some(if N::INTEGRAL { value.round() } else { value })
}

fn format_value<N: Numeric>(value: f64) -> String {
    if N::INTEGRAL {
        format!("{value:.0}")
    } else {
        let text = format!("{value:.6}");
        let text = text.trim_end_matches('0').trim_end_matches('.');
        if text == "-0" {
            "0".into()
        } else {
            text.into()
        }
    }
}

fn select_all(ui: &egui::Ui, id: egui::Id, text: &str) {
    let mut state = egui::TextEdit::load_state(ui.ctx(), id).unwrap_or_default();
    state
        .cursor
        .set_char_range(Some(egui::text::CCursorRange::two(
            egui::text::CCursor::default(),
            egui::text::CCursor::new(text.chars().count()),
        )));
    state.store(ui.ctx(), id);
}

#[cfg(test)]
mod tests;
