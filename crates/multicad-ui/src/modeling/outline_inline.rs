//! Inline segment entry while a wall path is drawn (CAD-427).
//!
//! After the first corner, two fields float beside the pending segment:
//! Length first, then Angle (Tab moves between them). An entered length
//! fixes the segment's length while the direction follows the pointer,
//! snapped to 0°, 90°, 180° or 270°; an entered angle fixes the direction
//! exactly (absolute, counter-clockwise from +X). Enter places the corner
//! and the next segment starts from it, so a room is typed as a chain of
//! lengths. Lengths accept explicit units (`3.5 m`, `350cm`, `12 in`);
//! unqualified numbers use the Settings length unit, and `=` starts a
//! parameter expression. Typing a digit while drawing enters the Length
//! field, so the fields never need a pointer click. Unfocused fields never
//! intercept the pointer, which keeps placing corners by clicking.
use super::*;

/// One inline field.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum SegmentField {
    Length,
    Angle,
}
impl SegmentField {
    pub(crate) const ALL: [Self; 2] = [Self::Length, Self::Angle];
    fn label(self) -> &'static str {
        match self {
            Self::Length => "L",
            Self::Angle => "∠",
        }
    }
}

#[derive(Clone, Debug, Default)]
pub(super) struct InlineSegment {
    pub(super) texts: [String; 2],
    errors: [Option<String>; 2],
    focus_request: Option<usize>,
    /// A field owned keyboard focus in the last frame.
    had_focus: bool,
    /// The fields were painted in the last frame.
    shown: bool,
    /// Enter in a field asks to place the shown corner.
    commit_requested: bool,
}

fn field_id(index: usize) -> egui::Id {
    egui::Id::new(("wall-inline-segment", index))
}

/// The unit direction of an absolute angle in degrees; multiples of 90° are
/// exact so typed rooms close without rounding.
pub(super) fn exact_direction(angle_deg: f64) -> [f64; 2] {
    let a = angle_deg.rem_euclid(360.);
    for (quarter, d) in [
        (0., [1., 0.]),
        (90., [0., 1.]),
        (180., [-1., 0.]),
        (270., [0., -1.]),
        (360., [1., 0.]),
    ] {
        if (a - quarter).abs() < 1e-9 {
            return d;
        }
    }
    let (s, c) = a.to_radians().portable_sin_cos();
    [c, s]
}

/// The axis direction (0°, 90°, 180° or 270°) closest to `toward - from`;
/// +X when the pointer sits on the start.
pub(super) fn axis_direction(from: [f64; 2], toward: [f64; 2]) -> [f64; 2] {
    let d = sub(toward, from);
    if d[0].abs() < 1e-9 && d[1].abs() < 1e-9 {
        return [1., 0.];
    }
    if d[0].abs() >= d[1].abs() {
        [d[0].signum(), 0.]
    } else {
        [0., d[1].signum()]
    }
}

/// The next corner from `last` with the entered quantities held fixed:
/// a length alone follows the pointer's axis direction, an angle alone
/// projects the pointer onto that direction, both fix the corner. `None`
/// when nothing is entered.
pub(super) fn typed_corner(
    last: [f64; 2],
    pointer: Option<[f64; 2]>,
    length_mm: Option<f64>,
    angle_deg: Option<f64>,
) -> Option<[f64; 2]> {
    if length_mm.is_none() && angle_deg.is_none() {
        return None;
    }
    let toward = pointer.unwrap_or([last[0] + 1., last[1]]);
    let d = angle_deg.map_or_else(|| axis_direction(last, toward), exact_direction);
    let l = length_mm.unwrap_or_else(|| {
        let p = sub(toward, last);
        (p[0] * d[0] + p[1] * d[1]).max(0.)
    });
    Some([last[0] + d[0] * l, last[1] + d[1] * l])
}

/// Parse an entered length into millimeters. Empty text is "not entered".
pub(super) fn parse_length_mm(
    document: &Document,
    parameters: Option<&EvaluatedParameters>,
    text: &str,
) -> Result<Option<f64>, String> {
    if text.trim().is_empty() {
        return Ok(None);
    }
    let mut field = mm(0.);
    field.text = text.trim().to_string();
    let dimension = field.resolve(document)?;
    let evaluated;
    let parameters = match parameters {
        Some(p) => p,
        None => {
            evaluated = multicad_core::parameters::evaluate(document).map_err(|e| e.to_string())?;
            &evaluated
        }
    };
    // Exact rational conversion: `12 in` is 304.8 mm, not its float product.
    let value = multicad_core::parameters::exact_length_mm(&dimension, parameters)?.as_f64();
    if !(value > 1e-6 && value.is_finite()) {
        return Err("Enter a positive length".into());
    }
    Ok(Some(value))
}

/// Parse an entered absolute angle in degrees (`90`, `-45`, `22,5°`).
pub(super) fn parse_angle_deg(text: &str) -> Result<Option<f64>, String> {
    let text = text.trim();
    if text.is_empty() {
        return Ok(None);
    }
    let number = text
        .trim_end_matches('°')
        .trim_end_matches("deg")
        .trim()
        .replace(',', ".");
    let value: f64 = number
        .parse()
        .map_err(|_| format!("'{text}' is not an angle in degrees"))?;
    if !value.is_finite() || value.abs() > 360. {
        return Err("Enter an angle between -360° and 360°".into());
    }
    Ok(Some(value))
}

impl OutlineEditor {
    /// Whether the inline fields belong to the current frame.
    fn inline_stage(&self) -> bool {
        self.inline_entry && self.mode == Mode::Drawing && !self.corners.is_empty()
    }
    /// Entered values; invalid entries record their message and count as
    /// not entered.
    pub(super) fn inline_values(&mut self, document: &Document) -> (Option<f64>, Option<f64>) {
        let parameters = self
            .source
            .as_ref()
            .filter(|s| std::ptr::eq(s.session().document(), document))
            .map(|s| s.session().parameters().clone());
        let length = parse_length_mm(document, parameters.as_ref(), &self.inline.texts[0]);
        let angle = parse_angle_deg(&self.inline.texts[1]);
        self.inline.errors = [length.clone().err(), angle.clone().err()];
        (length.ok().flatten(), angle.ok().flatten())
    }
    /// The first invalid entry, which blocks placing a typed corner.
    pub(super) fn inline_invalid(&self) -> Option<&str> {
        self.inline
            .errors
            .iter()
            .flatten()
            .next()
            .map(String::as_str)
    }
    /// True while one of the inline fields owns keyboard focus.
    #[cfg(test)]
    pub(crate) fn inline_input_active(&self) -> bool {
        self.inline.had_focus && self.inline_stage()
    }
    /// Whether the fields are shown, so Tab from the canvas enters them.
    pub(crate) fn inline_fields_visible(&self) -> bool {
        self.inline.shown && self.inline_stage()
    }
    /// Tab from the canvas enters Length (or, backwards, Angle).
    pub(crate) fn focus_inline_field(&mut self, forward: bool) {
        self.inline.focus_request = Some(if forward { 0 } else { 1 });
    }
    /// Type into one inline field, as the keyboard does.
    #[cfg(test)]
    pub(crate) fn set_inline_text(&mut self, field: SegmentField, text: &str) {
        self.inline.texts[field as usize] = text.into();
    }
    #[cfg(test)]
    pub(crate) fn inline_text(&self, field: SegmentField) -> &str {
        &self.inline.texts[field as usize]
    }
    /// Ask to place the typed corner, as Enter in a field does.
    #[cfg(test)]
    pub(crate) fn request_inline_commit(&mut self) {
        self.inline.commit_requested = true;
    }
    pub(super) fn take_inline_commit(&mut self) -> bool {
        std::mem::take(&mut self.inline.commit_requested)
    }
    pub(super) fn clear_inline(&mut self) {
        self.inline.texts = Default::default();
        self.inline.errors = Default::default();
    }
    /// Keyboard handling before the path keys: a typed digit enters Length,
    /// and Tab, Shift+Tab and Enter in a focused field are consumed here so
    /// they never finish the path or switch workbenches.
    pub(super) fn inline_keyboard(&mut self, ctx: &egui::Context) {
        if !self.inline_stage() {
            return;
        }
        let focused = (0..2).find(|i| ctx.memory(|m| m.has_focus(field_id(*i))));
        if focused.is_none() && !ctx.text_edit_focused() && self.inline.focus_request.is_none() {
            let typed = ctx.input_mut(|input| {
                let mut typed = String::new();
                input.events.retain(|event| match event {
                    egui::Event::Text(text)
                        if typed.is_empty()
                            && text
                                .chars()
                                .next()
                                .is_some_and(|c| c.is_ascii_digit() || ".=".contains(c)) =>
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
            if tab || back {
                self.inline.focus_request = Some(1 - i);
            }
            if enter {
                self.inline.commit_requested = true;
                ctx.request_repaint();
            }
        }
    }
    /// Paint the floating fields beside the pending segment `last → point`.
    pub(super) fn inline_fields_ui(
        &mut self,
        view: &TerminationViewport<'_>,
        z: f64,
        last: [f64; 2],
        point: [f64; 2],
    ) {
        if !self.inline_stage() {
            self.inline.had_focus = false;
            self.inline.shown = false;
            return;
        }
        let ctx = view.painter.ctx().clone();
        let focused = (0..2).find(|i| ctx.memory(|m| m.has_focus(field_id(*i))));
        let to_screen = |p: [f64; 2]| {
            let s = view
                .camera
                .project([p[0] as f32, p[1] as f32, z as f32], view.rect);
            egui::pos2(s[0], s[1])
        };
        let (s, e) = (to_screen(last), to_screen(point));
        let d = e - s;
        let normal = if d.length() > 1e-3 {
            d.normalized().rot90()
        } else {
            egui::vec2(0., -1.)
        };
        let anchor = s + d * 0.5 + normal * 30. - egui::vec2(66., 17.);
        let measured = sub(point, last);
        let hints = [
            format!("{:.0}", dist(point, last)),
            format!(
                "{:.1}",
                measured[1]
                    .portable_atan2(measured[0])
                    .to_degrees()
                    .rem_euclid(360.)
            ),
        ];
        // A field never sits under the pointer, which must keep placing.
        let size = egui::vec2(132., 34.);
        let pointer = view.pointer;
        let mut had_focus = false;
        for (i, field) in SegmentField::ALL.iter().enumerate() {
            let mut at = anchor + egui::vec2(0., i as f32 * 36.);
            if let Some(p) = pointer
                && egui::Rect::from_min_size(at, size).expand(6.).contains(p)
            {
                at = egui::pos2(at.x, p.y + 26. + i as f32 * 36.);
            }
            let request = self.inline.focus_request == Some(i);
            let interactive = request || focused == Some(i);
            let error = self.inline.errors[i].clone();
            let unit = match field {
                SegmentField::Length => crate::preferences::length_unit().symbol(),
                SegmentField::Angle => "°",
            };
            let text = &mut self.inline.texts[i];
            egui::Area::new(egui::Id::new(("wall-inline-area", i)))
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
                                    .hint_text(&hints[i]),
                            );
                            if request {
                                response.request_focus();
                            }
                            had_focus |= response.has_focus() || request;
                            response.on_hover_text(match field {
                                SegmentField::Length => {
                                    "Segment length, e.g. 3.5 m or 350cm; the direction snaps to 0°/90°/180°/270° toward the pointer. Tab: angle · Enter places the corner"
                                }
                                SegmentField::Angle => {
                                    "Absolute direction in degrees, counter-clockwise from +X. Enter places the corner"
                                }
                            });
                            ui.label(unit);
                        });
                        if let Some(error) = &error {
                            ui.colored_label(egui::Color32::LIGHT_RED, error);
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
    fn typed_segments_snap_to_axes_or_follow_an_exact_angle() {
        let last = [1000., 500.];
        // A length alone takes the pointer's nearest axis direction.
        assert_eq!(
            typed_corner(last, Some([1200., 2400.]), Some(3000.), None),
            Some([1000., 3500.])
        );
        assert_eq!(
            typed_corner(last, Some([-5., 400.]), Some(4000.), None),
            Some([-3000., 500.])
        );
        // No pointer, or a pointer on the start, draws along +X.
        assert_eq!(
            typed_corner(last, None, Some(250.), None),
            Some([1250., 500.])
        );
        // Quarter turns are exact; other angles use the angle itself.
        assert_eq!(
            typed_corner(last, Some([0., 0.]), Some(3000.), Some(-90.)),
            Some([1000., -2500.])
        );
        let p = typed_corner([0., 0.], None, Some(2000.), Some(30.)).unwrap();
        assert!((p[0].portable_hypot(p[1]) - 2000.).abs() < 1e-9 && p[0] > 1700.);
        assert!((p[1] - 1000.).abs() < 1e-9);
        // An angle alone projects the pointer onto that direction.
        assert_eq!(
            typed_corner([0., 0.], Some([300., 900.]), None, Some(90.)),
            Some([0., 900.])
        );
        assert_eq!(typed_corner(last, Some([0., 0.]), None, None), None);
    }
    #[test]
    fn lengths_keep_their_units_exactly_and_bad_entries_are_refused() {
        let document = Document::new("Inline wall").unwrap();
        let value = |text| parse_length_mm(&document, None, text);
        for (text, mm) in [
            ("3.5 m", 3500.),
            ("350cm", 3500.),
            ("4 m", 4000.),
            ("2,5 m", 2500.),
            ("1200 mm", 1200.),
            ("1 ft", 304.8),
            ("12 in", 304.8),
        ] {
            assert_eq!(value(text), Ok(Some(mm)), "{text}");
        }
        assert_eq!(value("  "), Ok(None));
        for text in ["-2 m", "0", "abc", "5 deg", "3 yards", "=$missing"] {
            assert!(value(text).is_err(), "{text}");
        }
        assert_eq!(parse_angle_deg("90"), Ok(Some(90.)));
        assert_eq!(parse_angle_deg("22,5°"), Ok(Some(22.5)));
        assert_eq!(parse_angle_deg("-45 deg"), Ok(Some(-45.)));
        assert_eq!(parse_angle_deg(""), Ok(None));
        assert!(parse_angle_deg("north").is_err());
        assert!(parse_angle_deg("720").is_err());
    }
}
