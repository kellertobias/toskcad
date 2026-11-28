//! CAD-416: layout and editing of the shared numeric Inspector field.
use super::*;
use egui::accesskit::Role;

const WIDTH: f32 = 236.;

struct Painted {
    /// Content rect of the panel the field was painted in.
    content: egui::Rect,
    nodes: Vec<(Role, Option<String>, Option<String>, egui::Rect)>,
    /// Text of the nodes each text input is labelled by.
    labelled_by: Vec<String>,
    texts: Vec<(String, egui::Rect)>,
}
impl Painted {
    fn node(&self, role: Role) -> (Option<&str>, Option<&str>, egui::Rect) {
        let (_, label, value, rect) = self
            .nodes
            .iter()
            .find(|(r, ..)| *r == role)
            .unwrap_or_else(|| panic!("{role:?} node"));
        (label.as_deref(), value.as_deref(), *rect)
    }
    fn text(&self, text: &str) -> Option<egui::Rect> {
        self.texts.iter().find(|(t, _)| t == text).map(|(_, r)| *r)
    }
}

struct Harness {
    ctx: egui::Context,
    time: f64,
}
impl Harness {
    fn new() -> Self {
        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        crate::theme::configure(&ctx);
        Self { ctx, time: 0. }
    }
    fn frame(&mut self, events: Vec<egui::Event>, body: impl FnMut(&mut egui::Ui)) -> Painted {
        let mut body = body;
        self.time += 1. / 60.;
        let mut content = egui::Rect::NOTHING;
        let out = self.ctx.run_ui(
            egui::RawInput {
                screen_rect: Some(egui::Rect::from_min_size(
                    egui::Pos2::ZERO,
                    egui::vec2(400., 300.),
                )),
                events,
                focused: true,
                time: Some(self.time),
                ..Default::default()
            },
            |ui| {
                egui::Frame::NONE.inner_margin(12).show(ui, |ui| {
                    ui.set_width(WIDTH);
                    content = ui.max_rect();
                    crate::controls::inspector_layout(ui);
                    body(ui);
                });
            },
        );
        let rect = |b: egui::accesskit::Rect| {
            egui::Rect::from_min_max(
                egui::pos2(b.x0 as f32, b.y0 as f32),
                egui::pos2(b.x1 as f32, b.y1 as f32),
            )
        };
        let nodes = out
            .platform_output
            .accesskit_update
            .as_ref()
            .map(|update| {
                update
                    .nodes
                    .iter()
                    .filter_map(|(_, node)| {
                        Some((
                            node.role(),
                            node.label().map(str::to_owned),
                            node.value().map(str::to_owned),
                            rect(node.bounds()?),
                        ))
                    })
                    .collect()
            })
            .unwrap_or_default();
        let labelled_by = out
            .platform_output
            .accesskit_update
            .as_ref()
            .map(|update| {
                let text = |id: &egui::accesskit::NodeId| {
                    update
                        .nodes
                        .iter()
                        .find(|(n, _)| n == id)
                        .and_then(|(_, node)| node.label().or(node.value()).map(str::to_owned))
                };
                update
                    .nodes
                    .iter()
                    .filter(|(_, node)| node.role() == Role::TextInput)
                    .flat_map(|(_, node)| node.labelled_by().iter().filter_map(text))
                    .collect()
            })
            .unwrap_or_default();
        let texts = out
            .shapes
            .iter()
            .filter_map(|clipped| match &clipped.shape {
                egui::Shape::Text(text) => Some((
                    text.galley.text().to_owned(),
                    text.galley.rect.translate(text.pos.to_vec2()),
                )),
                _ => None,
            })
            .collect();
        out.drop_without_applying_deltas();
        Painted {
            content,
            nodes,
            labelled_by,
            texts,
        }
    }
    fn click(&mut self, pos: egui::Pos2, body: &mut impl FnMut(&mut egui::Ui)) {
        for pressed in [true, false] {
            self.frame(
                vec![
                    egui::Event::PointerMoved(pos),
                    egui::Event::PointerButton {
                        pos,
                        button: egui::PointerButton::Primary,
                        pressed,
                        modifiers: Default::default(),
                    },
                ],
                &mut *body,
            );
        }
        self.frame(vec![], &mut *body);
    }
    fn key(&mut self, key: egui::Key, body: &mut impl FnMut(&mut egui::Ui)) {
        for pressed in [true, false] {
            self.frame(
                vec![egui::Event::Key {
                    key,
                    physical_key: None,
                    pressed,
                    repeat: false,
                    modifiers: Default::default(),
                }],
                &mut *body,
            );
        }
    }
}

fn assert_layout(painted: &Painted, label: &str, value: &str, unit: &str) {
    let (_, input_value, input) = painted.node(Role::TextInput);
    // The input spans the full panel content width.
    assert!(
        (input.width() - painted.content.width()).abs() < 0.5
            && (input.left() - painted.content.left()).abs() < 0.5,
        "input {input:?} spans content {:?}",
        painted.content
    );
    // The label is its own text above the input and names it.
    let label_rect = painted.text(label).expect("label painted");
    assert!(
        label_rect.bottom() <= input.top() + 0.5,
        "label above input"
    );
    assert!(label_rect.left() >= input.left() - 0.5);
    assert_eq!(painted.labelled_by, [label]);
    // The editable value is the number alone; the unit is painted inside the
    // input, right-aligned, as separate text.
    assert_eq!(input_value, Some(value));
    assert!(!input_value.unwrap().contains(unit));
    let unit_rect = painted.text(unit).expect("unit painted");
    assert!(
        input.contains_rect(unit_rect),
        "unit {unit_rect:?} inside {input:?}"
    );
    assert!(
        input.right() - unit_rect.right() < 12.,
        "unit right-aligned: {unit_rect:?} in {input:?}"
    );
    let value_rect = painted.text(value).expect("value painted");
    assert!(
        value_rect.right() < unit_rect.left(),
        "unit after the value"
    );
}

#[test]
fn text_field_is_full_width_with_label_above_and_unit_inside() {
    let mut h = Harness::new();
    let mut text = "12".to_string();
    let mut error: Option<String> = None;
    let mut body = |ui: &mut egui::Ui| {
        let shown = NumericField::new("Distance", "mm")
            .id(egui::Id::new("distance"))
            .hint("Value or =$parameter")
            .error(error.as_deref())
            .show(ui, &mut text);
        assert!(shown.label.rect.bottom() <= shown.input.rect.top());
    };
    h.frame(vec![], &mut body);
    let painted = h.frame(vec![], &mut body);
    assert_layout(&painted, "Distance", "12", "mm");
    // Clicking the unit focuses the input; the unit itself cannot be selected
    // into the value: select-all and typing replace only the number.
    let (_, _, input) = painted.node(Role::TextInput);
    let unit = painted.text("mm").unwrap();
    h.click(unit.center(), &mut body);
    assert!(h.ctx.text_edit_focused());
    h.frame(
        vec![egui::Event::Key {
            key: egui::Key::A,
            physical_key: None,
            pressed: true,
            repeat: false,
            modifiers: egui::Modifiers::COMMAND,
        }],
        &mut body,
    );
    h.frame(vec![egui::Event::Text("=$width / 2".into())], &mut body);
    drop(body);
    assert_eq!(text, "=$width / 2");
    // A validation error renders below the input.
    error = Some("Unknown parameter $width".into());
    let mut body = |ui: &mut egui::Ui| {
        NumericField::new("Distance", "mm")
            .id(egui::Id::new("distance"))
            .error(error.as_deref())
            .show(ui, &mut text);
    };
    let painted = h.frame(vec![], &mut body);
    let message = painted.text("Unknown parameter $width").expect("error");
    assert!(
        message.top() >= input.bottom() - 0.5,
        "error below the input"
    );
}

#[test]
fn angle_and_count_units_use_their_symbol_or_none() {
    assert_eq!(unit_symbol("deg"), "°");
    assert_eq!(unit_symbol("count"), "");
    assert_eq!(unit_symbol("ratio"), "");
    assert_eq!(unit_symbol("in"), "in");
    let mut h = Harness::new();
    let mut angle = 45.;
    let mut body = |ui: &mut egui::Ui| {
        NumericField::new("Angle", "deg")
            .id(egui::Id::new("angle"))
            .value(ui, &mut angle, -360.0..=360.0);
    };
    h.frame(vec![], &mut body);
    let painted = h.frame(vec![], &mut body);
    assert_layout(&painted, "Angle", "45", "°");
}

#[test]
fn value_field_commits_clamps_rounds_and_reports_unparseable_text() {
    let mut h = Harness::new();
    let mut count: u32 = 4;
    let mut body = |ui: &mut egui::Ui| {
        NumericField::new("Count", "count")
            .id(egui::Id::new("count"))
            .value(ui, &mut count, 2..=64);
    };
    h.frame(vec![], &mut body);
    let painted = h.frame(vec![], &mut body);
    let (_, value, input) = painted.node(Role::TextInput);
    assert_eq!(value, Some("4"));
    // Focus selects the whole value, so typing replaces it.
    h.click(input.center(), &mut body);
    h.frame(vec![egui::Event::Text("7.6".into())], &mut body);
    h.frame(vec![], &mut body);
    h.key(egui::Key::Enter, &mut body);
    drop(body);
    assert_eq!(count, 8, "integer fields round");

    let mut h = Harness::new();
    let mut distance = 10.;
    let mut body = |ui: &mut egui::Ui| {
        NumericField::new("Offset", "mm")
            .id(egui::Id::new("offset"))
            .value(ui, &mut distance, 0.0..=100.0);
    };
    h.frame(vec![], &mut body);
    let painted = h.frame(vec![], &mut body);
    let (_, _, input) = painted.node(Role::TextInput);
    h.click(input.center(), &mut body);
    h.frame(vec![egui::Event::Text("250".into())], &mut body);
    drop(body);
    assert_eq!(distance, 100., "clamped to the range while typing");
    let mut body = |ui: &mut egui::Ui| {
        NumericField::new("Offset", "mm")
            .id(egui::Id::new("offset"))
            .value(ui, &mut distance, 0.0..=100.0);
    };
    h.frame(
        vec![egui::Event::Key {
            key: egui::Key::A,
            physical_key: None,
            pressed: true,
            repeat: false,
            modifiers: egui::Modifiers::COMMAND,
        }],
        &mut body,
    );
    h.frame(vec![egui::Event::Text("abc".into())], &mut body);
    let painted = h.frame(vec![], &mut body);
    assert!(
        painted.text("Enter a number").is_some(),
        "error below field"
    );
    // Escape leaves the field; the last valid value stays and its text returns.
    h.key(egui::Key::Escape, &mut body);
    let painted = h.frame(vec![], &mut body);
    assert!(!h.ctx.text_edit_focused());
    assert_eq!(painted.node(Role::TextInput).1, Some("100"));
    assert!(painted.text("Enter a number").is_none());
    drop(body);
    assert_eq!(distance, 100.);
}
