use super::*;
use multicad_core::constraint_sketch::Dimension;

fn mm(value: f64) -> DimensionField {
    DimensionField::from_dimension(
        &Dimension::Literal {
            value,
            unit: "mm".into(),
        },
        "mm",
    )
    .unwrap()
}

fn washer(bore: f32, outside: f32) -> Diagram {
    Diagram::new()
        .rect(-outside / 2., 0., -bore / 2., 2.)
        .rect(bore / 2., 0., outside / 2., 2.)
        .line([0., -1.], [0., 3.], Line::Center)
        .dim(
            "bore",
            [-bore / 2., 2.],
            [bore / 2., 2.],
            16.,
            Notation::Diameter,
            "Bore",
        )
        .dim(
            "outside",
            [-outside / 2., 0.],
            [outside / 2., 0.],
            -16.,
            Notation::Diameter,
            "Outside",
        )
        .dim(
            "thickness",
            [outside / 2., 0.],
            [outside / 2., 2.],
            -14.,
            Notation::Length,
            "Thickness",
        )
}

fn input(width: f32, events: Vec<egui::Event>) -> egui::RawInput {
    egui::RawInput {
        screen_rect: Some(Rect::from_min_size(Pos2::ZERO, Vec2::new(width, 900.))),
        focused: true,
        events,
        ..Default::default()
    }
}

fn click(pos: Pos2) -> Vec<egui::Event> {
    vec![
        egui::Event::PointerMoved(pos),
        egui::Event::PointerButton {
            pos,
            button: egui::PointerButton::Primary,
            pressed: true,
            modifiers: Default::default(),
        },
        egui::Event::PointerButton {
            pos,
            button: egui::PointerButton::Primary,
            pressed: false,
            modifiers: Default::default(),
        },
    ]
}

struct Form {
    link: DiagramLink,
    fields: [DimensionField; 3],
}
const KEYS: [Key; 3] = ["bore", "outside", "thickness"];
impl Form {
    fn frame(
        &mut self,
        ctx: &egui::Context,
        width: f32,
        events: Vec<egui::Event>,
    ) -> (Painted, Vec<String>) {
        let mut painted = Painted::default();
        let output = ctx.run_ui(input(width, events), |ui| {
            let values: Vec<String> = self.fields.iter().map(field_value).collect();
            let value = |key: Key| values[KEYS.iter().position(|k| *k == key).unwrap()].clone();
            painted = self
                .link
                .diagram(ui, "test", &washer(3.2, 7.), 200., &value);
            for (i, key) in KEYS.into_iter().enumerate() {
                self.link
                    .field(ui, key, &mut self.fields[i], ("field", i), key, "help");
            }
        });
        let text = output
            .shapes
            .iter()
            .filter_map(|s| match &s.shape {
                egui::Shape::Text(t) => Some(t.galley.text().to_owned()),
                _ => None,
            })
            .collect();
        output.drop_without_applying_deltas();
        (painted, text)
    }
}

#[test]
fn callouts_show_notation_and_staged_values_unclipped_at_compact_and_ordinary_widths() {
    for width in [160., 320.] {
        let ctx = egui::Context::default();
        let mut form = Form {
            link: DiagramLink::default(),
            fields: [mm(3.2), mm(7.), mm(0.5)],
        };
        form.frame(&ctx, width, vec![]);
        let (painted, text) = form.frame(&ctx, width, vec![]);
        assert!(text.contains(&"Ø Bore 3.2 mm".to_owned()), "{text:?}");
        assert!(text.contains(&"Ø Outside 7 mm".to_owned()), "{text:?}");
        assert!(text.contains(&"Thickness 0.5 mm".to_owned()), "{text:?}");
        assert_eq!(painted.labels.len(), 3);
        for (key, label) in &painted.labels {
            assert!(
                painted.rect.expand(0.5).contains_rect(*label),
                "{key} label {label:?} clipped by {:?} at width {width}",
                painted.rect
            );
        }
        for (i, (_, a)) in painted.labels.iter().enumerate() {
            for (_, b) in &painted.labels[i + 1..] {
                assert!(!a.intersects(*b), "overlapping labels at width {width}");
            }
        }
    }
}

#[test]
fn clicking_a_callout_focuses_its_field_and_typing_updates_the_same_value() {
    let ctx = egui::Context::default();
    let mut form = Form {
        link: DiagramLink::default(),
        fields: [mm(3.2), mm(7.), mm(0.5)],
    };
    form.frame(&ctx, 300., vec![]);
    let (painted, _) = form.frame(&ctx, 300., vec![]);
    let outside = painted
        .labels
        .iter()
        .find(|(k, _)| *k == "outside")
        .unwrap()
        .1
        .center();
    let (painted, _) = form.frame(&ctx, 300., click(outside));
    assert_eq!(painted.clicked, Some("outside"));
    assert_eq!(form.link.active(), Some("outside"));
    form.frame(&ctx, 300., vec![]);
    assert!(ctx.memory(|m| m.focused()).is_some());
    assert_eq!(form.link.pending_focus(), None, "focus request consumed");
    // Typing goes into the same canonical field the callout names.
    form.frame(
        &ctx,
        300.,
        vec![
            egui::Event::Key {
                key: egui::Key::A,
                physical_key: None,
                pressed: true,
                repeat: false,
                modifiers: egui::Modifiers::COMMAND,
            },
            egui::Event::Text("9".into()),
        ],
    );
    assert_eq!(form.fields[1].text, "9");
    assert_eq!(form.fields[0].text, "3.2");
    let (_, text) = form.frame(&ctx, 300., vec![]);
    assert!(text.contains(&"Ø Outside 9 mm".to_owned()), "{text:?}");
}

#[test]
fn focusing_a_list_field_highlights_its_measurement_and_validation_marks_it() {
    let ctx = egui::Context::default();
    let mut form = Form {
        link: DiagramLink::default(),
        fields: [mm(3.2), mm(7.), mm(0.5)],
    };
    form.frame(&ctx, 300., vec![]);
    assert_eq!(form.link.active(), None);
    form.link.focus("thickness");
    form.frame(&ctx, 300., vec![]);
    form.frame(&ctx, 300., vec![]);
    assert_eq!(form.link.active(), Some("thickness"));
    assert_eq!(form.link.emphasis().active, Some("thickness"));
    let document = Document::new("Validation").unwrap();
    form.fields[1].text = "=".into();
    form.link.invalid = first_invalid(&document, KEYS.into_iter().zip(form.fields.iter()));
    assert_eq!(form.link.invalid, Some("outside"));
    assert!(form.link.reveals(&["outside"]));
    assert!(
        !form.link.reveals(&["thickness"]),
        "the consumed focus request reveals nothing"
    );
}

#[test]
fn field_values_show_units_or_the_authored_expression() {
    let mut field = mm(4.);
    assert_eq!(field_value(&field), "4 mm");
    field.text = "=$bore + 1".into();
    assert_eq!(field_value(&field), "=$bore + 1");
    field.text = "1 in".into();
    assert_eq!(field_value(&field), "1 in");
    field.text = " ".into();
    assert_eq!(field_value(&field), "?");
    let document = Document::new("Values").unwrap();
    let parameters = multicad_core::parameters::evaluate(&document).unwrap();
    field.text = "1 in".into();
    assert_eq!(value_of(&field, &document, &parameters), Some(25.4));
    field.text = "=".into();
    assert_eq!(value_of(&field, &document, &parameters), None);
    assert_eq!(legible(None, 3., 1., 5.), 3.);
    assert_eq!(legible(Some(-2.), 3., 1., 5.), 3.);
    assert_eq!(legible(Some(50.), 3., 1., 5.), 5.);
}

#[derive(Clone, Copy, PartialEq, Debug)]
enum Family {
    A,
    B,
    C,
}

fn grid_frame(
    ctx: &egui::Context,
    width: f32,
    selected: &mut Family,
    enabled: bool,
    events: Vec<egui::Event>,
) -> (bool, Vec<(String, Pos2)>) {
    let mut changed = false;
    let output = ctx.run_ui(input(width, events), |ui| {
        let art = |w: f32| Diagram::new().rect(0., 0., w, 1.);
        changed = choice_grid(
            ui,
            "grid",
            selected,
            &[
                Choice {
                    value: Family::A,
                    label: "Plain washer",
                    hint: "a",
                    art: art(1.),
                },
                Choice {
                    value: Family::B,
                    label: "Hex nut",
                    hint: "b",
                    art: art(2.),
                },
                Choice {
                    value: Family::C,
                    label: "Bearing mounting envelope",
                    hint: "c",
                    art: art(3.),
                },
            ],
            enabled,
        );
    });
    let labels = output
        .shapes
        .iter()
        .filter_map(|s| match &s.shape {
            egui::Shape::Text(t) => {
                Some((t.galley.text().to_owned(), t.pos + t.galley.size() / 2.))
            }
            _ => None,
        })
        .collect();
    output.drop_without_applying_deltas();
    (changed, labels)
}

#[test]
fn choice_grid_selects_by_pointer_and_keyboard_and_respects_disabled_state() {
    for width in [160., 320.] {
        let ctx = egui::Context::default();
        let mut selected = Family::A;
        grid_frame(&ctx, width, &mut selected, true, vec![]);
        let (_, labels) = grid_frame(&ctx, width, &mut selected, true, vec![]);
        for (label, pos) in &labels {
            assert!(
                pos.x > 0. && pos.x < width,
                "{label} outside the grid at {width}"
            );
        }
        let at = |name: &str| labels.iter().find(|(l, _)| l == name).unwrap().1;
        let (changed, _) = grid_frame(&ctx, width, &mut selected, true, click(at("Hex nut")));
        assert!(changed);
        assert_eq!(selected, Family::B);
        // The clicked tile has focus; the arrow key moves the selection on.
        let key = |key| {
            vec![egui::Event::Key {
                key,
                physical_key: None,
                pressed: true,
                repeat: false,
                modifiers: Default::default(),
            }]
        };
        grid_frame(&ctx, width, &mut selected, true, vec![]);
        let (changed, _) = grid_frame(&ctx, width, &mut selected, true, key(egui::Key::ArrowRight));
        assert!(changed);
        assert_eq!(selected, Family::C);
        grid_frame(&ctx, width, &mut selected, true, vec![]);
        grid_frame(&ctx, width, &mut selected, true, key(egui::Key::ArrowLeft));
        assert_eq!(selected, Family::B);
        // A disabled selector (editing an existing part) never changes.
        let (changed, _) = grid_frame(&ctx, width, &mut selected, false, click(at("Plain washer")));
        assert!(!changed);
        assert_eq!(selected, Family::B);
    }
}

#[test]
fn sections_start_closed_reveal_on_request_and_keep_their_values() {
    let ctx = egui::Context::default();
    let mut field = mm(12.);
    let frame = |reveal: bool, field: &mut DimensionField| {
        let mut shown = false;
        ctx.run_ui(input(300., vec![]), |ui| {
            shown = section(ui, "advanced", "Advanced placement", reveal, |ui| {
                field.show(ui, "Origin X", "help");
            })
            .is_some();
        })
        .drop_without_applying_deltas();
        shown
    };
    assert!(!frame(false, &mut field));
    assert!(!frame(false, &mut field), "closed by default");
    assert_eq!(field.text, "12", "collapsed fields keep their values");
    assert!(frame(true, &mut field));
    assert!(frame(true, &mut field));
    assert!(frame(false, &mut field), "a revealed group stays open");
}
