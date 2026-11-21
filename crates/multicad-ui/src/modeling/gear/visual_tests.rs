//! CAD-394 gear diagrams: every family is drawn as itself and its callouts
//! edit the same staged recipe fields.
use super::*;
use multicad_core::commands::Session;

fn session() -> Session {
    Session::new(Document::new("Gear diagrams").unwrap()).unwrap()
}

fn panel(s: &Session) -> GearPanel {
    GearPanel::new(
        GearRequest::new("Diagram gear", multicad_core::model::GearParameters::default()),
        s.parameters().clone(),
        false,
    )
    .unwrap()
}

fn schematic(p: &GearPanel, paired: bool) -> diagram::Diagram {
    p.schematic(paired)
}

#[test]
fn every_supported_family_draws_its_own_diagram_with_only_visible_callouts() {
    let s = session();
    let mut seen = vec![];
    for family in GearFamily::ALL.into_iter().filter(|f| f.supported()) {
        let mut p = panel(&s);
        p.select_family(family);
        let d = schematic(&p, false);
        let fields: Vec<Key> = p.linked_fields(false).into_iter().map(|(k, _)| k).collect();
        assert!(d.callouts.len() >= 3, "{family:?} dimensions its recipe");
        for key in d.keys() {
            assert!(fields.contains(&key));
        }
        assert!(d.keys().contains(&visual::FIELDS[0]), "{family:?} names its pitch input");
        assert!(!seen.contains(&d), "{family:?} reuses another family's drawing");
        seen.push(d);
        // A pair hides the per-member mate inputs, and their callouts.
        let paired = schematic(&p, true);
        for key in [visual::WORM_MATE, visual::BEVEL_MATE] {
            assert!(!paired.keys().contains(&key));
        }
    }
    let mut p = panel(&s);
    p.select_family(GearFamily::Helical);
    let right = schematic(&p, false);
    assert!(right.keys().contains(&visual::FIELDS[4]), "helix angle is dimensioned");
    p.hand = GearHand::Left;
    assert_ne!(schematic(&p, false), right, "hand changes the tooth traces");
    p.select_family(GearFamily::DoubleHelical);
    assert!(schematic(&p, false).keys().contains(&visual::FIELDS[5]), "center gap");
    p.select_family(GearFamily::Spur);
    let spur = schematic(&p, false);
    assert!(!spur.keys().contains(&visual::FIELDS[4]), "spur has no helix");
    p.select_pitch_input(GearPitchInput::DiametralPitch);
    let dp = schematic(&p, false);
    let label = |d: &diagram::Diagram| {
        d.callouts
            .iter()
            .find(|c| c.key == visual::FIELDS[0])
            .unwrap()
            .name
            .clone()
    };
    assert_eq!(label(&spur), "Module");
    assert_eq!(label(&dp), "Diametral pitch");
}

struct Painted {
    texts: Vec<(String, egui::Rect)>,
}
fn frame(ctx: &egui::Context, p: &mut GearPanel, s: &Session, width: f32, events: Vec<egui::Event>) -> Painted {
    let output = ctx.run_ui(
        egui::RawInput {
            screen_rect: Some(egui::Rect::from_min_size(egui::Pos2::ZERO, egui::vec2(width, 3600.))),
            focused: true,
            events,
            ..Default::default()
        },
        |ui| {
            let _ = p.show(ui, s.document(), true, PreviewState::Idle);
        },
    );
    let texts = output
        .shapes
        .iter()
        .filter_map(|shape| match &shape.shape {
            egui::Shape::Text(t) => Some((t.galley.text().to_owned(), egui::Rect::from_min_size(t.pos, t.galley.size()))),
            _ => None,
        })
        .collect();
    output.drop_without_applying_deltas();
    Painted { texts }
}
fn click(pos: egui::Pos2) -> Vec<egui::Event> {
    [true, false]
        .into_iter()
        .flat_map(|pressed| {
            [
                egui::Event::PointerMoved(pos),
                egui::Event::PointerButton {
                    pos,
                    button: egui::PointerButton::Primary,
                    pressed,
                    modifiers: Default::default(),
                },
            ]
        })
        .collect()
}
fn find(p: &Painted, prefix: &str) -> Option<egui::Rect> {
    p.texts.iter().find(|(t, _)| t.starts_with(prefix)).map(|(_, r)| *r)
}

#[test]
fn family_tiles_and_face_width_callout_edit_the_staged_recipe() {
    let s = session();
    for width in [220., 460.] {
        let ctx = egui::Context::default();
        let mut p = panel(&s);
        frame(&ctx, &mut p, &s, width, vec![]);
        let painted = frame(&ctx, &mut p, &s, width, vec![]);
        let helical = find(&painted, "Helical").expect("helical tile");
        frame(&ctx, &mut p, &s, width, click(helical.center()));
        assert_eq!(p.family, GearFamily::Helical);
        frame(&ctx, &mut p, &s, width, vec![]);
        let painted = frame(&ctx, &mut p, &s, width, vec![]);
        assert!(find(&painted, "∠ Helix").is_some(), "diagram follows the family");
        let face = find(&painted, "Face width ").expect("face callout");
        assert!(face.left() >= 0. && face.right() <= width, "unclipped at {width}");
        frame(&ctx, &mut p, &s, width, click(face.center()));
        frame(&ctx, &mut p, &s, width, vec![]);
        frame(
            &ctx,
            &mut p,
            &s,
            width,
            vec![
                egui::Event::Key {
                    key: egui::Key::A,
                    physical_key: None,
                    pressed: true,
                    repeat: false,
                    modifiers: egui::Modifiers::COMMAND,
                },
                egui::Event::Text("13".into()),
            ],
        );
        assert_eq!(p.fields[3].text, "13");
        let painted = frame(&ctx, &mut p, &s, width, vec![]);
        assert!(find(&painted, "Face width 13 mm").is_some());
    }
}

#[test]
fn tolerances_are_collapsed_retained_and_revealed_by_validation() {
    let s = session();
    let ctx = egui::Context::default();
    let mut p = panel(&s);
    p.fields[7].set_literal(0.04, "mm");
    let request = p.request(s.document(), true).unwrap();
    frame(&ctx, &mut p, &s, 360., vec![]);
    let painted = frame(&ctx, &mut p, &s, 360., vec![]);
    assert!(find(&painted, "Design backlash per gear").is_none());
    assert!(find(&painted, "Backlash, clearance and root · backlash 0.04 mm").is_some());
    assert_eq!(p.request(s.document(), true).unwrap(), request);
    p.fields[9].text = "=".into();
    // A revealed group animates open over a few frames.
    for _ in 0..10 {
        frame(&ctx, &mut p, &s, 360., vec![]);
    }
    let painted = frame(&ctx, &mut p, &s, 360., vec![]);
    assert_eq!(p.link.invalid, Some(visual::FIELDS[9]));
    assert!(find(&painted, "Root fillet / module").is_some());
}
