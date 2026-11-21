//! CAD-394 thread diagrams: variant-specific profiles linked to the recipe.
use super::*;
use crate::modeling::diagram::Mark;
use multicad_core::parameters;

fn request(mode: ThreadMode) -> ThreadRequest {
    let mut request = ThreadRequest::new(
        "Diagram thread",
        ThreadParameters {
            major_diameter: model::length(6.),
            pitch: model::length(1.),
            length: model::length(10.),
            ..Default::default()
        },
    );
    request.mode = mode;
    request
}

fn panel(document: &Document, mode: ThreadMode) -> ThreadPanel {
    ThreadPanel::new(
        request(mode),
        parameters::evaluate(document).unwrap(),
        vec![],
        &Default::default(),
        true,
    )
    .unwrap()
}

#[test]
fn every_family_and_mode_draws_its_own_profile_with_linked_dimensions() {
    let document = Document::new("Thread families").unwrap();
    let mut seen = vec![];
    for family in ThreadFamily::ALL {
        for mode in [ThreadMode::NewExternal, ThreadMode::CutInternal] {
            let mut p = panel(&document, mode);
            p.family = family;
            let schematic = visual::diagram(&p.variant(mode), &p.diagram_values(&document));
            let fields: Vec<Key> = p.linked_fields().into_iter().map(|(k, _)| k).collect();
            for key in schematic.keys() {
                assert!(
                    fields.contains(&key),
                    "{family:?} callout {key} has no field"
                );
            }
            for key in [visual::DIAMETER, visual::PITCH, visual::LENGTH] {
                assert!(schematic.keys().contains(&key), "{family:?} lacks {key}");
            }
            assert!(
                schematic.marks.iter().any(
                    |m| matches!(m, Mark::Note { text, .. } if text.starts_with(family.label()))
                ),
                "{family:?} diagram names the selected family"
            );
            assert!(
                !seen.contains(&schematic),
                "{family:?} {mode:?} reuses another diagram"
            );
            seen.push(schematic);
        }
    }
    // Profile forms are the actual family geometry, not one generic V.
    let f = |family| visual::form(family, None, 1.);
    assert_eq!(f(ThreadFamily::IsoMetric).lower, 30.);
    assert_eq!(f(ThreadFamily::Bsw).lower, 27.5);
    assert_eq!(f(ThreadFamily::Acme).lower, 14.5);
    assert_eq!(f(ThreadFamily::Trapezoidal).lower, 15.);
    let buttress = f(ThreadFamily::Buttress);
    assert_ne!(
        buttress.lower, buttress.upper,
        "buttress flanks are asymmetric"
    );
    assert!(f(ThreadFamily::Npt).taper > 0. && f(ThreadFamily::Bspt).taper > 0.);
    assert_eq!(f(ThreadFamily::Bspp).taper, 0.);
    // Custom flanks follow the authored angles and add their callouts.
    let mut p = panel(&document, ThreadMode::NewExternal);
    p.family = ThreadFamily::Custom;
    p.custom[2].text = "10".into();
    p.custom[3].text = "40".into();
    let custom = visual::form(ThreadFamily::Custom, p.diagram_values(&document).custom, 1.);
    assert_eq!((custom.lower, custom.upper), (10., 40.));
    let keys = visual::diagram(
        &p.variant(ThreadMode::NewExternal),
        &p.diagram_values(&document),
    )
    .keys();
    assert!(keys.contains(&visual::CUSTOM[2]) && keys.contains(&visual::CUSTOM[3]));
}

#[test]
fn end_treatments_and_hand_change_the_drawing() {
    let document = Document::new("Thread treatments").unwrap();
    let mut p = panel(&document, ThreadMode::NewExternal);
    let plain = visual::diagram(
        &p.variant(ThreadMode::NewExternal),
        &p.diagram_values(&document),
    );
    p.end_treatments_enabled = true;
    p.start_chamfer = true;
    p.end_relief = true;
    let treated = visual::diagram(
        &p.variant(ThreadMode::NewExternal),
        &p.diagram_values(&document),
    );
    for key in [
        visual::START_DEPTH,
        visual::RELIEF_LENGTH,
        visual::RELIEF_DEPTH,
    ] {
        assert!(treated.keys().contains(&key) && !plain.keys().contains(&key));
        assert!(p.linked_fields().iter().any(|(k, _)| *k == key));
    }
    p.hand = ThreadHand::Left;
    let left = visual::diagram(
        &p.variant(ThreadMode::NewExternal),
        &p.diagram_values(&document),
    );
    assert_ne!(left, treated);
}

struct Painted {
    texts: Vec<(String, egui::Rect)>,
}

fn frame(
    ctx: &egui::Context,
    p: &mut ThreadPanel,
    d: &Document,
    width: f32,
    events: Vec<egui::Event>,
) -> Painted {
    let output = ctx.run_ui(
        egui::RawInput {
            screen_rect: Some(egui::Rect::from_min_size(
                egui::Pos2::ZERO,
                egui::vec2(width, 3200.),
            )),
            focused: true,
            events,
            ..Default::default()
        },
        |ui| {
            let _ = p.show(ui, d, true, PreviewState::Idle);
        },
    );
    let texts = output
        .shapes
        .iter()
        .filter_map(|s| match &s.shape {
            egui::Shape::Text(t) => Some((
                t.galley.text().to_owned(),
                egui::Rect::from_min_size(t.pos, t.galley.size()),
            )),
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

#[test]
fn clicking_the_pitch_callout_edits_the_recipe_pitch_at_ordinary_and_compact_widths() {
    let document = Document::new("Painted thread").unwrap();
    for width in [200., 420.] {
        let ctx = egui::Context::default();
        let mut p = panel(&document, ThreadMode::NewExternal);
        frame(&ctx, &mut p, &document, width, vec![]);
        let painted = frame(&ctx, &mut p, &document, width, vec![]);
        let pitch = painted
            .texts
            .iter()
            .find(|(t, _)| t == "Pitch 1 mm")
            .map(|(_, r)| *r)
            .expect("pitch callout");
        assert!(
            pitch.left() >= 0. && pitch.right() <= width,
            "unclipped at {width}"
        );
        frame(&ctx, &mut p, &document, width, click(pitch.center()));
        frame(&ctx, &mut p, &document, width, vec![]);
        frame(
            &ctx,
            &mut p,
            &document,
            width,
            vec![
                egui::Event::Key {
                    key: egui::Key::A,
                    physical_key: None,
                    pressed: true,
                    repeat: false,
                    modifiers: egui::Modifiers::COMMAND,
                },
                egui::Event::Text("1.25".into()),
            ],
        );
        // The callout edited the recipe's own pitch input, nothing else.
        assert_eq!(p.pitch.resolve(&document).unwrap(), model::length(1.25));
        assert_eq!(p.diameter.resolve(&document).unwrap(), model::length(6.));
        assert_eq!(p.length.resolve(&document).unwrap(), model::length(10.));
        let painted = frame(&ctx, &mut p, &document, width, vec![]);
        assert!(painted.texts.iter().any(|(t, _)| t == "Pitch 1.25 mm"));
    }
}

#[test]
fn fabrication_allowance_is_collapsed_retained_and_revealed_by_validation() {
    let document = Document::new("Allowance").unwrap();
    let ctx = egui::Context::default();
    let mut p = panel(&document, ThreadMode::NewExternal);
    p.radial.set_literal(0.05, "mm");
    let recipe = p.recipe(&document, ThreadMode::NewExternal).unwrap();
    frame(&ctx, &mut p, &document, 360., vec![]);
    let painted = frame(&ctx, &mut p, &document, 360., vec![]);
    let has = |painted: &Painted, text: &str| painted.texts.iter().any(|(t, _)| t == text);
    assert!(
        !has(&painted, "Radial clearance (mm)"),
        "advanced allowance starts collapsed"
    );
    assert!(
        painted
            .texts
            .iter()
            .any(|(t, _)| t.starts_with("Fabrication allowance") && t.contains("0.05 mm"))
    );
    assert!(
        !has(&painted, "Base center X (mm)"),
        "axis placement starts collapsed"
    );
    assert_eq!(
        p.recipe(&document, ThreadMode::NewExternal).unwrap(),
        recipe
    );
    p.radial.text = "=".into();
    for _ in 0..10 {
        frame(&ctx, &mut p, &document, 360., vec![]);
    }
    let painted = frame(&ctx, &mut p, &document, 360., vec![]);
    assert_eq!(p.link.invalid, Some(visual::ALLOWANCE[0]));
    assert!(
        has(&painted, "Radial clearance (mm)"),
        "validation reveals its section"
    );
    p.origin[2].text = "=".into();
    p.radial.set_literal(0.05, "mm");
    // A revealed group animates open over a few frames.
    for _ in 0..10 {
        frame(&ctx, &mut p, &document, 360., vec![]);
    }
    let painted = frame(&ctx, &mut p, &document, 360., vec![]);
    assert_eq!(p.link.invalid, Some(visual::ORIGIN[2]));
    assert!(has(&painted, "Base center Z (mm)"));
}
