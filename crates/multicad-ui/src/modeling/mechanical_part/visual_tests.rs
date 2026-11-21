//! CAD-386 observable behaviour of the visual hardware form.
use super::*;
use crate::modeling::diagram::{Callout, Mark};
use multicad_core::commands::Session;

fn session() -> Session {
    Session::new(Document::new("Visual hardware").unwrap()).unwrap()
}

fn washer_request() -> MechanicalPartRequest {
    MechanicalPartRequest::new(
        "Washer",
        MechanicalPartParameters::Washer {
            bore_diameter: length(3.2),
            outside_diameter: length(7.),
            thickness: length(0.5),
            frame: Default::default(),
        },
    )
}

fn panel(s: &Session, request: MechanicalPartRequest, editing: bool) -> MechanicalPartPanel {
    MechanicalPartPanel::new(request, s.parameters().clone(), editing).unwrap()
}

struct Frame {
    texts: Vec<(String, egui::Rect)>,
    event: Option<ModelEvent>,
}

fn frame(
    ctx: &egui::Context,
    p: &mut MechanicalPartPanel,
    s: &Session,
    width: f32,
    events: Vec<egui::Event>,
) -> Frame {
    let mut event = None;
    let output = ctx.run_ui(
        egui::RawInput {
            screen_rect: Some(egui::Rect::from_min_size(
                egui::Pos2::ZERO,
                egui::vec2(width, 3000.),
            )),
            focused: true,
            events,
            ..Default::default()
        },
        |ui| {
            event = p.show(ui, s.document(), true, PreviewState::Idle);
        },
    );
    let texts = output
        .shapes
        .iter()
        .filter_map(|shape| match &shape.shape {
            egui::Shape::Text(t) => Some((
                t.galley.text().to_owned(),
                egui::Rect::from_min_size(t.pos, t.galley.size()),
            )),
            _ => None,
        })
        .collect();
    output.drop_without_applying_deltas();
    Frame { texts, event }
}

fn click_at(pos: egui::Pos2) -> Vec<egui::Event> {
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

fn find<'a>(f: &'a Frame, text: &str) -> Option<&'a egui::Rect> {
    f.texts.iter().find(|(t, _)| t == text).map(|(_, r)| r)
}

fn type_replace(text: &str) -> Vec<egui::Event> {
    vec![
        egui::Event::Key {
            key: egui::Key::A,
            physical_key: None,
            pressed: true,
            repeat: false,
            modifiers: egui::Modifiers::COMMAND,
        },
        egui::Event::Text(text.into()),
    ]
}

fn panel_of(s: &Session, kind: Kind) -> MechanicalPartPanel {
    let mut p = panel(s, washer_request(), false);
    p.change_kind(kind);
    p
}

#[test]
fn every_family_has_its_own_diagram_whose_callouts_name_only_its_visible_fields() {
    let s = session();
    let mut seen: Vec<diagram::Diagram> = vec![];
    for kind in Kind::ALL {
        let p = panel_of(&s, kind);
        let schematic = p.schematic(s.document());
        let fields: Vec<Key> = p.linked_fields().into_iter().map(|(k, _)| k).collect();
        assert!(!schematic.callouts.is_empty(), "{kind:?} has dimensions");
        for key in schematic.keys() {
            assert!(fields.contains(&key), "{kind:?} callout {key} has no field");
        }
        for other in &seen {
            assert_ne!(
                other, &schematic,
                "{kind:?} reuses another family's diagram"
            );
        }
        // Thumbnails depict the family too, not one shared picture.
        assert!(visual::art(kind).callouts.is_empty());
        seen.push(schematic);
    }
    let notation = |kind, key| {
        panel_of(&s, kind)
            .schematic(s.document())
            .callouts
            .into_iter()
            .find(|c: &Callout| c.key == key)
            .unwrap()
            .notation
    };
    use crate::modeling::diagram::Notation;
    assert_eq!(notation(Kind::Washer, "washer.bore"), Notation::Diameter);
    assert_eq!(notation(Kind::Washer, "washer.thickness"), Notation::Length);
    assert_eq!(notation(Kind::HexNut, visual::NUT_FLATS), Notation::Length);
    assert_eq!(
        notation(Kind::HexNut, visual::THREAD_DIAMETER),
        Notation::Diameter
    );
    assert_eq!(
        notation(Kind::MachineScrew, visual::SCREW_SHANK),
        Notation::Diameter
    );
}

fn length_callout(p: &MechanicalPartPanel, s: &Session) -> Callout {
    p.schematic(s.document())
        .callouts
        .into_iter()
        .find(|c| c.key == visual::SCREW_LENGTH)
        .unwrap()
}

fn head_top(p: &MechanicalPartPanel, s: &Session) -> f32 {
    p.schematic(s.document())
        .callouts
        .into_iter()
        .find(|c| c.key == visual::HEAD_HEIGHT)
        .map(|c| c.from[1].max(c.to[1]))
        .unwrap()
}

#[test]
fn screw_length_endpoints_follow_the_convention_and_the_actual_head() {
    let s = session();
    for head in [
        screw_shape::HeadKind::SocketCap,
        screw_shape::HeadKind::Countersunk,
    ] {
        let mut p = panel_of(&s, Kind::MachineScrew);
        p.screw_shape.head = head;
        p.screw_length_convention = model::ScrewLengthConvention::Underhead;
        let under = length_callout(&p, &s);
        assert_eq!(
            under.from[1], 0.,
            "underhead length starts at the neck datum"
        );
        assert_eq!(under.name, "Length under head");
        p.screw_length_convention = model::ScrewLengthConvention::Overall;
        let overall = length_callout(&p, &s);
        assert_eq!(overall.name, "Overall length");
        assert!(
            (overall.from[1] - head_top(&p, &s)).abs() < 1e-4,
            "{head:?} overall length includes the actual head height"
        );
        // The same typed 10 mm ends at the same tip only when the convention
        // says so: overall moves the tip up by the head height.
        assert!(overall.to[1] > under.to[1]);
    }
    let mut p = panel_of(&s, Kind::MachineScrew);
    p.screw_shape.head = screw_shape::HeadKind::Countersunk;
    let keys = p.schematic(s.document()).keys();
    assert!(
        keys.contains(&visual::HEAD_ANGLE),
        "countersunk angle is drawn"
    );
    p.screw_shape.head = screw_shape::HeadKind::ExternalHex;
    let keys = p.schematic(s.document()).keys();
    assert!(keys.contains(&visual::HEAD_FLATS) && !keys.contains(&visual::HEAD_DIAMETER));
}

#[test]
fn clicking_a_washer_dimension_edits_the_canonical_field_and_the_list_updates_the_diagram() {
    let s = session();
    for width in [180., 420.] {
        let ctx = egui::Context::default();
        let mut p = panel(&s, washer_request(), false);
        frame(&ctx, &mut p, &s, width, vec![]);
        let f = frame(&ctx, &mut p, &s, width, vec![]);
        let label = *find(&f, "Ø Outside 7 mm").expect("outside callout painted");
        assert!(
            label.right() <= width && label.left() >= 0.,
            "callout unclipped at {width}"
        );
        for (name, rect) in &f.texts {
            if [
                "Plain washer",
                "Hex nut",
                "Machine screw",
                "Wood screw envelope",
                "Bearing mounting envelope",
            ]
            .contains(&name.as_str())
            {
                assert!(
                    rect.left() >= 0. && rect.right() <= width + 0.5,
                    "{name} tile label clipped at {width}"
                );
            }
        }
        frame(&ctx, &mut p, &s, width, click_at(label.center()));
        frame(&ctx, &mut p, &s, width, vec![]);
        let f = frame(&ctx, &mut p, &s, width, type_replace("9"));
        assert!(
            matches!(f.event, Some(ModelEvent::Changed)),
            "an edit restarts the automatic preview"
        );
        let MechanicalPartParameters::Washer {
            outside_diameter,
            bore_diameter,
            ..
        } = p.request(s.document(), true).unwrap().parameters
        else {
            unreachable!()
        };
        assert_eq!(outside_diameter, length(9.));
        assert_eq!(bore_diameter, length(3.2));
        let f = frame(&ctx, &mut p, &s, width, vec![]);
        assert!(
            find(&f, "Ø Outside 9 mm").is_some(),
            "diagram shows the staged value"
        );
        // Editing the list field directly updates the same diagram label.
        p.fields[0].text = "=2 + 2".into();
        let f = frame(&ctx, &mut p, &s, width, vec![]);
        assert!(find(&f, "Ø Bore =2 + 2").is_some());
    }
}

#[test]
fn placement_is_collapsed_keeps_values_and_validation_reveals_and_marks_the_field() {
    let s = session();
    let ctx = egui::Context::default();
    let mut request = washer_request();
    if let MechanicalPartParameters::Washer { frame, .. } = &mut request.parameters {
        frame.origin[0] = length(12.);
    }
    let mut p = panel(&s, request.clone(), false);
    frame(&ctx, &mut p, &s, 320., vec![]);
    let f = frame(&ctx, &mut p, &s, 320., vec![]);
    assert!(
        find(&f, "Base center X (mm)").is_none(),
        "advanced placement starts collapsed"
    );
    assert!(
        f.texts
            .iter()
            .any(|(t, _)| t.starts_with("Placement") && t.contains("12 mm")),
        "collapsed summary keeps the authored value visible"
    );
    assert_eq!(
        p.request(s.document(), true).unwrap(),
        request,
        "collapsed values are retained"
    );
    p.origin[0].text = "=".into();
    frame(&ctx, &mut p, &s, 320., vec![]);
    let f = frame(&ctx, &mut p, &s, 320., vec![]);
    assert_eq!(p.link.invalid, Some("hardware.x"));
    assert!(
        find(&f, "Base center X (mm)").is_some(),
        "validation reveals the containing section"
    );
    assert!(p.request(s.document(), true).is_err());
    p.origin[0].text = "12".into();
    assert_eq!(p.request(s.document(), true).unwrap(), request);
    let f = frame(&ctx, &mut p, &s, 320., vec![]);
    assert_eq!(p.link.invalid, None);
    assert!(
        find(&f, "Base center X (mm)").is_some(),
        "a revealed section stays open"
    );
}

#[test]
fn family_tiles_switch_families_and_are_fixed_when_editing() {
    let s = session();
    let ctx = egui::Context::default();
    let mut p = panel(&s, washer_request(), false);
    frame(&ctx, &mut p, &s, 320., vec![]);
    let f = frame(&ctx, &mut p, &s, 320., vec![]);
    let nut = find(&f, "Hex nut").unwrap().center();
    frame(&ctx, &mut p, &s, 320., click_at(nut));
    assert!(p.kind == Kind::HexNut);
    let f = frame(&ctx, &mut p, &s, 320., vec![]);
    assert!(
        f.texts.iter().any(|(t, _)| t.starts_with("Across flats ")),
        "nut diagram follows the selection"
    );
    assert!(f.texts.iter().any(|(t, _)| t.starts_with("Height ")));
    let mut editing = panel(&s, washer_request(), true);
    frame(&ctx, &mut editing, &s, 320., vec![]);
    let f = frame(&ctx, &mut editing, &s, 320., vec![]);
    let nut = find(&f, "Hex nut").unwrap().center();
    frame(&ctx, &mut editing, &s, 320., click_at(nut));
    assert!(
        editing.kind == Kind::Washer,
        "an existing part keeps its family"
    );
}

#[test]
fn diagram_marks_follow_current_values_but_never_reach_the_request() {
    let s = session();
    let mut p = panel(&s, washer_request(), false);
    let before = p.schematic(s.document());
    let request = p.request(s.document(), true).unwrap();
    p.fields[1].text = "12".into();
    let after = p.schematic(s.document());
    let radius = |d: &diagram::Diagram| match &d.marks[0] {
        Mark::Circle { radius, .. } => *radius,
        _ => unreachable!(),
    };
    assert!(radius(&after) > radius(&before));
    p.fields[1].text = "7".into();
    assert_eq!(p.request(s.document(), true).unwrap(), request);
}
