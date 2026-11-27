//! CAD-395 wood box schematic: the envelope shows present and omitted
//! faces; callouts edit the staged outside size and placement.
use super::*;
use crate::modeling::diagram::{
    Mark,
    harness::{self, WIDTHS},
};

fn show<'a>(p: &'a mut WoodBoxPanel, d: &'a Document) -> impl FnMut(&mut egui::Ui) + 'a {
    move |ui| {
        let _ = p.show(ui, d, true, PreviewState::Idle);
    }
}

fn open_notes(diagram: &super::super::diagram::Diagram) -> Vec<String> {
    diagram
        .marks
        .iter()
        .filter_map(|m| match m {
            Mark::Note { text, .. } if text.ends_with(" open") => Some(text.clone()),
            _ => None,
        })
        .collect()
}

#[test]
fn the_envelope_marks_omitted_faces_and_names_only_staged_fields() {
    let (session, _) = super::super::panel_fixture::flat_panel();
    let d = session.document();
    let mut p = WoodBoxPanel::new(d, "Drawer".into()).unwrap();
    let open_top = visual::diagram(&p.diagram_values(d));
    assert_eq!(open_notes(&open_top), ["Top open"]);
    let fields: Vec<_> = p.linked_fields().into_iter().map(|(k, _)| k).collect();
    for key in open_top.keys() {
        assert!(fields.contains(&key), "{key} edits a staged field");
    }
    for key in visual::SIZE {
        assert!(open_top.keys().contains(&key));
    }
    assert!(
        open_top.marks.iter().any(
            |m| matches!(m, Mark::Note { text, .. } if text.contains("5 of 6 faces") && text.contains("mm panels"))
        ),
        "the stock thickness is noted"
    );
    p.source.definition.faces = WoodBoxFace::ALL.into_iter().collect();
    let closed = visual::diagram(&p.diagram_values(d));
    assert!(open_notes(&closed).is_empty());
    p.source.definition.faces.remove(&WoodBoxFace::Left);
    assert_eq!(
        open_notes(&visual::diagram(&p.diagram_values(d))),
        ["Left open"]
    );
}

#[test]
fn clicking_the_depth_callout_edits_the_outside_depth_at_both_widths() {
    let (session, _) = super::super::panel_fixture::flat_panel();
    let d = session.document();
    for width in WIDTHS {
        let ctx = egui::Context::default();
        let mut p = WoodBoxPanel::new(d, "Drawer".into()).unwrap();
        let painted =
            harness::edit_through_callout(&ctx, width, "Depth 400 mm", "450", &mut show(&mut p, d));
        assert!(painted.has("Depth 450 mm"));
        let WoodBoxChange::Set(request) = p.request(d, true).unwrap() else {
            panic!("set")
        };
        assert_eq!(request.definition.depth, length(450.));
        assert_eq!(request.definition.width, length(600.), "only depth changed");
        let schematic = visual::diagram(&p.diagram_values(d));
        harness::assert_callouts_apart(&schematic, width, 200., &|k| p.callout_value(k));
    }
}

#[test]
fn placement_is_collapsed_retained_and_revealed() {
    let (session, _) = super::super::panel_fixture::flat_panel();
    let d = session.document();
    let ctx = egui::Context::default();
    let mut p = WoodBoxPanel::new(d, "Drawer".into()).unwrap();
    p.origin[0].text = "100".into();
    p.rotation.text = "15".into();
    let painted = harness::settle(&ctx, 360., &mut show(&mut p, d));
    assert!(!painted.has("X"), "placement starts collapsed");
    assert!(painted.has_prefix("Placement · 100, 0, 0 mm · 15°"));
    let WoodBoxChange::Set(request) = p.request(d, true).unwrap() else {
        panic!("set")
    };
    assert_eq!(request.definition.origin[0], length(100.));
    let origin = painted.find("Origin 100, 0, 0 mm").expect("origin callout");
    harness::frame(
        &ctx,
        360.,
        harness::click(origin.center()),
        &mut show(&mut p, d),
    );
    let painted = harness::settle(&ctx, 360., &mut show(&mut p, d));
    assert!(painted.has("X"), "the origin callout opens placement");
    assert_eq!(p.link.active(), Some(visual::ORIGIN[0]));

    let ctx = egui::Context::default();
    let mut p = WoodBoxPanel::new(d, "Drawer".into()).unwrap();
    p.rotation.text = "=".into();
    let painted = harness::settle(&ctx, 360., &mut show(&mut p, d));
    assert_eq!(p.link.invalid, Some(visual::ROTATION));
    assert!(painted.has("Rotation"), "validation reveals placement");
}
