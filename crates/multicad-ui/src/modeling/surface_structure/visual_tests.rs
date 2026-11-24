//! CAD-395 surface relief schematic: the pattern tiles, the relief window
//! and the pattern section are linked to the form's staged fields.
use super::*;
use crate::modeling::diagram::harness::{self, WIDTHS};
use multicad_core::model::length;

fn show<'a>(p: &'a mut SurfaceStructurePanel, d: &'a Document) -> impl FnMut(&mut egui::Ui) + 'a {
    move |ui| {
        let _ = p.show(ui, d, true, PreviewState::Idle);
    }
}

fn document() -> Document {
    Document::new("Relief diagrams").unwrap()
}

#[test]
fn every_pattern_draws_its_own_relief_and_callouts_name_staged_fields() {
    let d = document();
    let mut p = SurfaceStructurePanel::new().unwrap();
    let mut seen = vec![];
    for pattern in [Pattern::Ribs, Pattern::Waves, Pattern::Stipple] {
        p.pattern = pattern;
        let schematic = visual::diagram(p.kind(), &p.diagram_values(&d), p.image());
        let fields: Vec<Key> = p.linked_fields().into_iter().map(|(k, _)| k).collect();
        for key in schematic.keys() {
            assert!(fields.contains(&key), "{key} edits a staged field");
        }
        assert!(schematic.keys().contains(&visual::PITCH));
        assert!(
            !seen.contains(&schematic),
            "{:?} reuses a drawing",
            p.kind()
        );
        seen.push(schematic);
    }
    // A recessed relief flips the section and names the depth.
    let raised = visual::diagram(p.kind(), &p.diagram_values(&d), None);
    p.direction = Direction::Recessed;
    assert_ne!(
        visual::diagram(p.kind(), &p.diagram_values(&d), None),
        raised
    );
    assert_eq!(p.callout_value(visual::OFFSET[0]), "0, 0 mm");
    // A height map has no pitch callout or field.
    p.pattern = Pattern::HeightMap {
        width: 2,
        height: 2,
        samples: vec![0, 255, 255, 0],
        repeat: [1, 1],
    };
    let map = visual::diagram(p.kind(), &p.diagram_values(&d), p.image());
    assert!(!map.keys().contains(&visual::PITCH));
    assert!(!p.linked_fields().iter().any(|(k, _)| *k == visual::PITCH));
}

#[test]
fn clicking_the_pitch_callout_edits_the_pattern_pitch_at_both_widths() {
    let d = document();
    for width in WIDTHS {
        let ctx = egui::Context::default();
        let mut p = SurfaceStructurePanel::new().unwrap();
        let painted =
            harness::edit_through_callout(&ctx, width, "Pitch 4 mm", "6", &mut show(&mut p, &d));
        assert!(painted.has("Pitch 6 mm"));
        assert_eq!(p.pitch.resolve(&d).unwrap(), length(6.));
        assert_eq!(p.amplitude.resolve(&d).unwrap(), length(0.8), "only pitch");
        for pattern in [Pattern::Ribs, Pattern::Waves, Pattern::Stipple] {
            p.pattern = pattern;
            let schematic = visual::diagram(p.kind(), &p.diagram_values(&d), None);
            harness::assert_callouts_apart(&schematic, width, 240., &|k| p.callout_value(k));
        }
    }
}

#[test]
fn pattern_tiles_choose_the_pattern() {
    let d = document();
    let ctx = egui::Context::default();
    let mut p = SurfaceStructurePanel::new().unwrap();
    let painted = harness::settle(&ctx, 360., &mut show(&mut p, &d));
    let tile = painted
        .texts
        .iter()
        .find(|(t, _)| t == "Waves")
        .map(|(_, r)| *r)
        .expect("waves tile");
    let before = p.stamp();
    harness::frame(
        &ctx,
        360.,
        harness::click(tile.center()),
        &mut show(&mut p, &d),
    );
    assert_eq!(p.pattern, Pattern::Waves);
    assert_ne!(p.stamp(), before, "choosing a tile stages a change");
}

#[test]
fn window_placement_is_collapsed_retained_and_revealed() {
    let d = document();
    let ctx = egui::Context::default();
    let mut p = SurfaceStructurePanel::new().unwrap();
    p.offset[0].text = "3".into();
    p.rotation.text = "30".into();
    let painted = harness::settle(&ctx, 360., &mut show(&mut p, &d));
    assert!(!painted.has("Offset U (mm)"), "placement starts collapsed");
    assert!(painted.has_prefix("Window placement · offset 3, 0 mm · 30° · border 1 mm"));
    assert_eq!(p.offset[0].resolve(&d).unwrap(), length(3.));
    // The rotation callout opens the section at its field.
    let callout = painted.find("∠ Rotation 30°").expect("rotation callout");
    harness::frame(
        &ctx,
        360.,
        harness::click(callout.center()),
        &mut show(&mut p, &d),
    );
    let painted = harness::settle(&ctx, 360., &mut show(&mut p, &d));
    assert!(painted.has("Rotation (deg)"));
    assert_eq!(p.link.active(), Some(visual::ROTATION));
    // Validation reveals a collapsed field in a fresh form.
    let ctx = egui::Context::default();
    let mut p = SurfaceStructurePanel::new().unwrap();
    p.border.text = "=".into();
    let painted = harness::settle(&ctx, 360., &mut show(&mut p, &d));
    assert_eq!(p.link.invalid, Some(visual::BORDER));
    assert!(painted.has("Face border (mm)"));
}
