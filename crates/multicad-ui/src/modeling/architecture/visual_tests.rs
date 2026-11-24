//! CAD-395 architecture schematics: stair, floor slab and pedestal callouts
//! edit the panel's own staged fields; secondary inputs are named sections.
use super::*;
use crate::modeling::diagram::{
    Mark,
    harness::{self, WIDTHS},
};
use multicad_core::{commands::Channel, model::Turn};

/// A document with one level (2750 mm floor-to-floor, 200 mm slab).
fn leveled() -> Session {
    let mut session = Session::new(Document::new("Architecture diagrams").unwrap()).unwrap();
    let panel = ArchitecturePanel::new(session.document(), "Ground".into()).unwrap();
    let prepared = panel
        .request(session.document(), true)
        .unwrap()
        .prepare(&session, Actor::local(Channel::Ui))
        .unwrap();
    session.commit(prepared.prepared).unwrap();
    session
}

fn panel(session: &Session, kind: ArchitectureKind) -> ArchitecturePanel {
    let mut panel = ArchitecturePanel::new(session.document(), kind.label().into()).unwrap();
    panel.kind = kind;
    panel
}

fn show<'a>(
    panel: &'a mut ArchitecturePanel,
    document: &'a Document,
) -> impl FnMut(&mut egui::Ui) + 'a {
    move |ui| {
        let _ = panel.show(ui, document, true, PreviewState::Idle);
    }
}

fn mm(field: &DimensionField, document: &Document) -> Dimension {
    field.resolve(document).unwrap()
}

#[test]
fn every_stair_form_draws_its_own_plan_with_callouts_on_linked_inputs() {
    let session = leveled();
    let document = session.document();
    let mut p = panel(&session, ArchitectureKind::Stair);
    let linked = [
        visual::STAIR_WIDTH,
        visual::STAIR_GOING,
        visual::STAIR_RISERS,
        visual::STAIR_DIRECTION,
        visual::STAIR_START,
    ];
    let mut seen = vec![];
    for form in visual::FORMS {
        p.form = form;
        let values = p.stair_values(document);
        assert!(
            (values.riser.unwrap() - 2750. / 16.).abs() < 1e-9,
            "the elevation uses the level's actual riser"
        );
        let schematic = visual::stair(form, &values);
        let mut keys = schematic.keys();
        keys.sort_unstable();
        let mut expected = linked.to_vec();
        expected.sort_unstable();
        assert_eq!(keys, expected, "{form:?}");
        assert!(schematic.marks.iter().any(
            |m| matches!(m, Mark::Note { text, .. } if text.starts_with(visual::form_label(form)))
        ));
        assert!(!seen.contains(&schematic), "{form:?} reuses another plan");
        seen.push(schematic);
    }
    // The tiles illustrate distinct flights and landings.
    let left = visual::form_art(StairForm::QuarterTurn { turn: Turn::Left });
    let right = visual::form_art(StairForm::QuarterTurn { turn: Turn::Right });
    assert_ne!(left, right);
    // Values only proportion the drawing; the request is unchanged by it.
    let before = p.request(document, true).unwrap();
    let _ = visual::stair(p.form, &p.stair_values(document));
    assert_eq!(
        format!("{:?}", p.request(document, true).unwrap()),
        format!("{before:?}")
    );
}

#[test]
fn clicking_the_going_callout_edits_the_stair_going_at_compact_and_ordinary_widths() {
    let session = leveled();
    let document = session.document();
    for width in WIDTHS {
        let ctx = egui::Context::default();
        let mut p = panel(&session, ArchitectureKind::Stair);
        let painted = harness::edit_through_callout(
            &ctx,
            width,
            "Going 270 mm",
            "300",
            &mut show(&mut p, document),
        );
        assert_eq!(mm(&p.going, document), length(300.));
        assert_eq!(mm(&p.width, document), length(1000.), "only going changed");
        assert!(painted.has("Going 300 mm"), "the list updates the drawing");
        let ArchitectureCommand::Stair { stair, .. } = p.request(document, true).unwrap() else {
            panic!("stair")
        };
        assert_eq!(stair.going, length(300.));
    }
}

#[test]
fn stair_position_is_collapsed_retained_and_revealed_by_its_callout() {
    let session = leveled();
    let document = session.document();
    let ctx = egui::Context::default();
    let mut p = panel(&session, ArchitectureKind::Stair);
    p.start = [1200., 300.];
    p.direction_deg = 90.;
    let painted = harness::settle(&ctx, 360., &mut show(&mut p, document));
    assert!(!painted.has("Direction °"), "position starts collapsed");
    assert!(painted.has_prefix("Position · 1200, 300 mm · 90°"));
    let ArchitectureCommand::Stair { stair, .. } = p.request(document, true).unwrap() else {
        panic!("stair")
    };
    assert_eq!(stair.start, [length(1200.), length(300.)]);
    assert_eq!(stair.direction_deg, 90.);
    let callout = painted.find("∠ Direction 90°").expect("direction callout");
    harness::frame(
        &ctx,
        360.,
        harness::click(callout.center()),
        &mut show(&mut p, document),
    );
    let painted = harness::settle(&ctx, 360., &mut show(&mut p, document));
    assert!(painted.has("Direction °"), "the callout opens its section");
    assert_eq!(p.link.active(), Some(visual::STAIR_DIRECTION));
    assert_eq!(p.direction_deg, 90.);
}

#[test]
fn clicking_the_pedestal_height_callout_edits_its_height_at_both_widths() {
    let session = leveled();
    let document = session.document();
    for width in WIDTHS {
        let ctx = egui::Context::default();
        let mut p = panel(&session, ArchitectureKind::Pedestal);
        let painted = harness::edit_through_callout(
            &ctx,
            width,
            "Height 200 mm",
            "450",
            &mut show(&mut p, document),
        );
        assert!(painted.has("Height 450 mm"));
        let ArchitectureCommand::Slab { slab, .. } = p.request(document, true).unwrap() else {
            panic!("slab")
        };
        assert_eq!(slab.kind, SlabKind::Pedestal);
        assert_eq!(slab.thickness, Some(length(450.)));
    }
}

#[test]
fn slab_thickness_override_is_collapsed_retained_and_revealed_by_validation() {
    let session = leveled();
    let document = session.document();
    let ctx = egui::Context::default();
    let mut p = panel(&session, ArchitectureKind::Slab);
    // Without an override the callout reads the level and opens the override.
    let painted = harness::settle(&ctx, 360., &mut show(&mut p, document));
    assert!(!painted.has("Override the level's slab thickness"));
    let callout = painted
        .find("Level thickness 200 mm")
        .expect("level callout");
    harness::frame(
        &ctx,
        360.,
        harness::click(callout.center()),
        &mut show(&mut p, document),
    );
    let painted = harness::settle(&ctx, 360., &mut show(&mut p, document));
    assert!(painted.has("Override the level's slab thickness"));
    assert_eq!(p.link.active(), Some(visual::SLAB_OVERRIDE));

    // A fresh form keeps a collapsed override in its request.
    let ctx = egui::Context::default();
    let mut p = panel(&session, ArchitectureKind::Slab);
    p.override_thickness = true;
    p.thickness.text = "250".into();
    let painted = harness::settle(&ctx, 360., &mut show(&mut p, document));
    assert!(
        !painted.has("Thickness (mm)"),
        "the section starts collapsed"
    );
    assert!(painted.has_prefix("Thickness · 250 mm"));
    assert!(
        painted.has("Thickness 250 mm"),
        "the callout shows the override"
    );
    let ArchitectureCommand::Slab { slab, .. } = p.request(document, true).unwrap() else {
        panic!("slab")
    };
    assert_eq!(slab.thickness, Some(length(250.)));
    p.thickness.text = "=".into();
    let painted = harness::settle(&ctx, 360., &mut show(&mut p, document));
    assert_eq!(p.link.invalid, Some(visual::SLAB_THICKNESS));
    assert!(
        painted.has("Thickness (mm)"),
        "validation reveals the field"
    );
}
