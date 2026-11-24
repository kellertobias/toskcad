//! CAD-395 FDM aid schematics: each kind draws its own elevation and plan,
//! callouts edit the panel's staged fields and position sections open on
//! demand or on validation.
use super::*;
use crate::modeling::diagram::{
    Key,
    harness::{self, WIDTHS},
};
use multicad_core::{
    commands::Channel,
    model::{PrimitiveFrame, PrimitiveParameters, PrimitiveRequest},
};

fn session() -> Session {
    let mut s = Session::new(Document::new("FDM diagrams").unwrap()).unwrap();
    let block = PrimitiveRequest::new(
        "Block",
        PrimitiveParameters::Box {
            size: [length(40.), length(20.), length(10.)],
            frame: PrimitiveFrame {
                origin: [length(0.), length(0.), length(15.)],
                ..Default::default()
            },
        },
    );
    let p = model::prepare_primitive(&s, 0, Actor::local(Channel::Ui), &block).unwrap();
    s.commit(p.prepared).unwrap();
    s
}

fn panel(s: &Session, kind: AidKind) -> FdmAidPanel {
    let mut p = FdmAidPanel::new(s.document(), visual::label(kind).into());
    p.choose_kind(kind).unwrap();
    // Coordinates are typed here, not picked in the viewer.
    p.use_coordinates();
    p
}

fn show<'a>(p: &'a mut FdmAidPanel, d: &'a Document) -> impl FnMut(&mut egui::Ui) + 'a {
    move |ui| {
        let _ = p.show(ui, d, true, PreviewState::Idle);
    }
}

fn shape(p: &FdmAidPanel, d: &Document) -> AidShape {
    p.draft_request(d, true).unwrap().request.aid.shape
}

#[test]
fn every_aid_kind_draws_its_own_schematic_with_only_linked_callouts() {
    let s = session();
    let d = s.document();
    let mut seen = vec![];
    for (kind, essentials) in [
        (
            AidKind::Tab,
            vec![
                visual::LENGTH,
                visual::OVERLAP,
                visual::WIDTH,
                visual::THICKNESS,
            ],
        ),
        (AidKind::Bridge, vec![visual::WIDTH, visual::THICKNESS]),
        (
            AidKind::Support,
            vec![visual::LENGTH, visual::WIDTH, visual::BASE, visual::GAP],
        ),
    ] {
        let p = panel(&s, kind);
        let schematic = visual::diagram(kind, &p.diagram_values(d));
        let fields: Vec<Key> = p.linked_fields().into_iter().map(|(k, _)| k).collect();
        for key in schematic.keys() {
            assert!(fields.contains(&key), "{kind:?} callout {key} has no field");
        }
        for key in essentials {
            assert!(schematic.keys().contains(&key), "{kind:?} lacks {key}");
        }
        assert!(!seen.contains(&schematic), "{kind:?} reuses a drawing");
        seen.push(schematic);
    }
    // The bridge span comes from its staged anchors.
    let p = panel(&s, AidKind::Bridge);
    assert_eq!(p.diagram_values(d).span, Some(20.));
}

#[test]
fn clicking_each_kinds_callout_edits_its_staged_field_at_both_widths() {
    let s = session();
    let d = s.document();
    for width in WIDTHS {
        for (kind, label, text, typed) in [
            (AidKind::Tab, "Reach 8 mm", "12", "Reach 12 mm"),
            (
                AidKind::Bridge,
                "Thickness 0.4 mm",
                "0.6",
                "Thickness 0.6 mm",
            ),
            (AidKind::Support, "Top gap 0.2 mm", "0.3", "Top gap 0.3 mm"),
        ] {
            let ctx = egui::Context::default();
            let mut p = panel(&s, kind);
            let painted =
                harness::edit_through_callout(&ctx, width, label, text, &mut show(&mut p, d));
            assert!(
                painted.has(typed),
                "{kind:?} at {width}: the drawing follows"
            );
            match shape(&p, d) {
                AidShape::Tab {
                    length: l,
                    width: w,
                    ..
                } => {
                    assert_eq!(l, length(12.));
                    assert_eq!(w, length(6.), "only the reach changed");
                }
                AidShape::Bridge {
                    thickness,
                    width: w,
                    ..
                } => {
                    assert_eq!(thickness, length(0.6));
                    assert_eq!(w, length(6.));
                }
                AidShape::Support { gap, base, .. } => {
                    assert_eq!(gap, length(0.3));
                    assert_eq!(base, length(0.));
                }
            }
        }
    }
}

#[test]
fn position_sections_are_collapsed_retained_and_revealed_by_callouts_and_validation() {
    let s = session();
    let d = s.document();
    let ctx = egui::Context::default();
    let mut p = panel(&s, AidKind::Tab);
    p.origin[0].text = "5".into();
    p.angle.text = "90".into();
    let painted = harness::settle(&ctx, 360., &mut show(&mut p, d));
    assert!(!painted.has("Contact X (mm)"), "contact starts collapsed");
    assert!(painted.has_prefix("Contact and direction · 5, 0, 0 mm · 90°"));
    let AidShape::Tab {
        origin, angle: a, ..
    } = shape(&p, d)
    else {
        panic!("tab")
    };
    assert_eq!(origin[0], length(5.));
    assert_eq!(a, angle(90.));
    p.origin[2].text = "=".into();
    let painted = harness::settle(&ctx, 360., &mut show(&mut p, d));
    assert_eq!(p.link.invalid, Some(visual::ORIGIN[2]));
    assert!(
        painted.has("Contact Z (mm)"),
        "validation reveals the contact"
    );

    // A bridge's End callout opens its anchors and focuses End X.
    let ctx = egui::Context::default();
    let mut p = panel(&s, AidKind::Bridge);
    let painted = harness::settle(&ctx, 360., &mut show(&mut p, d));
    assert!(!painted.has("End X (mm)"));
    let end = painted.find("End 20, 0, 0 mm").expect("end callout");
    harness::frame(
        &ctx,
        360.,
        harness::click(end.center()),
        &mut show(&mut p, d),
    );
    let painted = harness::settle(&ctx, 360., &mut show(&mut p, d));
    assert!(painted.has("End X (mm)"));
    assert_eq!(p.link.active(), Some(visual::END[0]));
}

#[test]
fn kind_tiles_switch_new_aids_and_stay_fixed_when_editing() {
    let s = session();
    let d = s.document();
    let ctx = egui::Context::default();
    let mut p = panel(&s, AidKind::Tab);
    let painted = harness::settle(&ctx, 360., &mut show(&mut p, d));
    let tile = painted
        .texts
        .iter()
        .find(|(t, _)| t == "Support")
        .map(|(_, r)| *r)
        .expect("support tile");
    harness::frame(
        &ctx,
        360.,
        harness::click(tile.center()),
        &mut show(&mut p, d),
    );
    assert_eq!(p.kind(), AidKind::Support);
    assert!(p.choose_kind(AidKind::Tab).unwrap());
    p.existing = Some(EntityId::new());
    harness::frame(
        &ctx,
        360.,
        harness::click(tile.center()),
        &mut show(&mut p, d),
    );
    assert_eq!(p.kind(), AidKind::Tab, "an existing aid keeps its type");
}
