//! CAD-430 panel inspector placement fields: they show where the Part is
//! now, follow the draft until the user types or picks here, and only then
//! turn a panel edit into a placement edit.
use super::tests::session;
use super::*;
use multicad_core::{
    commands::{Actor, Channel},
    model::{BodyInput, NumericTransform, PanelFacing, TransformRequest},
};

fn ui() -> Actor {
    Actor::local(Channel::Ui)
}
fn panel(s: &mut Session) -> (EntityId, EntityId) {
    let mut form = CabinetPanel::panel(s.document(), "Board".into());
    form.origin = [field(10., "mm"), field(20., "mm"), field(30., "mm")];
    let prepared = form
        .request(s.document(), true)
        .unwrap()
        .prepare(s, ui())
        .unwrap();
    let out = (prepared.outputs[0].feature, prepared.outputs[0].body);
    s.commit(prepared.prepared).unwrap();
    out
}
fn text(f: &DimensionField) -> &str {
    f.text.as_str()
}

#[test]
fn an_existing_panel_shows_its_placement_and_an_untouched_edit_stays_a_panel_edit() {
    let mut s = session();
    let (feature, _) = panel(&mut s);
    let form = CabinetPanel::existing(s.document(), feature).unwrap();
    let fields = form.placement.as_ref().expect("placement fields");
    assert_eq!(fields.position.each_ref().map(text), ["10", "20", "30"]);
    assert_eq!(fields.facing, Some(PanelFacing::PosZ));
    assert_eq!(text(&fields.rotation), "0");
    let Ok(CabinetCommand::Panel { placement, .. }) = form.request(s.document(), true) else {
        panic!("panel command")
    };
    assert!(placement.is_none());
    // A new panel has no placement fields: its plane and origin are drafted.
    assert!(
        CabinetPanel::panel(s.document(), "New".into())
            .placement
            .is_none()
    );
}

#[test]
fn the_fields_follow_a_move_made_elsewhere_and_the_draft_origin() {
    let mut s = session();
    let (feature, body) = panel(&mut s);
    let request = TransformRequest::new(
        BodyInput {
            body,
            producer: feature,
        },
        NumericTransform {
            translation: [length(100.), length(0.), length(0.)],
            ..Default::default()
        },
    );
    let p = model::prepare_transform(&s, s.revision(), ui(), &request).unwrap();
    s.commit(p.prepared).unwrap();
    let mut form = CabinetPanel::existing(s.document(), feature).unwrap();
    let fields = form.placement.as_mut().unwrap();
    assert_eq!(fields.position.each_ref().map(text), ["110", "20", "30"]);
    // The viewport draft moves the panel feature's origin: the fields show
    // where that puts the Part.
    fields.follow(Some([10., 20., 30.]), PanelPlane::Xy);
    fields.follow(Some([15., 20., 30.]), PanelPlane::Xy);
    assert_eq!(fields.position.each_ref().map(text), ["115", "20", "30"]);
    // Once the user types here, the typed value is kept.
    fields.edited = true;
    fields.position[0] = field(500., "mm");
    fields.follow(Some([16., 20., 30.]), PanelPlane::Xy);
    assert_eq!(text(&fields.position[0]), "500");
    let Ok(CabinetCommand::Panel {
        placement: Some(placement),
        ..
    }) = form.request(s.document(), true)
    else {
        panic!("a placement edit")
    };
    assert_eq!(placement.position[0], length(500.));
    assert_eq!(placement.facing, PanelFacing::PosZ);
}
