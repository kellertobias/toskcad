//! CAD-392 / CAD-393: the joint form's operand state machine and its corner
//! guide. Panels are picked on their visible faces in the production
//! scene: A then B fill their own inputs, either one can be replaced
//! deliberately, duplicates and objects that are not panels are explained,
//! and a complete pair never silently turns into a new one.
use super::*;
use crate::modeling::panel_fixture::{actor, corner, flat_panel};
use multicad_core::model::{
    PrimitiveFrame, PrimitiveParameters, PrimitiveRequest, length as mm_len,
};
use multicad_model::Compiler;

/// The side's top end, then the bottom's top face, seen from above.
const SIDE: [f32; 3] = [9., 200., 300.];
const BOTTOM: [f32; 3] = [300., 200., 18.];

fn click(form: &mut JointPanel, s: &Session, point: [f32; 3]) -> bool {
    crate::modeling::panel_fixture::click(s, point, |view| {
        form.viewport(s.document(), s.revision(), view)
    })
}

#[test]
fn a_then_b_fill_their_own_inputs_and_a_complete_pair_waits_for_a_choice() {
    let (s, side, bottom) = corner(false);
    let mut form = JointPanel::new("Corner".into());
    assert_eq!(form.active, Some((0, 0)));
    assert!(click(&mut form, &s, SIDE));
    assert_eq!(form.pairs, vec![[Some(side), None]]);
    assert_eq!(form.active, Some((0, 1)), "the form now waits for B");
    assert!(click(&mut form, &s, BOTTOM));
    assert_eq!(form.pairs, vec![[Some(side), Some(bottom)]]);
    assert_eq!(form.active, None);
    // Another click neither replaces an operand nor adds a hidden pair.
    assert!(!click(&mut form, &s, SIDE));
    assert_eq!(form.pairs, vec![[Some(side), Some(bottom)]]);
    let note = form.pick_note.clone().unwrap();
    assert!(note.contains("every pair is chosen"), "{note}");
    assert!(note.contains("Pick beside A or B"), "{note}");
    // The pair previews as a valid corner.
    let request = form.request(s.document(), true).unwrap();
    assert_eq!(
        request.joint.pairs,
        vec![PanelJointPair { a: side, b: bottom }]
    );
    assert_eq!(request.prepare(&s, actor()).unwrap().outputs.len(), 2);
    assert!(form.diagnoses(s.document()).iter().all(Option::is_none));
}

#[test]
fn either_operand_is_replaced_deliberately_and_duplicates_are_refused() {
    let (s, side, bottom) = corner(false);
    let mut form = JointPanel::new("Corner".into()).with_pair(Some(side), None);
    // A selected panel seeds A; the form waits for B.
    assert_eq!(form.active, Some((0, 1)));
    // The same panel cannot also be B.
    assert!(!click(&mut form, &s, SIDE));
    let note = form.pick_note.clone().unwrap();
    assert!(note.contains("Already panel A of pair 1"), "{note}");
    assert_eq!(form.pairs, vec![[Some(side), None]]);
    assert_eq!(form.active, Some((0, 1)), "B is still waiting");
    assert!(click(&mut form, &s, BOTTOM));
    // Choose A to replace: swapping the roles is two deliberate picks.
    form.active = Some((0, 0));
    assert!(!click(&mut form, &s, BOTTOM));
    assert!(
        form.pick_note
            .as_deref()
            .unwrap()
            .contains("Already panel B")
    );
    form.active = Some((0, 1));
    assert!(!click(&mut form, &s, BOTTOM), "B is already the bottom");
    assert_eq!(form.pairs, vec![[Some(side), Some(bottom)]]);
    // Add pair waits for the new pair's A, leaving pair 1 intact.
    form.pairs.push([None, None]);
    form.active = Some((1, 0));
    assert!(click(&mut form, &s, BOTTOM));
    assert_eq!(form.active, Some((1, 1)));
    assert_eq!(form.pairs[0], [Some(side), Some(bottom)]);
    assert_eq!(
        form.request(s.document(), true).unwrap_err(),
        "Choose panel B of pair 2"
    );
    // A pair whose inputs name one panel twice (via the Inspector) is refused.
    form.pairs[1] = [Some(bottom), Some(bottom)];
    let error = form.request(s.document(), true).unwrap_err();
    assert!(error.contains("uses Bottom as both A and B"), "{error}");
}

#[test]
fn a_part_that_is_not_a_panel_is_named_and_leaves_the_operand_waiting() {
    let (mut s, side, _) = corner(false);
    let request = PrimitiveRequest::new(
        "Knob",
        PrimitiveParameters::Box {
            size: [mm_len(40.), mm_len(40.), mm_len(40.)],
            frame: PrimitiveFrame {
                origin: [mm_len(700.), mm_len(180.), mm_len(0.)],
                ..Default::default()
            },
        },
    );
    let p = model::prepare_primitive(&s, s.revision(), actor(), &request).unwrap();
    s.commit(p.prepared).unwrap();
    let mut form = JointPanel::new("Corner".into());
    assert!(click(&mut form, &s, SIDE));
    assert!(!click(&mut form, &s, [720., 200., 40.]));
    let note = form.pick_note.clone().unwrap();
    assert!(note.contains("Knob is not a cabinet panel"), "{note}");
    assert_eq!(form.pairs, vec![[Some(side), None]]);
    assert_eq!(form.active, Some((0, 1)));
}

#[test]
fn a_butted_corner_is_diagnosed_in_the_form_and_the_named_fix_previews() {
    let (mut s, side, bottom) = corner(true);
    let mut form = JointPanel::new("Corner".into());
    form.kind = Kind::Miter;
    assert!(click(&mut form, &s, SIDE));
    assert!(click(&mut form, &s, [300., 200., 18.]));
    assert_eq!(form.pairs, vec![[Some(side), Some(bottom)]]);
    let guides = form.diagnoses(s.document());
    let guide = guides[0].as_ref().expect("the butted corner is explained");
    assert_eq!(guide.extensions.len(), 1);
    assert_eq!(guide.extensions[0].panel, bottom);
    assert!(
        guide
            .message
            .contains("lengthen Bottom (B) by 18 mm at its start (length 582 → 600 mm)"),
        "{}",
        guide.message
    );
    // Preview is refused atomically with the same guidance.
    let before = s.document().clone();
    let error = form
        .request(s.document(), true)
        .unwrap()
        .prepare(&s, actor())
        .unwrap_err();
    assert_eq!(error, guide.message);
    assert_eq!(s.document(), &before);
    // Loose tenons and dowels join face to edge and get no corner guide.
    form.kind = Kind::LooseTenon;
    assert!(form.diagnoses(s.document()).iter().all(Option::is_none));
    form.kind = Kind::Miter;
    // Apply the named change with the ordinary panel edit.
    let ext = guide.extensions[0].clone();
    let feature = model::panel_feature(s.document(), bottom).unwrap();
    let mut def = model::read_panel(s.document(), feature).unwrap();
    let model::PanelSource::Rectangle { stock, plane, .. } = def.source.clone() else {
        panic!("a rectangle panel");
    };
    def.source = model::PanelSource::Rectangle {
        stock,
        plane,
        origin: ext.new_origin_mm.map(length),
        length: length(ext.required_mm()),
        width: length(400.),
    };
    let p = model::prepare_panel_edit(&s, s.revision(), actor(), feature, &def).unwrap();
    s.commit(p.prepared).unwrap();
    assert!(form.diagnoses(s.document()).iter().all(Option::is_none));
    let prepared = form
        .request(s.document(), true)
        .unwrap()
        .prepare(&s, actor())
        .unwrap();
    let model = Compiler::default()
        .compile(prepared.prepared.document(), s.revision() + 1)
        .unwrap();
    let v = model
        .solid(model.snapshot().bodies[&bottom])
        .unwrap()
        .volume_mm3()
        .unwrap();
    assert!(
        (v - (600. * 400. * 18. - 18. * 18. * 400. / 2.)).abs() < 1e-3,
        "{v}"
    );
}

#[test]
fn removing_a_pair_keeps_the_waiting_operand_on_its_pair() {
    let (_s, panel) = flat_panel();
    let mut form = JointPanel::new("Corner".into());
    form.pairs = vec![[Some(panel), None], [None, None], [None, None]];
    form.active = Some((2, 0));
    form.remove_pair(1);
    assert_eq!(form.pairs.len(), 2);
    assert_eq!(form.active, Some((1, 0)));
    // Removing the waiting pair moves on to the first empty operand.
    form.remove_pair(1);
    assert_eq!(form.active, Some((0, 1)));
    assert_eq!(form.next_open_slot(0), Some((0, 1)));
}
