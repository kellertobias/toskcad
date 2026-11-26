//! CAD-408: Boolean in plain words. Join and Intersect take two or more
//! bodies with no roles; Cut names "Body to cut" and "Bodies to cut away";
//! each message appears once and says what to do next. The model is the
//! reported one: a box and the corner-circle segment extruded from its top
//! face as a second body.
use super::ui_harness::Harness;
use super::*;
use crate::modeling::native_fixture::corner_segment_session;
use multicad_core::model::{self, Operation};
use multicad_model::Compiler;
use std::collections::BTreeSet;

fn open(h: &mut Harness, operation: Operation, selection: BTreeSet<EntityId>) {
    h.app
        .current_mut()
        .sketch
        .set_tree_selection_for_tests(selection);
    h.app
        .invoke(Command::BooleanOperation(operation), &[], Channel::Ui)
        .unwrap();
    h.frame(vec![]);
}
fn select(h: &mut Harness, selection: BTreeSet<EntityId>) {
    h.app
        .current_mut()
        .sketch
        .set_tree_selection_for_tests(selection);
    h.frame(vec![]);
    h.frame(vec![]);
}
fn can_apply(app: &Application) -> bool {
    let p = app.current();
    p.sketch.model_can_apply(&p.session, &p.resources)
}
fn count(texts: &[String], text: &str) -> usize {
    texts.iter().filter(|t| t.as_str() == text).count()
}
/// No visible text speaks of tools or targets.
fn assert_plain(texts: &[String]) {
    for text in texts {
        let lower = text.to_lowercase();
        assert!(
            !lower.contains("tool bod")
                && !lower.starts_with("tools ·")
                && !lower.contains("target bod")
                && !lower.starts_with("targets ·"),
            "role jargon shown: {text:?}"
        );
        assert!(
            !text.contains("1 and 128"),
            "internal limit shown: {text:?}"
        );
    }
}
fn volume(document: &Document, revision: u64, body: EntityId) -> Option<f64> {
    let native = Compiler::default().compile(document, revision).ok()?;
    native.solid(body)?.volume_mm3().ok()
}

#[test]
fn preselected_join_of_box_and_corner_segment_previews_and_enter_merges_them() {
    let (session, stock, segment) = corner_segment_session(5.);
    let document = session.document().clone();
    let box_volume = volume(&document, session.revision(), stock.body).unwrap();
    let segment_volume = volume(&document, session.revision(), segment.body).unwrap();
    let mut h = Harness::new(document);
    open(
        &mut h,
        Operation::Join,
        BTreeSet::from([stock.body, segment.body]),
    );
    // Both pre-selected bodies fill the one list; no roles to assign.
    h.until("the Join preview", can_apply);
    let texts = h.texts();
    assert_plain(&texts);
    assert_eq!(count(&texts, "Bodies to join · 2"), 1, "{texts:?}");
    assert!(!texts.iter().any(|t| t.starts_with("Select")), "{texts:?}");
    // Enter applies exactly the previewed union (CAD-406 exact curved Join).
    h.key(egui::Key::Enter);
    h.until("Enter to apply the Join", |app| {
        app.current().session.revision() == 1
    });
    assert!(!h.modeling());
    let joined = h.app.current().session.document().clone();
    let (target, tool) = {
        let EntityData::Part { features, .. } = &joined
            .entities()
            .find(|e| {
                matches!(&e.data, EntityData::Part { features, .. }
                    if features.iter().any(|f| model::read_boolean(&joined, *f).is_ok()))
            })
            .expect("a Part with the Join")
            .data
        else {
            unreachable!()
        };
        let saved = model::read_boolean(&joined, *features.last().unwrap()).unwrap();
        assert_eq!(saved.operation, Operation::Join);
        assert!(!saved.keep_tools);
        assert_eq!(saved.tools.len(), 1);
        (saved.target.body, saved.tools[0].body)
    };
    assert_eq!(
        BTreeSet::from([target, tool]),
        BTreeSet::from([stock.body, segment.body])
    );
    // One body remains: the joined one holds both volumes and the other is
    // consumed (kept only as source history).
    let native = Compiler::default().compile(&joined, 1).unwrap();
    assert!(native.snapshot().consumed_bodies.contains(&tool));
    assert!(!native.snapshot().consumed_bodies.contains(&target));
    let joined_volume = native.solid(target).unwrap().volume_mm3().unwrap();
    let expected = box_volume + segment_volume;
    assert!(
        (joined_volume - expected).abs() < 1e-6 * expected,
        "{joined_volume} != {expected}"
    );
}

#[test]
fn join_selected_after_opening_follows_the_selection_with_one_hint() {
    let (session, stock, segment) = corner_segment_session(5.);
    let mut h = Harness::new(session.document().clone());
    open(&mut h, Operation::Join, BTreeSet::new());
    assert!(h.app.current().sketch.awaiting_model_targets());
    let texts = h.texts();
    assert_plain(&texts);
    assert_eq!(count(&texts, "Join · select bodies"), 1, "{texts:?}");
    // One body is not enough to open the Join.
    select(&mut h, BTreeSet::from([stock.body]));
    assert!(h.app.current().sketch.awaiting_model_targets());
    // The second body opens it with both filled in; no Continue click.
    select(&mut h, BTreeSet::from([stock.body, segment.body]));
    h.until("the Join preview", can_apply);
    assert!(!h.app.current().sketch.awaiting_model_targets());
    // Deselecting a body in the viewer or browser updates the list: one
    // plain hint, shown once, and Apply is disabled.
    select(&mut h, BTreeSet::from([stock.body]));
    let texts = h.texts();
    assert_plain(&texts);
    assert_eq!(
        count(&texts, "Select at least one more body to join"),
        1,
        "{texts:?}"
    );
    assert_eq!(count(&texts, "Bodies to join · 1"), 1, "{texts:?}");
    assert!(!can_apply(&h.app));
    // Selecting it again previews the union; Apply commits.
    select(&mut h, BTreeSet::from([stock.body, segment.body]));
    h.until("the Join preview", can_apply);
    h.click_label("Apply");
    h.until("Apply to commit the Join", |app| {
        app.current().session.revision() == 1
    });
}

#[test]
fn cut_names_its_roles_and_a_body_moves_to_bodies_to_cut_away() {
    let (session, stock, segment) = corner_segment_session(5.);
    let segment_name = session.document().get(segment.body).unwrap().name.clone();
    let mut h = Harness::new(session.document().clone());
    open(
        &mut h,
        Operation::Cut,
        BTreeSet::from([stock.body, segment.body]),
    );
    assert!(h.modeling());
    // Pre-selected bodies fill the first role; the hint names the next step.
    let texts = h.texts();
    assert_plain(&texts);
    assert_eq!(count(&texts, "Body to cut · 2"), 1, "{texts:?}");
    assert_eq!(count(&texts, "Bodies to cut away · 0"), 1, "{texts:?}");
    assert_eq!(
        count(&texts, "Select one or more bodies to cut away"),
        1,
        "{texts:?}"
    );
    assert!(!can_apply(&h.app));
    // Tick the segment in "Bodies to cut away": it leaves "Body to cut".
    h.click_nth_checkbox(&format!("{segment_name} ·"), 1);
    let texts = h.texts();
    assert_plain(&texts);
    assert_eq!(count(&texts, "Body to cut · 1"), 1, "{texts:?}");
    assert_eq!(count(&texts, "Bodies to cut away · 1"), 1, "{texts:?}");
    assert!(!texts.iter().any(|t| t.starts_with("Select")), "{texts:?}");
    assert_eq!(count(&texts, "Keep the bodies cut away"), 1, "{texts:?}");
    // The viewer selection still shows both bodies.
    h.frame(vec![]);
    assert_eq!(
        h.app.current().sketch.tree_selection_for_tests(),
        &BTreeSet::from([stock.body, segment.body])
    );
    h.key(egui::Key::Escape);
    assert!(!h.modeling());
    assert_eq!(h.revision(), 0);
}

/// The body a committed Join kept (its target).
fn join_target(document: &Document) -> EntityId {
    document
        .entities()
        .filter_map(|e| match &e.data {
            EntityData::Part { features, .. } => features
                .iter()
                .find_map(|f| model::read_boolean(document, *f).ok()),
            _ => None,
        })
        .next()
        .expect("a committed Join")
        .target
        .body
}

#[test]
fn join_keeps_the_first_selected_body_and_the_order_can_be_changed() {
    let (session, stock, segment) = corner_segment_session(5.);
    let document = session.document().clone();
    // Picked one at a time: the first pick keeps its name and identity,
    // although it is listed second alphabetically.
    let mut h = Harness::new(document.clone());
    open(&mut h, Operation::Join, BTreeSet::new());
    select(&mut h, BTreeSet::from([segment.body]));
    select(&mut h, BTreeSet::from([stock.body, segment.body]));
    h.until("the Join preview", can_apply);
    let texts = h.texts();
    let first = texts
        .iter()
        .position(|t| t.starts_with("1. "))
        .expect("an order list");
    let segment_name = document.get(segment.body).unwrap().name.clone();
    assert!(
        texts[first].starts_with(&format!("1. {segment_name}")),
        "{texts:?}"
    );
    assert_eq!(count(&texts, "keeps its name"), 1, "{texts:?}");
    h.key(egui::Key::Enter);
    h.until("the Join", |app| app.current().session.revision() == 1);
    assert_eq!(
        join_target(h.app.current().session.document()),
        segment.body
    );

    // Pre-selected together: list order, then "Move up" reorders it.
    let mut h = Harness::new(document);
    open(
        &mut h,
        Operation::Join,
        BTreeSet::from([stock.body, segment.body]),
    );
    h.until("the Join preview", can_apply);
    let texts = h.texts();
    assert!(texts.iter().any(|t| t.starts_with("1. Box")), "{texts:?}");
    h.click_label("Move up");
    let texts = h.texts();
    assert!(
        texts
            .iter()
            .any(|t| t.starts_with(&format!("1. {segment_name}"))),
        "{texts:?}"
    );
    h.until("the reordered Join preview", can_apply);
    h.key(egui::Key::Enter);
    h.until("the reordered Join", |app| {
        app.current().session.revision() == 1
    });
    assert_eq!(
        join_target(h.app.current().session.document()),
        segment.body
    );
}
