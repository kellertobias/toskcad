//! CAD-426: deleting a part or feature with dependents through real pointer
//! events in the model tree, its confirmation, Cancel, Undo and save/reopen.
use super::*;
use multicad_core::{
    constraint_sketch::source::Source,
    document::Mutation,
    model::{
        self, BodyInput, EdgeTreatmentParameters, EdgeTreatmentProfile, EdgeTreatmentRequest,
        ExtrudeGroupDraft, ExtrudeRequest, GroupedExtrudeRequest,
    },
    project_file::{self, ProjectFile},
    references::{GeometryIntent, Origin, TopologyKind},
    sketch::Datum,
    sketch_regions,
};

fn actor() -> Actor {
    Actor::local(Channel::Ui)
}

fn edge(source: EntityId) -> GeometryIntent {
    GeometryIntent {
        id: EntityId::new(),
        kind: TopologyKind::Edge,
        origins: BTreeSet::from([Origin {
            source,
            role: "selected-source-edge".into(),
        }]),
        boundary: None,
        attachment: None,
    }
}

fn rename(session: &mut Session, id: EntityId, name: &str) {
    let mut entity = session.document().get(id).unwrap().clone();
    entity.name = name.into();
    session
        .apply(multicad_core::commands::Request {
            expected_revision: session.revision(),
            actor: actor(),
            label: "Name fixture".into(),
            edits: vec![Mutation::Replace(entity)],
        })
        .unwrap();
}

pub(super) struct Treated {
    pub sketch: EntityId,
    pub part: EntityId,
    pub extrude: EntityId,
    pub fillet: EntityId,
    pub chamfer: EntityId,
}

/// The reported model: an extruded part with a fillet and then a chamfer.
pub(super) fn treated_part(session: &mut Session, name: &str, origin: [f64; 2]) -> Treated {
    let (sketch, edits) = Source::create(Datum::Xy, &format!("{name} sketch")).unwrap();
    let document = session.document().edited(&edits).unwrap();
    let (mut source, _, _) = Source::read(&document, sketch).unwrap();
    source
        .add_rectangle(origin, [origin[0] + 30., origin[1] + 20.])
        .unwrap();
    let mut edits = edits;
    let replaced = source.replace(&document, sketch).unwrap();
    edits.extend(replaced);
    session
        .apply(multicad_core::commands::Request {
            expected_revision: session.revision(),
            actor: actor(),
            label: "Sketch fixture".into(),
            edits,
        })
        .unwrap();
    let regions = sketch_regions::extract(&source)
        .regions
        .iter()
        .map(|region| region.intent())
        .collect();
    let request = GroupedExtrudeRequest {
        extrusion: ExtrudeRequest::new(sketch, regions, Default::default()),
        group: ExtrudeGroupDraft::new("New body extrusion", None),
    };
    let extrusion =
        model::prepare_grouped_extrude(session, session.revision(), actor(), &request).unwrap();
    let output = extrusion.outputs[0].clone();
    session.commit(extrusion.prepared).unwrap();
    let mut treat = |producer, profile| {
        let request = EdgeTreatmentRequest::new(
            BodyInput {
                body: output.body,
                producer,
            },
            EdgeTreatmentParameters {
                profile,
                size: model::length(1.),
            },
            vec![edge(producer)],
        );
        let prepared =
            model::prepare_edge_treatment(session, session.revision(), actor(), &request).unwrap();
        session.commit(prepared.prepared).unwrap();
        request.feature
    };
    let fillet = treat(output.feature, EdgeTreatmentProfile::Fillet);
    let chamfer = treat(fillet, EdgeTreatmentProfile::Chamfer);
    rename(session, output.body, name);
    rename(session, output.feature, &format!("{name} extrusion"));
    rename(session, fillet, &format!("Fillet of {name}"));
    rename(session, chamfer, &format!("Chamfer of {name}"));
    Treated {
        sketch,
        part: output.body,
        extrude: output.feature,
        fillet,
        chamfer,
    }
}

struct Harness {
    context: egui::Context,
    tree: ModelTree,
    session: Session,
    selected: BTreeSet<EntityId>,
}
impl Harness {
    fn new(session: Session) -> Self {
        let context = egui::Context::default();
        context.enable_accesskit();
        Self {
            context,
            tree: ModelTree::default(),
            session,
            selected: BTreeSet::new(),
        }
    }
    fn frame(&mut self, events: Vec<egui::Event>) -> egui::FullOutput {
        self.context.run_ui(
            egui::RawInput {
                focused: true,
                screen_rect: Some(egui::Rect::from_min_size(
                    egui::Pos2::ZERO,
                    egui::vec2(1100., 700.),
                )),
                events,
                ..Default::default()
            },
            |ui| {
                egui::Panel::left("cascade-tree")
                    .exact_size(300.)
                    .show(ui, |ui| {
                        let event = self.tree.show(
                            ui,
                            &mut self.session,
                            &self.selected,
                            TreeOptions::default(),
                        );
                        assert!(event.error.is_none(), "{:?}", event.error);
                        if let Some(selected) = event.selection {
                            self.selected = selected;
                        }
                    });
                egui::CentralPanel::default().show(ui, |_| {});
            },
        )
    }
    fn texts(&mut self) -> Vec<(String, egui::Rect)> {
        self.frame(vec![]).drop_without_applying_deltas();
        let output = self.frame(vec![]);
        fn painted(shape: &egui::Shape, found: &mut Vec<(String, egui::Rect)>) {
            match shape {
                egui::Shape::Text(text) => found.push((
                    text.galley.text().to_owned(),
                    egui::Rect::from_min_size(text.pos, text.galley.size()),
                )),
                egui::Shape::Vec(shapes) => {
                    for shape in shapes {
                        painted(shape, found);
                    }
                }
                _ => {}
            }
        }
        let mut found = vec![];
        for shape in &output.shapes {
            painted(&shape.shape, &mut found);
        }
        output.drop_without_applying_deltas();
        found
    }
    fn shown(&mut self, label: &str) -> bool {
        self.texts().iter().any(|(text, _)| text == label)
    }
    fn bounds(&mut self, label: &str) -> egui::Rect {
        self.texts()
            .into_iter()
            .find(|(text, _)| text == label)
            .unwrap_or_else(|| panic!("Missing painted {label}"))
            .1
    }
    fn pointer(&mut self, pos: egui::Pos2, button: egui::PointerButton) {
        for pressed in [true, false] {
            self.frame(vec![
                egui::Event::PointerMoved(pos),
                egui::Event::PointerButton {
                    pos,
                    button,
                    pressed,
                    modifiers: egui::Modifiers::NONE,
                },
            ])
            .drop_without_applying_deltas();
        }
    }
    fn click(&mut self, label: &str) {
        let pos = self.bounds(label).center();
        self.pointer(pos, egui::PointerButton::Primary);
    }
    fn context_delete(&mut self, row: &str) {
        let pos = self.bounds(row).center();
        self.pointer(pos, egui::PointerButton::Secondary);
        self.click("Delete");
    }
}

fn reopened(document: &Document) -> Document {
    let mut bytes = vec![];
    project_file::write(
        &ProjectFile::new(document.clone()).unwrap(),
        &mut bytes,
        Default::default(),
    )
    .unwrap();
    project_file::read(
        &mut bytes.as_slice(),
        Default::default(),
        &project_file::supported_capabilities(),
    )
    .unwrap()
    .project
    .document
}

#[test]
fn deleting_an_extruded_part_with_fillet_and_chamfer_confirms_and_removes_them_in_one_undo_step() {
    let mut session = Session::new(Document::new("Reported model").unwrap()).unwrap();
    let kept = treated_part(&mut session, "Spacer", [60., 0.]);
    let target = treated_part(&mut session, "Bracket", [0., 0.]);
    let before = session.document().clone();
    let revision = session.revision();
    let mut h = Harness::new(session);

    let prompt = "Delete Bracket and 2 dependent features (Fillet of Bracket, Chamfer of Bracket)?";
    h.context_delete("Bracket");
    assert!(
        h.shown(prompt),
        "Delete must ask before removing dependents; painted: {:?}",
        h.texts().into_iter().map(|(t, _)| t).collect::<Vec<_>>()
    );
    assert_eq!(h.session.document(), &before, "asking changes nothing");

    h.click("Cancel");
    assert!(!h.shown(prompt));
    assert_eq!(h.session.document(), &before, "Cancel changes nothing");
    assert_eq!(h.session.revision(), revision);

    h.context_delete("Bracket");
    h.click("Delete all");
    assert!(!h.shown(prompt));
    assert_eq!(h.session.revision(), revision + 1, "one transaction");
    for id in [target.part, target.extrude, target.fillet, target.chamfer] {
        assert!(h.session.document().get(id).is_none(), "{id} remains");
    }
    assert_eq!(
        h.session.document().get(target.sketch),
        before.get(target.sketch),
        "the consumed sketch is an input, not a dependent"
    );
    for id in [
        kept.sketch,
        kept.part,
        kept.extrude,
        kept.fillet,
        kept.chamfer,
    ] {
        assert_eq!(h.session.document().get(id), before.get(id), "unrelated");
    }
    assert!(h.session.document().reference_diagnostics().is_empty());
    assert!(h.selected.is_empty());
    let after = h.session.document().clone();
    assert_eq!(&reopened(&after), &after, "save/reopen after the delete");

    h.session.undo(h.session.revision(), actor()).unwrap();
    assert_eq!(h.session.document(), &before, "one Undo step restores all");
    h.session.redo(h.session.revision(), actor()).unwrap();
    assert_eq!(h.session.document(), &after);
}

#[test]
fn deleting_a_feature_from_part_history_removes_its_nested_dependents_and_keeps_the_part() {
    let mut session = Session::new(Document::new("Nested history").unwrap()).unwrap();
    let target = treated_part(&mut session, "Bracket", [0., 0.]);
    let before = session.document().clone();
    let mut h = Harness::new(session);
    // Feature history rows are listed under their part.
    h.tree.set_expanded(target.part, true);
    let prompt = "Delete Fillet of Bracket and 1 dependent feature (Chamfer of Bracket)?";
    h.context_delete("1. Fillet of Bracket");
    assert!(h.shown(prompt));
    h.click("Delete all");
    let document = h.session.document();
    assert!(document.get(target.fillet).is_none());
    assert!(document.get(target.chamfer).is_none());
    let Some(EntityData::Part { features, .. }) = document.get(target.part).map(|e| &e.data) else {
        panic!("the part stays");
    };
    assert_eq!(features, &vec![target.extrude]);
    assert_eq!(document.get(target.extrude), before.get(target.extrude));
    h.session.undo(h.session.revision(), actor()).unwrap();
    assert_eq!(h.session.document(), &before);
}

#[test]
fn deleting_the_last_feature_needs_no_confirmation() {
    let mut session = Session::new(Document::new("Leaf feature").unwrap()).unwrap();
    let target = treated_part(&mut session, "Bracket", [0., 0.]);
    let revision = session.revision();
    let mut h = Harness::new(session);
    h.tree.set_expanded(target.part, true);
    h.context_delete("2. Chamfer of Bracket");
    assert!(!h.shown("Delete all"), "no dependents, no extra prompt");
    assert_eq!(h.session.revision(), revision + 1);
    assert!(h.session.document().get(target.chamfer).is_none());
    assert!(h.session.document().get(target.fillet).is_some());
}
