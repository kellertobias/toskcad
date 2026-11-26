//! CAD-404: the automatic default operation rule and its canonical record.
use multicad_core::{
    commands::{Actor, Channel, Session},
    constraint_sketch::source::Source,
    document::{Document, EntityData, EntityId},
    model::*,
    project_file::{self, ProjectFile},
    sketch::Datum,
    sketch_regions,
};

fn actor() -> Actor {
    Actor::local(Channel::Ui)
}
fn input() -> BodyInput {
    BodyInput {
        body: EntityId::new(),
        producer: EntityId::new(),
    }
}
fn selection(targets: &[&BodyInput]) -> CutTargetSelection {
    CutTargetSelection {
        defaults: targets.iter().map(|t| (*t).clone()).collect(),
        targets: targets.iter().map(|t| (*t).clone()).collect(),
        not_intersecting: vec![],
    }
}

#[test]
fn decision_joins_one_touching_body_cuts_one_entered_body_and_never_guesses() {
    let (a, b) = (input(), input());
    let none = selection(&[]);
    assert_eq!(
        AutomaticOperation::decide(&none, &none),
        AutomaticOperation::NewBody
    );
    // Touching only: Join.
    let join = AutomaticOperation::decide(&selection(&[&a]), &none);
    assert_eq!(join, AutomaticOperation::Join(a.clone()));
    assert_eq!(join.operation(), Operation::Join);
    assert_eq!(join.target(), Some(&a));
    // The contact probe also reports the entered body: Cut.
    let cut = AutomaticOperation::decide(&selection(&[&a]), &selection(&[&a]));
    assert_eq!(cut, AutomaticOperation::Cut(a.clone()));
    assert_eq!(cut.operation(), Operation::Cut);
    // An entered body missing from the contact probe still counts once.
    assert_eq!(
        AutomaticOperation::decide(&none, &selection(&[&a])),
        AutomaticOperation::Cut(a.clone())
    );
    // Several bodies, touched or entered in any mix: listed, New body.
    for (contact, volume) in [
        (selection(&[&a, &b]), none.clone()),
        (selection(&[&a]), selection(&[&b])),
        (selection(&[&a, &b]), selection(&[&a, &b])),
    ] {
        let several = AutomaticOperation::decide(&contact, &volume);
        assert_eq!(
            several,
            AutomaticOperation::Several(vec![a.clone(), b.clone()])
        );
        assert_eq!(several.operation(), Operation::NewBody);
        assert_eq!(several.target(), None);
    }
}

struct Fixture {
    session: Session,
    sketch: EntityId,
    regions: Vec<sketch_regions::RegionIntent>,
    stock: BodyInput,
}
fn fixture() -> Fixture {
    let (sketch, edits) = Source::create(Datum::Xy, "Tools").unwrap();
    let mut document = Document::new("Automatic operation")
        .unwrap()
        .edited(&edits)
        .unwrap();
    let mut source = Source::read(&document, sketch).unwrap().0;
    source.add_rectangle([0., 0.], [2., 2.]).unwrap();
    let regions = sketch_regions::extract(&source)
        .regions
        .iter()
        .map(|r| r.intent())
        .collect();
    document = document
        .edited(&source.replace(&document, sketch).unwrap())
        .unwrap();
    let mut session = Session::new(document).unwrap();
    let request = PrimitiveRequest::new(
        "Box",
        PrimitiveParameters::Box {
            size: [length(10.), length(8.), length(4.)],
            frame: Default::default(),
        },
    );
    let prepared = prepare_primitive(&session, session.revision(), actor(), &request).unwrap();
    session.commit(prepared.prepared).unwrap();
    Fixture {
        session,
        sketch,
        regions,
        stock: BodyInput {
            body: request.body,
            producer: request.feature,
        },
    }
}
fn reopen(document: &Document) -> Document {
    let mut bytes = vec![];
    project_file::write(
        &ProjectFile::new(document.clone()).unwrap(),
        &mut bytes,
        Default::default(),
    )
    .unwrap();
    project_file::read(
        &mut &bytes[..],
        Default::default(),
        &project_file::supported_capabilities(),
    )
    .unwrap()
    .project
    .document
}
fn raw(document: &Document, group: EntityId) -> (u32, serde_json::Value) {
    let EntityData::Extension { definition, .. } = &document.get(group).unwrap().data else {
        panic!()
    };
    (definition.schema_version, definition.data.clone())
}

#[test]
fn automatic_boolean_is_stored_in_schema_two_and_survives_save_and_reopen() {
    let mut f = fixture();
    for automatic in [false, true] {
        let mut request = ExtrudeBooleanRequest::new(
            ExtrudeRequest::new(f.sketch, f.regions.clone(), Default::default()),
            Operation::Join,
            vec![f.stock.clone()],
            false,
        );
        let mut group = ExtrudeGroupDraft::new("Extrusion operation", None);
        group.automatic = automatic;
        let id = group.id;
        request.group = Some(group);
        let prepared =
            prepare_extrude_boolean(&f.session, f.session.revision(), actor(), &request).unwrap();
        let document = prepared.prepared.document().clone();
        let (schema, data) = raw(&document, id);
        if automatic {
            assert_eq!(schema, EXTRUDE_GROUP_V2_SCHEMA);
            assert_eq!(data["automatic"], serde_json::Value::Bool(true));
        } else {
            // Explicit choices keep the existing record unchanged.
            assert_eq!(schema, EXTRUDE_GROUP_SCHEMA);
            assert!(data.get("automatic").is_none());
        }
        let reopened = reopen(&document);
        let state = inspect_extrude_group(&reopened, id).unwrap();
        assert_eq!(state.definition.automatic, automatic);
        assert_eq!(state.operation().unwrap(), Operation::Join);
        assert_eq!(state.definition.results[0].body, f.stock.body);
    }
    // A New body operation records the flag too.
    let mut group = ExtrudeGroupDraft::new("Extrusion operation", None);
    group.automatic = true;
    let id = group.id;
    let prepared = prepare_grouped_extrude(
        &f.session,
        f.session.revision(),
        actor(),
        &GroupedExtrudeRequest {
            extrusion: ExtrudeRequest::new(f.sketch, f.regions.clone(), Default::default()),
            group,
        },
    )
    .unwrap();
    f.session.commit(prepared.prepared).unwrap();
    let state = inspect_extrude_group(&reopen(f.session.document()), id).unwrap();
    assert!(state.definition.automatic);
    assert_eq!(state.definition.kind, ExtrudeGroupKind::NewBody);
}

#[test]
fn editing_an_automatic_operation_keeps_the_flag_only_while_its_operation_is_kept() {
    let mut f = fixture();
    let mut request = ExtrudeBooleanRequest::new(
        ExtrudeRequest::new(f.sketch, f.regions.clone(), Default::default()),
        Operation::Join,
        vec![f.stock.clone()],
        false,
    );
    let mut group = ExtrudeGroupDraft::new("Extrusion operation", None);
    group.automatic = true;
    let id = group.id;
    request.group = Some(group);
    let prepared =
        prepare_extrude_boolean(&f.session, f.session.revision(), actor(), &request).unwrap();
    f.session.commit(prepared.prepared).unwrap();
    let state = inspect_extrude_group(f.session.document(), id).unwrap();
    // A distance edit with the same operation stays automatic.
    let mut edit = ExtrudeGroupEdit::from_state(&state).unwrap();
    edit.parameters.extent = Extent::OneSided {
        end: Termination::Distance {
            distance: length(3.),
        },
    };
    let kept =
        prepare_extrude_group_edit(&f.session, f.session.revision(), actor(), &edit).unwrap();
    let kept = inspect_extrude_group(kept.prepared.document(), id).unwrap();
    assert!(kept.definition.automatic);
    // Choosing Cut explicitly ends the automatic record.
    edit.operation = Operation::Cut;
    let chosen =
        prepare_extrude_group_edit(&f.session, f.session.revision(), actor(), &edit).unwrap();
    f.session.commit(chosen.prepared).unwrap();
    let chosen = inspect_extrude_group(f.session.document(), id).unwrap();
    assert!(!chosen.definition.automatic);
    assert_eq!(chosen.operation().unwrap(), Operation::Cut);
}
