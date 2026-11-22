//! CAD-377 panel names generated from the chosen stock.
use super::*;
use multicad_core::{
    commands::{Actor, Channel, Request},
    document::{Entity, EntityData, Mutation},
    hierarchy::{self, Action},
    materials::{Grain, Stock, prepare_library},
    project_file::{self, ProjectFile},
};
use multicad_model::Compiler;

fn ui() -> Actor {
    Actor::local(Channel::Ui)
}
/// Birch material with stocks "MPX 6", "Ply 18" and an unnamed editor
/// placeholder "New panel stock".
fn session() -> (Session, [EntityId; 3]) {
    let mut s = Session::new(Document::new("Panel names").unwrap()).unwrap();
    let birch = Entity::new(
        "Birch plywood",
        EntityData::Material {
            density_kg_m3: None,
            color_rgba: [210, 180, 140, 255],
            section_hatch: None,
        },
    );
    let stock = |t: f64, name: &str| {
        LibraryDefinition::Stock(Stock {
            material: birch.id,
            specification: StockKind::Panel {
                thickness_mm: t,
                length_mm: 2500.,
                width_mm: 1250.,
                grain: Grain::AlongLength,
            },
        })
        .entity(EntityId::new(), name)
        .unwrap()
    };
    let stocks = [
        stock(6., "MPX 6"),
        stock(18., "Ply 18"),
        stock(12., "New panel stock"),
    ];
    let ids = [stocks[0].id, stocks[1].id, stocks[2].id];
    let mut entities = vec![birch];
    entities.extend(stocks);
    s.commit(prepare_library(&s, 0, ui(), &entities).unwrap())
        .unwrap();
    (s, ids)
}
fn create(s: &mut Session, form: &CabinetPanel) -> (EntityId, EntityId) {
    let prepared = form
        .request(s.document(), true)
        .unwrap()
        .prepare(s, ui())
        .unwrap();
    let output = (prepared.outputs[0].feature, prepared.outputs[0].body);
    s.commit(prepared.prepared).unwrap();
    output
}
fn new_panel(s: &Session, stock: EntityId) -> CabinetPanel {
    let mut form = CabinetPanel::generated_panel(s.document());
    form.select_stock(stock, s.document());
    form
}
fn name(s: &Session, id: EntityId) -> String {
    s.document().get(id).unwrap().name.clone()
}
fn assemblies(document: &Document) -> usize {
    document
        .entities()
        .filter(|e| matches!(e.data, EntityData::Assembly { .. }))
        .count()
}

#[test]
fn successive_panels_are_named_from_their_stock_with_numbered_suffixes() {
    let (mut s, [mpx, ply, _]) = session();
    let components = assemblies(s.document());
    let first = new_panel(&s, mpx);
    assert_eq!(first.name, "MPX 6 - 1");
    let (feature, body) = create(&mut s, &first);
    assert_eq!(name(&s, body), "MPX 6 - 1");
    assert_eq!(name(&s, feature), "MPX 6 - 1 feature");
    let second = new_panel(&s, mpx);
    assert_eq!(second.name, "MPX 6 - 2");
    let (_, second_body) = create(&mut s, &second);
    assert_eq!(name(&s, second_body), "MPX 6 - 2");
    assert_ne!(body, second_body);
    // Another stock has its own sequence.
    assert_eq!(new_panel(&s, ply).name, "Ply 18 - 1");
    // The name is the Part itself: no extra component holds it.
    assert_eq!(assemblies(s.document()), components);
}

#[test]
fn generated_names_skip_saved_duplicated_and_producer_names() {
    let (mut s, [mpx, ..]) = session();
    let taken = [
        Entity::new(
            "MPX 6 - 1",
            EntityData::Part {
                features: vec![],
                material: None,
                tags: Default::default(),
            },
        ),
        Entity::new("MPX 6 - 2 feature", EntityData::Group {}),
    ];
    let revision = s.revision();
    let prepared = s
        .prepare(Request {
            expected_revision: revision,
            actor: ui(),
            label: "Saved names".into(),
            edits: taken.map(Mutation::Insert).to_vec(),
        })
        .unwrap();
    s.commit(prepared).unwrap();
    // Reopen a saved copy: names already in the file still count.
    let mut bytes = vec![];
    project_file::write(
        &ProjectFile::new(s.document().clone()).unwrap(),
        &mut bytes,
        Default::default(),
    )
    .unwrap();
    let opened = project_file::read(
        &mut bytes.as_slice(),
        Default::default(),
        &project_file::supported_capabilities(),
    )
    .unwrap();
    let reopened = Session::new(opened.project.document).unwrap();
    assert_eq!(new_panel(&reopened, mpx).name, "MPX 6 - 3");
}

#[test]
fn stock_choice_follows_until_the_user_types_a_name() {
    let (s, [mpx, ply, placeholder]) = session();
    let mut form = new_panel(&s, mpx);
    form.select_stock(ply, s.document());
    assert_eq!(form.name, "Ply 18 - 1");
    // Placeholder stock names fall back to the material name.
    form.select_stock(placeholder, s.document());
    assert_eq!(form.name, "Birch plywood - 1");
    // A typed name is kept through later stock changes.
    let ctx = egui::Context::default();
    ctx.enable_accesskit();
    let rect = egui::Rect::from_min_size(egui::Pos2::ZERO, egui::vec2(360., 900.));
    let frame = |form: &mut CabinetPanel, events| {
        let output = ctx.run_ui(
            egui::RawInput {
                screen_rect: Some(rect),
                focused: true,
                events,
                ..Default::default()
            },
            |ui| {
                form.show(ui, s.document(), true, PreviewState::Idle);
            },
        );
        let field = output
            .platform_output
            .accesskit_update
            .as_ref()
            .unwrap()
            .nodes
            .iter()
            .find(|(_, node)| node.value() == Some(form.name.as_str()))
            .and_then(|(_, node)| node.bounds());
        output.drop_without_applying_deltas();
        field
    };
    let bounds = frame(&mut form, vec![]).expect("name field");
    let pos = egui::pos2(
        ((bounds.x0 + bounds.x1) * 0.5) as f32,
        ((bounds.y0 + bounds.y1) * 0.5) as f32,
    );
    for pressed in [true, false] {
        frame(
            &mut form,
            vec![
                egui::Event::PointerMoved(pos),
                egui::Event::PointerButton {
                    pos,
                    button: egui::PointerButton::Primary,
                    pressed,
                    modifiers: egui::Modifiers::NONE,
                },
            ],
        );
    }
    frame(&mut form, vec![egui::Event::Text(" shelf".into())]);
    assert!(form.name.contains(" shelf"), "{}", form.name);
    assert!(!form.name_generated, "typing ends generated naming");
    let typed = form.name.clone();
    form.select_stock(mpx, s.document());
    assert_eq!(form.name, typed);
    // Without stock the existing readable fallback stays.
    let mut empty = Session::new(Document::new("No stock").unwrap()).unwrap();
    assert_eq!(
        CabinetPanel::generated_panel(empty.document()).name,
        "Panel 1"
    );
    let _ = &mut empty;
}

#[test]
fn a_renamed_panel_keeps_its_name_through_creation_stock_edit_rebuild_and_reopen() {
    let (mut s, [mpx, ply, _]) = session();
    let first = new_panel(&s, mpx);
    let (feature, body) = create(&mut s, &first);
    let rename = hierarchy::prepare(
        &s,
        s.revision(),
        ui(),
        Action::Rename {
            id: body,
            name: "Back panel".into(),
        },
    )
    .unwrap();
    s.commit(rename.prepared).unwrap();
    // A later creation does not reuse or touch the renamed Part.
    let next = new_panel(&s, mpx);
    let (_, next) = create(&mut s, &next);
    // The renamed Part's producer still carries "MPX 6 - 1 feature".
    assert_eq!(name(&s, next), "MPX 6 - 2");
    assert_eq!(name(&s, body), "Back panel");
    // Editing the stock of the renamed panel keeps the user's name and id.
    let mut edit = CabinetPanel::existing(s.document(), feature).unwrap();
    assert_eq!(edit.name, "Back panel");
    edit.select_stock(ply, s.document());
    assert_eq!(edit.name, "Back panel");
    let prepared = edit
        .request(s.document(), true)
        .unwrap()
        .prepare(&s, ui())
        .unwrap();
    assert_eq!(prepared.outputs[0].body, body);
    s.commit(prepared.prepared).unwrap();
    assert_eq!(name(&s, body), "Back panel");
    // Rebuild and reopen keep it.
    Compiler::default()
        .compile(s.document(), s.revision())
        .unwrap();
    let mut bytes = vec![];
    project_file::write(
        &ProjectFile::new(s.document().clone()).unwrap(),
        &mut bytes,
        Default::default(),
    )
    .unwrap();
    let opened = project_file::read(
        &mut bytes.as_slice(),
        Default::default(),
        &project_file::supported_capabilities(),
    )
    .unwrap();
    assert_eq!(
        opened.project.document.get(body).unwrap().name,
        "Back panel"
    );
}
