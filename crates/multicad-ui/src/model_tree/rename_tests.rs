//! CAD-376 inline tree rename through real pointer and key events.
use super::*;
use multicad_core::{
    document::{Entity, Mutation},
    hierarchy::Hierarchy,
    materials::{Assignment, prepare_assignment},
    project_file::{self, ProjectFile},
};

fn fixture() -> (Session, EntityId, EntityId, EntityId) {
    let component = Entity::new("Cabinet", EntityData::Assembly { instances: vec![] });
    let group = Entity::new("Panels", EntityData::Group {});
    let material = Entity::new(
        "Birch",
        EntityData::Material {
            density_kg_m3: None,
            color_rgba: [210, 180, 140, 255],
            section_hatch: None,
        },
    );
    let body = Entity::new(
        "Panel 1",
        EntityData::Part {
            features: vec![],
            material: None,
            tags: BTreeSet::new(),
        },
    );
    let document = Document::new("Workshop")
        .unwrap()
        .edited(&[
            Mutation::Insert(component.clone()),
            Mutation::Insert(group.clone()),
            Mutation::Insert(material.clone()),
            Mutation::Insert(body.clone()),
            Mutation::SetHierarchy(Hierarchy {
                members: BTreeMap::from([
                    (component.id, vec![group.id]),
                    (group.id, vec![body.id]),
                ]),
                active_component: None,
                ..Default::default()
            }),
        ])
        .unwrap();
    (
        Session::new(document).unwrap(),
        component.id,
        body.id,
        material.id,
    )
}

struct Harness {
    context: egui::Context,
    tree: ModelTree,
    session: Session,
    selected: BTreeSet<EntityId>,
    time: f64,
    last: TreeEvents,
}
impl Harness {
    fn new(session: Session) -> Self {
        let context = egui::Context::default();
        context.enable_accesskit();
        let mut harness = Self {
            context,
            tree: ModelTree::default(),
            session,
            selected: BTreeSet::new(),
            time: 0.,
            last: TreeEvents::default(),
        };
        harness.frame(vec![]);
        harness
    }
    fn frame(&mut self, events: Vec<egui::Event>) {
        self.output(events).drop_without_applying_deltas();
    }
    fn output(&mut self, events: Vec<egui::Event>) -> egui::FullOutput {
        self.time += 0.05;
        self.context.run_ui(
            egui::RawInput {
                focused: true,
                time: Some(self.time),
                screen_rect: Some(egui::Rect::from_min_size(
                    egui::Pos2::ZERO,
                    egui::vec2(360., 480.),
                )),
                events,
                ..Default::default()
            },
            |ui| {
                let events = self.tree.show(
                    ui,
                    &mut self.session,
                    &self.selected,
                    TreeOptions::default(),
                );
                if let Some(next) = &events.selection {
                    self.selected = next.clone();
                }
                self.last = events;
            },
        )
    }
    fn bounds(&mut self, label: &str) -> egui::Rect {
        let output = self.output(vec![]);
        let bounds = output
            .platform_output
            .accesskit_update
            .as_ref()
            .unwrap()
            .nodes
            .iter()
            .find(|(_, node)| node.label() == Some(label))
            .and_then(|(_, node)| node.bounds());
        output.drop_without_applying_deltas();
        let bounds = bounds.unwrap_or_else(|| panic!("{label} row"));
        egui::Rect::from_min_max(
            egui::pos2(bounds.x0 as f32, bounds.y0 as f32),
            egui::pos2(bounds.x1 as f32, bounds.y1 as f32),
        )
    }
    fn click(&mut self, pos: egui::Pos2) {
        for pressed in [true, false] {
            self.frame(vec![
                egui::Event::PointerMoved(pos),
                egui::Event::PointerButton {
                    pos,
                    button: egui::PointerButton::Primary,
                    pressed,
                    modifiers: egui::Modifiers::NONE,
                },
            ]);
        }
    }
    fn double_click(&mut self, label: &str) {
        // A pause keeps egui from counting a previous gesture as a triple click.
        self.time += 1.;
        let pos = self.bounds(label).center();
        self.click(pos);
        self.click(pos);
        self.frame(vec![]);
    }
    /// Press and release a key; `last` keeps the events of the press frame.
    fn key(&mut self, key: egui::Key) {
        let mut pressed_events = TreeEvents::default();
        for pressed in [true, false] {
            self.frame(vec![egui::Event::Key {
                key,
                physical_key: None,
                pressed,
                repeat: false,
                modifiers: egui::Modifiers::NONE,
            }]);
            if pressed {
                pressed_events = std::mem::take(&mut self.last);
            }
        }
        self.last = pressed_events;
    }
    fn select_all_and_type(&mut self, text: &str) {
        self.frame(vec![egui::Event::Key {
            key: egui::Key::A,
            physical_key: None,
            pressed: true,
            repeat: false,
            modifiers: egui::Modifiers::COMMAND,
        }]);
        self.frame(vec![egui::Event::Text(text.into())]);
    }
    fn name(&self, id: EntityId) -> String {
        self.session.document().get(id).unwrap().name.clone()
    }
}

#[test]
fn double_click_renames_a_part_inline_with_enter_and_keeps_identity() {
    let (session, component, body, material) = fixture();
    let mut h = Harness::new(session);
    let before = h.session.document().clone();
    h.double_click("Panel 1");
    assert_eq!(
        h.selected,
        BTreeSet::from([body]),
        "the first click selects"
    );
    assert_eq!(h.tree.inline_renaming(), Some(body));
    assert!(
        h.context.text_edit_focused(),
        "the inline editor takes focus"
    );
    assert_eq!(h.session.revision(), 0, "editing alone changes nothing");
    // The editor starts with the whole name selected, so typing replaces it.
    h.frame(vec![egui::Event::Text("Left side".into())]);
    h.key(egui::Key::Enter);
    assert_eq!(h.tree.inline_renaming(), None);
    assert_eq!(h.session.revision(), 1);
    assert_eq!(h.name(body), "Left side");
    assert!(h.last.change.is_some());
    assert_eq!(h.session.activity().last().unwrap().label, "Rename part");
    // Same entity, data and containment: only the name differs.
    let renamed = h.session.document().get(body).unwrap();
    assert_eq!(renamed.data, before.get(body).unwrap().data);
    assert_eq!(h.session.document().hierarchy(), before.hierarchy());
    h.bounds("Left side");

    // Undo / Redo and save / reopen keep the name.
    h.session
        .undo(h.session.revision(), Actor::local(Channel::Ui))
        .unwrap();
    assert_eq!(h.name(body), "Panel 1");
    h.session
        .redo(h.session.revision(), Actor::local(Channel::Ui))
        .unwrap();
    assert_eq!(h.name(body), "Left side");
    let mut bytes = vec![];
    project_file::write(
        &ProjectFile::new(h.session.document().clone()).unwrap(),
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
    assert_eq!(opened.project.document.get(body).unwrap().name, "Left side");

    // A later material assignment does not overwrite the supplied name.
    let revision = h.session.revision();
    let prepared = prepare_assignment(
        &h.session,
        revision,
        Actor::local(Channel::Ui),
        &[body],
        &Assignment {
            material: Some(material),
            stock: None,
            appearance: None,
            hatch: None,
        },
    )
    .unwrap();
    h.session.commit(prepared).unwrap();
    assert_eq!(h.name(body), "Left side");

    // Components rename the same way; activation stays in the context menu.
    h.double_click("Cabinet");
    assert_eq!(h.tree.inline_renaming(), Some(component));
    h.select_all_and_type("Base cabinet");
    h.key(egui::Key::Enter);
    assert_eq!(h.name(component), "Base cabinet");
    assert_eq!(h.session.document().hierarchy().active_component, None);
}

#[test]
fn escape_cancels_and_empty_names_are_refused_in_place() {
    let (session, _, body, _) = fixture();
    let mut h = Harness::new(session);
    let before = h.session.document().clone();
    h.double_click("Panel 1");
    h.frame(vec![egui::Event::Text("Discarded".into())]);
    h.key(egui::Key::Escape);
    assert_eq!(h.tree.inline_renaming(), None);
    assert_eq!(h.session.document(), &before);
    assert_eq!(h.session.revision(), 0);

    h.double_click("Panel 1");
    h.key(egui::Key::Delete);
    h.frame(vec![egui::Event::Text("   ".into())]);
    h.key(egui::Key::Enter);
    let error = h.last.error.clone().expect("blank names are refused");
    assert!(error.contains("cannot be empty"), "{error}");
    assert_eq!(
        h.tree.inline_renaming(),
        Some(body),
        "the editor stays open"
    );
    assert_eq!(h.session.document(), &before);
    h.frame(vec![]);
    assert!(
        h.context.text_edit_focused(),
        "focus returns for a correction"
    );
    h.select_all_and_type("Shelf");
    h.key(egui::Key::Enter);
    assert_eq!(h.name(body), "Shelf");
    assert_eq!(h.session.revision(), 1);
}

#[test]
fn single_click_selection_and_caret_expansion_still_work() {
    let (session, component, body, _) = fixture();
    let mut h = Harness::new(session);
    let row = h.bounds("Panel 1");
    h.click(row.center());
    assert_eq!(h.selected, BTreeSet::from([body]));
    assert_eq!(h.tree.inline_renaming(), None, "one click never renames");
    // The caret left of the component name collapses its children.
    let cabinet = h.bounds("Cabinet");
    h.click(egui::pos2(cabinet.left() - 30., cabinet.center().y));
    h.frame(vec![]);
    assert!(
        !h.tree
            .rows(&h.session)
            .unwrap()
            .iter()
            .any(|row| row.id == body),
        "caret collapsed the component"
    );
    assert_eq!(h.tree.inline_renaming(), None);
    assert!(h.selected.contains(&body) && !h.selected.contains(&component));
    assert_eq!(h.session.revision(), 0);
}
