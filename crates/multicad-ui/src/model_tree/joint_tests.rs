//! CAD-432: the model tree shows joined panels as one component with the
//! panels and the joint as children, selects each child on its own and
//! deletes the joint from its context menu, keeping the panels and keeping
//! or dissolving the group.
use super::*;
use multicad_core::{
    document::Entity,
    materials::{Grain, LibraryDefinition, Stock, StockKind, prepare_library},
    model::{self, PanelJointDefinition, PanelJointKind, PanelJointPair, PanelPlane},
};
use std::time::{Duration, Instant};

fn ui_actor() -> Actor {
    Actor::local(Channel::Ui)
}
fn panel(
    s: &mut Session,
    stock: EntityId,
    name: &str,
    plane: PanelPlane,
    origin: f64,
    size: [f64; 2],
) -> EntityId {
    let request = model::PanelRequest::new(
        name,
        model::PanelDefinition {
            output: model::BodyOutput {
                body: EntityId::new(),
            },
            source: model::PanelSource::Rectangle {
                stock,
                plane,
                origin: [origin, 0., 0.].map(model::length),
                length: model::length(size[0]),
                width: model::length(size[1]),
            },
            grain: model::GrainAxis::U,
            edges: Default::default(),
            oversize: None,
        },
    );
    let p = model::prepare_panel(s, s.revision(), ui_actor(), &request).unwrap();
    s.commit(p.prepared).unwrap();
    request.body
}
/// Panel 1 (base) and Panel 2 (side) joined by "Corner joint"; Panel 3 stands
/// alone at the other end.
fn fixture() -> (Session, [EntityId; 3], EntityId) {
    let mut s = Session::new(Document::new("Workshop").unwrap()).unwrap();
    let birch = Entity::new(
        "Birch",
        EntityData::Material {
            density_kg_m3: None,
            color_rgba: [210, 180, 140, 255],
            section_hatch: None,
        },
    );
    let stock = LibraryDefinition::Stock(Stock {
        material: birch.id,
        specification: StockKind::Panel {
            thickness_mm: 18.,
            length_mm: 2500.,
            width_mm: 1250.,
            grain: Grain::AlongLength,
        },
    })
    .entity(EntityId::new(), "Ply 18")
    .unwrap();
    let id = stock.id;
    s.commit(prepare_library(&s, 0, ui_actor(), &[birch, stock]).unwrap())
        .unwrap();
    let base = panel(&mut s, id, "Panel 1", PanelPlane::Xy, 0., [600., 400.]);
    let side = panel(&mut s, id, "Panel 2", PanelPlane::Yz, 0., [400., 300.]);
    let other = panel(&mut s, id, "Panel 3", PanelPlane::Yz, 582., [400., 300.]);
    let def = PanelJointDefinition {
        generated_by: None,
        kind: PanelJointKind::Finger {
            count: 5,
            fit: model::length(0.),
        },
        pairs: vec![PanelJointPair { a: base, b: side }],
    };
    let p = model::prepare_panel_joint(&s, s.revision(), ui_actor(), "Corner joint", &def).unwrap();
    s.commit(p.model.prepared).unwrap();
    (s, [base, side, other], p.joint)
}

struct Harness {
    context: egui::Context,
    tree: ModelTree,
    session: Session,
    selected: BTreeSet<EntityId>,
    time: f64,
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
            time: 0.,
        }
    }
    fn output(&mut self, events: Vec<egui::Event>) -> egui::FullOutput {
        self.time += 0.05;
        self.context.run_ui(
            egui::RawInput {
                focused: true,
                time: Some(self.time),
                screen_rect: Some(egui::Rect::from_min_size(
                    egui::Pos2::ZERO,
                    egui::vec2(420., 640.),
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
                assert!(events.error.is_none(), "{:?}", events.error);
                if let Some(next) = events.selection {
                    self.selected = next;
                }
            },
        )
    }
    fn bounds(&mut self, label: &str) -> Option<egui::Rect> {
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
        bounds.map(|b| {
            egui::Rect::from_min_max(
                egui::pos2(b.x0 as f32, b.y0 as f32),
                egui::pos2(b.x1 as f32, b.y1 as f32),
            )
        })
    }
    /// The control labelled `label`, waiting a bounded time for it to appear.
    fn wait_for(&mut self, label: &str) -> egui::Pos2 {
        let deadline = Instant::now() + Duration::from_secs(10);
        loop {
            if let Some(rect) = self.bounds(label) {
                return rect.center();
            }
            assert!(Instant::now() < deadline, "timed out waiting for {label}");
        }
    }
    fn click(&mut self, label: &str, button: egui::PointerButton) {
        self.time += 1.;
        let pos = self.wait_for(label);
        for pressed in [true, false] {
            self.output(vec![
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
    fn rows(&mut self) -> Vec<(usize, String, NodeKind)> {
        self.tree
            .rows(&self.session)
            .unwrap()
            .iter()
            .map(|r| (r.depth, r.name.clone(), r.kind))
            .collect()
    }
}

#[test]
fn joined_panels_show_as_one_component_with_selectable_children() {
    let (session, [base, side, other], joint) = fixture();
    let mut h = Harness::new(session);
    let rows = h.rows();
    let component = rows
        .iter()
        .position(|(_, name, kind)| name == "Joined panels 1" && *kind == NodeKind::Component)
        .unwrap_or_else(|| panic!("component row in {rows:?}"));
    let depth = rows[component].0;
    assert_eq!(
        rows[component + 1..component + 4],
        [
            (depth + 1, "Panel 1".to_owned(), NodeKind::Body),
            (depth + 1, "Panel 2".to_owned(), NodeKind::Body),
            (depth + 1, "Corner joint".to_owned(), NodeKind::Joint),
        ]
    );
    // The lone panel stays outside the component.
    assert!(
        rows.iter()
            .any(|(d, name, _)| *d == depth && name == "Panel 3")
    );
    // Selection works per child.
    h.click("Corner joint", egui::PointerButton::Primary);
    assert_eq!(h.selected, BTreeSet::from([joint]));
    h.click("Panel 2", egui::PointerButton::Primary);
    assert_eq!(h.selected, BTreeSet::from([side]));
    h.click("Panel 1", egui::PointerButton::Primary);
    assert_eq!(h.selected, BTreeSet::from([base]));
    let _ = other;
}

#[test]
fn deleting_a_joint_from_the_tree_keeps_or_ungroups_the_panels() {
    let (session, [base, side, _], joint) = fixture();
    let after_join = session.document().clone();
    let mut h = Harness::new(session);
    let revision = h.session.revision();
    h.click("Corner joint", egui::PointerButton::Secondary);
    h.click(JOINT_DELETE_KEEP, egui::PointerButton::Primary);
    assert_eq!(h.session.revision(), revision + 1);
    let d = h.session.document();
    assert!(d.get(joint).is_none());
    let index = Index::new(d).unwrap();
    let component = index.parent(base).unwrap().expect("still grouped");
    assert_eq!(index.children(Some(component)).unwrap(), vec![base, side]);
    assert!(
        h.rows()
            .iter()
            .any(|(_, name, kind)| name == "Joined panels 1" && *kind == NodeKind::Component)
    );

    // From the joined state again, ungrouping dissolves the component.
    let mut h = Harness::new(Session::new(after_join).unwrap());
    h.click("Corner joint", egui::PointerButton::Secondary);
    h.click(JOINT_DELETE_UNGROUP, egui::PointerButton::Primary);
    let d = h.session.document();
    assert!(d.get(joint).is_none());
    assert!(
        !d.entities()
            .any(|e| matches!(e.data, EntityData::Assembly { .. }))
    );
    let index = Index::new(d).unwrap();
    assert_eq!(index.parent(base).unwrap(), None);
    assert_eq!(index.parent(side).unwrap(), None);
    let rows = h.rows();
    assert!(
        rows.iter().any(|(_, name, _)| name == "Panel 1"),
        "{rows:?}"
    );
    assert!(
        !rows.iter().any(|(_, name, _)| name == "Corner joint"),
        "{rows:?}"
    );
}
