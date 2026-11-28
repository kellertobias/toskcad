//! CAD-432: joining panels groups them and the joint into one canonical
//! component. A further join adds to it; deleting the joint keeps the
//! panels and either keeps or dissolves the component; Undo, save/reopen
//! and the component model (scope, ownership, cut scope) stay consistent.
use multicad_core::{
    commands::{Actor, Channel, Prepared, Session},
    document::{Document, Entity, EntityData, EntityId},
    hierarchy::{self, Action, ContainerKind, Index},
    materials::*,
    model::{
        self, BodyOutput, GrainAxis, JointGroupOnDelete, PanelDefinition, PanelEdges,
        PanelJointDefinition, PanelJointKind, PanelJointPair, PanelPlane, PanelRequest,
        PanelSource, length,
    },
    project_file::{self, ProjectFile},
};
use multicad_model::Compiler;

fn actor() -> Actor {
    Actor::local(Channel::Ui)
}
fn checked(s: &mut Session, prepared: Prepared) {
    Compiler::default()
        .compile(prepared.document(), s.revision() + 1)
        .unwrap();
    s.commit(prepared).unwrap();
}
struct Fixture {
    s: Session,
    stock: EntityId,
    base: EntityId,
    left: EntityId,
    right: EntityId,
}
fn add_panel(
    s: &mut Session,
    stock: EntityId,
    name: &str,
    plane: PanelPlane,
    origin: [f64; 3],
    size: [f64; 2],
) -> EntityId {
    let request = PanelRequest::new(
        name,
        PanelDefinition {
            output: BodyOutput {
                body: EntityId::new(),
            },
            source: PanelSource::Rectangle {
                stock,
                plane,
                origin: origin.map(length),
                length: length(size[0]),
                width: length(size[1]),
            },
            grain: GrainAxis::U,
            edges: PanelEdges::default(),
            oversize: None,
        },
    );
    let p = model::prepare_panel(s, s.revision(), actor(), &request).unwrap();
    checked(s, p.prepared);
    request.body
}
/// A 600 × 400 base with an 18 mm side standing on each end.
fn fixture() -> Fixture {
    let mut s = Session::new(Document::new("Joined panels").unwrap()).unwrap();
    let birch = Entity::new(
        "Birch plywood",
        EntityData::Material {
            density_kg_m3: None,
            color_rgba: [214, 186, 140, 255],
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
    .entity(EntityId::new(), "Birch ply 18")
    .unwrap();
    let id = stock.id;
    s.commit(prepare_library(&s, 0, actor(), &[birch, stock]).unwrap())
        .unwrap();
    let base = add_panel(
        &mut s,
        id,
        "Panel 1",
        PanelPlane::Xy,
        [0., 0., 0.],
        [600., 400.],
    );
    let left = add_panel(
        &mut s,
        id,
        "Panel 2",
        PanelPlane::Yz,
        [0., 0., 0.],
        [400., 300.],
    );
    let right = add_panel(
        &mut s,
        id,
        "Panel 3",
        PanelPlane::Yz,
        [582., 0., 0.],
        [400., 300.],
    );
    Fixture {
        s,
        stock: id,
        base,
        left,
        right,
    }
}
fn join(s: &mut Session, name: &str, a: EntityId, b: EntityId) -> EntityId {
    let def = PanelJointDefinition {
        generated_by: None,
        kind: PanelJointKind::Finger {
            count: 5,
            fit: length(0.),
        },
        pairs: vec![PanelJointPair { a, b }],
    };
    let p = model::prepare_panel_joint(s, s.revision(), actor(), name, &def).unwrap();
    checked(s, p.model.prepared);
    p.joint
}
fn components(document: &Document) -> Vec<EntityId> {
    document
        .entities()
        .filter(|e| matches!(e.data, EntityData::Assembly { .. }))
        .map(|e| e.id)
        .collect()
}
/// Model-tree children: Part feature history is listed under its Part,
/// so features a panel's later cuts add to its container are left out.
fn children(document: &Document, container: Option<EntityId>) -> Vec<EntityId> {
    Index::new(document)
        .unwrap()
        .children(container)
        .unwrap()
        .into_iter()
        .filter(|id| !matches!(document.get(*id).unwrap().data, EntityData::Feature { .. }))
        .collect()
}
fn names(document: &Document, ids: &[EntityId]) -> Vec<String> {
    ids.iter()
        .map(|id| document.get(*id).unwrap().name.clone())
        .collect()
}

#[test]
fn joining_two_ungrouped_panels_makes_one_component_in_one_step() {
    let mut f = fixture();
    let before = f.s.document().clone();
    let joint = join(&mut f.s, "Joint", f.base, f.left);
    let d = f.s.document();
    let [component] = components(d)[..] else {
        panic!("one component")
    };
    assert_eq!(d.get(component).unwrap().name, "Joined panels 1");
    assert_eq!(children(d, Some(component)), vec![f.base, f.left, joint]);
    assert_eq!(
        names(d, &children(d, Some(component))),
        ["Panel 1", "Panel 2", "Joint"]
    );
    // The component is where the first panel was; the other panel is not in it.
    let root = children(d, None);
    assert!(
        root.contains(&component) && root.contains(&f.right),
        "{root:?}"
    );
    assert!(!root.contains(&f.base) && !root.contains(&joint));
    // It is an ordinary component: owning component, scope, cut scope.
    let index = Index::new(d).unwrap();
    assert_eq!(index.owning_component(f.left).unwrap(), Some(component));
    assert_eq!(index.default_cut_bodies(Some(component)).unwrap().len(), 2);
    let p = hierarchy::prepare(
        &f.s,
        f.s.revision(),
        actor(),
        Action::SetActive {
            component: Some(component),
        },
    )
    .unwrap();
    assert_eq!(
        p.prepared.document().hierarchy().active_component,
        Some(component)
    );
    // One Undo removes the joint, its cuts and the component together.
    f.s.undo(f.s.revision(), actor()).unwrap();
    assert_eq!(f.s.document(), &before);
    f.s.redo(f.s.revision(), actor()).unwrap();
    assert_eq!(components(f.s.document()), vec![component]);
}

#[test]
fn joining_a_further_panel_adds_it_and_its_joint_to_the_component() {
    let mut f = fixture();
    let first = join(&mut f.s, "Left joint", f.base, f.left);
    let second = join(&mut f.s, "Right joint", f.right, f.base);
    let d = f.s.document();
    let [component] = components(d)[..] else {
        panic!("still one component")
    };
    assert_eq!(
        children(d, Some(component)),
        vec![f.base, f.left, first, f.right, second]
    );
    // Save and reopen keep the component and its order.
    let mut bytes = vec![];
    project_file::write(
        &ProjectFile::new(d.clone()).unwrap(),
        &mut bytes,
        Default::default(),
    )
    .unwrap();
    let opened = project_file::read(
        &mut bytes.as_slice(),
        Default::default(),
        &project_file::supported_capabilities(),
    )
    .unwrap()
    .project;
    assert_eq!(&opened.document, d);
    assert_eq!(
        children(&opened.document, Some(component)),
        children(d, Some(component))
    );
    Compiler::default().compile(&opened.document, 0).unwrap();
}

#[test]
fn panels_of_another_component_stay_there() {
    let mut f = fixture();
    let p = hierarchy::prepare(
        &f.s,
        f.s.revision(),
        actor(),
        Action::Create {
            kind: ContainerKind::Component,
            name: "Drawer".into(),
            parent: None,
        },
    )
    .unwrap();
    let drawer = p.created.unwrap();
    f.s.commit(p.prepared).unwrap();
    let p = hierarchy::prepare(
        &f.s,
        f.s.revision(),
        actor(),
        Action::Move {
            entries: vec![f.left],
            parent: Some(drawer),
            before: None,
        },
    )
    .unwrap();
    f.s.commit(p.prepared).unwrap();
    // The base joins the Drawer the left side already belongs to.
    let joint = join(&mut f.s, "Joint", f.base, f.left);
    let d = f.s.document();
    assert_eq!(components(d), vec![drawer]);
    assert_eq!(children(d, Some(drawer)), vec![f.left, f.base, joint]);
    // A user-made second component keeps its own member.
    let p = hierarchy::prepare(
        &f.s,
        f.s.revision(),
        actor(),
        Action::Create {
            kind: ContainerKind::Component,
            name: "Side".into(),
            parent: None,
        },
    )
    .unwrap();
    let side = p.created.unwrap();
    f.s.commit(p.prepared).unwrap();
    let p = hierarchy::prepare(
        &f.s,
        f.s.revision(),
        actor(),
        Action::Move {
            entries: vec![f.right],
            parent: Some(side),
            before: None,
        },
    )
    .unwrap();
    f.s.commit(p.prepared).unwrap();
    let second = join(&mut f.s, "Right joint", f.base, f.right);
    let d = f.s.document();
    assert_eq!(children(d, Some(side)), vec![f.right]);
    assert_eq!(
        children(d, Some(drawer)),
        vec![f.left, f.base, joint, second]
    );
}

#[test]
fn deleting_the_joint_keeps_the_panels_and_by_default_the_component() {
    let mut f = fixture();
    let joint = join(&mut f.s, "Joint", f.base, f.left);
    let [component] = components(f.s.document())[..] else {
        panic!("one component")
    };
    let p = model::prepare_panel_joint_delete(&f.s, f.s.revision(), actor(), joint).unwrap();
    checked(&mut f.s, p);
    let d = f.s.document();
    assert!(d.get(joint).is_none());
    assert_eq!(children(d, Some(component)), vec![f.base, f.left]);
    // The panels lost their joint cuts.
    let m = Compiler::default().compile(d, 0).unwrap();
    let volume = |part| {
        m.solid(m.snapshot().bodies[&part])
            .unwrap()
            .volume_mm3()
            .unwrap()
    };
    assert!((volume(f.base) - 600. * 400. * 18.).abs() < 1e-3);
    // Save and reopen keep the now joint-less component.
    let mut bytes = vec![];
    project_file::write(
        &ProjectFile::new(d.clone()).unwrap(),
        &mut bytes,
        Default::default(),
    )
    .unwrap();
    let opened = project_file::read(
        &mut bytes.as_slice(),
        Default::default(),
        &project_file::supported_capabilities(),
    )
    .unwrap()
    .project;
    assert_eq!(
        children(&opened.document, Some(component)),
        vec![f.base, f.left]
    );
}

#[test]
fn deleting_the_last_joint_can_ungroup_the_panels_in_place() {
    let mut f = fixture();
    let first = join(&mut f.s, "Left joint", f.base, f.left);
    let second = join(&mut f.s, "Right joint", f.right, f.base);
    let [component] = components(f.s.document())[..] else {
        panic!("one component")
    };
    // Another joint still joins panels in it: ungrouping is refused, the
    // document is unchanged.
    let e = model::prepare_panel_joint_delete_with(
        &f.s,
        f.s.revision(),
        actor(),
        first,
        JointGroupOnDelete::Ungroup,
    )
    .unwrap_err();
    assert!(e.to_string().contains("Other joints"), "{e}");
    let p = model::prepare_panel_joint_delete(&f.s, f.s.revision(), actor(), first).unwrap();
    checked(&mut f.s, p);
    let before = f.s.document().clone();
    let p = model::prepare_panel_joint_delete_with(
        &f.s,
        f.s.revision(),
        actor(),
        second,
        JointGroupOnDelete::Ungroup,
    )
    .unwrap();
    checked(&mut f.s, p);
    let d = f.s.document();
    assert!(components(d).is_empty());
    assert!(d.get(component).is_none());
    let root = children(d, None);
    for panel in [f.base, f.left, f.right] {
        assert!(root.contains(&panel), "{panel} at root: {root:?}");
    }
    // One Undo brings the joint and the component back.
    f.s.undo(f.s.revision(), actor()).unwrap();
    assert_eq!(f.s.document(), &before);
    let _ = f.stock;
}
