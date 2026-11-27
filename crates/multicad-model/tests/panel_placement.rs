//! CAD-430: a panel's position (origin corner X/Y/Z) and orientation (face
//! direction preset plus a turn about it) read where the Part is now and edit
//! as one canonical transaction each. Placements that the panel feature can
//! hold stay there; others are carried by one trailing Move / rotate that
//! later edits change in place. Values reflect Move / rotate made elsewhere,
//! survive Undo/Redo and save/reopen, and joints follow.
use multicad_core::{
    commands::{Actor, Channel, Prepared, Session},
    document::{Document, Entity, EntityData, EntityId},
    materials::*,
    model::{
        self, BodyInput, BodyOutput, GrainAxis, NumericTransform, PanelDefinition, PanelEdges,
        PanelFacing, PanelJointDefinition, PanelJointKind, PanelJointPair, PanelPlacement,
        PanelPlacementRequest, PanelPlane, PanelRequest, PanelSource, TransformRequest, angle,
        length,
    },
    parameters,
    project_file::{self, ProjectFile},
};
use multicad_model::{Compiler, Model};

fn actor() -> Actor {
    Actor::local(Channel::Ui)
}
fn checked(s: &mut Session, prepared: Prepared) {
    Compiler::default()
        .compile(prepared.document(), s.revision() + 1)
        .unwrap();
    s.commit(prepared).unwrap();
}
fn compile(document: &Document) -> Model {
    Compiler::default().compile(document, 0).unwrap()
}
fn bounds(document: &Document, part: EntityId) -> [[i64; 3]; 2] {
    let m = compile(document);
    let b = m
        .solid(m.snapshot().bodies[&part])
        .unwrap()
        .bounds()
        .unwrap();
    let r = |p: [f64; 3]| p.map(|x| (x * 1000.).round() as i64);
    [r(b.min), r(b.max)]
}
fn mm(p: [f64; 3]) -> [i64; 3] {
    p.map(|x| (x * 1000.).round() as i64)
}

struct Fixture {
    s: Session,
    stock: EntityId,
    feature: EntityId,
    part: EntityId,
}
fn new_panel(
    s: &mut Session,
    stock: EntityId,
    name: &str,
    plane: PanelPlane,
    origin: [f64; 3],
    size: [f64; 2],
) -> (EntityId, EntityId) {
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
    (request.feature, request.body)
}
fn fixture() -> Fixture {
    let mut s = Session::new(Document::new("Panel placement").unwrap()).unwrap();
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
    let (feature, part) = new_panel(
        &mut s,
        id,
        "Panel 1",
        PanelPlane::Xy,
        [0., 0., 0.],
        [600., 400.],
    );
    Fixture {
        s,
        stock: id,
        feature,
        part,
    }
}
fn placement(s: &Session, part: EntityId) -> PanelPlacement {
    model::panel_placement(s.document(), s.parameters(), part).unwrap()
}
fn features(s: &Session, part: EntityId) -> Vec<EntityId> {
    let Some(EntityData::Part { features, .. }) = s.document().get(part).map(|e| &e.data) else {
        panic!("not a Part")
    };
    features.clone()
}
fn place(s: &mut Session, feature: EntityId, position: [f64; 3], facing: PanelFacing, deg: f64) {
    let panel = model::read_panel(s.document(), feature).unwrap();
    let request = PanelPlacementRequest {
        position: position.map(length),
        facing,
        rotation: angle(deg),
    };
    let p = model::prepare_panel_placement(s, s.revision(), actor(), feature, &panel, &request)
        .unwrap();
    checked(s, p.prepared);
}

#[test]
fn a_new_panel_reads_its_feature_placement() {
    let f = fixture();
    assert_eq!(
        placement(&f.s, f.part),
        PanelPlacement {
            origin_mm: [0.; 3],
            facing: Some(PanelFacing::PosZ),
            rotation_deg: 0.,
        }
    );
    assert_eq!(
        bounds(f.s.document(), f.part),
        [[0; 3], mm([600., 400., 18.])]
    );
}

#[test]
fn position_and_principal_presets_edit_the_panel_feature_in_one_step() {
    let mut f = fixture();
    let before = f.s.document().clone();
    place(&mut f.s, f.feature, [100., 50., 25.], PanelFacing::PosZ, 0.);
    // Stays in the panel feature: no new history.
    assert_eq!(features(&f.s, f.part), vec![f.feature]);
    let PanelSource::Rectangle { origin, plane, .. } =
        model::read_panel(f.s.document(), f.feature).unwrap().source
    else {
        unreachable!()
    };
    assert_eq!(origin, [100., 50., 25.].map(length));
    assert_eq!(plane, PanelPlane::Xy);
    assert_eq!(
        bounds(f.s.document(), f.part),
        [mm([100., 50., 25.]), mm([700., 450., 43.])]
    );
    // +X is the YZ panel plane: thickness toward +X.
    place(&mut f.s, f.feature, [100., 50., 25.], PanelFacing::PosX, 0.);
    assert_eq!(features(&f.s, f.part), vec![f.feature]);
    let PanelSource::Rectangle { plane, .. } =
        model::read_panel(f.s.document(), f.feature).unwrap().source
    else {
        unreachable!()
    };
    assert_eq!(plane, PanelPlane::Yz);
    assert_eq!(
        bounds(f.s.document(), f.part),
        [mm([100., 50., 25.]), mm([118., 650., 425.])]
    );
    assert_eq!(placement(&f.s, f.part).facing, Some(PanelFacing::PosX));
    // One Undo step per edit.
    f.s.undo(f.s.revision(), actor()).unwrap();
    f.s.undo(f.s.revision(), actor()).unwrap();
    assert_eq!(f.s.document(), &before);
}

#[test]
fn other_presets_and_turns_use_one_trailing_move_that_later_edits_change() {
    let mut f = fixture();
    // Face down: the panel hangs below its origin plane, spanning −Y.
    place(&mut f.s, f.feature, [0., 0., 0.], PanelFacing::NegZ, 0.);
    let history = features(&f.s, f.part);
    assert_eq!(history.len(), 2);
    assert_eq!(
        f.s.document().get(history[1]).unwrap().name,
        "Panel placement"
    );
    assert_eq!(
        bounds(f.s.document(), f.part),
        [mm([0., -400., -18.]), mm([600., 0., 0.])]
    );
    assert_eq!(
        placement(&f.s, f.part),
        PanelPlacement {
            origin_mm: [0.; 3],
            facing: Some(PanelFacing::NegZ),
            rotation_deg: 0.,
        }
    );
    // A turn about the face direction edits the same Move / rotate.
    let undo_to = f.s.document().clone();
    place(&mut f.s, f.feature, [10., 20., 30.], PanelFacing::PosZ, 90.);
    assert_eq!(features(&f.s, f.part), history);
    let read = placement(&f.s, f.part);
    assert_eq!(read.facing, Some(PanelFacing::PosZ));
    assert!((read.rotation_deg - 90.).abs() < 1e-9, "{read:?}");
    assert_eq!(mm(read.origin_mm), mm([10., 20., 30.]));
    // Length now runs along +Y, width along −X.
    assert_eq!(
        bounds(f.s.document(), f.part),
        [mm([-390., 20., 30.]), mm([10., 620., 48.])]
    );
    // A free turn reads back exactly.
    place(&mut f.s, f.feature, [10., 20., 30.], PanelFacing::PosY, 30.);
    let read = placement(&f.s, f.part);
    assert_eq!(read.facing, Some(PanelFacing::PosY));
    assert!((read.rotation_deg - 30.).abs() < 1e-9, "{read:?}");
    assert_eq!(features(&f.s, f.part), history);
    f.s.undo(f.s.revision(), actor()).unwrap();
    f.s.undo(f.s.revision(), actor()).unwrap();
    assert_eq!(f.s.document(), &undo_to);
}

#[test]
fn values_follow_a_move_made_elsewhere_and_survive_save_and_reopen() {
    let mut f = fixture();
    // A Move / rotate made by the Move tool or MCP.
    let request = TransformRequest::new(
        BodyInput {
            body: f.part,
            producer: f.feature,
        },
        NumericTransform {
            translation: [length(5.), length(6.), length(7.)],
            axis: [0., 0., 1.],
            angle: angle(180.),
            ..Default::default()
        },
    );
    let move_id = request.feature;
    let p = model::prepare_transform(&f.s, f.s.revision(), actor(), &request).unwrap();
    checked(&mut f.s, p.prepared);
    let read = placement(&f.s, f.part);
    assert_eq!(mm(read.origin_mm), mm([5., 6., 7.]));
    assert_eq!(read.facing, Some(PanelFacing::PosZ));
    assert!((read.rotation_deg - 180.).abs() < 1e-9, "{read:?}");
    // A position edit now edits that literal trailing Move / rotate.
    place(
        &mut f.s,
        f.feature,
        [50., 60., 70.],
        PanelFacing::PosZ,
        180.,
    );
    assert_eq!(features(&f.s, f.part), vec![f.feature, move_id]);
    assert_eq!(
        bounds(f.s.document(), f.part),
        [mm([-550., -340., 70.]), mm([50., 60., 88.])]
    );
    let mut bytes = vec![];
    project_file::write(
        &ProjectFile::new(f.s.document().clone()).unwrap(),
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
    assert_eq!(&opened.document, f.s.document());
    let parameters = parameters::evaluate(&opened.document).unwrap();
    assert_eq!(
        model::panel_placement(&opened.document, &parameters, f.part).unwrap(),
        placement(&f.s, f.part)
    );
}

#[test]
fn a_parameter_bound_move_is_kept_and_a_new_one_appended() {
    let mut f = fixture();
    let shift = Entity::new(
        "shift",
        EntityData::Parameter {
            expression: "10".into(),
            unit: "mm".into(),
            exposed: true,
            inputs: vec![],
        },
    );
    let prepared =
        f.s.prepare(multicad_core::commands::Request {
            expected_revision: f.s.revision(),
            actor: actor(),
            label: "Add parameter".into(),
            edits: vec![multicad_core::document::Mutation::Insert(shift)],
        })
        .unwrap();
    f.s.commit(prepared).unwrap();
    let expression =
        multicad_core::parameters::BoundExpression::bind(f.s.document(), "=$shift", "mm", None)
            .unwrap();
    let request = TransformRequest::new(
        BodyInput {
            body: f.part,
            producer: f.feature,
        },
        NumericTransform {
            translation: [
                multicad_core::constraint_sketch::Dimension::Expression { expression },
                length(0.),
                length(0.),
            ],
            ..Default::default()
        },
    );
    let p = model::prepare_transform(&f.s, f.s.revision(), actor(), &request).unwrap();
    checked(&mut f.s, p.prepared);
    assert_eq!(mm(placement(&f.s, f.part).origin_mm), mm([10., 0., 0.]));
    // The bound Move / rotate stays as authored; a literal one is appended.
    place(&mut f.s, f.feature, [0., 0., 0.], PanelFacing::PosZ, 0.);
    let history = features(&f.s, f.part);
    assert_eq!(history.len(), 3);
    assert_eq!(history[1], request.feature);
    assert_eq!(
        model::read_transform(f.s.document(), request.feature)
            .unwrap()
            .parameters,
        request.parameters
    );
    assert_eq!(mm(placement(&f.s, f.part).origin_mm), [0; 3]);
}
#[test]
fn a_joint_follows_a_placement_edit() {
    let mut f = fixture();
    let (_, upright) = new_panel(
        &mut f.s,
        f.stock,
        "Panel 2",
        PanelPlane::Yz,
        [291., 0., 0.],
        [400., 300.],
    );
    // A 400 × 300 upright standing on the centre line of the base.
    let def = PanelJointDefinition {
        generated_by: None,
        kind: PanelJointKind::Finger {
            count: 5,
            fit: length(0.),
        },
        pairs: vec![PanelJointPair {
            a: f.part,
            b: upright,
        }],
    };
    let p = model::prepare_panel_joint(&f.s, f.s.revision(), actor(), "Joint", &def).unwrap();
    checked(&mut f.s, p.model.prepared);
    let volume = |s: &Session| {
        let m = compile(s.document());
        m.solid(m.snapshot().bodies[&f.part])
            .unwrap()
            .volume_mm3()
            .unwrap()
    };
    let joined = volume(&f.s);
    // Slide the base along X with the inspector: the joint still fits.
    place(&mut f.s, f.feature, [-100., 0., 0.], PanelFacing::PosZ, 0.);
    // The joint cut moved with the base: the same material is removed.
    assert!((volume(&f.s) - joined).abs() < 1e-3);
    assert!(joined < 600. * 400. * 18. - 1.);
    // Lifting the base off the upright is refused with the joint named.
    let panel = model::read_panel(f.s.document(), f.feature).unwrap();
    let lifted = model::prepare_panel_placement(
        &f.s,
        f.s.revision(),
        actor(),
        f.feature,
        &panel,
        &PanelPlacementRequest {
            position: [-100., 0., 500.].map(length),
            facing: PanelFacing::PosZ,
            rotation: angle(0.),
        },
    )
    .unwrap();
    let Err(e) = Compiler::default().compile(lifted.prepared.document(), 0) else {
        panic!("a placement that separates joined panels must not compile")
    };
    assert!(e.to_string().contains("Joint"), "{e}");
}
