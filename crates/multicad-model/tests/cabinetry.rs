//! Cabinetry panels: a parametric carcass generating named panel Parts with
//! finished and core sizes, parameter edits and Undo, shelf and back
//! changes, explicit banding and grain on rectangle panels, sketch-region
//! panels, face-to-face stacks, and save/reopen. Every size is measured from
//! the native solids.
use multicad_core::{
    commands::{Actor, Channel, Prepared, Request, Session},
    constraint_sketch::{Curve, CurveKind, Dimension, source::Source},
    document::{Document, Entity, EntityData, EntityId, Mutation},
    manufacturing::PlanarEdge,
    materials::*,
    model::{
        self, Back, Band, BodyOutput, CarcassDefinition, CarcassRequest, Construction, GrainAxis,
        MachiningKind, PanelDefinition, PanelEdges, PanelJointKind, PanelPlane, PanelRequest,
        PanelRole, PanelSource, ProfileSource, ShelfLayout, Shelves, angle, length,
    },
    parameters::BoundExpression,
    project_file::{self, ProjectFile},
    references::{self, AnalyticGeometry, Attachment, GeometryIntent, Resolution},
    sketch::Datum,
    sketch_regions,
};
use multicad_kernel::{GeometryPolicy, MeshOptions, PlaneFrame, Solid};
use multicad_model::{Compiler, CutRecord, GeometryBackend, Model};

fn actor() -> Actor {
    Actor::local(Channel::Ui)
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() <= 1e-6, "{a} != {b}");
}
fn same_cut_record(actual: &CutRecord, expected: &CutRecord) {
    // Native integration can differ in its last bits after recompilation.
    near(actual.volume_mm3, expected.volume_mm3);
    let metadata = CutRecord {
        volume_mm3: expected.volume_mm3,
        ..actual.clone()
    };
    assert_eq!(&metadata, expected);
}
fn same_cut_records(actual: &[CutRecord], expected: &[CutRecord]) {
    assert_eq!(actual.len(), expected.len());
    for (actual, expected) in actual.iter().zip(expected) {
        same_cut_record(actual, expected);
    }
}
fn pair(a: [f64; 2], b: [f64; 2]) {
    near(a[0], b[0]);
    near(a[1], b[1]);
}
fn compile(s: &Session) -> Model {
    Compiler::default()
        .compile(s.document(), s.revision())
        .unwrap()
}
/// Evaluate the candidate natively, then commit it.
fn checked(s: &mut Session, prepared: Prepared) {
    Compiler::default()
        .compile(prepared.document(), s.revision() + 1)
        .unwrap();
    s.commit(prepared).unwrap();
}
fn insert(s: &mut Session, edits: Vec<Mutation>) {
    s.apply(Request {
        expected_revision: s.revision(),
        actor: actor(),
        label: "Fixture".into(),
        edits,
    })
    .unwrap();
}
fn record<'a>(list: &'a [CutRecord], name: &str) -> &'a CutRecord {
    list.iter()
        .find(|r| r.name == name)
        .unwrap_or_else(|| panic!("no cut record {name}"))
}

struct Fixture {
    s: Session,
    ply: EntityId,
    ply_material: EntityId,
    hdf: EntityId,
    width: EntityId,
}
fn material(name: &str, color: [u8; 4]) -> Entity {
    Entity::new(
        name,
        EntityData::Material {
            density_kg_m3: Some(650.),
            color_rgba: color,
            section_hatch: None,
        },
    )
}
fn stock(material: EntityId, thickness: f64, grain: Grain, name: &str) -> Entity {
    LibraryDefinition::Stock(Stock {
        material,
        specification: StockKind::Panel {
            thickness_mm: thickness,
            length_mm: 2500.,
            width_mm: 1250.,
            grain,
        },
    })
    .entity(EntityId::new(), name)
    .unwrap()
}
fn fixture() -> Fixture {
    let mut s = Session::new(Document::new("Cabinetry").unwrap()).unwrap();
    let birch = material("Birch plywood", [214, 186, 140, 255]);
    let fibre = material("HDF", [120, 96, 70, 255]);
    let ply = stock(birch.id, 18., Grain::AlongLength, "Birch ply 18");
    let hdf = stock(fibre.id, 6., Grain::None, "HDF 6");
    let mut width = Entity::new(
        "cab_width",
        EntityData::Parameter {
            expression: "600".into(),
            unit: "mm".into(),
            exposed: true,
            inputs: vec![],
        },
    );
    width.description = "Outside cabinet width".into();
    let ids = (ply.id, birch.id, hdf.id, width.id);
    s.commit(prepare_library(&s, 0, actor(), &[birch, fibre]).unwrap())
        .unwrap();
    s.commit(prepare_library(&s, s.revision(), actor(), &[ply, hdf]).unwrap())
        .unwrap();
    insert(&mut s, vec![Mutation::Insert(width)]);
    Fixture {
        s,
        ply: ids.0,
        ply_material: ids.1,
        hdf: ids.2,
        width: ids.3,
    }
}
fn set_parameter(s: &mut Session, id: EntityId, expression: &str) {
    let mut entity = s.document().get(id).unwrap().clone();
    let EntityData::Parameter { expression: e, .. } = &mut entity.data else {
        unreachable!()
    };
    *e = expression.into();
    let prepared = s
        .prepare(Request {
            expected_revision: s.revision(),
            actor: actor(),
            label: "Edit parameter".into(),
            edits: vec![Mutation::Replace(entity)],
        })
        .unwrap();
    checked(s, prepared);
}
fn band(thickness: f64, allowance: f64) -> Band {
    Band {
        thickness: length(thickness),
        allowance: Some(length(allowance)),
    }
}
fn cabinet(f: &Fixture) -> CarcassDefinition {
    CarcassDefinition {
        cabinet: None,
        origin: [length(0.), length(0.), length(0.)],
        rotation: angle(0.),
        width: Dimension::Expression {
            expression: BoundExpression::bind(f.s.document(), "=$cab_width", "mm", None).unwrap(),
        },
        height: length(720.),
        depth: length(560.),
        stock: f.ply,
        construction: Construction::SidesFullHeight,
        back: Back::Inset {
            stock: f.hdf,
            recess: length(10.),
        },
        shelves: Shelves {
            layout: ShelfLayout::Even { count: 2 },
            setback: length(20.),
            clearance: length(1.),
        },
        front_band: Some(band(2., 0.5)),
        group: None,
    }
}
fn bounds(m: &Model, s: &Session, part: EntityId) -> multicad_kernel::Bounds {
    let head = m.snapshot().bodies[&part];
    let _ = s;
    m.solid(head).unwrap().bounds().unwrap()
}

fn explicit_cabinet(f: &Fixture) -> CarcassDefinition {
    let mut definition = cabinet(f);
    definition.width = length(800.);
    definition.height = length(2000.);
    definition.depth = length(600.);
    definition.cabinet = Some(model::CabinetOptions {
        toe_clearance: length(100.),
        shelf_stock: f.ply,
        shelves: vec![
            model::CabinetShelf {
                key: 0,
                underside: length(500.),
                depth: length(420.),
                stock: None,
                support: model::ShelfSupport::Fixed,
            },
            model::CabinetShelf {
                key: 3,
                underside: length(1053.5),
                depth: length(350.),
                stock: Some(f.hdf),
                support: model::ShelfSupport::System32 { level: 28 },
            },
        ],
        joinery: Default::default(),
        pins: Default::default(),
    });
    definition
}

#[test]
fn explicit_cabinet_floor_datum_independent_shelves_and_stable_rows_are_native() {
    let mut f = fixture();
    let definition = explicit_cabinet(&f);
    let created = model::prepare_carcass(
        &f.s,
        f.s.revision(),
        actor(),
        &CarcassRequest::new("Tall cabinet", definition.clone()),
    )
    .unwrap();
    let carcass = created.carcass;
    checked(&mut f.s, created.model.prepared);
    let native = compile(&f.s);
    let cuts = native.cut_list(f.s.document()).unwrap();
    let bottom = record(&cuts, "Tall cabinet bottom");
    near(bounds(&native, &f.s, bottom.part).min[2], 100.);
    let side = record(&cuts, "Tall cabinet left side");
    pair(side.cut.finished_mm, [2000., 600.]);
    near(bounds(&native, &f.s, side.part).min[2], 0.);
    let fixed = record(&cuts, "Tall cabinet shelf 1");
    pair(fixed.cut.finished_mm, [764., 420.]);
    near(fixed.cut.thickness_mm, 18.);
    let adjustable = record(&cuts, "Tall cabinet shelf 4");
    pair(adjustable.cut.finished_mm, [762., 350.]);
    near(adjustable.cut.thickness_mm, 6.);
    assert_eq!(adjustable.stock, f.hdf);
    near(bounds(&native, &f.s, adjustable.part).min[2], 1053.5);

    // Remove the first row and move the cabinet; the retained fourth shelf
    // still owns its original Part, stock, depth and local floor dimension.
    let mut edited = definition;
    edited.cabinet.as_mut().unwrap().shelves.remove(0);
    edited.origin = [length(200.), length(100.), length(50.)];
    edited.rotation = angle(90.);
    edited.back = Back::None;
    let p = model::prepare_carcass_edit(&f.s, f.s.revision(), actor(), carcass, &edited).unwrap();
    checked(&mut f.s, p.model.prepared);
    let result = compile(&f.s);
    let changed = result.cut_list(f.s.document()).unwrap();
    let retained = record(&changed, "Tall cabinet shelf 4");
    assert_eq!(retained.part, adjustable.part);
    assert_eq!(retained.stock, f.hdf);
    near(bounds(&result, &f.s, retained.part).min[2], 1103.5);
    near(bounds(&result, &f.s, bottom.part).min[2], 150.);
    assert_eq!(changed.len(), 5);
    f.s.undo(f.s.revision(), actor()).unwrap();
    same_cut_records(&compile(&f.s).cut_list(f.s.document()).unwrap(), &cuts);
    f.s.redo(f.s.revision(), actor()).unwrap();
    same_cut_records(&compile(&f.s).cut_list(f.s.document()).unwrap(), &changed);

    let project = ProjectFile::new(f.s.document().clone()).unwrap();
    assert!(project.capabilities().contains(model::CABINET_CAPABILITY));
    let mut bytes = vec![];
    project_file::write(&project, &mut bytes, Default::default()).unwrap();
    let mut old_reader = project_file::supported_capabilities();
    old_reader.remove(model::CABINET_CAPABILITY);
    assert!(project_file::read(&mut bytes.as_slice(), Default::default(), &old_reader).is_err());
    let opened = project_file::read(
        &mut bytes.as_slice(),
        Default::default(),
        &project_file::supported_capabilities(),
    )
    .unwrap()
    .project;
    assert_eq!(&opened.document, f.s.document());
    same_cut_records(
        &Compiler::default()
            .compile(&opened.document, 0)
            .unwrap()
            .cut_list(&opened.document)
            .unwrap(),
        &changed,
    );
}

#[test]
fn explicit_cabinet_bad_sizes_collisions_keys_and_support_levels_are_source_neutral() {
    let f = fixture();
    let original = f.s.document().clone();
    let definition = explicit_cabinet(&f);
    let mut invalid = vec![];
    let mut d = definition.clone();
    d.cabinet.as_mut().unwrap().toe_clearance = length(-1.);
    invalid.push((d, "Toe clearance"));
    let mut d = definition.clone();
    d.cabinet.as_mut().unwrap().shelves[0].depth = length(600.);
    invalid.push((d, "requested depth"));
    let mut d = definition.clone();
    d.cabinet.as_mut().unwrap().shelves[1].underside = length(501.);
    d.cabinet.as_mut().unwrap().shelves[1].support = model::ShelfSupport::Fixed;
    invalid.push((d, "overlaps"));
    let mut d = definition.clone();
    d.cabinet.as_mut().unwrap().shelves[1].support = model::ShelfSupport::System32 { level: 10000 };
    invalid.push((d, "does not fit"));
    let mut d = definition.clone();
    d.cabinet.as_mut().unwrap().shelves[1].key = 0;
    invalid.push((d, "distinct stable slots"));
    let mut d = definition;
    d.cabinet.as_mut().unwrap().shelves[0].underside = length(105.);
    invalid.push((d, "does not fit"));
    for (d, message) in invalid {
        let error = model::prepare_carcass(
            &f.s,
            f.s.revision(),
            actor(),
            &CarcassRequest::new("Invalid", d),
        )
        .expect_err("invalid cabinet refused");
        assert!(error.to_string().contains(message), "{error}");
        assert_eq!(f.s.document(), &original);
    }
}

#[test]
fn explicit_cabinet_default_stock_updates_inheriting_shelves_only() {
    let mut f = fixture();
    let definition = explicit_cabinet(&f);
    let created = model::prepare_carcass(
        &f.s,
        f.s.revision(),
        actor(),
        &CarcassRequest::new("Stock cabinet", definition.clone()),
    )
    .unwrap();
    let carcass = created.carcass;
    checked(&mut f.s, created.model.prepared);
    let before = compile(&f.s).cut_list(f.s.document()).unwrap();
    let thicker = stock(f.ply_material, 24., Grain::AlongLength, "Birch ply 24");
    let thicker_id = thicker.id;
    let p = prepare_library(&f.s, f.s.revision(), actor(), &[thicker]).unwrap();
    f.s.commit(p).unwrap();
    let original = f.s.document().clone();
    let mut changed = definition;
    changed.cabinet.as_mut().unwrap().shelf_stock = thicker_id;
    let p = model::prepare_carcass_edit(&f.s, f.s.revision(), actor(), carcass, &changed).unwrap();
    // Native preview is a candidate, not a source mutation.
    Compiler::default()
        .compile(p.model.prepared.document(), f.s.revision() + 1)
        .unwrap();
    assert_eq!(f.s.document(), &original);
    checked(&mut f.s, p.model.prepared);
    let after = compile(&f.s).cut_list(f.s.document()).unwrap();
    for name in [
        "Stock cabinet shelf 1",
        "Stock cabinet shelf 4",
        "Stock cabinet left side",
    ] {
        assert_eq!(record(&before, name).part, record(&after, name).part);
    }
    assert_eq!(record(&after, "Stock cabinet shelf 1").stock, thicker_id);
    near(
        record(&after, "Stock cabinet shelf 1").cut.thickness_mm,
        24.,
    );
    for name in ["Stock cabinet shelf 4", "Stock cabinet left side"] {
        same_cut_record(record(&after, name), record(&before, name));
    }
}

#[test]
fn cabinet_connections_follow_parameters_levels_joinery_and_user_owned_dependencies() {
    let mut f = fixture();
    let mut definition = explicit_cabinet(&f);
    definition.width = length(300.);
    definition.height = length(400.);
    definition.depth = length(220.);
    definition.back = Back::None;
    let depth = Dimension::Expression {
        expression: BoundExpression::bind(f.s.document(), "=$cab_width / 4", "mm", None).unwrap(),
    };
    let c = definition.cabinet.as_mut().unwrap();
    c.toe_clearance = length(20.);
    c.shelves[0].underside = length(200.);
    c.shelves[0].depth = depth.clone();
    c.shelves[1].depth = depth;
    c.shelves[1].support = model::ShelfSupport::System32 { level: 3 };
    let created = model::prepare_carcass(
        &f.s,
        f.s.revision(),
        actor(),
        &CarcassRequest::new("Connected", definition.clone()),
    )
    .unwrap();
    let id = created.carcass;
    checked(&mut f.s, created.model.prepared);
    let joints = |s: &Session| {
        model::panel_joints(s.document())
            .filter(|e| {
                model::read_panel_joint(s.document(), e.id)
                    .unwrap()
                    .generated_by
                    == Some(id)
            })
            .map(|e| e.id)
            .collect::<std::collections::BTreeSet<_>>()
    };
    let operations = |s: &Session| {
        model::machinings(s.document())
            .filter(|e| {
                model::read_machining(s.document(), e.id)
                    .unwrap()
                    .generated_by
                    == Some(id)
            })
            .map(|e| e.id)
            .collect::<Vec<_>>()
    };
    let original_joints = joints(&f.s);
    assert_eq!(original_joints.len(), 6);
    let pins = operations(&f.s);
    assert_eq!(pins.len(), 1);
    let list = compile(&f.s).cut_list(f.s.document()).unwrap();
    let adjustable = record(&list, "Connected shelf 4").part;
    let fixed = record(&list, "Connected shelf 1").part;
    let side = record(&list, "Connected left side").part;
    near(bounds(&compile(&f.s), &f.s, adjustable).min[2], 173.5);
    let holes = |s: &Session| {
        compile(s)
            .machining(s.document())
            .unwrap()
            .into_iter()
            .find(|r| r.operation == pins[0] && r.part == side)
            .unwrap()
            .holes
    };
    assert!(
        holes(&f.s)
            .iter()
            .any(|h| (h.center[1] - 133.).abs() < 1e-6)
    );
    set_parameter(&mut f.s, f.width, "640");
    assert_eq!(joints(&f.s), original_joints);
    assert_eq!(operations(&f.s), pins);
    assert!(
        holes(&f.s)
            .iter()
            .any(|h| (h.center[1] - 143.).abs() < 1e-6)
    );
    let current = compile(&f.s).cut_list(f.s.document()).unwrap();
    pair(
        record(&current, "Connected shelf 4").cut.finished_mm,
        [262., 160.],
    );
    // Moving the bottom datum keeps the selected adjustable level and moves
    // its pin support together, while the fixed shelf retains 200 mm.
    definition.cabinet.as_mut().unwrap().toe_clearance = length(30.);
    let before = f.s.document().clone();
    let p = model::prepare_carcass_edit(&f.s, f.s.revision(), actor(), id, &definition).unwrap();
    Compiler::default()
        .compile(p.model.prepared.document(), f.s.revision() + 1)
        .unwrap();
    assert_eq!(f.s.document(), &before);
    checked(&mut f.s, p.model.prepared);
    let m = compile(&f.s);
    near(bounds(&m, &f.s, adjustable).min[2], 183.5);
    near(bounds(&m, &f.s, fixed).min[2], 200.);
    assert_eq!(joints(&f.s), original_joints);
    assert_eq!(operations(&f.s), pins);
    let snapshot = f.s.document().clone();
    assert!(
        model::prepare_panel_joint_delete(
            &f.s,
            f.s.revision(),
            actor(),
            *original_joints.first().unwrap()
        )
        .is_err()
    );
    assert!(model::prepare_machining_delete(&f.s, f.s.revision(), actor(), pins[0]).is_err());
    assert_eq!(f.s.document(), &snapshot);
    // Switching connection family replaces only owned connection geometry.
    definition.cabinet.as_mut().unwrap().joinery =
        model::CabinetJoinery::Domino(Box::new(model::LooseTenonSpec {
            thickness: length(6.),
            width: length(20.),
            length: length(40.),
            count: 2,
            inset: length(60.),
            offset: None,
            face_depth: length(12.),
            edge_depth: length(30.),
            face_play: length(0.),
            edge_play: length(0.),
            fit: length(0.1),
        }));
    let p = model::prepare_carcass_edit(&f.s, f.s.revision(), actor(), id, &definition).unwrap();
    checked(&mut f.s, p.model.prepared);
    assert_eq!(joints(&f.s).len(), 6);
    assert!(joints(&f.s).is_disjoint(&original_joints));
    let m = compile(&f.s);
    let records = m.machining(f.s.document()).unwrap();
    assert_eq!(
        records
            .iter()
            .filter(|r| r.kind == "loose-tenon mortise")
            .count(),
        12
    );
    assert_eq!(records.iter().filter(|r| r.kind == "shelf pins").count(), 2);
    assert!(
        records
            .iter()
            .filter(|r| r.kind == "loose-tenon mortise")
            .all(|r| r.tenon.is_some())
    );
    // Making the fixed shelf adjustable removes its two fixed joints only.
    definition.cabinet.as_mut().unwrap().shelves[0].support =
        model::ShelfSupport::System32 { level: 0 };
    let p = model::prepare_carcass_edit(&f.s, f.s.revision(), actor(), id, &definition).unwrap();
    checked(&mut f.s, p.model.prepared);
    assert_eq!(joints(&f.s).len(), 4);
    assert_eq!(operations(&f.s), pins);
    assert!(
        !compile(&f.s)
            .machining(f.s.document())
            .unwrap()
            .iter()
            .any(|r| r.part == fixed && r.kind == "loose-tenon mortise")
    );
    let user = model::prepare_machining(
        &f.s,
        f.s.revision(),
        actor(),
        "User cable hole",
        &model::MachiningDefinition {
            generated_by: None,
            kind: model::MachiningKind::Drill {
                holes: vec![[length(20.), length(20.)]],
                diameter: length(3.),
                depth: None,
            },
            face: model::MachiningFace::Top,
            targets: vec![adjustable],
            suppressed: false,
        },
    )
    .unwrap();
    let user_id = user.operation;
    checked(&mut f.s, user.model.prepared);
    let mut remove = definition.clone();
    remove.cabinet.as_mut().unwrap().shelves.remove(1);
    let snapshot = f.s.document().clone();
    let error = model::prepare_carcass_edit(&f.s, f.s.revision(), actor(), id, &remove)
        .expect_err("user machining blocks panel removal");
    assert!(error.to_string().contains("User cable hole"), "{error}");
    assert_eq!(f.s.document(), &snapshot);
    let deleted = model::prepare_machining_delete(&f.s, f.s.revision(), actor(), user_id).unwrap();
    checked(&mut f.s, deleted);
    let p = model::prepare_carcass_edit(&f.s, f.s.revision(), actor(), id, &remove).unwrap();
    checked(&mut f.s, p.model.prepared);
    assert_eq!(joints(&f.s).len(), 4);
    assert_eq!(operations(&f.s), pins);
    assert!(f.s.document().get(adjustable).is_none());
    assert!(f.s.document().get(fixed).is_some());
}

#[test]
fn cabinet_receiving_cut_collisions_are_refused_before_apply_and_after_parameter_edits() {
    let mut f = fixture();
    let mut definition = explicit_cabinet(&f);
    definition.width = length(300.);
    definition.height = length(400.);
    definition.depth = length(220.);
    definition.back = Back::None;
    let c = definition.cabinet.as_mut().unwrap();
    c.toe_clearance = length(20.);
    c.pins.start = length(43.);
    c.shelves[0].underside = length(200.);
    c.shelves[0].depth = length(150.);
    c.shelves[1].depth = length(150.);
    c.shelves[1].support = model::ShelfSupport::System32 { level: 3 };
    let model::CabinetJoinery::Dowels(spec) = &mut c.joinery else {
        unreachable!()
    };
    spec.inset = Dimension::Expression {
        expression: BoundExpression::bind(f.s.document(), "=$cab_width / 10", "mm", None).unwrap(),
    };
    let created = model::prepare_carcass(
        &f.s,
        f.s.revision(),
        actor(),
        &CarcassRequest::new("Collision cabinet", definition.clone()),
    )
    .unwrap();
    let id = created.carcass;
    let upstream = created
        .panels
        .iter()
        .find(|(r, _)| *r == PanelRole::LeftSide)
        .unwrap()
        .1
        .feature;
    checked(&mut f.s, created.model.prepared);
    let original = f.s.document().clone();
    let mut invalid = definition.clone();
    let model::CabinetJoinery::Dowels(spec) = &mut invalid.cabinet.as_mut().unwrap().joinery else {
        unreachable!()
    };
    spec.inset = length(37.);
    let error = model::prepare_carcass_edit(&f.s, f.s.revision(), actor(), id, &invalid)
        .expect_err("dowel and pin collision refused");
    assert!(error.to_string().contains("receiving cuts"), "{error}");
    assert_eq!(f.s.document(), &original);
    invalid.width = length(70.);
    invalid.cabinet.as_mut().unwrap().shelves.clear();
    let error = model::prepare_carcass_edit(&f.s, f.s.revision(), actor(), id, &invalid)
        .expect_err("opposed blind depths refused");
    assert!(error.to_string().contains("overlap"), "{error}");
    let mut parameter = f.s.document().get(f.width).unwrap().clone();
    let EntityData::Parameter { expression, .. } = &mut parameter.data else {
        unreachable!()
    };
    *expression = "370".into();
    let candidate =
        f.s.prepare(Request {
            expected_revision: f.s.revision(),
            actor: actor(),
            label: "Colliding parameter edit".into(),
            edits: vec![Mutation::Replace(parameter)],
        })
        .unwrap();
    let error = Compiler::default()
        .compile(candidate.document(), f.s.revision() + 1)
        .err()
        .expect("parameter-only collision refused");
    assert!(error.to_string().contains("receiving cuts"), "{error}");
    assert_eq!(f.s.document(), &original);
    // A downstream failure still allows upstream face/sketch repair queries.
    Compiler::default()
        .compile_sources(candidate.document(), f.s.revision() + 1, &[upstream])
        .unwrap();
    let mut suppressed = vec![];
    for entity in model::machining_cuts(candidate.document()) {
        let mut entity = entity.clone();
        if let EntityData::Feature { suppressed, .. } = &mut entity.data {
            *suppressed = true;
        }
        suppressed.push(Mutation::Replace(entity));
    }
    let muted = candidate.document().edited(&suppressed).unwrap();
    Compiler::default()
        .compile(&muted, f.s.revision() + 1)
        .unwrap();
    assert!(
        model::check_cabinet_connections(f.s.document(), f.s.parameters(), id, Some(&|| true))
            .is_err()
    );
}

#[test]
fn carcass_generates_named_panels_with_finished_and_core_sizes() {
    let mut f = fixture();
    let created = model::prepare_carcass(
        &f.s,
        f.s.revision(),
        actor(),
        &CarcassRequest::new("Base cabinet", cabinet(&f)),
    )
    .unwrap();
    let carcass = created.carcass;
    assert_eq!(created.panels.len(), 7);
    checked(&mut f.s, created.model.prepared);
    let m = compile(&f.s);
    let list = m.cut_list(f.s.document()).unwrap();
    assert_eq!(list.len(), 7);
    for r in &list {
        assert!(r.issues.is_empty(), "{}: {:?}", r.name, r.issues);
        assert!(r.fits_stock, "{}", r.name);
        assert_eq!(r.carcass, Some(carcass));
    }

    // Sides: full height, grain up, front edge banded 2 mm with 0.5 allowance.
    let side = record(&list, "Base cabinet left side");
    assert_eq!(side.role, Some(PanelRole::LeftSide));
    assert_eq!(side.material, Some(f.ply_material));
    assert_eq!(side.stock, f.ply);
    pair(side.cut.finished_mm, [720., 560.]);
    pair(side.cut.core_mm, [720., 558.5]);
    near(side.cut.thickness_mm, 18.);
    near(side.volume_mm3, 18. * 560. * 720.);
    assert_eq!(side.bands.len(), 1);
    assert_eq!(side.bands[0].0, "u_start");
    // Top and bottom fit between the sides, grain along the width.
    let top = record(&list, "Base cabinet top");
    pair(top.cut.finished_mm, [564., 560.]);
    pair(top.cut.core_mm, [564., 558.5]);
    // The inset back sits 10 mm in from the rear, captured on all four sides.
    let back = record(&list, "Base cabinet back");
    assert_eq!(back.stock, f.hdf);
    pair(back.cut.finished_mm, [684., 564.]);
    near(back.cut.thickness_mm, 6.);
    assert!(back.bands.is_empty());
    let b = bounds(&m, &f.s, back.part);
    near(b.min[1], 544.);
    near(b.max[1], 550.);
    // Two evenly spaced shelves, 1 mm clear of each side, from the setback to the back.
    let shelf = record(&list, "Base cabinet shelf 1");
    pair(shelf.cut.finished_mm, [562., 524.]);
    pair(shelf.cut.core_mm, [562., 522.5]);
    near(bounds(&m, &f.s, shelf.part).min[2], 18. + 216.);
    let shelf2 = record(&list, "Base cabinet shelf 2");
    near(bounds(&m, &f.s, shelf2.part).min[2], 18. + 2. * 216. + 18.);
    // The panels are held by a group named after the carcass.
    let group = model::read_carcass(f.s.document(), carcass)
        .unwrap()
        .group
        .unwrap();
    let members = &f.s.document().hierarchy().members[&group];
    assert!(
        created
            .panels
            .iter()
            .all(|(_, o)| members.contains(&o.body))
    );

    // A parameter edit resizes the panels in place; Undo restores them.
    let before = f.s.document().clone();
    set_parameter(&mut f.s, f.width, "800");
    let wide = compile(&f.s).cut_list(f.s.document()).unwrap();
    pair(
        record(&wide, "Base cabinet top").cut.finished_mm,
        [764., 560.],
    );
    pair(
        record(&wide, "Base cabinet shelf 2").cut.core_mm,
        [762., 522.5],
    );
    pair(
        record(&wide, "Base cabinet back").cut.finished_mm,
        [684., 764.],
    );
    assert_eq!(record(&wide, "Base cabinet top").part, top.part);
    f.s.undo(f.s.revision(), actor()).unwrap();
    assert_eq!(f.s.document(), &before);

    // Explicit shelf offsets add a third shelf; existing shelves keep their Parts.
    let mut edited = cabinet(&f);
    edited.shelves.layout = ShelfLayout::Offsets {
        offsets: vec![length(100.), length(300.), length(500.)],
    };
    let p = model::prepare_carcass_edit(&f.s, f.s.revision(), actor(), carcass, &edited).unwrap();
    checked(&mut f.s, p.model.prepared);
    let m = compile(&f.s);
    let list = m.cut_list(f.s.document()).unwrap();
    assert_eq!(list.len(), 8);
    assert_eq!(record(&list, "Base cabinet shelf 1").part, shelf.part);
    let third = record(&list, "Base cabinet shelf 3");
    near(bounds(&m, &f.s, third.part).min[2], 518.);
    assert!(f.s.document().hierarchy().members[&group].contains(&third.part));
    // One shelf removes the others with their properties; Undo brings them back.
    let mut one = cabinet(&f);
    one.shelves.layout = ShelfLayout::Even { count: 1 };
    let with_three = f.s.document().clone();
    let p = model::prepare_carcass_edit(&f.s, f.s.revision(), actor(), carcass, &one).unwrap();
    checked(&mut f.s, p.model.prepared);
    assert!(f.s.document().get(third.part).is_none());
    assert!(
        part_properties(f.s.document(), third.part)
            .unwrap()
            .is_none()
    );
    let list = compile(&f.s).cut_list(f.s.document()).unwrap();
    assert_eq!(list.len(), 6);
    near(
        bounds(&compile(&f.s), &f.s, shelf.part).min[2],
        18. + (684. - 18.) / 2.,
    );
    f.s.undo(f.s.revision(), actor()).unwrap();
    assert_eq!(f.s.document(), &with_three);

    // An overlay back covers the whole rear; sides, top and shelves get shallower.
    let mut overlay = cabinet(&f);
    overlay.back = Back::Overlay { stock: f.hdf };
    overlay.construction = Construction::TopBottomFullWidth;
    let p = model::prepare_carcass_edit(&f.s, f.s.revision(), actor(), carcass, &overlay).unwrap();
    checked(&mut f.s, p.model.prepared);
    let list = compile(&f.s).cut_list(f.s.document()).unwrap();
    pair(
        record(&list, "Base cabinet back").cut.finished_mm,
        [720., 600.],
    );
    pair(
        record(&list, "Base cabinet left side").cut.finished_mm,
        [684., 554.],
    );
    pair(
        record(&list, "Base cabinet top").cut.finished_mm,
        [600., 554.],
    );
    pair(
        record(&list, "Base cabinet shelf 1").cut.finished_mm,
        [562., 534.],
    );

    // Save and reopen keep the editable carcass and identical cut records.
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
    assert!(opened.capabilities().contains(model::CARCASS_CAPABILITY));
    assert!(opened.capabilities().contains(model::PANEL_CAPABILITY));
    assert_eq!(&opened.document, f.s.document());
    let reopened = Compiler::default()
        .compile(&opened.document, 0)
        .unwrap()
        .cut_list(&opened.document)
        .unwrap();
    same_cut_records(&reopened, &list);

    // Sizes that cannot be built are refused before anything is staged.
    let mut crowded = cabinet(&f);
    crowded.shelves.layout = ShelfLayout::Even { count: 40 };
    let e =
        model::prepare_carcass_edit(&f.s, f.s.revision(), actor(), carcass, &crowded).unwrap_err();
    assert!(e.message.contains("do not fit"), "{}", e.message);
    let mut narrow = cabinet(&f);
    narrow.width = length(30.);
    assert!(model::prepare_carcass_edit(&f.s, f.s.revision(), actor(), carcass, &narrow).is_err());
    let mut overlapping = cabinet(&f);
    overlapping.shelves.layout = ShelfLayout::Offsets {
        offsets: vec![length(100.), length(110.)],
    };
    let e = model::prepare_carcass_edit(&f.s, f.s.revision(), actor(), carcass, &overlapping)
        .unwrap_err();
    assert!(e.message.contains("overlaps"), "{}", e.message);
}

#[test]
fn rotated_carcass_places_panels_about_its_origin() {
    let mut f = fixture();
    let mut def = cabinet(&f);
    def.origin = [length(1000.), length(500.), length(100.)];
    def.rotation = angle(90.);
    def.shelves.layout = ShelfLayout::Even { count: 0 };
    let created = model::prepare_carcass(
        &f.s,
        f.s.revision(),
        actor(),
        &CarcassRequest {
            group: false,
            ..CarcassRequest::new("Tall", def)
        },
    )
    .unwrap();
    assert_eq!(created.panels.len(), 5);
    checked(&mut f.s, created.model.prepared);
    let m = compile(&f.s);
    let left = created.panels[0].1.body;
    // Width now runs along +Y, depth along −X.
    let b = bounds(&m, &f.s, left);
    near(b.min[0], 1000. - 560.);
    near(b.max[0], 1000.);
    near(b.min[1], 500.);
    near(b.max[1], 518.);
    near(b.min[2], 100.);
    near(b.max[2], 820.);
    let list = m.cut_list(f.s.document()).unwrap();
    pair(
        record(&list, "Tall left side").cut.finished_mm,
        [720., 560.],
    );
}

#[test]
fn rectangle_panels_band_edges_follow_grain_and_stack_face_to_face() {
    let mut f = fixture();
    let rect = |stock: EntityId, z: f64, l: f64, w: f64, grain: GrainAxis| PanelDefinition {
        output: BodyOutput {
            body: EntityId::new(),
        },
        source: PanelSource::Rectangle {
            stock,
            plane: PanelPlane::Xy,
            origin: [length(0.), length(0.), length(z)],
            length: length(l),
            width: length(w),
        },
        grain,
        edges: PanelEdges::default(),
        oversize: None,
    };
    let mut top = rect(f.ply, 0., 800., 400., GrainAxis::U);
    top.edges = PanelEdges {
        u_start: Some(band(1., 0.)),
        u_end: Some(band(1., 0.)),
        v_start: Some(band(2., 1.)),
        v_end: None,
    };
    top.oversize = Some(length(2.));
    let a = PanelRequest::new("Desk top", top);
    let p = model::prepare_panel(&f.s, f.s.revision(), actor(), &a).unwrap();
    checked(&mut f.s, p.prepared);
    let b = PanelRequest::new("Desk lining", rect(f.hdf, 18., 800., 400., GrainAxis::Free));
    let p = model::prepare_panel(&f.s, f.s.revision(), actor(), &b).unwrap();
    checked(&mut f.s, p.prepared);
    let list = compile(&f.s).cut_list(f.s.document()).unwrap();
    let desk = record(&list, "Desk top");
    pair(desk.cut.finished_mm, [800., 400.]);
    // u loses two 1 mm bands; v loses 2 mm and keeps 1 mm allowance; both get 2 mm oversize.
    pair(desk.cut.core_mm, [800., 401.]);
    assert_eq!(desk.bands.len(), 3);
    assert_eq!(desk.carcass, None);

    // Grain decides whether a long panel fits the 2500 × 1250 stock with grain along its length.
    let along = PanelRequest::new("Rail", rect(f.ply, 100., 1300., 300., GrainAxis::U));
    let p = model::prepare_panel(&f.s, f.s.revision(), actor(), &along).unwrap();
    checked(&mut f.s, p.prepared);
    let mut across = rect(f.ply, 100., 1300., 300., GrainAxis::V);
    across.output.body = along.body;
    let p =
        model::prepare_panel_edit(&f.s, f.s.revision(), actor(), along.feature, &across).unwrap();
    let list = Compiler::default()
        .compile(p.prepared.document(), 0)
        .unwrap()
        .cut_list(p.prepared.document())
        .unwrap();
    assert!(!record(&list, "Rail").fits_stock);
    pair(record(&list, "Rail").cut.core_mm, [300., 1300.]);
    assert!(record(&compile(&f.s).cut_list(f.s.document()).unwrap(), "Rail").fits_stock);

    // A face-to-face stack keeps each layer's own Part, material and stock.
    let (stack, p) = model::prepare_stack(
        &f.s,
        f.s.revision(),
        actor(),
        "Desk laminate",
        &[a.body, b.body],
        Some("PVA".into()),
        true,
    )
    .unwrap();
    checked(&mut f.s, p);
    let m = compile(&f.s);
    let stacks = m.stacks(f.s.document()).unwrap();
    assert_eq!(stacks.len(), 1);
    let s0 = &stacks[0];
    assert!(s0.face_to_face, "{:?}", s0.issues);
    near(s0.total_thickness_mm, 24.);
    assert_eq!(
        s0.layers.iter().map(|l| l.stock).collect::<Vec<_>>(),
        [f.ply, f.hdf]
    );
    assert_ne!(s0.layers[0].material, s0.layers[1].material);
    assert_eq!(s0.adhesive.as_deref(), Some("PVA"));
    let list = m.cut_list(f.s.document()).unwrap();
    assert_eq!(record(&list, "Desk top").stack, Some(stack));
    assert_eq!(record(&list, "Desk lining").stack, Some(stack));
    let group = model::read_stack(f.s.document(), stack)
        .unwrap()
        .group
        .unwrap();
    assert_eq!(
        f.s.document().hierarchy().members[&group],
        vec![a.body, b.body]
    );

    // Moving a layer away is reported; the stack is not silently regrouped.
    let mut lifted = rect(f.hdf, 20., 800., 400., GrainAxis::Free);
    lifted.output.body = b.body;
    let p = model::prepare_panel_edit(&f.s, f.s.revision(), actor(), b.feature, &lifted).unwrap();
    checked(&mut f.s, p.prepared);
    let stacks = compile(&f.s).stacks(f.s.document()).unwrap();
    assert!(!stacks[0].face_to_face);
    assert!(
        stacks[0].issues[0].contains("2.000 mm apart"),
        "{:?}",
        stacks[0].issues
    );

    // Only panel Parts stack, and a layer belongs to one stack.
    let e = model::prepare_stack(
        &f.s,
        f.s.revision(),
        actor(),
        "Again",
        &[a.body, along.body],
        None,
        false,
    )
    .unwrap_err();
    assert!(e.message.contains("already a layer"), "{}", e.message);
    let e = model::prepare_stack(
        &f.s,
        f.s.revision(),
        actor(),
        "Odd",
        &[along.body, f.ply],
        None,
        false,
    )
    .unwrap_err();
    assert!(e.message.contains("not a panel"), "{}", e.message);
    // Carcass panels cannot be edited directly.
    let carcass_panel = PanelDefinition {
        output: BodyOutput { body: a.body },
        source: PanelSource::Carcass {
            carcass: EntityId::new(),
            role: PanelRole::Top,
        },
        grain: GrainAxis::U,
        edges: PanelEdges::default(),
        oversize: None,
    };
    assert!(
        model::prepare_panel_edit(&f.s, f.s.revision(), actor(), a.feature, &carcass_panel)
            .is_err()
    );
}

#[test]
fn sketch_region_panel_takes_stock_thickness_and_follows_stock_edits() {
    let mut f = fixture();
    let (sketch, edits) = Source::create(Datum::Xz, "Side profile").unwrap();
    insert(&mut f.s, edits);
    let mut source = Source::read(f.s.document(), sketch).unwrap().0;
    source.add_rectangle([0., 0.], [450., 700.]).unwrap();
    let regions = sketch_regions::extract(&source);
    assert_eq!(regions.regions.len(), 1);
    let region = regions.regions[0].intent();
    let edits = source.replace(f.s.document(), sketch).unwrap();
    insert(&mut f.s, edits);
    let request = PanelRequest::new(
        "Shaped side",
        PanelDefinition {
            output: BodyOutput {
                body: EntityId::new(),
            },
            source: PanelSource::Profile {
                stock: f.ply,
                profile: ProfileSource { sketch, region },
            },
            grain: GrainAxis::V,
            edges: PanelEdges::default(),
            oversize: None,
        },
    );
    let p = model::prepare_panel(&f.s, f.s.revision(), actor(), &request).unwrap();
    checked(&mut f.s, p.prepared);
    let list = compile(&f.s).cut_list(f.s.document()).unwrap();
    let side = record(&list, "Shaped side");
    pair(side.cut.finished_mm, [700., 450.]);
    near(side.cut.thickness_mm, 18.);
    near(side.volume_mm3, 450. * 700. * 18.);

    // Banding is refused on sketch-region panels.
    let mut banded = model::read_panel(f.s.document(), request.feature).unwrap();
    banded.edges.u_start = Some(band(2., 0.));
    assert!(
        model::prepare_panel_edit(&f.s, f.s.revision(), actor(), request.feature, &banded).is_err()
    );

    // Editing the stock thickness rebuilds the panel; the assignment shows review.
    let stock = f.s.document().get(f.ply).unwrap().clone();
    let LibraryDefinition::Stock(mut spec) = read(&stock).unwrap().unwrap() else {
        unreachable!()
    };
    let StockKind::Panel { thickness_mm, .. } = &mut spec.specification else {
        unreachable!()
    };
    *thickness_mm = 19.;
    let replaced = LibraryDefinition::Stock(spec)
        .entity(f.ply, stock.name.clone())
        .unwrap();
    let p = prepare_library(&f.s, f.s.revision(), actor(), &[replaced]).unwrap();
    checked(&mut f.s, p);
    let list = compile(&f.s).cut_list(f.s.document()).unwrap();
    near(record(&list, "Shaped side").cut.thickness_mm, 19.);
    near(record(&list, "Shaped side").volume_mm3, 450. * 700. * 19.);
}

type Point = [f64; 3];
fn dot3(a: Point, b: Point) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
/// The panel feature at the root of a stock-source chain of joint and
/// machining cuts.
fn source_panel(document: &Document, mut id: EntityId) -> EntityId {
    loop {
        if let Ok(cut) = model::read_panel_joint_cut(document, id) {
            id = cut.input.producer;
        } else if let Ok(cut) = model::read_machining_cut(document, id) {
            id = cut.input.producer;
        } else {
            model::read_panel(document, id).unwrap();
            return id;
        }
    }
}
/// The finished panel box resolved from the saved panel definition.
fn panel_blank(
    document: &Document,
    parameters: &multicad_core::parameters::EvaluatedParameters,
    panel: EntityId,
) -> model::PanelBox {
    let definition = model::read_panel(document, panel).unwrap();
    match model::resolve_panel(document, parameters, &definition)
        .unwrap()
        .shape
    {
        model::PanelShape::Box(blank) => blank,
        model::PanelShape::Profile { .. } => panic!("Expected a rectangular stock panel"),
    }
}
/// The saved sketch region of a profile panel, resolved by core.
fn panel_region(
    document: &Document,
    profile: &ProfileSource,
) -> multicad_core::sketch_regions::Region {
    let source = Source::read(document, profile.sketch).unwrap().0;
    sketch_regions::extract(&source)
        .resolve(&profile.region)
        .expect("profile panel region resolves")
        .clone()
}
/// Extents of a region's outer boundary in its sketch coordinates.
fn region_extents(region: &multicad_core::sketch_regions::Region) -> [[f64; 2]; 2] {
    use multicad_core::sketch_regions::ExactCurve;
    let mut lo = [f64::INFINITY; 2];
    let mut hi = [f64::NEG_INFINITY; 2];
    let mut add = |p: [f64; 2]| {
        for i in 0..2 {
            lo[i] = lo[i].min(p[i]);
            hi[i] = hi[i].max(p[i]);
        }
    };
    for edge in &region.outer.edges {
        match &edge.curve {
            ExactCurve::Line { a, b } => {
                add(*a);
                add(*b);
            }
            ExactCurve::Circle { center, radius_mm } => {
                add([center[0] - radius_mm, center[1] - radius_mm]);
                add([center[0] + radius_mm, center[1] + radius_mm]);
            }
            other => panic!("Unexpected bundled panel boundary {other:?}"),
        }
    }
    [lo, hi]
}

fn blank_corners(blank: &model::PanelBox) -> Vec<Point> {
    let n = blank.normal();
    (0..8)
        .map(|k| {
            let [a, b, c] = [
                f64::from(k & 1) * blank.length_mm,
                f64::from((k >> 1) & 1) * blank.width_mm,
                f64::from((k >> 2) & 1) * blank.thickness_mm,
            ];
            std::array::from_fn(|i| blank.origin[i] + a * blank.u[i] + b * blank.v[i] + c * n[i])
        })
        .collect()
}
fn polygon_area(polygon: &[[f64; 2]]) -> f64 {
    (0..polygon.len())
        .map(|i| {
            let [a, b] = [polygon[i], polygon[(i + 1) % polygon.len()]];
            a[0] * b[1] - b[0] * a[1]
        })
        .sum::<f64>()
        .abs()
        / 2.
}
/// Clip a polygon to an axis-aligned rectangle [x0, x1] x [y0, y1].
fn clip_rect(polygon: &[[f64; 2]], [x0, x1, y0, y1]: [f64; 4]) -> Vec<[f64; 2]> {
    let mut out = polygon.to_vec();
    for (axis, bound, keep_above) in [(0, x0, true), (0, x1, false), (1, y0, true), (1, y1, false)]
    {
        let inside = |p: [f64; 2]| {
            if keep_above {
                p[axis] >= bound
            } else {
                p[axis] <= bound
            }
        };
        let input = std::mem::take(&mut out);
        for i in 0..input.len() {
            let [p, q] = [input[i], input[(i + 1) % input.len()]];
            if inside(p) {
                out.push(p);
            }
            if inside(p) != inside(q) {
                let t = (bound - p[axis]) / (q[axis] - p[axis]);
                out.push([p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])]);
            }
        }
    }
    out
}
/// Material a resolved cut prism removes from a rectangular blank whose
/// edges are parallel to the prism's U, V and sweep.
fn prism_removal(prism: &model::CutPrism, blank: &model::PanelBox) -> f64 {
    let depth = dot3(prism.sweep, prism.sweep).sqrt();
    let axis = prism.sweep.map(|x| x / depth);
    let local: Vec<[f64; 3]> = blank_corners(blank)
        .into_iter()
        .map(|p| {
            let d = std::array::from_fn(|i| p[i] - prism.origin[i]);
            [dot3(d, prism.u), dot3(d, prism.v), dot3(d, axis)]
        })
        .collect();
    let range = |k: usize| {
        local
            .iter()
            .fold([f64::INFINITY, f64::NEG_INFINITY], |[lo, hi], p| {
                [lo.min(p[k]), hi.max(p[k])]
            })
    };
    let [a, b, s] = [range(0), range(1), range(2)];
    for edge in [blank.u, blank.v, blank.normal()] {
        assert!(
            [prism.u, prism.v, axis]
                .iter()
                .any(|d| (dot3(*d, edge).abs() - 1.).abs() < 1e-9),
            "Cut prism is not aligned with its panel"
        );
    }
    let overlap = (depth.min(s[1]) - s[0].max(0.)).max(0.);
    let area = if let Some((center, radius)) = prism.circle {
        assert!(center[0] - radius >= a[0] - 1e-9 && center[0] + radius <= a[1] + 1e-9);
        assert!(center[1] - radius >= b[0] - 1e-9 && center[1] + radius <= b[1] + 1e-9);
        std::f64::consts::PI * radius * radius
    } else if let Some((p, q, radius)) = prism.slot {
        2. * radius * (q[0] - p[0]).hypot(q[1] - p[1]) + std::f64::consts::PI * radius * radius
    } else {
        polygon_area(&clip_rect(&prism.polygon, [a[0], a[1], b[0], b[1]]))
    };
    area * overlap
}

/// Material removed from a rectangular blank by a set of resolved prisms.
/// Rectangular prisms may overlap one another (a groove running through a
/// finger), so their clipped boxes are united exactly by coordinate
/// compression; round holes and slots (clipped by [`prism_removal`]) must
/// stay clear of every box.
fn removed_material(prisms: &[model::CutPrism], blank: &model::PanelBox) -> f64 {
    let n = blank.normal();
    let local = |p: Point| {
        let d: Point = std::array::from_fn(|i| p[i] - blank.origin[i]);
        [dot3(d, blank.u), dot3(d, blank.v), dot3(d, n)]
    };
    let limits = [blank.length_mm, blank.width_mm, blank.thickness_mm];
    let mut boxes: Vec<[[f64; 2]; 3]> = Vec::new();
    let mut round = 0.;
    let mut round_boxes = Vec::new();
    for prism in prisms {
        let depth = dot3(prism.sweep, prism.sweep).sqrt();
        let axis = prism.sweep.map(|x| x / depth);
        let outline: Vec<[f64; 2]> = if let Some((c, r)) = prism.circle {
            vec![[c[0] - r, c[1] - r], [c[0] + r, c[1] + r]]
        } else if let Some((a, b, r)) = prism.slot {
            vec![
                [a[0].min(b[0]) - r, a[1].min(b[1]) - r],
                [a[0].max(b[0]) + r, a[1].max(b[1]) + r],
            ]
        } else {
            prism.polygon.clone()
        };
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for [x, y] in &outline {
            for t in [0., depth] {
                let p = local(std::array::from_fn(|i| {
                    prism.origin[i] + x * prism.u[i] + y * prism.v[i] + t * axis[i]
                }));
                for k in 0..3 {
                    lo[k] = lo[k].min(p[k]);
                    hi[k] = hi[k].max(p[k]);
                }
            }
        }
        let clipped: [[f64; 2]; 3] = std::array::from_fn(|k| [lo[k].max(0.), hi[k].min(limits[k])]);
        if prism.circle.is_some() || prism.slot.is_some() {
            round += prism_removal(prism, blank);
            round_boxes.push(clipped);
        } else {
            let span = |k: usize| (hi[k] - lo[k]).max(0.);
            let section = polygon_area(&prism.polygon) * depth;
            assert!(
                (section - span(0) * span(1) * span(2)).abs() < 1e-6,
                "Joint and recess prisms here are aligned rectangular boxes"
            );
            if (0..3).all(|k| clipped[k][1] > clipped[k][0]) {
                boxes.push(clipped);
            }
        }
    }
    for r in &round_boxes {
        for b in &boxes {
            assert!(
                (0..3).any(|k| r[k][1] <= b[k][0] + 1e-9 || b[k][1] <= r[k][0] + 1e-9),
                "Round holes do not cross rectangular cuts"
            );
        }
    }
    let mut cuts: [Vec<f64>; 3] = std::array::from_fn(|k| {
        let mut values: Vec<f64> = boxes.iter().flat_map(|b| b[k]).collect();
        values.sort_by(f64::total_cmp);
        values.dedup();
        values
    });
    let [xs, ys, zs] = std::mem::take(&mut cuts);
    let mut union = 0.;
    for x in xs.windows(2) {
        for y in ys.windows(2) {
            for z in zs.windows(2) {
                let mid = [(x[0] + x[1]) / 2., (y[0] + y[1]) / 2., (z[0] + z[1]) / 2.];
                if boxes
                    .iter()
                    .any(|b| (0..3).all(|k| b[k][0] < mid[k] && mid[k] < b[k][1]))
                {
                    union += (x[1] - x[0]) * (y[1] - y[0]) * (z[1] - z[0]);
                }
            }
        }
    }
    union + round
}

/// Every resolved joint and machining prism in the cut chain ending at `id`.
fn chain_prisms(
    document: &Document,
    parameters: &multicad_core::parameters::EvaluatedParameters,
    mut id: EntityId,
) -> Vec<model::CutPrism> {
    let mut prisms = Vec::new();
    loop {
        if let Ok(cut) = model::read_panel_joint_cut(document, id) {
            let joint = model::read_panel_joint(document, cut.joint).unwrap();
            let corner =
                model::resolve_panel_joint_pair(document, parameters, &joint, cut.pair as usize)
                    .unwrap();
            prisms.extend(match cut.side {
                model::PanelJointSide::A => corner.cut_a,
                model::PanelJointSide::B => corner.cut_b,
            });
            id = cut.input.producer;
        } else if let Ok(cut) = model::read_machining_cut(document, id) {
            let definition = model::read_machining(document, cut.operation).unwrap();
            prisms.extend(
                model::resolve_machining(document, parameters, &definition, cut.output.body)
                    .unwrap()
                    .prisms,
            );
            id = cut.input.producer;
        } else {
            return prisms;
        }
    }
}
/// Independent expectation of a bundled stock source from core resolution
/// only: the saved panel blank (or sketch region times thickness) less the
/// material of every resolved joint and machining prism in its chain, and
/// the blank's box. Returns (volume, min, max).
fn expected_stock_source(
    document: &Document,
    parameters: &multicad_core::parameters::EvaluatedParameters,
    model: &Model,
    id: EntityId,
) -> (f64, Point, Point) {
    let panel = source_panel(document, id);
    let definition = model::read_panel(document, panel).unwrap();
    let prisms = chain_prisms(document, parameters, id);
    let resolved = model::resolve_panel(document, parameters, &definition).unwrap();
    let (volume, corners) = match resolved.shape {
        model::PanelShape::Box(blank) => (
            blank.volume_mm3() - removed_material(&prisms, &blank),
            blank_corners(&blank),
        ),
        model::PanelShape::Profile {
            profile,
            thickness_mm,
        } => {
            assert!(prisms.is_empty(), "Profile panels here carry no machining");
            let region = panel_region(document, &profile);
            let [lo, hi] = region_extents(&region);
            let frame = model.snapshot().panel_frames[&panel];
            let n = [
                frame.u[1] * frame.v[2] - frame.u[2] * frame.v[1],
                frame.u[2] * frame.v[0] - frame.u[0] * frame.v[2],
                frame.u[0] * frame.v[1] - frame.u[1] * frame.v[0],
            ];
            let corners = (0..8)
                .map(|k| {
                    let x = if k & 1 == 0 { lo[0] } else { hi[0] };
                    let y = if k & 2 == 0 { lo[1] } else { hi[1] };
                    let t = f64::from((k >> 2) & 1) * thickness_mm;
                    std::array::from_fn(|i| {
                        frame.origin[i] + x * frame.u[i] + y * frame.v[i] + t * n[i]
                    })
                })
                .collect();
            (region.area_mm2 * thickness_mm, corners)
        }
    };
    let min = std::array::from_fn(|i| corners.iter().map(|c| c[i]).fold(f64::INFINITY, f64::min));
    let max = std::array::from_fn(|i| {
        corners
            .iter()
            .map(|c| c[i])
            .fold(f64::NEG_INFINITY, f64::max)
    });
    (volume, min, max)
}
/// Assert a Rust stock source against its independent expectation at the
/// former cross-backend budget; returns its solid count.
fn check_stock_source(
    document: &Document,
    parameters: &multicad_core::parameters::EvaluatedParameters,
    name: &str,
    id: EntityId,
    rust: &Model,
) -> u32 {
    assert!(rust.solid(id).unwrap().uses_rust_geometry(), "{name}");
    let shape = &rust.snapshot().shapes[&id];
    assert_eq!(shape.solid_count, 1, "{name}");
    let (volume, min, max) = expected_stock_source(document, parameters, rust, id);
    assert!(
        (shape.volume_mm3 - volume).abs() <= 1e-7 + volume * 1e-10,
        "{name}: {} != {volume}",
        shape.volume_mm3
    );
    for (a, b) in shape
        .bounds
        .min
        .into_iter()
        .chain(shape.bounds.max)
        .zip(min.into_iter().chain(max))
    {
        near(a, b);
    }
    shape.solid_count
}
/// Loose tenon inserts: one solid per resolved tenon, bounded by each
/// stadium (width x thickness) between its face and edge insertion depths.
fn check_tenon_source(
    document: &Document,
    parameters: &multicad_core::parameters::EvaluatedParameters,
    name: &str,
    id: EntityId,
    rust: &Model,
) -> u32 {
    assert!(rust.solid(id).unwrap().uses_rust_geometry(), "{name}");
    let shape = &rust.snapshot().shapes[&id];
    let layout = model::resolve_loose_tenon_feature(document, parameters, id).unwrap();
    let mut min = [f64::INFINITY; 3];
    let mut max = [f64::NEG_INFINITY; 3];
    for tenon in &layout {
        let side = [
            tenon.axis[1] * tenon.long[2] - tenon.axis[2] * tenon.long[1],
            tenon.axis[2] * tenon.long[0] - tenon.axis[0] * tenon.long[2],
            tenon.axis[0] * tenon.long[1] - tenon.axis[1] * tenon.long[0],
        ];
        for i in 0..3 {
            let ends = [
                tenon.center[i] - tenon.into_face_mm * tenon.axis[i],
                tenon.center[i] + tenon.into_edge_mm * tenon.axis[i],
            ];
            let spread =
                tenon.width_mm / 2. * tenon.long[i].abs() + tenon.thickness_mm / 2. * side[i].abs();
            min[i] = min[i].min(ends[0].min(ends[1]) - spread);
            max[i] = max[i].max(ends[0].max(ends[1]) + spread);
        }
    }
    for (a, b) in shape
        .bounds
        .min
        .into_iter()
        .chain(shape.bounds.max)
        .zip(min.into_iter().chain(max))
    {
        near(a, b);
    }
    shape.solid_count
}
/// The cut list records from core: the resolved cut of each panel's blank
/// (or profile extents), with its stock, grain and Part.
fn check_cut_list(
    document: &Document,
    parameters: &multicad_core::parameters::EvaluatedParameters,
    rust: &Model,
) -> usize {
    let cuts = rust.cut_list(document).unwrap();
    for record in &cuts {
        let definition = model::read_panel(document, record.feature).unwrap();
        let resolved = model::resolve_panel(document, parameters, &definition).unwrap();
        let [u, v, t] = match &resolved.shape {
            model::PanelShape::Box(blank) => [blank.length_mm, blank.width_mm, blank.thickness_mm],
            model::PanelShape::Profile {
                profile,
                thickness_mm,
            } => {
                let [lo, hi] = region_extents(&panel_region(document, profile));
                [hi[0] - lo[0], hi[1] - lo[1], *thickness_mm]
            }
        };
        let expected = resolved.cut(u, v, t).unwrap();
        assert_eq!(
            (record.part, record.stock, record.grain),
            (definition.output.body, resolved.stock, resolved.grain),
            "{}",
            record.name
        );
        pair(record.cut.finished_mm, expected.finished_mm);
        near(record.cut.thickness_mm, expected.thickness_mm);
    }
    cuts.len()
}

type Expected = (f64, Point, Point);
fn hull(points: impl IntoIterator<Item = Point>) -> (Point, Point) {
    let mut min = [f64::INFINITY; 3];
    let mut max = [f64::NEG_INFINITY; 3];
    for p in points {
        for i in 0..3 {
            min[i] = min[i].min(p[i]);
            max[i] = max[i].max(p[i]);
        }
    }
    (min, max)
}
fn resolved_member(
    document: &Document,
    parameters: &multicad_core::parameters::EvaluatedParameters,
    feature: EntityId,
) -> model::ResolvedMember {
    let definition = model::read_member(document, feature).unwrap();
    let path = model::point_path(&definition.path, parameters)
        .unwrap()
        .expect("bundled members use point paths");
    model::resolve_member(document, parameters, &definition, path).unwrap()
}
/// A straight weldment member from its authored stock, path and the plane
/// end cuts of its resolved connections. A prism cut by planes has the
/// volume of its section area times the axial length through the section
/// centroid; its box is the hull of the section corners on both end planes.
fn expected_member(
    document: &Document,
    parameters: &multicad_core::parameters::EvaluatedParameters,
    part: EntityId,
    feature: EntityId,
) -> Expected {
    let member = resolved_member(document, parameters, feature);
    let section = &member.section;
    assert_eq!(
        section.shape,
        multicad_core::materials::ProfileShape::Rectangular,
        "{}",
        section.designation
    );
    let mut cuts: [Option<([f64; 3], [f64; 3])>; 2] = [None, None];
    for (_, connection) in model::member_connections(document, part).unwrap() {
        let side = usize::from(connection.members[1] == part);
        let members = connection.members.map(|m| {
            resolved_member(
                document,
                parameters,
                model::member_feature(document, m).unwrap(),
            )
        });
        let resolved = model::resolve_connection(
            connection.kind,
            &members[0],
            &members[1],
            connection.gap_mm(parameters).unwrap(),
        )
        .unwrap();
        let Some(treatment) = resolved.treatments[side] else {
            continue;
        };
        let model::EndCut::Plane { point, normal } = treatment.cut else {
            panic!("Bundled workbench joints are plane cuts");
        };
        let end = usize::from(treatment.end == model::MemberEnd::End);
        assert!(cuts[end].is_none(), "One cut per member end");
        cuts[end] = Some((point, normal));
    }
    let axis = member.frame.axis;
    let station = |p: [f64; 2], end: usize| {
        let base = member.point(p, 0.);
        match cuts[end] {
            Some((point, normal)) => {
                let d = std::array::from_fn(|i| point[i] - base[i]);
                dot3(d, normal) / dot3(axis, normal)
            }
            None => member.start_offset_mm() + if end == 1 { member.length_mm() } else { 0. },
        }
    };
    let (w, h) = (section.width_mm, section.height_mm);
    let centroid = [w / 2., h / 2.];
    let volume = section.area_mm2() * (station(centroid, 1) - station(centroid, 0));
    let corners = [[0., 0.], [w, 0.], [w, h], [0., h]];
    let (min, max) = hull(
        corners
            .into_iter()
            .flat_map(|c| [0, 1].map(|end| member.point(c, station(c, end)))),
    );
    (volume, min, max)
}
fn cylinder_bounds(
    radius: f64,
    height: f64,
    frame: &multicad_core::sketch::Frame,
) -> (Point, Point) {
    let mut min = [0.; 3];
    let mut max = [0.; 3];
    for i in 0..3 {
        let ends = [frame.origin[i], frame.origin[i] + height * frame.n[i]];
        let rim = radius * frame.u[i].hypot(frame.v[i]);
        min[i] = ends[0].min(ends[1]) - rim;
        max[i] = ends[0].max(ends[1]) + rim;
    }
    (min, max)
}
fn primitive(
    document: &Document,
    parameters: &multicad_core::parameters::EvaluatedParameters,
    feature: EntityId,
) -> model::ResolvedPrimitive {
    model::read_primitive(document, feature)
        .unwrap()
        .parameters
        .resolve(parameters)
        .unwrap()
}
fn expected_primitive(resolved: &model::ResolvedPrimitive) -> Expected {
    match resolved {
        model::ResolvedPrimitive::Box { size_mm, frame } => {
            let (min, max) = hull((0..8).map(|k| {
                let [a, b, c] = [0, 1, 2].map(|j| f64::from((k >> j) & 1) * size_mm[j]);
                std::array::from_fn(|i| {
                    frame.origin[i] + a * frame.u[i] + b * frame.v[i] + c * frame.n[i]
                })
            }));
            (size_mm.iter().product(), min, max)
        }
        model::ResolvedPrimitive::Cylinder {
            radius_mm,
            height_mm,
            frame,
        } => {
            let (min, max) = cylinder_bounds(*radius_mm, *height_mm, frame);
            (
                std::f64::consts::PI * radius_mm * radius_mm * height_mm,
                min,
                max,
            )
        }
        other => panic!("Unexpected bundled primitive {other:?}"),
    }
}
/// Material shared by a cylinder of radius `a` and a perpendicular cylinder
/// of radius `big` whose axes intersect, with `a` passing fully through:
/// 4 a^2 * integral of cos^2(t) sqrt(big^2 - a^2 sin^2(t)) over -pi/2..pi/2,
/// a smooth periodic integrand summed by the exponentially convergent
/// trapezoid rule.
fn crossing_cylinders(a: f64, big: f64) -> f64 {
    let n = 4096;
    let step = std::f64::consts::PI / f64::from(n);
    (0..n)
        .map(|k| {
            let t = -std::f64::consts::FRAC_PI_2 + (f64::from(k) + 0.5) * step;
            t.cos().powi(2) * (big * big - a * a * t.sin().powi(2)).sqrt()
        })
        .sum::<f64>()
        * step
        * 4.
        * a
        * a
}
/// The bundled vise Boolean bodies, from their authored primitives:
/// a box cut by disjoint cylinders parallel to one box edge, and a coaxial
/// screw-and-hub join cut by a perpendicular cylinder crossing the hub axis.
fn expected_boolean(
    document: &Document,
    parameters: &multicad_core::parameters::EvaluatedParameters,
    feature: EntityId,
) -> Expected {
    use model::ResolvedPrimitive::{Box as Block, Cylinder};
    let cut = model::read_boolean(document, feature).unwrap();
    assert_eq!(cut.operation, model::Operation::Cut);
    let tools: Vec<_> = cut
        .tools
        .iter()
        .map(|tool| primitive(document, parameters, tool.producer))
        .collect();
    if let Ok(join) = model::read_boolean(document, cut.target.producer) {
        assert_eq!(join.operation, model::Operation::Join);
        let parts = std::iter::once(join.target.producer)
            .chain(join.tools.iter().map(|t| t.producer))
            .map(|id| primitive(document, parameters, id))
            .collect::<Vec<_>>();
        let [
            Cylinder {
                radius_mm: r1,
                height_mm: h1,
                frame: f1,
            },
            Cylinder {
                radius_mm: r2,
                height_mm: h2,
                frame: f2,
            },
        ] = parts.as_slice()
        else {
            panic!("The joined screw is two cylinders");
        };
        let axis = f1.n;
        assert!((dot3(axis, f2.n) - 1.).abs() < 1e-12);
        let along = |p: Point| dot3(p, axis);
        let offset: Point = std::array::from_fn(|i| f2.origin[i] - f1.origin[i]);
        assert!(
            dot3(offset, offset) - along(offset).powi(2) < 1e-18,
            "coaxial"
        );
        let spans = [
            [along(f1.origin), along(f1.origin) + h1],
            [along(f2.origin), along(f2.origin) + h2],
        ];
        let overlap = (spans[0][1].min(spans[1][1]) - spans[0][0].max(spans[1][0])).max(0.);
        let pi = std::f64::consts::PI;
        let union = pi * r1 * r1 * h1 + pi * r2 * r2 * h2 - pi * r1.min(*r2).powi(2) * overlap;
        let [
            Cylinder {
                radius_mm: a,
                height_mm: ha,
                frame: fa,
            },
        ] = tools.as_slice()
        else {
            panic!("One crossing bore");
        };
        assert!(dot3(fa.n, axis).abs() < 1e-12, "perpendicular bore");
        // The bore crosses only the larger hub, which it spans completely.
        let (big, span) = if r2 > r1 {
            (*r2, spans[1])
        } else {
            (*r1, spans[0])
        };
        let at = along(fa.origin);
        assert!(at - a >= span[0] && at + a <= span[1]);
        let from_axis: Point = std::array::from_fn(|i| fa.origin[i] - f1.origin[i]);
        let radial: Point = std::array::from_fn(|i| from_axis[i] - dot3(from_axis, axis) * axis[i]);
        assert!(
            dot3(radial, fa.n).abs() < 1e-12
                || dot3(radial, radial) - dot3(radial, fa.n).powi(2) < 1e-18
        );
        let reach = [dot3(fa.origin, fa.n), dot3(fa.origin, fa.n) + ha];
        let hub_center = dot3(f1.origin, fa.n);
        assert!(reach[0] <= hub_center - big && reach[1] >= hub_center + big);
        let (min1, max1) = cylinder_bounds(*r1, *h1, f1);
        let (min2, max2) = cylinder_bounds(*r2, *h2, f2);
        let min = std::array::from_fn(|i| min1[i].min(min2[i]));
        let max = std::array::from_fn(|i| max1[i].max(max2[i]));
        return (union - crossing_cylinders(*a, big), min, max);
    }
    let block = primitive(document, parameters, cut.target.producer);
    let Block { size_mm, frame } = &block else {
        panic!("The bored housing is a box");
    };
    let (volume, min, max) = expected_primitive(&block);
    let axes = [frame.u, frame.v, frame.n];
    let mut removed = 0.;
    for tool in &tools {
        let Cylinder {
            radius_mm,
            height_mm,
            frame: f,
        } = tool
        else {
            panic!("Housing bores are cylinders");
        };
        let k = (0..3)
            .find(|&k| (dot3(axes[k], f.n).abs() - 1.).abs() < 1e-12)
            .expect("bore parallel to a housing edge");
        let local = |p: Point| {
            let d: Point = std::array::from_fn(|i| p[i] - frame.origin[i]);
            [0, 1, 2].map(|j| dot3(d, axes[j]))
        };
        let c = local(f.origin);
        for j in (0..3).filter(|&j| j != k) {
            assert!(c[j] - radius_mm >= 0. && c[j] + radius_mm <= size_mm[j]);
        }
        let ends = [c[k], c[k] + height_mm * dot3(f.n, axes[k])];
        let inside = ends[0].max(ends[1]).min(size_mm[k]) - ends[0].min(ends[1]).max(0.);
        removed += std::f64::consts::PI * radius_mm * radius_mm * inside.max(0.);
    }
    for (i, a) in tools.iter().enumerate() {
        for b in &tools[i + 1..] {
            let (
                Cylinder {
                    radius_mm: ra,
                    frame: fa,
                    ..
                },
                Cylinder {
                    radius_mm: rb,
                    frame: fb,
                    ..
                },
            ) = (a, b)
            else {
                unreachable!()
            };
            let d: Point = std::array::from_fn(|i| fb.origin[i] - fa.origin[i]);
            let lateral = dot3(d, d) - dot3(d, fa.n).powi(2);
            assert!(lateral.sqrt() >= ra + rb, "Housing bores are disjoint");
        }
    }
    (volume - removed, min, max)
}
/// Independent expectation of every final bundled workbench body.
fn expected_workbench_body(
    document: &Document,
    parameters: &multicad_core::parameters::EvaluatedParameters,
    rust: &Model,
    part: EntityId,
    feature: EntityId,
) -> Expected {
    if model::read_member(document, feature).is_ok() {
        expected_member(document, parameters, part, feature)
    } else if model::read_primitive(document, feature).is_ok() {
        expected_primitive(&primitive(document, parameters, feature))
    } else if model::read_boolean(document, feature).is_ok() {
        expected_boolean(document, parameters, feature)
    } else {
        expected_stock_source(document, parameters, rust, feature)
    }
}

fn read_bundled_project(mut bytes: &[u8]) -> ProjectFile {
    project_file::read(
        &mut bytes,
        Default::default(),
        &project_file::supported_capabilities(),
    )
    .unwrap()
    .project
}

#[test]
fn rust_bundled_workbench_rebuilds_all_parts_with_through_and_blind_drilling() {
    let project = read_bundled_project(include_bytes!("../../../examples/welded-workbench.tcd"));
    let document = &project.document;
    let compile = |backend| {
        Compiler::default()
            .compile_with_backend(document, 0, backend, None)
            .unwrap_or_else(|error| panic!("{backend:?}: {error:?}"))
    };
    let parameters = multicad_core::parameters::evaluate(document).unwrap();
    let rust = compile(GeometryBackend::Rust);
    // 19 welded members, the worktop and the 12 under-bench vise parts that
    // are not consumed as Boolean tools.
    assert_eq!(rust.snapshot().bodies.len(), 32);
    for (&part, &feature) in &rust.snapshot().bodies {
        let name = &document.get(part).unwrap().name;
        assert!(rust.solid(feature).unwrap().uses_rust_geometry(), "{name}");
        let shape = &rust.snapshot().shapes[&feature];
        assert_eq!(shape.solid_count, 1, "{name}");
        let (volume, min, max) =
            expected_workbench_body(document, &parameters, &rust, part, feature);
        assert!(
            (shape.volume_mm3 - volume).abs() <= 1e-7 + volume * 1e-10,
            "{name}: {} != {volume}",
            shape.volume_mm3
        );
        for (a, b) in shape
            .bounds
            .min
            .into_iter()
            .chain(shape.bounds.max)
            .zip(min.into_iter().chain(max))
        {
            near(a, b);
        }
    }

    // Machining records follow the saved operations resolved by core: the
    // worktop grid and tab pilots and the four vise jaw bores.
    let resolved_holes: Vec<_> = model::machining_cuts(document)
        .flat_map(|entity| {
            let cut = model::read_machining_cut(document, entity.id).unwrap();
            let definition = model::read_machining(document, cut.operation).unwrap();
            model::resolve_machining(document, &parameters, &definition, cut.output.body)
                .unwrap()
                .holes
        })
        .collect();
    let records = rust.machining(document).unwrap();
    assert_eq!(records.len(), model::machinings(document).count());
    assert_eq!(records.len(), 6);
    let through = resolved_holes.iter().filter(|hole| hole.through).count();
    assert_eq!(through, 106);
    assert_eq!(
        records
            .iter()
            .flat_map(|record| &record.holes)
            .filter(|hole| hole.through)
            .count(),
        through
    );
    assert_eq!(
        records
            .iter()
            .flat_map(|record| &record.holes)
            .filter(|hole| !hole.through)
            .count(),
        resolved_holes.len() - through
    );
    let worktop = rust.cut_list(document).unwrap().pop().unwrap().part;
    let feature = rust.snapshot().bodies[&worktop];
    let solid = rust.solid(feature).unwrap();
    let geometry = solid.topology_snapshot(feature, 0).unwrap();
    let rims = geometry
        .elements
        .iter()
        .filter(|element| matches!(element.geometry, Some(AnalyticGeometry::CircularArc { .. })));
    assert_eq!(rims.count(), 208);
    let receipt = solid
        .manufacturing_mesh_with_receipt(MeshOptions::default())
        .unwrap();
    assert!(receipt.edge_pairs > 0 && receipt.vertex_links > 0);
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("worktop.step");
    solid.write_step(&path).unwrap();
    let policy = GeometryPolicy::default();
    let readbacks = [
        Solid::read_step(&path, policy).unwrap(),
        Solid::import_rust_step(&std::fs::read_to_string(&path).unwrap(), feature, policy).unwrap(),
    ];
    for imported in readbacks {
        assert_eq!(imported.solid_count().unwrap(), 1);
        near(imported.volume_mm3().unwrap(), solid.volume_mm3().unwrap());
    }
}

fn compile_stock_sources(
    document: &Document,
    roots: &[EntityId],
    backend: GeometryBackend,
) -> Model {
    Compiler::default()
        .compile_sources_with_backend(document, 0, roots, backend, None)
        .unwrap_or_else(|error| {
            let name = error
                .entity
                .and_then(|id| document.get(id))
                .map(|e| &e.name);
            panic!("{backend:?} failed for {name:?}: {error}")
        })
}

/// Two Rust evaluations of the same saved panels report the same cuts.
fn same_stock_cut_lists(document: &Document, rust: &Model, rebuilt: &Model) -> usize {
    let rust_cuts = rust.cut_list(document).unwrap();
    let rebuilt_cuts = rebuilt.cut_list(document).unwrap();
    assert_eq!(rust_cuts.len(), rebuilt_cuts.len());
    for (a, b) in rust_cuts.iter().zip(rebuilt_cuts) {
        assert_eq!((a.part, a.stock, a.grain), (b.part, b.stock, b.grain));
        pair(a.cut.finished_mm, b.cut.finished_mm);
        near(a.cut.thickness_mm, b.cut.thickness_mm);
    }
    rust_cuts.len()
}

#[test]
fn rust_rebuilds_the_bundled_cabinet_panels_finger_joints_and_drawer_grooves() {
    let project = read_bundled_project(include_bytes!("../../../examples/cabinet.tcd"));
    let document = &project.document;
    let panels: Vec<_> = model::panels(document).map(|entity| entity.id).collect();
    let finger_cuts: Vec<_> = model::panel_joint_cuts(document)
        .filter(|entity| {
            let cut = model::read_panel_joint_cut(document, entity.id).unwrap();
            let joint = model::read_panel_joint(document, cut.joint).unwrap();
            matches!(joint.kind, PanelJointKind::Finger { .. })
        })
        .collect();
    let grooves: Vec<_> = model::machining_cuts(document)
        .filter(|entity| {
            let cut = model::read_machining_cut(document, entity.id).unwrap();
            let source = model::read_machining(document, cut.operation).unwrap();
            matches!(source.kind, MachiningKind::Groove { .. })
        })
        .map(|entity| entity.id)
        .collect();
    assert_eq!(
        (panels.len(), finger_cuts.len(), grooves.len()),
        (24, 16, 8)
    );
    let roots: Vec<_> = panels
        .iter()
        .copied()
        .chain(finger_cuts.iter().map(|entity| entity.id))
        .chain(grooves)
        .collect();
    let parameters = multicad_core::parameters::evaluate(document).unwrap();
    let rust = compile_stock_sources(document, &roots, GeometryBackend::Rust);
    for id in roots {
        let name = &document.get(id).unwrap().name;
        check_stock_source(document, &parameters, name, id, &rust);
        if panels.contains(&id) {
            // The panel frame and both cap attachments are the resolved
            // finished panel box: its corner, U and V, and the far face.
            let blank = panel_blank(document, &parameters, id);
            assert_eq!(
                rust.snapshot().panel_frames[&id],
                PlaneFrame {
                    origin: blank.origin,
                    u: blank.u,
                    v: blank.v,
                }
            );
            let a = &rust.snapshot().shapes[&id];
            let normal = blank.normal();
            for (role, depth) in [("cap.start", 0.), ("cap.end", blank.thickness_mm)] {
                let face = a
                    .topology
                    .elements
                    .iter()
                    .find(|e| e.frames.iter().any(|f| f.origin.role == role))
                    .unwrap_or_else(|| panic!("{name} has no {role} attachment"));
                let anchor = face.frames.iter().find(|f| f.origin.role == role).unwrap();
                let mut intent = GeometryIntent::capture(&a.topology, face.key).unwrap();
                intent.attachment = Some(Attachment {
                    anchor: anchor.origin.clone(),
                    offset_mm: [0.; 3],
                    angle_rad: 0.,
                });
                let Resolution::Resolved {
                    frame: Some(frame), ..
                } = references::resolve(&intent, &a.topology, 0)
                else {
                    panic!("{name} lost authored face attachment {:?}", anchor.origin);
                };
                let origin: Point = std::array::from_fn(|i| blank.origin[i] + depth * normal[i]);
                for (a, b) in frame
                    .origin_mm
                    .into_iter()
                    .chain(frame.u)
                    .chain(frame.v)
                    .zip(origin.into_iter().chain(blank.u).chain(blank.v))
                {
                    near(a, b);
                }
            }
        }
    }
    assert!(check_cut_list(document, &parameters, &rust) >= 8);
}

#[test]
fn rust_rebuilds_the_whole_bundled_cabinet_with_exact_tenon_inserts_and_mortises() {
    let project = read_bundled_project(include_bytes!("../../../examples/cabinet.tcd"));
    let document = &project.document;
    let parameters = multicad_core::parameters::evaluate(document).unwrap();
    let mortises: Vec<_> = model::panel_joint_cuts(document)
        .filter_map(|entity| {
            let cut = model::read_panel_joint_cut(document, entity.id).unwrap();
            let joint = model::read_panel_joint(document, cut.joint).unwrap();
            matches!(joint.kind, PanelJointKind::LooseTenon(_)).then_some((entity, cut, joint))
        })
        .collect();
    let inserts: Vec<_> = model::loose_tenons(document).collect();
    assert_eq!((mortises.len(), inserts.len()), (16, 8));
    let roots: Vec<_> = mortises
        .iter()
        .map(|(entity, _, _)| entity.id)
        .chain(inserts.iter().map(|entity| entity.id))
        .collect();
    let compile = |backend| {
        Compiler::default()
            .compile_with_backend(document, 0, backend, None)
            .unwrap_or_else(|error| panic!("{backend:?} failed: {error}"))
    };
    let rust = compile(GeometryBackend::Rust);
    let a = rust.snapshot();
    assert_eq!((a.shapes.len(), a.bodies.len()), (137, 19));
    for (&id, shape) in &a.shapes {
        assert!(rust.solid(id).unwrap().uses_rust_geometry(), "{id}");
        assert!(!shape.mesh.indices.is_empty(), "{id}");
    }
    let mut kinds = [0; 3];
    for (&body, &id) in &a.bodies {
        let name = &document.get(body).unwrap().name;
        if model::resolve_loose_tenon_feature(document, &parameters, id).is_ok() {
            check_tenon_source(document, &parameters, name, id, &rust);
            kinds[1] += 1;
        } else if let Ok(primitive) = model::read_primitive(document, id) {
            // The one shelf support pin is a placed cylinder primitive.
            let model::ResolvedPrimitive::Cylinder {
                radius_mm,
                height_mm,
                frame,
            } = primitive.parameters.resolve(&parameters).unwrap()
            else {
                panic!("{name} is not a cylinder");
            };
            let shape = &a.shapes[&id];
            assert_eq!(shape.solid_count, 1);
            near(
                shape.volume_mm3,
                std::f64::consts::PI * radius_mm * radius_mm * height_mm,
            );
            for i in 0..3 {
                let ends = [frame.origin[i], frame.origin[i] + height_mm * frame.n[i]];
                let rim = radius_mm * frame.u[i].hypot(frame.v[i]);
                near(shape.bounds.min[i], ends[0].min(ends[1]) - rim);
                near(shape.bounds.max[i], ends[0].max(ends[1]) + rim);
            }
            kinds[2] += 1;
        } else {
            check_stock_source(document, &parameters, name, id, &rust);
            kinds[0] += 1;
        }
    }
    assert_eq!(kinds, [10, 8, 1]);
    let stadium_area =
        |span: f64, radius: f64| 2. * radius * span + std::f64::consts::PI * radius.powi(2);
    for (entity, cut, joint) in mortises {
        let id = entity.id;
        check_stock_source(document, &parameters, &entity.name, id, &rust);
        let PanelJointKind::LooseTenon(spec) = joint.kind else {
            unreachable!();
        };
        let (corner, layout) = model::resolve_loose_tenon(
            document,
            &parameters,
            &spec,
            &joint.pairs[cut.pair as usize],
        )
        .unwrap();
        let tools = match cut.side {
            model::PanelJointSide::A => corner.cut_a,
            model::PanelJointSide::B => corner.cut_b,
        };
        let recesses = if cut.side == layout.face_side {
            layout.face_mortises
        } else {
            layout.edge_mortises
        };
        assert_eq!(tools.len(), recesses.len());
        let solid = rust.solid(id).unwrap();
        let faces = solid.faces().unwrap();
        let mut removed = 0.;
        for (tool, recess) in tools.iter().zip(recesses) {
            let (a, b, radius) = tool.slot.unwrap();
            let area = stadium_area((b[0] - a[0]).hypot(b[1] - a[1]), radius);
            removed += area * recess.depth_mm;
            let center: [f64; 3] = std::array::from_fn(|i| {
                tool.origin[i]
                    + tool.sweep[i]
                    + tool.u[i] * (a[0] + b[0]) / 2.
                    + tool.v[i] * (a[1] + b[1]) / 2.
            });
            let floor = faces
                .iter()
                .find(|face| {
                    face.plane.is_some()
                        && face
                            .center_mm
                            .into_iter()
                            .zip(center)
                            .all(|(a, b)| (a - b).abs() < 1e-6)
                })
                .expect("mortise has a planar floor at its resolved center");
            near(floor.area_mm2, area);
        }
        near(
            rust.solid(cut.input.producer)
                .unwrap()
                .volume_mm3()
                .unwrap()
                - solid.volume_mm3().unwrap(),
            removed,
        );
    }
    for entity in inserts {
        let id = entity.id;
        let layout = model::resolve_loose_tenon_feature(document, &parameters, id).unwrap();
        assert_eq!(
            check_tenon_source(document, &parameters, &entity.name, id, &rust) as usize,
            layout.len()
        );
        let volume = layout
            .iter()
            .map(|tenon| {
                let radius = tenon.thickness_mm / 2.;
                stadium_area((tenon.width_mm - 2. * radius).max(0.), radius)
                    * (tenon.into_face_mm + tenon.into_edge_mm)
            })
            .sum();
        near(rust.solid(id).unwrap().volume_mm3().unwrap(), volume);
    }
    for id in roots {
        let name = &document.get(id).unwrap().name;
        multicad_model::export::mesh::tessellate_manufacturing(
            rust.solid(id).unwrap(),
            &Default::default(),
            name,
        )
        .unwrap_or_else(|error| panic!("{name}: {error}"));
    }
}

#[test]
fn rust_rebuilds_both_hinge_boring_groups_on_each_bundled_cabinet_door() {
    let project = read_bundled_project(include_bytes!("../../../examples/cabinet.tcd"));
    let document = &project.document;
    let parameters = multicad_core::parameters::evaluate(document).unwrap();
    let operations: Vec<_> = model::machining_cuts(document)
        .filter_map(|entity| {
            let cut = model::read_machining_cut(document, entity.id).unwrap();
            document
                .get(cut.operation)
                .unwrap()
                .name
                .contains("door hinge boring")
                .then_some((entity, cut))
        })
        .collect();
    assert_eq!(operations.len(), 4);
    let roots: Vec<_> = operations.iter().map(|(entity, _)| entity.id).collect();
    let rust = compile_stock_sources(document, &roots, GeometryBackend::Rust);
    let mut surviving_floors = 0;
    for (entity, cut) in operations {
        let id = entity.id;
        check_stock_source(document, &parameters, &entity.name, id, &rust);
        let definition = model::read_machining(document, cut.operation).unwrap();
        let resolved =
            model::resolve_machining(document, &parameters, &definition, cut.output.body).unwrap();
        assert_eq!(resolved.holes.len(), 3);
        let removed: f64 = resolved
            .holes
            .iter()
            .map(|hole| {
                assert!(!hole.through);
                std::f64::consts::PI * (hole.diameter_mm / 2.).powi(2) * hole.depth_mm
            })
            .sum();
        let solid = rust.solid(id).unwrap();
        near(
            rust.solid(cut.input.producer)
                .unwrap()
                .volume_mm3()
                .unwrap()
                - solid.volume_mm3().unwrap(),
            removed,
        );
        let snapshot = &rust.snapshot().shapes[&id].topology;
        let first_hole = resolved.holes[0];
        let expected: [f64; 3] = std::array::from_fn(|i| {
            first_hole.center[i] + first_hole.axis[i] * first_hole.depth_mm
        });
        let floor = solid
            .faces()
            .unwrap()
            .into_iter()
            .find(|face| {
                face.plane.is_some()
                    && face
                        .center_mm
                        .into_iter()
                        .zip(expected)
                        .all(|(a, b)| (a - b).abs() < 1e-6)
            })
            .expect("hinge cup has a planar floor at its resolved center");
        near(
            floor.area_mm2,
            std::f64::consts::PI * (first_hole.diameter_mm / 2.).powi(2),
        );
        near(
            floor
                .plane
                .unwrap()
                .1
                .into_iter()
                .zip(first_hole.axis)
                .map(|(a, b)| a * b)
                .sum(),
            -1.,
        );
        let entry = snapshot
            .elements
            .iter()
            .find(|item| {
                matches!(item.geometry, Some(AnalyticGeometry::Plane { .. }))
                    && item.key.ordinal == floor.index
            })
            .unwrap();
        let anchor = entry
            .frames
            .iter()
            .find(|frame| frame.origin.source == id)
            .unwrap();
        let mut intent = GeometryIntent::capture(snapshot, entry.key).unwrap();
        intent.attachment = Some(Attachment {
            anchor: anchor.origin.clone(),
            offset_mm: [0.; 3],
            angle_rad: 0.,
        });
        if let Some(next) = roots.iter().copied().find(|next| {
            model::read_machining_cut(document, *next)
                .unwrap()
                .input
                .producer
                == id
        }) {
            let terminal = &rust.snapshot().shapes[&next].topology;
            let Resolution::Resolved {
                frame: Some(frame), ..
            } = references::resolve(&intent, terminal, 0)
            else {
                panic!("the second hinge group lost the first floor attachment");
            };
            assert_eq!(frame, anchor.frame);
            surviving_floors += 1;
        }
    }
    assert_eq!(surviving_floors, 2);
}

#[test]
fn rust_drills_bundled_drawer_runner_holes_after_finger_joints_and_grooves() {
    let project = read_bundled_project(include_bytes!("../../../examples/cabinet.tcd"));
    let document = &project.document;
    let parameters = multicad_core::parameters::evaluate(document).unwrap();
    let operations: Vec<_> = model::machining_cuts(document)
        .filter_map(|entity| {
            let cut = model::read_machining_cut(document, entity.id).unwrap();
            let name = &document.get(cut.operation).unwrap().name;
            (name.contains("runner screws") && name.ends_with("box side")).then_some((entity, cut))
        })
        .collect();
    assert_eq!(operations.len(), 4);
    let roots: Vec<_> = operations.iter().map(|(entity, _)| entity.id).collect();
    let rust = compile_stock_sources(document, &roots, GeometryBackend::Rust);
    let mut saved_floor_intents = Vec::new();
    let dot = |a: [f64; 3], b: [f64; 3]| a.into_iter().zip(b).map(|(x, y)| x * y).sum::<f64>();
    let same_point =
        |a: [f64; 3], b: [f64; 3]| a.into_iter().zip(b).all(|(x, y)| (x - y).abs() < 1e-6);
    let attach = |topology: &references::TopologySnapshot, face: u32, source: EntityId| {
        let entry = topology
            .elements
            .iter()
            .find(|e| {
                e.key.ordinal == face && matches!(e.geometry, Some(AnalyticGeometry::Plane { .. }))
            })
            .unwrap();
        let anchor = entry
            .frames
            .iter()
            .find(|a| a.origin.source == source)
            .unwrap();
        let mut intent = GeometryIntent::capture(topology, entry.key).unwrap();
        intent.attachment = Some(Attachment {
            anchor: anchor.origin.clone(),
            offset_mm: [0.; 3],
            angle_rad: 0.,
        });
        (intent, anchor.frame)
    };
    for (entity, cut) in operations {
        let id = entity.id;
        let groove_id = cut.input.producer;
        let groove_cut = model::read_machining_cut(document, groove_id).unwrap();
        let groove_definition = model::read_machining(document, groove_cut.operation).unwrap();
        assert!(matches!(
            groove_definition.kind,
            MachiningKind::Groove { .. }
        ));
        let last_joint = model::read_panel_joint_cut(document, groove_cut.input.producer).unwrap();
        let joint = model::read_panel_joint(document, last_joint.joint).unwrap();
        assert!(matches!(joint.kind, PanelJointKind::Finger { .. }));
        check_stock_source(document, &parameters, &entity.name, id, &rust);
        check_stock_source(
            document,
            &parameters,
            "grooved runner stock",
            groove_id,
            &rust,
        );

        let definition = model::read_machining(document, cut.operation).unwrap();
        let resolved =
            model::resolve_machining(document, &parameters, &definition, cut.output.body).unwrap();
        assert_eq!(resolved.holes.len(), 4);
        let solid = rust.solid(id).unwrap();
        let grooved = rust.solid(groove_id).unwrap();
        let faces = solid.faces().unwrap();
        let groove_faces = grooved.faces().unwrap();
        // Four blind bores add a wall and a floor each; every other face of
        // the grooved stock is unchanged or keeps its plane and loses
        // exactly the four bore mouths it is drilled through.
        assert_eq!(faces.len(), groove_faces.len() + 8);
        let mouth = std::f64::consts::PI * 1.25_f64.powi(2);
        let mut unchanged = 0;
        let mut pierced = 0.;
        for face in &groove_faces {
            if faces.iter().any(|f| {
                f.plane.is_some() == face.plane.is_some()
                    && same_point(f.center_mm, face.center_mm)
                    && (f.area_mm2 - face.area_mm2).abs() < 1e-6
            }) {
                unchanged += 1;
            } else {
                pierced += face.area_mm2;
            }
        }
        assert_eq!(unchanged, groove_faces.len() - 1, "{}", entity.name);
        let entry = faces
            .iter()
            .filter(|f| {
                !groove_faces.iter().any(|g| {
                    same_point(f.center_mm, g.center_mm) && (f.area_mm2 - g.area_mm2).abs() < 1e-6
                })
            })
            .filter(|f| f.plane.is_some() && (f.area_mm2 - mouth).abs() > 1e-6)
            .map(|f| f.area_mm2)
            .sum::<f64>();
        near(entry, pierced - 4. * mouth);

        let topology = &rust.snapshot().shapes[&id].topology;
        let groove_topology = &rust.snapshot().shapes[&groove_id].topology;
        let groove = model::resolve_machining(
            document,
            &parameters,
            &groove_definition,
            groove_cut.output.body,
        )
        .unwrap();
        let depth = groove.recesses[0].depth_mm;
        let groove_floor = groove_faces
            .iter()
            .find(|face| {
                face.plane.is_some_and(|(point, normal)| {
                    let delta = std::array::from_fn(|i| point[i] - groove.frame.origin[i]);
                    (dot(delta, groove.frame.inward) - depth).abs() < 1e-6
                        && (dot(normal, groove.frame.inward) + 1.).abs() < 1e-6
                })
            })
            .expect("the rectangular groove leaves a planar floor");
        let (groove_intent, groove_frame) = attach(groove_topology, groove_floor.index, groove_id);
        let Resolution::Resolved {
            key,
            frame: Some(frame),
            ..
        } = references::resolve(&groove_intent, topology, 0)
        else {
            panic!("runner drilling lost the earlier groove floor attachment");
        };
        assert_eq!(frame, groove_frame);
        let preserved_floor = faces.iter().find(|face| face.index == key.ordinal).unwrap();
        near(preserved_floor.area_mm2, groove_floor.area_mm2);
        saved_floor_intents.push((id, groove_intent, groove_frame));

        let receipt = solid
            .manufacturing_mesh_with_receipt(MeshOptions::default())
            .unwrap();
        assert!(receipt.edge_pairs > 0 && receipt.vertex_links > 0);
        assert_eq!(receipt.maximum_displacement_mm, 0.);
        let mut removed = 0.;
        for hole in resolved.holes {
            assert!(!hole.through);
            near(hole.diameter_mm, 2.5);
            near(hole.depth_mm, 8.);
            let area = std::f64::consts::PI * (hole.diameter_mm / 2.).powi(2);
            removed += area * hole.depth_mm;
            let center = std::array::from_fn(|i| hole.center[i] + hole.depth_mm * hole.axis[i]);
            let floor = faces
                .iter()
                .find(|face| face.plane.is_some() && same_point(face.center_mm, center))
                .expect("runner screw bore retains a floor at its authored depth");
            near(floor.area_mm2, area);
            near(dot(floor.plane.unwrap().1, hole.axis), -1.);
            assert!(receipt.mesh.face_indices.contains(&floor.index));
            let (intent, frame) = attach(topology, floor.index, id);
            let contour = solid
                .planar_face_contours(
                    floor.index,
                    PlaneFrame {
                        origin: frame.origin_mm,
                        u: frame.u,
                        v: frame.v,
                    },
                )
                .unwrap();
            assert_eq!(contour.loops.len(), 1);
            assert!(matches!(
                contour.loops[0].edges.as_slice(),
                [PlanarEdge::Arc { sweep_rad, .. }]
                    if (sweep_rad.abs() - std::f64::consts::TAU).abs() < 1e-7
            ));
            near(contour.area_mm2(), area);
            saved_floor_intents.push((id, intent, frame));
        }
        near(
            grooved.volume_mm3().unwrap() - solid.volume_mm3().unwrap(),
            removed,
        );
    }
    assert_eq!(check_cut_list(document, &parameters, &rust), 4);

    let mut saved = Vec::new();
    project_file::write(&project, &mut saved, Default::default()).unwrap();
    let reopened = read_bundled_project(&saved);
    let rebuilt = compile_stock_sources(&reopened.document, &roots, GeometryBackend::Rust);
    for (feature, intent, expected) in saved_floor_intents {
        let topology = &rebuilt.snapshot().shapes[&feature].topology;
        let Resolution::Resolved {
            frame: Some(frame), ..
        } = references::resolve(&intent, topology, 0)
        else {
            panic!("runner or groove floor attachment did not survive reopen");
        };
        assert_eq!(frame, expected);
    }
    assert_eq!(same_stock_cut_lists(document, &rust, &rebuilt), 4);
}

#[test]
fn rust_rebuilds_the_bundled_toolcase_stock_panels_with_exact_hole_outlines() {
    let project = read_bundled_project(include_bytes!("../../../examples/toolcase.tcd"));
    let document = &project.document;
    let panels: Vec<_> = model::panels(document).collect();
    assert_eq!(panels.len(), 7);
    let roots: Vec<_> = panels.iter().map(|panel| panel.id).collect();
    let parameters = multicad_core::parameters::evaluate(document).unwrap();
    let rust = compile_stock_sources(document, &roots, GeometryBackend::Rust);
    for panel in panels {
        let id = panel.id;
        check_stock_source(document, &parameters, &panel.name, id, &rust);
        let snap = rust.snapshot();
        let start = snap.shapes[&id]
            .topology
            .elements
            .iter()
            .find(|element| element.frames.iter().any(|f| f.origin.role == "cap.start"))
            .expect("panel start cap has an authored attachment");
        let actual = rust
            .solid(id)
            .unwrap()
            .planar_face_contours(start.key.ordinal, snap.panel_frames[&id])
            .unwrap();
        // Independent expectation: the saved sketch region resolved by core.
        let model::PanelShape::Profile { profile, .. } = model::resolve_panel(
            document,
            &parameters,
            &model::read_panel(document, id).unwrap(),
        )
        .unwrap()
        .shape
        else {
            panic!("{} is a sketch-region panel", panel.name);
        };
        let region = panel_region(document, &profile);
        near(actual.area_mm2(), region.area_mm2);
        assert_eq!(actual.loops.len(), 1 + region.holes.len(), "{}", panel.name);
        let circles = actual
            .loops
            .iter()
            .filter(|wire| {
                !wire.outer
                    && matches!(wire.edges.as_slice(), [PlanarEdge::Arc { sweep_rad, .. }]
                        if (sweep_rad.abs() - std::f64::consts::TAU).abs() < 1e-7)
            })
            .count();
        let expected_circles = match panel.name.as_str() {
            "Left side feature" | "Right side feature" => 10,
            _ => 0,
        };
        assert_eq!(circles, expected_circles, "{}", panel.name);
    }
    // The two side panels have later edge treatments outside these source roots.
    assert_eq!(check_cut_list(document, &parameters, &rust), 5);
}

#[test]
fn rust_stock_backed_profile_panel_keeps_its_round_through_hole_and_cut_dimensions() {
    let mut f = fixture();
    let (sketch, edits) = Source::create(Datum::Xz, "Drilled panel profile").unwrap();
    insert(&mut f.s, edits);
    let mut source = Source::read(f.s.document(), sketch).unwrap().0;
    source.add_rectangle([0., 0.], [120., 70.]).unwrap();
    let center = source.add_point([27., 33.]);
    source.sketch.curves.push(Curve {
        id: EntityId::new(),
        construction: false,
        kind: CurveKind::Circle {
            center,
            radius_mm: 6.,
        },
    });
    let region = sketch_regions::extract(&source)
        .regions
        .into_iter()
        .find(|region| region.holes.len() == 1)
        .unwrap()
        .intent();
    let edits = source.replace(f.s.document(), sketch).unwrap();
    insert(&mut f.s, edits);
    let panel = PanelRequest::new(
        "Drilled panel",
        PanelDefinition {
            output: BodyOutput {
                body: EntityId::new(),
            },
            source: PanelSource::Profile {
                stock: f.ply,
                profile: ProfileSource { sketch, region },
            },
            grain: GrainAxis::V,
            edges: PanelEdges::default(),
            oversize: None,
        },
    );
    let prepared = model::prepare_panel(&f.s, f.s.revision(), actor(), &panel).unwrap();
    checked(&mut f.s, prepared.prepared);
    let rust = Compiler::default()
        .compile_with_backend(f.s.document(), f.s.revision(), GeometryBackend::Rust, None)
        .unwrap();
    assert!(rust.solid(panel.feature).unwrap().uses_rust_geometry());
    let list = rust.cut_list(f.s.document()).unwrap();
    let cut = record(&list, "Drilled panel");
    assert_eq!(cut.stock, f.ply);
    pair(cut.cut.finished_mm, [70., 120.]);
    near(cut.cut.thickness_mm, 18.);
    near(
        cut.volume_mm3,
        (120. * 70. - std::f64::consts::PI * 6. * 6.) * 18.,
    );
    same_cut_records(&list, &compile(&f.s).cut_list(f.s.document()).unwrap());
    let topology = &rust.snapshot().shapes[&panel.feature].topology;
    let round_edges = topology.elements.iter().filter(|e| {
        matches!(e.geometry, Some(AnalyticGeometry::CircularArc { radius_mm, .. }) if radius_mm == 6.)
    }).count();
    assert_eq!(round_edges, 2);
}
