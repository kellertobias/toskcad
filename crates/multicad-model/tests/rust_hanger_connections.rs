//! Normal Rust connection admission, independent from source-only pose tests.
use multicad_core::{
    commands::{Actor, Channel, Session},
    document::EntityId,
    model::*,
    project_file::{self, ProjectFile},
};
use multicad_model::Compiler;
#[path = "support/custom_screw_connections.rs"]
mod support;
fn actor() -> Actor {
    Actor::local(Channel::Mcp)
}
fn compiler() -> Compiler {
    Compiler {
        omit_display: true,
        ..Default::default()
    }
}
#[test]
fn normal_connection_manufactures_actual_members_and_live_edit_history_reopens() {
    let (doc, r) = support::fixture(ScrewConnectionMode::HangerBolt);
    let mut session = Session::new(doc).unwrap();
    let candidate = prepare_screw_connection(&session, 0, actor(), &r).unwrap();
    let model = Compiler::default()
        .compile(candidate.prepared.document(), 1)
        .unwrap();
    let check = &model.snapshot().screw_connection_checks[0];
    assert!(check.native_component_collision_checked);
    assert_eq!(check.mode, ScrewConnectionMode::HangerBolt);
    assert!(
        check.actual_wood_engagement_mm.unwrap() > 11.99
            && check.actual_engagement_mm.unwrap() > 2.39
    );
    for root in r
        .hardware
        .iter()
        .map(|h| OccurrencePath {
            root: h.instance,
            members: vec![],
        })
        .chain(r.receivers.iter().map(|h| h.source.occurrence.clone()))
    {
        let shape = model.leaf_solid(&root).unwrap();
        let manufacturing = shape
            .manufacturing_mesh_with_receipt(Default::default())
            .unwrap();
        assert!(manufacturing.edge_pairs > 0 && manufacturing.vertex_links > 0);
    }
    session.commit(candidate.prepared).unwrap();
    let row = r
        .hardware
        .iter()
        .find(|h| h.role == ScrewConnectionHardwareRole::NutWasher)
        .unwrap();
    let ScrewConnectionHardwareSource::Generated { part } = &row.source else {
        panic!()
    };
    let mut edit = read_mechanical_part(session.document(), part.feature).unwrap();
    let MechanicalPartParameters::Washer { thickness, .. } = &mut edit.parameters else {
        panic!()
    };
    *thickness = length(0.27);
    let before = session.document().clone();
    let candidate = prepare_mechanical_part_edit(
        &session,
        session.revision(),
        actor(),
        part.feature,
        &edit.parameters,
    )
    .unwrap();
    compiler()
        .compile(candidate.prepared.document(), session.revision() + 1)
        .unwrap();
    session.commit(candidate.prepared).unwrap();
    let after = session.document().clone();
    let mut bytes = vec![];
    project_file::write(
        &ProjectFile::new(after.clone()).unwrap(),
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
    assert!(
        compiler()
            .compile(&opened.project.document, 0)
            .unwrap()
            .snapshot()
            .screw_connection_checks[0]
            .native_component_collision_checked
    );
    session.undo(session.revision(), actor()).unwrap();
    assert_eq!(*session.document(), before);
    compiler()
        .compile(session.document(), session.revision())
        .unwrap();
    session.redo(session.revision(), actor()).unwrap();
    assert_eq!(*session.document(), after);
}
#[test]
fn actual_prefix_floor_and_other_component_failures_leave_the_canonical_document_unchanged() {
    for case in 0..5 {
        let (doc, mut r) = support::fixture(ScrewConnectionMode::HangerBolt);
        let mut session = Session::new(doc).unwrap();
        if case >= 3 {
            let parameters = if case == 3 {
                PrimitiveParameters::Box {
                    size: [length(1.), length(1.), length(1.)],
                    frame: PrimitiveFrame {
                        origin: [length(9.5), length(9.5), length(-8.)],
                        ..Default::default()
                    },
                }
            } else {
                PrimitiveParameters::Cylinder {
                    radius: length(1.),
                    height: length(2.),
                    frame: PrimitiveFrame {
                        origin: [length(9.5), length(9.5), length(-8.)],
                        ..Default::default()
                    },
                }
            };
            let obstacle = PrimitiveRequest::new("Other physical material", parameters);
            let candidate =
                prepare_primitive(&session, session.revision(), actor(), &obstacle).unwrap();
            session.commit(candidate.prepared).unwrap();
        } else if case == 0 {
            let row = r
                .hardware
                .iter_mut()
                .find(|h| h.role == ScrewConnectionHardwareRole::HangerBolt)
                .unwrap();
            let ScrewConnectionHardwareSource::Generated { part } = &mut row.source else {
                panic!()
            };
            let MechanicalPartParameters::CustomConnection { parameters } = &mut part.parameters
            else {
                panic!()
            };
            let CustomConnectionHardwareParameters::HangerBolt { plain_diameter, .. } =
                parameters.as_mut()
            else {
                panic!()
            };
            *plain_diameter = length(3.6);
        } else if case == 1 {
            r.parameters.blind_tip_clearance = length(5.);
        } else {
            let row = r
                .hardware
                .iter_mut()
                .find(|h| h.role == ScrewConnectionHardwareRole::Nut)
                .unwrap();
            let ScrewConnectionHardwareSource::Generated { part } = &mut row.source else {
                panic!()
            };
            let MechanicalPartParameters::HexNut { thread, .. } = &mut part.parameters else {
                panic!()
            };
            thread.hand = ThreadHand::Left;
        }
        let before = session.document().clone();
        match prepare_screw_connection(&session, session.revision(), actor(), &r) {
            Ok(candidate) => {
                let result =
                    compiler().compile(candidate.prepared.document(), session.revision() + 1);
                assert!(result.is_err(), "case{case}");
                if case >= 3 {
                    assert!(
                        result
                            .err()
                            .unwrap()
                            .to_string()
                            .contains("Clearance from another component")
                    );
                }
            }
            Err(_) => assert!(
                case < 3,
                "other actual source must reach material validation"
            ),
        }
        assert_eq!(*session.document(), before);
    }
}
#[test]
fn far_actual_primitive_material_earns_strict_clearance() {
    let frame = PrimitiveFrame {
        origin: [length(100.), length(100.), length(100.)],
        u: [0.8, 0.6, 0.],
        v: [-0.36, 0.48, 0.8],
    };
    for parameters in [
        PrimitiveParameters::Box {
            size: std::array::from_fn(|_| length(2.)),
            frame: frame.clone(),
        },
        PrimitiveParameters::Cylinder {
            radius: length(2.),
            height: length(3.),
            frame: frame.clone(),
        },
        PrimitiveParameters::Sphere {
            radius: length(2.),
            frame,
        },
    ] {
        let (doc, r) = support::fixture(ScrewConnectionMode::HangerBolt);
        let mut session = Session::new(doc).unwrap();
        let other = PrimitiveRequest::new("Other actual stock", parameters);
        let candidate = prepare_primitive(&session, 0, actor(), &other).unwrap();
        session.commit(candidate.prepared).unwrap();
        let candidate =
            prepare_screw_connection(&session, session.revision(), actor(), &r).unwrap();
        assert!(
            compiler()
                .compile(candidate.prepared.document(), session.revision() + 1)
                .unwrap()
                .snapshot()
                .screw_connection_checks[0]
                .native_component_collision_checked
        );
    }
}

#[test]
fn finite_scene_query_budget_refuses_atomically_after_real_source_checks() {
    let (doc, mut r) = support::fixture(ScrewConnectionMode::HangerBolt);
    r.hardware
        .retain(|h| h.role != ScrewConnectionHardwareRole::NutWasher);
    for _ in 0..62 {
        r.hardware.push(ScrewConnectionHardwareRequest {
            role_id: EntityId::new(),
            role: ScrewConnectionHardwareRole::NutWasher,
            name: "Actual thin ring".into(),
            source: ScrewConnectionHardwareSource::Generated {
                part: Box::new(MechanicalPartRequest::new(
                    "Ring",
                    MechanicalPartParameters::Washer {
                        outside_diameter: length(7.),
                        bore_diameter: length(3.4),
                        thickness: length(0.01),
                        frame: Default::default(),
                    },
                )),
            },
            instance: EntityId::new(),
            driver: EntityId::new(),
        });
    }
    let mut session = Session::new(doc).unwrap();
    for i in 0..65 {
        let request = PrimitiveRequest::new(
            "Actual separated stock",
            PrimitiveParameters::Box {
                size: [length(1.), length(1.), length(1.)],
                frame: PrimitiveFrame {
                    origin: [length(100. + 2. * f64::from(i)), length(100.), length(100.)],
                    ..Default::default()
                },
            },
        );
        let p = prepare_primitive(&session, session.revision(), actor(), &request).unwrap();
        session.commit(p.prepared).unwrap();
    }
    let before = session.document().clone();
    let p = prepare_screw_connection(&session, session.revision(), actor(), &r).unwrap();
    let error = compiler()
        .compile(p.prepared.document(), session.revision() + 1)
        .err()
        .unwrap();
    assert!(
        error.to_string().contains("4096 pair query budget"),
        "{error}"
    );
    assert_eq!(*session.document(), before);
}
