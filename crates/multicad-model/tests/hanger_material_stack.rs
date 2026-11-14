use multicad_core::{
    commands::{Actor, Channel, Session},
    document::Document,
    model::*,
    project_file::{self, ProjectFile},
};
use multicad_kernel::rust_kernel::wood_thread_contact::{
    HangerMatingRequirements, HangerStack, PassageRequirements, Requirements, StackRequirements,
};
use multicad_model::Compiler;
#[path = "support/custom_screw_connections.rs"]
mod support;
#[test]
fn actual_canonical_hanger_stack_has_no_related_material_intrusion() {
    let (doc, r) = support::fixture(ScrewConnectionMode::HangerBolt);
    let session = Session::new(doc).unwrap();
    let prepared = prepare_screw_connection(&session, 0, Actor::local(Channel::Mcp), &r).unwrap();
    let model = compile_sources(prepared.prepared.document(), 1, &r);
    check(&model, &r);
}
fn compile_sources(
    doc: &Document,
    revision: u64,
    r: &ScrewConnectionRequest,
) -> multicad_model::Model {
    let roots: Vec<_> = r
        .hardware
        .iter()
        .map(|h| h.instance)
        .chain(r.receivers.iter().map(|h| h.source.output.body))
        .collect();
    // Source-only compilation is diagnostic setup, not connection admission.
    Compiler {
        omit_display: true,
        ..Default::default()
    }
    .compile_sources(doc, revision, &roots)
    .unwrap()
}
fn check(model: &multicad_model::Model, r: &ScrewConnectionRequest) -> HangerStack {
    attempt(model, r).unwrap()
}
fn attempt(
    model: &multicad_model::Model,
    r: &ScrewConnectionRequest,
) -> Result<HangerStack, multicad_kernel::GeometryError> {
    let member = |role| {
        let row = r.hardware.iter().find(|h| h.role == role).unwrap();
        model
            .leaf_solid(&OccurrencePath {
                root: row.instance,
                members: vec![],
            })
            .unwrap()
    };
    let stock = |index: usize| {
        model
            .leaf_solid(&r.receivers[index].source.occurrence)
            .unwrap()
    };
    let hanger = member(ScrewConnectionHardwareRole::HangerBolt);
    let nut = member(ScrewConnectionHardwareRole::Nut);
    let washer = member(ScrewConnectionHardwareRole::NutWasher);
    let req = StackRequirements {
        wood: Requirements {
            minimum_wall_mm: 1.,
            minimum_floor_mm: 1.,
            minimum_engagement_mm: 2.,
            tip_clearance_mm: 1.,
            linear_tolerance_mm: 1e-7,
            minimum_feature_mm: 1e-5,
        },
        minimum_clamp_wall_mm: 1.,
        passage: PassageRequirements {
            minimum_radial_clearance_mm: 1e-5,
            linear_tolerance_mm: 1e-7,
            minimum_feature_mm: 1e-5,
        },
        mating: HangerMatingRequirements {
            minimum_engagement_mm: 2.,
            ..Default::default()
        },
        minimum_bearing_land_mm: 1e-5,
        maximum_bearing_gap_mm: 1e-7,
    };
    stock(1)
        .rust_hanger_wood_contact(hanger, req.wood, None)
        .unwrap_or_else(|e| {
            panic!(
                "wood proof: {e:?}; actual hanger {:?}",
                hanger.rust_hanger_data().unwrap().frame
            )
        });
    stock(0)
        .rust_hanger_through_clearance(hanger, req.minimum_clamp_wall_mm, req.passage, None)
        .unwrap_or_else(|e| panic!("clamp passage: {e:?}"));
    washer
        .rust_hanger_annulus_clearance(hanger, req.passage, None)
        .unwrap_or_else(|e| panic!("washer passage: {e:?}"));
    hanger
        .rust_hanger_nut_mating(nut, req.mating, None)
        .unwrap_or_else(|e| panic!("nut mating: {e:?}"));
    hanger.rust_hanger_stack(stock(1), stock(0), nut, &[washer], req, None)
}

#[test]
fn actual_handed_edited_and_reopened_sources_keep_floor_mating_and_material_proofs() {
    let (doc, mut r) = support::fixture(ScrewConnectionMode::HangerBolt);
    for row in &mut r.hardware {
        let ScrewConnectionHardwareSource::Generated { part } = &mut row.source else {
            panic!()
        };
        match &mut part.parameters {
            MechanicalPartParameters::CustomConnection { parameters } => {
                let CustomConnectionHardwareParameters::HangerBolt {
                    wood_thread,
                    metric_thread,
                    ..
                } = parameters.as_mut()
                else {
                    panic!()
                };
                wood_thread.hand = ThreadHand::Left;
                metric_thread.hand = ThreadHand::Left;
            }
            MechanicalPartParameters::HexNut { thread, .. } => thread.hand = ThreadHand::Left,
            _ => {}
        }
    }
    let mut session = Session::new(doc).unwrap();
    let prepared = prepare_screw_connection(&session, 0, Actor::local(Channel::Mcp), &r).unwrap();
    check(&compile_sources(prepared.prepared.document(), 1, &r), &r);
    assert!(
        Compiler {
            omit_display: true,
            ..Default::default()
        }
        .compile(prepared.prepared.document(), 1)
        .unwrap()
        .snapshot()
        .screw_connection_checks[0]
            .native_component_collision_checked
    );
    session.commit(prepared.prepared).unwrap();
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
        Actor::local(Channel::Mcp),
        part.feature,
        &edit.parameters,
    )
    .unwrap();
    let edited = compile_sources(candidate.prepared.document(), session.revision() + 1, &r);
    assert!(
        Compiler {
            omit_display: true,
            ..Default::default()
        }
        .compile(candidate.prepared.document(), session.revision() + 1)
        .unwrap()
        .snapshot()
        .screw_connection_checks[0]
            .native_component_collision_checked
    );
    let proof = check(&edited, &r);
    assert!(proof.wood.tip_clearance_mm[0] >= 1. && proof.wood.engagement_mm[0] >= 2.);
    let hanger = r
        .hardware
        .iter()
        .find(|h| h.role == ScrewConnectionHardwareRole::HangerBolt)
        .unwrap();
    let shape = edited
        .leaf_solid(&OccurrencePath {
            root: hanger.instance,
            members: vec![],
        })
        .unwrap();
    assert!(
        shape
            .manufacturing_mesh_with_receipt(Default::default())
            .unwrap()
            .edge_pairs
            > 0
    );
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
    check(&compile_sources(&opened.project.document, 0, &r), &r);
    session
        .undo(session.revision(), Actor::local(Channel::Mcp))
        .unwrap();
    assert_eq!(*session.document(), before);
    check(
        &compile_sources(session.document(), session.revision(), &r),
        &r,
    );
    session
        .redo(session.revision(), Actor::local(Channel::Mcp))
        .unwrap();
    assert_eq!(*session.document(), after);
}

#[test]
fn actual_proper_and_large_source_frames_need_no_compiler_geometry_relocation() {
    for case in 0..3 {
        let posed = case != 1;
        let (doc, r) = support::fixture(ScrewConnectionMode::HangerBolt);
        let mut session = Session::new(doc).unwrap();
        let u = if posed { [0.8, 0.6, 0.] } else { [1., 0., 0.] };
        let v = if posed {
            [-0.36, 0.48, 0.8]
        } else {
            [0., 1., 0.]
        };
        // The carrier uses this literal cross. Replacing it by the nearby
        // ideal decimal .64 really changes the authored stock contact plane.
        let n = if case == 2 {
            [0.48, -0.64, 0.6]
        } else {
            [
                u[1] * v[2] - u[2] * v[1],
                u[2] * v[0] - u[0] * v[2],
                u[0] * v[1] - u[1] * v[0],
            ]
        };
        for (i, row) in r.receivers.iter().enumerate() {
            let feature = row.source.input.producer;
            let mut primitive = read_primitive(session.document(), feature).unwrap();
            let PrimitiveParameters::Box { frame, size } = &mut primitive.parameters else {
                panic!()
            };
            let origin = if posed {
                if i == 1 {
                    size[2] = length(32.);
                    n.map(|x| -32. * x)
                } else {
                    [0.; 3]
                }
            } else {
                [0., 0., 1e6 - if i == 1 { 20. } else { 0. }]
            };
            *frame = PrimitiveFrame {
                origin: origin.map(length),
                u,
                v,
            };
            let prepared = prepare_primitive_edit(
                &session,
                session.revision(),
                Actor::local(Channel::Mcp),
                feature,
                &primitive.parameters,
            )
            .unwrap();
            session.commit(prepared.prepared).unwrap();
        }
        let prepared =
            prepare_screw_connection(&session, session.revision(), Actor::local(Channel::Mcp), &r)
                .unwrap();
        let compiled = compile_sources(prepared.prepared.document(), session.revision() + 1, &r);
        let normal = Compiler {
            omit_display: true,
            ..Default::default()
        }
        .compile(prepared.prepared.document(), session.revision() + 1);
        if case == 2 {
            assert!(
                normal.is_err(),
                "original source stock overlap must not be waived"
            );
            // The exact original-stock gap is about -2.2737e-15mm. Hardware
            // authoring must never relocate those stock bodies or waive it.
            assert_eq!(
                attempt(&compiled, &r).unwrap_err().kind,
                multicad_kernel::GeometryErrorKind::Unsupported
            );
        } else {
            check(&compiled, &r);
            assert!(
                normal.unwrap().snapshot().screw_connection_checks[0]
                    .native_component_collision_checked
            );
        }
    }
}
