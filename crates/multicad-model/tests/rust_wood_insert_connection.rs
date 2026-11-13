//! A real insert connection must pass before negative candidates are meaningful.
#[path = "support/custom_screw_connections.rs"]
mod support;
use multicad_core::{
    commands::{Actor, Channel, Session},
    model::*,
};
use multicad_model::{Compiler, GeometryBackend};
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
fn retained_insert_connection_has_distinct_wood_and_metric_engagement_and_history() {
    let (document, request) = support::fixture(ScrewConnectionMode::WoodInsert);
    let mut session = Session::new(document.clone()).unwrap();
    let prepared = prepare_screw_connection(&session, 0, actor(), &request).unwrap();
    let model = compiler()
        .compile_with_backend(prepared.prepared.document(), 1, GeometryBackend::Rust, None)
        .unwrap();
    let [check] = model.snapshot().screw_connection_checks.as_slice() else {
        panic!("one actual connection")
    };
    assert!((check.actual_wood_engagement_mm.unwrap() - 4.).abs() < 1e-7);
    assert!((check.actual_engagement_mm.unwrap() - 3.).abs() < 1e-7);
    assert!(check.native_component_collision_checked);
    assert_eq!(*session.document(), document);
    session.commit(prepared.prepared).unwrap();
    let committed = session.document().clone();
    session.undo(session.revision(), actor()).unwrap();
    assert_eq!(*session.document(), document);
    session.redo(session.revision(), actor()).unwrap();
    assert_eq!(*session.document(), committed);
}
#[test]
fn actual_insert_refuses_invalid_pilot_engagement_and_wrong_thread_without_committing() {
    for case in 0..6 {
        let (document, mut request) = support::fixture(ScrewConnectionMode::WoodInsert);
        let session = Session::new(document.clone()).unwrap();
        match case {
            0 => {
                request.receivers[0].source.cut = ConnectionReceiverCut::Round {
                    diameter: length(7.5),
                }
            }
            1 => {
                request.receivers[0].source.cut = ConnectionReceiverCut::Round {
                    diameter: length(10.5),
                }
            }
            2 => {
                request.receivers[0].source.extent = ConnectionReceiverExtent::Blind {
                    depth: length(4.5),
                    floor: length(1.),
                }
            }
            3 => request.parameters.minimum_engagement = length(3.5),
            4 => request.parameters.blind_tip_clearance = length(2.5),
            _ => {
                let row = request
                    .hardware
                    .iter_mut()
                    .find(|h| h.role == ScrewConnectionHardwareRole::Screw)
                    .unwrap();
                let ScrewConnectionHardwareSource::Generated { part } = &mut row.source else {
                    panic!()
                };
                let MechanicalPartParameters::MachineScrew { thread, .. } = &mut part.parameters
                else {
                    panic!()
                };
                thread.hand = ThreadHand::Left;
            }
        }
        let prepared = prepare_screw_connection(&session, 0, actor(), &request).unwrap();
        assert!(
            compiler()
                .compile_with_backend(prepared.prepared.document(), 1, GeometryBackend::Rust, None)
                .is_err(),
            "case{case}"
        );
        assert_eq!(*session.document(), document);
    }
}

#[test]
fn screw_can_continue_beyond_through_insert_when_actual_pilot_floor_still_clears() {
    let (document, mut request) = support::fixture(ScrewConnectionMode::WoodInsert);
    let session = Session::new(document.clone()).unwrap();
    let row = request
        .hardware
        .iter_mut()
        .find(|h| h.role == ScrewConnectionHardwareRole::Screw)
        .unwrap();
    let ScrewConnectionHardwareSource::Generated { part } = &mut row.source else {
        panic!()
    };
    let MechanicalPartParameters::MachineScrew {
        underhead_length,
        thread,
        ..
    } = &mut part.parameters
    else {
        panic!()
    };
    *underhead_length = Some(length(4.5));
    thread.length = length(4.5);
    thread.metric_fit_engagement = Some(MetricFitEngagement::NormalGroupAssumed {});
    let prepared = prepare_screw_connection(&session, 0, actor(), &request).unwrap();
    let model = compiler()
        .compile_with_backend(prepared.prepared.document(), 1, GeometryBackend::Rust, None)
        .unwrap();
    let check = &model.snapshot().screw_connection_checks[0];
    assert!((check.actual_engagement_mm.unwrap() - 4.).abs() < 1e-7);
    assert!(check.native_component_collision_checked);
    assert_eq!(*session.document(), document);
}
