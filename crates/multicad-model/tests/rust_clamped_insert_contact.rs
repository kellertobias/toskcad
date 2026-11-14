//! Real separate clamp/wood inputs; source-only staging never claims a validated connection.
#[path = "support/custom_screw_connections.rs"]
mod custom;
#[path = "../../multicad-mcp/tests/support/screw_connection.rs"]
mod metric;
use multicad_core::{
    commands::{Actor, Channel, Session},
    model::*,
};
use multicad_model::Compiler;
fn request() -> (Session, ScrewConnectionRequest) {
    let (document, mut r) = custom::fixture(ScrewConnectionMode::HangerBolt);
    let (_, insert) = custom::fixture(ScrewConnectionMode::WoodInsert);
    r.mode = ScrewConnectionMode::WoodInsert;
    r.hardware = insert.hardware;
    let (_, screw, _) = metric::fixture();
    *r.hardware
        .iter_mut()
        .find(|row| row.role == ScrewConnectionHardwareRole::Screw)
        .unwrap() = screw
        .hardware
        .into_iter()
        .find(|row| row.role == ScrewConnectionHardwareRole::Screw)
        .unwrap();
    let pilot = r.receivers.last_mut().unwrap();
    pilot.source.cut = ConnectionReceiverCut::Round {
        diameter: length(8.),
    };
    pilot.source.extent = ConnectionReceiverExtent::Blind {
        depth: length(6.),
        floor: length(1.),
    };
    r.receivers[0].source.seat = Some(ConnectionReceiverSeat::Counterbore {
        diameter: length(7.),
        depth: length(1.),
    });
    (Session::new(document).unwrap(), r)
}
#[test]
fn actual_two_stock_sources_expose_independent_clock_and_clearance_contracts() {
    let (session, r) = request();
    let pending = prepare_screw_connection(&session, 0, Actor::local(Channel::Mcp), &r).unwrap();
    let roots: Vec<_> = r
        .receivers
        .iter()
        .map(|row| row.source.output.body)
        .chain(r.hardware.iter().map(|row| row.instance))
        .collect();
    let model = Compiler {
        omit_display: true,
        ..Default::default()
    }
    .compile_sources(pending.prepared.document(), 1, &roots)
    .unwrap();
    assert!(model.snapshot().screw_connection_checks.is_empty());
    let hardware = |role| {
        let row = r.hardware.iter().find(|row| row.role == role).unwrap();
        model
            .leaf_solid(&OccurrencePath {
                root: row.instance,
                members: vec![],
            })
            .unwrap()
    };
    let insert = hardware(ScrewConnectionHardwareRole::Insert);
    let screw = hardware(ScrewConnectionHardwareRole::Screw);
    let data = insert.rust_wood_insert_data().unwrap();
    let (actual, screw_frame) = screw.rust_machine_screw_geometry().unwrap();
    eprintln!(
        "INSERT frame={:?}; SCREW frame={:?}, source={:?}",
        data.frame,
        screw_frame,
        actual.thread().spec()
    );
    assert_eq!(data.source.full_length_mm(), 4.);
    assert_eq!(actual.thread().spec().length_mm, 8.);
    let wood = model.leaf_solid(&r.receivers[1].source.occurrence).unwrap();
    let policy = wood.geometry_policy();
    let req = multicad_kernel::rust_kernel::wood_thread_contact::Requirements {
        minimum_wall_mm: 1.,
        minimum_floor_mm: 1.,
        minimum_engagement_mm: 2.,
        tip_clearance_mm: 1.,
        linear_tolerance_mm: policy.linear_tolerance_mm,
        minimum_feature_mm: policy.minimum_feature_mm,
    };
    let contact = wood.rust_wood_insert_contact(insert, req, None);
    let mating = screw.rust_wood_insert_mating(insert, 2., None);
    let clearance = wood.rust_wood_insert_screw_clearance(screw, 1., None);
    eprintln!(
        "actual contact={:?}; mating={:?}; screw-floor={:?}",
        contact, mating, clearance
    );
    let frame =
        |f: multicad_kernel::PlaneFrame| multicad_kernel::rust_kernel::wood_thread_contact::Frame {
            origin_mm: f.origin,
            axes: [f.u, f.v, f.normal()],
        };
    let new = multicad_kernel::rust_kernel::wood_thread_contact::certify_clocked_insert_mating(
        &data.source,
        frame(data.frame),
        actual,
        frame(screw_frame),
        multicad_kernel::rust_kernel::wood_thread_contact::ClockedMatingOptions {
            requirements: multicad_kernel::rust_kernel::wood_thread_contact::MatingRequirements {
                minimum_engagement_mm: 2.,
                minimum_clearance_mm: 4. * policy.linear_tolerance_mm,
            },
            maximum_cells: 4096,
            linear_tolerance_mm: policy.linear_tolerance_mm,
        },
        &|| false,
    );
    match &new {
        Ok(v) => eprintln!(
            "actual CLOCKED cells={}, radial={}, axial={}, angle={}, mating={:?}",
            v.cells_used, v.radial_map_error_mm, v.axial_map_error_mm, v.angle_error_rad, v.mating
        ),
        Err(e) => eprintln!("actual CLOCKED refusal={e:?}"),
    };
    assert!(new.is_ok());
    let general_floor = wood.rust_clocked_wood_insert_screw_clearance(screw, 1., None);
    eprintln!("actual CLOCKED floor={general_floor:?}");
    assert!(general_floor.is_ok());
    assert!(contact.is_ok());
    assert!(mating.is_err());
    assert!(clearance.is_err());
}
