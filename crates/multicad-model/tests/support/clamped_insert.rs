//! Actual separate counterbored clamp and blind wood insert inputs.
#[path = "custom_screw_connections.rs"]
mod custom;
#[path = "../../../multicad-mcp/tests/support/screw_connection.rs"]
mod metric;
use multicad_core::{commands::Session, model::*};
pub fn request() -> (Session, ScrewConnectionRequest) {
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
