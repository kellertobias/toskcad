//! Actual finite head/cavity bearing, not metadata or proximity to a proposal.
use super::*;
use multicad_kernel::{GeometryPolicy, rust_kernel::conical_screw_seat as source};

pub(super) struct Proof {
    station_mm: f64,
    actual: source::BetweenFrameCertificate,
}
impl Proof {
    pub(super) fn station_mm(&self) -> f64 {
        self.station_mm
    }
    pub(super) fn maximum_gap_mm(&self) -> f64 {
        self.actual.maximum_physical_gap_mm
    }
}

pub(super) fn prove(
    screw: &Solid,
    receiver: &Solid,
    receiver_feature: EntityId,
    bore: rust_through::ThroughReceiverProof,
    policy: GeometryPolicy,
    connection: EntityId,
) -> Result<Proof, ModelError> {
    let fail = |message: &str| error(Some(connection), message);
    let (head, frame) = screw
        .rust_profiled_screw_geometry()
        .ok_or_else(|| fail("Conical bearing requires the actual retained screw head"))?;
    let requested = bore
        .head_seat
        .filter(|seat| seat.conical)
        .ok_or_else(|| fail("Conical bearing requires a current conical receiver"))?;
    let seats = receiver
        .rust_countersink_geometry()
        .ok_or_else(|| fail("The receiving part has no retained finite cone"))?;
    let expected_origin = connection_receivers::origin(receiver_feature, "seat");
    let matches: Vec<_> = seats
        .iter()
        .filter(|seat| seat.origin == expected_origin)
        .collect();
    let [seat] = matches.as_slice() else {
        return Err(fail(
            "The physical receiver has no unique current head-seat support",
        ));
    };
    let options = source::Options {
        tolerance_mm: policy.linear_tolerance_mm,
        minimum_land_mm: 5.
            * policy
                .minimum_feature_mm
                .max(16. * policy.linear_tolerance_mm),
    };
    let actual = source::certify_between_frames(
        head.profile(),
        source::Seat {
            bore_radius_mm: seat.bore_radius_mm,
            mouth_radius_mm: seat.mouth_radius_mm,
            depth_mm: seat.depth_mm,
        },
        options,
        source::Placement {
            origin_mm: frame.origin,
            axes: [frame.u, frame.v, frame.normal()],
            axial_offset_mm: -head.spec().shank_length_mm,
        },
        source::Frame {
            origin_mm: seat.frame.origin,
            axes: [seat.frame.u, seat.frame.v, seat.frame.normal()],
            offset_mm: seat.offset_mm,
            source_error_mm: seat.source_error_mm,
        },
    )
    .map_err(|cause| {
        fail(&format!(
            "Actual finite conical bearing cannot be certified: {cause:?}"
        ))
    })?;
    // The generated pose uses this requested engineering chart. Physical
    // clearance above was independently checked against the actual receiver.
    let nominal = source::certify(
        head.profile(),
        source::Seat {
            bore_radius_mm: bore.radius_mm,
            mouth_radius_mm: requested.radius_mm,
            depth_mm: requested.depth_mm,
        },
        options,
    )
    .map_err(|cause| {
        fail(&format!(
            "Generated head datum cannot be certified: {cause:?}"
        ))
    })?;
    Ok(Proof {
        station_mm: nominal.safe_neck_z_mm,
        actual,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::connection_receivers::rust_through::HeadSeatProof;
    use multicad_kernel::rust_kernel::{head_profile::Head, profiled_screw::Spec};
    use multicad_kernel::{CircularPocketTool, HelicalThread, RigidTransform};
    #[test]
    fn actual_head_and_receiving_material_qualify_after_pose_and_refuse_wrong_or_floating_supports()
    {
        let policy = GeometryPolicy::default();
        let receiver_id = EntityId::new();
        let axis = PlaneFrame {
            origin: [10., 10., 4.1],
            u: [1., 0., 0.],
            v: [0., -1., 0.],
        };
        let stock = Solid::box_solid([20., 20., 4.1], policy).unwrap();
        let tool = CircularPocketTool {
            frame: PlaneFrame {
                origin: [10., 10., 4.10001],
                ..axis
            },
            center_uv: [0.; 2],
            radius_mm: 1.7,
            displacement_mm: [0., 0., -4.10002],
            origin: connection_receivers::origin(receiver_id, "bore"),
        };
        let (receiver, mouths_mm) = stock
            .rust_connection_countersink(
                &tool,
                3.5,
                1.8,
                1.,
                connection_receivers::origin(receiver_id, "seat"),
                policy,
            )
            .unwrap();
        let bore = rust_through::ThroughReceiverProof {
            mouths_mm,
            radius_mm: 1.7,
            protected_radius_mm: 2.7,
            head_seat: Some(HeadSeatProof {
                conical: true,
                radius_mm: 3.5,
                depth_mm: 1.8,
            }),
        };
        let screw = Solid::rust_profiled_machine_screw(
            HelicalThread {
                major_radius_mm: 1.5,
                pitch_mm: 0.5,
                length_mm: 6.,
                depth_mm: 0.3,
                crest_width_mm: 0.0625,
                lower_flank_deg: 30.,
                upper_flank_deg: 30.,
                radial_taper: 0.,
                starts: 1,
                left_hand: false,
            },
            None,
            Spec {
                head: Head::Countersunk {
                    radius_mm: 3.,
                    height_mm: 1.65,
                    included_angle_deg: 90.,
                    edge_chamfer_mm: 0.1,
                    neck_radius_mm: 0.1,
                },
                shank_length_mm: 2.,
                socket_across_flats_mm: 2.5,
                socket_depth_mm: 0.8,
                minimum_wall_mm: 0.1,
            },
            EntityId::new(),
            PlaneFrame::default(),
            policy,
        )
        .unwrap();
        let proposal = source::certify(
            screw.rust_profiled_screw_geometry().unwrap().0.profile(),
            source::Seat {
                bore_radius_mm: 1.7,
                mouth_radius_mm: 3.5,
                depth_mm: 1.8,
            },
            source::Options {
                tolerance_mm: policy.linear_tolerance_mm,
                minimum_land_mm: policy.minimum_feature_mm,
            },
        )
        .unwrap();
        let screw = screw
            .transformed(RigidTransform {
                rotation: [[1., 0., 0.], [0., -1., 0.], [0., 0., -1.]],
                translation_mm: [10., 10., 4.1 - proposal.safe_neck_z_mm],
            })
            .unwrap();
        for pose in [
            RigidTransform::default(),
            RigidTransform {
                rotation: [[0.8, -0.36, 0.48], [0.6, 0.48, -0.64], [0., 0.8, 0.6]],
                translation_mm: [1e6, -2e6, 3e6],
            },
        ] {
            let head = screw.transformed(pose).unwrap();
            let cavity = receiver.transformed(pose).unwrap();
            let proof = prove(&head, &cavity, receiver_id, bore, policy, EntityId::new()).unwrap();
            assert!(
                proof.maximum_gap_mm() > 0. && proof.maximum_gap_mm() <= policy.linear_tolerance_mm
            );
            assert!((proof.station_mm() - proposal.safe_neck_z_mm).abs() < 1e-12);
            assert!(
                prove(
                    &head,
                    &cavity,
                    EntityId::new(),
                    bore,
                    policy,
                    EntityId::new()
                )
                .is_err()
            );
            let moved = head
                .transformed(RigidTransform {
                    translation_mm: [0.001, 0., 0.],
                    ..Default::default()
                })
                .unwrap();
            assert!(prove(&moved, &cavity, receiver_id, bore, policy, EntityId::new()).is_err());
        }
        for translation_mm in [[0., 0., -0.001], [0., 0., 0.001]] {
            let moved = screw
                .transformed(RigidTransform {
                    translation_mm,
                    ..Default::default()
                })
                .unwrap();
            assert!(
                prove(
                    &moved,
                    &receiver,
                    receiver_id,
                    bore,
                    policy,
                    EntityId::new()
                )
                .is_err()
            );
        }
        assert!(prove(&screw, &stock, receiver_id, bore, policy, EntityId::new()).is_err());
    }
}
