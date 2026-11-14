use super::*;
use crate::{commands::Session, document::Document};
fn fixture(
    left: bool,
    m6: bool,
) -> (
    Session,
    Vec<model::ScrewConnectionHardwareRequest>,
    Vec<ResolvedMechanicalPart>,
) {
    let session = Session::new(Document::new("Actual canonical seating sources").unwrap()).unwrap();
    let mut hanger = model::default_custom_connection_hanger_bolt();
    let model::MechanicalPartParameters::CustomConnection { parameters } = &mut hanger else {
        panic!()
    };
    let model::CustomConnectionHardwareParameters::HangerBolt {
        wood_thread,
        metric_thread,
        plain_diameter,
        ..
    } = parameters.as_mut()
    else {
        panic!()
    };
    let hand = if left {
        model::ThreadHand::Left
    } else {
        model::ThreadHand::Right
    };
    wood_thread.hand = hand;
    metric_thread.hand = hand;
    if m6 {
        wood_thread.major_diameter = length(12.);
        *plain_diameter = length(6.);
        metric_thread.major_diameter = length(6.);
        metric_thread.pitch = length(1.);
    }
    let nut = model::MechanicalPartParameters::HexNut {
        across_flats: length(if m6 { 10. } else { 5.5 }),
        thread: Box::new(model::ThreadParameters {
            major_diameter: length(if m6 { 6. } else { 3. }),
            pitch: length(if m6 { 1. } else { 0.5 }),
            length: length(if m6 { 4. } else { 2.4 }),
            hand,
            fit: model::ThreadFit::Metric6HInternal,
            ..Default::default()
        }),
    };
    let washer = model::MechanicalPartParameters::Washer {
        bore_diameter: length(if m6 { 6.4 } else { 3.4 }),
        outside_diameter: length(if m6 { 12. } else { 7. }),
        thickness: length(0.2),
        frame: Default::default(),
    };
    let mut parts = Vec::new();
    let mut rows = Vec::new();
    for (role, parameters) in [
        (Role::HangerBolt, hanger),
        (Role::Nut, nut),
        (Role::NutWasher, washer),
    ] {
        parts.push(parameters.resolve(session.parameters()).unwrap());
        rows.push(model::ScrewConnectionHardwareRequest {
            role_id: EntityId::new(),
            role,
            name: "Retained canonical hardware".into(),
            source: HardwareSource::Generated {
                part: Box::new(model::MechanicalPartRequest::new("Hardware", parameters)),
            },
            instance: EntityId::new(),
            driver: EntityId::new(),
        });
    }
    (session, rows, parts)
}

#[test]
fn actual_nonquarter_pose_roll_tracks_handed_lead_and_large_world_rounding() {
    for (left, m6, origin, posed) in [
        (false, false, [10., 10., 4.1], false),
        (true, false, [1e6, -1e6, 1e6], false),
        (true, false, [11., -7., 5.], true),
        (false, true, [10., 10., 4.1], false),
    ] {
        let (session, rows, parts) = fixture(left, m6);
        let axis = RigidPose {
            rotation: if posed {
                [[0.8, -0.36, 0.48], [0.6, 0.48, -0.64], [0., 0.8, 0.6]]
            } else {
                [[1., 0., 0.], [0., -1., 0.], [0., 0., -1.]]
            },
            translation_mm: origin,
        };
        let (authored_hanger_station, hanger) = pose(
            &rows,
            &parts,
            rows[0].role_id,
            4.1,
            axis,
            session.parameters(),
            Options::canonical(),
        )
        .unwrap();
        let (authored_nut_station, nut) = pose(
            &rows,
            &parts,
            rows[1].role_id,
            4.1,
            axis,
            session.parameters(),
            Options::canonical(),
        )
        .unwrap();
        let (_, ring) = pose(
            &rows,
            &parts,
            rows[2].role_id,
            4.1,
            axis,
            session.parameters(),
            Options::canonical(),
        )
        .unwrap();
        let inv = axis.inverse().unwrap();
        let hanger_at = inv.point(hanger.point(parts[0].mount_pose().translation_mm))[2];
        let nut_at = inv.point(nut.point(parts[1].mount_pose().translation_mm))[2];
        let reference = datum(&rows, &parts, authored_hanger_station).unwrap();
        assert_eq!(reference.1, if m6 { 1. } else { 0.5 });
        let expected =
            signed_turn(reference.2, authored_nut_station - reference.0, reference.1).unwrap();
        let local = inv.vector(nut.vector(parts[1].mount_pose().vector([1., 0., 0.])));
        let (s, c) = expected.sin_cos();
        assert!((local[0] - c).abs() < 1e-10 && (local[1] - s).abs() < 1e-10);
        assert!(nut_at < if m6 { -4.2 } else { -2.6 });
        assert!(inv.point(ring.translation_mm)[2] < -0.2);
        // Exact cardinal wood datum remains unchanged: no new floor intrusion.
        if !posed && origin[2] == 4.1 {
            assert_eq!(hanger.translation_mm[2], origin[2] - 4.1);
        } else {
            assert!((4.1..=4.1000001).contains(&hanger_at));
        }
    }
}

#[test]
fn numerical_budget_precision_and_wrong_typed_sources_refuse_without_mutation() {
    let (session, rows, parts) = fixture(false, false);
    let before_rows = rows.clone();
    let before_parts = parts.clone();
    let axis = RigidPose {
        translation_mm: [10., 10., 4.1],
        rotation: [[1., 0., 0.], [0., -1., 0.], [0., 0., -1.]],
    };
    for options in [
        Options {
            maximum_attempts: 1,
            ..Options::canonical()
        },
        Options {
            maximum_gap_mm: 1e-20,
            ..Options::canonical()
        },
        Options {
            maximum_attempts: 0,
            ..Options::canonical()
        },
        Options {
            maximum_attempts: 5,
            ..Options::canonical()
        },
        Options {
            maximum_gap_mm: f64::NAN,
            ..Options::canonical()
        },
    ] {
        assert!(
            pose(
                &rows,
                &parts,
                rows[1].role_id,
                4.1,
                axis,
                session.parameters(),
                options
            )
            .is_err()
        );
    }
    assert!(
        pose(
            &rows,
            &parts,
            EntityId::new(),
            4.1,
            axis,
            session.parameters(),
            Options::canonical()
        )
        .is_err()
    );
    assert!(
        pose(
            &rows,
            &parts[..2],
            rows[1].role_id,
            4.1,
            axis,
            session.parameters(),
            Options::canonical()
        )
        .is_err()
    );
    let mut wrong = rows.clone();
    wrong[2].role = Role::HeadWasher;
    assert!(
        pose(
            &wrong,
            &parts,
            wrong[1].role_id,
            4.1,
            axis,
            session.parameters(),
            Options::canonical()
        )
        .is_err()
    );
    assert_eq!(rows, before_rows);
    assert_eq!(parts, before_parts);
    assert!(
        extent(
            &model::default_machine_screw()
                .resolve(session.parameters())
                .unwrap()
        )
        .is_err()
    );
    assert!(
        extent(
            &model::default_custom_connection_wood_insert()
                .resolve(session.parameters())
                .unwrap()
        )
        .is_err()
    );
    assert!(datum(&rows[..1], &parts[1..2], 0.).is_err());
    let mut wrong_parts = parts.clone();
    wrong_parts[0] = model::default_custom_connection_wood_insert()
        .resolve(session.parameters())
        .unwrap();
    assert!(datum(&rows, &wrong_parts, 0.).is_err());
    assert!(
        pose(
            &rows,
            &parts,
            rows[0].role_id,
            4.1,
            RigidPose {
                rotation: [[0.; 3]; 3],
                ..Default::default()
            },
            session.parameters(),
            Options::canonical()
        )
        .is_err()
    );
}
#[test]
fn directed_projection_arithmetic_covers_signed_corners_and_exact_identities() {
    for x in [-3., 0., 1., 2.6] {
        assert_eq!(I::p(x).sub(I::p(x)).abs(), 0.);
        assert_eq!(I::p(0.).mul(I::p(x)).abs(), 0.);
        assert_eq!(I::p(1.).mul(I::p(x)).midpoint(), x);
        assert_eq!(I::p(-1.).mul(I::p(x)).midpoint(), -x);
        assert_eq!(I::p(x).mul(I::p(-1.)).midpoint(), -x);
    }
    let a = I { lo: -0.3, hi: 0.4 };
    let b = I { lo: -2., hi: 3. };
    let product = a.mul(b);
    for x in [a.lo, a.hi] {
        for y in [b.lo, b.hi] {
            assert!(product.lo <= x * y && product.hi >= x * y);
        }
    }
    let denominator = I { lo: 0.9, hi: 1.1 };
    let ratio = a.axis_div(denominator);
    for x in [a.lo, a.hi] {
        for y in [denominator.lo, denominator.hi] {
            assert!(ratio.lo <= x / y && ratio.hi >= x / y);
        }
    }
    let axis = RigidPose {
        translation_mm: [1e6, 1e6, 1e6],
        ..Default::default()
    };
    let pose = RigidPose {
        translation_mm: [1e6, 1e6, 1e6 + 2.6],
        ..Default::default()
    };
    axis.validate().unwrap();
    pose.validate().unwrap();
    let measured = residual(axis, pose, RigidPose::default(), I::p(2.6), 3., 12.);
    assert!(measured > 0. && measured < 1e-9);
}

#[test]
fn public_bound_driver_uses_authored_spacing_and_keeps_exact_wood_datum() {
    use crate::{
        commands::{Actor, Channel},
        references::{
            AnalyticGeometry, AnchoredFrame, Origin, ReferenceFrame, TopologyElement, TopologyKey,
            TopologyKind, TopologySnapshot,
        },
    };
    let (mut session, rows, _) = fixture(false, false);
    let actor = || Actor::local(Channel::Mcp);
    let clamp = model::PrimitiveRequest::new(
        "Clamp",
        model::PrimitiveParameters::Box {
            size: [length(20.), length(20.), length(4.1)],
            frame: Default::default(),
        },
    );
    let wood = model::PrimitiveRequest::new(
        "Wood",
        model::PrimitiveParameters::Box {
            size: std::array::from_fn(|_| length(20.)),
            frame: model::PrimitiveFrame {
                origin: [length(0.), length(0.), length(-20.)],
                ..Default::default()
            },
        },
    );
    for part in [&clamp, &wood] {
        let q = model::prepare_primitive(&session, session.revision(), actor(), part).unwrap();
        session.commit(q.prepared).unwrap();
    }
    // This explicit synthetic snapshot tests Core reference/driver routing;
    // the separate Model test obtains its topology from actual Native stock.
    let frame = ReferenceFrame {
        origin_mm: [0., 0., 4.1],
        u: [1., 0., 0.],
        v: [0., 1., 0.],
    };
    let origin = Origin {
        source: clamp.feature,
        role: "cap.end".into(),
    };
    let key = TopologyKey {
        kind: TopologyKind::Face,
        ordinal: 0,
    };
    let mut topology = TopologySnapshot {
        owner: clamp.feature,
        revision: session.revision() + 1,
        elements: vec![TopologyElement {
            key,
            origins: BTreeSet::from([origin.clone()]),
            boundary: Default::default(),
            frames: vec![AnchoredFrame { origin, frame }],
            geometry: Some(AnalyticGeometry::Plane { frame }),
        }],
    };
    let face = model::capture_profile_face(&topology, key, None).unwrap();
    let axis = EntityId::new();
    let mut source = model::ScrewConnectionAxisDefinition::new(
        model::BodyInput {
            body: clamp.body,
            producer: clamp.feature,
        },
        model::OccurrencePath {
            root: clamp.body,
            members: vec![],
        },
        face.intent.id,
    );
    source.reverse = true;
    source.offset = [length(10.), length(10.), length(0.)];
    let request = ScrewConnectionRequest {
        connection: EntityId::new(),
        name: "Actual canonical Hanger sources".into(),
        mode: Mode::HangerBolt,
        axis: model::ScrewConnectionAxisRequest {
            id: axis,
            source,
            face,
        },
        hardware: rows,
        receivers: [(&clamp, false), (&wood, true)]
            .into_iter()
            .map(|(part, blind)| model::ConnectionReceiverRequest {
                feature: EntityId::new(),
                name: "Retained receiver".into(),
                source: model::ConnectionReceiverDefinition {
                    output: model::BodyOutput { body: part.body },
                    input: model::BodyInput {
                        body: part.body,
                        producer: part.feature,
                    },
                    occurrence: model::OccurrencePath {
                        root: part.body,
                        members: vec![],
                    },
                    axis,
                    station: length(if blind { 4.1 } else { 0. }),
                    entry_reference: None,
                    cut: model::ConnectionReceiverCut::Round {
                        diameter: length(if blind { 4. } else { 3.4 }),
                    },
                    extent: if blind {
                        model::ConnectionReceiverExtent::Blind {
                            depth: length(16.),
                            floor: length(1.),
                        }
                    } else {
                        model::ConnectionReceiverExtent::Through {}
                    },
                    minimum_wall: length(1.),
                    seat: None,
                },
            })
            .collect(),
        parameters: model::ScrewConnectionParameters {
            grip: length(4.1),
            end_reference: None,
            minimum_engagement: length(2.),
            blind_tip_clearance: length(1.),
        },
    };
    let q =
        model::prepare_screw_connection(&session, session.revision(), actor(), &request).unwrap();
    session.commit(q.prepared).unwrap();
    topology.revision = session.revision();
    let snapshots = BTreeMap::from([(clamp.feature, topology)]);
    for role in [Role::HangerBolt, Role::Nut] {
        let row = request.hardware.iter().find(|h| h.role == role).unwrap();
        let actual = model::resolve_instance_pose_with_snapshots(
            session.document(),
            row.instance,
            session.parameters(),
            session.revision(),
            &snapshots,
        )
        .unwrap();
        if role == Role::HangerBolt {
            assert_eq!(actual.translation_mm[2], 0.);
        } else {
            assert!(actual.translation_mm[2] > 6.7);
        }
    }
}
