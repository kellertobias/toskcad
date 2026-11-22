//! Branch controls for the Rust rotary/rack validator: each accepted mesh
//! carries its witness, and wrong phase, separation and unsupported member
//! combinations fail with their own diagnostics.
use super::*;
use multicad_core::{
    commands::{Actor, Channel, Session},
    model::*,
};

fn count(value: f64) -> Dimension {
    Dimension::Literal {
        value,
        unit: "count".into(),
    }
}

fn try_session(request: &GearPairRequest) -> Option<Session> {
    let mut session = Session::new(Document::new("Rust pair branches").unwrap()).unwrap();
    let prepared = prepare_gear_pair(&session, 0, Actor::local(Channel::Ui), request).ok()?;
    session.commit(prepared.prepared).unwrap();
    Some(session)
}

fn shifted(family: GearFamily) -> GearPairRequest {
    let mut recipes = std::array::from_fn(|i| GearParameters {
        family,
        hand: if family == GearFamily::Spur {
            GearHand::None
        } else if i == 0 {
            GearHand::Right
        } else {
            GearHand::Left
        },
        helix_angle: angle(if family == GearFamily::Spur { 0. } else { 25. }),
        center_gap: length(if family == GearFamily::DoubleHelical {
            0.5
        } else {
            0.
        }),
        module: length(1.),
        teeth: count([24., 36.][i]),
        profile_shift: ratio([0.2, 0.3][i]),
        bore_diameter: length(2.),
        face_width: length(2.),
        minimum_wall: length(0.5),
        hub: Some(GearHub {
            diameter: length(8.),
            axial_start: length(-1.),
            axial_end: length(3.),
        }),
        ..Default::default()
    });
    couple_shifted_gear_pair_recipes(&mut recipes).unwrap();
    GearPairRequest::new("Pair", recipes.map(|p| GearRequest::new("Member", p)))
}

fn rack_request() -> GearPairRequest {
    let shared = GearParameters {
        module: length(1.),
        face_width: length(2.),
        ..Default::default()
    };
    let rack = GearParameters {
        family: GearFamily::Rack,
        rack: Some(RackDimensions {
            backing_thickness: length(3.),
        }),
        teeth: count(20.),
        bore_diameter: length(0.),
        ..shared.clone()
    };
    GearPairRequest::new(
        "Rack",
        [
            GearRequest::new("Pinion", shared),
            GearRequest::new("Rack", rack),
        ],
    )
}

struct Fixture {
    session: Session,
    source: GearPairDefinition,
    resolved: model::ResolvedGearPair,
    solids: [Solid; 2],
    pair: EntityId,
}

fn fixture(request: &GearPairRequest) -> Fixture {
    let session = try_session(request).expect("valid pair");
    let document = session.document();
    let source = read_gear_pair(document, request.pair).unwrap();
    let resolved =
        resolve_gear_pair_parameters(document, session.parameters(), &source.parameters).unwrap();
    let solids = source.parameters.gears.map(|feature| {
        let recipe = read_gear(document, feature).unwrap();
        crate::gears::build_rust(
            feature,
            &recipe,
            session.parameters(),
            GeometryPolicy::default(),
        )
        .unwrap()
    });
    Fixture {
        session,
        source,
        resolved,
        solids,
        pair: request.pair,
    }
}

impl Fixture {
    /// Unit vector from the first member's bounds centre to the second's.
    fn apart(&self) -> [f64; 3] {
        let [a, b] = self.placed(0., [RigidTransform::default(); 2]);
        let center = |s: &Solid| {
            let bounds = s.bounds().unwrap();
            std::array::from_fn::<f64, 3, _>(|i| (bounds.min[i] + bounds.max[i]) / 2.)
        };
        let [p, q] = [center(&a), center(&b)];
        let delta: [f64; 3] = std::array::from_fn(|i| q[i] - p[i]);
        let norm = delta.iter().map(|v| v * v).sum::<f64>().sqrt();
        delta.map(|v| v / norm)
    }

    fn pose(&self) -> Pose<'static> {
        let (mesh, engagement) = classify(
            self.session.document(),
            self.pair,
            &self.source,
            &self.resolved,
        )
        .unwrap();
        Pose {
            pair: self.pair,
            mesh,
            engagement,
            policy: GeometryPolicy::default(),
            interrupt: None,
        }
    }

    fn placed(&self, drive: f64, extra: [RigidTransform; 2]) -> [Solid; 2] {
        self.members([drive; 2], extra)
    }

    /// Place each member at its own sampled drive, so unequal drives put
    /// the pair out of phase.
    fn members(&self, drives: [f64; 2], extra: [RigidTransform; 2]) -> [Solid; 2] {
        std::array::from_fn(|i| {
            let sample =
                sample_document(self.session.document(), self.pair, &self.source, drives[i])
                    .unwrap();
            let pose = model::resolve_instance_pose(
                &sample,
                self.source.outputs.instances[i],
                self.session.parameters(),
            )
            .unwrap();
            placed(&self.solids[i], pose, self.pair)
                .unwrap()
                .transformed(extra[i])
                .unwrap()
        })
    }
}

fn shift(by: [f64; 3]) -> RigidTransform {
    RigidTransform {
        translation_mm: by,
        ..Default::default()
    }
}

#[test]
fn shifted_external_families_carry_named_flank_witnesses_and_refuse_bad_poses() {
    for (family, kind, chevron) in [
        (GearFamily::Spur, "shifted_spur_flank_faces", false),
        (GearFamily::Helical, "shifted_helical_flank_faces", false),
        (GearFamily::Herringbone, "shifted_chevron_both_bands", true),
        (
            GearFamily::DoubleHelical,
            "shifted_chevron_both_bands",
            true,
        ),
    ] {
        let fixture = fixture(&shifted(family));
        let pose = fixture.pose();
        assert!(pose.mesh == Mesh::External(family));
        let engagement = pose.engagement.expect("shifted pairs need engagement");
        assert_eq!(engagement.kind, kind);
        assert_eq!(engagement.chevron, chevron);
        assert!((engagement.clearance_mm - 0.1).abs() < 1e-9);
        assert!(scope(pose.mesh, true).contains("profile-shifted"));
        let none = [RigidTransform::default(); 2];
        for drive in [0., TAU / 24. / 2.] {
            let [a, b] = fixture.placed(drive, none);
            let witness = pose.check([&a, &b], drive).unwrap().expect("witness");
            assert_eq!(witness.kind, kind);
            assert!(witness.faces.is_some());
            assert!(witness.distance_mm <= 0.1 + 1e-6, "{family:?}");
        }
        // Moving the second member along its line of centres by a module
        // keeps it clear but removes all flank proximity.
        let [a, b] = fixture.placed(0., [RigidTransform::default(), shift(fixture.apart())]);
        let failure = pose.check([&a, &b], 0.).unwrap_err();
        assert!(
            failure.message.contains("flank proximity"),
            "{}",
            failure.message
        );
        // Half a tooth of the second member out of phase shares material.
        let [a, b] = fixture.members([0., TAU / 24. / 2.], [RigidTransform::default(); 2]);
        let failure = pose.check([&a, &b], 0.).unwrap_err();
        refused_pose(&failure);
    }
}

#[test]
fn unshifted_pairs_need_no_witness_and_wrong_phase_is_an_overlap() {
    let mut request = shifted(GearFamily::Spur);
    for gear in &mut request.gears {
        gear.parameters.profile_shift = ratio(0.);
        gear.parameters.shift_mate = None;
    }
    let fixture = fixture(&request);
    let pose = fixture.pose();
    assert!(pose.engagement.is_none());
    assert!(scope(pose.mesh, false).contains("unshifted external spur"));
    let [a, b] = fixture.placed(0., [RigidTransform::default(); 2]);
    assert!(pose.check([&a, &b], 0.).unwrap().is_none());
    let [a, b] = fixture.members([0., TAU / 24. / 2.], [RigidTransform::default(); 2]);
    let failure = pose.check([&a, &b], 0.).unwrap_err();
    refused_pose(&failure);
}

#[test]
fn rack_pairs_need_rack_flank_proximity_at_samples_and_endpoints() {
    let fixture = fixture(&rack_request());
    let pose = fixture.pose();
    assert!(pose.mesh == Mesh::Rack);
    let engagement = pose.engagement.expect("racks need engagement");
    assert_eq!(engagement.kind, "rack_pinion_flanks");
    assert!(scope(Mesh::Rack, true).contains("straight rack"));
    for drive in cycle_phases(&fixture.resolved) {
        let [a, b] = fixture.placed(drive, [RigidTransform::default(); 2]);
        let witness = pose.check([&a, &b], drive).unwrap().expect("rack witness");
        assert_eq!(witness.kind, "rack_pinion_flanks");
        assert!(witness.distance_mm <= engagement.clearance_mm + 1e-6);
    }
    // Lifting the pinion away from the rack clears but has no engagement.
    let away = fixture.apart().map(|v| -v);
    let [a, b] = fixture.placed(0., [shift(away), RigidTransform::default()]);
    let failure = pose.check([&a, &b], 0.).unwrap_err();
    assert!(
        failure.message.contains("flank proximity"),
        "{}",
        failure.message
    );
    // Half a pinion tooth puts a pinion tooth onto a rack tooth.
    let [a, b] = fixture.members([PI / 24., 0.], [RigidTransform::default(); 2]);
    let failure = pose.check([&a, &b], 0.).unwrap_err();
    refused_pose(&failure);
}

#[test]
fn unsupported_member_combinations_keep_a_typed_refusal() {
    // A helical pinion on a straight rack has no Rust rack certificate.
    let mut request = rack_request();
    request.gears[0].parameters.family = GearFamily::Helical;
    request.gears[0].parameters.hand = GearHand::Right;
    request.gears[0].parameters.helix_angle = angle(20.);
    let rack = refused(&request, "straight rack with a spur pinion");
    // An internal ring pinion with a hub needs a certificate not provided.
    let ring = GearParameters {
        family: GearFamily::InternalRing,
        module: length(1.),
        face_width: length(2.),
        teeth: count(48.),
        bore_diameter: length(0.),
        outer_diameter: Some(length(56.)),
        ..Default::default()
    };
    let pinion = GearParameters {
        module: length(1.),
        face_width: length(2.),
        bore_diameter: length(2.),
        hub: Some(GearHub {
            diameter: length(8.),
            axial_start: length(0.),
            axial_end: length(4.),
        }),
        ..Default::default()
    };
    let request = GearPairRequest::new(
        "Hubbed internal",
        [
            GearRequest::new("Ring", ring),
            GearRequest::new("Pinion", pinion),
        ],
    );
    let ring = refused(&request, "hubless spur pinion");
    // Shared source validation already rejects a helical rack pinion; the
    // hubbed ring pinion is accepted there and must meet the Rust refusal.
    assert!(!rack && ring);
}

/// A recipe the shared source validation already rejects never reaches the
/// Rust validator; one it accepts must meet the typed refusal.
fn refused(request: &GearPairRequest, diagnostic: &str) -> bool {
    let Some(session) = try_session(request) else {
        return false;
    };
    let document = session.document();
    let source = read_gear_pair(document, request.pair).unwrap();
    let resolved =
        resolve_gear_pair_parameters(document, session.parameters(), &source.parameters).unwrap();
    let failure = classify(document, request.pair, &source, &resolved)
        .err()
        .expect("unsupported member combination is refused");
    assert!(failure.message.contains(diagnostic), "{}", failure.message);
    true
}

/// An out-of-phase pose is either a certified overlap or lacks a positive
/// clearance certificate; both refuse the pose.
fn refused_pose(failure: &ModelError) {
    assert!(
        failure.message.starts_with("Rust gear pair ")
            && (failure.message.contains("overlaps")
                || failure.message.contains("certified positive clearance")),
        "{}",
        failure.message
    );
}
