use super::*;
use crate::commands::Session;
#[test]
fn actual_two_start_composition_projection_uses_advance_per_revolution() {
    let session =
        Session::new(Document::new("Actual multi-start connection phase").unwrap()).unwrap();
    for hand in [ThreadHand::Right, ThreadHand::Left] {
        let mut parameters = default_machine_screw();
        let MechanicalPartParameters::MachineScrew { thread, .. } = &mut parameters else {
            panic!()
        };
        thread.family = ThreadFamily::Custom;
        thread.fit = ThreadFit::NominalBasic;
        thread.metric_fit_engagement = None;
        thread.starts = 2;
        thread.pitch = length(0.5);
        thread.custom = Some(CustomThreadProfile {
            depth: length(0.2),
            crest_width: length(0.05),
            lower_flank: angle(30.),
            upper_flank: angle(30.),
            ..Default::default()
        });
        thread.hand = hand;
        let source = parameters.resolve(session.parameters()).unwrap();
        let ResolvedMechanicalPart::MachineScrew { thread, .. } = &source else {
            panic!()
        };
        let (start, advance, resolved_hand) = composition_thread_phase(&source, -0.2, 4.1).unwrap();
        assert_eq!(advance, thread.lead_mm);
        assert_eq!(advance, 1.);
        assert_eq!(resolved_hand, hand);
        let station = start + 0.3;
        let sign = if hand == ThreadHand::Right { 1. } else { -1. };
        let roll = sign * std::f64::consts::TAU * (station - start) / advance;
        assert!((roll - sign * std::f64::consts::TAU * 0.3).abs() < 1e-12);
    }
}

#[test]
fn resolved_hanger_single_start_and_nonmetric_hardware_keep_their_phase_domains() {
    let session = Session::new(Document::new("Distinct source phase domains").unwrap()).unwrap();
    let hanger = default_custom_connection_hanger_bolt()
        .resolve(session.parameters())
        .unwrap();
    let ResolvedMechanicalPart::CustomConnection { parameters } = &hanger else {
        panic!()
    };
    let ResolvedCustomConnectionHardware::HangerBolt {
        metric_thread,
        plain_length_mm,
        ..
    } = parameters.as_ref()
    else {
        panic!()
    };
    let (start, lead, hand) = composition_thread_phase(&hanger, 0., 4.1).unwrap();
    assert_eq!(start, 4.1 - plain_length_mm - metric_thread.length_mm);
    assert_eq!(lead, metric_thread.lead_mm);
    assert_eq!(hand, metric_thread.hand);
    for parameters in [
        default_custom_connection_wood_screw(),
        default_custom_connection_wood_insert(),
    ] {
        assert!(
            composition_thread_phase(&parameters.resolve(session.parameters()).unwrap(), 0., 4.1)
                .is_none()
        );
    }
    let washer = MechanicalPartParameters::Washer {
        bore_diameter: length(3.4),
        outside_diameter: length(7.),
        thickness: length(0.2),
        frame: Default::default(),
    }
    .resolve(session.parameters())
    .unwrap();
    assert!(composition_thread_phase(&washer, 0., 4.1).is_none());
}
