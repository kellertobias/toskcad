use super::*;
use multicad_rust_kernel::{
    authored_flat_thread::{AuthoredThread, SourceOptions},
    authored_hanger_bolt::{Hanger, Spec as HangerSpec},
    authored_metric_thread::{AuthoredMetricThread, Spec as MetricSpec},
};
fn fixture(
    left: bool,
    nonbinary: bool,
    policy: GeometryPolicy,
) -> (Solid, Solid, Solid, Requirements) {
    let (stock, _, mut req) = super::super::tests::fixture();
    let major = (2.874 + 2.980) / 4.;
    let pitch_mid = (2.580 + 2.655) / 2.;
    let source = Hanger::new(
        HangerSpec {
            wood: AuthoredThread::new(
                multicad_rust_kernel::helical_thread::Spec {
                    major_radius_mm: 3.,
                    pitch_mm: 2.,
                    length_mm: if nonbinary { 12.3 } else { 12. },
                    depth_mm: 1.,
                    crest_width_mm: 0.4,
                    lower_flank_deg: 30.,
                    upper_flank_deg: 30.,
                    starts: 1,
                    left_hand: left,
                    root_radius_mm: None,
                },
                SourceOptions {
                    position_limit_mm: 1e-7,
                    normal_angle_limit_rad: 1e-7,
                },
            )
            .unwrap(),
            metric: AuthoredMetricThread::new(MetricSpec {
                major_radius_mm: major,
                pitch_mm: 0.5,
                crest_width_mm: 0.25 - (2. * major - pitch_mid) / 3_f64.sqrt(),
                length_mm: if nonbinary { 12.7 } else { 12. },
                starts: 1,
                left_hand: left,
            })
            .unwrap(),
            plain_length_mm: if nonbinary { 2.15 } else { 2. },
            plain_radius_mm: 1.5,
        },
        None,
    )
    .unwrap();
    let body_frame = PlaneFrame {
        origin: [10., 10., 20.],
        u: [1., 0., 0.],
        v: [0., -1., 0.],
    };
    let hanger = Solid::rust_hanger(source, EntityId::new(), body_frame, policy, None).unwrap();
    let station = -6.7;
    let start = -source.spec().plain_length_mm - source.metric().full_length_mm();
    let angle = std::f64::consts::TAU * (station - start) / 0.5 * if left { -1. } else { 1. };
    let (s, c) = angle.sin_cos();
    let nut = Solid::rust_hex_nut(
        crate::HelicalThread {
            major_radius_mm: 1.5,
            pitch_mm: 0.5,
            length_mm: 2.4,
            depth_mm: (3. - (2.459 + 2.599) / 2.) / 2.,
            crest_width_mm: 0.25 - (3. - (2.675 + 2.775) / 2.) / 3_f64.sqrt(),
            lower_flank_deg: 30.,
            upper_flank_deg: 30.,
            radial_taper: 0.,
            starts: 1,
            left_hand: left,
        },
        5.5,
        0.2,
        EntityId::new(),
        PlaneFrame {
            origin: [10., 10., 20. - station],
            u: [c, -s, 0.],
            v: [-s, -c, 0.],
        },
        policy,
    )
    .unwrap();
    req.linear_tolerance_mm = policy.linear_tolerance_mm;
    req.minimum_feature_mm = policy.minimum_feature_mm;
    (stock, hanger, nut, req)
}
fn mating(policy: GeometryPolicy) -> HangerMatingRequirements {
    HangerMatingRequirements {
        minimum_engagement_mm: 2.,
        linear_tolerance_mm: policy.linear_tolerance_mm,
        ..Default::default()
    }
}
#[test]
fn actual_native_hanger_pilot_and_nonquarter_retained_nut_certify_both_hands_and_placed_sources() {
    for (left, nb, posed) in [
        (false, false, false),
        (true, true, false),
        (false, true, true),
        (true, false, true),
    ] {
        let policy = GeometryPolicy::default();
        let (mut stock, mut hanger, mut nut, req) = fixture(left, nb, policy);
        if posed {
            let pose = RigidTransform {
                rotation: [[0.8, -0.36, 0.48], [0.6, 0.48, -0.64], [0., 0.8, 0.6]],
                translation_mm: [11., -7., 5.],
            };
            stock = stock.transformed(pose).unwrap();
            hanger = hanger.transformed(pose).unwrap();
            nut = nut.transformed(pose).unwrap();
        }
        let retained = hanger.rust_hanger_data().unwrap();
        let wood = stock.rust_hanger_wood_contact(&hanger, req, None).unwrap();
        let threads = hanger
            .rust_hanger_nut_mating(&nut, mating(policy), None)
            .unwrap();
        assert!(wood.engagement_mm[0] >= 8. && wood.witness_volume_lower_mm3 > 0.);
        assert!(
            threads.engagement_mm[0] >= 2. && threads.minimum_radial_clearance_lower_mm >= 1e-7
        );
        assert_eq!(hanger.rust_hanger_data().unwrap(), retained);
    }
}
#[test]
fn actual_source_policies_invalid_roles_offsets_and_cancellation_refuse_atomically() {
    let policy = GeometryPolicy::default();
    let (stock, hanger, nut, req) = fixture(false, false, policy);
    let retained = hanger.rust_hanger_data();
    let unknown = Solid::box_solid([2.; 3], policy).unwrap();
    for result in [
        unknown.rust_hanger_wood_contact(&hanger, req, None),
        stock.rust_hanger_wood_contact(&unknown, req, None),
        stock.rust_hanger_wood_contact(&nut, req, None),
    ] {
        assert_eq!(result.unwrap_err().kind, GeometryErrorKind::Unsupported);
    }
    assert_eq!(
        unknown
            .rust_hanger_nut_mating(&nut, mating(policy), None)
            .unwrap_err()
            .kind,
        GeometryErrorKind::Unsupported
    );
    assert_eq!(
        hanger
            .rust_hanger_nut_mating(&unknown, mating(policy), None)
            .unwrap_err()
            .kind,
        GeometryErrorKind::Unsupported
    );
    for result in [
        stock
            .rust_hanger_wood_contact(&hanger, req, Some(&|| true))
            .map(|_| ()),
        hanger
            .rust_hanger_nut_mating(&nut, mating(policy), Some(&|| true))
            .map(|_| ()),
    ] {
        assert_eq!(result.unwrap_err().kind, GeometryErrorKind::Cancelled);
    }
    let shifted = nut
        .transformed(RigidTransform {
            translation_mm: [0.001, 0., 0.],
            ..Default::default()
        })
        .unwrap();
    assert_eq!(
        hanger
            .rust_hanger_nut_mating(&shifted, mating(policy), None)
            .unwrap_err()
            .kind,
        GeometryErrorKind::Unsupported
    );
    let mut bad = req;
    bad.linear_tolerance_mm = 1e-6;
    assert_eq!(
        stock
            .rust_hanger_wood_contact(&hanger, bad, None)
            .unwrap_err()
            .kind,
        GeometryErrorKind::Input
    );
    let mut bad = mating(policy);
    bad.linear_tolerance_mm = 1e-6;
    assert_eq!(
        hanger
            .rust_hanger_nut_mating(&nut, bad, None)
            .unwrap_err()
            .kind,
        GeometryErrorKind::Input
    );
    let mut other_policy = policy;
    other_policy.linear_tolerance_mm = 1e-6;
    let (_, other_hanger, other_nut, _) = fixture(false, false, other_policy);
    assert_eq!(
        stock
            .rust_hanger_wood_contact(&other_hanger, req, None)
            .unwrap_err()
            .kind,
        GeometryErrorKind::Input
    );
    assert_eq!(
        hanger
            .rust_hanger_nut_mating(&other_nut, mating(policy), None)
            .unwrap_err()
            .kind,
        GeometryErrorKind::Input
    );
    assert_eq!(hanger.rust_hanger_data(), retained);
}
