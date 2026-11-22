//! A short centerline cannot certify the radial extent of a wide thread tool.
use multicad_core::document::EntityId;
use multicad_kernel::{
    BlindThreadCut, GeometryPolicy, HelicalThread, PlaneFrame, Solid, rust_kernel::Rational,
};
use multicad_precision::PortableMath as _;
#[test]
fn shallow_wide_blind_thread_cannot_project_a_proper_tilted_tool_to_the_stock_cap() {
    let policy = GeometryPolicy::default();
    let owner = EntityId::new();
    let stock = Solid::rust_box(
        owner.uuid_bits(),
        [
            Rational::new(10_000_000, 1).unwrap(),
            Rational::new(10_000_000, 1).unwrap(),
            Rational::new(20, 1).unwrap(),
        ],
        PlaneFrame::default(),
        policy,
    )
    .unwrap();
    let before = stock.topology_snapshot(owner, 0).unwrap();
    let tilt: f64 = 1e-6;
    let profile = HelicalThread {
        major_radius_mm: 1000.,
        pitch_mm: 0.01,
        length_mm: 0.01,
        depth_mm: 0.002,
        crest_width_mm: 0.001,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        radial_taper: 0.,
        starts: 1,
        left_hand: false,
    };
    assert!(profile.length_mm * tilt < policy.linear_tolerance_mm);
    assert!(profile.major_radius_mm * tilt > 1000. * policy.linear_tolerance_mm);
    let cut = BlindThreadCut {
        profile,
        entry: PlaneFrame {
            origin: [2_000_000., 2_000_000., 0.],
            u: [tilt.portable_cos(), 0., tilt.portable_sin()],
            v: [0., 1., 0.],
        },
        receiving_radius_mm: 1000.,
        minimum_wall_mm: 1.,
        minimum_floor_mm: policy.minimum_feature_mm,
    };
    let failure = stock
        .rust_cut_blind_internal_thread(cut, EntityId::new(), policy)
        .err()
        .expect("the full radial tool extent must reject a genuine tilt");
    assert!(
        failure.message.contains("requested tool placement"),
        "{failure}"
    );
    assert_eq!(stock.topology_snapshot(owner, 0).unwrap(), before);
}

#[test]
fn full_tool_frame_guard_keeps_exact_and_resolvable_near_cardinal_sources() {
    let policy = GeometryPolicy::default();
    let owner = EntityId::new();
    let stock = Solid::rust_box(
        owner.uuid_bits(),
        [
            Rational::new(10_000_000, 1).unwrap(),
            Rational::new(10_000_000, 1).unwrap(),
            Rational::new(20, 1).unwrap(),
        ],
        PlaneFrame::default(),
        policy,
    )
    .unwrap();
    let before = stock.topology_snapshot(owner, 0).unwrap();
    for tilt in [0_f64, 1e-12] {
        let profile = HelicalThread {
            major_radius_mm: 1000.,
            pitch_mm: 0.01,
            length_mm: 0.01,
            depth_mm: 0.002,
            crest_width_mm: 0.001,
            lower_flank_deg: 30.,
            upper_flank_deg: 30.,
            radial_taper: 0.,
            starts: 1,
            left_hand: false,
        };
        let cut = stock
            .rust_cut_blind_internal_thread(
                BlindThreadCut {
                    profile,
                    entry: PlaneFrame {
                        origin: [2_000_000., 2_000_000., 0.],
                        u: [tilt.portable_cos(), 0., tilt.portable_sin()],
                        v: [0., 1., 0.],
                    },
                    receiving_radius_mm: 1000.,
                    minimum_wall_mm: 1.,
                    minimum_floor_mm: policy.minimum_feature_mm,
                },
                EntityId::new(),
                policy,
            )
            .unwrap();
        let floor = cut
            .faces()
            .unwrap()
            .into_iter()
            .find(|face| face.index == 1)
            .unwrap();
        assert!((floor.center_mm[2] - profile.length_mm).abs() < policy.linear_tolerance_mm);
        assert_eq!(stock.topology_snapshot(owner, 0).unwrap(), before);
    }
}
