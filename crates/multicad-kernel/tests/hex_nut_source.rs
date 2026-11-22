use multicad_core::document::EntityId;
use multicad_kernel::{
    GeometryPolicy, HelicalThread, PlaneFrame, RigidTransform, Solid, rust_kernel::Rational,
};
fn spec(left_hand: bool) -> HelicalThread {
    HelicalThread {
        major_radius_mm: 1.55,
        pitch_mm: 0.5,
        length_mm: 2.4,
        depth_mm: 0.28,
        crest_width_mm: 0.08,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        radial_taper: 0.,
        starts: 1,
        left_hand,
    }
}
#[test]
fn retained_nut_getter_owns_actual_hex_material_and_current_proper_pose() {
    let policy = GeometryPolicy::default();
    let frame = PlaneFrame {
        origin: [5., -3., 2.],
        u: [0.8, 0.6, 0.],
        v: [-0.36, 0.48, 0.8],
    };
    let pose = RigidTransform {
        rotation: [[0., 0., 1.], [1., 0., 0.], [0., 1., 0.]],
        translation_mm: [100., -30., 50.],
    };
    for left in [false, true] {
        let solid =
            Solid::rust_hex_nut(spec(left), 5.5, 0.2, EntityId::new(), frame, policy).unwrap();
        let (original, source_frame) = solid.rust_hex_nut_geometry().unwrap();
        let original = *original;
        assert_eq!(source_frame, frame);
        assert_eq!(original.across_flats_mm(), 5.5);
        assert_eq!(original.thread().spec().left_hand, left);
        let moved = solid.transformed(pose).unwrap();
        let (actual, placed) = moved.rust_hex_nut_geometry().unwrap();
        assert_eq!(*actual, original);
        for i in 0..3 {
            let expected = pose.translation_mm[i]
                + pose.rotation[i]
                    .into_iter()
                    .zip(frame.origin)
                    .map(|(r, p)| r * p)
                    .sum::<f64>();
            assert!((placed.origin[i] - expected).abs() < 1e-12);
        }
        assert_eq!(solid.rust_hex_nut_geometry().unwrap(), (&original, frame));
        let bare =
            Solid::rust_helical_thread(spec(left), None, EntityId::new(), frame, policy).unwrap();
        assert!(bare.rust_hex_nut_geometry().is_none());
    }
    let stock = Solid::rust_box(
        EntityId::new().uuid_bits(),
        [Rational::new(20, 1).unwrap(); 3],
        PlaneFrame::default(),
        policy,
    )
    .unwrap();
    assert!(stock.rust_hex_nut_geometry().is_none());
}
