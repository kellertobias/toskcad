use super::*;
use multicad_rust_kernel::{
    authored_flat_thread::{AuthoredThread, SourceOptions},
    helical_thread::Spec,
    wood_insert::Insert,
};
fn fixture() -> (Solid, Solid, Requirements) {
    let (stock, _, mut req) = super::super::tests::fixture();
    let source = |radius, pitch, depth, crest| {
        AuthoredThread::new(
            Spec {
                major_radius_mm: radius,
                pitch_mm: pitch,
                length_mm: 4.,
                depth_mm: depth,
                crest_width_mm: crest,
                lower_flank_deg: 30.,
                upper_flank_deg: 30.,
                starts: 1,
                left_hand: false,
                root_radius_mm: None,
            },
            SourceOptions {
                position_limit_mm: 1e-8,
                normal_angle_limit_rad: 1e-6,
            },
        )
        .unwrap()
    };
    let insert = Insert::new(
        source(3., 2., 1., 0.4),
        source(1., 0.5, 0.25, 0.05),
        0.5,
        None,
    )
    .unwrap();
    let insert = Solid::rust_wood_insert(
        insert,
        EntityId::new(),
        PlaneFrame {
            origin: [10., 10., 20.],
            u: [1., 0., 0.],
            v: [0., -1., 0.],
        },
        GeometryPolicy::default(),
        None,
    )
    .unwrap();
    req.minimum_engagement_mm = 2.;
    (stock, insert, req)
}
#[test]
fn actual_native_insert_and_drilled_source_prove_positive_contact() {
    let (stock, insert, req) = fixture();
    assert!(
        stock
            .rust_wood_insert_contact(&insert, req, None)
            .unwrap()
            .witness_volume_lower_mm3
            > 0.
    );
}
#[test]
fn actual_source_policy_and_foreign_geometry_refusals_are_atomic() {
    let (stock, insert, req) = fixture();
    let original = insert.rust_wood_insert_data().unwrap().source;
    assert_eq!(
        stock
            .rust_wood_insert_contact(&insert, req, Some(&|| true))
            .unwrap_err()
            .kind,
        GeometryErrorKind::Cancelled
    );
    let mut bad = req;
    bad.linear_tolerance_mm *= 2.;
    assert_eq!(
        stock
            .rust_wood_insert_contact(&insert, bad, None)
            .unwrap_err()
            .kind,
        GeometryErrorKind::Input
    );
    let mut bad = req;
    bad.minimum_wall_mm = 100.;
    assert!(stock.rust_wood_insert_contact(&insert, bad, None).is_err());
    let plain = Solid::box_solid([20.; 3], GeometryPolicy::default()).unwrap();
    assert_eq!(
        plain
            .rust_wood_insert_contact(&insert, req, None)
            .unwrap_err()
            .kind,
        GeometryErrorKind::Unsupported
    );
    assert_eq!(
        insert
            .rust_wood_insert_contact(&insert, req, None)
            .unwrap_err()
            .kind,
        GeometryErrorKind::Unsupported
    );
    assert_eq!(
        stock
            .rust_wood_insert_contact(&plain, req, None)
            .unwrap_err()
            .kind,
        GeometryErrorKind::Unsupported
    );
    assert_eq!(insert.rust_wood_insert_data().unwrap().source, original);
    assert!(stock.rust_wood_insert_contact(&insert, req, None).is_ok());
}
