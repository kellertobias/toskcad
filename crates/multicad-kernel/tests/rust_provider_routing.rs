use multicad_kernel::{BooleanOperation, GeometryErrorKind, GeometryPolicy, PlaneFrame, Solid};

#[test]
fn oblique_cylinder_refuses_unsupported_geometry_and_preserves_cardinal_cut_retry() {
    let policy = GeometryPolicy::default();
    let stock = Solid::box_solid([10.; 3], policy).unwrap();
    let before = stock.bounds().unwrap();
    let s = 2_f64.sqrt().recip();
    let t = 6_f64.sqrt().recip();
    let oblique = Solid::cylinder(
        PlaneFrame {
            origin: [0.; 3],
            u: [s, -s, 0.],
            v: [t, t, -2. * t],
        },
        0.1,
        20.,
        policy,
    )
    .unwrap();
    let failure = stock
        .boolean(&oblique, BooleanOperation::Subtract, policy)
        .err()
        .unwrap();
    assert_eq!(failure.kind, GeometryErrorKind::Unsupported);
    assert!(failure.message.contains("Oblique cylindrical subtraction"));
    assert_eq!(stock.bounds().unwrap(), before);
    assert_eq!(stock.volume_mm3().unwrap(), 1000.);
    let tangent = Solid::cylinder(
        PlaneFrame {
            origin: [0., 5., -1.],
            ..Default::default()
        },
        0.1,
        12.,
        policy,
    )
    .unwrap();
    let failure = stock
        .boolean(&tangent, BooleanOperation::Subtract, policy)
        .err()
        .unwrap();
    assert_eq!(failure.kind, GeometryErrorKind::Input);
    let cardinal = Solid::cylinder(
        PlaneFrame {
            origin: [5., 5., -1.],
            ..Default::default()
        },
        0.1,
        12.,
        policy,
    )
    .unwrap();
    let cut = stock
        .boolean(&cardinal, BooleanOperation::Subtract, policy)
        .unwrap()
        .unwrap();
    assert!((cut.volume_mm3().unwrap() - (1000. - std::f64::consts::PI * 0.01 * 10.)).abs() < 1e-9);
    assert_eq!(stock.bounds().unwrap(), before);
    assert_eq!(stock.volume_mm3().unwrap(), 1000.);
}
