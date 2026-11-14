use super::*;
use crate::display::{Tessellation, axial_profile_stack};
use crate::manufacturing::certify;
fn spec() -> RectangularBlindCounterbore {
    RectangularBlindCounterbore {
        min_mm: [-10.; 2],
        max_mm: [10.; 2],
        height_mm: 20.,
        bore_radius_mm: 2.,
        bore_depth_mm: 16.,
        seat_radius_mm: 6.,
        seat_depth_mm: 1.,
        minimum_wall_mm: 1.,
        minimum_floor_mm: 1.,
    }
}
#[test]
fn blind_counterbore_retains_real_seat_shoulder_pilot_and_floor() {
    let s = spec();
    let body = ProfileStack::blind_counterbore(s, 1e-7, 1e-5).unwrap();
    assert_eq!(
        (
            body.vertices().len(),
            body.edges().len(),
            body.faces().len()
        ),
        (24, 44, 22)
    );
    let mut uses = vec![[0; 2]; body.edges().len()];
    for face in body.faces() {
        for edge in &face.boundary {
            uses[edge.edge][usize::from(edge.forward)] += 1;
        }
    }
    assert!(uses.iter().all(|u| *u == [1, 1]));
    let expected = 8000. - PI * (4. * 16. + 32.);
    assert!((body.volume_mm3() - expected).abs() < 1e-9);
    let mesh = axial_profile_stack(
        &body,
        Tessellation {
            deflection: 0.001,
            angle: 0.1,
            triangle_limit: 20000,
        },
    )
    .unwrap();
    let proof = certify(&mesh).unwrap();
    assert!((proof.signed_volume_mm3 - expected).abs() < 1.);
    assert_eq!(
        body.faces()
            .iter()
            .filter(|f| matches!(f.support, Support::InwardCylinder { .. }))
            .count(),
        8
    );
    let area: f64 = body.faces().iter().map(|f| f.area_mm2).sum();
    let expected_area = 2400. - PI * 36. + PI * 32. + PI * 4. + 2. * PI * (6. + 2. * 15.);
    assert!(
        (area - expected_area).abs() < 1e-8,
        "{area} != {expected_area}"
    );
}
#[test]
fn blind_counterbore_asymmetric_stock_preserves_entry_area_and_centroid() {
    let mut s = spec();
    s.min_mm = [-8., -9.];
    s.max_mm = [12., 15.];
    let body = ProfileStack::blind_counterbore(s, 1e-7, 1e-5).unwrap();
    let entry: Vec<_> = body.faces().iter().filter(|f| f.sources == [5]).collect();
    let area: f64 = entry.iter().map(|f| f.area_mm2).sum();
    assert!((area - (480. - 36. * PI)).abs() < 1e-9);
    for (axis, expected) in [(0, 960.), (1, 1440.)] {
        let moment: f64 = entry.iter().map(|f| f.area_mm2 * f.center_mm[axis]).sum();
        assert!((moment - expected).abs() < 1e-9);
    }
    let mesh = axial_profile_stack(
        &body,
        Tessellation {
            deflection: 0.001,
            angle: 0.1,
            triangle_limit: 20000,
        },
    )
    .unwrap();
    assert!(certify(&mesh).is_ok());
}
#[test]
fn blind_counterbore_rejects_invalid_or_unprotected_source_stations() {
    for s in [
        RectangularBlindCounterbore {
            height_mm: f64::NAN,
            ..spec()
        },
        RectangularBlindCounterbore {
            bore_radius_mm: 0.,
            ..spec()
        },
        RectangularBlindCounterbore {
            seat_radius_mm: 1.,
            ..spec()
        },
        RectangularBlindCounterbore {
            seat_depth_mm: 16.,
            ..spec()
        },
        RectangularBlindCounterbore {
            min_mm: [10.; 2],
            ..spec()
        },
    ] {
        assert!(matches!(
            ProfileStack::blind_counterbore(s, 1e-7, 1e-5),
            Err(Error::Invalid)
        ));
    }
    for s in [
        RectangularBlindCounterbore {
            minimum_wall_mm: 5.,
            ..spec()
        },
        RectangularBlindCounterbore {
            minimum_floor_mm: 5.,
            ..spec()
        },
        RectangularBlindCounterbore {
            min_mm: [-6., -10.],
            ..spec()
        },
    ] {
        assert!(matches!(
            ProfileStack::blind_counterbore(s, 1e-7, 1e-5),
            Err(Error::Unsupported)
        ));
    }
    for (t, m) in [
        (-1., 1e-5),
        (f64::NAN, 1e-5),
        (1e-5, 1e-5),
        (0., f64::INFINITY),
    ] {
        assert!(matches!(
            ProfileStack::blind_counterbore(spec(), t, m),
            Err(Error::Invalid)
        ));
    }
    assert!(sum_leq(6., 4., 10.));
    assert!(!sum_leq(1., 2f64.powi(-54), 1.));
    assert!(!sum_leq(6., 5., 10.));
    let boundary = RectangularBlindCounterbore {
        minimum_wall_mm: 4.,
        minimum_floor_mm: 4.,
        ..spec()
    };
    assert!(ProfileStack::blind_counterbore(boundary, 1e-7, 1e-5).is_ok());
}
#[test]
fn blind_counterbore_step_is_independently_decoded_as_one_material_solid() {
    let body = ProfileStack::blind_counterbore(spec(), 1e-7, 1e-5).unwrap();
    let text =
        crate::step::write_axial_profile_stack(&body, [0.; 3], [1., 0., 0.], [0., 1., 0.]).unwrap();
    let read = crate::step_read::read_analytic_solid(&text, 1e-7).unwrap();
    assert_eq!(read.solid_count, 1);
    assert_eq!(read.face_count, 22);
    assert!((read.volume_mm3 - body.volume_mm3()).abs() < 1e-7);
}
