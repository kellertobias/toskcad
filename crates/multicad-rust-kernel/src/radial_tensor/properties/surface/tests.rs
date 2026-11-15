use super::*;
use crate::axial_pipe_thread::authored::{AuthoredSection, AuthoredSpec, AuthoredThread};
fn actual_bspt(external: bool) -> RadialBoundary {
    let p = 25.4 / 28.;
    let taper = if external { 1. / 32. } else { -1. / 32. };
    let raw = AuthoredThread::new(
        AuthoredSpec {
            major_radius_mm: 9.728 / 2. - 4. / 32.,
            pitch_mm: p,
            length_mm: 6.5,
            radial_taper: taper,
            section: AuthoredSection::BsptRounded {
                depth_mm: p / (3. * 27.5_f64.to_radians().tan()),
            },
            crest_radius_limit_mm: None,
        },
        1e-7,
        1e-6,
    )
    .unwrap()
    .raw_boundary_ownership(1e-7, None)
    .unwrap();
    raw.ownership
        .bind_surfaces(1e-7, None)
        .unwrap()
        .certify_axial(1e-7, 1e-6, None)
        .unwrap()
        .certify_caps(None)
        .unwrap()
        .certify_ordering(None)
        .unwrap()
        .certify_material(1e-7, None)
        .unwrap()
}
#[test]
fn actual_canonical_bspt_faces_fit_the_existing_work_cap_without_inertia() {
    for external in [true, false] {
        let body = actual_bspt(external);
        let original = body.input().clone();
        let now = std::time::Instant::now();
        let value = body.face_properties_with_interrupt(1e-4, None);
        eprintln!(
            "actual BSPT surface external{external} {:?}: {:?}",
            now.elapsed(),
            value.as_ref().map(|v| v.len())
        );
        let value = value.unwrap();
        let total = value.iter().fold(I::point(0.), |sum, face| {
            sum + I::new(face.area_mm2[0], face.area_mm2[1])
        });
        eprintln!(
            "actual BSPT summed area external{external}: {total:?}; width {}",
            total.width()
        );
        assert_eq!(value.len(), body.input().faces.len());
        assert!(value.iter().all(|f| {
            f.area_mm2[0] > 0.
                && f.area_mm2[1] >= f.area_mm2[0]
                && f.centroid_mm
                    .iter()
                    .all(|r| r[0].is_finite() && r[1].is_finite() && r[0] <= r[1])
        }));
        assert_eq!(body.input(), &original);
    }
}
fn cam() -> RadialBoundary {
    BoundaryOwnership::new(crate::radial_tensor::tests::cam::wedge(&[1]), 1e-7, None)
        .unwrap()
        .bind_surfaces(1e-7, None)
        .unwrap()
        .certify_axial(1e-7, 1e-6, None)
        .unwrap()
        .certify_caps(None)
        .unwrap()
        .certify_ordering(None)
        .unwrap()
        .certify_material(1e-7, None)
        .unwrap()
}
fn contains(range: [f64; 2], expected: f64) {
    assert!(
        range[0] <= expected && expected <= range[1],
        "{range:?} misses {expected}"
    );
}
#[test]
fn independent_cam_surface_areas_centroids_and_rigid_world_transport() {
    let body = cam();
    let phi = 4_f64.atan2(3.);
    let angle = 2. * (std::f64::consts::PI - phi);
    let area = 25. * (std::f64::consts::PI - phi) + 30.;
    let expected = [
        (area, [(106. / 3.) / area, 0., 0.]),
        (area, [(106. / 3.) / area, 0., 4.]),
        (4. * 13_f64.sqrt(), [4.5, -3., 2.]),
        (16., [6., 0., 2.]),
        (4. * 13_f64.sqrt(), [4.5, 3., 2.]),
        (20. * angle, [-8. / angle, 0., 2.]),
    ];
    let full = body.properties_with_interrupt(1e-5, None).unwrap();
    let actual = body.face_properties_with_interrupt(1e-5, None).unwrap();
    assert_eq!(actual.len(), expected.len());
    for ((face, (area, centroid)), original) in actual.iter().zip(expected).zip(&full.faces) {
        contains(face.area_mm2, area);
        assert_eq!(face.face_index, original.face_index);
        assert_eq!(face.planar_support, original.planar_support);
        assert!(
            face.area_mm2[0] <= original.area_mm2[1] && original.area_mm2[0] <= face.area_mm2[1]
        );
        for (range, value) in face.centroid_mm.into_iter().zip(centroid) {
            contains(range, value);
        }
    }
    let transformed = body
        .transformed_with_interrupt(
            [[20., -10., 30.], [0., 1., 0.], [0., 0., 1.]],
            1e-7,
            1e-6,
            None,
        )
        .unwrap();
    let posed = transformed
        .face_properties_with_interrupt(1e-5, None)
        .unwrap();
    for (face, (area, [x, y, z])) in posed.iter().zip(expected) {
        contains(face.area_mm2, area);
        for (range, value) in face
            .centroid_mm
            .into_iter()
            .zip([20. + z, -10. + x, 30. + y])
        {
            contains(range, value);
        }
    }
}
#[test]
fn face_queries_refuse_invalid_tight_work_and_cancellation_then_retry() {
    let body = cam();
    let before = body.input().clone();
    for accuracy in [0., -1., 1., f64::NAN, f64::INFINITY] {
        assert_eq!(
            body.face_properties_with_interrupt(accuracy, None).err(),
            Some(Error::Invalid)
        );
    }
    for budget in [0, 1, 1000] {
        assert_eq!(
            body.face_properties_with_budget(1e-5, None, budget).err(),
            Some(Error::Limit)
        );
    }
    assert_eq!(
        body.face_properties_with_interrupt(1e-30, None).err(),
        Some(Error::Unresolved)
    );
    let count = std::cell::Cell::new(0usize);
    let counting = || {
        count.set(count.get() + 1);
        false
    };
    let expected = body
        .face_properties_with_interrupt(1e-5, Some(&counting))
        .unwrap();
    let end = count.get();
    for stop in [1, 2, end / 2, end - 1, end] {
        count.set(0);
        let cancel = || {
            count.set(count.get() + 1);
            count.get() >= stop
        };
        assert_eq!(
            body.face_properties_with_interrupt(1e-5, Some(&cancel))
                .err(),
            Some(Error::Cancelled)
        );
    }
    assert_eq!(body.input(), &before);
    assert_eq!(
        body.face_properties_with_interrupt(1e-5, None).unwrap(),
        expected
    );
}
