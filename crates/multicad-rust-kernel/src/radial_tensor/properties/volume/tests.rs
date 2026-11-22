use super::*;
fn material(input: BoundaryInput) -> RadialBoundary {
    BoundaryOwnership::new(input, 1e-7, None)
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
fn cam() -> RadialBoundary {
    material(crate::radial_tensor::tests::cam::wedge(&[1, 2]))
}
#[test]
fn volume_only_encloses_the_same_actual_finite_source_as_full_properties() {
    let body = cam();
    let full = body.properties_with_interrupt(1e-5, None).unwrap();
    let volume = body.volume_with_interrupt(1e-5, None).unwrap();
    assert!(volume.volume_mm3[0] <= full.volume_mm3[1]);
    assert!(volume.volume_mm3[1] >= full.volume_mm3[0]);
    assert_eq!(
        volume.characteristic_length_mm,
        full.characteristic_length_mm
    );
    assert!(
        volume.volume_mm3[1] - volume.volume_mm3[0]
            <= 1e-5 * volume.characteristic_length_mm.powi(3)
    );
    let before = body.input().clone();
    for accuracy in [0., -1., 1., f64::NAN, f64::INFINITY] {
        assert_eq!(
            body.volume_with_interrupt(accuracy, None).err(),
            Some(Error::Invalid)
        );
    }
    assert_eq!(
        body.volume_with_interrupt(1e-5, Some(&|| true)).err(),
        Some(Error::Cancelled)
    );
    assert_eq!(
        body.volume_with_budget(1e-5, None, 0).err(),
        Some(Error::Limit)
    );
    assert_eq!(
        body.volume_with_interrupt(1e-30, None).err(),
        Some(Error::Unresolved)
    );
    let count = std::cell::Cell::new(0usize);
    let counting = || {
        count.set(count.get() + 1);
        false
    };
    let measured = body.volume_with_interrupt(1e-5, Some(&counting)).unwrap();
    assert_eq!(measured, volume);
    let end = count.get();
    for stop in [2, end - 1, end] {
        count.set(0);
        let cancel = || {
            count.set(count.get() + 1);
            count.get() >= stop
        };
        assert_eq!(
            body.volume_with_interrupt(1e-5, Some(&cancel)).err(),
            Some(Error::Cancelled)
        );
    }
    assert_eq!(body.input(), &before);
    assert_eq!(body.volume_with_interrupt(1e-5, None).unwrap(), volume);
}
#[test]
fn authored_npt_volume_query_keeps_analytic_volume_and_avoids_area_integration() {
    use crate::axial_pipe_thread::authored::{AuthoredSection, AuthoredSpec, AuthoredThread};
    let p = 25.4 / 27.;
    let k = 1. / 32.;
    let m = 3_f64.sqrt();
    let slopes = 1. / (m - k) + 1. / (m + k);
    let d = p / slopes - 2. * 0.0645 * p;
    let w = (p - d * slopes) / 2.;
    let spec = AuthoredSpec {
        major_radius_mm: (0.37360 * 25.4 - 0.1615 * 25.4 / 16. + d) / 2.,
        pitch_mm: p,
        length_mm: 0.2639 * 25.4,
        radial_taper: k,
        section: AuthoredSection::NptFlat {
            depth_mm: d,
            flat_width_mm: w,
        },
        crest_radius_limit_mm: None,
    };
    let raw = AuthoredThread::new(spec, 1e-7, 1e-6)
        .unwrap()
        .raw_boundary_ownership(1e-7, None)
        .unwrap();
    let body = raw
        .ownership
        .bind_surfaces(1e-7, None)
        .unwrap()
        .certify_axial(1e-7, 1e-6, None)
        .unwrap()
        .certify_caps(None)
        .unwrap()
        .certify_ordering(None)
        .unwrap()
        .certify_material(1e-7, None)
        .unwrap();
    let full_started = std::time::Instant::now();
    let full = body.properties_with_interrupt(1e-5, None).unwrap();
    let full_elapsed = full_started.elapsed();
    let started = std::time::Instant::now();
    let volume = body.volume_with_interrupt(1e-5, None).unwrap();
    let elapsed = started.elapsed();
    let root_width = p - w - d * slopes;
    let area = |z: f64| {
        let r = spec.major_radius_mm + k * z;
        let root = r - d;
        std::f64::consts::PI / p
            * (w * r * r + root_width * root * root + slopes * (r.powi(3) - root.powi(3)) / 3.)
    };
    let half = spec.length_mm / 2.;
    let expected = half * (area(half - half / 3_f64.sqrt()) + area(half + half / 3_f64.sqrt()));
    assert!(volume.volume_mm3[0] <= expected && expected <= volume.volume_mm3[1]);
    assert!(
        volume.volume_mm3[0] <= full.volume_mm3[1] && full.volume_mm3[0] <= volume.volume_mm3[1]
    );
    eprintln!(
        "actual NPT full mass={full_elapsed:?}, volume only={elapsed:?}, volume={:?}, independent={expected}",
        volume.volume_mm3
    );
}
