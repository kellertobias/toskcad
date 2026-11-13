use super::*;
use crate::{head_profile::TransportError, helical_thread};
#[test]
fn actual_authored_source_volumes_are_paired_with_cavity_and_profile_precision() {
    for index in 0..3 {
        for plain in [0., 0.2, 1.] {
            for left in [false, true] {
                let s = super::super::tests::fixture(index, plain, left);
                let bounds = s.source_volume_bounds().unwrap();
                let [lo, hi] = bounds.volume_mm3();
                assert!(lo > 0. && lo <= s.nominal_volume_mm3() && s.nominal_volume_mm3() <= hi);
                assert!(hi - lo < 1e-4);
                assert!(bounds.profile_precision_error_mm3() > 0.);
                assert!(bounds.cavity_source_width_mm3() <= 1e-5);
            }
        }
    }
}
#[test]
fn actual_fitted_rounded_roots_use_the_original_literal_circle_and_interval() {
    let s = super::super::tests::fixture(2, 0., false);
    for left in [false, true] {
        let thread = Thread::new(helical_thread::Spec {
            major_radius_mm: 1.4635,
            pitch_mm: 0.5,
            length_mm: 8.,
            depth_mm: 0.29908756729740676,
            crest_width_mm: 0.07131009168581062,
            lower_flank_deg: 30.,
            upper_flank_deg: 30.,
            starts: 1,
            left_hand: left,
            root_radius_mm: Some(0.07216878364870323),
        })
        .unwrap();
        let original = thread;
        let s = PhillipsScrew::new(
            thread,
            *s.profile(),
            s.spec(),
            s.cavity().clone(),
            100_000,
            &|| false,
        )
        .unwrap();
        let interval = thread_volume(s.thread(), s.visible_start_mm).unwrap();
        let n = 32768;
        let mut integral = 0.;
        for band in s.thread().profile_breaks().windows(2) {
            let [a, b] = [band[0], band[1]];
            let mut value =
                s.thread().radius_at_phase(a).powi(2) + s.thread().radius_at_phase(b).powi(2);
            let mut correction = 0.;
            for i in 1..n {
                let phase = a + (b - a) * i as f64 / n as f64;
                let term =
                    (if i % 2 == 0 { 2. } else { 4. }) * s.thread().radius_at_phase(phase).powi(2);
                let corrected = term - correction;
                let next = value + corrected;
                correction = (next - value) - corrected;
                value = next;
            }
            integral += value * (b - a) / (3. * n as f64);
        }
        let numeric = PI * (s.thread().spec().length_mm - s.visible_start_mm)
            / s.thread().spec().pitch_mm
            * integral;
        assert!(
            interval.lo <= numeric && numeric <= interval.hi,
            "{interval:?}/{numeric}"
        );
        let bounds = s.source_volume_bounds().unwrap().volume_mm3();
        assert!(bounds[0] <= s.nominal_volume_mm3() && s.nominal_volume_mm3() <= bounds[1]);
        assert_eq!(s.thread(), &original);
        assert!(root_primitive(s.thread().root_circle().unwrap(), 100.).is_none());
        assert!(thread_volume(s.thread(), 9.).is_none());
    }
}
#[test]
fn unresolved_transported_head_precision_cannot_be_hidden_by_nominal_mass() {
    let s = super::super::tests::fixture(1, 1., false);
    let segments: Vec<_> = s.profile().segments().collect();
    let p = Profile::from_transported_segments(
        &segments,
        TransportError {
            numeric_conversion_mm: 0.005,
            arithmetic_mm: 0.,
        },
    )
    .unwrap();
    let s = PhillipsScrew::new(
        *s.thread(),
        p,
        SourceSpec {
            head: None,
            ..s.spec()
        },
        s.cavity().clone(),
        100_000,
        &|| false,
    )
    .unwrap();
    assert_eq!(s.source_volume_bounds(), Err(Error::Unresolved));
}
#[test]
fn individually_positive_component_bounds_can_exhaust_total_material_precision() {
    use crate::head_profile::Head;
    let s = super::super::tests::fixture(1, 0., false);
    let p = Profile::new(
        Head::Cylindrical {
            radius_mm: 2.8,
            height_mm: 2.4,
        },
        1.5,
        0.,
    )
    .unwrap();
    let segments: Vec<_> = p.segments().collect();
    let profile = |error| {
        Profile::from_transported_segments(
            &segments,
            TransportError {
                numeric_conversion_mm: error,
                arithmetic_mm: 0.,
            },
        )
        .unwrap()
    };
    // A conservative declared transport reserve can leave the head's own
    // lower volume positive while consuming the body-minus-cavity budget.
    // Find that real precision boundary, without changing any literal support.
    let (mut low, mut high) = (1e-6, 0.05);
    for _ in 0..64 {
        let mid = low + (high - low) / 2.;
        let enough = crate::head_profile::source_volume::volume_bounds(&profile(mid))
            .is_some_and(|v| v.lo > 0.5);
        if enough {
            low = mid;
        } else {
            high = mid;
        }
    }
    let p = profile(low);
    assert!(
        crate::head_profile::source_volume::volume_bounds(&p)
            .unwrap()
            .lo
            > 0.
    );
    let mut spec = s.thread().spec();
    spec.length_mm = 0.5;
    let t = Thread::new(spec).unwrap();
    let s = PhillipsScrew::new(
        t,
        p,
        SourceSpec {
            head: None,
            ..s.spec()
        },
        s.cavity().clone(),
        100_000,
        &|| false,
    )
    .unwrap();
    assert!(s.nominal_volume_mm3() > 0.);
    assert_eq!(s.source_volume_bounds(), Err(Error::Unresolved));
}
