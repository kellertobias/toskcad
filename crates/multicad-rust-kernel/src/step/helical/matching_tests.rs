use super::*;
fn source(rounded: bool, left: bool) -> Thread {
    let root = rounded.then_some(1. / (4. * 3f64.sqrt()));
    Thread::new(crate::helical_thread::Spec {
        major_radius_mm: 3.,
        pitch_mm: 1.,
        length_mm: 3.7,
        depth_mm: root.map_or(5. * 3f64.sqrt() / 16., |r| 0.875 * 3f64.sqrt() / 2. - r),
        crest_width_mm: 0.125,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts: 1,
        left_hand: left,
        root_radius_mm: root,
    })
    .unwrap()
}
fn grid(value: &ThreadSplines) -> RetainedReferenceGrid {
    RetainedReferenceGrid {
        radial_poles: std::array::from_fn(|i| value.bands[i].poles.len()),
        axial_poles: value.axial.poles.len(),
        axial_entry: value.visible_start,
    }
}
#[test]
fn matching_reference_keeps_identical_real_factors_both_roots_hands_and_trim_prefixes() {
    for rounded in [false, true] {
        for left in [false, true] {
            for visible in [0., 0.37] {
                let thread = source(rounded, left);
                let original = prepare_retained_reference(&thread, visible, 1e-7).unwrap();
                let matched =
                    prepare_retained_reference_matching(&thread, visible, 1e-7, grid(&original))
                        .unwrap();
                assert_eq!(matched.visible_start, original.visible_start);
                assert_eq!(matched.heights, original.heights);
                assert_eq!(matched.axial.knots, original.axial.knots);
                assert_eq!(matched.axial.poles, original.axial.poles);
                for (a, b) in matched.bands.iter().zip(&original.bands) {
                    assert_eq!(a.knots, b.knots);
                    assert_eq!(a.poles, b.poles);
                }
            }
        }
    }
}
#[test]
fn actual_reference_dimension_entry_and_invalid_domain_mismatches_refuse() {
    for rounded in [false, true] {
        for visible in [0., 0.37] {
            let thread = source(rounded, false);
            let reference = prepare_retained_reference(&thread, visible, 1e-7).unwrap();
            let good = grid(&reference);
            let mut radial = good;
            radial.radial_poles[0] += 11;
            let mut zero = good;
            zero.radial_poles[0] = 0;
            for wrong in [
                radial,
                zero,
                RetainedReferenceGrid {
                    axial_poles: good.axial_poles + 11,
                    ..good
                },
                RetainedReferenceGrid {
                    axial_entry: good.axial_entry + 11,
                    ..good
                },
                RetainedReferenceGrid {
                    axial_poles: usize::MAX,
                    ..good
                },
                RetainedReferenceGrid {
                    axial_entry: usize::MAX,
                    ..good
                },
            ] {
                assert!(matches!(
                    prepare_retained_reference_matching(&thread, visible, 1e-7, wrong),
                    Err(StepError::Unsupported)
                ));
            }
            for (v, t) in [
                (-1., 1e-7),
                (thread.spec().length_mm, 1e-7),
                (visible, 0.),
                (f64::NAN, 1e-7),
            ] {
                assert!(matches!(
                    prepare_retained_reference_matching(&thread, v, t, good),
                    Err(StepError::InvalidGeometry)
                ));
            }
            let mut spec = thread.spec();
            spec.starts = 2;
            let different = Thread::new(spec).unwrap();
            assert!(matches!(
                prepare_retained_reference_matching(&different, visible, 1e-7, good),
                Err(StepError::Unsupported)
            ));
        }
    }
}
#[test]
fn actual_circle_factor_checks_count_before_any_jet_generation() {
    let calls = std::cell::Cell::new(0usize);
    let result = Spline::with_jets_matching(0., 1., 0.25, Some(12), |t| {
        calls.set(calls.get() + 1);
        affine_jet(t, 1., 0., 1.)
    });
    assert!(matches!(result, Err(StepError::Unsupported)));
    assert_eq!(calls.get(), 0);
    let matched = Spline::with_jets_matching(0., 1., 0.25, Some(45), |t| {
        calls.set(calls.get() + 1);
        affine_jet(t, 1., 0., 1.)
    })
    .unwrap();
    assert_eq!(matched.poles.len(), 45);
    assert!(calls.get() > 0);
}

#[test]
fn actual_rounded_count_mismatch_and_excess_accuracy_refuse_without_source_mutation() {
    let thread = source(true, false);
    let before = thread;
    let reference = prepare_retained_reference(&thread, 0., 1e-7).unwrap();
    let mut wrong_root = grid(&reference);
    wrong_root.radial_poles[1] += 11;
    assert!(matches!(
        prepare_retained_reference_matching(&thread, 0., 1e-7, wrong_root),
        Err(StepError::Unsupported)
    ));
    // The complete actual reference exceeds the fixed1024-span bound at
    // this accuracy request; no caller quota is increased.
    assert!(matches!(
        prepare_retained_reference(&thread, 0., 1e-36),
        Err(StepError::Limit)
    ));
    assert_eq!(thread, before);
    assert!(prepare_retained_reference(&thread, 0., 1e-7).is_ok());
}
#[test]
fn existing_affine_factor_entry_retains_identical_unrestricted_coefficients() {
    let original = Spline::new(0., 1., [3., 2.5], -1., 1e-8).unwrap();
    let delegated = Spline::new_matching(0., 1., [3., 2.5], -1., 1e-8, None).unwrap();
    assert_eq!(original.knots, delegated.knots);
    assert_eq!(original.poles, delegated.poles);
}
