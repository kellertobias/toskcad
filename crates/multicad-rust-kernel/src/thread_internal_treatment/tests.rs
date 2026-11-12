use super::*;
use crate::helical_thread::Spec;
fn tool(left: bool, starts: u32) -> Thread {
    Thread::new(Spec {
        major_radius_mm: 1.5,
        pitch_mm: 0.5,
        length_mm: 3.3,
        depth_mm: 0.2355,
        crest_width_mm: 0.5 / 2. - (3. - 2.725) / 3_f64.sqrt(),
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts,
        left_hand: left,
        root_radius_mm: None,
    })
    .unwrap()
}
fn both() -> Intent {
    Intent {
        entry_depth_mm: Some(0.4),
        terminal: Some(Terminal {
            width_mm: 0.3,
            radial_depth_mm: 0.05,
        }),
    }
}
#[test]
fn actual_internal_max_source_keeps_phase_hands_starts_endpoints_and_largest_envelope() {
    for left in [false, true] {
        for starts in [1, 2, 8] {
            let original = tool(left, starts);
            let source = Cavity::new(original, both()).unwrap();
            assert_eq!(source.source_thread(), &original);
            assert_eq!(source.intent(), both());
            assert_eq!(source.full_length_mm(), 3.3);
            assert_eq!(
                source.entry_radius_mm(),
                Some(original.root_radius_mm() + 0.4)
            );
            assert_eq!(
                source.entry_cone_radius_mm(0.4),
                Some(original.root_radius_mm())
            );
            assert_eq!(source.terminal_radius_mm(), Some(1.55));
            assert_eq!(
                source.maximum_radius_bound_mm(),
                source.entry_radius_mm().unwrap()
            );
            let tail = source.terminal_station_mm().unwrap();
            assert_eq!(tail, 3.3 - 0.3);
            for theta in [0., 0.43, 2.1] {
                for z in [0., 0.01, 0.2, 0.4, 0.71, 2.5, tail, 3.3] {
                    let radius = source.radius_mm(theta, z);
                    if z >= tail {
                        assert_eq!(radius, 1.55);
                    } else if z <= 0.4 {
                        assert_eq!(
                            radius,
                            original
                                .radius(theta, z)
                                .max(source.entry_cone_radius_mm(z).unwrap())
                        );
                    } else {
                        assert_eq!(radius, original.radius(theta, z));
                    }
                    assert!(source.contains_void(
                        [
                            radius * theta.cos() * 0.999,
                            radius * theta.sin() * 0.999,
                            z
                        ],
                        0.
                    ));
                    assert!(!source.contains_void(
                        [
                            radius * theta.cos() * 1.001,
                            radius * theta.sin() * 1.001,
                            z
                        ],
                        0.
                    ));
                }
            }
            let knots = original.profile_breaks();
            assert!(source.entry_seam_station_mm(knots[2]).unwrap() > 0.399999999);
            assert!(source.entry_seam_station_mm(0.).unwrap() > 0.16);
            assert!(source.entry_cone_radius_mm(-0.001).is_none());
            assert!(source.entry_cone_radius_mm(f64::NAN).is_none());
            assert!(source.entry_seam_station_mm(f64::NAN).is_none());
            assert!(source.radius_mm(0., -0.1).is_nan());
            assert!(source.radius_mm(f64::INFINITY, 0.).is_nan());
            assert!(!source.contains_void([0., 0., -0.1], 0.));
            assert!(!source.contains_void([0., 0., 3.4], 0.));
            assert!(!source.contains_void([0.; 3], -1.));
            assert!(!source.contains_void([f64::NAN, 0., 0.], 0.));
        }
    }
}
#[test]
fn terminal_only_has_independent_analytic_volume_and_disabled_entry_keeps_original_sections() {
    let thread = tool(false, 1);
    let source = Cavity::new(
        thread,
        Intent {
            entry_depth_mm: None,
            terminal: both().terminal,
        },
    )
    .unwrap();
    assert!(source.entry_radius_mm().is_none());
    assert!(source.entry_cone_radius_mm(0.).is_none());
    assert!(source.entry_seam_station_mm(0.).is_none());
    let f = source.terminal_station_mm().unwrap();
    assert_eq!(source.radius_mm(0.73, 0.19), thread.radius(0.73, 0.19));
    assert_eq!(source.maximum_radius_bound_mm(), 1.55);
    let volume = source.removed_volume_mm3(Default::default()).unwrap();
    let independent = thread.volume_mm3() / thread.spec().length_mm * f
        + std::f64::consts::PI * 1.55_f64.powi(2) * (thread.spec().length_mm - f);
    assert!(
        volume.lower <= independent && independent <= volume.upper,
        "{volume:?} {independent}"
    );
}
#[test]
fn combined_volume_encloses_independent_dense_source_integration_and_is_hand_invariant() {
    let source = Cavity::new(tool(false, 1), both()).unwrap();
    let volume = source.removed_volume_mm3(Default::default()).unwrap();
    let other = Cavity::new(tool(true, 2), both())
        .unwrap()
        .removed_volume_mm3(Default::default())
        .unwrap();
    assert_eq!(volume, other);
    // Independent midpoint integration of the actual source over theta and z.
    // Its explicit coarse/fine convergence band is a numerical test oracle,
    // never the production source-volume certificate.
    let integrate = |n: usize| {
        let mut sum = 0.;
        let dt = std::f64::consts::TAU / n as f64;
        let dz = source.full_length_mm() / n as f64;
        for iz in 0..n {
            for it in 0..n {
                sum += source
                    .radius_mm((it as f64 + 0.5) * dt, (iz as f64 + 0.5) * dz)
                    .powi(2);
            }
        }
        sum * dt * dz / 2.
    };
    let coarse = integrate(512);
    let fine = integrate(1024);
    assert!((fine - (volume.lower + volume.upper) / 2.).abs() < 2. * (coarse - fine).abs() + 1e-5);
    assert!(volume.upper - volume.lower <= 1e-7);
}
#[test]
fn invalid_sources_and_query_budgets_do_not_change_the_original_thread() {
    let source = tool(false, 1);
    let before = source;
    for intent in [
        Intent {
            entry_depth_mm: None,
            terminal: None,
        },
        Intent {
            entry_depth_mm: Some(f64::NAN),
            terminal: None,
        },
        Intent {
            entry_depth_mm: Some(0.0001),
            terminal: None,
        },
        Intent {
            entry_depth_mm: Some(3.),
            terminal: both().terminal,
        },
        Intent {
            entry_depth_mm: None,
            terminal: Some(Terminal {
                width_mm: 0.3,
                radial_depth_mm: -1.,
            }),
        },
        Intent {
            entry_depth_mm: Some(1e8),
            terminal: None,
        },
        Intent {
            entry_depth_mm: None,
            terminal: Some(Terminal {
                width_mm: 0.3,
                radial_depth_mm: 1e8,
            }),
        },
    ] {
        assert!(Cavity::new(source, intent).is_err());
        assert_eq!(source, before);
    }
    let cavity = Cavity::new(source, both()).unwrap();
    assert_eq!(
        cavity.removed_volume_mm3(QueryOptions {
            tolerance_mm3: 0.,
            ..Default::default()
        }),
        Err(Error::Invalid)
    );
    for maximum_cells in [0, 1, 1_000_001] {
        assert_eq!(
            cavity.removed_volume_mm3(QueryOptions {
                maximum_cells,
                ..Default::default()
            }),
            Err(Error::Limit)
        );
    }
    assert_eq!(
        cavity.removed_volume_with_interrupt(Default::default(), Some(&|| true)),
        Err(Error::Cancelled)
    );
    assert!(
        cavity
            .removed_volume_mm3(QueryOptions {
                tolerance_mm3: 1e-300,
                maximum_cells: 1024
            })
            .is_err()
    );
    assert!(cavity.removed_volume_mm3(Default::default()).is_ok());
    assert_eq!(source, before);
}

#[test]
fn entry_only_scaling_and_unresolved_literal_root_keep_source_precision_explicit() {
    let original = tool(false, 1);
    let mut normalized = Vec::new();
    for scale in [0.1, 1., 1000.] {
        let spec = original.spec();
        let source = Thread::new(Spec {
            major_radius_mm: spec.major_radius_mm * scale,
            pitch_mm: spec.pitch_mm * scale,
            length_mm: spec.length_mm * scale,
            depth_mm: spec.depth_mm * scale,
            crest_width_mm: spec.crest_width_mm * scale,
            ..spec
        })
        .unwrap();
        let cavity = Cavity::new(
            source,
            Intent {
                entry_depth_mm: Some(0.4 * scale),
                terminal: None,
            },
        )
        .unwrap();
        assert!(cavity.terminal_station_mm().is_none());
        assert!(cavity.terminal_radius_mm().is_none());
        let volume = cavity
            .removed_volume_mm3(QueryOptions {
                tolerance_mm3: 1e-7 * scale.powi(3),
                ..Default::default()
            })
            .unwrap();
        normalized.push((volume.lower / scale.powi(3), volume.upper / scale.powi(3)));
        assert_eq!(
            cavity.radius_mm(0.3, 2. * scale),
            source.radius(0.3, 2. * scale)
        );
    }
    let intersection = (
        normalized
            .iter()
            .map(|v| v.0)
            .fold(f64::NEG_INFINITY, f64::max),
        normalized.iter().map(|v| v.1).fold(f64::INFINITY, f64::min),
    );
    assert!(intersection.0 <= intersection.1);
    let tiny = Thread::new(Spec {
        major_radius_mm: 1.,
        depth_mm: 1_f64.next_down(),
        pitch_mm: 4.,
        length_mm: 20.,
        crest_width_mm: 1.,
        ..original.spec()
    })
    .unwrap();
    assert_eq!(
        Cavity::new(
            tiny,
            Intent {
                entry_depth_mm: Some(0.1),
                terminal: None
            }
        ),
        Err(Error::Unresolved)
    );
    let spec = original.spec();
    let root = 0.05;
    let circle = Thread::new(Spec {
        root_radius_mm: Some(root),
        depth_mm: (spec.pitch_mm - spec.crest_width_mm) / (2. * 30_f64.to_radians().tan()) - root,
        ..spec
    })
    .unwrap();
    assert_eq!(Cavity::new(circle, both()), Err(Error::Invalid));
}

#[test]
fn source_maximum_encloses_actual_flank_coefficients_and_mid_query_cancel_is_retryable() {
    use std::cell::Cell;
    let source = Cavity::new(
        tool(false, 1),
        Intent {
            entry_depth_mm: Some(0.01),
            terminal: None,
        },
    )
    .unwrap();
    let knots = source.source_thread().profile_breaks();
    let mut maximum = source.source_thread().spec().major_radius_mm;
    for band in [1, 3] {
        let anchor = if band == 1 { knots[1] } else { knots[4] };
        let slope = source.source_thread().slope_at_phase(anchor, band);
        for q in [knots[band], knots[band + 1]] {
            let r = slope.mul_add(q - anchor, source.source_thread().spec().major_radius_mm);
            maximum = maximum.max(r);
        }
    }
    assert!(source.maximum_radius_bound_mm() >= maximum);
    let expected = source.removed_volume_mm3(Default::default()).unwrap();
    let calls = Cell::new(0);
    let cancel = || {
        calls.set(calls.get() + 1);
        calls.get() >= 12
    };
    assert_eq!(
        source.removed_volume_with_interrupt(Default::default(), Some(&cancel)),
        Err(Error::Cancelled)
    );
    assert_eq!(calls.get(), 12);
    assert_eq!(
        source.removed_volume_mm3(Default::default()).unwrap(),
        expected
    );
}

#[test]
fn large_finite_terminal_cavity_reports_global_precision_floor_after_small_cells_converge() {
    let source = Cavity::new(
        tool(false, 1),
        Intent {
            entry_depth_mm: None,
            terminal: Some(Terminal {
                width_mm: 0.3,
                radial_depth_mm: 1e7 - 1.5,
            }),
        },
    )
    .unwrap();
    assert_eq!(source.terminal_radius_mm(), Some(1e7));
    assert_eq!(
        source.removed_volume_mm3(Default::default()),
        Err(Error::Unresolved)
    );
    assert!(
        source
            .removed_volume_mm3(QueryOptions {
                tolerance_mm3: 1.,
                ..Default::default()
            })
            .is_ok()
    );
}
