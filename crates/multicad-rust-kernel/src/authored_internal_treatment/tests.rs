use super::*;
use crate::authored_flat_thread::SourceOptions;
use crate::helical_thread::Spec;
fn thread(left_hand: bool, starts: u32, scale: f64) -> AuthoredThread {
    let pitch_mid = (2.675 + 2.775) / 2.;
    let minor_mid = (2.459 + 2.599) / 2.;
    AuthoredThread::new(
        Spec {
            major_radius_mm: 1.5 * scale,
            pitch_mm: 0.5 * scale,
            length_mm: 3.3 * scale,
            depth_mm: (3. - minor_mid) / 2. * scale,
            crest_width_mm: (0.25 - (3. - pitch_mid) / 3_f64.sqrt()) * scale,
            lower_flank_deg: 30.,
            upper_flank_deg: 30.,
            starts,
            left_hand,
            root_radius_mm: None,
        },
        SourceOptions {
            position_limit_mm: 1e-7 * scale,
            normal_angle_limit_rad: 1e-7,
        },
    )
    .unwrap()
}
fn intent(scale: f64) -> Intent {
    Intent {
        entry_depth_mm: Some(0.4 * scale),
        terminal: Some(Terminal {
            width_mm: 0.3 * scale,
            radial_depth_mm: 0.05 * scale,
        }),
    }
}
#[test]
fn actual_shared_max_cavity_retains_full_source_and_common_band_seams() {
    for left in [false, true] {
        for starts in [1, 2, 3, 8] {
            let original = thread(left, starts, 1.);
            let source = Cavity::new(original, intent(1.)).unwrap();
            assert_eq!(source.source_thread(), &original);
            assert_eq!(source.intent(), intent(1.));
            assert_eq!(source.full_length_mm(), 3.3);
            assert_eq!(
                source.entry_radius_mm(),
                Some(original.minimum_radius_mm() + 0.4)
            );
            assert_eq!(source.terminal_radius_mm(), Some(1.55));
            assert_eq!(
                source.maximum_radius_mm(),
                source.entry_radius_mm().unwrap()
            );
            assert_eq!(source.entry_cone_radius_mm(0.), source.entry_radius_mm());
            assert_eq!(
                source.entry_cone_radius_mm(0.4),
                Some(original.minimum_radius_mm())
            );
            for band in 0..4 {
                assert_eq!(
                    source.entry_seam_station_in_band(band, 1.).unwrap(),
                    source.entry_seam_station_in_band(band + 1, 0.).unwrap()
                );
            }
            assert_eq!(
                source.entry_seam_station_in_band(2, 0.5).unwrap(),
                Some(0.4)
            );
            for theta in [0., 0.43, 2.1] {
                for z in [0., 0.01, 0.2, 0.4, 0.71, 2.5, 3.1, 3.3] {
                    let radius = source.radius_mm(theta, z);
                    assert!(source.contains_void(
                        [
                            radius * theta.portable_cos() * 0.999,
                            radius * theta.portable_sin() * 0.999,
                            z
                        ],
                        0.
                    ));
                    assert!(!source.contains_void(
                        [
                            radius * theta.portable_cos() * 1.001,
                            radius * theta.portable_sin() * 1.001,
                            z
                        ],
                        0.
                    ));
                    if z >= source.terminal_station_mm().unwrap() {
                        assert_eq!(radius, 1.55);
                    }
                }
            }
            assert!(source.entry_cone_radius_mm(-0.01).is_none());
            assert!(source.entry_cone_radius_mm(f64::NAN).is_none());
            assert!(source.radius_mm(f64::NAN, 0.).is_nan());
            assert!(source.radius_mm(0., -1.).is_nan());
            assert!(!source.contains_void([0., 0., -1.], 0.));
            assert!(!source.contains_void([0., 0., 4.], 0.));
            assert!(!source.contains_void([f64::INFINITY, 0., 0.], 0.));
            assert!(!source.contains_void([0.; 3], -1.));
            assert_eq!(
                source.entry_seam_station_in_band(5, 0.),
                Err(Error::Invalid)
            );
            assert_eq!(
                source.entry_seam_station_in_band(0, -1.),
                Err(Error::Invalid)
            );
        }
    }
}
fn independent_closed_form(source: &Cavity) -> f64 {
    let thread = source.source_thread();
    let p = thread.spec().pitch_mm;
    let f = source
        .terminal_station_mm()
        .unwrap_or(source.full_length_mm());
    let a = source.entry_radius_mm();
    let c = source.intent().entry_depth_mm;
    let integral = thread
        .bands()
        .into_iter()
        .map(|b| {
            let w = b.phase_mm[1] - b.phase_mm[0];
            let [r, s] = b.radius_mm;
            let second = w * (r * r + r * s + s * s) / 3.;
            let third = w * (r.powi(3) + r * r * s + r * s * s + s.powi(3)) / 4.;
            if let (Some(a), Some(c)) = (a, c) {
                // This fixture covers the whole profile with the entry cone, so
                // its MAX column is a cubic on each affine band (no clipping root).
                assert!(a > thread.maximum_radius_mm());
                let k = (a - thread.minimum_radius_mm()) / c;
                second * f + (a.powi(3) * w - 3. * a * second + 2. * third) / (3. * k)
            } else {
                second * f
            }
        })
        .sum::<f64>();
    std::f64::consts::PI / p * integral
        + source.terminal_radius_mm().map_or(0., |r| {
            std::f64::consts::PI * r * r * (source.full_length_mm() - f)
        })
}
#[test]
fn exact_turn_volume_encloses_independent_cubic_integrals_for_both_hands_and_scales() {
    let mut normalized = Vec::new();
    for scale in [0.1, 1., 1000.] {
        for left in [false, true] {
            let source = Cavity::new(thread(left, 3, scale), intent(scale)).unwrap();
            let volume = source
                .removed_volume_mm3(QueryOptions {
                    tolerance_mm3: 1e-7 * scale.powi(3),
                    ..Default::default()
                })
                .unwrap();
            let exact = independent_closed_form(&source);
            assert!(
                volume.lower <= exact && exact <= volume.upper,
                "{volume:?} {exact}"
            );
            normalized.push((volume.lower / scale.powi(3), volume.upper / scale.powi(3)));
        }
    }
    assert!(
        normalized
            .iter()
            .map(|v| v.0)
            .fold(f64::NEG_INFINITY, f64::max)
            <= normalized.iter().map(|v| v.1).fold(f64::INFINITY, f64::min)
    );
}
#[test]
fn independently_owned_entry_only_small_entry_and_terminal_only_sources_remain_distinct() {
    let thread = thread(false, 1, 1.);
    let tail = Cavity::new(
        thread,
        Intent {
            entry_depth_mm: None,
            terminal: intent(1.).terminal,
        },
    )
    .unwrap();
    assert!(tail.entry_radius_mm().is_none());
    assert!(tail.entry_cone_radius_mm(0.).is_none());
    assert_eq!(tail.entry_seam_station_in_band(0, 0.5).unwrap(), None);
    let volume = tail.removed_volume_mm3(Default::default()).unwrap();
    let exact = independent_closed_form(&tail);
    assert!(volume.lower <= exact && exact <= volume.upper);
    for depth in [0.01, 0.4] {
        let entry = Cavity::new(
            thread,
            Intent {
                entry_depth_mm: Some(depth),
                terminal: None,
            },
        )
        .unwrap();
        assert!(entry.terminal_radius_mm().is_none());
        assert!(entry.terminal_station_mm().is_none());
        assert!(entry.entry_seam_station_in_band(1, 0.5).unwrap().is_some());
        let volume = entry.removed_volume_mm3(Default::default()).unwrap();
        if depth == 0.4 {
            let exact = independent_closed_form(&entry);
            assert!(volume.lower <= exact && exact <= volume.upper);
        } else {
            assert_eq!(entry.maximum_radius_mm(), thread.maximum_radius_mm());
            let integrate = |n: usize| {
                let mut sum = 0.;
                let dz = entry.full_length_mm() / n as f64;
                let da = std::f64::consts::TAU / n as f64;
                for z in 0..n {
                    for a in 0..n {
                        sum += entry
                            .radius_mm((a as f64 + 0.5) * da, (z as f64 + 0.5) * dz)
                            .powi(2);
                    }
                }
                sum * dz * da / 2.
            };
            let coarse = integrate(512);
            let fine = integrate(1024);
            assert!(
                (fine - (volume.lower + volume.upper) / 2.).abs()
                    < 2. * (fine - coarse).abs() + 1e-5
            );
        }
    }
}
#[test]
fn invalid_source_cancel_work_and_global_precision_floor_never_publish_a_partial_cavity() {
    use std::cell::Cell;
    let thread = thread(false, 1, 1.);
    for bad in [
        Intent {
            entry_depth_mm: None,
            terminal: None,
        },
        Intent {
            entry_depth_mm: Some(f64::NAN),
            terminal: None,
        },
        Intent {
            entry_depth_mm: Some(3.),
            terminal: intent(1.).terminal,
        },
        Intent {
            entry_depth_mm: None,
            terminal: Some(Terminal {
                width_mm: 0.3,
                radial_depth_mm: 1e8,
            }),
        },
    ] {
        assert!(Cavity::new(thread, bad).is_err());
    }
    let cavity = Cavity::new(thread, intent(1.)).unwrap();
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
        cavity.removed_volume_mm3(QueryOptions {
            tolerance_mm3: 0.,
            ..Default::default()
        }),
        Err(Error::Invalid)
    );
    assert!(
        cavity
            .removed_volume_mm3(QueryOptions {
                tolerance_mm3: 1e-300,
                maximum_cells: 1024
            })
            .is_err()
    );
    assert_eq!(
        cavity.removed_volume_with_interrupt(Default::default(), Some(&|| true)),
        Err(Error::Cancelled)
    );
    let calls = Cell::new(0);
    let cancel = || {
        calls.set(calls.get() + 1);
        calls.get() >= 12
    };
    assert_eq!(
        cavity.removed_volume_with_interrupt(Default::default(), Some(&cancel)),
        Err(Error::Cancelled)
    );
    assert_eq!(calls.get(), 12);
    assert!(cavity.removed_volume_mm3(Default::default()).is_ok());
    let large = Cavity::new(
        thread,
        Intent {
            entry_depth_mm: None,
            terminal: Some(Terminal {
                width_mm: 0.3,
                radial_depth_mm: 1e7 - 1.5,
            }),
        },
    )
    .unwrap();
    assert_eq!(
        large.removed_volume_mm3(Default::default()),
        Err(Error::Unresolved)
    );
    assert!(
        large
            .removed_volume_mm3(QueryOptions {
                tolerance_mm3: 1.,
                ..Default::default()
            })
            .is_ok()
    );
}
