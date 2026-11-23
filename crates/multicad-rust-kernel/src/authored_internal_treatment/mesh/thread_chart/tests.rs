use super::*;
use std::{cell::Cell, collections::BTreeMap};
fn source(wood: bool, left: bool, starts: u32) -> AuthoredThread {
    if !wood {
        return *super::super::tests::source(left, starts, Some(0.4), true).source_thread();
    }
    AuthoredThread::new(
        crate::helical_thread::Spec {
            major_radius_mm: 5.,
            pitch_mm: 2.,
            length_mm: 4.,
            depth_mm: 1.,
            crest_width_mm: 0.4,
            lower_flank_deg: 30.,
            upper_flank_deg: 30.,
            starts,
            left_hand: left,
            root_radius_mm: None,
        },
        crate::authored_flat_thread::SourceOptions {
            position_limit_mm: 1e-7,
            normal_angle_limit_rad: 1e-7,
        },
    )
    .unwrap()
}
fn options() -> Tessellation {
    super::super::tests::options()
}
#[test]
fn actual_bare_thread_shell_has_two_owned_open_rims_and_whole_source_receipts() {
    for wood in [false, true] {
        for (left, starts) in [(false, 1), (true, 2)] {
            let source = source(wood, left, starts);
            let timer = std::time::Instant::now();
            let candidate = source
                .tessellate_shell(options(), 16_000_000, None)
                .unwrap();
            eprintln!(
                "bare shell wood{wood}/left{left}/N{starts}: {}facets {:?} {}work",
                candidate.mesh().face_indices.len(),
                timer.elapsed(),
                candidate.work_used()
            );
            assert_eq!(candidate.source(), &source);
            assert!(candidate.source_chord_bound_mm() <= options().deflection);
            assert!(candidate.source_normal_angle_bound_rad() <= options().angle);
            assert!(candidate.retained_bytes() < 64 * 1024 * 1024);
            assert!(candidate.work_used() > 0);
            let mesh = candidate.mesh();
            assert!(
                crate::manufacturing::certify(mesh).is_err(),
                "uncapped shell cannot be manufacturing material"
            );
            let mut edges = BTreeMap::<[usize; 2], (usize, i32)>::new();
            for tri in mesh.indices.as_chunks::<3>().0 {
                for [a, b] in [[tri[0], tri[1]], [tri[1], tri[2]], [tri[2], tri[0]]] {
                    let (a, b) = (a as usize, b as usize);
                    let key = [a.min(b), a.max(b)];
                    let e = edges.entry(key).or_default();
                    e.0 += 1;
                    e.1 += if a < b { 1 } else { -1 };
                }
            }
            let mut rim = BTreeMap::new();
            for end in [ShellEnd::Start, ShellEnd::End] {
                assert!(candidate.rim_count(end) > 2);
                assert_eq!(candidate.rim_edge(end, usize::MAX), None);
                let mut degrees = BTreeMap::<usize, usize>::new();
                for i in 0..candidate.rim_count(end) {
                    let e = candidate.rim_edge(end, i).unwrap();
                    let [a, b] = e.vertex_indices;
                    assert_ne!(a, b);
                    assert_eq!(
                        mesh.vertices[a * 3 + 2],
                        if end == ShellEnd::Start {
                            0.
                        } else {
                            source.full_length_mm()
                        }
                    );
                    assert!(
                        e.point_error_mm
                            .into_iter()
                            .all(|x| x.is_finite() && x >= 0.)
                    );
                    assert!(
                        e.phase_slope[0] <= e.phase_slope[1]
                            && e.phase_intercept_mm[0] <= e.phase_intercept_mm[1]
                    );
                    assert!(e.angular_turns.into_iter().all(|p| p[0] <= p[1]));
                    assert!(rim.insert([a.min(b), a.max(b)], end).is_none());
                    *degrees.entry(a).or_default() += 1;
                    *degrees.entry(b).or_default() += 1;
                }
                assert!(degrees.values().all(|d| *d == 2));
            }
            for (key, (count, sign)) in edges {
                if rim.contains_key(&key) {
                    assert_eq!(count, 1);
                } else {
                    assert_eq!((count, sign), (2, 0));
                }
            }
        }
    }
}
#[test]
fn source_shell_cancel_limits_invalid_quality_and_retry_are_atomic() {
    let source = source(false, false, 1);
    let before = source;
    let good = source
        .tessellate_shell(options(), 16_000_000, None)
        .unwrap();
    let calls = Cell::new(0);
    let count = || {
        calls.set(calls.get() + 1);
        false
    };
    source
        .tessellate_shell(options(), 16_000_000, Some(&count))
        .unwrap();
    let total = calls.get();
    for stop in [1, 2, 100, total] {
        let n = Cell::new(0);
        let cancel = || {
            n.set(n.get() + 1);
            n.get() == stop
        };
        assert!(matches!(
            source.tessellate_shell(options(), 16_000_000, Some(&cancel)),
            Err(SourceError::Cancelled)
        ));
    }
    assert!(matches!(
        source.tessellate_shell(options(), 0, None),
        Err(SourceError::Limit)
    ));
    assert!(matches!(
        source.tessellate_shell(options(), good.work_used() - 1, None),
        Err(SourceError::Limit)
    ));
    assert!(matches!(
        source.tessellate_shell(
            Tessellation {
                deflection: 0.,
                ..options()
            },
            16_000_000,
            None
        ),
        Err(SourceError::Invalid)
    ));
    assert!(matches!(
        source.tessellate_shell(
            Tessellation {
                triangle_limit: 1,
                ..options()
            },
            16_000_000,
            None
        ),
        Err(SourceError::Limit)
    ));
    assert!(matches!(
        source.tessellate_shell(
            Tessellation {
                deflection: 1e-12,
                ..options()
            },
            1_000_000,
            None
        ),
        Err(SourceError::Limit)
    ));
    assert_eq!(source, before);
    let retry = source
        .tessellate_shell(options(), 16_000_000, None)
        .unwrap();
    assert_eq!(retry.mesh().vertices, good.mesh().vertices);
    assert_eq!(retry.mesh().indices, good.mesh().indices);
    assert_eq!(retry.mesh().face_indices, good.mesh().face_indices);
    assert_eq!(failure(Error::Unresolved), SourceError::Unresolved);
}

#[test]
fn real_nominal_source_repartitions_without_changing_phase_or_admitting_noop_treatment() {
    let thread = AuthoredThread::new(
        crate::helical_thread::Spec {
            major_radius_mm: 1.5,
            pitch_mm: 0.5,
            length_mm: 4.,
            depth_mm: 0.27063293868263705,
            crest_width_mm: 0.0625,
            lower_flank_deg: 30.,
            upper_flank_deg: 30.,
            starts: 1,
            left_hand: false,
            root_radius_mm: None,
        },
        crate::authored_flat_thread::SourceOptions {
            position_limit_mm: 1e-7,
            normal_angle_limit_rad: 1e-7,
        },
    )
    .unwrap();
    // The lower flank knot is pitch - half - depth tan(30 deg), with depth
    // 0.15625 sqrt 3: the source law itself, within ulps of row 40 of the
    // 64-row chart. The upper knot sits within ulps of row 24 and misses it,
    // which forces the 64-row repartition asserted below.
    let spec = thread.spec();
    let flank = |deg: f64| spec.depth_mm * deg.to_radians().portable_tan();
    let half = spec.crest_width_mm * 0.5;
    let bands = thread.bands();
    assert_eq!(bands[1].phase_mm[1], half + flank(spec.upper_flank_deg));
    assert_eq!(
        bands[2].phase_mm[1],
        (spec.pitch_mm - half) - flank(spec.lower_flank_deg)
    );
    for (band, row) in [(1, 24.), (2, 40.)] {
        let exact = spec.pitch_mm * row / 64.;
        assert!((bands[band].phase_mm[1] - exact).abs() <= 4. * f64::EPSILON * exact);
    }
    assert_ne!(bands[1].phase_mm[1], spec.pitch_mm * 24. / 64.);
    assert!(
        Cavity::new(
            thread,
            Intent {
                entry_depth_mm: None,
                terminal: None
            }
        )
        .is_err()
    );
    let mut work = Work {
        remaining: 16_000_000,
        used: 0,
        interrupt: None,
    };
    assert!(matches!(
        chart::traced_chart(
            ThreadChart::bare(&thread),
            options(),
            (1, 64, false),
            &mut work
        )
        .unwrap(),
        Attempt::Repartition
    ));
    let candidate = thread
        .tessellate_shell(options(), 16_000_000, None)
        .unwrap();
    assert_eq!(*candidate.source(), thread);
    // Neither the refused treatment, the rejected chart nor the retiled shell
    // changes a source phase.
    assert_eq!(thread.bands(), bands);
    assert_eq!(candidate.source().bands(), bands);
    assert!(candidate.source_chord_bound_mm() <= options().deflection);
    assert!(candidate.source_normal_angle_bound_rad() <= options().angle);
    assert!(candidate.rim_count(ShellEnd::Start) > 2 && candidate.rim_count(ShellEnd::End) > 2);
}

#[test]
fn actual_shifted_source_footpoints_compare_by_exact_physical_turns_for_both_hands() {
    use num_rational::BigRational;
    use num_traits::{One, Zero};
    for left in [false, true] {
        let spec = source(false, left, 1).spec();
        let thread = AuthoredThread::new(
            spec,
            crate::authored_flat_thread::SourceOptions {
                position_limit_mm: 1e-7,
                normal_angle_limit_rad: 1e-7,
            },
        )
        .unwrap();
        let mut w = Work {
            remaining: 16_000_000,
            used: 0,
            interrupt: None,
        };
        let mut grid = Grid::new();
        grid.shifted = true;
        let shell = loop {
            match chart::traced_chart(
                ThreadChart::bare(&thread),
                options(),
                (grid.subdivisions, grid.rows(), grid.shifted),
                &mut w,
            )
            .unwrap()
            {
                Attempt::Complete((shell, trace)) => {
                    break ThreadShell {
                        thread,
                        shell,
                        trace,
                        work_used: w.used,
                        end_station_mm: thread.full_length_mm(),
                    };
                }
                Attempt::Refine => grid.advance(false).unwrap(),
                Attempt::Repartition => grid.advance(true).unwrap(),
            }
        };
        let before = *shell.source();
        let mut found_after_one = false;
        let mut selected = None;
        for end in [ShellEnd::Start, ShellEnd::End] {
            for i in 0..shell.rim_count(end) {
                let vertex = shell.rim_edge(end, i).unwrap().vertex_indices[0];
                let original = shell
                    .trace
                    .parameter(vertex, &mut w)
                    .unwrap()
                    .rational(&mut w)
                    .unwrap();
                found_after_one |= original >= BigRational::one();
                let (physical, cost) = shell.exact_rim_turns(vertex, 16_000_000, None).unwrap();
                assert!(
                    cost > 0 && physical >= BigRational::zero() && physical < BigRational::one()
                );
                let mut expected = if original >= BigRational::one() {
                    original - BigRational::one()
                } else {
                    original
                };
                if left && !expected.is_zero() {
                    expected = BigRational::one() - expected;
                }
                assert_eq!(physical, expected);
                selected = Some(vertex);
            }
        }
        assert!(
            found_after_one,
            "actual source has a footpoint at the full-period boundary"
        );
        let vertex = selected.unwrap();
        assert!(matches!(
            shell.exact_rim_turns(usize::MAX, 16_000_000, None),
            Err(SourceError::Unresolved)
        ));
        assert!(matches!(
            shell.exact_rim_turns(vertex, 0, None),
            Err(SourceError::Limit)
        ));
        let calls = Cell::new(0);
        let count = || {
            calls.set(calls.get() + 1);
            false
        };
        shell
            .exact_rim_turns(vertex, 16_000_000, Some(&count))
            .unwrap();
        let total = calls.get();
        for stop in 1..=total {
            let n = Cell::new(0);
            let cancel = || {
                n.set(n.get() + 1);
                n.get() == stop
            };
            assert!(matches!(
                shell.exact_rim_turns(vertex, 16_000_000, Some(&cancel)),
                Err(SourceError::Cancelled)
            ));
        }
        assert_eq!(*shell.source(), before);
    }
}
#[test]
fn actual_finite_prefix_keeps_original_source_phase_and_earns_both_rims() {
    let source = source(false, false, 1);
    let end = source.full_length_mm() - 0.001;
    let candidate = source
        .tessellate_prefix_shell(end, options(), 16_000_000, None)
        .unwrap();
    assert_eq!(candidate.source(), &source);
    assert_eq!(candidate.end_station_mm(), end);
    assert!(candidate.source_chord_bound_mm() <= options().deflection);
    assert!(candidate.source_normal_angle_bound_rad() <= options().angle);
    for end_kind in [ShellEnd::Start, ShellEnd::End] {
        let station = if end_kind == ShellEnd::Start { 0. } else { end };
        for i in 0..candidate.rim_count(end_kind) {
            let rim = candidate.rim_edge(end_kind, i).unwrap();
            for v in rim.vertex_indices {
                assert_eq!(candidate.mesh().vertices[3 * v + 2], station)
            }
        }
    }
    assert!(crate::manufacturing::certify(candidate.mesh()).is_err());
    assert!(
        source
            .tessellate_prefix_shell(f64::NAN, options(), 16_000_000, None)
            .is_err()
    );
    assert!(
        source
            .tessellate_prefix_shell(source.spec().pitch_mm / 2., options(), 16_000_000, None)
            .is_err()
    );
    assert!(
        source
            .tessellate_prefix_shell(source.full_length_mm() + 0.01, options(), 16_000_000, None)
            .is_err()
    );
    assert!(
        source
            .tessellate_prefix_shell(end, options(), candidate.work_used() - 1, None)
            .is_err()
    );
    assert_eq!(
        source
            .tessellate_prefix_shell(end, options(), 16_000_000, Some(&|| true))
            .err(),
        Some(crate::authored_flat_thread::Error::Cancelled)
    );
}
#[test]
fn alternate_prefix_partition_retains_original_geometry_and_actual_work_frontiers() {
    let thread = source(false, false, 1);
    let end = thread.full_length_mm() - 0.001;
    let candidate = thread
        .tessellate_prefix_shell_shifted(end, options(), 16_000_000, None)
        .unwrap();
    assert_eq!(candidate.source(), &thread);
    assert_eq!(candidate.end_station_mm(), end);
    for side in [ShellEnd::Start, ShellEnd::End] {
        let z = if side == ShellEnd::Start { 0. } else { end };
        for i in 0..candidate.rim_count(side) {
            for vertex in candidate.rim_edge(side, i).unwrap().vertex_indices {
                assert_eq!(candidate.mesh().vertices[3 * vertex + 2], z)
            }
        }
    }
    assert!(
        thread
            .tessellate_prefix_shell_shifted(end, options(), candidate.work_used() - 1, None)
            .is_err()
    );
    assert_eq!(
        thread
            .tessellate_prefix_shell_shifted(f64::NAN, options(), 16_000_000, None)
            .err(),
        Some(SourceError::Invalid)
    );
    assert_eq!(
        thread
            .tessellate_prefix_shell_shifted(end, options(), 16_000_000, Some(&|| true))
            .err(),
        Some(SourceError::Cancelled)
    );
}
