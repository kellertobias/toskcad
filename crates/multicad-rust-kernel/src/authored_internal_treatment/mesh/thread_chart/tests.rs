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
    assert_eq!(thread.bands()[2].phase_mm[1], 0.31250000000000006);
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
    assert!(candidate.source_chord_bound_mm() <= options().deflection);
    assert!(candidate.source_normal_angle_bound_rad() <= options().angle);
    assert!(candidate.rim_count(ShellEnd::Start) > 2 && candidate.rim_count(ShellEnd::End) > 2);
}
