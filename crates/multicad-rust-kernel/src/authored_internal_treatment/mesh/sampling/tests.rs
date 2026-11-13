use super::*;
use crate::authored_flat_thread::{AuthoredThread, SourceOptions};
fn source(nonbinary: bool) -> Cavity {
    let thread = AuthoredThread::new(
        crate::helical_thread::Spec {
            major_radius_mm: 1.5,
            pitch_mm: if nonbinary { 0.35 } else { 0.5 },
            length_mm: 4.,
            depth_mm: if nonbinary { 0.15 } else { 0.27063293868263705 },
            crest_width_mm: if nonbinary { 0.08 } else { 0.0625 },
            lower_flank_deg: 30.,
            upper_flank_deg: 30.,
            starts: if nonbinary { 3 } else { 1 },
            left_hand: nonbinary,
            root_radius_mm: None,
        },
        SourceOptions {
            position_limit_mm: 1e-7,
            normal_angle_limit_rad: 1e-7,
        },
    )
    .unwrap();
    Cavity::new(
        thread,
        Intent {
            entry_depth_mm: Some(0.5),
            terminal: Some(Terminal {
                width_mm: 0.5,
                radial_depth_mm: 0.125,
            }),
        },
    )
    .unwrap()
}
fn quality() -> Tessellation {
    Tessellation {
        deflection: 0.1,
        angle: 0.3,
        triangle_limit: 1_000_000,
    }
}
fn work() -> Work<'static> {
    Work {
        remaining: 16_000_000,
        used: 0,
        interrupt: None,
    }
}
#[test]
fn real_nominal_partition_rejects_then_retiles_without_altering_original_source() {
    let body = source(false);
    let original = body;
    assert_eq!(body.thread.bands()[1].phase_mm[1], 0.18749999999999994);
    assert_eq!(body.thread.bands()[2].phase_mm[1], 0.31250000000000006);
    let mut w = work();
    assert!(matches!(
        chart::shell_attempt(&body, quality(), 1, 64, &mut w),
        Ok(Attempt::Repartition)
    ));
    // The existing Option helper preserves the same rejection for proof callers.
    assert!(
        chart::shell(&body, quality(), 1, 64, &mut w)
            .unwrap()
            .is_none()
    );
    let mut grid = Grid::new();
    grid.advance(true).unwrap();
    assert_eq!((grid.subdivisions, grid.rows()), (1, 67));
    assert!(matches!(
        chart::shell_attempt(&body, quality(), 1, 67, &mut w),
        Ok(Attempt::Complete(_))
    ));
    let candidate = tessellate_candidate(&body, quality(), None).unwrap();
    assert!(
        candidate.work_used <= 16_000_000
            && candidate.source_chord_bound_mm <= 0.1
            && candidate.source_normal_angle_bound_rad <= 0.3
    );
    assert_eq!(body, original);
}
#[test]
fn real_alternate_failure_returns_to_original_refinement_and_never_retries_errors() {
    let body = source(true);
    let mut w = work();
    let mut grid = Grid::new();
    assert!(matches!(
        chart::shell_attempt(&body, quality(), 1, 64, &mut w),
        Ok(Attempt::Repartition)
    ));
    grid.advance(true).unwrap();
    assert!(matches!(
        chart::shell_attempt(&body, quality(), 1, 67, &mut w),
        Ok(Attempt::Repartition)
    ));
    grid.advance(true).unwrap();
    assert_eq!((grid.subdivisions, grid.rows()), (2, 128));
    let mut cancelled = Work {
        remaining: 16_000_000,
        used: 0,
        interrupt: Some(&|| true),
    };
    assert!(matches!(
        chart::shell_attempt(&body, quality(), 1, 64, &mut cancelled),
        Err(Error::Cancelled)
    ));
    let mut limited = work();
    assert!(matches!(
        chart::shell_attempt(
            &source(false),
            Tessellation {
                triangle_limit: 8,
                ..quality()
            },
            1,
            64,
            &mut limited
        ),
        Err(Error::Limit)
    ));
    let mut refine = Grid::new();
    refine.advance(false).unwrap();
    assert_eq!((refine.subdivisions, refine.rows()), (2, 128));
    loop {
        if refine.advance(false).is_err() {
            break;
        }
    }
    assert_eq!(
        (refine.subdivisions, refine.rows()),
        (MAX_SUBDIVISIONS, MAX_ROWS)
    );
    assert!(matches!(refine.advance(true), Err(Error::Limit)));
}
#[test]
fn cancellation_at_failed_partition_and_new_allocation_never_publishes_a_partial_candidate() {
    let body = source(false);
    let calls = std::cell::Cell::new(0usize);
    let count = || {
        calls.set(calls.get() + 1);
        false
    };
    let mut w = Work {
        remaining: 16_000_000,
        used: 0,
        interrupt: Some(&count),
    };
    assert!(matches!(
        chart::shell_attempt(&body, quality(), 1, 64, &mut w),
        Ok(Attempt::Repartition)
    ));
    let first = calls.get();
    calls.set(0);
    tessellate_candidate(&body, quality(), Some(&count)).unwrap();
    let total = calls.get();
    for stop in [0, first - 1, first, first + 1, total - 1] {
        calls.set(0);
        let cancel = || {
            calls.set(calls.get() + 1);
            calls.get() > stop
        };
        assert!(matches!(
            tessellate_candidate(&body, quality(), Some(&cancel)),
            Err(Error::Cancelled)
        ));
    }
}

#[test]
fn a_real_large_terminal_refines_after_qualified_side_sources_and_strict_targets_refuse() {
    let base = source(false);
    let thread = AuthoredThread::new(
        crate::helical_thread::Spec {
            depth_mm: 0.25,
            crest_width_mm: 0.125,
            ..base.thread.spec()
        },
        SourceOptions {
            position_limit_mm: 1e-7,
            normal_angle_limit_rad: 1e-7,
        },
    )
    .unwrap();
    let body = Cavity::new(
        thread,
        Intent {
            terminal: Some(Terminal {
                width_mm: 0.5,
                radial_depth_mm: 10.,
            }),
            ..base.intent()
        },
    )
    .unwrap();
    let mut w = work();
    let Attempt::Complete(shell) = chart::shell_attempt(&body, quality(), 1, 64, &mut w).unwrap()
    else {
        panic!("actual side sources should already qualify")
    };
    assert!(
        close::close(shell, &body, quality(), &mut w)
            .unwrap()
            .is_none(),
        "the larger terminal has its own source chord proof"
    );
    let result = tessellate_candidate(&body, quality(), None);
    match result {
        Ok(candidate) => assert!(candidate.source_chord_bound_mm <= 0.1),
        Err(error) => assert_eq!(error, Error::Limit),
    }
    assert!(
        chart::traced_shell(
            &body,
            Tessellation {
                deflection: 1e-12,
                ..quality()
            },
            1,
            64,
            &mut work()
        )
        .unwrap()
        .is_none()
    );
}
