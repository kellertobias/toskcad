use super::*;
use crate::authored_flat_thread::{AuthoredThread, SourceOptions};
use multicad_precision::PortableMath as _;
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
    // The flank knots are the source law itself, half + depth tan(30 deg) and
    // pitch - half - depth tan(30 deg). The depth is 0.15625 sqrt 3, so both
    // lie within ulps of rows 24 and 40 of the 64-row chart; the upper knot
    // must miss its row, which is what forces the 64-row repartition below.
    let spec = body.thread.spec();
    let flank = |deg: f64| spec.depth_mm * deg.to_radians().portable_tan();
    let half = spec.crest_width_mm * 0.5;
    let knots = [
        half + flank(spec.upper_flank_deg),
        (spec.pitch_mm - half) - flank(spec.lower_flank_deg),
    ];
    for (band, row, knot) in [(1, 24., knots[0]), (2, 40., knots[1])] {
        let phase = body.thread.bands()[band].phase_mm[1];
        assert_eq!(phase, knot);
        let exact = spec.pitch_mm * row / 64.;
        assert!((phase - exact).abs() <= 4. * f64::EPSILON * exact);
    }
    assert_ne!(
        body.thread.bands()[1].phase_mm[1],
        spec.pitch_mm * 24. / 64.
    );
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
    assert_eq!(
        (grid.subdivisions, grid.rows(), grid.shifted),
        (1, 64, true)
    );
    grid.advance(true).unwrap();
    assert_eq!((grid.subdivisions, grid.rows()), (1, 67));
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
    let mut angular = Grid::new();
    while angular.rows() < MAX_ROWS {
        angular.refine_angular().unwrap();
    }
    assert_eq!(angular.subdivisions, 1);
    assert!(matches!(angular.refine_angular(), Err(Error::Limit)));
    assert_eq!(angular.rows(), MAX_ROWS);
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
    let closed = close::close(shell, &body, quality(), &mut w)
        .unwrap()
        .unwrap();
    assert!(closed.source_chord_bound_mm <= quality().deflection);
    crate::manufacturing::certify(&closed.mesh).unwrap();
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

#[test]
fn actual_large_terminal_uses_local_conforming_shoulders_without_refining_the_thread() {
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
        panic!("side source must qualify")
    };
    let before = w.used;
    let tail = body.tail.unwrap();
    let chord = shell
        .top
        .iter()
        .map(|&e| {
            certificate::chord(
                &close::field(&body, e, tail.station, close::Surface::Shoulder),
                [0.; 3],
            )
        })
        .fold(0_f64, f64::max);
    eprintln!(
        "actual shoulder curvature-only bound={chord}, requested={}, sidework={before}, sidefacets={}",
        quality().deflection,
        shell.mesh.face_indices.len()
    );
    assert!(chord > quality().deflection);
    let start = std::time::Instant::now();
    let candidate = tessellate_candidate(&body, quality(), None);
    eprintln!(
        "before actual full candidate={:?}, elapsed={:?}",
        candidate.as_ref().map(|x| (
            x.mesh.face_indices.len(),
            x.work_used,
            x.source_chord_bound_mm
        )),
        start.elapsed()
    );
    let candidate = candidate.unwrap();
    assert!(candidate.source_chord_bound_mm <= quality().deflection);
    assert!(candidate.source_normal_angle_bound_rad <= quality().angle);
    let proof = crate::manufacturing::certify(&candidate.mesh).unwrap();
    assert!(proof.signed_volume_mm3 > 0. && proof.vertex_links > 0);
    assert!(candidate.work_used < 4_000_000);
}

#[test]
fn shifted_sampling_seam_preserves_actual_fine_pitch_source_and_closed_material() {
    use crate::{Circle, Rational, cylinder::AxialCylinder};
    let body = source(true);
    let original = body;
    let start = std::time::Instant::now();
    let candidate = tessellate_candidate(&body, quality(), None).unwrap();
    eprintln!(
        "fine-pitch cavity {} facets, {} work, {:?}",
        candidate.mesh.face_indices.len(),
        candidate.work_used,
        start.elapsed()
    );
    assert!(candidate.source_chord_bound_mm <= quality().deflection);
    assert!(candidate.source_normal_angle_bound_rad <= quality().angle);
    assert!(candidate.work_used < 16_000_000);
    assert!(
        crate::manufacturing::certify(&candidate.mesh)
            .unwrap()
            .signed_volume_mm3
            > 0.
    );
    assert_eq!(body, original);
    let zero = Rational::new(0, 1).unwrap();
    let stock = AxialCylinder::new(
        7,
        Circle::new([zero; 2], Rational::new(5, 2).unwrap()).unwrap(),
        zero,
        Rational::new(4, 1).unwrap(),
    )
    .unwrap();
    let material = super::super::super::Material::new(stock, body, 0.5).unwrap();
    let start = std::time::Instant::now();
    let candidate = material.tessellate_candidate(quality(), None).unwrap();
    eprintln!(
        "fine-pitch material {} facets, {} work, {:?}",
        candidate.mesh.face_indices.len(),
        candidate.work_used,
        start.elapsed()
    );
    assert!(candidate.source_chord_bound_mm <= quality().deflection);
    assert!(candidate.source_normal_angle_bound_rad <= quality().angle);
    assert!(
        crate::manufacturing::certify(&candidate.mesh)
            .unwrap()
            .signed_volume_mm3
            > 0.
    );
    assert_eq!(*material.cavity(), original);
}

#[test]
fn actual_shift_fallback_cancels_before_and_during_new_source_sampling() {
    let body = source(true);
    let original = body;
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
    let mut grid = Grid::new();
    for _ in 0..2 {
        assert!(matches!(
            chart::grid_attempt(&body, quality(), &grid, &mut w),
            Ok(Attempt::Repartition)
        ));
        grid.advance(true).unwrap();
    }
    assert!(grid.shifted);
    let boundary = calls.get();
    for stop in [boundary - 1, boundary, boundary + 1, boundary + 10_000] {
        calls.set(0);
        let cancel = || {
            calls.set(calls.get() + 1);
            calls.get() > stop
        };
        assert_eq!(
            tessellate_candidate(&body, quality(), Some(&cancel)).unwrap_err(),
            Error::Cancelled
        );
        assert_eq!(body, original);
    }
}

#[test]
fn actual_extreme_terminal_cap_requests_angular_refinement_and_keeps_real_work_limit() {
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
                radial_depth_mm: 100.,
            }),
            ..base.intent()
        },
    )
    .unwrap();
    let unchanged = body;
    let mut w = work();
    let Attempt::Complete(shell) = chart::shell_attempt(&body, quality(), 1, 64, &mut w).unwrap()
    else {
        panic!("actual side should qualify");
    };
    assert!(
        close::close(shell, &body, quality(), &mut w)
            .unwrap()
            .is_none()
    );
    assert_eq!(
        tessellate_candidate(&body, quality(), None).unwrap_err(),
        Error::Limit
    );
    assert_eq!(body, unchanged);
}

#[test]
fn actual_traced_retry_cancels_before_its_angular_table_allocation() {
    let body = source(true);
    let mut w = Work {
        remaining: 16_000_000,
        used: 0,
        interrupt: Some(&|| true),
    };
    assert!(matches!(
        chart::traced_grid_attempt(&body, quality(), &Grid::new(), &mut w),
        Err(Error::Cancelled)
    ));
}

#[test]
fn only_edges_far_below_facet_scale_propose_another_chart() {
    // The rail/row near-coincidence at a mouth: about 1e-15 mm at r = 1.63.
    let a = [-1.1521364980995048, -1.1521364980995044, 0.];
    let b = [-1.1521364980995041, -1.1521364980995048, 0.];
    assert!(conditioning_edge([a, b]));
    assert!(conditioning_edge([a, a]));
    // A micrometre edge, far below any print deflection, is still a facet.
    assert!(!conditioning_edge([a, [a[0] + 1e-6, a[1], a[2]]]));
    // The bound scales with the coordinates, not with an absolute length.
    let far = [1e6, 0., 0.];
    assert!(conditioning_edge([far, [1e6 + 1e-10, 0., 0.]]));
    assert!(!conditioning_edge([[0., 0., 0.], [1e-10, 0., 0.]]));
}
