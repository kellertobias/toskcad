//! Private stage measurements. Production builds contain no probe hooks.
use super::*;
use std::{
    cell::{Cell, RefCell},
    time::{Duration, Instant},
};
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Stage {
    Prepare,
    Emit,
    Preflight,
    Check,
    Cells,
    Embedding,
}
#[derive(Default)]
struct Stats {
    active: Option<Stage>,
    time: [Duration; 6],
    work: [usize; 6],
    calls: [usize; 6],
    finished: usize,
}
thread_local! {static STATS:RefCell<Option<Stats>>=const{RefCell::new(None)};}
pub(super) type Stamp = Option<(Stage, Instant, usize)>;
pub(super) fn start(stage: Stage, work: &Work<'_>) -> Stamp {
    STATS.with(|s| {
        s.borrow_mut().as_mut().map(|stats| {
            stats.active = Some(stage);
            (stage, Instant::now(), work.remaining)
        })
    })
}
pub(super) fn end(stamp: Stamp, work: &Work<'_>) {
    if let Some((stage, at, before)) = stamp {
        STATS.with(|s| {
            if let Some(s) = s.borrow_mut().as_mut() {
                s.active = None;
                let i = stage as usize;
                s.time[i] += at.elapsed();
                s.work[i] += before - work.remaining;
                s.calls[i] += 1;
            }
        });
    }
}
pub(super) fn finished(work: usize) {
    STATS.with(|s| {
        if let Some(s) = s.borrow_mut().as_mut() {
            s.finished = work;
        }
    });
}
fn fixture() -> RadialBoundary {
    fixture_length(6.5, 1e-6)
}
fn fixture_length(length: f64, minimum_feature: f64) -> RadialBoundary {
    use crate::axial_pipe_thread::authored::{AuthoredSection, AuthoredSpec, AuthoredThread};
    let p = 25.4 / 28.;
    let source = AuthoredThread::new(
        AuthoredSpec {
            major_radius_mm: 9.728 / 2. - 4. / 32.,
            pitch_mm: p,
            length_mm: length,
            radial_taper: 1. / 32.,
            section: AuthoredSection::BsptRounded {
                depth_mm: p / (3. * 27.5_f64.to_radians().tan()),
            },
            crest_radius_limit_mm: None,
        },
        1e-7,
        minimum_feature,
    )
    .unwrap();
    source
        .raw_boundary_ownership(1e-7, None)
        .unwrap()
        .ownership
        .bind_surfaces(1e-7, None)
        .unwrap()
        .certify_axial(1e-7, minimum_feature, None)
        .unwrap()
        .certify_caps(None)
        .unwrap()
        .certify_ordering(None)
        .unwrap()
        .certify_material(1e-7, None)
        .unwrap()
}
fn hash(mesh: &DisplayMesh) -> u64 {
    let mut value = 0xcbf29ce484222325u64;
    for b in mesh
        .vertices
        .iter()
        .flat_map(|v| v.to_bits().to_le_bytes())
        .chain(mesh.indices.iter().flat_map(|v| v.to_le_bytes()))
        .chain(mesh.face_indices.iter().flat_map(|v| v.to_le_bytes()))
    {
        value = (value ^ u64::from(b)).wrapping_mul(0x100000001b3);
    }
    value
}
thread_local! {pub(super) static WORKERS:Cell<Option<usize>>=const{Cell::new(None)}; pub(super) static ACTIVE:Cell<bool>=const{Cell::new(false)}; pub(super) static MIN_FACETS:Cell<Option<usize>>=const{Cell::new(None)};}

#[test]
fn paired_parallel_bspt_preserves_every_emitted_bit_and_callback_thread() {
    let body = fixture();
    let original = body.input().clone();
    let coordinator = std::thread::current().id();
    let callbacks = RefCell::new(0usize);
    let interrupt = || {
        assert_eq!(std::thread::current().id(), coordinator);
        *callbacks.borrow_mut() += 1;
        false
    };
    let options = Tessellation {
        deflection: 0.1,
        angle: 0.3,
        triangle_limit: 1_000_000,
    };
    let mut meshes = Vec::new();
    for (workers, threshold) in [(0, 8192), (8, 8192)] {
        MIN_FACETS.with(|t| t.set(Some(threshold)));
        WORKERS.with(|w| w.set(Some(workers)));
        STATS.with(|s| *s.borrow_mut() = Some(Stats::default()));
        let at = Instant::now();
        let result = body.tessellate_with_interrupt(options, Some(&interrupt));
        WORKERS.with(|w| w.set(None));
        MIN_FACETS.with(|t| t.set(None));
        let stats = STATS.with(|s| s.borrow_mut().take().unwrap());
        let mesh = result.unwrap();
        eprintln!(
            "PAIRED workers={workers} threshold_facets={threshold} wall={:?} work={} cells_wall={:?} cells_work={} hash={:016x}",
            at.elapsed(),
            stats.finished,
            stats.time[Stage::Cells as usize],
            stats.work[Stage::Cells as usize],
            hash(&mesh)
        );
        for stage in [
            Stage::Prepare,
            Stage::Emit,
            Stage::Preflight,
            Stage::Check,
            Stage::Cells,
            Stage::Embedding,
        ] {
            let i = stage as usize;
            eprintln!(
                "PAIRED_STAGE workers={workers} threshold_facets={threshold} stage={stage:?} wall={:?} work={}",
                stats.time[i], stats.work[i]
            );
        }
        meshes.push(mesh);
    }
    for mesh in &meshes[1..] {
        assert_eq!(meshes[0], *mesh);
    }
    // Width-cohort planning adds valid equal-class subdivisions. Both paths
    // independently pass all source/facet and embedding gates before this
    // new deterministic literal mesh golden is frozen.
    assert_eq!(hash(&meshes[0]), 0x11756b273029d41c);
    assert_eq!(body.input(), &original);
    assert!(*callbacks.borrow() > 0);
}

#[test]
fn actual_parallel_callback_cancel_and_budget_exhaustion_publish_no_mesh() {
    let body = fixture();
    let original = body.input().clone();
    let coordinator = std::thread::current().id();
    let reached = Cell::new(false);
    let interrupt = || {
        assert_eq!(std::thread::current().id(), coordinator);
        let active = ACTIVE.with(Cell::get);
        reached.set(reached.get() || active);
        active
    };
    let options = Tessellation {
        deflection: 0.1,
        angle: 0.3,
        triangle_limit: 1_000_000,
    };
    {
        let (workers, threshold) = (8, 8192);
        MIN_FACETS.with(|t| t.set(Some(threshold)));
        WORKERS.with(|w| w.set(Some(workers)));
        let cancelled = build(&body, options, Some(&interrupt));
        WORKERS.with(|w| w.set(None));
        assert_eq!(cancelled, Err(Error::Cancelled));
        assert!(reached.get());
        reached.set(false);
        assert!(!ACTIVE.with(Cell::get));
        WORKERS.with(|w| w.set(Some(workers)));
        let limited = build_with_limits(
            &body,
            options,
            None,
            Limits {
                work: 230_000_000,
                ..DEFAULT_LIMITS
            },
        );
        WORKERS.with(|w| w.set(None));
        assert_eq!(limited, Err(Error::Limit));
        MIN_FACETS.with(|t| t.set(None));
    }
    assert_eq!(body.input(), &original);
    assert!(!ACTIVE.with(Cell::get));
}

pub(super) fn multi_u_cam() -> RadialBoundary {
    let mut input = crate::radial_tensor::tests::five_bands::graph(
        &[[5., 0.], [2., 5.], [-4., 3.], [-4., -3.], [2., -5.]],
        false,
    );
    for side in 0..5 {
        let poles = spline_poles(&input.curves[side]).expect("actual spline cam rim");
        let mut controls = poles.to_vec();
        let direction = [
            controls[0][0] + controls[9][0],
            controls[0][1] + controls[9][1],
        ];
        for (i, p) in controls.iter_mut().enumerate() {
            let c = (i * (9 - i)) as f64 / 72.;
            for k in 0..2 {
                p[k] += direction[k] * 0.4 * c;
            }
        }
        let mut left = vec![controls[0]];
        let mut right = vec![controls[9]];
        while controls.len() > 1 {
            controls = controls
                .windows(2)
                .map(|p| std::array::from_fn(|k| (p[0][k] + p[1][k]) * 0.5))
                .collect();
            left.push(controls[0]);
            right.push(*controls.last().unwrap());
        }
        right.reverse();
        left.extend_from_slice(&right[1..]);
        input.curves[side] = Curve::Spline {
            spans: 2,
            poles: left.clone(),
        };
        input.curves[5 + side] = Curve::Spline {
            spans: 2,
            poles: left.iter().map(|p| [p[0], p[1], 4.]).collect(),
        };
        let poles = left
            .iter()
            .flat_map(|p| (0..=9).map(move |j| [p[0], p[1], 4. * j as f64 / 9.]))
            .collect();
        input.surfaces[2 + side] = Surface::Tensor(RetainedTensor::new([2, 1], poles).unwrap());
    }
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

#[test]
fn general_multi_u_cam_parallel_cache_zero_tiny_and_budget_receipts() {
    let body = multi_u_cam();
    let original = body.input().clone();
    let options = Tessellation {
        deflection: 0.025,
        angle: 0.15,
        triangle_limit: 100_000,
    };
    for strips in [0, 1, 2, 128] {
        let mut baseline = None;
        let mut serial_work = 0;
        for workers in [0, 8] {
            MIN_FACETS.with(|t| t.set(Some(0)));
            WORKERS.with(|w| w.set(Some(workers)));
            STATS.with(|s| *s.borrow_mut() = Some(Stats::default()));
            let result = build_with_limits(
                &body,
                options,
                None,
                Limits {
                    strips,
                    ..DEFAULT_LIMITS
                },
            );
            WORKERS.with(|w| w.set(None));
            MIN_FACETS.with(|t| t.set(None));
            let stats = STATS.with(|s| s.borrow_mut().take().unwrap());
            let mesh = result.unwrap();
            if workers == 0 {
                serial_work = stats.finished;
            }
            eprintln!(
                "CAM_CACHE strips={strips} workers={workers} work={} facets={} hash={:016x}",
                stats.finished,
                mesh.indices.len() / 3,
                hash(&mesh)
            );
            if let Some(previous) = &baseline {
                assert_eq!(&mesh, previous);
            } else {
                baseline = Some(mesh);
            }
        }
        MIN_FACETS.with(|t| t.set(Some(0)));
        WORKERS.with(|w| w.set(Some(8)));
        let bounded = build_with_limits(
            &body,
            options,
            None,
            Limits {
                strips,
                work: serial_work + 1000,
                ..DEFAULT_LIMITS
            },
        );
        WORKERS.with(|w| w.set(None));
        MIN_FACETS.with(|t| t.set(None));
        assert_eq!(bounded.unwrap(), baseline.unwrap());
    }
    assert_eq!(body.input(), &original);
}

#[test]
fn external_callback_unwind_restores_parallel_activity_state() {
    let body = multi_u_cam();
    let original = body.input().clone();
    let coordinator = std::thread::current().id();
    let interrupt = || {
        assert_eq!(std::thread::current().id(), coordinator);
        assert!(!ACTIVE.with(Cell::get), "injected external callback panic");
        false
    };
    MIN_FACETS.with(|t| t.set(Some(0)));
    WORKERS.with(|w| w.set(Some(8)));
    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        build(
            &body,
            Tessellation {
                deflection: 0.025,
                angle: 0.15,
                triangle_limit: 100_000,
            },
            Some(&interrupt),
        )
    }));
    MIN_FACETS.with(|t| t.set(None));
    WORKERS.with(|w| w.set(None));
    assert!(result.is_err());
    assert!(!ACTIVE.with(Cell::get));
    assert_eq!(body.input(), &original);
}

#[test]
fn actual_bspt_per_owner_static_cache_eligibility_keeps_certified_mesh() {
    fn assert_send_sync<T: Send + Sync>() {}
    assert_send_sync::<RadialBoundary>();
    assert_send_sync::<Prepared<'static>>();
    assert_send_sync::<tensor::Tensor>();
    let body = fixture();
    MIN_FACETS.with(|t| t.set(None));
    WORKERS.with(|w| w.set(Some(8)));
    STATS.with(|s| *s.borrow_mut() = Some(Stats::default()));
    let at = Instant::now();
    let result = build(
        &body,
        Tessellation {
            deflection: 0.1,
            angle: 0.3,
            triangle_limit: 1_000_000,
        },
        None,
    );
    WORKERS.with(|w| w.set(None));
    MIN_FACETS.with(|t| t.set(None));
    let stats = STATS.with(|s| s.borrow_mut().take().unwrap());
    let mesh = result.unwrap();
    eprintln!(
        "PER_OWNER_BSPT wall={:?} work={} hash={:016x}",
        at.elapsed(),
        stats.finished,
        hash(&mesh)
    );
    assert_eq!(mesh.indices.len() / 3, 186576);
    assert!(
        crate::manufacturing::certify(&mesh)
            .unwrap()
            .signed_volume_mm3
            > 0.
    );
    assert_eq!(hash(&mesh), 0x11756b273029d41c);
}

fn spline_poles(curve: &Curve) -> Option<&[Point]> {
    match curve {
        Curve::Spline { poles, .. } => Some(poles),
        _ => None,
    }
}
#[test]
fn cam_fixture_refuses_non_spline_rims() {
    assert!(
        spline_poles(&Curve::Circle {
            axis: Axis {
                origin: [0.; 3],
                normal: [0., 0., 1.],
                reference: [1., 0., 0.]
            },
            radius_mm: 1.
        })
        .is_none()
    );
}

#[test]
fn cancellation_inside_actual_large_grid_normal_preflight_preserves_source() {
    let body = fixture();
    let original = body.input().clone();
    STATS.with(|stats| *stats.borrow_mut() = Some(Stats::default()));
    let cancelled = Cell::new(false);
    let interrupt = || {
        let preflight =
            STATS.with(|stats| stats.borrow().as_ref().unwrap().active == Some(Stage::Preflight));
        cancelled.set(cancelled.get() || preflight);
        preflight
    };
    let result = build(
        &body,
        Tessellation {
            deflection: 0.1,
            angle: 0.3,
            triangle_limit: 1_000_000,
        },
        Some(&interrupt),
    );
    let stats = STATS.with(|stats| stats.borrow_mut().take().unwrap());
    assert_eq!(result, Err(Error::Cancelled));
    assert!(cancelled.get());
    assert_eq!(stats.active, Some(Stage::Preflight));
    assert!(stats.calls[Stage::Emit as usize] > 0);
    assert_eq!(body.input(), &original);
}

#[test]
fn actual_canonical_eight_mm_bspt_mesh_reaches_manufacturing_quality() {
    let body = fixture_length(8., 1e-5);
    let original = body.input().clone();
    STATS.with(|s| *s.borrow_mut() = Some(Stats::default()));
    let at = Instant::now();
    let result = body.tessellate_with_interrupt(
        Tessellation {
            deflection: 0.1,
            angle: 0.3,
            triangle_limit: 1_000_000,
        },
        None,
    );
    let stats = STATS.with(|s| s.borrow_mut().take().unwrap());
    eprintln!(
        "canonical8 BSPT elapsed={:?}, error={:?}, unfinished={:?}, work={}",
        at.elapsed(),
        result.as_ref().err(),
        stats.active,
        stats.finished
    );
    for stage in [
        Stage::Prepare,
        Stage::Emit,
        Stage::Preflight,
        Stage::Check,
        Stage::Cells,
        Stage::Embedding,
    ] {
        let i = stage as usize;
        eprintln!(
            "{stage:?} calls={} work={} time={:?}",
            stats.calls[i], stats.work[i], stats.time[i]
        );
    }
    let mesh = result.expect("actual canonical8 mesh");
    assert!(
        crate::manufacturing::certify(&mesh)
            .unwrap()
            .signed_volume_mm3
            > 0.
    );
    assert_eq!(body.input(), &original);
}

#[test]
fn actual_failed_transverse_cell_is_not_recovered_by_longitudinal_splits() {
    let body = fixture_length(8., 1e-5);
    let options = Tessellation {
        deflection: 0.1,
        angle: 0.3,
        triangle_limit: 1_000_000,
    };
    let mut w = Work {
        remaining: DEFAULT_LIMITS.work,
        interrupt: None,
    };
    let p = Prepared::new(&body, options, DEFAULT_LIMITS, &mut w).unwrap();
    let original: [[f64; 2]; 2] = [[0., 0.0625], [0.2676056338028169, 0.2680457746478873]];
    for axis in [2, 0, 1] {
        for part in 0..if axis == 2 { 1 } else { 2 } {
            let mut range = original;
            if axis < 2 {
                let m = range[axis][0].midpoint(range[axis][1]);
                range[axis][usize::from(part == 0)] = m;
            }
            let uv = [
                [range[0][0], range[1][0]],
                [range[0][1], range[1][0]],
                [range[0][1], range[1][1]],
                [range[0][0], range[1][1]],
            ];
            let vertices = uv.map(|uv| p.point(0, uv[0], uv[1], None, &mut w).unwrap());
            let actual = quad_triangles(vertices);
            let Side::Tensor(t) = &p.sides[0] else {
                panic!("tensor")
            };
            let mut cache = strip::Cache::new(t, DEFAULT_LIMITS.strips);
            cache.advance(range[1][0], &mut w).unwrap();
            let result = p.cell_check(
                CandidateCell {
                    side: 0,
                    range,
                    vertices,
                    actual,
                },
                None,
                Some(&mut cache),
                &mut w,
            );
            let expected = if axis == 0 {
                CellCheck::Match
            } else {
                CellCheck::Normal
            };
            assert!(result == Ok(expected));
        }
    }
}
