use super::*;
use crate::{
    head_profile::{Head, Profile},
    helical_thread::{self, Thread},
    phillips_drive::{Cavity, Cutter, PhillipsScrew},
    profiled_screw::SourceSpec,
};
fn text(index: usize, posed: bool) -> String {
    let (s, d) = crate::phillips_drive::test_fixture(index);
    let cavity = Cavity::new(Cutter::new(s, d).unwrap(), 1e-5, 100_000).unwrap();
    let (rs, pitch, r, h, crown, land, neck) = [
        (1., 0.4, 2., 1.6, 3.2, 1.25, 0.1),
        (1.5, 0.5, 2.8, 2.4, 3., 2.1, 0.1),
        (1.5, 0.5, 2.75, 1.65, 0., 0., 0.1),
        (2., 0.7, 4., 3.1, 6.5, 3., 0.2),
        (3., 1., 6., 4.6, 10., 4.2, 0.25),
        (4., 1.25, 8., 6., 13., 5.4, 0.4),
    ][index];
    let head = if index == 2 {
        Head::Countersunk {
            radius_mm: r,
            height_mm: h,
            included_angle_deg: 90.,
            edge_chamfer_mm: 0.05,
            neck_radius_mm: neck,
        }
    } else {
        Head::Crowned {
            radius_mm: r,
            height_mm: h,
            crown_radius_mm: crown,
            top_land_radius_mm: land,
            neck_radius_mm: neck,
        }
    };
    let thread = Thread::new(helical_thread::Spec {
        major_radius_mm: rs - 0.03,
        pitch_mm: pitch,
        length_mm: 2.,
        depth_mm: 0.6 * pitch,
        crest_width_mm: 0.125 * pitch,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts: 1,
        left_hand: posed,
        root_radius_mm: None,
    })
    .unwrap();
    let screw = PhillipsScrew::new(
        thread,
        Profile::new_for_thread_union(head, rs, 0.).unwrap(),
        SourceSpec {
            head: Some(head),
            shank_length_mm: 0.,
            minimum_wall_mm: 1e-5,
        },
        cavity,
        100_000,
        &|| false,
    )
    .unwrap();
    let (origin, u, v) = if posed {
        ([105., -20., 60.], [0.8, 0.6, 0.], [0.36, -0.48, -0.8])
    } else {
        ([0.; 3], [1., 0., 0.], [0., 1., 0.])
    };
    crate::step::write_phillips_screw(&screw, origin, u, v, 1e-7).unwrap()
}
#[test]
fn actual_metric_file_declared_cubic_branches_are_qualified_without_a_source_recipe() {
    let text = text(0, false);
    let b = read_phillips_screw_boundary(&text, 1e-7, 1_000_000).unwrap();
    let c0 = b
        .cubic_support_bounds(
            CubicSupportOptions {
                maximum_distance_mm: 1e-7,
                work_limit: 10_000_000,
            },
            None,
        )
        .unwrap();
    let tube = c0.intersection_tubes(Options::default(), None).unwrap();
    assert!(std::ptr::eq(tube.boundary(), &b));
    assert_eq!(tube.edges().len(), 32);
    let mut largest = 0_f64;
    let mut derivative = 0_f64;
    for edge in tube.edges() {
        assert_eq!(edge.edge_index(), edge.edge);
        assert_eq!(edge.section_range(), edge.sections);
        assert_eq!(edge.section_direction(), edge.functional.direction);
        assert_eq!(edge.pivot_coordinate(), edge.functional.pivot);
        assert!(matches!(
            b.edges()[edge.edge].curve,
            Curve::CubicIntersection { .. }
        ));
        assert_eq!(edge.functional.direction[edge.functional.pivot], 1.);
        assert!(!edge.sections.is_empty());
        let mut previous = None;
        for cell in &tube.sections()[edge.sections.clone()] {
            assert_eq!(cell.span_index(), cell.span);
            assert_eq!(cell.parameter_range(), cell.proof.parameter);
            assert_eq!(cell.branch_coordinate_bounds_mm(), cell.proof.source_box);
            assert_eq!(cell.position_error_bounds_mm(), cell.proof.correction_mm);
            assert_eq!(
                cell.derivative_error_bounds_mm(),
                cell.proof.derivative_correction_mm
            );
            assert_eq!(cell.contraction_bound(), cell.proof.contraction);
            assert!(cell.proof.parameter[0] < cell.proof.parameter[1]);
            assert!(
                cell.proof
                    .source_box
                    .into_iter()
                    .all(|[a, b]| a.is_finite() && a <= b)
            );
            assert!(cell.proof.contraction < 0.5);
            largest = largest.max(
                cell.proof
                    .correction_mm
                    .into_iter()
                    .map(|v| v * v)
                    .sum::<f64>()
                    .sqrt(),
            );
            derivative = derivative.max(
                cell.proof
                    .derivative_correction_mm
                    .into_iter()
                    .fold(0., f64::max),
            );
            if let Some((span, end)) = previous {
                if span == cell.span {
                    assert_eq!(end, cell.proof.parameter[0]);
                } else {
                    assert_eq!(end, 1.);
                    assert_eq!(cell.proof.parameter[0], 0.);
                }
            }
            previous = Some((cell.span, cell.proof.parameter[1]));
        }
        assert_eq!(previous.unwrap().1, 1.);
    }
    assert!(largest <= 1e-7 && derivative.is_finite());
    assert!(tube.work_used() <= 16_000_000);
    assert!(tube.retained_capacity_bytes() > std::mem::size_of_val(&tube));
    assert!(tube.retained_capacity_bytes() < BYTES);
    eprintln!(
        "actualPhi tubes={} sections={} work={} correction={largest} derivative={derivative}",
        tube.edges().len(),
        tube.sections().len(),
        tube.work_used()
    );
}

#[test]
fn all_six_metric_files_and_proper_poses_earn_independent_declared_edge_tubes() {
    for index in 0..6 {
        for posed in [false, true] {
            let b = read_phillips_screw_boundary(&text(index, posed), 1e-7, 1_000_000).unwrap();
            let c0 = b
                .cubic_support_bounds(
                    CubicSupportOptions {
                        maximum_distance_mm: 1e-7,
                        work_limit: 10_000_000,
                    },
                    None,
                )
                .unwrap();
            let tube = c0
                .intersection_tubes(Options::default(), None)
                .unwrap_or_else(|e| panic!("metric{index} posed{posed} refused{e:?}"));
            assert_eq!(tube.edges().len(), 32);
            assert!(tube.sections().iter().all(|c| c.proof.contraction < 0.5));
            eprintln!(
                "metric{index} posed{posed}: sections{} work{}",
                tube.sections().len(),
                tube.work_used()
            );
        }
    }
}

#[test]
fn source_bound_tube_queries_keep_real_work_bytes_cell_precision_and_cancellation_limits() {
    let input = text(0, false);
    let b = read_phillips_screw_boundary(&input, 1e-7, 1_000_000).unwrap();
    let options = Options::default();
    let query = |o: Options, interrupt: Option<&dyn Fn() -> bool>| {
        b.cubic_support_bounds(
            CubicSupportOptions {
                maximum_distance_mm: 1e-7,
                work_limit: 10_000_000,
            },
            None,
        )
        .unwrap()
        .intersection_tubes(o, interrupt)
    };
    for o in [
        Options {
            maximum_position_error_mm: 0.,
            ..options
        },
        Options {
            maximum_position_error_mm: f64::NAN,
            ..options
        },
        Options {
            maximum_cells: 0,
            ..options
        },
        Options {
            maximum_cells: CELLS + 1,
            ..options
        },
        Options {
            maximum_work: 0,
            ..options
        },
        Options {
            maximum_work: WORK + 1,
            ..options
        },
        Options {
            maximum_bytes: 0,
            ..options
        },
        Options {
            maximum_bytes: BYTES + 1,
            ..options
        },
    ] {
        assert!(matches!(query(o, None), Err(PhillipsReadError::Invalid)));
    }
    for o in [
        Options {
            maximum_cells: 1,
            ..options
        },
        Options {
            maximum_work: 1,
            ..options
        },
        Options {
            maximum_bytes: 1,
            ..options
        },
    ] {
        assert!(matches!(query(o, None), Err(PhillipsReadError::Limit)));
    }
    assert!(
        query(
            Options {
                maximum_position_error_mm: 1e-15,
                ..options
            },
            None
        )
        .is_err()
    );
    let calls = std::cell::Cell::new(0usize);
    let count = || {
        calls.set(calls.get() + 1);
        false
    };
    let complete = query(options, Some(&count)).unwrap();
    let total = calls.get();
    for stop in [0, 1, 2, 3, 31, 32, total / 2, total - 1] {
        calls.set(0);
        let cancel = || {
            calls.set(calls.get() + 1);
            calls.get() > stop
        };
        assert!(matches!(
            query(options, Some(&cancel)),
            Err(PhillipsReadError::Cancelled)
        ));
    }
    assert_eq!(
        query(options, None).unwrap().sections().len(),
        complete.sections().len()
    );
    assert!(std::ptr::eq(complete.boundary(), &b));
}

#[test]
fn an_actual_off_support_pole_cannot_be_promoted_from_a_loose_residual_receipt() {
    let text = text(0, false);
    let original = read_phillips_screw_boundary(&text, 1e-7, 1_000_000).unwrap();
    let (id, point) = original
        .edges()
        .iter()
        .find_map(|e| {
            if let Curve::CubicIntersection {
                pole_entities,
                spans,
                ..
            } = &e.curve
            {
                Some((pole_entities[1], spans[0].poles_mm[1]))
            } else {
                None
            }
        })
        .unwrap();
    let prefix = format!("#{id}=");
    let changed = text
        .lines()
        .map(|line| {
            if line.starts_with(&prefix) {
                format!(
                    "#{id}=CARTESIAN_POINT('',({:.16E},{:.16E},{:.16E}));",
                    point[0] + 0.01,
                    point[1],
                    point[2]
                )
            } else {
                line.to_owned()
            }
        })
        .collect::<Vec<_>>()
        .join("\n");
    assert_ne!(changed, text);
    let actual = read_phillips_screw_boundary(&changed, 1e-7, 1_000_000).unwrap();
    let c0 = actual
        .cubic_support_bounds(
            CubicSupportOptions {
                maximum_distance_mm: 0.1,
                work_limit: 10_000_000,
            },
            None,
        )
        .unwrap();
    assert!(
        c0.edges()
            .iter()
            .any(|e| e.owner_distances_mm.iter().any(|v| *v > 1e-7))
    );
    assert!(c0.intersection_tubes(Options::default(), None).is_err());
    assert_eq!(original.edges().len(), actual.edges().len());
}

#[test]
fn actual_scratch_capacity_and_cancellation_are_checked_before_cell_processing() {
    let mut stack = Vec::<Cell>::new();
    stack.try_reserve_exact(300).unwrap();
    let mut budget = Budget {
        remaining: WORK,
        interrupt: None,
    };
    assert!(matches!(
        scratch(allocation(&stack), &mut budget),
        Err(PhillipsReadError::Limit)
    ));
    let cancel = || true;
    budget.interrupt = Some(&cancel);
    assert!(matches!(
        scratch(0, &mut budget),
        Err(PhillipsReadError::Cancelled)
    ));
}
