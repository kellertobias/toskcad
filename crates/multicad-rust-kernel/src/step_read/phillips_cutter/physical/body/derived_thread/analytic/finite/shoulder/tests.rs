use super::*;
#[test]
fn actual_m2_shared_crest_shoulder_has_original_complement_plane_and_bounded_transfer() {
    degree11::tests::with_head(&mut |head| {
        let field = head
            .common_chart_periodic_thread_repair(Default::default(), None)
            .unwrap();
        let finite = head
            .periodic_thread_finite_curves(&field, Default::default(), None)
            .unwrap();
        let start = std::time::Instant::now();
        let join = head
            .periodic_shared_crest_shoulder(&finite, Default::default(), None)
            .unwrap();
        assert!(std::ptr::eq(join.finite_curves(), &finite));
        assert!(std::ptr::eq(join.head_source(), head.source()));
        assert_eq!(join.head_bottom_station_mm(), [1., 0.]);
        let c = join.original_complement();
        assert_eq!(
            (
                c.face_entity,
                c.edge_entity,
                c.curve_entity,
                join.original_plane_entity()
            ),
            (37880, 3463, 3462, 37879)
        );
        assert_eq!(c.vertex_entities, [2, 8]);
        assert_eq!(c.kind, PeriodicFiniteKind::Circle);
        assert!(
            c.distance_mm < 1e-7 && c.cell_derivative_mm < 1e-7 && join.plane_distance_mm() < 1e-7
        );
        assert_eq!(join.inherited_work_used(), finite.total_work_used());
        assert_eq!(
            join.total_work_used(),
            join.inherited_work_used() + join.new_work_used()
        );
        assert!(join.retained_budget_bytes() + SCRATCH < BYTES);
        assert_eq!(
            join.periodic_source().phase_node_expressions_mm(),
            field.source().phase_node_expressions_mm()
        );
        let ledger = || Work {
            remaining: WORK,
            interrupt: None,
            cancelled: Cell::new(false),
        };
        let source = *field.source();
        let stations = finite.cap_stations_mm();
        let bottom = join.head_bottom_station_mm();
        for (b, z) in [
            ([2., bottom[1]], stations),
            ([bottom[0], stations[0] + 0.01], stations),
            (bottom, [f64::NAN, stations[1]]),
            (bottom, [stations[0], f64::INFINITY]),
            (bottom, [stations[0], stations[0]]),
            (bottom, [stations[0], source.length_mm() + 1.]),
        ] {
            assert!(domain::Domain::new(source, b, z, &mut ledger()).is_err());
        }
        let mut work = ledger();
        work.remaining = 0;
        assert!(domain::Domain::new(source, bottom, stations, &mut work).is_err());
        let r = join
            .radial_bounds_mm(1, [0., 1.], [0., 0.], 128, None)
            .unwrap();
        let root = source.radius_nodes_mm()[1];
        assert!(r[0] <= root && r[1] >= root);
        let r = join
            .radial_bounds_mm(1, [0., 1.], [1., 1.], 128, None)
            .unwrap();
        assert!(r[0] <= bottom[0] && r[1] >= bottom[0]);
        for (band, u, t) in [
            (3, [0., 1.], [0., 1.]),
            (0, [f64::NAN, 0.], [0., 1.]),
            (0, [-0.1, 1.], [0., 1.]),
            (0, [0., 1.1], [0., 1.]),
            (0, [0., 1.], [f64::NAN, 0.]),
            (0, [0., 1.], [-0.1, 1.]),
            (0, [0., 1.], [0., 1.1]),
        ] {
            assert!(matches!(
                join.radial_bounds_mm(band, u, t, 128, None),
                Err(PhillipsReadError::Invalid)
            ));
        }
        assert!(matches!(
            join.radial_bounds_mm(0, [0., 1.], [0., 1.], WORK + 1, None),
            Err(PhillipsReadError::Invalid)
        ));
        assert!(matches!(
            join.radial_bounds_mm(0, [0., 1.], [0., 1.], 127, None),
            Err(PhillipsReadError::Limit)
        ));
        assert!(matches!(
            join.radial_bounds_mm(0, [0., 1.], [0., 1.], 128, Some(&|| true)),
            Err(PhillipsReadError::Cancelled)
        ));

        eprintln!(
            "SHOULDER qualified {:?} new{} inherited{} total{} retained{} plane{} C0{} D1{} incidence{}",
            start.elapsed(),
            join.new_work_used(),
            join.inherited_work_used(),
            join.total_work_used(),
            join.retained_budget_bytes(),
            join.plane_distance_mm(),
            c.distance_mm,
            c.cell_derivative_mm,
            c.endpoint_incidence_reserve_mm
        );
    });
}
#[test]
fn whole_disk_plane_uses_normalized_numeric_family_and_honest_refusals() {
    let ledger = || Work {
        remaining: WORK,
        interrupt: None,
        cancelled: Cell::new(false),
    };
    let p = I::point;
    let o = [p(0.); 3];
    let n = [p(0.), p(0.), p(2.)];
    let u = [p(3.), p(0.), p(0.)];
    let chart = CommonMouthChart {
        columns: [[0.; 3], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
        outer_top_z: 0.,
    };
    assert!(plane::disk_transfer(o, n, u, chart, [2., 0.], &mut ledger()).unwrap() < 1e-12);
    let mut shifted = chart;
    shifted.columns[0][2] = 0.5;
    let d = plane::disk_transfer(o, n, u, shifted, [2., 0.], &mut ledger()).unwrap();
    assert!((0.5..0.5000000001).contains(&d));
    for bottom in [[0., 0.], [-1., 0.], [f64::NAN, 0.], [1., f64::INFINITY]] {
        assert!(plane::disk_transfer(o, n, u, chart, bottom, &mut ledger()).is_err());
    }
    let mut bad = chart;
    bad.outer_top_z = f64::NAN;
    assert!(plane::disk_transfer(o, n, u, bad, [1., 0.], &mut ledger()).is_err());
    bad = chart;
    bad.columns[1][0] = f64::INFINITY;
    assert!(plane::disk_transfer(o, n, u, bad, [1., 0.], &mut ledger()).is_err());
    assert!(plane::disk_transfer(o, o, u, chart, [1., 0.], &mut ledger()).is_err());
    assert!(plane::disk_transfer(o, n.map(|x| -x), u, chart, [1., 0.], &mut ledger()).is_err());
    bad = chart;
    bad.columns[0][2] = f64::MAX;
    bad.outer_top_z = -f64::MAX;
    assert!(plane::disk_transfer(o, n, u, bad, [1., f64::MAX], &mut ledger()).is_err());
    let mut w = ledger();
    w.remaining = 0;
    assert!(plane::disk_transfer(o, n, u, chart, [1., 0.], &mut w).is_err());
}
#[test]
fn actual_shoulder_public_frontiers_cover_whole_stage_without_restarting_inherited_quota() {
    degree11::tests::with_head(&mut |head| {
        let field = head
            .common_chart_periodic_thread_repair(Default::default(), None)
            .unwrap();
        let finite = head
            .periodic_thread_finite_curves(&field, Default::default(), None)
            .unwrap();
        let o = PeriodicShoulderOptions::default();
        let good = head
            .periodic_shared_crest_shoulder(&finite, o, None)
            .unwrap();
        let bytes = good.retained_budget_bytes() + SCRATCH;
        for bad in [
            PeriodicShoulderOptions {
                maximum_distance_mm: 0.,
                ..o
            },
            PeriodicShoulderOptions {
                maximum_distance_mm: f64::NAN,
                ..o
            },
            PeriodicShoulderOptions {
                maximum_distance_mm: 1e-6,
                ..o
            },
            PeriodicShoulderOptions {
                maximum_cell_derivative_mm: 0.,
                ..o
            },
            PeriodicShoulderOptions {
                maximum_cell_derivative_mm: f64::NAN,
                ..o
            },
            PeriodicShoulderOptions {
                maximum_cell_derivative_mm: 1e-6,
                ..o
            },
            PeriodicShoulderOptions {
                maximum_work: WORK + 1,
                ..o
            },
            PeriodicShoulderOptions {
                maximum_total_work: TOTAL + 1,
                ..o
            },
            PeriodicShoulderOptions {
                maximum_bytes: BYTES + 1,
                ..o
            },
        ] {
            assert!(matches!(
                head.periodic_shared_crest_shoulder(&finite, bad, None),
                Err(PhillipsReadError::Invalid)
            ));
        }
        for bad in [
            PeriodicShoulderOptions {
                maximum_work: 0,
                ..o
            },
            PeriodicShoulderOptions {
                maximum_work: 1,
                ..o
            },
            PeriodicShoulderOptions {
                maximum_work: good.new_work_used() - 1,
                ..o
            },
            PeriodicShoulderOptions {
                maximum_total_work: finite.total_work_used() - 1,
                ..o
            },
            PeriodicShoulderOptions {
                maximum_total_work: finite.total_work_used(),
                ..o
            },
            PeriodicShoulderOptions {
                maximum_total_work: good.total_work_used() - 1,
                ..o
            },
            PeriodicShoulderOptions {
                maximum_bytes: bytes - 1,
                ..o
            },
        ] {
            assert!(matches!(
                head.periodic_shared_crest_shoulder(&finite, bad, None),
                Err(PhillipsReadError::Limit)
            ));
        }
        for bad in [
            PeriodicShoulderOptions {
                maximum_distance_mm: good.original_complement().distance_mm / 2.,
                ..o
            },
            PeriodicShoulderOptions {
                maximum_cell_derivative_mm: good.original_complement().cell_derivative_mm / 2.,
                ..o
            },
        ] {
            assert!(matches!(
                head.periodic_shared_crest_shoulder(&finite, bad, None),
                Err(PhillipsReadError::Unsupported)
            ));
        }
        assert!(
            head.periodic_shared_crest_shoulder(
                &finite,
                PeriodicShoulderOptions {
                    maximum_work: good.new_work_used(),
                    maximum_total_work: good.total_work_used(),
                    maximum_bytes: bytes,
                    ..o
                },
                None
            )
            .is_ok()
        );
        degree11::tests::with_head(&mut |other| {
            assert!(matches!(
                other.periodic_shared_crest_shoulder(&finite, o, None),
                Err(PhillipsReadError::Invalid)
            ));
        });
        let calls = Cell::new(0);
        let count = || {
            calls.set(calls.get() + 1);
            false
        };
        head.periodic_shared_crest_shoulder(&finite, o, Some(&count))
            .unwrap();
        let complete = calls.get();
        assert!(complete > 20);
        for stop in 1..=complete {
            calls.set(0);
            let cancel = || {
                calls.set(calls.get() + 1);
                calls.get() == stop
            };
            assert!(
                matches!(
                    head.periodic_shared_crest_shoulder(&finite, o, Some(&cancel)),
                    Err(PhillipsReadError::Cancelled)
                ),
                "callback{stop}/{complete}"
            );
        }
        assert!(
            head.periodic_shared_crest_shoulder(&finite, o, None)
                .is_ok()
        );
    });
}
#[test]
fn actual_off_nominal_annular_shoulder_is_an_explicit_separate_family() {
    let source = crate::step_read::phillips_cutter::physical::tests::fixture(0, false, false, 1);
    let text =
        crate::step::write_phillips_screw(&source, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7)
            .unwrap();
    crate::step_read::phillips_cutter::physical::body::common_mouth::tests::derived_thread_probe::with_head(&text,&mut |head| {
  let field=head.common_chart_periodic_thread_repair(Default::default(),None).unwrap();
  let finite=head.periodic_thread_finite_curves(&field,Default::default(),None).unwrap();
  assert_ne!(head.virtual_bottom_station()[0],field.source().radius_nodes_mm()[0]);
  assert!(matches!(head.periodic_shared_crest_shoulder(&finite,Default::default(),None),Err(PhillipsReadError::Unsupported)));
 });
}
#[test]
fn actual_five_nominal_metric_heads_and_posed_left_hand_join_matrix() {
    for index in [0, 1, 3, 4, 5] {
        for left in [false, true] {
            let source =
                crate::step_read::phillips_cutter::physical::tests::fixture(index, left, true, 1);
            let (o, u, v) = if left {
                ([105., -20., 60.], [0.8, 0.6, 0.], [0.36, -0.48, -0.8])
            } else {
                ([0.; 3], [1., 0., 0.], [0., 1., 0.])
            };
            let text = crate::step::write_phillips_screw(&source, o, u, v, 1e-7).unwrap();
            crate::step_read::phillips_cutter::physical::body::common_mouth::tests::derived_thread_probe::with_head(&text,&mut |head| {
   let field=head.common_chart_periodic_thread_repair(Default::default(),None).unwrap();
   let finite=head.periodic_thread_finite_curves(&field,Default::default(),None).unwrap();
   let join=head.periodic_shared_crest_shoulder(&finite,Default::default(),None);
   eprintln!("NOMINAL SHOULDER index{index} left{left} bottom{:?} threadR{} accepted{}",head.virtual_bottom_station(),field.source().radius_nodes_mm()[0],join.is_ok());
   let join=join.unwrap();assert!(join.original_complement().distance_mm<1e-7&&join.plane_distance_mm()<1e-7);
  });
        }
    }
}
#[test]
fn legal_thread_with_depth_lost_at_large_radius_does_not_earn_noncollapsed_annulus() {
    use crate::helical_thread::{Spec, Thread};
    let thread = Thread::new(Spec {
        major_radius_mm: 1e6,
        pitch_mm: 0.01,
        length_mm: 1.,
        depth_mm: 1e-12,
        crest_width_mm: 0.002,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts: 1,
        left_hand: false,
        root_radius_mm: None,
    })
    .unwrap();
    assert_eq!(thread.root_radius_mm(), thread.spec().major_radius_mm);
    let mut work = Work {
        remaining: WORK,
        interrupt: None,
        cancelled: Cell::new(false),
    };
    let source = PeriodicThreadSource::new(thread, &mut work).unwrap();
    assert!(matches!(
        domain::Domain::new(source, [1e6, 0.], [0., 1.], &mut work),
        Err(PhillipsReadError::Unsupported)
    ));
}
