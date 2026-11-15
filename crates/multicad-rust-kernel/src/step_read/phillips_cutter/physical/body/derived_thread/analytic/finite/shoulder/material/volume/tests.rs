use super::*;
#[test]
fn actual_m2_periodic_repair_source_volume_composes_actual_components() {
    degree11::tests::with_head(&mut |head| {
        let field = head
            .common_chart_periodic_thread_repair(Default::default(), None)
            .unwrap();
        let finite = head
            .periodic_thread_finite_curves(&field, Default::default(), None)
            .unwrap();
        let join = head
            .periodic_shared_crest_shoulder(&finite, Default::default(), None)
            .unwrap();
        let material = head
            .periodic_closed_material(&join, Default::default(), None)
            .unwrap();
        let start = std::time::Instant::now();
        let volume = material.volume_bounds(Default::default(), None).unwrap();
        assert!(std::ptr::eq(volume.source(), &material));
        let local = volume.local_volume_mm3();
        let world = volume.world_volume_mm3();
        assert!(local[0] > 0. && local[1] - local[0] <= 1e-6 && world[1] - world[0] <= 1e-6);
        assert_eq!(volume.inherited_work_used(), material.total_work_used());
        assert_eq!(
            volume.total_work_used(),
            volume.inherited_work_used() + volume.new_work_used()
        );
        assert!(volume.retained_budget_bytes() < BYTES);
        assert!(
            volume.head_cells() > 0
                && volume.cavity_evaluations() > 0
                && volume.head_work_used() > 0
                && volume.thread_work_used() > 0
        );
        eprintln!(
            "PERIODIC VOLUME {:?} local{local:?} world{world:?} components{:?} cells{} evaluations{} headWork{} threadWork{} new{} total{}",
            start.elapsed(),
            volume.component_volumes_mm3(),
            volume.head_cells(),
            volume.cavity_evaluations(),
            volume.head_work_used(),
            volume.thread_work_used(),
            volume.new_work_used(),
            volume.total_work_used()
        );
    });
}
#[test]
fn periodic_section_measure_matches_independent_exact_simpson_and_machin_oracle() {
    use num_bigint::BigInt as Z;
    use num_rational::BigRational as Q;
    let q = |x: f64| Q::from_float(x).unwrap();
    let atan = |d: i64| {
        let mut sum = Q::from_integer(Z::from(0));
        for k in 0..32_u32 {
            let term = Q::new(Z::from(1), Z::from(2 * k + 1) * Z::from(d).pow(2 * k + 1));
            if k % 2 == 0 {
                sum += term;
            } else {
                sum -= term;
            }
        }
        let upper = &sum + Q::new(Z::from(1), Z::from(65) * Z::from(d).pow(65));
        [sum, upper]
    };
    let a = atan(5);
    let b = atan(239);
    let pi = [
        &a[0] * Z::from(16) - &b[1] * Z::from(4),
        &a[1] * Z::from(16) - &b[0] * Z::from(4),
    ];
    for left in [false, true] {
        let thread = crate::helical_thread::Thread::new(crate::helical_thread::Spec {
            major_radius_mm: 1.,
            pitch_mm: 1.,
            length_mm: 1.,
            depth_mm: 0.2,
            crest_width_mm: 0.25,
            lower_flank_deg: 30.,
            upper_flank_deg: 30.,
            starts: 1,
            left_hand: left,
            root_radius_mm: None,
        })
        .unwrap();
        let mut work = Work {
            remaining: WORK,
            interrupt: None,
            cancelled: Cell::new(false),
        };
        let source = PeriodicThreadSource::new(thread, &mut work).unwrap();
        let phase = source.phase_node_expressions_mm().map(|[a, b]| q(a) + q(b));
        let radius = source.radius_nodes_mm().map(q);
        let mut integral = Q::from_integer(Z::from(0));
        for i in 0..4 {
            let mid = (&radius[i] + &radius[i + 1]) / Z::from(2);
            integral += (&phase[i + 1] - &phase[i])
                * (&radius[i] * &radius[i]
                    + &mid * &mid * Z::from(4)
                    + &radius[i + 1] * &radius[i + 1])
                / Z::from(6);
        }
        let stations = [0.125, 0.875];
        let exact_factor = integral / q(source.pitch_mm()) * (q(stations[1]) - q(stations[0]));
        let actual = periodic_volume(&source, stations, &mut work).unwrap();
        assert!(q(actual.lo) <= &pi[0] * &exact_factor && q(actual.hi) >= &pi[1] * &exact_factor);
        for limit in [0, 1024, 5120, 9216, 9472, 9728, 9984] {
            let mut w = Work {
                remaining: limit,
                interrupt: None,
                cancelled: Cell::new(false),
            };
            assert!(matches!(
                periodic_volume(&source, stations, &mut w),
                Err(PhillipsReadError::Limit)
            ));
        }
    }
}
#[test]
fn actual_periodic_volume_resource_and_cancellation_frontiers_are_whole_stage() {
    degree11::tests::with_head(&mut |head| {
        let field = head
            .common_chart_periodic_thread_repair(Default::default(), None)
            .unwrap();
        let finite = head
            .periodic_thread_finite_curves(&field, Default::default(), None)
            .unwrap();
        let join = head
            .periodic_shared_crest_shoulder(&finite, Default::default(), None)
            .unwrap();
        let material = head
            .periodic_closed_material(&join, Default::default(), None)
            .unwrap();
        let o = PeriodicVolumeOptions::default();
        let good = material.volume_bounds(o, None).unwrap();
        let minimum_bytes = good.retained_budget_bytes()
            + BASE
            + STREAM
                .max(RepairedHeadVolume::scratch_budget_bytes(good.head_cells()))
                .max(Cutter::volume_scratch_budget_bytes(
                    good.cavity_evaluations(),
                ));
        for bad in [
            PeriodicVolumeOptions {
                maximum_width_mm3: 0.,
                ..o
            },
            PeriodicVolumeOptions {
                maximum_width_mm3: f64::NAN,
                ..o
            },
            PeriodicVolumeOptions {
                maximum_head_cells: 0,
                ..o
            },
            PeriodicVolumeOptions {
                maximum_head_cells: 65537,
                ..o
            },
            PeriodicVolumeOptions {
                maximum_cavity_evaluations: 0,
                ..o
            },
            PeriodicVolumeOptions {
                maximum_cavity_evaluations: 1_000_001,
                ..o
            },
            PeriodicVolumeOptions {
                maximum_work: WORK + 1,
                ..o
            },
            PeriodicVolumeOptions {
                maximum_total_work: TOTAL + 1,
                ..o
            },
            PeriodicVolumeOptions {
                maximum_bytes: BYTES + 1,
                ..o
            },
        ] {
            assert!(matches!(
                material.volume_bounds(bad, None),
                Err(PhillipsReadError::Invalid)
            ));
        }
        for bad in [
            PeriodicVolumeOptions {
                maximum_work: 0,
                ..o
            },
            PeriodicVolumeOptions {
                maximum_work: 1,
                ..o
            },
            PeriodicVolumeOptions {
                maximum_work: META,
                ..o
            },
            PeriodicVolumeOptions {
                maximum_work: META + 1,
                ..o
            },
            PeriodicVolumeOptions {
                maximum_work: 1 + META + good.thread_work_used() + good.head_work_used(),
                ..o
            },
            PeriodicVolumeOptions {
                maximum_work: good.new_work_used() - 1,
                ..o
            },
            PeriodicVolumeOptions {
                maximum_total_work: material.total_work_used() - 1,
                ..o
            },
            PeriodicVolumeOptions {
                maximum_total_work: material.total_work_used(),
                ..o
            },
            PeriodicVolumeOptions {
                maximum_total_work: good.total_work_used() - 1,
                ..o
            },
            PeriodicVolumeOptions {
                maximum_bytes: good.retained_budget_bytes() + BASE + STREAM - 1,
                ..o
            },
            PeriodicVolumeOptions {
                maximum_bytes: minimum_bytes - 1,
                ..o
            },
            PeriodicVolumeOptions {
                maximum_head_cells: 1,
                ..o
            },
            PeriodicVolumeOptions {
                maximum_cavity_evaluations: 1,
                ..o
            },
        ] {
            assert!(matches!(
                material.volume_bounds(bad, None),
                Err(PhillipsReadError::Limit)
            ));
        }
        assert!(matches!(
            material.volume_bounds(
                PeriodicVolumeOptions {
                    maximum_width_mm3: 1e-30,
                    ..o
                },
                None
            ),
            Err(PhillipsReadError::Unsupported)
        ));
        assert!(
            material
                .volume_bounds(
                    PeriodicVolumeOptions {
                        maximum_work: good.new_work_used(),
                        maximum_total_work: good.total_work_used(),
                        maximum_bytes: minimum_bytes,
                        ..o
                    },
                    None
                )
                .is_ok()
        );
        // Calibrate child frontiers through their real public queries. This is only
        // callback accounting; the independent exact-Q test is the area oracle.
        let thread = interval(good.component_volumes_mm3()[1]);
        let determinant = interval(material.determinant_bounds());
        let source = head.source();
        let mut radius = 0_f64;
        for i in 0..source.segment_count() {
            radius = radius.max(source.certified_interval_bounds(i, [0., 1.]).position_mm[0][1]);
        }
        let pi = I::new(
            std::f64::consts::PI.next_down(),
            std::f64::consts::PI.next_up(),
        );
        let outer = pi
            * I::point(radius).square()
            * (I::point(source.station(source.segment_count()).unwrap()[1])
                - I::point(source.station(0).unwrap()[1]))
            + I::point(thread.hi);
        let reserve = ((I::point(determinant.hi) - I::point(determinant.lo)) * outer).hi;
        let width = (I::point((I::point(o.maximum_width_mm3) - I::point(reserve)).lo)
            / (I::point(8.) * I::point(determinant.hi)))
        .lo;
        let calls = Cell::new(0);
        let count = || {
            calls.set(calls.get() + 1);
            false
        };
        source
            .volume_bounds(
                RepairedHeadVolumeOptions {
                    maximum_width_mm3: width,
                    ..Default::default()
                },
                Some(&count),
            )
            .unwrap();
        let head_calls = calls.get();
        calls.set(0);
        material
            .containment()
            .source()
            .volume_bounds_with_cancel(width, o.maximum_cavity_evaluations, &count)
            .unwrap();
        let cavity_calls = calls.get();
        calls.set(0);
        material.volume_bounds(o, Some(&count)).unwrap();
        let complete = calls.get();
        let before_head = 9;
        assert_eq!(complete, before_head + head_calls + 1 + cavity_calls + 2);
        let stops = (1..=10).chain([
            before_head + head_calls,
            before_head + head_calls + 1,
            before_head + head_calls + 2,
            complete - 2,
            complete - 1,
            complete,
        ]);
        for stop in stops {
            calls.set(0);
            let cancel = || {
                calls.set(calls.get() + 1);
                calls.get() == stop
            };
            assert!(
                matches!(
                    material.volume_bounds(o, Some(&cancel)),
                    Err(PhillipsReadError::Cancelled)
                ),
                "callback{stop}/{complete}"
            );
        }
        assert!(material.volume_bounds(o, None).is_ok());
    });
}
#[test]
fn actual_large_head_preserves_affine_precision_reserve_in_source_mass() {
    use crate::{
        head_profile::{Head, Profile},
        phillips_drive::{Cavity, Cutter, PhillipsScrew},
        profiled_screw::SourceSpec,
    };
    let (s, d) = crate::phillips_drive::test_fixture(0);
    let cavity = Cavity::new(Cutter::new(s, d).unwrap(), 1e-5, 100_000).unwrap();
    let head = Head::Crowned {
        radius_mm: 2000.,
        height_mm: 1600.,
        crown_radius_mm: 3200.,
        top_land_radius_mm: 1250.,
        neck_radius_mm: 0.1,
    };
    let base = crate::step_read::phillips_cutter::physical::tests::fixture(0, false, true, 1);
    let source = PhillipsScrew::new(
        *base.thread(),
        Profile::new_for_thread_union(head, 1., 0.).unwrap(),
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
    let text =
        crate::step::write_phillips_screw(&source, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7)
            .unwrap();
    crate::step_read::phillips_cutter::physical::body::common_mouth::tests::derived_thread_probe::with_head(&text,&mut |head| {
  let field=head.common_chart_periodic_thread_repair(Default::default(),None).unwrap();let finite=head.periodic_thread_finite_curves(&field,Default::default(),None).unwrap();let join=head.periodic_shared_crest_shoulder(&finite,Default::default(),None).unwrap();let material=head.periodic_closed_material(&join,Default::default(),None).unwrap();
  assert!(matches!(material.volume_bounds(Default::default(),None),Err(PhillipsReadError::Unsupported)));
 });
}
