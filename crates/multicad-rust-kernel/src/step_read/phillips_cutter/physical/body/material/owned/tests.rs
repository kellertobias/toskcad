use super::super::super::common_mouth::tests::derived_thread_probe::with_head;
use super::*;
fn text(scale: f64, posed: bool) -> String {
    super::super::super::super::branches::representative::scaled_text(scale, posed)
}
#[test]
fn actual_owned_sources_and_local_volume_survive_all_original_owners() {
    for (scale, posed) in [(1., false), (1.137, true)] {
        let (owned, samples, faces, source_world) = {
            let input = text(scale, posed);
            let mut result = None;
            with_head(&input, &mut |head| {
                let thread = head
                    .derived_thread_carrier(Default::default(), None)
                    .unwrap()
                    .certify_occupied(1e-7, 1e-5, None)
                    .unwrap();
                let source = head
                    .closed_repaired_material(&thread, Default::default(), None)
                    .unwrap();
                let samples = (0..source.head().segment_count())
                    .map(|i| {
                        source
                            .head()
                            .evaluate(i, [0.125, 0.875])
                            .unwrap()
                            .position_mm
                    })
                    .collect::<Vec<_>>();
                let faces = (0..source.face_count())
                    .map(|i| source.face_source(i))
                    .collect::<Vec<_>>();
                let value = source.volume_bounds(Default::default(), None).unwrap();
                let world = value.world_volume_mm3();
                let local = value.local_volume_mm3();
                let components = value.component_volumes_mm3();
                let chart = source.chart();
                let determinant = source.determinant_bounds();
                let join = source.join_radial_margin_mm();
                let clearances = [
                    source.containment().minimum_radial_clearance_mm(),
                    source.containment().axial_floor_clearance_mm(),
                ];
                let owned = value.into_owned(Default::default(), None).unwrap();
                assert_eq!(owned.original_input(), input);
                assert_eq!(owned.local_volume_mm3(), local);
                assert_eq!(owned.component_volumes_mm3(), components);
                assert_eq!(owned.chart(), chart);
                assert_eq!(owned.determinant_bounds(), determinant);
                assert_eq!(owned.join_radial_margin_mm(), join);
                assert_eq!(owned.containment_clearances_mm(), clearances);
                assert_eq!(owned.join_inner_coedges(), source.join_inner_coedges());
                assert_eq!(owned.cavity().spec(), source.cavity().spec());
                assert_eq!(owned.cavity().domain(), source.cavity().domain());
                assert_eq!(owned.thread().input(), source.thread().input());
                assert!(owned.work_used() > source.work_used());
                assert!(owned.retained_budget_bytes() < BYTES);
                result = Some((owned, samples, faces, world));
            });
            result.unwrap()
        };
        // Parser AST, input String, all source/volume receipts are now gone.
        let head = owned.head().view();
        for (i, expected) in samples.into_iter().enumerate() {
            assert_eq!(
                head.evaluate(i, [0.125, 0.875]).unwrap().position_mm,
                expected
            );
            assert!(head.segment_certificate(i).is_some());
        }
        assert_eq!(owned.face_count(), faces.len());
        for (i, expected) in faces.into_iter().enumerate() {
            assert_eq!(owned.face_source(i), expected);
        }
        assert_eq!(owned.face_source(owned.face_count()), None);
        assert_eq!(owned.world_volume_mm3(), source_world);
        let thread = owned
            .thread()
            .volume_with_interrupt(1e-5, None)
            .unwrap()
            .volume_mm3;
        let cavity = owned
            .cavity()
            .volume_bounds_mm3(1e-5, 100_000)
            .unwrap()
            .volume_mm3;
        let components = owned.component_volumes_mm3();
        assert!(thread[0] <= components[1][0] && thread[1] >= components[1][1]);
        assert!(cavity[0] <= components[2][0] && cavity[1] >= components[2][1]);

        assert!(owned.original_input().contains("ISO-10303-21"));
    }
}
#[test]
fn extraction_preflights_resources_and_cancels_without_mutating_original() {
    let input = text(1., false);
    with_head(&input, &mut |head| {
        let thread = head
            .derived_thread_carrier(Default::default(), None)
            .unwrap()
            .certify_occupied(1e-7, 1e-5, None)
            .unwrap();
        let source = head
            .closed_repaired_material(&thread, Default::default(), None)
            .unwrap();
        let run = |options, cancel: Option<&dyn Fn() -> bool>| {
            source
                .volume_bounds(Default::default(), None)
                .unwrap()
                .into_owned(options, cancel)
        };
        let owned = run(Default::default(), None).unwrap();
        let exact_peak = source
            .volume_bounds(Default::default(), None)
            .unwrap()
            .owned_peak_budget_bytes();
        assert!(
            run(
                OwnedMaterialOptions {
                    maximum_bytes: exact_peak,
                    ..Default::default()
                },
                None
            )
            .is_ok()
        );
        assert!(matches!(
            run(
                OwnedMaterialOptions {
                    maximum_bytes: exact_peak - 1,
                    ..Default::default()
                },
                None
            ),
            Err(PhillipsReadError::Limit)
        ));

        for options in [
            OwnedMaterialOptions {
                maximum_work: 512_000_001,
                ..Default::default()
            },
            OwnedMaterialOptions {
                maximum_bytes: BYTES + 1,
                ..Default::default()
            },
        ] {
            assert!(matches!(
                run(options, None),
                Err(PhillipsReadError::Invalid)
            ));
        }
        for options in [
            OwnedMaterialOptions {
                maximum_work: 0,
                ..Default::default()
            },
            OwnedMaterialOptions {
                maximum_bytes: 0,
                ..Default::default()
            },
            OwnedMaterialOptions {
                maximum_work: owned.work_used() - 1,
                ..Default::default()
            },
        ] {
            assert!(matches!(run(options, None), Err(PhillipsReadError::Limit)));
        }
        assert!(
            run(
                OwnedMaterialOptions {
                    maximum_work: owned.work_used(),
                    ..Default::default()
                },
                None
            )
            .is_ok()
        );
        for cutoff in 0..5 {
            let n = Cell::new(0);
            let cancel = || {
                let k = n.get();
                n.set(k + 1);
                k >= cutoff
            };
            assert!(matches!(
                run(Default::default(), Some(&cancel)),
                Err(PhillipsReadError::Cancelled)
            ));
        }
        // Exercise the actual-capacity publication guard with a genuine
        // independently extracted candidate and a lower destination budget.
        let candidate = run(Default::default(), None).unwrap();
        let options = OwnedMaterialOptions {
            maximum_bytes: candidate.retained_budget_bytes() - 1,
            ..Default::default()
        };
        let mut work = Work {
            remaining: options.maximum_work,
            interrupt: None,
            cancelled: Cell::new(false),
        };
        assert!(matches!(
            publish(candidate, 0, options, &mut work),
            Err(PhillipsReadError::Limit)
        ));
        let retry = run(Default::default(), Some(&|| false)).unwrap();
        assert_eq!(retry.local_volume_mm3(), owned.local_volume_mm3());
        assert_eq!(source.body().literal_input(), input);
    });
}
