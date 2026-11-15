use super::*;
#[test]
fn direct_authored_metric_factory_does_not_need_a_step_reader_or_decoded_owner() {
    let source = crate::step_read::phillips_cutter::physical::tests::fixture(0, false, true, 1);
    assert!(OwnedPeriodicMaterial::supports_authored_family(
        source.thread(),
        source.profile(),
        source.spec()
    ));
    let start = std::time::Instant::now();
    let owned = OwnedPeriodicMaterial::from_authored(&source, Default::default(), None).unwrap();
    assert!(matches!(
        owned.provenance(),
        PeriodicMaterialProvenance::Authored { .. }
    ));
    assert!(owned.original_input().is_none());
    assert!(owned.original_curve_owners().is_none());
    assert!(owned.thread_transfers().is_none());
    assert_eq!(owned.head_z_shift_mm(), 0.);
    assert_eq!(owned.inherited_work_used(), 0);
    assert!(owned.local_volume_mm3()[0] > 0.);
    assert_eq!(owned.determinant_bounds(), [1., 1.]);
    eprintln!(
        "AUTHORED PERIODIC {:?} bytes{} work{} volume{:?} caps{:?}",
        start.elapsed(),
        owned.retained_budget_bytes(),
        owned.total_work_used(),
        owned.world_volume_mm3(),
        owned.cap_expressions_mm()
    );
}
fn compare_original_source(source: &PhillipsScrew, label: &str) {
    let expected = source.source_volume_bounds().unwrap();
    let start = std::time::Instant::now();
    let actual = OwnedPeriodicMaterial::from_authored(source, Default::default(), None).unwrap();
    let volume = actual.world_volume_mm3();
    let oracle = expected.volume_mm3();
    assert!(
        volume[0] <= oracle[1] && volume[1] >= oracle[0],
        "{label}: new{volume:?} original{oracle:?}"
    );
    assert_eq!(actual.head_z_shift_mm(), source.spec().shank_length_mm);
    assert_eq!(
        actual.join_station_expression_mm(),
        [source.profile().end()[1], -source.spec().shank_length_mm]
    );
    assert_eq!(
        actual.cap_expressions_mm(),
        [
            [source.profile().end()[1], -source.spec().shank_length_mm],
            [source.thread().spec().length_mm, 0.]
        ]
    );
    assert_eq!(
        actual.cavity_mouth_expression_mm(),
        [source.profile().top()[1], -source.spec().shank_length_mm]
    );
    eprintln!(
        "AUTHORED ORACLE {label} elapsed{:?} thread{:?} head{:?} plain{} cutter{:?} new{volume:?} original{oracle:?} headShift{} caps{:?}",
        start.elapsed(),
        source.thread().spec(),
        source.spec().head,
        source.spec().shank_length_mm,
        source.cavity().source().spec(),
        actual.head_z_shift_mm(),
        actual.cap_expressions_mm()
    );
}
#[test]
fn identical_actual_normal_m2_and_m3_include_nonzero_plain_shank_in_independent_oracle() {
    let normal = degree11::tests::fixture::actual_source();
    assert_eq!(normal.thread().spec().length_mm, 11.5);
    assert_eq!(normal.spec().shank_length_mm, 0.5);
    compare_original_source(&normal, "actual normal ISO7045 M2 length11.5 plain0.5");
    for index in [0, 1] {
        let base =
            crate::step_read::phillips_cutter::physical::tests::fixture(index, false, true, 1);
        compare_original_source(
            &base,
            if index == 0 {
                "short M2 length2 plain0"
            } else {
                "short M3 length2 plain0"
            },
        );
        for plain in [0.03, 0.3] {
            let mut spec = base.spec();
            spec.shank_length_mm = plain;
            let profile = Profile::new_for_thread_union(
                spec.head.unwrap(),
                base.thread().spec().major_radius_mm,
                plain,
            )
            .unwrap();
            let source = PhillipsScrew::new(
                *base.thread(),
                profile,
                spec,
                base.cavity().clone(),
                100_000,
                &|| false,
            )
            .unwrap();
            compare_original_source(
                &source,
                &format!("metric index{index} length2 plain{plain}"),
            );
        }
    }
}
#[test]
fn direct_actual_metric_matrix_retains_hand_and_no_decoded_provenance() {
    for index in [0, 1, 3, 4, 5] {
        for left in [false, true] {
            let source =
                crate::step_read::phillips_cutter::physical::tests::fixture(index, left, true, 1);
            let owned = OwnedPeriodicMaterial::from_authored_parts(
                source.thread(),
                source.profile(),
                source.cavity().source(),
                source.spec(),
                Default::default(),
                None,
            )
            .unwrap();
            assert_eq!(owned.thread().left_hand(), left);
            assert!(owned.original_input().is_none());
            assert!(owned.original_curve_owners().is_none());
            assert!(owned.thread_transfers().is_none());
            assert_eq!(
                owned.head().view().station(0),
                source.profile().segments().next().map(|s| s.start)
            );
            let volume = owned.world_volume_mm3();
            let oracle = source.source_volume_bounds().unwrap().volume_mm3();
            assert!(volume[0] <= oracle[1] && volume[1] >= oracle[0]);
        }
    }
}
#[test]
fn direct_authored_parts_validate_requested_geometry_and_operation_frontiers() {
    let source = crate::step_read::phillips_cutter::physical::tests::fixture(0, false, true, 1);
    let run = |o, interrupt: Option<&dyn Fn() -> bool>| {
        OwnedPeriodicMaterial::from_authored(&source, o, interrupt)
    };
    let o = AuthoredPeriodicOptions::default();
    let good = run(o, None).unwrap();
    for bad in [
        AuthoredPeriodicOptions {
            minimum_clearance_mm: 0.,
            ..o
        },
        AuthoredPeriodicOptions {
            minimum_clearance_mm: f64::NAN,
            ..o
        },
        AuthoredPeriodicOptions {
            maximum_width_mm3: 0.,
            ..o
        },
        AuthoredPeriodicOptions {
            maximum_width_mm3: f64::NAN,
            ..o
        },
        AuthoredPeriodicOptions {
            maximum_containment_cells: 0,
            ..o
        },
        AuthoredPeriodicOptions {
            maximum_containment_cells: 65537,
            ..o
        },
        AuthoredPeriodicOptions {
            maximum_head_cells: 0,
            ..o
        },
        AuthoredPeriodicOptions {
            maximum_head_cells: 65537,
            ..o
        },
        AuthoredPeriodicOptions {
            maximum_cavity_evaluations: 0,
            ..o
        },
        AuthoredPeriodicOptions {
            maximum_cavity_evaluations: 1_000_001,
            ..o
        },
        AuthoredPeriodicOptions {
            maximum_work: TOTAL + 1,
            ..o
        },
        AuthoredPeriodicOptions {
            maximum_bytes: BYTES + 1,
            ..o
        },
    ] {
        assert!(matches!(run(bad, None), Err(PhillipsReadError::Invalid)));
    }
    for bad in [
        AuthoredPeriodicOptions {
            maximum_work: 0,
            ..o
        },
        AuthoredPeriodicOptions {
            maximum_work: 32768,
            ..o
        },
        AuthoredPeriodicOptions {
            maximum_work: 32769,
            ..o
        },
        AuthoredPeriodicOptions {
            maximum_work: good.new_work_used() - 1,
            ..o
        },
        AuthoredPeriodicOptions {
            maximum_bytes: FIXED + ClosedBoundaryCertificate::scratch_budget_bytes() - 1,
            ..o
        },
        AuthoredPeriodicOptions {
            maximum_containment_cells: 1,
            ..o
        },
        AuthoredPeriodicOptions {
            maximum_head_cells: 1,
            ..o
        },
        AuthoredPeriodicOptions {
            maximum_cavity_evaluations: 1,
            ..o
        },
    ] {
        assert!(matches!(run(bad, None), Err(PhillipsReadError::Limit)));
    }
    assert!(
        run(
            AuthoredPeriodicOptions {
                maximum_work: good.new_work_used(),
                ..o
            },
            None
        )
        .is_ok()
    );
    assert!(matches!(
        run(
            AuthoredPeriodicOptions {
                minimum_clearance_mm: 1000.,
                ..o
            },
            None
        ),
        Err(PhillipsReadError::Unsupported)
    ));
    assert!(matches!(
        run(
            AuthoredPeriodicOptions {
                maximum_width_mm3: 1e-30,
                ..o
            },
            None
        ),
        Err(PhillipsReadError::Unsupported)
    ));
    let mut spec = source.spec();
    spec.shank_length_mm = -1.;
    assert!(matches!(
        OwnedPeriodicMaterial::from_authored_parts(
            source.thread(),
            source.profile(),
            source.cavity().source(),
            spec,
            o,
            None
        ),
        Err(PhillipsReadError::Invalid)
    ));
    spec = source.spec();
    spec.minimum_wall_mm = f64::NAN;
    assert!(matches!(
        OwnedPeriodicMaterial::from_authored_parts(
            source.thread(),
            source.profile(),
            source.cavity().source(),
            spec,
            o,
            None
        ),
        Err(PhillipsReadError::Invalid)
    ));
    let mut thread = source.thread().spec();
    thread.starts = 2;
    let multi = Thread::new(thread).unwrap();
    assert!(matches!(
        OwnedPeriodicMaterial::from_authored_parts(
            &multi,
            source.profile(),
            source.cavity().source(),
            source.spec(),
            o,
            None
        ),
        Err(PhillipsReadError::Unsupported)
    ));
    thread.starts = 1;
    thread.major_radius_mm = 2.;
    let wide = Thread::new(thread).unwrap();
    assert!(matches!(
        OwnedPeriodicMaterial::from_authored_parts(
            &wide,
            source.profile(),
            source.cavity().source(),
            source.spec(),
            o,
            None
        ),
        Err(PhillipsReadError::Unsupported)
    ));
    let annular = crate::step_read::phillips_cutter::physical::tests::fixture(0, false, false, 1);
    assert!(matches!(
        OwnedPeriodicMaterial::from_authored(&annular, o, None),
        Err(PhillipsReadError::Unsupported)
    ));
    let countersunk =
        crate::step_read::phillips_cutter::physical::tests::fixture(2, false, true, 1);
    assert!(!OwnedPeriodicMaterial::supports_authored_family(
        countersunk.thread(),
        countersunk.profile(),
        countersunk.spec()
    ));
    assert!(!OwnedPeriodicMaterial::supports_authored_family(
        annular.thread(),
        annular.profile(),
        annular.spec()
    ));
    assert!(!OwnedPeriodicMaterial::supports_authored_family(
        &multi,
        source.profile(),
        source.spec()
    ));
    assert!(matches!(
        OwnedPeriodicMaterial::from_authored(&countersunk, o, None),
        Err(PhillipsReadError::Unsupported)
    ));
    // Count real public component callback paths, then stop at each handoff. The
    // constructors remain unchanged and each numerical child retains its cap.
    let calls = Cell::new(0);
    let count = || {
        calls.set(calls.get() + 1);
        false
    };
    let head_proof = source
        .profile()
        .closed_repaired_meridian(RepairedHeadOptions::default(), Some(&count))
        .unwrap();
    let head_calls = calls.get();
    calls.set(0);
    let local_proof = source
        .cavity()
        .source()
        .closed_boundary_certificate(Default::default(), Some(&count))
        .unwrap();
    let cavity_calls = calls.get();
    calls.set(0);
    run(o, Some(&count)).unwrap();
    let complete = calls.get();
    calls.set(0);
    let containment = local_proof
        .certify_repaired_head(
            &head_proof,
            RepairedHeadContainmentOptions {
                minimum_clearance_mm: o.minimum_clearance_mm.max(source.spec().minimum_wall_mm),
                ..Default::default()
            },
            Some(&count),
        )
        .unwrap();
    let containment_calls = calls.get();
    let paid = head_proof.work_used() + local_proof.work_used();
    let before_mass_work = 1
        + 32768
        + head_proof.work_used()
        + local_proof.work_used()
        + 64
        + (containment.work_used() - paid);
    assert!(matches!(
        run(
            AuthoredPeriodicOptions {
                maximum_work: before_mass_work,
                ..o
            },
            None
        ),
        Err(PhillipsReadError::Limit)
    ));
    let before_head = 2;
    let after_head = before_head + head_calls;
    let after_cavity = after_head + 1 + cavity_calls;
    let after_containment = after_cavity + 2 + containment_calls;
    for stop in (1..=4).chain([
        after_containment + 1,
        after_containment + 2,
        after_containment + 3,
        after_head,
        after_head + 1,
        after_cavity,
        after_cavity + 1,
        after_cavity + 2,
        complete - 2,
        complete - 1,
        complete,
    ]) {
        calls.set(0);
        let cancel = || {
            calls.set(calls.get() + 1);
            calls.get() == stop
        };
        assert!(
            matches!(run(o, Some(&cancel)), Err(PhillipsReadError::Cancelled)),
            "callback{stop}/{complete}"
        );
    }
    assert!(run(o, None).is_ok());
}
