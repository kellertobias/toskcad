use super::super::super::super::tests::actual_m3_parts;
use super::*;
use crate::step_read::phillips_cutter::physical::body::common_mouth::tests::derived_thread_probe::with_head;
#[test]
fn actual_fitted_m3_supplied_author_matches_every_closed_material_face() {
    let (author, frame) = actual_m3_parts(true, false);
    let raw =
        crate::step::write_phillips_screw(&author, frame[0], frame[1], frame[2], 1e-7).unwrap();
    with_head(&raw, &mut |head| {
        let occupied = head
            .degree11_thread_carrier(Default::default(), None)
            .unwrap()
            .certify_occupied(Default::default(), None)
            .unwrap();
        let material = head
            .closed_degree11_material(&occupied, Default::default(), None)
            .unwrap();
        let volume = material.volume_bounds(Default::default(), None).unwrap();
        let a = author.cavity().source().surface_patches().unwrap();
        let b = material.cavity().surface_patches().unwrap();
        let matched = volume
            .match_supplied_author(&author, frame, Default::default(), None)
            .unwrap();
        let accuracy = matched.cavity_accuracy_mm();
        let direct = a
            .source_correspondence(
                &b,
                crate::phillips_drive::PatchCorrespondenceOptions {
                    maximum_position_mm: accuracy[0],
                    maximum_derivative_mm: accuracy[1],
                    maximum_cells: 262144,
                },
                None,
            )
            .unwrap();
        assert_eq!(matched.cavity_cells_used(), direct.cells_used);
        assert!(std::ptr::eq(matched.author(), &author));
        assert!(std::ptr::eq(matched.decoded_volume(), &volume));
        assert!(std::ptr::eq(
            matched.supplied_thread().body(),
            material.body()
        ));
        assert_eq!(matched.supplied_thread().entry_rim_transfers().len(), 4);
        assert!(matched.maximum_world_distance_mm() <= 1e-7);
        assert!(matched.maximum_world_chart_derivative_mm() <= 1e-7);
        assert!(
            matched
                .independent_volume_difference_mm3()
                .into_iter()
                .all(|x| x.abs() <= 1e-5)
        );
        assert!(matched.total_work_used() <= WORK_LIMIT);
        assert!(matched.accounted_peak_bytes() <= MAX_BYTES);
        eprintln!(
            "D11_SUPPLIED_MATERIAL distance={} derivative={} mass={:?} new={} total={} cells={} retained={} peak={} proposal={:?}",
            matched.maximum_world_distance_mm(),
            matched.maximum_world_chart_derivative_mm(),
            matched.independent_volume_difference_mm3(),
            matched.new_work_used(),
            matched.total_work_used(),
            matched.cavity_cells_used(),
            matched.retained_budget_bytes(),
            matched.accounted_peak_bytes(),
            matched.original_head_proposal_bounds_mm()
        );
        let options = SuppliedDegree11MaterialOptions::default();
        let exact = matched.total_work_used();
        let exact_cells = matched.cavity_cells_used();
        assert!(
            volume
                .match_supplied_author(
                    &author,
                    frame,
                    SuppliedDegree11MaterialOptions {
                        maximum_work: exact,
                        maximum_cavity_cells: exact_cells,
                        ..options
                    },
                    None
                )
                .is_ok()
        );
        for o in [
            SuppliedDegree11MaterialOptions {
                maximum_work: exact - 1,
                ..options
            },
            SuppliedDegree11MaterialOptions {
                maximum_work: volume.work_used() - 1,
                ..options
            },
            SuppliedDegree11MaterialOptions {
                maximum_work: volume.work_used(),
                ..options
            },
            // Interior real workload frontiers: author head construction,
            // proposal transfer, and the first paid paired Taylor attempt.
            SuppliedDegree11MaterialOptions {
                maximum_work: 194_344_000,
                ..options
            },
            SuppliedDegree11MaterialOptions {
                maximum_work: 194_350_000,
                ..options
            },
            SuppliedDegree11MaterialOptions {
                maximum_work: 222_800_000,
                ..options
            },
            SuppliedDegree11MaterialOptions {
                maximum_bytes: 0,
                ..options
            },
            SuppliedDegree11MaterialOptions {
                maximum_cavity_cells: exact_cells - 1,
                ..options
            },
        ] {
            assert!(matches!(
                volume.match_supplied_author(&author, frame, o, None),
                Err(PhillipsReadError::Limit)
            ));
        }
        for o in [
            SuppliedDegree11MaterialOptions {
                maximum_distance_mm: 0.,
                ..options
            },
            SuppliedDegree11MaterialOptions {
                maximum_derivative_mm: f64::NAN,
                ..options
            },
            SuppliedDegree11MaterialOptions {
                maximum_distance_mm: 1.01e-7,
                ..options
            },
            SuppliedDegree11MaterialOptions {
                maximum_volume_discrepancy_mm3: 0.,
                ..options
            },
            SuppliedDegree11MaterialOptions {
                maximum_volume_discrepancy_mm3: 1.01e-5,
                ..options
            },
            SuppliedDegree11MaterialOptions {
                maximum_work: WORK_LIMIT + 1,
                ..options
            },
            SuppliedDegree11MaterialOptions {
                maximum_bytes: MAX_BYTES + 1,
                ..options
            },
            SuppliedDegree11MaterialOptions {
                maximum_cavity_cells: 0,
                ..options
            },
            SuppliedDegree11MaterialOptions {
                maximum_cavity_cells: 262145,
                ..options
            },
        ] {
            assert!(matches!(
                volume.match_supplied_author(&author, frame, o, None),
                Err(PhillipsReadError::Invalid)
            ));
        }
        for o in [
            SuppliedDegree11MaterialOptions {
                maximum_distance_mm: matched.maximum_world_distance_mm() / 2.,
                ..options
            },
            SuppliedDegree11MaterialOptions {
                maximum_derivative_mm: matched.maximum_world_chart_derivative_mm() / 2.,
                ..options
            },
            SuppliedDegree11MaterialOptions {
                maximum_volume_discrepancy_mm3: 1e-12,
                ..options
            },
            SuppliedDegree11MaterialOptions {
                maximum_distance_mm: 1e-300,
                ..options
            },
        ] {
            assert!(matches!(
                volume.match_supplied_author(&author, frame, o, None),
                Err(PhillipsReadError::Unsupported)
            ));
        }
        let (different, _) = actual_m3_parts(false, false);
        assert!(matches!(
            volume.match_supplied_author(&different, frame, options, None),
            Err(PhillipsReadError::Unsupported)
        ));
        let cylindrical = crate::head_profile::Head::Cylindrical {
            radius_mm: 2.8,
            height_mm: 2.4,
        };
        let profile = crate::head_profile::Profile::new(
            cylindrical,
            author.profile().end()[0],
            author.spec().shank_length_mm,
        )
        .unwrap();
        let mut source_spec = author.spec();
        source_spec.head = Some(cylindrical);
        let different_head = PhillipsScrew::new(
            *author.thread(),
            profile,
            source_spec,
            author.cavity().clone(),
            100_000,
            &|| false,
        )
        .unwrap();
        // The actual SAME thread has an earned source/file relation, yet that
        // relation must never grant a different physical head a material match.
        assert!(
            head.supplied_thread_correspondence(&different_head, frame, Default::default(), None)
                .is_ok()
        );
        assert!(matches!(
            volume.match_supplied_author(&different_head, frame, options, None),
            Err(PhillipsReadError::Unsupported)
        ));
        let mut wrong = frame;
        wrong[0][0] += 1.;
        assert!(matches!(
            volume.match_supplied_author(&author, wrong, options, None),
            Err(PhillipsReadError::Unsupported)
        ));
        wrong[1] = [0.; 3];
        assert!(matches!(
            volume.match_supplied_author(&author, wrong, options, None),
            Err(PhillipsReadError::Invalid)
        ));
        assert!(matches!(
            volume.match_supplied_author(&author, frame, options, Some(&|| true)),
            Err(PhillipsReadError::Cancelled)
        ));
        let calls = Cell::new(0usize);
        volume
            .match_supplied_author(
                &author,
                frame,
                options,
                Some(&|| {
                    calls.set(calls.get() + 1);
                    false
                }),
            )
            .unwrap();
        let end = calls.get();
        calls.set(0);
        assert!(matches!(
            volume.match_supplied_author(
                &author,
                frame,
                options,
                Some(&|| {
                    calls.set(calls.get() + 1);
                    calls.get() == end
                })
            ),
            Err(PhillipsReadError::Cancelled)
        ));
        assert!(
            volume
                .match_supplied_author(&author, frame, options, None)
                .is_ok()
        );
        assert_eq!(material.body().literal_input(), raw);
    });
}

#[test]
fn actual_affine_gram_inverse_and_symbolic_origin_have_independent_oracles() {
    let frame = [[2., 3., 5.], [1., 0., 0.], [0., 1., 0.]];
    let target = CommonMouthChart {
        columns: [[2., 3., 7.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
        outer_top_z: 2.,
    };
    let exact = chart::compare(frame, target).unwrap();
    assert!(exact.translation.hi < 1e-12 && exact.linear.hi < 1e-12);
    assert!(exact.author_norm.hi >= 1. && exact.author_norm.hi < 1.00001);
    let identity = CommonMouthChart {
        columns: [[0.; 3], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
        outer_top_z: 0.,
    };
    let shear = chart::compare([[0.; 3], [1., 0., 0.], [1., 1., 0.]], identity).unwrap();
    let golden = (1. + 5f64.sqrt()) / 2.;
    assert!(shear.author_norm.hi >= golden && shear.author_norm.hi < 2.);
    assert!(shear.inverse_author_norm.hi >= golden && shear.inverse_author_norm.hi < 2.);
    let scale = chart::compare([[0.; 3], [2., 0., 0.], [0., 3., 0.]], identity).unwrap();
    assert!(scale.author_norm.hi >= 6. && scale.author_norm.hi < 6.00001);
    assert!(scale.inverse_author_norm.hi >= 0.5 && scale.inverse_author_norm.hi < 0.50001);
    assert!(scale.author_determinant.lo <= 36. && scale.author_determinant.hi >= 36.);
    assert!(chart::compare([[0.; 3], [0.; 3], [0., 1., 0.]], identity).is_err());
    let mut invalid = identity;
    invalid.columns[1][1] = f64::NAN;
    assert!(chart::compare(frame, invalid).is_err());
    assert!(!maximum([I::point(1.), I::point(f64::NAN)]).valid());
}
#[test]
fn normalized_ray_derivatives_cover_independent_vectors_and_opposite_points() {
    let bounds = chart::ray(
        I::point(2f64.sqrt()),
        I::point(18f64.sqrt()),
        I::point(5f64.sqrt()),
        I::point(5.),
    )
    .unwrap();
    // a=(3,4),b=(4,3),a'=(2,-1),b'=(-1,2). Their exact normalized
    // derivative difference is(.616,-.616), independently projected here.
    assert!(bounds[0].hi >= 0.08f64.sqrt());
    assert!(bounds[1].hi >= 2f64.sqrt() * 0.616);
    let opposite = chart::ray(I::point(2.), I::point(2.), I::point(1.), I::point(1.)).unwrap();
    assert!(opposite[0].hi >= 2. && opposite[1].hi >= 2.);
    assert!(chart::ray(I::point(0.), I::point(0.), I::point(1.), I::point(0.)).is_err());
    assert!(chart::ray(I::point(f64::NAN), I::point(1.), I::point(1.), I::point(1.)).is_err());
    assert!(chart::ray(I::point(1.), I::point(1.), I::point(1.), I::point(1e-300)).is_err());
}

#[test]
fn actual_h0_rounded_cavity_cannot_bypass_its_requested_comparison_cell_quota() {
    // Independent real H0 source: rounded DRIVE cavity, not a relabeled fit or
    // a change to the screw's canonical flat-root/shank parameters.
    let author = crate::step_read::phillips_cutter::physical::tests::fixture(0, false, false, 1);
    assert!(matches!(
        author.cavity().source().spec().inner,
        crate::phillips_drive::Inner::Rounded { .. }
    ));
    let frame = [[0.; 3], [1., 0., 0.], [0., 1., 0.]];
    let raw =
        crate::step::write_phillips_screw(&author, frame[0], frame[1], frame[2], 1e-7).unwrap();
    with_head(&raw, &mut |head| {
        let occupied = head
            .degree11_thread_carrier(Default::default(), None)
            .unwrap()
            .certify_occupied(Default::default(), None)
            .unwrap();
        let material = head
            .closed_degree11_material(&occupied, Default::default(), None)
            .unwrap();
        let volume = material.volume_bounds(Default::default(), None).unwrap();
        assert!(matches!(
            volume.match_supplied_author(
                &author,
                frame,
                SuppliedDegree11MaterialOptions {
                    maximum_cavity_cells: 1,
                    ..Default::default()
                },
                None
            ),
            Err(PhillipsReadError::Limit)
        ));
        assert_eq!(material.body().literal_input(), raw);
    });
}

#[test]
fn allocated_pair_accuracy_satisfies_independent_complete_annulus_inequalities() {
    let identity = CommonMouthChart {
        columns: [[0.; 3], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
        outer_top_z: 0.,
    };
    for scale in [1., 2.] {
        let frame = [[0.; 3], [scale, 0., 0.], [0., scale, 0.]];
        let target = CommonMouthChart {
            columns: [[0.; 3], frame[1], frame[2], [0., 0., scale * scale]],
            outer_top_z: 0.,
        };
        let world = chart::compare(frame, target).unwrap();
        let budget = CavityBudget {
            world: &world,
            extent: I::point(10.),
            other_target_speed: I::point(5.),
            options: Default::default(),
        };
        let [ep, ed] = cavity_allocation(
            &budget,
            I::point(0.5),
            I::point(2.),
            I::point(1e-12),
            I::point(1e-12),
            [I::point(3.), I::point(4.)],
        )
        .unwrap();
        // Direct independent real scalar ray/annulus inequalities, multiplied
        // by the actual affine operator upper bound, not a copied allocation.
        let position = 2e-12 + 8. * ep;
        let du = 6e-12 + (2. + 1e-12) * ed / 0.5 + 5. * 2. * 3. * ep / (0.5 * 0.5);
        let dv = 1e-12 + 9. * ep;
        assert!(world.author_norm.hi * position < 1e-7);
        assert!(world.author_norm.hi * du.max(dv).max(ed) < 1e-7);
        assert!(ep > 0. && ed > 0. && ep <= 1e-7 && ed <= 1e-7);
    }
    let mut world = chart::compare([[0.; 3], [1., 0., 0.], [0., 1., 0.]], identity).unwrap();
    world.translation = I::point(2e-7);
    let budget = CavityBudget {
        world: &world,
        extent: I::point(10.),
        other_target_speed: I::point(5.),
        options: Default::default(),
    };
    assert!(
        cavity_allocation(
            &budget,
            I::point(0.5),
            I::point(2.),
            I::point(0.),
            I::point(0.),
            [I::point(3.), I::point(4.)]
        )
        .is_err()
    );
    assert!(
        cavity_allocation(
            &budget,
            I::point(0.),
            I::point(2.),
            I::point(0.),
            I::point(0.),
            [I::point(3.), I::point(4.)]
        )
        .is_err()
    );
}
