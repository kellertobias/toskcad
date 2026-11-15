use super::*;
fn identity() -> Frame {
    Frame {
        origin_mm: [0.; 3],
        axes: [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
    }
}
fn options(engagement: f64) -> ClockedMatingOptions {
    ClockedMatingOptions {
        requirements: MatingRequirements {
            minimum_engagement_mm: engagement,
            minimum_clearance_mm: 4e-7,
        },
        maximum_cells: 4096,
        linear_tolerance_mm: 1e-7,
    }
}
#[test]
fn shared_clocked_receiver_matches_established_insert_fixture_before_migration() {
    for left in [false, true] {
        let (insert, screw, frame) = super::super::tests::fixture(left);
        let old =
            established_insert_reference(&insert, identity(), &screw, frame, options(4.), &|| {
                false
            })
            .unwrap();
        let shared = certify_clocked_receiver(
            &Receiver::insert(&insert, identity()).unwrap(),
            MatingExternal {
                source: MatingScrew::Regular(&screw),
                frame,
            },
            options(4.),
            &|| false,
        )
        .unwrap();
        let a = &old.mating;
        let b = &shared.bounds;
        assert_eq!(a.engagement_mm, b.engagement_mm);
        assert_eq!(
            a.minimum_radial_clearance_lower_mm,
            b.minimum_radial_clearance_lower_mm
        );
        assert_eq!(a.source_phase_offset_mm, b.source_phase_offset_mm);
        assert_eq!(a.prefix_entry_clearance_mm, b.prefix_entry_clearance_mm);
        assert_eq!(a.thread_exit_clearance_mm, b.thread_exit_clearance_mm);
        assert_eq!(
            a.entire_thread_within_insert,
            b.entire_thread_within_receiver
        );
        assert_eq!(
            a.external_indexed_phase_error_mm,
            b.external_indexed_phase_error_mm
        );
        assert_eq!(
            a.external_evaluator_phase_error_mm,
            b.external_evaluator_phase_error_mm
        );
        assert_eq!(
            a.inner_whole_period_envelope_used,
            b.inner_whole_period_envelope_used
        );
        assert_eq!(
            a.outer_whole_source_envelope_used,
            b.outer_whole_source_envelope_used
        );
        assert_eq!(old.radial_map_error_mm, shared.radial_map_error_mm);
        assert_eq!(old.axial_map_error_mm, shared.axial_map_error_mm);
        assert_eq!(old.angle_error_rad, shared.angle_error_rad);
        assert_eq!(old.cells_used, shared.cells_used);
    }
}
fn builder_fixture() -> (
    crate::line8_nut::Nut,
    crate::profiled_screw::Screw,
    Frame,
    Frame,
) {
    use crate::{
        authored_flat_thread::{AuthoredThread, SourceOptions},
        head_profile::Head,
        helical_thread::{Spec, Thread},
        line8_nut::Relief,
        linear_sweep::Sweep,
        profiled_screw::Spec as ScrewSpec,
    };
    // Exact original Builder::screw/Builder::nut construction inputs. Derive
    // stations using their original expressions, rather than rounded .65/3.85.
    let pitch = 1.25_f64;
    let depth = 5. * 3_f64.sqrt() / 16. * pitch;
    let height = 4.4;
    let underhead = 14. - height;
    let runout = 1.25;
    let neck = -5. + height;
    let external_start = neck + runout;
    let seat = 4.5;
    let station: f64 = seat - external_start;
    let angle = std::f64::consts::TAU * station.rem_euclid(pitch) / pitch;
    let (s, c) = angle.sin_cos();
    let raw = s.atan2(c);
    let roll = (raw + std::f64::consts::TAU / 2.).rem_euclid(std::f64::consts::TAU)
        - std::f64::consts::TAU / 2.;
    let (s, c) = roll.sin_cos();
    let external = Thread::new(Spec {
        major_radius_mm: 4. - 0.08,
        pitch_mm: pitch,
        length_mm: underhead - runout,
        depth_mm: depth,
        crest_width_mm: pitch / 8. - 0.04,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts: 1,
        left_hand: false,
        root_radius_mm: None,
    })
    .unwrap();
    let screw = crate::profiled_screw::Screw::new_with_shank_radius(
        external,
        4.,
        ScrewSpec {
            head: Head::Countersunk {
                radius_mm: 8.,
                height_mm: height,
                included_angle_deg: 90.,
                edge_chamfer_mm: 0.05,
                neck_radius_mm: 0.4,
            },
            shank_length_mm: underhead - external.spec().length_mm,
            socket_across_flats_mm: 5.08,
            socket_depth_mm: 3.2,
            minimum_wall_mm: 1e-5,
        },
    )
    .unwrap();
    let internal = AuthoredThread::new(
        Spec {
            major_radius_mm: 4. + 0.08,
            length_mm: 4.1,
            crest_width_mm: pitch / 8. + 0.04,
            ..external.spec()
        },
        SourceOptions {
            position_limit_mm: 1e-8,
            normal_angle_limit_rad: 1e-8,
        },
    )
    .unwrap();
    let section = vec![
        [-6.9, 0.],
        [6.9, 0.],
        [6.5, 3.],
        [5.4, 5.],
        [3.3, 7.3],
        [-3.3, 7.3],
        [-5.4, 5.],
        [-6.5, 3.],
    ];
    let stock = Sweep::new(vec![section.clone()], vec![section], 23.).unwrap();
    let nut = crate::line8_nut::Nut::new(
        stock,
        7.5,
        internal,
        [[c, -s], [s, c]],
        Relief {
            start_mm: 4.1 - 0.001,
            radius_mm: 4.10,
            tool_length_mm: 7.3 - 4.1 + 0.02,
        },
        None,
    )
    .unwrap();
    let stock_frame = Frame {
        origin_mm: [-7.5, 0., seat],
        ..identity()
    };
    let screw_frame = Frame {
        origin_mm: [0., 0., external_start],
        ..identity()
    };
    (nut, screw, stock_frame, screw_frame)
}
#[test]
fn original_builder_m8_nut_stock_relief_frames_and_screw_qualify_independently() {
    let (nut, screw, stock_frame, screw_frame) = builder_fixture();
    let placement = ClockedLine8Placement {
        source: &nut,
        stock_frame,
        outer_frame: identity(),
        minimum_wall_mm: 0.5,
    };
    let proof = certify_clocked_line8_mating(
        placement,
        MatingExternal {
            source: MatingScrew::Profiled(&screw),
            frame: screw_frame,
        },
        options(2. * 1.25),
        &|| false,
    )
    .unwrap();
    assert!(std::ptr::eq(proof.nut().source, &nut));
    assert!(matches!(proof.external().source,MatingScrew::Profiled(s)if std::ptr::eq(s,&screw)));
    assert_eq!(nut.thread().spec().length_mm, 4.1);
    assert_eq!(nut.relief().start_mm, 4.1 - 0.001);
    assert!(proof.complete_turns_lower >= 2.);
    assert!(proof.relief_clearance_lower_mm > 0.1);
    assert!(proof.stock_wall_lower_mm > 0.5);
}

#[test]
fn original_line8_proof_preserves_symbolic_outer_placement_and_rejects_real_interference() {
    let (nut, screw, stock_frame, screw_frame) = builder_fixture();
    let mut outer = identity();
    outer.origin_mm = [37., -11., 103.];
    outer.axes = [[0., 1., 0.], [-1., 0., 0.], [0., 0., 1.]];
    let mut placed_screw = screw_frame;
    placed_screw.origin_mm = [37., -11., 103. + screw_frame.origin_mm[2]];
    placed_screw.axes = outer.axes;
    let placement = ClockedLine8Placement {
        source: &nut,
        stock_frame,
        outer_frame: outer,
        minimum_wall_mm: 0.5,
    };
    let external = MatingExternal {
        source: MatingScrew::Profiled(&screw),
        frame: placed_screw,
    };
    let proof = certify_clocked_line8_mating(placement, external, options(2.5), &|| false).unwrap();
    assert!(proof.complete_turns_lower > 3.);
    assert!(proof.mating.minimum_radial_clearance_lower_mm > 4e-7);
    let mut interference = placed_screw;
    interference.origin_mm[0] += 0.1;
    assert!(matches!(
        certify_clocked_line8_mating(
            placement,
            MatingExternal {
                frame: interference,
                ..external
            },
            options(2.5),
            &|| false
        ),
        Err(Error::Unsupported)
    ));
    assert!(matches!(
        certify_clocked_line8_mating(
            ClockedLine8Placement {
                minimum_wall_mm: 2.,
                ..placement
            },
            external,
            options(2.5),
            &|| false
        ),
        Err(Error::Unsupported)
    ));
    assert!(matches!(
        certify_clocked_line8_mating(placement, external, options(2.5), &|| true),
        Err(Error::Cancelled)
    ));
    assert!(matches!(
        certify_clocked_line8_mating(
            placement,
            external,
            ClockedMatingOptions {
                maximum_cells: 1,
                ..options(2.5)
            },
            &|| false
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        certify_clocked_line8_mating(
            placement,
            external,
            ClockedMatingOptions {
                maximum_cells: 0,
                ..options(2.5)
            },
            &|| false
        ),
        Err(Error::Invalid)
    ));
    let calls = std::cell::Cell::new(0usize);
    assert!(matches!(
        certify_clocked_line8_mating(placement, external, options(2.5), &|| {
            calls.set(calls.get() + 1);
            calls.get() > 10
        }),
        Err(Error::Cancelled)
    ));
}

#[test]
fn shared_clocked_insert_preserves_established_public_failure_categories() {
    for left in [false, true] {
        let (insert, screw, frame) = super::super::tests::fixture(left);
        for n in 0..8 {
            let mut request = options(4.);
            let mut placed = frame;
            match n {
                0 => request.maximum_cells = 0,
                1 => request.maximum_cells = 1,
                2 => request.requirements.minimum_engagement_mm = 4_f64.next_up(),
                3 => placed.origin_mm[2] = 0.01,
                4 => placed.origin_mm[0] = 0.01,
                5 => placed.axes[2] = [0., 0., -1.],
                6 => request.linear_tolerance_mm = f64::NAN,
                _ => request.requirements.minimum_clearance_mm = 1.,
            }
            let previous =
                established_insert_reference(&insert, identity(), &screw, placed, request, &|| {
                    false
                })
                .err()
                .unwrap();
            let shared = certify_clocked_insert_mating(
                &insert,
                identity(),
                &screw,
                placed,
                request,
                &|| false,
            )
            .err()
            .unwrap();
            assert_eq!(previous, shared, "hand {left}, case {n}");
        }
        let wrong = crate::helical_thread::Thread::new(crate::helical_thread::Spec {
            pitch_mm: 0.6,
            root_radius_mm: None,
            ..screw.thread().spec()
        })
        .unwrap();
        let wrong = crate::headed_screw::Screw::new(wrong, screw.spec()).unwrap();
        assert_eq!(
            established_insert_reference(&insert, identity(), &wrong, frame, options(4.), &|| {
                false
            })
            .err(),
            certify_clocked_insert_mating(&insert, identity(), &wrong, frame, options(4.), &|| {
                false
            })
            .err()
        );
        let short = crate::helical_thread::Thread::new(crate::helical_thread::Spec {
            length_mm: 4.,
            ..screw.thread().spec()
        })
        .unwrap();
        let short = crate::headed_screw::Screw::new(short, screw.spec()).unwrap();
        assert_eq!(
            established_insert_reference(&insert, identity(), &short, frame, options(4.), &|| {
                false
            })
            .err(),
            certify_clocked_insert_mating(&insert, identity(), &short, frame, options(4.), &|| {
                false
            })
            .err()
        );
    }
}

#[test]
fn shared_clocked_source_boxes_cover_actual_phillips_head_and_keep_line8_admission_typed() {
    use crate::{
        head_profile::{Head, Profile},
        phillips_drive::{self, Cavity, Cutter, PhillipsScrew},
        profiled_screw::SourceSpec,
    };
    let (insert, regular, placed) = super::super::tests::fixture(false);
    let head = Head::Countersunk {
        radius_mm: 2.75,
        height_mm: 1.65,
        included_angle_deg: 90.,
        edge_chamfer_mm: 0.05,
        neck_radius_mm: 0.1,
    };
    let (shape, datum) = phillips_drive::test_fixture(2);
    let cavity = Cavity::new(Cutter::new(shape, datum).unwrap(), 1e-5, 100_000).unwrap();
    let screw = PhillipsScrew::new(
        *regular.thread(),
        Profile::new(head, 1.5, 1.).unwrap(),
        SourceSpec {
            head: Some(head),
            shank_length_mm: 1.,
            minimum_wall_mm: 1e-5,
        },
        cavity,
        100_000,
        &|| false,
    )
    .unwrap();
    let receiver = Receiver::insert(&insert, identity()).unwrap();
    let proof = certify_clocked_receiver(
        &receiver,
        MatingExternal {
            source: MatingScrew::Phillips(&screw),
            frame: placed,
        },
        options(4.),
        &|| false,
    )
    .unwrap();
    assert!(proof.bounds.minimum_radial_clearance_lower_mm > 0.007);
    let (nut, _, stock_frame, _) = builder_fixture();
    let placement = ClockedLine8Placement {
        source: &nut,
        stock_frame,
        outer_frame: identity(),
        minimum_wall_mm: 0.5,
    };
    assert!(matches!(
        certify_clocked_line8_mating(
            placement,
            MatingExternal {
                source: MatingScrew::Regular(&regular),
                frame: placed
            },
            options(2.5),
            &|| false
        ),
        Err(Error::Unsupported)
    ));
    assert!(matches!(
        certify_clocked_line8_mating(
            placement,
            MatingExternal {
                source: MatingScrew::Phillips(&screw),
                frame: placed
            },
            options(2.5),
            &|| false
        ),
        Err(Error::Unsupported)
    ));
}

#[test]
fn finite_original_prefix_must_earn_two_actual_turns_and_physical_engagement() {
    use crate::authored_flat_thread::{AuthoredThread, SourceOptions};
    let (original, screw, stock_frame, screw_frame) = builder_fixture();
    let source_options = SourceOptions {
        position_limit_mm: 1e-8,
        normal_angle_limit_rad: 1e-8,
    };
    let pitch = 2.;
    let male = crate::helical_thread::Thread::new(crate::helical_thread::Spec {
        pitch_mm: pitch,
        ..screw.thread().spec()
    })
    .unwrap();
    let external = crate::profiled_screw::Screw::new_with_recess(
        male,
        4.,
        screw.spec().head.unwrap(),
        screw.spec().shank_length_mm,
        screw.spec().minimum_wall_mm,
        screw.recess().clone(),
    )
    .unwrap();
    let internal = AuthoredThread::new(
        crate::helical_thread::Spec {
            pitch_mm: pitch,
            ..original.thread().spec()
        },
        source_options,
    )
    .unwrap();
    let station = stock_frame.origin_mm[2] - screw_frame.origin_mm[2];
    let angle = std::f64::consts::TAU * station.rem_euclid(pitch) / pitch;
    let (s, c) = angle.sin_cos();
    let raw = s.atan2(c);
    let roll = (raw + std::f64::consts::TAU / 2.).rem_euclid(std::f64::consts::TAU)
        - std::f64::consts::TAU / 2.;
    let (s, c) = roll.sin_cos();
    let nut = crate::line8_nut::Nut::new(
        original.stock().clone(),
        original.bore_from_start_mm(),
        internal,
        [[c, -s], [s, c]],
        crate::line8_nut::Relief {
            start_mm: 3.5,
            tool_length_mm: 7.3 - 3.5 + 0.02,
            ..original.relief()
        },
        None,
    )
    .unwrap();
    let placement = ClockedLine8Placement {
        source: &nut,
        stock_frame,
        outer_frame: identity(),
        minimum_wall_mm: 0.5,
    };
    let receiver = Receiver {
        source: nut.thread(),
        length: nut.relief().start_mm,
        frame: SourceFrame::literal(stock_frame)
            .unwrap()
            .translated(0, nut.bore_from_start_mm())
            .unwrap()
            .phase(nut.phase_rotation())
            .unwrap(),
    };
    let request = options(1.);
    let actual = MatingExternal {
        source: MatingScrew::Profiled(&external),
        frame: screw_frame,
    };
    let core = certify_clocked_receiver(&receiver, actual, request, &|| false).unwrap();
    assert!(core.external_span.hi / pitch < 2.);
    assert!(matches!(
        certify_clocked_line8_mating(placement, actual, request, &|| false),
        Err(Error::Unsupported)
    ));
    let original_placement = ClockedLine8Placement {
        source: &original,
        ..placement
    };
    assert!(matches!(
        certify_clocked_line8_mating(
            original_placement,
            MatingExternal {
                source: MatingScrew::Profiled(&screw),
                frame: screw_frame
            },
            options(original.relief().start_mm),
            &|| false
        ),
        Err(Error::Unsupported)
    ));
}
