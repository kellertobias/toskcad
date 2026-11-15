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
