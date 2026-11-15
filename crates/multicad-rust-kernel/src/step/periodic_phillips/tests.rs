use super::*;
fn source(left: bool, extra: f64) -> Material {
    use crate::{
        head_profile::{Head, Profile},
        helical_thread::{Spec as ThreadSpec, Thread},
        phillips_drive::{Cutter, Domain, Inner, Spec},
        profiled_screw::SourceSpec,
    };
    // The actual ISO7045 M2 normal factory packet, independently qualified by
    // the public authored source constructor. No fixture values enter decoding.
    let head = Head::Crowned {
        radius_mm: 2.,
        height_mm: 1.6,
        crown_radius_mm: 3.2,
        top_land_radius_mm: 1.235645633015642,
        neck_radius_mm: 0.1,
    };
    let profile = Profile::new_for_thread_union(head, 1., 0.5).unwrap();
    let thread = Thread::new(ThreadSpec {
        major_radius_mm: 1.,
        pitch_mm: 0.4,
        length_mm: 11.5 + extra,
        depth_mm: 0.21650635094610965,
        crest_width_mm: 0.05,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts: 1,
        left_hand: left,
        root_radius_mm: None,
    })
    .unwrap();
    let cutter = Cutter::new(
        Spec {
            mouth_diameter_mm: 1.9,
            mouth_radius_mm: 0.3,
            b_width_mm: 0.61,
            g_diameter_mm: 0.81,
            tilt_angle_deg: 7.,
            outer_flank_angle_deg: 26.5,
            floor_angle_deg: 28.,
            outer_corner_angle_deg: 92.,
            inner: Inner::Rounded { radius_mm: 0.25 },
        },
        Domain {
            corner_start_depth_mm: -3.8,
            corner_end_depth_mm: 3.208443210824065,
            corner_reach_mm: 16.,
        },
    )
    .unwrap();
    Material::from_authored_parts(
        &thread,
        &profile,
        &cutter,
        SourceSpec {
            head: Some(head),
            shank_length_mm: 0.5,
            minimum_wall_mm: 1e-5,
        },
        Default::default(),
        None,
    )
    .unwrap()
}
fn occurrence(posed: bool) -> [[f64; 3]; 4] {
    if posed {
        [
            [105., -20., 60.],
            [0.8, 0.6, 0.],
            [0.36, -0.48, -0.8],
            [-0.48, 0.64, -0.6],
        ]
    } else {
        [[0.; 3], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]
    }
}
#[test]
fn actual_normal_metric_export_independently_decodes_matches_and_reexports() {
    for (left, posed) in [(false, false), (true, true)] {
        let source = source(left, 0.);
        let frame = occurrence(posed);
        let now = std::time::Instant::now();
        let file = export(&source, frame, Default::default(), None).unwrap();
        let verified = file.qualification().verification();
        assert!(file.as_str().starts_with("ISO-10303-21;"));
        assert!(file.proposal_work_used() > 0);
        assert!(verified.maximum_world_distance_mm <= 1e-7);
        assert!(verified.maximum_world_chart_derivative_mm <= 1e-7);
        assert!(verified.independent_volume_discrepancy_mm3 <= 1e-5);
        assert!(verified.reserved_reader_work > 16_000_000);
        assert!(verified.decoded_material_proof_work > 0);
        assert!(verified.correspondence_fixed_work > 0 && verified.correspondence_cavity_cells > 0);
        let decoded = file.qualification().material();
        assert!(matches!(
            decoded.provenance(),
            step_read::PhillipsPeriodicMaterialProvenance::Decoded { .. }
        ));
        assert!(decoded.original_input().unwrap() == file.as_str());
        assert!(
            decoded
                .original_curve_owners()
                .unwrap()
                .iter()
                .all(|c| c.edge_entity > 0)
        );
        eprintln!(
            "qualified actual metric STEP LH{left} pose{posed} {}bytes {:?} C0{} D{} mass{} fixed{} cells{}",
            file.as_str().len(),
            now.elapsed(),
            verified.maximum_world_distance_mm,
            verified.maximum_world_chart_derivative_mm,
            verified.independent_volume_discrepancy_mm3,
            verified.correspondence_fixed_work,
            verified.correspondence_cavity_cells
        );
        let retained = qualify_actual_file(&source, file.as_str(), frame, Default::default(), None)
            .unwrap()
            .into_material();
        // Reexport the independently decoded geometry with its actual stored
        // construction chart. The authored-only proposal getter stays guarded.
        assert!(
            retained
                .authored_mesh_approximation(16_000_000, None)
                .is_err()
        );
        export(&retained, occurrence(false), Default::default(), None).unwrap();
    }
}
#[test]
fn actual_stale_and_corrupted_fields_limits_and_cancellation_refuse() {
    let expected = source(false, 0.);
    let frame = occurrence(false);
    let file = export(&expected, frame, Default::default(), None).unwrap();
    let stale = export(&source(false, 0.4), frame, Default::default(), None).unwrap();
    assert!(matches!(
        qualify_actual_file(&expected, stale.as_str(), frame, Default::default(), None),
        Err(Error::Unsupported)
    ));
    let tampered = file.as_str().replacen("CARTESIAN_POINT(", "DIRECTION(", 1);
    assert_ne!(tampered, file.as_str());
    assert!(matches!(
        qualify_actual_file(&expected, &tampered, frame, Default::default(), None),
        Err(Error::Reader(_))
    ));
    for cap in [0, 64 * 1024 * 1024 + 1] {
        assert!(matches!(
            qualify_actual_file(
                &expected,
                file.as_str(),
                frame,
                Options {
                    maximum_output_bytes: cap,
                    ..Default::default()
                },
                None
            ),
            Err(Error::Invalid)
        ));
    }
    assert!(matches!(
        qualify_actual_file(
            &expected,
            file.as_str(),
            frame,
            Options {
                maximum_output_bytes: file.as_str().len() - 1,
                ..Default::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert!(
        qualify_actual_file(
            &expected,
            file.as_str(),
            frame,
            Options {
                maximum_output_bytes: file.as_str().len(),
                ..Default::default()
            },
            None
        )
        .is_ok()
    );
    assert!(matches!(
        export(
            &expected,
            frame,
            Options {
                maximum_proposal_work: 0,
                ..Default::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
    for boundary in [1, 10, 100, 1000] {
        let calls = std::cell::Cell::new(0);
        let stop = || {
            calls.set(calls.get() + 1);
            calls.get() >= boundary
        };
        assert!(matches!(
            export(&expected, frame, Default::default(), Some(&stop)),
            Err(Error::Cancelled)
                | Err(Error::Reader(step_read::PhillipsPeriodicReadError {
                    cause: PhillipsReadError::Cancelled,
                    ..
                }))
        ));
    }
}
