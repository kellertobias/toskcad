use super::super::common_mouth::tests::derived_thread_probe::with_head;
use super::*;
fn text(scale: f64, posed: bool) -> String {
    super::super::super::branches::representative::scaled_text(scale, posed)
}
#[test]
fn actual_same_body_sources_form_a_closed_repaired_material_with_real_rim_and_faces() {
    for (scale, posed) in [(1., false), (1.137, true)] {
        let input = text(scale, posed);
        with_head(&input, &mut |head| {
            let carrier = head
                .derived_thread_carrier(Default::default(), None)
                .unwrap();
            let work = carrier.work_used();
            let radius = carrier.entry_radius_upper_mm();
            let thread = carrier.certify_occupied(1e-7, 1e-5, None).unwrap();
            assert_eq!(thread.construction_work_used(), work);
            assert_eq!(thread.entry_radius_upper_mm(), radius);
            assert!(thread.retained_budget_bytes() > head.body().retained_budget_bytes());
            let material = head
                .closed_repaired_material(&thread, Default::default(), None)
                .unwrap();
            assert!(std::ptr::eq(material.body(), head.body()));
            assert!(std::ptr::eq(material.head(), head.source()));
            assert!(std::ptr::eq(material.thread(), thread.boundary()));
            assert!(std::ptr::eq(
                material.cavity(),
                head.common_mouth()
                    .closed_cavity()
                    .local_source_certificate()
                    .source()
            ));
            assert_eq!(material.chart(), head.chart());
            assert!(material.determinant_bounds()[0] > 0.);
            assert!(material.join_radial_margin_mm() > 1e-7);
            assert!(material.containment().axial_floor_clearance_mm() > 1e-7);
            assert!(material.work_used() < WORK);
            assert!(material.retained_budget_bytes() > thread.retained_budget_bytes());
            assert_eq!(
                material.join_inner_coedges(),
                &thread.boundary().input().faces[0].outer
            );
            let faces = (0..material.face_count())
                .map(|i| material.face_source(i).unwrap())
                .collect::<Vec<_>>();
            assert!(
                faces
                    .iter()
                    .any(|f| matches!(f, MaterialFaceSource::HeadMeridian { .. }))
            );
            assert!(
                faces
                    .iter()
                    .any(|f| matches!(f, MaterialFaceSource::ThreadSide { .. }))
            );
            assert!(faces.contains(&MaterialFaceSource::ThreadTip));
            assert!(faces.contains(&MaterialFaceSource::TopAnnulus));
            assert!(faces.contains(&MaterialFaceSource::JoinAnnulus));
            assert_eq!(
                faces
                    .iter()
                    .filter(|f| matches!(f, MaterialFaceSource::ReversedCavity { .. }))
                    .count(),
                head.common_mouth()
                    .closed_cavity()
                    .local_source_certificate()
                    .topology()
                    .faces
                    - 1
            );
            assert!(material.face_source(material.face_count()).is_none());
            assert!(material.face_source(usize::MAX).is_none());
            assert_eq!(head.body().literal_input(), input);
        });
    }
}
#[test]
fn actual_material_queries_enforce_aggregate_budgets_and_atomic_cancellation() {
    let input = text(1., false);
    with_head(&input, &mut |head| {
        let thread = head
            .derived_thread_carrier(Default::default(), None)
            .unwrap()
            .certify_occupied(1e-7, 1e-5, None)
            .unwrap();
        let proof = head
            .closed_repaired_material(&thread, Default::default(), None)
            .unwrap();
        let peak =
            proof.retained_budget_bytes() + RepairedHeadContainment::scratch_budget_bytes(CELLS);
        let options = ClosedMaterialOptions {
            maximum_work: proof.work_used(),
            maximum_bytes: peak,
            ..Default::default()
        };
        assert!(
            head.closed_repaired_material(&thread, options, None)
                .is_ok()
        );
        for options in [
            ClosedMaterialOptions {
                maximum_work: 0,
                ..options
            },
            ClosedMaterialOptions {
                maximum_work: 1,
                ..options
            },
            ClosedMaterialOptions {
                maximum_work: 1
                    + head.work_used()
                    + thread.construction_work_used()
                    + METADATA_WORK
                    - 1,
                ..options
            },
            ClosedMaterialOptions {
                maximum_work: options.maximum_work - 1,
                ..options
            },
            ClosedMaterialOptions {
                maximum_bytes: 0,
                ..options
            },
            ClosedMaterialOptions {
                maximum_bytes: peak - 1,
                ..options
            },
            ClosedMaterialOptions {
                maximum_cells: 1,
                ..options
            },
            ClosedMaterialOptions {
                maximum_work: 1
                    + head.work_used()
                    + thread.construction_work_used()
                    + METADATA_WORK
                    + head
                        .common_mouth()
                        .closed_cavity()
                        .local_source_certificate()
                        .work_used()
                    + head.source().work_used()
                    - 1,
                ..options
            },
        ] {
            assert!(matches!(
                head.closed_repaired_material(&thread, options, None),
                Err(PhillipsReadError::Limit)
            ));
        }
        for options in [
            ClosedMaterialOptions {
                minimum_clearance_mm: f64::NAN,
                ..options
            },
            ClosedMaterialOptions {
                minimum_clearance_mm: 0.,
                ..options
            },
            ClosedMaterialOptions {
                maximum_cells: 0,
                ..options
            },
            ClosedMaterialOptions {
                maximum_cells: CELLS + 1,
                ..options
            },
            ClosedMaterialOptions {
                maximum_work: WORK + 1,
                ..options
            },
            ClosedMaterialOptions {
                maximum_bytes: BYTES + 1,
                ..options
            },
        ] {
            assert!(matches!(
                head.closed_repaired_material(&thread, options, None),
                Err(PhillipsReadError::Invalid)
            ));
        }
        assert!(matches!(
            head.closed_repaired_material(
                &thread,
                ClosedMaterialOptions {
                    minimum_clearance_mm: 1.,
                    ..Default::default()
                },
                None
            ),
            Err(PhillipsReadError::Unsupported)
        ));
        let n = std::cell::Cell::new(0);
        head.closed_repaired_material(
            &thread,
            options,
            Some(&|| {
                n.set(n.get() + 1);
                false
            }),
        )
        .unwrap();
        for at in 0..n.get() {
            let calls = std::cell::Cell::new(0);
            assert!(matches!(
                head.closed_repaired_material(
                    &thread,
                    options,
                    Some(&|| {
                        let old = calls.get();
                        calls.set(old + 1);
                        old == at
                    })
                ),
                Err(PhillipsReadError::Cancelled)
            ));
        }
        assert!(
            head.closed_repaired_material(&thread, options, None)
                .is_ok()
        );
    });
}
#[test]
fn distinct_real_parsed_arenas_cannot_be_mixed_even_for_equal_input_geometry() {
    let input = text(1., false);
    with_head(&input, &mut |head| {
        let thread = head
            .derived_thread_carrier(Default::default(), None)
            .unwrap()
            .certify_occupied(1e-7, 1e-5, None)
            .unwrap();
        with_head(&input, &mut |other| {
            assert!(!std::ptr::eq(head.body(), other.body()));
            assert!(matches!(
                other.closed_repaired_material(&thread, Default::default(), None),
                Err(PhillipsReadError::Unsupported)
            ));
        });
        assert!(
            head.closed_repaired_material(&thread, Default::default(), None)
                .is_ok()
        );
    });
}
fn shallow_text() -> String {
    use crate::{
        head_profile::{Head, Profile},
        helical_thread::{Spec, Thread},
        phillips_drive::{Cavity, PhillipsScrew},
        profiled_screw::SourceSpec,
    };
    let (s, d) = crate::phillips_drive::test_fixture(0);
    let cutter = Cutter::new(s, d).unwrap();
    let height = cutter.envelope().tip_depth_mm_bounds()[1] + 0.001;
    let cavity = Cavity::new(cutter, 1e-5, 100_000).unwrap();
    let head = Head::Cylindrical {
        radius_mm: 2.,
        height_mm: height,
    };
    let thread = Thread::new(Spec {
        major_radius_mm: 0.97,
        pitch_mm: 0.4,
        length_mm: 2.,
        depth_mm: 0.24,
        crest_width_mm: 0.05,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts: 1,
        left_hand: false,
        root_radius_mm: None,
    })
    .unwrap();
    let screw = PhillipsScrew::new(
        thread,
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
    crate::step::write_phillips_screw(&screw, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7).unwrap()
}
#[test]
fn real_arbitrary_shallow_head_requires_actual_cavity_floor_clearance() {
    let input = shallow_text();
    with_head(&input, &mut |head| {
        let thread = head
            .derived_thread_carrier(Default::default(), None)
            .unwrap()
            .certify_occupied(1e-7, 1e-5, None)
            .unwrap();
        let accepted = head
            .closed_repaired_material(&thread, Default::default(), None)
            .unwrap();
        assert!(accepted.containment().axial_floor_clearance_mm() > 0.0009);
        assert!(accepted.containment().axial_floor_clearance_mm() < 0.0011);
        assert!(accepted.join_radial_margin_mm() > 0.01);
        assert!(matches!(
            head.closed_repaired_material(
                &thread,
                ClosedMaterialOptions {
                    minimum_clearance_mm: 0.005,
                    ..Default::default()
                },
                None
            ),
            Err(PhillipsReadError::Unsupported)
        ));
        assert!(
            head.closed_repaired_material(&thread, Default::default(), None)
                .is_ok()
        );
    });
}
