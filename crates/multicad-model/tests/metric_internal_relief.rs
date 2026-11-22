#[path = "metric_single_internal_support/mod.rs"]
mod support;
#[test]
fn normal_metric_tail_only_receiving_cut_builds_for_both_hands_and_real_canonical_lengths() {
    support::require_positive(false);
}

#[test]
fn actual_tail_only_source_manufactures_edits_clears_saves_and_reopens_without_step_output() {
    use multicad_core::{
        commands::{Actor, Channel, Session},
        constraint_sketch::Dimension,
        model::{self, ThreadEndTreatment, ThreadStartTreatment, length},
        parameters::{self, BoundExpression},
        project_file::{self, ProjectFile},
    };
    use multicad_kernel::{GeometryErrorKind, MeshOptions, RigidTransform};
    use multicad_model::Compiler;
    for left in [false, true] {
        for expression in [false, true] {
            let height = if expression {
                Dimension::Expression {
                    expression: BoundExpression {
                        text: "=1.1*3".into(),
                        unit: "mm".into(),
                        bindings: vec![],
                    },
                }
            } else {
                length(3.3000000000000003)
            };
            let (mut session, thread) = support::source_dimension(false, left, height);
            let compile = |s: &Session| {
                (Compiler {
                    omit_display: true,
                    ..Default::default()
                })
                .compile(s.document(), s.revision())
                .unwrap()
            };
            let built = compile(&session);
            let body = built.solid(thread.feature).unwrap();
            let data = body.rust_internal_thread_treatment_data().unwrap();
            let end = parameters::exact_length_mm(&thread.parameters.length, session.parameters())
                .unwrap();
            assert_eq!(data.material.stock().height(), end);
            assert_eq!(data.material.cavity().exact_receiving_extent(), Some(end));
            assert_eq!(data.material.cavity().entry_radius_mm(), None);
            assert_eq!(
                data.material.cavity().source_thread().spec().left_hand,
                left
            );
            let before_volume = data.material.volume_bounds(Default::default()).unwrap();
            let scalar = body.volume_mm3().unwrap();
            assert!(before_volume.lower <= scalar && scalar <= before_volume.upper);
            let source_before = data.material;
            let moved = body
                .transformed(RigidTransform {
                    rotation: [[0.6, 0., 0.8], [0., 1., 0.], [-0.8, 0., 0.6]],
                    translation_mm: [11., -7., 5.],
                })
                .unwrap();
            let quality = MeshOptions {
                linear_deflection_mm: 0.02,
                angular_deflection_rad: 0.3,
                maximum_triangles: 500_000,
            };
            for source in [body, &moved] {
                let mfg = source.manufacturing_mesh_with_receipt(quality).unwrap();
                assert!(mfg.edge_pairs > 0 && mfg.vertex_links > 0);
                assert_eq!(
                    source
                        .rust_internal_thread_treatment_data()
                        .unwrap()
                        .material,
                    source_before
                );
            }
            assert_eq!(
                body.manufacturing_mesh(MeshOptions {
                    maximum_triangles: 16,
                    ..quality
                })
                .err()
                .unwrap()
                .kind,
                GeometryErrorKind::Limit
            );
            let dir = tempfile::tempdir().unwrap();
            let destination = dir.path().join("tail.step");
            std::fs::write(&destination, b"preserve original manufacturing file").unwrap();
            assert!(body.write_step(&destination).is_err());
            assert_eq!(
                std::fs::read(&destination).unwrap(),
                b"preserve original manufacturing file"
            );
            assert_eq!(std::fs::read_dir(dir.path()).unwrap().count(), 1);
            let original = model::read_thread(session.document(), thread.feature)
                .unwrap()
                .parameters;
            let mut edited = original.clone();
            let ends = edited.end_treatments.as_mut().unwrap();
            assert_eq!(ends.start, ThreadStartTreatment::FlatCut);
            ends.end = ThreadEndTreatment::CylindricalRelief {
                length: length(0.35),
                radial_depth: length(0.075),
            };
            let actor = || Actor::local(Channel::Ui);
            session
                .commit(
                    model::prepare_thread_edit(
                        &session,
                        session.revision(),
                        actor(),
                        thread.feature,
                        &edited,
                    )
                    .unwrap()
                    .prepared,
                )
                .unwrap();
            let changed = compile(&session);
            let changed_body = changed.solid(thread.feature).unwrap();
            let changed_data = changed_body.rust_internal_thread_treatment_data().unwrap();
            assert_eq!(changed_data.material.stock(), data.material.stock());
            assert_eq!(
                changed_data.material.cavity().source_thread(),
                data.material.cavity().source_thread()
            );
            assert!(
                changed_data
                    .material
                    .volume_bounds(Default::default())
                    .unwrap()
                    .upper
                    < before_volume.lower
            );
            changed_body
                .manufacturing_mesh_with_receipt(quality)
                .unwrap();
            let mut clear = edited.clone();
            clear.end_treatments = None;
            session
                .commit(
                    model::prepare_thread_edit(
                        &session,
                        session.revision(),
                        actor(),
                        thread.feature,
                        &clear,
                    )
                    .unwrap()
                    .prepared,
                )
                .unwrap();
            let cleared = compile(&session);
            let plain = cleared.solid(thread.feature).unwrap();
            assert!(plain.rust_internal_thread_treatment_data().is_none());
            assert!(
                plain
                    .rust_thread_contains([
                        1.525,
                        0.,
                        data.material.cavity().source_thread().full_length_mm() - 0.1
                    ])
                    .unwrap()
            );
            session
                .commit(
                    model::prepare_thread_edit(
                        &session,
                        session.revision(),
                        actor(),
                        thread.feature,
                        &edited,
                    )
                    .unwrap()
                    .prepared,
                )
                .unwrap();
            let restored = compile(&session);
            assert_eq!(
                restored
                    .solid(thread.feature)
                    .unwrap()
                    .rust_internal_thread_treatment_data()
                    .unwrap()
                    .material,
                changed_data.material
            );
            let expected = changed_data.material;
            let mut bytes = Vec::new();
            project_file::write(
                &ProjectFile::new(session.document().clone()).unwrap(),
                &mut bytes,
                Default::default(),
            )
            .unwrap();
            drop(moved);
            drop(built);
            drop(changed);
            drop(cleared);
            drop(restored);
            let document = project_file::read(
                &mut std::io::Cursor::new(bytes),
                Default::default(),
                &project_file::supported_capabilities(),
            )
            .unwrap()
            .project
            .document;
            let cold = Session::new(document).unwrap();
            let reopened = compile(&cold);
            let actual = reopened.solid(thread.feature).unwrap();
            assert_eq!(
                actual
                    .rust_internal_thread_treatment_data()
                    .unwrap()
                    .material,
                expected
            );
            actual.manufacturing_mesh_with_receipt(quality).unwrap();
        }
    }
}
