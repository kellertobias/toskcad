use super::*;
use multicad_kernel::RigidTransform;
pub(super) fn fixture(
    case: usize,
    variant: ProfileConnectionNut,
) -> (
    Session,
    [MemberRequest; 2],
    ProfileConnectionRequest,
    Vec<model::CreatedBody>,
) {
    let posed = |[x, y, z]: [f64; 3]| {
        if case == 1 {
            [11. + 0.8 * x - 0.6 * y, -7. + 0.6 * x + 0.8 * y, 5. + z]
        } else {
            [x, y, z]
        }
    };
    let (mut session, members) = prepare_members([
        (posed([0., 0., 0.]), posed([100., 0., 0.])),
        (posed([0., 40., 0.]), posed([100., 40., 0.])),
    ]);
    let mut sites = [
        groove(&members[0], ProfileGrooveSide::MaxY, 40.),
        groove(&members[1], ProfileGrooveSide::MaxY, 40.),
    ];
    if case == 2 {
        for i in 0..2 {
            let placement = model::ComponentPlacementRequest::new(
                members[i].member.output.body,
                model::NumericTransform {
                    translation: [length(11.), length(-7.), length(5.)],
                    axis: [0.3, 1., 0.4],
                    angle: model::angle(0.37),
                    ..Default::default()
                },
            );
            session
                .commit(
                    model::prepare_place_component(
                        &session,
                        session.revision(),
                        actor(),
                        &placement,
                    )
                    .unwrap()
                    .prepared,
                )
                .unwrap();
            sites[i].occurrence = OccurrencePath {
                root: placement.instance,
                members: vec![],
            };
        }
    }
    let mut request = ProfileConnectionRequest::new(
        "Modeled Rust Line8 nuts",
        ProfileConnectionKind::JoiningPlate,
        sites[0].clone(),
        sites[1].clone(),
    );
    request.nut = variant;
    let draft =
        model::prepare_profile_connection(&session, session.revision(), actor(), &request).unwrap();
    let outputs = draft.outputs.clone();
    session.commit(draft.prepared).unwrap();
    (session, members, request, outputs)
}
pub(super) fn compiler() -> Compiler {
    Compiler {
        omit_display: true,
        ..Default::default()
    }
}
#[test]
fn rust_line8_nuts_compile_whole_projects_reopen_and_manufacture_in_actual_occurrences() {
    for (case, variant) in [
        (0, ProfileConnectionNut::ZincPlated),
        (1, ProfileConnectionNut::Stainless),
        (2, ProfileConnectionNut::ZincPlated),
    ] {
        let (session, members, request, outputs) = fixture(case, variant);
        let before = session.document().clone();
        let revision = session.revision();
        let mut bytes = Vec::new();
        project_file::write(
            &ProjectFile::new(before.clone()).unwrap(),
            &mut bytes,
            Default::default(),
        )
        .unwrap();
        let reopened = project_file::read(
            &mut bytes.as_slice(),
            Default::default(),
            &project_file::supported_capabilities(),
        )
        .unwrap();
        assert_eq!(reopened.project.document, before);
        let mut original = Vec::new();
        for (round, document) in [&before, &reopened.project.document]
            .into_iter()
            .enumerate()
        {
            let started = std::time::Instant::now();
            let model = compiler()
                .compile_with_backend(document, revision, GeometryBackend::Rust, None)
                .unwrap_or_else(|e| panic!("case {case}, round {round}: {e}"));
            eprintln!(
                "Line8 whole project case={case} round={round}: {:?}",
                started.elapsed()
            );
            let parameters = multicad_core::parameters::evaluate(document).unwrap();
            let snapshots = model
                .snapshot()
                .shapes
                .iter()
                .map(|(id, shape)| (*id, shape.topology.clone()))
                .collect();
            for (nut_index, screw_index) in [(2, 1), (4, 3)] {
                let output = &outputs[nut_index];
                let solid = model.solid(output.feature).unwrap();
                let source = solid.rust_line8_nut_source().unwrap();
                assert_eq!(source.thread().spec().length_mm, 4.1);
                assert_eq!(source.relief().start_mm, 4.1 - 0.001);
                assert_eq!(source.bore_from_start_mm(), 7.5);
                assert_eq!(source.stock().height_mm(), 23.);
                assert_eq!(source.stock().base()[0][0], [-6.9, 0.]);
                assert_eq!(solid.rust_line8_nut_requested_wall_mm(), Some(0.5));
                assert!(solid.uses_rust_geometry());
                assert_eq!(solid.solid_count().unwrap(), 1);
                assert_eq!(
                    model.snapshot().shapes[&output.feature].producer,
                    output.feature
                );
                assert_eq!(model.snapshot().bodies[&output.body], output.feature);
                solid
                    .topology_snapshot(output.feature, revision)
                    .unwrap()
                    .validate()
                    .unwrap();
                let screw = model.solid(outputs[screw_index].feature).unwrap();
                let mating = screw.rust_line8_nut_mating(solid, None).unwrap();
                assert!(mating.complete_turns_lower >= 2.);
                assert!(mating.relief_clearance_lower_mm > 0.1);
                for (site, member) in [
                    (&request.first, &members[0]),
                    (&request.second, &members[1]),
                ] {
                    let pose = model::resolve_occurrence_pose(
                        document,
                        &site.occurrence,
                        site.member.body,
                        &parameters,
                        revision,
                        &snapshots,
                    )
                    .unwrap();
                    let receipt = solid
                        .rust_line8_receiver_clearance(
                            model.solid(member.feature).unwrap(),
                            RigidTransform {
                                rotation: pose.rotation,
                                translation_mm: pose.translation_mm,
                            },
                            0.001,
                            None,
                        )
                        .unwrap();
                    assert!(receipt.overlap_volume_upper_mm3 <= 0.001);
                }
                let current = (solid.volume_mm3().unwrap(), solid.bounds().unwrap());
                if round == 0 {
                    original.push(current);
                } else {
                    assert_eq!(current, original[nut_index / 2 - 1]);
                    let receipt = solid
                        .manufacturing_mesh_with_receipt(Default::default())
                        .unwrap_or_else(|e| {
                            panic!("manufacturing case {case}, nut {nut_index}: {e}")
                        });
                    assert!(!receipt.mesh.indices.is_empty());
                    assert_eq!(receipt.edge_pairs as usize * 2, receipt.mesh.indices.len());
                    assert_eq!(
                        receipt.vertex_links as usize,
                        receipt.mesh.vertices.len() / 3
                    );
                    assert_eq!(receipt.maximum_displacement_mm, 0.);
                }
            }
        }
        assert_eq!(session.document(), &before);
    }
}
#[test]
fn rust_line8_wrong_owners_and_changed_occurrence_fail_atomically_and_undo_recovers() {
    let (mut session, _members, request, outputs) = fixture(2, ProfileConnectionNut::ZincPlated);
    let original = session.document().clone();
    let mut wrong = request.clone();
    wrong.first.member.producer = request.second.member.producer;
    assert!(
        model::prepare_profile_connection(&session, session.revision(), actor(), &wrong).is_err()
    );
    assert_eq!(session.document(), &original);
    let moved = model::InstanceMoveRequest::new(
        request.second.occurrence.root,
        model::NumericTransform {
            translation: [length(0.5), length(0.), length(0.)],
            ..Default::default()
        },
    );
    session
        .commit(
            model::prepare_instance_moves(&session, session.revision(), actor(), &[moved])
                .unwrap()
                .prepared,
        )
        .unwrap();
    let changed = session.document().clone();
    let error = compiler()
        .compile_sources_with_backend(
            &changed,
            session.revision(),
            &[outputs[2].feature],
            GeometryBackend::Rust,
            None,
        )
        .err()
        .unwrap();
    assert!(
        error.to_string().contains("grooves")
            || error.to_string().contains("station")
            || error.to_string().contains("pose"),
        "{error}"
    );
    assert_eq!(session.document(), &changed);
    session.undo(session.revision(), actor()).unwrap();
    assert_eq!(session.document(), &original);
    compiler()
        .compile_sources_with_backend(
            session.document(),
            session.revision(),
            &[outputs[2].feature],
            GeometryBackend::Rust,
            None,
        )
        .unwrap();
    assert!(
        compiler()
            .compile_sources_with_backend(
                session.document(),
                session.revision(),
                &[outputs[2].feature],
                GeometryBackend::Rust,
                Some(&|| true)
            )
            .is_err()
    );
    assert_eq!(session.document(), &original);
}

#[test]
fn rust_line8_normal_viewport_builds_both_actual_nut_meshes() {
    let (session, _members, _request, outputs) = fixture(2, ProfileConnectionNut::Stainless);
    let start = std::time::Instant::now();
    let model = Compiler::default()
        .compile_with_backend(
            session.document(),
            session.revision(),
            GeometryBackend::Rust,
            None,
        )
        .unwrap();
    eprintln!("Line8 normal viewport build: {:?}", start.elapsed());
    for i in [2, 4] {
        let shape = &model.snapshot().shapes[&outputs[i].feature];
        assert!(!shape.mesh.indices.is_empty());
        assert_eq!(shape.mesh.face_indices.len() * 3, shape.mesh.indices.len());
        assert!(
            model
                .solid(outputs[i].feature)
                .unwrap()
                .rust_line8_nut_source()
                .is_some()
        );
    }
}
