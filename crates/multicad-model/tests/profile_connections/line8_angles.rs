use super::*;
use multicad_kernel::RigidTransform;
fn fixture(
    case: usize,
) -> (
    Session,
    [MemberRequest; 2],
    ProfileConnectionRequest,
    Vec<model::CreatedBody>,
) {
    let pose = |[x, y, z]: [f64; 3]| {
        if case == 1 {
            [11. + 0.8 * x - 0.6 * y, -7. + 0.6 * x + 0.8 * y, 5. + z]
        } else {
            [x, y, z]
        }
    };
    let (mut session, mut members) = prepare_members([
        (pose([0., 0., 0.]), pose([100., 0., 0.])),
        (pose([-20., 0., 20.]), pose([-20., 0., 120.])),
    ]);
    if case == 1 {
        // Rotating only a vertical path leaves its section roll unchanged.
        // Author the actual section rotation too, so BOTH member solids have
        // the intended common physical pose rather than just moved endpoints.
        let mut member = members[1].member.clone();
        member.rotation = Some(model::angle(0.6_f64.atan2(0.8).to_degrees()));
        session
            .commit(
                model::prepare_member_edit(
                    &session,
                    session.revision(),
                    actor(),
                    members[1].feature,
                    &member,
                )
                .unwrap()
                .prepared,
            )
            .unwrap();
        members[1].member = member;
    }
    let mut sites = [
        groove(&members[0], ProfileGrooveSide::MaxY, 20.),
        groove(&members[1], ProfileGrooveSide::MaxX, 20.),
    ];
    if case == 2 {
        for i in 0..2 {
            let placement = model::ComponentPlacementRequest::new(
                members[i].member.output.body,
                model::NumericTransform {
                    translation: [length(11.), length(-7.), length(5.)],
                    axis: [0.3, 1., 0.4],
                    angle: model::angle(0.37_f64.to_degrees()),
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
    let request = ProfileConnectionRequest::new(
        "Rust retained angle bracket",
        ProfileConnectionKind::AngleBracket,
        sites[0].clone(),
        sites[1].clone(),
    );
    let draft =
        model::prepare_profile_connection(&session, session.revision(), actor(), &request).unwrap();
    let outputs = draft.outputs.clone();
    session.commit(draft.prepared).unwrap();
    (session, members, request, outputs)
}
#[test]
fn rust_actual_angle_projects_reopen_with_all_five_physical_parts_and_both_receivers() {
    for case in 0..3 {
        let (session, members, request, outputs) = fixture(case);
        let original = session.document().clone();
        let revision = session.revision();
        let mut bytes = Vec::new();
        project_file::write(
            &ProjectFile::new(original.clone()).unwrap(),
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
        assert_eq!(reopened.project.document, original);
        for (round, document) in [&original, &reopened.project.document]
            .into_iter()
            .enumerate()
        {
            let now = std::time::Instant::now();
            let compiled = super::line8_nuts::compiler()
                .compile_with_backend(document, revision, GeometryBackend::Rust, None)
                .unwrap_or_else(|e| panic!("Angle case={case} round={round}: {e}"));
            eprintln!(
                "Actual angle project case={case} round={round}: {:?}",
                now.elapsed()
            );
            assert_eq!(
                compiled.snapshot().leaf_occurrences.len(),
                if case == 2 { 9 } else { 7 }
            );
            let parameters = multicad_core::parameters::evaluate(document).unwrap();
            let snapshots = compiled
                .snapshot()
                .shapes
                .iter()
                .map(|(id, s)| (*id, s.topology.clone()))
                .collect();
            let bracket = compiled.solid(outputs[0].feature).unwrap();
            let (source, _) = bracket.rust_orthogonal_countersink_source().unwrap();
            assert_eq!(source.bores().len(), 2);
            assert_eq!(source.seats().len(), 2);
            assert!(source.bores().iter().all(|b| b.radius_mm == 4.9));
            assert!(
                source.seats().iter().all(
                    |s| (s.radius_mm - 8.02).abs() < 1e-12 && (s.depth_mm - 3.12).abs() < 1e-12
                )
            );
            for (index, output) in outputs.iter().enumerate() {
                let shape = compiled.solid(output.feature).unwrap();
                assert_eq!(shape.solid_count().unwrap(), 1);
                assert_eq!(
                    compiled.snapshot().shapes[&output.feature].producer,
                    output.feature
                );
                assert_eq!(compiled.snapshot().bodies[&output.body], output.feature);
                shape
                    .topology_snapshot(output.feature, revision)
                    .unwrap()
                    .validate()
                    .unwrap();
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
                    let receiver = compiled.solid(member.feature).unwrap();
                    let transform = RigidTransform {
                        rotation: pose.rotation,
                        translation_mm: pose.translation_mm,
                    };
                    let upper = match index {
                        0 => {
                            shape
                                .rust_seated_stock_receiver_clearance(
                                    receiver, transform, 0.001, None,
                                )
                                .unwrap()
                                .overlap_volume_upper_mm3
                        }
                        1 | 3 => {
                            shape
                                .rust_profiled_screw_receiver_clearance(
                                    receiver, transform, 0.001, None,
                                )
                                .unwrap()
                                .overlap_volume_upper_mm3
                        }
                        2 | 4 => {
                            shape
                                .rust_line8_receiver_clearance(receiver, transform, 0.001, None)
                                .unwrap()
                                .overlap_volume_upper_mm3
                        }
                        _ => unreachable!(),
                    };
                    assert!(upper <= 0.001);
                }
                if round == 1 {
                    let proof = shape
                        .manufacturing_mesh_with_receipt(Default::default())
                        .unwrap_or_else(|e| panic!("Angle MFG case={case} part={index}: {e}"));
                    assert_eq!(proof.edge_pairs as usize * 2, proof.mesh.indices.len());
                    assert_eq!(proof.vertex_links as usize, proof.mesh.vertices.len() / 3);
                    assert_eq!(proof.maximum_displacement_mm, 0.);
                }
            }
            for (screw, nut) in [(1, 2), (3, 4)] {
                let screw = compiled.solid(outputs[screw].feature).unwrap();
                let nut = compiled.solid(outputs[nut].feature).unwrap();
                assert_eq!(
                    nut.rust_line8_nut_source()
                        .unwrap()
                        .thread()
                        .spec()
                        .length_mm,
                    4.1
                );
                assert_eq!(
                    nut.rust_line8_nut_source().unwrap().relief().start_mm,
                    4.1 - 0.001
                );
                assert!(
                    screw
                        .rust_line8_nut_mating(nut, None)
                        .unwrap()
                        .complete_turns_lower
                        >= 2.
                );
            }
        }
        assert_eq!(session.document(), &original);
    }
}
#[test]
fn rust_angle_warm_viewport_reuses_all_hardware_and_invalidates_both_actual_pose_paths() {
    let (mut session, members, request, outputs) = fixture(2);
    let mut cache = multicad_model::EvaluationCache::default();
    let compiler = Compiler::default();
    let compile = |session: &Session, cache: &mut multicad_model::EvaluationCache| {
        compiler.compile_cached_with_backend(
            cache,
            session.document(),
            session.revision(),
            GeometryBackend::Rust,
            None,
        )
    };
    let cold = compile(&session, &mut cache).unwrap();
    let now = std::time::Instant::now();
    let hot = compiler
        .compile_cached_with_backend(
            &mut cache,
            session.document(),
            session.revision() + 1,
            GeometryBackend::Rust,
            None,
        )
        .unwrap();
    eprintln!("Actual Angle warm whole viewport: {:?}", now.elapsed());
    for o in &outputs {
        assert!(cache.last_compile().reused.contains(&o.feature));
        assert!(!cache.last_compile().evaluated.contains(&o.feature));
        assert!(std::sync::Arc::ptr_eq(
            &cold.snapshot().shapes[&o.feature].mesh,
            &hot.snapshot().shapes[&o.feature].mesh
        ));
    }
    let moves = [
        request.first.occurrence.root,
        request.second.occurrence.root,
    ]
    .map(|id| {
        model::InstanceMoveRequest::new(
            id,
            model::NumericTransform {
                translation: [length(1.), length(0.), length(0.)],
                ..Default::default()
            },
        )
    });
    session
        .commit(
            model::prepare_instance_moves(&session, session.revision(), actor(), &moves)
                .unwrap()
                .prepared,
        )
        .unwrap();
    let moved = compile(&session, &mut cache).unwrap();
    for o in &outputs {
        assert!(cache.last_compile().evaluated.contains(&o.feature));
        assert!(!cache.last_compile().reused.contains(&o.feature));
        assert_ne!(
            cold.solid(o.feature).unwrap().bounds().unwrap(),
            moved.solid(o.feature).unwrap().bounds().unwrap()
        );
    }
    let mut member = members[0].member.clone();
    member.path = MemberPath::Points {
        start: [length(0.), length(0.), length(0.)],
        end: [length(120.), length(0.), length(0.)],
    };
    session
        .commit(
            model::prepare_member_edit(
                &session,
                session.revision(),
                actor(),
                members[0].feature,
                &member,
            )
            .unwrap()
            .prepared,
        )
        .unwrap();
    let edited = compile(&session, &mut cache).unwrap();
    for o in &outputs {
        assert!(cache.last_compile().evaluated.contains(&o.feature));
        assert!(!cache.last_compile().reused.contains(&o.feature));
        assert_eq!(
            edited.solid(o.feature).unwrap().bounds().unwrap(),
            moved.solid(o.feature).unwrap().bounds().unwrap()
        );
    }
    let valid = session.document().clone();
    let bad = model::InstanceMoveRequest::new(
        request.second.occurrence.root,
        model::NumericTransform {
            translation: [length(0.5), length(0.), length(0.)],
            ..Default::default()
        },
    );
    session
        .commit(
            model::prepare_instance_moves(&session, session.revision(), actor(), &[bad])
                .unwrap()
                .prepared,
        )
        .unwrap();
    let before = session.document().clone();
    assert!(compile(&session, &mut cache).is_err());
    assert_eq!(session.document(), &before);
    session.undo(session.revision(), actor()).unwrap();
    assert_eq!(session.document(), &valid);
    let recovered = compile(&session, &mut cache).unwrap();
    for o in &outputs {
        assert_eq!(
            recovered.solid(o.feature).unwrap().bounds().unwrap(),
            moved.solid(o.feature).unwrap().bounds().unwrap()
        );
    }
}
