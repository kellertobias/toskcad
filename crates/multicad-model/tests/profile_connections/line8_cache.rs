use super::line8_nuts::{compiler, fixture};
use super::*;
use multicad_model::EvaluationCache;
#[test]
fn rust_line8_normal_viewport_build_reuses_both_nuts_on_the_next_compile() {
    let (session, _members, _request, outputs) = fixture(2, ProfileConnectionNut::Stainless);
    let mut cache = multicad_model::EvaluationCache::default();
    let compiler = Compiler::default();
    let start = std::time::Instant::now();
    let cold = compiler
        .compile_cached_with_backend(
            &mut cache,
            session.document(),
            session.revision(),
            GeometryBackend::Rust,
            None,
        )
        .unwrap();
    let cold_time = start.elapsed();
    for i in [2, 4] {
        let shape = &cold.snapshot().shapes[&outputs[i].feature];
        assert!(!shape.mesh.indices.is_empty());
        assert_eq!(shape.mesh.face_indices.len() * 3, shape.mesh.indices.len());
    }
    let start = std::time::Instant::now();
    let hot = compiler
        .compile_cached_with_backend(
            &mut cache,
            session.document(),
            session.revision() + 1,
            GeometryBackend::Rust,
            None,
        )
        .unwrap();
    eprintln!(
        "Line8 normal viewport cold={cold_time:?}, cached={:?}",
        start.elapsed()
    );
    for i in [2, 4] {
        let feature = outputs[i].feature;
        assert!(cache.last_compile().reused.contains(&feature));
        assert!(!cache.last_compile().evaluated.contains(&feature));
        assert!(std::sync::Arc::ptr_eq(
            &cold.snapshot().shapes[&feature].mesh,
            &hot.snapshot().shapes[&feature].mesh
        ));
        assert_eq!(
            cold.solid(feature).unwrap().volume_mm3().unwrap(),
            hot.solid(feature).unwrap().volume_mm3().unwrap()
        );
    }
}

fn cached(session: &Session, cache: &mut EvaluationCache) -> Model {
    compiler()
        .compile_cached_with_backend(
            cache,
            session.document(),
            session.revision(),
            GeometryBackend::Rust,
            None,
        )
        .unwrap()
}
fn nuts(outputs: &[model::CreatedBody]) -> [EntityId; 2] {
    [outputs[2].feature, outputs[4].feature]
}
fn assert_reused(cache: &EvaluationCache, ids: [EntityId; 2]) {
    for id in ids {
        assert!(cache.last_compile().reused.contains(&id));
        assert!(!cache.last_compile().evaluated.contains(&id));
    }
}
fn assert_evaluated(cache: &EvaluationCache, ids: [EntityId; 2]) {
    for id in ids {
        assert!(cache.last_compile().evaluated.contains(&id));
        assert!(!cache.last_compile().reused.contains(&id));
    }
}
fn parameter(name: &str, value: &str) -> Entity {
    Entity::new(
        name,
        EntityData::Parameter {
            expression: value.into(),
            unit: "mm".into(),
            exposed: true,
            inputs: vec![],
        },
    )
}
fn apply(session: &mut Session, label: &str, edits: Vec<Mutation>) {
    session
        .apply(Request {
            expected_revision: session.revision(),
            actor: actor(),
            label: label.into(),
            edits,
        })
        .unwrap();
}
#[test]
fn rust_line8_cache_tracks_bound_parameters_and_actual_member_source_but_not_unrelated_parameters()
{
    use multicad_core::{constraint_sketch::Dimension, parameters::BoundExpression};
    let (mut session, members, mut request, outputs) = fixture(0, ProfileConnectionNut::ZincPlated);
    let station = parameter("mount_station", "=40 mm");
    let station_id = station.id;
    apply(
        &mut session,
        "Add station parameter",
        vec![Mutation::Insert(station)],
    );
    let bound = Dimension::Expression {
        expression: BoundExpression::bind(session.document(), "=$mount_station", "mm", None)
            .unwrap(),
    };
    for site in [&mut request.first, &mut request.second] {
        site.anchor = ProfileConnectorAnchor::Groove {
            side: ProfileGrooveSide::MaxY,
            station: bound.clone(),
        };
    }
    session
        .commit(
            model::prepare_profile_connection(&session, session.revision(), actor(), &request)
                .unwrap()
                .prepared,
        )
        .unwrap();
    let mut cache = EvaluationCache::default();
    let original = cached(&session, &mut cache);
    cached(&session, &mut cache);
    assert_reused(&cache, nuts(&outputs));
    apply(
        &mut session,
        "Unrelated parameter",
        vec![Mutation::Insert(parameter("unrelated", "=7 mm"))],
    );
    cached(&session, &mut cache);
    assert_reused(&cache, nuts(&outputs));
    let mut edit = session.document().get(station_id).unwrap().clone();
    let EntityData::Parameter { expression, .. } = &mut edit.data else {
        unreachable!()
    };
    *expression = "=45 mm".into();
    apply(
        &mut session,
        "Move actual bound station",
        vec![Mutation::Replace(edit)],
    );
    let moved = cached(&session, &mut cache);
    assert_evaluated(&cache, nuts(&outputs));
    for id in nuts(&outputs) {
        assert_ne!(
            moved.solid(id).unwrap().bounds().unwrap(),
            original.solid(id).unwrap().bounds().unwrap()
        );
    }
    let mut member = members[0].member.clone();
    member.path = MemberPath::Points {
        start: std::array::from_fn(|_| length(0.)),
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
    cached(&session, &mut cache);
    assert_evaluated(&cache, nuts(&outputs));
}
#[test]
fn rust_line8_cache_rejects_changed_component_binding_and_bad_pose_then_reuses_after_undo() {
    let (mut session, members, request, outputs) = fixture(2, ProfileConnectionNut::ZincPlated);
    let mut cache = EvaluationCache::default();
    cached(&session, &mut cache);
    cached(&session, &mut cache);
    assert_reused(&cache, nuts(&outputs));
    let original = session.document().clone();
    let mut instance = session
        .document()
        .get(request.first.occurrence.root)
        .unwrap()
        .clone();
    let EntityData::Instance { component, .. } = &mut instance.data else {
        unreachable!()
    };
    *component = members[1].member.output.body;
    apply(
        &mut session,
        "Replace occurrence component",
        vec![Mutation::Replace(instance)],
    );
    let changed = session.document().clone();
    assert!(
        compiler()
            .compile_sources_cached_with_backend(
                &mut cache,
                &changed,
                session.revision(),
                &nuts(&outputs),
                GeometryBackend::Rust,
                None
            )
            .is_err()
    );
    assert_eq!(session.document(), &changed);
    session.undo(session.revision(), actor()).unwrap();
    assert_eq!(session.document(), &original);
    cached(&session, &mut cache);
    assert_reused(&cache, nuts(&outputs));
    let move_ = model::InstanceMoveRequest::new(
        request.second.occurrence.root,
        model::NumericTransform {
            translation: [length(0.5), length(0.), length(0.)],
            ..Default::default()
        },
    );
    session
        .commit(
            model::prepare_instance_moves(&session, session.revision(), actor(), &[move_])
                .unwrap()
                .prepared,
        )
        .unwrap();
    let changed = session.document().clone();
    assert!(
        compiler()
            .compile_sources_cached_with_backend(
                &mut cache,
                &changed,
                session.revision(),
                &nuts(&outputs),
                GeometryBackend::Rust,
                None
            )
            .is_err()
    );
    assert_eq!(session.document(), &changed);
    session.undo(session.revision(), actor()).unwrap();
    cached(&session, &mut cache);
    assert_reused(&cache, nuts(&outputs));
}
#[test]
fn rust_line8_cache_includes_common_parent_absolute_pose_and_membership() {
    let (mut session, members, mut request, outputs) = fixture(0, ProfileConnectionNut::ZincPlated);
    let mut instances = Vec::new();
    for member in &members {
        let instance =
            model::ComponentPlacementRequest::new(member.member.output.body, Default::default());
        session
            .commit(
                model::prepare_place_component(&session, session.revision(), actor(), &instance)
                    .unwrap()
                    .prepared,
            )
            .unwrap();
        instances.push(instance.instance);
    }
    let group = model::AssemblyRequest::new(
        "Receivers",
        instances.clone(),
        members.iter().map(|m| m.member.output.body).collect(),
    );
    session
        .commit(
            model::prepare_assembly(&session, session.revision(), actor(), &group)
                .unwrap()
                .prepared,
        )
        .unwrap();
    let parent = model::ComponentPlacementRequest::new(
        group.assembly,
        model::NumericTransform {
            translation: [length(11.), length(-7.), length(5.)],
            angle: model::angle(0.37),
            ..Default::default()
        },
    );
    session
        .commit(
            model::prepare_place_component(&session, session.revision(), actor(), &parent)
                .unwrap()
                .prepared,
        )
        .unwrap();
    for (site, instance) in [&mut request.first, &mut request.second]
        .into_iter()
        .zip(&instances)
    {
        site.occurrence = OccurrencePath {
            root: parent.instance,
            members: vec![*instance],
        };
    }
    session
        .commit(
            model::prepare_profile_connection(&session, session.revision(), actor(), &request)
                .unwrap()
                .prepared,
        )
        .unwrap();
    let mut cache = EvaluationCache::default();
    let old = cached(&session, &mut cache);
    cached(&session, &mut cache);
    assert_reused(&cache, nuts(&outputs));
    let move_ = model::InstanceMoveRequest::new(
        parent.instance,
        model::NumericTransform {
            translation: [length(2.), length(3.), length(4.)],
            ..Default::default()
        },
    );
    session
        .commit(
            model::prepare_instance_moves(&session, session.revision(), actor(), &[move_])
                .unwrap()
                .prepared,
        )
        .unwrap();
    let changed = cached(&session, &mut cache);
    assert_evaluated(&cache, nuts(&outputs));
    for id in nuts(&outputs) {
        assert_ne!(
            old.solid(id).unwrap().bounds().unwrap(),
            changed.solid(id).unwrap().bounds().unwrap()
        );
    }
    let mut reordered = group.clone();
    reordered.members.reverse();
    session
        .commit(
            model::prepare_assembly_edit(&session, session.revision(), actor(), &reordered)
                .unwrap()
                .prepared,
        )
        .unwrap();
    cached(&session, &mut cache);
    assert_evaluated(&cache, nuts(&outputs));
}

fn align_instances_to_geometry(
    session: &mut Session,
    request: &ProfileConnectionRequest,
    input: model::AlignmentInput,
) {
    let align = [&request.first, &request.second].map(|site| {
        model::InstanceAlignmentRequest::new(
            site.occurrence.root,
            model::AlignmentParameters {
                target: model::AlignmentAnchor::Geometry {
                    reference: input.intent.id,
                },
                mode: model::AlignmentMode::Position,
                ..Default::default()
            },
            vec![input.clone()],
        )
    });
    session
        .commit(
            model::prepare_instance_alignments(session, session.revision(), actor(), &align)
                .unwrap()
                .prepared,
        )
        .unwrap();
}
fn support_capture(model: &Model, source: EntityId) -> model::AlignmentInput {
    use multicad_core::references::{AnalyticGeometry, Attachment, GeometryIntent, TopologyKind};
    let snapshot = &model.snapshot().shapes[&source].topology;
    let face=snapshot.elements.iter().find(|element|element.key.kind==TopologyKind::Face && matches!(element.geometry,Some(AnalyticGeometry::Plane{frame})if frame.normal()[2]>0.9) && !element.frames.is_empty()).unwrap();
    let mut intent = GeometryIntent::capture(snapshot, face.key).unwrap();
    intent.attachment = Some(Attachment {
        anchor: face.frames[0].origin.clone(),
        offset_mm: [0.; 3],
        angle_rad: 0.,
    });
    model::AlignmentInput { source, intent }
}
#[test]
fn rust_line8_cache_rechecks_actual_geometry_support_and_reference_binding() {
    let (mut session, _members, request, outputs) = fixture(2, ProfileConnectionNut::ZincPlated);
    let mut support = model::PrimitiveRequest::new(
        "Actual alignment support",
        model::PrimitiveParameters::Box {
            size: [length(10.), length(20.), length(30.)],
            frame: model::PrimitiveFrame {
                origin: [length(100.), length(200.), length(300.)],
                ..Default::default()
            },
        },
    );
    session
        .commit(
            model::prepare_primitive(&session, session.revision(), actor(), &support)
                .unwrap()
                .prepared,
        )
        .unwrap();
    let source = compiler()
        .compile_with_backend(
            session.document(),
            session.revision(),
            GeometryBackend::Rust,
            None,
        )
        .unwrap();
    align_instances_to_geometry(
        &mut session,
        &request,
        support_capture(&source, support.feature),
    );
    let mut cache = EvaluationCache::default();
    let old = cached(&session, &mut cache);
    cached(&session, &mut cache);
    assert_reused(&cache, nuts(&outputs));
    let model::PrimitiveParameters::Box { size, .. } = &mut support.parameters else {
        unreachable!()
    };
    size[2] = length(40.);
    session
        .commit(
            model::prepare_primitive_edit(
                &session,
                session.revision(),
                actor(),
                support.feature,
                &support.parameters,
            )
            .unwrap()
            .prepared,
        )
        .unwrap();
    let changed = cached(&session, &mut cache);
    assert_evaluated(&cache, nuts(&outputs));
    for id in nuts(&outputs) {
        assert_ne!(
            old.solid(id).unwrap().bounds().unwrap(),
            changed.solid(id).unwrap().bounds().unwrap()
        );
    }
    cached(&session, &mut cache);
    assert_reused(&cache, nuts(&outputs));
    // A new capture is a new explicit reference binding, even when it resolves
    // to the same unchanged supporting face and the resulting pose is equal.
    align_instances_to_geometry(
        &mut session,
        &request,
        support_capture(&changed, support.feature),
    );
    let rebound = cached(&session, &mut cache);
    assert_evaluated(&cache, nuts(&outputs));
    for id in nuts(&outputs) {
        assert_eq!(
            changed.solid(id).unwrap().bounds().unwrap(),
            rebound.solid(id).unwrap().bounds().unwrap()
        );
    }
}
