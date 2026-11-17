mod common;
use common::*;
use multicad_core::{
    commands::Session,
    constraint_sketch::{Dimension, source::Source},
    document::{Document, Entity, EntityData, Mutation},
    model::{
        self, BodyInput, Direction, Extent, ExtrudeParameters, Operation, PrimitiveFrame,
        PrimitiveParameters, PrimitiveRequest, RevolveParameters, RevolveRequest, Termination,
        angle, length,
    },
    parameters::BoundExpression,
    project_file::{self, ContentHash, ProjectFile},
    references::{Origin, TopologyKind},
    sketch::{self, CreationPlane, PlaneOptions},
    sketch_regions,
};
use multicad_kernel::Solid;
use multicad_model::{Compiler, GeometryBackend, Model, ShapeSnapshot};

fn shape(model: &Model, body: multicad_core::document::EntityId) -> &ShapeSnapshot {
    let producer = model.snapshot().bodies.get(&body).expect("Part producer");
    model
        .snapshot()
        .shapes
        .get(producer)
        .expect("Exact shape snapshot")
}
#[test]
fn one_shared_sketch_compiles_two_real_solids_without_mutating_canonical_source() {
    let fixture = SharedProfiles::new(Default::default());
    let before = fixture.session.document().clone();
    let hash = ContentHash::document(&before).unwrap();
    let model = Compiler::default()
        .compile(&before, fixture.session.revision())
        .unwrap();
    let snapshot = model.snapshot();
    assert_eq!(snapshot.stamp.document, before.id());
    assert_eq!(snapshot.stamp.revision, fixture.session.revision());
    assert_eq!(snapshot.stamp.source_hash, hash);
    assert_eq!(snapshot.bodies.len(), 2);
    assert_eq!(snapshot.shapes.len(), 2);
    let outer = shape(&model, fixture.outer().body);
    let inner = shape(&model, fixture.inner().body);
    near(outer.volume_mm3, (WIDTH * HEIGHT - INNER_AREA) * 10.);
    near(inner.volume_mm3, INNER_AREA * 10.);
    near(outer.volume_mm3 + inner.volume_mm3, WIDTH * HEIGHT * 10.);
    point_near(outer.bounds.min, [0., 0., 0.]);
    point_near(outer.bounds.max, [WIDTH, HEIGHT, 10.]);
    point_near(inner.bounds.min, [14., 14., 0.]);
    point_near(inner.bounds.max, [26., 26., 10.]);
    for output in &fixture.request.outputs {
        let snapshot = shape(&model, output.body);
        assert_eq!(snapshot.producer, output.feature);
        assert_eq!(snapshot.solid_count, 1);
        assert!(!snapshot.mesh.vertices.is_empty());
        assert!(!snapshot.mesh.indices.is_empty());
        assert!(!snapshot.topology.elements.is_empty());
        let exact = model.solid(output.feature).expect("Owned native B-rep");
        assert_eq!(exact.solid_count().unwrap(), 1);
        near(exact.volume_mm3().unwrap(), snapshot.volume_mm3);
    }
    assert_eq!(fixture.session.document(), &before);
    assert_eq!(
        ContentHash::document(fixture.session.document()).unwrap(),
        hash
    );
}
#[test]
fn sketch_edit_rebuilds_profiles_with_the_same_feature_and_body_ids() {
    let mut fixture = SharedProfiles::new(Default::default());
    let compiler = Compiler::default();
    let first = compiler
        .compile(fixture.session.document(), fixture.session.revision())
        .unwrap();
    let identities = first.snapshot().bodies.clone();
    let authored_ids: Vec<_> = fixture
        .source()
        .sketch
        .curves
        .iter()
        .map(|c| c.id)
        .collect();
    fixture.set_width(100.);
    let changed = fixture.session.document().clone();
    let second = compiler
        .compile(&changed, fixture.session.revision())
        .unwrap();
    assert_eq!(second.snapshot().bodies, identities);
    assert_eq!(
        fixture
            .source()
            .sketch
            .curves
            .iter()
            .map(|c| c.id)
            .collect::<Vec<_>>(),
        authored_ids
    );
    near(
        shape(&second, fixture.outer().body).volume_mm3,
        (100. * HEIGHT - INNER_AREA) * 10.,
    );
    near(
        shape(&second, fixture.inner().body).volume_mm3,
        INNER_AREA * 10.,
    );
    point_near(
        shape(&second, fixture.outer().body).bounds.max,
        [100., HEIGHT, 10.],
    );
    assert_ne!(first.snapshot().stamp, second.snapshot().stamp);
    assert_eq!(fixture.session.document(), &changed);
}
#[test]
fn changed_holes_and_invalid_boundary_intent_fail_without_adopting_a_different_region() {
    let mut fixture = SharedProfiles::new(Default::default());
    let compiler = Compiler::default();
    let previous = compiler
        .compile(fixture.session.document(), fixture.session.revision())
        .unwrap();
    let outer = fixture.outer().feature;
    fixture.add_unexpected_hole();
    let changed = fixture.session.document().clone();
    let error = compiler
        .compile(&changed, fixture.session.revision())
        .err()
        .expect("Stale selected holes must fail");
    assert_eq!(error.entity, Some(outer));
    assert!(error.message.contains("Selected profile changed"));
    assert_eq!(fixture.session.document(), &changed);
    // A caller may retain this old display result, but the failed compile has
    // not relabeled it as current or silently generated a newly holed part.
    assert_ne!(
        previous.snapshot().stamp.revision,
        fixture.session.revision()
    );

    let mut fixture = SharedProfiles::new(Default::default());
    let feature = fixture.outer().feature;
    let mut source = model::read_extrude(fixture.session.document(), feature).unwrap();
    let model::ModelProfile::Sketch(profile) = &mut source.profile else {
        panic!("Sketch profile")
    };
    profile.region.outer.curves = profile
        .region
        .outer
        .curves
        .iter()
        .copied()
        .take(1)
        .collect();
    profile.region.holes.clear();
    let mut entity = fixture.session.document().get(feature).unwrap().clone();
    let EntityData::Feature { definition, .. } = &mut entity.data else {
        panic!("Feature")
    };
    *definition = source.definition().unwrap();
    fixture.apply(vec![Mutation::Replace(entity)]);
    let error = compiler
        .compile(fixture.session.document(), fixture.session.revision())
        .err()
        .expect("One open line is not the selected area");
    assert_eq!(error.entity, Some(feature));
    assert!(error.message.contains("Selected profile changed"));
}
#[test]
fn offsets_tilt_and_symmetric_total_have_real_sweep_bounds_and_volumes() {
    let s = std::f64::consts::FRAC_1_SQRT_2;
    let cases = [
        (
            90.,
            false,
            Extent::OneSided {
                end: Termination::Distance {
                    distance: length(10.),
                },
            },
            [0., 0., 7.],
            [WIDTH, HEIGHT, 17.],
            10.,
        ),
        (
            45.,
            false,
            Extent::OneSided {
                end: Termination::Distance {
                    distance: length(10.),
                },
            },
            [0., 0., 7.],
            [WIDTH + 10. * s, HEIGHT, 7. + 10. * s],
            10. * s,
        ),
        (
            135.,
            false,
            Extent::OneSided {
                end: Termination::Distance {
                    distance: length(10.),
                },
            },
            [-10. * s, 0., 7.],
            [WIDTH, HEIGHT, 7. + 10. * s],
            10. * s,
        ),
        (
            45.,
            true,
            Extent::OneSided {
                end: Termination::Distance {
                    distance: length(10.),
                },
            },
            [-10. * s, 0., 7. - 10. * s],
            [WIDTH, HEIGHT, 7.],
            10. * s,
        ),
        // A signed negative one-sided distance sweeps the same magnitude in
        // the opposite direction, measured from the offset start plane.
        (
            90.,
            false,
            Extent::OneSided {
                end: Termination::Distance {
                    distance: length(-25.),
                },
            },
            [0., 0., -18.],
            [WIDTH, HEIGHT, 7.],
            25.,
        ),
        (
            45.,
            false,
            Extent::OneSided {
                end: Termination::Distance {
                    distance: length(-10.),
                },
            },
            [-10. * s, 0., 7. - 10. * s],
            [WIDTH, HEIGHT, 7.],
            10. * s,
        ),
        (
            45.,
            false,
            Extent::Symmetric { total: length(10.) },
            [-5. * s, 0., 7. - 5. * s],
            [WIDTH + 5. * s, HEIGHT, 7. + 5. * s],
            10. * s,
        ),
        (
            90.,
            false,
            Extent::TwoSided {
                forward: Termination::Distance {
                    distance: length(10.),
                },
                backward: Termination::Distance {
                    distance: length(4.),
                },
            },
            [0., 0., 3.],
            [WIDTH, HEIGHT, 17.],
            14.,
        ),
    ];
    for (elevation, reverse, extent, min, max, normal_height) in cases {
        let parameters = ExtrudeParameters {
            direction: Direction {
                elevation: angle(elevation),
                azimuth: angle(0.),
                reverse,
            },
            start_offset: length(7.),
            extent,
            ..Default::default()
        };
        let fixture = SharedProfiles::new(parameters);
        let before = fixture.session.document().clone();
        let model = Compiler::default()
            .compile(&before, fixture.session.revision())
            .unwrap();
        let outer = shape(&model, fixture.outer().body);
        point_near(outer.bounds.min, min);
        point_near(outer.bounds.max, max);
        near(
            outer.volume_mm3,
            (WIDTH * HEIGHT - INNER_AREA) * normal_height,
        );
        near(
            shape(&model, fixture.inner().body).volume_mm3,
            INNER_AREA * normal_height,
        );
        assert_eq!(fixture.session.document(), &before);
    }
}
#[test]
fn unsupported_extrusion_operation_is_an_error_instead_of_default_geometry() {
    let mut fixture = SharedProfiles::new(Default::default());
    let feature = fixture.outer().feature;
    let mut source = model::read_extrude(fixture.session.document(), feature).unwrap();
    source.operation = Operation::Cut;
    source.targets = vec![BodyInput {
        body: fixture.inner().body,
        producer: fixture.inner().feature,
    }];
    let mut entity = fixture.session.document().get(feature).unwrap().clone();
    let EntityData::Feature {
        definition, inputs, ..
    } = &mut entity.data
    else {
        panic!()
    };
    *definition = source.definition().unwrap();
    *inputs = source.inputs(&[]).unwrap();
    fixture.apply(vec![Mutation::Replace(entity)]);
    let before = fixture.session.document().clone();
    let error = Compiler::default()
        .compile(&before, fixture.session.revision())
        .err()
        .expect("Unsupported mode must not produce default extrusion");
    assert_eq!(error.entity, Some(feature));
    assert_eq!(fixture.session.document(), &before);
}
#[test]
fn exact_compiled_body_writes_and_reads_step_without_using_the_display_mesh() {
    let fixture = SharedProfiles::new(Default::default());
    let model = Compiler::default()
        .compile(fixture.session.document(), fixture.session.revision())
        .unwrap();
    let exact = model.solid(fixture.inner().feature).unwrap();
    let path = std::env::temp_dir().join(format!(
        "multicad-model-{}.step",
        multicad_core::document::EntityId::new()
    ));
    exact.write_step(&path).unwrap();
    let imported = Solid::read_step(&path, Default::default()).unwrap();
    std::fs::remove_file(&path).unwrap();
    near(imported.volume_mm3().unwrap(), INNER_AREA * 10.);
    assert_eq!(imported.solid_count().unwrap(), 1);
    let bounds = imported.bounds().unwrap();
    point_near(bounds.min, [14., 14., 0.]);
    point_near(bounds.max, [26., 26., 10.]);
}
#[test]
fn face_attached_child_body_follows_an_edited_parent_cap_without_rewriting_its_sketch() {
    let mut fixture = SharedProfiles::new(Default::default());
    let compiler = Compiler::default();
    let parent = fixture.outer().feature;
    let initial = compiler
        .compile(fixture.session.document(), fixture.session.revision())
        .unwrap();
    let snapshot = &initial.snapshot().shapes[&parent].topology;
    let cap = Origin {
        source: parent,
        role: "cap.end".into(),
    };
    let face = snapshot
        .elements
        .iter()
        .find(|element| element.key.kind == TopologyKind::Face && element.origins.contains(&cap))
        .expect("Named generated end cap");
    let created = sketch::prepare_new(
        &fixture.session,
        fixture.session.revision(),
        actor(),
        "Attached cap sketch",
        CreationPlane::Face {
            snapshot,
            key: face.key,
            anchor: Some(&cap),
        },
        PlaneOptions::default(),
    )
    .unwrap();
    let child_sketch = created.sketch;
    let child_plane = created.plane;
    fixture.session.commit(created.prepared).unwrap();
    let EntityData::Sketch { definition, .. } =
        &fixture.session.document().get(child_sketch).unwrap().data
    else {
        panic!()
    };
    let mut source = Source::from_definition(definition).unwrap();
    source.add_rectangle([5., 5.], [15., 15.]).unwrap();
    let region = sketch_regions::extract(&source).regions[0].intent();
    let edits = source
        .replace(fixture.session.document(), child_sketch)
        .unwrap();
    fixture.apply(edits);
    let child_request = model::ExtrudeRequest::new(
        child_sketch,
        vec![region],
        ExtrudeParameters {
            extent: Extent::OneSided {
                end: Termination::Distance {
                    distance: length(3.),
                },
            },
            ..Default::default()
        },
    );
    let prepared = model::prepare_extrude(
        &fixture.session,
        fixture.session.revision(),
        actor(),
        &child_request,
    )
    .unwrap();
    fixture.session.commit(prepared.prepared).unwrap();
    let before = fixture.session.document().clone();
    let child = child_request.outputs[0].body;
    let first = compiler
        .compile(&before, fixture.session.revision())
        .unwrap();
    assert_eq!(first.snapshot().bodies.len(), 3);
    let child_shape = shape(&first, child);
    near(child_shape.volume_mm3, 300.);
    point_near(child_shape.bounds.min, [5., 5., 10.]);
    point_near(child_shape.bounds.max, [15., 15., 13.]);
    let identities = first.snapshot().bodies.clone();

    let mut parent_parameters = fixture.request.parameters.clone();
    parent_parameters.extent = Extent::OneSided {
        end: Termination::Distance {
            distance: length(20.),
        },
    };
    let changed = model::prepare_extrude_edit(
        &fixture.session,
        fixture.session.revision(),
        actor(),
        parent,
        &parent_parameters,
    )
    .unwrap();
    fixture.session.commit(changed.prepared).unwrap();
    let canonical = fixture.session.document().clone();
    assert_eq!(canonical.get(child_sketch), before.get(child_sketch));
    assert_eq!(canonical.get(child_plane), before.get(child_plane));
    let second = compiler
        .compile(&canonical, fixture.session.revision())
        .unwrap();
    assert_eq!(second.snapshot().bodies, identities);
    near(shape(&second, child).volume_mm3, 300.);
    point_near(shape(&second, child).bounds.min, [5., 5., 20.]);
    point_near(shape(&second, child).bounds.max, [15., 15., 23.]);
    let frame = sketch::inspect_plane(
        &canonical,
        child_plane,
        fixture.session.revision(),
        Some(&second.snapshot().references),
    )
    .unwrap();
    point_near(frame.origin, [0., 0., 20.]);
    assert_eq!(fixture.session.document(), &canonical);
    assert!(
        sketch::prepare_new(
            &fixture.session,
            fixture.session.revision(),
            actor(),
            "Stale face",
            CreationPlane::Face {
                snapshot,
                key: face.key,
                anchor: Some(&cap)
            },
            PlaneOptions::default()
        )
        .is_err()
    );
}
#[test]
fn parameter_bound_extrusion_uses_source_units_and_stable_identity_after_rename_and_edit() {
    let mut fixture = SharedProfiles::new(Default::default());
    let feature = fixture.outer().feature;
    let depth = Entity::new(
        "depth",
        EntityData::Parameter {
            expression: "2".into(),
            unit: "cm".into(),
            exposed: true,
            inputs: vec![],
        },
    );
    fixture.apply(vec![Mutation::Insert(depth.clone())]);
    let expression =
        BoundExpression::bind(fixture.session.document(), "=$depth", "mm", None).unwrap();
    let parameters = ExtrudeParameters {
        extent: Extent::OneSided {
            end: Termination::Distance {
                distance: Dimension::Expression {
                    expression: expression.clone(),
                },
            },
        },
        ..Default::default()
    };
    let prepared = model::prepare_extrude_edit(
        &fixture.session,
        fixture.session.revision(),
        actor(),
        feature,
        &parameters,
    )
    .unwrap();
    fixture.session.commit(prepared.prepared).unwrap();
    let compiler = Compiler::default();
    let first = compiler
        .compile(fixture.session.document(), fixture.session.revision())
        .unwrap();
    near(
        shape(&first, fixture.outer().body).volume_mm3,
        (WIDTH * HEIGHT - INNER_AREA) * 20.,
    );
    near(
        shape(&first, fixture.inner().body).volume_mm3,
        INNER_AREA * 10.,
    );
    let bound_source = model::read_extrude(fixture.session.document(), feature).unwrap();
    let mut renamed = depth.clone();
    renamed.name = "renamed_depth".into();
    let EntityData::Parameter { expression, .. } = &mut renamed.data else {
        panic!()
    };
    *expression = "3".into();
    fixture.apply(vec![Mutation::Replace(renamed)]);
    let canonical = fixture.session.document().clone();
    let second = compiler
        .compile(&canonical, fixture.session.revision())
        .unwrap();
    near(
        shape(&second, fixture.outer().body).volume_mm3,
        (WIDTH * HEIGHT - INNER_AREA) * 30.,
    );
    point_near(
        shape(&second, fixture.outer().body).bounds.max,
        [WIDTH, HEIGHT, 30.],
    );
    assert_eq!(second.snapshot().bodies, first.snapshot().bodies);
    assert_eq!(
        model::read_extrude(&canonical, feature).unwrap(),
        bound_source
    );
    assert_eq!(fixture.session.document(), &canonical);
    let EntityData::Parameter { unit, .. } = &canonical.get(depth.id).unwrap().data else {
        panic!()
    };
    assert_eq!(unit, "cm");
    fixture.apply(vec![
        Mutation::Remove(depth.id),
        Mutation::Insert(Entity::new(
            "depth",
            EntityData::Parameter {
                expression: "90".into(),
                unit: "mm".into(),
                exposed: true,
                inputs: vec![],
            },
        )),
    ]);
    let error = compiler
        .compile(fixture.session.document(), fixture.session.revision())
        .err()
        .expect("Deleted binding must not rebind to a reused name");
    assert_eq!(error.entity, Some(depth.id));
}
#[test]
fn rotated_box_and_cylinder_primitives_compile_named_exact_solids_in_their_frames() {
    let mut session = Session::new(Document::new("Placed primitives").unwrap()).unwrap();
    let box_request = PrimitiveRequest::new(
        "Rotated box",
        PrimitiveParameters::Box {
            size: [length(10.), length(20.), length(30.)],
            frame: PrimitiveFrame {
                origin: [length(100.), length(-10.), length(5.)],
                u: [0., 1., 0.],
                v: [0., 0., 1.],
            },
        },
    );
    let cylinder = PrimitiveRequest::new(
        "Rotated cylinder",
        PrimitiveParameters::Cylinder {
            radius: length(3.),
            height: length(12.),
            frame: PrimitiveFrame {
                origin: [length(-20.), length(15.), length(8.)],
                u: [1., 0., 0.],
                v: [0., 0., -1.],
            },
        },
    );
    for request in [&box_request, &cylinder] {
        let prepared =
            model::prepare_primitive(&session, session.revision(), actor(), request).unwrap();
        session.commit(prepared.prepared).unwrap();
    }
    let original = session.document().clone();
    for backend in [GeometryBackend::OpenCascade, GeometryBackend::Rust] {
        let model = Compiler::default()
            .compile_with_backend(&original, session.revision(), backend, None)
            .unwrap();
        assert_eq!(model.snapshot().bodies.len(), 2);
        let box_shape = shape(&model, box_request.body);
        near(box_shape.volume_mm3, 6000.);
        point_near(box_shape.bounds.min, [100., -10., 5.]);
        point_near(box_shape.bounds.max, [130., 0., 25.]);
        let cylinder_shape = shape(&model, cylinder.body);
        near(cylinder_shape.volume_mm3, 108. * std::f64::consts::PI);
        point_near(cylinder_shape.bounds.min, [-23., 15., 5.]);
        point_near(cylinder_shape.bounds.max, [-17., 27., 11.]);
        for request in [&box_request, &cylinder] {
            assert_eq!(
                model.solid(request.feature).unwrap().uses_rust_geometry(),
                backend == GeometryBackend::Rust
            );
            let shape = shape(&model, request.body);
            assert_eq!(shape.solid_count, 1);
            for role in ["cap.start", "cap.end"] {
                let cap = shape
                    .topology
                    .elements
                    .iter()
                    .find(|element| {
                        element.key.kind == TopologyKind::Face
                            && element.origins.contains(&Origin {
                                source: request.feature,
                                role: role.into(),
                            })
                    })
                    .expect("Named primitive cap");
                assert!(
                    !cap.frames.is_empty(),
                    "Primitive planar cap has no truthful source frame"
                );
            }
        }
    }
    assert_eq!(session.document(), &original);
}
#[test]
fn sphere_primitive_compiles_edits_and_reopens_as_an_exact_named_solid() {
    let mut session = Session::new(Document::new("Sphere primitive").unwrap()).unwrap();
    let center = [100., -10., 5.];
    let frame = PrimitiveFrame {
        origin: center.map(length),
        u: [0., 1., 0.],
        v: [0., 0., 1.],
    };
    let sphere = |radius| PrimitiveParameters::Sphere {
        radius: length(radius),
        frame: frame.clone(),
    };
    let request = PrimitiveRequest::new("Sphere", sphere(10.));
    let initial = session.document().clone();
    let prepared = model::prepare_primitive(&session, 0, actor(), &request).unwrap();
    assert_eq!(session.document(), &initial);
    session.commit(prepared.prepared).unwrap();

    let compiler = Compiler::default();
    let check = |document: &Document, revision, radius: f64| {
        let model = compiler
            .compile_with_backend(document, revision, GeometryBackend::Rust, None)
            .unwrap();
        let sphere = shape(&model, request.body);
        assert_eq!(sphere.producer, request.feature);
        assert_eq!(sphere.solid_count, 1);
        near(
            sphere.volume_mm3,
            4. / 3. * std::f64::consts::PI * radius.powi(3),
        );
        point_near(sphere.bounds.min, center.map(|value| value - radius));
        point_near(sphere.bounds.max, center.map(|value| value + radius));
        assert!(!sphere.mesh.indices.is_empty());
        assert!(sphere.topology.elements.iter().any(|element| {
            element.key.kind == TopologyKind::Face
                && element.frames.is_empty()
                && element
                    .origins
                    .iter()
                    .any(|origin| origin.source == request.feature)
        }));
        sphere.topology.validate().unwrap();
        let exact = model.solid(request.feature).unwrap();
        let faces = exact.faces().unwrap();
        assert_eq!(faces.len(), 1);
        near(
            faces[0].area_mm2,
            4. * std::f64::consts::PI * radius.powi(2),
        );
        point_near(faces[0].center_mm, center);
        assert_eq!(faces[0].plane, None);
        let mesh = exact
            .manufacturing_mesh_with_receipt(Default::default())
            .unwrap();
        assert!(mesh.edge_pairs > 0 && mesh.vertex_links > 0);
        assert!(mesh.mesh.face_indices.iter().all(|&face| face == 0));
        assert_eq!(exact.display_edges(0.1).unwrap().len(), 1);
        let poles = exact.display_vertices().unwrap();
        assert_eq!(poles.len(), 2);
        point_near(
            poles[0].position_mm,
            [center[0] - radius, center[1], center[2]],
        );
        point_near(
            poles[1].position_mm,
            [center[0] + radius, center[1], center[2]],
        );
        let translated = exact
            .transformed(multicad_kernel::RigidTransform {
                rotation: [[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]],
                translation_mm: [20., 30., 40.],
            })
            .unwrap();
        let moved_center = [30., 130., 45.];
        let moved_poles = translated.display_vertices().unwrap();
        point_near(moved_poles[0].position_mm, [30., 130. - radius, 45.]);
        point_near(moved_poles[1].position_mm, [30., 130. + radius, 45.]);
        point_near(
            translated.bounds().unwrap().min,
            moved_center.map(|x| x - radius),
        );
        assert_eq!(
            exact
                .topology_snapshot(request.body, revision)
                .unwrap()
                .elements
                .into_iter()
                .map(|e| e.origins)
                .collect::<Vec<_>>(),
            translated
                .topology_snapshot(request.body, revision)
                .unwrap()
                .elements
                .into_iter()
                .map(|e| e.origins)
                .collect::<Vec<_>>()
        );
    };
    check(session.document(), session.revision(), 10.);
    assert!(
        model::prepare_primitive_edit(
            &session,
            session.revision(),
            actor(),
            request.feature,
            &sphere(0.),
        )
        .is_err()
    );
    assert_eq!(session.revision(), 1);
    let edit = model::prepare_primitive_edit(
        &session,
        session.revision(),
        actor(),
        request.feature,
        &sphere(4.),
    )
    .unwrap();
    session.commit(edit.prepared).unwrap();
    check(session.document(), session.revision(), 4.);
    session.undo(session.revision(), actor()).unwrap();
    check(session.document(), session.revision(), 10.);
    session.redo(session.revision(), actor()).unwrap();
    check(session.document(), session.revision(), 4.);

    let project = ProjectFile::new(session.document().clone()).unwrap();
    assert!(project.capabilities().contains(model::SPHERE_CAPABILITY));
    let mut bytes = vec![];
    project_file::write(&project, &mut bytes, Default::default()).unwrap();
    let supported = project_file::supported_capabilities();
    let opened = project_file::read(&mut bytes.as_slice(), Default::default(), &supported)
        .unwrap()
        .project;
    assert_eq!(opened.document, project.document);
    assert_eq!(
        model::read_primitive(&opened.document, request.feature)
            .unwrap()
            .parameters,
        sphere(4.)
    );
    check(&opened.document, 0, 4.);
    let mut legacy = supported;
    legacy.remove(model::SPHERE_CAPABILITY);
    assert!(project_file::read(&mut bytes.as_slice(), Default::default(), &legacy).is_err());
}

fn revolved_fixture(degrees: f64) -> (SharedProfiles, RevolveRequest) {
    let mut fixture = SharedProfiles::new(Default::default());
    let regions = fixture
        .request
        .outputs
        .iter()
        .map(|output| output.region.clone())
        .collect();
    let removals = fixture
        .request
        .outputs
        .iter()
        .flat_map(|output| {
            [
                Mutation::Remove(output.feature),
                Mutation::Remove(output.body),
            ]
        })
        .collect();
    fixture.apply(removals);
    let request = RevolveRequest::new(
        fixture.sketch,
        regions,
        RevolveParameters {
            axis_origin: [length(-10.), length(0.)],
            axis_angle: angle(90.),
            angle: angle(degrees),
        },
    );
    let prepared = model::prepare_revolve(
        &fixture.session,
        fixture.session.revision(),
        actor(),
        &request,
    )
    .unwrap();
    fixture.session.commit(prepared.prepared).unwrap();
    (fixture, request)
}
#[test]
fn shared_profiles_revolve_into_signed_partial_and_full_real_solids() {
    for degrees in [90_f64, -90., 360., -360.] {
        let (fixture, request) = revolved_fixture(degrees);
        let original = fixture.session.document().clone();
        let model = Compiler::default()
            .compile(&original, fixture.session.revision())
            .unwrap();
        assert_eq!(model.snapshot().bodies.len(), 2);
        let outer = request
            .outputs
            .iter()
            .find(|output| output.region.holes.len() == 1)
            .unwrap();
        let inner = request
            .outputs
            .iter()
            .find(|output| output.region.holes.is_empty())
            .unwrap();
        let outer_shape = shape(&model, outer.body);
        let inner_shape = shape(&model, inner.body);
        let sweep = degrees.to_radians().abs();
        near(
            outer_shape.volume_mm3,
            (WIDTH * HEIGHT * 50. - INNER_AREA * 30.) * sweep,
        );
        near(inner_shape.volume_mm3, INNER_AREA * 30. * sweep);
        near(
            outer_shape.volume_mm3 + inner_shape.volume_mm3,
            240000. * sweep,
        );
        if degrees.abs() == 360. {
            point_near(outer_shape.bounds.min, [-100., 0., -90.]);
            point_near(outer_shape.bounds.max, [80., HEIGHT, 90.]);
            point_near(inner_shape.bounds.min, [-46., 14., -36.]);
            point_near(inner_shape.bounds.max, [26., 26., 36.]);
        } else if degrees > 0. {
            point_near(outer_shape.bounds.min, [-10., 0., -90.]);
            point_near(outer_shape.bounds.max, [80., HEIGHT, 0.]);
            point_near(inner_shape.bounds.min, [-10., 14., -36.]);
            point_near(inner_shape.bounds.max, [26., 26., 0.]);
        } else {
            point_near(outer_shape.bounds.min, [-10., 0., 0.]);
            point_near(outer_shape.bounds.max, [80., HEIGHT, 90.]);
            point_near(inner_shape.bounds.min, [-10., 14., 0.]);
            point_near(inner_shape.bounds.max, [26., 26., 36.]);
        }
        for output in &request.outputs {
            let result = shape(&model, output.body);
            assert_eq!(result.solid_count, 1);
            let caps = result
                .topology
                .elements
                .iter()
                .filter(|element| {
                    element.key.kind == TopologyKind::Face
                        && element.origins.iter().any(|origin| {
                            origin.source == output.feature
                                && ["cap.start", "cap.end"].contains(&origin.role.as_str())
                        })
                })
                .count();
            assert_eq!(caps, if degrees.abs() == 360. { 0 } else { 2 });
        }
        assert_eq!(fixture.session.document(), &original);
        assert_eq!(
            model.snapshot().stamp.source_hash,
            ContentHash::document(&original).unwrap()
        );
    }
}
#[test]
fn revolve_angle_edits_preserve_identity_and_crossing_axis_is_an_explicit_error() {
    let (mut fixture, request) = revolved_fixture(90.);
    let compiler = Compiler::default();
    let first = compiler
        .compile(fixture.session.document(), fixture.session.revision())
        .unwrap();
    let original_sketch = fixture
        .session
        .document()
        .get(fixture.sketch)
        .unwrap()
        .clone();
    let mut changed = request.parameters.clone();
    changed.angle = angle(-90.);
    for output in &request.outputs {
        let prepared = model::prepare_revolve_edit(
            &fixture.session,
            fixture.session.revision(),
            actor(),
            output.feature,
            &changed,
        )
        .unwrap();
        fixture.session.commit(prepared.prepared).unwrap();
    }
    let second = compiler
        .compile(fixture.session.document(), fixture.session.revision())
        .unwrap();
    assert_eq!(second.snapshot().bodies, first.snapshot().bodies);
    assert_eq!(
        fixture.session.document().get(fixture.sketch).unwrap(),
        &original_sketch
    );
    let outer = request
        .outputs
        .iter()
        .find(|output| output.region.holes.len() == 1)
        .unwrap();
    near(
        shape(&second, outer.body).volume_mm3,
        shape(&first, outer.body).volume_mm3,
    );
    point_near(shape(&second, outer.body).bounds.max, [80., HEIGHT, 90.]);
    let mut crossing = changed;
    crossing.axis_origin = [length(40.), length(0.)];
    let prepared = model::prepare_revolve_edit(
        &fixture.session,
        fixture.session.revision(),
        actor(),
        outer.feature,
        &crossing,
    )
    .unwrap();
    fixture.session.commit(prepared.prepared).unwrap();
    let original = fixture.session.document().clone();
    let error = compiler
        .compile(&original, fixture.session.revision())
        .err()
        .expect("A section crossing the axis must not produce a self-overlapping solid");
    assert_eq!(error.entity, Some(outer.feature));
    assert!(error.message.contains("crosses"));
    assert_eq!(fixture.session.document(), &original);
    let mut source = model::read_revolve(&original, outer.feature).unwrap();
    source.parameters.angle = angle(0.);
    let mut feature = original.get(outer.feature).unwrap().clone();
    let EntityData::Feature { definition, .. } = &mut feature.data else {
        panic!()
    };
    *definition = source.definition().unwrap();
    fixture.apply(vec![Mutation::Replace(feature)]);
    let error = compiler
        .compile(fixture.session.document(), fixture.session.revision())
        .err()
        .expect("Zero angle cannot default to a full turn");
    assert_eq!(error.entity, Some(outer.feature));
    assert!(error.message.contains("nonzero"));
}

/// Closed, consistently oriented manufacturing triangles whose enclosed volume
/// matches the exact solid within the requested chordal accuracy.
fn assert_closed_accurate_mesh(solid: &Solid, exact: f64) {
    let deflection = 0.01;
    let receipt = solid
        .manufacturing_mesh_with_receipt(multicad_kernel::MeshOptions {
            linear_deflection_mm: deflection,
            maximum_triangles: 4_000_000,
            ..Default::default()
        })
        .unwrap();
    let mesh = &receipt.mesh;
    let point = |i: u32| -> [f64; 3] { std::array::from_fn(|k| mesh.vertices[3 * i as usize + k]) };
    let mut directed = std::collections::BTreeMap::<(u32, u32), i32>::new();
    let mut volume = 0.;
    for t in mesh.indices.chunks(3) {
        for k in 0..3 {
            let (a, b) = (t[k], t[(k + 1) % 3]);
            *directed.entry((a.min(b), a.max(b))).or_default() += if a < b { 1 } else { -1 };
        }
        let [a, b, c] = [point(t[0]), point(t[1]), point(t[2])];
        volume += (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0]))
            / 6.;
    }
    assert!(directed.values().all(|&uses| uses == 0));
    assert!(receipt.maximum_displacement_mm <= deflection);
    let area: f64 = solid.faces().unwrap().iter().map(|f| f.area_mm2).sum();
    assert!(
        (volume - exact).abs() <= deflection * area,
        "{volume} != {exact}"
    );
}
#[test]
fn curved_revolution_keeps_identity_through_sketch_and_angle_edits_undo_redo_and_cold_reopen() {
    let (mut fixture, request) = revolved_fixture(90.);
    let outer = request
        .outputs
        .iter()
        .find(|output| output.region.holes.len() == 1)
        .unwrap();
    let (feature, body) = (outer.feature, outer.body);
    // Pappus: the plate's first moment about x = -10 minus the disc's.
    let expected = |width: f64, degrees: f64| {
        (width * HEIGHT * (width / 2. + 10.) - INNER_AREA * 30.) * degrees.to_radians().abs()
    };
    let check = |document: &Document, revision: u64, width: f64| {
        let degrees = 90.;
        let model = Compiler::default().compile(document, revision).unwrap();
        let result = shape(&model, body);
        assert_eq!(result.producer, feature);
        assert_eq!(result.solid_count, 1);
        near(result.volume_mm3, expected(width, degrees));
        // A positive quarter turn sweeps toward -z from the sketch plane.
        point_near(result.bounds.max, [width, HEIGHT, 0.]);
        result.topology.validate().unwrap();
        assert_closed_accurate_mesh(model.solid(feature).unwrap(), expected(width, degrees));
        result.topology.elements.clone()
    };
    let first = check(
        fixture.session.document(),
        fixture.session.revision(),
        WIDTH,
    );

    // A sketch edit moves the outer corners; the feature and Part persist.
    fixture.set_width(70.);
    let narrowed = check(fixture.session.document(), fixture.session.revision(), 70.);
    let faces = |elements: &[multicad_core::references::TopologyElement]| {
        elements
            .iter()
            .filter(|e| e.key.kind == TopologyKind::Face)
            .map(|e| (e.key, e.origins.clone()))
            .collect::<Vec<_>>()
    };
    assert_eq!(faces(&narrowed), faces(&first));
    fixture
        .session
        .undo(fixture.session.revision(), actor())
        .unwrap();
    assert_eq!(
        check(
            fixture.session.document(),
            fixture.session.revision(),
            WIDTH
        ),
        first
    );
    fixture
        .session
        .redo(fixture.session.revision(), actor())
        .unwrap();
    assert_eq!(
        check(fixture.session.document(), fixture.session.revision(), 70.),
        narrowed
    );

    // A full reverse turn closes the hole into a void and drops both caps.
    let mut full = request.parameters.clone();
    full.angle = angle(-360.);
    let prepared = model::prepare_revolve_edit(
        &fixture.session,
        fixture.session.revision(),
        actor(),
        feature,
        &full,
    )
    .unwrap();
    fixture.session.commit(prepared.prepared).unwrap();
    let model = Compiler::default()
        .compile(fixture.session.document(), fixture.session.revision())
        .unwrap();
    let result = shape(&model, body);
    near(result.volume_mm3, expected(70., 360.));
    point_near(result.bounds.min, [-90., 0., -80.]);
    assert!(!result.topology.elements.iter().any(|e| {
        e.origins
            .iter()
            .any(|o| o.source == feature && (o.role == "cap.start" || o.role == "cap.end"))
    }));
    assert_closed_accurate_mesh(model.solid(feature).unwrap(), expected(70., 360.));
    let turned = result.topology.elements.clone();

    let project = ProjectFile::new(fixture.session.document().clone()).unwrap();
    let mut bytes = vec![];
    project_file::write(&project, &mut bytes, Default::default()).unwrap();
    let supported = project_file::supported_capabilities();
    let opened = project_file::read(&mut bytes.as_slice(), Default::default(), &supported)
        .unwrap()
        .project;
    assert_eq!(opened.document, project.document);
    assert_eq!(
        model::read_revolve(&opened.document, feature)
            .unwrap()
            .parameters,
        full
    );
    let reopened = Compiler::default().compile(&opened.document, 0).unwrap();
    let again = shape(&reopened, body);
    near(again.volume_mm3, expected(70., 360.));
    assert_eq!(again.topology.elements, turned);
}
