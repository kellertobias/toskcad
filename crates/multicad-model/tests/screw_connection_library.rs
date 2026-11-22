//! Real pinned mechanical Library roles remain associative connection sources.
use multicad_precision::PortableMath as _;
#[path = "support/custom_screw_connections.rs"]
mod custom_fixture;
#[path = "../../multicad-mcp/tests/support/screw_connection.rs"]
mod fixture;
use multicad_core::{
    commands::{Actor, Channel, Request, Session},
    constraint_sketch::Dimension,
    document::{Document, Entity, EntityData, EntityId, Mutation, Placement},
    model::*,
    parameters::BoundExpression,
    project_file::{self, ProjectFile},
    symbols::{
        self, SymbolAuthorRequest, SymbolDefinition, SymbolEditRequest, SymbolKind,
        SymbolPlaceRequest, SymbolValue,
    },
};
use multicad_model::{Compiler, EvaluationCache, Model};
use std::collections::BTreeMap;
fn actor() -> Actor {
    Actor::local(Channel::Mcp)
}

#[test]
fn actual_pinned_hanger_nut_and_dynamic_ring_keep_sources_through_edit_and_reopen() {
    for (left, posed) in [(false, false), (true, false), (false, true), (true, true)] {
        let (doc, mut request) = custom_fixture::fixture(ScrewConnectionMode::HangerBolt);
        let mut ring_parameter = None;
        for row in &mut request.hardware {
            let ScrewConnectionHardwareSource::Generated { part } = &row.source else {
                panic!()
            };
            let mut parameters = part.parameters.clone();
            if left {
                match &mut parameters {
                    MechanicalPartParameters::CustomConnection { parameters } => {
                        let CustomConnectionHardwareParameters::HangerBolt {
                            wood_thread,
                            metric_thread,
                            ..
                        } = parameters.as_mut()
                        else {
                            panic!()
                        };
                        wood_thread.hand = ThreadHand::Left;
                        metric_thread.hand = ThreadHand::Left;
                    }
                    MechanicalPartParameters::HexNut { thread, .. } => {
                        thread.hand = ThreadHand::Left
                    }
                    _ => {}
                }
            }
            let dynamic = row.role == ScrewConnectionHardwareRole::NutWasher;
            let (_, definition, parameter) = author(parameters, dynamic);
            let pinned = library_row(&definition, row.role, "Pinned actual hardware");
            if dynamic {
                ring_parameter = Some((pinned.instance, parameter.unwrap()));
            }
            *row = pinned;
        }
        let mut session = Session::new(doc).unwrap();
        if posed {
            let u = [0.8, 0.6, 0.];
            let v = [-0.36, 0.48, 0.8];
            let n = [
                u[1] * v[2] - u[2] * v[1],
                u[2] * v[0] - u[0] * v[2],
                u[0] * v[1] - u[1] * v[0],
            ];
            for (i, row) in request.receivers.iter().enumerate() {
                let feature = row.source.input.producer;
                let mut primitive = read_primitive(session.document(), feature).unwrap();
                let PrimitiveParameters::Box { frame, size } = &mut primitive.parameters else {
                    panic!()
                };
                let origin = if i == 1 {
                    size[2] = length(32.);
                    n.map(|x| -32. * x)
                } else {
                    [0.; 3]
                };
                *frame = PrimitiveFrame {
                    origin: origin.map(length),
                    u,
                    v,
                };
                let p = prepare_primitive_edit(
                    &session,
                    session.revision(),
                    actor(),
                    feature,
                    &primitive.parameters,
                )
                .unwrap();
                session.commit(p.prepared).unwrap();
            }
        }
        let prepared =
            prepare_screw_connection(&session, session.revision(), actor(), &request).unwrap();
        let mut cache = EvaluationCache::default();
        let native = compiler()
            .compile_cached(
                &mut cache,
                prepared.prepared.document(),
                session.revision() + 1,
                None,
            )
            .unwrap();
        assert!(native.snapshot().screw_connection_checks[0].native_component_collision_checked);
        for row in &request.hardware {
            let current = leaf(&native, row.instance);
            assert!(
                !current.path.members.is_empty(),
                "pin must retain mapped member ownership"
            );
            let solid = native.leaf_solid(&current.path).unwrap();
            assert!(
                solid
                    .manufacturing_mesh_with_receipt(Default::default())
                    .unwrap()
                    .edge_pairs
                    > 0
            );
        }
        drop(native);
        session.commit(prepared.prepared).unwrap();
        let before = session.document().clone();
        let (record, parameter) = ring_parameter.unwrap();
        let edit = SymbolEditRequest {
            record,
            overrides: BTreeMap::from([(
                parameter,
                SymbolValue {
                    value: 0.27,
                    unit: "mm".into(),
                },
            )]),
        };
        let p = symbols::prepare_edit(&session, session.revision(), actor(), &edit).unwrap();
        let native = compiler()
            .compile_cached(
                &mut cache,
                p.prepared.document(),
                session.revision() + 1,
                None,
            )
            .unwrap();
        assert!(native.snapshot().screw_connection_checks[0].native_component_collision_checked);
        drop(native);
        session.commit(p.prepared).unwrap();
        let after = session.document().clone();
        let opened = reopened(&after);
        assert!(
            compiler()
                .compile(&opened, 0)
                .unwrap()
                .snapshot()
                .screw_connection_checks[0]
                .native_component_collision_checked
        );
        session.undo(session.revision(), actor()).unwrap();
        assert_eq!(*session.document(), before);
        session.redo(session.revision(), actor()).unwrap();
        assert_eq!(*session.document(), after);
        let invalid = SymbolEditRequest {
            record,
            overrides: BTreeMap::from([(
                parameter,
                SymbolValue {
                    value: 100.,
                    unit: "mm".into(),
                },
            )]),
        };
        if let Ok(p) = symbols::prepare_edit(&session, session.revision(), actor(), &invalid) {
            assert!(
                compiler()
                    .compile(p.prepared.document(), session.revision() + 1)
                    .is_err()
            );
        }
        assert_eq!(*session.document(), after);
        let mut stale = request.clone();
        let ScrewConnectionHardwareSource::Library { definition, .. } =
            &mut stale.hardware[0].source
        else {
            panic!()
        };
        definition.name.push_str(" stale pin");
        assert!(prepare_screw_connection(&session, session.revision(), actor(), &stale).is_err());
        assert_eq!(*session.document(), after);
    }
}
fn compiler() -> Compiler {
    Compiler {
        omit_display: true,
        ..Default::default()
    }
}
fn author(
    mut parameters: MechanicalPartParameters,
    dynamic: bool,
) -> (Document, SymbolDefinition, Option<EntityId>) {
    let mut s = Session::new(Document::new("Authored connection Library source").unwrap()).unwrap();
    let parameter = if dynamic {
        let e = Entity::new(
            "washer_thickness",
            EntityData::Parameter {
                expression: "0.2".into(),
                unit: "mm".into(),
                exposed: true,
                inputs: vec![],
            },
        );
        let id = e.id;
        s.apply(Request {
            expected_revision: 0,
            actor: actor(),
            label: "Author washer input".into(),
            edits: vec![Mutation::Insert(e)],
        })
        .unwrap();
        let MechanicalPartParameters::Washer { thickness, .. } = &mut parameters else {
            unreachable!()
        };
        *thickness = Dimension::Expression {
            expression: BoundExpression::bind(s.document(), "=$washer_thickness", "mm", None)
                .unwrap(),
        };
        Some(id)
    } else {
        None
    };
    let part = MechanicalPartRequest::new("Authored mechanical provider", parameters);
    let q = prepare_mechanical_part(&s, s.revision(), actor(), &part).unwrap();
    s.commit(q.prepared).unwrap();
    let request = SymbolAuthorRequest::new(
        "Test-authored exact mechanical pin",
        part.body,
        if dynamic {
            SymbolKind::Dynamic
        } else {
            SymbolKind::Static
        },
        parameter.into_iter().collect(),
        None,
    );
    let q = symbols::prepare_author(&s, s.revision(), actor(), &request).unwrap();
    s.commit(q.prepared).unwrap();
    let definition = symbols::read_definition(s.document(), request.definition).unwrap();
    (s.document().clone(), definition, parameter)
}
fn library_row(
    definition: &SymbolDefinition,
    role: ScrewConnectionHardwareRole,
    name: &str,
) -> ScrewConnectionHardwareRequest {
    let placement = SymbolPlaceRequest::from_definition(
        EntityId::new(),
        definition,
        name,
        BTreeMap::new(),
        Placement::default(),
    )
    .unwrap();
    ScrewConnectionHardwareRequest {
        role_id: EntityId::new(),
        role,
        name: name.into(),
        instance: placement.instance,
        driver: EntityId::new(),
        source: ScrewConnectionHardwareSource::Library {
            definition: Box::new(definition.clone()),
            placement: Box::new(placement),
        },
    }
}
fn placement(row: &ScrewConnectionHardwareRequest) -> &SymbolPlaceRequest {
    let ScrewConnectionHardwareSource::Library { placement, .. } = &row.source else {
        unreachable!()
    };
    placement
}
fn leaf(native: &Model, instance: EntityId) -> &multicad_model::LeafOccurrenceSnapshot {
    let leaves = native
        .snapshot()
        .leaf_occurrences
        .values()
        .filter(|leaf| leaf.path.root == instance)
        .collect::<Vec<_>>();
    assert_eq!(
        leaves.len(),
        1,
        "A Library role is one actual physical mechanical leaf"
    );
    leaves[0]
}
fn assert_native(native: &Model, request: &ScrewConnectionRequest, washer: f64) -> RigidPose {
    assert_eq!(native.snapshot().leaf_occurrences.len(), 4);
    let [check] = native.snapshot().screw_connection_checks.as_slice() else {
        panic!("One full native connection receipt")
    };
    assert_eq!(check.connection, request.connection);
    assert!(check.native_component_collision_checked);
    assert!((check.actual_engagement_mm.unwrap() - 2.4).abs() < 1e-7);
    assert!((check.grip_mm - 4.1).abs() < 1e-7);
    let screw = leaf(native, request.hardware[0].instance);
    let ring = leaf(native, request.hardware[1].instance);
    let nut = leaf(native, request.hardware[2].instance);
    assert!((screw.pose.translation_mm[2] - (4.1 + washer)).abs() < 1e-7);
    assert!((ring.pose.translation_mm[2] - (4.1 + washer)).abs() < 1e-7);
    assert!(nut.pose.translation_mm[2].abs() < 1e-7);
    // The former arbitrary CSG cube/solid intersection witness is unavailable
    // for Rust helices. Check the entire retained mating profiles and actual
    // head/shank separation, then independently manufacture and reopen STEP.
    let screw_solid = native.leaf_solid(&screw.path).unwrap();
    let nut_solid = native.leaf_solid(&nut.path).unwrap();
    let [(external, a), (internal, b)] = screw_solid.rust_screw_nut_threads(nut_solid).unwrap();
    let dot = |a: [f64; 3], b: [f64; 3]| a.into_iter().zip(b).map(|(a, b)| a * b).sum::<f64>();
    let alpha = dot(b.u, a.v).portable_atan2(dot(b.u, a.u));
    let station = dot(
        std::array::from_fn(|i| b.origin[i] - a.origin[i]),
        a.normal(),
    );
    assert!(
        external
            .minimum_radial_clearance(&internal, external.helical_rate() * alpha - station, false)
            .unwrap()
            > 1e-6
    );
    let (actual, frame) = screw_solid.rust_machine_screw_geometry().unwrap();
    let head_station = dot(
        std::array::from_fn(|i| {
            frame.origin[i] + frame.normal()[i] * actual.head_seat_z() - b.origin[i]
        }),
        b.normal(),
    );
    assert!(head_station + actual.spec().shank_length_mm < 0.);
    let directory = tempfile::tempdir().unwrap();
    for (index, leaf) in [screw, ring, nut].into_iter().enumerate() {
        let solid = native.leaf_solid(&leaf.path).unwrap();
        let receipt = solid
            .manufacturing_mesh_with_receipt(Default::default())
            .unwrap();
        assert!(receipt.edge_pairs > 0 && receipt.vertex_links > 0);
        let path = directory.path().join(format!("hardware-{index}.step"));
        solid.write_step(&path).unwrap();
        let back = multicad_kernel::Solid::read_step(&path, Default::default()).unwrap();
        assert_eq!(back.solid_count().unwrap(), 1);
        assert!((back.volume_mm3().unwrap() - solid.volume_mm3().unwrap()).abs() < 1e-5);
    }

    nut.pose
}
fn reopened(document: &Document) -> Document {
    let mut bytes = vec![];
    project_file::write(
        &ProjectFile::new(document.clone()).unwrap(),
        &mut bytes,
        Default::default(),
    )
    .unwrap();
    project_file::read(
        &mut bytes.as_slice(),
        Default::default(),
        &project_file::supported_capabilities(),
    )
    .unwrap()
    .project
    .document
}

#[test]
fn replacing_saved_hardware_source_preserves_the_single_physical_occurrence() {
    let (document, mut request, _) = fixture::fixture();
    let parameters = MechanicalPartParameters::Washer {
        bore_diameter: length(3.2),
        outside_diameter: length(7.),
        thickness: length(0.5),
        frame: Default::default(),
    };
    let part = MechanicalPartRequest::new("Generated washer", parameters.clone());
    let generated = ScrewConnectionHardwareRequest {
        role_id: EntityId::new(),
        role: ScrewConnectionHardwareRole::HeadWasher,
        name: part.name.clone(),
        instance: EntityId::new(),
        driver: EntityId::new(),
        source: ScrewConnectionHardwareSource::Generated {
            part: Box::new(part.clone()),
        },
    };
    request.hardware.insert(1, generated.clone());
    let mut session = Session::new(document.clone()).unwrap();
    let initial = prepare_screw_connection(&session, 0, actor(), &request).unwrap();
    session.commit(initial.prepared).unwrap();
    let before = reopened(session.document());
    let mut session = Session::new(before.clone()).unwrap();
    let source = project_file::read(
        &mut include_bytes!("../../../library/washer-normal-m3.tcd").as_slice(),
        Default::default(),
        &project_file::supported_capabilities(),
    )
    .unwrap();
    let mut definitions = symbols::definitions(&source.project.document).unwrap();
    let (_, shipped) = definitions.pop().unwrap();
    assert!(definitions.is_empty());
    let (_, alternate, _) = author(parameters, false);
    let mut retired = vec![part.body, part.feature];
    for definition in [shipped, alternate] {
        let mut incoming = SymbolPlaceRequest::from_definition(
            EntityId::new(),
            &definition,
            &definition.name,
            BTreeMap::new(),
            Placement::default(),
        )
        .unwrap();
        incoming.instance = generated.instance;
        request.hardware[1].source = ScrewConnectionHardwareSource::Library {
            definition: Box::new(definition.clone()),
            placement: Box::new(incoming.clone()),
        };
        let original = session.document().clone();
        let edit =
            prepare_screw_connection(&session, session.revision(), actor(), &request).unwrap();
        assert_eq!(session.document(), &original);
        let candidate = edit.prepared.document();
        assert!(retired.iter().all(|id| candidate.get(*id).is_none()));
        let saved = read_screw_connection_request(candidate, request.connection).unwrap();
        assert_eq!(saved.hardware[1].instance, generated.instance);
        assert_eq!(saved.hardware[1].role_id, generated.role_id);
        assert_eq!(saved.hardware[1].driver, generated.driver);
        assert_eq!(saved, request);
        let pin = symbols::read_instance(candidate, incoming.record).unwrap();
        assert_eq!(pin.pin, definition.pin().unwrap());
        retired.extend([incoming.record, incoming.wrapper]);
        session.commit(edit.prepared).unwrap();
        assert_eq!(reopened(session.document()), *session.document());
        session.undo(session.revision(), actor()).unwrap();
        assert_eq!(session.document(), &original);
        session.redo(session.revision(), actor()).unwrap();
    }
    let native = compiler()
        .compile(session.document(), session.revision())
        .unwrap();
    assert_native(&native, &request, 0.5);

    request.hardware[1] = generated;
    let edit = prepare_screw_connection(&session, session.revision(), actor(), &request).unwrap();
    let candidate = edit.prepared.document();
    assert!(retired[2..].iter().all(|id| candidate.get(*id).is_none()));
    assert_eq!(
        read_screw_connection_request(candidate, request.connection).unwrap(),
        request
    );
}

#[test]
fn pinned_metric_nut_and_dynamic_library_washer_follow_overrides_with_native_proof_and_reopen() {
    let (document, mut request, _) = fixture::fixture();
    let ScrewConnectionHardwareSource::Generated { part } = &request.hardware[1].source else {
        unreachable!()
    };
    let nut_parameters = part.parameters.clone();
    let (nut_author, nut_definition, _) = author(nut_parameters.clone(), false);
    let (washer_author, washer_definition, thickness) = author(
        MechanicalPartParameters::Washer {
            bore_diameter: length(3.2),
            outside_diameter: length(7.),
            thickness: length(0.2),
            frame: Default::default(),
        },
        true,
    );
    let thickness = thickness.unwrap();
    let nut_author_bytes = reopened(&nut_author);
    let washer_author_bytes = reopened(&washer_author);
    request.hardware[1] = library_row(
        &nut_definition,
        ScrewConnectionHardwareRole::Nut,
        "Pinned 6H nut",
    );
    request.hardware.insert(
        1,
        library_row(
            &washer_definition,
            ScrewConnectionHardwareRole::HeadWasher,
            "Pinned parametric washer",
        ),
    );
    let ring = placement(&request.hardware[1]).clone();
    let nut = placement(&request.hardware[2]).clone();
    let mut s = Session::new(document).unwrap();
    let before = s.document().clone();
    let q = prepare_screw_connection(&s, 0, actor(), &request).unwrap();
    let native = compiler().compile(q.prepared.document(), 1).unwrap();
    let initial_pose = assert_native(&native, &request, 0.2);
    let initial_hash = native.snapshot().stamp.source_hash;
    drop(native);
    assert_eq!(s.document(), &before, "Preview is source-neutral");
    s.commit(q.prepared).unwrap();
    let source = s.document().get(request.connection).unwrap().clone();
    let initial = s.document().clone();
    let ring_before = symbols::read_instance(&initial, ring.instance).unwrap();
    let nut_before = symbols::read_instance(&initial, nut.instance).unwrap();
    assert_eq!(ring_before.pin, washer_definition.pin().unwrap());
    assert_eq!(nut_before.pin, nut_definition.pin().unwrap());
    assert_eq!(ring_before.ids, ring.ids);
    assert_eq!(nut_before.ids, nut.ids);
    assert_ne!(ring.ids[&thickness], thickness);
    assert!(
        symbols::owned_entities(&initial)
            .unwrap()
            .contains(&ring.ids[&thickness])
    );
    let reread = read_screw_connection_request(&initial, request.connection).unwrap();
    let ResolvedMechanicalPart::HexNut { thread, .. } =
        resolve_screw_connection_hardware(&initial, s.parameters(), &reread.hardware[2]).unwrap()
    else {
        unreachable!()
    };
    assert!(thread.metric_fit.is_some());
    let EntityData::Part { features, .. } =
        &initial.get(nut.ids[&nut_definition.root]).unwrap().data
    else {
        unreachable!()
    };
    let MechanicalPartParameters::HexNut {
        thread: authored, ..
    } = read_mechanical_part(&initial, features[0])
        .unwrap()
        .parameters
    else {
        unreachable!()
    };
    assert_eq!(authored.fit, ThreadFit::Metric6HInternal);
    assert_eq!(thread.representation, ThreadRepresentation::Modeled);
    let on_disk = reopened(&initial);
    assert_eq!(on_disk, initial);
    let mut cache = EvaluationCache::default();
    assert_native(
        &compiler()
            .compile_cached(&mut cache, &on_disk, s.revision(), None)
            .unwrap(),
        &request,
        0.2,
    );
    let edit = SymbolEditRequest {
        record: ring.instance,
        overrides: BTreeMap::from([(
            thickness,
            SymbolValue {
                value: 0.3,
                unit: "mm".into(),
            },
        )]),
    };
    let q = symbols::prepare_edit(&s, s.revision(), actor(), &edit).unwrap();
    let after_native = compiler()
        .compile_cached(&mut cache, q.prepared.document(), s.revision() + 1, None)
        .unwrap();
    assert_ne!(
        assert_native(&after_native, &request, 0.3),
        initial_pose,
        "Current washer thickness changes the actual mating phase"
    );
    assert_ne!(after_native.snapshot().stamp.source_hash, initial_hash);
    drop(after_native);
    assert_eq!(s.document(), &initial);
    s.commit(q.prepared).unwrap();
    assert_eq!(
        s.document().get(request.connection),
        Some(&source),
        "No connection recomposition"
    );
    let ring_after = symbols::read_instance(s.document(), ring.instance).unwrap();
    assert_eq!(ring_after.pin, ring_before.pin);
    assert_eq!(ring_after.ids, ring_before.ids);
    assert_eq!(
        symbols::read_instance(s.document(), nut.instance).unwrap(),
        nut_before
    );
    let changed = initial
        .entities()
        .filter(|e| s.document().get(e.id) != Some(*e))
        .map(|e| e.id)
        .collect::<Vec<_>>();
    assert_eq!(
        changed.len(),
        2,
        "Only override record and mapped scalar may change: {changed:?}"
    );
    assert!(changed.contains(&ring.record));
    assert!(changed.contains(&ring.ids[&thickness]));
    assert_eq!(initial.entities().count(), s.document().entities().count());
    let final_document = reopened(s.document());
    assert_eq!(&final_document, s.document());
    assert_native(
        &compiler()
            .compile_cached(&mut cache, &final_document, s.revision(), None)
            .unwrap(),
        &request,
        0.3,
    );
    // A private mapped ID is not a public exposed source key.
    let bad_private = SymbolEditRequest {
        record: ring.instance,
        overrides: BTreeMap::from([(
            ring.ids[&thickness],
            SymbolValue {
                value: 0.4,
                unit: "mm".into(),
            },
        )]),
    };
    assert!(symbols::prepare_edit(&s, s.revision(), actor(), &bad_private).is_err());
    // A well-typed override can still make native geometry invalid; never commit that candidate.
    let bad_geometry = SymbolEditRequest {
        record: ring.instance,
        overrides: BTreeMap::from([(
            thickness,
            SymbolValue {
                value: -0.1,
                unit: "mm".into(),
            },
        )]),
    };
    let invalid = symbols::prepare_edit(&s, s.revision(), actor(), &bad_geometry).unwrap();
    assert!(
        compiler()
            .compile(invalid.prepared.document(), s.revision() + 1)
            .is_err()
    );
    assert_eq!(s.document(), &final_document);
    let mut incompatible_parameters = nut_parameters;
    let MechanicalPartParameters::HexNut {
        thread,
        across_flats,
    } = &mut incompatible_parameters
    else {
        unreachable!()
    };
    thread.major_diameter = length(6.);
    thread.pitch = length(1.);
    thread.length = length(4.);
    thread.fit = ThreadFit::Metric6HInternal;
    thread.metric_fit_engagement = Some(MetricFitEngagement::NormalGroupAssumed {});
    *across_flats = length(10.);
    let (_, incompatible_definition, _) = author(incompatible_parameters, false);
    let mut incompatible = request.clone();
    incompatible.hardware[2] = library_row(
        &incompatible_definition,
        ScrewConnectionHardwareRole::Nut,
        "Incompatible supported M6 pin",
    );
    let scratch = Session::new(before.clone()).unwrap();
    let candidate = prepare_screw_connection(&scratch, 0, actor(), &incompatible).unwrap();
    let failure = compiler()
        .compile(candidate.prepared.document(), 1)
        .err()
        .expect("Different sourced nominal/pitch must not mate");
    assert!(
        failure.to_string().contains("mating")
            || failure.to_string().contains("nominal")
            || failure.to_string().contains("pitch")
            // A larger incompatible nut can fail actual bearing before the
            // nominal/pitch gate. Both physical refusals preserve this intent.
            || failure.to_string().contains("annular land"),
        "{failure}"
    );
    assert_eq!(scratch.document(), &before);
    assert_eq!(s.document(), &final_document);
    assert_eq!(nut_author, nut_author_bytes);
    assert_eq!(washer_author, washer_author_bytes);
    assert_eq!(symbols::read_definition(&nut_author,nut_author.entities().find(|e| matches!(&e.data,EntityData::Extension{definition,..} if definition.type_id==symbols::DEFINITION_TYPE)).unwrap().id).unwrap(),nut_definition);
}

fn pinned_case(left: bool, posed: bool) -> (Session, ScrewConnectionRequest) {
    let (document, mut request, _) = fixture::fixture();
    let hand = if left {
        ThreadHand::Left
    } else {
        ThreadHand::Right
    };
    for row in &mut request.hardware {
        let ScrewConnectionHardwareSource::Generated { part } = &mut row.source else {
            unreachable!()
        };
        match &mut part.parameters {
            MechanicalPartParameters::MachineScrew { thread, .. }
            | MechanicalPartParameters::HexNut { thread, .. } => thread.hand = hand,
            _ => unreachable!(),
        }
    }
    let ScrewConnectionHardwareSource::Generated { part } = &request.hardware[1].source else {
        unreachable!()
    };
    let (_, nut, _) = author(part.parameters.clone(), false);
    request.hardware[1] = library_row(
        &nut,
        ScrewConnectionHardwareRole::Nut,
        "Pinned current metric nut",
    );
    let (_, washer, _) = author(
        MechanicalPartParameters::Washer {
            bore_diameter: length(3.2),
            outside_diameter: length(7.),
            thickness: length(0.2),
            frame: Default::default(),
        },
        true,
    );
    request.hardware.insert(
        1,
        library_row(
            &washer,
            ScrewConnectionHardwareRole::HeadWasher,
            "Pinned live washer",
        ),
    );
    let mut session = Session::new(document).unwrap();
    if posed {
        let stock = request.receivers[0].source.input.producer;
        let mut primitive = read_primitive(session.document(), stock).unwrap();
        let PrimitiveParameters::Box { frame, .. } = &mut primitive.parameters else {
            unreachable!()
        };
        *frame = PrimitiveFrame {
            origin: [length(100.), length(-30.), length(50.)],
            u: [0.8, 0.6, 0.],
            v: [-0.36, 0.48, 0.8],
        };
        let q = prepare_primitive_edit(
            &session,
            session.revision(),
            actor(),
            stock,
            &primitive.parameters,
        )
        .unwrap();
        session.commit(q.prepared).unwrap();
    }
    (session, request)
}
#[test]
fn current_pinned_metric_material_is_qualified_in_both_hands_and_proper_poses() {
    for left in [false, true] {
        for posed in [false, true] {
            let (session, request) = pinned_case(left, posed);
            let before = session.document().clone();
            let q =
                prepare_screw_connection(&session, session.revision(), actor(), &request).unwrap();
            let model = compiler()
                .compile(q.prepared.document(), session.revision() + 1)
                .unwrap();
            let [check] = model.snapshot().screw_connection_checks.as_slice() else {
                panic!("one actual connection")
            };
            assert!(check.native_component_collision_checked);
            assert!((check.actual_engagement_mm.unwrap() - 2.4).abs() < 1e-7);
            for row in &request.hardware {
                let leaf = leaf(&model, row.instance);
                let solid = model.leaf_solid(&leaf.path).unwrap();
                let receipt = solid
                    .manufacturing_mesh_with_receipt(Default::default())
                    .unwrap();
                assert!(receipt.edge_pairs > 0 && receipt.vertex_links > 0);
                match row.role {
                    ScrewConnectionHardwareRole::Nut => {
                        let (actual, _) = solid.rust_hex_nut_geometry().unwrap();
                        assert_eq!(actual.thread().spec().left_hand, left);
                        assert_eq!(actual.across_flats_mm(), 5.5);
                    }
                    ScrewConnectionHardwareRole::HeadWasher => {
                        let (actual, _) = solid.rust_annulus_data().unwrap();
                        assert_eq!(actual.height_mm(), 0.2);
                    }
                    _ => assert!(solid.rust_machine_screw_geometry().is_some()),
                }
            }
            assert_eq!(session.document(), &before);
        }
    }
}
#[test]
fn pinned_hardware_cannot_waive_actual_bearing_pose_or_current_identity() {
    let (mut session, request) = pinned_case(false, true);
    let q = prepare_screw_connection(&session, session.revision(), actor(), &request).unwrap();
    compiler()
        .compile(q.prepared.document(), session.revision() + 1)
        .unwrap();
    session.commit(q.prepared).unwrap();
    let before = session.document().clone();
    for (index, tilted) in [(1, false), (1, true), (2, false), (2, true)] {
        let row = &request.hardware[index];
        let mut driver = read_screw_hardware_placement(&before, row.driver).unwrap();
        driver.connection = None;
        driver.station =
            length(if index == 1 { -0.2 } else { 4.1 } + if tilted { 0. } else { 0.01 });
        if tilted {
            let a: f64 = 0.01;
            driver.local_mount.rotation = [
                [1., 0., 0.],
                [0., a.portable_cos(), -a.portable_sin()],
                [0., a.portable_sin(), a.portable_cos()],
            ];
        }
        let mut entity = before.get(row.driver).unwrap().clone();
        entity.data = EntityData::Extension {
            definition: driver.definition().unwrap(),
            inputs: driver.inputs().unwrap(),
        };
        let q = session
            .prepare(Request {
                expected_revision: session.revision(),
                actor: actor(),
                label: "Displace actual pinned bearing".into(),
                edits: vec![Mutation::Replace(entity)],
            })
            .unwrap();
        let failure = compiler()
            .compile(q.document(), session.revision() + 1)
            .err()
            .expect("pin cannot replace actual bearing proof");
        assert!(
            failure.to_string().contains("bearing") || failure.to_string().contains("station"),
            "{failure}"
        );
        assert_eq!(session.document(), &before);
    }
    let q = prepare_screw_connection(&session, session.revision(), actor(), &request).unwrap();
    let failure = compiler()
        .compile_with_backend(
            q.prepared.document(),
            session.revision() + 1,
            multicad_model::GeometryBackend::Rust,
            Some(&|| true),
        )
        .err()
        .unwrap();
    assert!(multicad_model::is_superseded(&failure));
    assert_eq!(session.document(), &before);
    let mut mismatched = request.clone();
    let ScrewConnectionHardwareSource::Library { definition, .. } =
        &mut mismatched.hardware[2].source
    else {
        unreachable!()
    };
    definition.name.push_str(" stale replacement");
    assert!(prepare_screw_connection(&session, session.revision(), actor(), &mismatched).is_err());
    assert_eq!(session.document(), &before);
}
