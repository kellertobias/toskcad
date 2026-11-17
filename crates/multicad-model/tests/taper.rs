use multicad_core::{
    commands::{Actor, Channel, Request, Session},
    constraint_sketch::{Curve, CurveKind, source::Source},
    document::{Document, Entity, EntityData, EntityId, Mutation},
    model::{
        self, Direction, Extent, ExtrudeParameters, ExtrudeRequest, FaceInput, PrimitiveParameters,
        PrimitiveRequest, Termination, angle, length,
    },
    references::{GeometryIntent, Origin, TopologyKind},
    sketch::Datum,
    sketch_regions,
};
use multicad_model::{Compiler, Model, ShapeSnapshot};
use std::f64::consts::PI;
fn actor() -> Actor {
    Actor::local(Channel::Mcp)
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 4e-5, "{a} != {b}");
}
fn circle(source: &mut Source, radius: f64) {
    let center = source.add_point([0., 0.]);
    source.sketch.curves.push(Curve {
        id: EntityId::new(),
        construction: false,
        kind: CurveKind::Circle {
            center,
            radius_mm: radius,
        },
    });
}
struct Fixture {
    session: Session,
    request: ExtrudeRequest,
}
impl Fixture {
    fn new(build: impl FnOnce(&mut Source)) -> Self {
        let (sketch, edits) = Source::create(Datum::Xy, "Exact tapered profile").unwrap();
        let mut doc = Document::new("Native taper")
            .unwrap()
            .edited(&edits)
            .unwrap();
        let mut source = Source::read(&doc, sketch).unwrap().0;
        build(&mut source);
        let extraction = sketch_regions::extract(&source);
        assert!(extraction.diagnostics.is_empty());
        let region = extraction
            .regions
            .iter()
            .max_by_key(|r| r.intent().holes.len())
            .unwrap()
            .intent();
        doc = doc.edited(&source.replace(&doc, sketch).unwrap()).unwrap();
        Self {
            session: Session::new(doc).unwrap(),
            request: ExtrudeRequest::new(sketch, vec![region], Default::default()),
        }
    }
    fn apply(&mut self, edits: Vec<Mutation>) {
        self.session
            .apply(Request {
                expected_revision: self.session.revision(),
                actor: actor(),
                label: "Taper fixture".into(),
                edits,
            })
            .unwrap();
    }
    fn commit(&mut self, parameters: ExtrudeParameters) {
        self.request.parameters = parameters;
        let before = self.session.document().clone();
        let p = model::prepare_extrude(
            &self.session,
            self.session.revision(),
            actor(),
            &self.request,
        )
        .unwrap();
        assert_eq!(self.session.document(), &before);
        self.session.commit(p.prepared).unwrap();
    }
    fn compile(&self) -> Model {
        Compiler::default()
            .compile(self.session.document(), self.session.revision())
            .unwrap()
    }
    fn shape<'a>(&self, model: &'a Model) -> &'a ShapeSnapshot {
        &model.snapshot().shapes[&self.request.outputs[0].feature]
    }
    fn parameters(delta: f64) -> ExtrudeParameters {
        ExtrudeParameters {
            taper: angle((delta / 10.).atan().to_degrees()),
            ..Default::default()
        }
    }
}
fn cap(shape: &ShapeSnapshot, role: &str) -> u32 {
    let o = Origin {
        source: shape.producer,
        role: role.into(),
    };
    let caps: Vec<_> = shape
        .topology
        .elements
        .iter()
        .filter(|e| e.key.kind == TopologyKind::Face && e.origins.contains(&o))
        .collect();
    assert_eq!(
        caps.len(),
        1,
        "Logical {role} must identify one genuine cap"
    );
    caps[0].key.ordinal
}
fn radii_at(shape: &ShapeSnapshot, center: [f64; 3]) -> (f64, f64) {
    let radii: Vec<_> = shape
        .edges
        .iter()
        .filter(|e| {
            e.points_mm.len() > 2 && e.points_mm.iter().all(|p| (p[2] - center[2]).abs() < 1e-6)
        })
        .flat_map(|e| &e.points_mm)
        .map(|p| (p[0] - center[0]).hypot(p[1] - center[1]))
        .collect();
    assert!(!radii.is_empty(), "Expected native circular cap edges");
    (
        radii.iter().copied().fold(f64::INFINITY, f64::min),
        radii.iter().copied().fold(f64::NEG_INFINITY, f64::max),
    )
}
#[test]
fn signed_convex_taper_has_exact_edge_offsets_volume_and_source_neutral_rebuild() {
    for delta in [-1., 1.] {
        let mut f = Fixture::new(|s| {
            s.add_rectangle([0., 0.], [10., 6.]).unwrap();
        });
        f.commit(Fixture::parameters(delta));
        let original = f.session.document().clone();
        let model = f.compile();
        let shape = f.shape(&model);
        near(
            shape.volume_mm3,
            10. * (60. + 16. * delta + 4. * delta * delta / 3.),
        );
        near(shape.bounds.min[0], -delta.max(0.));
        near(shape.bounds.max[0], 10. + delta.max(0.));
        assert_eq!(shape.solid_count, 1);
        let end = cap(shape, "cap.end");
        let start = cap(shape, "cap.start");
        let faces = model.solid(shape.producer).unwrap().faces().unwrap();
        near(
            faces
                .iter()
                .find(|face| face.index == end)
                .unwrap()
                .area_mm2,
            (10. + 2. * delta) * (6. + 2. * delta),
        );
        near(
            faces
                .iter()
                .find(|face| face.index == start)
                .unwrap()
                .area_mm2,
            60.,
        );
        let second = f.compile();
        assert_eq!(model.snapshot().stamp, second.snapshot().stamp);
        assert_eq!(shape.topology, f.shape(&second).topology);
        assert_eq!(f.session.document(), &original);
    }
}
#[test]
fn circular_hole_taper_expands_material_and_shrinks_hole_for_both_angle_signs() {
    for delta in [-1., 1.] {
        let mut f = Fixture::new(|s| {
            circle(s, 5.);
            circle(s, 2.);
        });
        f.commit(Fixture::parameters(delta));
        let model = f.compile();
        let shape = f.shape(&model);
        near(shape.volume_mm3, 10. * PI * (21. + 7. * delta));
        let end = cap(shape, "cap.end");
        let faces = model.solid(shape.producer).unwrap().faces().unwrap();
        near(
            faces
                .iter()
                .find(|face| face.index == end)
                .unwrap()
                .area_mm2,
            PI * ((5. + delta).powi(2) - (2. - delta).powi(2)),
        );
        let (inner, outer) = radii_at(shape, [0., 0., 10.]);
        near(inner, 2. - delta);
        near(outer, 5. + delta);
        assert_eq!(shape.solid_count, 1);
        near(shape.bounds.max[0], 5. + delta.max(0.));
    }
}
#[test]
fn tilted_reverse_taper_keeps_normal_start_offset_independent_of_sweep_length() {
    for reverse in [false, true] {
        for delta in [-1., 1.] {
            let mut f = Fixture::new(|s| circle(s, 3.));
            let p = ExtrudeParameters {
                direction: Direction {
                    elevation: angle(0.8_f64.asin().to_degrees()),
                    azimuth: angle(0.),
                    reverse,
                },
                start_offset: length(7.),
                ..Fixture::parameters(delta)
            };
            f.commit(p);
            let model = f.compile();
            let shape = f.shape(&model);
            near(
                shape.volume_mm3,
                8. * PI * (9. + 3. * delta + delta * delta / 3.),
            );
            near(shape.bounds.min[2], if reverse { -1. } else { 7. });
            near(shape.bounds.max[2], if reverse { 7. } else { 15. });
            let end = cap(shape, "cap.end");
            let faces = model.solid(shape.producer).unwrap().faces().unwrap();
            near(
                faces
                    .iter()
                    .find(|face| face.index == end)
                    .unwrap()
                    .area_mm2,
                PI * (3. + delta).powi(2),
            );
            let center = [
                if reverse { -6. } else { 6. },
                0.,
                if reverse { -1. } else { 15. },
            ];
            let (min, max) = radii_at(shape, center);
            near(min, 3. + delta);
            near(max, 3. + delta);
        }
    }
}
#[test]
fn symmetric_and_independent_ends_both_draft_from_original_start_section() {
    for symmetric in [false, true] {
        let forward = if symmetric { 5. } else { 10. };
        let backward = 5.;
        let extent = if symmetric {
            Extent::Symmetric { total: length(10.) }
        } else {
            Extent::TwoSided {
                forward: Termination::Distance {
                    distance: length(forward),
                },
                backward: Termination::Distance {
                    distance: length(backward),
                },
            }
        };
        let mut f = Fixture::new(|s| circle(s, 3.));
        f.commit(ExtrudeParameters {
            extent,
            start_offset: length(3.),
            ..Fixture::parameters(1.)
        });
        let model = f.compile();
        let shape = f.shape(&model);
        let volume = |l: f64| {
            let d = l * 0.1;
            l * PI * (9. + 3. * d + d * d / 3.)
        };
        near(shape.volume_mm3, volume(forward) + volume(backward));
        near(shape.bounds.min[2], 3. - backward);
        near(shape.bounds.max[2], 3. + forward);
        let solid = model.solid(shape.producer).unwrap();
        let faces = solid.faces().unwrap();
        for (role, z, l) in [
            ("cap.start", 3. - backward, backward),
            ("cap.end", 3. + forward, forward),
        ] {
            let id = cap(shape, role);
            let face = faces.iter().find(|face| face.index == id).unwrap();
            near(face.center_mm[2], z);
            near(face.area_mm2, PI * (3. + l * 0.1).powi(2));
            let (min, max) = radii_at(shape, [0., 0., z]);
            near(min, 3. + l * 0.1);
            near(max, 3. + l * 0.1);
        }
        // Native B-rep edge at the meeting plane retains the authored R3;
        // there is no interior cap or a shifted/re-tapered starting section.
        let waist: Vec<_> = shape
            .edges
            .iter()
            .filter(|e| e.points_mm.iter().all(|p| (p[2] - 3.).abs() < 1e-6))
            .collect();
        assert!(!waist.is_empty());
        for p in waist.iter().flat_map(|e| &e.points_mm) {
            near(p[0].hypot(p[1]), 3.);
        }
        assert!(
            !faces
                .iter()
                .any(|face| face.plane.is_some() && (face.center_mm[2] - 3.).abs() < 1e-6)
        );
    }
}
#[test]
fn unsupported_concavity_and_collapsing_outer_or_hole_never_return_default_geometry() {
    for case in 0..3 {
        let mut f = Fixture::new(|s| match case {
            0 => {
                let ids: Vec<_> = [[0., 0.], [4., 0.], [4., 4.], [2., 2.], [0., 4.]]
                    .into_iter()
                    .map(|p| s.add_point(p))
                    .collect();
                for i in 0..ids.len() {
                    s.add_line(ids[i], ids[(i + 1) % ids.len()]);
                }
            }
            1 => circle(s, 3.),
            _ => {
                circle(s, 5.);
                circle(s, 2.);
            }
        });
        let delta = match case {
            0 => 1.,
            1 => -4.,
            _ => 2.,
        };
        f.commit(Fixture::parameters(delta));
        let original = f.session.document().clone();
        let error = Compiler::default()
            .compile(&original, f.session.revision())
            .err()
            .expect("Unsupported taper must fail");
        assert_eq!(error.entity, Some(f.request.outputs[0].feature));
        assert!(!error.message.is_empty());
        assert_eq!(f.session.document(), &original);
    }
}
#[test]
fn taper_with_plane_or_face_termination_fails_explicitly_before_native_sweep() {
    for face_target in [false, true] {
        let mut f = Fixture::new(|s| {
            s.add_rectangle([0., 0.], [10., 6.]).unwrap();
        });
        let end = if face_target {
            let target = PrimitiveRequest::new(
                "Target",
                PrimitiveParameters::Box {
                    size: [length(10.), length(6.), length(20.)],
                    frame: Default::default(),
                },
            );
            let prepared = model::prepare_primitive(&f.session, 0, actor(), &target).unwrap();
            f.session.commit(prepared.prepared).unwrap();
            let initial = f.compile();
            let shape = &initial.snapshot().shapes[&target.feature];
            let key = shape
                .topology
                .elements
                .iter()
                .find(|e| {
                    e.key.kind == TopologyKind::Face
                        && e.origins.contains(&Origin {
                            source: target.feature,
                            role: "cap.end".into(),
                        })
                })
                .unwrap()
                .key;
            let intent = GeometryIntent::capture(&shape.topology, key).unwrap();
            let reference = intent.id;
            f.request.faces = vec![FaceInput {
                source: target.feature,
                intent,
            }];
            Termination::ToFace {
                reference,
                offset: length(0.),
            }
        } else {
            let plane = Entity::new(
                "End",
                EntityData::Plane {
                    origin_mm: [0., 0., 20.],
                    normal: [0., 0., 1.],
                    anchor: None,
                },
            );
            let id = plane.id;
            f.apply(vec![Mutation::Insert(plane)]);
            Termination::ToPlane {
                plane: id,
                offset: length(0.),
            }
        };
        f.commit(ExtrudeParameters {
            extent: Extent::OneSided { end },
            ..Fixture::parameters(1.)
        });
        let original = f.session.document().clone();
        let error = Compiler::default()
            .compile(&original, f.session.revision())
            .err()
            .expect("Referenced taper must fail");
        assert_eq!(error.entity, Some(f.request.outputs[0].feature));
        assert!(
            error.message.contains("Taper with an associative"),
            "{}",
            error.message
        );
        assert_eq!(f.session.document(), &original);
    }
}
/// Closed, consistently oriented manufacturing triangles whose enclosed volume
/// matches the exact solid within the requested chordal accuracy.
fn assert_closed_accurate_mesh(solid: &multicad_kernel::Solid, exact: f64, surface: f64) {
    let deflection = 0.01;
    let receipt = solid
        .manufacturing_mesh_with_receipt(multicad_kernel::MeshOptions {
            linear_deflection_mm: deflection,
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
    // Every edge has one use in each direction: closed and oriented.
    assert!(directed.values().all(|&uses| uses == 0));
    assert!(receipt.maximum_displacement_mm <= deflection);
    assert!(
        (volume - exact).abs() <= deflection * surface,
        "{volume} != {exact}"
    );
}
#[test]
fn curved_oblique_taper_keeps_identity_through_edits_undo_redo_and_cold_reopen() {
    let mut f = Fixture::new(|s| {
        circle(s, 5.);
        circle(s, 2.);
    });
    let tilted = |delta: f64| ExtrudeParameters {
        direction: Direction {
            elevation: angle(0.8_f64.asin().to_degrees()),
            azimuth: angle(0.),
            reverse: false,
        },
        ..Fixture::parameters(delta)
    };
    // Oblique frustum with a circular hole, normal height eight.
    let expected = |delta: f64| {
        let (outer, hole) = (5. + delta, 2. - delta);
        8. * PI / 3. * ((25. + 5. * outer + outer * outer) - (4. + 2. * hole + hole * hole))
    };
    f.commit(tilted(1.));
    let feature = f.request.outputs[0].feature;
    let check = |document: &Document, revision: u64, delta: f64| {
        let model = Compiler::default().compile(document, revision).unwrap();
        let shape = &model.snapshot().shapes[&feature];
        assert_eq!(shape.producer, feature);
        assert_eq!(shape.solid_count, 1);
        near(shape.volume_mm3, expected(delta));
        let end = cap(shape, "cap.end");
        cap(shape, "cap.start");
        let solid = model.solid(feature).unwrap();
        let faces = solid.faces().unwrap();
        let cap_face = faces.iter().find(|face| face.index == end).unwrap();
        near(
            cap_face.area_mm2,
            PI * ((5. + delta).powi(2) - (2. - delta).powi(2)),
        );
        // The end section sits on the sloped generator: center (6, 0, 8).
        for (axis, value) in [6., 0., 8.].into_iter().enumerate() {
            near(cap_face.center_mm[axis], value);
        }
        let (inner, outer) = radii_at(shape, [6., 0., 8.]);
        near(inner, 2. - delta);
        near(outer, 5. + delta);
        shape.topology.validate().unwrap();
        let area: f64 = faces.iter().map(|face| face.area_mm2).sum();
        assert_closed_accurate_mesh(solid, expected(delta), area);
        shape.topology.elements.clone()
    };
    let first = check(f.session.document(), f.session.revision(), 1.);

    // A collapsing hole is refused by the compiler without publishing.
    let revision = f.session.revision();
    let invalid =
        model::prepare_extrude_edit(&f.session, revision, actor(), feature, &tilted(3.)).unwrap();
    f.session.commit(invalid.prepared).unwrap();
    let error = Compiler::default()
        .compile(f.session.document(), f.session.revision())
        .err()
        .expect("A collapsed hole must not produce geometry");
    assert_eq!(error.entity, Some(feature));
    f.session.undo(f.session.revision(), actor()).unwrap();
    assert_eq!(check(f.session.document(), f.session.revision(), 1.), first);

    let revision = f.session.revision();
    let edit =
        model::prepare_extrude_edit(&f.session, revision, actor(), feature, &tilted(-1.)).unwrap();
    f.session.commit(edit.prepared).unwrap();
    let narrowed = check(f.session.document(), f.session.revision(), -1.);
    let keys = |elements: &[multicad_core::references::TopologyElement]| {
        elements
            .iter()
            .filter(|e| e.key.kind == TopologyKind::Face)
            .map(|e| (e.key, e.origins.clone()))
            .collect::<Vec<_>>()
    };
    // Face identities survive the dimension change.
    assert_eq!(keys(&narrowed), keys(&first));
    f.session.undo(f.session.revision(), actor()).unwrap();
    assert_eq!(check(f.session.document(), f.session.revision(), 1.), first);
    f.session.redo(f.session.revision(), actor()).unwrap();
    assert_eq!(
        check(f.session.document(), f.session.revision(), -1.),
        narrowed
    );

    let project =
        multicad_core::project_file::ProjectFile::new(f.session.document().clone()).unwrap();
    let mut bytes = vec![];
    multicad_core::project_file::write(&project, &mut bytes, Default::default()).unwrap();
    let supported = multicad_core::project_file::supported_capabilities();
    let opened =
        multicad_core::project_file::read(&mut bytes.as_slice(), Default::default(), &supported)
            .unwrap()
            .project;
    assert_eq!(opened.document, project.document);
    assert_eq!(
        model::read_extrude(&opened.document, feature)
            .unwrap()
            .parameters,
        tilted(-1.)
    );
    assert_eq!(check(&opened.document, 0, -1.), narrowed);
}
