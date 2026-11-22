use multicad_core::{
    commands::{Actor, Channel, Request, Session},
    constraint_sketch::{Curve, CurveKind, source::Source},
    document::{Document, Entity, EntityData, EntityId, Mutation},
    model::{
        self, BodyInput, BooleanRequest, Direction, Extent, ExtrudeParameters, ExtrudeRequest,
        FaceInput, Operation, PrimitiveFrame, PrimitiveParameters, PrimitiveRequest, Termination,
        angle, length,
    },
    references::{Attachment, GeometryIntent, Origin, TopologyKind},
    sketch::Datum,
    sketch_regions,
};
use multicad_kernel::{GeometryPolicy, Solid};
use multicad_model::{Compiler, GeometryBackend, Model, ModelError, ShapeSnapshot};
use multicad_precision::PortableMath as _;
fn actor() -> Actor {
    Actor::local(Channel::Mcp)
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 3e-5, "{a} != {b}");
}
struct Fixture {
    session: Session,
    request: ExtrudeRequest,
}
impl Fixture {
    fn new() -> Self {
        Self::rectangle(false)
    }
    fn rectangle(with_hole: bool) -> Self {
        Self::profile(usize::from(with_hole), |source| {
            source.add_rectangle([0., 0.], [10., 6.]).unwrap();
            if with_hole {
                source.add_rectangle([3., 2.], [7., 4.]).unwrap();
            }
        })
    }
    fn circle() -> Self {
        Self::profile(0, |source| {
            let center = source.add_point([5., 4.]);
            source.sketch.curves.push(Curve {
                id: EntityId::new(),
                construction: false,
                kind: CurveKind::Circle {
                    center,
                    radius_mm: 2.,
                },
            });
        })
    }
    fn profile(holes: usize, build: impl FnOnce(&mut Source)) -> Self {
        let (sketch, edits) = Source::create(Datum::Xy, "Terminated profile").unwrap();
        let mut doc = Document::new("Native terminations")
            .unwrap()
            .edited(&edits)
            .unwrap();
        let mut source = Source::read(&doc, sketch).unwrap().0;
        build(&mut source);
        let region = sketch_regions::extract(&source)
            .regions
            .into_iter()
            .find(|region| region.holes.len() == holes)
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
                label: "Termination fixture edit".into(),
                edits,
            })
            .unwrap();
    }
    fn plane(&mut self, height: f64, slope: f64, reverse: bool) -> EntityId {
        let scale = 1_f64.portable_hypot(slope);
        let sign = if reverse { -1. } else { 1. };
        let entity = Entity::new(
            "Target plane",
            EntityData::Plane {
                origin_mm: [0., 0., height],
                normal: [-slope / scale * sign, 0., 1. / scale * sign],
                anchor: None,
            },
        );
        let id = entity.id;
        self.apply(vec![Mutation::Insert(entity)]);
        id
    }
    fn commit(&mut self) {
        let prepared = model::prepare_extrude(
            &self.session,
            self.session.revision(),
            actor(),
            &self.request,
        )
        .unwrap();
        self.session.commit(prepared.prepared).unwrap();
    }
    fn compile(&self) -> Model {
        self.compile_rust().unwrap()
    }
    fn compile_rust(&self) -> Result<Model, ModelError> {
        self.compile_with_backend(GeometryBackend::Rust)
    }
    fn compile_with_backend(&self, backend: GeometryBackend) -> Result<Model, ModelError> {
        Compiler::default().compile_with_backend(
            self.session.document(),
            self.session.revision(),
            backend,
            None,
        )
    }
    /// Check one Rust result against an analytic expectation of its whole
    /// box, its volume and the first (support) frame of each logical cap.
    fn verify(&self, model: &Model, expected: &Expected) {
        let feature = self.request.outputs[0].feature;
        let solid = model.solid(feature).unwrap();
        assert!(solid.uses_rust_geometry());
        assert_eq!(solid.solid_count().unwrap(), 1);
        assert_eq!(
            model.snapshot().bodies.get(&self.request.outputs[0].body),
            Some(&feature)
        );
        let actual = self.shape(model);
        near(actual.volume_mm3, expected.volume);
        for i in 0..3 {
            near(actual.bounds.min[i], expected.min[i]);
            near(actual.bounds.max[i], expected.max[i]);
        }
        for (role, cap_expected) in ["cap.start", "cap.end"].into_iter().zip(&expected.caps) {
            let frame = cap(actual, role).frames[0].frame;
            for (a, b) in frame
                .origin_mm
                .into_iter()
                .chain(frame.u)
                .chain(frame.normal())
                .zip(
                    cap_expected
                        .origin
                        .into_iter()
                        .chain(cap_expected.u)
                        .chain(cap_expected.normal),
                )
            {
                near(a, b);
            }
        }
    }
    fn shape<'a>(&self, model: &'a Model) -> &'a ShapeSnapshot {
        &model.snapshot().shapes[&self.request.outputs[0].feature]
    }
    fn primitive(&mut self, parameters: PrimitiveParameters) -> PrimitiveRequest {
        let request = PrimitiveRequest::new("Target solid", parameters);
        let prepared =
            model::prepare_primitive(&self.session, self.session.revision(), actor(), &request)
                .unwrap();
        self.session.commit(prepared.prepared).unwrap();
        request
    }
}

/// Analytic first frame of one logical cap: its origin, in-plane U and
/// normal. A cap keeps the frame of its generator, so a sketch-side cap has
/// the working sketch frame and a plane cap has the target plane's own frame,
/// whichever way the material lies.
struct CapFrame {
    origin: [f64; 3],
    u: [f64; 3],
    normal: [f64; 3],
}
struct Expected {
    volume: f64,
    min: [f64; 3],
    max: [f64; 3],
    caps: [CapFrame; 2],
}
fn unit(v: [f64; 3]) -> [f64; 3] {
    let length = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    v.map(|x| x / length)
}
/// The fixture's target-plane normal, as authored in [`Fixture::plane`].
fn plane_normal(slope: f64, reverse: bool) -> [f64; 3] {
    let sign = if reverse { -1. } else { 1. };
    unit([-slope * sign, 0., sign])
}
/// A construction plane has no intrinsic roll: its documented U is the
/// horizontal (-ny, nx, 0) direction, or projected +X for a horizontal plane.
fn plane_u(normal: [f64; 3]) -> [f64; 3] {
    if normal[0].portable_hypot(normal[1]) > 1e-7 {
        unit([-normal[1], normal[0], 0.])
    } else {
        unit([
            1. - normal[0] * normal[0],
            -normal[0] * normal[1],
            -normal[0] * normal[2],
        ])
    }
}
/// A cap from the XY working sketch (start offset or distance end) keeps the
/// sketch frame translated along +Z.
fn sketch_cap(z: f64) -> CapFrame {
    CapFrame {
        origin: [0., 0., z],
        u: [1., 0., 0.],
        normal: [0., 0., 1.],
    }
}
/// A cap on a fixture target plane through (0, 0, height), offset along its
/// own normal.
fn plane_cap(height: f64, slope: f64, reverse: bool, offset: f64) -> CapFrame {
    let normal = plane_normal(slope, reverse);
    CapFrame {
        origin: [normal[0] * offset, 0., height + normal[2] * offset],
        u: plane_u(normal),
        normal,
    }
}

fn check_curved_step(solid: &Solid, feature: EntityId, volume: f64) {
    let folder = tempfile::tempdir().unwrap();
    let path = folder.path().join("termination.step");
    solid.write_step(&path).unwrap();
    let policy = GeometryPolicy::default();
    let check = |readback: Solid| {
        assert_eq!(readback.solid_count().unwrap(), 1);
        near(readback.volume_mm3().unwrap(), volume);
        let expected = solid.bounds().unwrap();
        let actual = readback.bounds().unwrap();
        for axis in 0..3 {
            near(actual.min[axis], expected.min[axis]);
            near(actual.max[axis], expected.max[axis]);
        }
    };
    check(
        Solid::import_rust_step(&std::fs::read_to_string(&path).unwrap(), feature, policy).unwrap(),
    );
}

#[test]
fn rust_circular_extrusion_to_face_follows_the_reference_and_keeps_a_through_cut() {
    let mut f = Fixture::circle();
    let target = f.primitive(PrimitiveParameters::Box {
        size: [length(12.), length(10.), length(8.)],
        frame: Default::default(),
    });
    let initial = f.compile_rust().unwrap();
    let target_shape = &initial.snapshot().shapes[&target.feature];
    let target_cap = cap(target_shape, "cap.end");
    let original_anchor = &target_cap.frames[0];
    let mut intent = GeometryIntent::capture(&target_shape.topology, target_cap.key).unwrap();
    intent.attachment = Some(Attachment {
        anchor: original_anchor.origin.clone(),
        offset_mm: [1.5, 2., 0.],
        angle_rad: std::f64::consts::FRAC_PI_2,
    });
    f.request.parameters.extent = Extent::OneSided {
        end: Termination::ToFace {
            reference: intent.id,
            offset: length(0.),
        },
    };
    f.request.faces = vec![FaceInput {
        source: target.feature,
        intent,
    }];
    f.commit();
    let tool = &f.request.outputs[0];
    let boolean = BooleanRequest::new(
        Operation::Cut,
        BodyInput {
            body: target.body,
            producer: target.feature,
        },
        vec![BodyInput {
            body: tool.body,
            producer: tool.feature,
        }],
    );
    let feature = tool.feature;
    let prepared =
        model::prepare_boolean(&f.session, f.session.revision(), actor(), &boolean).unwrap();
    f.session.commit(prepared.prepared).unwrap();
    for height in [8., 11.] {
        if height != 8. {
            let edit = model::prepare_primitive_edit(
                &f.session,
                f.session.revision(),
                actor(),
                target.feature,
                &PrimitiveParameters::Box {
                    size: [length(12.), length(10.), length(height)],
                    frame: Default::default(),
                },
            )
            .unwrap();
            f.session.commit(edit.prepared).unwrap();
        }
        let rust = f.compile_rust().unwrap();
        let solid = rust.solid(feature).unwrap();
        assert!(solid.uses_rust_geometry());
        let volume = std::f64::consts::PI * 4. * height;
        near(solid.volume_mm3().unwrap(), volume);
        near(f.shape(&rust).bounds.max[2], height);
        let anchor = cap(f.shape(&rust), "cap.end").frames[0].frame;
        let basis = original_anchor.frame;
        for i in 0..3 {
            near(
                anchor.origin_mm[i],
                basis.origin_mm[i]
                    + 1.5 * basis.u[i]
                    + 2. * basis.v[i]
                    + if i == 2 { height - 8. } else { 0. },
            );
            near(anchor.u[i], basis.v[i]);
            near(anchor.v[i], -basis.u[i]);
        }
        near(
            rust.solid(boolean.feature).unwrap().volume_mm3().unwrap(),
            120. * height - volume,
        );
        check_curved_step(solid, feature, volume);
    }
}

#[test]
fn rust_circular_extrusion_to_plane_carries_its_cap_anchor_through_rigid_placement() {
    let mut f = Fixture::circle();
    let plane = f.plane(-9., 0., true);
    let mut entity = f.session.document().get(plane).unwrap().clone();
    let EntityData::Plane { origin_mm, .. } = &mut entity.data else {
        panic!("Expected a construction plane");
    };
    *origin_mm = [3., 6., -9.];
    f.apply(vec![Mutation::Replace(entity)]);
    f.request.parameters.direction.reverse = true;
    f.request.parameters.extent = Extent::OneSided {
        end: Termination::ToPlane {
            plane,
            offset: length(0.),
        },
    };
    f.commit();
    let rust = f.compile_rust().unwrap();
    let shape = f.shape(&rust);
    let cap_end = cap(shape, "cap.end");
    let anchor = cap_end.frames[0].frame;
    assert_eq!(anchor.origin_mm, [3., 6., -9.]);
    assert_eq!(anchor.u, [1., 0., 0.]);
    assert_eq!(anchor.v, [0., -1., 0.]);
    let placed = rust
        .solid(shape.producer)
        .unwrap()
        .transformed(multicad_kernel::RigidTransform {
            rotation: [[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]],
            translation_mm: [7., 8., 9.],
        })
        .unwrap();
    let topology = placed.topology_snapshot(shape.producer, 0).unwrap();
    let anchor = topology
        .elements
        .iter()
        .find(|element| element.key == cap_end.key)
        .unwrap()
        .frames[0]
        .frame;
    assert_eq!(anchor.origin_mm, [1., 11., 0.]);
    assert_eq!(anchor.u, [0., 1., 0.]);
    assert_eq!(anchor.v, [1., 0., 0.]);
}

#[test]
fn rust_arc_extrusion_between_parallel_caps_preserves_reverse_sides_and_step() {
    let mut f = Fixture::profile(0, |source| {
        let a = source.add_point([0., 0.]);
        let b = source.add_point([8., 0.]);
        let c = source.add_point([8., 6.]);
        let d = source.add_point([0., 6.]);
        let center = source.add_point([8., 3.]);
        source.add_line(a, b);
        source.sketch.curves.push(Curve {
            id: EntityId::new(),
            construction: false,
            kind: CurveKind::Arc {
                center,
                start: b,
                end: c,
                clockwise: false,
            },
        });
        source.add_line(c, d);
        source.add_line(d, a);
    });
    let plane = f.plane(-8., 0., true);
    f.request.parameters.start_offset = length(1.);
    f.request.parameters.direction.reverse = true;
    f.request.parameters.extent = Extent::TwoSided {
        forward: Termination::ToPlane {
            plane,
            offset: length(2.),
        },
        backward: Termination::Distance {
            distance: length(3.),
        },
    };
    f.commit();
    let rust = f.compile_rust().unwrap();
    let shape = f.shape(&rust);
    let volume = (48. + 4.5 * std::f64::consts::PI) * 14.;
    near(shape.volume_mm3, volume);
    near(shape.bounds.min[2], -10.);
    near(shape.bounds.max[2], 4.);
    near(cap(shape, "cap.start").frames[0].frame.origin_mm[2], 4.);
    near(cap(shape, "cap.end").frames[0].frame.origin_mm[2], -10.);
    check_curved_step(rust.solid(shape.producer).unwrap(), shape.producer, volume);
}

#[test]
fn rust_curved_termination_refuses_sloped_caps_oblique_travel_or_reversed_extent() {
    for (height, slope, elevation, message) in [
        (8., 1e-9, 90., "parallel caps"),
        (8., 0., 80., "normal travel"),
        (-8., 0., 90., "sufficient material"),
    ] {
        let mut f = Fixture::circle();
        let plane = f.plane(height, slope, false);
        f.request.parameters.direction.elevation = angle(elevation);
        f.request.parameters.extent = Extent::OneSided {
            end: Termination::ToPlane {
                plane,
                offset: length(0.),
            },
        };
        f.commit();
        let error = f.compile_rust().err().expect("Unsupported curved cap");
        assert_eq!(error.entity, Some(f.request.outputs[0].feature));
        assert!(error.message.contains(message), "{error}");
    }
}
fn cap<'a>(shape: &'a ShapeSnapshot, role: &str) -> &'a multicad_core::references::TopologyElement {
    let origin = Origin {
        source: shape.producer,
        role: role.into(),
    };
    let caps: Vec<_> = shape
        .topology
        .elements
        .iter()
        .filter(|e| e.key.kind == TopologyKind::Face && e.origins.contains(&origin))
        .collect();
    assert_eq!(caps.len(), 1, "Cap {role} is missing or ambiguous");
    caps[0]
}
#[test]
fn nonparallel_plane_offset_and_later_plane_move_rebuild_real_geometry() {
    let mut f = Fixture::new();
    let plane = f.plane(12., 0.5, false);
    f.request.parameters.start_offset = length(3.);
    f.request.parameters.extent = Extent::OneSided {
        end: Termination::ToPlane {
            plane,
            offset: length(2.),
        },
    };
    f.commit();
    let original = f.session.document().clone();
    let source = model::read_extrude(&original, f.request.outputs[0].feature).unwrap();
    let first = f.compile();
    let shape = f.shape(&first);
    let intercept = 12. + 2. * 1_f64.portable_hypot(0.5);
    near(shape.volume_mm3, 60. * (intercept + 2.5 - 3.));
    near(shape.bounds.min[2], 3.);
    near(shape.bounds.max[2], intercept + 5.);
    let end = cap(shape, "cap.end");
    assert!(!end.frames.is_empty());
    for frame in &end.frames {
        near(
            frame.frame.origin_mm[2] - 0.5 * frame.frame.origin_mm[0],
            intercept,
        );
    }
    let sloped = |height: f64| {
        let intercept = height + 2. * 1_f64.portable_hypot(0.5);
        Expected {
            volume: 60. * (intercept + 2.5 - 3.),
            min: [0., 0., 3.],
            max: [10., 6., intercept + 5.],
            caps: [sketch_cap(3.), plane_cap(height, 0.5, false, 2.)],
        }
    };
    f.verify(&first, &sloped(12.));
    assert_eq!(f.session.document(), &original);
    let mut entity = original.get(plane).unwrap().clone();
    let EntityData::Plane { origin_mm, .. } = &mut entity.data else {
        panic!()
    };
    origin_mm[2] = 20.;
    f.apply(vec![Mutation::Replace(entity)]);
    let second = f.compile();
    assert_eq!(second.snapshot().bodies, first.snapshot().bodies);
    near(
        f.shape(&second).volume_mm3,
        60. * (20. + 2. * 1_f64.portable_hypot(0.5) + 2.5 - 3.),
    );
    f.verify(&second, &sloped(20.));
    assert_eq!(
        model::read_extrude(f.session.document(), f.request.outputs[0].feature).unwrap(),
        source
    );
}
#[test]
fn target_normal_offset_and_tilt_are_independent_of_start_plane_offset() {
    let mut f = Fixture::new();
    let plane = f.plane(20., 0., true);
    f.request.parameters.extent = Extent::OneSided {
        end: Termination::ToPlane {
            plane,
            offset: length(2.),
        },
    };
    f.commit();
    let model = f.compile();
    near(f.shape(&model).volume_mm3, 60. * 18.);
    near(f.shape(&model).bounds.max[2], 18.);
    f.verify(
        &model,
        &Expected {
            volume: 60. * 18.,
            min: [0.; 3],
            max: [10., 6., 18.],
            caps: [sketch_cap(0.), plane_cap(20., 0., true, 2.)],
        },
    );
    let mut f = Fixture::new();
    let plane = f.plane(20., 0., false);
    f.request.parameters = ExtrudeParameters {
        direction: Direction {
            elevation: angle(45.),
            azimuth: angle(0.),
            reverse: false,
        },
        start_offset: length(3.),
        extent: Extent::OneSided {
            end: Termination::ToPlane {
                plane,
                offset: length(0.),
            },
        },
        ..Default::default()
    };
    f.commit();
    let model = f.compile();
    near(f.shape(&model).volume_mm3, 60. * 17.);
    near(f.shape(&model).bounds.min[0], 0.);
    near(f.shape(&model).bounds.max[0], 27.);
    near(f.shape(&model).bounds.min[2], 3.);
    near(f.shape(&model).bounds.max[2], 20.);
    f.verify(
        &model,
        &Expected {
            volume: 60. * 17.,
            min: [0., 0., 3.],
            max: [27., 6., 20.],
            caps: [sketch_cap(3.), plane_cap(20., 0., false, 0.)],
        },
    );
}
#[test]
fn behind_parallel_crossing_and_zero_height_targets_fail_instead_of_guessing_distance() {
    for (height, slope, parallel) in [
        (-1., 0., false),
        (0., 0., false),
        (-5., 1., false),
        (10., 0., true),
    ] {
        let mut f = Fixture::new();
        let plane = f.plane(height, slope, false);
        if parallel {
            let mut entity = f.session.document().get(plane).unwrap().clone();
            let EntityData::Plane { normal, .. } = &mut entity.data else {
                panic!()
            };
            *normal = [1., 0., 0.];
            f.apply(vec![Mutation::Replace(entity)]);
        }
        f.request.parameters.extent = Extent::OneSided {
            end: Termination::ToPlane {
                plane,
                offset: length(0.),
            },
        };
        f.commit();
        let original = f.session.document().clone();
        let error = Compiler::default()
            .compile(&original, f.session.revision())
            .err()
            .expect("Unreachable termination must not produce a default solid");
        assert_eq!(error.entity, Some(f.request.outputs[0].feature));
        let rust_error = f
            .compile_rust()
            .err()
            .expect("Rust must reject the same plane");
        assert_eq!(rust_error.entity, error.entity);
        assert_eq!(f.session.document(), &original);
    }
}
#[test]
fn independent_plane_and_distance_ends_have_one_truthful_logical_cap_per_side() {
    for case in 0..3 {
        let mut f = Fixture::new();
        let forward = f.plane(15., 0., false);
        let backward = f.plane(-5., 0., false);
        f.request.parameters.start_offset = length(2.);
        let (a, b, min, max, caps) = match case {
            0 => (
                Termination::ToPlane {
                    plane: forward,
                    offset: length(0.),
                },
                Termination::ToPlane {
                    plane: backward,
                    offset: length(0.),
                },
                -5.,
                15.,
                [plane_cap(-5., 0., false, 0.), plane_cap(15., 0., false, 0.)],
            ),
            1 => (
                Termination::Distance {
                    distance: length(7.),
                },
                Termination::ToPlane {
                    plane: backward,
                    offset: length(0.),
                },
                -5.,
                9.,
                [plane_cap(-5., 0., false, 0.), sketch_cap(9.)],
            ),
            _ => (
                Termination::ToPlane {
                    plane: forward,
                    offset: length(0.),
                },
                Termination::Distance {
                    distance: length(4.),
                },
                -2.,
                15.,
                [sketch_cap(-2.), plane_cap(15., 0., false, 0.)],
            ),
        };
        f.request.parameters.extent = Extent::TwoSided {
            forward: a,
            backward: b,
        };
        f.commit();
        let model = f.compile();
        let shape = f.shape(&model);
        near(shape.volume_mm3, 60. * (max - min));
        near(shape.bounds.min[2], min);
        near(shape.bounds.max[2], max);
        for (role, z) in [("cap.start", min), ("cap.end", max)] {
            let end = cap(shape, role);
            assert!(!end.frames.is_empty());
            for frame in &end.frames {
                near(frame.frame.origin_mm[2], z);
            }
        }
        f.verify(
            &model,
            &Expected {
                volume: 60. * (max - min),
                min: [0., 0., min],
                max: [10., 6., max],
                caps,
            },
        );
    }
}

#[test]
fn oblique_two_sided_plane_sweep_preserves_the_hole_on_both_sloped_caps() {
    let mut f = Fixture::rectangle(true);
    let forward = f.plane(15., 0.2, false);
    let backward = f.plane(-8., -0.1, true);
    f.request.parameters.direction = Direction {
        elevation: angle(45.),
        azimuth: angle(0.),
        reverse: false,
    };
    f.request.parameters.start_offset = length(2.);
    f.request.parameters.extent = Extent::TwoSided {
        forward: Termination::ToPlane {
            plane: forward,
            offset: length(0.),
        },
        backward: Termination::ToPlane {
            plane: backward,
            offset: length(0.),
        },
    };
    f.commit();
    let rust = f.compile();
    let shape = f.shape(&rust);
    // Integrate the affine travel gap over the source face. Its centroid
    // stays at (5, 3) because the inner rectangle is concentric.
    let front_height = (15. + 0.2 * 5. - 2.) / (1. - 0.2);
    let back_height = (-8. - 0.1 * 5. - 2.) / (1. + 0.1);
    near(shape.volume_mm3, 52. * (front_height - back_height));
    // Generators leave the z = 2 section along (1, 0, 1)/sqrt(2); the rise
    // to z - 0.2 x = 15 and the fall to z + 0.1 x = -8 are affine in x, so the
    // outer corners at x = 0 and x = 10 bound the sheared prism.
    let mut min = [f64::INFINITY, 0., f64::INFINITY];
    let mut max = [f64::NEG_INFINITY, 6., f64::NEG_INFINITY];
    for x in [0., 10.] {
        for rise in [(13. + 0.2 * x) / 0.8, -(10. + 0.1 * x) / 1.1] {
            for (i, value) in [(0, x + rise), (2, 2. + rise)] {
                min[i] = min[i].min(value);
                max[i] = max[i].max(value);
            }
        }
    }
    f.verify(
        &rust,
        &Expected {
            volume: 52. * (front_height - back_height),
            min,
            max,
            caps: [
                plane_cap(-8., -0.1, true, 0.),
                plane_cap(15., 0.2, false, 0.),
            ],
        },
    );
    let solid = rust.solid(shape.producer).unwrap();
    for role in ["cap.start", "cap.end"] {
        let face = cap(shape, role);
        let frame = face.frames[0].frame;
        let contours = solid
            .planar_face_contours(
                face.key.ordinal,
                multicad_kernel::PlaneFrame {
                    origin: frame.origin_mm,
                    u: frame.u,
                    v: frame.v,
                },
            )
            .unwrap();
        assert_eq!(contours.loops.len(), 2);
    }
}

#[test]
fn two_sided_plane_at_the_working_sketch_cannot_become_a_one_sided_sweep() {
    let mut f = Fixture::new();
    let backward = f.plane(2., 0., false);
    f.request.parameters.start_offset = length(2.);
    f.request.parameters.extent = Extent::TwoSided {
        forward: Termination::Distance {
            distance: length(7.),
        },
        backward: Termination::ToPlane {
            plane: backward,
            offset: length(0.),
        },
    };
    f.commit();
    let error = f
        .compile_rust()
        .err()
        .expect("Zero backward span is invalid");
    assert_eq!(error.entity, Some(f.request.outputs[0].feature));
}

#[test]
fn face_termination_follows_a_real_primitive_cap_height_and_rejects_missing_reference() {
    let mut f = Fixture::new();
    let target = f.primitive(PrimitiveParameters::Box {
        size: [length(30.), length(30.), length(20.)],
        frame: Default::default(),
    });
    let initial = f.compile();
    let target_shape = &initial.snapshot().shapes[&target.feature];
    let key = cap(target_shape, "cap.end").key;
    let intent = GeometryIntent::capture(&target_shape.topology, key).unwrap();
    f.request.parameters.extent = Extent::OneSided {
        end: Termination::ToFace {
            reference: intent.id,
            offset: length(0.),
        },
    };
    f.request.faces = vec![FaceInput {
        source: target.feature,
        intent,
    }];
    f.commit();
    let original = f.session.document().clone();
    let first = f.compile();
    near(f.shape(&first).volume_mm3, 1200.);
    near(f.shape(&first).bounds.max[2], 20.);
    let face_end = |height: f64| Expected {
        volume: 60. * height,
        min: [0.; 3],
        max: [10., 6., height],
        caps: [
            sketch_cap(0.),
            CapFrame {
                origin: [0., 0., height],
                u: [1., 0., 0.],
                normal: [0., 0., 1.],
            },
        ],
    };
    f.verify(&first, &face_end(20.));
    assert_eq!(f.session.document(), &original);
    let parameters = PrimitiveParameters::Box {
        size: [length(30.), length(30.), length(35.)],
        frame: Default::default(),
    };
    let prepared = model::prepare_primitive_edit(
        &f.session,
        f.session.revision(),
        actor(),
        target.feature,
        &parameters,
    )
    .unwrap();
    f.session.commit(prepared.prepared).unwrap();
    let updated = f.compile();
    near(f.shape(&updated).volume_mm3, 2100.);
    near(f.shape(&updated).bounds.max[2], 35.);
    f.verify(&updated, &face_end(35.));
    assert_eq!(updated.snapshot().bodies, first.snapshot().bodies);
    assert_eq!(
        f.session.document().get(f.request.sketch),
        original.get(f.request.sketch)
    );
    let feature = f.request.outputs[0].feature;
    let mut entity = f.session.document().get(feature).unwrap().clone();
    let EntityData::Feature { inputs, .. } = &mut entity.data else {
        panic!()
    };
    inputs
        .iter_mut()
        .find(|l| l.geometry.is_some())
        .unwrap()
        .geometry
        .as_mut()
        .unwrap()
        .origins = std::collections::BTreeSet::from([Origin {
        source: target.feature,
        role: "missing.cap".into(),
    }]);
    f.apply(vec![Mutation::Replace(entity)]);
    let error = Compiler::default()
        .compile(f.session.document(), f.session.revision())
        .err()
        .expect("Missing face reference cannot fall back to a saved plane");
    assert_eq!(error.entity, Some(feature));
    assert!(error.message.contains("unresolved"));
    let rust_error = f
        .compile_rust()
        .err()
        .expect("Rust must reject the missing reference");
    assert_eq!(rust_error.entity, Some(feature));
    assert!(rust_error.message.contains("unresolved"));
}
#[test]
fn genuine_curved_face_is_rejected_as_a_planar_termination_target() {
    {
        let backend = GeometryBackend::Rust;
        let mut f = Fixture::new();
        let target = f.primitive(PrimitiveParameters::Cylinder {
            radius: length(20.),
            height: length(20.),
            frame: Default::default(),
        });
        let initial = f.compile_with_backend(backend).unwrap();
        let snapshot = &initial.snapshot().shapes[&target.feature].topology;
        let solid = initial.solid(target.feature).unwrap();
        let faces = solid.faces().unwrap();
        let key = snapshot
            .elements
            .iter()
            .find(|e| {
                e.key.kind == TopologyKind::Face
                    && faces
                        .iter()
                        .any(|face| face.index == e.key.ordinal && face.plane.is_none())
            })
            .unwrap()
            .key;
        let intent = GeometryIntent::capture(snapshot, key).unwrap();
        f.request.parameters.extent = Extent::OneSided {
            end: Termination::ToFace {
                reference: intent.id,
                offset: length(0.),
            },
        };
        f.request.faces = vec![FaceInput {
            source: target.feature,
            intent,
        }];
        f.commit();
        let original = f.session.document().clone();
        let error = f
            .compile_with_backend(backend)
            .err()
            .expect("Curved face must not be approximated by a tangent or center plane");
        assert_eq!(error.entity, Some(f.request.outputs[0].feature));
        assert!(error.message.contains("curved"));
        assert_eq!(f.session.document(), &original);
    }
}
#[test]
fn signed_face_offset_uses_native_outward_normal_with_opposing_inherited_frames() {
    for negative_first in [false, true] {
        let mut f = Fixture::new();
        let ids = [
            EntityId::parse("00000000-0000-4000-8000-000000000010").unwrap(),
            EntityId::parse("00000000-0000-4000-8000-000000000020").unwrap(),
        ];
        let mut positive = PrimitiveRequest::new(
            "Positive-frame box",
            PrimitiveParameters::Box {
                size: [length(10.), length(6.), length(10.)],
                frame: Default::default(),
            },
        );
        let mut negative = PrimitiveRequest::new(
            "Negative-frame same box",
            PrimitiveParameters::Box {
                size: [length(10.), length(6.), length(10.)],
                frame: PrimitiveFrame {
                    origin: [length(0.), length(6.), length(10.)],
                    u: [1., 0., 0.],
                    v: [0., -1., 0.],
                },
            },
        );
        positive.feature = ids[usize::from(negative_first)];
        negative.feature = ids[usize::from(!negative_first)];
        for request in [&positive, &negative] {
            let prepared =
                model::prepare_primitive(&f.session, f.session.revision(), actor(), request)
                    .unwrap();
            f.session.commit(prepared.prepared).unwrap();
        }
        let boolean = BooleanRequest::new(
            Operation::Join,
            BodyInput {
                body: positive.body,
                producer: positive.feature,
            },
            vec![BodyInput {
                body: negative.body,
                producer: negative.feature,
            }],
        );
        let prepared =
            model::prepare_boolean(&f.session, f.session.revision(), actor(), &boolean).unwrap();
        f.session.commit(prepared.prepared).unwrap();
        let initial = f.compile();
        let solid = initial.solid(boolean.feature).unwrap();
        let face = solid
            .faces()
            .unwrap()
            .into_iter()
            .find(|face| {
                face.plane
                    .is_some_and(|(p, n)| (p[2] - 10.).abs() < 1e-6 && n[2] > 0.99)
            })
            .expect("Native outward top cap");
        let topology = &initial.snapshot().shapes[&boolean.feature].topology;
        let element = topology
            .elements
            .iter()
            .find(|e| e.key.kind == TopologyKind::Face && e.key.ordinal == face.index)
            .unwrap();
        assert!(
            element
                .frames
                .iter()
                .any(|frame| frame.frame.normal()[2] > 0.99)
        );
        let opposite = element
            .frames
            .iter()
            .find(|frame| frame.frame.normal()[2] < -0.99)
            .expect("Coincident Boolean inherited the opposite generator frame");
        assert_eq!(opposite.origin.source, negative.feature);
        assert_eq!(
            element.frames[0].frame.normal()[2] < 0.,
            negative_first,
            "The fixture did not reverse the stored frame ordering"
        );
        let captured = GeometryIntent::capture(topology, element.key).unwrap();
        let original_sketch = f.session.document().get(f.request.sketch).unwrap().clone();
        let mut created = false;
        for anchored in [false, true] {
            for offset in [-2., 2.] {
                let mut intent = captured.clone();
                if anchored {
                    intent.attachment = Some(Attachment {
                        anchor: opposite.origin.clone(),
                        offset_mm: [0.; 3],
                        angle_rad: 0.,
                    });
                }
                let faces = [FaceInput {
                    source: boolean.feature,
                    intent,
                }];
                let parameters = ExtrudeParameters {
                    extent: Extent::OneSided {
                        end: Termination::ToFace {
                            reference: captured.id,
                            offset: length(offset),
                        },
                    },
                    ..Default::default()
                };
                if created {
                    let prepared = model::prepare_extrude_edit_with_faces(
                        &f.session,
                        f.session.revision(),
                        actor(),
                        f.request.outputs[0].feature,
                        &parameters,
                        &faces,
                    )
                    .unwrap();
                    f.session.commit(prepared.prepared).unwrap();
                } else {
                    f.request.parameters = parameters;
                    f.request.faces = faces.to_vec();
                    f.commit();
                    created = true;
                }
                let original = f.session.document().clone();
                let current = f.compile();
                let shape = f.shape(&current);
                near(shape.volume_mm3, 60. * (10. + offset));
                near(shape.bounds.max[2], 10. + offset);
                let end = cap(shape, "cap.end");
                for frame in &end.frames {
                    assert!(frame.frame.normal()[2] > 0.99);
                    near(frame.frame.origin_mm[2], 10. + offset);
                }
                assert_eq!(f.session.document(), &original);
                assert_eq!(original.get(f.request.sketch).unwrap(), &original_sketch);
            }
        }
    }
}
