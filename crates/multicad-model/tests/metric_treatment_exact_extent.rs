//! Exact receiving end and actual treatment source lifecycle; generic pair
//! intersection remains an independent unsupported gate in metric_fit_threads.
use multicad_core::{
    commands::{Actor, Channel, Session},
    constraint_sketch::Dimension,
    document::Document,
    model::{self, *},
    parameters::BoundExpression,
    project_file::{self, ProjectFile},
};
use multicad_kernel::{
    MeshOptions, RigidTransform, rust_kernel::authored_internal_treatment as source,
};
use multicad_model::Compiler;
fn actor() -> Actor {
    Actor::local(Channel::Ui)
}
fn compile(s: &Session) -> multicad_model::Model {
    Compiler {
        omit_display: true,
        ..Default::default()
    }
    .compile(s.document(), s.revision())
    .unwrap()
}
fn pair(
    d: f64,
    pitch: f64,
    hand: ThreadHand,
    treated: bool,
    wall_outer: Option<f64>,
) -> (Session, ThreadRequest, ThreadRequest) {
    pair_with_length(d, pitch, hand, treated, wall_outer, length(1.1 * d))
}
fn pair_with_length(
    d: f64,
    pitch: f64,
    hand: ThreadHand,
    treated: bool,
    wall_outer: Option<f64>,
    extent: multicad_core::constraint_sketch::Dimension,
) -> (Session, ThreadRequest, ThreadRequest) {
    let mut s = Session::new(Document::new("Metric fit metrology").unwrap()).unwrap();
    let mut p = ThreadParameters {
        major_diameter: length(d),
        pitch: length(pitch),
        length: extent,
        hand,
        fit: ThreadFit::Metric6gExternal,
        ..Default::default()
    };
    if treated {
        p.end_treatments = Some(ThreadEndTreatments {
            start: ThreadStartTreatment::Chamfer45 {
                radial_depth: length(0.8 * pitch),
            },
            end: ThreadEndTreatment::CylindricalRelief {
                length: length(0.6 * pitch),
                radial_depth: length(0.1 * pitch),
            },
        });
    }
    let male = ThreadRequest::new("6g shaft", p.clone());
    let q = model::prepare_thread(&s, s.revision(), actor(), &male).unwrap();
    s.commit(q.prepared).unwrap();
    let blank = PrimitiveRequest::new(
        "Receiver",
        PrimitiveParameters::Cylinder {
            radius: length(wall_outer.unwrap_or(d / 2. + 2.)),
            height: p.length.clone(),
            frame: Default::default(),
        },
    );
    let q = model::prepare_primitive(&s, s.revision(), actor(), &blank).unwrap();
    s.commit(q.prepared).unwrap();
    p.fit = ThreadFit::Metric6HInternal;
    let female = ThreadRequest::receiving(
        "6H internal",
        ThreadMode::CutInternal,
        BodyInput {
            body: blank.body,
            producer: blank.feature,
        },
        p,
    )
    .unwrap();
    let q = model::prepare_thread(&s, s.revision(), actor(), &female).unwrap();
    s.commit(q.prepared).unwrap();
    (s, male, female)
}

#[test]
fn canonical_metric_receiving_end_is_shared_by_both_hands_and_survives_clear_edit_and_reopen() {
    for hand in [ThreadHand::Right, ThreadHand::Left] {
        for expression in [false, true] {
            let dimension = if expression {
                Dimension::Expression {
                    expression: BoundExpression {
                        text: "=1.1*3".into(),
                        unit: "mm".into(),
                        bindings: vec![],
                    },
                }
            } else {
                length(1.1 * 3.)
            };
            let (mut session, male, female) =
                pair_with_length(3., 0.5, hand, true, None, dimension);
            let built = compile(&session);
            let a = built.solid(male.feature).unwrap();
            let b = built.solid(female.feature).unwrap();
            let canonical = multicad_core::parameters::exact_length_mm(
                &female.parameters.length,
                session.parameters(),
            )
            .unwrap();
            let data = b.rust_internal_thread_treatment_data().unwrap();
            assert_eq!(data.material.stock().height(), canonical);
            assert_eq!(
                data.material.cavity().exact_receiving_extent(),
                Some(canonical)
            );
            let boundary = data
                .material
                .boundary_with_interrupt(Default::default(), None)
                .unwrap();
            assert_eq!(boundary.face_loops(1, None).unwrap().len(), 2);
            let r = female
                .parameters
                .resolve(session.parameters(), ThreadMode::CutInternal)
                .unwrap();
            let ends = r.end_treatments.as_ref().unwrap();
            // Actual-source section witnesses are kept separate from the still-open
            // whole-pair Boolean/separation assertion; they cannot replace that gate.
            for (z, radius) in [
                (
                    0.01,
                    r.major_radius_mm - r.profile.depth_mm + ends.start_chamfer_radial_depth_mm
                        - 0.01,
                ),
                (
                    r.length_mm - 0.1,
                    r.major_radius_mm + ends.end_relief_radial_depth_mm,
                ),
            ] {
                for phi in [0f64, 0.43, 2.1] {
                    for (delta, void) in [(-0.001, true), (0.001, false)] {
                        let radial = radius + delta;
                        assert_eq!(
                            data.material
                                .cavity()
                                .contains_void([radial * phi.cos(), radial * phi.sin(), z], 0.),
                            void
                        );
                    }
                }
            }
            let male_resolved = male
                .parameters
                .resolve(session.parameters(), ThreadMode::NewExternal)
                .unwrap();
            let treatment = male_resolved.end_treatments.as_ref().unwrap();
            for (z, radius) in [
                (
                    0.01,
                    male_resolved.major_radius_mm - treatment.start_chamfer_radial_depth_mm + 0.01,
                ),
                (
                    male_resolved.length_mm - 0.1,
                    male_resolved.major_radius_mm
                        - male_resolved.profile.depth_mm
                        - treatment.end_relief_radial_depth_mm,
                ),
            ] {
                for phi in [0f64, 0.43, 2.1] {
                    for (delta, inside) in [(-0.001, true), (0.001, false)] {
                        let radius = radius + delta;
                        assert_eq!(
                            a.rust_thread_contains([radius * phi.cos(), radius * phi.sin(), z])
                                .unwrap(),
                            inside
                        );
                    }
                }
            }
            let cap = b
                .rust_internal_thread_treatment_world_face_properties_with_interrupt(
                    source::FaceKey::EndCap,
                    Default::default(),
                    None,
                )
                .unwrap()
                .unwrap();
            assert!(cap.world.area_mm2.lower > 0.);
            let moved = b
                .transformed(RigidTransform {
                    rotation: [[0.6, 0., 0.8], [0., 1., 0.], [-0.8, 0., 0.6]],
                    translation_mm: [11., -7., 5.],
                })
                .unwrap();
            assert_eq!(
                moved
                    .rust_internal_thread_treatment_data()
                    .unwrap()
                    .material,
                data.material
            );
            moved
                .rust_internal_thread_treatment_world_face_properties_with_interrupt(
                    source::FaceKey::EndCap,
                    Default::default(),
                    None,
                )
                .unwrap()
                .unwrap();

            let original = model::read_thread(session.document(), female.feature).unwrap();
            let mut clear = original.parameters.clone();
            clear.end_treatments = None;
            session
                .commit(
                    model::prepare_thread_edit(
                        &session,
                        session.revision(),
                        actor(),
                        female.feature,
                        &clear,
                    )
                    .unwrap()
                    .prepared,
                )
                .unwrap();
            let cleared = compile(&session);
            let plain = cleared.solid(female.feature).unwrap();
            assert!(plain.rust_internal_thread_treatment_data().is_none());
            assert!(
                plain.rust_thread_contains([1.6, 0., 0.01]).unwrap(),
                "clearing entry treatment restores actual wall"
            );
            session
                .commit(
                    model::prepare_thread_edit(
                        &session,
                        session.revision(),
                        actor(),
                        female.feature,
                        &original.parameters,
                    )
                    .unwrap()
                    .prepared,
                )
                .unwrap();
            let restored = compile(&session);
            assert_eq!(
                restored
                    .solid(female.feature)
                    .unwrap()
                    .rust_internal_thread_treatment_data()
                    .unwrap()
                    .material,
                data.material
            );
            let mut bytes = Vec::new();
            project_file::write(
                &ProjectFile::new(session.document().clone()).unwrap(),
                &mut bytes,
                Default::default(),
            )
            .unwrap();
            drop(built);
            drop(cleared);
            drop(restored);
            let reopened = project_file::read(
                &mut std::io::Cursor::new(bytes),
                Default::default(),
                &project_file::supported_capabilities(),
            )
            .unwrap()
            .project;
            let fresh = Session::new(reopened.document).unwrap();
            let rebuilt = compile(&fresh);
            assert_eq!(
                model::read_thread(fresh.document(), female.feature)
                    .unwrap()
                    .parameters,
                original.parameters
            );
            let solid = rebuilt.solid(female.feature).unwrap();
            assert_eq!(
                solid
                    .rust_internal_thread_treatment_data()
                    .unwrap()
                    .material,
                data.material
            );

            eprintln!(
                "METRIC_EXACT_RECEIVING hand={hand:?} expression={expression} end={:?}",
                canonical.ratio()
            );
        }
    }
}
#[test]
fn exact_shared_end_keeps_largest_wall_refusal() {
    let (session, _, _) = pair(3., 0.5, ThreadHand::Right, true, Some(2.55));
    assert!(
        Compiler {
            omit_display: true,
            ..Default::default()
        }
        .compile(session.document(), session.revision())
        .is_err()
    );
}

// This is a distinct real manufacturing gate, not a replacement for source-end
// checks or for the whole-pair intersection assertion in metric_fit_threads.
#[test]
fn actual_metric_receiving_source_requires_fresh_literal_manufacturing_embedding() {
    let mut failures = Vec::new();
    for hand in [ThreadHand::Right, ThreadHand::Left] {
        for expression in [false, true] {
            let dimension = if expression {
                Dimension::Expression {
                    expression: BoundExpression {
                        text: "=1.1*3".into(),
                        unit: "mm".into(),
                        bindings: vec![],
                    },
                }
            } else {
                length(1.1 * 3.)
            };
            let (session, _, female) = pair_with_length(3., 0.5, hand, true, None, dimension);
            let model = compile(&session);
            let body = model.solid(female.feature).unwrap();
            let moved = body
                .transformed(RigidTransform {
                    rotation: [[0.6, 0., 0.8], [0., 1., 0.], [-0.8, 0., 0.6]],
                    translation_mm: [11., -7., 5.],
                })
                .unwrap();
            for (placed, body) in [(false, body), (true, &moved)] {
                let result = body.manufacturing_mesh(MeshOptions {
                    linear_deflection_mm: 0.02,
                    angular_deflection_rad: 0.3,
                    maximum_triangles: 500_000,
                });
                if let Err(e) = result {
                    failures.push(format!(
                        "{hand:?} expression={expression} placed={placed}: {e:?}"
                    ));
                }
            }
        }
    }
    assert!(
        failures.is_empty(),
        "Actual metric manufacturing gates: {failures:#?}"
    );
}
