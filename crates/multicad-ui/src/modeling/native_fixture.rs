//! Tests use the production native model worker and detached geometry receipt.
use crate::viewport::{self, Scene, SceneBackend, SceneReply, SceneStamp};
use multicad_core::{
    document::Document,
    model::{EvaluationOptions, ModelEvaluationInfo},
    project_file::ProjectResources,
};
use multicad_model::{Compiler, ModelError, ModelSnapshot, Runtime};
use std::sync::Arc;

#[path = "../../../multicad-mcp/tests/support/screw_connection.rs"]
pub(crate) mod screw_connection;
/// Native aid inspections as scene review data.
fn aid_reviews(
    inspections: &Result<Vec<multicad_model::AidInspection>, String>,
) -> Result<Vec<viewport::AidReview>, String> {
    let list = inspections.as_ref().map_err(Clone::clone)?;
    Ok(list
        .iter()
        .map(|a| viewport::AidReview {
            feature: a.feature,
            part: a.part,
            target: a.target.part,
            touches_target: a.target.relation != multicad_kernel::SolidPairRelation::Separated,
            target_distance_mm: a.target.distance_mm,
            crowded: a.crowded.clone(),
            issues: a.issues.clone(),
            contact_issues: a.contact_issues.clone(),
            bounds_mm: a.bounds_mm,
            exact_contact: a.exact_contact.map(|c| viewport::ExactContact {
                point_mm: c.point_mm,
                angle_deg: c.angle_deg,
            }),
        })
        .collect())
}
pub(crate) fn scene(
    _document: &Document,
    model: &ModelSnapshot,
    generation: u64,
) -> Result<Arc<Scene>, ModelError> {
    let make = || -> Result<Scene, String> {
        let geometries = model
            .shapes
            .iter()
            .map(|(id, shape)| {
                viewport::Geometry::new(
                    *id,
                    shape.mesh.vertices.as_chunks::<3>().0.to_vec(),
                    shape.mesh.indices.clone(),
                    shape.mesh.face_indices.clone(),
                    shape.topology.elements.clone(),
                )?
                .with_topology_display(
                    shape
                        .edges
                        .iter()
                        .map(|edge| viewport::DisplayEdge {
                            index: edge.index,
                            points_mm: edge.points_mm.clone(),
                        })
                        .collect(),
                    shape
                        .vertices
                        .iter()
                        .map(|point| viewport::DisplayPoint {
                            index: point.index,
                            position_mm: point.position_mm,
                        })
                        .collect(),
                    0.1,
                )
                .map(Arc::new)
            })
            .collect::<Result<Vec<_>, String>>()?;
        let mut instances: Vec<_> = model
            .bodies
            .iter()
            .map(|(body, producer)| viewport::Instance {
                id: *body,
                geometry: *producer,
                transform: Default::default(),
                color: [128, 154, 204, 255],
            })
            .collect();
        for occurrence in model.occurrences.values() {
            instances.push(viewport::Instance {
                id: occurrence.instance,
                geometry: occurrence.producer,
                transform: viewport::Transform::from_rigid_pose(occurrence.pose)?,
                color: [128, 154, 204, 255],
            });
        }
        let mut native_bounds: std::collections::BTreeMap<_, _> = model
            .shapes
            .iter()
            .map(|(id, shape)| (*id, [shape.bounds.min, shape.bounds.max]))
            .collect();
        for (body, producer) in &model.body_sources {
            let bounds = model
                .shapes
                .get(producer)
                .ok_or("Native Part bounds have no current producer")?
                .bounds;
            native_bounds.insert(*body, [bounds.min, bounds.max]);
        }
        for occurrence in model.occurrences.values() {
            native_bounds.insert(
                occurrence.instance,
                [occurrence.bounds.min, occurrence.bounds.max],
            );
        }
        let scene = Scene::new(
            SceneStamp {
                document: model.stamp.document,
                revision: model.stamp.revision,
                generation,
            },
            geometries,
            instances,
        )?
        .with_reference_aliases(model.body_sources.clone())?;
        let scene = match &model.scope {
            multicad_model::EvaluationScope::Full => scene,
            multicad_model::EvaluationScope::Sources(roots) => {
                scene.with_source_scope(roots.clone())?
            }
        };
        let scene = scene
            .with_model_evaluation(Arc::new(ModelEvaluationInfo {
                screw_connection_checks: model.screw_connection_checks.clone(),
                source_hash: model.stamp.source_hash,
                consumed_bodies: model.consumed_bodies.clone(),
                cut_targets: model.cut_targets.clone(),
            }))?
            .with_native_reference_report(model.stamp.source_hash, model.references.clone())?
            .with_native_bounds(model.stamp.source_hash, native_bounds)?;
        if !scene.is_full_model() {
            return Ok(scene);
        }
        let scene = scene.with_assembly_evaluation(Arc::new(
            multicad_core::model::AssemblyEvaluationInfo {
                source_hash: model.stamp.source_hash,
                leaves: model
                    .leaf_occurrences
                    .values()
                    .map(|leaf| multicad_core::model::AssemblyLeafReceipt {
                        path: leaf.path.clone(),
                        component: leaf.component,
                        producer: leaf.producer,
                        pose: leaf.pose,
                    })
                    .collect(),
            },
        ))?;
        scene
            .with_transfer_evaluation(Arc::new(multicad_core::model::TransferEvaluationInfo {
                source_hash: model.stamp.source_hash,
                counts: model.transfer_counts.clone(),
            }))?
            .with_fdm_reviews(model.stamp.source_hash, aid_reviews(&model.fdm_aids))
    };
    make().map(Arc::new).map_err(|message| ModelError {
        entity: None,
        message,
    })
}
pub(crate) struct NativeBackend {
    runtime: Runtime<Arc<Scene>>,
    expected: Option<(u64, SceneStamp)>,
}
impl NativeBackend {
    pub fn new() -> Self {
        Self {
            runtime: Runtime::with_publisher(Compiler::default(), scene).unwrap(),
            expected: None,
        }
    }
}
impl SceneBackend for NativeBackend {
    fn request(&mut self, document: &Document, stamp: SceneStamp) -> Result<(), String> {
        self.request_evaluation(document, &Default::default(), stamp, &Default::default())
    }
    fn request_evaluation(
        &mut self,
        document: &Document,
        _resources: &ProjectResources,
        stamp: SceneStamp,
        options: &EvaluationOptions,
    ) -> Result<(), String> {
        let generation =
            self.runtime
                .request_evaluation(document, stamp.revision, stamp.generation, options)?;
        self.expected = Some((generation, stamp));
        Ok(())
    }
    fn request_sources(
        &mut self,
        document: &Document,
        _resources: &ProjectResources,
        stamp: SceneStamp,
        roots: &[multicad_core::document::EntityId],
    ) -> Result<(), String> {
        let generation =
            self.runtime
                .request_sources(document, stamp.revision, stamp.generation, roots)?;
        self.expected = Some((generation, stamp));
        Ok(())
    }
    fn cancel(&mut self) {
        self.expected = None;
        self.runtime.cancel();
    }
    fn poll(&mut self) -> Option<SceneReply> {
        let (expected, stamp) = self.expected?;
        let reply = self.runtime.poll()?;
        if reply.generation != expected {
            return None;
        }
        self.expected = None;
        Some(SceneReply {
            stamp,
            result: reply.result.map_err(|e| e.to_string()),
        })
    }
}

/// CAD-405/409/411: a 20 × 20 × 10 box with a face sketch circle around its
/// corner, the clipped segment extruded `height` as a second, new body.
/// Returns the session, the box and the segment body.
pub(crate) fn corner_segment_session(
    height: f64,
) -> (
    multicad_core::commands::Session,
    multicad_core::model::BodyInput,
    multicad_core::model::BodyInput,
) {
    use multicad_core::{
        commands::{Actor, Channel, Request, Session},
        constraint_sketch::{Curve, CurveKind, source::Source},
        document::EntityId,
        model,
        references::AnalyticGeometry,
        references::TopologyKind,
        sketch::{self, CreationPlane, PlaneOptions},
    };
    let actor = || Actor::local(Channel::Ui);
    let compile = |s: &Session| {
        Compiler::default()
            .compile(s.document(), s.revision())
            .unwrap()
    };
    let mut session = Session::new(Document::new("Corner segment").unwrap()).unwrap();
    let stock = model::PrimitiveRequest::new(
        "Box",
        model::PrimitiveParameters::Box {
            size: [model::length(20.), model::length(20.), model::length(10.)],
            frame: Default::default(),
        },
    );
    session
        .commit(
            model::prepare_primitive(&session, session.revision(), actor(), &stock)
                .unwrap()
                .prepared,
        )
        .unwrap();
    let compiled = compile(&session);
    let topology = &compiled.snapshot().shapes[&stock.feature].topology;
    let top = topology
        .elements
        .iter()
        .find(|e| {
            e.key.kind == TopologyKind::Face
                && matches!(e.geometry, Some(AnalyticGeometry::Plane { frame })
                    if frame.normal()[2] > 0.9 && (frame.origin_mm[2] - 10.).abs() < 1e-9)
        })
        .unwrap();
    let anchor = top.frames.first().unwrap().origin.clone();
    let created = sketch::prepare_new(
        &session,
        session.revision(),
        actor(),
        "Corner circle",
        CreationPlane::Face {
            snapshot: topology,
            key: top.key,
            anchor: Some(&anchor),
        },
        PlaneOptions::default(),
    )
    .unwrap();
    let (sketch, frame) = (created.sketch, created.frame);
    session.commit(created.prepared).unwrap();
    let report = compile(&session).snapshot().references.clone();
    let mut source = Source::read_for_evaluation(
        session.document(),
        sketch,
        session.revision(),
        Some(&report),
    )
    .unwrap()
    .0;
    let center = source.add_point(frame.coordinates([2., 3., 10.]));
    source.sketch.curves.push(Curve {
        id: EntityId::new(),
        construction: false,
        kind: CurveKind::Circle {
            center,
            radius_mm: 8.,
        },
    });
    session
        .apply(Request {
            expected_revision: session.revision(),
            actor: actor(),
            label: "Circle around the face corner".into(),
            edits: source.replace(session.document(), sketch).unwrap(),
        })
        .unwrap();
    let report = compile(&session).snapshot().references.clone();
    let queries = multicad_kernel::NativeRegionQueries::analytic(Default::default()).unwrap();
    let semantic = multicad_solver::evaluate_document_sketch(
        session.document(),
        session.revision(),
        sketch,
        Some(&report),
    )
    .unwrap();
    let inside = semantic.frame.coordinates([1., 1., 10.]);
    let region = semantic
        .regions_with_queries(&queries)
        .unwrap()
        .regions
        .into_iter()
        .find(|r| r.contains_with_queries(inside, &queries).unwrap())
        .unwrap();
    let request = model::ExtrudeRequest::new(
        sketch,
        vec![region.intent()],
        model::ExtrudeParameters {
            extent: model::Extent::OneSided {
                end: model::Termination::Distance {
                    distance: model::length(height),
                },
            },
            ..Default::default()
        },
    );
    let body = model::BodyInput {
        body: request.outputs[0].body,
        producer: request.outputs[0].feature,
    };
    session
        .commit(
            model::prepare_extrude(&session, session.revision(), actor(), &request)
                .unwrap()
                .prepared,
        )
        .unwrap();
    let stock = model::BodyInput {
        body: stock.body,
        producer: stock.feature,
    };
    (session, stock, body)
}

/// One closed sketch loop: corner points, and for each edge from point `i`
/// to point `i + 1` an optional arc center with its direction.
pub(crate) type SketchLoop = (Vec<[f64; 2]>, Vec<Option<([f64; 2], bool)>>);

/// CAD-421: a top-plane sketch of `loops` (outer first), whose region with
/// every inner loop as a hole is extruded `height` as a new body.
pub(crate) fn extruded_loops_session(
    name: &str,
    loops: &[SketchLoop],
    height: f64,
) -> (
    multicad_core::commands::Session,
    multicad_core::model::BodyInput,
) {
    use multicad_core::{
        commands::{Actor, Channel, Session},
        constraint_sketch::{Curve, CurveKind, source::Source},
        document::EntityId,
        model,
        sketch::Datum,
        sketch_regions,
    };
    let (sketch, edits) = Source::create(Datum::Xy, name).unwrap();
    let document = Document::new(name).unwrap().edited(&edits).unwrap();
    let mut profile = Source::read(&document, sketch).unwrap().0;
    for (points, arcs) in loops {
        let ids: Vec<_> = points.iter().map(|&p| profile.add_point(p)).collect();
        for (i, arc) in arcs.iter().enumerate() {
            let (a, b) = (ids[i], ids[(i + 1) % ids.len()]);
            match *arc {
                None => {
                    profile.add_line(a, b);
                }
                Some((center, clockwise)) => {
                    let center = profile.add_point(center);
                    profile.sketch.curves.push(Curve {
                        id: EntityId::new(),
                        construction: false,
                        kind: CurveKind::Arc {
                            center,
                            start: a,
                            end: b,
                            clockwise,
                        },
                    });
                }
            }
        }
    }
    let region = sketch_regions::extract(&profile)
        .regions
        .into_iter()
        .find(|r| r.holes.len() == loops.len() - 1)
        .unwrap()
        .intent();
    let document = document
        .edited(&profile.replace(&document, sketch).unwrap())
        .unwrap();
    let mut session = Session::new(document).unwrap();
    let request = model::ExtrudeRequest::new(
        sketch,
        vec![region],
        model::ExtrudeParameters {
            extent: model::Extent::OneSided {
                end: model::Termination::Distance {
                    distance: model::length(height),
                },
            },
            ..Default::default()
        },
    );
    let body = model::BodyInput {
        body: request.outputs[0].body,
        producer: request.outputs[0].feature,
    };
    session
        .commit(
            model::prepare_extrude(
                &session,
                session.revision(),
                Actor::local(Channel::Ui),
                &request,
            )
            .unwrap()
            .prepared,
        )
        .unwrap();
    (session, body)
}

/// CAD-421: the 20 × 20 L with 8 mm arms whose inner corner is rounded
/// with radius 3; its inner chain runs from (20, 8) to (11, 8), around
/// (11, 11) and from (8, 11) to (8, 20).
pub(crate) fn rounded_l_loops() -> Vec<SketchLoop> {
    vec![(
        vec![
            [0., 0.],
            [20., 0.],
            [20., 8.],
            [11., 8.],
            [8., 11.],
            [8., 20.],
            [0., 20.],
        ],
        vec![None, None, None, Some(([11., 11.], true)), None, None, None],
    )]
}

/// CAD-421: a 30 × 20 plate with a pocket from (10, 6) to (20, 14) whose
/// corners are rounded with radius 2.
pub(crate) fn pocket_plate_loops() -> Vec<SketchLoop> {
    let outer = (
        vec![[0., 0.], [30., 0.], [30., 20.], [0., 20.]],
        vec![None; 4],
    );
    let pocket = (
        vec![
            [12., 6.],
            [18., 6.],
            [20., 8.],
            [20., 12.],
            [18., 14.],
            [12., 14.],
            [10., 12.],
            [10., 8.],
        ],
        vec![
            None,
            Some(([18., 8.], false)),
            None,
            Some(([18., 12.], false)),
            None,
            Some(([12., 12.], false)),
            None,
            Some(([12., 8.], false)),
        ],
    );
    vec![outer, pocket]
}
