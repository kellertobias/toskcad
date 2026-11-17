//! Normal canonical single-mode receiving cuts; no injected Native bodies.
use multicad_core::{
    commands::{Actor, Channel, Session},
    document::Document,
    model::{self, *},
};
use multicad_model::Compiler;
pub fn source(entry: bool, left: bool, height: f64) -> (Session, ThreadRequest) {
    source_dimension(entry, left, length(height))
}
pub fn source_dimension(
    entry: bool,
    left: bool,
    height: multicad_core::constraint_sketch::Dimension,
) -> (Session, ThreadRequest) {
    source_dimension_with_fit(entry, left, height, ThreadFit::Metric6HInternal)
}
pub fn source_dimension_with_fit(
    entry: bool,
    left: bool,
    height: multicad_core::constraint_sketch::Dimension,
    fit: ThreadFit,
) -> (Session, ThreadRequest) {
    let actor = Actor::local(Channel::Ui);
    let mut session =
        Session::new(Document::new("Single metric receiving treatment").unwrap()).unwrap();
    let stock = PrimitiveRequest::new(
        "Actual cylinder stock",
        PrimitiveParameters::Cylinder {
            radius: length(3.5),
            height: height.clone(),
            frame: Default::default(),
        },
    );
    let p = model::prepare_primitive(&session, session.revision(), actor.clone(), &stock).unwrap();
    session.commit(p.prepared).unwrap();
    let thread = ThreadRequest::receiving(
        "M3 single receiving treatment",
        ThreadMode::CutInternal,
        BodyInput {
            body: stock.body,
            producer: stock.feature,
        },
        ThreadParameters {
            major_diameter: length(3.),
            pitch: length(0.5),
            length: height,
            fit,
            hand: if left {
                ThreadHand::Left
            } else {
                ThreadHand::Right
            },
            end_treatments: Some(ThreadEndTreatments {
                start: if entry {
                    ThreadStartTreatment::Chamfer45 {
                        radial_depth: length(0.4),
                    }
                } else {
                    ThreadStartTreatment::FlatCut
                },
                end: if entry {
                    ThreadEndTreatment::FlatCut
                } else {
                    ThreadEndTreatment::CylindricalRelief {
                        length: length(0.3),
                        radial_depth: length(0.05),
                    }
                },
            }),
            ..Default::default()
        },
    )
    .unwrap();
    let p = model::prepare_thread(&session, session.revision(), actor, &thread).unwrap();
    session.commit(p.prepared).unwrap();
    (session, thread)
}
pub fn require_positive(entry: bool) {
    let mut failures = Vec::new();
    for left in [false, true] {
        for height in [4., 3.3000000000000003] {
            let (session, thread) = source(entry, left, height);
            match (Compiler {
                omit_display: true,
                ..Default::default()
            })
            .compile(session.document(), session.revision())
            {
                Ok(model) => {
                    assert_eq!(
                        model.solid(thread.feature).unwrap().solid_count().unwrap(),
                        1
                    );
                }
                Err(e) => failures.push(format!(
                    "entry={entry} left={left} height={height:?}: {e:?}"
                )),
            }
        }
    }
    assert!(
        failures.is_empty(),
        "Normal single receiving modes: {failures:#?}"
    );
}
