//! CAD-382: face sketches face outward from the face's real material side,
//! while attachments saved before the `outward` flag keep their meaning.
use multicad_core::{
    commands::*,
    document::*,
    references::*,
    sketch::{self, CreationPlane, PlaneOptions},
};
use std::collections::BTreeSet;

fn near(a: [f64; 3], b: [f64; 3]) -> bool {
    (0..3).all(|i| (a[i] - b[i]).abs() < 1e-12)
}

/// A prism start cap: its semantic frame is the profile frame (normal +Z,
/// into the part) while the B-rep face normal leaves the material (-Z).
fn bottom(owner: EntityId, revision: u64, geometry: bool) -> TopologySnapshot {
    let start = Origin {
        source: owner,
        role: "cap.start".into(),
    };
    TopologySnapshot {
        owner,
        revision,
        elements: vec![TopologyElement {
            key: TopologyKey {
                kind: TopologyKind::Face,
                ordinal: 0,
            },
            origins: BTreeSet::from([start.clone()]),
            boundary: BTreeSet::new(),
            frames: vec![AnchoredFrame {
                origin: start,
                frame: ReferenceFrame {
                    origin_mm: [0.; 3],
                    u: [1., 0., 0.],
                    v: [0., 1., 0.],
                },
            }],
            geometry: geometry.then_some(AnalyticGeometry::Plane {
                frame: ReferenceFrame {
                    origin_mm: [0.; 3],
                    u: [1., 0., 0.],
                    v: [0., -1., 0.],
                },
            }),
        }]
        .into(),
    }
}

fn resolved(intent: &GeometryIntent, snapshot: &TopologySnapshot) -> ReferenceFrame {
    match resolve(intent, snapshot, snapshot.revision) {
        Resolution::Resolved {
            frame: Some(frame), ..
        } => frame,
        other => panic!("{other:?}"),
    }
}

#[test]
fn saved_attachments_without_the_outward_flag_keep_the_source_frame_and_round_trip() {
    let owner = EntityId::new();
    let snapshot = bottom(owner, 0, true);
    // A face reference as written before CAD-382: no `outward` key.
    let saved = serde_json::json!({
        "id": EntityId::new(),
        "kind": "Face",
        "origins": [{"source": owner, "role": "cap.start"}],
        "boundary": null,
        "attachment": {
            "anchor": {"source": owner, "role": "cap.start"},
            "offset_mm": [0., 0., 2.],
            "angle_rad": 0.
        }
    });
    let intent: GeometryIntent = serde_json::from_value(saved.clone()).unwrap();
    assert!(!intent.attachment.as_ref().unwrap().outward);
    let frame = resolved(&intent, &snapshot);
    // Unchanged meaning: inward normal, and the offset follows it.
    assert!(near(frame.normal(), [0., 0., 1.]));
    assert!(near(frame.origin_mm, [0., 0., 2.]));
    assert!(near(frame.u, [1., 0., 0.]));
    // Saving again writes the same document; no new key appears.
    assert_eq!(serde_json::to_value(&intent).unwrap(), saved);
}

#[test]
fn outward_attachments_reverse_only_inward_frames_before_offset_and_angle() {
    let owner = EntityId::new();
    let snapshot = bottom(owner, 0, true);
    let mut intent = GeometryIntent::capture(&snapshot, snapshot.elements[0].key).unwrap();
    intent.attachment = Some(Attachment {
        anchor: snapshot.elements[0].frames[0].origin.clone(),
        offset_mm: [0., 0., 2.],
        angle_rad: std::f64::consts::FRAC_PI_2,
        outward: true,
    });
    let frame = resolved(&intent, &snapshot);
    assert!(near(frame.normal(), [0., 0., -1.]));
    // A positive normal offset now leaves the material.
    assert!(near(frame.origin_mm, [0., 0., -2.]));
    // The angle turns about the outward normal: u goes from +X toward -Y.
    assert!(near(frame.u, [0., -1., 0.]));
    let value = serde_json::to_value(&intent).unwrap();
    assert_eq!(value["attachment"]["outward"], serde_json::json!(true));
    let reopened: GeometryIntent = serde_json::from_value(value).unwrap();
    assert_eq!(reopened, intent);

    // A frame that already faces outward is left alone.
    let mut top = snapshot.clone();
    top.elements[0].geometry = Some(AnalyticGeometry::Plane {
        frame: top.elements[0].frames[0].frame,
    });
    assert!(near(resolved(&intent, &top).normal(), [0., 0., 1.]));
    // Without real face geometry nothing is guessed.
    let unknown = bottom(owner, 0, false);
    assert!(near(resolved(&intent, &unknown).normal(), [0., 0., 1.]));
}

#[test]
fn new_face_sketches_save_an_outward_attachment_on_a_start_cap() {
    let feature = Entity::new(
        "Panel",
        EntityData::Feature {
            definition: Definition {
                type_id: "test.feature".into(),
                schema_version: 1,
                data: serde_json::json!({}),
            },
            inputs: vec![],
            suppressed: false,
        },
    );
    let owner = feature.id;
    let mut document = Document::new("Outward").unwrap();
    document.insert(feature).unwrap();
    let session = Session::new(document).unwrap();
    let snapshot = bottom(owner, session.revision(), true);
    let created = sketch::prepare_new(
        &session,
        session.revision(),
        Actor::local(Channel::Ui),
        "Bottom",
        CreationPlane::Face {
            snapshot: &snapshot,
            key: snapshot.elements[0].key,
            anchor: None,
        },
        PlaneOptions {
            normal_offset_mm: 1.,
            angle_rad: 0.,
        },
    )
    .unwrap();
    assert!(near(created.frame.n, [0., 0., -1.]));
    assert!(near(created.frame.origin, [0., 0., -1.]));
    let mut session = session;
    session.commit(created.prepared).unwrap();
    let Some(EntityData::Plane {
        normal,
        anchor: Some(link),
        ..
    }) = session.document().get(created.plane).map(|e| &e.data)
    else {
        panic!("face plane")
    };
    assert!(near(*normal, [0., 0., -1.]));
    let attachment = link.geometry.as_ref().unwrap().attachment.as_ref().unwrap();
    assert!(attachment.outward);
}
