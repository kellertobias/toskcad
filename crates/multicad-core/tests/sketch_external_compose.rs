//! CAD-383 snap-time projection: several captures of one sketch, prepared
//! against the same base document, compose into one transaction.
use multicad_core::{
    constraint_sketch::source::Source,
    document::{Definition, Document, Entity, EntityData, EntityId, Mutation},
    references::{
        AnalyticGeometry, GeometryIntent, Origin, TopologyElement, TopologyKey, TopologyKind,
        TopologySnapshot,
    },
    sketch::Datum,
    sketch_external::{self, Behavior, CaptureInput, CaptureRequest, Method},
};
use std::collections::{BTreeMap, BTreeSet};

struct Fixture {
    document: Document,
    sketch: EntityId,
    owner: EntityId,
    map: BTreeMap<EntityId, TopologySnapshot>,
}
fn fixture() -> Fixture {
    let (sketch, edits) = Source::create(Datum::Xy, "Face sketch").unwrap();
    let mut document = Document::new("Compose captures")
        .unwrap()
        .edited(&edits)
        .unwrap();
    let source = Entity::new(
        "Analytic corners",
        EntityData::Feature {
            definition: Definition {
                type_id: "test.analytic".into(),
                schema_version: 1,
                data: serde_json::json!({}),
            },
            inputs: vec![],
            suppressed: false,
        },
    );
    let owner = source.id;
    document.insert(source).unwrap();
    let corner = |ordinal: u32, position_mm: [f64; 3]| TopologyElement {
        key: TopologyKey {
            kind: TopologyKind::Vertex,
            ordinal,
        },
        origins: BTreeSet::from([Origin {
            source: owner,
            role: format!("test.corner.{ordinal}"),
        }]),
        boundary: BTreeSet::new(),
        frames: vec![],
        geometry: Some(AnalyticGeometry::Point { position_mm }),
    };
    let snapshot = TopologySnapshot {
        owner,
        revision: 0,
        elements: vec![corner(0, [0., 0., 0.]), corner(1, [30., 20., 0.])].into(),
    };
    Fixture {
        document,
        sketch,
        owner,
        map: BTreeMap::from([(owner, snapshot)]),
    }
}
fn capture(f: &Fixture, ordinal: u32) -> (Document, CaptureRequest) {
    let intent = GeometryIntent::capture(
        &f.map[&f.owner],
        TopologyKey {
            kind: TopologyKind::Vertex,
            ordinal,
        },
    )
    .unwrap();
    let request = CaptureRequest::new(
        &f.document,
        0,
        CaptureInput {
            sketch: f.sketch,
            source: f.owner,
            intent,
            method: Method::Coplanar,
            behavior: Behavior::ReferenceOnly,
        },
        &f.map,
    )
    .unwrap();
    let preview = f
        .document
        .edited(&sketch_external::capture_edits(&f.document, 0, &request, &f.map).unwrap())
        .unwrap();
    (preview, request)
}
fn source(document: &Document, sketch: EntityId) -> Source {
    let Some(EntityData::Sketch { definition, .. }) = document.get(sketch).map(|e| &e.data) else {
        panic!("sketch")
    };
    Source::from_definition(definition).unwrap()
}

#[test]
fn two_corner_captures_compose_into_one_sketch_replacement_with_both_references() {
    let f = fixture();
    let (first, a) = capture(&f, 0);
    let (second, b) = capture(&f, 1);
    let edits =
        sketch_external::merge_capture_previews(&f.document, f.sketch, &[&first, &second]).unwrap();
    assert_eq!(edits.len(), 3, "two reference inserts and one sketch edit");
    assert!(matches!(&edits[2], Mutation::Replace(e) if e.id == f.sketch));
    let composed = f.document.edited(&edits).unwrap();
    let bindings = sketch_external::bindings(&composed, f.sketch).unwrap();
    let references: BTreeSet<_> = bindings.iter().map(|b| b.reference).collect();
    assert_eq!(references, BTreeSet::from([a.reference, b.reference]));
    let points = source(&composed, f.sketch).sketch.points;
    for (request, expected) in [(&a, [0., 0.]), (&b, [30., 20.])] {
        let id = request.output.points()[0];
        let point = points.iter().find(|p| p.id == id).unwrap();
        assert_eq!(point.position_mm, expected, "exact projected corner");
    }
    // One capture alone is exactly that capture.
    let single = sketch_external::merge_capture_previews(&f.document, f.sketch, &[&first]).unwrap();
    assert_eq!(f.document.edited(&single).unwrap(), first);
}

#[test]
fn captures_that_change_existing_source_or_repeat_an_identity_are_refused() {
    let f = fixture();
    let (first, _) = capture(&f, 0);
    assert!(
        sketch_external::merge_capture_previews(&f.document, f.sketch, &[&first, &first]).is_err(),
        "the same capture twice shares its identities"
    );
    assert!(sketch_external::merge_capture_previews(&f.document, f.sketch, &[]).is_err());
    // A preview that also renames an unrelated entity is not a pure capture.
    let mut renamed = f.document.get(f.owner).unwrap().clone();
    renamed.name = "Changed".into();
    let tampered = first.edited(&[Mutation::Replace(renamed)]).unwrap();
    let error = sketch_external::merge_capture_previews(&f.document, f.sketch, &[&tampered])
        .unwrap_err()
        .to_string();
    assert!(
        error.contains("changed existing project entities"),
        "{error}"
    );
    // A preview whose sketch edits an existing element is refused too.
    let base = source(&f.document, f.sketch);
    let mut grown = base.clone();
    grown
        .sketch
        .points
        .push(multicad_core::constraint_sketch::Point {
            id: EntityId::new(),
            position_mm: [1., 1.],
        });
    let with_point = f
        .document
        .edited(&grown.replace(&f.document, f.sketch).unwrap())
        .unwrap();
    let (capture_on_grown, _) = {
        let g = Fixture {
            document: with_point.clone(),
            sketch: f.sketch,
            owner: f.owner,
            map: f.map.clone(),
        };
        capture(&g, 1)
    };
    let mut moved = source(&capture_on_grown, f.sketch);
    moved.sketch.points[0].position_mm = [2., 2.];
    let edited = capture_on_grown
        .edited(&moved.replace(&capture_on_grown, f.sketch).unwrap())
        .unwrap();
    let error = sketch_external::merge_capture_previews(&with_point, f.sketch, &[&edited])
        .unwrap_err()
        .to_string();
    assert!(error.contains("edited existing sketch source"), "{error}");
}
