//! Fillets and chamfers as first-class shapes (CAD-414): recorded by the
//! corner treatment, linked to the shape they modify, sized by their
//! canonical dimension, labelled with that size, pruned, validated, saved
//! and reopened under their own capability.
use multicad_core::{
    constraint_sketch::{
        ConstraintKind, Dimension, Shape, ShapeKind, Sketch,
        binding::DimensionInput,
        shape_edit::{self, ShapeSize},
        source::{self, Source},
    },
    document::Document,
    project_file::{self, ProjectFile},
    sketch::Datum,
    sketch_corner::{CornerTreatment, corner_draft},
    sketch_draw::{self, DrawInput, DrawRequest, PointInput},
};

fn rectangle() -> Source {
    sketch_draw::draft(
        &Source {
            u: [1., 0., 0.],
            sketch: Sketch::default(),
        },
        &DrawRequest {
            input: DrawInput::Rectangle {
                a: PointInput::Position([0., 0.]),
                b: PointInput::Position([100., 60.]),
            },
            construction: false,
        },
    )
    .unwrap()
    .source
}
fn corner(source: &Source, at: [f64; 2]) -> multicad_core::document::EntityId {
    source
        .sketch
        .points
        .iter()
        .find(|p| p.position_mm == at)
        .unwrap()
        .id
}
/// A rectangle with a 3 mm fillet at (100, 60) and a 2 mm chamfer at (0, 0).
fn treated() -> (Source, Shape, Shape, Shape) {
    let source = rectangle();
    let rectangle = source.sketch.shapes[0].clone();
    let filleted = corner_draft(
        &source,
        corner(&source, [100., 60.]),
        CornerTreatment::Fillet { radius_mm: 3. },
    )
    .unwrap();
    let fillet = filleted.source.sketch.shapes[1].clone();
    assert_eq!(fillet.curves, vec![filleted.created_curves[0]]);
    let chamfered = corner_draft(
        &filleted.source,
        corner(&filleted.source, [0., 0.]),
        CornerTreatment::Chamfer { distance_mm: 2. },
    )
    .unwrap();
    let chamfer = chamfered.source.sketch.shapes[2].clone();
    assert_eq!(chamfer.curves, vec![chamfered.created_curves[0]]);
    (chamfered.source, rectangle, fillet, chamfer)
}
fn literal(value: f64) -> DimensionInput {
    DimensionInput::Literal {
        value,
        unit: "mm".into(),
    }
}
fn values(source: &Source, constraints: &[multicad_core::document::EntityId]) -> Vec<f64> {
    constraints
        .iter()
        .map(|id| {
            let c = source
                .sketch
                .constraints
                .iter()
                .find(|c| c.id == *id)
                .unwrap();
            match c.kind.dimension() {
                Some(Dimension::Literal { value, .. }) => *value,
                other => panic!("{other:?}"),
            }
        })
        .collect()
}

#[test]
fn corner_treatments_record_typed_shapes_linked_to_the_modified_shape() {
    let (source, rectangle, fillet, chamfer) = treated();
    assert_eq!(source.sketch.shapes.len(), 3);
    assert_eq!(fillet.kind, ShapeKind::Fillet);
    assert_eq!(fillet.name, "Fillet 1");
    assert_eq!(fillet.modifies, Some(rectangle.id));
    assert_eq!(chamfer.kind, ShapeKind::Chamfer);
    assert_eq!(chamfer.name, "Chamfer 1");
    assert_eq!(chamfer.modifies, Some(rectangle.id));
    assert_eq!(rectangle.modifies, None);
    // The rectangle keeps its four (trimmed) edges; helpers belong to none.
    assert_eq!(source.sketch.shapes[0].curves, rectangle.curves);
    let helpers = source.sketch.curves.iter().filter(|c| c.construction);
    for helper in helpers {
        assert!(source.sketch.shape_of(helper.id).is_none());
    }
    assert_eq!(ShapeKind::Fillet.label(), "Fillet");
    assert_eq!(ShapeKind::Chamfer.label(), "Chamfer");
}

#[test]
fn sizes_and_labels_come_from_the_canonical_dimensions() {
    let (source, _, fillet, chamfer) = treated();
    assert_eq!(shape_edit::sizes(&fillet), &[ShapeSize::Radius]);
    assert_eq!(shape_edit::sizes(&chamfer), &[ShapeSize::Distance]);
    let radius = shape_edit::corner_constraints(&source.sketch, &fillet);
    assert_eq!(radius.len(), 1);
    assert!(matches!(
        source.sketch.constraints.iter().find(|c| c.id == radius[0]).unwrap().kind,
        ConstraintKind::Radius { curve, .. } if curve == fillet.curves[0]
    ));
    // Both chamfer sides carry the setback.
    let setbacks = shape_edit::corner_constraints(&source.sketch, &chamfer);
    assert_eq!(values(&source, &setbacks), vec![2., 2.]);
    assert_eq!(
        shape_edit::corner_size(&source.sketch, &fillet),
        Some(Dimension::Literal {
            value: 3.,
            unit: "mm".into()
        })
    );
    assert_eq!(
        shape_edit::corner_label(&source.sketch, &fillet),
        "Fillet R3"
    );
    assert_eq!(
        shape_edit::corner_label(&source.sketch, &chamfer),
        "Chamfer 2"
    );
    let rectangle = &source.sketch.shapes[0];
    assert_eq!(
        shape_edit::corner_label(&source.sketch, rectangle),
        "Rectangle 1"
    );
}

#[test]
fn set_size_edits_the_existing_radius_and_every_chamfer_setback() {
    let (source, _, fillet, chamfer) = treated();
    let document = Document::new("Corner sizes").unwrap();
    let before = source.sketch.constraints.len();
    let resized = shape_edit::set_size(
        &document,
        &source,
        fillet.id,
        ShapeSize::Radius,
        &literal(10.),
    )
    .unwrap();
    assert_eq!(resized.sketch.constraints.len(), before, "edited in place");
    let radius = shape_edit::corner_constraints(&resized.sketch, &fillet);
    assert_eq!(values(&resized, &radius), vec![10.]);
    assert_eq!(
        shape_edit::corner_label(&resized.sketch, &fillet),
        "Fillet R10"
    );
    let resized = shape_edit::set_size(
        &document,
        &resized,
        chamfer.id,
        ShapeSize::Distance,
        &literal(4.5),
    )
    .unwrap();
    let setbacks = shape_edit::corner_constraints(&resized.sketch, &chamfer);
    assert_eq!(values(&resized, &setbacks), vec![4.5, 4.5]);
    assert_eq!(
        shape_edit::corner_label(&resized.sketch, &chamfer),
        "Chamfer 4.5"
    );
    // A size a fillet does not have is refused.
    assert!(
        shape_edit::set_size(
            &document,
            &resized,
            fillet.id,
            ShapeSize::Width,
            &literal(1.)
        )
        .is_err()
    );
}

#[test]
fn pruning_unlinks_a_treatment_whose_shape_is_gone_and_validation_guards_the_link() {
    let (source, rectangle, fillet, _) = treated();
    let mut pruned = source.clone();
    pruned
        .sketch
        .shapes
        .retain(|shape| shape.id != rectangle.id);
    pruned.sketch.prune_shapes();
    assert!(pruned.sketch.shapes.iter().all(|s| s.modifies.is_none()));
    // Deleting the fillet arc removes the fillet shape on write.
    let mut deleted = source.clone();
    deleted.sketch.curves.retain(|c| c.id != fillet.curves[0]);
    deleted.sketch.prune_shapes();
    assert!(deleted.sketch.shapes.iter().all(|s| s.id != fillet.id));
    // Only a fillet or chamfer links to another shape, never to itself.
    let mut wrong = source.clone();
    wrong.sketch.shapes[0].modifies = Some(fillet.id);
    assert!(wrong.validate().is_err());
    let mut own = source.clone();
    own.sketch.shapes[1].modifies = Some(fillet.id);
    assert!(own.validate().is_err());
    let mut two = source.clone();
    two.sketch.shapes[1].curves.push(two.sketch.curves[0].id);
    two.sketch.shapes[0]
        .curves
        .retain(|c| *c != two.sketch.curves[0].id);
    assert!(two.validate().is_err());
}

#[test]
fn corner_shapes_survive_save_and_reopen_and_require_their_capability() {
    let (id, edits) = Source::create(Datum::Xy, "Sketch 1").unwrap();
    let document = Document::new("Corner shapes")
        .unwrap()
        .edited(&edits)
        .unwrap();
    let plain = document
        .edited(&rectangle().replace(&document, id).unwrap())
        .unwrap();
    assert!(source::requires_shapes(&plain));
    assert!(
        !source::requires_corner_shapes(&plain),
        "old files unchanged"
    );
    let (treated, rectangle, fillet, chamfer) = treated();
    let document = document
        .edited(&treated.replace(&document, id).unwrap())
        .unwrap();
    assert!(source::requires_corner_shapes(&document));
    let file = ProjectFile::new(document.clone()).unwrap();
    assert!(
        file.capabilities()
            .contains(source::CORNER_SHAPES_CAPABILITY)
    );
    let mut bytes = vec![];
    project_file::write(&file, &mut bytes, Default::default()).unwrap();
    let opened = project_file::read(
        &mut std::io::Cursor::new(bytes),
        Default::default(),
        &project_file::supported_capabilities(),
    )
    .unwrap()
    .project
    .document;
    assert_eq!(opened, document);
    let (reopened, _, _) = Source::read(&opened, id).unwrap();
    assert_eq!(
        reopened.sketch.shapes,
        vec![rectangle, fillet, chamfer],
        "kinds and links reopen unchanged"
    );
    // A shape without the field still reads, with no link.
    let json = serde_json::json!({
        "id": multicad_core::document::EntityId::new(),
        "kind": "line",
        "name": "Line 1",
        "curves": [multicad_core::document::EntityId::new()],
    });
    let shape: Shape = serde_json::from_value(json).unwrap();
    assert_eq!(shape.modifies, None);
    assert!(
        serde_json::to_value(&shape)
            .unwrap()
            .get("modifies")
            .is_none()
    );
}
