//! CAD-434 profile stock picker: a section picture for every type, a
//! dimensioned diagram whose callouts focus and highlight their fields,
//! plain names beside catalogue names, the bar length after the section,
//! and a validated tube corner radius that the member geometry uses.
use super::*;
use crate::materials::profile_picker as picker;
use crate::modeling::diagram::harness::{self, WIDTHS};
use multicad_core::{document::Document, project_file::ProjectResources};
use std::time::{Duration, Instant};

#[allow(dead_code)]
#[path = "../../../../multicad-app/examples/shared_materials_review.rs"]
mod shared_fixture;

const PI: f64 = std::f64::consts::PI;

fn profile(shape: ProfileShape, w: f64, h: f64, wall: Option<f64>, r: Option<f64>) -> StockKind {
    StockKind::Profile {
        designation: "Catalogue 0.0.1".into(),
        length_mm: 6000.,
        width_mm: w,
        height_mm: h,
        wall_mm: wall,
        shape,
        section: None,
        corner_radius_mm: r,
    }
}

/// One full egui pass with accessibility, returning the output.
fn run(
    ctx: &egui::Context,
    size: egui::Vec2,
    events: Vec<egui::Event>,
    show: &mut dyn FnMut(&mut egui::Ui),
) -> egui::FullOutput {
    let mut output = ctx.run_ui(
        egui::RawInput {
            screen_rect: Some(egui::Rect::from_min_size(egui::Pos2::ZERO, size)),
            events,
            focused: true,
            ..Default::default()
        },
        |ui| show(ui),
    );
    // Shapes are inspected, never painted: the texture upload is dropped.
    output.textures_delta.clear();
    output
}
fn node_center(
    output: &egui::FullOutput,
    matches: impl Fn(&egui::accesskit::Node) -> bool,
) -> Option<egui::Pos2> {
    output
        .platform_output
        .accesskit_update
        .as_ref()?
        .nodes
        .iter()
        .find_map(|(_, node)| matches(node).then(|| node.bounds()).flatten())
        .map(|b| egui::pos2(((b.x0 + b.x1) / 2.) as f32, ((b.y0 + b.y1) / 2.) as f32))
}
fn texts(output: &egui::FullOutput) -> Vec<(String, egui::Rect)> {
    output
        .shapes
        .iter()
        .filter_map(|s| match &s.shape {
            egui::Shape::Text(t) => Some((
                t.galley.text().to_owned(),
                egui::Rect::from_min_size(t.pos, t.galley.size()),
            )),
            _ => None,
        })
        .collect()
}
fn y_of(texts: &[(String, egui::Rect)], prefix: &str) -> f32 {
    texts
        .iter()
        .find(|(t, _)| t.starts_with(prefix))
        .unwrap_or_else(|| {
            panic!(
                "{prefix:?} painted; texts: {:?}",
                texts.iter().map(|t| &t.0).collect::<Vec<_>>()
            )
        })
        .1
        .top()
}
fn stock_form<'a>(spec: &'a mut StockKind) -> impl FnMut(&mut egui::Ui) + 'a {
    move |ui| {
        let mut request = None;
        stock_controls(ui, spec, false, &mut request);
    }
}

#[test]
fn every_profile_type_has_a_section_picture_from_its_canonical_section() {
    // Tube and bar variants of rectangular and round, angle, slot profile
    // and custom sketch each draw a picture; hollow types show their bore.
    for (shape, hollow, loops) in [
        (ProfileShape::Rectangular, true, 2),
        (ProfileShape::Rectangular, false, 1),
        (ProfileShape::Round, true, 2),
        (ProfileShape::Round, false, 1),
        (ProfileShape::Angle, true, 1),
        (ProfileShape::ItemLine8_40x40, false, 2),
    ] {
        let section = picker::sample(shape, hollow);
        assert_eq!(picker::loops(&section).len(), loops, "{shape:?} {hollow}");
        assert_eq!(
            picker::loops(&section).len(),
            1 + section.holes().len(),
            "the picture is the canonical outside plus bores"
        );
        assert!(!picker::family_art(shape, hollow).marks.is_empty());
    }
    assert!(
        !picker::family_art(ProfileShape::Sketch, true)
            .marks
            .is_empty()
    );
    // The rounded sample is drawn with its arcs, not as a sharp rectangle.
    let rounded = picker::loops(&picker::sample(ProfileShape::Rectangular, true));
    assert!(rounded[0].len() > 8, "arcs are sampled into the outline");

    // Open the picker in the stock form: every family is listed with its
    // picture cells (two for the families with tube and bar variants).
    let ctx = egui::Context::default();
    ctx.enable_accesskit();
    let mut spec = profile(ProfileShape::Rectangular, 40., 40., Some(2.), None);
    let size = egui::vec2(420., 1200.);
    let output = run(&ctx, size, vec![], &mut stock_form(&mut spec));
    let current_pictures = output
        .shapes
        .iter()
        .filter(
            |s| matches!(&s.shape, egui::Shape::Rect(r) if r.rect.size() == egui::vec2(40., 30.)),
        )
        .count();
    assert_eq!(current_pictures, 1, "the chosen type shows its picture");
    let combo = node_center(&output, |n| n.value() == Some("Rectangular")).expect("type picker");
    for events in harness::click(combo).chunks(2) {
        run(&ctx, size, events.to_vec(), &mut stock_form(&mut spec));
    }
    let output = run(&ctx, size, vec![], &mut stock_form(&mut spec));
    let painted = texts(&output);
    for (family, label, hint) in picker::FAMILIES {
        assert!(painted.iter().any(|(t, _)| t == label), "{family:?} listed");
        assert!(
            painted.iter().any(|(t, _)| t == hint),
            "{family:?} explained"
        );
    }
    let cells = output
        .shapes
        .iter()
        .filter(
            |s| matches!(&s.shape, egui::Shape::Rect(r) if r.rect.size() == egui::vec2(34., 26.)),
        )
        .count();
    assert_eq!(
        cells,
        picker::FAMILIES.len() + 2,
        "tube and bar pictures for rectangular and round"
    );
}

#[test]
fn every_editable_field_is_a_linked_callout_for_each_type() {
    for spec in [
        profile(ProfileShape::Rectangular, 40., 40., Some(2.), Some(4.)),
        profile(ProfileShape::Rectangular, 40., 8., None, None),
        profile(ProfileShape::Round, 33.7, 33.7, Some(2.6), None),
        profile(ProfileShape::Round, 20., 20., None, None),
        profile(ProfileShape::Angle, 40., 30., Some(4.), None),
    ] {
        let diagram = picker::diagram(&spec);
        let mut drawn = diagram.keys();
        let mut fields = picker::keys(&spec);
        drawn.sort_unstable();
        fields.sort_unstable();
        assert_eq!(drawn, fields, "{spec:?}");
        assert!(!fields.is_empty());
        for width in WIDTHS {
            harness::assert_callouts_apart(&diagram, width, 200., &|k| {
                picker::callout_value(&spec, k)
            });
        }
        // Each field is painted with its label in the form.
        let ctx = egui::Context::default();
        let mut staged = spec.clone();
        let painted = harness::frame(&ctx, 420., vec![], &mut stock_form(&mut staged));
        for label in match fields.len() {
            _ if fields.contains(&picker::DIAMETER) => vec!["Outside diameter"],
            _ => vec!["Width", "Height"],
        } {
            assert!(painted.has(label), "{label} field for {spec:?}");
        }
    }
    // The angle calls its wall the leg thickness on the drawing too.
    let angle = picker::diagram(&profile(ProfileShape::Angle, 40., 30., Some(4.), None));
    assert!(angle.callouts.iter().any(|c| c.name == "Leg thickness"));
}

#[test]
fn clicking_a_callout_focuses_highlights_and_edits_its_field_at_both_widths() {
    for width in WIDTHS {
        for (callout, typed, check) in [
            (
                "R Corner radius sharp",
                "4",
                (|s: &StockKind| matches!(s, StockKind::Profile { corner_radius_mm: Some(r), .. } if *r == 4.))
                    as fn(&StockKind) -> bool,
            ),
            (
                "Width 40 mm",
                "60",
                |s| matches!(s, StockKind::Profile { width_mm, height_mm, .. } if *width_mm == 60. && *height_mm == 40.),
            ),
            (
                "Height 40 mm",
                "30",
                |s| matches!(s, StockKind::Profile { width_mm, height_mm, .. } if *width_mm == 40. && *height_mm == 30.),
            ),
            (
                "Wall 2 mm",
                "3",
                |s| matches!(s, StockKind::Profile { wall_mm: Some(t), .. } if *t == 3.),
            ),
        ] {
            let ctx = egui::Context::default();
            let mut spec = profile(ProfileShape::Rectangular, 40., 40., Some(2.), None);
            harness::edit_through_callout(&ctx, width, callout, typed, &mut stock_form(&mut spec));
            assert!(check(&spec), "{callout} → {typed} at {width}: {spec:?}");
            // The field being edited is the highlighted callout and label.
            let link = ctx
                .data(|d| {
                    d.get_temp::<crate::modeling::diagram::DiagramLink>(egui::Id::new(
                        "profile-stock-diagram-link",
                    ))
                })
                .unwrap();
            let key = match callout.split(' ').next().unwrap() {
                "R" => picker::RADIUS,
                "Width" => picker::WIDTH,
                "Height" => picker::HEIGHT,
                _ => picker::WALL,
            };
            assert_eq!(link.emphasis().active, Some(key), "{callout} at {width}");
            let output = run(
                &ctx,
                egui::vec2(width, 4000.),
                vec![],
                &mut stock_form(&mut spec),
            );
            let selection = ctx.global_style().visuals.selection.stroke.color;
            let name = match key {
                picker::RADIUS => "Corner radius",
                picker::WIDTH => "Width",
                picker::HEIGHT => "Height",
                _ => "Wall",
            };
            let highlighted = output.shapes.iter().any(|s| {
                matches!(&s.shape,
                egui::Shape::Text(t) if t.galley.text() == name
                    && t.galley.job.sections.iter().any(|sec| sec.format.color == selection))
            });
            assert!(highlighted, "{name} label is highlighted at {width}");
        }
    }
}

#[test]
fn catalogue_names_get_a_plain_description_and_picture() {
    // Model tool stock rows: picture, plain description, name secondary.
    let line8 = profile(ProfileShape::ItemLine8_40x40, 40., 40., None, None);
    let ctx = egui::Context::default();
    let painted = harness::frame(&ctx, 420., vec![], &mut |ui| {
        picker::stock_row(ui, false, &line8, "item 0.0.026.03");
    });
    assert!(painted.has("Aluminium slot profile 40×40, 8 mm slot (item Line 8 system)"));
    assert!(
        painted.has("item 0.0.026.03"),
        "the system name stays as secondary text"
    );
    // The form names the chosen type in plain words above its diagram.
    let mut spec = profile(ProfileShape::Rectangular, 60., 40., Some(3.), Some(6.));
    let painted = harness::frame(&ctx, 420., vec![], &mut stock_form(&mut spec));
    assert!(painted.has("Rectangular tube 60×40×3 mm, corner radius R6 mm"));
    let mut line8_form = line8.clone();
    let painted = harness::frame(&ctx, 420., vec![], &mut stock_form(&mut line8_form));
    assert!(painted.has("Aluminium slot profile 40×40, 8 mm slot (item Line 8 system)"));
    assert!(
        painted.has("Slot profile 40×40 (item Line 8)"),
        "the picker names the family"
    );

    // Project → Materials tree rows: the tube stock row says what it is.
    let mut fixture = shared_fixture::fixture().unwrap();
    let mut editor = LibraryEditor::default();
    editor.set_backend(Box::new(
        crate::modeling::native_fixture::NativeBackend::new(),
    ));
    editor.open(&fixture.session);
    // Selecting the tube opens its Profile stock group in the tree.
    let tube = fixture
        .session
        .document()
        .get(fixture.stocks[2])
        .unwrap()
        .clone();
    editor.select(&fixture.session, tube);
    let ctx = egui::Context::default();
    let size = egui::vec2(1100., 900.);
    let mut texts_seen = vec![];
    for _ in 0..3 {
        let output = run(&ctx, size, vec![], &mut |ui| {
            editor.show_inline(ui, &mut fixture.session, &Default::default(), true)
        });
        texts_seen = texts(&output);
    }
    assert!(
        texts_seen
            .iter()
            .any(|(t, _)| t == "Rectangular tube 30×20×2 mm"),
        "{:?}",
        texts_seen.iter().map(|t| &t.0).collect::<Vec<_>>()
    );
    assert!(texts_seen.iter().any(|(t, _)| t == "Steel tube 30x20x2"));
}

#[test]
fn stock_length_follows_the_type_and_section_and_does_not_size_members() {
    let ctx = egui::Context::default();
    let mut spec = profile(ProfileShape::Rectangular, 40., 40., Some(2.), None);
    let painted = harness::frame(&ctx, 420., vec![], &mut stock_form(&mut spec));
    let order = [
        y_of(&painted.texts, "Section type"),
        y_of(&painted.texts, "Width 40 mm"),
        y_of(&painted.texts, "Designation"),
        y_of(&painted.texts, "Stock length"),
    ];
    assert!(order.windows(2).all(|w| w[0] < w[1]), "{order:?}");
    assert!(painted.has_prefix("The bar length you buy"));

    // A new profile stock starts with its type, not a length.
    let session = Session::new(Document::new("New stock").unwrap()).unwrap();
    let editor = LibraryEditor::for_model_stock(&session, ToolStock::Profile);
    let mut draft = draft_stock(&editor).specification;
    assert!(matches!(draft, StockKind::Profile { .. }));
    let painted = harness::frame(&ctx, 420., vec![], &mut stock_form(&mut draft));
    assert!(y_of(&painted.texts, "Section type") < y_of(&painted.texts, "Stock length"));

    // Members take their length from their path: changing the bar length
    // keeps the member Part's volume.
    let mut fixture = shared_fixture::fixture().unwrap();
    let mut editor = LibraryEditor::default();
    editor.set_backend(Box::new(
        crate::modeling::native_fixture::NativeBackend::new(),
    ));
    editor.open(&fixture.session);
    let stock = fixture
        .session
        .document()
        .get(fixture.stocks[2])
        .unwrap()
        .clone();
    editor.select(&fixture.session, stock);
    let before = part_volume(&fixture.session, fixture.parts[2]);
    set_stock(&mut editor, |stock| {
        if let StockKind::Profile { length_mm, .. } = &mut stock.specification {
            *length_mm = 2500.;
        }
    });
    apply(&mut editor, &mut fixture.session).unwrap();
    assert!((part_volume(&fixture.session, fixture.parts[2]) - before).abs() < 1e-6);
}

#[test]
fn tube_corner_radius_is_validated_and_shapes_the_member_exactly() {
    let mut fixture = shared_fixture::fixture().unwrap();
    let mut editor = LibraryEditor::default();
    editor.set_backend(Box::new(
        crate::modeling::native_fixture::NativeBackend::new(),
    ));
    editor.open(&fixture.session);
    let stock = fixture
        .session
        .document()
        .get(fixture.stocks[2])
        .unwrap()
        .clone();
    editor.select(&fixture.session, stock);
    let sharp = part_volume(&fixture.session, fixture.parts[2]);
    let length = sharp / (30. * 20. - 26. * 16.);
    let ctx = egui::Context::default();
    let size = egui::vec2(1100., 1400.);
    let radius = |editor: &mut LibraryEditor, r: f64| {
        set_stock(editor, |stock| {
            if let StockKind::Profile {
                corner_radius_mm, ..
            } = &mut stock.specification
            {
                *corner_radius_mm = Some(r);
            }
        });
    };
    let painted = |editor: &mut LibraryEditor, s: &mut Session| {
        let mut seen = vec![];
        for _ in 0..2 {
            seen = texts(&run(&ctx, size, vec![], &mut |ui| {
                editor.show_inline(ui, s, &Default::default(), true)
            }));
        }
        seen
    };
    // Too large for the 20 mm side: refused in the form and at Apply.
    radius(&mut editor, 12.);
    let seen = painted(&mut editor, &mut fixture.session);
    assert!(
        seen.iter()
            .any(|(t, _)| t == "The corner radius can be at most half the smaller side (10 mm)"),
        "{:?}",
        seen.iter().map(|t| &t.0).collect::<Vec<_>>()
    );
    assert!(apply(&mut editor, &mut fixture.session).is_err());
    assert!((part_volume(&fixture.session, fixture.parts[2]) - sharp).abs() < 1e-6);
    // R3 outside on a 2 mm wall: inside R1, shown and built exactly.
    radius(&mut editor, 3.);
    let seen = painted(&mut editor, &mut fixture.session);
    assert!(
        seen.iter()
            .any(|(t, _)| t == "Inside radius R1 mm (corner radius − wall)")
    );
    apply(&mut editor, &mut fixture.session).unwrap();
    let area = (600. - (4. - PI) * 9.) - (26. * 16. - (4. - PI) * 1.);
    let rounded = part_volume(&fixture.session, fixture.parts[2]);
    assert!(
        (rounded - area * length).abs() < 1e-6 * rounded,
        "{rounded} vs {}",
        area * length
    );
    // R1.5 is below the wall: the bore keeps sharp corners.
    radius(&mut editor, 1.5);
    let seen = painted(&mut editor, &mut fixture.session);
    assert!(
        seen.iter()
            .any(|(t, _)| t == "Inside corners stay sharp (corner radius ≤ wall)")
    );
}

fn draft_stock(editor: &LibraryEditor) -> Stock {
    let Some(LibraryDefinition::Stock(stock)) =
        read(&editor.draft.as_ref().unwrap().entity).unwrap()
    else {
        panic!("stock draft")
    };
    stock
}
fn set_stock(editor: &mut LibraryEditor, f: impl FnOnce(&mut Stock)) {
    let d = editor.draft.as_mut().unwrap();
    let EntityData::Extension {
        definition: wire,
        inputs,
    } = &mut d.entity.data
    else {
        panic!("stock draft");
    };
    let LibraryDefinition::Stock(mut stock) = serde_json::from_value(wire.data.clone()).unwrap()
    else {
        panic!("stock draft");
    };
    f(&mut stock);
    let definition = LibraryDefinition::Stock(stock);
    wire.data = serde_json::to_value(&definition).unwrap();
    *inputs = definition.links();
}
fn apply(editor: &mut LibraryEditor, s: &mut Session) -> Result<(), String> {
    let resources = ProjectResources::default();
    editor.refresh_preview(s, &resources, Channel::Ui)?;
    let deadline = Instant::now() + Duration::from_secs(15);
    while editor.pending() {
        assert!(
            Instant::now() < deadline,
            "native definition preview timed out"
        );
        editor.poll(s, &resources);
        std::thread::sleep(Duration::from_millis(1));
    }
    if let Some(error) = &editor.error {
        return Err(error.clone());
    }
    editor.apply(s, &resources, Channel::Ui)
}
fn part_volume(s: &Session, part: EntityId) -> f64 {
    multicad_model::Compiler::default()
        .compile(s.document(), s.revision())
        .unwrap()
        .solid(part)
        .unwrap()
        .volume_mm3()
        .unwrap()
}
