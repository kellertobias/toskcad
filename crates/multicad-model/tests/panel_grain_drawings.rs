//! CAD-365: faint panel grain in 2D drawings follows each panel's canonical
//! grain axis (including after a rotate), is clipped at the perimeter and a
//! through hole, is a persisted page setting that never touches geometry or
//! cut-list grain, and reaches the PDF only while enabled.
use multicad_core::{
    commands::{Actor, Channel, Prepared, Session},
    document::{Document, Entity, EntityData, EntityId},
    drawings::{
        self, Annotation, AnnotationKind, DrawingEvaluation, DrawingVector, DrawingView, Page,
        PageAction, Scale, VectorCurve, VectorRole, ViewEvaluation, ViewKind, ViewSource,
    },
    materials::*,
    model::*,
    project_file::{self, ProjectFile},
    references::ReferenceFrame,
};
use multicad_model::{Compiler, DrawingProviders, GRAIN_RGB};
use multicad_precision::PortableMath as _;

fn actor() -> Actor {
    Actor::local(Channel::Ui)
}
fn checked(s: &mut Session, prepared: Prepared) {
    Compiler::default()
        .compile(prepared.document(), s.revision() + 1)
        .unwrap();
    s.commit(prepared).unwrap();
}
fn panel(stock: EntityId, name: &str, origin: [f64; 3], grain: GrainAxis) -> PanelRequest {
    PanelRequest::new(
        name,
        PanelDefinition {
            output: BodyOutput {
                body: EntityId::new(),
            },
            source: PanelSource::Rectangle {
                stock,
                plane: PanelPlane::Xy,
                origin: origin.map(length),
                length: length(600.),
                width: length(400.),
            },
            grain,
            edges: PanelEdges::default(),
            oversize: None,
        },
    )
}

struct Fixture {
    s: Session,
    panels: [PanelRequest; 4],
    book: EntityId,
    page: EntityId,
    view: EntityId,
}
fn fixture() -> Fixture {
    let mut s = Session::new(Document::new("Grain drawing").unwrap()).unwrap();
    let birch = Entity::new(
        "Birch",
        EntityData::Material {
            density_kg_m3: None,
            color_rgba: [210, 180, 140, 255],
            section_hatch: None,
        },
    );
    let stock = LibraryDefinition::Stock(Stock {
        material: birch.id,
        specification: StockKind::Panel {
            thickness_mm: 18.,
            length_mm: 2500.,
            width_mm: 1250.,
            grain: Grain::AlongLength,
        },
    })
    .entity(EntityId::new(), "Ply 18")
    .unwrap();
    let stock_id = stock.id;
    s.commit(prepare_library(&s, 0, actor(), &[birch, stock]).unwrap())
        .unwrap();
    // Two panels share one stock but keep their own grain.
    let panels = [
        panel(stock_id, "Along length", [0., 0., 0.], GrainAxis::U),
        panel(stock_id, "Along width", [800., 0., 0.], GrainAxis::V),
        panel(stock_id, "Rotated", [0., 700., 0.], GrainAxis::U),
        panel(stock_id, "Any", [800., 700., 0.], GrainAxis::Free),
    ];
    for p in &panels {
        let prepared = prepare_panel(&s, s.revision(), actor(), p).unwrap();
        checked(&mut s, prepared.prepared);
    }
    // A Ø40 through hole in the first panel.
    let hole = prepare_machining(
        &s,
        s.revision(),
        actor(),
        "Cable hole",
        &MachiningDefinition {
            generated_by: None,
            suppressed: false,
            kind: MachiningKind::Drill {
                holes: vec![[length(300.), length(200.)]],
                diameter: length(40.),
                depth: None,
            },
            face: MachiningFace::Top,
            targets: vec![panels[0].body],
        },
    )
    .unwrap();
    checked(&mut s, hole.model.prepared);
    // Rotate the third panel 30° about +Z through its center.
    let rotated = panels[2].body;
    let producer = match &s.document().get(rotated).unwrap().data {
        EntityData::Part { features, .. } => *features.last().unwrap(),
        _ => unreachable!(),
    };
    let transform = TransformRequest::new(
        BodyInput {
            body: rotated,
            producer,
        },
        NumericTransform {
            translation: [length(0.), length(0.), length(0.)],
            pivot: [length(300.), length(900.), length(0.)],
            axis: [0., 0., 1.],
            angle: angle(30.),
        },
    );
    let prepared = prepare_transform(&s, s.revision(), actor(), &transform).unwrap();
    checked(&mut s, prepared.prepared);

    let view = DrawingView {
        id: EntityId::new(),
        name: "Top".into(),
        source: ViewSource::Model {
            entities: panels.iter().map(|p| p.body).collect(),
        },
        frame: ReferenceFrame {
            origin_mm: [0.; 3],
            u: [1., 0., 0.],
            v: [0., 1., 0.],
        },
        kind: ViewKind::Projection,
        position_mm: [60., 150.],
        scale: Scale {
            numerator: 1,
            denominator: 10,
        },
        hidden_lines: false,
        tangent_edges: false,
        paper_chord_tolerance_mm: 0.01,
        guide: None,
        equipment_categories: None,
    };
    let view_id = view.id;
    let mut page = Page::new("1");
    page.views.push(view);
    page.annotations.push(Annotation {
        id: EntityId::new(),
        kind: AnnotationKind::DraftDimension {
            a: [40., 30.],
            b: [100., 30.],
            offset_mm: -6.,
            precision: 0,
        },
    });
    let page_id = page.id;
    let (book, p) =
        drawings::prepare_create(&s, s.revision(), actor(), "Panels", Some(page)).unwrap();
    s.commit(p).unwrap();
    Fixture {
        s,
        panels,
        book,
        page: page_id,
        view: view_id,
    }
}
fn evaluate(s: &Session) -> DrawingEvaluation {
    Compiler::default()
        .compile(s.document(), s.revision())
        .unwrap()
        .drawings(s.document(), s.revision(), &DrawingProviders::default())
        .unwrap()
}
fn vectors(e: &DrawingEvaluation, view: EntityId) -> &[DrawingVector] {
    match &e.views[&view] {
        ViewEvaluation::Ready { vectors, .. } => vectors,
        other => panic!("{other:?}"),
    }
}
fn set_grain(f: &mut Fixture, on: bool) {
    let mut page = drawings::read(f.s.document(), f.book).unwrap().pages[0].clone();
    page.settings.panel_grain = on;
    let prepared = drawings::prepare(
        &f.s,
        f.s.revision(),
        actor(),
        f.book,
        PageAction::Replace(page),
    )
    .unwrap();
    f.s.commit(prepared).unwrap();
}
fn line(v: &DrawingVector) -> ([f64; 2], [f64; 2]) {
    match v.curve {
        VectorCurve::Line { a, b } => (a, b),
        _ => panic!("grain is straight"),
    }
}
const GRAY_RG: &str = "0.769 0.769 0.769 RG";

#[test]
fn page_grain_follows_each_panel_clips_holes_and_reaches_pdf_only_when_enabled() {
    let mut f = fixture();
    // Off by default: no grain anywhere, and the PDF has no grain stroke.
    let off = evaluate(&f.s);
    assert!(
        vectors(&off, f.view)
            .iter()
            .all(|v| v.role != VectorRole::Grain)
    );
    let pdf_off =
        drawings::pdf::book_pdf(f.s.document(), f.s.revision(), &[(f.book, f.page)], &off).unwrap();
    assert!(!String::from_utf8_lossy(&pdf_off).contains(GRAY_RG));

    let model_before: Vec<_> =
        f.s.document()
            .entities()
            .filter(|e| e.id != f.book)
            .cloned()
            .collect();
    let cut_before = Compiler::default()
        .compile(f.s.document(), f.s.revision())
        .unwrap()
        .cut_list(f.s.document())
        .unwrap();
    set_grain(&mut f, true);
    let on = evaluate(&f.s);
    let all = vectors(&on, f.view);
    let grain: Vec<_> = all.iter().filter(|v| v.role == VectorRole::Grain).collect();
    assert!(!grain.is_empty());
    // Grain draws first, beneath every outline.
    let first_outline = all
        .iter()
        .position(|v| v.role != VectorRole::Grain)
        .unwrap();
    assert!(
        all[..first_outline]
            .iter()
            .all(|v| v.role == VectorRole::Grain)
    );
    assert!(
        grain
            .iter()
            .all(|v| v.color_rgb == Some(GRAIN_RGB) && v.line_width_mm < 0.1)
    );

    let of = |i: usize| -> Vec<_> {
        grain
            .iter()
            .filter(|v| v.source == f.panels[i].body)
            .map(|v| line(v))
            .collect()
    };
    let dir = |(a, b): ([f64; 2], [f64; 2])| {
        let (dx, dy) = (b[0] - a[0], b[1] - a[1]);
        let l = dx.portable_hypot(dy);
        [dx / l, dy / l]
    };
    // Along length: horizontal, inside 0..600 × 0..400, never inside the hole.
    let along = of(0);
    assert!(along.len() > 10);
    for &(a, b) in &along {
        assert!(dir((a, b))[1].abs() < 1e-9);
        for p in [a, b] {
            assert!((-1e-6..=600. + 1e-6).contains(&p[0]) && (-1e-6..=400. + 1e-6).contains(&p[1]));
        }
        let y = a[1] - 200.;
        if y.abs() < 19. {
            // The chord through the hole is split around it.
            let half = (20f64.powi(2) - y * y).sqrt();
            let (lo, hi) = (a[0].min(b[0]), a[0].max(b[0]));
            assert!(
                hi <= 300. - half + 0.05 || lo >= 300. + half - 0.05,
                "grain crosses the hole: {a:?} {b:?}"
            );
        }
    }
    assert!(
        along
            .iter()
            .any(|(a, b)| (a[1] - 200.).abs() < 19. && a[0].max(b[0]) < 300.)
            && along
                .iter()
                .any(|(a, b)| (a[1] - 200.).abs() < 19. && a[0].min(b[0]) > 300.),
        "grain resumes on both sides of the hole"
    );
    // Along width: vertical.
    let across = of(1);
    assert!(!across.is_empty());
    assert!(across.iter().all(|&l| dir(l)[0].abs() < 1e-9));
    // Rotated 30° about +Z: the length grain turned with the panel.
    let rotated = of(2);
    assert!(!rotated.is_empty());
    for &l in &rotated {
        let d = dir(l);
        let along =
            d[0] * 30f64.to_radians().portable_cos() + d[1] * 30f64.to_radians().portable_sin();
        assert!(along.abs() > 1. - 1e-9, "{d:?}");
    }
    // Any: no fixed direction, no grain.
    assert!(of(3).is_empty());

    // Display only: geometry, panels and cut-list grain are unchanged.
    let model_after: Vec<_> =
        f.s.document()
            .entities()
            .filter(|e| e.id != f.book)
            .cloned()
            .collect();
    assert_eq!(model_before, model_after);
    let cut_after = Compiler::default()
        .compile(f.s.document(), f.s.revision())
        .unwrap()
        .cut_list(f.s.document())
        .unwrap();
    assert_eq!(cut_before, cut_after);
    let grain_of = |part| cut_after.iter().find(|r| r.part == part).unwrap().grain;
    assert_eq!(
        f.panels.each_ref().map(|p| grain_of(p.body)),
        [GrainAxis::U, GrainAxis::V, GrainAxis::U, GrainAxis::Free]
    );

    // PDF carries the faint gray grain together with the draft dimension.
    let render =
        drawings::render::render_page(f.s.document(), f.s.revision(), f.book, f.page, &on).unwrap();
    assert!(
        render
            .polylines
            .iter()
            .any(|p| p.layer == drawings::render::Layer::Grain)
    );
    assert!(
        render
            .polylines
            .iter()
            .any(|p| p.layer == drawings::render::Layer::Dimension)
            || !render.texts.is_empty()
    );
    let pdf_on =
        drawings::pdf::book_pdf(f.s.document(), f.s.revision(), &[(f.book, f.page)], &on).unwrap();
    assert!(String::from_utf8_lossy(&pdf_on).contains(GRAY_RG));
    // Visual review: MULTICAD_GRAIN_REVIEW=<dir> writes both PDFs.
    if let Ok(directory) = std::env::var("MULTICAD_GRAIN_REVIEW") {
        std::fs::create_dir_all(&directory).unwrap();
        let dir = std::path::Path::new(&directory);
        std::fs::write(dir.join("grain-on.pdf"), &pdf_on).unwrap();
        std::fs::write(dir.join("grain-off.pdf"), &pdf_off).unwrap();
    }

    // Save/reopen keeps the setting and declares its capability; Undo
    // restores the page without grain.
    let mut bytes = vec![];
    let file = ProjectFile::new(f.s.document().clone()).unwrap();
    assert!(
        file.capabilities()
            .contains(drawings::PANEL_GRAIN_CAPABILITY)
    );
    project_file::write(&file, &mut bytes, Default::default()).unwrap();
    let opened = project_file::read(
        &mut std::io::Cursor::new(bytes),
        Default::default(),
        &project_file::supported_capabilities(),
    )
    .unwrap();
    let reopened = drawings::read(&opened.project.document, f.book).unwrap();
    assert!(reopened.pages[0].settings.panel_grain);
    f.s.undo(f.s.revision(), actor()).unwrap();
    let book = drawings::read(f.s.document(), f.book).unwrap();
    assert!(!book.pages[0].settings.panel_grain);
    let json = serde_json::to_string(&book).unwrap();
    assert!(!json.contains("panel_grain"), "default stays omitted");
}

#[test]
fn pages_without_the_setting_open_and_unknown_values_are_rejected() {
    let mut settings = serde_json::to_value(drawings::PageSettings::default()).unwrap();
    assert!(settings.get("panel_grain").is_none());
    let parsed: drawings::PageSettings = serde_json::from_value(settings.clone()).unwrap();
    assert!(!parsed.panel_grain);
    settings["panel_grain"] = serde_json::json!("faint");
    assert!(serde_json::from_value::<drawings::PageSettings>(settings).is_err());
    assert!(serde_json::from_value::<GrainAxis>(serde_json::json!("diagonal")).is_err());
}
