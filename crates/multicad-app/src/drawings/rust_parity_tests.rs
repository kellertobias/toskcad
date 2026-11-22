//! CAD-346 drawing parity through the normal desktop drawing worker: a drilled
//! plate and a straight spur gear drawn from saved model features, followed
//! through source edits, stale requests, cancellation and project reopening.
use super::tests::{add_box, preview, projection_page};
use super::*;
use multicad_core::{
    commands::{Actor, Channel, Session},
    document::Document,
    drawings::{self, VectorCurve, VectorRole, ViewEvaluation},
    model::{self, PrimitiveFrame, PrimitiveParameters, PrimitiveRequest, length},
};
use std::time::{Duration, Instant};

struct Fixture {
    session: Session,
    plate: EntityId,
    drill: PrimitiveRequest,
    gear: EntityId,
    page: drawings::Page,
    views: [EntityId; 5],
    findings: drawings::Page,
    book: EntityId,
}

fn drill(radius: f64) -> PrimitiveRequest {
    PrimitiveRequest::new(
        "Drill",
        PrimitiveParameters::Cylinder {
            radius: length(radius),
            height: length(9.),
            frame: PrimitiveFrame {
                origin: [length(10.), length(5.), length(-2.)],
                ..Default::default()
            },
        },
    )
}

fn fixture() -> Fixture {
    let actor = Actor::local(Channel::Ui);
    let mut session = Session::new(Document::new("Drawing parity").unwrap()).unwrap();
    // 20 × 10 × 5 plate with a radius-2 through bore at (10, 5).
    let plate = add_box(&mut session, "Plate", 0.);
    let drill = drill(2.);
    let prepared =
        model::prepare_primitive(&session, session.revision(), actor.clone(), &drill).unwrap();
    session.commit(prepared.prepared).unwrap();
    let cut = model::BooleanRequest::new(
        model::Operation::Cut,
        plate.clone(),
        vec![model::BodyInput {
            body: drill.body,
            producer: drill.feature,
        }],
    );
    let prepared =
        model::prepare_boolean(&session, session.revision(), actor.clone(), &cut).unwrap();
    session.commit(prepared.prepared).unwrap();
    let gear = model::GearRequest::new(
        "Spur",
        model::GearParameters {
            family: model::GearFamily::Spur,
            face_width: length(4.),
            bore_diameter: length(6.),
            frame: PrimitiveFrame {
                origin: [length(60.), length(0.), length(0.)],
                ..Default::default()
            },
            ..Default::default()
        },
    );
    let prepared = model::prepare_gear(&session, session.revision(), actor.clone(), &gear).unwrap();
    session.commit(prepared.prepared).unwrap();

    let mut page = projection_page(plate.body);
    let axial = page.views[0].id;
    page.views[0].hidden_lines = true;
    let section = |name: &str, source, frame, offset_mm, range| drawings::DrawingView {
        id: EntityId::new(),
        name: name.into(),
        source: drawings::ViewSource::Model {
            entities: vec![source],
        },
        frame,
        kind: drawings::ViewKind::Section {
            offset_mm,
            range,
            window: Default::default(),
        },
        hidden_lines: false,
        ..page.views[0].clone()
    };
    let top = multicad_core::references::ReferenceFrame {
        origin_mm: [0.; 3],
        u: [1., 0., 0.],
        v: [0., 1., 0.],
    };
    let front = multicad_core::references::ReferenceFrame {
        origin_mm: [0., 5., 0.],
        u: [1., 0., 0.],
        v: [0., 0., 1.],
    };
    let views = [
        section(
            "Plate mid-depth section",
            plate.body,
            top,
            2.5,
            drawings::SectionRange::Depth { depth_mm: 1. },
        ),
        section(
            "Plate section through the bore axis",
            plate.body,
            front,
            0.,
            drawings::SectionRange::Unlimited,
        ),
        section(
            "Gear mid-face section",
            gear.body,
            top,
            2.,
            drawings::SectionRange::SectionOnly,
        ),
        section(
            "Gear retained section",
            gear.body,
            top,
            2.,
            drawings::SectionRange::Unlimited,
        ),
    ];
    let mut hidden = section(
        "Plate section with hidden lines",
        plate.body,
        front,
        0.,
        drawings::SectionRange::Unlimited,
    );
    hidden.hidden_lines = true;
    let ids = [axial, views[0].id, views[1].id, views[2].id, views[3].id];
    page.views.extend(views);
    // Unresolved views are never cached, so the finding lives on its own page.
    let mut findings = drawings::Page::new("2");
    findings.views.push(hidden);
    let mut books = Vec::new();
    for (name, page) in [("Parity drawing", &page), ("Parity findings", &findings)] {
        let (book, prepared) = drawings::prepare_create(
            &session,
            session.revision(),
            actor.clone(),
            name,
            Some(page.clone()),
        )
        .unwrap();
        session.commit(prepared).unwrap();
        books.push(book);
    }
    Fixture {
        session,
        plate: plate.body,
        drill,
        gear: gear.body,
        page,
        views: ids,
        findings,
        book: books[0],
    }
}

fn vectors(result: &DrawingEvaluation, view: EntityId) -> &[drawings::DrawingVector] {
    let ViewEvaluation::Ready { vectors, .. } = &result.views[&view] else {
        panic!("{:?}", result.views[&view]);
    };
    vectors
}

fn bore_radius(vectors: &[drawings::DrawingVector], role: VectorRole) -> Vec<f64> {
    vectors
        .iter()
        .filter(|v| v.role == role)
        .filter_map(|v| match v.curve {
            VectorCurve::Ellipse {
                major_mm,
                minor_mm,
                start_rad,
                end_rad,
                ..
            } if (major_mm - minor_mm).abs() < 1e-9
                && ((end_rad - start_rad).abs() - std::f64::consts::TAU).abs() < 1e-9 =>
            {
                Some(major_mm)
            }
            _ => None,
        })
        .collect()
}

fn check(result: &DrawingEvaluation, f: &Fixture, radius: f64) {
    let [axial, depth, front, gear_only, gear_retained] = f.views;
    let axial = vectors(result, axial);
    assert!(axial.iter().all(|v| v.source == f.plate));
    assert_eq!(bore_radius(axial, VectorRole::Visible), [radius]);
    let depth = vectors(result, depth);
    assert_eq!(bore_radius(depth, VectorRole::Section), [radius]);
    // The plate has no material, so the model publishes no hatch; hatch
    // regions are covered by the kernel regressions.
    let front = vectors(result, front);
    let mut cut_edges: Vec<_> = front
        .iter()
        .filter(|v| v.role == VectorRole::Section)
        .filter_map(|v| match v.curve {
            VectorCurve::Line { a, b } if (a[0] - b[0]).abs() < 1e-9 => Some(a[0]),
            _ => None,
        })
        .collect();
    cut_edges.sort_by(f64::total_cmp);
    cut_edges.dedup_by(|a, b| (*a - *b).abs() < 1e-9);
    assert_eq!(cut_edges, [0., 10. - radius, 10. + radius, 20.]);
    let gear_only = vectors(result, gear_only);
    assert!(gear_only.iter().all(|v| v.source == f.gear));
    assert!(gear_only.iter().all(|v| v.role == VectorRole::Section));
    assert!(
        gear_only
            .iter()
            .any(|v| matches!(v.curve, VectorCurve::BSpline { degree: 3, .. }))
    );
    assert_eq!(bore_radius(gear_only, VectorRole::Section), [3.]);
    let retained = vectors(result, gear_retained);
    assert_eq!(
        retained
            .iter()
            .filter(|v| v.role == VectorRole::Visible)
            .count(),
        retained
            .iter()
            .filter(|v| v.role == VectorRole::Section)
            .count()
    );
}

fn stamp(session: &Session) -> DrawingStamp {
    DrawingStamp {
        document: session.document().id(),
        revision: session.revision(),
        source_hash: multicad_core::project_file::ContentHash::document(session.document())
            .unwrap(),
    }
}

#[test]
fn rust_worker_draws_bored_stock_and_gear_sections_through_edits_cancellation_and_reopen() {
    use multicad_core::project_file::ProjectFile;
    use multicad_ui::{application::DocumentStorage, storage::SaveMode};

    let directory = tempfile::tempdir().unwrap();
    let cache = cache::Cache::isolated(directory.path().join("drawing-cache"));
    let mut f = fixture();
    let mut worker = Backend::with_cache_and_backend(Some(cache.clone()), GeometryBackend::Rust);
    let started = Instant::now();
    let original = preview(&mut worker, &f.session, f.page.id);
    let cold = started.elapsed();
    check(&original, &f, 2.);
    let started = Instant::now();
    assert_eq!(
        preview(&mut worker, &f.session, f.page.id).views,
        original.views
    );
    let warm = started.elapsed();
    // Hidden edges behind a curved retained section draw through the scene
    // engine (CAD-346). The plane y = 5 through the bore axis retains the
    // half plate behind it: its rear face and the two half rims of the bore
    // project onto the outline of the cut face and the cap lines across the
    // open bore, and the bore's outline generators lie in the cut plane.
    // Every rear edge therefore coincides with a visible one and the view
    // has no hidden line.
    let finding = preview(&mut worker, &f.session, f.findings.id);
    let hidden = vectors(&finding, f.findings.views[0].id);
    for role in [VectorRole::Section, VectorRole::Visible] {
        assert!(hidden.iter().any(|vector| vector.role == role), "{role:?}");
    }
    assert!(
        hidden
            .iter()
            .all(|vector| vector.role != VectorRole::Hidden)
    );

    // A cancelled request publishes nothing, and the worker accepts the next one.
    worker
        .request_page(f.session.document(), stamp(&f.session), f.page.id)
        .unwrap();
    worker.cancel();
    let deadline = Instant::now() + Duration::from_millis(500);
    while Instant::now() < deadline {
        assert!(worker.poll().is_none(), "a cancelled page must stay silent");
        std::thread::sleep(Duration::from_millis(5));
    }

    // A stale request superseded by an edit is never delivered.
    worker
        .request_page(f.session.document(), stamp(&f.session), f.page.id)
        .unwrap();
    let stale = stamp(&f.session);
    let prepared = model::prepare_primitive_edit(
        &f.session,
        f.session.revision(),
        Actor::local(Channel::Ui),
        f.drill.feature,
        &drill(3.).parameters,
    )
    .unwrap();
    f.session.commit(prepared.prepared).unwrap();
    let edited = preview(&mut worker, &f.session, f.page.id);
    assert_ne!(edited.stamp, stale);
    assert_eq!(edited.stamp.revision, f.session.revision());
    assert_ne!(edited.stamp.source_hash, original.stamp.source_hash);
    check(&edited, &f, 3.);
    // Only the plate changed; the gear views are reproduced exactly.
    for view in &f.views[3..] {
        assert_eq!(edited.views[view], original.views[view]);
    }
    drop(worker);

    // Save, reopen from disk and restore the edited page from the cache.
    let path = directory.path().join("parity.tcd");
    let mut store = crate::files::ProjectStore::default();
    store
        .save_project(
            &path,
            &ProjectFile::new(f.session.document().clone()).unwrap(),
            SaveMode::Create,
        )
        .unwrap();
    let reopened = Session::new(store.open_project(&path).unwrap().file.document).unwrap();
    assert_eq!(reopened.document(), f.session.document());
    let mut cached = Backend::with_cache_and_backend(Some(cache), GeometryBackend::Rust);
    let started = Instant::now();
    let restored = preview(&mut cached, &reopened, f.page.id);
    let reopen = started.elapsed();
    assert_eq!(restored.views, edited.views);
    assert!(
        cached
            .runtime
            .as_ref()
            .unwrap()
            .last_compile_stats()
            .is_none()
    );
    // A fresh worker without a cache rebuilds the identical page.
    let mut fresh = Backend::with_cache_and_backend(None, GeometryBackend::Rust);
    let started = Instant::now();
    assert_eq!(
        preview(&mut fresh, &reopened, f.page.id).views,
        edited.views
    );
    let rebuilt = started.elapsed();
    println!(
        "{{\"case\":\"CAD-346 parity page\",\"cold_ms\":{},\"warm_ms\":{},\"cached_reopen_ms\":{},\"uncached_reopen_ms\":{}}}",
        cold.as_secs_f64() * 1e3,
        warm.as_secs_f64() * 1e3,
        reopen.as_secs_f64() * 1e3,
        rebuilt.as_secs_f64() * 1e3,
    );
    exports(&reopened, &f, &restored);
}

/// Render the reopened page to the SVG, PDF and DXF the export workspace
/// writes, then read the actual files back for scale and curve content.
fn exports(session: &Session, f: &Fixture, evaluation: &DrawingEvaluation) {
    use multicad_core::drawings::{dxf, pdf, render, svg};
    let page = render::render_page(
        session.document(),
        session.revision(),
        f.book,
        f.page.id,
        evaluation,
    )
    .unwrap();
    // View geometry is published as its preview polyline within the view's
    // paper chord tolerance; only annotations keep exact paper arcs. At 1:1
    // the edited bore is a closed ring of points 3 mm from its center.
    let tolerance = f.page.views[0].paper_chord_tolerance_mm + 1e-9;
    let is_bore = |points: &[[f64; 2]]| {
        // A closed preview repeats its first point; count it once.
        let points = match points {
            [first, rest @ .., last] if first == last && !rest.is_empty() => {
                &points[..points.len() - 1]
            }
            _ => points,
        };
        let n = points.len() as f64;
        let center = [0, 1].map(|i| points.iter().map(|p| p[i]).sum::<f64>() / n);
        points.len() > 8
            && points
                .iter()
                .all(|p| ((p[0] - center[0]).hypot(p[1] - center[1]) - 3.).abs() <= tolerance)
    };
    assert!(page.polylines.iter().any(|line| is_bore(&line.points)));
    let directory = tempfile::tempdir().unwrap();
    let svg_path = directory.path().join("parity.svg");
    let pdf_path = directory.path().join("parity.pdf");
    let dxf_path = directory.path().join("parity.dxf");
    std::fs::write(&svg_path, svg::render_svg(&page)).unwrap();
    std::fs::write(
        &pdf_path,
        pdf::render_pdf(std::slice::from_ref(&page), "Parity"),
    )
    .unwrap();
    std::fs::write(&dxf_path, dxf::render_dxf(&page, 1., 4)).unwrap();
    if let Some(keep) = std::env::var_os("CAD346_ARTIFACTS") {
        std::fs::create_dir_all(&keep).unwrap();
        for path in [&svg_path, &pdf_path, &dxf_path] {
            std::fs::copy(
                path,
                std::path::Path::new(&keep).join(path.file_name().unwrap()),
            )
            .unwrap();
        }
    }

    let svg = std::fs::read_to_string(&svg_path).unwrap();
    let [width, height] = page.size_mm;
    assert!(svg.contains(&format!(r#"width="{width}mm" height="{height}mm""#)));
    let svg_bore = svg.split("points=\"").skip(1).any(|rest| {
        let points: Vec<[f64; 2]> = rest
            .split('"')
            .next()
            .unwrap()
            .split_whitespace()
            .filter_map(|pair| {
                let (x, y) = pair.split_once(',')?;
                Some([x.parse().ok()?, y.parse().ok()?])
            })
            .collect();
        is_bore(&points)
    });
    assert!(svg_bore, "the SVG file keeps the 3 mm bore");

    let dxf = std::fs::read_to_string(&dxf_path).unwrap();
    let pairs: Vec<_> = dxf.lines().map(str::trim).collect();
    assert!(dxf.contains("$INSUNITS"));
    let lines: Vec<[[f64; 2]; 2]> = pairs
        .windows(12)
        .filter(|w| w[0] == "0" && w[1] == "LINE")
        .filter_map(|w| {
            let value = |code: &str| {
                w.chunks(2)
                    .find(|pair| pair[0] == code)
                    .and_then(|pair| pair[1].parse::<f64>().ok())
            };
            Some([[value("10")?, value("20")?], [value("11")?, value("21")?]])
        })
        .collect();
    let segments: usize = page
        .polylines
        .iter()
        .map(|line| line.points.len() - 1 + usize::from(line.closed && line.points.len() > 2))
        .sum();
    assert_eq!(lines.len(), segments, "every paper segment is one DXF LINE");
    // Millimetres at 1:1: the plate keeps its exact 20 mm edges.
    assert!(
        lines
            .iter()
            .any(|[a, b]| ((a[0] - b[0]).hypot(a[1] - b[1]) - 20.).abs() < 1e-6),
        "plate edge length"
    );
    // DXF keeps millimetres with paper y upward: every paper segment of
    // the bore ring is one LINE at 1:1 within the file's 1e-6 rounding.
    let to_dxf = |q: [f64; 2]| [q[0], height - q[1]];
    let bore = page
        .polylines
        .iter()
        .find(|line| is_bore(&line.points))
        .expect("bore ring");
    for pair in bore.points.windows(2) {
        let [a, b] = [to_dxf(pair[0]), to_dxf(pair[1])];
        assert!(
            lines.iter().any(|[p, q]| {
                (p[0] - a[0]).abs() < 1e-5
                    && (p[1] - a[1]).abs() < 1e-5
                    && (q[0] - b[0]).abs() < 1e-5
                    && (q[1] - b[1]).abs() < 1e-5
            }),
            "the DXF file keeps the 3 mm bore segment {a:?} → {b:?}"
        );
    }

    let pdf = std::fs::read(&pdf_path).unwrap();
    assert!(pdf.starts_with(b"%PDF-"));
    let text = String::from_utf8_lossy(&pdf);
    let to_mm = |pt: f64| pt * 25.4 / 72.;
    let media: Vec<f64> = text
        .split("/MediaBox [")
        .nth(1)
        .and_then(|rest| rest.split(']').next())
        .unwrap()
        .split_whitespace()
        .map(|v| v.parse().unwrap())
        .collect();
    assert!(
        (to_mm(media[2]) - width).abs() < 1e-3 && (to_mm(media[3]) - height).abs() < 1e-3,
        "PDF page size {media:?} matches the {width} × {height} mm sheet"
    );
    // Each stroked path is `x y m x y l … s|S` in points with y upward.
    let pdf_bore = text.lines().any(|line| {
        let tokens: Vec<_> = line.split_whitespace().collect();
        let mut points = Vec::new();
        for (i, token) in tokens.iter().enumerate() {
            if (*token == "m" || *token == "l")
                && i >= 2
                && let (Ok(x), Ok(y)) = (tokens[i - 2].parse::<f64>(), tokens[i - 1].parse::<f64>())
            {
                points.push([to_mm(x), height - to_mm(y)]);
            }
        }
        is_bore(&points)
    });
    assert!(pdf_bore, "the PDF file strokes the 3 mm bore");
}
