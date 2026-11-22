//! CAD-346 acceptance through external control: the body families that draw
//! through the exact-curve scene engine are authored with the normal MCP
//! tools, drawn in a side elevation and a section through every axis by
//! drawing guides, edited and regenerated, refused when stale or cancelled,
//! saved and reopened. The oblique bore, its cube and the stock of the
//! internal thread come from the opened project, since MCP has no posed
//! primitive tools; everything else is created and edited over MCP.
use multicad_core::{
    commands::{Actor, Channel, Session},
    document::{Document, EntityId},
    drawings::{self, DrawingVector, VectorRole, ViewEvaluation},
    model::{self, PrimitiveFrame, PrimitiveParameters, PrimitiveRequest, length},
    project_file::{self, ProjectFile},
};
use multicad_mcp::{
    Headless, Storage, Workspace,
    tools::{Call, dispatch},
};
use multicad_model::{Compiler, GeometryBackend};
use multicad_precision::PortableMath as _;
use serde_json::{Value, json};
use std::{
    collections::BTreeMap,
    sync::{Arc, atomic::AtomicBool},
};

fn call(cancelled: bool) -> Call {
    Call {
        actor: Actor {
            id: "mcp:scene-families".into(),
            channel: Channel::Mcp,
        },
        cancel: Arc::new(AtomicBool::new(cancelled)),
    }
}

fn mm(v: f64) -> Value {
    json!({ "value": v, "unit": "mm" })
}

fn count(v: f64) -> Value {
    json!({ "value": v, "unit": "count" })
}

fn at(x: f64) -> Value {
    json!([mm(x), mm(0.), mm(0.)])
}

fn run(host: &Headless, name: &str, mut args: Value) -> Value {
    args["expected_revision"] = json!(revision(host));
    dispatch(host, &call(false), name, args.clone())
        .unwrap_or_else(|e| panic!("{name} {args}: {e:?}"))
}

fn revision(host: &Headless) -> u64 {
    dispatch(host, &call(false), "document_info", json!({})).unwrap()["revision"]
        .as_u64()
        .unwrap()
}

fn id(value: &Value) -> EntityId {
    serde_json::from_value(value.clone()).unwrap()
}

struct Store;
impl Storage for Store {
    fn open(&mut self, path: &std::path::Path) -> Result<ProjectFile, String> {
        project_file::read(
            &mut std::fs::File::open(path).map_err(|e| e.to_string())?,
            Default::default(),
            &project_file::supported_capabilities(),
        )
        .map(|p| p.project)
    }
    fn save(&mut self, path: &std::path::Path, p: &ProjectFile) -> Result<(), String> {
        project_file::write(
            p,
            &mut std::fs::File::create(path).map_err(|e| e.to_string())?,
            Default::default(),
        )
        .map_err(|e| e.to_string())
    }
}

/// The saved project MCP opens: a block overlapping the spur gear MCP will
/// create at x = 100, a 10 mm cube at x = 900 with a 0.5 mm bore along its
/// main diagonal, and cylinder stock at x = 1400 for the internal thread.
/// Returns the block, cube and stock Parts.
fn seed(path: &std::path::Path) -> [EntityId; 3] {
    let actor = Actor::local(Channel::Ui);
    let mut session = Session::new(Document::new("MCP scene families").unwrap()).unwrap();
    let add = |session: &mut Session, request: &PrimitiveRequest| {
        let prepared =
            model::prepare_primitive(session, session.revision(), actor.clone(), request).unwrap();
        session.commit(prepared.prepared).unwrap();
        model::BodyInput {
            body: request.body,
            producer: request.feature,
        }
    };
    let corner = |x: f64, y: f64, z: f64| PrimitiveFrame {
        origin: [length(x), length(y), length(z)],
        ..Default::default()
    };
    let block = add(
        &mut session,
        &PrimitiveRequest::new(
            "Overlapping block",
            PrimitiveParameters::Box {
                size: [length(20.), length(10.), length(6.)],
                frame: corner(118., -5., -1.),
            },
        ),
    );
    let cube = add(
        &mut session,
        &PrimitiveRequest::new(
            "Cube",
            PrimitiveParameters::Box {
                size: [length(10.), length(10.), length(10.)],
                frame: corner(900., -5., 0.),
            },
        ),
    );
    let (sx, ty) = (0.5_f64.sqrt(), (1_f64 / 6.).sqrt());
    let start = 1. / 3. - 3_f64.sqrt().recip();
    let drill = add(
        &mut session,
        &PrimitiveRequest::new(
            "Diagonal drill",
            PrimitiveParameters::Cylinder {
                radius: length(0.5),
                height: length(20.),
                frame: PrimitiveFrame {
                    origin: [length(900. + start), length(start - 5.), length(start)],
                    u: [sx, -sx, 0.],
                    v: [ty, ty, -2. * ty],
                },
            },
        ),
    );
    let cut = model::BooleanRequest::new(model::Operation::Cut, cube.clone(), vec![drill]);
    let prepared =
        model::prepare_boolean(&session, session.revision(), actor.clone(), &cut).unwrap();
    session.commit(prepared.prepared).unwrap();
    let stock = add(
        &mut session,
        &PrimitiveRequest::new(
            "Cylinder stock",
            PrimitiveParameters::Cylinder {
                radius: length(2.5),
                height: length(4.),
                frame: corner(1400., 0., 0.),
            },
        ),
    );
    Store
        .save(path, &ProjectFile::new(session.document().clone()).unwrap())
        .unwrap();
    [block.body, cube.body, stock.body]
}

fn open(path: &std::path::Path) -> Headless {
    Headless::with_geometry_backend(
        Workspace::project(path.to_owned(), Box::new(Store)).unwrap(),
        Some(path.parent().unwrap().into()),
        false,
        GeometryBackend::Rust,
    )
}

/// Evaluate the host's current drawings with the Rust backend, as the
/// guide tool does before it commits.
fn views(host: &Headless) -> BTreeMap<EntityId, ViewEvaluation> {
    let document = host.document();
    Compiler {
        omit_display: true,
        ..Default::default()
    }
    .compile_with_backend(&document, 0, GeometryBackend::Rust, None)
    .unwrap()
    .drawings(&document, 0, &multicad_model::DrawingProviders::default())
    .unwrap()
    .views
}

fn ready(views: &BTreeMap<EntityId, ViewEvaluation>, view: EntityId) -> &[DrawingVector] {
    match &views[&view] {
        ViewEvaluation::Ready { vectors, .. } => vectors,
        other => panic!("{other:?}"),
    }
}

fn of(vectors: &[DrawingVector], source: EntityId) -> Vec<DrawingVector> {
    vectors
        .iter()
        .filter(|v| v.source == source)
        .cloned()
        .collect()
}

/// u extents of a source's preview points.
fn span(vectors: &[DrawingVector]) -> [f64; 2] {
    vectors
        .iter()
        .flat_map(|v| &v.preview_mm)
        .fold([f64::INFINITY, f64::NEG_INFINITY], |s, p| {
            [s[0].min(p[0]), s[1].max(p[0])]
        })
}

/// Section points on the wall `u`, as their v span.
fn wall(vectors: &[DrawingVector], u: f64) -> f64 {
    let points: Vec<_> = vectors
        .iter()
        .filter(|v| v.role == VectorRole::Section)
        .flat_map(|v| &v.preview_mm)
        .filter(|p| (p[0] - u).abs() < 1e-6)
        .map(|p| p[1])
        .collect();
    points.iter().copied().fold(f64::NEG_INFINITY, f64::max)
        - points.iter().copied().fold(f64::INFINITY, f64::min)
}

/// A drawn family: its Part, axis x and exact side-elevation half-width
/// about the axis where it has a closed form.
struct Family {
    name: &'static str,
    part: EntityId,
    x: f64,
    half_width: Option<f64>,
}

fn guide(
    host: &Headless,
    name: &str,
    view: Value,
    sources: &[String],
    guide: Option<&str>,
) -> Value {
    let mut args = json!({ "view": view, "sources": sources,
        "scale_denominator": 1, "generate_pages": true });
    match guide {
        Some(guide) => args["guide"] = json!(guide),
        None => args["name"] = json!(name),
    }
    run(host, "set_drawing_guide", args)
}

#[test]
fn mcp_guides_draw_edit_and_reopen_scene_engine_families() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("scene-families.tcd");
    let [block, cube, stock] = seed(&path);
    let host = open(&path);
    let part = |value: &Value| id(&value["part"]);
    let gear = |host: &Headless, name: &str, args: Value| {
        let mut args = args;
        args["name"] = json!(name);
        run(host, "set_gear", args)
    };
    let mut families = Vec::new();
    let twisted = |family: &str, x: f64| {
        json!({ "family": family, "face_width": mm(8.), "bore_diameter": mm(6.),
            "helix_angle": { "value": 15., "unit": "deg" }, "hand": "right",
            "center_gap": mm(if family == "double_helical" { 2. } else { 0. }),
            "origin": at(x) })
    };
    let spur = gear(
        &host,
        "Spur",
        json!({ "family": "spur", "face_width": mm(4.), "bore_diameter": mm(6.), "origin": at(100.) }),
    );
    families.push(Family {
        name: "spur gear",
        part: part(&spur),
        x: 100.,
        half_width: Some(26.),
    });
    let ring = gear(
        &host,
        "Ring",
        json!({ "family": "internal_ring", "module": mm(1.), "teeth": count(48.),
            "bore_diameter": mm(0.), "outer_diameter": mm(56.), "face_width": mm(3.),
            "origin": at(200.) }),
    );
    families.push(Family {
        name: "internal ring",
        part: part(&ring),
        x: 200.,
        half_width: Some(28.),
    });
    let tip = 24. / 15_f64.to_radians().portable_cos() + 2.;
    let helical = gear(&host, "Helical", twisted("helical", 300.));
    families.push(Family {
        name: "helical gear",
        part: part(&helical),
        x: 300.,
        half_width: Some(tip),
    });
    for (name, family, x) in [
        ("herringbone gear", "herringbone", 400.),
        ("double-helical gear", "double_helical", 500.),
    ] {
        let created = gear(&host, name, twisted(family, x));
        families.push(Family {
            name,
            part: part(&created),
            x,
            half_width: Some(tip),
        });
    }
    let bevel = gear(
        &host,
        "Bevel",
        json!({ "family": "bevel", "teeth": count(24.), "face_width": mm(6.),
            "mate_teeth": count(32.), "shaft_angle": { "value": 90., "unit": "deg" },
            "origin": at(600.) }),
    );
    families.push(Family {
        name: "bevel gear",
        part: part(&bevel),
        x: 600.,
        half_width: Some(25.6),
    });
    let worm = gear(
        &host,
        "Worm",
        json!({ "family": "worm", "worm_role": "worm", "module": mm(1.),
            "mate_teeth": count(32.), "face_width": mm(20.), "bore_diameter": mm(2.),
            "origin": at(700.) }),
    );
    families.push(Family {
        name: "za worm",
        part: part(&worm),
        x: 700.,
        half_width: Some(6.),
    });
    let wheel = gear(
        &host,
        "Worm wheel",
        json!({ "family": "worm", "worm_role": "wheel", "module": mm(1.),
            "teeth": count(32.), "mate_teeth": count(1.), "face_width": mm(2.),
            "bore_diameter": mm(6.), "origin": at(800.) }),
    );
    families.push(Family {
        name: "worm wheel",
        part: part(&wheel),
        x: 800.,
        half_width: None,
    });
    families.push(Family {
        name: "oblique bore",
        part: cube,
        x: 905.,
        half_width: Some(5.),
    });
    let ends = |start: Value, end: Value| json!({ "start": start, "end": end });
    let relief = |width: f64| {
        ends(
            json!("flat_cut"),
            json!({ "cylindrical_relief": { "length": mm(width), "radial_depth": mm(0.1) } }),
        )
    };
    let metric = |name: &str, x: f64, treatments: Value| {
        json!({ "name": name, "family": "iso_metric", "fit": "metric_6g_external",
            "major_diameter": mm(6.), "pitch": mm(1.), "length": mm(6.6),
            "end_treatments": treatments, "origin": at(x) })
    };
    let relieved = run(
        &host,
        "set_thread",
        metric("M6 relieved", 1000., relief(0.6)),
    );
    families.push(Family {
        name: "thread end relief",
        part: part(&relieved),
        x: 1000.,
        half_width: Some(2.942),
    });
    let bsw = model::whitworth_thread_presets()
        .into_iter()
        .find(|p| p.family == model::ThreadFamily::Bsw)
        .unwrap();
    let rounded = run(
        &host,
        "set_thread",
        json!({ "name": "BSW rounded", "family": "bsw", "fit": "whitworth_medium_external",
            "major_diameter": { "value": bsw.nominal_diameter_inch.unwrap(), "unit": "in" },
            "pitch": { "value": 1. / f64::from(bsw.threads_per_inch), "unit": "in" },
            "length": mm(bsw.pitch_mm * 2.37), "origin": at(1100.) }),
    );
    families.push(Family {
        name: "rounded thread",
        part: part(&rounded),
        x: 1100.,
        half_width: None,
    });
    let chamfered = run(
        &host,
        "set_thread",
        metric(
            "M6 conical entry",
            1200.,
            ends(
                json!({ "chamfer45": { "radial_depth": mm(0.3) } }),
                json!("flat_cut"),
            ),
        ),
    );
    families.push(Family {
        name: "thread entry chamfer",
        part: part(&chamfered),
        x: 1200.,
        half_width: Some(2.942),
    });
    let screw = run(
        &host,
        "set_mechanical_part",
        json!({ "name": "Machine screw", "recipe": { "kind": "machine_screw",
            "underhead_length": mm(6.), "shank_diameter": mm(3.),
            "head": { "kind": "socket_cap", "diameter": mm(5.68), "height": mm(3.) },
            "drive": { "kind": "hex_socket", "across_flats": mm(2.55), "depth": mm(1.3) },
            "thread": { "major_diameter": mm(3.), "pitch": mm(0.5), "length": mm(6.),
                "origin": at(1300.) } } }),
    );
    families.push(Family {
        name: "machine screw",
        part: part(&screw),
        x: 1300.,
        half_width: Some(2.84),
    });
    run(
        &host,
        "set_thread",
        json!({ "name": "Treated M3 internal", "mode": "cut_internal",
            "target_part": stock, "family": "iso_metric", "fit": "metric_6h_internal",
            "major_diameter": mm(3.), "pitch": mm(0.5), "length": mm(4.),
            "minimum_wall": mm(0.5), "origin": at(1400.),
            "end_treatments": ends(
                json!({ "chamfer45": { "radial_depth": mm(0.5) } }),
                json!({ "cylindrical_relief": { "length": mm(0.5), "radial_depth": mm(0.125) } }),
            ) }),
    );
    families.push(Family {
        name: "internal thread treatment",
        part: stock,
        x: 1400.,
        half_width: Some(2.5),
    });
    // ISO 7 R 1/8; posed pipe constructions are not certified, so it stays
    // at the origin.
    let pipe = model::PipeThreadIntent {
        form: model::PipeThreadForm::R,
        size: model::PipeThreadSize::OneEighth,
        gauge: model::PipeGaugeConvention::Iso7Nominal1994,
    };
    let pipe_thread = run(
        &host,
        "set_thread",
        json!({ "name": "R 1/8", "family": "bspt", "fit": "bspt_nominal", "pipe": pipe,
            "major_diameter": mm(pipe.base_major_mm()), "pitch": mm(pipe.preset().pitch_mm),
            "length": mm(pipe.default_length_mm()), "representation": "modeled" }),
    );
    families.push(Family {
        name: "axial pipe thread",
        part: part(&pipe_thread),
        x: 0.,
        half_width: None,
    });

    // One side elevation and one section along y = 0 through every axis.
    let text = |part: EntityId| json!(part).as_str().unwrap().to_owned();
    let mut sources: Vec<String> = families.iter().map(|f| text(f.part)).collect();
    sources.push(text(block));
    let elevation_view = json!({ "kind": "elevation", "direction_deg": 90. });
    let section_view = json!({ "kind": "section", "start_mm": [-50., 0.], "end_mm": [1500., 0.] });
    let elevation = guide(&host, "Side", elevation_view.clone(), &sources, None);
    let section = guide(
        &host,
        "Through the axes",
        section_view.clone(),
        &sources,
        None,
    );
    let mut generated = Vec::new();
    for created in [&elevation, &section] {
        let generation = &created["generation"];
        assert!(
            generation["unresolved"].as_array().unwrap().is_empty(),
            "{created}"
        );
        let [view] = generation["created"].as_array().unwrap().as_slice() else {
            panic!("{created}");
        };
        generated.push(id(view));
    }
    let [side, cut] = [generated[0], generated[1]];
    let guides = [
        elevation["guide"].as_str().unwrap().to_owned(),
        section["guide"].as_str().unwrap().to_owned(),
    ];

    let block_owner = families[0].part;
    let check = |evaluated: &BTreeMap<EntityId, ViewEvaluation>, helical_bore: f64| {
        let side = ready(evaluated, side);
        let cut = ready(evaluated, cut);
        assert!(side.iter().all(|v| v.role == VectorRole::Visible));
        for family in &families {
            let elevation = of(side, family.part);
            let section = of(cut, family.part);
            assert!(!elevation.is_empty(), "{} elevation", family.name);
            assert!(
                section.iter().any(|v| v.role == VectorRole::Section),
                "{} section",
                family.name
            );
            let [lo, hi] = span(&elevation);
            if let Some(half) = family.half_width {
                // The overlapping block hides the spur gear's right flank.
                let hi = if family.part == block_owner {
                    family.x + half
                } else {
                    hi
                };
                assert!(
                    (lo - (family.x - half)).abs() < 0.011
                        && (hi - (family.x + half)).abs() < 0.011,
                    "{} elevation [{lo}, {hi}] vs {} ± {half}",
                    family.name,
                    family.x
                );
            }
            // The section frame starts at x = -50 on the plan line.
            let [lo, hi] = span(&section);
            assert!(
                lo >= family.x + 50. - 30. && hi <= family.x + 50. + 30.,
                "{} section [{lo}, {hi}]",
                family.name
            );
        }
        // The block overlaps the spur gear in both views.
        assert!(!of(side, block).is_empty() && !of(cut, block).is_empty());
        // The twisted gears' bore walls are cut along their whole 8 mm face.
        for family in &families[2..5] {
            let bore = if family.name == "helical gear" {
                helical_bore
            } else {
                3.
            };
            let section = of(cut, family.part);
            for u in [family.x + 50. - bore, family.x + 50. + bore] {
                assert!(
                    (wall(&section, u) - 8.).abs() < 1e-6,
                    "{} wall {u}",
                    family.name
                );
            }
        }
    };
    let original = views(&host);
    check(&original, 3.);

    // Edits regenerate both guided views; only the edited Parts change.
    let helical_part = families[2].part;
    let relief_part = families[9].part;
    run(
        &host,
        "set_gear",
        json!({ "gear": helical["gear"], "bore_diameter": mm(8.) }),
    );
    run(
        &host,
        "set_thread",
        json!({ "thread": relieved["thread"], "end_treatments": relief(1.1) }),
    );
    for (guide_id, kind, view) in [
        (&guides[0], &elevation_view, side),
        (&guides[1], &section_view, cut),
    ] {
        let updated = self::guide(&host, "", kind.clone(), &sources, Some(guide_id));
        let generation = &updated["generation"];
        assert_eq!(
            generation["updated"].as_array().unwrap().as_slice(),
            [json!(view)],
            "{updated}"
        );
        assert!(generation["created"].as_array().unwrap().is_empty());
        assert!(generation["unresolved"].as_array().unwrap().is_empty());
    }
    let edited = views(&host);
    check(&edited, 4.);
    for family in &families {
        for view in [side, cut] {
            let before = of(ready(&original, view), family.part);
            let after = of(ready(&edited, view), family.part);
            // The bore is invisible from the side; its section changes.
            let changed =
                family.part == relief_part || (family.part == helical_part && view == cut);
            assert_eq!(before != after, changed, "{} {view}", family.name);
        }
    }

    // Stale and cancelled requests leave the document untouched.
    let before = host.document();
    let current = revision(&host);
    let stale = dispatch(
        &host,
        &call(false),
        "set_drawing_guide",
        json!({ "expected_revision": current - 1, "guide": guides[1], "view": section_view,
            "sources": sources, "scale_denominator": 1, "generate_pages": true }),
    )
    .unwrap_err();
    assert_eq!(stale.code, "stale_revision", "{stale:?}");
    let cancelled = dispatch(
        &host,
        &call(true),
        "set_gear",
        json!({ "expected_revision": current, "gear": helical["gear"], "bore_diameter": mm(10.) }),
    )
    .unwrap_err();
    assert_eq!(cancelled.code, "cancelled");
    let cancelled = dispatch(
        &host,
        &call(true),
        "set_drawing_guide",
        json!({ "expected_revision": current, "guide": guides[0], "view": elevation_view,
            "sources": sources, "scale_denominator": 1, "generate_pages": true }),
    );
    assert!(cancelled.is_err());
    assert_eq!(host.document(), before);

    // Save over MCP, reopen in an independent host and redraw identically.
    run(&host, "save_project", json!({}));
    drop(host);
    let reopened = open(&path);
    assert_eq!(reopened.document(), before);
    assert_eq!(views(&reopened), edited);
    let book = id(&section["generation"]["book"]);
    let drawing = drawings::read(&reopened.document(), book).unwrap();
    assert!(
        drawing
            .pages
            .iter()
            .flat_map(|p| &p.views)
            .any(|v| v.id == cut)
    );
}
