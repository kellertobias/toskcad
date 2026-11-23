//! CAD-364: author-defined fixed plan symbols of generic library parts.
//! Plan projections draw the authored sketch at each instance's pose, composed
//! with the definition's explicit alignment; elevations and 3D keep the model.
#[allow(dead_code)]
#[path = "../../multicad-app/examples/symbol_review.rs"]
mod review;
use multicad_core::{
    commands::Session,
    document::{Document, EntityId, Placement},
    drawings::{self, DrawingView, Scale, ViewEvaluation, ViewKind, ViewSource},
    model::{self, PrimitiveParameters, PrimitiveRequest, length},
    references::ReferenceFrame,
    symbols::{self, PlanAlignment, SymbolPlaceRequest, SymbolPlan, SymbolUpdateRequest},
};
use multicad_model::{Compiler, DrawingProviders};
use std::collections::BTreeMap;

type Bounds = [f64; 4];

fn view(entities: Vec<EntityId>, u: [f64; 3], v: [f64; 3], denominator: u32) -> DrawingView {
    DrawingView {
        id: EntityId::new(),
        name: "Plan".into(),
        source: ViewSource::Model { entities },
        frame: ReferenceFrame {
            origin_mm: [0.; 3],
            u,
            v,
        },
        kind: ViewKind::Projection,
        position_mm: [150., 120.],
        scale: Scale {
            numerator: 1,
            denominator,
        },
        hidden_lines: false,
        tangent_edges: false,
        paper_chord_tolerance_mm: 0.01,
        guide: None,
        equipment_categories: None,
    }
}
const TOP: ([f64; 3], [f64; 3]) = ([1., 0., 0.], [0., 1., 0.]);
const FRONT: ([f64; 3], [f64; 3]) = ([1., 0., 0.], [0., 0., 1.]);

/// Per-source bounds of a plan's vectors in model millimeters.
fn evaluate(session: &Session, view: &DrawingView) -> BTreeMap<EntityId, (Bounds, usize)> {
    let native = Compiler::default()
        .compile(session.document(), session.revision())
        .unwrap();
    let ViewEvaluation::Ready { vectors, .. } = native
        .drawing_view(
            session.document(),
            session.revision(),
            view,
            &DrawingProviders::default(),
        )
        .unwrap()
    else {
        panic!("plan view must evaluate")
    };
    let mut out: BTreeMap<EntityId, (Bounds, usize)> = BTreeMap::new();
    for vector in vectors {
        let entry = out.entry(vector.source).or_insert((
            [
                f64::INFINITY,
                f64::INFINITY,
                f64::NEG_INFINITY,
                f64::NEG_INFINITY,
            ],
            0,
        ));
        for p in &vector.preview_mm {
            entry.0 = [
                entry.0[0].min(p[0]),
                entry.0[1].min(p[1]),
                entry.0[2].max(p[0]),
                entry.0[3].max(p[1]),
            ];
        }
        entry.1 += 1;
    }
    out
}
fn near(a: Bounds, b: Bounds) {
    for i in 0..4 {
        assert!((a[i] - b[i]).abs() < 1e-6, "{a:?} != {b:?}");
    }
}
/// Bounds of `local` after a quarter-turn count about Z and a translation.
fn place(local: Bounds, quarter_turns: u32, translation: [f64; 2]) -> Bounds {
    let corners = [
        [local[0], local[1]],
        [local[2], local[1]],
        [local[2], local[3]],
        [local[0], local[3]],
    ];
    let mut b = [
        f64::INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NEG_INFINITY,
    ];
    for mut p in corners {
        for _ in 0..quarter_turns {
            p = [-p[1], p[0]];
        }
        let p = [p[0] + translation[0], p[1] + translation[1]];
        b = [
            b[0].min(p[0]),
            b[1].min(p[1]),
            b[2].max(p[0]),
            b[3].max(p[1]),
        ];
    }
    b
}
fn plan_edit(session: &mut Session, definition: EntityId, plan: Option<SymbolPlan>) {
    let current = symbols::read_definition(session.document(), definition).unwrap();
    let revision =
        symbols::next_revision(session.document(), current.symbol, &current.revision).unwrap();
    let request =
        SymbolUpdateRequest::plan_edit(session.document(), definition, revision, plan).unwrap();
    let prepared =
        symbols::prepare_update(session, session.revision(), review::actor(), &request).unwrap();
    review::checked(session, prepared.prepared).unwrap();
}

#[test]
fn authored_plan_symbols_replace_plan_projection_at_each_pose_and_survive_edits_history_and_files()
{
    let mut f = review::seed(false).unwrap();
    let [a, b] = [f.placed[0].instance, f.placed[1].instance];
    // An ordinary Part without an authored plan stays a geometry projection.
    let box_request = PrimitiveRequest::new(
        "Ordinary crate 333 × 77",
        PrimitiveParameters::Box {
            size: [length(333.), length(77.), length(40.)],
            frame: Default::default(),
        },
    );
    let prepared = model::prepare_primitive(
        &f.session,
        f.session.revision(),
        review::actor(),
        &box_request,
    )
    .unwrap();
    review::checked(&mut f.session, prepared.prepared).unwrap();
    let crate_body = box_request.body;
    let sources = vec![a, b, crate_body];
    let top = view(sources.clone(), TOP.0, TOP.1, 10);

    // Identity alignment: the authored plan is the panel profile, drawn once
    // per instance at its own translation, rotation and exposed size.
    let initial = evaluate(&f.session, &top);
    let plan_sketch = symbols::read_definition(f.session.document(), f.definition)
        .unwrap()
        .plan
        .unwrap()
        .sketch;
    let native = Compiler::default()
        .compile(f.session.document(), f.session.revision())
        .unwrap();
    let mut local = BTreeMap::new();
    for (instance, turns, translation) in [(a, 0, [250., 50.]), (b, 1, [700., 50.])] {
        let leaf = native
            .snapshot()
            .leaf_occurrences
            .values()
            .find(|leaf| leaf.path.root == instance)
            .unwrap();
        let (bounds, count) = initial[&instance];
        near(
            bounds,
            [
                leaf.bounds.min[0],
                leaf.bounds.min[1],
                leaf.bounds.max[0],
                leaf.bounds.max[1],
            ],
        );
        // Only the authored sketch's four edges, not the projected solid.
        assert_eq!(count, 4, "authored linework replaces the projection");
        let mut back = [
            bounds[0] - translation[0],
            bounds[1] - translation[1],
            bounds[2] - translation[0],
            bounds[3] - translation[1],
        ];
        for _ in 0..(4 - turns) % 4 {
            back = place(back, 1, [0., 0.]);
        }
        local.insert(instance, (back, turns, translation));
    }
    near(initial[&crate_body].0, [0., 0., 333., 77.]);
    assert!((local[&a].0[2] - local[&a].0[0] - 150.).abs() < 1e-6);
    assert!((local[&b].0[2] - local[&b].0[0] - 200.).abs() < 1e-6);
    let front = view(sources.clone(), FRONT.0, FRONT.1, 10);
    let elevation = evaluate(&f.session, &front);
    assert!(
        elevation[&a].0[3] - elevation[&a].0[1] > 1.,
        "elevations keep the physical thickness"
    );

    // Edit the alignment: rotate a quarter turn, then shift 300 mm along -Y.
    let alignment = PlanAlignment {
        offset_mm: [0., -300.],
        rotation_deg: 90.,
    };
    let records: Vec<_> = f.placed.iter().map(|p| p.record).collect();
    let before_ids: Vec<_> = records
        .iter()
        .map(|r| symbols::read_instance(f.session.document(), *r).unwrap())
        .collect();
    let before_edit = f.session.document().clone();
    // Invalid alignment and a non-template sketch are refused without change.
    let current = symbols::read_definition(f.session.document(), f.definition).unwrap();
    let bad = SymbolPlan::aligned(
        plan_sketch,
        PlanAlignment {
            offset_mm: [f64::NAN, 0.],
            rotation_deg: 0.,
        },
    );
    assert!(
        SymbolUpdateRequest::plan_edit(f.session.document(), f.definition, "9", Some(bad)).is_err()
    );
    assert!(
        SymbolUpdateRequest::plan_edit(
            f.session.document(),
            f.definition,
            "9",
            Some(SymbolPlan::new(EntityId::new()))
        )
        .is_err()
    );
    assert!(
        SymbolUpdateRequest::plan_edit(
            f.session.document(),
            f.definition,
            "9",
            current.plan.clone()
        )
        .is_err(),
        "an unchanged plan is not a new revision"
    );
    assert_eq!(f.session.document(), &before_edit);
    plan_edit(
        &mut f.session,
        f.definition,
        Some(SymbolPlan::aligned(plan_sketch, alignment)),
    );
    let after_edit = f.session.document().clone();
    let edited = evaluate(&f.session, &top);
    let expected = |instance: EntityId| {
        let (bounds, turns, translation) = local[&instance];
        let aligned = place(bounds, 1, [0., -300.]);
        place(aligned, turns, translation)
    };
    for instance in [a, b] {
        near(edited[&instance].0, expected(instance));
        assert_eq!(edited[&instance].1, 4);
    }
    near(edited[&crate_body].0, [0., 0., 333., 77.]);
    for (record, before) in records.iter().zip(&before_ids) {
        let after = symbols::read_instance(f.session.document(), *record).unwrap();
        assert_eq!(after.instance, before.instance);
        assert_eq!(after.ids, before.ids);
        assert_eq!(after.overrides, before.overrides);
        assert_ne!(after.pin, before.pin);
    }
    let stored = symbols::read_definition(f.session.document(), f.definition).unwrap();
    assert_eq!(stored.plan.unwrap().alignment, alignment);
    assert_eq!(stored.revision, "2");
    // Elevations still project the detailed model, unchanged by the plan.
    let elevation_after = evaluate(&f.session, &front);
    near(elevation_after[&a].0, elevation[&a].0);

    // Undo/Redo restore and reapply the representation with the source.
    f.session
        .undo(f.session.revision(), review::actor())
        .unwrap();
    assert_eq!(f.session.document(), &before_edit);
    near(evaluate(&f.session, &top)[&a].0, initial[&a].0);
    f.session
        .redo(f.session.revision(), review::actor())
        .unwrap();
    assert_eq!(f.session.document(), &after_edit);

    // Save/reopen keeps the embedded definition and its alignment offline.
    let reopened = review::reopen(&review::archive(f.session.document()).unwrap()).unwrap();
    assert_eq!(reopened.document(), f.session.document());
    let reopened_plan = evaluate(&reopened, &top);
    for instance in [a, b] {
        near(reopened_plan[&instance].0, expected(instance));
    }

    // Library placement in another project carries the same plan.
    let definition = symbols::read_definition(f.session.document(), f.definition).unwrap();
    let mut other = Session::new(Document::new("Another project").unwrap()).unwrap();
    let id = symbols::library_placement_definition(other.document(), &definition)
        .unwrap()
        .unwrap_or_default();
    let request = SymbolPlaceRequest::from_definition(
        id,
        &definition,
        "Panel C · half turn",
        review::values(&f.library, 150., 90., 6.),
        Placement {
            translation_mm: [-500., 400., 0.],
            rotation_xyzw: [0., 0., 1., 0.],
        },
    )
    .unwrap();
    let prepared = symbols::prepare_library_place(
        &other,
        other.revision(),
        review::actor(),
        &definition,
        &request,
        None,
    )
    .unwrap();
    review::checked(&mut other, prepared.prepared).unwrap();
    let placed = evaluate(&other, &view(vec![request.instance], TOP.0, TOP.1, 10));
    let (bounds, _, _) = local[&a];
    near(
        placed[&request.instance].0,
        place(place(bounds, 1, [0., -300.]), 2, [-500., 400.]),
    );

    // Removing the plan restores the ordinary projection.
    plan_edit(&mut f.session, f.definition, None);
    let projected = evaluate(&f.session, &top);
    near(projected[&a].0, initial[&a].0);
    assert!(projected[&a].1 >= 4);
    f.session
        .undo(f.session.revision(), review::actor())
        .unwrap();
    assert_eq!(f.session.document(), &after_edit);
}

#[test]
fn exported_plan_pdf_draws_authored_symbols_beside_ordinary_geometry_and_notes() {
    let mut f = review::seed(false).unwrap();
    let plan_sketch = symbols::read_definition(f.session.document(), f.definition)
        .unwrap()
        .plan
        .unwrap()
        .sketch;
    plan_edit(
        &mut f.session,
        f.definition,
        Some(SymbolPlan::aligned(
            plan_sketch,
            PlanAlignment {
                offset_mm: [0., -300.],
                rotation_deg: 90.,
            },
        )),
    );
    let box_request = PrimitiveRequest::new(
        "Ordinary crate 333 × 77",
        PrimitiveParameters::Box {
            size: [length(333.), length(77.), length(40.)],
            frame: Default::default(),
        },
    );
    let prepared = model::prepare_primitive(
        &f.session,
        f.session.revision(),
        review::actor(),
        &box_request,
    )
    .unwrap();
    review::checked(&mut f.session, prepared.prepared).unwrap();
    let plan = view(
        vec![f.placed[0].instance, f.placed[1].instance, box_request.body],
        TOP.0,
        TOP.1,
        10,
    );
    let view_id = plan.id;
    let mut page = drawings::Page::new("1");
    page.views.push(plan);
    page.notes.push(drawings::TextNote {
        id: EntityId::new(),
        position_mm: [30., 30.],
        text: "Authored symbols at 1:10".into(),
    });
    let page_id = page.id;
    let (book, prepared) = drawings::prepare_create(
        &f.session,
        f.session.revision(),
        review::actor(),
        "Plan symbols",
        Some(page),
    )
    .unwrap();
    review::checked(&mut f.session, prepared).unwrap();
    let native = Compiler::default()
        .compile(f.session.document(), f.session.revision())
        .unwrap();
    let evaluation = native
        .drawings(
            f.session.document(),
            f.session.revision(),
            &DrawingProviders::default(),
        )
        .unwrap();
    let ViewEvaluation::Ready { vectors, .. } = &evaluation.views[&view_id] else {
        panic!("linked plan view must evaluate")
    };
    let symbol_vectors = vectors
        .iter()
        .filter(|v| v.source == f.placed[0].instance)
        .count();
    assert_eq!(symbol_vectors, 4);
    assert!(vectors.iter().any(|v| v.source == box_request.body));
    let render = drawings::render::render_page(
        f.session.document(),
        f.session.revision(),
        book,
        page_id,
        &evaluation,
    )
    .unwrap();
    assert!(
        render
            .texts
            .iter()
            .any(|t| t.text == "Authored symbols at 1:10")
    );
    // The rotated panel A plan (90 × 150 mm in model space) is drawn at 1:10.
    let segments: Vec<f64> = render
        .polylines
        .iter()
        .flat_map(|p| {
            p.points.windows(2).map(|w| {
                let d = [w[1][0] - w[0][0], w[1][1] - w[0][1]];
                (d[0] * d[0] + d[1] * d[1]).sqrt()
            })
        })
        .collect();
    for expected in [9., 15., 33.3, 7.7] {
        assert!(
            segments.iter().any(|l| (l - expected).abs() < 1e-6),
            "paper segment {expected} mm missing from {segments:?}"
        );
    }
    let pdf = drawings::pdf::book_pdf(
        f.session.document(),
        f.session.revision(),
        &[(book, page_id)],
        &evaluation,
    )
    .unwrap();
    assert!(pdf.starts_with(b"%PDF"));
    // Every rendered vertex of the plan, symbols included, is stroked in the PDF.
    let stroked = pdf_points(&pdf);
    let height = render.size_mm[1];
    let mut checked = 0;
    for polyline in &render.polylines {
        for p in &polyline.points {
            assert!(
                stroked.iter().any(|q| (q[0] - p[0]).abs() < 0.01
                    && ((q[1] - p[1]).abs() < 0.01 || (q[1] - (height - p[1])).abs() < 0.01)),
                "PDF lacks rendered vertex {p:?}"
            );
            checked += 1;
        }
    }
    assert!(checked >= 16);
}
fn pdf_points(bytes: &[u8]) -> Vec<[f64; 2]> {
    let text = String::from_utf8_lossy(bytes);
    let mut points = vec![];
    let mut rest = &text[..];
    while let Some(begin) = rest.find("stream\n") {
        let stream = &rest[begin + 7..];
        let end = stream.find("endstream").unwrap();
        for line in stream[..end].lines() {
            let tokens: Vec<_> = line.split_whitespace().collect();
            for (i, t) in tokens.iter().enumerate() {
                if (*t == "m" || *t == "l")
                    && i >= 2
                    && let (Ok(x), Ok(y)) =
                        (tokens[i - 2].parse::<f64>(), tokens[i - 1].parse::<f64>())
                {
                    points.push([x / drawings::pdf::PT_PER_MM, y / drawings::pdf::PT_PER_MM]);
                }
            }
        }
        rest = &stream[end + "endstream".len()..];
    }
    points
}
