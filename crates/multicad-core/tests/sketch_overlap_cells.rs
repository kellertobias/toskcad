//! CAD-400: ordinary sketches whose closed shapes cross, touch or share
//! straight edges are split into exact cells; authored curves stay intact.
use multicad_core::{
    constraint_sketch::{Curve, CurveKind, Sketch, source::Source},
    document::EntityId,
    sketch_regions::{
        self, AnalyticQueries, Extraction, Problem, extract_sketch_with_queries_and_interrupt,
    },
};

fn source() -> Source {
    Source {
        u: [1., 0., 0.],
        sketch: Sketch::default(),
    }
}
fn line(s: &mut Source, a: [f64; 2], b: [f64; 2]) -> EntityId {
    let a = s.add_point(a);
    let b = s.add_point(b);
    s.add_line(a, b)
}
fn circle(s: &mut Source, p: [f64; 2], radius_mm: f64) -> EntityId {
    let center = s.add_point(p);
    let id = EntityId::new();
    s.sketch.curves.push(Curve {
        id,
        construction: false,
        kind: CurveKind::Circle { center, radius_mm },
    });
    id
}
fn cells(s: &Source) -> Extraction {
    extract_sketch_with_queries_and_interrupt(
        &s.sketch.points,
        &s.sketch.curves,
        &[],
        &AnalyticQueries,
        None,
    )
    .unwrap()
}
fn areas(result: &Extraction) -> Vec<f64> {
    let mut areas: Vec<_> = result.regions.iter().map(|r| r.area_mm2).collect();
    areas.sort_by(f64::total_cmp);
    areas
}
fn near(a: &[f64], b: &[f64]) {
    assert_eq!(a.len(), b.len(), "{a:?} != {b:?}");
    for (a, b) in a.iter().zip(b) {
        assert!((a - b).abs() < 1e-7, "{a} != {b}");
    }
}
/// Rectangle A (0,0)-(100,60) and B from (50,30) to A's corner (100,0).
fn reported() -> (Source, EntityId) {
    let mut s = source();
    s.add_rectangle([0., 0.], [100., 60.]).unwrap();
    let corner = s
        .sketch
        .points
        .iter()
        .find(|p| p.position_mm == [100., 0.])
        .unwrap()
        .id;
    let ids = [
        s.add_point([50., 30.]),
        s.add_point([100., 30.]),
        corner,
        s.add_point([50., 0.]),
    ];
    for i in 0..4 {
        s.add_line(ids[i], ids[(i + 1) % 4]);
    }
    (s, corner)
}

#[test]
fn reported_overlapping_rectangles_split_into_the_l_and_the_quarter() {
    let (s, _) = reported();
    // The authored-loop extraction refuses it, as before.
    let authored = sketch_regions::extract(&s);
    assert!(authored.regions.is_empty());
    let result = cells(&s);
    near(&areas(&result), &[1500., 4500.]);
    assert!(result.diagnostics.is_empty(), "{:?}", result.diagnostics);
    let intents: Vec<_> = result.regions.iter().map(|r| r.intent()).collect();
    assert_ne!(intents[0], intents[1]);
    for intent in &intents {
        assert!(result.resolve(intent).is_some());
    }
    // Source order and curve direction do not change the selection intents.
    let mut reordered = s.clone();
    reordered.sketch.curves.reverse();
    for curve in &mut reordered.sketch.curves {
        if let CurveKind::Line { start, end } = &mut curve.kind {
            std::mem::swap(start, end);
        }
    }
    let again = cells(&reordered);
    for intent in &intents {
        assert!(again.resolve(intent).is_some());
    }
    // The authored curves are never trimmed.
    assert_eq!(s.sketch.curves.len(), 8);
}

#[test]
fn clean_sketches_keep_their_authored_extraction() {
    let mut s = source();
    s.add_rectangle([0., 0.], [100., 60.]).unwrap();
    circle(&mut s, [50., 30.], 5.);
    assert_eq!(cells(&s), sketch_regions::extract(&s));
}

#[test]
fn crossing_line_partitions_and_open_ends_stay_diagnosed() {
    let mut s = source();
    s.add_rectangle([0., 0.], [20., 10.]).unwrap();
    // Touches the rectangle at (20, 5) and ends in empty space.
    let stray = line(&mut s, [20., 5.], [30., 5.]);
    let divider = line(&mut s, [10., 0.], [10., 10.]);
    let result = cells(&s);
    near(&areas(&result), &[100., 100.]);
    assert!(
        result
            .regions
            .iter()
            .all(|r| r.outer.intent.curves.contains(&divider))
    );
    let open: Vec<_> = result
        .diagnostics
        .iter()
        .filter(|d| d.problem == Problem::OpenBoundary)
        .collect();
    assert_eq!(open.len(), 1, "{:?}", result.diagnostics);
    assert_eq!(open[0].curves, vec![stray]);
    assert_eq!(open[0].at_mm, Some([30., 5.]));
}

#[test]
fn curved_overlaps_and_unsplittable_circles_are_refused_explicitly() {
    let mut s = source();
    circle(&mut s, [0., 0.], 10.);
    circle(&mut s, [0., 0.], 10.);
    let result = cells(&s);
    assert!(result.regions.is_empty());
    assert!(
        result
            .diagnostics
            .iter()
            .any(|d| d.message.contains("could not be split")),
        "{:?}",
        result.diagnostics
    );
    // Analytic queries cannot split a circle; the refusal says so and the
    // unrelated rectangle stays available.
    let mut s = source();
    s.add_rectangle([20., 20.], [30., 30.]).unwrap();
    circle(&mut s, [0., 0.], 10.);
    line(&mut s, [-10., 0.], [10., 0.]);
    let result = cells(&s);
    near(&areas(&result), &[100.]);
    assert!(
        result
            .diagnostics
            .iter()
            .any(|d| d.message.contains("could not be split") && d.message.contains("native")),
        "{:?}",
        result.diagnostics
    );
}

#[test]
fn degenerate_and_unrelated_points_do_not_block_cells() {
    let (mut s, _) = reported();
    let zero = line(&mut s, [500., 0.], [500., 0.]);
    let result = cells(&s);
    near(&areas(&result), &[1500., 4500.]);
    assert!(
        result
            .diagnostics
            .iter()
            .any(|d| d.problem == Problem::DegenerateCurve && d.curves == vec![zero])
    );
}
