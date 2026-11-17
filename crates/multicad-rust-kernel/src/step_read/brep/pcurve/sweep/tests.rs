use super::*;
use crate::brep::{Frame, NurbsCurve};
use crate::step_read::brep::Source;
use crate::step_read::{ReadError, reference};
const TOL: f64 = 1e-7;
const INPUT: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/src/step_read/tests/fixtures/occt-spline-prism.step"
));
fn line(origin: [f64; 3], direction: [f64; 3]) -> Curve {
    Curve::Line { origin, direction }
}
fn actual(input: &str, pc: u32, edge_id: u32) -> (Surface, Curve, Edge) {
    let source = Source::parse(input, TOL).unwrap();
    let mut prep = Budget::new(512_000_000, None);
    let def = super::super::uv::definition(&source, pc, &mut prep).unwrap();
    let args = source.args(edge_id, "EDGE_CURVE").unwrap();
    let vertices = [&args[1], &args[2]].map(|v| {
        let args = source.args(reference(v).unwrap(), "VERTEX_POINT").unwrap();
        source.point(reference(&args[1]).unwrap()).unwrap()
    });
    let curve = source.curve(reference(&args[3]).unwrap(), 0).unwrap();
    let edge = Edge::new(curve, [0, 1], vertices, TOL).unwrap();
    (source.surface(def.surface).unwrap(), def.geometry, edge)
}
fn certified(surface: &Surface, uv: &Curve, edge: &Edge) -> usize {
    let mut budget = Budget::new(512_000_000, None);
    match certify(surface, uv, edge, TOL, &mut budget).unwrap() {
        Proof::Certified(upper) => assert!(upper.is_finite() && (0. ..=TOL).contains(&upper)),
        Proof::Disagreement(lower) => panic!("unexpected physical disagreement {lower}"),
        Proof::Unsupported => panic!("actual isocurve unproved"),
    }
    budget.used()
}
fn spline(points: Vec<[f64; 3]>, knots: [f64; 2]) -> Curve {
    Curve::nurbs(
        NurbsCurve::new(3, points, None, knots.to_vec(), vec![4, 4]).unwrap(),
        TOL,
    )
}
fn finite(curve: Curve) -> Edge {
    let points = match &curve {
        Curve::Nurbs { curve, .. } => [curve.points()[0], *curve.points().last().unwrap()],
        _ => [curve.evaluate(0.).0, curve.evaluate(1.).0],
    };
    Edge::new(curve, [0, 1], points, TOL).unwrap()
}
#[test]
fn actual_occt_extrusion_all_four_unchanged_original_isocurves_are_certified_after_parser_drop() {
    for (pc, edge) in [(43, 21), (102, 87), (147, 140), (171, 164)] {
        let (surface, uv, edge) = actual(INPUT, pc, edge);
        certified(&surface, &uv, &edge);
    }
    // The whole decoded translated cubic changes, not just its endpoints.
    let changed = INPUT.replace(
        "#150 = CARTESIAN_POINT('',(0.,-6.));",
        "#150 = CARTESIAN_POINT('',(0.,-5.9));",
    );
    assert_ne!(changed, INPUT);
    let (surface, uv, edge) = actual(&changed, 147, 140);
    let mut budget = Budget::new(512_000_000, None);
    assert!(!matches!(
        certify(&surface, &uv, &edge, TOL, &mut budget).unwrap(),
        Proof::Certified(_)
    ));
    let (surface, uv, edge) = actual(INPUT, 147, 140);
    certified(&surface, &uv, &edge);
}
#[test]
fn extrusion_generator_parameter_maps_and_line_scales_do_not_need_to_match() {
    let (surface, _, original) = actual(INPUT, 171, 164);
    let Surface::Extrusion {
        curve: Curve::Nurbs { curve, .. },
        ..
    } = &surface
    else {
        panic!()
    };
    let mut poles = curve.points().to_vec();
    poles.reverse();
    let opposite = finite(spline(poles, [3., 9.]));
    certified(&surface, &line([14., 0., 0.], [-7., 0., 0.]), &opposite);
    certified(&surface, &line([-31., 0., 0.], [3., 0., 0.]), &original);
    let mut opposite_budget = Budget::new(512_000_000, None);
    assert!(matches!(
        certify(
            &surface,
            &line([0.; 3], [-1., 0., 0.]),
            &original,
            TOL,
            &mut opposite_budget
        )
        .unwrap(),
        Proof::Unsupported
    ));
    let (_, uv, _) = actual(INPUT, 43, 21);
    let edge = finite(line([0., 0., 12.], [0., 0., 7.]));
    certified(&surface, &uv, &edge);
}
#[test]
fn revolution_original_generators_and_full_periodic_isocurves_have_whole_support_proofs() {
    let poles = vec![[3., 0., 0.], [5., 0., 1.], [5., 0., 3.], [3., 0., 4.]];
    let generator = spline(poles.clone(), [0., 1.]);
    let surface = Surface::Revolution {
        curve: generator,
        origin: [0.; 3],
        axis: [0., 0., 1.],
    };
    let rotated = finite(spline(
        poles.iter().rev().map(|p| [0., p[0], p[2]]).collect(),
        [2., 8.],
    ));
    certified(
        &surface,
        &line([std::f64::consts::FRAC_PI_2, 20., 0.], [0., -3., 0.]),
        &rotated,
    );
    let frame = Frame {
        origin: [0.; 3],
        x: [0., 1., 0.],
        y: [1., 0., 0.],
        z: [0., 0., -1.],
    };
    let circle = Curve::Circle { frame, radius: 3. };
    let edge = Edge::new(circle, [0, 0], [[0., 3., 0.]; 2], TOL).unwrap();
    certified(&surface, &line([21., 0., 0.], [-4., 0., 0.]), &edge);
    let bad = finite(spline(
        poles.iter().map(|p| [0., p[0] + 0.01, p[2]]).collect(),
        [2., 8.],
    ));
    let mut budget = Budget::new(512_000_000, None);
    assert!(!matches!(
        certify(
            &surface,
            &line([std::f64::consts::FRAC_PI_2, 0., 0.], [0., 1., 0.]),
            &bad,
            TOL,
            &mut budget
        )
        .unwrap(),
        Proof::Certified(_)
    ));
}
#[test]
fn bounded_generator_domain_and_non_isocurve_inputs_are_not_repaired() {
    let (surface, _, edge) = actual(INPUT, 43, 21);
    let mut budget = Budget::new(512_000_000, None);
    assert!(matches!(
        certify(
            &surface,
            &line([2., 0., 0.], [0., 1., 0.]),
            &edge,
            TOL,
            &mut budget
        ),
        Err(ReadError::Unsupported) | Ok(Proof::Unsupported)
    ));
    for uv in [
        line([0.; 3], [1., 1., 0.]),
        line([0.; 3], [0.; 3]),
        line([0., 0., 1.], [1., 0., 0.]),
    ] {
        let mut budget = Budget::new(512_000_000, None);
        assert!(matches!(
            certify(&surface, &uv, &edge, TOL, &mut budget).unwrap(),
            Proof::Unsupported
        ));
    }
}
#[test]
fn actual_curved_isocurve_work_frontier_and_observed_cancellation_preserve_retry() {
    use std::cell::Cell;
    let (surface, uv, edge) = actual(INPUT, 147, 140);
    let calls = Cell::new(0);
    let callback = || {
        calls.set(calls.get() + 1);
        false
    };
    let mut first = Budget::new(512_000_000, Some(&callback));
    assert!(matches!(
        certify(&surface, &uv, &edge, TOL, &mut first).unwrap(),
        Proof::Certified(_)
    ));
    let used = first.used();
    let polls = calls.get();
    assert!(used > 2048 && polls > 4);
    let mut below = Budget::new(used - 1, None);
    assert!(matches!(
        certify(&surface, &uv, &edge, TOL, &mut below),
        Err(ReadError::Limit)
    ));
    let mut exact = Budget::new(used, None);
    assert!(matches!(
        certify(&surface, &uv, &edge, TOL, &mut exact).unwrap(),
        Proof::Certified(_)
    ));
    for stop in [1, polls / 2, polls] {
        let seen = Cell::new(0);
        let callback = || {
            seen.set(seen.get() + 1);
            seen.get() >= stop
        };
        let mut interrupted = Budget::new(512_000_000, Some(&callback));
        assert!(matches!(
            certify(&surface, &uv, &edge, TOL, &mut interrupted),
            Err(ReadError::Limit)
        ));
        assert!(seen.get() >= stop);
        assert_eq!(certified(&surface, &uv, &edge), used);
    }
}
