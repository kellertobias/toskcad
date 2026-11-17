use super::super::Value;
use super::*;
use crate::step_read::brep::pcurve::{Proof, uv};

const CYLINDER: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../multicad-app/tests/fixtures/exchange/corpus/occt-cylinder.step"
));
const CONE: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../multicad-app/tests/fixtures/exchange/corpus/occt-cone-frustum.step"
));
const TORUS: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../multicad-app/tests/fixtures/exchange/corpus/occt-torus.step"
));
fn pairs(source: &Source<'_>) -> Vec<(u32, uv::Definition, Surface, Edge)> {
    let mut result = Vec::new();
    let mut budget = Budget::new(512_000_000, None);
    for (&id, parts) in &source.entities {
        let [entity] = parts.as_slice() else { continue };
        if entity.name != "EDGE_CURVE" {
            continue;
        }
        let [_, start, end, curve_id, same_sense] = entity.args.as_slice() else {
            panic!()
        };
        let curve_id = reference(curve_id).unwrap();
        let [_, _, associated, _] = source.simple(curve_id).unwrap().args.as_slice() else {
            panic!()
        };
        let curve = source.curve(curve_id, 0).unwrap();
        if !matches!(
            curve,
            Curve::Line { .. } | Curve::Circle { .. } | Curve::Ellipse { .. }
        ) {
            continue;
        }
        let mut vertices = [reference(start).unwrap(), reference(end).unwrap()];
        if !crate::step_read::sense(same_sense).unwrap() {
            vertices.swap(0, 1);
        }
        let points = vertices.map(|vertex| {
            let [_, point] = source.args(vertex, "VERTEX_POINT").unwrap() else {
                panic!()
            };
            source.point(reference(point).unwrap()).unwrap()
        });
        let edge = Edge::new(curve, vertices.map(|v| v as usize), points, 1e-7).unwrap();
        for pc in crate::step_read::refs(associated).unwrap() {
            let definition = uv::definition(source, pc, &mut budget).unwrap();
            let surface = source.surface(definition.surface).unwrap();
            if matches!(
                surface,
                Surface::Plane(_)
                    | Surface::Cylinder { .. }
                    | Surface::Cone { .. }
                    | Surface::Sphere { .. }
                    | Surface::Torus { .. }
            ) && matches!(
                definition.geometry,
                Curve::Line { .. } | Curve::Circle { .. } | Curve::Ellipse { .. }
            ) {
                result.push((id, definition, surface, edge.clone()));
            }
        }
    }
    result
}
fn curve_id(source: &Source<'_>, edge: u32) -> u32 {
    let [_, _, _, curve, _] = source.args(edge, "EDGE_CURVE").unwrap() else {
        panic!()
    };
    reference(curve).unwrap()
}
#[test]
fn actual_occt_analytic_pairs_keep_original_real_bounds_and_positive_sense() {
    let mut counts = Vec::new();
    for input in [CYLINDER, CONE, TORUS] {
        let source = Source::parse(input, 1e-7).unwrap();
        let pairs = pairs(&source);
        let mut budget = Budget::new(512_000_000, None);
        assert!(!pairs.is_empty());
        for (id, definition, surface, edge) in &pairs {
            let proof =
                retained::certify(surface, &definition.geometry, edge, 1e-7, &mut budget).unwrap();
            assert!(
                matches!(proof, Proof::Certified(_)),
                "actual analytic pair edge #{id}"
            );
            let a = curve_reserve(&source, curve_id(&source, *id), edge, &mut budget).unwrap();
            let b = image_reserve(&source, definition, surface, edge, &mut budget).unwrap();
            assert!(a >= 0. && b >= 0. && a + b < 1e-7, "edge #{id}: {a}+{b}");
            sense(
                &source,
                curve_id(&source, *id),
                definition,
                surface,
                edge,
                &mut budget,
            )
            .unwrap();
        }
        counts.push(pairs.len());
    }
    assert_eq!(counts, [6, 6, 4]);
}
#[test]
fn original_nonzero_underflow_cannot_turn_diagonal_uv_into_an_isocurve() {
    let source = Source::parse(CYLINDER, 1e-7).unwrap();
    let (_, definition, _, _) = pairs(&source)
        .into_iter()
        .find(|(_, d, s, _)| {
            matches!(s, Surface::Cylinder { .. })
                && matches!(d.geometry,Curve::Line{direction:[_,v,_],..} if v==0.)
        })
        .unwrap();
    let [_, _, vector_id] = source.args(definition.curve_id, "LINE").unwrap() else {
        panic!()
    };
    let [_, direction, _] = source
        .args(reference(vector_id).unwrap(), "VECTOR")
        .unwrap()
    else {
        panic!()
    };
    let direction = reference(direction).unwrap();
    let entity = source.simple(direction).unwrap();
    let [name, values] = entity.args.as_slice() else {
        panic!()
    };
    assert!(matches!(name, Value::String(_)));
    let values = list(values).unwrap();
    let Value::Number(x) = values[0] else {
        panic!()
    };
    let replacement = format!("DIRECTION('',({x:.17E},1E-1000))");
    let original = entity.source.unwrap();
    let changed = CYLINDER.replacen(original, &replacement, 1);
    assert_ne!(changed, CYLINDER);
    let source = Source::parse(&changed, 1e-7).unwrap();
    let (_, definition, surface, edge) = pairs(&source)
        .into_iter()
        .find(|(_, d, s, _)| {
            matches!(s, Surface::Cylinder { .. })
                && matches!(d.geometry,Curve::Line{direction:[_,v,_],..} if v==0.)
        })
        .unwrap();
    let mut budget = Budget::new(512_000_000, None);
    assert!(matches!(
        retained::certify(&surface, &definition.geometry, &edge, 1e-7, &mut budget).unwrap(),
        Proof::Certified(_)
    ));
    assert_eq!(
        image_reserve(&source, &definition, &surface, &edge, &mut budget),
        Err(ReadError::Unsupported)
    );
}
#[test]
fn line_lower_transfer_covers_raw_closest_parameters_as_well_as_candidate_map() {
    let reference = Image::Line {
        origin: cached([0., 0., 0.]),
        direction: cached([1., 0., 0.]),
    };
    // A physical point well away from the reference line can have its RAW
    // closest parameter outside the retained candidate interval [-.1,.1].
    let raw = Image::Line {
        origin: cached([0., 0., 0.]),
        direction: cached([1., 1e-3, 0.]),
    };
    let edge = Edge {
        curve: Curve::Line {
            origin: [0., 1000., 0.],
            direction: [1., 0., 0.],
        },
        range: [-0.1, 0.1],
        vertices: [0, 1],
    };
    let mut budget = Budget::new(512_000_000, None);
    let error = line_image_reserve(raw, reference, &edge, &mut budget).unwrap();
    assert!(
        error > 0.99,
        "a candidate-only reserve would be merely .0001: {error}"
    );
    let reverse = Image::Line {
        origin: cached([0., 0., 1.]),
        direction: cached([-1., 0., 0.]),
    };
    let flat = Edge {
        curve: Curve::Line {
            origin: [0.; 3],
            direction: [1., 0., 0.],
        },
        range: [-1., 1.],
        vertices: [0, 1],
    };
    assert!(
        line_image_reserve(reverse, reverse, &flat, &mut budget).unwrap() < 1e-7,
        "a real shifted-line discrepancy reserve must not require same sense"
    );
    let tiny = Image::Line {
        origin: cached([0.; 3]),
        direction: [I::new(-1e-10, 1e-10), I::point(0.), I::point(0.)],
    };
    assert_eq!(
        line_image_reserve(tiny, reference, &edge, &mut budget),
        Err(ReadError::Unsupported)
    );
}
#[test]
fn large_raw_angle_interval_encloses_the_true_extremum_without_quotient_search() {
    let center = 100_000_000_000f64 * std::f64::consts::TAU;
    let interval = I::new(center.next_down(), center.next_up());
    let cosine = trig(interval, true).unwrap();
    assert_eq!(cosine.hi, 1.);
    for i in 0..=32 {
        let t = interval.lo + (interval.hi - interval.lo) * f64::from(i) / 32.;
        assert!(cosine.lo <= t.cos() && t.cos() <= cosine.hi);
    }
}
#[test]
fn actual_raw_analytic_reserve_work_cancel_and_retry_are_bounded() {
    let source = Source::parse(CYLINDER, 1e-7).unwrap();
    let (_, definition, surface, edge) = pairs(&source)
        .into_iter()
        .find(|(_, d, s, _)| {
            matches!(s, Surface::Cylinder { .. })
                && matches!(d.geometry,Curve::Line{direction:[_,v,_],..} if v==0.)
        })
        .unwrap();
    let initial = 512_000_000;
    let mut budget = Budget::new(initial, None);
    let value = image_reserve(&source, &definition, &surface, &edge, &mut budget).unwrap();
    let used = initial - budget.remaining;
    let mut exact = Budget::new(used, None);
    assert_eq!(
        image_reserve(&source, &definition, &surface, &edge, &mut exact).unwrap(),
        value
    );
    let mut limited = Budget::new(used - 1, None);
    assert_eq!(
        image_reserve(&source, &definition, &surface, &edge, &mut limited),
        Err(ReadError::Limit)
    );
    let polls = std::cell::Cell::new(0);
    let stop = || {
        polls.set(polls.get() + 1);
        polls.get() == 12
    };
    let mut cancelled = Budget::new(initial, Some(&stop));
    assert_eq!(
        image_reserve(&source, &definition, &surface, &edge, &mut cancelled),
        Err(ReadError::Limit)
    );
    assert!(cancelled.cancelled.get());
    let mut retry = Budget::new(initial, None);
    assert_eq!(
        image_reserve(&source, &definition, &surface, &edge, &mut retry).unwrap(),
        value
    );
    assert_eq!(retry.remaining, budget.remaining);
}

#[test]
fn original_decimal_cardinal_plane_translation_cancels_before_enclosure() {
    // Authentic additional geometric records in a real source file. This
    // exercises the source-bound point operator, not material ownership.
    let extra = "\n#900001=CARTESIAN_POINT('',(1000000000000.,0.,5.));\n#900002=AXIS2_PLACEMENT_3D('',#900001,#900004,#900005);\n#900003=PLANE('',#900002);\n#900004=DIRECTION('',(0.,0.,1.));\n#900005=DIRECTION('',(1.,0.,0.));\n#900006=CARTESIAN_POINT('',(-999999999996.85840734641020676154,2.71828182845904523536));\n";
    let end = CYLINDER.rfind("ENDSEC;").unwrap();
    let text = format!("{}{}{}", &CYLINDER[..end], extra, &CYLINDER[end..]);
    let source = Source::parse(&text, 1e-7).unwrap();
    let mut budget = Budget::new(512_000_000, None);
    let p = plane_point(&source, 900003, 900006, &mut budget).unwrap();
    let denominator = num_bigint::BigInt::from(10).pow(20);
    for (bound, numerator) in p.into_iter().zip([
        "314159265358979323846",
        "271828182845904523536",
        "500000000000000000000",
    ]) {
        let exact = num_rational::BigRational::new(
            num_bigint::BigInt::parse_bytes(numerator.as_bytes(), 10).unwrap(),
            denominator.clone(),
        );
        let lo = num_rational::BigRational::from_float(bound.lo).unwrap();
        let hi = num_rational::BigRational::from_float(bound.hi).unwrap();
        assert!(lo <= exact && exact <= hi);
        assert!(
            bound.width() < 1e-12,
            "lost original affine cancellation: {bound:?}"
        );
    }
}

#[test]
fn actual_plane_packet_reuse_matches_independent_original_queries_with_less_work() {
    let mut checked = 0;
    for text in [CYLINDER, CONE, TORUS] {
        let source = Source::parse(text, 1e-7).unwrap();
        for (edge_id, definition, surface, edge) in pairs(&source) {
            if !matches!(surface, Surface::Plane(_)) {
                continue;
            }
            let id = curve_id(&source, edge_id);
            let mut independent = Budget::new(512_000_000, None);
            let image = actual_image(&source, &definition, &surface, &mut independent).unwrap();
            let reserve = curve_reserve(&source, id, &edge, &mut independent).unwrap();
            let Proof::Certified(upper) =
                retained::certify_image(image, &edge, 1e-7, &mut independent).unwrap()
            else {
                panic!("actual plane association lacks full-field agreement")
            };
            sense(&source, id, &definition, &surface, &edge, &mut independent).unwrap();
            let expected = (I::point(upper) + I::point(reserve)).hi;
            let mut reused = Budget::new(512_000_000, None);
            let actual = super::super::plane_analytic_certificate(
                &source,
                id,
                &definition,
                &surface,
                &edge,
                1e-7,
                &mut reused,
            )
            .unwrap();
            assert!(matches!(actual,Some(Proof::Certified(value))if value==expected));
            assert!(reused.used() < independent.used());
            checked += 1;
        }
    }
    assert!(checked >= 4);
}

#[test]
fn actual_reused_plane_packets_keep_work_frontiers_cancellation_and_retry() {
    let source = Source::parse(CYLINDER, 1e-7).unwrap();
    let (edge_id, definition, surface, edge) = pairs(&source)
        .into_iter()
        .find(|(_, _, s, _)| matches!(s, Surface::Plane(_)))
        .unwrap();
    let id = curve_id(&source, edge_id);
    let run = |b: &mut Budget<'_>| {
        super::super::plane_analytic_certificate(&source, id, &definition, &surface, &edge, 1e-7, b)
    };
    let mut full = Budget::new(512_000_000, None);
    let Some(Proof::Certified(first)) = run(&mut full).unwrap() else {
        panic!()
    };
    let used = full.used();
    let mut exact = Budget::new(used, None);
    assert!(matches!(run(&mut exact).unwrap(),Some(Proof::Certified(v))if v==first));
    assert_eq!(exact.remaining, 0);
    assert!(matches!(
        run(&mut Budget::new(used - 1, None)),
        Err(ReadError::Limit)
    ));
    let polls = std::cell::Cell::new(0);
    let stop = || {
        polls.set(polls.get() + 1);
        polls.get() == 20
    };
    let mut cancelled = Budget::new(512_000_000, Some(&stop));
    assert!(matches!(run(&mut cancelled), Err(ReadError::Limit)));
    assert!(cancelled.cancelled.get());
    let mut retry = Budget::new(512_000_000, None);
    assert!(matches!(run(&mut retry).unwrap(),Some(Proof::Certified(v))if v==first));
    assert_eq!(retry.used(), used);
}
