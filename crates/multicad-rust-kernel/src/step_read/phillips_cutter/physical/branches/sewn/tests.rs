use super::*;
fn prerequisites(b: &Boundary) -> (FiniteLifts<'_>, CubicSupportBounds<'_>) {
    let endpoints = b
        .endpoint_coupling_bounds(
            EndpointCouplingOptions {
                maximum_distance_mm: 1e-7,
                work_limit: 1_000_000,
            },
            None,
        )
        .unwrap();
    let analytic = b
        .analytic_support_bounds(
            &endpoints,
            AnalyticSupportOptions {
                maximum_distance_mm: 1e-7,
                work_limit: 10_000_000,
            },
            None,
        )
        .unwrap();
    let lifts = analytic
        .finite_edge_lifts(endpoints, LiftOptions::default(), None)
        .unwrap();
    let cubics = b
        .cubic_support_bounds(
            CubicSupportOptions {
                maximum_distance_mm: 1e-7,
                work_limit: 10_000_000,
            },
            None,
        )
        .unwrap();
    (lifts, cubics)
}
#[test]
fn actual_metric_sources_share_original_point_leaves_with_explicit_derived_transfers() {
    for index in 0..6 {
        for posed in [false, true] {
            let text = super::super::tests::text(index, posed);
            let body = read_phillips_body(&text, 1e-7, 1_000_000).unwrap();
            let b = body.cavity();
            let (lifts, cubics) = prerequisites(b);
            let sewn = lifts
                .sewn_edge_sources(&body, &cubics, SewnOptions::default(), None)
                .unwrap();
            assert!(std::ptr::eq(sewn.boundary(), b));
            assert!(std::ptr::eq(sewn.body(), &body));
            assert_eq!(sewn.edges().len(), 76);
            assert!(sewn.work_used() < WORK);
            assert!(sewn.retained_capacity_bytes() > std::mem::size_of_val(&sewn));
            let mut max_position = 0_f64;
            let mut max_derivative = 0_f64;
            for record in sewn.edges() {
                let edge = &b.edges()[record.edge_index()];
                assert_eq!(
                    record.source_endpoint_point_references(),
                    edge.vertices.map(|i| b.vertices()[i].point_entity)
                );
                max_position = max_position.max(record.position_transfer_mm());
                max_derivative = max_derivative.max(record.derivative_transfer_mm());
                assert!(
                    record
                        .raw_owner_distances_mm()
                        .into_iter()
                        .all(|x| (0. ..1e-7).contains(&x))
                );
                let start = sewn
                    .evaluate(record.edge_index(), 0, [0.; 2], None)
                    .unwrap();
                let end = sewn
                    .evaluate(record.edge_index(), record.span_count() - 1, [1.; 2], None)
                    .unwrap();
                assert_eq!(start.original_point_reference, Some(record.endpoints[0]));
                assert_eq!(end.original_point_reference, Some(record.endpoints[1]));
                for span in 0..record.span_count() {
                    let middle = sewn
                        .evaluate(record.edge_index(), span, [0.25, 0.75], None)
                        .unwrap();
                    assert_eq!(middle.original_point_reference, None);
                    assert!(
                        middle
                            .position_bounds_mm
                            .into_iter()
                            .chain(middle.derivative_bounds_mm)
                            .all(|[a, b]| a.is_finite() && b.is_finite() && a <= b)
                    );
                    if span + 1 < record.span_count() {
                        let a = sewn
                            .evaluate(record.edge_index(), span, [1.; 2], None)
                            .unwrap();
                        let c = sewn
                            .evaluate(record.edge_index(), span + 1, [0.; 2], None)
                            .unwrap();
                        assert_eq!(a.original_point_reference, c.original_point_reference);
                        assert_eq!(a.position_bounds_mm, c.position_bounds_mm);
                    }
                }
                match edge.curve {
                    Curve::CubicIntersection { .. } => {
                        assert_eq!(record.position_transfer_mm(), 0.);
                        assert_eq!(record.derivative_transfer_mm(), 0.);
                        assert!(
                            matches!(start.parameterization, SewnParameter::UnitBezierSpan { raw_curve_forward, .. } if raw_curve_forward == edge.same_sense)
                        );
                    }
                    _ => {
                        assert!(record.position_transfer_mm() > 0.);
                        assert_eq!(start.parameterization, SewnParameter::UnitAnalyticEdge);
                    }
                }
            }
            eprintln!(
                "derived edges H{index} posed{posed}: C0={max_position}mm C1={max_derivative}mm/unit, work={}",
                sewn.work_used()
            );
        }
    }
}

fn decimal(token: &str) -> num_rational::BigRational {
    use num_bigint::BigInt;
    use num_rational::BigRational as Q;
    let (mantissa, exponent) = token.split_once('E').unwrap();
    let (integer, fractional) = mantissa.split_once('.').unwrap();
    let numerator: BigInt = format!("{integer}{fractional}").parse().unwrap();
    let power = exponent.parse::<i32>().unwrap() - fractional.len() as i32;
    if power >= 0 {
        Q::from_integer(numerator * BigInt::from(10).pow(power as u32))
    } else {
        Q::new(numerator, BigInt::from(10).pow((-power) as u32))
    }
}
#[test]
fn actual_derived_lines_match_exact_original_decimal_vertex_interpolation() {
    use num_rational::BigRational as Q;
    for posed in [false, true] {
        let text = super::super::tests::text(0, posed);
        let body = read_phillips_body(&text, 1e-7, 1_000_000).unwrap();
        let (lifts, cubics) = prerequisites(body.cavity());
        let sewn = lifts
            .sewn_edge_sources(&body, &cubics, SewnOptions::default(), None)
            .unwrap();
        for edge in sewn.edges() {
            if !matches!(
                body.cavity().edges()[edge.edge_index()].curve,
                Curve::Line { .. }
            ) {
                continue;
            }
            let endpoints = edge.source_endpoint_point_references().map(|entity| {
                let leaf = body.exact_point_leaf(entity, 1_000_000, None).unwrap();
                assert!(std::ptr::eq(leaf.body(), &body));
                assert_eq!(leaf.entity(), entity);
                leaf.coordinate_spellings().map(decimal)
            });
            for (s, q) in [
                (0.25, Q::new(1.into(), 4.into())),
                (0.5, Q::new(1.into(), 2.into())),
                (0.75, Q::new(3.into(), 4.into())),
            ] {
                let actual = sewn.evaluate(edge.edge_index(), 0, [s; 2], None);
                let actual = actual.unwrap();
                for (k, (a, b)) in endpoints[0].iter().zip(&endpoints[1]).enumerate() {
                    let expected = (a * (Q::from_integer(1.into()) - &q)) + b * &q;
                    let derivative = b - a;
                    let p = actual.position_bounds_mm[k].map(|x| Q::from_float(x).unwrap());
                    let d = actual.derivative_bounds_mm[k].map(|x| Q::from_float(x).unwrap());
                    assert!(p[0] <= expected && expected <= p[1]);
                    assert!(d[0] <= derivative && derivative <= d[1]);
                }
            }
        }
    }
}

#[test]
fn derived_edges_retain_the_actual_positive_literal_mouth_support_gap() {
    use crate::step_read::part21::{self, Value};
    use num_rational::BigRational as Q;
    use num_traits::Signed;
    let text = super::super::tests::text(0, false);
    let body = read_phillips_body(&text, 1e-7, 1_000_000).unwrap();
    let b = body.cavity();
    let entities = part21::parse(&text).unwrap();
    let reference = |value: &Value<'_>| {
        if let Value::Reference(id) = value {
            *id
        } else {
            panic!()
        }
    };
    let mouth = b
        .faces()
        .iter()
        .position(|f| f.material_entity.is_none())
        .unwrap();
    let torus = b
        .faces()
        .iter()
        .position(|f| matches!(f.support, Surface::Torus { .. }))
        .unwrap();
    let original_origin = |face: usize| {
        let support = &entities[&b.faces()[face].surface_entity][0];
        let axis = &entities[&reference(&support.args[1])][0];
        let point = body
            .exact_point_leaf(reference(&axis.args[1]), 1_000_000, None)
            .unwrap();
        let direction = &entities[&reference(&axis.args[2])][0];
        let words = part21::numeric_spellings(direction).unwrap();
        assert_eq!(
            words[..2].iter().map(|x| decimal(x)).collect::<Vec<_>>(),
            [Q::from_integer(0.into()), Q::from_integer(0.into())]
        );
        assert_eq!(decimal(words[2]).abs(), Q::from_integer(1.into()));
        point.coordinate_spellings().map(decimal)
    };
    let p = original_origin(mouth);
    let o = original_origin(torus);
    let support = &entities[&b.faces()[torus].surface_entity][0];
    let radii = part21::numeric_spellings(support).unwrap();
    let radius = decimal(radii[1]);
    let gap = (&p[2] - &o[2]).abs() - radius;
    assert_eq!(gap, Q::new(11.into(), 100000000000000000_i64.into()));
    let (lifts, cubics) = prerequisites(b);
    let sewn = lifts
        .sewn_edge_sources(&body, &cubics, SewnOptions::default(), None)
        .unwrap();
    let mut count = 0;
    for record in sewn.edges() {
        let edge = &b.edges()[record.edge_index()];
        if edge.owners.contains(&mouth) && edge.owners.contains(&torus) {
            let owner = edge.owners.iter().position(|&f| f == torus).unwrap();
            assert!(Q::from_float(record.raw_owner_distances_mm()[owner]).unwrap() >= gap);
            assert!(record.position_transfer_mm() > 0.);
            count += 1;
        }
    }
    // One retained torus face owns this arc. Other fragments keep separate
    // face/entity ownership even when their complete support laws coincide.
    assert_eq!(count, 1);
    assert_eq!(sewn.body().literal_input(), text);
}

#[test]
fn actual_equal_coordinates_do_not_replace_the_cubic_endpoint_reference_proof() {
    use crate::step_read::part21;
    let text = super::super::tests::text(0, false);
    let body = read_phillips_body(&text, 1e-7, 1_000_000).unwrap();
    let entities = part21::parse(&text).unwrap();
    let (spline, point) = body
        .cavity()
        .edges()
        .iter()
        .find_map(|e| {
            if let Curve::CubicIntersection {
                spline_entity,
                pole_entities,
                ..
            } = &e.curve
            {
                Some((*spline_entity, pole_entities[0]))
            } else {
                None
            }
        })
        .unwrap();
    let next = entities.keys().last().unwrap() + 1;
    let original_point = text
        .lines()
        .find(|l| l.starts_with(&format!("#{point}=")))
        .unwrap();
    let duplicate = original_point.replacen(&format!("#{point}="), &format!("#{next}="), 1);
    let original_spline = text
        .lines()
        .find(|l| l.starts_with(&format!("#{spline}=")))
        .unwrap();
    let changed_spline = original_spline.replacen(&format!("(#{point},"), &format!("(#{next},"), 1);
    assert_ne!(changed_spline, original_spline);
    let mut changed = text.replacen(original_spline, &changed_spline, 1);
    let end = changed.rfind("ENDSEC;").unwrap();
    changed.insert_str(end, &format!("{duplicate}\n"));
    assert!(read_phillips_body(&changed, 1e-7, 1_000_000).is_err());
}

#[test]
fn source_binding_numeric_limits_and_stage_spanning_cancellation_publish_no_partial_receipt() {
    let text = super::super::tests::text(0, false);
    let body = read_phillips_body(&text, 1e-7, 1_000_000).unwrap();
    let foreign = read_phillips_body(&text, 1e-7, 1_000_000).unwrap();
    let (lifts, cubics) = prerequisites(body.cavity());
    let (_, foreign_cubics) = prerequisites(foreign.cavity());
    let options = SewnOptions::default();
    assert!(matches!(
        lifts.sewn_edge_sources(&foreign, &cubics, options, None),
        Err(PhillipsReadError::Invalid)
    ));
    assert!(matches!(
        lifts.sewn_edge_sources(&body, &foreign_cubics, options, None),
        Err(PhillipsReadError::Invalid)
    ));
    for o in [
        SewnOptions {
            maximum_position_transfer_mm: 0.,
            ..options
        },
        SewnOptions {
            maximum_position_transfer_mm: f64::NAN,
            ..options
        },
        SewnOptions {
            maximum_derivative_transfer_mm: 0.,
            ..options
        },
        SewnOptions {
            maximum_derivative_transfer_mm: f64::INFINITY,
            ..options
        },
        SewnOptions {
            maximum_work: 0,
            ..options
        },
        SewnOptions {
            maximum_work: WORK + 1,
            ..options
        },
        SewnOptions {
            maximum_bytes: 0,
            ..options
        },
        SewnOptions {
            maximum_bytes: BYTES + 1,
            ..options
        },
    ] {
        assert!(matches!(
            lifts.sewn_edge_sources(&body, &cubics, o, None),
            Err(PhillipsReadError::Invalid)
        ));
    }
    for o in [
        SewnOptions {
            maximum_work: 1,
            ..options
        },
        SewnOptions {
            maximum_bytes: 1,
            ..options
        },
    ] {
        assert!(matches!(
            lifts.sewn_edge_sources(&body, &cubics, o, None),
            Err(PhillipsReadError::Limit)
        ));
    }
    for o in [
        SewnOptions {
            maximum_position_transfer_mm: 1e-20,
            ..options
        },
        SewnOptions {
            maximum_derivative_transfer_mm: 1e-20,
            ..options
        },
    ] {
        assert!(matches!(
            lifts.sewn_edge_sources(&body, &cubics, o, None),
            Err(PhillipsReadError::Unsupported)
        ));
    }
    let calls = std::cell::Cell::new(0usize);
    let count = || {
        calls.set(calls.get() + 1);
        false
    };
    let complete = lifts
        .sewn_edge_sources(&body, &cubics, options, Some(&count))
        .unwrap();
    let total = calls.get();
    let mut stops = (0..2048.min(total)).collect::<Vec<_>>();
    stops.extend([total / 4, total / 2, total - 2, total - 1]);
    stops.sort_unstable();
    stops.dedup();
    for stop in stops {
        calls.set(0);
        let cancel = || {
            calls.set(calls.get() + 1);
            calls.get() > stop
        };
        assert!(matches!(
            lifts.sewn_edge_sources(&body, &cubics, options, Some(&cancel)),
            Err(PhillipsReadError::Cancelled)
        ));
    }
    assert_eq!(
        lifts
            .sewn_edge_sources(&body, &cubics, options, None)
            .unwrap()
            .edges()
            .len(),
        complete.edges().len()
    );
    for (edge, span, parameter) in [
        (usize::MAX, 0, [0.; 2]),
        (0, usize::MAX, [0.; 2]),
        (0, 0, [f64::NAN, 1.]),
        (0, 0, [1., 0.]),
        (0, 0, [-1., 0.]),
        (0, 0, [0., 2.]),
    ] {
        assert!(matches!(
            complete.evaluate(edge, span, parameter, None),
            Err(PhillipsReadError::Invalid)
        ));
    }
    assert!(matches!(
        complete.evaluate(0, 0, [0.5; 2], Some(&|| true)),
        Err(PhillipsReadError::Cancelled)
    ));
    calls.set(0);
    let late = || {
        calls.set(calls.get() + 1);
        calls.get() > 1
    };
    assert!(matches!(
        complete.evaluate(0, 0, [0.5; 2], Some(&late)),
        Err(PhillipsReadError::Cancelled)
    ));
    eprintln!(
        "derived source constructor cancellation checkpoints={total}, bounded stage-spanning matrix qualified"
    );
}

#[test]
fn original_cubic_discrete_basis_cannot_be_claimed_from_rounded_real_tokens() {
    let text = super::super::tests::text(0, false);
    let body = read_phillips_body(&text, 1e-7, 1_000_000).unwrap();
    let spline = body
        .cavity()
        .edges()
        .iter()
        .find_map(|e| {
            if let Curve::CubicIntersection { spline_entity, .. } = e.curve {
                Some(spline_entity)
            } else {
                None
            }
        })
        .unwrap();
    let line = text
        .lines()
        .find(|l| l.starts_with(&format!("#{spline}=")))
        .unwrap();
    for changed_line in [
        line.replacen(
            "B_SPLINE_CURVE_WITH_KNOTS('',3,",
            "B_SPLINE_CURVE_WITH_KNOTS('',3.0000000000000001,",
            1,
        ),
        line.replacen(",(4,3,", ",(4,3.0000000000000001,", 1),
    ] {
        assert_ne!(changed_line, line);
        let changed = text.replacen(line, &changed_line, 1);
        // A structural packet may retain unsupported metadata, but a derived
        // manufacturing-source receipt must refuse its rounded basis claim.
        if let Ok(actual) = read_phillips_body(&changed, 1e-7, 1_000_000) {
            let (lifts, cubics) = prerequisites(actual.cavity());
            assert!(matches!(
                lifts.sewn_edge_sources(&actual, &cubics, SewnOptions::default(), None),
                Err(PhillipsReadError::Unsupported)
            ));
        }
    }
}

#[test]
fn physically_reversed_cubic_source_preserves_reference_identity_and_reverses_local_derivative() {
    use num_rational::BigRational as Q;
    let text = super::super::tests::text(0, false);
    let original = read_phillips_body(&text, 1e-7, 1_000_000).unwrap();
    let edge = original
        .cavity()
        .edges()
        .iter()
        .find(|e| matches!(e.curve, Curve::CubicIntersection { .. }))
        .unwrap();
    assert!(edge.same_sense);
    let [a, b] = edge
        .vertices
        .map(|i| original.cavity().vertices()[i].entity);
    let old = text
        .lines()
        .find(|l| l.starts_with(&format!("#{id}=", id = edge.entity)))
        .unwrap();
    let mut changed = text.replacen(
        old,
        &format!(
            "#{}=EDGE_CURVE('',#{b},#{a},#{},.F.);",
            edge.entity, edge.geometry
        ),
        1,
    );
    for line in text
        .lines()
        .filter(|l| l.contains("=ORIENTED_EDGE(") && l.contains(&format!(",#{},", edge.entity)))
    {
        let flipped = if line.ends_with(".T.);") {
            line.replace(".T.);", ".F.);")
        } else {
            line.replace(".F.);", ".T.);")
        };
        changed = changed.replacen(line, &flipped, 1);
    }
    let body = read_phillips_body(&changed, 1e-7, 1_000_000).unwrap();
    let index = body
        .cavity()
        .edges()
        .iter()
        .position(|e| e.entity == edge.entity)
        .unwrap();
    let actual = &body.cavity().edges()[index];
    assert!(!actual.same_sense);
    let (lifts, cubics) = prerequisites(body.cavity());
    let sewn = lifts
        .sewn_edge_sources(&body, &cubics, SewnOptions::default(), None)
        .unwrap();
    let Curve::CubicIntersection {
        pole_entities,
        spans,
        ..
    } = &actual.curve
    else {
        panic!()
    };
    let result = sewn.evaluate(index, 0, [0.; 2], None).unwrap();
    assert_eq!(
        result.parameterization,
        SewnParameter::UnitBezierSpan {
            original_span: spans.len() - 1,
            raw_curve_forward: false
        }
    );
    assert_eq!(
        result.original_point_reference,
        pole_entities.last().copied()
    );
    let last = pole_entities.len() - 1;
    let points = [last, last - 1].map(|i| {
        body.exact_point_leaf(pole_entities[i], 1_000_000, None)
            .unwrap()
            .coordinate_spellings()
            .map(decimal)
    });
    for (k, (a, b)) in points[0].iter().zip(&points[1]).enumerate() {
        let expected = Q::from_integer(3.into()) * (b - a);
        let bounds = result.derivative_bounds_mm[k].map(|x| Q::from_float(x).unwrap());
        assert!(bounds[0] <= expected && expected <= bounds[1]);
    }
}
