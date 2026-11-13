use super::*;
#[test]
fn actual_all_metric_cutter_files_bind_every_finite_source_support_and_carrier() {
    let u = [0.36, 0.48, 0.8];
    let v = [-0.8, 0.6, 0.];
    for metric in 0..=5 {
        let (spec, domain) = crate::phillips_drive::test_fixture(metric);
        let source = Cutter::new(spec, domain).unwrap();
        for (origin, u, v) in [
            ([0.; 3], [1., 0., 0.], [0., 1., 0.]),
            ([13., -7., 4.], u, v),
            ([1e6, -2e6, 3e6], u, v),
        ] {
            let output = crate::step::write_phillips_cutter(&source, origin, u, v, 1e-7).unwrap();
            let receipt = verify_phillips_cutter_source_known(&output, &source, origin, u, v, 1e-7)
                .unwrap_or_else(|e| panic!("metric {metric} at {origin:?}: {e:?}"));
            assert_eq!(receipt.vertex_count(), 49);
            assert_eq!(receipt.edge_count(), 76);
            assert_eq!(receipt.face_count(), 29);
            assert!(receipt.certified_span_count() >= 32);
            assert!(receipt.maximum_vertex_error_mm() <= 1e-7 / 16.);
            assert!(receipt.maximum_analytic_support_error_mm() <= 1e-7 / 16.);
            assert!(receipt.maximum_polynomial_error_mm() <= 1e-7 / 16.);
            assert!(receipt.maximum_polynomial_support_error_mm() <= 1e-7 / 8.);
        }
    }
}
fn fixture() -> (Cutter, String) {
    let (s, d) = crate::phillips_drive::test_fixture(0);
    let source = Cutter::new(s, d).unwrap();
    let text =
        crate::step::write_phillips_cutter(&source, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7)
            .unwrap();
    (source, text)
}
fn value(v: &Value<'_>) -> String {
    match v {
        Value::String(s) => format!("'{s}'"),
        Value::Enum(s) => format!(".{s}."),
        Value::Reference(id) => format!("#{id}"),
        Value::Number(n) => format!("{n:.17E}"),
        Value::List(values) => format!(
            "({})",
            values.iter().map(value).collect::<Vec<_>>().join(",")
        ),
        Value::Omitted => "$".into(),
        Value::Derived => "*".into(),
        Value::Typed(e) => format!(
            "{}({})",
            e.name,
            e.args.iter().map(value).collect::<Vec<_>>().join(",")
        ),
    }
}
fn id(text: &str, name: &str) -> u32 {
    part21::parse(text)
        .unwrap()
        .into_iter()
        .find(|(_, e)| e.len() == 1 && e[0].name == name)
        .unwrap()
        .0
}
fn rewrite<'a>(text: &'a str, id: u32, change: impl FnOnce(&mut Entity<'a>)) -> String {
    let mut entities = part21::parse(text).unwrap();
    let e = entities.get_mut(&id).unwrap();
    assert_eq!(e.len(), 1);
    change(&mut e[0]);
    let body = format!(
        "#{id}={}({});",
        e[0].name,
        e[0].args.iter().map(value).collect::<Vec<_>>().join(",")
    );
    text.lines()
        .map(|line| {
            if line.starts_with(&format!("#{id}=")) {
                body.as_str()
            } else {
                line
            }
        })
        .collect::<Vec<_>>()
        .join("\n")
        + "\n"
}
fn refuse(text: &str, source: &Cutter) {
    assert!(
        verify_phillips_cutter_source_known(
            text,
            source,
            [0.; 3],
            [1., 0., 0.],
            [0., 1., 0.],
            1e-7
        )
        .is_err()
    );
}
#[test]
fn source_known_reader_checks_actual_supports_poles_knots_and_finite_senses() {
    let (source, text) = fixture();
    for name in [
        "CYLINDRICAL_SURFACE",
        "TOROIDAL_SURFACE",
        "CONICAL_SURFACE",
        "CIRCLE",
        "ELLIPSE",
    ] {
        let target = id(&text, name);
        let altered = rewrite(&text, target, |e| {
            if let Value::Number(n) = &mut e.args[2] {
                *n += 0.01;
            } else {
                panic!("dimension")
            }
        });
        refuse(&altered, &source);
        let malformed = rewrite(&text, target, |e| {
            e.args.pop();
        });
        refuse(&malformed, &source);
    }
    let cylinder = id(&text, "CYLINDRICAL_SURFACE");
    let entities = part21::parse(&text).unwrap();
    let axis_id = reference(&entities[&cylinder][0].args[1]).unwrap();
    let origin = reference(&entities[&axis_id][0].args[1]).unwrap();
    let normal = reference(&entities[&axis_id][0].args[2]).unwrap();
    refuse(
        &rewrite(&text, origin, |e| {
            e.args[1] = Value::List((0..3).map(|_| Value::Number(1e308)).collect())
        }),
        &source,
    );
    refuse(
        &rewrite(&text, normal, |e| {
            e.args[1] = Value::List((0..3).map(|_| Value::Number(0.)).collect())
        }),
        &source,
    );
    let spline = id(&text, "B_SPLINE_CURVE_WITH_KNOTS");
    let entities = part21::parse(&text).unwrap();
    let poles = refs(&entities[&spline][0].args[2]).unwrap();
    let shell = Shell::new(&entities, 1e-7, true);
    let start_y = shell.point(poles[0]).unwrap()[1];
    refuse(
        &rewrite(&text, poles[1], |e| {
            if let Value::List(p) = &mut e.args[1] {
                p[1] = Value::Number(start_y + 0.01);
            }
        }),
        &source,
    );
    refuse(
        &rewrite(&text, poles[1], |e| {
            if let Value::List(p) = &mut e.args[1] {
                p[2] = Value::Number(0.5);
            }
        }),
        &source,
    );
    refuse(
        &rewrite(&text, spline, |e| e.args[1] = Value::Number(4.)),
        &source,
    );
    refuse(
        &rewrite(&text, spline, |e| {
            if let Value::List(p) = &mut e.args[2] {
                p[0] = Value::Reference(poles[1]);
            }
        }),
        &source,
    );
    refuse(
        &rewrite(&text, spline, |e| {
            if let Value::List(p) = &mut e.args[6] {
                p[0] = Value::Number(3.);
            }
        }),
        &source,
    );
    refuse(
        &rewrite(&text, spline, |e| {
            if let Value::List(p) = &mut e.args[7] {
                p[1] = Value::Number(0.);
            }
        }),
        &source,
    );
    refuse(
        &rewrite(&text, spline, |e| {
            e.args[7] = Value::List((0..4098).map(|_| Value::Number(0.)).collect())
        }),
        &source,
    );
    let intersection = id(&text, "INTERSECTION_CURVE");
    refuse(
        &rewrite(&text, intersection, |e| {
            e.args[3] = Value::Enum("PCURVE_S1")
        }),
        &source,
    );
    refuse(
        &rewrite(&text, intersection, |e| {
            e.args[2] = Value::List((0..2).map(|_| Value::Reference(cylinder)).collect())
        }),
        &source,
    );
    let line = id(&text, "LINE");
    let entities = part21::parse(&text).unwrap();
    let vector = reference(&entities[&line][0].args[2]).unwrap();
    let direction = reference(&entities[&vector][0].args[1]).unwrap();
    refuse(
        &rewrite(&text, vector, |e| e.args[2] = Value::Number(-1.)),
        &source,
    );
    refuse(
        &rewrite(&text, direction, |e| {
            e.args[1] = Value::List((0..3).map(|_| Value::Number(0.)).collect())
        }),
        &source,
    );
    let use_id = id(&text, "ORIENTED_EDGE");
    refuse(
        &rewrite(&text, use_id, |e| {
            e.args[4] = Value::Enum(if sense(&e.args[4]).unwrap() { "F" } else { "T" })
        }),
        &source,
    );
    let face = id(&text, "ADVANCED_FACE");
    refuse(
        &rewrite(&text, face, |e| {
            e.args[3] = Value::Enum(if sense(&e.args[3]).unwrap() { "F" } else { "T" })
        }),
        &source,
    );
    let edge = id(&text, "EDGE_CURVE");
    refuse(
        &rewrite(&text, edge, |e| {
            e.args[4] = Value::Enum(if sense(&e.args[4]).unwrap() { "F" } else { "T" })
        }),
        &source,
    );
}
#[test]
fn bounded_source_known_residual_does_not_claim_exact_literal_incidence() {
    let (source, text) = fixture();
    let base = verify_phillips_cutter_source_known(
        &text,
        &source,
        [0.; 3],
        [1., 0., 0.],
        [0., 1., 0.],
        1e-7,
    )
    .unwrap();
    let changed = rewrite(&text, id(&text, "CYLINDRICAL_SURFACE"), |e| {
        if let Value::Number(r) = &mut e.args[2] {
            *r += 1e-9;
        }
    });
    let receipt = verify_phillips_cutter_source_known(
        &changed,
        &source,
        [0.; 3],
        [1., 0., 0.],
        [0., 1., 0.],
        1e-7,
    )
    .unwrap();
    assert!(
        receipt.maximum_analytic_support_error_mm()
            > base.maximum_analytic_support_error_mm() + 1e-9
    );
    assert!(
        receipt.maximum_polynomial_support_error_mm() > base.maximum_polynomial_support_error_mm()
    );
}
#[test]
fn literal_graph_reordering_keeps_geometry_and_does_not_repair_wrong_incidence() {
    let (source, mut text) = fixture();
    let shell = id(&text, "CLOSED_SHELL");
    text = rewrite(&text, shell, |e| {
        if let Value::List(v) = &mut e.args[1] {
            v.reverse();
        }
    });
    let loops = part21::parse(&text)
        .unwrap()
        .into_iter()
        .filter(|(_, e)| e.len() == 1 && e[0].name == "EDGE_LOOP")
        .map(|(i, _)| i)
        .collect::<Vec<_>>();
    for target in loops {
        text = rewrite(&text, target, |e| {
            if let Value::List(v) = &mut e.args[1] {
                v.rotate_left(1);
            }
        });
    }
    let use_id = id(&text, "ORIENTED_EDGE");
    let entities = part21::parse(&text).unwrap();
    let edge_id = reference(&entities[&use_id][0].args[3]).unwrap();
    let uses = entities
        .iter()
        .filter(|(_, e)| {
            e.len() == 1
                && e[0].name == "ORIENTED_EDGE"
                && reference(&e[0].args[3]).unwrap() == edge_id
        })
        .map(|(&id, _)| id)
        .collect::<Vec<_>>();
    text = rewrite(&text, edge_id, |e| {
        e.args.swap(1, 2);
        e.args[4] = Value::Enum(if sense(&e.args[4]).unwrap() { "F" } else { "T" });
    });
    for target in uses {
        text = rewrite(&text, target, |e| {
            e.args[4] = Value::Enum(if sense(&e.args[4]).unwrap() { "F" } else { "T" })
        });
    }
    assert!(
        verify_phillips_cutter_source_known(
            &text,
            &source,
            [0.; 3],
            [1., 0., 0.],
            [0., 1., 0.],
            1e-7
        )
        .is_ok()
    );
}
#[test]
fn source_known_invalid_precision_cancel_and_wrong_geometry_leave_source_unchanged() {
    let (source, text) = fixture();
    let before = source.bounding_box_mm();
    for tolerance in [0., f64::NAN] {
        assert_eq!(
            verify_phillips_cutter_source_known(
                &text,
                &source,
                [0.; 3],
                [1., 0., 0.],
                [0., 1., 0.],
                tolerance
            ),
            Err(PhillipsReadError::Invalid)
        );
    }
    assert_eq!(
        verify_phillips_cutter_source_known(&text, &source, [0.; 3], [0.; 3], [0.; 3], 1e-7),
        Err(PhillipsReadError::Invalid)
    );
    refuse(
        &text,
        &Cutter::new(
            crate::phillips_drive::test_fixture(1).0,
            crate::phillips_drive::test_fixture(1).1,
        )
        .unwrap(),
    );
    for at in [0, 1, 20, 150, 500] {
        let calls = Cell::new(0);
        assert_eq!(
            verify_phillips_cutter_source_known_with_interrupt(
                &text,
                &source,
                [0.; 3],
                [1., 0., 0.],
                [0., 1., 0.],
                1e-7,
                Some(&|| {
                    let n = calls.get();
                    calls.set(n + 1);
                    n >= at
                })
            ),
            Err(PhillipsReadError::Cancelled)
        );
    }
    assert_eq!(source.bounding_box_mm(), before);
    assert!(
        verify_phillips_cutter_source_known(
            &text,
            &source,
            [0.; 3],
            [1., 0., 0.],
            [0., 1., 0.],
            1e-7
        )
        .is_ok()
    );
}
#[test]
fn malformed_axis_records_curve_lists_and_cell_records_are_never_coerced() {
    let (source, text) = fixture();
    let entities = part21::parse(&text).unwrap();
    let cone = id(&text, "CONICAL_SURFACE");
    refuse(
        &rewrite(&text, cone, |e| e.args[3] = Value::Number(1e-300)),
        &source,
    );
    let circle = id(&text, "CIRCLE");
    refuse(
        &rewrite(&text, circle, |e| e.args[2] = Value::Number(1e300)),
        &source,
    );
    refuse(
        &rewrite(&text, circle, |e| e.args[2] = Value::Number(0.)),
        &source,
    );
    let fake_ellipse = rewrite(&text, circle, |e| {
        e.name = "ELLIPSE";
        let Value::Number(r) = e.args[2] else {
            panic!("radius")
        };
        e.args.push(Value::Number(r));
    });
    refuse(&fake_ellipse, &source);
    let line = id(&text, "LINE");
    let edge = entities
        .iter()
        .find(|(_, e)| {
            e.len() == 1 && e[0].name == "EDGE_CURVE" && reference(&e[0].args[3]).unwrap() == line
        })
        .unwrap()
        .0;
    refuse(
        &rewrite(&text, *edge, |e| e.args[3] = Value::Reference(circle)),
        &source,
    );
    let circle_edge = entities
        .iter()
        .find(|(_, e)| {
            e.len() == 1 && e[0].name == "EDGE_CURVE" && reference(&e[0].args[3]).unwrap() == circle
        })
        .unwrap()
        .0;
    refuse(
        &rewrite(&text, *circle_edge, |e| {
            e.args[4] = Value::Enum(if sense(&e.args[4]).unwrap() { "F" } else { "T" })
        }),
        &source,
    );
    let cylinder = id(&text, "CYLINDRICAL_SURFACE");
    let axis_id = reference(&entities[&cylinder][0].args[1]).unwrap();
    let normal = reference(&entities[&axis_id][0].args[2]).unwrap();
    for target in [axis_id, normal, *edge, id(&text, "VERTEX_POINT")] {
        refuse(
            &rewrite(&text, target, |e| {
                e.args.pop();
            }),
            &source,
        );
    }
    let face = id(&text, "ADVANCED_FACE");
    refuse(
        &rewrite(&text, face, |e| e.args[1] = Value::List(vec![])),
        &source,
    );
    let vector = reference(&entities[&line][0].args[2]).unwrap();
    let direction = reference(&entities[&vector][0].args[1]).unwrap();
    for target in [vector, direction] {
        refuse(
            &rewrite(&text, target, |e| {
                e.args.pop();
            }),
            &source,
        );
    }
    let spline = id(&text, "B_SPLINE_CURVE_WITH_KNOTS");
    for field in [6, 7] {
        refuse(
            &rewrite(&text, spline, |e| e.args[field] = Value::Number(0.)),
            &source,
        );
    }
    refuse(
        &rewrite(&text, spline, |e| {
            if let Value::List(v) = &mut e.args[7] {
                v[1] = Value::Reference(circle);
            }
        }),
        &source,
    );
}
#[test]
fn real_layer_failures_keep_their_typed_error_categories() {
    let (s, d) = crate::phillips_drive::test_fixture(1);
    let mut bad = s;
    bad.mouth_diameter_mm = f64::NAN;
    assert_eq!(
        cutter_error(Cutter::new(bad, d).unwrap_err()),
        ReadError::Invalid
    );
    let mut clipped = d;
    clipped.corner_reach_mm = 0.001;
    let source = Cutter::new(s, clipped).unwrap();
    assert_eq!(
        cutter_error(source.boundary_topology().unwrap_err()),
        ReadError::Unsupported
    );
    let source = Cutter::new(s, d).unwrap();
    let cancelled = crate::step::write_phillips_cutter_with_interrupt(
        &source,
        [0.; 3],
        [1., 0., 0.],
        [0., 1., 0.],
        1e-7,
        Some(&|| true),
    )
    .unwrap_err();
    assert_eq!(source_error(cancelled), ReadError::Limit);
}

#[test]
fn actual_changed_corner_join_roles_remain_a_finite_source_known_boundary() {
    let (mut spec, domain) = crate::phillips_drive::test_fixture(1);
    spec.inner = crate::phillips_drive::Inner::ShallowV {
        angle_deg: 138.,
        e_width_mm: 0.55,
    };
    let source = Cutter::new(spec, domain).unwrap();
    let output =
        crate::step::write_phillips_cutter(&source, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7)
            .unwrap();
    assert!(
        verify_phillips_cutter_source_known(
            &output,
            &source,
            [0.; 3],
            [1., 0., 0.],
            [0., 1., 0.],
            1e-7
        )
        .is_ok()
    );
}
