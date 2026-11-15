use super::*;
#[test]
fn actual_m2_all_original_finite_curves_have_measured_transfer_and_graph_ownership() {
    degree11::tests::with_head(&mut |head| {
        let field = head
            .common_chart_periodic_thread_repair(Default::default(), None)
            .unwrap();
        let start = std::time::Instant::now();
        let packet = head
            .periodic_thread_finite_curves(&field, Default::default(), None)
            .unwrap();
        assert!(std::ptr::eq(packet.thread(), &field));
        assert_eq!(
            packet.cap_stations_mm(),
            head.body().inferred_visible_interval_mm()
        );
        assert_eq!(packet.curves().len(), 16);
        assert_eq!(
            packet
                .curves()
                .iter()
                .filter(|c| c.kind == PeriodicFiniteKind::Circle)
                .count(),
            1
        );
        assert_eq!(packet.inherited_work_used(), field.total_work_used());
        assert_eq!(
            packet.total_work_used(),
            packet.new_work_used() + field.total_work_used()
        );
        assert!(packet.retained_budget_bytes() < BYTES);
        let mut uses = std::collections::BTreeMap::<u32, Vec<bool>>::new();
        for face in packet.curves().as_chunks::<4>().0 {
            for (i, c) in face.iter().enumerate() {
                assert_eq!(c.role, i);
                assert!(c.face_entity > 0 && c.curve_entity > 0);
                assert_eq!(c.forward_in_face, i < 2);
                assert!(
                    c.distance_mm < 1e-7
                        && c.cell_derivative_mm < 1e-7
                        && c.endpoint_incidence_reserve_mm < 1e-7
                );
                if c.kind == PeriodicFiniteKind::Polynomial {
                    assert_eq!(
                        c.original_curve_controls.len(),
                        c.original_surface_controls.len()
                    );
                    assert!(!c.original_curve_controls.is_empty());
                }
                let ends = if c.forward_in_face {
                    c.vertex_entities
                } else {
                    [c.vertex_entities[1], c.vertex_entities[0]]
                };
                let next = &face[(i + 1) % 4];
                let start = next.vertex_entities[usize::from(!next.forward_in_face)];
                assert_eq!(ends[1], start);
                uses.entry(c.edge_entity)
                    .or_default()
                    .push(c.forward_in_face);
            }
        }
        assert_eq!(
            uses.values()
                .filter(|v| v.len() == 2 && v[0] != v[1])
                .count(),
            4
        );
        eprintln!(
            "FINITE M2 {:?} new{} inherited{} total{} retained{} transfers{:?}",
            start.elapsed(),
            packet.new_work_used(),
            packet.inherited_work_used(),
            packet.total_work_used(),
            packet.retained_budget_bytes(),
            packet
                .curves()
                .iter()
                .map(|c| (
                    c.edge_entity,
                    c.kind,
                    c.distance_mm,
                    c.cell_derivative_mm,
                    c.endpoint_incidence_reserve_mm
                ))
                .collect::<Vec<_>>()
        );
    });
}
#[test]
fn actual_finite_curve_public_frontiers_preserve_caller_total_and_body_identity() {
    degree11::tests::with_head(&mut |head| {
        let field = head
            .common_chart_periodic_thread_repair(Default::default(), None)
            .unwrap();
        let o = PeriodicFiniteOptions::default();
        let good = head.periodic_thread_finite_curves(&field, o, None).unwrap();
        let used = good.new_work_used();
        let total = good.total_work_used();
        let bytes = good.retained_budget_bytes() + FINITE_SCRATCH;
        for bad in [
            PeriodicFiniteOptions {
                maximum_distance_mm: 0.,
                ..o
            },
            PeriodicFiniteOptions {
                maximum_cell_derivative_mm: f64::NAN,
                ..o
            },
            PeriodicFiniteOptions {
                maximum_distance_mm: 1e-6,
                ..o
            },
            PeriodicFiniteOptions {
                maximum_work: WORK + 1,
                ..o
            },
            PeriodicFiniteOptions {
                maximum_total_work: FINITE_TOTAL + 1,
                ..o
            },
            PeriodicFiniteOptions {
                maximum_bytes: BYTES + 1,
                ..o
            },
            PeriodicFiniteOptions {
                maximum_work: 0,
                ..o
            },
            PeriodicFiniteOptions {
                maximum_work: 1,
                ..o
            },
            PeriodicFiniteOptions {
                maximum_work: used - 1,
                ..o
            },
            PeriodicFiniteOptions {
                maximum_total_work: field.total_work_used() - 1,
                ..o
            },
            PeriodicFiniteOptions {
                maximum_total_work: field.total_work_used(),
                ..o
            },
            PeriodicFiniteOptions {
                maximum_total_work: total - 1,
                ..o
            },
            PeriodicFiniteOptions {
                maximum_bytes: bytes - 1,
                ..o
            },
            PeriodicFiniteOptions {
                maximum_distance_mm: 1e-12,
                ..o
            },
            PeriodicFiniteOptions {
                maximum_cell_derivative_mm: 1e-12,
                ..o
            },
        ] {
            assert!(
                head.periodic_thread_finite_curves(&field, bad, None)
                    .is_err()
            );
        }
        let exact = head
            .periodic_thread_finite_curves(
                &field,
                PeriodicFiniteOptions {
                    maximum_work: used,
                    maximum_total_work: total,
                    maximum_bytes: bytes,
                    ..o
                },
                None,
            )
            .unwrap();
        assert_eq!(exact.new_work_used(), used);
        let count = Cell::new(0);
        let observe = || {
            count.set(count.get() + 1);
            false
        };
        head.periodic_thread_finite_curves(&field, o, Some(&observe))
            .unwrap();
        let total = count.get();
        for frontier in [1, 2, total / 2, total] {
            let n = Cell::new(0);
            let cancel = || {
                n.set(n.get() + 1);
                n.get() == frontier
            };
            assert!(matches!(
                head.periodic_thread_finite_curves(&field, o, Some(&cancel)),
                Err(PhillipsReadError::Cancelled)
            ));
        }
        degree11::tests::with_head(&mut |other| {
            assert!(matches!(
                other.periodic_thread_finite_curves(&field, o, None),
                Err(PhillipsReadError::Invalid)
            ));
        });
        assert_eq!(
            head.periodic_thread_finite_curves(&field, o, None)
                .unwrap()
                .cap_stations_mm(),
            good.cap_stations_mm()
        );
    });
}
#[test]
fn actual_file_decimal_curve_alias_retains_its_own_point_owner_and_leaf() {
    let mut changed = None;
    degree11::tests::with_head(&mut |head| {
        let field = head
            .common_chart_periodic_thread_repair(Default::default(), None)
            .unwrap();
        let packet = head
            .periodic_thread_finite_curves(&field, Default::default(), None)
            .unwrap();
        let curve = &packet.curves()[0];
        let old = curve.original_curve_controls[1];
        let body = head.body();
        let leaf = body.exact_point_leaf(old, 100_000, None).unwrap();
        let tokens = leaf.coordinate_spellings();
        let axis = tokens
            .iter()
            .position(|s| s.parse::<f64>().unwrap() != 0.)
            .unwrap();
        let (mantissa, exponent) = tokens[axis].split_once('E').unwrap();
        let replacement = format!("{mantissa}000000000000000000000000001E{exponent}");
        assert_eq!(
            replacement.parse::<f64>().unwrap(),
            tokens[axis].parse::<f64>().unwrap()
        );
        assert_ne!(replacement, tokens[axis]);
        let new = body.entities.keys().max().unwrap() + 1;
        let record = body.entities[&old][0]
            .source
            .unwrap()
            .replacen(tokens[axis], &replacement, 1);
        let line = body
            .literal_input()
            .lines()
            .find(|l| l.starts_with(&format!("#{}=", curve.curve_entity)))
            .unwrap();
        let next = line.replacen(&format!("#{old},"), &format!("#{new},"), 1);
        assert_ne!(line, next);
        let mut text = body.literal_input().replacen(line, &next, 1);
        let end = text.rfind("ENDSEC;").unwrap();
        text.insert_str(end, &format!("#{new}={};\n", record.trim_end_matches(';')));
        changed = Some((text, old, new, axis, replacement));
    });
    let (text, old, new, axis, replacement) = changed.unwrap();
    super::super::super::super::common_mouth::tests::derived_thread_probe::with_head(
        &text,
        &mut |head| {
            let field = head
                .common_chart_periodic_thread_repair(Default::default(), None)
                .unwrap();
            let packet = head
                .periodic_thread_finite_curves(&field, Default::default(), None)
                .unwrap();
            let curve = &packet.curves()[0];
            assert_eq!(curve.original_curve_controls[1], new);
            assert_eq!(curve.original_surface_controls[1], old);
            assert_ne!(new, old);
            assert_eq!(
                head.body()
                    .exact_point_leaf(new, 100_000, None)
                    .unwrap()
                    .coordinate_spellings()[axis],
                replacement
            );
            assert!(curve.distance_mm < 1e-7 && curve.cell_derivative_mm < 1e-7);
        },
    );
}
#[test]
fn original_subnormal_circle_direction_refuses_at_the_earliest_precision_provider() {
    let mut changed = None;
    degree11::tests::with_head(&mut |head| {
        let field = head
            .common_chart_periodic_thread_repair(Default::default(), None)
            .unwrap();
        let packet = head
            .periodic_thread_finite_curves(&field, Default::default(), None)
            .unwrap();
        let circle = packet
            .curves()
            .iter()
            .find(|c| c.kind == PeriodicFiniteKind::Circle)
            .unwrap();
        let body = head.body();
        let shell = Shell::new(&body.entities, body.tolerance, true);
        let a = shell.args(circle.curve_entity, "CIRCLE").unwrap();
        let axis = shell
            .args(reference(&a[1]).unwrap(), "AXIS2_PLACEMENT_3D")
            .unwrap();
        let origin = reference(&axis[1]).unwrap();
        let direction = reference(&axis[3]).unwrap();
        let next = body.entities.keys().max().unwrap() + 1;
        let axis_id = next + 1;
        let old = reference(&a[1]).unwrap();
        let line = body
            .literal_input()
            .lines()
            .find(|l| l.starts_with(&format!("#{}=", circle.curve_entity)))
            .unwrap();
        let new = line.replacen(&format!("#{old},"), &format!("#{axis_id},"), 1);
        assert_ne!(new, line);
        let mut text = body.literal_input().replacen(line, &new, 1);
        let end = text.rfind("ENDSEC;").unwrap();
        // Exactly the same negative Z direction, with an independently owned
        // original decimal scale. The adjacent f64 evaluation box includes zero,
        // so the numeric family cannot certify a uniform normalization.
        text.insert_str(end,&format!("#{next}=DIRECTION('',(0.,0.,-5.0E-324));\n#{axis_id}=AXIS2_PLACEMENT_3D('',#{origin},#{next},#{direction});\n"));
        changed = Some(text);
    });
    let text = changed.unwrap();
    let body = read_phillips_body(&text, 1e-7, 1_000_000).unwrap();
    // The inherited actual outer-support provider already protects this
    // unresolved normalization. Do not fabricate a Head to bypass that guard.
    assert!(matches!(
        body.outer_carrier_bounds(
            OuterCarrierOptions {
                maximum_distance_mm: 1e-7,
                work_limit: 10_000_000
            },
            None
        ),
        Err(PhillipsReadError::Unsupported)
    ));
}
#[test]
fn actual_circle_numeric_qualification_cancellation_uses_its_real_public_callback_frontier() {
    degree11::tests::with_head(&mut |head| {
        let field = head
            .common_chart_periodic_thread_repair(Default::default(), None)
            .unwrap();
        let options = PeriodicFiniteOptions::default();
        let packet = head
            .periodic_thread_finite_curves(&field, options, None)
            .unwrap();
        let body = head.body();
        let shell = Shell::new(&body.entities, body.tolerance, true);
        let ordinal = packet
            .curves()
            .iter()
            .position(|c| c.kind == PeriodicFiniteKind::Circle)
            .unwrap();
        let circle = &packet.curves()[ordinal];
        let a = shell.args(circle.curve_entity, "CIRCLE").unwrap();
        let axis = shell
            .args(reference(&a[1]).unwrap(), "AXIS2_PLACEMENT_3D")
            .unwrap();
        let point_calls = |id| {
            let count = Cell::new(0);
            let observe = || {
                count.set(count.get() + 1);
                false
            };
            body.exact_point_leaf(id, WORK, Some(&observe)).unwrap();
            count.get()
        };
        let numeric_calls = |id| {
            let count = Cell::new(0);
            let observe = || {
                count.set(count.get() + 1);
                false
            };
            body.original_numeric_leaves(id, WORK, Some(&observe))
                .unwrap();
            count.get()
        };
        let p = point_calls(circle.vertex_point_entities[0]) + 1;
        let leaf_calls = point_calls(reference(&axis[1]).unwrap())
            + 1
            + numeric_calls(reference(&axis[2]).unwrap())
            + 1
            + numeric_calls(reference(&axis[3]).unwrap())
            + 1
            + numeric_calls(circle.curve_entity)
            + 1;
        let count = Cell::new(0);
        let observe = || {
            count.set(count.get() + 1);
            false
        };
        field
            .source()
            .bounds(
                ordinal / 4,
                [0.; 2],
                [packet.cap_stations_mm()[0]; 2],
                WORK,
                Some(&observe),
            )
            .unwrap();
        let query_calls = count.get();
        // Calibrate SAME public source helpers, actual original control counts and
        // final complete callback count; no forged Head or production instrumentation.
        let prefix = 1 + 3 * field.transfers().len() + 1;
        let polynomial_calls =
            |c: &PeriodicFiniteCurve| 1 + c.original_curve_controls.len() * (1 + 2 * p) + 4 * p;
        let before = prefix
            + packet.curves()[..ordinal]
                .iter()
                .map(polynomial_calls)
                .sum::<usize>();
        let target_constructor = before + 1 + 1;
        let raw_constructor = target_constructor + leaf_calls + 1;
        let numeric_constructor = raw_constructor + 2 * p + 1;
        let circle_calls = 1 + 1 + leaf_calls + 1 + 2 * p + 1 + 2 * query_calls + 2;
        let expected = prefix
            + packet
                .curves()
                .iter()
                .filter(|c| c.kind == PeriodicFiniteKind::Polynomial)
                .map(polynomial_calls)
                .sum::<usize>()
            + circle_calls
            + 1;
        let count = Cell::new(0);
        let observe = || {
            count.set(count.get() + 1);
            false
        };
        head.periodic_thread_finite_curves(&field, options, Some(&observe))
            .unwrap();
        assert_eq!(count.get(), expected);
        for frontier in [target_constructor, raw_constructor, numeric_constructor] {
            let count = Cell::new(0);
            let cancel = || {
                count.set(count.get() + 1);
                count.get() == frontier
            };
            assert!(matches!(
                head.periodic_thread_finite_curves(&field, options, Some(&cancel)),
                Err(PhillipsReadError::Cancelled)
            ));
            assert_eq!(count.get(), frontier);
        }
    });
}
mod metric_matrix;
