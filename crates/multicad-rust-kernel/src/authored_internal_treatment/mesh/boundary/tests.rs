use super::*;
#[test]
fn actual_unsampled_cells_pair_by_entire_exact_source_curve_laws() {
    for (left, starts, entry, tail) in [
        (false, 1, Some(0.4), true),
        (true, 2, Some(0.4), true),
        (false, 1, None, true),
        (true, 2, Some(0.4), false),
    ] {
        let body = super::super::tests::source(left, starts, entry, tail);
        let graph = Graph::raw(body, BoundaryOptions::default(), None).unwrap();
        let mut internal = 0;
        let mut physical = 0;
        let mut rims = [0; 2];
        for curve in &graph.curves {
            assert!(curve.endpoints[0] < curve.endpoints[1]);
            assert!(!curve.exact.is_empty());
            assert!(curve.base.iter().chain(&curve.delta).all(|v| v.valid()));
            let a = curve.uses[0].unwrap();
            if let Some(b) = curve.uses[1] {
                assert_ne!(a.forward, b.forward);
                if a.face == b.face {
                    internal += 1;
                    assert!(curve.station.is_none());
                } else {
                    physical += 1;
                }
            } else {
                rims[match curve.station.unwrap() {
                    0 => 0,
                    2 => 1,
                    _ => panic!("unpaired actual interior support"),
                }] += 1;
            }
        }
        assert!(internal > 0 && physical > 0 && rims.into_iter().all(|n| n > 0));
        assert!(graph.vertices.iter().all(|v| !v.exact.is_empty()
            && v.rzv.iter().all(|b| b.valid())
            && matches!(v.origin, Origin::Parameter {cell,corner} if cell<graph.cells.len() && corner<graph.cells[cell].count)));
        assert!(graph.cells.iter().all(|c| {
            c.count >= 3
                && c.count <= 6
                && c.band < 5
                && c.period >= -i64::from(starts)
                && c.face >= 2
                && c.parameters[..c.count]
                    .iter()
                    .flatten()
                    .all(|i| i.lo.is_finite() && i.hi.is_finite())
                && matches!(c.kind, Kind::Cone | Kind::Wall)
        }));
        let w = &mut Work {
            remaining: 16_000_000,
            used: 0,
            interrupt: None,
        };
        for face in 2..(2 + 4 * starts + u32::from(entry.is_some())) {
            let loops = graph.face_loops(face, w).unwrap();
            assert_eq!(loops.len(), if face == 2 + 4 * starts { 2 } else { 1 });
        }
        assert!(!graph.rim(0, w).unwrap().is_empty());
        assert!(!graph.rim(1, w).unwrap().is_empty());
        eprintln!(
            "sourcecells {left}/{starts}/{entry:?}/{tail}: {}cells,{}nodes,{}curves,{internal}internal,{physical}creases,{rims:?}rims",
            graph.cells.len(),
            graph.vertices.len(),
            graph.curves.len()
        );
    }
}

#[test]
fn actual_closed_tool_source_has_all_cap_tail_face_loops_and_opposite_coedges() {
    for (left, starts, entry, tail) in [
        (false, 1, Some(0.4), true),
        (true, 2, Some(0.4), true),
        (false, 1, None, true),
        (true, 2, Some(0.4), false),
    ] {
        let body = super::super::tests::source(left, starts, entry, tail);
        let mut w = Work {
            remaining: 16_000_000,
            used: 0,
            interrupt: None,
        };
        let graph = Graph::raw(body, BoundaryOptions::default(), None)
            .unwrap()
            .closed(&mut w)
            .unwrap();
        let base = 2 + 4 * starts + u32::from(entry.is_some());
        for face in 0..base + 2 * u32::from(tail) {
            let loops = graph.face_loops(face, &mut w).unwrap();
            let expected = if (entry.is_some() && face == 2 + 4 * starts) || (tail && face >= base)
            {
                2
            } else {
                1
            };
            assert_eq!(loops.len(), expected, "face{face}");
        }
        assert!(graph.curves.iter().all(|c| c.uses[1].is_some()));
        let (vertices, edges) = graph.qualify(&mut w).unwrap();
        eprintln!("closedsource {left}/{starts}/{entry:?}/{tail}: {vertices}vertices,{edges}edges");

        for node in &graph.vertices {
            if let Origin::Tail { source, end } = node.origin {
                assert!((source as usize) < graph.vertices.len() && end < 2);
            }
        }
    }
}
#[test]
fn actual_stock_minus_void_graph_has_closed_annuli_outer_wall_and_genus_one_links() {
    use crate::{Circle, Rational, cylinder::AxialCylinder};
    for (left, starts, entry, tail) in [
        (false, 1, Some(0.4), true),
        (true, 2, Some(0.4), true),
        (false, 1, None, true),
        (true, 2, Some(0.4), false),
    ] {
        let body = super::super::tests::source(left, starts, entry, tail);
        let zero = Rational::new(0, 1).unwrap();
        // Deliberately actual binary64 3.3, not decimal33/10. The Material
        // constructor must retain its exact stock/tool height equality proof.
        let stock = AxialCylinder::new(
            7,
            Circle::new([zero; 2], Rational::new(5, 2).unwrap()).unwrap(),
            zero,
            Rational::new(3_715_469_692_580_659, 1_125_899_906_842_624).unwrap(),
        )
        .unwrap();
        let material = Material::new(stock, body, 0.5).unwrap();
        let mut w = Work {
            remaining: 16_000_000,
            used: 0,
            interrupt: None,
        };
        let graph = Graph::raw(body, BoundaryOptions::default(), None)
            .unwrap()
            .closed(&mut w)
            .unwrap()
            .material(material, &mut w)
            .unwrap();
        let (vertices, edges) = graph.qualify(&mut w).unwrap();
        let published = material
            .boundary_with_interrupt(BoundaryOptions::default(), None)
            .unwrap();
        assert_eq!(published.material(), Some(&material));
        assert_eq!(published.topology_counts(), (vertices, edges));
        assert!(published.retained_bytes() < 200_000);
        for face in 0..published.face_count() {
            assert!(!published.face_loops(face, None).unwrap().is_empty());
        }
        for node in &graph.vertices {
            if let Origin::Stock { end } = node.origin {
                assert!(end < 2);
            }
        }
        eprintln!("actualmaterial {left}/{starts}: {vertices}vertices,{edges}edges");
    }
}
#[test]
fn immutable_boundary_drops_build_tables_and_keeps_original_domains_and_actual_curve_points() {
    let body = super::super::tests::source(true, 2, Some(0.4), true);
    let boundary = body
        .boundary_with_interrupt(BoundaryOptions::default(), None)
        .unwrap();
    assert_eq!(boundary.cavity(), &body);
    assert!(boundary.material().is_none());
    assert_eq!(boundary.topology_counts().0, boundary.vertices().count());
    assert_eq!(boundary.topology_counts().1, boundary.curves().count());
    assert!(boundary.retained_bytes() < 200_000);
    assert!(boundary.work_used() > 0 && boundary.work_used() <= 16_000_000);
    assert!(boundary.graph.point_owners.is_none() && boundary.graph.curve_owners.is_none());
    for face in 0..boundary.face_count() {
        assert!(!boundary.face_loops(face, None).unwrap().is_empty());
    }
    for (index, cell) in boundary.cells().iter().enumerate() {
        assert!(cell.period() >= -2 && cell.band() < 5 && cell.face() < boundary.face_count());
        assert_eq!(
            cell.is_entry_cone(),
            cell.face() == 2 + 4 * body.source_thread().spec().starts
        );
        for corner in 0..cell.points().len() {
            let uv = boundary.cell_parameter_bounds(index, corner).unwrap();
            assert!(uv[1].lower >= 0. && uv[1].upper <= 1.);
            assert!(
                boundary.points()[cell.points()[corner] as usize]
                    .source_cell()
                    .is_some()
            );
        }
    }
    for (id, curve) in boundary.curves() {
        assert!(curve.faces().into_iter().all(|f| f < boundary.face_count()));
        assert!(
            curve
                .affine_rzv_bounds()
                .into_iter()
                .flatten()
                .all(|v| v.lower <= v.upper)
        );
        for t in [0., 0.25, 0.5, 0.75, 1.] {
            let sample = boundary.curve_point(id, t, None).unwrap();
            assert!(sample.error_mm.is_finite() && sample.point_mm.into_iter().all(f64::is_finite));
        }
        for (t, end) in [(0., 0), (1., 1)] {
            let point = &boundary.points()[curve.endpoints()[end] as usize];
            assert!(point.is_vertex());
            let source = point.position_bounds();
            let emitted = boundary.curve_point(id, t, None).unwrap();
            assert!((0..3).all(|i| source[i].lower <= emitted.point_mm[i]
                && emitted.point_mm[i] <= source[i].upper));
        }
    }
}
#[test]
fn declared_lower_budgets_cancellation_and_invalid_queries_never_publish_partial_boundaries() {
    let body = super::super::tests::source(false, 1, Some(0.4), true);
    let original = body;
    let normal = BoundaryOptions::default();
    for limits in [
        BoundaryOptions { cells: 1, ..normal },
        BoundaryOptions {
            vertices: 1,
            ..normal
        },
        BoundaryOptions {
            curves: 1,
            ..normal
        },
        BoundaryOptions { work: 1, ..normal },
        BoundaryOptions { bytes: 1, ..normal },
    ] {
        assert!(matches!(
            body.boundary_with_interrupt(limits, None),
            Err(Error::Limit)
        ));
    }
    for limits in [
        BoundaryOptions { cells: 0, ..normal },
        BoundaryOptions {
            vertices: 8193,
            ..normal
        },
        BoundaryOptions {
            curves: 8193,
            ..normal
        },
        BoundaryOptions {
            work: 16_000_001,
            ..normal
        },
        BoundaryOptions {
            bytes: BYTES + 1,
            ..normal
        },
    ] {
        assert!(matches!(
            body.boundary_with_interrupt(limits, None),
            Err(Error::Invalid)
        ));
    }
    let counter = std::cell::Cell::new(0usize);
    let count = || {
        counter.set(counter.get() + 1);
        false
    };
    body.boundary_with_interrupt(normal, Some(&count)).unwrap();
    let total = counter.get();
    for stop in [0, 1, total / 4, total / 2, total - 2, total - 1] {
        let calls = std::cell::Cell::new(0usize);
        let cancel = || {
            calls.set(calls.get() + 1);
            calls.get() > stop
        };
        assert!(matches!(
            body.boundary_with_interrupt(normal, Some(&cancel)),
            Err(Error::Cancelled)
        ));
    }
    assert_eq!(body, original);
    let small = Cavity::new(
        *body.source_thread(),
        Intent {
            entry_depth_mm: Some(0.001),
            terminal: body.intent().terminal,
        },
    )
    .unwrap();
    assert!(matches!(
        small.boundary_with_interrupt(normal, None),
        Err(Error::Unresolved)
    ));
    let boundary = body.boundary_with_interrupt(normal, None).unwrap();
    assert!(boundary.points().iter().any(|p| !p.is_vertex()));
    assert!(boundary.points().iter().any(|p| p.source_cell().is_none()));
    assert!(matches!(
        boundary.face_loops(boundary.face_count(), None),
        Err(Error::Invalid)
    ));
    assert!(matches!(
        boundary.face_loops(0, Some(&|| true)),
        Err(Error::Cancelled)
    ));
    assert!(matches!(
        boundary.cell_parameter_bounds(usize::MAX, 0),
        Err(Error::Invalid)
    ));
    assert!(matches!(
        boundary.cell_parameter_bounds(0, usize::MAX),
        Err(Error::Invalid)
    ));
    assert!(matches!(
        boundary.curve_point(u32::MAX, 0., None),
        Err(Error::Invalid)
    ));
    let id = boundary.curves().next().unwrap().0;
    for value in [-1., 2., f64::NAN] {
        assert!(matches!(
            boundary.curve_point(id, value, None),
            Err(Error::Invalid)
        ));
    }
    assert!(matches!(
        boundary.curve_point(id, 0.5, Some(&|| true)),
        Err(Error::Cancelled)
    ));
    let calls = std::cell::Cell::new(0usize);
    let cancel = || {
        calls.set(calls.get() + 1);
        calls.get() > 1
    };
    assert!(matches!(
        boundary.curve_point(id, 0.5, Some(&cancel)),
        Err(Error::Cancelled)
    ));
}
#[test]
fn actual_trim_wire_records_refuse_missing_duplicate_reversed_and_foreign_indices() {
    let body = super::super::tests::source(false, 1, Some(0.4), true);
    let mut w = Work {
        remaining: 16_000_000,
        used: 0,
        interrupt: None,
    };
    let graph = Graph::raw(body, BoundaryOptions::default(), None)
        .unwrap()
        .closed(&mut w)
        .unwrap();
    let mut wire = graph.face_loops(0, &mut w).unwrap().remove(0);
    let removed = wire.pop().unwrap();
    assert!(matches!(graph.trace(&wire, &mut w), Err(Error::Unresolved)));
    wire.push(removed);
    wire.push(removed);
    assert!(matches!(graph.trace(&wire, &mut w), Err(Error::Unresolved)));
    wire.pop();
    wire[0].forward = !wire[0].forward;
    assert!(matches!(graph.trace(&wire, &mut w), Err(Error::Unresolved)));
    wire[0].forward = !wire[0].forward;
    wire[0].curve = usize::MAX;
    assert!(matches!(graph.trace(&wire, &mut w), Err(Error::Unresolved)));
    assert!(matches!(graph.trace(&[], &mut w), Err(Error::Unresolved)));
}
#[test]
fn nonbinary_pitch_and_odd_or_eight_starts_preserve_actual_periodic_domains() {
    use crate::authored_flat_thread::{AuthoredThread, SourceOptions};
    let base = super::super::tests::source(false, 1, Some(0.4), true);
    for starts in [3, 8] {
        for left_hand in [false, true] {
            let source = AuthoredThread::new(
                crate::helical_thread::Spec {
                    pitch_mm: 0.35,
                    depth_mm: 0.15,
                    crest_width_mm: 0.08,
                    starts,
                    left_hand,
                    ..base.source_thread().spec()
                },
                SourceOptions {
                    position_limit_mm: 1e-7,
                    normal_angle_limit_rad: 1e-7,
                },
            )
            .unwrap();
            let body = Cavity::new(source, base.intent()).unwrap();
            let boundary = body
                .boundary_with_interrupt(BoundaryOptions::default(), None)
                .unwrap();
            assert_eq!(boundary.face_count(), 2 + 4 * starts + 3);
            assert!(boundary.retained_bytes() < 400_000);
            for face in 0..boundary.face_count() {
                assert!(!boundary.face_loops(face, None).unwrap().is_empty());
            }
            for (index, cell) in boundary.cells().iter().enumerate() {
                for corner in 0..cell.points().len() {
                    let uv = boundary.cell_parameter_bounds(index, corner).unwrap();
                    assert!(uv[1].lower >= 0. && uv[1].upper <= 1.);
                }
            }
        }
    }
}

#[test]
fn original_constraints_and_exact_laws_remain_owned_after_construction_tables_drop() {
    let cavity = super::super::tests::source(false, 2, Some(0.4), true);
    let zero = crate::Rational::new(0, 1).unwrap();
    let stock = crate::cylinder::AxialCylinder::new(
        7,
        crate::Circle::new([zero; 2], crate::Rational::new(5, 2).unwrap()).unwrap(),
        zero,
        crate::Rational::new(3715469692580659, 1125899906842624).unwrap(),
    )
    .unwrap();
    let material = Material::new(stock, cavity, 0.1).unwrap();
    let boundary = material
        .boundary_with_interrupt(BoundaryOptions::default(), None)
        .unwrap();
    let mut tail = 0;
    let mut stock = 0;
    for node in boundary.points() {
        assert!(!node.exact_parameter_key().is_empty());
        match node.role() {
            BoundaryPointRole::Source(key) => {
                assert_eq!(node.source_key, key);
                assert!(key.0.iter().any(|c| matches!(
                    c,
                    BoundaryConstraint::Rail(_) | BoundaryConstraint::Row(_)
                )));
            }
            BoundaryPointRole::Tail { original, end } => {
                tail += 1;
                assert!(end < 2);
                let Origin::Tail { source, .. } = node.origin else {
                    panic!("actual tail role")
                };
                assert_eq!(
                    boundary.points()[source as usize].role(),
                    BoundaryPointRole::Source(original)
                );
            }
            BoundaryPointRole::StockRim { end } => {
                stock += 1;
                assert!(end < 2);
            }
        }
    }
    assert!(tail > 0);
    assert_eq!(stock, 2);
    for (_, curve) in boundary.curves() {
        assert!(!curve.exact_affine_key().is_empty());
    }
}
#[test]
fn smaller_combined_memory_budget_is_checked_before_owner_allocation_or_exact_work() {
    let body = super::super::tests::source(false, 1, Some(0.4), true);
    let mut w = Work {
        remaining: 16_000_000,
        used: 0,
        interrupt: None,
    };
    let required = Owners::requested_bytes(8192).unwrap();
    assert!(matches!(
        Owners::new_reserved(8192, EXACT_SCRATCH, required + EXACT_SCRATCH - 1, &mut w),
        Err(Error::Limit)
    ));
    // Only the entry checkpoint runs: allocation-table work is not charged.
    assert_eq!(w.used, 1);
    let opts = BoundaryOptions {
        bytes: EXACT_SCRATCH - 1,
        ..BoundaryOptions::default()
    };
    assert!(matches!(Graph::new(body, opts, &mut w), Err(Error::Limit)));
    let normal = body
        .boundary_with_interrupt(BoundaryOptions::default(), None)
        .unwrap();
    assert!(normal.graph.temporary == 0 && normal.graph.point_owners.is_none());
    assert!(normal.retained_bytes() < EXACT_SCRATCH);
}

#[test]
fn whole_source_polylines_keep_canonical_ends_and_independent_stock_circle_chords() {
    use crate::{Circle, Rational, cylinder::AxialCylinder};
    for left in [false, true] {
        let cavity = super::super::tests::source(left, 2, Some(0.4), true);
        let zero = Rational::new(0, 1).unwrap();
        let stock = AxialCylinder::new(
            7,
            Circle::new([zero; 2], Rational::new(5, 2).unwrap()).unwrap(),
            zero,
            Rational::new(3715469692580659, 1125899906842624).unwrap(),
        )
        .unwrap();
        let material = Material::new(stock, cavity, 0.1).unwrap();
        let boundary = material
            .boundary_with_interrupt(BoundaryOptions::default(), None)
            .unwrap();
        let mut stock_rims = 0;
        for (id, curve) in boundary.curves() {
            let polyline = boundary.curve_polyline(id, 0.01, 500_000, None).unwrap();
            assert!(polyline.maximum_error_mm <= 0.01);
            assert_eq!(
                polyline.points_mm.first().copied(),
                Some(boundary.curve_point(id, 0., None).unwrap().point_mm)
            );
            assert_eq!(
                polyline.points_mm.last().copied(),
                Some(boundary.curve_point(id, 1., None).unwrap().point_mm)
            );
            if matches!(
                boundary.points()[curve.endpoints()[0] as usize].role(),
                BoundaryPointRole::StockRim { .. }
            ) {
                stock_rims += 1;
                // Independent exact cardinal source points, unrelated to any
                // sampled manufacturing/display mesh or the polyline formula.
                for (t, expected) in [
                    (0., [2.5, 0.]),
                    (0.25, [0., if left { -2.5 } else { 2.5 }]),
                    (0.5, [-2.5, 0.]),
                    (0.75, [0., if left { 2.5 } else { -2.5 }]),
                    (1., [2.5, 0.]),
                ] {
                    let p = boundary.curve_point(id, t, None).unwrap();
                    assert!(
                        (p.point_mm[0] - expected[0]).hypot(p.point_mm[1] - expected[1])
                            <= p.error_mm
                    );
                }
                for pair in polyline.points_mm.windows(2) {
                    let midpoint = [
                        (pair[0][0] + pair[1][0]) * 0.5,
                        (pair[0][1] + pair[1][1]) * 0.5,
                    ];
                    let sagitta = 2.5 - midpoint[0].hypot(midpoint[1]);
                    assert!(sagitta <= polyline.maximum_error_mm);
                }
                assert!(matches!(
                    boundary.curve_polyline(id, 0.01, 2, None),
                    Err(Error::Limit)
                ));
                assert!(matches!(
                    boundary.curve_polyline(id, f64::MIN_POSITIVE, 500_000, None),
                    Err(Error::Limit)
                ));
            }
        }
        assert_eq!(stock_rims, 2);
        let id = boundary.curves().next().unwrap().0;
        for deflection in [0., -1., f64::NAN] {
            assert!(matches!(
                boundary.curve_polyline(id, deflection, 1000, None),
                Err(Error::Invalid)
            ));
        }
        for limit in [0, 1, 500_001] {
            assert!(matches!(
                boundary.curve_polyline(id, 0.1, limit, None),
                Err(Error::Invalid)
            ));
        }
        assert!(matches!(
            boundary.curve_polyline(u32::MAX, 0.1, 1000, None),
            Err(Error::Invalid)
        ));
        assert!(matches!(
            boundary.curve_polyline(id, 0.1, 1000, Some(&|| true)),
            Err(Error::Cancelled)
        ));
        let count = std::cell::Cell::new(0usize);
        let counter = || {
            count.set(count.get() + 1);
            false
        };
        boundary
            .curve_polyline(id, 0.01, 500_000, Some(&counter))
            .unwrap();
        for stop in [1, count.get() / 2, count.get() - 1] {
            let calls = std::cell::Cell::new(0usize);
            let cancel = || {
                calls.set(calls.get() + 1);
                calls.get() > stop
            };
            assert!(matches!(
                boundary.curve_polyline(id, 0.01, 500_000, Some(&cancel)),
                Err(Error::Cancelled)
            ));
        }
    }
}

#[test]
fn actual_closed_incidence_records_reject_duplicate_attachments_and_missing_or_foreign_uses() {
    let source = super::super::tests::source(false, 1, Some(0.4), true);
    let fresh = || Graph::raw(source, BoundaryOptions::default(), None).unwrap();
    let work = || Work {
        remaining: 16_000_000,
        used: 0,
        interrupt: None,
    };
    let mut w = work();
    let mut graph = fresh();
    let rim = graph.rim(0, &mut w).unwrap();
    let edge = rim[0];
    assert_eq!(
        graph.attach(edge.curve, 0, edge.forward, &mut w),
        Err(Error::Unresolved)
    );
    graph.attach(edge.curve, 0, !edge.forward, &mut w).unwrap();
    assert_eq!(
        graph.attach(edge.curve, 0, !edge.forward, &mut w),
        Err(Error::Unresolved)
    );
    assert_eq!(
        graph.attach(usize::MAX, 0, false, &mut w),
        Err(Error::Unresolved)
    );
    let mut graph = fresh().closed(&mut w).unwrap();
    let curve = graph
        .curves
        .iter()
        .position(|c| c.uses[0].unwrap().face != c.uses[1].unwrap().face)
        .unwrap();
    graph.curves[curve].uses[1] = None;
    assert_eq!(graph.qualify(&mut w), Err(Error::Unresolved));
    let mut graph = fresh().closed(&mut w).unwrap();
    graph.curves[curve].uses[1].as_mut().unwrap().face = u32::MAX;
    assert_eq!(graph.qualify(&mut w), Err(Error::Unresolved));
    let mut graph = fresh().closed(&mut w).unwrap();
    graph.curves[curve].uses[1].as_mut().unwrap().forward =
        graph.curves[curve].uses[0].unwrap().forward;
    assert_eq!(graph.qualify(&mut w), Err(Error::Unresolved));
    // These are real source coedges: reversing the outgoing neighbor gives two
    // incoming germs at the same vertex with distinct starting vertices.
    let graph = fresh().closed(&mut w).unwrap();
    let wire = graph.face_loops(0, &mut w).unwrap().remove(0);
    assert!(wire.len() > 2);
    let mut next = wire[1];
    next.forward = !next.forward;
    assert!(matches!(
        graph.trace(&[wire[0], next], &mut w),
        Err(Error::Unresolved)
    ));
}

#[test]
fn actual_source_ownership_refuses_a_third_coedge_even_with_identical_whole_curve_law() {
    let source = super::super::tests::source(false, 1, Some(0.4), true);
    let mut w = Work {
        remaining: 16_000_000,
        used: 0,
        interrupt: None,
    };
    let mut graph = Graph::raw(source, BoundaryOptions::default(), None)
        .unwrap()
        .closed(&mut w)
        .unwrap();
    let curve = &graph.curves[0];
    let q = |i, w: &mut Work<'_>| close::scalar(&curve.exact, i, w).unwrap();
    let base = [q(0, &mut w), q(1, &mut w), q(2, &mut w)];
    let end = std::array::from_fn(|i| base[i].add(&q(3 + i, &mut w), &mut w).unwrap());
    let ids = curve.endpoints;
    let side = curve.uses[0].unwrap();
    let station = curve.station;
    assert_eq!(
        graph.register(ids, [base, end], station, side, &mut w),
        Err(Error::Unresolved)
    );
}

#[test]
fn valid_sources_reject_unclosed_or_mismatched_material_handoffs() {
    use crate::{Circle, Rational, cylinder::AxialCylinder};
    let source = super::super::tests::source(false, 1, Some(0.4), true);
    let zero = Rational::new(0, 1).unwrap();
    let stock = AxialCylinder::new(
        7,
        Circle::new([zero; 2], Rational::new(5, 2).unwrap()).unwrap(),
        zero,
        Rational::new(3715469692580659, 1125899906842624).unwrap(),
    )
    .unwrap();
    let material = Material::new(stock, source, 0.1).unwrap();
    let other = super::super::tests::source(true, 1, Some(0.4), true);
    let other_material = Material::new(stock, other, 0.1).unwrap();
    let mut w = Work {
        remaining: 16_000_000,
        used: 0,
        interrupt: None,
    };
    let fresh = || Graph::raw(source, BoundaryOptions::default(), None).unwrap();
    assert!(matches!(
        fresh().material(material, &mut w),
        Err(Error::Unresolved)
    ));
    assert!(matches!(
        fresh()
            .closed(&mut w)
            .unwrap()
            .material(other_material, &mut w),
        Err(Error::Invalid)
    ));
    let mut graph = fresh()
        .closed(&mut w)
        .unwrap()
        .material(material, &mut w)
        .unwrap();
    // Three disjoint real rim cycles exceed the source face contract's two.
    let mut edges = graph
        .face_loops(0, &mut w)
        .unwrap()
        .into_iter()
        .flatten()
        .collect::<Vec<_>>();
    let upper = graph.face_loops(1, &mut w).unwrap();
    let circle = upper.into_iter().find(|chain| chain.len() == 1).unwrap();
    edges.extend(circle);
    assert!(matches!(
        graph.trace(&edges, &mut w),
        Err(Error::Unresolved)
    ));
    // Reusing an occupied result is never an alternate source admission route.
    graph.material = Some(material);
    assert!(matches!(
        graph.material(material, &mut w),
        Err(Error::Invalid)
    ));
}
#[test]
fn exact_unsampled_domains_agree_with_independent_qz_face_area_and_semantic_ownership() {
    use crate::authored_internal_treatment::{FaceKey, FacePropertyOptions};
    let tau = I::new(
        std::f64::consts::TAU.next_down(),
        std::f64::consts::TAU.next_up(),
    );
    for (left, starts) in [(false, 1), (true, 2), (false, 3)] {
        let body = super::super::tests::source(left, starts, Some(0.4), true);
        let boundary = body
            .boundary_with_interrupt(BoundaryOptions::default(), None)
            .unwrap();
        let mut areas = vec![I::point(0.); body.face_count() as usize];
        let lead = I::point(body.thread.spec().pitch_mm) * I::point(f64::from(starts));
        let bounds = |e: Enclosure| I::new(e.lower, e.upper);
        for (index, cell) in boundary.cells().iter().enumerate() {
            let uv = (0..cell.points().len())
                .map(|corner| {
                    boundary
                        .cell_parameter_bounds(index, corner)
                        .unwrap()
                        .map(bounds)
                })
                .collect::<Vec<_>>();
            let mut signed_area = I::point(0.);
            let mut qrange = I::new(f64::INFINITY, f64::NEG_INFINITY);
            let mut vrange = qrange;
            for (i, p) in uv.iter().enumerate() {
                let next = uv[(i + 1) % uv.len()];
                signed_area = signed_area + p[0] * next[1] - next[0] * p[1];
                qrange = I::new(qrange.lo.min(p[0].lo), qrange.hi.max(p[0].hi));
                vrange = I::new(vrange.lo.min(p[1].lo), vrange.hi.max(p[1].hi));
            }
            let area = signed_area / I::point(2.);
            assert!(area.hi > 0.);
            // Exact clipping supplies a positive domain; an interval shoelace
            // lower bound may conservatively include zero, never erase area.
            let area = I::new(area.lo.max(0.), area.hi);
            let (key, jacobian) = if cell.is_entry_cone() {
                let entry = body.entry.unwrap();
                let k = (I::point(entry.radius) - I::point(body.thread.minimum_radius_mm()))
                    / I::point(entry.depth);
                let z = qrange + lead * vrange;
                let r = I::point(entry.radius) - k * z;
                // Actual clipped cone lies between its mouth radius and the
                // shared thread seam; this bound is independently source-owned.
                let r = I::new(
                    r.lo.max(body.thread.minimum_radius_mm()),
                    r.hi.min(entry.radius),
                );
                (
                    FaceKey::EntryCone,
                    tau * r * (I::point(1.) + k.square()).sqrt(),
                )
            } else {
                let key = body.thread_face_key(cell.period(), cell.band()).unwrap();
                let band = body.thread.bands()[cell.band()];
                let m = (I::point(band.radius_mm[1]) - I::point(band.radius_mm[0]))
                    / (I::point(band.phase_mm[1]) - I::point(band.phase_mm[0]));
                let q0 = I::point(band.phase_mm[0])
                    + I::point(cell.period() as f64) * I::point(body.thread.spec().pitch_mm);
                let r = I::point(band.radius_mm[0]) + m * (qrange - q0);
                (
                    key,
                    ((tau * r).square() * (I::point(1.) + m.square()) + (m * lead).square()).sqrt(),
                )
            };
            assert_eq!(body.face_index(key), Some(cell.face()));
            areas[cell.face() as usize] = areas[cell.face() as usize] + area * jacobian;
        }
        for (face, area) in areas
            .into_iter()
            .enumerate()
            .filter(|(_, area)| area.hi > 0.)
        {
            let key = body.face_key(face as u32).unwrap();
            let qz = body
                .local_face_properties(key, FacePropertyOptions::default())
                .unwrap();
            assert!(
                area.lo <= qz.area_mm2.upper && qz.area_mm2.lower <= area.hi,
                "actual unsampled q/v versus independent q/z {left}/{starts}/{key:?}: {area:?} vs {:?}",
                qz.area_mm2
            );
            assert!(area.hi.is_finite() && qz.area_mm2.lower > 0.);
        }
    }
}

#[test]
fn real_narrow_crest_curves_refuse_an_unrepresentable_endpoint_error_budget() {
    use crate::authored_flat_thread::{AuthoredThread, SourceOptions};
    let base = super::super::tests::source(false, 1, Some(0.4), true);
    let source = AuthoredThread::new(
        crate::helical_thread::Spec {
            crest_width_mm: 1e-6,
            ..base.thread.spec()
        },
        SourceOptions {
            position_limit_mm: 1e-7,
            normal_angle_limit_rad: 1e-7,
        },
    )
    .unwrap();
    let body = Cavity::new(source, base.intent()).unwrap();
    let boundary = body
        .boundary_with_interrupt(BoundaryOptions::default(), None)
        .unwrap();
    let mut refused_precision = false;
    for (id, _) in boundary.curves() {
        let sample = boundary.curve_point(id, 0., None).unwrap();
        if matches!(
            boundary.curve_polyline(id, sample.error_mm * 0.125, 500_000, None),
            Err(Error::Unresolved)
        ) {
            refused_precision = true;
        }
    }
    assert!(refused_precision);
}

#[test]
fn source_boundary_stage_matrix_cancels_before_publishing_any_geometry() {
    use crate::authored_flat_thread::{AuthoredThread, SourceOptions};
    let base = super::super::tests::source(false, 1, Some(0.4), true);
    let thread = AuthoredThread::new(
        crate::helical_thread::Spec {
            length_mm: 1.5,
            ..base.thread.spec()
        },
        SourceOptions {
            position_limit_mm: 1e-7,
            normal_angle_limit_rad: 1e-7,
        },
    )
    .unwrap();
    let body = Cavity::new(thread, base.intent()).unwrap();
    let options = BoundaryOptions {
        cells: 128,
        vertices: 128,
        curves: 256,
        ..BoundaryOptions::default()
    };
    let calls = std::cell::Cell::new(0usize);
    let counter = || {
        calls.set(calls.get() + 1);
        false
    };
    body.boundary_with_interrupt(options, Some(&counter))
        .unwrap();
    let count = calls.get();
    eprintln!("actual complete source boundary cancellation checkpoints: {count}");
    // The exhaustive qualification run checked all 28,159 positions. Routine
    // regression retains allocation/publication windows and sampled body paths.
    let mut stops = std::collections::BTreeSet::new();
    stops.extend((0..count).step_by(count.div_ceil(512)));
    stops.extend(0..count.min(64));
    stops.extend(count.saturating_sub(64)..count);
    // Locate real phase boundaries without adding production stage logging.
    calls.set(0);
    let mut work = Work {
        remaining: options.work,
        used: 0,
        interrupt: Some(&counter),
    };
    let graph = Graph::raw_in(body, options, &mut work).unwrap();
    let raw = calls.get() + 1;
    let graph = graph.closed(&mut work).unwrap();
    let closed = calls.get() + 1;
    graph.qualify(&mut work).unwrap();
    let qualified = calls.get() + 1;
    for end in [raw, closed, qualified] {
        stops.extend(end.saturating_sub(64)..end.saturating_add(64).min(count));
    }
    assert!(stops.len() <= 1024);
    // Repeat the separately recorded exhaustive qualification when explicitly
    // requested; ordinary regression always retains the bounded stage matrix.
    if std::env::var("MULTICAD_BOUNDARY_EXHAUSTIVE_CANCEL").as_deref() == Ok("1") {
        stops.extend(0..count);
    }
    for stop in stops {
        calls.set(0);
        let cancel = || {
            calls.set(calls.get() + 1);
            calls.get() > stop
        };
        assert!(
            matches!(
                body.boundary_with_interrupt(options, Some(&cancel)),
                Err(Error::Cancelled)
            ),
            "checkpoint {stop}/{count}"
        );
    }
    assert_eq!(body.source_thread(), &thread);
    assert!(body.boundary_with_interrupt(options, None).is_ok());
}

#[test]
fn distinct_actual_cap_vertices_cannot_be_merged_into_a_pinched_source_boundary() {
    let source = super::super::tests::source(false, 1, Some(0.4), true);
    let mut w = Work {
        remaining: 16_000_000,
        used: 0,
        interrupt: None,
    };
    let mut graph = Graph::raw(source, BoundaryOptions::default(), None)
        .unwrap()
        .closed(&mut w)
        .unwrap();
    let mouth = graph
        .curves
        .iter()
        .find(|c| c.station == Some(0))
        .unwrap()
        .endpoints[0];
    let far = graph
        .vertices
        .iter()
        .position(|p| matches!(p.origin, Origin::Tail { end: 1, .. }))
        .unwrap() as u32;
    assert_ne!(mouth, far);
    // Tamper only the proof's incidence records: preserve all actual source
    // curve laws, face/coedge identities and opposite orientations. Joining
    // disjoint cap stars must fail the connected vertex-link certificate.
    for curve in &mut graph.curves {
        for endpoint in &mut curve.endpoints {
            if *endpoint == far {
                *endpoint = mouth;
            }
        }
    }
    assert_eq!(graph.qualify(&mut w), Err(Error::Unresolved));
}

#[test]
fn actual_source_rim_receipts_refuse_multiple_cycles_reversed_direction_and_wrong_winding() {
    let source = super::super::tests::source(false, 1, Some(0.4), true);
    let fresh = || Graph::raw(source, BoundaryOptions::default(), None).unwrap();
    let mut w = Work {
        remaining: 16_000_000,
        used: 0,
        interrupt: None,
    };
    let mut graph = fresh();
    for curve in &mut graph.curves {
        if curve.station == Some(2) {
            curve.station = Some(0);
        }
    }
    assert!(matches!(graph.rim(0, &mut w), Err(Error::Unresolved)));
    let mut graph = fresh();
    let rim = graph.rim(0, &mut w).unwrap();
    graph.curves[rim[0].curve].v_direction = graph.curves[rim[0].curve].v_direction.reverse();
    assert!(matches!(graph.rim(0, &mut w), Err(Error::Unresolved)));
    let mut graph = fresh();
    let rim = graph.rim(0, &mut w).unwrap();
    let curve = rim
        .iter()
        .find(|e| graph.curves[e.curve].wrap_jump != 0)
        .unwrap()
        .curve;
    graph.curves[curve].wrap_jump = 0;
    assert!(matches!(graph.rim(0, &mut w), Err(Error::Unresolved)));
}

#[test]
fn actual_shell_receipts_retain_connectivity_genus_and_face_loop_contracts() {
    let source = super::super::tests::source(false, 1, Some(0.4), true);
    let mut w = Work {
        remaining: 16_000_000,
        used: 0,
        interrupt: None,
    };
    let mut graph = Graph::raw(source, BoundaryOptions::default(), None)
        .unwrap()
        .closed(&mut w)
        .unwrap();
    let (vertices, edges) = graph.qualify(&mut w).unwrap();
    let faces = source.face_count();
    let mut chi = 0;
    for face in 0..faces {
        chi += 2 - graph.face_loops(face, &mut w).unwrap().len() as isize;
    }
    let mut reached = vec![true; faces as usize];
    assert_eq!(
        qualify::shell_receipt(false, vertices, edges, chi, &reached),
        Ok(())
    );
    assert_eq!(
        qualify::shell_receipt(true, vertices, edges, chi, &reached),
        Err(Error::Unresolved)
    );
    reached[faces as usize - 1] = false;
    assert_eq!(
        qualify::shell_receipt(false, vertices, edges, chi, &reached),
        Err(Error::Unresolved)
    );
    // A real complete cone mouth cycle moved to the far-cap receipt produces
    // two simple disjoint cap loops. It must not masquerade as one tool disk.
    let cone = source
        .face_index(crate::authored_internal_treatment::FaceKey::EntryCone)
        .unwrap();
    for curve in &mut graph.curves {
        if curve.station == Some(0) {
            for side in curve.uses.iter_mut().flatten() {
                if side.face == cone {
                    side.face = 1;
                }
            }
        }
    }
    assert_eq!(graph.qualify(&mut w), Err(Error::Unresolved));
}

#[test]
fn valid_source_cannot_exceed_a_lower_requested_finite_cell_budget() {
    use crate::authored_flat_thread::{AuthoredThread, SourceOptions};
    let base = super::super::tests::source(false, 1, Some(0.4), true);
    let thread = AuthoredThread::new(
        crate::helical_thread::Spec {
            length_mm: 10.,
            ..base.thread.spec()
        },
        SourceOptions {
            position_limit_mm: 1e-7,
            normal_angle_limit_rad: 1e-7,
        },
    )
    .unwrap();
    let body = Cavity::new(thread, base.intent()).unwrap();
    assert!(matches!(
        body.boundary_with_interrupt(
            BoundaryOptions {
                cells: 32,
                ..BoundaryOptions::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert_eq!(body.full_length_mm(), 10.);
}

#[test]
fn actual_material_handoff_cancels_each_exact_stock_circle_allocation_and_publication_checkpoint() {
    use crate::{Circle, Rational, cylinder::AxialCylinder};
    let source = super::super::tests::source(false, 1, Some(0.4), true);
    let zero = Rational::new(0, 1).unwrap();
    let stock = AxialCylinder::new(
        7,
        Circle::new([zero; 2], Rational::new(5, 2).unwrap()).unwrap(),
        zero,
        Rational::new(3715469692580659, 1125899906842624).unwrap(),
    )
    .unwrap();
    let material = Material::new(stock, source, 0.1).unwrap();
    let fresh = || {
        let mut w = Work {
            remaining: 16_000_000,
            used: 0,
            interrupt: None,
        };
        Graph::raw(source, BoundaryOptions::default(), None)
            .unwrap()
            .closed(&mut w)
            .unwrap()
    };
    let calls = std::cell::Cell::new(0usize);
    let counter = || {
        calls.set(calls.get() + 1);
        false
    };
    fresh()
        .material(
            material,
            &mut Work {
                remaining: 16_000_000,
                used: 0,
                interrupt: Some(&counter),
            },
        )
        .unwrap();
    let count = calls.get();
    for stop in 0..count {
        calls.set(0);
        let cancel = || {
            calls.set(calls.get() + 1);
            calls.get() > stop
        };
        assert!(matches!(
            fresh().material(
                material,
                &mut Work {
                    remaining: 16_000_000,
                    used: 0,
                    interrupt: Some(&cancel)
                }
            ),
            Err(Error::Cancelled)
        ));
    }
}
