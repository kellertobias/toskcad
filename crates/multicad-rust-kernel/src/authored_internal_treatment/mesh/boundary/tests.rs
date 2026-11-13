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
        assert!(graph.cells.iter().all(|c| c.count >= 3
            && c.count <= 6
            && c.band < 5
            && c.period >= -i64::from(starts)
            && c.face >= 2
            && c.wraps < 64
            && matches!(c.kind, Kind::Cone | Kind::Wall)));
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
    let calls = std::cell::Cell::new(0);
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
