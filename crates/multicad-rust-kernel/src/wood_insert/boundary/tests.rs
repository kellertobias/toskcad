use super::*;
use crate::{authored_flat_thread::SourceOptions, helical_thread::Spec};
fn source(outer_left: bool, inner_left: bool, outer_starts: u32, inner_starts: u32) -> Insert {
    let options = SourceOptions {
        position_limit_mm: 1e-7,
        normal_angle_limit_rad: 1e-7,
    };
    let spec = Spec {
        major_radius_mm: 5.,
        pitch_mm: 1.9,
        length_mm: 4.7,
        depth_mm: 1.,
        crest_width_mm: 0.4,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts: outer_starts,
        left_hand: outer_left,
        root_radius_mm: None,
    };
    let outer = AuthoredThread::new(spec, options).unwrap();
    let inner = AuthoredThread::new(
        Spec {
            major_radius_mm: 1.5,
            pitch_mm: 0.5,
            depth_mm: 0.3,
            crest_width_mm: 0.1,
            starts: inner_starts,
            left_hand: inner_left,
            ..spec
        },
        options,
    )
    .unwrap();
    Insert::new(outer, inner, 0.5, None).unwrap()
}
#[test]
fn real_paired_sources_have_two_annuli_and_one_closed_oriented_torus() {
    for a in [false, true] {
        for b in [false, true] {
            for (outer, inner) in [(1, 1), (2, 3)] {
                let source = source(a, b, outer, inner);
                let graph = source
                    .boundary_with_interrupt(BoundaryOptions::default(), None)
                    .unwrap();
                assert_eq!(graph.source(), &source);
                assert_eq!(graph.source_boundary(Side::Outer).source(), source.outer());
                assert_eq!(graph.source_boundary(Side::Inner).source(), source.inner());
                assert_eq!(graph.face_count(), source.face_count());
                for cap in [FaceKey::StartCap, FaceKey::EndCap] {
                    assert_eq!(graph.face_loops(cap).unwrap().len(), 2);
                }
                let (v, e) = graph.topology_counts();
                assert_eq!(v as isize - e as isize + graph.face_count() as isize - 2, 0);
                assert_eq!(graph.vertices().count(), v);
                assert_eq!(graph.curves().count(), e);
                assert!(graph.work_used() <= BoundaryOptions::default().work);
                assert!(graph.retained_bytes() <= BoundaryOptions::default().bytes);
                for side in [Side::Outer, Side::Inner] {
                    assert_eq!(
                        graph.cell_surface_same_sense(side),
                        graph.source_boundary(side).cell_surface_same_sense()
                            != (side == Side::Inner)
                    );
                    for (i, _) in graph.cells(side).iter().enumerate() {
                        assert!(graph.cell_face(side, i).is_some());
                        assert!(graph.cell_coedges(side, i).unwrap().count() >= 3);
                    }
                }
            }
        }
    }
}

fn requalify(graph: &Boundary) -> Result<(usize, usize), Error> {
    graph.qualify(
        &mut Memory {
            used: graph.retained_bytes(),
            limit: BoundaryOptions::default().bytes,
        },
        &mut Work {
            remaining: 16_000_000,
            used: 0,
            interrupt: None,
        },
    )
}
#[test]
fn actual_mapped_source_use_directions_exclude_a_globally_reversed_cavity() {
    let source = source(false, true, 2, 3);
    let mut graph = source
        .boundary_with_interrupt(BoundaryOptions::default(), None)
        .unwrap();
    for face in &mut graph.loops {
        for wire in face {
            if wire
                .first()
                .is_some_and(|edge| edge.curve.owner() == Side::Inner)
            {
                wire.reverse();
                for edge in wire {
                    edge.forward = !edge.forward;
                }
            }
        }
    }
    assert_eq!(requalify(&graph), Err(Error::Unresolved));
    assert_eq!(graph.source(), &source);
}
#[test]
fn public_namespace_queries_retain_original_source_objects_and_reject_foreign_keys() {
    use std::collections::{BTreeSet, HashSet};
    let source = source(true, false, 2, 3);
    let graph = source
        .boundary_with_interrupt(BoundaryOptions::default(), None)
        .unwrap();
    let mut vertices = BTreeSet::new();
    let mut curves = HashSet::new();
    for (id, point) in graph.vertices() {
        assert!(vertices.insert(id));
        assert!(std::ptr::eq(graph.point(id).unwrap(), point));
        let _ = format!("{id:?}");
    }
    for (id, curve) in graph.curves() {
        assert!(curves.insert(id));
        assert!(std::ptr::eq(graph.curve(id).unwrap(), curve));
        let bare = graph.source_boundary(id.owner());
        assert!(std::ptr::eq(
            bare.source_curve(id.source_index()).unwrap(),
            curve
        ));
        let line = graph.curve_polyline(id, 0.05, 10000, None).unwrap();
        assert!(line.maximum_error_mm <= 0.05);
        assert_eq!(
            line.points_mm.first().copied(),
            Some(graph.curve_point(id, 0., None).unwrap().point_mm)
        );
        assert_eq!(
            line.points_mm.last().copied(),
            Some(graph.curve_point(id, 1., None).unwrap().point_mm)
        );
        assert_eq!(
            graph.curve_point(id, 0.37, None).unwrap(),
            bare.curve_point(id.source_index(), 0.37, None).unwrap()
        );
        assert_eq!(
            graph.curve_endpoints(id).unwrap(),
            curve.endpoints().map(|i| vertex_key(id.owner(), i))
        );
        assert!(graph.curve_faces(id).is_some());
        let _ = format!("{id:?}");
    }
    assert_ne!(BoundaryVertexKey::Outer(0), BoundaryVertexKey::Inner(0));
    assert!(graph.point(BoundaryVertexKey::Outer(u32::MAX)).is_none());
    assert!(graph.point(BoundaryVertexKey::Inner(u32::MAX)).is_none());
    let bad = BoundaryCurveKey::Inner(u32::MAX);
    assert!(graph.curve(bad).is_none());
    assert!(graph.curve_endpoints(bad).is_none());
    assert!(graph.curve_faces(bad).is_none());
    assert_eq!(graph.curve_point(bad, 0., None), Err(Error::Invalid));
    assert_eq!(
        graph.curve_point(BoundaryCurveKey::Outer(0), f64::NAN, None),
        Err(Error::Invalid)
    );
    assert_eq!(
        graph.curve_point(BoundaryCurveKey::Outer(0), 0., Some(&|| true)),
        Err(Error::Cancelled)
    );
    assert!(
        graph
            .face_loops(FaceKey::InnerThread {
                start: u32::MAX,
                band: crate::authored_internal_treatment::ThreadBandRole::Root
            })
            .is_none()
    );
    for side in [Side::Outer, Side::Inner] {
        assert!(graph.cell_face(side, usize::MAX).is_none());
        assert!(matches!(
            graph.cell_coedges(side, usize::MAX),
            Err(Error::Invalid)
        ));
        for i in 0..graph.cells(side).len() {
            let bare = graph
                .source_boundary(side)
                .cell_coedges(i)
                .unwrap()
                .collect::<Vec<_>>();
            let actual = graph.cell_coedges(side, i).unwrap().collect::<Vec<_>>();
            for (j, edge) in actual.iter().enumerate() {
                let original = bare[if side == Side::Inner {
                    bare.len() - 1 - j
                } else {
                    j
                }];
                assert_eq!(edge.curve, curve_key(side, original.curve as u32));
                assert_eq!(edge.forward, original.forward != (side == Side::Inner));
            }
        }
    }
}

#[test]
fn aggregate_packet_limits_are_not_reset_for_each_component_or_topology_query() {
    let source = source(false, true, 1, 1);
    let options = BoundaryOptions::default();
    assert!(matches!(
        source.boundary_with_interrupt(options, Some(&|| true)),
        Err(Error::Cancelled)
    ));
    for invalid in [
        BoundaryOptions {
            cells: 0,
            ..options
        },
        BoundaryOptions {
            vertices: 8193,
            ..options
        },
        BoundaryOptions {
            curves: 0,
            ..options
        },
        BoundaryOptions { work: 0, ..options },
        BoundaryOptions {
            bytes: 0,
            ..options
        },
    ] {
        assert!(matches!(
            source.boundary_with_interrupt(invalid, None),
            Err(Error::Invalid)
        ));
    }
    for budget in [1, 2, 20] {
        assert!(matches!(
            source.boundary_with_interrupt(
                BoundaryOptions {
                    work: budget,
                    ..options
                },
                None
            ),
            Err(Error::Limit)
        ));
    }
    for budget in [1, 1024, 8 * 1024 * 1024] {
        assert!(matches!(
            source.boundary_with_interrupt(
                BoundaryOptions {
                    bytes: budget,
                    ..options
                },
                None
            ),
            Err(Error::Limit)
        ));
    }
    let outer = source
        .outer()
        .boundary_with_interrupt(options, None)
        .unwrap();
    let inner = source
        .inner()
        .boundary_with_interrupt(options, None)
        .unwrap();
    let single_work = outer.work_used().max(inner.work_used());
    for component in [source.outer(), source.inner()] {
        assert!(
            component
                .boundary_with_interrupt(
                    BoundaryOptions {
                        work: single_work,
                        ..options
                    },
                    None
                )
                .is_ok()
        );
    }
    assert!(matches!(
        source.boundary_with_interrupt(
            BoundaryOptions {
                work: single_work,
                ..options
            },
            None
        ),
        Err(Error::Limit)
    ));
    let baseline = source.boundary_with_interrupt(options, None).unwrap();
    assert!(matches!(
        source.boundary_with_interrupt(
            BoundaryOptions {
                work: baseline.work_used() - 1,
                ..options
            },
            None
        ),
        Err(Error::Limit)
    ));
    for requested in [
        BoundaryOptions {
            cells: 1,
            ..options
        },
        BoundaryOptions {
            vertices: 1,
            ..options
        },
        BoundaryOptions {
            curves: 1,
            ..options
        },
    ] {
        assert!(matches!(
            source.boundary_with_interrupt(requested, None),
            Err(Error::Limit)
        ));
    }
    // Search only a resource threshold, with the exact same real source. Both
    // individually fitting components still need aggregate packet ownership.
    let mut low = 1;
    let mut high = options.bytes;
    while low < high {
        let mid = low + (high - low) / 2;
        let fits = [source.outer(), source.inner()]
            .into_iter()
            .all(|component| {
                component
                    .boundary_with_interrupt(
                        BoundaryOptions {
                            bytes: mid,
                            ..options
                        },
                        None,
                    )
                    .is_ok()
            });
        if fits {
            high = mid;
        } else {
            low = mid + 1;
        }
    }
    assert!(
        [source.outer(), source.inner()]
            .into_iter()
            .all(|component| component
                .boundary_with_interrupt(
                    BoundaryOptions {
                        bytes: low,
                        ..options
                    },
                    None
                )
                .is_ok())
    );
    assert!(matches!(
        source.boundary_with_interrupt(
            BoundaryOptions {
                bytes: low,
                ..options
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert_eq!(baseline.source(), &source);
}
#[test]
fn cancellation_across_both_real_sources_annular_sewing_and_final_qualification_is_atomic() {
    let source = source(true, false, 1, 2);
    let calls = std::cell::Cell::new(0usize);
    let count = || {
        calls.set(calls.get() + 1);
        false
    };
    let graph = source
        .boundary_with_interrupt(BoundaryOptions::default(), Some(&count))
        .unwrap();
    let total = calls.get();
    assert!(total > 100);
    for at in [
        1,
        2,
        10,
        100,
        total / 4,
        total / 2,
        3 * total / 4,
        total - 2,
        total - 1,
        total,
    ] {
        calls.set(0);
        let cancel = || {
            calls.set(calls.get() + 1);
            calls.get() >= at
        };
        assert!(matches!(
            source.boundary_with_interrupt(BoundaryOptions::default(), Some(&cancel)),
            Err(Error::Cancelled)
        ));
    }
    assert_eq!(graph.source(), &source);
    assert!(
        source
            .boundary_with_interrupt(BoundaryOptions::default(), None)
            .is_ok()
    );
}
#[test]
fn missing_foreign_duplicate_and_misordered_actual_wire_uses_never_qualify() {
    let source = source(false, false, 1, 1);
    let fresh = || {
        source
            .boundary_with_interrupt(BoundaryOptions::default(), None)
            .unwrap()
    };
    let mut graph = fresh();
    graph.loops.pop();
    assert_eq!(requalify(&graph), Err(Error::Unresolved));
    let mut graph = fresh();
    graph.loops[0].pop();
    assert_eq!(requalify(&graph), Err(Error::Unresolved));
    let mut graph = fresh();
    graph.loops[0][0].clear();
    assert_eq!(requalify(&graph), Err(Error::Unresolved));
    let mut graph = fresh();
    graph.loops[0][0][0].curve = BoundaryCurveKey::Outer(u32::MAX);
    assert_eq!(requalify(&graph), Err(Error::Unresolved));
    let mut graph = fresh();
    graph.loops[0][0][0].curve = graph.loops[1][0][0].curve;
    assert_eq!(requalify(&graph), Err(Error::Unresolved));
    let mut graph = fresh();
    let edge = graph.loops[0][0][0];
    graph.loops[0][0].push(edge);
    assert_eq!(requalify(&graph), Err(Error::Unresolved));
    let mut graph = fresh();
    graph.loops[0][0].swap(0, 1);
    assert_eq!(requalify(&graph), Err(Error::Unresolved));
    let mut graph = fresh();
    graph.loops[0][0][0].forward = !graph.loops[0][0][0].forward;
    assert_eq!(requalify(&graph), Err(Error::Unresolved));
}

#[test]
fn allocation_accounting_and_incidence_only_checks_refuse_overflow_and_nonmanifold_data() {
    use super::topology::{Corner, Use, edge_receipt, shell_receipt, vertex_links};
    let mut memory = Memory { used: 1, limit: 2 };
    assert_eq!(memory.add(2), Err(Error::Limit));
    assert!(matches!(
        memory.allocate::<u64>(
            2,
            &mut Work {
                remaining: 100,
                used: 0,
                interrupt: None
            }
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        memory.retain(Vec::<u64>::with_capacity(2)),
        Err(Error::Limit)
    ));
    let mut unlimited = Work {
        remaining: usize::MAX,
        used: 0,
        interrupt: None,
    };
    assert!(matches!(
        memory.allocate::<u8>(usize::MAX, &mut unlimited),
        Err(Error::Limit)
    ));
    let mut unlimited = Work {
        remaining: usize::MAX,
        used: 0,
        interrupt: None,
    };
    let mut huge = Memory {
        used: 0,
        limit: usize::MAX,
    };
    assert!(matches!(
        huge.allocate::<u8>(isize::MAX as usize + 1, &mut unlimited),
        Err(Error::Limit)
    ));
    assert_eq!(
        super::topology::edge_count([(true, [None, None])], &mut unlimited),
        Err(Error::Unresolved)
    );
    let mut ordinary = Memory {
        used: 0,
        limit: 100,
    };
    assert!(
        ordinary
            .allocate::<u8>(2, &mut unlimited)
            .unwrap()
            .capacity()
            >= 2
    );
    assert!(
        ordinary
            .allocate::<u64>(2, &mut unlimited)
            .unwrap()
            .capacity()
            >= 2
    );
    let incidence = |uses| super::topology::Incidence {
        uses,
        heads: vec![],
        corners: vec![],
        seen: vec![],
        reached: vec![false],
        queue: vec![],
        chi: 0,
    };
    assert_eq!(
        incidence(vec![[None, None]]).finish([true], (0, 0), &mut unlimited),
        Err(Error::Unresolved)
    );
    assert_eq!(
        incidence(vec![]).finish([], (1, 0), &mut unlimited),
        Err(Error::Unresolved)
    );
    let use0 = Use {
        face: 0,
        forward: false,
    };
    assert!(edge_receipt(true, [Some(use0), None]).is_err());
    assert!(edge_receipt(false, [Some(use0), None]).is_err());
    let corner = |a, b, next| Corner {
        germs: [a, b],
        next,
    };
    let cases = [
        (vec![0], vec![corner(0, 1, usize::MAX)]),
        (
            vec![0],
            vec![corner(0, 1, 1), corner(0, 1, 2), corner(0, 1, usize::MAX)],
        ),
        (vec![0, 0], vec![corner(0, 1, 1), corner(0, 1, usize::MAX)]),
        (
            vec![0],
            vec![
                corner(0, 1, 1),
                corner(0, 1, 2),
                corner(2, 3, 3),
                corner(2, 3, usize::MAX),
            ],
        ),
    ];
    for (heads, corners) in cases {
        assert_eq!(
            vertex_links(
                &heads,
                &corners,
                &mut [false; 8],
                &mut Work {
                    remaining: 1000,
                    used: 0,
                    interrupt: None
                }
            ),
            Err(Error::Unresolved)
        );
    }
    assert_eq!(
        shell_receipt((1, 1), (1, 2), 0, &[true]),
        Err(Error::Unresolved)
    );
    assert_eq!(
        shell_receipt((1, 1), (1, 1), 0, &[false]),
        Err(Error::Unresolved)
    );
    assert_eq!(
        shell_receipt((1, 1), (1, 1), 1, &[true]),
        Err(Error::Unresolved)
    );
}
#[test]
fn hidden_cell_seams_and_duplicate_cycles_cannot_be_promoted_to_material_wires() {
    let source = source(false, false, 1, 1);
    let mut graph = source
        .boundary_with_interrupt(BoundaryOptions::default(), None)
        .unwrap();
    let hidden = graph
        .outer
        .source_curves()
        .find(|(_, curve)| curve.faces()[0] == curve.faces()[1] && curve.faces()[0] >= 2)
        .map(|(i, c)| (i, c.faces()[0]))
        .unwrap();
    let (id, face) = hidden;
    let forward = graph.outer.source_curve(id).unwrap().face_uses()[0].1;
    graph.loops[face as usize][0] = vec![BoundaryCoedge {
        curve: BoundaryCurveKey::Outer(id),
        forward,
    }];
    assert_eq!(requalify(&graph), Err(Error::Unresolved));
    let mut graph = source
        .boundary_with_interrupt(BoundaryOptions::default(), None)
        .unwrap();
    let copy = graph.loops[0][0].clone();
    graph.loops[0][0].extend(copy);
    assert_eq!(requalify(&graph), Err(Error::Unresolved));
}
#[test]
fn original_curve_precision_and_parent_query_bytes_remain_real_limits() {
    let source = source(false, true, 1, 1);
    let outer = source
        .outer()
        .boundary_with_interrupt(BoundaryOptions::default(), None)
        .unwrap();
    assert!(matches!(
        source.boundary_with_interrupt(
            BoundaryOptions {
                work: outer.work_used() + 1,
                ..BoundaryOptions::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
    let limits = BoundaryOptions {
        cells: 256,
        vertices: 1024,
        curves: 1024,
        bytes: 9 * 1024 * 1024,
        ..BoundaryOptions::default()
    };
    let graph = source.boundary_with_interrupt(limits, None).unwrap();
    let key = graph
        .outer
        .curves()
        .max_by(|(_, a), (_, b)| {
            let score = |c: &BoundaryCurve| {
                let p = c.affine_rzv_bounds();
                p[0][0].upper * p[1][2].upper.abs().powi(2)
            };
            score(a).total_cmp(&score(b))
        })
        .map(|(id, _)| BoundaryCurveKey::Outer(id))
        .unwrap();
    assert!(matches!(
        graph.curve_polyline(key, 2e-10, 500000, None),
        Err(Error::Limit)
    ));
    assert!(matches!(
        graph.curve_polyline(key, 0.05, 1000, Some(&|| true)),
        Err(Error::Cancelled)
    ));
    let original = source.outer().spec();
    let options = SourceOptions {
        position_limit_mm: 1e-7,
        normal_angle_limit_rad: 1e-7,
    };
    let tiny_source = Insert::new(
        AuthoredThread::new(
            Spec {
                length_mm: 4.,
                ..original
            },
            options,
        )
        .unwrap(),
        AuthoredThread::new(
            Spec {
                length_mm: 4.,
                ..source.inner().spec()
            },
            options,
        )
        .unwrap(),
        0.5,
        None,
    )
    .unwrap();
    let tiny = tiny_source
        .boundary_with_interrupt(BoundaryOptions::default(), None)
        .unwrap();
    let key = tiny
        .outer
        .curves()
        .min_by(|(_, a), (_, b)| {
            a.affine_rzv_bounds()[1][2]
                .upper
                .abs()
                .total_cmp(&b.affine_rzv_bounds()[1][2].upper.abs())
        })
        .map(|(id, _)| BoundaryCurveKey::Outer(id))
        .unwrap();
    assert!(matches!(
        tiny.curve_polyline(key, 1e-20, 1000, None),
        Err(Error::Unresolved)
    ));
}
