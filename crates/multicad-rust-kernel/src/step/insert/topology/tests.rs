use super::*;
#[test]
fn actual_local_tensor_shell_has_annular_caps_and_closed_manifold_genus_one() {
    for kind in 0..4 {
        let body = if kind == 0 {
            super::super::tests::metric_graph()
        } else {
            super::super::tests::graph(kind == 2, if kind == 3 { 2 } else { 1 })
        };
        let d = DomainPlan::prepare(&body, DomainOptions::default(), None).unwrap();
        let s = SurfacePlan::prepare(&d, CurveOptions::default(), None).unwrap();
        let t = TopologyPlan::prepare(&s, 16_000_000, 64 * 1024 * 1024, None).unwrap();
        assert!(std::ptr::eq(t.surface(), &s));
        assert_eq!(t.faces().len(), s.patches().len() + 2);
        assert_eq!(t.used_vertices(), s.vertices().len());
        assert_eq!(
            t.used_vertices() as isize - t.edges().len() as isize + t.faces().len() as isize - 2,
            0
        );
        for cap in &t.faces()[s.patches().len()..] {
            assert!(matches!(cap.geometry, FaceGeometry::Cap(_)));
            assert_eq!(cap.loops.len(), 2);
        }
        let winding = t
            .winding_bounds(16_000_000, 64 * 1024 * 1024, None)
            .unwrap();
        assert!(winding.certified_segments > 0);
        assert!(t.work_used() <= 16_000_000);
        assert!(t.retained_bytes() <= 64 * 1024 * 1024);
        eprintln!(
            "TOPO {kind} V{} E{} F{} work{}",
            t.used_vertices(),
            t.edges().len(),
            t.faces().len(),
            t.work_used()
        );
    }
}
#[test]
fn graph_and_winding_aggregate_ledgers_refuse_before_publication() {
    let body = super::super::tests::metric_graph();
    let d = DomainPlan::prepare(&body, DomainOptions::default(), None).unwrap();
    let s = SurfacePlan::prepare(&d, CurveOptions::default(), None).unwrap();
    for bytes in [0, 64 * 1024 * 1024 + 1] {
        assert!(matches!(
            TopologyPlan::prepare(&s, 16_000_000, bytes, None),
            Err(Error::Invalid)
        ));
    }
    assert!(matches!(
        TopologyPlan::prepare(&s, 1, 64 * 1024 * 1024, None),
        Err(Error::Limit)
    ));
    assert!(matches!(
        TopologyPlan::prepare(&s, 16_000_000, s.retained_bytes(), None),
        Err(Error::Limit)
    ));
    assert!(matches!(
        TopologyPlan::prepare(&s, 16_000_000, 64 * 1024 * 1024, Some(&|| true)),
        Err(Error::Cancelled)
    ));
    let t = TopologyPlan::prepare(&s, 16_000_000, 64 * 1024 * 1024, None).unwrap();
    for bytes in [0, 64 * 1024 * 1024 + 1] {
        assert!(matches!(
            t.winding_bounds(16_000_000, bytes, None),
            Err(Error::Invalid)
        ));
    }
    assert!(matches!(
        t.winding_bounds(16_000_000, t.retained_bytes(), None),
        Err(Error::Limit)
    ));
    assert!(matches!(
        t.winding_bounds(1, 64 * 1024 * 1024, None),
        Err(Error::Limit)
    ));
    assert!(matches!(
        t.winding_bounds(16_000_000, 64 * 1024 * 1024, Some(&|| true)),
        Err(Error::Cancelled)
    ));
}
#[test]
fn literal_graph_mutations_refuse_open_misoriented_and_nonmanifold_incidence() {
    let body = super::super::tests::metric_graph();
    let d = DomainPlan::prepare(&body, DomainOptions::default(), None).unwrap();
    let s = SurfacePlan::prepare(&d, CurveOptions::default(), None).unwrap();
    for kind in 0..9 {
        let mut t = TopologyPlan::prepare(&s, 16_000_000, 64 * 1024 * 1024, None).unwrap();
        match kind {
            0 => {
                t.edges[0].ends[1] = t.edges[0].ends[0];
            }
            1 => {
                t.edges[0].ends[1] = s.vertices().len();
            }
            2 => {
                t.faces[0].loops.clear();
            }
            3 => {
                t.faces[0].loops[0].truncate(2);
            }
            4 => {
                t.faces[0].loops[0][0].edge = t.edges.len();
            }
            5 => {
                t.faces[0].loops[0][0].forward ^= true;
            }
            6 => {
                let c = t.faces[0].loops[0][0];
                t.faces[0].loops[0].push(c);
            }
            7 => {
                t.faces.pop();
            }
            _ => {
                let c = t.faces[0].loops[0][0];
                t.faces[0].loops[0][1] = c;
            }
        }
        let mut w = Work::new(16_000_000, None).unwrap();
        let mut memory = Memory::new(t.retained_bytes(), 64 * 1024 * 1024).unwrap();
        assert!(
            matches!(
                verify(s.vertices().len(), &t.edges, &t.faces, &mut memory, &mut w),
                Err(Error::Unresolved)
            ),
            "mutation{kind}"
        );
    }
}
#[test]
fn actual_closed_graph_rejects_disconnected_pinched_and_wrong_genus_shells() {
    let body = super::super::tests::metric_graph();
    let d = DomainPlan::prepare(&body, DomainOptions::default(), None).unwrap();
    let s = SurfacePlan::prepare(&d, CurveOptions::default(), None).unwrap();
    let t = TopologyPlan::prepare(&s, 16_000_000, 64 * 1024 * 1024, None).unwrap();
    let clone_faces = || {
        t.faces
            .iter()
            .map(|f| ShellFace {
                geometry: f.geometry,
                same_sense: f.same_sense,
                loops: f.loops.clone(),
            })
            .collect::<Vec<_>>()
    };
    let run = |edges: &[ShellEdge], faces: &[ShellFace], vertices| {
        let mut w = Work::new(16_000_000, None).unwrap();
        let mut m = Memory::new(0, 64 * 1024 * 1024).unwrap();
        verify(vertices, edges, faces, &mut m, &mut w)
    };
    let mut faces = clone_faces();
    faces[0].loops[0].reverse();
    for c in &mut faces[0].loops[0] {
        c.forward ^= true;
    }
    assert!(matches!(
        run(&t.edges, &faces, s.vertices().len()),
        Err(Error::Unresolved)
    ));
    let mut faces = clone_faces();
    let f = &t.faces[0];
    faces.push(ShellFace {
        geometry: f.geometry,
        same_sense: f.same_sense,
        loops: f.loops.clone(),
    });
    assert!(matches!(
        run(&t.edges, &faces, s.vertices().len()),
        Err(Error::Unresolved)
    ));
    let mut faces = clone_faces();
    faces[0].loops[0].last_mut().unwrap().edge = usize::MAX;
    assert!(matches!(
        run(&t.edges, &faces, s.vertices().len()),
        Err(Error::Unresolved)
    ));
    let mut edges = t.edges.clone();
    let mut faces = clone_faces();
    edges.extend(t.edges.iter().map(|e| ShellEdge {
        ends: e.ends.map(|v| v + s.vertices().len()),
        geometry: e.geometry,
    }));
    faces.extend(t.faces.iter().map(|f| {
        ShellFace {
            geometry: f.geometry,
            same_sense: f.same_sense,
            loops: f
                .loops
                .iter()
                .map(|wire| {
                    wire.iter()
                        .map(|c| ShellCoedge {
                            edge: c.edge + t.edges.len(),
                            forward: c.forward,
                        })
                        .collect()
                })
                .collect(),
        }
    }));
    assert!(matches!(
        run(&edges, &faces, 2 * s.vertices().len()),
        Err(Error::Unresolved)
    ));
    // Identify two actual, nonadjacent vertices: edge pairing and face
    // connectivity survive, but the resulting pinched vertex has TWO link fans.
    let a = t.edges[0].ends[0];
    let b = (0..s.vertices().len())
        .find(|&b| {
            b != a
                && !t
                    .edges
                    .iter()
                    .any(|e| e.ends.contains(&a) && e.ends.contains(&b))
        })
        .unwrap();
    let edges = t
        .edges
        .iter()
        .map(|e| ShellEdge {
            ends: e.ends.map(|v| if v == b { a } else { v }),
            geometry: e.geometry,
        })
        .collect::<Vec<_>>();
    assert!(matches!(
        run(&edges, &t.faces, s.vertices().len()),
        Err(Error::Unresolved)
    ));
    // A real closed tetrahedron is manifold and connected, but has genus zero.
    let mut edges = Vec::<ShellEdge>::new();
    let mut faces = Vec::new();
    for tri in [[0, 2, 1], [0, 1, 3], [1, 2, 3], [2, 0, 3]] {
        let mut wire = Vec::new();
        for i in 0..3 {
            let pair = [tri[i], tri[(i + 1) % 3]];
            let ends = if pair[0] < pair[1] {
                pair
            } else {
                [pair[1], pair[0]]
            };
            let edge = edges
                .iter()
                .position(|e| e.ends == ends)
                .unwrap_or_else(|| {
                    edges.push(ShellEdge {
                        ends,
                        geometry: EdgeGeometry::Ruling,
                    });
                    edges.len() - 1
                });
            wire.push(ShellCoedge {
                edge,
                forward: pair == ends,
            });
        }
        faces.push(ShellFace {
            geometry: FaceGeometry::Cap(0),
            same_sense: true,
            loops: vec![wire],
        });
    }
    assert!(matches!(run(&edges, &faces, 4), Err(Error::Unresolved)));
}
#[test]
fn winding_refuses_modified_literal_cap_ownership_and_reversed_single_loop() {
    let body = super::super::tests::metric_graph();
    let d = DomainPlan::prepare(&body, DomainOptions::default(), None).unwrap();
    let s = SurfacePlan::prepare(&d, CurveOptions::default(), None).unwrap();
    let mut t = TopologyPlan::prepare(&s, 16_000_000, 64 * 1024 * 1024, None).unwrap();
    let cap = s.patches().len();
    let row = t.faces[cap].loops[0][0].edge;
    let geometry = t.edges[row].geometry;
    t.edges[row].geometry = EdgeGeometry::Ruling;
    assert!(matches!(
        t.winding_bounds(16_000_000, 64 * 1024 * 1024, None),
        Err(Error::Unresolved)
    ));
    t.edges[row].geometry = geometry;
    t.faces[cap].loops[0].reverse();
    for c in &mut t.faces[cap].loops[0] {
        c.forward ^= true;
    }
    assert!(matches!(
        t.winding_bounds(16_000_000, 64 * 1024 * 1024, None),
        Err(Error::Unresolved)
    ));
}
#[test]
fn actual_cap_identity_and_shared_ruling_cache_reject_foreign_endpoint_aliases() {
    let body = super::super::tests::metric_graph();
    let d = DomainPlan::prepare(&body, DomainOptions::default(), None).unwrap();
    let s = SurfacePlan::prepare(&d, CurveOptions::default(), None).unwrap();
    let cap = body
        .face_loops(crate::wood_insert::FaceKey::StartCap)
        .unwrap();
    let original = body.curve_endpoints(cap[0][0].curve).unwrap();
    let foreign = body.curve_endpoints(cap[1][0].curve).unwrap();
    assert!(matches!(
        cap_forward(foreign.map(VertexOwner::Original), original, true),
        Err(Error::Unresolved)
    ));
    let p = s
        .patches()
        .iter()
        .find(|p| {
            let [a, b] = p.rows.map(|i| &s.rows()[i]);
            a.endpoints()[0] != b.endpoints()[0]
        })
        .unwrap();
    let slab = &d.slabs()[p.slab];
    let [a, b] = p.rows.map(|i| &s.rows()[i]);
    let ends = [a.endpoints()[0], b.endpoints()[0]];
    let mut keys = Vec::new();
    let mut edges = Vec::new();
    let mut w = Work::new(16_000_000, None).unwrap();
    let e = ruling(&s, slab, 0, ends, &mut keys, &mut edges, &mut w).unwrap();
    edges[e.edge].ends.reverse();
    assert!(matches!(
        ruling(&s, slab, 0, ends, &mut keys, &mut edges, &mut w),
        Err(Error::Unresolved)
    ));
}
#[test]
fn cancellation_at_every_actual_small_shell_checkpoint_is_atomic() {
    use crate::{
        authored_flat_thread::{AuthoredThread, SourceOptions},
        helical_thread::Spec,
        wood_insert::{BoundaryOptions, Insert},
    };
    let spec = Spec {
        major_radius_mm: 5.,
        pitch_mm: 2.,
        length_mm: 2.,
        depth_mm: 1.,
        crest_width_mm: 0.4,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts: 1,
        left_hand: false,
        root_radius_mm: None,
    };
    let options = SourceOptions {
        position_limit_mm: 1e-7,
        normal_angle_limit_rad: 1e-7,
    };
    let outer = AuthoredThread::new(spec, options).unwrap();
    let inner = AuthoredThread::new(
        Spec {
            major_radius_mm: 1.5,
            pitch_mm: 2.,
            depth_mm: 0.3,
            crest_width_mm: 0.1,
            ..spec
        },
        options,
    )
    .unwrap();
    let body = Insert::new(outer, inner, 0.5, None)
        .unwrap()
        .boundary_with_interrupt(BoundaryOptions::default(), None)
        .unwrap();
    let d = DomainPlan::prepare(&body, DomainOptions::default(), None).unwrap();
    let s = SurfacePlan::prepare(&d, CurveOptions::default(), None).unwrap();
    let polls = std::cell::Cell::new(0);
    let count = || {
        polls.set(polls.get() + 1);
        false
    };
    TopologyPlan::prepare(&s, 16_000_000, 64 * 1024 * 1024, Some(&count)).unwrap();
    let total = polls.get();
    for stop in 1..=total {
        polls.set(0);
        let cancel = || {
            polls.set(polls.get() + 1);
            polls.get() == stop
        };
        assert!(
            matches!(
                TopologyPlan::prepare(&s, 16_000_000, 64 * 1024 * 1024, Some(&cancel)),
                Err(Error::Cancelled)
            ),
            "checkpoint{stop}/{total}"
        );
    }
    eprintln!("actual small shell cancellation checkpoints {total}");
}
