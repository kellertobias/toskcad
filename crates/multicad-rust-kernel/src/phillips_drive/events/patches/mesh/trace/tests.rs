use super::*;
use std::cell::Cell;
fn options() -> Tessellation {
    Tessellation {
        deflection: 0.1,
        angle: 0.3,
        triangle_limit: 300_000,
    }
}
#[test]
fn actual_shaped_mouth_trace_preserves_source_keys_mesh_and_directed_loop() {
    for (index, scale) in [(0, 1.), (2, 1.137), (5, 1.)] {
        let (mut s, mut d) = crate::phillips_drive::test_fixture(index);
        s.mouth_diameter_mm *= scale;
        s.mouth_radius_mm *= scale;
        s.b_width_mm *= scale;
        s.g_diameter_mm *= scale;
        s.inner = match s.inner {
            Inner::Rounded { radius_mm } => Inner::Rounded {
                radius_mm: radius_mm * scale,
            },
            Inner::ShallowV {
                angle_deg,
                e_width_mm,
            } => Inner::ShallowV {
                angle_deg,
                e_width_mm: e_width_mm * scale,
            },
        };
        d.corner_start_depth_mm *= scale;
        d.corner_end_depth_mm *= scale;
        d.corner_reach_mm *= scale;
        let source = Cutter::new(s, d).unwrap();
        let atlas = source.surface_patches().unwrap();
        let old = atlas.boundary_mesh(options()).unwrap();
        let traced = atlas
            .boundary_mesh_with_mouth_trace(options(), Default::default(), &|| false)
            .unwrap();
        assert!(std::ptr::eq(traced.source(), &source));
        assert_eq!(old.mesh(), traced.boundary_mesh().mesh());
        assert_eq!(
            old.source_distance_mm(),
            traced.boundary_mesh().source_distance_mm()
        );
        assert_eq!(old.normal_chord(), traced.boundary_mesh().normal_chord());
        assert_eq!(old.refinement(), traced.boundary_mesh().refinement());
        let mut arcs = 0;
        let mut traces = 0;
        let mut vertices = std::collections::BTreeSet::new();
        for p in traced.mouth() {
            assert!(vertices.insert(p.vertex));
            let at = p.vertex as usize * 3;
            for k in 0..3 {
                assert_eq!(
                    traced.boundary_mesh().mesh().vertices[at + k],
                    b(p.source_position_mm[k]).midpoint()
                );
            }
            match p.curve {
                MeshCurveKey::Circle(MeshRow::Mouth, _) => arcs += 1,
                MeshCurveKey::Trace(..) => traces += 1,
                _ => {}
            }
            assert!((0. ..=1.).contains(&p.parameter));
        }
        assert!(arcs > 0 && traces > 0);
        // Exact indexed mesh incidence independently identifies the entire cap
        // boundary. No XY sorting, circle replacement or proximity matching.
        let cap = atlas
            .topology
            .faces()
            .iter()
            .position(|f| f.source_identity() == FaceKey::Mouth)
            .unwrap() as u32;
        let mut uses = std::collections::BTreeMap::<(u32, u32), i32>::new();
        for (tri, &face) in old
            .mesh()
            .indices
            .as_chunks::<3>()
            .0
            .iter()
            .zip(&old.mesh().face_indices)
        {
            if face != cap {
                continue;
            }
            for (a, z) in [(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])] {
                *uses.entry((a.min(z), a.max(z))).or_default() += if a < z { 1 } else { -1 };
            }
        }
        let boundary = uses
            .into_iter()
            .filter(|(_, s)| *s != 0)
            .map(|((a, z), s)| if s > 0 { (a, z) } else { (z, a) })
            .collect::<std::collections::BTreeSet<_>>();
        let ring = traced.mouth();
        let actual = (0..ring.len())
            .map(|i| (ring[i].vertex, ring[(i + 1) % ring.len()].vertex))
            .collect::<std::collections::BTreeSet<_>>();
        assert_eq!(boundary, actual);
        assert!(traced.retained_metadata_bytes() <= traced.peak_metadata_bytes());
        assert!(traced.work_used() < WORK);
    }
}
#[test]
fn actual_trace_work_bytes_and_final_cancellation_are_atomic() {
    let (s, d) = crate::phillips_drive::test_fixture(0);
    let source = Cutter::new(s, d).unwrap();
    let atlas = source.surface_patches().unwrap();
    let untraced_calls = Cell::new(0);
    let untraced_poll = || {
        untraced_calls.set(untraced_calls.get() + 1);
        false
    };
    let untraced = atlas
        .boundary_mesh_with_cancel(options(), &untraced_poll)
        .unwrap();
    let untraced_count = untraced_calls.get();
    let calls = Cell::new(0);
    let poll = || {
        calls.set(calls.get() + 1);
        false
    };
    let result = atlas
        .boundary_mesh_with_mouth_trace(options(), Default::default(), &poll)
        .unwrap();
    let polls = calls.get();
    let limits = MouthTraceOptions {
        maximum_work: result.work_used(),
        maximum_bytes: result.peak_metadata_bytes(),
    };
    assert!(
        atlas
            .boundary_mesh_with_mouth_trace(options(), limits, &|| false)
            .is_ok()
    );
    for trace in [
        MouthTraceOptions {
            maximum_work: 0,
            ..limits
        },
        MouthTraceOptions {
            maximum_work: limits.maximum_work - 1,
            ..limits
        },
        MouthTraceOptions {
            maximum_work: limits.maximum_work - result.mouth().len() - 2,
            ..limits
        },
        MouthTraceOptions {
            maximum_bytes: 0,
            ..limits
        },
        MouthTraceOptions {
            maximum_bytes: limits.maximum_bytes - 1,
            ..limits
        },
    ] {
        assert!(matches!(
            atlas.boundary_mesh_with_mouth_trace(options(), trace, &|| false),
            Err(Error::Limit)
        ));
    }
    for trace in [
        MouthTraceOptions {
            maximum_work: WORK + 1,
            ..limits
        },
        MouthTraceOptions {
            maximum_bytes: BYTES + 1,
            ..limits
        },
    ] {
        assert!(matches!(
            atlas.boundary_mesh_with_mouth_trace(options(), trace, &|| false),
            Err(Error::Invalid)
        ));
    }
    // The tail has one poll per actual mouth point, then capacity/publication
    // polls. Before it are the sort allocation/query and finite wire gates.
    // This bounded stage matrix uses actual returned source keys/counts.
    let points = result.mouth().len();
    let curves = result
        .mouth()
        .iter()
        .map(|p| p.curve)
        .collect::<std::collections::HashSet<_>>()
        .len();
    assert!(polls > points + curves + 8);
    for cutoff in [
        0,
        1,
        polls - 1,
        polls - 2,
        polls - points - 3,
        polls - points - 4,
        polls - points - curves - 7,
        polls - points - curves - 8,
    ] {
        let calls = Cell::new(0);
        let cancel = || {
            let n = calls.get();
            calls.set(n + 1);
            n >= cutoff
        };
        assert!(matches!(
            atlas.boundary_mesh_with_mouth_trace(options(), limits, &cancel),
            Err(Error::Cancelled)
        ));
    }
    let retry = atlas
        .boundary_mesh_with_mouth_trace(options(), limits, &|| false)
        .unwrap();
    assert_eq!(retry.boundary_mesh().mesh(), result.boundary_mesh().mesh());
    untraced_calls.set(0);
    let unchanged = atlas
        .boundary_mesh_with_cancel(options(), &untraced_poll)
        .unwrap();
    assert_eq!(unchanged.mesh(), untraced.mesh());
    assert_eq!(untraced_calls.get(), untraced_count);
}
#[test]
fn actual_mesh_wire_candidates_reject_missing_reversed_branch_and_second_cycles() {
    let (s, d) = crate::phillips_drive::test_fixture(0);
    let source = Cutter::new(s, d).unwrap();
    let atlas = source.surface_patches().unwrap();
    let mesh = atlas.boundary_mesh(options()).unwrap();
    let cap = atlas
        .topology
        .faces()
        .iter()
        .position(|f| f.source_identity() == FaceKey::Mouth)
        .unwrap() as u32;
    let mut sums = std::collections::BTreeMap::<(usize, usize), i32>::new();
    for (tri, &face) in mesh
        .mesh()
        .indices
        .as_chunks::<3>()
        .0
        .iter()
        .zip(&mesh.mesh().face_indices)
    {
        if face != cap {
            continue;
        }
        for (a, z) in [
            (tri[0] as usize, tri[1] as usize),
            (tri[1] as usize, tri[2] as usize),
            (tri[2] as usize, tri[0] as usize),
        ] {
            *sums.entry((a.min(z), a.max(z))).or_default() += if a < z { 1 } else { -1 };
        }
    }
    let ends = sums
        .into_iter()
        .filter(|(_, s)| *s != 0)
        .map(|((a, z), s)| if s > 0 { [a, z] } else { [z, a] })
        .collect::<Vec<_>>();
    let uses = (0..ends.len()).map(|i| (i, true)).collect::<Vec<_>>();
    let run = |ends: &[[usize; 2]], uses: &[(usize, bool)]| {
        let mut work = Counted {
            remaining: WORK,
            maximum: WORK,
            bytes: BYTES,
            peak: 0,
            cancel: &|| false,
        };
        ordered_wire(ends, uses, &mut work, 0)
    };
    assert_eq!(run(&ends, &uses).unwrap().len(), ends.len());
    assert!(matches!(run(&ends, &[]), Err(Error::Unresolved)));
    assert!(matches!(run(&ends, &uses[..1]), Err(Error::Unresolved)));
    assert!(matches!(
        run(&ends, &uses[..uses.len() - 1]),
        Err(Error::Unresolved)
    ));
    let mut reversed = uses.clone();
    reversed[1].1 = false;
    assert!(matches!(run(&ends, &reversed), Err(Error::Unresolved)));
    let mut branch = uses.clone();
    branch.push(uses[0]);
    assert!(matches!(run(&ends, &branch), Err(Error::Unresolved)));
    // A genuine non-mouth interior facet contributes a disjoint closed cycle.
    // It is an invalid proposed mouth wire, never a fabricated source solid.
    let nodes = ends
        .iter()
        .flatten()
        .copied()
        .collect::<std::collections::BTreeSet<_>>();
    let other = mesh
        .mesh()
        .indices
        .as_chunks::<3>()
        .0
        .iter()
        .find(|t| t.iter().all(|v| !nodes.contains(&(*v as usize))))
        .unwrap();
    let mut extra = ends.clone();
    extra.extend([
        [other[0] as usize, other[1] as usize],
        [other[1] as usize, other[2] as usize],
        [other[2] as usize, other[0] as usize],
    ]);
    let extra_uses = (0..extra.len()).map(|i| (i, true)).collect::<Vec<_>>();
    assert!(matches!(run(&extra, &extra_uses), Err(Error::Unresolved)));
}
