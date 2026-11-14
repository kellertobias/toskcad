use super::*;
use crate::authored_metric_thread::Spec;
fn source(left: bool, nb: bool) -> AuthoredMetricThread {
    let p = 0.5;
    let major = (2.874 + 2.980) / 4.;
    let pd = (2.580 + 2.655) / 2.;
    AuthoredMetricThread::new(Spec {
        major_radius_mm: major,
        pitch_mm: p,
        length_mm: if nb { 12.7 } else { 12. },
        crest_width_mm: p / 2. - (2. * major - pd) / 3_f64.sqrt(),
        starts: 1,
        left_hand: left,
    })
    .unwrap()
}
#[test]
fn actual_canonical_rounded_metric_shell_is_source_qualified() {
    for (left, nb) in [(false, false), (true, true)] {
        let source = source(left, nb);
        let timer = std::time::Instant::now();
        let shell = metric_shell(
            &source,
            Tessellation {
                deflection: 0.1,
                angle: 0.3,
                triangle_limit: 1_000_000,
            },
            16_000_000,
            None,
        )
        .unwrap();
        eprintln!(
            "actualmetricLH{left}nb{nb} {}facets {}work {:?} chord{}normal{}",
            shell.mesh.face_indices.len(),
            shell.work,
            timer.elapsed(),
            shell.chord,
            shell.normal
        );
        assert_eq!(*shell.source(), source);
        assert!(shell.source_chord_bound_mm() <= 0.1);
        assert!(shell.source_normal_angle_bound_rad() <= 0.3);
        assert!(shell.rim_count(ShellEnd::Start) > 3 && shell.rim_count(ShellEnd::End) > 3);
    }
}
#[test]
fn clipped_root_cell_fan_diagnostic() {
    let source = source(false, false);
    let mut w = Work {
        remaining: 16_000_000,
        used: 0,
        interrupt: None,
    };
    let grid = Grid {
        counts: [1, 1, 16, 1, 1],
        rows: 64,
        odd: false,
        shifted: false,
    };
    let patch = Patch::new(&source, 23, 2, grid, &mut w).unwrap();
    let strip = patch.strip(2, &mut w).unwrap();
    let polygon = patch.cell(&strip, 37, &mut w).unwrap();
    let mut points = Vec::new();
    let mut errors = Vec::new();
    for p in &polygon {
        let point = patch.point(p, &strip, &mut w).unwrap();
        let literal = point.map(I::midpoint);
        errors.push(
            point
                .into_iter()
                .zip(literal)
                .fold(I::point(0.), |s, (p, l)| {
                    s + I::point((p - I::point(l)).magnitude())
                })
                .hi,
        );
        points.push(literal);
    }
    assert_eq!(polygon.len(), 5);
    let field = {
        let _old = OldRoot::new();
        patch.field(&polygon, &strip, &mut w).unwrap()
    };
    assert!(matches!(
        qualify_attempt(&field, &points, &errors, [0, 2, 3], options(), &mut w).unwrap(),
        Attempt::Refine
    ));
    let proofs = match qualify_fan(&field, &points, &errors, options(), &mut w).unwrap() {
        Attempt::Complete(p) => p,
        _ => panic!("actual alternative fan must qualify"),
    };
    assert_eq!(
        proofs
            .into_iter()
            .flatten()
            .map(|p| p.ids)
            .collect::<Vec<_>>(),
        vec![[2, 3, 4], [2, 4, 0], [2, 0, 1]]
    );
    assert!(matches!(
        qualify_fan(
            &field,
            &points,
            &errors,
            Tessellation {
                angle: 0.01,
                ..options()
            },
            &mut w
        )
        .unwrap(),
        Attempt::Refine
    ));
    let collapsed = vec![points[0]; points.len()];
    assert!(matches!(
        qualify_fan(&field, &collapsed, &errors, options(), &mut w).unwrap(),
        Attempt::Repartition
    ));
}
thread_local! { pub(super) static OLD_ROOT_NORMAL:std::cell::Cell<bool> = const { std::cell::Cell::new(false) }; }
thread_local! { pub(super) static FORCE_CLIP: std::cell::Cell<bool> = const { std::cell::Cell::new(false) }; }
struct FullClip;
impl FullClip {
    fn new() -> Self {
        FORCE_CLIP.set(true);
        Self
    }
}
impl Drop for FullClip {
    fn drop(&mut self) {
        FORCE_CLIP.set(false);
    }
}
fn options() -> Tessellation {
    Tessellation {
        deflection: 0.1,
        angle: 0.3,
        triangle_limit: 1_000_000,
    }
}
pub(super) fn small_source() -> AuthoredMetricThread {
    let mut spec = source(false, false).spec();
    spec.length_mm = 1.;
    AuthoredMetricThread::new(spec).unwrap()
}
#[test]
fn strict_affine_preflight_matches_forced_exact_clip_and_owned_rims() {
    let source = small_source();
    let fast = source
        .tessellate_shell(options(), 16_000_000, None)
        .unwrap();
    let _forced = FullClip::new();
    let exact = source
        .tessellate_shell(options(), 16_000_000, None)
        .unwrap();
    assert_eq!(fast.mesh().vertices, exact.mesh().vertices);
    assert_eq!(fast.mesh().indices, exact.mesh().indices);
    assert_eq!(fast.mesh().face_indices, exact.mesh().face_indices);
    for end in [ShellEnd::Start, ShellEnd::End] {
        assert_eq!(fast.rim_count(end), exact.rim_count(end));
        for i in 0..fast.rim_count(end) {
            let a = fast.rim_edge(end, i).unwrap();
            assert_eq!(Some(a), exact.rim_edge(end, i));
            for v in a.vertex_indices {
                assert_eq!(
                    fast.exact_rim_turns(v, 10000, None).unwrap().0,
                    exact.exact_rim_turns(v, 10000, None).unwrap().0
                );
            }
        }
        assert!(fast.rim_edge(end, usize::MAX).is_none());
    }
    assert!(fast.retained_bytes() < BYTES);
    assert!(fast.work_used() < exact.work_used());
}
#[test]
fn owned_shell_has_only_its_two_actual_boundary_cycles() {
    let shell = small_source()
        .tessellate_shell(options(), 16_000_000, None)
        .unwrap();
    let mut edges = std::collections::BTreeMap::<(u32, u32), Vec<(u32, u32)>>::new();
    for tri in shell.mesh().indices.as_chunks::<3>().0 {
        for (a, b) in [(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])] {
            edges.entry((a.min(b), a.max(b))).or_default().push((a, b));
        }
    }
    let mut boundary = std::collections::BTreeSet::new();
    for (edge, uses) in edges {
        if uses.len() == 1 {
            boundary.insert(edge);
        } else {
            assert_eq!(uses.len(), 2);
            assert_eq!(uses[0], (uses[1].1, uses[1].0));
        }
    }
    let mut expected = std::collections::BTreeSet::new();
    for end in [ShellEnd::Start, ShellEnd::End] {
        let mut degree = std::collections::BTreeMap::new();
        for i in 0..shell.rim_count(end) {
            let [a, b] = shell
                .rim_edge(end, i)
                .unwrap()
                .vertex_indices
                .map(|x| x as u32);
            assert!(expected.insert((a.min(b), a.max(b))));
            *degree.entry(a).or_insert(0) += 1;
            *degree.entry(b).or_insert(0) += 1;
        }
        assert!(degree.values().all(|&d| d == 2));
    }
    assert_eq!(boundary, expected);
}
#[test]
fn cancelled_invalid_and_resource_refusals_do_not_publish_or_mutate_source() {
    let source = small_source();
    let before = source;
    assert!(matches!(
        source.tessellate_shell(
            Tessellation {
                angle: 0.,
                ..options()
            },
            16_000_000,
            None
        ),
        Err(crate::authored_metric_thread::Error::Invalid)
    ));
    assert!(matches!(
        source.tessellate_shell(options(), 0, None),
        Err(crate::authored_metric_thread::Error::Limit)
    ));
    assert!(matches!(
        source.tessellate_shell(options(), 16_000_000, Some(&|| true)),
        Err(crate::authored_metric_thread::Error::Cancelled)
    ));
    let calls = std::cell::Cell::new(0);
    assert!(matches!(
        source.tessellate_shell(
            options(),
            16_000_000,
            Some(&|| {
                calls.set(calls.get() + 1);
                calls.get() > 5000
            })
        ),
        Err(crate::authored_metric_thread::Error::Cancelled)
    ));
    assert_eq!(source, before);
    for limit in [1, 100, 1000] {
        assert!(
            source
                .tessellate_shell(
                    Tessellation {
                        triangle_limit: limit,
                        ..options()
                    },
                    16_000_000,
                    None
                )
                .is_err()
        );
    }
    let shell = source
        .tessellate_shell(options(), 16_000_000, None)
        .unwrap();
    assert!(shell.exact_rim_turns(usize::MAX, 10000, None).is_err());
    assert!(shell.exact_rim_turns(0, 0, None).is_err());
    assert!(shell.exact_rim_turns(0, 10000, Some(&|| true)).is_err());
}
#[test]
fn trace_layout_preflight_counts_the_exact_owned_layout() {
    let maximum = 131072;
    let bytes = trace::Trace::requested_bytes(maximum).unwrap();
    let mut w = Work {
        remaining: 16_000_000,
        used: 0,
        interrupt: None,
    };
    let trace = trace::Trace::at_origin(maximum, 0., &mut w).unwrap();
    assert_eq!(trace.memory_bytes(), bytes);
    assert!(storage::bytes(&[(bytes, 1)], 0, bytes).is_ok());
    assert_eq!(
        storage::bytes(&[(bytes, 1)], 0, bytes - 1),
        Err(Error::Limit)
    );
    assert_eq!(
        trace::Trace::requested_bytes(maximum + 1),
        Err(Error::Limit)
    );
    assert_eq!(trace::Trace::requested_bytes(usize::MAX), Err(Error::Limit));
}

#[test]
fn bounded_grid_planning_refuses_extreme_public_quality_without_panicking() {
    let source = small_source();
    // Public options are checked before Grid::new: the root count is bounded
    // by ceil((2π/3)/.01)=210 before next_power_of_two is reached.
    assert!(matches!(
        source.tessellate_shell(
            Tessellation {
                angle: f64::from_bits(1),
                ..options()
            },
            16_000_000,
            None
        ),
        Err(crate::authored_metric_thread::Error::Invalid)
    ));
    assert!(matches!(
        source.tessellate_shell(
            Tessellation {
                angle: 0.01,
                ..options()
            },
            16_000_000,
            None
        ),
        Err(crate::authored_metric_thread::Error::Limit)
    ));

    assert!(matches!(
        source.tessellate_shell(
            Tessellation {
                deflection: f64::from_bits(1),
                ..options()
            },
            16_000_000,
            None
        ),
        Err(crate::authored_metric_thread::Error::Limit)
    ));
    let mut grid = Grid::new(&source, options()).unwrap();
    grid.advance(Reason::Repartition).unwrap();
    assert!(grid.odd);
    grid.advance(Reason::Repartition).unwrap();
    assert!(grid.shifted && !grid.odd);
    grid.advance(Reason::Repartition).unwrap();
    assert!(grid.odd);
    grid.advance(Reason::Repartition).unwrap();
    assert!(!grid.odd);
    grid.rows = sampling::MAX_ROWS;
    assert!(grid.advance(Reason::Repartition).is_err());
    grid.counts[2] = sampling::MAX_SUBDIVISIONS;
    assert!(grid.advance(Reason::Refine(2)).is_err());
    let mut w = Work {
        remaining: 16_000_000,
        used: 0,
        interrupt: None,
    };
    assert!(attempt(&source, options(), grid, &mut w).is_err());
    grid.rows = 64;
    grid.counts[2] = 0;
    assert!(attempt(&source, options(), grid, &mut w).is_err());
    assert!(reserve::<u8>(usize::MAX).is_err());
    for (from, to) in [
        (
            crate::authored_metric_thread::Error::Invalid,
            Error::Invalid,
        ),
        (crate::authored_metric_thread::Error::Limit, Error::Limit),
        (
            crate::authored_metric_thread::Error::Unresolved,
            Error::Unresolved,
        ),
        (
            crate::authored_metric_thread::Error::Cancelled,
            Error::Cancelled,
        ),
    ] {
        assert_eq!(source_error(from), to);
    }
}
#[test]
fn actual_capacity_and_literal_source_point_publication_have_separate_guards() {
    let mut v = Vec::<u32>::with_capacity(3);
    assert_eq!(room(v.len(), 3, v.capacity()), Ok(()));
    v.extend([1, 2, 3]);
    assert_eq!(room(v.len(), 1, v.capacity()), Err(Error::Limit));
    let source = small_source();
    let mut w = Work {
        remaining: 16_000_000,
        used: 0,
        interrupt: None,
    };
    let grid = Grid::new(&source, options()).unwrap();
    let patch = Patch::new(&source, 0, 2, grid, &mut w).unwrap();
    let strip = patch.strip(1, &mut w).unwrap();
    let polygon = patch.cell(&strip, 1, &mut w).unwrap();
    let point = patch.point(&polygon[0], &strip, &mut w).unwrap();
    let (_, error) = literal_point(point, 0.1).unwrap();
    assert!(error > 0.);
    assert_eq!(literal_point(point, error), Err(Error::Unresolved));
}
#[test]
fn actual_shifted_closed_phase_and_left_hand_rim_tokens_remain_source_owned() {
    let mut spec = small_source().spec();
    spec.left_hand = true;
    let source = AuthoredMetricThread::new(spec).unwrap();
    let grid = Grid {
        counts: [1, 1, 16, 1, 1],
        rows: 64,
        odd: false,
        shifted: true,
    };
    let mut w = Work {
        remaining: 16_000_000,
        used: 0,
        interrupt: None,
    };
    let shell = match attempt(&source, options(), grid, &mut w).unwrap() {
        Ok(shell) => shell,
        Err(_) => panic!("shifted actual source must qualify"),
    };
    for end in [ShellEnd::Start, ShellEnd::End] {
        for i in 0..shell.rim_count(end) {
            for vertex in shell.rim_edge(end, i).unwrap().vertex_indices {
                let (turn, used) = shell.exact_rim_turns(vertex, 10000, None).unwrap();
                assert!(
                    turn >= num_rational::BigRational::from_integer(0.into())
                        && turn < num_rational::BigRational::from_integer(1.into())
                );
                assert!(used > 0);
            }
        }
    }
}
#[test]
fn genuine_near_grid_cap_repartitions_without_dropping_its_positive_domain() {
    let original = source(false, false);
    let grid = Grid {
        counts: [1, 1, 16, 1, 1],
        rows: 64,
        odd: false,
        shifted: false,
    };
    let mut w = Work {
        remaining: 16_000_000,
        used: 0,
        interrupt: None,
    };
    let patch = Patch::new(&original, 1, 2, grid, &mut w).unwrap();
    let strip = patch.strip(2, &mut w).unwrap();
    let polygon = patch.cell(&strip, 37, &mut w).unwrap();
    let length = patch.point(&polygon[0], &strip, &mut w).unwrap()[2].midpoint();
    assert_eq!(length, 0.988_031_463_692_017_1);
    let mut spec = original.spec();
    spec.length_mm = length;
    let source = AuthoredMetricThread::new(spec).unwrap();
    let mut work = Work {
        remaining: 16_000_000,
        used: 0,
        interrupt: None,
    };
    assert!(matches!(
        attempt(&source, options(), grid, &mut work).unwrap(),
        Err(Reason::Repartition)
    ));
    let shell = source
        .tessellate_shell(options(), 16_000_000, None)
        .unwrap();
    assert_eq!(*shell.source(), source);
    assert!(shell.chord <= 0.1 && shell.normal <= 0.3);
}

#[test]
fn source_circle_scaling_qualifies_real_root8_facet_without_quality_change() {
    let source = source(false, false);
    let grid = Grid {
        counts: [1, 1, 8, 1, 1],
        rows: 64,
        odd: false,
        shifted: false,
    };
    let mut w = Work {
        remaining: 16_000_000,
        used: 0,
        interrupt: None,
    };
    let patch = Patch::new(&source, -1, 2, grid, &mut w).unwrap();
    let strip = patch.strip(0, &mut w).unwrap();
    let polygon = patch.cell(&strip, 39, &mut w).unwrap();
    let mut points = Vec::new();
    let mut errors = Vec::new();
    for p in &polygon {
        let (point, error) = literal_point(patch.point(p, &strip, &mut w).unwrap(), 0.1).unwrap();
        points.push(point);
        errors.push(error);
    }
    let field = patch.field(&polygon, &strip, &mut w).unwrap();
    let next = qualify_fan(&field, &points, &errors, options(), &mut w).unwrap();
    let old_field = {
        let _old = OldRoot::new();
        patch.field(&polygon, &strip, &mut w).unwrap()
    };
    let old = qualify_fan(&old_field, &points, &errors, options(), &mut w).unwrap();
    eprintln!(
        "scaledroot8new{} old{}",
        matches!(next, Attempt::Complete(_)),
        matches!(old, Attempt::Complete(_))
    );
    assert!(matches!(old, Attempt::Refine));
    assert!(matches!(next, Attempt::Complete(_)));
    assert!(matches!(
        qualify_fan(
            &field,
            &points,
            &errors,
            Tessellation {
                angle: 0.01,
                ..options()
            },
            &mut w
        )
        .unwrap(),
        Attempt::Refine
    ));
    let mut reverse = points.clone();
    reverse.reverse();
    let mut reverse_errors = errors.clone();
    reverse_errors.reverse();
    let proof = qualify_fan(&field, &reverse, &reverse_errors, options(), &mut w).unwrap();
    assert!(matches!(proof, Attempt::Complete(_)));
}

struct OldRoot(bool);
impl OldRoot {
    fn new() -> Self {
        Self(OLD_ROOT_NORMAL.replace(true))
    }
}
impl Drop for OldRoot {
    fn drop(&mut self) {
        OLD_ROOT_NORMAL.set(self.0);
    }
}
#[test]
fn paired_actual_shell_normal_conditioning_and_v32_proposal() {
    for (left, nb) in [(false, false), (true, true)] {
        let source = source(left, nb);
        let now = std::time::Instant::now();
        let old = {
            let _old = OldRoot::new();
            metric_shell(&source, options(), 16_000_000, None).unwrap()
        };
        eprintln!(
            "oldcircleLH{left}nb{nb}:{}facets,{}work,{:?}",
            old.mesh.face_indices.len(),
            old.work,
            now.elapsed()
        );
        let now = std::time::Instant::now();
        let new = metric_shell(&source, options(), 16_000_000, None).unwrap();
        eprintln!(
            "scaledcircleLH{left}nb{nb}:{}facets,{}work,{:?}",
            new.mesh.face_indices.len(),
            new.work,
            now.elapsed()
        );
        assert_eq!(*old.source(), *new.source());
        assert!(new.mesh.face_indices.len() < old.mesh.face_indices.len());
        assert!(new.work < old.work);
        let grid = Grid {
            counts: [1, 1, 8, 1, 1],
            rows: 32,
            odd: false,
            shifted: false,
        };
        let mut w = Work {
            remaining: 16_000_000,
            used: 0,
            interrupt: None,
        };
        let now = std::time::Instant::now();
        // Coarser angular sampling remains a real rejected proposal; positive
        // normal conditioning does not waive the requested facet certificate.
        assert!(matches!(
            attempt(&source, options(), grid, &mut w).unwrap(),
            Err(Reason::Refine(2))
        ));
        eprintln!(
            "V32scaledLH{left}nb{nb}:actual root refinement after{}work,{:?}",
            w.used,
            now.elapsed()
        );
    }
}

#[test]
fn consuming_shell_keeps_actual_earned_rim_errors_and_detached_mesh() {
    let source = small_source();
    let shell = source
        .tessellate_shell(options(), 16_000_000, None)
        .unwrap();
    for end in [ShellEnd::Start, ShellEnd::End] {
        for i in 0..shell.rim_count(end) {
            let edge = shell.rim_edge(end, i).unwrap();
            assert!(
                edge.point_error_mm
                    .into_iter()
                    .all(|e| e.is_finite() && (0. ..0.1).contains(&e))
            );
        }
    }
    let count = shell.mesh().face_indices.len();
    assert_eq!(shell.into_mesh().face_indices.len(), count);
}
