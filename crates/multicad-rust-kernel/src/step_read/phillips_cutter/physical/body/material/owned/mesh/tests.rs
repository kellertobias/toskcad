use super::*;
#[test]
fn actual_cap_source_coarse_head_proposals_refine_and_refuse_nonprogress() {
    let source = super::super::mesh_probe::owned_neck(1., false, false, 0.4, 0.125);
    let quality = JoinedMeshOptions::default().tessellation;
    let trace = source
        .thread()
        .tessellate_with_cap_trace(quality, Default::default(), None)
        .unwrap();
    actual_trace_owner_chart_and_curve_family_are_earned_before_head_sampling(&source, &trace);
    let mut seed_work = Work {
        remaining: MESH_WORK,
        held: source.retained_budget_bytes()
            + trace.retained_metadata_bytes()
            + mesh_bytes(trace.mesh()),
        peak: 0,
        maximum_bytes: MESH_BYTES,
        maximum_triangles: quality.triangle_limit,
        delegated: Default::default(),
        interrupt: None,
    };
    let mut angular = angular_from_trace(source.thread(), &trace, &mut seed_work).unwrap();
    angular.retain(|a| a.range[0] == 0.);
    for a in &mut angular {
        a.range = [0., 1.];
    }
    assert_eq!(angular.len(), source.thread().side_faces().len());
    let initial = angular.clone();
    let counts = vec![1; source.head().view().segment_count()];
    let mut work = Work {
        remaining: MESH_WORK - trace.work_used(),
        held: source.retained_budget_bytes()
            + trace.retained_metadata_bytes()
            + mesh_bytes(trace.mesh())
            + angular.capacity() * std::mem::size_of::<Angular<'_>>(),
        peak: 0,
        maximum_bytes: MESH_BYTES,
        maximum_triangles: quality.triangle_limit,
        delegated: DelegatedWorkspaceBounds::default(),
        interrupt: None,
    };
    let coarse_count = angular.len();
    let head = make_head(&source, &mut angular, quality, &mut work).unwrap();
    assert!(angular.len() > coarse_count);
    assert!(head.distance <= quality.deflection);
    work.held += mesh_bytes(&head.mesh) + (head.top.capacity() + head.bottom.capacity()) * 4;
    actual_plane_proposals(&source, &trace, &head, &mut work);
    assembly_cancellation_matrix(&source, &trace, quality);
    let polls = std::cell::Cell::new(0usize);
    let observe = || {
        polls.set(polls.get() + 1);
        false
    };
    let mut probe = Work {
        remaining: MESH_WORK - trace.work_used(),
        held: work.held,
        peak: 0,
        maximum_bytes: MESH_BYTES,
        maximum_triangles: quality.triangle_limit,
        delegated: Default::default(),
        interrupt: Some(&observe),
    };
    assert!(
        head_attempt(
            &source,
            &initial,
            &vec![1; counts.len()],
            quality,
            &mut probe
        )
        .is_ok()
    );
    for stop in 0..polls.get() {
        let counter = std::cell::Cell::new(0usize);
        let cancel = || {
            let n = counter.get();
            counter.set(n + 1);
            n == stop
        };
        let mut probe = Work {
            remaining: MESH_WORK - trace.work_used(),
            held: work.held,
            peak: 0,
            maximum_bytes: MESH_BYTES,
            maximum_triangles: quality.triangle_limit,
            delegated: Default::default(),
            interrupt: Some(&cancel),
        };
        assert!(matches!(
            head_attempt(
                &source,
                &initial,
                &vec![1; counts.len()],
                quality,
                &mut probe
            ),
            Err(PhillipsReadError::Cancelled)
        ));
    }
    let mut tiny = quality;
    tiny.triangle_limit = 1;
    assert!(matches!(
        head_attempt(&source, &angular, &counts, tiny, &mut work),
        Err(PhillipsReadError::Limit)
    ));
    let first: P = std::array::from_fn(|k| head.mesh.vertices[k]);
    assert!(first.iter().all(|v| v.is_finite()));
    let mut adjacent = vec![Angular {
        range: [1f64.next_down(), 1.],
        ..angular[0]
    }];
    work.held += adjacent.capacity() * std::mem::size_of::<Angular<'_>>();
    let before = adjacent[0].range;
    assert!(matches!(
        refine_angular(&mut adjacent, &mut work),
        Err(PhillipsReadError::Unsupported)
    ));
    assert_eq!(adjacent[0].range, before);
    // Actual literal head facets provide positive/opposite normal controls.
    let ids: [u32; 3] = head.mesh.indices[..3].try_into().unwrap();
    let points: [P; 3] =
        ids.map(|id| std::array::from_fn(|k| head.mesh.vertices[3 * id as usize + k]));
    let p = points.map(|p| p.map(I::point));
    let n = geometry::cross(
        std::array::from_fn(|k| p[1][k] - p[0][k]),
        std::array::from_fn(|k| p[2][k] - p[0][k]),
    );
    assert!(normal_match(points, n, 0.3));
    assert!(!normal_match(points, n.map(|v| -v), 0.3));
    let mut lifted = DisplayMesh {
        vertices: Vec::new(),
        indices: Vec::new(),
        face_indices: Vec::new(),
    };
    let offset = append_mesh(&mut lifted, &head.mesh, &mut Some, 1e7, &mut work).unwrap();
    assert!(matches!(
        translation_transfer(&head.mesh, &lifted, offset, 1e7, 0., 0.01, &mut work),
        Err(PhillipsReadError::Unsupported)
    ));
    assert!(matches!(
        translation_transfer(&head.mesh, &lifted, offset, 1e7, 0.01, 1e-20, &mut work),
        Err(PhillipsReadError::Unsupported)
    ));
}

pub(in super::super) fn publication_frontiers(earned: &JoinedMaterialCandidate<'_>) {
    let owned = earned.clone();
    let earned = &owned;
    let prefix = earned.work_used();
    let bytes = earned.source().retained_budget_bytes() + mesh_bytes(earned.mesh()) + MESH_SCRATCH;
    for (remaining, maximum_bytes, success) in
        [(2, bytes, true), (1, bytes, false), (2, bytes - 1, false)]
    {
        let work = Work {
            remaining,
            held: earned.source().retained_budget_bytes(),
            peak: 0,
            maximum_bytes,
            maximum_triangles: earned.mesh().face_indices.len(),
            delegated: earned.delegated_workspace_bounds(),
            interrupt: None,
        };
        let result = publish(earned.clone(), prefix + remaining, work);
        if success {
            let result = result.unwrap();
            assert_eq!(result.mesh(), earned.mesh());
            assert_eq!(result.remaining_work(), 0);
            assert_eq!(result.work_used(), prefix + 2);
        } else {
            assert!(matches!(result, Err(PhillipsReadError::Limit)));
        }
    }
    let work = Work {
        remaining: 2,
        held: earned.source().retained_budget_bytes(),
        peak: 0,
        maximum_bytes: bytes,
        maximum_triangles: earned.mesh().face_indices.len() - 1,
        delegated: earned.delegated_workspace_bounds(),
        interrupt: None,
    };
    assert!(matches!(
        publish(earned.clone(), prefix + 2, work),
        Err(PhillipsReadError::Limit)
    ));
    for stop in [0, 1] {
        let calls = std::cell::Cell::new(0usize);
        let cancel = || {
            let n = calls.get();
            calls.set(n + 1);
            n == stop
        };
        let work = Work {
            remaining: 2,
            held: earned.source().retained_budget_bytes(),
            peak: 0,
            maximum_bytes: bytes,
            maximum_triangles: earned.mesh().face_indices.len(),
            delegated: earned.delegated_workspace_bounds(),
            interrupt: Some(&cancel),
        };
        assert!(matches!(
            publish(earned.clone(), prefix + 2, work),
            Err(PhillipsReadError::Cancelled)
        ));
    }
}

#[test]
fn actual_head_circle_near_semicircle_uses_conservative_sag_and_finite_fields() {
    let source = super::super::mesh_probe::owned(1., false);
    let radius = source.head().view().station(0).unwrap()[0];
    let alpha = std::f64::consts::FRAC_PI_2 - 1e-6;
    let a = [I::point(alpha).cos(), -I::point(alpha).sin(), I::point(0.)];
    let b = [I::point(alpha).cos(), I::point(alpha).sin(), I::point(0.)];
    let theta = I::new(-alpha, alpha);
    let gamma = [
        I::point(radius) * theta.cos(),
        I::point(radius) * theta.sin(),
        I::point(0.),
    ];
    let (sag, d) = geometry::minor_arc_sag(a, b, gamma).unwrap();
    assert!(sag.hi >= 1. - alpha.portable_cos() && sag.hi <= 1.);
    assert!(d[0].lo <= 1. && d[0].hi >= 1.);
    // This genuine actual-circle enclosure is deliberately too broad to earn
    // the minor branch. The production policy asks for U refinement, without
    // claiming that this oracle is an admitted retained-thread trace.
    let rejection = geometry::minor_arc_sag(a, b, [I::new(-radius, radius); 3]).unwrap_err();
    assert!(matches!(angular_rejection(rejection), Ok(Err(Split::U))));
    // Conservative fields for DISTINCT points of this same actual head circle
    // near its positive-Y cardinal direction. Enlarging only their upper-Y
    // bounds retains true unit rays, but cannot earn a finite dot-product proof.
    let beta = I::point(0.1);
    let a = [beta.sin(), I::new(beta.cos().lo, 1e169), I::point(0.)];
    let b = [-beta.sin(), I::new(beta.cos().lo, 1e169), I::point(0.)];
    let gamma = [
        I::new(-radius * beta.sin().hi, radius * beta.sin().hi),
        I::new(radius * beta.cos().lo, radius),
        I::point(0.),
    ];
    assert!(matches!(
        geometry::minor_arc_sag(a, b, gamma),
        Err(PhillipsReadError::Unsupported)
    ));
}
fn actual_trace_owner_chart_and_curve_family_are_earned_before_head_sampling(
    source: &OwnedDerivedMaterial,
    trace: &crate::radial_tensor::TracedMesh<'_>,
) {
    let quality = JoinedMeshOptions::default().tessellation;
    let mut work = Work {
        remaining: MESH_WORK,
        held: source.retained_budget_bytes(),
        peak: 0,
        maximum_bytes: MESH_BYTES,
        maximum_triangles: quality.triangle_limit,
        delegated: Default::default(),
        interrupt: None,
    };
    let clone = source.thread().clone();
    assert!(matches!(
        angular_from_trace(&clone, trace, &mut work),
        Err(PhillipsReadError::Unsupported)
    ));
    let angles = angular_from_trace(source.thread(), trace, &mut work).unwrap();
    let mut whole = angles[0].gamma([0., 1.], &mut work).unwrap();
    for a in angles.iter().filter(|a| a.range[0] == 0.) {
        let p = a.gamma([0., 1.], &mut work).unwrap();
        whole = std::array::from_fn(|k| geometry::hull(whole[k], p[k]));
    }
    assert!(geometry::unit_direction(whole).is_err());
    assert!(angles[0].gamma([-1., -0.5], &mut work).is_err());
    let circular = crate::radial_tensor::BoundaryOwnership::new(
        crate::radial_tensor::tests::cam::wedge(&[1, 2]),
        1e-7,
        None,
    )
    .unwrap()
    .bind_surfaces(1e-7, None)
    .unwrap()
    .certify_axial(1e-7, 1e-6, None)
    .unwrap()
    .certify_caps(None)
    .unwrap()
    .certify_ordering(None)
    .unwrap()
    .certify_material(1e-7, None)
    .unwrap();
    let trace = circular
        .tessellate_with_cap_trace(quality, Default::default(), None)
        .unwrap();
    assert!(matches!(
        angular_from_trace(&circular, &trace, &mut work),
        Err(PhillipsReadError::Unsupported)
    ));
    let placed = circular
        .transformed_with_interrupt(
            [[10., 20., 30.], [0.8, 0.6, 0.], [-0.36, 0.48, 0.8]],
            1e-6,
            1e-5,
            None,
        )
        .unwrap();
    let trace = placed
        .tessellate_with_cap_trace(quality, Default::default(), None)
        .unwrap();
    assert!(matches!(
        angular_from_trace(&placed, &trace, &mut work),
        Err(PhillipsReadError::Unsupported)
    ));
}

fn actual_plane_proposals(
    source: &OwnedDerivedMaterial,
    trace: &crate::radial_tensor::TracedMesh<'_>,
    head: &Head,
    work: &mut Work<'_>,
) {
    let mut target = DisplayMesh {
        vertices: Vec::new(),
        indices: Vec::new(),
        face_indices: Vec::new(),
    };
    let outer_offset = append_mesh(&mut target, &head.mesh, &mut Some, 0., work).unwrap();
    let inner_offset = append_mesh(
        &mut target,
        trace.mesh(),
        &mut |f| if f == 0 { None } else { Some(f) },
        0.,
        work,
    )
    .unwrap();
    let outer = head
        .bottom
        .iter()
        .map(|&i| i + outer_offset)
        .collect::<Vec<_>>();
    let inner = trace.cap_rims()[0]
        .iter()
        .map(|p| p.vertex + inner_offset)
        .collect::<Vec<_>>();
    let station = source
        .head()
        .view()
        .station(source.head().view().segment_count())
        .unwrap()[1];
    let delta = trace.cap_rims()[0]
        .iter()
        .map(|p| p.outgoing_chord_error_mm)
        .fold(head.distance, f64::max);
    let options = JoinedMeshOptions::default().tessellation;
    let saved = work.maximum_triangles;
    work.maximum_triangles = target.face_indices.len();
    assert!(matches!(
        append_mesh(&mut target, trace.mesh(), &mut Some, 0., work),
        Err(PhillipsReadError::Limit)
    ));
    assert!(matches!(
        plane(
            &mut target,
            &outer,
            &inner,
            false,
            PlaneProof {
                face: 0,
                positive: true,
                source_z: station,
                boundary_error: delta,
                options
            },
            work
        ),
        Err(PhillipsReadError::Limit)
    ));
    work.maximum_triangles = saved;
    assert!(matches!(
        plane(
            &mut target,
            &outer,
            &inner,
            true,
            PlaneProof {
                face: 0,
                positive: true,
                source_z: station,
                boundary_error: delta,
                options
            },
            work
        ),
        Err(PhillipsReadError::Unsupported)
    ));
    let mut tight = options;
    tight.deflection = delta / 2.;
    assert!(matches!(
        plane(
            &mut target,
            &outer,
            &inner,
            false,
            PlaneProof {
                face: 0,
                positive: true,
                source_z: station,
                boundary_error: delta,
                options: tight
            },
            work
        ),
        Err(PhillipsReadError::Unsupported)
    ));
    // Actual finite candidate vertices are deliberately tilted while their
    // source stays fixed: distance still fits, but the plane normal does not.
    for p in target.vertices.as_chunks_mut::<3>().0 {
        p[2] += 0.05 * p[0];
    }
    tight = options;
    tight.angle = 0.01;
    assert!(matches!(
        plane(
            &mut target,
            &outer,
            &inner,
            false,
            PlaneProof {
                face: 0,
                positive: true,
                source_z: station,
                boundary_error: delta,
                options: tight
            },
            work
        ),
        Err(PhillipsReadError::Unsupported)
    ));
}

// Calibrated once against real scalar assembly callbacks, with no production
// stage instrumentation. Relative locations keep cancellation portable when
// an equivalent platform source changes callback totals. The independent
// source/normal/manifold and recovery checks still execute on every host.
fn assembly_cancellation_matrix(
    source: &OwnedDerivedMaterial,
    trace: &crate::radial_tensor::TracedMesh<'_>,
    options: Tessellation,
) {
    let mut prefix = Work {
        remaining: MESH_WORK - trace.work_used(),
        held: source.retained_budget_bytes()
            + trace.retained_metadata_bytes()
            + mesh_bytes(trace.mesh()),
        peak: 0,
        maximum_bytes: MESH_BYTES,
        maximum_triangles: options.triangle_limit,
        delegated: Default::default(),
        interrupt: None,
    };
    let q = source_quality(options).unwrap();
    let wrong = source_quality(Tessellation {
        angle: options.angle / 2.,
        ..options
    })
    .unwrap();
    assert!(matches!(
        components(source, trace, wrong, &mut prefix),
        Err(PhillipsReadError::Unsupported)
    ));
    let other = super::super::mesh_probe::owned_neck(1., true, false, 0.4, 0.125);
    assert!(matches!(
        components(&other, trace, q, &mut prefix),
        Err(PhillipsReadError::Unsupported)
    ));
    let parts = components(source, trace, q, &mut prefix).unwrap();
    let polls = std::cell::Cell::new(0usize);
    let observe = || {
        polls.set(polls.get() + 1);
        false
    };
    let mut w = prefix.clone();
    w.interrupt = Some(&observe);
    let baseline = assemble(&parts, MESH_WORK, w).unwrap();
    let total = polls.get();
    // Actual canonical checkpoints span owner/array preparation, source head,
    // each physical component copy, literal transfer, both joining faces,
    // each plane's preallocation, compaction and final publication.
    for mark in [
        0usize, 90, 91, 86076, 86077, 103241, 108611, 152079, 152080, 195546, 238266, 238267,
        238269, 238447, 238449, 239667, 372005, 372006, 372007,
    ] {
        let stop = (mark as u128 * total as u128 / 372008) as usize;
        let counter = std::cell::Cell::new(0usize);
        let cancel = || {
            let n = counter.get();
            counter.set(n + 1);
            n == stop
        };
        let mut w = prefix.clone();
        w.interrupt = Some(&cancel);
        assert!(matches!(
            assemble(&parts, MESH_WORK, w),
            Err(PhillipsReadError::Cancelled)
        ));
    }
    let retry = assemble(&parts, MESH_WORK, prefix).unwrap();
    assert_eq!(baseline.mesh(), retry.mesh());
}

#[test]
fn actual_cap_three_span_parameter_uses_literal_rounded_knots() {
    use num_rational::BigRational as Q;
    let source = super::super::mesh_probe::owned(1., false);
    let input = source.thread().input();
    // Coalesce THREE consecutive actual retained cap Bernstein pieces for a
    // numerical parameterization control. All original coefficients/shared
    // endpoints remain unchanged; no Owned material or new admission is forged.
    let mut poles = Vec::new();
    'edges: for &(index, forward) in &input.faces[0].outer {
        let edge = input.edges[index];
        let crate::radial_tensor::Curve::Spline {
            spans,
            poles: actual,
        } = &input.curves[edge.curve]
        else {
            panic!("actual derived cap is spline-owned")
        };
        let mut actual = actual.clone();
        if edge.forward != forward {
            actual.reverse();
        }
        for span in 0..*spans {
            let panel = &actual[9 * span..=9 * span + 9];
            if poles.is_empty() {
                poles.extend(panel);
            } else {
                assert_eq!(poles.last(), panel.first());
                poles.extend(&panel[1..]);
            }
            if poles.len() == 28 {
                break 'edges;
            }
        }
    }
    assert_eq!(poles.len(), 28);
    let knots: [f64; 4] = [0., 1. / 3., 2. / 3., 1.];
    let angular = Angular {
        poles: &poles,
        spans: 3,
        reverse: false,
        range: [0., 1.],
    };
    let mut work = Work {
        remaining: MESH_WORK,
        held: source.retained_budget_bytes(),
        peak: 0,
        maximum_bytes: MESH_BYTES,
        maximum_triangles: 1_000_000,
        delegated: Default::default(),
        interrupt: None,
    };
    let exact = |x: f64| Q::from_float(x).unwrap();
    for u in [
        knots[1].next_down(),
        knots[1],
        knots[1].next_up(),
        knots[2].next_down(),
        knots[2],
        knots[2].next_up(),
    ] {
        let span = knots.partition_point(|&x| x <= u).saturating_sub(1).min(2);
        let t = (exact(u) - exact(knots[span])) / (exact(knots[span + 1]) - exact(knots[span]));
        let bounds = angular.gamma([u, u], &mut work).unwrap();
        for (axis, bound) in bounds.iter().enumerate() {
            let mut values: Vec<Q> = poles[9 * span..=9 * span + 9]
                .iter()
                .map(|p| exact(p[axis]))
                .collect();
            for level in (1..10).rev() {
                for i in 0..level {
                    values[i] = (Q::from_integer(1.into()) - &t) * &values[i] + &t * &values[i + 1];
                }
            }
            assert!(exact(bound.lo) <= values[0] && values[0] <= exact(bound.hi));
        }
    }
}

#[test]
fn actual_public_producer_cancellation_tracks_component_array_frontiers() {
    // Derive locations from THIS host's actual certified components. The old
    // fraction of one historical callback total drifts when a source mesher
    // changes its facet count; no synthetic body or production stage hook is
    // needed to exercise the public caller's real copy/join refusal paths.
    let source = super::super::mesh_probe::owned_neck(1., false, false, 0.4, 0.125);
    let original = source.original_input().to_owned();
    let query = JoinedMeshOptions {
        tessellation: Tessellation {
            deflection: 0.2,
            angle: 0.5,
            ..JoinedMeshOptions::default().tessellation
        },
        ..Default::default()
    };
    let polls = std::cell::Cell::new(0usize);
    let observe = || {
        polls.set(polls.get() + 1);
        false
    };
    // Reproduce only the public PREPARATION contract to retain its earned
    // components. Budget and callback conditions are identical to the caller;
    // equal final meshes alone do not imply equal construction callbacks.
    let mut prefix = Work {
        remaining: MESH_WORK,
        held: source.retained_budget_bytes(),
        peak: 0,
        maximum_bytes: MESH_BYTES,
        maximum_triangles: query.tessellation.triangle_limit,
        delegated: Default::default(),
        interrupt: Some(&observe),
    };
    prefix.storage(0).unwrap();
    prefix
        .charge(source.thread().input().surfaces.len())
        .unwrap();
    prefix.delegated.retained_thread_bytes =
        delegated::retained(source.thread(), query.tessellation.triangle_limit).unwrap();
    let trace = source
        .thread()
        .tessellate_with_cap_trace_budget(
            query.tessellation,
            Default::default(),
            prefix.remaining,
            prefix.interrupt,
        )
        .unwrap();
    prefix.charge(trace.work_used()).unwrap();
    prefix.held += mesh_bytes(trace.mesh()) + trace.retained_metadata_bytes();
    prefix.storage(0).unwrap();
    let parts = components(
        &source,
        &trace,
        source_quality(query.tessellation).unwrap(),
        &mut prefix,
    )
    .unwrap();
    let public_prefix = polls.get();
    let assembly = assemble(&parts, MESH_WORK, prefix.clone()).unwrap();
    let total_manual = polls.get();
    let assembly_polls = total_manual - public_prefix;
    polls.set(0);
    let public = source.joined_mesh_candidate(query, Some(&observe)).unwrap();
    let public_polls = polls.get();
    assert_eq!(public.mesh(), assembly.mesh());
    assert_eq!(public_polls, total_manual);
    // Observe the real angular/head preparation once. The remaining copy and
    // transfer routines poll once per actual vertex/facet and twice per copy
    // reservation; those source-array counts avoid guessed global ratios.
    polls.set(0);
    let mut preparation = prefix.clone();
    preparation.interrupt = Some(&observe);
    let mut angular = angular_from_trace(source.thread(), &trace, &mut preparation).unwrap();
    preparation.held += angular.capacity() * std::mem::size_of::<Angular<'_>>();
    preparation.storage(0).unwrap();
    let head = make_head(&source, &mut angular, query.tessellation, &mut preparation).unwrap();
    preparation.held += mesh_bytes(&head.mesh) + (head.top.capacity() + head.bottom.capacity()) * 4;
    preparation.storage(0).unwrap();
    let before_copies = public_prefix + polls.get();
    let copy_polls = |mesh: &DisplayMesh| 2 + mesh.vertices.len() / 3 + mesh.face_indices.len();
    let before_cavity = before_copies + copy_polls(trace.mesh()) + copy_polls(&head.mesh);
    let cavity = parts.cavity.boundary_mesh().mesh();
    let after_cavity = before_cavity + copy_polls(cavity);
    let transfer_polls = cavity.vertices.len() / 3 + cavity.face_indices.len();
    let mouth_faces = parts.atlas.topology().faces();
    let kept_thread = trace
        .mesh()
        .face_indices
        .iter()
        .filter(|&&f| f != 0)
        .count();
    let kept_cavity = cavity
        .face_indices
        .iter()
        .filter(|&&f| mouth_faces[f as usize].source_identity() != CavityFace::Mouth)
        .count();
    let reversed_polls = kept_thread + head.mesh.face_indices.len() + kept_cavity;
    let before_wire_storage = after_cavity + 1 + transfer_polls + reversed_polls;
    let before_join = before_wire_storage + 1;
    let join_face =
        (source.head().view().segment_count() + source.thread().side_faces().len() + 1) as u32;
    let join_triangles = public
        .mesh()
        .face_indices
        .iter()
        .filter(|&&f| f == join_face)
        .count();
    assert!(join_triangles > 0);
    // Plane: one source-boundary charge, one pre-CDT allocation check, then
    // the actual post-CDT reservation. Its triangle count comes from the
    // successful PUBLIC mesh, not a formula that assumes a polygon shape.
    let before_top = before_join + 4 + join_triangles;
    let marks = [
        ("cavity copy", before_cavity),
        ("cavity copy work handoff", after_cavity),
        ("literal cavity translation", after_cavity + 1),
        ("joining wire storage", before_wire_storage),
        ("actual join post-CDT storage", before_join + 2),
        ("actual top plane", before_top),
    ];
    for (phase, stop) in marks {
        assert!(stop < public_polls);
        let calls = std::cell::Cell::new(0usize);
        let cancel = || {
            let at = calls.get();
            calls.set(at + 1);
            at == stop
        };
        assert!(
            matches!(
                source.joined_mesh_candidate(query, Some(&cancel)),
                Err(PhillipsReadError::Cancelled)
            ),
            "{phase}"
        );
        assert_eq!(calls.get(), stop + 1, "{phase}");
        assert_eq!(source.original_input(), original, "{phase}");
    }
    let recovered = source.joined_mesh_candidate(query, None).unwrap();
    assert!(std::ptr::eq(recovered.source(), &source));
    assert_eq!(public.mesh(), recovered.mesh());
    assert_eq!(
        public.certified_tessellation().deflection,
        recovered.certified_tessellation().deflection
    );
    eprintln!(
        "PUBLIC JOIN FRONTIERS public{public_polls} assembly{assembly_polls} prefix{public_prefix} marks{marks:?}"
    );
}

#[test]
fn exact_planar_disk_and_annulus_keep_distinct_hole_counts_and_orientation() {
    for hole in [false, true] {
        for positive in [false, true] {
            let mut mesh = DisplayMesh {
                vertices: vec![-2., -2., 0., 2., -2., 0., 2., 2., 0., -2., 2., 0.],
                indices: vec![],
                face_indices: vec![],
            };
            if hole {
                mesh.vertices
                    .extend([-1., -1., 0., 1., -1., 0., 1., 1., 0., -1., 1., 0.]);
            }
            let mut work = Work {
                remaining: MESH_WORK,
                held: 0,
                peak: 0,
                maximum_bytes: MESH_BYTES,
                maximum_triangles: 1_000_000,
                delegated: Default::default(),
                interrupt: None,
            };
            let inner: &[u32] = if hole { &[4, 5, 6, 7] } else { &[] };
            let error = plane(
                &mut mesh,
                &[0, 1, 2, 3],
                inner,
                false,
                PlaneProof {
                    face: 7,
                    positive,
                    source_z: 0.,
                    boundary_error: 0.,
                    options: JoinedMeshOptions::default().tessellation,
                },
                &mut work,
            )
            .unwrap();
            assert!(error <= 0.1);
            assert_eq!(mesh.face_indices.len(), if hole { 8 } else { 2 });
            assert!(mesh.face_indices.iter().all(|&f| f == 7));
            let area = mesh
                .indices
                .as_chunks::<3>()
                .0
                .iter()
                .map(|ids| {
                    let p = ids.map(|i| {
                        [
                            mesh.vertices[3 * i as usize],
                            mesh.vertices[3 * i as usize + 1],
                        ]
                    });
                    ((p[1][0] - p[0][0]) * (p[2][1] - p[0][1])
                        - (p[1][1] - p[0][1]) * (p[2][0] - p[0][0]))
                        * 0.5
                })
                .sum::<f64>();
            assert_eq!(
                area,
                (if hole { 12. } else { 16. }) * (if positive { 1. } else { -1. })
            );
        }
    }
}
