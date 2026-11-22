use super::*;
use crate::step_read::phillips_cutter::physical::body::{
    common_mouth::tests::derived_thread_probe::with_head, derived_thread::actual_m3,
};
pub(super) fn actual_owned() -> OwnedDegree11Material {
    let text = actual_m3(true, false);
    let mut owned = None;
    with_head(&text, &mut |head| {
        let thread = head
            .degree11_thread_carrier(Default::default(), None)
            .unwrap()
            .certify_occupied(Default::default(), None)
            .unwrap();
        let material = head
            .closed_degree11_material(&thread, Default::default(), None)
            .unwrap();
        let volume = material.volume_bounds(Default::default(), None).unwrap();
        owned = Some(volume.owned_snapshot(Default::default(), None).unwrap());
    });
    owned.unwrap()
}
#[test]
fn actual_rounded_m3_owned_thread_head_and_cavity_have_coarse_source_certified_facets() {
    let source = actual_owned();
    let options = JoinedMeshOptions {
        tessellation: Tessellation {
            deflection: 0.5,
            angle: 0.7,
            triangle_limit: 1_000_000,
        },
        ..Default::default()
    };
    for bad in [
        JoinedMeshOptions {
            maximum_work: MESH_WORK + 1,
            ..options
        },
        JoinedMeshOptions {
            maximum_bytes: MESH_BYTES + 1,
            ..options
        },
    ] {
        assert!(matches!(
            source.joined_mesh_candidate(bad, None),
            Err(PhillipsReadError::Invalid)
        ));
    }
    for limited in [
        JoinedMeshOptions {
            maximum_work: 1,
            ..options
        },
        JoinedMeshOptions {
            maximum_bytes: source.retained_budget_bytes() - 1,
            ..options
        },
    ] {
        assert!(matches!(
            source.joined_mesh_candidate(limited, None),
            Err(PhillipsReadError::Limit)
        ));
    }
    assert!(matches!(
        source.joined_mesh_candidate(options, Some(&|| true)),
        Err(PhillipsReadError::Cancelled)
    ));
    let polls = std::cell::Cell::new(0usize);
    let cancel = || {
        polls.set(polls.get() + 1);
        polls.get() >= 32
    };
    assert!(matches!(
        source.joined_mesh_candidate(options, Some(&cancel)),
        Err(PhillipsReadError::Cancelled)
    ));
    let mesh = source.joined_mesh_candidate(options, None).unwrap();
    eprintln!(
        "D11_JOIN facets={} work={} peak={} combined={} error={}",
        mesh.mesh().face_indices.len(),
        mesh.work_used(),
        mesh.accounted_peak_bytes(),
        mesh.combined_live_upper_bound_bytes(),
        mesh.new_surface_error_mm()
    );
    assert!(std::ptr::eq(mesh.source(), &source));
    assert!(mesh.new_surface_error_mm() <= 0.5);
    assert!(mesh.work_used() <= MESH_WORK);
    assert!(mesh.accounted_peak_bytes() <= MESH_BYTES);
    assert!(!mesh.mesh().face_indices.is_empty());
    assert!(mesh.mesh().vertices.iter().all(|v| v.is_finite()));
    let proof = crate::manufacturing::certify(mesh.mesh()).unwrap();
    assert_eq!(proof, mesh.connectivity());
    for face in 0..source.face_count() as u32 {
        assert!(
            mesh.mesh().face_indices.contains(&face),
            "missing actual face{face}"
        );
    }
}

#[test]
fn actual_helical_chord_cancellation_prioritizes_v_without_sample_admission() {
    use crate::tapered_thread::NormalPolynomial11;
    let source = actual_owned();
    let root = source.thread_faces()[1].source().cell([0, 0]).unwrap();
    let field = NormalPolynomial11::new(&root).unwrap();
    let mut leaf = root.subdivide(0).unwrap()[0].clone();
    for _ in 0..3 {
        leaf = leaf.subdivide(1).unwrap()[0].clone();
    }
    assert_eq!(leaf.domain(), [[0., 0.5], [0., 0.125]]);
    let mut work = Work {
        remaining: MESH_WORK,
        held: source.retained_budget_bytes(),
        peak: 0,
        maximum_bytes: MESH_BYTES,
        maximum_triangles: 1_000_000,
        delegated: Default::default(),
        interrupt: None,
    };
    let axis = thread::preferred_axis(&leaf, &field, &mut work).unwrap();
    eprintln!(
        "D11_CHORD preferred_axis={axis} work={}",
        MESH_WORK - work.remaining
    );
    assert_eq!(axis, 1);
    let reduced = std::sync::Arc::new(
        crate::tapered_thread::ReducedNormal11::new(std::sync::Arc::clone(
            source.thread_normal(1, [0, 0]).unwrap(),
        ))
        .unwrap(),
    );
    let mut normal = reduced.root_cell().subdivide(0).unwrap()[0].clone();
    for _ in 0..3 {
        normal = normal.subdivide(1).unwrap()[0].clone();
    }
    for level in [3, 4] {
        let corners = leaf
            .corner_bounds()
            .map(|p| p.map(|[a, b]| I::new(a, b).midpoint()));
        let mut refusals = 0;
        for tri in [[0, 1, 2], [0, 2, 3]] {
            let p = tri.map(|i| corners[i].map(I::point));
            let n = super::super::geometry::cross(
                std::array::from_fn(|k| p[1][k] - p[0][k]),
                std::array::from_fn(|k| p[2][k] - p[0][k]),
            );
            let cone = normal
                .certifies_normal_cone(n.map(|v| [v.lo, v.hi]), 0.3)
                .unwrap();
            if level == 4 {
                assert!(cone);
            } else {
                refusals += usize::from(!cone);
            }
        }
        if level == 3 {
            assert!(refusals > 0);
        }
        leaf = leaf.subdivide(1).unwrap()[0].clone();
        normal = normal.subdivide(1).unwrap()[0].clone();
    }
}

#[test]
fn actual_rounded_m3_default_source_certified_facets_after_original_lattice() {
    let source = actual_owned();
    let candidate = source
        .joined_mesh_candidate(Default::default(), None)
        .unwrap();
    eprintln!(
        "D11_ORIGINAL_LATTICE_DEFAULT facets={} work={} remaining={} peak={} combined={} error={}",
        candidate.mesh().face_indices.len(),
        candidate.work_used(),
        candidate.remaining_work(),
        candidate.accounted_peak_bytes(),
        candidate.combined_live_upper_bound_bytes(),
        candidate.new_surface_error_mm()
    );
    assert!(candidate.new_surface_error_mm() <= 0.1);
    assert_eq!(candidate.certified_tessellation().angle, 0.3);
    assert_eq!(
        candidate.work_used() + candidate.remaining_work(),
        MESH_WORK
    );
    assert!(candidate.accounted_peak_bytes() <= MESH_BYTES);
    assert_eq!(
        crate::manufacturing::certify(candidate.mesh()).unwrap(),
        candidate.connectivity()
    );
    for limited in [
        JoinedMeshOptions {
            maximum_work: candidate.work_used() - 1,
            ..Default::default()
        },
        JoinedMeshOptions {
            maximum_bytes: candidate.accounted_peak_bytes() - 1,
            ..Default::default()
        },
        JoinedMeshOptions {
            tessellation: Tessellation {
                triangle_limit: candidate.mesh().face_indices.len() - 1,
                ..candidate.certified_tessellation()
            },
            ..Default::default()
        },
    ] {
        assert!(matches!(
            source.joined_mesh_candidate(limited, None),
            Err(PhillipsReadError::Limit)
        ));
    }
    let polls = std::cell::Cell::new(0usize);
    let cancel = || {
        polls.set(polls.get() + 1);
        polls.get() >= 4096
    };
    assert!(matches!(
        source.joined_mesh_candidate(Default::default(), Some(&cancel)),
        Err(PhillipsReadError::Cancelled)
    ));
    let retry = source
        .joined_mesh_candidate(Default::default(), None)
        .unwrap();
    assert_eq!(
        candidate
            .mesh()
            .vertices
            .iter()
            .map(|v| v.to_bits())
            .collect::<Vec<_>>(),
        retry
            .mesh()
            .vertices
            .iter()
            .map(|v| v.to_bits())
            .collect::<Vec<_>>()
    );
    assert_eq!(candidate.mesh().indices, retry.mesh().indices);
    assert_eq!(candidate.mesh().face_indices, retry.mesh().face_indices);
    assert_eq!(candidate.work_used(), retry.work_used());
}

#[test]
fn actual_prepared_head_matches_all_uncached_literal_bits_and_refusals() {
    let source = actual_owned();
    let quality = Tessellation {
        deflection: 0.5,
        angle: 0.7,
        triangle_limit: 1_000_000,
    };
    let budget = |held| Work {
        remaining: MESH_WORK,
        held,
        peak: 0,
        maximum_bytes: MESH_BYTES,
        maximum_triangles: 1_000_000,
        delegated: Default::default(),
        interrupt: None,
    };
    let mut initial = budget(source.retained_budget_bytes());
    let thread = thread::build(&source, quality, &mut initial).unwrap();
    let mut angular = Vec::new();
    for (curve, &divisions) in thread.curves.iter().zip(&thread.divisions) {
        for k in 0..divisions {
            angular.push(AngularCurve::<12> {
                poles: curve,
                spans: 1,
                reverse: false,
                range: [
                    k as f64 / divisions as f64,
                    (k + 1) as f64 / divisions as f64,
                ],
            });
        }
    }
    let mut counts = vec![1usize; source.head().view().segment_count()];
    let base = source.retained_budget_bytes()
        + mesh_bytes(&thread.mesh)
        + (thread.entry.capacity() + thread.tip.capacity()) * 4
        + thread.curves.capacity() * std::mem::size_of::<[P; 12]>()
        + thread.divisions.capacity() * std::mem::size_of::<usize>();
    for mode in 0..3 {
        let mut poor = budget(base);
        match mode {
            0 => poor.remaining = 1,
            1 => poor.maximum_bytes = base,
            2 => poor.interrupt = Some(&|| true),
            _ => unreachable!(),
        }
        let error = head_attempt_curve_cached(source.head(), &angular, &counts, quality, &mut poor);
        assert!(matches!(
            (mode, error),
            (0 | 1, Err(PhillipsReadError::Limit)) | (2, Err(PhillipsReadError::Cancelled))
        ));
    }
    for _ in 0..64 {
        let resident = base
            + angular.capacity() * std::mem::size_of::<AngularCurve<'_, 12>>()
            + counts.capacity() * std::mem::size_of::<usize>();
        let mut old_work = budget(resident);
        let old =
            head_attempt_curve(source.head(), &angular, &counts, quality, &mut old_work).unwrap();
        let extra = match &old {
            Ok(head) => mesh_bytes(&head.mesh) + (head.top.capacity() + head.bottom.capacity()) * 4,
            Err(_) => 0,
        };
        let mut new_work = budget(resident + extra);
        let new =
            head_attempt_curve_cached(source.head(), &angular, &counts, quality, &mut new_work)
                .unwrap();
        match (old, new) {
            (Ok(old), Ok(new)) => {
                assert_eq!(
                    old.mesh
                        .vertices
                        .iter()
                        .map(|x| x.to_bits())
                        .collect::<Vec<_>>(),
                    new.mesh
                        .vertices
                        .iter()
                        .map(|x| x.to_bits())
                        .collect::<Vec<_>>()
                );
                assert_eq!(old.mesh.indices, new.mesh.indices);
                assert_eq!(old.mesh.face_indices, new.mesh.face_indices);
                assert_eq!(old.top, new.top);
                assert_eq!(old.bottom, new.bottom);
                assert_eq!(old.distance.to_bits(), new.distance.to_bits());
                assert!(new_work.remaining > old_work.remaining);
                return;
            }
            (Err(Split::T(a)), Err(Split::T(b))) => {
                assert_eq!(a, b);
                counts[a] *= 2;
            }
            (Err(Split::U), Err(Split::U)) => {
                angular = angular
                    .iter()
                    .flat_map(|&a| {
                        let m = a.range[0] + (a.range[1] - a.range[0]) / 2.;
                        [
                            AngularCurve {
                                range: [a.range[0], m],
                                ..a
                            },
                            AngularCurve {
                                range: [m, a.range[1]],
                                ..a
                            },
                        ]
                    })
                    .collect();
            }
            _ => panic!("prepared source must preserve actual head refinement decision"),
        }
    }
    panic!("actual head must converge within the production attempt bound");
}
