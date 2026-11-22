use super::{INVALID, MAX_TRIANGLES, MAX_VERTICES, certify, shape_budget};
use crate::display::{DisplayMesh, MeshError};

fn tetrahedron(origin: [f64; 3], size: f64, outward: bool) -> DisplayMesh {
    let points = [[0., 0., 0.], [size, 0., 0.], [0., size, 0.], [0., 0., size]];
    let mut triangles = [[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]];
    if !outward {
        for triangle in &mut triangles {
            triangle.swap(1, 2);
        }
    }
    DisplayMesh {
        vertices: points
            .into_iter()
            .flat_map(|p| std::array::from_fn::<_, 3, _>(|i| p[i] + origin[i]))
            .collect(),
        indices: triangles.into_iter().flatten().collect(),
        face_indices: vec![0, 1, 2, 3],
    }
}

fn append(target: &mut DisplayMesh, shell: DisplayMesh) {
    let offset = target.vertices.len() as u32 / 3;
    target.vertices.extend(shell.vertices);
    target
        .indices
        .extend(shell.indices.into_iter().map(|i| i + offset));
    target.face_indices.extend(shell.face_indices);
}

#[test]
fn separate_components_and_inward_cavity_preserve_signed_volume() {
    let mut mesh = tetrahedron([0.; 3], 8., true);
    append(&mut mesh, tetrahedron([1.; 3], 1., false));
    append(&mut mesh, tetrahedron([1e12, -1e12, 1e12], 4., true));
    let proof = certify(&mesh).unwrap();
    assert_eq!(proof.edge_pairs, 18);
    assert_eq!(proof.vertex_links, 12);
    assert!((proof.signed_volume_mm3 - (512. - 1. + 64.) / 6.).abs() < 1e-12);
}

#[test]
fn coincident_points_keep_separate_shell_indices_and_shared_fans_refuse() {
    let mut mesh = tetrahedron([0.; 3], 2., true);
    let mut touching = tetrahedron([0.; 3], 2., true);
    // Reflect through the origin and restore the outward winding.
    for point in touching.vertices.as_chunks_mut::<3>().0 {
        point[0] = -point[0];
        point[1] = -point[1];
        point[2] = -point[2];
    }
    for triangle in touching.indices.as_chunks_mut::<3>().0 {
        triangle.swap(1, 2);
    }
    append(&mut mesh, touching);
    let proof = certify(&mesh).unwrap();
    assert_eq!(proof.vertex_links, 8);
    assert!((proof.signed_volume_mm3 - 16. / 6.).abs() < 1e-12);
    let counted = super::certify_counted(&mesh, &mut |_, _| Ok::<_, ()>(())).unwrap();
    assert_eq!(counted, proof);
    for index in &mut mesh.indices {
        if *index == 4 {
            *index = 0;
        }
    }
    assert_eq!(certify(&mesh), Err(INVALID));
    // Both shells retain valid opposite edge uses, but identifying only the
    // touching vertex creates two fan cycles. Exercise the counted cycle
    // validator as well as the independent union-based implementation.
    let original = mesh.clone();
    assert_eq!(
        super::certify_counted(&mesh, &mut |_, _| Ok::<_, ()>(())),
        Err(super::CountedError::Mesh(INVALID))
    );
    assert_eq!(mesh, original);
}

#[test]
fn open_overused_and_misoriented_edges_are_rejected() {
    let original = tetrahedron([0.; 3], 1., true);
    let mut open = original.clone();
    open.indices.truncate(9);
    open.face_indices.truncate(3);
    assert_eq!(certify(&open), Err(INVALID));

    let mut overused = original.clone();
    overused.indices.extend_from_slice(&original.indices[..3]);
    overused.face_indices.push(4);
    assert_eq!(certify(&overused), Err(INVALID));

    let mut inconsistent = original;
    inconsistent.indices.swap(0, 1);
    assert_eq!(certify(&inconsistent), Err(INVALID));
}

#[test]
fn mirrored_duplicate_triangles_cannot_hide_in_an_outward_mesh() {
    let mut mesh = tetrahedron([0.; 3], 1., true);
    mesh.vertices.extend([2., 0., 0., 3., 0., 0., 2., 1., 0.]);
    mesh.indices.extend([4, 5, 6, 6, 5, 4]);
    mesh.face_indices.extend([4, 5]);
    assert_eq!(certify(&mesh), Err(INVALID));
}

#[test]
fn malformed_and_degenerate_meshes_fail_before_volume_integration() {
    let original = tetrahedron([0.; 3], 1., true);
    let mut bad = original.clone();
    bad.vertices.push(1.);
    assert_eq!(certify(&bad), Err(INVALID));
    bad = original.clone();
    bad.indices.push(1);
    assert_eq!(certify(&bad), Err(INVALID));
    bad = original.clone();
    bad.indices.clear();
    assert_eq!(certify(&bad), Err(INVALID));
    bad = original.clone();
    bad.face_indices.pop();
    assert_eq!(certify(&bad), Err(INVALID));
    bad = original.clone();
    bad.vertices[0] = f64::NAN;
    assert_eq!(certify(&bad), Err(INVALID));
    bad = original.clone();
    bad.indices[0] = 4;
    assert_eq!(certify(&bad), Err(INVALID));
    bad = original.clone();
    bad.indices[1] = bad.indices[0];
    assert_eq!(certify(&bad), Err(INVALID));
    bad = original;
    bad.vertices[0] = f64::MAX;
    bad.vertices[3] = -f64::MAX;
    assert_eq!(certify(&bad), Err(INVALID));
}

#[test]
fn a_nonpositive_or_unrepresentable_total_volume_is_rejected() {
    assert_eq!(certify(&tetrahedron([0.; 3], 2., false)), Err(INVALID));
    assert_eq!(certify(&tetrahedron([0.; 3], 1e120, true)), Err(INVALID));
}

#[test]
fn counts_are_bounded_before_any_certification_allocation() {
    assert_eq!(shape_budget(MAX_VERTICES, MAX_TRIANGLES), Ok(()));
    for (vertices, triangles) in [(MAX_VERTICES + 1, 1), (1, MAX_TRIANGLES + 1)] {
        assert_eq!(
            shape_budget(vertices, triangles),
            Err(MeshError::BudgetExceeded)
        );
    }
}

#[test]
fn counted_connectivity_preserves_valid_invalid_and_multi_shell_proofs() {
    let valid = tetrahedron([0.; 3], 8., true);
    let mut multi = valid.clone();
    append(&mut multi, tetrahedron([1.; 3], 1., false));
    append(&mut multi, tetrahedron([1e12, -1e12, 1e12], 4., true));
    let mut open = valid.clone();
    open.indices.truncate(9);
    open.face_indices.truncate(3);
    let mut wrong = valid.clone();
    wrong.indices.swap(0, 1);
    let mut degenerate = valid.clone();
    degenerate.indices[0] = degenerate.indices[1];
    let mut nonfinite = valid.clone();
    nonfinite.vertices[0] = f64::INFINITY;
    let mut badindex = valid.clone();
    badindex.indices[0] = u32::MAX;
    for mesh in [valid, multi, open, wrong, degenerate, nonfinite, badindex] {
        let original = mesh.clone();
        let mut work = 0;
        let mut peak = 0;
        let counted = super::certify_counted(&mesh, &mut |n, bytes| {
            work += n;
            peak = peak.max(bytes.unwrap_or(0));
            Ok::<_, ()>(())
        })
        .map_err(|error| match error {
            super::CountedError::Mesh(error) => error,
            super::CountedError::Interrupted(()) => panic!("ample budget"),
        });
        assert_eq!(counted, certify(&mesh));
        assert_eq!(mesh, original);
        assert!(work > 0);
        if counted.is_ok() {
            assert!(peak > 0);
        }
    }
}
#[test]
fn counted_connectivity_budget_cancel_storage_failure_and_clean_retry() {
    let mesh = tetrahedron([0.; 3], 8., true);
    let expected = certify(&mesh).unwrap();
    let mut total = 0;
    super::certify_counted(&mesh, &mut |n, _| {
        total += n;
        Ok::<_, ()>(())
    })
    .unwrap();
    for limit in [0, total / 2, total - 1] {
        let mut left = limit;
        let result = super::certify_counted(&mesh, &mut |n, _| {
            left = left.checked_sub(n).ok_or("budget")?;
            Ok(())
        });
        assert_eq!(result, Err(super::CountedError::Interrupted("budget")));
    }
    let mut calls = 0;
    super::certify_counted(&mesh, &mut |_, _| {
        calls += 1;
        Ok::<_, ()>(())
    })
    .unwrap();
    for checkpoint in [1, calls / 2, calls] {
        let mut seen = 0;
        assert_eq!(
            super::certify_counted(&mesh, &mut |_, _| {
                seen += 1;
                if seen == checkpoint {
                    Err("cancel")
                } else {
                    Ok(())
                }
            }),
            Err(super::CountedError::Interrupted("cancel"))
        );
    }
    assert_eq!(
        super::certify_counted(&mesh, &mut |_, storage| if storage.is_some() {
            Err("bytes")
        } else {
            Ok(())
        }),
        Err(super::CountedError::Interrupted("bytes"))
    );
    assert_eq!(
        super::certify_counted(&mesh, &mut |_, _| Ok::<_, ()>(())),
        Ok(expected)
    );
}

#[test]
fn public_uncounted_specialization_never_invokes_budget_callbacks() {
    let mesh = tetrahedron([0.; 3], 8., true);
    let mut never = |_: usize, _: Option<usize>| -> Result<(), ()> { panic!("uncounted callback") };
    assert_eq!(
        super::certify_impl::<_, false>(&mesh, &mut super::Budget(&mut never)),
        Ok(certify(&mesh).unwrap())
    );
}

#[test]
fn counted_vertex_cycles_match_independent_shells_cavities_and_unused_vertices() {
    let mut mesh = tetrahedron([0.; 3], 4., true);
    append(&mut mesh, tetrahedron([0.125; 3], 0.25, false));
    for i in 1..128 {
        let x = i as f64 * 8.;
        append(&mut mesh, tetrahedron([x, 0., 0.], 4., true));
        append(
            &mut mesh,
            tetrahedron([x + 0.125, 0.125, 0.125], 0.25, false),
        );
    }
    mesh.vertices.extend([9999.; 3]);
    let original = mesh.clone();
    let expected = certify(&mesh).unwrap();
    assert_eq!(expected.vertex_links, 128 * 8);
    let mut calls = 0;
    let mut last_storage = 0;
    let actual = super::certify_counted(&mesh, &mut |_, bytes| {
        calls += 1;
        if bytes.is_some() {
            last_storage = calls;
        }
        Ok::<_, ()>(())
    })
    .unwrap();
    assert_eq!(actual, expected);
    // The last allocation is the vertex-link table; interrupt after its
    // allocation and within traversal, after private links have changed.
    let mut seen = 0;
    assert_eq!(
        super::certify_counted(&mesh, &mut |_, _| {
            seen += 1;
            if seen == last_storage + 4 {
                Err("cycle cancel")
            } else {
                Ok(())
            }
        }),
        Err(super::CountedError::Interrupted("cycle cancel"))
    );
    assert_eq!(mesh, original);
    assert_eq!(
        super::certify_counted(&mesh, &mut |_, _| Ok::<_, ()>(())),
        Ok(expected)
    );
}
