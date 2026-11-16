use super::*;
use crate::step_read::phillips_cutter::physical::body::{
    common_mouth::tests::derived_thread_probe::with_head, derived_thread::actual_m3,
};
fn actual_owned() -> OwnedDegree11Material {
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
