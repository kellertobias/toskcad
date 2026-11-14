use super::*;
use std::hash::{Hash, Hasher};
fn source(left: bool, nonbinary: bool) -> Insert {
    let original = super::super::super::tests::fixture(left, 1, nonbinary);
    let mut outer = original.outer().spec();
    outer.crest_width_mm = 0.4;
    outer.lower_flank_deg = 30.;
    outer.upper_flank_deg = 30.;
    let mut inner = original.inner().spec();
    inner.left_hand = left;
    let options = crate::authored_flat_thread::SourceOptions {
        position_limit_mm: 1e-7,
        normal_angle_limit_rad: 1e-7,
    };
    Insert::new(
        AuthoredThread::new(outer, options).unwrap(),
        AuthoredThread::new(inner, options).unwrap(),
        0.5,
        None,
    )
    .unwrap()
}
fn options() -> Tessellation {
    Tessellation {
        deflection: 0.1,
        angle: 0.3,
        triangle_limit: 1_000_000,
    }
}
#[test]
fn actual_default_insert_candidate_has_closed_annular_material_and_paired_volume() {
    let body = source(false, false);
    let before = body;
    let time = std::time::Instant::now();
    let c = body.tessellate_candidate(options(), None).unwrap();
    eprintln!(
        "actual wood insert {}facets {:?} {}work source{:?} literal{:?}",
        c.mesh().face_indices.len(),
        time.elapsed(),
        c.work_used(),
        c.source_volume_mm3(),
        c.literal_volume_mm3()
    );
    let mut hash = std::collections::hash_map::DefaultHasher::new();
    for v in &c.mesh().vertices {
        v.to_bits().hash(&mut hash);
    }
    c.mesh().indices.hash(&mut hash);
    c.mesh().face_indices.hash(&mut hash);
    eprintln!("default actualmesh hash {:016x}", hash.finish());
    eprintln!(
        "candidate earned chord {} normal {}",
        c.source_chord_bound_mm(),
        c.source_normal_angle_bound_rad()
    );
    assert_eq!(*c.source(), before);
    assert!(c.source_chord_bound_mm() <= options().deflection);
    assert!(c.source_normal_angle_bound_rad() <= options().angle);
    let proof = crate::manufacturing::certify(c.mesh()).unwrap();
    assert!(proof.edge_pairs > 0 && proof.vertex_links > 0 && proof.signed_volume_mm3 > 0.);
    let v = c.literal_volume_mm3();
    assert!(v[0] <= proof.signed_volume_mm3 && proof.signed_volume_mm3 <= v[1]);
    let source = c.source_volume_mm3();
    assert!((proof.signed_volume_mm3 - source[0]).abs() < source[0] * 0.01);
    assert_eq!(body, before);
    let moved = c.into_mesh();
    assert!(crate::manufacturing::certify(&moved).is_ok());
}

#[test]
fn actual_left_hand_and_nonbinary_insert_materials_keep_source_and_closed_manifold() {
    for (left, nonbinary) in [(true, false), (false, true)] {
        let body = source(left, nonbinary);
        let timer = std::time::Instant::now();
        let c = body.tessellate_candidate(options(), None).unwrap();
        let proof = crate::manufacturing::certify(c.mesh()).unwrap();
        let literal = c.literal_volume_mm3();
        assert!(literal[0] <= proof.signed_volume_mm3 && proof.signed_volume_mm3 <= literal[1]);
        assert_eq!(*c.source(), body);
        assert!(c.source_chord_bound_mm() <= 0.1 && c.source_normal_angle_bound_rad() <= 0.3);
        eprintln!(
            "insert LH{left}/nonbinary{nonbinary}: {}facets {:?} {}work",
            c.mesh().face_indices.len(),
            timer.elapsed(),
            c.work_used()
        );
    }
}

fn budget() -> Work<'static> {
    Work {
        used: 0,
        limit: WORK,
        bytes: BYTES,
        interrupt: None,
    }
}
fn graph(angles: &[f64], links: &[[usize; 2]]) -> Vec<Option<Node>> {
    // Pure index/parameter protocol inputs; these are not geometry/source seals.
    let mut w = budget();
    angles
        .iter()
        .zip(links)
        .map(|(a, links)| {
            Some(Node {
                neighbors: *links,
                degree: 2,
                angle: Angle::encode(&binary(*a), &mut w).unwrap(),
            })
        })
        .collect()
}
#[test]
fn rim_index_parameter_protocol_rejects_branches_disconnected_and_wrong_winding() {
    let mut w = budget();
    let mut memory = Memory { used: SCRATCH };
    assert_eq!(endpoints([0, 0], 3), Err(Error::Unresolved));
    assert_eq!(endpoints([0, 3], 3), Err(Error::Unresolved));
    assert_eq!(endpoints([1, 2], 3), Ok([1, 2]));
    let mut n = Node {
        neighbors: [usize::MAX; 2],
        degree: 0,
        angle: Angle::encode(&Q::zero(), &mut w).unwrap(),
    };
    neighbor(&mut n, 1).unwrap();
    assert_eq!(neighbor(&mut n, 1), Err(Error::Unresolved));
    neighbor(&mut n, 2).unwrap();
    assert_eq!(neighbor(&mut n, 3), Err(Error::Unresolved));
    assert!(matches!(
        order_nodes(vec![], 0, &mut memory, &mut w),
        Err(Error::Unresolved)
    ));
    let mut chain = graph(&[0., 0.25, 0.5], &[[1, 2], [0, 2], [0, 1]]);
    chain[0].as_mut().unwrap().degree = 1;
    assert!(matches!(
        order_nodes(chain, 0, &mut memory, &mut w),
        Err(Error::Unresolved)
    ));
    let disconnected = graph(
        &[0., 0.125, 0.25, 0.5, 0.625, 0.75],
        &[[1, 2], [0, 2], [0, 1], [4, 5], [3, 5], [3, 4]],
    );
    assert!(matches!(
        order_nodes(disconnected, 0, &mut memory, &mut w),
        Err(Error::Unresolved)
    ));
    let reversed = graph(&[0., 0.5, 0.25, 0.75], &[[1, 3], [0, 2], [1, 3], [0, 2]]);
    assert!(matches!(
        order_nodes(reversed, 0, &mut memory, &mut w),
        Err(Error::Unresolved)
    ));
    let broken = graph(&[0., 0.25, 0.5, 0.75], &[[1, 2], [2, 3], [0, 1], [1, 2]]);
    assert!(matches!(
        order_nodes(broken, 0, &mut memory, &mut w),
        Err(Error::Unresolved)
    ));
    let good = graph(&[0., 0.25, 0.5], &[[1, 2], [0, 2], [0, 1]]);
    assert_eq!(
        order_nodes(good, 4, &mut memory, &mut w)
            .unwrap()
            .vertices
            .iter()
            .map(|n| n.index)
            .collect::<Vec<_>>(),
        vec![4, 5, 6]
    );
}
#[test]
fn finite_storage_preflight_and_typed_source_errors_preserve_failures() {
    let mut w = budget();
    let mut memory = Memory { used: BYTES - 1 };
    assert_eq!(memory.check(2), Err(Error::Limit));
    assert!(matches!(
        allocate::<u64>(usize::MAX, &mut memory, &mut w),
        Err(Error::Limit)
    ));
    for (e, expected) in [
        (SourceError::Invalid, Error::Invalid),
        (SourceError::Unsupported, Error::Unsupported),
        (SourceError::Unresolved, Error::Unresolved),
        (SourceError::Limit, Error::Limit),
        (SourceError::Cancelled, Error::Cancelled),
    ] {
        assert_eq!(failure(e), expected);
    }
    let body = source(false, false);
    let before = body;
    assert!(matches!(
        body.tessellate_candidate(
            Tessellation {
                deflection: 0.,
                ..options()
            },
            None
        ),
        Err(Error::Invalid)
    ));
    assert!(matches!(
        body.tessellate_candidate(
            Tessellation {
                triangle_limit: 1,
                ..options()
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        body.tessellate_candidate(
            Tessellation {
                triangle_limit: 6000,
                ..options()
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert_eq!(body, before);
}
#[test]
fn actual_full_insert_cancellation_never_publishes_partial_material() {
    use std::cell::Cell;
    let body = source(false, false);
    let before = body;
    let calls = Cell::new(0);
    let count = || {
        calls.set(calls.get() + 1);
        false
    };
    let c = body.tessellate_candidate(options(), Some(&count)).unwrap();
    let total = calls.get();
    let triangles = c.mesh().face_indices.len();
    for stop in [
        1,
        1000,
        total / 2,
        total - triangles - 14,
        total - triangles - 2,
        total,
    ] {
        let n = Cell::new(0);
        let cancel = || {
            n.set(n.get() + 1);
            n.get() == stop
        };
        assert!(
            matches!(
                body.tessellate_candidate(options(), Some(&cancel)),
                Err(Error::Cancelled)
            ),
            "checkpoint{stop}/{total}"
        );
        assert_eq!(body, before);
    }
    let mut inward = c.into_mesh();
    for t in inward.indices.as_chunks_mut::<3>().0 {
        t.swap(1, 2);
    }
    assert_eq!(
        literal_volume(&inward, &mut budget()),
        Err(Error::Unresolved)
    );
}

#[test]
fn actual_source_caps_reject_literal_tampering_and_accept_independent_periodic_cuts() {
    let body = source(false, false);
    let cap_source = body
        .cap_domains(
            CapEnd::Start,
            CapOptions {
                maximum_work: WORK,
                maximum_bytes: 1024 * 1024,
            },
            None,
        )
        .unwrap();
    let outer = body
        .outer()
        .tessellate_shell(options(), WORK, None)
        .unwrap();
    let inner = body
        .inner()
        .tessellate_shell(options(), WORK, None)
        .unwrap();
    let mut w = budget();
    let mut memory = Memory { used: SCRATCH };
    let mut outer_ring = ring(&outer, ShellEnd::Start, 0, &mut memory, &mut w).unwrap();
    let inner_ring = ring(
        &inner,
        ShellEnd::Start,
        outer.mesh().vertices.len() / 3,
        &mut memory,
        &mut w,
    )
    .unwrap();
    let mut vertices = outer.mesh().vertices.clone();
    vertices.extend(&inner.mesh().vertices);
    let empty = || DisplayMesh {
        vertices: vertices.clone(),
        indices: vec![],
        face_indices: vec![],
    };
    let mut wrong_plane = empty();
    wrong_plane.vertices[outer_ring.vertices[0].index * 3 + 2] += 0.125;
    assert_eq!(
        cap(
            &mut wrong_plane,
            &cap_source,
            &outer_ring,
            &inner_ring,
            options(),
            &mut w
        ),
        Err(Error::Unresolved)
    );
    let mut collapsed = empty();
    let a = outer_ring.vertices[0].index;
    let b = inner_ring.vertices[0].index;
    collapsed.vertices[b * 3] = collapsed.vertices[a * 3];
    collapsed.vertices[b * 3 + 1] = collapsed.vertices[a * 3 + 1];
    assert_eq!(
        cap(
            &mut collapsed,
            &cap_source,
            &outer_ring,
            &inner_ring,
            options(),
            &mut w
        ),
        Err(Error::Unresolved)
    );
    let mut limited = empty();
    assert_eq!(
        cap(
            &mut limited,
            &cap_source,
            &outer_ring,
            &inner_ring,
            Tessellation {
                triangle_limit: 1,
                ..options()
            },
            &mut w
        ),
        Err(Error::Limit)
    );
    // Same actual boundary/source points, a different legitimate periodic cut.
    // Wrapped vertices keep their identities; only exact proof parameters gain1.
    outer_ring.vertices.rotate_left(1);
    let last = outer_ring.vertices.len() - 1;
    let q = outer_ring.vertices[last].angle.decode(&mut w).unwrap() + Q::one();
    outer_ring.vertices[last].angle = Angle::encode(&q, &mut w).unwrap();
    let mut mesh = empty();
    let bound = cap(
        &mut mesh,
        &cap_source,
        &outer_ring,
        &inner_ring,
        options(),
        &mut w,
    )
    .unwrap();
    assert!(bound <= 0.1);
    assert_eq!(mesh.face_indices.len(), outer_ring.len() + inner_ring.len());
}

#[test]
fn actual_literal_plane_proof_cancels_every_coordinate_transfer_checkpoint() {
    use std::cell::Cell;
    let body = source(false, false);
    let mut w = budget();
    let cap = body
        .cap_domains(
            CapEnd::Start,
            CapOptions {
                maximum_work: WORK,
                maximum_bytes: 1024 * 1024,
            },
            None,
        )
        .unwrap();
    let uv = [
        [Q::zero(), Q::one()],
        [Q::zero(), Q::zero()],
        [Q::new(1.into(), 64.into()), Q::one()],
    ];
    let p = uv.each_ref().map(|v| {
        super::super::affine_probe::literal(cap.owner(), &Q::zero(), &v[0], &v[1], &mut w)
    });
    let calls = Cell::new(0);
    let count = || {
        calls.set(calls.get() + 1);
        false
    };
    let mut w = Work {
        interrupt: Some(&count),
        ..budget()
    };
    assert!(!plane_orientation(p, 0., &mut w).unwrap());
    let total = calls.get();
    for stop in 1..=total {
        let n = Cell::new(0);
        let cancel = || {
            n.set(n.get() + 1);
            n.get() == stop
        };
        let mut w = Work {
            interrupt: Some(&cancel),
            ..budget()
        };
        assert_eq!(plane_orientation(p, 0., &mut w), Err(Error::Cancelled));
    }
}
