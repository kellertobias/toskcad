use super::*;
use std::cell::Cell;

fn owner() -> Rc<Placement> {
    let solid = Solid::rust_retained_helical_thread_with_interrupt(
        crate::HelicalThread {
            major_radius_mm: 4.,
            pitch_mm: 1.25,
            length_mm: 1.25,
            depth_mm: 0.5,
            crest_width_mm: 0.25,
            lower_flank_deg: 12.,
            upper_flank_deg: 18.,
            radial_taper: 0.,
            starts: 1,
            left_hand: false,
        },
        EntityId::new(),
        PlaneFrame::default(),
        GeometryPolicy::default(),
        None,
    )
    .unwrap();
    Rc::clone(&solid.radial_boundary_body().unwrap().0.source)
}
fn sibling(source: &Rc<Placement>) -> Rc<Placement> {
    Rc::new(Placement {
        boundary: Rc::clone(&source.boundary),
        properties: RefCell::new(None),
        volume: RefCell::new(None),
        bounds: RefCell::new(None),
    })
}
fn options() -> Tessellation {
    Tessellation {
        deflection: 0.1,
        angle: 0.3,
        triangle_limit: 1_000_000,
    }
}
// Opaque detached arrays exercise storage only. These private cache controls
// do not manufacture or certify a solid from this triangle.
fn triangle() -> Mesh {
    Mesh {
        vertices: vec![0., 0., 0., 1., 0., 0., 0., 1., 0.],
        indices: vec![0, 1, 2],
        face_indices: vec![0],
    }
}

fn bounded(payload: usize, entries: usize, receipt: usize) -> Cache {
    Cache::new(payload + entries * size_of::<Entry>(), entries, receipt)
}

#[test]
fn cache_uses_live_source_identity_all_quality_bits_and_detached_results() {
    let source = owner();
    let other = sibling(&source);
    let mesh = triangle();
    let key = Quality::from(options());
    let mut cache = Cache::default();
    cache.remember(&source, key, &mesh);
    assert_eq!(cache.get(&Rc::clone(&source), key), Some(mesh.clone()));
    assert!(cache.get(&other, key).is_none());
    for altered in [
        Tessellation {
            deflection: 0.2,
            ..options()
        },
        Tessellation {
            angle: 0.4,
            ..options()
        },
        Tessellation {
            triangle_limit: 999,
            ..options()
        },
    ] {
        assert!(cache.get(&source, Quality::from(altered)).is_none());
    }
    assert_ne!(
        Quality::from(Tessellation {
            deflection: 0.,
            ..options()
        }),
        Quality::from(Tessellation {
            deflection: -0.,
            ..options()
        })
    );
    assert_ne!(
        Quality::from(Tessellation {
            angle: 0.,
            ..options()
        }),
        Quality::from(Tessellation {
            angle: -0.,
            ..options()
        })
    );
    let mut detached = cache.get(&source, key).unwrap();
    detached.vertices[0] = 123.;
    assert_eq!(cache.get(&source, key), Some(mesh.clone()));
    let replacement = Quality::from(Tessellation {
        angle: 0.4,
        ..options()
    });
    cache.remember(&source, replacement, &mesh);
    assert_eq!(cache.entries.len(), 1);
    assert!(cache.get(&source, key).is_none());
    assert_eq!(cache.get(&source, replacement), Some(mesh));
}

#[test]
fn aggregate_byte_entry_and_receipt_limits_evict_lru_and_release_expired_sources() {
    let a = owner();
    let b = sibling(&a);
    let c = sibling(&a);
    let mesh = triangle();
    let size = payload_bytes(&mesh).unwrap();
    let key = Quality::from(options());
    let mut cache = bounded(size * 2, 10, size);
    let queue = cache.bytes;
    cache.remember(&a, key, &mesh);
    cache.remember(&b, key, &mesh);
    assert!(cache.get(&a, key).is_some());
    cache.remember(&c, key, &mesh);
    assert!(cache.get(&b, key).is_none());
    assert_eq!(cache.bytes, queue + size * 2);
    assert!(cache.get(&a, key).is_some());
    drop(c);
    cache.expire();
    assert_eq!(cache.entries.len(), 1);
    assert_eq!(cache.bytes, queue + size);
    let mut entries = bounded(size * 10, 1, size);
    let queue = entries.bytes;
    entries.remember(&a, key, &mesh);
    entries.remember(&b, key, &mesh);
    assert!(entries.get(&a, key).is_none());
    assert_eq!(entries.bytes, queue + size);
    for limits in [
        (size - 1, 10, size),
        (size * 10, 10, size - 1),
        (size * 10, 0, size),
    ] {
        let mut bounded = bounded(limits.0, limits.1, limits.2);
        let queue = bounded.bytes;
        bounded.remember(&a, key, &mesh);
        bounded.remember_owned(&a, key, mesh.clone());
        assert_eq!(bounded.bytes, queue);
    }
    let mut large_capacity = mesh.clone();
    large_capacity.vertices.reserve(128);
    let mut bounded = bounded(size * 2, 10, size);
    bounded.remember(&a, key, &large_capacity);
    assert!(bounded.entries.is_empty());
    assert!(payload_capacities([usize::MAX, 0, 0]).is_none());
    assert!(payload_capacities([0, usize::MAX, 0]).is_none());
    assert!(payload_capacities([0, 0, usize::MAX]).is_none());
    assert!(payload_capacities([usize::MAX / 8, usize::MAX / 4, 0]).is_none());
}

#[test]
fn manufacturing_always_builds_and_failed_or_invalid_work_does_not_replace_a_display() {
    DISPLAY.with(|cache| *cache.borrow_mut() = Cache::default());
    let source = owner();
    let builds = Cell::new(0);
    let build = || {
        builds.set(builds.get() + 1);
        Ok(triangle())
    };
    mesh(&source, options(), MeshUse::Display, None, build).unwrap();
    mesh(&source, options(), MeshUse::Display, None, build).unwrap();
    assert_eq!(builds.get(), 1);
    mesh(&source, options(), MeshUse::Manufacturing, None, build).unwrap();
    mesh(&source, options(), MeshUse::Manufacturing, None, build).unwrap();
    assert_eq!(builds.get(), 3);
    let alternate = Tessellation {
        angle: 0.4,
        ..options()
    };
    let cancelled = mesh(&source, alternate, MeshUse::Display, None, || {
        Err(error("test", GeometryErrorKind::Cancelled, "superseded"))
    })
    .unwrap_err();
    assert_eq!(cancelled.kind, GeometryErrorKind::Cancelled);
    let invalid = mesh(&source, alternate, MeshUse::Display, None, || {
        Ok(Mesh {
            indices: vec![99, 99, 99],
            ..triangle()
        })
    })
    .unwrap_err();
    assert_eq!(invalid.kind, GeometryErrorKind::InvalidResult);
    mesh(&source, options(), MeshUse::Display, None, build).unwrap();
    assert_eq!(builds.get(), 3);
    mesh(&source, alternate, MeshUse::Display, None, build).unwrap();
    assert_eq!(builds.get(), 4);
}

#[test]
fn actual_source_meshes_keep_pose_quality_and_manufacturing_proofs_on_cached_requests() {
    let feature = EntityId::new();
    let build = |radius| {
        Solid::rust_retained_helical_thread_with_interrupt(
            crate::HelicalThread {
                major_radius_mm: radius,
                pitch_mm: 1.25,
                length_mm: 1.25,
                depth_mm: 0.5,
                crest_width_mm: 0.25,
                lower_flank_deg: 12.,
                upper_flank_deg: 18.,
                radial_taper: 0.,
                starts: 1,
                left_hand: false,
            },
            feature,
            PlaneFrame::default(),
            GeometryPolicy::default(),
            None,
        )
        .unwrap()
    };
    let source = build(4.);
    let options = MeshOptions {
        linear_deflection_mm: 0.2,
        angular_deflection_rad: 0.6,
        maximum_triangles: 1_000_000,
    };
    let first = source.mesh(options).unwrap();
    assert_eq!(source.share().unwrap().mesh(options).unwrap(), first);
    let rejected = source
        .mesh(MeshOptions {
            maximum_triangles: 1,
            ..options
        })
        .unwrap_err();
    assert_eq!(rejected.kind, GeometryErrorKind::Limit);
    assert_eq!(source.mesh(options).unwrap(), first);
    let manufacturing = source.manufacturing_mesh_with_receipt(options).unwrap();
    assert_eq!(manufacturing.mesh, first);
    assert!(manufacturing.edge_pairs > 0 && manufacturing.vertex_links > 0);
    let moved = source
        .transformed(RigidTransform {
            rotation: [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
            translation_mm: [10., 0., 0.],
        })
        .unwrap()
        .mesh(options)
        .unwrap();
    assert_ne!(moved.vertices, first.vertices);
    assert_ne!(build(4.1).mesh(options).unwrap().vertices, first.vertices);
}

#[test]
fn cancellation_checkpoints_cover_hot_clones_and_new_results_without_hiding_errors() {
    DISPLAY.with(|cache| *cache.borrow_mut() = Cache::default());
    let source = owner();
    mesh(
        &source,
        options(),
        MeshUse::Display,
        None,
        || Ok(triangle()),
    )
    .unwrap();
    let never = || panic!("a cancelled warm display must not build");
    assert_eq!(
        mesh(&source, options(), MeshUse::Display, Some(&|| true), never)
            .unwrap_err()
            .kind,
        GeometryErrorKind::Cancelled
    );
    let calls = Cell::new(0);
    let after_clone = || {
        let n = calls.get() + 1;
        calls.set(n);
        n > 1
    };
    assert_eq!(
        mesh(
            &source,
            options(),
            MeshUse::Display,
            Some(&after_clone),
            never
        )
        .unwrap_err()
        .kind,
        GeometryErrorKind::Cancelled
    );
    let other = sibling(&source);
    calls.set(0);
    assert_eq!(
        mesh(
            &other,
            options(),
            MeshUse::Display,
            Some(&after_clone),
            || Ok(triangle())
        )
        .unwrap_err()
        .kind,
        GeometryErrorKind::Cancelled
    );
    assert!(
        DISPLAY
            .with(|cache| cache.borrow_mut().get(&other, Quality::from(options())))
            .is_none()
    );
    calls.set(0);
    let invalid = mesh(
        &other,
        options(),
        MeshUse::Display,
        Some(&after_clone),
        || Err(error("test", GeometryErrorKind::Input, "original failure")),
    )
    .unwrap_err();
    assert_eq!(invalid.kind, GeometryErrorKind::Input);
    assert!(
        DISPLAY
            .with(|cache| cache.borrow_mut().get(&source, Quality::from(options())))
            .is_some()
    );
}
