//! Actual joined source-candidate controls; Native admission remains separate.
use super::super::super::common_mouth::tests::derived_thread_probe::with_head;
use super::*;
pub(super) fn owned(scale: f64, posed: bool) -> OwnedDerivedMaterial {
    owned_hand(scale, posed, false)
}
fn owned_hand(scale: f64, posed: bool, left_hand: bool) -> OwnedDerivedMaterial {
    owned_length(scale, posed, left_hand, 2. * scale)
}
fn owned_length(scale: f64, posed: bool, left_hand: bool, length: f64) -> OwnedDerivedMaterial {
    owned_neck(scale, posed, left_hand, length, 0.1 * scale)
}
pub(super) fn owned_neck(
    scale: f64,
    posed: bool,
    left_hand: bool,
    length: f64,
    neck: f64,
) -> OwnedDerivedMaterial {
    let input = super::super::super::super::branches::representative::scaled_text_with_neck(
        scale, posed, left_hand, length, neck,
    );
    let mut result = None;
    with_head(&input, &mut |head| {
        let thread = head
            .derived_thread_carrier(Default::default(), None)
            .unwrap()
            .certify_occupied(1e-7, 1e-5, None)
            .unwrap();
        let source = head
            .closed_repaired_material(&thread, Default::default(), None)
            .unwrap();
        result = Some(
            source
                .volume_bounds(Default::default(), None)
                .unwrap()
                .into_owned(Default::default(), None)
                .unwrap(),
        );
    });
    result.unwrap()
}

fn byte_hash(mesh: &crate::display::DisplayMesh) -> u64 {
    // Golden checksum of the EXACT independently embedded prototype bytes,
    // not a geometry/mass certificate or an import source.
    let mut bytes = Vec::new();
    bytes.extend(*b"CAD352M1");
    for n in [
        mesh.vertices.len(),
        mesh.indices.len(),
        mesh.face_indices.len(),
    ] {
        bytes.extend((n as u64).to_le_bytes());
    }
    for &v in &mesh.vertices {
        bytes.extend(v.to_bits().to_le_bytes());
    }
    for &i in &mesh.indices {
        bytes.extend(i.to_le_bytes());
    }
    for &i in &mesh.face_indices {
        bytes.extend(i.to_le_bytes());
    }
    bytes
        .into_iter()
        .fold(14_695_981_039_346_656_037u64, |h, b| {
            (h ^ u64::from(b)).wrapping_mul(1_099_511_628_211)
        })
}
#[test]
fn actual_production_joined_h0_candidate_uses_owned_geometry_and_bounded_ledger() {
    for (scale, posed) in [(1., false), (1.137, true)] {
        let source = owned(scale, posed);
        let start = std::time::Instant::now();
        let candidate = source
            .joined_mesh_candidate(JoinedMeshOptions::default(), Some(&|| false))
            .unwrap();
        assert!(std::ptr::eq(candidate.source(), &source));
        // The exact canonical Mac prototype's SHA/embedding proof is recorded
        // externally. Runtime equality/recovery checks compare this host's
        // actual source, avoiding a regenerated-platform-trig golden fixture.
        assert!(
            candidate
                .mesh()
                .face_indices
                .iter()
                .all(|&id| (id as usize) < source.face_count())
        );
        assert!(candidate.new_surface_error_mm() <= 0.1);
        assert_eq!(candidate.certified_tessellation().angle, 0.3);
        assert_eq!(candidate.certified_tessellation().deflection, 0.1);
        assert!(candidate.accounted_peak_bytes() <= 64 * 1024 * 1024);
        let bounds = candidate.delegated_workspace_bounds();
        assert!(bounds.retained_thread_bytes > bounds.connectivity_bytes);
        assert!(bounds.cutter_bytes > bounds.planar_bytes);
        assert!(candidate.combined_live_upper_bound_bytes() > candidate.accounted_peak_bytes());
        assert!(candidate.work_used() < 512_000_000);
        assert_eq!(
            candidate.work_used() + candidate.remaining_work(),
            512_000_000
        );
        let checked = crate::manufacturing::certify(candidate.mesh()).unwrap();
        assert_eq!(candidate.connectivity(), checked);
        assert_eq!(
            checked.vertex_links as usize,
            candidate.mesh().vertices.len() / 3
        );
        eprintln!(
            "PRODUCTION_JOINED scale{scale} posed{posed} facets{} work{} accountedpeak{} error{} time{:?}",
            candidate.mesh().indices.len() / 3,
            candidate.work_used(),
            candidate.accounted_peak_bytes(),
            candidate.new_surface_error_mm(),
            start.elapsed()
        );
    }
}
#[test]
fn actual_joined_mesh_resource_frontiers_and_cancellation_preserve_owned_source() {
    let source = owned(1., false);
    let defaults = JoinedMeshOptions::default();
    let original = source.joined_mesh_candidate(defaults, None).unwrap();
    let expected = byte_hash(original.mesh());
    // The exhaustive whole-workflow exact/one-less frontiers passed in
    // /tmp/multicad-joined-producer-first-coverage.log (605.87s unoptimized).
    // Routine checks exercise the earned candidate's real publication phase;
    // they conservatively include all original proof work as a prerequisite.
    super::mesh::tests::publication_frontiers(&original);
    for query in [
        JoinedMeshOptions {
            maximum_work: 512_000_001,
            ..defaults
        },
        JoinedMeshOptions {
            maximum_bytes: 64 * 1024 * 1024 + 1,
            ..defaults
        },
        JoinedMeshOptions {
            tessellation: crate::display::Tessellation {
                deflection: f64::NAN,
                ..defaults.tessellation
            },
            ..defaults
        },
        JoinedMeshOptions {
            tessellation: crate::display::Tessellation {
                triangle_limit: 0,
                ..defaults.tessellation
            },
            ..defaults
        },
    ] {
        assert!(matches!(
            source.joined_mesh_candidate(query, None),
            Err(PhillipsReadError::Invalid)
        ));
    }
    assert!(matches!(
        source.joined_mesh_candidate(
            JoinedMeshOptions {
                tessellation: crate::display::Tessellation {
                    angle: 0.01,
                    ..defaults.tessellation
                },
                ..defaults
            },
            None
        ),
        Err(PhillipsReadError::Unsupported)
    ));
    assert!(matches!(
        source.joined_mesh_candidate(defaults, Some(&|| true)),
        Err(PhillipsReadError::Cancelled)
    ));
    let floor = source.thread().boundary_error_mm();
    assert!(matches!(
        source.joined_mesh_candidate(
            JoinedMeshOptions {
                tessellation: crate::display::Tessellation {
                    deflection: 2. * floor,
                    ..defaults.tessellation
                },
                ..defaults
            },
            None
        ),
        Err(PhillipsReadError::Unsupported)
    ));
    let calls = std::cell::Cell::new(0usize);
    assert!(matches!(
        source.joined_mesh_candidate(
            defaults,
            Some(&|| {
                let n = calls.get() + 1;
                calls.set(n);
                n >= 3
            })
        ),
        Err(PhillipsReadError::Cancelled)
    ));
    // Coarse, actual source budgets fail without changing the owner, and an
    // ordinary retry recreates exactly the independently checked candidate.
    assert!(matches!(
        source.joined_mesh_candidate(
            JoinedMeshOptions {
                tessellation: crate::display::Tessellation {
                    triangle_limit: 16,
                    ..defaults.tessellation
                },
                ..defaults
            },
            None
        ),
        Err(PhillipsReadError::Limit)
    ));
    assert_eq!(
        byte_hash(source.joined_mesh_candidate(defaults, None).unwrap().mesh()),
        expected
    );
}

#[test]
fn actual_left_handed_joined_source_keeps_original_cap_orientation() {
    let source = owned_hand(1.137, true, true);
    let candidate = source
        .joined_mesh_candidate(Default::default(), None)
        .unwrap();
    assert_eq!(
        candidate.connectivity().vertex_links as usize,
        candidate.mesh().vertices.len() / 3
    );
    assert!(candidate.connectivity().signed_volume_mm3 > 0.);
    assert!(candidate.new_surface_error_mm() <= 0.1);
    assert!(source.original_input().contains("ISO-10303-21"));
    eprintln!(
        "ACTUAL_LH_JOIN facets{} work{} accounted{} combinedupper{}",
        candidate.mesh().face_indices.len(),
        candidate.work_used(),
        candidate.accounted_peak_bytes(),
        candidate.combined_live_upper_bound_bytes()
    );
}

#[test]
fn actual_short_thread_head_uses_independent_angular_grid_when_needed() {
    for length in [0.4, 0.5] {
        let source = owned_length(1., false, false, length);
        let start = std::time::Instant::now();
        let candidate = source
            .joined_mesh_candidate(Default::default(), None)
            .unwrap();
        assert!(candidate.new_surface_error_mm() <= 0.1);
        assert!(candidate.connectivity().signed_volume_mm3 > 0.);
        eprintln!(
            "SHORT_JOIN length{length} facets{} work{} time{:?}",
            candidate.mesh().face_indices.len(),
            candidate.work_used(),
            start.elapsed()
        );
    }
}
#[test]
fn actual_joined_source_at_native_transfer_quality_keeps_full_source_certificates() {
    let source = owned(1., false);
    for (deflection, angle) in [(0.1, 0.297), (0.05, 0.297)] {
        let start = std::time::Instant::now();
        let result = source.joined_mesh_candidate(
            JoinedMeshOptions {
                tessellation: crate::display::Tessellation {
                    deflection,
                    angle,
                    triangle_limit: 1_000_000,
                },
                ..Default::default()
            },
            None,
        );
        let result = result.unwrap();
        eprintln!(
            "TRANSFER_QUALITY D{deflection} A{angle} PASS{}facets work{} time{:?}",
            result.mesh().face_indices.len(),
            result.work_used(),
            start.elapsed()
        );
        assert!(result.new_surface_error_mm() <= deflection);
    }
}
