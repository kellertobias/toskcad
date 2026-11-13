use super::*;
use crate::{
    head_profile::{Head, Profile},
    helical_thread,
    phillips_drive::{Cavity, Cutter},
    profiled_screw::SourceSpec,
};
fn fixture(index: usize) -> PhillipsScrew {
    let (spec, domain) = crate::phillips_drive::test_fixture(index);
    let cavity = Cavity::new(Cutter::new(spec, domain).unwrap(), 1e-5, 100_000).unwrap();
    let (rs, pitch, r, h, crown, land, neck) = [
        (1., 0.4, 2., 1.6, 3.2, 1.25, 0.1),
        (1.5, 0.5, 2.8, 2.4, 3., 2.1, 0.1),
        (1.5, 0.5, 2.75, 1.65, 0., 0., 0.1),
        (2., 0.7, 4., 3.1, 6.5, 3., 0.2),
        (3., 1., 6., 4.6, 10., 4.2, 0.25),
        (4., 1.25, 8., 6., 13., 5.4, 0.4),
    ][index];
    let head = if index == 2 {
        Head::Countersunk {
            radius_mm: r,
            height_mm: h,
            included_angle_deg: 90.,
            edge_chamfer_mm: 0.05,
            neck_radius_mm: neck,
        }
    } else {
        Head::Crowned {
            radius_mm: r,
            height_mm: h,
            crown_radius_mm: crown,
            top_land_radius_mm: land,
            neck_radius_mm: neck,
        }
    };
    let thread = crate::helical_thread::Thread::new(helical_thread::Spec {
        major_radius_mm: rs - 0.03,
        pitch_mm: pitch,
        length_mm: 8.,
        depth_mm: 0.6 * pitch,
        crest_width_mm: 0.125 * pitch,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts: 1,
        left_hand: false,
        root_radius_mm: None,
    })
    .unwrap();
    PhillipsScrew::new(
        thread,
        Profile::new_for_thread_union(head, rs, 0.).unwrap(),
        SourceSpec {
            head: Some(head),
            shank_length_mm: 0.,
            minimum_wall_mm: 1e-5,
        },
        cavity,
        100_000,
        &|| false,
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
fn actual_metric_cavity_sews_to_original_thread_head_with_shared_closed_indices() {
    for index in 0..6 {
        let source = fixture(index);
        let thread = *source.thread();
        let now = std::time::Instant::now();
        let candidate = tessellate_phillips_screw(&source, options(), &|| false).unwrap();
        let mesh = candidate.mesh();
        eprintln!(
            "Phillips head {index}: {} triangles {:?}, volume {} nominal {}",
            mesh.indices.len() / 3,
            now.elapsed(),
            candidate.connectivity().signed_volume_mm3,
            source.nominal_volume_mm3()
        );
        assert_eq!(source.thread(), &thread);
        assert_eq!(candidate.faces(), source.faces());
        assert!(candidate.connectivity().signed_volume_mm3 > 0.);
        let literal = candidate.literal_volume_mm3();
        assert!(
            literal[0] <= candidate.connectivity().signed_volume_mm3
                && candidate.connectivity().signed_volume_mm3 <= literal[1]
        );
        let source_bounds = candidate.source_volume().volume_mm3();
        assert!(
            source_bounds[0] <= source.nominal_volume_mm3()
                && source.nominal_volume_mm3() <= source_bounds[1]
        );
        assert!(
            candidate.volume_difference_bound_mm3()
                >= (candidate.connectivity().signed_volume_mm3 - source.nominal_volume_mm3()).abs()
        );
        // Diagnostic agreement only. A precise source-to-mesh volume error
        // budget and global 3D embedding remain separate admission gates.
        let nominal = source.nominal_volume_mm3();
        assert!((candidate.connectivity().signed_volume_mm3 - nominal).abs() < nominal * 0.05);
        assert!(
            !candidate
                .faces()
                .contains(&ScrewFace::Cavity(FaceKey::Mouth))
        );
        for id in 0..candidate.faces().len() {
            assert!(
                mesh.face_indices.contains(&(id as u32)),
                "missing {index} {:?}",
                candidate.faces()[id]
            );
        }
        assert_eq!(
            manufacturing::certify(mesh).unwrap().edge_pairs,
            candidate.connectivity().edge_pairs
        );
    }
}
#[test]
fn directed_literal_volume_survives_large_exact_translation_without_geometry_changes() {
    // Binary-rational cube: its exact literal volume remains1/64 even at a
    // distant origin. This tests arithmetic, not any CAD feature admission.
    for origin in [[0.; 3], [1e12, -2e12, 3e12]] {
        let mut mesh = DisplayMesh {
            vertices: Vec::new(),
            indices: vec![
                0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 1, 2, 6, 1, 6, 5, 2, 3, 7, 2,
                7, 6, 3, 0, 4, 3, 4, 7,
            ],
            face_indices: vec![0; 12],
        };
        for p in [
            [0., 0., 0.],
            [0.25, 0., 0.],
            [0.25, 0.25, 0.],
            [0., 0.25, 0.],
            [0., 0., 0.25],
            [0.25, 0., 0.25],
            [0.25, 0.25, 0.25],
            [0., 0.25, 0.25],
        ] {
            mesh.vertices
                .extend(std::array::from_fn::<_, 3, _>(|i| origin[i] + p[i]));
        }
        manufacturing::certify(&mesh).unwrap();
        let interval = literal_volume(&mesh);
        assert!(interval.lo <= 1. / 64. && 1. / 64. <= interval.hi);
        assert!(interval.width() < 1e-12);
    }
}
#[test]
fn cancellation_precision_and_resource_failures_return_no_candidate() {
    let source = fixture(1);
    assert_eq!(
        tessellate_phillips_screw(&source, options(), &|| true).unwrap_err(),
        PhillipsMeshError::Cancelled
    );
    let mut invalid = options();
    invalid.deflection = 0.;
    assert_eq!(
        tessellate_phillips_screw(&source, invalid, &|| false).unwrap_err(),
        PhillipsMeshError::Geometry(MeshError::InvalidOptions)
    );
    invalid = options();
    invalid.deflection = 1e-20;
    assert_eq!(
        tessellate_phillips_screw(&source, invalid, &|| false).unwrap_err(),
        PhillipsMeshError::Geometry(MeshError::UnresolvedBoundary)
    );
    invalid = options();
    invalid.triangle_limit = 10;
    assert!(tessellate_phillips_screw(&source, invalid, &|| false).is_err());
    let checks = std::cell::Cell::new(0);
    assert_eq!(
        tessellate_phillips_screw(&source, options(), &|| {
            checks.set(checks.get() + 1);
            checks.get() > 1
        })
        .unwrap_err(),
        PhillipsMeshError::Cancelled
    );
    invalid = options();
    invalid.triangle_limit = 3000;
    assert_eq!(
        tessellate_phillips_screw(&source, invalid, &|| false).unwrap_err(),
        PhillipsMeshError::Geometry(MeshError::BudgetExceeded)
    );
    // Keep the boundary algorithm's defensive typed failure contract explicit;
    // an invalid/unresolved support is never converted into a partial mesh.
    assert_eq!(
        source_error(Error::Invalid),
        PhillipsMeshError::Geometry(MeshError::UnresolvedBoundary)
    );
    assert_eq!(
        source_error(Error::Unresolved),
        PhillipsMeshError::Geometry(MeshError::UnresolvedBoundary)
    );
}
#[test]
fn actual_cap_uses_and_literal_plane_order_cannot_be_repaired_or_welded() {
    let source = fixture(1);
    let topology = source.cavity().source().boundary_topology().unwrap();
    let cap = topology
        .faces()
        .iter()
        .position(|f| f.source_identity() == FaceKey::Mouth)
        .unwrap() as u32;
    let candidate = source
        .cavity()
        .source()
        .surface_patches()
        .unwrap()
        .boundary_mesh(Tessellation {
            angle: 0.5,
            ..options()
        })
        .unwrap();
    let raw = candidate.mesh();
    let proof = source.cavity().mouth_sweep();
    let ring = mouth_ring(raw, cap, proof).unwrap();
    let first = raw.face_indices.iter().position(|&f| f == cap).unwrap();
    let mut changed = raw.clone();
    changed.indices.swap(3 * first + 1, 3 * first + 2);
    assert!(mouth_ring(&changed, cap, proof).is_err());
    changed = raw.clone();
    changed
        .indices
        .extend_from_slice(&raw.indices[3 * first..3 * first + 3]);
    changed.face_indices.push(cap);
    assert!(mouth_ring(&changed, cap, proof).is_err());
    // Two literal cap regions pinched at the same vertex have ambiguous
    // outgoing boundary uses even though no edge is multiply covered.
    let n = (raw.vertices.len() / 3) as u32;
    changed = raw.clone();
    let root = ring[0];
    changed
        .vertices
        .extend([[4., 0., 0.], [4., 1., 0.]].into_iter().flatten());
    changed.indices.extend([root, n, n + 1]);
    changed.face_indices.push(cap);
    assert!(mouth_ring(&changed, cap, proof).is_err());
    changed = raw.clone();
    changed.face_indices.fill(cap + 100);
    assert!(mouth_ring(&changed, cap, proof).is_err());
    changed = raw.clone();
    for p in changed.vertices.as_chunks_mut::<3>().0 {
        let [x, y, z] = *p;
        *p = [
            x * 0.01f64.cos() - y * 0.01f64.sin(),
            x * 0.01f64.sin() + y * 0.01f64.cos(),
            z,
        ];
    }
    assert!(mouth_ring(&changed, cap, proof).is_err());
    changed = raw.clone();
    changed.vertices[3 * ring[0] as usize + 2] = 1e-9;
    assert!(mouth_ring(&changed, cap, proof).is_err());
    changed = raw.clone();
    changed.vertices[3 * ring[0] as usize] = 0.;
    changed.vertices[3 * ring[0] as usize + 1] = 0.;
    assert!(mouth_ring(&changed, cap, proof).is_err());
    // A second disconnected copy has valid internal opposite uses. Its cap
    // must still refuse the required single source mouth cycle.
    changed = raw.clone();
    let offset = (raw.vertices.len() / 3) as u32;
    changed.vertices.extend(
        raw.vertices
            .as_chunks::<3>()
            .0
            .iter()
            .flat_map(|&[x, y, z]| [x + 10., y, z]),
    );
    for (tri, &face) in raw.indices.as_chunks::<3>().0.iter().zip(&raw.face_indices) {
        if face == cap {
            changed.indices.extend(tri.map(|id| id + offset));
            changed.face_indices.push(face);
        }
    }
    assert!(mouth_ring(&changed, cap, proof).is_err());
}
#[test]
fn real_outer_shell_limits_late_cancel_and_unresolved_thin_annulus_refuse() {
    let source = fixture(1);
    for budget in [7500, 8000, 8500, 9000] {
        assert_eq!(
            tessellate_phillips_screw(
                &source,
                Tessellation {
                    triangle_limit: budget,
                    ..options()
                },
                &|| false
            )
            .unwrap_err(),
            PhillipsMeshError::Geometry(MeshError::BudgetExceeded)
        );
    }
    let checks = std::cell::Cell::new(0);
    let candidate = tessellate_phillips_screw(&source, options(), &|| {
        checks.set(checks.get() + 1);
        false
    })
    .unwrap();
    let total = checks.get();
    checks.set(0);
    assert_eq!(
        tessellate_phillips_screw(&source, options(), &|| {
            checks.set(checks.get() + 1);
            checks.get() + 10 > total
        })
        .unwrap_err(),
        PhillipsMeshError::Cancelled
    );
    assert_eq!(source.thread().spec().length_mm, 8.);
    assert!(candidate.connectivity().signed_volume_mm3 > 0.);
    let radius = source.cavity().mouth_sweep().maximum_radius_mm() + 1e-4;
    let head = Head::Cylindrical {
        radius_mm: radius,
        height_mm: 2.4,
    };
    let profile = Profile::new_for_thread_union(head, 1.5, 0.).unwrap();
    let thin = PhillipsScrew::new(
        *source.thread(),
        profile,
        SourceSpec {
            head: Some(head),
            ..source.spec()
        },
        source.cavity().clone(),
        100_000,
        &|| false,
    )
    .unwrap();
    assert!(thin.host_envelope().minimum_radial_clearance_mm() > thin.spec().minimum_wall_mm);
    // The true source fits; its coarse literal outer chords do not enclose
    // this very thin positive annulus. A planar constraint conflict refuses.
    assert_eq!(
        tessellate_phillips_screw(&thin, options(), &|| false).unwrap_err(),
        PhillipsMeshError::Geometry(MeshError::UnresolvedBoundary)
    );
}
#[test]
fn real_negative_zipper_cell_and_nonplanar_cap_are_refused() {
    let source = fixture(0);
    let work = options();
    let strips = phase_strips(source.thread(), 0.1, work).unwrap();
    let mut builder = Builder::new(source.thread(), 64, strips, work);
    let mut shell = HeadMesh {
        builder: &mut builder,
        angles: Vec::new(),
    };
    // Captured from the initial angular-zipper failure on the actual H0
    // source. The constrained planar region replaces that invalid cell.
    let tri = [
        [0.6204759418910268, 0.9511575102589167, -1.6],
        [0.7609517862609008, 0.9916916753640439, -1.6],
        [0.5609499586990121, 0.8932422448752108, -1.6],
    ]
    .map(|p| shell.push(p));
    assert_eq!(
        top_triangle(&mut shell, tri, 0),
        Err(MeshError::UnresolvedBoundary)
    );
    let tri = [[0., 0., 0.], [1., 0., 0.], [0., 1., 1e-9]].map(|p| shell.push(p));
    assert_eq!(
        top_triangle(&mut shell, tri, 0),
        Err(MeshError::UnresolvedBoundary)
    );
}
