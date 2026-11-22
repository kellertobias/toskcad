use super::*;
use multicad_precision::PortableMath as _;
#[test]
fn metric_source_mesh_has_shared_closed_indices_and_actual_facet_bounds() {
    for fixture in 0..=5 {
        let (s, d) = super::super::super::super::tests::fixture(fixture);
        let source = Cutter::new(s, d).unwrap();
        let atlas = source.surface_patches().unwrap();
        let start = std::time::Instant::now();
        let result = atlas
            .boundary_mesh(Tessellation {
                deflection: 0.02,
                angle: 0.5,
                triangle_limit: 300_000,
            })
            .unwrap();
        eprintln!(
            "Phillips {fixture}: {} triangles {:?}, distance={}, normal={}",
            result.mesh().indices.len() / 3,
            start.elapsed(),
            result.source_distance_mm(),
            result.normal_chord()
        );
        assert_eq!(
            result.refinement(),
            TopologyCertificate {
                vertices: 66,
                edges: 152,
                faces: 88,
                euler_characteristic: 2
            }
        );
        assert_eq!(
            result.connectivity().vertex_links as usize,
            result.mesh().vertices.len() / 3
        );
        assert!(result.connectivity().signed_volume_mm3 > 0.);
        assert!(result.source_distance_mm() <= 0.02);
        assert!(result.normal_chord() <= 2. * (0.5_f64 / 2.).portable_sin());
        let source_volume = source.volume_bounds_mm3(0.001, 100_000).unwrap();
        assert!(
            (result.connectivity().signed_volume_mm3 - source_volume.volume_mm3[0]).abs()
                < source_volume.volume_mm3[1] * 0.03
        );
    }
}
#[test]
fn candidate_cancel_invalid_precision_and_work_limits_publish_no_partial_shell() {
    use std::cell::Cell;
    let (s, d) = super::super::super::super::tests::fixture(0);
    let source = Cutter::new(s, d).unwrap();
    let atlas = source.surface_patches().unwrap();
    let valid = Tessellation {
        deflection: 0.02,
        angle: 0.5,
        triangle_limit: 300_000,
    };
    for options in [
        Tessellation {
            deflection: f64::NAN,
            ..valid
        },
        Tessellation { angle: 0., ..valid },
        Tessellation {
            triangle_limit: 0,
            ..valid
        },
    ] {
        assert_eq!(atlas.boundary_mesh(options).unwrap_err(), Error::Invalid);
    }
    for limit in [1, 200, 500] {
        assert_eq!(
            atlas
                .boundary_mesh(Tessellation {
                    triangle_limit: limit,
                    ..valid
                })
                .unwrap_err(),
            Error::Limit
        );
    }
    assert_eq!(
        atlas
            .boundary_mesh(Tessellation {
                deflection: 1e-30,
                triangle_limit: 500,
                ..valid
            })
            .unwrap_err(),
        Error::Limit
    );
    for at in [0, 1, 10, 100, 300] {
        let calls = Cell::new(0);
        assert_eq!(
            atlas
                .boundary_mesh_with_cancel(valid, &|| {
                    let n = calls.get();
                    calls.set(n + 1);
                    n >= at
                })
                .unwrap_err(),
            Error::Cancelled
        );
        assert_eq!(calls.get(), at + 1);
    }
    let mut broken = atlas.clone();
    broken.patches.pop();
    assert!(matches!(
        Builder::new(&broken, &|| false),
        Err(Error::Unresolved)
    ));
    assert_eq!(
        budget([usize::MAX; 8], &atlas.patches, usize::MAX),
        Err(Error::Limit)
    );
    assert_eq!(mesh_error(MeshError::BudgetExceeded), Error::Limit);
    let invalid = DisplayMesh {
        vertices: vec![],
        indices: vec![],
        face_indices: vec![],
    };
    assert_eq!(
        manufacturing::certify(&invalid).map_err(mesh_error),
        Err(Error::Unresolved)
    );
}
#[test]
fn literal_facets_must_be_nondegenerate_oriented_and_match_the_true_normal() {
    let normal = [I::point(0.), I::point(0.), I::point(1.)];
    let valid = [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]];
    assert!(facet_chord(valid, normal, 0., Law::Floor(0)).unwrap() < 1e-12);
    assert!(facet_chord(valid, normal, 0.3, Law::Floor(0)).unwrap() >= 0.3);
    assert_eq!(
        facet_chord([[0.; 3]; 3], normal, 0., Law::Floor(0)),
        Err(Error::Unresolved)
    );
    assert_eq!(
        facet_chord([valid[0], valid[2], valid[1]], normal, 0., Law::Floor(0)),
        Err(Error::Unresolved)
    );
    assert_eq!(
        facet_chord(
            [[0.; 3], [f64::MAX, 0., 0.], [0., f64::MAX, 0.]],
            normal,
            0.,
            Law::Floor(0)
        ),
        Err(Error::Unresolved)
    );
    assert_eq!(
        facet_chord(valid, normal, f64::INFINITY, Law::Floor(0)),
        Err(Error::Unresolved)
    );
    let tilted = [[0., 0., 0.], [1., 0., 1.], [0., 1., 0.]];
    assert!(facet_chord(tilted, normal, 0., Law::Floor(0)).unwrap() > 0.7);
}

#[test]
fn real_contact_fold_requests_angular_refinement_without_winding_repair() {
    let (s, d) = super::super::super::super::tests::fixture(0);
    let source = Cutter::new(s, d).unwrap();
    let atlas = source.surface_patches().unwrap();
    let index = atlas
        .patches
        .iter()
        .position(|p| matches!(p.law, Law::Envelope(3, 0)))
        .unwrap();
    let law = atlas.patches[index].law;
    let mut builder = Builder::new(&atlas, &|| false).unwrap();
    let uv = [[0.875, 0.], [1., 0.], [1., 0.25], [0.875, 0.25]];
    let ids = uv.map(|p| builder.point(index, p).unwrap());
    let literal = builder.literal(ids);
    let center = atlas
        .state(index, [I::point(0.9375), I::point(0.125)])
        .unwrap()
        .1;
    let bad = facet_chord([literal[0], literal[3], literal[2]], center, 0., law);
    assert_eq!(bad, Err(Error::Unresolved));
    let mut refine = [false; 8];
    assert_eq!(resolve_facet(bad, law, &mut refine, [3, 4]), None);
    assert!(!refine[3] && refine[4]);
    let mut refine = [false; 8];
    assert_eq!(
        resolve_facet(
            Err(Error::Unresolved),
            Law::Wall(0, 0, false),
            &mut refine,
            [0, 5]
        ),
        None
    );
    assert!(refine[0] && refine[5]);
}

#[test]
fn default_h0_quality_refines_the_real_torus_normal_direction() {
    let (s, d) = super::super::super::super::tests::fixture(0);
    let source = Cutter::new(s, d).unwrap();
    let atlas = source.surface_patches().unwrap();
    let mesh = atlas
        .boundary_mesh(Tessellation {
            deflection: 0.1,
            angle: 0.3,
            triangle_limit: 100_000,
        })
        .unwrap();
    assert!(mesh.source_distance_mm() <= 0.1);
    assert!(mesh.normal_chord() <= 2. * (0.3_f64 / 2.).portable_sin());
    assert_eq!(
        atlas
            .boundary_mesh(Tessellation {
                deflection: 0.0001,
                angle: 0.5,
                triangle_limit: 5000
            })
            .unwrap_err(),
        Error::Limit
    );
}
#[test]
fn apex_preview_and_reordered_real_charts_preserve_the_same_closed_source() {
    let (s, d) = super::super::super::super::tests::fixture(0);
    let source = Cutter::new(s, d).unwrap();
    let mut atlas = source.surface_patches().unwrap();
    let fan = atlas
        .patches
        .iter()
        .position(|p| matches!(p.law, Law::Floor(0)))
        .unwrap();
    let mut preview = Builder::new(&atlas, &|| false).unwrap();
    assert!(quad_normal_score(&mut preview, fan, [[0., 1.], [0., 0.5]]).unwrap() > 0.);
    let first = atlas
        .patches
        .iter()
        .position(|p| matches!(p.law, Law::Envelope(3, 0)))
        .unwrap();
    atlas.patches.rotate_left(first);
    let result = atlas
        .boundary_mesh(Tessellation {
            deflection: 0.02,
            angle: 1.,
            triangle_limit: 100_000,
        })
        .unwrap();
    assert!(result.source_distance_mm() <= 0.02);
    assert!(result.connectivity().signed_volume_mm3 > 0.);
    assert_eq!(result.refinement().euler_characteristic, 2);
}
#[test]
fn real_folded_contact_grid_recovers_by_refinement_without_geometry_or_winding_repair() {
    let (s, d) = super::super::super::super::tests::fixture(0);
    let source = Cutter::new(s, d).unwrap();
    let atlas = source.surface_patches().unwrap();
    let options = Tessellation {
        deflection: 100.,
        angle: 1.,
        triangle_limit: 100_000,
    };
    for grid in [[0; 8], [3; 8]] {
        assert_eq!(
            atlas
                .refine_from_grid(options, &|| false, grid)
                .unwrap_err(),
            Error::Invalid
        );
    }
    // This is the independently checked actual contact fold from the coarse
    // triangle control above, now exercising the complete shared-grid recovery.
    let candidate = atlas
        .refine_from_grid(options, &|| false, [8, 8, 2, 8, 4, 4, 2, 4])
        .unwrap();
    assert!(candidate.connectivity().signed_volume_mm3 > 0.);
    assert!(candidate.source_distance_mm() <= options.deflection);
    assert!(candidate.normal_chord() <= 2. * (options.angle / 2.).portable_sin());
}

#[test]
fn thumbnail_contact_diagonal_earns_the_original_source_quality() {
    let (spec, _) = super::super::super::super::tests::fixture(0);
    let source = Cutter::new(
        spec,
        Domain {
            corner_start_depth_mm: -3.8,
            corner_end_depth_mm: 3.208443210824065,
            corner_reach_mm: 16.,
        },
    )
    .unwrap();
    let atlas = source.surface_patches().unwrap();
    let options = Tessellation {
        deflection: 0.024999999999707968,
        angle: 0.2591813939211579,
        triangle_limit: 2_000_000,
    };
    let mesh = atlas.boundary_mesh(options).unwrap();
    eprintln!(
        "M2 source mesh: {} triangles",
        mesh.mesh().indices.len() / 3
    );
    assert!(mesh.source_distance_mm() <= options.deflection);
    assert!(mesh.normal_chord() <= 2. * (options.angle / 2.).portable_sin());
    assert!(mesh.connectivity().signed_volume_mm3 > 0.);
    let index = atlas
        .patches
        .iter()
        .position(|p| matches!(p.law, Law::Envelope(3, 0)))
        .unwrap();
    let mut builder = Builder::new(&atlas, &|| false).unwrap();
    let uv = [[0.75, 1.], [0.75, 0.875]];
    let pairs = [
        [uv[0][0], uv[1][0]],
        [uv[0][1], uv[1][0]],
        [uv[0][1], uv[1][1]],
        [uv[0][0], uv[1][1]],
    ];
    let keys = pairs.map(|p| builder.point(index, p).unwrap());
    let literal = builder.literal(keys);
    let bound = atlas.certify_quad(index, uv, literal).unwrap();
    let center = atlas
        .state(index, uv.map(|[a, b]| I::point(a + (b - a) / 2.)))
        .unwrap()
        .1;
    let original = [[0, 2, 1], [0, 3, 2]].map(|tri| {
        facet_chord(
            tri.map(|i| literal[i]),
            center,
            bound.normal_radius_chord,
            atlas.patches[index].law,
        )
        .unwrap()
    });
    let allowed = 2. * (options.angle / 2.).portable_sin();
    assert!(original.into_iter().any(|score| score > allowed));
    let selected = choose_facets(
        keys,
        literal,
        center,
        bound.normal_radius_chord,
        atlas.patches[index].law,
        -1,
    );
    assert_eq!(selected[0].corners, [0, 3, 1]);
    assert!(
        selected
            .into_iter()
            .all(|f| f.chord.unwrap().unwrap() <= allowed)
    );
}

#[test]
fn diagonal_selection_does_not_erase_a_collapsed_literal_cell() {
    let literal = [[0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [1., 0., 0.]];
    let facets = choose_facets(
        [0, 1, 2, 1],
        literal,
        [I::point(0.), I::point(0.), I::point(1.)],
        0.,
        Law::Envelope(3, 0),
        1,
    );
    assert_eq!(facets[0].corners, [0, 1, 2]);
    assert!(
        facets
            .into_iter()
            .any(|facet| facet.chord == Some(Err(Error::Unresolved)))
    );
    let empty = choose_facets(
        [0; 4],
        [[0.; 3]; 4],
        [I::point(0.); 3],
        0.,
        Law::Envelope(3, 0),
        1,
    );
    assert_eq!(empty[0].corners, [0, 1, 2]);
    assert!(empty.into_iter().all(|facet| facet.chord.is_none()));
}

#[test]
fn actual_cell_and_preview_work_frontiers_refuse_atomically() {
    let (spec, domain) = super::super::super::super::tests::fixture(0);
    let source = Cutter::new(spec, domain).unwrap();
    let atlas = source.surface_patches().unwrap();
    let options = Tessellation {
        deflection: 100.,
        angle: 0.01,
        triangle_limit: 300_000,
    };
    let initial = 4 * atlas.patches.len();
    // The complete source graph precedes the first Envelope cell (two
    // diagonal candidates), then its eight prospective child candidates.
    for maximum_work in [initial, initial + 2] {
        let error = atlas
            .boundary_mesh_with_mouth_trace(
                options,
                MouthTraceOptions {
                    maximum_work,
                    ..Default::default()
                },
                &|| false,
            )
            .err();
        assert_eq!(error, Some(Error::Limit));
    }
}
