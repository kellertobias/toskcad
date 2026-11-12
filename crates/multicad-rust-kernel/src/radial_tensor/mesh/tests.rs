use super::*;
use std::cell::Cell;
pub(super) fn cam() -> RadialBoundary {
    BoundaryOwnership::new(crate::radial_tensor::tests::cam::wedge(&[1, 2]), 1e-7, None)
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
        .unwrap()
}
pub(super) fn options() -> Tessellation {
    Tessellation {
        deflection: 0.05,
        angle: 0.3,
        triangle_limit: 50000,
    }
}
#[test]
fn restricted_source_chords_keep_actual_facets_ranges_and_owners_independent() {
    let thread = crate::tapered_thread::Thread::new(
        crate::helical_thread::Spec {
            major_radius_mm: 4.,
            pitch_mm: 1.25,
            length_mm: 1.25,
            depth_mm: 0.5,
            crest_width_mm: 0.25,
            lower_flank_deg: 12.,
            upper_flank_deg: 18.,
            starts: 1,
            left_hand: false,
            root_radius_mm: None,
        },
        0.03125,
        1e-7,
        1e-6,
    )
    .unwrap();
    let body = thread
        .retained_boundary_with_interrupt([[0.; 3], [1., 0., 0.], [0., 1., 0.]], 1e-7, 1e-6, None)
        .unwrap()
        .boundary;
    let mut work = Work {
        remaining: 100_000_000,
        interrupt: None,
    };
    let prepared = Prepared::new(&body, options(), DEFAULT_LIMITS, &mut work).unwrap();
    let side = prepared
        .sides
        .iter()
        .position(|s| matches!(s, Side::Tensor(_)))
        .unwrap();
    let domains = if let Side::Tensor(t) = &prepared.sides[side] {
        t.patches([[0., 1.]; 2], &mut work)
            .unwrap()
            .iter()
            .map(|p| p.domain)
            .collect::<Vec<_>>()
    } else {
        Vec::new()
    };
    assert!(!domains.is_empty());
    for fraction in [[0.11, 0.13], [0.61, 0.63]] {
        let domain = domains[0];
        let range = domain.map(|d| fraction.map(|f| d[0] + (d[1] - d[0]) * f));
        let uv = [
            [range[0][0], range[1][0]],
            [range[0][1], range[1][0]],
            [range[0][1], range[1][1]],
        ];
        let mut vertices = [[0.; 3]; 3];
        for i in 0..3 {
            vertices[i] = prepared
                .point(side, uv[i][0], uv[i][1], None, &mut work)
                .unwrap();
        }
        let patches = if let Side::Tensor(t) = &prepared.sides[side] {
            t.patches(range, &mut work).unwrap()
        } else {
            Vec::new()
        };
        assert_eq!(patches.len(), 1);
        let patch = patches[0];
        let mut source = patch.controls.clone();
        let cell = certificate::SourceCell::new(&source, domain, range, &mut work).unwrap();
        let strip = certificate::SourceStrip::new(&source, domain, range[0], &mut work).unwrap();
        let cached = strip.cell(range[1], &mut work).unwrap();
        let affine = certificate::Affine::new(uv, vertices).unwrap();
        let bound = cell.chord_error(&affine, &mut work).unwrap();
        assert_eq!(
            cached.chord_error(&affine, &mut work).unwrap().to_bits(),
            bound.to_bits()
        );
        assert!(bound < 0.1);
        // Independently evaluate the actual residual polynomial. Subtracting
        // before evaluation avoids treating correlated world coordinates as
        // unrelated interval boxes near a zero residual.
        let residual = source
            .iter()
            .enumerate()
            .map(|(i, row)| {
                row.iter()
                    .enumerate()
                    .map(|(j, p)| {
                        let local = [
                            I::point(i as f64) / I::point((source.len() - 1) as f64),
                            I::point(j as f64) / I::point((row.len() - 1) as f64),
                        ];
                        let global = std::array::from_fn::<_, 2, _>(|k| {
                            I::point(domain[k][0])
                                + (I::point(domain[k][1]) - I::point(domain[k][0])) * local[k]
                        });
                        let at = std::array::from_fn::<_, 2, _>(|k| {
                            (global[k] - I::point(range[k][0]))
                                / (I::point(range[k][1]) - I::point(range[k][0]))
                        });
                        std::array::from_fn(|k| {
                            p[k] - (I::point(vertices[0][k])
                                + (I::point(vertices[1][k]) - I::point(vertices[0][k])) * at[0]
                                + (I::point(vertices[2][k]) - I::point(vertices[1][k])) * at[1])
                        })
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        for du in [0., 0.37, 1.] {
            for dv in [0., 0.63, 1.] {
                let at = [du, dv];
                let local = std::array::from_fn::<_, 2, _>(|k| {
                    let global = I::point(range[k][0])
                        + (I::point(range[k][1]) - I::point(range[k][0])) * I::point(at[k]);
                    (global - I::point(domain[k][0]))
                        / (I::point(domain[k][1]) - I::point(domain[k][0]))
                });
                let error = evaluate_source(&residual, local)
                    .into_iter()
                    .map(I::magnitude)
                    .sum::<f64>();
                assert!(error <= bound, "residual={error}, bound={bound}");
            }
        }
        for sign in [-1., 1.] {
            let shifted = vertices.map(|mut p| {
                p[0] += sign;
                p
            });
            let other = certificate::Affine::new(uv, shifted).unwrap();
            assert!(cell.chord_error(&other, &mut work).unwrap() > 0.9);
        }
        for p in source.iter_mut().flatten() {
            p[0] = p[0] + I::point(1.);
        }
        let other = certificate::SourceCell::new(&source, domain, range, &mut work).unwrap();
        assert!(other.chord_error(&affine, &mut work).unwrap() > 0.9);
        assert_eq!(cell.chord_error(&affine, &mut work).unwrap(), bound);
    }
}
fn evaluate_source(poles: &[Vec<[I; 3]>], at: [I; 2]) -> [I; 3] {
    let reduce = |mut row: Vec<[I; 3]>, t: I| {
        while row.len() > 1 {
            row = row
                .windows(2)
                .map(|p| std::array::from_fn(|k| (I::point(1.) - t) * p[0][k] + t * p[1][k]))
                .collect();
        }
        row[0]
    };
    reduce(
        poles.iter().map(|row| reduce(row.clone(), at[1])).collect(),
        at[0],
    )
}
#[test]
fn real_source_mesh_is_cache_independent_and_failed_requests_are_atomic() {
    let body = cam();
    let original = body.input().clone();
    let edge = body
        .sample_edge_with_interrupt(0, options().deflection, 1000, None)
        .unwrap();
    let endpoints = original.edges[0].vertices.map(|i| original.vertices[i]);
    assert_eq!(edge.first(), Some(&endpoints[0]));
    assert_eq!(edge.last(), Some(&endpoints[1]));
    let regular = body.tessellate_with_interrupt(options(), None).unwrap();
    let uncached = build_with_limits(
        &body,
        options(),
        None,
        Limits {
            points: 0,
            certificates: 0,
            strips: 0,
            ..DEFAULT_LIMITS
        },
    )
    .unwrap();
    assert_eq!(regular, uncached);
    let tiny = build_with_limits(
        &body,
        options(),
        None,
        Limits {
            points: 1,
            certificates: 1,
            strips: 1,
            ..DEFAULT_LIMITS
        },
    )
    .unwrap();
    assert_eq!(regular, tiny);
    assert!(
        crate::manufacturing::certify(&regular)
            .unwrap()
            .signed_volume_mm3
            > 0.
    );
    for limits in [
        Limits {
            work: 0,
            ..DEFAULT_LIMITS
        },
        Limits {
            passes: 0,
            ..DEFAULT_LIMITS
        },
        Limits {
            work: 100,
            ..DEFAULT_LIMITS
        },
    ] {
        assert_eq!(
            build_with_limits(&body, options(), None, limits).err(),
            Some(Error::Limit)
        );
    }
    let count = Cell::new(0);
    let probe = || {
        count.set(count.get() + 1);
        false
    };
    build(&body, options(), Some(&probe)).unwrap();
    let total = count.get();
    assert!(total > 100);
    for at in [1, total / 2, total] {
        count.set(0);
        let cancel = || {
            count.set(count.get() + 1);
            count.get() >= at
        };
        assert_eq!(
            build(&body, options(), Some(&cancel)).err(),
            Some(Error::Cancelled)
        );
        assert_eq!(count.get(), at);
    }
    assert_eq!(body.input(), &original);
    assert_eq!(build(&body, options(), None).unwrap(), regular);
}
#[test]
fn actual_proper_and_distant_source_frames_keep_certified_closed_meshes() {
    let body = cam();
    for origin in [[100., -30., 50.], [1e6, -2e6, 3e6]] {
        let frame = [origin, [0.8, 0.6, 0.], [-0.36, 0.48, 0.8]];
        let placed = body
            .transformed_with_interrupt(frame, 1e-6, 1e-5, None)
            .unwrap();
        let mesh = build(&placed, options(), None).unwrap();
        assert!(
            crate::manufacturing::certify(&mesh)
                .unwrap()
                .signed_volume_mm3
                > 0.
        );
    }
}
#[test]
fn invalid_precision_and_allocation_limits_refuse_before_publication() {
    let body = cam();
    for invalid in [0., -1., f64::NAN, f64::INFINITY] {
        assert_eq!(
            build(
                &body,
                Tessellation {
                    deflection: invalid,
                    ..options()
                },
                None
            )
            .err(),
            Some(Error::Invalid)
        );
        assert_eq!(
            build(
                &body,
                Tessellation {
                    angle: invalid,
                    ..options()
                },
                None
            )
            .err(),
            Some(Error::Invalid)
        );
    }
    assert_eq!(
        build(
            &body,
            Tessellation {
                angle: 4.,
                ..options()
            },
            None
        )
        .err(),
        Some(Error::Invalid)
    );
    assert_eq!(
        build(
            &body,
            Tessellation {
                triangle_limit: 0,
                ..options()
            },
            None
        )
        .err(),
        Some(Error::Invalid)
    );
    for limit in [1, usize::MAX] {
        assert_eq!(
            build(
                &body,
                Tessellation {
                    triangle_limit: limit,
                    ..options()
                },
                None
            )
            .err(),
            Some(Error::Limit)
        );
    }
    assert_eq!(
        build(
            &body,
            Tessellation {
                deflection: body.error_mm,
                ..options()
            },
            None
        )
        .err(),
        Some(Error::Unresolved)
    );
    let mut work = Work {
        remaining: 1000,
        interrupt: None,
    };
    let a = certificate::Affine::new(
        [[0., 0.], [1., 0.], [0., 1.]],
        [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]],
    )
    .unwrap();
    let degenerate = certificate::Affine::new([[0., 0.]; 3], [[0.; 3]; 3]);
    assert_eq!(degenerate.err(), Some(Error::Unresolved));
    let linear = vec![vec![[I::point(0.); 3]; 2]; 2];
    assert!(
        certificate::chord_error(&linear, [[0., 1.]; 2], [[0., 1.]; 2], &a, &mut work).unwrap()
            >= 1.
    );
    assert_eq!(uncertain(Err(Error::Cancelled)), Err(Error::Cancelled));
    assert!(!uncertain(Err(Error::Unresolved)).unwrap());
}
#[test]
fn actual_piecewise_curve_knots_and_facet_normal_failures_remain_explicit() {
    let poles = (0..=18)
        .map(|i| {
            let y = if i <= 9 {
                i as f64 / 9.
            } else {
                1. + 2. * (i - 9) as f64 / 9.
            };
            [I::point(4.), I::point(y), I::point(i as f64 / 9.)]
        })
        .collect();
    let rail = spline::Spline::new(poles);
    let mut work = Work {
        remaining: 100000,
        interrupt: None,
    };
    let jump = rail.angular_jumps([0., 1.], &mut work).unwrap();
    assert!(jump.hi >= 8. / 17. && jump.hi < 8. / 17. + 1e-12);
    assert_eq!(rail.angular_jumps([0., 0.5], &mut work).unwrap().hi, 0.);
    let normal = vec![vec![[I::point(0.), I::point(0.), I::point(1.)]]];
    for (facet, angle, expected) in [
        ([0., 0., 1.], 0.1, true),
        ([0., 0., -1.], 0.1, false),
        ([1., 0., 0.1], 0.1, false),
        ([1., 0., 0.1], 2., true),
    ] {
        assert_eq!(
            certificate::normal_matches(
                normal.clone(),
                [[0., 1.]; 2],
                [[0., 1.]; 2],
                facet.map(I::point),
                angle,
                &mut work
            )
            .unwrap(),
            expected
        );
    }
    let body = cam();
    assert_eq!(
        build(
            &body,
            Tessellation {
                angle: 1e-14,
                ..options()
            },
            None
        )
        .err(),
        Some(Error::Unresolved)
    );
    assert!(crate::manufacturing::certify(&build(&body, options(), None).unwrap()).is_ok());
}
fn candidate(ring: usize, turns: usize) -> DisplayMesh {
    let mut mesh = DisplayMesh {
        vertices: Vec::new(),
        indices: Vec::new(),
        face_indices: Vec::new(),
    };
    for end in 0..2 {
        for i in 0..ring {
            let a = std::f64::consts::TAU * turns as f64 * i as f64 / ring as f64;
            mesh.vertices
                .extend([5. * a.cos(), 5. * a.sin(), 4. * end as f64]);
        }
    }
    mesh.vertices.extend([0., 0., -0.01, 0., 0., 4.01]);
    for i in 0..ring {
        for t in [
            [i, (i + 1) % ring, ring + (i + 1) % ring],
            [i, ring + (i + 1) % ring, ring + i],
        ] {
            mesh.triangle(t, 2);
        }
    }
    for i in 0..ring {
        mesh.triangle([2 * ring, (i + 1) % ring, i], 0);
        mesh.triangle([2 * ring + 1, ring + i, ring + (i + 1) % ring], 1);
    }
    mesh
}
#[test]
fn literal_embedding_rejects_self_wrapping_folded_and_unseparated_candidates() {
    let mut work = Work {
        remaining: 100000,
        interrupt: None,
    };
    let good = candidate(16, 1);
    assert!(matches!(
        embedding::certify(&good, 16, 2, &mut work).unwrap(),
        embedding::Verdict::Certified
    ));
    let mut changed = good.clone();
    changed.vertices[3 * 32 + 2] = 0.;
    assert!(matches!(
        embedding::certify(&changed, 16, 2, &mut work).unwrap(),
        embedding::Verdict::HeightUncertain
    ));
    let mut changed = good.clone();
    changed.vertices[3 * 33 + 2] = 4.;
    assert!(matches!(
        embedding::certify(&changed, 16, 2, &mut work).unwrap(),
        embedding::Verdict::HeightUncertain
    ));
    let mut changed = good.clone();
    changed.vertices[2] = 3.;
    changed.vertices[3 * 16 + 2] = 2.;
    assert!(matches!(
        embedding::certify(&changed, 16, 2, &mut work).unwrap(),
        embedding::Verdict::HeightUncertain
    ));
    let mut changed = good.clone();
    changed.vertices[3..6].copy_from_slice(&good.vertices[..3]);
    assert!(matches!(
        embedding::certify(&changed, 16, 2, &mut work).unwrap(),
        embedding::Verdict::AngularUncertain(_)
    ));
    let mut changed = good;
    changed.indices.swap(0, 1);
    assert!(matches!(
        embedding::certify(&changed, 16, 2, &mut work).unwrap(),
        embedding::Verdict::RadialUncertain(_)
    ));
    assert!(matches!(
        embedding::certify(&candidate(16, 2), 16, 2, &mut work),
        Err(Error::Unresolved)
    ));
}
#[test]
fn source_cap_seams_endpoint_ownership_and_precision_receipts_are_explicit() {
    let body = cam();
    let mut work = Work {
        remaining: 10000000,
        interrupt: None,
    };
    let prepared = Prepared::new(&body, options(), DEFAULT_LIMITS, &mut work).unwrap();
    for end in 0..2 {
        for side in 0..prepared.sides.len() {
            let p = prepared.caps[end]
                .point(&body, end, side, 1., &mut work)
                .unwrap();
            let edge = *body.ordered.caps.rims[end][side].last().unwrap();
            let vertex = body.input().vertices[caps::ends(body.input(), edge).1];
            for k in 0..3 {
                assert!(p[k].lo <= vertex[k] && p[k].hi >= vertex[k]);
            }
            if body.ordered.caps.rims[end][side].len() == 2 {
                assert_eq!(
                    prepared.caps[end]
                        .curve_error(side, [0.25, 0.75], [[0.; 3]; 2], &mut work)
                        .err(),
                    Some(Error::Invalid)
                );
            }
        }
    }
    let closure = prepared
        .sides
        .iter()
        .filter_map(|s| {
            if let Side::Cone(c) = s {
                Some(c.closure_error)
            } else {
                None
            }
        })
        .fold(0_f64, f64::max);
    let floor = (I::point(4.) * I::point(body.error_mm)).hi;
    assert!(closure > floor);
    let deflection = floor + (closure - floor) * 0.5;
    assert_eq!(
        build(
            &body,
            Tessellation {
                deflection,
                ..options()
            },
            None
        )
        .err(),
        Some(Error::Unresolved)
    );
    let grid = Grid {
        u: vec![vec![0., f64::from_bits(1)]; prepared.sides.len()],
        v: vec![0., 1.],
    };
    let mut changes = Refinement::new(prepared.sides.len());
    assert_eq!(
        changes
            .mark(&prepared, &grid, [0, 0, 0], 0, &mut work)
            .err(),
        Some(Error::Unresolved)
    );
    let normal = vec![vec![[I::point(0.), I::point(0.), I::point(1.)]]];
    assert_eq!(
        certificate::normal_matches(
            normal,
            [[0., 1.]; 2],
            [[2., 3.]; 2],
            [I::point(0.); 3],
            0.1,
            &mut work
        )
        .err(),
        Some(Error::Invalid)
    );
    let mut added = BTreeSet::new();
    added.insert(0.5_f64.to_bits());
    assert_eq!(
        merge(&[0., 0.5, 1.], added, &mut work).unwrap(),
        vec![0., 0.5, 1.]
    );
}
fn certify(input: BoundaryInput) -> RadialBoundary {
    BoundaryOwnership::new(input, 1e-7, None)
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
        .unwrap()
}
#[test]
fn actual_cap_sewing_collar_requires_height_retry_and_refuses_unproved_closure() {
    let mut input = crate::radial_tensor::tests::cam::wedge(&[1]);
    let point = |p: Point| {
        let mut q = p.map(|v| v * 1e-3);
        if q[2] == 0. {
            q[2] -= 3e-9;
        }
        q
    };
    for vertex in &mut input.vertices {
        *vertex = point(*vertex);
    }
    for curve in &mut input.curves {
        match curve {
            Curve::Spline { poles, .. } => {
                for p in poles {
                    *p = point(*p);
                }
            }
            Curve::Circle { axis, radius_mm } => {
                axis.origin = axis.origin.map(|v| v * 1e-3);
                *radius_mm *= 1e-3;
            }
        }
    }
    for surface in &mut input.surfaces {
        match surface {
            Surface::Tensor(t) => {
                *t = RetainedTensor::new(t.spans(), t.poles().iter().copied().map(point).collect())
                    .unwrap()
            }
            Surface::Plane(axis) => axis.origin = axis.origin.map(|v| v * 1e-3),
            Surface::Cone {
                axis, radius_mm, ..
            } => {
                axis.origin = axis.origin.map(|v| v * 1e-3);
                *radius_mm *= 1e-3;
            }
        }
    }
    // The actual retained cap curves are within the certified sewing
    // collar, but below their exact support plane. They are not snapped.
    let body = certify(input);
    let mesh = build(
        &body,
        Tessellation {
            deflection: 4e-8,
            ..options()
        },
        None,
    )
    .unwrap();
    let count = mesh.vertices.len() / 3;
    // The first two requested center shifts leave the center above the
    // actual bottom ring. A later certified shift must be retained.
    assert!(mesh.vertices[3 * (count - 2) + 2] < -3e-9);
    assert!(crate::manufacturing::certify(&mesh).is_ok());
}
#[test]
fn oppositely_stored_cap_support_keeps_outward_material_orientation() {
    let mut input = crate::radial_tensor::tests::cam::wedge(&[1]);
    let surface = input.faces[0].surface;
    let Surface::Plane(axis) = &mut input.surfaces[surface] else {
        panic!()
    };
    axis.normal = axis.normal.map(|v| -v);
    input.faces[0].forward = false;
    let body = certify(input);
    assert!(
        crate::manufacturing::certify(&build(&body, options(), None).unwrap())
            .unwrap()
            .signed_volume_mm3
            > 0.
    );
}
#[test]
fn direct_signed_taper_and_multistart_sources_get_actual_manifold_meshes() {
    for (starts, taper, left_hand) in [
        (1, 0.03125, false),
        (8, 0.03125, false),
        (2, -0.03125, true),
    ] {
        let spec = crate::helical_thread::Spec {
            major_radius_mm: 4.,
            pitch_mm: 1.25,
            length_mm: 1.25,
            depth_mm: 0.5,
            crest_width_mm: 0.25,
            lower_flank_deg: 12.,
            upper_flank_deg: 18.,
            starts,
            left_hand,
            root_radius_mm: None,
        };
        let source = crate::tapered_thread::Thread::new(spec, taper, 1e-7, 1e-6).unwrap();
        let authored = source
            .retained_boundary_with_interrupt(
                [[0.; 3], [1., 0., 0.], [0., 1., 0.]],
                1e-7,
                1e-6,
                None,
            )
            .unwrap();
        let mesh = build(
            &authored.boundary,
            Tessellation {
                deflection: 0.1,
                angle: 0.3,
                triangle_limit: 200000,
            },
            None,
        )
        .unwrap();
        assert!(
            crate::manufacturing::certify(&mesh)
                .unwrap()
                .signed_volume_mm3
                > 0.
        );
        assert_eq!(
            mesh.face_indices
                .iter()
                .copied()
                .collect::<BTreeSet<_>>()
                .len(),
            authored.boundary.input().faces.len()
        );
    }
}
#[test]
fn failed_actual_candidate_facet_checks_and_mid_proof_cancellation_cannot_publish() {
    let body = cam();
    let mut work = Work {
        remaining: 10000000,
        interrupt: None,
    };
    let mut source = Prepared::new(&body, options(), DEFAULT_LIMITS, &mut work).unwrap();
    let side = source
        .sides
        .iter()
        .position(|s| matches!(s, Side::Tensor(_)))
        .unwrap();
    let range = [[0., 1.]; 2];
    let mut corners = [[0.; 3]; 4];
    for (index, [u, v]) in [[0., 0.], [1., 0.], [1., 1.], [0., 1.]]
        .into_iter()
        .enumerate()
    {
        corners[index] = source.point(side, u, v, None, &mut work).unwrap();
    }
    let owner = match &source.sides[side] {
        Side::Tensor(t) => t,
        Side::Cone(_) => unreachable!(),
    };
    let mut strip_cache = strip::Cache::new(owner, DEFAULT_LIMITS.strips);
    assert!(
        !source
            .cell_matches(
                side,
                range,
                corners,
                [[2, 1, 0], [3, 2, 0]],
                None,
                Some(&mut strip_cache),
                &mut work
            )
            .unwrap()
    );
    // The first actual facet may qualify against a source enclosure reused by
    // the second. That enclosure cannot certify the second facet's reversed
    // literal orientation.
    assert!(
        !source
            .cell_matches(
                side,
                range,
                corners,
                [[0, 1, 2], [3, 2, 0]],
                None,
                Some(&mut strip_cache),
                &mut work
            )
            .unwrap()
    );
    let count = Cell::new(0);
    let observe = || {
        count.set(count.get() + 1);
        false
    };
    let mut counted = Work {
        remaining: 10000000,
        interrupt: Some(&observe),
    };
    assert!(
        source
            .cell_matches(
                side,
                range,
                corners,
                [[0, 1, 2], [0, 2, 3]],
                None,
                Some(&mut strip_cache),
                &mut counted
            )
            .unwrap()
    );
    let total = count.get();
    assert!(total > 2);
    for at in 1..=total {
        count.set(0);
        let cancel = || {
            count.set(count.get() + 1);
            count.get() >= at
        };
        let mut cut = Work {
            remaining: 10000000,
            interrupt: Some(&cancel),
        };
        assert_eq!(
            source
                .cell_matches(
                    side,
                    range,
                    corners,
                    [[0, 1, 2], [0, 2, 3]],
                    None,
                    Some(&mut strip_cache),
                    &mut cut
                )
                .err(),
            Some(Error::Cancelled)
        );
    }
    let grid = Grid {
        u: source
            .sides
            .iter()
            .map(|s| {
                if matches!(s, Side::Cone(_)) {
                    vec![0., 0.5, 1.]
                } else {
                    vec![0., 1.]
                }
            })
            .collect(),
        v: vec![0., 1.],
    };
    let ring = Ring::new(&grid);
    let candidate = source.emit(&grid, &ring, 0.001, 50000, &mut work).unwrap();
    let mut changes = Refinement::new(grid.u.len());
    let mut cut = Work {
        remaining: 1000000,
        interrupt: Some(&|| true),
    };
    assert_eq!(
        source
            .check(&candidate, &grid, &ring, 0, &mut changes, &mut cut)
            .err(),
        Some(Error::Cancelled)
    );
    assert!(changes.empty());
    assert!(source.hint(side, range, 1, &mut work).unwrap() <= 1);
    let cone = source
        .sides
        .iter()
        .position(|s| matches!(s, Side::Cone(_)))
        .unwrap();
    assert_eq!(
        source.hint(cone, range, 0, &mut cut).err(),
        Some(Error::Cancelled)
    );
}
#[test]
fn recovery_planner_uses_actual_rejected_candidates_without_publishing_them() {
    let body = cam();
    let mut work = Work {
        remaining: 10000000,
        interrupt: None,
    };
    let source = Prepared::new(&body, options(), DEFAULT_LIMITS, &mut work).unwrap();
    let grid = Grid {
        u: vec![vec![0., 0.5, 1.]; source.sides.len()],
        v: vec![0., 1.],
    };
    let ring = Ring::new(&grid);
    let n = ring.entries.len();
    assert_eq!(n, 16);
    let good = candidate(n, 1);
    let mut changes = Refinement::new(grid.u.len());
    assert!(
        changes
            .embedding(&source, &grid, &ring, &good, &mut work)
            .unwrap()
    );
    assert!(changes.empty());
    let mut angular = good.clone();
    angular.vertices[3..6].copy_from_slice(&good.vertices[..3]);
    assert!(
        !changes
            .embedding(&source, &grid, &ring, &angular, &mut work)
            .unwrap()
    );
    assert!(!changes.empty());
    let mut radial = good;
    radial.indices.swap(0, 1);
    let mut changes = Refinement::new(grid.u.len());
    assert!(
        !changes
            .embedding(&source, &grid, &ring, &radial, &mut work)
            .unwrap()
    );
    assert!(!changes.v.is_empty());
    assert!(changes.u.iter().any(|u| !u.is_empty()));
    for pass in [0, 1] {
        assert_eq!(refinement_axis([0., 0.], pass), pass);
    }
    assert_eq!(refinement_axis([2., 1.], 0), 0);
    assert_eq!(refinement_axis([1., 2.], 0), 1);
    for slope in [I::point(-10000.), I::point(10000.)] {
        let angle = cone::angular_reserve(I::point(1e-6), slope, 1e-9).unwrap();
        assert!(angle >= std::f64::consts::PI / (1_f64 + 1e8).sqrt());
        assert!(angle < 0.001);
    }
    let rail = spline::Spline::new(
        (0..=18)
            .map(|i| [I::point(4.), I::point(i as f64), I::point(i as f64)])
            .collect(),
    );
    assert_eq!(
        rail.bounds([0.8, 0.2], &mut work).err(),
        Some(Error::Unresolved)
    );
    assert!(rail.bounds([0.1, 0.9], &mut work).is_ok());
    for side in &source.sides {
        if let Side::Cone(cone) = side {
            let split = cone.split_at([0.1, 0.8]);
            assert!(split > 0.1 && split < 0.8);
        }
    }
    let one = certify(crate::radial_tensor::tests::cam::wedge(&[1]));
    let source = Prepared::new(&one, options(), DEFAULT_LIMITS, &mut work).unwrap();
    let cone = source
        .sides
        .iter()
        .find_map(|side| {
            if let Side::Cone(c) = side {
                Some(c)
            } else {
                None
            }
        })
        .unwrap();
    let rows = cone.rows([0., 0.1], &mut work).unwrap();
    assert_eq!(
        cone.cell_error(&one, [0., 1.], &rows, [[0.; 3]; 4], false, &mut work)
            .err(),
        Some(Error::Unresolved)
    );
}

fn every_cancellation<T>(mut run: impl FnMut(&mut Work<'_>) -> Result<T, Error>) {
    // The immutable source jet has its own cold-construction cancellation
    // tests. Count every actual candidate checkpoint with a stable warm
    // source cache, so the observation itself cannot change the next run.
    let mut warm = Work {
        remaining: 10000000,
        interrupt: None,
    };
    assert!(!matches!(
        run(&mut warm),
        Err(Error::Cancelled | Error::Limit)
    ));
    let count = Cell::new(0);
    let observe = || {
        count.set(count.get() + 1);
        false
    };
    let mut work = Work {
        remaining: 10000000,
        interrupt: Some(&observe),
    };
    assert!(!matches!(
        run(&mut work),
        Err(Error::Cancelled | Error::Limit)
    ));
    let total = count.get();
    assert!(total > 0);
    for at in 1..=total {
        count.set(0);
        let cancel = || {
            count.set(count.get() + 1);
            count.get() >= at
        };
        let mut work = Work {
            remaining: 10000000,
            interrupt: Some(&cancel),
        };
        assert!(
            matches!(run(&mut work), Err(Error::Cancelled)),
            "{at}/{total}"
        );
        assert_eq!(count.get(), at);
    }
}

#[test]
fn actual_candidate_stages_cancel_at_every_checkpoint_and_refuse_displaced_caps() {
    let body = cam();
    let mut work = Work {
        remaining: 10000000,
        interrupt: None,
    };
    let mut source = Prepared::new(&body, options(), DEFAULT_LIMITS, &mut work).unwrap();
    let tensor = source
        .sides
        .iter()
        .position(|s| matches!(s, Side::Tensor(_)))
        .unwrap();
    let cone = source
        .sides
        .iter()
        .position(|s| matches!(s, Side::Cone(_)))
        .unwrap();
    let Side::Cone(c) = &source.sides[cone] else {
        unreachable!()
    };
    let row = c.row(I::point(0.25), &mut work).unwrap();
    every_cancellation(|work| source.point(cone, 0.25, 0.25, Some(&row), work));
    for side in [tensor, cone] {
        let range = [[0.1, 0.15], [0.1, 0.15]];
        let rows = match &source.sides[side] {
            Side::Cone(c) => Some(c.rows(range[1], &mut work).unwrap()),
            Side::Tensor(_) => None,
        };
        let mut corners = [[0.; 3]; 4];
        for (i, [u, v]) in [[0, 0], [1, 0], [1, 1], [0, 1]].into_iter().enumerate() {
            let row = match &source.sides[side] {
                Side::Cone(c) => Some(c.row(I::point(range[1][v]), &mut work).unwrap()),
                Side::Tensor(_) => None,
            };
            corners[i] = source
                .point(side, range[0][u], range[1][v], row.as_ref(), &mut work)
                .unwrap();
        }
        for actual in [[[0, 1, 2], [0, 2, 3]], [[2, 1, 0], [3, 2, 0]]] {
            every_cancellation(|work| {
                source.cell_matches(side, range, corners, actual, rows.as_ref(), None, work)
            });
        }
    }
    let grid = Grid {
        u: source
            .sides
            .iter()
            .map(|s| {
                if matches!(s, Side::Cone(_)) {
                    vec![0., 0.5, 1.]
                } else {
                    vec![0., 1.]
                }
            })
            .collect(),
        v: vec![0., 1.],
    };
    let ring = Ring::new(&grid);
    let candidate = source.emit(&grid, &ring, 0.001, 50000, &mut work).unwrap();
    every_cancellation(|work| {
        source.certificates.clear();
        let mut changes = Refinement::new(grid.u.len());
        source.check(&candidate, &grid, &ring, 0, &mut changes, work)
    });
    let mut displaced = candidate;
    let center = ring.entries.len() * grid.v.len();
    displaced.vertices[3 * center] += 1.;
    let mut changes = Refinement::new(grid.u.len());
    assert_eq!(
        source.check(&displaced, &grid, &ring, 0, &mut changes, &mut work),
        Err(Error::Unresolved)
    );
    let rail = Spline::new(
        (0..=18)
            .map(|i| [I::point(4.), I::point(i as f64), I::point(i as f64)])
            .collect(),
    );
    assert_eq!(
        rail.bounds([0.4, 0.3], &mut work).err(),
        Some(Error::Unresolved)
    );
    every_cancellation(|work| rail.bounds([0.1, 0.9], work));

    let calls = Cell::new(0);
    let observe = || {
        calls.set(calls.get() + 1);
        false
    };
    let mut measured = Work {
        remaining: 10000000,
        interrupt: Some(&observe),
    };
    Prepared::new(&body, options(), DEFAULT_LIMITS, &mut measured).unwrap();
    // The build owns an initial checkpoint and a candidate-attempt
    // checkpoint before the emitter's first checkpoint. Cancel there,
    // after successful preparation, and require a clean full retry.
    let at = calls.get() + 3;
    calls.set(0);
    let cancel = || {
        calls.set(calls.get() + 1);
        calls.get() >= at
    };
    assert_eq!(
        build(&body, options(), Some(&cancel)).err(),
        Some(Error::Cancelled)
    );
    assert!(build(&body, options(), None).is_ok());
}

#[test]
fn loose_actual_twisted_sources_still_require_literal_embedding() {
    for length_mm in [1.25, 1.65, 1.875, 2.5] {
        let thread = crate::tapered_thread::Thread::new(
            crate::helical_thread::Spec {
                major_radius_mm: 4.,
                pitch_mm: 1.25,
                length_mm,
                depth_mm: 0.5,
                crest_width_mm: 0.25,
                lower_flank_deg: 12.,
                upper_flank_deg: 18.,
                starts: 1,
                left_hand: false,
                root_radius_mm: None,
            },
            0.03125,
            1e-7,
            1e-6,
        )
        .unwrap();
        let body = thread
            .retained_boundary_with_interrupt(
                [[0.; 3], [1., 0., 0.], [0., 1., 0.]],
                1e-7,
                1e-6,
                None,
            )
            .unwrap()
            .boundary;
        let mut work = Work {
            remaining: 10000000,
            interrupt: None,
        };
        let source = Prepared::new(&body, options(), DEFAULT_LIMITS, &mut work).unwrap();
        let side = source
            .sides
            .iter()
            .position(|s| matches!(s, Side::Cone(_)))
            .unwrap();
        assert!(source.hint(side, [[0., 1.]; 2], 0, &mut work).unwrap() <= 1);
        let mesh = build(
            &body,
            Tessellation {
                deflection: 100.,
                angle: std::f64::consts::PI,
                triangle_limit: 200000,
            },
            None,
        )
        .unwrap();
        assert!(crate::manufacturing::certify(&mesh).is_ok());
    }
}
