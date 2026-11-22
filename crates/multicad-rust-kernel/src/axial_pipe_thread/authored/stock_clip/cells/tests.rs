use super::*;
#[test]
fn actual_stock_source_angular_cells_preserve_every_finite_event() {
    let mut spec = super::super::super::tests::source(1. / 32.).spec();
    spec.crest_radius_limit_mm = Some(0.405 * 25.4 / 2.);
    let source = AuthoredThread::new(spec, 1e-7, 1e-5).unwrap();
    let atlas = source
        .stock_clip_plan(StockClipOptions::default(), None)
        .unwrap();
    let plan = atlas
        .cell_plan(64, 128, 512_000_000, 64 * 1024 * 1024, None)
        .unwrap();
    eprintln!(
        "NPT CELLS V{} C{} cap{} work{} bytes{}",
        plan.vertex_count(),
        plan.cell_count(),
        plan.cap_edge_count(),
        plan.work_used(),
        plan.retained_budget_bytes()
    );
    assert_eq!(*plan.source().source(), source);
}

#[test]
fn actual_clipped_angular_cells_close_the_manufacturing_candidate() {
    let mut spec = super::super::super::tests::source(1. / 32.).spec();
    spec.crest_radius_limit_mm = Some(0.405 * 25.4 / 2.);
    let source = AuthoredThread::new(spec, 1e-7, 1e-5).unwrap();
    let atlas = source
        .stock_clip_plan(StockClipOptions::default(), None)
        .unwrap();
    let plan = atlas
        .cell_plan(128, 128, 512_000_000, 64 * 1024 * 1024, None)
        .unwrap();
    let quality = crate::display::Tessellation {
        deflection: 0.1,
        angle: 0.3,
        triangle_limit: 1_000_000,
    };
    let candidate = plan
        .mesh_candidate(quality, 512_000_000, 64 * 1024 * 1024, None)
        .unwrap();
    eprintln!(
        "NPT CLOSED V{} T{} work{} bytes{} C0{} normal{} volumes {:?} {:?}",
        candidate.mesh().vertices.len() / 3,
        candidate.mesh().indices.len() / 3,
        candidate.work_used(),
        candidate.retained_budget_bytes(),
        candidate.source_distance_mm(),
        candidate.normal_chord(),
        candidate.source_volume_mm3(),
        candidate.literal_volume_mm3()
    );
    assert_eq!(
        candidate.connectivity().vertex_links as usize,
        candidate.mesh().vertices.len() / 3
    );
    assert!(candidate.source_distance_mm() <= quality.deflection);
}
#[test]
fn actual_both_pipe_sizes_and_scales_keep_closed_source_quality_stock_and_volume() {
    for quarter in [false, true] {
        for scale in [0.1, 1., 3.] {
            let source = super::super::tests::canonical(quarter, scale);
            let atlas = source.stock_clip_plan(Default::default(), None).unwrap();
            for (deflection, angle, n) in [(0.1, 0.3, 128), (0.02, 0.5, 256)] {
                let plan = atlas
                    .cell_plan(n, 128, 512_000_000, 64 * 1024 * 1024, None)
                    .unwrap();
                let candidate = plan
                    .mesh_candidate(
                        crate::display::Tessellation {
                            deflection,
                            angle,
                            triangle_limit: 1_000_000,
                        },
                        512_000_000,
                        64 * 1024 * 1024,
                        None,
                    )
                    .unwrap();
                let data = candidate.mesh();
                let count = data.vertices.len() / 3;
                let t = data.indices.len() / 3;
                assert_eq!(candidate.connectivity().vertex_links as usize, count);
                assert_eq!(t, 2 * count - 4); // one closed genus-zero indexed shell
                assert!(candidate.source_distance_mm() <= deflection);
                assert!(candidate.normal_chord() <= 2. * (angle / 2.).sin());
                assert!(
                    candidate.work_used() <= 512_000_000
                        && candidate.retained_budget_bytes() <= 64 * 1024 * 1024
                );
                for point in data.vertices.as_chunks::<3>().0 {
                    let radial = (I::point(point[0]).square() + I::point(point[1]).square()).sqrt();
                    assert!(radial.hi <= atlas.stock_radius_limit_mm());
                }
                let v = candidate.literal_volume_mm3();
                assert!(
                    v[0] > 0.
                        && v[0] <= candidate.connectivity().signed_volume_mm3
                        && candidate.connectivity().signed_volume_mm3 <= v[1]
                );
                assert_eq!(*plan.source().source(), source);
                eprintln!(
                    "NPT MATRIX {quarter}/{scale}/{deflection}/{angle}:V{count} T{t} work{} bytes{}",
                    candidate.work_used(),
                    candidate.retained_budget_bytes()
                );
            }
        }
    }
}
#[test]
fn cell_and_mesh_public_resources_cancel_and_precision_refuse_without_changing_source() {
    let source = super::super::tests::canonical(false, 1.);
    let atlas = source.stock_clip_plan(Default::default(), None).unwrap();
    let plan = atlas
        .cell_plan(128, 128, 512_000_000, 64 * 1024 * 1024, None)
        .unwrap();
    let q = crate::display::Tessellation {
        deflection: 0.1,
        angle: 0.3,
        triangle_limit: 1_000_000,
    };
    let candidate = plan
        .mesh_candidate(q, 512_000_000, 64 * 1024 * 1024, None)
        .unwrap();
    assert!(
        atlas
            .cell_plan(
                128,
                128,
                plan.work_used(),
                plan.retained_budget_bytes(),
                None
            )
            .is_ok()
    );
    assert!(matches!(
        atlas.cell_plan(128, 128, plan.work_used() - 1, 64 * 1024 * 1024, None),
        Err(Error::Limit)
    ));
    assert!(matches!(
        atlas.cell_plan(
            128,
            128,
            512_000_000,
            plan.retained_budget_bytes() - 1,
            None
        ),
        Err(Error::Limit)
    ));
    for (n, m, work, bytes) in [
        (3, 128, 512_000_000, 64 * 1024 * 1024),
        (128, 3, 512_000_000, 64 * 1024 * 1024),
        (128, 128, 0, 64 * 1024 * 1024),
        (128, 128, 512_000_001, 64 * 1024 * 1024),
        (128, 128, 512_000_000, 0),
        (128, 128, 512_000_000, 64 * 1024 * 1024 + 1),
    ] {
        assert!(matches!(
            atlas.cell_plan(n, m, work, bytes, None),
            Err(Error::Invalid)
        ));
    }
    assert!(matches!(
        atlas.cell_plan(128, 128, 512_000_000, 64 * 1024 * 1024, Some(&|| true)),
        Err(Error::Cancelled)
    ));
    assert!(
        plan.mesh_candidate(
            q,
            candidate.work_used(),
            candidate.retained_budget_bytes(),
            None
        )
        .is_ok()
    );
    assert!(matches!(
        plan.mesh_candidate(q, candidate.work_used() - 1, 64 * 1024 * 1024, None),
        Err(Error::Limit)
    ));
    assert!(matches!(
        plan.mesh_candidate(q, 512_000_000, candidate.retained_budget_bytes() - 1, None),
        Err(Error::Limit)
    ));
    assert!(matches!(
        plan.mesh_candidate(q, 512_000_000, 64 * 1024 * 1024, Some(&|| true)),
        Err(Error::Cancelled)
    ));
    assert!(matches!(
        plan.mesh_candidate(
            crate::display::Tessellation {
                triangle_limit: 1,
                ..q
            },
            512_000_000,
            64 * 1024 * 1024,
            None
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        plan.mesh_candidate(
            crate::display::Tessellation {
                deflection: f64::from_bits(1),
                ..q
            },
            512_000_000,
            64 * 1024 * 1024,
            None
        ),
        Err(Error::Unresolved)
    ));
    assert_eq!(*atlas.source(), source);
    assert_eq!(
        plan.mesh_candidate(q, 512_000_000, 64 * 1024 * 1024, None)
            .unwrap()
            .mesh(),
        candidate.mesh()
    );
}
#[test]
fn exact_contact_and_root_grid_coincidences_retain_the_original_event() {
    let source = AuthoredThread::new(
        AuthoredSpec {
            major_radius_mm: 4.,
            pitch_mm: 1.,
            length_mm: 4.,
            radial_taper: 1. / 32.,
            section: AuthoredSection::NptFlat {
                depth_mm: 0.5,
                flat_width_mm: 0.1,
            },
            crest_radius_limit_mm: Some(4.0625),
        },
        1e-7,
        1e-5,
    )
    .unwrap();
    let atlas = source.stock_clip_plan(Default::default(), None).unwrap();
    let event = atlas.contact_height_bounds_mm();
    assert!(event[0] <= 2. && event[1] >= 2.);
    let plan = atlas
        .cell_plan(128, 128, 512_000_000, 64 * 1024 * 1024, None)
        .unwrap();
    let mesh = plan
        .mesh_candidate(
            crate::display::Tessellation {
                deflection: 0.1,
                angle: 0.3,
                triangle_limit: 1_000_000,
            },
            512_000_000,
            64 * 1024 * 1024,
            None,
        )
        .unwrap();
    assert_eq!(
        mesh.connectivity().vertex_links as usize,
        mesh.mesh().vertices.len() / 3
    );
    assert_eq!(atlas.contact_height_bounds_mm(), event);
    assert_eq!(*plan.source().source(), source);
}
#[test]
fn source_cells_keep_exact_owners_closed_edges_and_original_radial_field() {
    use std::collections::{BTreeMap, BTreeSet};
    let source = super::super::tests::canonical(false, 1.);
    let atlas = source.stock_clip_plan(Default::default(), None).unwrap();
    let plan = atlas
        .cell_plan(128, 128, 512_000_000, 64 * 1024 * 1024, None)
        .unwrap();
    let mut uses: BTreeMap<(u32, u32), Vec<(u32, u32)>> = BTreeMap::new();
    let mut originals = BTreeSet::new();
    for index in 0..plan.vertex_count() {
        if let StockClipPointOwner::Vertex(key) = plan.vertex_owner(index).unwrap() {
            assert!(originals.insert(key));
            let original = atlas.vertex_bounds(key, Default::default(), None).unwrap();
            let actual = plan.vertices[index].position;
            for (i, coordinate) in actual.iter().enumerate() {
                let expected = original.position_mm()[i];
                assert!(coordinate.lo <= expected[1] && expected[0] <= coordinate.hi);
            }
        }
    }
    assert_eq!(originals, atlas.vertices().into_iter().collect());
    assert!(plan.vertex_owner(plan.vertex_count()).is_none());
    assert!(plan.cell_grid(plan.cell_count()).is_none());
    assert!(plan.cell_patch(plan.cell_count()).is_none());
    assert!(plan.cell_vertices(plan.cell_count()).is_none());
    for index in 0..plan.cell_count() {
        let cell = &plan.cells[index];
        assert_eq!(plan.cell_patch(index), Some(cell.patch));
        let (col, row) = plan.cell_grid(index).unwrap();
        let [angular, axial] = plan.subdivisions();
        let ids = plan.cell_vertices(index).unwrap();
        assert_eq!(ids.len(), cell.count);
        for i in 0..cell.count {
            let a = ids[i];
            let b = ids[(i + 1) % cell.count];
            uses.entry((a.min(b), a.max(b))).or_default().push((b, a));
        }
        // The SAME clipped cell owns the original unwrapped phase, including
        // root-period aliases; compare its interior to the independent old radial
        // source evaluator, rather than only to the new field implementation.
        if index % 37 == 0 {
            let s = cell.parameters[..cell.count]
                .iter()
                .map(|p| (p[0].lo + p[0].hi) * 0.5)
                .sum::<f64>()
                / cell.count as f64;
            let z = cell.parameters[..cell.count]
                .iter()
                .map(|p| (p[1].lo + p[1].hi) * 0.5)
                .sum::<f64>()
                / cell.count as f64;
            let phase = (s + z).rem_euclid(source.spec.pitch_mm);
            let radius = source.radial_bounds([phase, phase], [z, z]).unwrap();
            let actual = ids
                .iter()
                .map(|id| {
                    let r = plan.vertices[*id as usize].radius;
                    (r.lo + r.hi) * 0.5
                })
                .sum::<f64>()
                / cell.count as f64;
            assert!((actual - (radius[0] + radius[1]) * 0.5).abs() < 1e-9);
            assert!(
                s >= (atlas.q0
                    + I::point(col as f64) * I::point(source.spec.pitch_mm)
                        / I::point(angular as f64))
                .lo - 1e-12
            );
            assert!(z >= row as f64 * source.spec.length_mm / axial as f64 - 1e-12);
        }
    }
    let mut successor = [BTreeMap::new(), BTreeMap::new()];
    for edge in &plan.caps {
        let [a, b] = edge.ids;
        uses.entry((a.min(b), a.max(b))).or_default().push((a, b));
        let side = usize::from(edge.station == StockClipStation::End);
        assert!(successor[side].insert(a, b).is_none());
    }
    for chain in &successor {
        let start = *chain.keys().next().unwrap();
        let mut next = start;
        let mut seen = BTreeSet::new();
        loop {
            assert!(seen.insert(next));
            next = chain[&next];
            if next == start {
                break;
            }
        }
        assert_eq!(seen.len(), chain.len());
    }
    for pair in uses.values() {
        assert_eq!(pair.len(), 2);
        assert_eq!(pair[0], (pair[1].1, pair[1].0));
    }
    assert_eq!(
        plan.vertex_count() as isize - uses.len() as isize + plan.cell_count() as isize + 2,
        2
    );
    assert_eq!(*plan.source().source(), source);
}

#[test]
fn actual_public_construction_and_final_volume_stages_share_one_budget() {
    let source = super::super::tests::canonical(false, 1.);
    let atlas = source.stock_clip_plan(Default::default(), None).unwrap();
    for work in [1, 640, 641, 768, 769] {
        assert!(matches!(
            atlas.cell_plan(4, 4, work, 64 * 1024 * 1024, None),
            Err(Error::Limit)
        ));
    }
    assert!(matches!(
        atlas.cell_plan(128, 128, 512_000_000, 1, None),
        Err(Error::Limit)
    ));
    let plan = atlas
        .cell_plan(128, 128, 512_000_000, 64 * 1024 * 1024, None)
        .unwrap();
    let q = crate::display::Tessellation {
        deflection: 0.1,
        angle: 0.3,
        triangle_limit: 1_000_000,
    };
    let candidate = plan
        .mesh_candidate(q, 512_000_000, 64 * 1024 * 1024, None)
        .unwrap();
    for delta in [1, 2, 3, 512] {
        assert!(matches!(
            plan.mesh_candidate(q, candidate.work_used() - delta, 64 * 1024 * 1024, None),
            Err(Error::Limit)
        ));
    }
    for (work, bytes) in [
        (0, 64 * 1024 * 1024),
        (512_000_001, 64 * 1024 * 1024),
        (512_000_000, 0),
        (512_000_000, 64 * 1024 * 1024 + 1),
    ] {
        assert!(matches!(
            plan.mesh_candidate(q, work, bytes, None),
            Err(Error::Invalid)
        ));
    }
    assert!(matches!(
        plan.mesh_candidate(q, 512_000_000, 1, None),
        Err(Error::Limit)
    ));
    assert!(matches!(
        plan.mesh_candidate(
            crate::display::Tessellation {
                deflection: 0.,
                ..q
            },
            512_000_000,
            64 * 1024 * 1024,
            None
        ),
        Err(Error::Invalid)
    ));
    assert!(
        plan.mesh_candidate(
            crate::display::Tessellation {
                deflection: f64::MAX,
                ..q
            },
            512_000_000,
            64 * 1024 * 1024,
            None
        )
        .is_ok()
    );
    assert_eq!(*plan.source().source(), source);
}
#[test]
fn direct_rail_inverses_keep_tiny_slabs_and_original_curve_tokens() {
    let mut spec = super::super::tests::canonical(false, 1.).spec();
    spec.crest_radius_limit_mm = Some(spec.major_radius_mm.next_up());
    let source = AuthoredThread::new(spec, 1e-7, 1e-5).unwrap();
    let atlas = source.stock_clip_plan(Default::default(), None).unwrap();
    let plan = atlas
        .cell_plan(128, 128, 512_000_000, 64 * 1024 * 1024, None)
        .unwrap();
    assert!(
        plan.cells
            .iter()
            .any(|c| c.patch.slab == StockClipSlab::Before)
    );
    assert!(atlas.contact_height_bounds_mm()[0] > 0.);
    assert_eq!(*plan.source().source(), source);
    for curve in atlas.curves() {
        if let StockClipCurve::Axial { slab, rail } = curve {
            let (a, rate) = rail_law(&atlas, slab, rail, false);
            assert!(rate.lo >= 0.98 && rate.hi <= 1.02);
            for t in [0., 0.3, 1.] {
                let p = atlas
                    .curve_bounds(curve, [t, t], Default::default(), None)
                    .unwrap();
                let z = b(p.axial_mm());
                let q = a - rate * z + z;
                let actual = b(p.phase_mm());
                assert!(q.lo <= actual.hi && actual.lo <= q.hi);
            }
        }
    }
}
#[test]
fn real_auxiliary_cut_uncertainty_refuses_then_retries_the_unchanged_source() {
    let k = 1. / 32.;
    let depth = (3_f64.sqrt() - k) * 9. / 128.;
    let source = AuthoredThread::new(
        AuthoredSpec {
            major_radius_mm: 4.,
            pitch_mm: 1.,
            length_mm: 4.,
            radial_taper: k,
            section: AuthoredSection::NptFlat {
                depth_mm: depth,
                flat_width_mm: 0.1,
            },
            crest_radius_limit_mm: Some(4.0625),
        },
        1e-7,
        1e-5,
    )
    .unwrap();
    let atlas = source.stock_clip_plan(Default::default(), None).unwrap();
    // The actual irrational source station is within the arithmetic enclosure
    // of this optional angular cut; do not weld/snap them based on overlap.
    assert!(matches!(
        atlas.cell_plan(128, 128, 512_000_000, 64 * 1024 * 1024, None),
        Err(Error::Unresolved)
    ));
    let retry = atlas
        .cell_plan(64, 128, 512_000_000, 64 * 1024 * 1024, None)
        .unwrap();
    assert_eq!(*retry.source().source(), source);
    assert_eq!(atlas.source().spec(), source.spec());
}
#[test]
fn late_source_and_literal_cancellation_stop_before_publication() {
    use std::sync::atomic::{AtomicUsize, Ordering};
    let source = super::super::tests::canonical(false, 1.);
    let atlas = source.stock_clip_plan(Default::default(), None).unwrap();
    let calls = AtomicUsize::new(0);
    let poll = || {
        calls.fetch_add(1, Ordering::Relaxed);
        false
    };
    let plan = atlas
        .cell_plan(128, 128, 512_000_000, 64 * 1024 * 1024, Some(&poll))
        .unwrap();
    let total = calls.load(Ordering::Relaxed);
    assert!(total > 1000);
    calls.store(0, Ordering::Relaxed);
    let late = || calls.fetch_add(1, Ordering::Relaxed) >= total - 1;
    assert!(matches!(
        atlas.cell_plan(128, 128, 512_000_000, 64 * 1024 * 1024, Some(&late)),
        Err(Error::Cancelled)
    ));
    calls.store(0, Ordering::Relaxed);
    let quality = crate::display::Tessellation {
        deflection: 0.1,
        angle: 0.3,
        triangle_limit: 1_000_000,
    };
    let candidate = plan
        .mesh_candidate(quality, 512_000_000, 64 * 1024 * 1024, Some(&poll))
        .unwrap();
    let total = calls.load(Ordering::Relaxed);
    calls.store(0, Ordering::Relaxed);
    let late = || calls.fetch_add(1, Ordering::Relaxed) >= total - 1;
    assert!(matches!(
        plan.mesh_candidate(quality, 512_000_000, 64 * 1024 * 1024, Some(&late)),
        Err(Error::Cancelled)
    ));
    assert_eq!(
        plan.mesh_candidate(quality, 512_000_000, 64 * 1024 * 1024, None)
            .unwrap()
            .mesh(),
        candidate.mesh()
    );
    assert_eq!(*plan.source().source(), source);
}
#[test]
fn exact_product_witness_requires_the_full_dyadic_tail_not_rounded_zero() {
    let tiny = f64::from_bits(1);
    assert_eq!(exact_product(tiny, 0.5), None);
    assert_eq!(exact_product(tiny, 1.5), None);
    assert!(!same_products([tiny, 0.5], [tiny, 0.5])); // two failures are not equality
    assert_eq!(exact_product(tiny, 1.), Some([tiny, 0.]));
    assert!(same_products([tiny, 1.], [tiny, 1.]));
    assert_eq!(
        exact_product(f64::MIN_POSITIVE, 0.5),
        Some([f64::MIN_POSITIVE * 0.5, 0.])
    );
    assert_eq!(exact_product(0., tiny), Some([0., 0.]));
    assert_eq!(exact_product(f64::MAX, 2.), None);
    assert_eq!(exact_product(f64::NAN, 0.), None);
    let step = 2f64.powi(-52);
    assert_eq!(
        exact_product(1. + step, 1. - step),
        Some([1., -2f64.powi(-104)])
    );
    // The smallest admitted source pitch/length, largest grid divisors and
    // taper factor stay at least 996 binary orders above the excluded tail.
    for value in [0.001, 1., 1e7] {
        for divisor in [4., 512., 2048.] {
            let a = value / divisor / 32.;
            assert!(lowest_bit(a) >= -78);
            for index in [-66560., -1., 0., 1., 66560.] {
                assert!(exact_product(a, index).is_some());
            }
        }
    }
}
#[test]
fn real_npt_allowance_edits_certify_whole_alternate_fans_and_same_owned_arrays() {
    for quarter in [false, true] {
        let mut spec = super::super::tests::canonical(quarter, 1.).spec();
        spec.major_radius_mm -= 0.001;
        let source = AuthoredThread::new(spec, 1e-7, 1e-5).unwrap();
        let atlas = source.stock_clip_plan(Default::default(), None).unwrap();
        let plan = atlas
            .cell_plan(128, 128, 512_000_000, 64 * 1024 * 1024, None)
            .unwrap();
        let q = crate::display::Tessellation {
            deflection: 0.1,
            angle: 15f64.to_radians(),
            triangle_limit: 2_000_000,
        };
        let c = plan
            .mesh_candidate(q, 512_000_000, 64 * 1024 * 1024, None)
            .unwrap();
        assert!(c.source_distance_mm() <= q.deflection);
        assert!(c.normal_chord() <= 2. * (q.angle / 2.).sin());
        assert!(c.connectivity().edge_pairs > 0 && c.connectivity().vertex_links > 0);
        assert_eq!(*plan.source().source(), source);
        // Real public work frontier: prerequisite, literal vertices and the
        // first aspect proposal fit; the first complete source/facet proof
        // does not. A failed fan never receives a renewed ledger.
        let first_count = plan.cell_vertices(0).unwrap().len();
        let triangles = (0..plan.cell_count())
            .map(|i| plan.cell_vertices(i).unwrap().len() - 2)
            .sum::<usize>()
            + plan.cap_edge_count();
        let prefix = plan.work_used()
            + (plan.vertex_count() + 2 + triangles) * 32
            + plan.vertex_count() * 256
            + first_count * first_count * 64;
        assert!(matches!(
            plan.mesh_candidate(q, prefix + 1023, 64 * 1024 * 1024, None),
            Err(Error::Limit)
        ));
        let exact = c.work_used();
        assert!(matches!(
            plan.mesh_candidate(q, exact - 1, 64 * 1024 * 1024, None),
            Err(Error::Limit)
        ));
        let retry = plan
            .mesh_candidate(q, exact, 64 * 1024 * 1024, None)
            .unwrap();
        assert_eq!(c.mesh().vertices, retry.mesh().vertices);
        assert_eq!(c.mesh().indices, retry.mesh().indices);
        assert_eq!(c.mesh().face_indices, retry.mesh().face_indices);
        let pointers = (
            c.mesh().vertices.as_ptr(),
            c.mesh().indices.as_ptr(),
            c.mesh().face_indices.as_ptr(),
        );
        let mesh = c.into_mesh();
        assert_eq!(
            pointers,
            (
                mesh.vertices.as_ptr(),
                mesh.indices.as_ptr(),
                mesh.face_indices.as_ptr()
            )
        );
        let proof = crate::manufacturing::certify(&mesh).unwrap();
        assert_eq!(proof.vertex_links as usize, mesh.vertices.len() / 3);
        assert!(matches!(
            plan.mesh_candidate(
                crate::display::Tessellation { angle: 0.01, ..q },
                512_000_000,
                64 * 1024 * 1024,
                None
            ),
            Err(Error::Unresolved)
        ));
    }
}
