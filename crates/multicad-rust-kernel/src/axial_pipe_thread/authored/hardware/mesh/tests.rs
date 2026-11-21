use super::*;
fn quality() -> Tessellation {
    Tessellation {
        deflection: 0.05,
        angle: 0.3,
        triangle_limit: 1_000_000,
    }
}
fn options() -> MeshOptions {
    MeshOptions {
        tessellation: Tessellation {
            deflection: 0.1,
            angle: 0.4,
            ..quality()
        },
        maximum_work: 128_000_000,
        maximum_bytes: 64 * 1024 * 1024,
        additional_resident_bytes: 0,
    }
}
fn verify(source: ThreadBoundaryMesh<'_, '_>, kind: Kind) {
    let result = join(source, kind, options(), None).unwrap();
    let mesh = result.mesh();
    let vertices = mesh.vertices.len() / 3;
    let triangles = mesh.indices.len() / 3;
    let nut = matches!(kind, Kind::HexNut { .. });
    assert_eq!(triangles, 2 * vertices - if nut { 0 } else { 4 });
    assert_eq!(result.connectivity().edge_pairs as usize, 3 * triangles / 2);
    assert_eq!(result.connectivity().vertex_links as usize, vertices);
    assert!(result.source_distance_mm() <= options().tessellation.deflection);
    assert!(result.joining_work_used() <= options().maximum_work);
    assert_eq!(result.inherited_source_work_used(), source.work());
    assert!(result.accounted_peak_bytes() <= options().maximum_bytes);
    let inherited = match result.source() {
        ThreadBoundaryMesh::Clipped(t) => t.work_used(),
        ThreadBoundaryMesh::Retained(t) => t.work_used(),
    };
    assert_eq!(inherited, source.work());
    // Resource refusal follows the actual counted algorithm, rather than
    // assuming the former128-per-index upper estimate consumes64M.
    let required = result.joining_work_used();
    assert!(matches!(
        join(
            source,
            kind,
            MeshOptions {
                maximum_work: required - 1,
                ..options()
            },
            None
        ),
        Err(Error::Limit)
    ));
    let exact = join(
        source,
        kind,
        MeshOptions {
            maximum_work: required,
            ..options()
        },
        None,
    )
    .unwrap();
    assert_eq!(exact.mesh(), mesh);
    assert_eq!(exact.connectivity(), result.connectivity());
    let volume = mesh
        .indices
        .as_chunks::<3>()
        .0
        .iter()
        .map(|ids| {
            let [a, b, c] = ids.map(|id| &mesh.vertices[3 * id as usize..][..3]);
            (a[0] * (b[1] * c[2] - b[2] * c[1])
                + a[1] * (b[2] * c[0] - b[0] * c[2])
                + a[2] * (b[0] * c[1] - b[1] * c[0]))
                / 6.
        })
        .sum::<f64>();
    assert!(volume > 0. && volume.is_finite());
    eprintln!(
        "PIPE HEX JOIN nut{nut} V{vertices} T{triangles} volume{volume} work{} bytes{}",
        result.joining_work_used(),
        result.accounted_peak_bytes()
    );
    assert_eq!(result.into_mesh().indices.len(), 3 * triangles);
}
fn clipped(quarter: bool) {
    let source = super::super::super::stock_clip::tests::canonical(quarter, 1.);
    let plan = source
        .stock_clip_plan(Default::default(), None)
        .unwrap()
        .cell_plan(256, 128, 512_000_000, 64 * 1024 * 1024, None)
        .unwrap();
    let trace = plan
        .mesh_candidate_with_cap_trace(quality(), 512_000_000, 64 * 1024 * 1024, None)
        .unwrap();
    verify(
        ThreadBoundaryMesh::Clipped(&trace),
        Kind::HexHead {
            across_flats_mm: 16.,
            height_mm: 3.,
        },
    );
}
#[test]
fn actual_eighth_clipped_thread_joins_to_closed_hex_head() {
    clipped(false);
}
#[test]
fn actual_quarter_clipped_thread_joins_to_closed_hex_head() {
    clipped(true);
}
fn retained(external: bool) {
    let pitch = 25.4 / 28.;
    let source = AuthoredThread::new(
        AuthoredSpec {
            major_radius_mm: 9.728 / 2. - 4. / 32.,
            pitch_mm: pitch,
            length_mm: 6.5,
            radial_taper: if external { 1. / 32. } else { -1. / 32. },
            section: AuthoredSection::BsptRounded {
                depth_mm: pitch / (3. * 27.5_f64.to_radians().tan()),
            },
            crest_radius_limit_mm: None,
        },
        1e-7,
        1e-6,
    )
    .unwrap();
    let raw = source.raw_boundary_ownership(1e-6, None).unwrap();
    let body = raw
        .ownership
        .bind_surfaces(1e-6, None)
        .and_then(|b| b.certify_axial(1e-6, 1e-6, None))
        .and_then(|b| b.certify_caps(None))
        .and_then(|b| b.certify_ordering(None))
        .and_then(|b| b.certify_material(1e-6, None))
        .unwrap();
    let trace = body
        .tessellate_with_cap_trace(quality(), Default::default(), None)
        .unwrap();
    verify(
        ThreadBoundaryMesh::Retained(&trace),
        if external {
            Kind::HexHead {
                across_flats_mm: 16.,
                height_mm: 3.,
            }
        } else {
            Kind::HexNut {
                across_flats_mm: 16.,
            }
        },
    );
}
#[test]
fn actual_original_rounded_external_thread_joins_to_closed_hex_head() {
    retained(true);
}
#[test]
fn actual_original_rounded_internal_thread_joins_to_closed_through_hex_nut() {
    retained(false);
}
#[test]
fn actual_join_refuses_contact_precision_and_budgets_without_mutating_original_mesh() {
    let source = super::super::super::stock_clip::tests::canonical(false, 1.);
    let plan = source
        .stock_clip_plan(Default::default(), None)
        .unwrap()
        .cell_plan(256, 128, 512_000_000, 64 * 1024 * 1024, None)
        .unwrap();
    let trace = plan
        .mesh_candidate_with_cap_trace(quality(), 512_000_000, 64 * 1024 * 1024, None)
        .unwrap();
    let source = ThreadBoundaryMesh::Clipped(&trace);
    let original = trace.candidate().mesh().clone();
    let kind = Kind::HexHead {
        across_flats_mm: 16.,
        height_mm: 3.,
    };
    for bad in [
        Kind::HexHead {
            across_flats_mm: f64::NAN,
            height_mm: 3.,
        },
        Kind::HexHead {
            across_flats_mm: 16.,
            height_mm: 0.,
        },
        Kind::HexNut {
            across_flats_mm: -1.,
        },
    ] {
        assert!(matches!(
            join(source, bad, options(), None),
            Err(Error::Invalid)
        ));
    }
    for bad in [
        Kind::HexHead {
            across_flats_mm: 9.,
            height_mm: 3.,
        },
        Kind::HexHead {
            across_flats_mm: 16.,
            height_mm: f64::MIN_POSITIVE,
        },
        Kind::HexHead {
            across_flats_mm: 16.,
            height_mm: 1e308,
        },
        Kind::HexHead {
            across_flats_mm: 1e308,
            height_mm: 3.,
        },
    ] {
        assert!(matches!(
            join(source, bad, options(), None),
            Err(Error::Unsupported) | Err(Error::Unresolved)
        ));
    }
    for query in [
        MeshOptions {
            maximum_work: 1,
            ..options()
        },
        MeshOptions {
            maximum_bytes: 1,
            ..options()
        },
        MeshOptions {
            additional_resident_bytes: usize::MAX,
            ..options()
        },
        MeshOptions {
            tessellation: Tessellation {
                triangle_limit: 1,
                ..options().tessellation
            },
            ..options()
        },
    ] {
        assert!(matches!(join(source, kind, query, None), Err(Error::Limit)));
    }
    for query in [
        MeshOptions {
            maximum_work: 0,
            ..options()
        },
        MeshOptions {
            maximum_work: 512_000_001,
            ..options()
        },
        MeshOptions {
            maximum_bytes: 0,
            ..options()
        },
        MeshOptions {
            maximum_bytes: 64 * 1024 * 1024 + 1,
            ..options()
        },
    ] {
        assert!(matches!(
            join(source, kind, query, None),
            Err(Error::Invalid)
        ));
    }
    let bad = MeshOptions {
        tessellation: Tessellation {
            deflection: f64::NAN,
            ..options().tessellation
        },
        ..options()
    };
    assert!(matches!(join(source, kind, bad, None), Err(Error::Invalid)));
    let too_precise = MeshOptions {
        tessellation: Tessellation {
            deflection: 1e-8,
            ..options().tessellation
        },
        ..options()
    };
    assert!(matches!(
        join(source, kind, too_precise, None),
        Err(Error::Unsupported)
    ));
    assert!(matches!(
        join(source, kind, options(), Some(&|| true)),
        Err(Error::Cancelled)
    ));
    let calls = std::cell::Cell::new(0usize);
    let cancel = || {
        calls.set(calls.get() + 1);
        calls.get() > 50
    };
    assert!(matches!(
        join(source, kind, options(), Some(&cancel)),
        Err(Error::Cancelled)
    ));
    assert_eq!(calls.get(), 51);
    assert_eq!(trace.candidate().mesh(), &original);
    assert!(join(source, kind, options(), None).is_ok());
}

#[test]
fn through_bore_joins_to_closed_cylindrical_stock_with_exact_annular_caps() {
    use crate::display::Tessellation;
    let cavity = crate::step::radial_boundary_tests::internal_npt_cavity();
    let tessellation = Tessellation {
        deflection: 0.05,
        angle: 0.3,
        triangle_limit: 1_000_000,
    };
    let trace = cavity
        .tessellate_with_cap_trace(tessellation, Default::default(), None)
        .unwrap();
    let options = MeshOptions {
        tessellation,
        maximum_work: 512_000_000,
        maximum_bytes: 64 * 1024 * 1024,
        additional_resident_bytes: 0,
    };
    let radius = 7.;
    let joined = join_sleeve(&trace, radius, options, None).unwrap();
    let mesh = joined.mesh();
    let n = sleeve_segments(radius, tessellation).unwrap();
    assert_eq!(
        mesh.face_indices.len(),
        trace.mesh().face_indices.len() + 4 * n
    );
    // Closed, consistently oriented and positive: stock minus bore.
    let proof = crate::manufacturing::certify(mesh).unwrap();
    let bore = cavity.volume_with_interrupt(1e-6, None).unwrap().volume_mm3;
    let polygon = 0.5 * n as f64 * radius * radius * (std::f64::consts::TAU / n as f64).sin();
    let expected = polygon * 8. - (bore[0] + bore[1]) / 2.;
    let source_area: f64 = cavity
        .face_properties_with_interrupt(1e-4, None)
        .unwrap()
        .iter()
        .map(|f| f.area_mm2[1])
        .sum();
    assert!(
        (proof.signed_volume_mm3 - expected).abs() <= tessellation.deflection * source_area,
        "{} {expected}",
        proof.signed_volume_mm3
    );
    // Wall vertices lie on the stock circle, at the cap heights only.
    let wall = trace.mesh().face_indices.iter().max().unwrap() + 1;
    for (t, &face) in mesh
        .indices
        .as_chunks::<3>()
        .0
        .iter()
        .zip(&mesh.face_indices)
    {
        if face == wall {
            for &v in t {
                let p = &mesh.vertices[3 * v as usize..3 * v as usize + 3];
                assert!((p[0].hypot(p[1]) - radius).abs() < 1e-12);
                assert!(p[2] == 0. || p[2] == 8.);
            }
        }
    }
    assert!(joined.source_distance_mm() <= tessellation.deflection);
    // Stock that does not clear the bore, and invalid stock, are refused.
    assert_eq!(
        join_sleeve(&trace, 4.9, options, None).err(),
        Some(Error::Unsupported)
    );
    assert_eq!(
        join_sleeve(&trace, f64::NAN, options, None).err(),
        Some(Error::Invalid)
    );
    assert_eq!(
        join_sleeve(&trace, radius, options, Some(&|| true)).err(),
        Some(Error::Cancelled)
    );
}
