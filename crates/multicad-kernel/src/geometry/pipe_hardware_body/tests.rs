use super::*;
fn spec(rounded: bool, external: bool, quarter: bool) -> AuthoredSpec {
    let k = if external { 1. / 32. } else { -1. / 32. };
    if rounded {
        let p = 25.4 / 28.;
        AuthoredSpec {
            major_radius_mm: 9.728 / 2. - 4. / 32.,
            pitch_mm: p,
            length_mm: 6.5,
            radial_taper: k,
            section:
                multicad_rust_kernel::axial_pipe_thread::authored::AuthoredSection::BsptRounded {
                    depth_mm: p / (3. * 27.5_f64.to_radians().tan()),
                },
            crest_radius_limit_mm: None,
        }
    } else {
        let (tpi, e1, l1, l2, od, max) = if quarter {
            (18., 0.49163, 0.2278, 0.4018, 0.540, 0.088)
        } else {
            (27., 0.37360, 0.1615, 0.2639, 0.405, 0.096)
        };
        let p = 25.4 / tpi;
        let m = 3_f64.sqrt();
        let s = 1. / (m - k) + 1. / (m + k);
        let d = p / s - (0.033 + max) * p;
        AuthoredSpec {
            major_radius_mm: (e1 * 25.4 - if external { 2. * k * l1 * 25.4 } else { 0. } + d) / 2.,
            pitch_mm: p,
            length_mm: l2 * 25.4,
            radial_taper: k,
            section: multicad_rust_kernel::axial_pipe_thread::authored::AuthoredSection::NptFlat {
                depth_mm: d,
                flat_width_mm: (p - d * s) / 2.,
            },
            crest_radius_limit_mm: external.then_some(od * 25.4 / 2.),
        }
    }
}
fn frame(posed: bool) -> PlaneFrame {
    if posed {
        PlaneFrame {
            origin: [10., 20., 30.],
            u: [0.8, 0.6, 0.],
            v: [-0.36, 0.48, 0.8],
        }
    } else {
        PlaneFrame::default()
    }
}
fn run(rounded: bool, external: bool, quarter: bool, posed: bool, options: MeshOptions) {
    let spec = spec(rounded, external, quarter);
    let kind = if external {
        PipeHardwareKind::HexHead {
            across_flats_mm: 16.,
            height_mm: 3.,
        }
    } else {
        PipeHardwareKind::HexNut {
            across_flats_mm: 16.,
        }
    };
    let frame = frame(posed);
    let policy = GeometryPolicy::default();
    let shape = Solid::rust_authored_axial_pipe_hardware_with_interrupt(
        spec,
        kind,
        0.2,
        EntityId::new(),
        frame,
        policy,
        None,
    )
    .unwrap();
    assert_eq!(
        shape.rust_pipe_hardware_intent().unwrap().thread().spec(),
        spec
    );
    assert!(shape.volume_mm3().unwrap() > 0.);
    let faces = shape.faces().unwrap();
    assert!(
        faces
            .iter()
            .all(|f| f.area_mm2 > 0. && f.center_mm.iter().all(|x| x.is_finite()))
    );
    let display = shape
        .topology_with_display(EntityId::new(), 0, Some(0.1))
        .unwrap();
    assert_eq!(
        display
            .topology
            .elements
            .iter()
            .filter(|e| e.key.kind == TopologyKind::Face)
            .count(),
        faces.len()
    );
    assert!(display.edges.len() > 18);
    assert!(display.vertices.len() >= 12);
    assert!(display.edges.iter().all(|e| e.points_mm.len() >= 2));
    let receipt = shape.manufacturing_mesh_with_receipt(options).unwrap();
    let mesh = &receipt.mesh;
    let v = mesh.vertices.len() / 3;
    let t = mesh.indices.len() / 3;
    assert_eq!(t, 2 * v - if external { 4 } else { 0 });
    assert_eq!(receipt.edge_pairs as usize, 3 * t / 2);
    assert_eq!(receipt.vertex_links as usize, v);
    assert!(
        mesh.face_indices
            .iter()
            .all(|i| (*i as usize) < faces.len())
    );
    let b = shape.bounds().unwrap();
    assert!(
        mesh.vertices
            .as_chunks::<3>()
            .0
            .iter()
            .all(|p| (0..3).all(|k| b.min[k] <= p[k] && p[k] <= b.max[k]))
    );
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
    assert!(volume > 0.);
    // Canonical neck is retained directly: the head shoulder's axis origin
    // must not acquire a separately rounded tip+length datum.
    if external {
        let base = shape.pipe_hardware_body().unwrap().source.faces.len() - 2;
        let shoulder = faces[base].plane.unwrap().0;
        assert!(
            (0..3).all(|k| (shoulder[k] - frame.origin[k]).abs() <= policy.linear_tolerance_mm)
        );
    }
    let local = if external {
        [0., 0., spec.length_mm / 2.]
    } else {
        [7., 0., spec.length_mm / 2.]
    };
    let source = shape.pipe_hardware_body().unwrap();
    let point = source.placement.point(local).unwrap().0;
    assert!(shape.rust_pipe_hardware_contains(point).unwrap());
    assert!(
        !shape
            .rust_pipe_hardware_contains(
                source
                    .placement
                    .point([100., 0., spec.length_mm / 2.])
                    .unwrap()
                    .0
            )
            .unwrap()
    );
    let first = shape.mesh(options).unwrap();
    let builds = shape.rust_pipe_hardware_mesh_build_count().unwrap();
    let second = shape.mesh(options).unwrap();
    assert_eq!(first.vertices, second.vertices);
    assert_eq!(first.indices, second.indices);
    assert_eq!(shape.rust_pipe_hardware_mesh_build_count(), Some(builds));
    assert!(
        shape
            .mesh(MeshOptions {
                maximum_triangles: 1,
                ..options
            })
            .is_err()
    );
    let failed = shape.rust_pipe_hardware_mesh_build_count().unwrap();
    assert_eq!(shape.mesh(options).unwrap().vertices, first.vertices);
    assert_eq!(shape.rust_pipe_hardware_mesh_build_count(), Some(failed));
    assert!(!shape.supports_step_export());
}
#[test]
fn actual_eighth_npt_head_has_coupled_neck_and_closed_native_manufacturing() {
    run(false, true, false, false, MeshOptions::default());
}
#[test]
fn actual_quarter_npt_head_has_proper_world_normals_and_closed_native_manufacturing() {
    run(false, true, true, true, MeshOptions::default());
}
#[test]
fn actual_internal_npt_hex_nut_uses_retained_flat_source_and_closed_native_manufacturing() {
    run(false, false, false, true, MeshOptions::default());
}
#[test]
fn actual_bspt_hex_head_has_proper_world_normals_and_closed_native_manufacturing() {
    run(true, true, false, true, MeshOptions::default());
}
#[test]
fn actual_bspt_hex_nut_has_retained_material_and_closed_native_manufacturing() {
    run(true, false, false, false, MeshOptions::default());
}
#[test]
fn actual_canonical_posed_npt_internal_nut_keeps_strict_manufacturing_work_ceiling() {
    use multicad_core::{
        document::Document,
        model::{
            NptGaugeConvention, NptThreadForm, NptThreadIntent, NptThreadSize, NptTruncationTarget,
            ThreadMode,
        },
        parameters,
    };
    let document = Document::new("canonical posed nut").unwrap();
    let evaluated = parameters::evaluate(&document).unwrap();
    let parameters = NptThreadIntent {
        form: NptThreadForm::Internal,
        size: NptThreadSize::OneEighth,
        gauge: NptGaugeConvention::Nbs1959PlainFace,
        truncation: NptTruncationTarget::MidpointConeParallel,
    }
    .parameters();
    let resolved = parameters
        .resolve(&evaluated, ThreadMode::CutInternal)
        .unwrap();
    let npt = resolved.npt.as_ref().unwrap();
    let spec = AuthoredSpec {
        major_radius_mm: resolved.major_radius_mm,
        pitch_mm: resolved.pitch_mm,
        length_mm: resolved.length_mm,
        radial_taper: resolved.profile.radial_taper,
        section: multicad_rust_kernel::axial_pipe_thread::authored::AuthoredSection::NptFlat {
            depth_mm: npt.depth_mm,
            flat_width_mm: npt.flat_width_mm,
        },
        crest_radius_limit_mm: npt.stock_crest_radius_limit_mm,
    };
    eprintln!("CANONICAL_NPT_NUT {spec:?}");
    let start = std::time::Instant::now();
    let shape = Solid::rust_authored_axial_pipe_hardware_with_interrupt(
        spec,
        PipeHardwareKind::HexNut {
            across_flats_mm: 18.,
        },
        resolved.minimum_wall_mm,
        EntityId::new(),
        PlaneFrame {
            origin: [11., -7., 5.],
            u: [0., 1., 0.],
            v: [0., 0., 1.],
        },
        GeometryPolicy::default(),
        None,
    )
    .unwrap();
    eprintln!("CANONICAL_NPT_NUT constructor {:?}", start.elapsed());
    let start = std::time::Instant::now();
    let receipt = shape
        .manufacturing_mesh_with_receipt(MeshOptions {
            linear_deflection_mm: 0.05,
            angular_deflection_rad: 0.3,
            maximum_triangles: 500_000,
        })
        .unwrap();
    eprintln!(
        "CANONICAL_NPT_NUT fresh {:?} facets {}",
        start.elapsed(),
        receipt.mesh.indices.len() / 3
    );
    assert!(receipt.edge_pairs > 0 && receipt.vertex_links > 0);
}

#[test]
fn actual_rotated_bspt_head_keeps_strict_model_manufacturing_quality_and_work_ceiling() {
    run(
        true,
        true,
        false,
        true,
        MeshOptions {
            linear_deflection_mm: 0.05,
            angular_deflection_rad: 0.3,
            maximum_triangles: 500_000,
        },
    );
}

fn canonical_rc_eighth_nut(posed: bool) {
    use multicad_core::{
        document::Document,
        model::{
            PipeGaugeConvention, PipeThreadForm, PipeThreadIntent, PipeThreadSize, ThreadMode,
        },
        parameters,
    };
    let document = Document::new("canonical Rc nut").unwrap();
    let evaluated = parameters::evaluate(&document).unwrap();
    let parameters = PipeThreadIntent {
        form: PipeThreadForm::Rc,
        size: PipeThreadSize::OneEighth,
        gauge: PipeGaugeConvention::Iso7Nominal1994,
    }
    .parameters();
    let resolved = parameters
        .resolve(&evaluated, ThreadMode::CutInternal)
        .unwrap();
    let spec = AuthoredSpec {
        major_radius_mm: resolved.major_radius_mm,
        pitch_mm: resolved.pitch_mm,
        length_mm: resolved.length_mm,
        radial_taper: resolved.profile.radial_taper,
        section: multicad_rust_kernel::axial_pipe_thread::authored::AuthoredSection::BsptRounded {
            depth_mm: resolved.profile.depth_mm,
        },
        crest_radius_limit_mm: None,
    };
    eprintln!(
        "CANONICAL_RC_EIGHTH {spec:?} wall={}",
        resolved.minimum_wall_mm
    );
    let shape = Solid::rust_authored_axial_pipe_hardware_with_interrupt(
        spec,
        PipeHardwareKind::HexNut {
            across_flats_mm: 18.,
        },
        resolved.minimum_wall_mm,
        EntityId::new(),
        frame(posed),
        GeometryPolicy::default(),
        None,
    )
    .unwrap();
    let receipt = shape
        .manufacturing_mesh_with_receipt(MeshOptions {
            linear_deflection_mm: 0.05,
            angular_deflection_rad: 0.3,
            maximum_triangles: 500_000,
        })
        .unwrap();
    assert!(receipt.edge_pairs > 0 && receipt.vertex_links > 0);
    let hash = receipt
        .mesh
        .vertices
        .iter()
        .flat_map(|v| v.to_bits().to_le_bytes())
        .chain(
            receipt
                .mesh
                .indices
                .iter()
                .chain(&receipt.mesh.face_indices)
                .flat_map(|v| v.to_le_bytes()),
        )
        .fold(0xcbf29ce484222325u64, |h, b| {
            (h ^ u64::from(b)).wrapping_mul(0x100000001b3)
        });
    eprintln!("CANONICAL_RC_MESH posed={posed} hash={hash:016x}");
}

#[test]
fn actual_canonical_rc_eighth_nut_keeps_strict_model_manufacturing_work_ceiling() {
    canonical_rc_eighth_nut(false);
}
#[test]
fn actual_canonical_rc_eighth_nut_arbitrary_pose_keeps_strict_model_work_ceiling() {
    canonical_rc_eighth_nut(true);
}

#[test]
fn actual_curved_face_descriptors_are_lazy_exact_and_shared_across_stock_edits() {
    let spec = spec(true, false, false);
    let owner = EntityId::new();
    let policy = GeometryPolicy::default();
    let mut factory = PipeThreadFactory::new(32 * 1024 * 1024, 16);
    factory.begin_compile();
    let started = std::time::Instant::now();
    let shape = factory
        .build(
            spec,
            PipeHardwareKind::HexNut {
                across_flats_mm: 18.,
            },
            1.,
            owner,
            frame(false),
            policy,
            None,
        )
        .unwrap();
    eprintln!(
        "LAZY_BSPT_COLD_NATIVE ms={:.6}",
        started.elapsed().as_secs_f64() * 1000.
    );
    factory.finish_compile(true);
    let body = shape.pipe_hardware_body().unwrap();
    let prepared = &body.source._prepared;
    assert_eq!(prepared.descriptor_queries.get(), 0);
    assert!(prepared.lazy_measured.borrow().is_none());
    assert!(shape.volume_mm3().unwrap() > 0.);
    shape.bounds().unwrap();
    shape.topology_snapshot(owner, 0).unwrap();
    assert_eq!(
        prepared.descriptor_queries.get(),
        0,
        "construction, bounds, volume and topology cannot materialize descriptors"
    );
    let deferred = body
        .source
        .faces
        .iter()
        .find_map(|f| match f.measure {
            FaceMeasure::Retained(i) => Some(i),
            _ => None,
        })
        .unwrap();
    assert_eq!(
        prepared.measure(deferred, Some(&|| true)).unwrap_err().kind,
        GeometryErrorKind::Cancelled
    );
    assert!(prepared.lazy_measured.borrow().is_none());
    let (boundary, _) = prepared.thread.rust_radial_boundary().unwrap();
    let calls = std::cell::Cell::new(0usize);
    let observe = || {
        calls.set(calls.get() + 1);
        false
    };
    let independent = boundary
        .face_properties_with_interrupt(source::MEASURE_ACCURACY, Some(&observe))
        .unwrap();
    let stop = calls.get() + 2; // private lazy query checkpoint, full query, publication checkpoint
    calls.set(0);
    let cancel = || {
        calls.set(calls.get() + 1);
        calls.get() >= stop
    };
    assert_eq!(
        prepared.measure(deferred, Some(&cancel)).unwrap_err().kind,
        GeometryErrorKind::Cancelled
    );
    assert!(
        prepared.lazy_measured.borrow().is_none(),
        "cancelled descriptor result cannot publish"
    );
    let faces = shape.faces().unwrap();
    assert_eq!(prepared.descriptor_queries.get(), 2);
    for (index, local) in body.source.faces.iter().enumerate() {
        if let FaceMeasure::Retained(old) = local.measure {
            let expected = independent
                .iter()
                .find(|f| f.face_index == old as usize)
                .unwrap();
            let (area, center) = prepared.measure(old, None).unwrap();
            assert_eq!([area.lo, area.hi], expected.area_mm2);
            assert_eq!(center.map(|v| [v.lo, v.hi]), expected.centroid_mm);
            assert!(faces[index].area_mm2 > 0.);
        }
    }
    assert_eq!(shape.faces().unwrap(), faces);
    assert_eq!(prepared.descriptor_queries.get(), 2);
    factory.begin_compile();
    let posed = factory
        .build(
            spec,
            PipeHardwareKind::HexNut {
                across_flats_mm: 20.,
            },
            1.,
            owner,
            frame(true),
            policy,
            None,
        )
        .unwrap();
    factory.finish_compile(true);
    assert_eq!(factory.stats().hits, 1);
    let other = &posed.pipe_hardware_body().unwrap().source._prepared;
    assert!(Rc::ptr_eq(prepared, other));
    posed.faces().unwrap();
    assert_eq!(prepared.descriptor_queries.get(), 2);
    let mut invalid_frame = frame(false);
    invalid_frame.origin = [policy.maximum_coordinate_mm; 3];
    assert!(
        factory
            .build(
                spec,
                PipeHardwareKind::HexNut {
                    across_flats_mm: 20.
                },
                1.,
                owner,
                invalid_frame,
                policy,
                None
            )
            .is_err()
    );
}
