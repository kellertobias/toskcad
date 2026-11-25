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
                    depth_mm: p / (3. * 27.5_f64.to_radians().portable_tan()),
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
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("hardware.step");
    std::fs::write(&path, b"existing output").unwrap();
    // Clipped NPT hex heads are one mapped item requalified from their bytes;
    // BSPT heads and every nut keep their retained thread tables beside exact
    // hex stock, decoded back bit for bit.
    let mapped = !rounded && external;
    assert!(shape.supports_step_export());
    let written = shape
        .write_pipe_hardware_step_with_verification(&path, 1e-7, None)
        .unwrap()
        .unwrap();
    let actual = std::fs::read_to_string(&path).unwrap();
    assert_eq!(actual.contains("MAPPED_ITEM"), mapped);
    assert_eq!(
        written.method,
        if mapped {
            super::PipeHardwareStepMethod::MappedClippedHead
        } else {
            super::PipeHardwareStepMethod::RetainedHexStock
        }
    );
    let fresh = shape
        .verify_pipe_hardware_step_source(&actual, 1e-7, None)
        .unwrap()
        .unwrap();
    assert_eq!(fresh.world_sewn_volume_mm3, written.world_sewn_volume_mm3);
    assert!(fresh.world_sewn_reference_position_error_mm <= 1e-7);
    assert!(fresh.head_author_corner_error_mm <= 1e-7);
    assert!(fresh.factor_chain_discrepancy_mm <= 1e-7);
    let [lo, hi] = fresh.world_sewn_volume_mm3;
    assert!(0. < lo && lo <= hi);
    if !mapped {
        let volume = shape.volume_mm3().unwrap();
        assert!(lo <= volume && volume <= hi);
    }
    // One micrometre of extra occurrence is a different source, a bad
    // tolerance is refused, and a failed write keeps the existing destination.
    let moved = shape
        .transformed(crate::RigidTransform {
            translation_mm: [1e-3, 0., 0.],
            ..Default::default()
        })
        .unwrap();
    assert!(
        moved
            .verify_pipe_hardware_step_source(&actual, 1e-7, None)
            .unwrap()
            .is_err()
    );
    assert!(
        shape
            .verify_pipe_hardware_step_source(&actual, 0., None)
            .unwrap()
            .is_err()
    );
    let kept = dir.path().join("kept.step");
    std::fs::write(&kept, b"existing output").unwrap();
    assert!(
        shape
            .write_pipe_hardware_step_with_verification(&kept, 1e-7, Some(&|| true))
            .unwrap()
            .is_err()
    );
    assert_eq!(std::fs::read(&kept).unwrap(), b"existing output");
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

// Independent numerical Green integrals over the actual retained rim. Degree-9
// spline positions/derivatives use the explicit Bernstein formula, and circles
// use their literal axes and endpoint sweep. This is a test oracle, not a
// geometry receipt; step doubling must converge before comparison to the
// independently certified cap enclosures.
fn cap_green_oracle(
    input: &multicad_rust_kernel::radial_tensor::BoundaryInput,
    cap: usize,
    steps: usize,
) -> (f64, [f64; 3]) {
    use multicad_rust_kernel::radial_tensor::{Curve, Surface};
    fn bernstein(p: &[[f64; 3]], t: f64) -> [f64; 3] {
        let weights: &[f64] = match p.len() {
            10 => &[1., 9., 36., 84., 126., 126., 84., 36., 9., 1.],
            9 => &[1., 8., 28., 56., 70., 56., 28., 8., 1.],
            _ => panic!("actual retained cap degree"),
        };
        std::array::from_fn(|k| {
            p.iter()
                .enumerate()
                .map(|(i, q)| {
                    weights[i] * t.powi(i as i32) * (1. - t).powi((p.len() - 1 - i) as i32) * q[k]
                })
                .sum()
        })
    }
    let mut sums = [0.; 3];
    let mut integrate = |evaluate: &dyn Fn(f64) -> ([f64; 3], [f64; 3])| {
        for i in 0..=steps {
            let (p, d) = evaluate(i as f64 / steps as f64);
            let weight = if i == 0 || i == steps {
                1.
            } else if i % 2 == 0 {
                2.
            } else {
                4.
            };
            for (sum, value) in sums.iter_mut().zip([
                (p[0] * d[1] - p[1] * d[0]) / 2.,
                p[0] * p[0] * d[1] / 2.,
                -p[1] * p[1] * d[0] / 2.,
            ]) {
                *sum += weight * value / (3. * steps as f64);
            }
        }
    };
    for &(edge, direction) in &input.faces[cap].outer {
        let edge = input.edges[edge];
        match &input.curves[edge.curve] {
            Curve::Spline { spans, poles } => {
                let mut poles = poles.clone();
                if edge.forward != direction {
                    poles.reverse();
                }
                for span in 0..*spans {
                    let p = &poles[9 * span..=9 * span + 9];
                    let derivative = p
                        .windows(2)
                        .map(|p| std::array::from_fn(|k| 9. * (p[1][k] - p[0][k])))
                        .collect::<Vec<_>>();
                    integrate(&|t| (bernstein(p, t), bernstein(&derivative, t)));
                }
            }
            Curve::Circle { axis, radius_mm } => {
                assert!(axis.normal[0] == 0. && axis.normal[1] == 0. && axis.normal[2].abs() == 1.);
                assert_eq!(axis.reference[2], 0.);
                let norm = axis.reference[0].portable_hypot(axis.reference[1]);
                let u = [axis.reference[0] / norm, axis.reference[1] / norm, 0.];
                let v = [-axis.normal[2] * u[1], axis.normal[2] * u[0], 0.];
                let [a, b] = if direction {
                    edge.vertices
                } else {
                    [edge.vertices[1], edge.vertices[0]]
                };
                let phase = |index: usize| {
                    let p = input.vertices[index];
                    let x = p[0] - axis.origin[0];
                    let y = p[1] - axis.origin[1];
                    (x * v[0] + y * v[1]).portable_atan2(x * u[0] + y * u[1])
                };
                let start = phase(a);
                let positive = (phase(b) - start).rem_euclid(std::f64::consts::TAU);
                let sweep = if edge.forward == direction {
                    positive
                } else {
                    positive - std::f64::consts::TAU
                };
                integrate(&|t| {
                    let angle = start + sweep * t;
                    (
                        std::array::from_fn(|k| {
                            axis.origin[k]
                                + radius_mm
                                    * (u[k] * angle.portable_cos() + v[k] * angle.portable_sin())
                        }),
                        std::array::from_fn(|k| {
                            radius_mm
                                * sweep
                                * (-u[k] * angle.portable_sin() + v[k] * angle.portable_cos())
                        }),
                    )
                });
            }
        }
    }
    let Surface::Plane(axis) = input.surfaces[input.faces[cap].surface] else {
        panic!("actual cap plane")
    };
    (
        sums[0].abs(),
        [sums[1] / sums[0], sums[2] / sums[0], axis.origin[2]],
    )
}

#[test]
fn actual_external_pipe_head_lazy_ready_descriptors_match_independent_cap_and_hex_integrals() {
    let source_spec = spec(true, true, false);
    let policy = GeometryPolicy::default();
    let owner = EntityId::new();
    let across = 18.;
    let height = 3.;
    let shape = Solid::rust_authored_axial_pipe_hardware_with_interrupt(
        source_spec,
        PipeHardwareKind::HexHead {
            across_flats_mm: across,
            height_mm: height,
        },
        1.,
        owner,
        frame(false),
        policy,
        None,
    )
    .unwrap();
    let body = shape.pipe_hardware_body().unwrap();
    let prepared = &body.source._prepared;
    assert_eq!(prepared.descriptor_queries.get(), 0);
    let (boundary, pose) = prepared.thread.rust_radial_boundary().unwrap();
    assert_eq!(pose, PlaneFrame::default());
    let original = boundary.input().clone();
    let full = boundary
        .face_properties_with_interrupt(source::MEASURE_ACCURACY, None)
        .unwrap();
    let contains = |range: I, value: f64| {
        assert!(
            range.lo <= value && value <= range.hi,
            "{range:?} misses {value}"
        )
    };
    let oracles = prepared.caps.map(|cap| {
        let coarse = cap_green_oracle(&original, cap as usize, 256);
        let fine = cap_green_oracle(&original, cap as usize, 512);
        assert!((coarse.0 - fine.0).abs() < 1e-8);
        for (a, b) in coarse.1.into_iter().zip(fine.1) {
            assert!((a - b).abs() < 1e-9);
        }
        let actual = prepared.measure(cap, None).unwrap();
        let old = full.iter().find(|f| f.face_index == cap as usize).unwrap();
        assert_eq!([actual.0.lo, actual.0.hi], old.area_mm2);
        assert_eq!(actual.1.map(|x| [x.lo, x.hi]), old.centroid_mm);
        contains(actual.0, fine.0);
        for (range, value) in actual.1.into_iter().zip(fine.1) {
            contains(range, value);
        }
        fine
    });
    assert_eq!(
        prepared.descriptor_queries.get(),
        0,
        "actual ready cap reads cannot force curved integration"
    );
    let faces = shape.faces().unwrap();
    assert_eq!(prepared.descriptor_queries.get(), 1);
    assert_eq!(shape.faces().unwrap(), faces);
    assert_eq!(prepared.descriptor_queries.get(), 1);
    let hex_area = 3_f64.sqrt() * across * across / 2.;
    let last = body.source.faces.len() - 1;
    let FaceMeasure::Ready { area, center } = body.source.faces[last].measure else {
        panic!("head top is actual planar stock")
    };
    contains(area, hex_area);
    for (range, value) in center
        .into_iter()
        .zip([0., 0., source_spec.length_mm + height])
    {
        contains(range, value);
    }
    assert!((faces[last].area_mm2 - hex_area).abs() < 1e-10);
    assert!((faces[last].center_mm[2] + height).abs() < 1e-12);
    let shoulder_area = hex_area - oracles[1].0;
    let expected = [
        -oracles[1].0 * oracles[1].1[0] / shoulder_area,
        -oracles[1].0 * oracles[1].1[1] / shoulder_area,
        source_spec.length_mm,
    ];
    let FaceMeasure::Ready { area, center } = body.source.faces[last - 1].measure else {
        panic!("actual cap subtraction")
    };
    contains(area, shoulder_area);
    for (range, value) in center.into_iter().zip(expected) {
        contains(range, value);
    }
    let bottom = body.source.face_map[&prepared.caps[0]] as usize;
    assert!((faces[bottom].area_mm2 - oracles[0].0).abs() < 1e-4);
    assert!((faces[bottom].center_mm[2] - source_spec.length_mm).abs() < 1e-12);
    let lateral_area = across / 3_f64.sqrt() * height;
    for (i, face) in faces[last - 7..last - 1].iter().enumerate() {
        assert!((face.area_mm2 - lateral_area).abs() < 1e-10);
        let angles = [i as f64, (i + 1) as f64].map(|i| i * std::f64::consts::PI / 3.);
        let radius = across / 3_f64.sqrt();
        let expected = [
            radius * (angles[0].portable_cos() + angles[1].portable_cos()) / 2.,
            -radius * (angles[0].portable_sin() + angles[1].portable_sin()) / 2.,
            -height / 2.,
        ];
        for (actual, expected) in face.center_mm.into_iter().zip(expected) {
            assert!((actual - expected).abs() < 1e-11);
        }
    }
    assert_eq!(boundary.input(), &original);
}

#[test]
fn actual_pipe_affine_centroid_shift_refuses_before_lazy_queries_and_valid_retry_keeps_source() {
    let source_spec = spec(true, false, false);
    let owner = EntityId::new();
    let policy = GeometryPolicy::default();
    let mut factory = PipeThreadFactory::new(32 * 1024 * 1024, 16);
    let kind = PipeHardwareKind::HexNut {
        across_flats_mm: 18.,
    };
    factory.begin_compile();
    let shape = factory
        .build(source_spec, kind, 1., owner, frame(false), policy, None)
        .unwrap();
    factory.finish_compile(true);
    let body = shape.pipe_hardware_body().unwrap();
    let prepared = &body.source._prepared;
    let (boundary, _) = prepared.thread.rust_radial_boundary().unwrap();
    let original = boundary.input().clone();
    let transform = |epsilon| RigidTransform {
        rotation: [[1. + epsilon, 0., 0.], [0., 1. - epsilon, 0.], [0., 0., 1.]],
        translation_mm: [0.; 3],
    };
    let rejected = transform(4.9e-11);
    rejected.validate(policy).unwrap(); // reach the actual provider, not the public rigid-input guard
    let p = body.placement.transformed(rejected, policy).unwrap();
    let span = body.source.bounds.iter().map(|b| b[1] - b[0]).sum::<f64>();
    let scale = p.area_scale();
    let delta = scale.hi / scale.lo - 1.;
    assert!(p.sigma[1] * span * delta / (1. - delta) > policy.linear_tolerance_mm / 8.);
    assert_eq!(
        shape.transformed(rejected).err().unwrap().kind,
        GeometryErrorKind::Unsupported
    );
    assert_eq!(prepared.descriptor_queries.get(), 0);
    assert!(prepared.lazy_measured.borrow().is_none());
    let valid = shape.transformed(transform(1e-13)).unwrap();
    assert!(Rc::ptr_eq(
        prepared,
        &valid.pipe_hardware_body().unwrap().source._prepared
    ));
    assert!((valid.volume_mm3().unwrap() / shape.volume_mm3().unwrap() - 1.).abs() < 1e-12);
    factory.begin_compile();
    let retry = factory
        .build(source_spec, kind, 1., owner, frame(true), policy, None)
        .unwrap();
    factory.finish_compile(true);
    assert_eq!(factory.stats().hits, 1);
    assert_eq!(prepared.descriptor_queries.get(), 0);
    assert!(Rc::ptr_eq(
        prepared,
        &retry.pipe_hardware_body().unwrap().source._prepared
    ));
    retry.bounds().unwrap();
    assert_eq!(boundary.input(), &original);
}
