use super::*;
use multicad_rust_kernel::{
    head_profile::{Head, Profile},
    helical_thread::{Spec as ThreadSpec, Thread},
    phillips_drive::{Cavity, Cutter, Domain, Inner, PhillipsScrew, Spec},
    profiled_screw::SourceSpec,
    step, step_read,
};
pub(super) fn source(scale: f64, left: bool, posed: bool) -> Rc<OwnedDerivedMaterial> {
    // Real physical STEP producer data. These values only author the fixture;
    // the import chain earns all support/net/trim fits from the actual text.
    let spec = Spec {
        mouth_diameter_mm: 1.9 * scale,
        mouth_radius_mm: 0.3 * scale,
        b_width_mm: 0.61 * scale,
        g_diameter_mm: 0.81 * scale,
        tilt_angle_deg: 7.,
        outer_flank_angle_deg: 26.5,
        floor_angle_deg: 28.,
        outer_corner_angle_deg: 92.,
        inner: Inner::Rounded {
            radius_mm: 0.25 * scale,
        },
    };
    let tip = (spec.mouth_diameter_mm - spec.g_diameter_mm) / (2. * 26.5_f64.to_radians().tan())
        + spec.g_diameter_mm / 2. * 28_f64.to_radians().tan();
    let cutter = Cutter::new(
        spec,
        Domain {
            corner_start_depth_mm: -2. * spec.mouth_diameter_mm,
            corner_end_depth_mm: tip + spec.mouth_diameter_mm,
            corner_reach_mm: 100. * scale,
        },
    )
    .unwrap();
    let cavity = Cavity::new(cutter, 1e-5, 100_000).unwrap();
    let head = Head::Crowned {
        radius_mm: 2. * scale,
        height_mm: 1.6 * scale,
        crown_radius_mm: 3.2 * scale,
        top_land_radius_mm: 1.25 * scale,
        neck_radius_mm: 0.1 * scale,
    };
    let thread = Thread::new(ThreadSpec {
        major_radius_mm: 0.97 * scale,
        pitch_mm: 0.4 * scale,
        length_mm: 2. * scale,
        depth_mm: 0.24 * scale,
        crest_width_mm: 0.05 * scale,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts: 1,
        left_hand: left,
        root_radius_mm: None,
    })
    .unwrap();
    let screw = PhillipsScrew::new(
        thread,
        Profile::new_for_thread_union(head, scale, 0.).unwrap(),
        SourceSpec {
            head: Some(head),
            shank_length_mm: 0.,
            minimum_wall_mm: 1e-5,
        },
        cavity,
        100_000,
        &|| false,
    )
    .unwrap();
    let (origin, u, v) = if posed {
        ([105., -20., 60.], [0.8, 0.6, 0.], [0.36, -0.48, -0.8])
    } else {
        ([0.; 3], [1., 0., 0.], [0., 1., 0.])
    };
    let text = step::write_phillips_screw(&screw, origin, u, v, 1e-7).unwrap();
    let body = step_read::read_phillips_body(&text, 1e-7, 1_000_000).unwrap();
    let outer = body
        .outer_carrier_bounds(
            step_read::PhillipsOuterCarrierOptions {
                maximum_distance_mm: 1e-7,
                work_limit: 10_000_000,
            },
            None,
        )
        .unwrap();
    let cavity = body.cavity();
    let ends = cavity
        .endpoint_coupling_bounds(
            step_read::PhillipsEndpointCouplingOptions {
                maximum_distance_mm: 1e-7,
                work_limit: 1_000_000,
            },
            None,
        )
        .unwrap();
    let analytic = cavity
        .analytic_support_bounds(
            &ends,
            step_read::PhillipsAnalyticSupportOptions {
                maximum_distance_mm: 1e-7,
                work_limit: 10_000_000,
            },
            None,
        )
        .unwrap();
    let lifts = analytic
        .finite_edge_lifts(ends, Default::default(), None)
        .unwrap();
    let cubics = cavity
        .cubic_support_bounds(
            step_read::PhillipsCubicSupportOptions {
                maximum_distance_mm: 1e-7,
                work_limit: 10_000_000,
            },
            None,
        )
        .unwrap();
    let sewn = lifts
        .sewn_edge_sources(&body, &cubics, Default::default(), None)
        .unwrap();
    let derived = body
        .derived_cavity_representative(&sewn, Default::default(), None)
        .unwrap();
    let closed = derived
        .closed_boundary_certificate(Default::default(), None)
        .unwrap();
    let common = closed
        .common_mouth_head_envelope(&outer, Default::default(), None)
        .unwrap();
    let head = common
        .closed_repaired_head(Default::default(), None)
        .unwrap();
    let thread = head
        .derived_thread_carrier(Default::default(), None)
        .unwrap()
        .certify_occupied(1e-7, 1e-5, None)
        .unwrap();
    let material = head
        .closed_repaired_material(&thread, Default::default(), None)
        .unwrap();
    let owned = material
        .volume_bounds(Default::default(), None)
        .unwrap()
        .into_owned(Default::default(), None)
        .unwrap();
    assert_eq!(owned.original_input(), text);
    Rc::new(owned)
}
fn feature() -> EntityId {
    EntityId::from_uuid_bits(0x00ca_d352_0000_0000_0000_0000_0000_0001)
}
#[test]
fn actual_repaired_native_default_proves_source_pose_and_remaining_embedding() {
    for (scale, left, raw_pose) in [(1., false, false), (1.137, true, true)] {
        let source = source(scale, left, raw_pose);
        let solid = Solid::rust_repaired_phillips_material(
            Rc::clone(&source),
            feature(),
            PlaneFrame::default(),
            GeometryPolicy::default(),
            None,
        )
        .unwrap();
        assert!(solid.uses_rust_geometry());
        assert!(std::ptr::eq(
            solid.rust_repaired_phillips_source().unwrap().0,
            source.as_ref()
        ));
        let bounds = solid.bounds().unwrap();
        // Independent semantic station witness. At a source cap centre the
        // common chart uses z-topZ; volume/self-consistent bounds alone cannot
        // detect accidentally applying the mouth columns to unshifted z.
        let chart = source.chart();
        let top = source.head().view().station(0).unwrap();
        assert_eq!(top[1], chart.outer_top_z);
        let placement = solid
            .repaired_phillips_body()
            .unwrap()
            .0
            .placement(PlaneFrame::default())
            .unwrap();
        let (mouth_rim, measured) = placement.point([top[0], 0., top[1]]).unwrap();
        for (i, point) in mouth_rim.iter().enumerate() {
            let expected =
                I::point(chart.columns[0][i]) + I::point(chart.columns[1][i]) * I::point(top[0]);
            assert!((I::point(*point) - expected).magnitude() <= measured);
        }

        let tip = source.thread().cap_axes()[1].origin;
        let world_tip: [f64; 3] = std::array::from_fn(|i| {
            chart.columns[0][i]
                + chart.columns[1][i] * tip[0]
                + chart.columns[2][i] * tip[1]
                + chart.columns[3][i] * (tip[2] - chart.outer_top_z)
        });
        for (i, point) in world_tip.iter().enumerate() {
            assert!(bounds.min[i] <= *point && *point <= bounds.max[i]);
        }

        let volume = solid.rust_repaired_phillips_volume_bounds().unwrap();
        assert!(
            volume[0] <= solid.volume_mm3().unwrap() && solid.volume_mm3().unwrap() <= volume[1]
        );
        let start = std::time::Instant::now();
        let (mesh, receipt) = solid
            .rust_repaired_phillips_mesh_with_receipt(Default::default(), Default::default(), None)
            .unwrap();
        eprintln!(
            "REPAIRED_NATIVE scale{scale} LH{left} rawPose{raw_pose} {}facets work{} remaining{} peak{} time{:?}",
            mesh.indices.len() / 3,
            receipt.query_work_used(),
            receipt.remaining_work(),
            receipt.consumer_accounted_peak_bytes(),
            start.elapsed()
        );
        eprintln!(
            "REPAIRED_NATIVE phases[source,placement,manifold,embedding+publication] {:?}",
            receipt.timings()
        );
        assert_eq!(
            receipt.proof_work_used() + receipt.proof_remaining_work(),
            WORK
        );
        assert_eq!(receipt.proof_work_used() + receipt.remaining_work(), WORK);
        assert!(receipt.source_work_used() < receipt.proof_work_used());
        assert!(receipt.embedding_pairs() > 0);
        assert!(!receipt.reused_display());
        assert_eq!(receipt.new_source_volume_mm3(), volume);
        assert!(std::ptr::eq(receipt.source(), source.as_ref()));
        assert_eq!(receipt.quality(), MeshOptions::default());
        assert_eq!(receipt.occurrence(), PlaneFrame::default());
        assert!(
            receipt.point_transfer_bound_mm() <= MeshOptions::default().linear_deflection_mm / 4.
        );
        assert!(receipt.consumer_accounted_peak_bytes() <= BYTES);
        assert_eq!(
            receipt.validator_workspace_upper_bound_bytes(),
            workspace::validators(mesh.vertices.len() / 3, mesh.indices.len() / 3).unwrap()
        );
        assert_eq!(
            receipt.combined_live_upper_bound_bytes(),
            receipt.producer_combined_upper_bound_bytes()
                + receipt.consumer_accounted_peak_bytes()
                + receipt.validator_workspace_upper_bound_bytes()
        );
        assert!(
            receipt.combined_live_upper_bound_bytes()
                >= receipt.producer_combined_upper_bound_bytes()
        );
        assert!(
            receipt.combined_live_upper_bound_bytes() >= receipt.consumer_accounted_peak_bytes()
        );
        let tip_face = (0..source.face_count())
            .find(|&i| {
                source.face_source(i) == Some(step_read::PhillipsMaterialFaceSource::ThreadTip)
            })
            .unwrap() as u32;
        let top_face = (0..source.face_count())
            .find(|&i| {
                source.face_source(i) == Some(step_read::PhillipsMaterialFaceSource::TopAnnulus)
            })
            .unwrap() as u32;
        let mouth_normal = cross(chart.columns[1], chart.columns[2]);
        let mouth_norm = dot(mouth_normal, mouth_normal).sqrt();
        let mut nearest_tip = f64::INFINITY;
        let mut top_triangles = 0;
        for (ids, &face) in mesh
            .indices
            .as_chunks::<3>()
            .0
            .iter()
            .zip(&mesh.face_indices)
        {
            for &id in ids {
                let point: Point = std::array::from_fn(|i| mesh.vertices[3 * id as usize + i]);
                if face == tip_face {
                    let delta: Point = std::array::from_fn(|i| point[i] - world_tip[i]);
                    nearest_tip = nearest_tip.min(dot(delta, delta).sqrt());
                }
                if face == top_face {
                    let delta = std::array::from_fn(|i| point[i] - chart.columns[0][i]);
                    assert!(dot(mouth_normal, delta).abs() / mouth_norm <= 0.1);
                }
            }
            if face == top_face {
                top_triangles += 1;
            }
        }
        assert!(nearest_tip <= 0.1 && top_triangles > 0);

        assert!(receipt.literal_signed_volume_mm3() > 0.);
        let inward = multicad_rust_kernel::display::DisplayMesh {
            vertices: mesh.vertices.clone(),
            indices: mesh
                .indices
                .as_chunks::<3>()
                .0
                .iter()
                .flat_map(|t| [t[0], t[2], t[1]])
                .collect(),
            face_indices: mesh.face_indices.clone(),
        };
        assert_eq!(
            mesh::placed_connectivity(&inward).unwrap_err().kind,
            GeometryErrorKind::InvalidResult
        );
        for p in mesh.vertices.as_chunks::<3>().0 {
            for (i, point) in p.iter().enumerate() {
                assert!(bounds.min[i] - 0.1 <= *point && *point <= bounds.max[i] + 0.1);
            }
        }
    }
}

#[test]
fn repeated_native_display_reuses_proof_but_manufacturing_remains_fresh() {
    let source = source(1., false, false);
    let solid = Solid::rust_repaired_phillips_material(
        Rc::clone(&source),
        feature(),
        PlaneFrame::default(),
        Default::default(),
        None,
    )
    .unwrap();
    let (body, _) = solid.repaired_phillips_body().unwrap();
    let start = std::time::Instant::now();
    let (first, cold) = solid
        .rust_repaired_phillips_display_mesh_with_receipt(
            Default::default(),
            Default::default(),
            None,
        )
        .unwrap();
    let cold_time = start.elapsed();
    assert!(!cold.reused_display());
    assert_eq!(body.builds.get(), 1);
    let start = std::time::Instant::now();
    let (warm, reused) = solid
        .rust_repaired_phillips_display_mesh_with_receipt(
            Default::default(),
            Default::default(),
            None,
        )
        .unwrap();
    let warm_time = start.elapsed();
    assert_eq!(first, warm);
    assert!(reused.reused_display());
    assert_eq!(body.builds.get(), 1);
    assert_eq!(reused.proof_work_used(), cold.proof_work_used());
    assert_eq!(reused.proof_remaining_work(), cold.proof_remaining_work());
    assert!(reused.query_work_used() < reused.proof_work_used());
    assert_eq!(reused.query_work_used() + reused.remaining_work(), WORK);
    assert_eq!(solid.mesh(Default::default()).unwrap(), first);
    assert_eq!(body.builds.get(), 1);
    eprintln!(
        "REPAIRED_DISPLAY cold{cold_time:?} warm{warm_time:?} phases{:?}",
        cold.timings()
    );
    let polls = std::cell::Cell::new(0usize);
    let observe = || {
        polls.set(polls.get() + 1);
        false
    };
    solid
        .rust_repaired_phillips_display_mesh_with_receipt(
            Default::default(),
            Default::default(),
            Some(&observe),
        )
        .unwrap();
    for stop in 0..polls.get() {
        let n = std::cell::Cell::new(0);
        let cancel = || {
            let k = n.get();
            n.set(k + 1);
            k == stop
        };
        let error = solid
            .rust_repaired_phillips_display_mesh_with_receipt(
                Default::default(),
                Default::default(),
                Some(&cancel),
            )
            .err()
            .unwrap();
        assert_eq!(error.kind, GeometryErrorKind::Cancelled);
        assert_eq!(body.builds.get(), 1);
    }
    for options in [
        MeshOptions {
            linear_deflection_mm: f64::NAN,
            ..Default::default()
        },
        MeshOptions {
            angular_deflection_rad: 0.,
            ..Default::default()
        },
        MeshOptions {
            maximum_triangles: 0,
            ..Default::default()
        },
    ] {
        assert_eq!(
            solid
                .rust_repaired_phillips_display_mesh_with_receipt(options, Default::default(), None)
                .err()
                .unwrap()
                .kind,
            GeometryErrorKind::Input
        );
    }
    for budget in [
        RepairedPhillipsMeshBudget {
            maximum_work: 0,
            ..Default::default()
        },
        RepairedPhillipsMeshBudget {
            maximum_bytes: 0,
            ..Default::default()
        },
        RepairedPhillipsMeshBudget {
            maximum_work: WORK + 1,
            ..Default::default()
        },
        RepairedPhillipsMeshBudget {
            maximum_bytes: BYTES + 1,
            ..Default::default()
        },
    ] {
        assert!(
            solid
                .rust_repaired_phillips_display_mesh_with_receipt(Default::default(), budget, None)
                .is_err()
        );
        assert_eq!(body.builds.get(), 1);
    }
    assert_eq!(solid.mesh(Default::default()).unwrap(), first);
    let manufactured = solid
        .manufacturing_mesh_with_receipt(Default::default())
        .unwrap();
    assert_eq!(manufactured.mesh, first);
    assert_eq!(body.builds.get(), 2);
    assert!(manufactured.vertex_links > 0 && manufactured.edge_pairs > 0);
    // Manufacturing did not replace the still-valid display receipt.
    assert_eq!(solid.mesh(Default::default()).unwrap(), first);
    assert_eq!(body.builds.get(), 2);
    // A genuine source-provider refusal/cancellation after cache miss must
    // preserve the old earned entry, not only reject bad lookup arguments.
    let tight = MeshOptions {
        angular_deflection_rad: 0.01,
        ..Default::default()
    };
    assert_eq!(
        solid
            .rust_repaired_phillips_display_mesh_with_receipt(tight, Default::default(), None)
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Unsupported
    );
    let limited = RepairedPhillipsMeshBudget {
        maximum_work: 1025,
        ..Default::default()
    };
    assert_eq!(
        solid
            .rust_repaired_phillips_display_mesh_with_receipt(Default::default(), limited, None)
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Limit
    );
    let n = std::cell::Cell::new(0usize);
    let cancel = || {
        let k = n.get();
        n.set(k + 1);
        k == 3
    };
    assert_eq!(
        solid
            .rust_repaired_phillips_mesh_with_receipt(
                Default::default(),
                Default::default(),
                Some(&cancel)
            )
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Cancelled
    );
    // Initial Native reservation fits exactly, but the delegated provider
    // still needs its own scratch. The old cached display remains resident.
    let frontier = body.resident_bytes() + mesh::SCRATCH;
    let lower = RepairedPhillipsMeshBudget {
        maximum_bytes: frontier,
        ..Default::default()
    };
    let frontier_polls = std::cell::Cell::new(0usize);
    let observe_frontier = || {
        frontier_polls.set(frontier_polls.get() + 1);
        false
    };
    assert_eq!(
        solid
            .rust_repaired_phillips_display_mesh_with_receipt(
                Default::default(),
                lower,
                Some(&observe_frontier)
            )
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Limit
    );
    for stop in 0..frontier_polls.get() {
        let count = std::cell::Cell::new(0usize);
        let cancel = || {
            let n = count.get();
            count.set(n + 1);
            n == stop
        };
        assert_eq!(
            solid
                .rust_repaired_phillips_display_mesh_with_receipt(
                    Default::default(),
                    lower,
                    Some(&cancel)
                )
                .err()
                .unwrap()
                .kind,
            GeometryErrorKind::Cancelled
        );
    }
    let expected_builds = body.builds.get();
    assert_eq!(solid.mesh(Default::default()).unwrap(), first);
    assert_eq!(body.builds.get(), expected_builds);
    let changed = MeshOptions {
        linear_deflection_mm: 0.2,
        ..Default::default()
    };
    let (_, receipt) = solid
        .rust_repaired_phillips_display_mesh_with_receipt(changed, Default::default(), None)
        .unwrap();
    assert!(!receipt.reused_display());
    assert_eq!(receipt.quality(), changed);
    assert_eq!(body.builds.get(), expected_builds + 1);
    let posed = solid
        .transformed(RigidTransform {
            rotation: [[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]],
            translation_mm: [10., 20., 30.],
        })
        .unwrap();
    let (mesh, receipt) = posed
        .rust_repaired_phillips_display_mesh_with_receipt(changed, Default::default(), None)
        .unwrap();
    assert!(!receipt.reused_display());
    assert_eq!(body.builds.get(), expected_builds + 2);
    assert_eq!(posed.mesh(changed).unwrap(), mesh);
    assert_eq!(body.builds.get(), expected_builds + 2);
    let calls = std::cell::Cell::new(0usize);
    let nested = || {
        let n = calls.get();
        calls.set(n + 1);
        if n == 4 {
            let (_, proof) = posed
                .rust_repaired_phillips_display_mesh_with_receipt(
                    MeshOptions {
                        linear_deflection_mm: 0.3,
                        ..changed
                    },
                    Default::default(),
                    None,
                )
                .unwrap();
            assert!(!proof.reused_display());
        }
        false
    };
    let (again, proof) = posed
        .rust_repaired_phillips_display_mesh_with_receipt(
            changed,
            Default::default(),
            Some(&nested),
        )
        .unwrap();
    assert_eq!(again, mesh);
    assert!(proof.reused_display());
    assert_eq!(body.builds.get(), expected_builds + 3); // nested build, no outer source reproof
    assert_eq!(posed.mesh(changed).unwrap(), mesh);
    assert_eq!(body.builds.get(), expected_builds + 3);
}

#[test]
fn actual_source_factory_failures_and_unearned_apis_preserve_ownership() {
    let source = source(1., false, false);
    let original = source.original_input().to_owned();
    let refs = Rc::strong_count(&source);
    let retained = source.retained_budget_bytes()
        + std::mem::size_of::<RepairedPhillipsBody>()
        + 8 * std::mem::size_of::<usize>();
    assert!(check_source_storage(&source, retained).is_ok());
    assert_eq!(
        check_source_storage(&source, retained - 1)
            .unwrap_err()
            .kind,
        GeometryErrorKind::Limit
    );

    let make = |feature, frame, policy, cancel: Option<&dyn Fn() -> bool>| {
        Solid::rust_repaired_phillips_material(Rc::clone(&source), feature, frame, policy, cancel)
    };
    assert_eq!(
        make(
            EntityId::nil(),
            Default::default(),
            Default::default(),
            None
        )
        .err()
        .unwrap()
        .kind,
        GeometryErrorKind::Input
    );
    assert_eq!(
        make(
            feature(),
            Default::default(),
            GeometryPolicy {
                linear_tolerance_mm: 0.,
                ..Default::default()
            },
            None
        )
        .err()
        .unwrap()
        .kind,
        GeometryErrorKind::Input
    );
    assert_eq!(
        make(
            feature(),
            Default::default(),
            GeometryPolicy {
                boolean_fuzzy_mm: 1e-8,
                ..Default::default()
            },
            None
        )
        .err()
        .unwrap()
        .kind,
        GeometryErrorKind::Unsupported
    );
    for frame in [
        PlaneFrame {
            origin: [f64::NAN, 0., 0.],
            ..Default::default()
        },
        PlaneFrame {
            u: [0.; 3],
            ..Default::default()
        },
        PlaneFrame {
            v: [1e-7, 1., 0.],
            ..Default::default()
        },
        PlaneFrame {
            origin: [1e8, 0., 0.],
            ..Default::default()
        },
    ] {
        assert!(make(feature(), frame, Default::default(), None).is_err());
        assert_eq!(Rc::strong_count(&source), refs);
    }
    for stop in 0..3 {
        let n = std::cell::Cell::new(0);
        let cancel = || {
            let k = n.get();
            n.set(k + 1);
            k == stop
        };
        assert_eq!(
            make(
                feature(),
                Default::default(),
                Default::default(),
                Some(&cancel)
            )
            .err()
            .unwrap()
            .kind,
            GeometryErrorKind::Cancelled
        );
        assert_eq!(Rc::strong_count(&source), refs);
    }
    let solid = make(
        feature(),
        Default::default(),
        Default::default(),
        Some(&|| false),
    )
    .unwrap();
    assert_eq!(source.original_input(), original);
    let shared = solid.share().unwrap();
    assert!(std::ptr::eq(
        shared.rust_repaired_phillips_source().unwrap().0,
        source.as_ref()
    ));
    assert_eq!(shared.bounds().unwrap(), solid.bounds().unwrap());
    assert_eq!(
        shared.rust_repaired_phillips_volume_bounds().unwrap(),
        solid.rust_repaired_phillips_volume_bounds().unwrap()
    );
    assert!(solid.has_certified_source_mesh());

    let constrained = RepairedPhillipsMeshBudget {
        maximum_work: 1025,
        ..Default::default()
    };
    assert_eq!(
        solid
            .rust_repaired_phillips_mesh_with_receipt(Default::default(), constrained, None)
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Limit
    );
    let n = std::cell::Cell::new(0usize);
    let cancel_source = || {
        let k = n.get();
        n.set(k + 1);
        k == 3
    };
    assert_eq!(
        solid
            .rust_repaired_phillips_mesh_with_receipt(
                Default::default(),
                Default::default(),
                Some(&cancel_source)
            )
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Cancelled
    );
    let tight_angle = MeshOptions {
        angular_deflection_rad: 0.01,
        ..Default::default()
    };
    assert_eq!(
        solid
            .rust_repaired_phillips_mesh_with_receipt(tight_angle, Default::default(), None)
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Unsupported
    );

    assert_eq!(
        solid.faces().unwrap_err().kind,
        GeometryErrorKind::Unsupported
    );
    assert_eq!(
        solid.topology_snapshot(feature(), 1).unwrap_err().kind,
        GeometryErrorKind::Unsupported
    );
    assert_eq!(
        solid.display_edges(0.1).unwrap_err().kind,
        GeometryErrorKind::Unsupported
    );
    assert_eq!(
        solid
            .manufacturing_mesh_with_receipt(MeshOptions {
                maximum_triangles: 1_000_001, // Valid generic budget; exceeds this provider's hard cap.
                ..Default::default()
            })
            .unwrap_err()
            .kind,
        GeometryErrorKind::Input
    );
    let clocked = solid
        .transformed(RigidTransform {
            rotation: [[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]],
            translation_mm: [0.; 3],
        })
        .unwrap();
    assert_eq!(
        clocked
            .rust_repaired_phillips_mesh_with_receipt(tight_angle, Default::default(), None)
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Unsupported
    );
    assert!(!solid.supports_step_export());
    let path = tempfile::NamedTempFile::new().unwrap();
    std::fs::write(path.path(), b"existing output").unwrap();
    assert_eq!(
        solid.write_step(path.path()).unwrap_err().kind,
        GeometryErrorKind::Unsupported
    );
    assert_eq!(std::fs::read(path.path()).unwrap(), b"existing output");
    assert!(solid.append_step(&mut step::Writer::new()).is_err());
    let cube = Solid::box_solid([1.; 3], Default::default()).unwrap();
    assert!(cube.rust_repaired_phillips_source().is_none());
    assert_eq!(
        cube.rust_repaired_phillips_volume_bounds()
            .unwrap_err()
            .kind,
        GeometryErrorKind::Unsupported
    );
    assert_eq!(
        cube.rust_repaired_phillips_mesh_with_receipt(Default::default(), Default::default(), None)
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Unsupported
    );
    assert!(cube.volume_mm3().unwrap() > 0.);
    assert!(cube.manufacturing_mesh(Default::default()).is_ok());
}

#[test]
fn full_common_chart_keeps_symbolic_station_and_independent_fourth_column() {
    let source = source(1., false, false);
    let chart = source.chart();
    let identity = Placement::new(chart, Default::default()).unwrap();
    assert!(identity.is_identity());
    let mut nearby = chart;
    nearby.columns[0][2] = chart.outer_top_z.next_up();
    assert!(
        !Placement::new(nearby, Default::default())
            .unwrap()
            .is_identity()
    );
    let p = [0.25, 0.5, 0.4];
    let (same, error) = identity.point(p).unwrap();
    assert!(same.iter().zip(p).all(|(a, b)| (*a - b).abs() <= error));
    let mut sheared = chart;
    let epsilon = 2_f64.powi(-30);
    sheared.columns[3][0] = epsilon;
    let mapped = Placement::new(sheared, Default::default()).unwrap();
    assert!(!mapped.is_identity());
    let (actual, measured) = mapped.point(p).unwrap();
    let exact_x = I::point(sheared.columns[0][0])
        + I::point(p[0])
        + I::point(epsilon) * (I::point(p[2]) - I::point(chart.outer_top_z));
    assert!((I::point(actual[0]) - exact_x).magnitude() <= measured);
    assert_ne!(actual[0], same[0]); // a cross(U,V) substitution would lose this.
    let triangle = [p, [0.5, 0.5, 0.9], [0.25, 0.75, 0.4]];
    let world = triangle.map(|p| mapped.point(p).unwrap().0);
    assert!(mapped.normal_transfer(triangle, world, 0.001));
    assert!(mapped.point_within(p, measured).is_ok());
    assert!(mapped.point_within(p, 0.).is_err());
    assert!(mapped.point_within(p, f64::NAN).is_err());
    assert!(mapped.require_normal(triangle, world, 0.001).is_ok());
    assert!(
        mapped
            .require_normal(triangle, [world[0], world[2], world[1]], 0.001)
            .is_err()
    );
    let mut contracting = chart;
    contracting.columns[2][1] = 1. - 1e-8;
    let transport = Placement::new(contracting, Default::default()).unwrap();
    let subnormal = f64::from_bits(1);
    assert!(
        transport
            .volume_interval([subnormal, 2. * subnormal])
            .is_err()
    );
    assert!(transport.volume_interval([1., 2.]).unwrap()[0] < 1.);
    assert!(!mapped.normal_transfer(triangle, [world[0], world[2], world[1]], 0.001));
    assert!(!mapped.normal_transfer(triangle, [world[0]; 3], 0.001));
    for mode in 0..4 {
        let mut invalid = chart;
        match mode {
            0 => invalid.outer_top_z = f64::NAN,
            1 => invalid.columns[0][0] = f64::INFINITY,
            2 => invalid.columns[3] = [0., 0., -1.],
            _ => invalid.columns[3][0] = 0.1,
        }
        assert!(Placement::new(invalid, Default::default()).is_err());
    }
    assert!(mapped.point([f64::INFINITY, 0., 0.]).is_err());
    let enclosure = source.local_enclosure_bounds();
    let bounds = Bounds {
        min: enclosure.map(|b| b[0]),
        max: enclosure.map(|b| b[1]),
    };
    let mut distant = chart;
    distant.columns[0] = [f64::MAX; 3];
    assert!(
        Placement::new(distant, Default::default())
            .unwrap()
            .constructed_bounds(bounds)
            .is_err()
    );
    let occurrence = PlaneFrame {
        origin: [f64::MAX; 3],
        ..Default::default()
    };
    assert!(
        Placement::new(chart, occurrence)
            .unwrap()
            .bounds(bounds)
            .is_err()
    );
}

#[test]
fn coarse_actual_source_cancellation_spans_transfer_and_publication_stages() {
    use multicad_rust_kernel::{display::Tessellation, step_read::PhillipsJoinedMeshOptions};
    let source = source(1., false, false);
    let options = MeshOptions {
        linear_deflection_mm: 1.,
        angular_deflection_rad: 1.,
        ..Default::default()
    };
    for frame in [
        PlaneFrame::default(),
        PlaneFrame {
            origin: [10., 20., 30.],
            u: [0., 1., 0.],
            v: [-1., 0., 0.],
        },
    ] {
        let solid = Solid::rust_repaired_phillips_material(
            Rc::clone(&source),
            feature(),
            frame,
            Default::default(),
            None,
        )
        .unwrap();
        let body = solid.repaired_phillips_body().unwrap().0;
        let local = mesh::local_quality(
            Tessellation {
                deflection: 1.,
                angle: 1.,
                triangle_limit: 1_000_000,
            },
            body.placement(frame).unwrap(),
        )
        .unwrap();
        let polls = std::cell::Cell::new(0usize);
        let observe = || {
            polls.set(polls.get() + 1);
            false
        };
        let candidate = source
            .joined_mesh_candidate(
                PhillipsJoinedMeshOptions {
                    tessellation: local,
                    maximum_work: WORK - 1025,
                    maximum_bytes: BYTES
                        - (body.resident_bytes() - source.retained_budget_bytes())
                        - mesh::SCRATCH,
                },
                Some(&observe),
            )
            .unwrap();
        let vertices = candidate.mesh().vertices.len() / 3;
        let triangles = candidate.mesh().indices.len() / 3;
        let source_polls = polls.get();
        drop(candidate);
        polls.set(0);
        let (actual, receipt) = solid
            .rust_repaired_phillips_mesh_with_receipt(options, Default::default(), Some(&observe))
            .unwrap();
        let total = polls.get();
        let used = receipt.query_work_used();
        assert_eq!(actual.indices.len() / 3, triangles);
        eprintln!(
            "COARSE_NATIVE vertices{vertices} triangles{triangles} sourcepolls{source_polls} totalpolls{total} work{used}"
        );
        drop(actual);
        // Each source provider poll is measured independently on the SAME
        // immutable source/quality. Three Native setup polls precede it.
        let mut stops = vec![
            source_polls + 3,
            source_polls + 4,
            source_polls + 5,
            total - 1,
            total - 2,
            total - 3,
        ];
        if frame != PlaneFrame::default() {
            stops.extend([
                source_polls + 6,
                source_polls + 6 + vertices,
                source_polls + 6 + vertices + triangles,
            ]);
        }
        stops.sort_unstable();
        stops.dedup();
        for stop in stops {
            polls.set(0);
            let cancel = || {
                let n = polls.get();
                polls.set(n + 1);
                n == stop
            };
            let failure = solid
                .rust_repaired_phillips_mesh_with_receipt(
                    options,
                    Default::default(),
                    Some(&cancel),
                )
                .err()
                .unwrap();
            assert_eq!(
                failure.kind,
                GeometryErrorKind::Cancelled,
                "stop {stop} actual {}",
                polls.get()
            );
        }
        if frame == PlaneFrame::default() {
            let (_, exact) = solid
                .rust_repaired_phillips_mesh_with_receipt(
                    options,
                    RepairedPhillipsMeshBudget {
                        maximum_work: used,
                        ..Default::default()
                    },
                    None,
                )
                .unwrap();
            assert_eq!(exact.remaining_work(), 0);
            assert_eq!(
                solid
                    .rust_repaired_phillips_mesh_with_receipt(
                        options,
                        RepairedPhillipsMeshBudget {
                            maximum_work: used - 1,
                            ..Default::default()
                        },
                        None
                    )
                    .err()
                    .unwrap()
                    .kind,
                GeometryErrorKind::Limit
            );
        }
    }
}
