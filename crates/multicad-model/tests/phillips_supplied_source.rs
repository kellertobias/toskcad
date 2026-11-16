//! Actual Core-fit/Model authors versus independent raw-cell and closed-source
//! contracts. These receipts do not enable the still-unqualified native export.
use multicad_core::{
    commands::{Actor, Channel, Session},
    document::Document,
    model::{
        self, MechanicalPartParameters, MechanicalPartRequest, PhillipsHParameters, PhillipsHSize,
        PhillipsInnerForm, ResolvedMechanicalPart, ScrewDrive, ScrewHead, ThreadFit, ThreadHand,
        ThreadRepresentation, angle, length,
    },
    project_file::{self, ProjectFile},
};
use multicad_kernel::{
    RigidTransform,
    rust_kernel::{step, step_read::*},
};
use multicad_model::{Compiler, GeometryBackend};
fn parameters(nominal: f64, pitch: f64, left: bool) -> MechanicalPartParameters {
    let mut p = model::default_machine_screw();
    let MechanicalPartParameters::MachineScrew {
        underhead_length,
        shank_diameter,
        head,
        drive,
        thread,
        ..
    } = &mut p
    else {
        unreachable!()
    };
    let small = nominal == 3.;
    *underhead_length = Some(length(if small { 9. } else { 13. }));
    *shank_diameter = length(nominal);
    *head = ScrewHead::Pan {
        profile: Box::new(model::CrownedScrewHeadParameters {
            diameter: length(if small { 5.6 } else { 12. }),
            height: length(if small { 2.4 } else { 4.6 }),
            crown_radius: length(if small { 5. } else { 10. }),
            top_land_diameter: length(if small { 4.1 } else { 8.4 }),
            neck_radius: length(if small { 0.1 } else { 0.25 }),
        }),
    };
    *drive = ScrewDrive::PhillipsH {
        profile: Box::new(PhillipsHParameters {
            size: if small {
                PhillipsHSize::H1
            } else {
                PhillipsHSize::H3
            },
            mouth_diameter: length(if small { 3.2 } else { 6.9 }),
            mouth_radius: length(if small { 0.5 } else { 0.8 }),
            b_width: length(if small { 0.955 } else { 2.395 }),
            g_diameter: length(if small { 1.27 } else { 3.81 }),
            tilt_angle: angle(if small { 7. } else { 5.75 }),
            outer_flank_angle: angle(26.5),
            floor_angle: angle(28.),
            outer_corner_angle: angle(92.),
            inner: PhillipsInnerForm::ShallowV {
                angle: angle(if small { 138. } else { 146. }),
                e_width: length(if small { 0.435 } else { 2.005 }),
            },
        }),
    };
    thread.major_diameter = length(nominal);
    thread.pitch = length(pitch);
    thread.length = length(if small { 8. } else { 12. });
    thread.fit = ThreadFit::Metric6gExternal;
    thread.metric_fit_engagement = Some(model::MetricFitEngagement::IntendedLength {
        length: length(nominal),
    });
    thread.hand = if left {
        ThreadHand::Left
    } else {
        ThreadHand::Right
    };
    thread.representation = ThreadRepresentation::Modeled;
    p
}
fn compiler() -> Compiler {
    Compiler {
        omit_display: true,
        ..Default::default()
    }
}
fn verify(nominal: f64, pitch: f64, left: bool) {
    let mut session =
        Session::new(Document::new("Actual metric supplied source").unwrap()).unwrap();
    let request = MechanicalPartRequest::new(
        format!("M{nominal} actual 6g"),
        parameters(nominal, pitch, left),
    );
    session
        .commit(
            model::prepare_mechanical_part(
                &session,
                session.revision(),
                Actor::local(Channel::Mcp),
                &request,
            )
            .unwrap()
            .prepared,
        )
        .unwrap();
    let ResolvedMechanicalPart::MachineScrew {
        thread: resolved,
        shank_radius_mm,
        underhead_length_mm,
        ..
    } = request.parameters.resolve(session.parameters()).unwrap()
    else {
        unreachable!()
    };
    let root = resolved
        .metric_fit
        .as_ref()
        .unwrap()
        .external_root_radius_mm
        .unwrap();
    assert_eq!(shank_radius_mm, nominal / 2.);
    assert!(resolved.major_radius_mm < shank_radius_mm && root > 0.);
    let compiled = Compiler::default()
        .compile_with_backend(
            session.document(),
            session.revision(),
            GeometryBackend::Rust,
            None,
        )
        .unwrap();
    let display = &compiled.snapshot().shapes[&request.feature].mesh;
    assert!(!display.indices.is_empty());
    assert!(display.vertices.iter().all(|v| v.is_finite()));
    let manufacturing = compiled
        .solid(request.feature)
        .unwrap()
        .manufacturing_mesh_with_receipt(Default::default())
        .unwrap();
    assert_eq!(
        manufacturing.edge_pairs as usize * 2,
        manufacturing.mesh.indices.len()
    );
    let pose = if left {
        RigidTransform {
            translation_mm: [105., -20., 60.],
            rotation: [[0.8, 0.36, -0.48], [0.6, -0.48, 0.64], [0., -0.8, -0.6]],
        }
    } else {
        Default::default()
    };
    let solid = compiled
        .solid(request.feature)
        .unwrap()
        .transformed(pose)
        .unwrap();
    let (source, frame) = solid
        .rust_phillips_screw_geometry()
        .expect("actual compiled fitted source, not nominal substitute");
    assert_eq!(source.thread().spec().root_radius_mm, Some(root));
    assert_eq!(
        source.thread().spec().major_radius_mm,
        resolved.major_radius_mm
    );
    assert_eq!(source.profile().end()[0], shank_radius_mm);
    assert_eq!(
        source.spec().shank_length_mm,
        underhead_length_mm - resolved.length_mm
    );
    assert_eq!(source.thread().spec().left_hand, left);
    let extent = source
        .profile()
        .support([1., 0., 0.])
        .max(resolved.length_mm)
        .max(1.);
    let properties_options = multicad_kernel::rust_kernel::phillips_drive::FacePropertyOptions {
        area_tolerance_mm2: 1e-7 * extent / 16.,
        first_moment_tolerance_mm3: 1e-7 * extent * extent / 16.,
        centroid_tolerance_mm: 1e-7 / 16.,
        maximum_cells: 65536,
    };
    let properties = source
        .face_properties_with_cancel(properties_options, &|| false)
        .unwrap();
    for face in properties {
        assert!(face.area_mm2[1] - face.area_mm2[0] <= properties_options.area_tolerance_mm2);
        for value in face.first_moment_mm3 {
            assert!(value[1] - value[0] <= properties_options.first_moment_tolerance_mm3);
        }
        for value in face.centroid_mm {
            assert!(value[1] - value[0] <= properties_options.centroid_tolerance_mm);
        }
    }
    let actual_frame = [frame.origin, frame.u, frame.v];
    let raw = step::write_phillips_screw(source, frame.origin, frame.u, frame.v, 1e-7)
        .expect("actual low-level source writer");
    let temp = tempfile::tempdir().unwrap();
    let folder = std::env::var_os("MULTICAD_SUPPLIED_MODEL_OUTPUT")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|| temp.path().to_path_buf())
        .join(format!(
            "m{}-{}",
            nominal as u32,
            if left { "lh" } else { "rh" }
        ));
    std::fs::create_dir_all(&folder).unwrap();
    std::fs::write(folder.join("actual-source.step"), &raw).unwrap();
    let source_known =
        verify_phillips_screw_source_known(&raw, source, frame.origin, frame.u, frame.v, 1e-7)
            .expect("supplied source-known carrier/cavity incidence contract");
    assert_eq!(source_known.material_cavity_face_count(), 28);
    let body = read_phillips_body(&raw, 1e-7, 1_000_000).expect("independent body reader");
    let outer = body
        .outer_carrier_bounds(
            PhillipsOuterCarrierOptions {
                maximum_distance_mm: 1e-7,
                work_limit: 10_000_000,
            },
            None,
        )
        .unwrap();
    let b = body.cavity();
    let endpoints = b
        .endpoint_coupling_bounds(
            PhillipsEndpointCouplingOptions {
                maximum_distance_mm: 1e-7,
                work_limit: 1_000_000,
            },
            None,
        )
        .unwrap();
    let analytic = b
        .analytic_support_bounds(
            &endpoints,
            PhillipsAnalyticSupportOptions {
                maximum_distance_mm: 1e-7,
                work_limit: 10_000_000,
            },
            None,
        )
        .unwrap();
    let lifts = analytic
        .finite_edge_lifts(endpoints, Default::default(), None)
        .unwrap();
    let cubics = b
        .cubic_support_bounds(
            PhillipsCubicSupportOptions {
                maximum_distance_mm: 1e-7,
                work_limit: 10_000_000,
            },
            None,
        )
        .unwrap();
    let sewn = lifts
        .sewn_edge_sources(&body, &cubics, Default::default(), None)
        .unwrap();
    let representative = body
        .derived_cavity_representative(&sewn, Default::default(), None)
        .unwrap();
    let closed = representative
        .closed_boundary_certificate(Default::default(), None)
        .unwrap();
    let common = closed
        .common_mouth_head_envelope(&outer, Default::default(), None)
        .unwrap();
    let head = common
        .closed_repaired_head(Default::default(), None)
        .unwrap();
    let inferred = body.inferred_thread_spec();
    let correspondence = head
        .supplied_thread_correspondence(source, actual_frame, Default::default(), None)
        .expect("actual compiled author/frame whole-cell correspondence");
    assert!(std::ptr::eq(correspondence.source(), source));
    assert_eq!(body.inferred_thread_spec(), inferred);
    let mut bounds = [0f64; 3];
    let mut cells = 0usize;
    for t in correspondence.transfers() {
        bounds[0] = bounds[0].max(t.distance_mm);
        bounds[1] = bounds[1].max(t.derivative_mm[0]);
        bounds[2] = bounds[2].max(t.derivative_mm[1]);
        cells += t.original_cells[0] * t.original_cells[1];
    }
    assert!(bounds.into_iter().all(|b| b <= 1e-7));
    let occupied = head
        .degree11_thread_carrier(Default::default(), None)
        .unwrap()
        .certify_occupied(Default::default(), None)
        .expect("independent original D11 occupancy");
    let material = head
        .closed_degree11_material(&occupied, Default::default(), None)
        .expect("same original source closed material");
    let volume = material
        .volume_bounds(Default::default(), None)
        .expect("independent original D11 signed flux mass");
    let matched = volume
        .match_supplied_author(source, actual_frame, Default::default(), None)
        .expect("complete original author to independently closed D11 material correspondence");
    assert!(std::ptr::eq(matched.author(), source));
    assert!(std::ptr::eq(matched.decoded_volume(), &volume));
    assert!(matched.maximum_world_distance_mm() <= 1e-7);
    assert!(matched.maximum_world_chart_derivative_mm() <= 1e-7);
    assert!(
        matched
            .independent_volume_difference_mm3()
            .into_iter()
            .all(|v| v.abs() <= 1e-5)
    );
    assert!(matched.total_work_used() <= 512_000_000);
    assert!(matched.accounted_peak_bytes() <= 64 * 1024 * 1024);
    assert_eq!(matched.supplied_thread().entry_rim_transfers().len(), 4);
    eprintln!(
        "MODEL_SUPPLIED_MATERIAL nominal={nominal} left={left} c0={} d1={} mass={:?} total_work={} fresh_work={} cavity_cells={} retained_bytes={} peak_bytes={} cavity_accuracy={:?}",
        matched.maximum_world_distance_mm(),
        matched.maximum_world_chart_derivative_mm(),
        matched.independent_volume_difference_mm3(),
        matched.total_work_used(),
        matched.new_work_used(),
        matched.cavity_cells_used(),
        matched.retained_budget_bytes(),
        matched.accounted_peak_bytes(),
        matched.cavity_accuracy_mm()
    );
    let author_receipt = source.source_volume_bounds().unwrap();
    let author_volume = author_receipt.volume_mm3();
    let decoded_volume = volume.world_volume_mm3();
    let determinant = correspondence.frame_determinant_bounds();
    let author_volume = [
        (author_volume[0] * determinant[0]).next_down(),
        (author_volume[1] * determinant[1]).next_up(),
    ];
    // These independently certified intervals describe DISTINCT sources.
    // Report their earned scalar difference; whole-material correspondence
    // (head/cavity/annuli plus sewing) is still a separate unearned contract.
    let difference = [
        (author_volume[0] - decoded_volume[1]).next_down(),
        (author_volume[1] - decoded_volume[0]).next_up(),
    ];
    assert!(difference.into_iter().all(|v| v.is_finite()));
    let original_head = source
        .profile()
        .closed_repaired_meridian(Default::default(), None)
        .unwrap();
    let original_head_mass = original_head
        .volume_bounds(
            multicad_kernel::rust_kernel::head_profile::RepairedHeadVolumeOptions {
                maximum_width_mm3: 1e-9,
                ..Default::default()
            },
            None,
        )
        .unwrap()
        .volume_mm3();
    eprintln!(
        "MODEL_COMPONENTS nominal={nominal} left={left} source_head_repaired={original_head_mass:?} source_head_nominal_debug={} source_thread_nominal_debug={} source_cavity={:?} decoded_components={:?} source_precision={} source_phase_error={} source_cavity_width={} world_difference={difference:?} whole_material_match=EARNED_UNDER_POLICY source_spec={:?} decoded_spec={:?}",
        source.profile().volume_mm3(),
        source.thread().volume_mm3(),
        source.cavity().volume().volume_mm3,
        volume.component_volumes_mm3(),
        author_receipt.profile_precision_error_mm3(),
        author_receipt.thread_phase_domain_error_mm3(),
        author_receipt.cavity_source_width_mm3(),
        source.cavity().source().spec(),
        material.cavity().spec()
    );
    let owned = volume.owned_snapshot(Default::default(), None).unwrap();
    assert_eq!(owned.original_input(), raw);
    eprintln!(
        "MODEL_SUPPLIED nominal={nominal} left={left} actual_thread={:?} nominal_shank={shank_radius_mm} bytes={} cells={cells} c0_du_dv={bounds:?} source_work={} author_mass={author_volume:?} decoded_mass={decoded_volume:?} material_work={} owned_bytes={}",
        source.thread().spec(),
        raw.len(),
        correspondence.total_work_used(),
        volume.work_used(),
        owned.retained_budget_bytes()
    );
    // Explicit actual-owner Native handoff. Cheap plane selection must not
    // request the unrelated complete curved-surface measurement inventory.
    let owned = std::rc::Rc::new(owned);
    let native = multicad_kernel::Solid::rust_degree11_repaired_phillips_material(
        owned.clone(),
        request.feature,
        Default::default(),
        Default::default(),
        None,
    )
    .unwrap();
    let (retained, occurrence) = native.rust_degree11_repaired_phillips_source().unwrap();
    assert!(std::ptr::eq(retained, owned.as_ref()));
    assert_eq!(occurrence, multicad_kernel::PlaneFrame::default());
    let mut planes = 0;
    for index in 0..owned.face_count() {
        if native.face_plane(index as u32).unwrap().is_some() {
            planes += 1;
        }
    }
    let topology = native.topology_snapshot(request.feature, 0).unwrap();
    use multicad_core::references::TopologyKind;
    let faces = topology
        .elements
        .iter()
        .filter(|e| e.key.kind == TopologyKind::Face)
        .collect::<Vec<_>>();
    let edges = topology
        .elements
        .iter()
        .filter(|e| e.key.kind == TopologyKind::Edge)
        .collect::<Vec<_>>();
    let vertices = topology
        .elements
        .iter()
        .filter(|e| e.key.kind == TopologyKind::Vertex)
        .collect::<Vec<_>>();
    let n = owned.thread_faces().len();
    let h = owned.head().view().segment_count();
    let cavity = owned.cavity().boundary_topology().unwrap();
    assert_eq!(faces.len(), owned.face_count());
    assert_eq!(edges.len(), 3 * n + h + 1 + cavity.edges().len());
    assert_eq!(vertices.len(), 2 * n + cavity.vertices().len());
    for edge in &edges {
        assert_eq!(
            faces
                .iter()
                .filter(|face| face.origins.iter().any(|o| edge.boundary.contains(o)))
                .count(),
            2
        );
    }
    for index in 0..n {
        let role = format!("thread/source/{index}");
        let face = faces
            .iter()
            .find(|f| f.origins.iter().any(|o| o.role == role))
            .unwrap();
        assert_eq!(face.boundary.len(), 4);
        let entry = format!("thread/cap/0/rim/{index}");
        let edge = edges
            .iter()
            .find(|e| e.origins.iter().any(|o| o.role == entry))
            .unwrap();
        assert!(edge.boundary.iter().any(|o| o.role == "screw.thread/start"));
    }
    assert!(edges.iter().any(|e| {
        e.origins
            .iter()
            .any(|o| o.role == "screw.head/thread.junction")
    }));
    assert!(planes >= 3);
    let far = [1e6, 1e6, 1e6];
    assert!(!native.rust_thread_contains(far).unwrap());
    let query = native
        .rust_degree11_repaired_phillips_classify_world_box(
            far.map(|v| [v, v]),
            Default::default(),
            None,
        )
        .unwrap();
    assert!(std::ptr::eq(query.source(), owned.as_ref()));
    assert_eq!(
        query.location(),
        multicad_kernel::rust_kernel::phillips_drive::Classification::Outside
    );
    assert!(query.accounted_peak_bytes() <= 64 * 1024 * 1024);

    assert!(native.face_plane(owned.face_count() as u32).is_err());
    let curved = (0..owned.face_count())
        .find(|&i| {
            matches!(
                owned.face_source(i),
                Some(PhillipsDegree11MaterialFace::ThreadSide { .. })
            )
        })
        .unwrap();
    assert!(native.face_plane(curved as u32).unwrap().is_none());
    let tip = (0..owned.face_count())
        .find(|&i| {
            matches!(
                owned.face_source(i),
                Some(PhillipsDegree11MaterialFace::ThreadTip { .. })
            )
        })
        .unwrap();
    let selected = native.face(tip as u32).unwrap();
    assert_eq!(selected.index, tip as u32);
    let repeated = native.face(tip as u32).unwrap();
    assert_eq!(repeated.center_mm, selected.center_mm);
    assert_eq!(repeated.area_mm2.to_bits(), selected.area_mm2.to_bits());
    let translated = native
        .transformed(RigidTransform {
            translation_mm: [3., -2., 5.],
            ..Default::default()
        })
        .unwrap();
    let moved = translated.face(tip as u32).unwrap();
    for k in 0..3 {
        assert!((moved.center_mm[k] - selected.center_mm[k] - [3., -2., 5.][k]).abs() <= 1e-7);
    }
    assert_eq!(moved.area_mm2.to_bits(), selected.area_mm2.to_bits());

    assert!(selected.area_mm2 > 0. && selected.center_mm.into_iter().all(f64::is_finite));
    assert_eq!(
        native
            .rust_degree11_repaired_phillips_source()
            .unwrap()
            .0
            .original_input(),
        raw
    );
    // Saved/cold canonical reconstruction must retain the actual fit and source.
    let mut bytes = Vec::new();
    project_file::write(
        &ProjectFile::new(session.document().clone()).unwrap(),
        &mut bytes,
        Default::default(),
    )
    .unwrap();
    std::fs::write(folder.join("canonical-project.toscut"), &bytes).unwrap();
    let reopened = project_file::read(
        &mut bytes.as_slice(),
        Default::default(),
        &project_file::supported_capabilities(),
    )
    .unwrap()
    .project;
    let cold = compiler()
        .compile_with_backend(
            &reopened.document,
            session.revision(),
            GeometryBackend::Rust,
            None,
        )
        .unwrap();
    let restored = cold
        .solid(request.feature)
        .unwrap()
        .transformed(pose)
        .unwrap();
    let (cold_source, cold_frame) = restored.rust_phillips_screw_geometry().unwrap();
    assert_eq!(cold_source.thread().spec(), source.thread().spec());
    assert_eq!(cold_source.profile(), source.profile());
    assert_eq!(cold_source.spec(), source.spec());
    assert_eq!(
        step::write_phillips_screw(
            cold_source,
            cold_frame.origin,
            cold_frame.u,
            cold_frame.v,
            1e-7
        )
        .unwrap(),
        raw
    );
    // This bridge must not silently broaden the normal Native STEP gate.
    assert!(!solid.supports_step_export());
    let sentinel = folder.join("normal-export-unqualified.step");
    std::fs::write(&sentinel, b"existing accepted file").unwrap();
    assert!(solid.write_step(&sentinel).is_err());
    assert_eq!(std::fs::read(sentinel).unwrap(), b"existing accepted file");
}
#[test]
fn core_m3_right_actual_supplied_source_bridge() {
    verify(3., 0.5, false)
}
#[test]
fn core_m3_left_posed_actual_supplied_source_bridge() {
    verify(3., 0.5, true)
}
#[test]
fn core_m6_right_actual_supplied_source_bridge() {
    verify(6., 1., false)
}
#[test]
fn core_m6_left_posed_actual_supplied_source_bridge() {
    verify(6., 1., true)
}

#[test]
fn actual_core_m3_owned_native_coarse_candidate_earns_fresh_world_embedding() {
    let mut session = Session::new(Document::new("Native actual M3 mesh").unwrap()).unwrap();
    let request = MechanicalPartRequest::new("Actual M3", parameters(3., 0.5, false));
    let prepared = model::prepare_mechanical_part(
        &session,
        session.revision(),
        Actor::local(Channel::Mcp),
        &request,
    )
    .unwrap()
    .prepared;
    session.commit(prepared).unwrap();
    let compiled = compiler()
        .compile_with_backend(
            session.document(),
            session.revision(),
            GeometryBackend::Rust,
            None,
        )
        .unwrap();
    let (author, frame) = compiled
        .solid(request.feature)
        .unwrap()
        .rust_phillips_screw_geometry()
        .unwrap();
    let text = step::write_phillips_screw(author, frame.origin, frame.u, frame.v, 1e-7).unwrap();
    let owned = std::rc::Rc::new(
        read_degree11_phillips_material(&text, Default::default(), None)
            .unwrap()
            .into_material(),
    );
    let native = multicad_kernel::Solid::rust_degree11_repaired_phillips_material(
        owned.clone(),
        request.feature,
        Default::default(),
        Default::default(),
        None,
    )
    .unwrap();
    let options = multicad_kernel::MeshOptions {
        // A nonidentity actual source frame receives half the caller's
        // distance and99% of its angle before Native transport validation.
        // This explicitly coarse integration fixture leaves the qualified
        // 0.5mm/0.7rad source request intact; it is not a default-quality gate.
        linear_deflection_mm: 1.1,
        angular_deflection_rad: 0.72,
        maximum_triangles: 1_000_000,
    };
    let (mesh, receipt) = native
        .rust_degree11_repaired_phillips_mesh_with_receipt(options, Default::default(), None)
        .unwrap();
    assert!(std::ptr::eq(receipt.source(), owned.as_ref()));
    assert!(!mesh.indices.is_empty());
    assert!(receipt.literal_signed_volume_mm3() > 0.);
    assert!(receipt.embedding_pairs() > 0);
    assert!(receipt.query_work_used() <= 512_000_000);
    assert!(receipt.consumer_accounted_peak_bytes() <= 64 * 1024 * 1024);
    assert!(receipt.local_source_surface_error_mm() <= options.linear_deflection_mm);
    assert!(!receipt.reused_display_receipt());
}
