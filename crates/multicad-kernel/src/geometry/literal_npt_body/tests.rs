use super::*;
fn stock(quarter: bool) -> multicad_rust_kernel::axial_pipe_thread::authored::StockClipPlan {
    use multicad_rust_kernel::axial_pipe_thread::authored::{
        AuthoredSection, AuthoredSpec, AuthoredThread,
    };
    let (tpi, e1, l1, l2, od, max) = if quarter {
        (18., 0.49163, 0.2278, 0.4018, 0.540, 0.088)
    } else {
        (27., 0.37360, 0.1615, 0.2639, 0.405, 0.096)
    };
    let pitch = 25.4 / tpi;
    let k = 1. / 32.;
    let m = 3f64.sqrt();
    let s = 1. / (m - k) + 1. / (m + k);
    let depth = pitch / s - (0.033 + max) * pitch;
    let spec = AuthoredSpec {
        major_radius_mm: (e1 * 25.4 - 2. * k * l1 * 25.4 + depth) / 2.,
        pitch_mm: pitch,
        length_mm: l2 * 25.4,
        radial_taper: k,
        section: AuthoredSection::NptFlat {
            depth_mm: depth,
            flat_width_mm: (pitch - depth * s) / 2.,
        },
        crest_radius_limit_mm: Some(od * 25.4 / 2.),
    };
    AuthoredThread::new(spec, 1e-7, 1e-6)
        .unwrap()
        .stock_clip_plan(Default::default(), None)
        .unwrap()
}
fn text(quarter: bool) -> String {
    use multicad_rust_kernel::step::npt_stock as step;
    let source = stock(quarter);
    let curves = step::CurvePlan::new(&source, Default::default(), None).unwrap();
    let surfaces = step::SurfacePlan::new(curves, Default::default(), None).unwrap();
    step::write_source_candidate(&surfaces, 16_000_000, 64 * 1024 * 1024, None)
        .unwrap()
        .into_string()
}
fn mapped(text: &str) -> String {
    let doc = multicad_rust_kernel::step_read::owned_insert::CompoundDocument::parse(
        text,
        Default::default(),
        None,
    )
    .unwrap();
    let solid = doc.represented_members()[0].item_entity();
    let old = text
        .lines()
        .find(|s| s.contains("=ADVANCED_BREP_SHAPE_REPRESENTATION("))
        .unwrap();
    let context = old.rsplit('#').next().unwrap().split(')').next().unwrap();
    let changed = old.replacen(&format!("#{solid}"), "#900020", 1);
    let extra = format!(
        "#900010=CARTESIAN_POINT('',(0.,0.,0.));\n#900011=AXIS2_PLACEMENT_3D('',#900010,$,$);\n#900012=SHAPE_REPRESENTATION('',(#{solid},#900011),#{context});\n#900013=REPRESENTATION_MAP(#900011,#900012);\n#900014=CARTESIAN_POINT('',(11.,-5.,3.));\n#900015=DIRECTION('',(0.,2.,0.));\n#900016=DIRECTION('',(-3.,0.,0.));\n#900017=DIRECTION('',(0.,-4.,3.));\n#900018=(CARTESIAN_TRANSFORMATION_OPERATOR(#900015,#900016,#900014,1.)CARTESIAN_TRANSFORMATION_OPERATOR_3D(#900017)FUNCTIONALLY_DEFINED_TRANSFORMATION('','')GEOMETRIC_REPRESENTATION_ITEM()REPRESENTATION_ITEM(''));\n#900020=MAPPED_ITEM('',#900013,#900018);\n"
    );
    text.replacen(old, &changed, 1).replacen(
        "ENDSEC;\nEND-ISO-10303-21;",
        &format!("{extra}ENDSEC;\nEND-ISO-10303-21;"),
        1,
    )
}
fn source(text: &str) -> read::OwnedMaterial {
    let doc = multicad_rust_kernel::step_read::owned_insert::CompoundDocument::parse(
        text,
        Default::default(),
        None,
    )
    .unwrap();
    read::Selected::parse(
        &doc,
        doc.represented_members()[0].item_entity(),
        Default::default(),
        None,
    )
    .unwrap()
    .retain(Default::default(), None)
    .unwrap()
}
#[test]
fn actual_literal_native_owns_source_topology_mass_queries_and_manufacturing() {
    for (quarter, placed) in [(false, false), (true, false), (false, true)] {
        let text = text(quarter);
        let text = if placed { mapped(&text) } else { text };
        let material = source(&text);
        let expected = material.world_volume_mm3();
        let inside = material
            .map_box_mm([[0., 0.], [0., 0.], [2., 2.]])
            .unwrap()
            .map(|v| v[0].midpoint(v[1]));
        let owner = EntityId::new();
        let solid = Solid::import_rust_step(&text, owner, Default::default()).unwrap();
        assert!(solid.rust_clipped_npt_source().is_none());
        assert!(solid.rust_authored_axial_pipe_thread_source().is_none());
        let actual = solid.volume_mm3().unwrap();
        assert!(actual >= expected[0] && actual <= expected[1]);
        let top = solid.topology_snapshot(owner, 0).unwrap();
        let combined = solid.topology_with_display(owner, 0, Some(0.05)).unwrap();
        assert_eq!(combined.topology, top);
        assert_eq!(combined.edges.len(), 20);
        assert_eq!(combined.vertices.len(), 12);
        assert_eq!(solid.display_edges(0.05).unwrap().len(), 20);

        for (kind, count) in [
            (TopologyKind::Face, 10),
            (TopologyKind::Edge, 20),
            (TopologyKind::Vertex, 12),
        ] {
            assert_eq!(
                top.elements.iter().filter(|e| e.key.kind == kind).count(),
                count
            );
        }
        let faces = solid.faces().unwrap();
        assert_eq!(faces.len(), 10);
        assert!(
            faces
                .iter()
                .all(|f| f.area_mm2 > 0. && f.area_mm2.is_finite())
        );
        let cold = solid
            .mesh(Default::default())
            .unwrap_or_else(|e| panic!("quarter={quarter} placed={placed}: {e:?}"));
        let evidence = solid.rust_literal_npt_runtime_evidence().unwrap();
        let hot = solid.mesh(Default::default()).unwrap();
        let repeated = solid.rust_literal_npt_runtime_evidence().unwrap();
        assert_eq!(cold.vertices, hot.vertices);
        assert_eq!(evidence.mesh_builds, repeated.mesh_builds);
        assert_eq!(evidence.last_mesh_work, repeated.last_mesh_work);
        assert!(repeated.last_mesh_actual_work.unwrap() < evidence.last_mesh_actual_work.unwrap());
        assert!(
            repeated.resident_budget_bytes < BYTES
                && repeated.last_mesh_peak_bytes.unwrap() <= BYTES
        );
        let mesh = solid
            .manufacturing_mesh_with_receipt(MeshOptions {
                linear_deflection_mm: 0.1,
                angular_deflection_rad: 0.3,
                maximum_triangles: 100_000,
            })
            .unwrap();
        assert!(mesh.edge_pairs > 0 && mesh.vertex_links > 0);
        assert!(mesh.mesh.face_indices.iter().all(|f| *f < 10));
        assert!(solid.rust_literal_npt_contains(inside).unwrap());
        assert!(!solid.rust_literal_npt_contains([100., 0., 2.]).unwrap());
        assert_eq!(
            solid.share().unwrap().topology_snapshot(owner, 0).unwrap(),
            top
        );
        let moved = solid
            .transformed(RigidTransform {
                translation_mm: [1., 0., 0.],
                ..Default::default()
            })
            .unwrap();
        assert_eq!(
            moved
                .rust_literal_npt_occurrence()
                .unwrap()
                .transforms()
                .len(),
            1
        );
        assert_eq!(
            moved
                .rust_literal_npt_source()
                .unwrap()
                .original_mapping_identity(),
            solid
                .rust_literal_npt_source()
                .unwrap()
                .original_mapping_identity()
        );
        assert!(
            moved
                .rust_literal_npt_contains([inside[0] + 1., inside[1], inside[2]])
                .unwrap()
        );
        assert!(
            moved
                .topology_with_display(owner, 0, Some(0.05))
                .unwrap()
                .edges
                .len()
                == 20
        );
    }
}
#[test]
fn actual_normal_step_import_keeps_mixed_literal_members_separate() {
    use multicad_rust_kernel::step::{Writer, npt_stock as step};
    let a = stock(false);
    let ac = step::CurvePlan::new(&a, Default::default(), None).unwrap();
    let ap = step::SurfacePlan::new(ac, Default::default(), None).unwrap();
    let b = stock(true);
    let bc = step::CurvePlan::new(&b, Default::default(), None).unwrap();
    let bp = step::SurfacePlan::new(bc, Default::default(), None).unwrap();
    let av = ap.verify_candidate(Default::default(), None).unwrap();
    let bv = bp.verify_candidate(Default::default(), None).unwrap();
    let mut writer = Writer::new();
    writer
        .append_verified_source_npt_candidate(&av, Default::default(), None)
        .unwrap();
    writer
        .append_box(
            [2, 3, 4].map(|v| Rational::new(v, 1).unwrap()),
            [-20., 0., 0.],
            [1., 0., 0.],
            [0., 1., 0.],
        )
        .unwrap();
    writer
        .append_verified_source_npt_candidate(&bv, Default::default(), None)
        .unwrap();
    let text = writer.finish().unwrap();
    let owner = EntityId::new();
    let solid = Solid::import_rust_step(&text, owner, Default::default()).unwrap();
    let compound = solid.rust_compound().unwrap();
    assert_eq!(compound.members().count(), 3);
    assert_eq!(
        compound
            .members()
            .filter(|m| m.rust_literal_npt_source().is_some())
            .count(),
        2
    );
    let display = solid.topology_with_display(owner, 0, Some(0.1)).unwrap();
    assert_eq!(display.edges.len(), 52);
    assert_eq!(display.vertices.len(), 32);
    let volume = solid.volume_mm3().unwrap();
    let expected = av.literal_material().world_volume_mm3()[0]
        + bv.literal_material().world_volume_mm3()[0]
        + 24.;
    assert!((volume - expected).abs() < 1e-4);
}
#[test]
fn actual_literal_native_preserves_validation_cancellation_and_atomic_refusal() {
    let material = source(&text(false));
    let signature = source_signature(&material).unwrap();
    assert_eq!(signature, source_signature(&material.clone()).unwrap());
    assert_eq!(
        Solid::from_literal_npt_source(material.clone(), EntityId::nil(), Default::default(), None)
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Input
    );
    assert_eq!(
        Solid::from_literal_npt_source(
            material.clone(),
            EntityId::new(),
            Default::default(),
            Some(&|| true)
        )
        .err()
        .unwrap()
        .kind,
        GeometryErrorKind::Cancelled
    );
    let solid = Solid::from_literal_npt_source(material, EntityId::new(), Default::default(), None)
        .unwrap();
    assert_eq!(
        solid
            .mesh(MeshOptions {
                maximum_triangles: 1,
                ..Default::default()
            })
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Limit
    );
    let evidence = solid.rust_literal_npt_runtime_evidence().unwrap();
    assert_eq!(evidence.mesh_builds, 0);
    assert!(evidence.last_mesh_work.is_none());
    assert!(
        solid
            .mesh(MeshOptions {
                angular_deflection_rad: f64::NAN,
                ..Default::default()
            })
            .is_err()
    );
    assert!(solid.rust_literal_npt_contains([f64::NAN, 0., 0.]).is_err());
    assert!(!solid.supports_step_export());
}
#[test]
fn actual_normal_mixed_torx_and_literal_npt_preserves_both_orders() {
    use multicad_rust_kernel::{
        head_profile::Head,
        helical_thread::{Spec, Thread},
        hexalobular_drive::{Cone, Drive, TangentSpec},
        profiled_screw::Screw,
        screw_recess::Recess,
        step::{Writer, npt_stock as step},
    };
    let drive = Drive::from_authored_tangent(
        TangentSpec {
            major_diameter_mm: 2.8,
            minor_diameter_mm: 2.05,
            outer_lobe_radius_mm: 0.28,
        },
        Cone {
            tip_depth_mm: 1.1 + 1.4 / 3f64.sqrt(),
            depth_per_radius: 1. / 3f64.sqrt(),
        },
    )
    .unwrap();
    let recess = Recess::conical_hexalobular(drive, 1e-7, 100_000).unwrap();
    let thread = Thread::new(Spec {
        major_radius_mm: 1.5,
        pitch_mm: 0.5,
        length_mm: 3.,
        depth_mm: 0.2,
        crest_width_mm: 0.0625,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts: 1,
        left_hand: false,
        root_radius_mm: None,
    })
    .unwrap();
    let torx = Screw::new_with_recess(
        thread,
        1.5,
        Head::Crowned {
            radius_mm: 2.8,
            height_mm: 2.4,
            crown_radius_mm: 5.,
            top_land_radius_mm: 2.05,
            neck_radius_mm: 0.125,
        },
        0.,
        0.1,
        recess,
    )
    .unwrap();
    let source = stock(false);
    let curves = step::CurvePlan::new(&source, Default::default(), None).unwrap();
    let surfaces = step::SurfacePlan::new(curves, Default::default(), None).unwrap();
    let verified = surfaces.verify_candidate(Default::default(), None).unwrap();
    for torx_first in [true, false] {
        let mut writer = Writer::new();
        for append_torx in [torx_first, !torx_first] {
            if append_torx {
                writer
                    .append_profiled_screw(&torx, [-20., 0., 0.], [1., 0., 0.], [0., 1., 0.], 1e-7)
                    .unwrap();
            } else {
                writer
                    .append_verified_source_npt_candidate(&verified, Default::default(), None)
                    .unwrap();
            }
        }
        let text = writer.finish().unwrap();
        let owner = EntityId::new();
        let solid = Solid::import_rust_step(&text, owner, Default::default()).unwrap();
        let members: Vec<_> = solid.rust_compound().unwrap().members().collect();
        assert_eq!(members.len(), 2);
        assert_eq!(
            members
                .iter()
                .filter(|m| m.rust_literal_npt_source().is_some())
                .count(),
            1
        );
        let physical: Vec<_> = members
            .iter()
            .filter_map(|m| m.rust_profiled_screw_geometry())
            .collect();
        assert_eq!(physical.len(), 1);
        let Recess::ConicalHexalobular(recess) = physical[0].0.recess() else {
            panic!("actual physical Torx recess");
        };
        assert!(recess.finite_binding().is_some());
        assert!(
            solid
                .topology_with_display(owner, 0, Some(0.1))
                .unwrap()
                .edges
                .len()
                > 20
        );
        assert!(solid.volume_mm3().unwrap() > verified.literal_material().world_volume_mm3()[1]);
    }
}
#[test]
fn actual_original_scaled_mapping_keeps_usable_cap_frames_and_exact_identity() {
    let placed = mapped(&text(false));
    for scale in ["0.1", "3.", "0.3333333333333333"] {
        let scaled = placed.replace("#900014,1.)", &format!("#900014,{scale})"));
        assert_ne!(scaled, placed);
        let material = source(&scaled);
        let solid = Solid::from_literal_npt_source(
            material.clone(),
            EntityId::new(),
            Default::default(),
            None,
        )
        .unwrap();
        let body = solid.literal_npt_body().unwrap();
        for (index, face) in material.faces().iter().enumerate() {
            if matches!(face.support(), read::Support::Plane { .. }) {
                let frame = body.plane(index).expect("earned usable scaled cap frame");
                assert!((frame.u.into_iter().map(|x| x * x).sum::<f64>() - 1.).abs() < 1e-12);
                assert!((frame.v.into_iter().map(|x| x * x).sum::<f64>() - 1.).abs() < 1e-12);
                assert!(
                    solid
                        .topology_snapshot(EntityId::new(), 0)
                        .unwrap()
                        .elements[index]
                        .geometry
                        .is_some()
                );
            }
        }
    }
    let original_mapping = source(&placed);
    let mut rebased = placed.clone();
    for id in 900010..=900020 {
        rebased = rebased.replace(&format!("#{id}"), &format!("#{}", id + 1_000_000));
    }
    let rebased_mapping = source(&rebased);
    assert_eq!(
        original_mapping.original_mapping_identity(),
        rebased_mapping.original_mapping_identity()
    );
    assert_eq!(
        source_signature(&original_mapping).unwrap(),
        source_signature(&rebased_mapping).unwrap()
    );
    let a = source(&placed.replace(
        "DIRECTION('',(0.,2.,0.))",
        "DIRECTION('',(0.2500000000000000001,2.,0.))",
    ));
    let b = source(&placed.replace(
        "DIRECTION('',(0.,2.,0.))",
        "DIRECTION('',(0.2500000000000000002,2.,0.))",
    ));
    assert_eq!(a.linear_columns_bounds(), b.linear_columns_bounds());
    assert_eq!(a.origin_mm(), b.origin_mm());
    assert_ne!(a.original_mapping_identity(), b.original_mapping_identity());
    assert_ne!(source_signature(&a).unwrap(), source_signature(&b).unwrap());
    // This reader deliberately admits cap anchors only on local X=Y=0.
    // Preserve an authentic translated-origin refusal rather than fabricating
    // a retained source outside that sealed original-file contract.
    let original = text(false);
    let material = source(&original);
    let cap = material
        .faces()
        .iter()
        .find_map(|f| {
            if let read::Support::Plane { height } = f.support() {
                Some(height)
            } else {
                None
            }
        })
        .unwrap();
    let point = material.point(cap).unwrap();
    let key = format!("#{}=", point.entity());
    let line = original.lines().find(|l| l.starts_with(&key)).unwrap();
    let replacement = format!(
        "{key}CARTESIAN_POINT('',(1.,2.,{}));",
        point.coordinates()[2].spelling()
    );
    let changed = original.replacen(line, &replacement, 1);
    assert_eq!(
        Solid::import_rust_step(&changed, EntityId::new(), Default::default())
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Unsupported
    );
}
#[test]
fn actual_warm_display_survives_failed_quality_and_fresh_manufacturing_queries() {
    let solid = Solid::from_literal_npt_source(
        source(&text(false)),
        EntityId::new(),
        Default::default(),
        None,
    )
    .unwrap();
    let original = solid.mesh(Default::default()).unwrap();
    let warm = solid.rust_literal_npt_runtime_evidence().unwrap();
    let body = solid.literal_npt_body().unwrap();
    let cache = || {
        let cache = body.display.borrow();
        let cache = cache.as_ref().unwrap();
        (
            cache.key,
            LiteralNptBody::mesh_bytes(&cache.mesh),
            cache.work,
            cache.peak,
        )
    };
    let old_cache = cache();
    assert_eq!(
        solid
            .mesh(MeshOptions {
                maximum_triangles: 1,
                ..Default::default()
            })
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Limit
    );
    let refused = solid.rust_literal_npt_runtime_evidence().unwrap();
    assert_eq!(cache(), old_cache);
    assert_eq!(refused.mesh_builds, warm.mesh_builds);
    assert_eq!(refused.resident_budget_bytes, warm.resident_budget_bytes);
    assert_eq!(refused.last_mesh_work, warm.last_mesh_work);
    assert_eq!(refused.last_mesh_actual_work, warm.last_mesh_actual_work);
    assert_eq!(refused.last_mesh_peak_bytes, warm.last_mesh_peak_bytes);
    assert_eq!(
        solid.mesh(Default::default()).unwrap().vertices,
        original.vertices
    );
    assert_eq!(
        solid
            .rust_literal_npt_runtime_evidence()
            .unwrap()
            .mesh_builds,
        warm.mesh_builds
    );
    let changed = MeshOptions {
        linear_deflection_mm: 0.11,
        ..Default::default()
    };
    solid.mesh(changed).unwrap();
    let replaced = cache();
    assert_ne!(replaced.0, old_cache.0);
    let new_warm = solid.rust_literal_npt_runtime_evidence().unwrap();
    assert_eq!(new_warm.mesh_builds, warm.mesh_builds + 1);
    solid
        .manufacturing_mesh_with_receipt(Default::default())
        .unwrap();
    let manufactured = solid.rust_literal_npt_runtime_evidence().unwrap();
    assert_eq!(cache(), replaced);
    assert_eq!(manufactured.mesh_builds, new_warm.mesh_builds + 1);
    assert_eq!(
        manufactured.resident_budget_bytes,
        new_warm.resident_budget_bytes
    );
    solid.mesh(changed).unwrap();
    let hot = solid.rust_literal_npt_runtime_evidence().unwrap();
    assert_eq!(hot.mesh_builds, manufactured.mesh_builds);
    assert!(hot.last_mesh_actual_work.unwrap() < manufactured.last_mesh_actual_work.unwrap());
    assert!(
        hot.last_mesh_peak_bytes.unwrap() <= BYTES
            && manufactured.last_mesh_peak_bytes.unwrap() <= BYTES
    );
}

#[test]
fn actual_native_occurrences_keep_source_and_fresh_mesh_face_query_quality() {
    for placed in [false, true] {
        let original = text(false);
        let original = if placed { mapped(&original) } else { original };
        let solid =
            Solid::import_rust_step(&original, EntityId::new(), Default::default()).unwrap();
        let old_faces = solid.faces().unwrap();
        let rotate = RigidTransform {
            rotation: [[0.8, -0.36, 0.48], [0.6, 0.48, -0.64], [0., 0.8, 0.6]],
            translation_mm: [17., -8., 13.],
        };
        let first = solid.transformed(rotate).unwrap();
        let moved = first
            .transformed(RigidTransform {
                rotation: [[0., 0., 1.], [0., 1., 0.], [-1., 0., 0.]],
                translation_mm: [2., 3., 4.],
            })
            .unwrap();
        assert!(Rc::ptr_eq(
            &solid.literal_npt_body().unwrap().plan,
            &moved.literal_npt_body().unwrap().plan
        ));
        assert_eq!(
            moved
                .rust_literal_npt_occurrence()
                .unwrap()
                .transforms()
                .len(),
            2
        );
        assert_eq!(
            solid.rust_literal_npt_source().unwrap().point_count(),
            moved.rust_literal_npt_source().unwrap().point_count()
        );
        let faces = moved.faces().unwrap();
        assert_eq!(faces.len(), old_faces.len());
        let pose = moved.rust_literal_npt_occurrence().unwrap();
        for (old, new) in old_faces.iter().zip(&faces) {
            let center = pose.map_box_mm(old.center_mm.map(|v| [v, v])).unwrap();
            assert!(
                new.center_mm
                    .into_iter()
                    .zip(center)
                    .all(|(x, b)| x >= b[0] - 1e-7 && x <= b[1] + 1e-7)
            );
            assert!((new.area_mm2 - old.area_mm2).abs() < 1e-6);
        }
        let local = solid
            .rust_literal_npt_source()
            .unwrap()
            .map_box_mm([[0., 0.], [0., 0.], [2., 2.]])
            .unwrap();
        let world = pose.map_box_mm(local).unwrap().map(|b| b[0].midpoint(b[1]));
        assert!(moved.rust_literal_npt_contains(world).unwrap());
        assert!(!moved.rust_literal_npt_contains([100., 100., 100.]).unwrap());
        let top = moved
            .topology_with_display(EntityId::new(), 0, Some(0.05))
            .unwrap();
        assert_eq!(top.edges.len(), 20);
        assert_eq!(top.vertices.len(), 12);
        let mesh = moved.mesh(Default::default()).unwrap();
        assert!(!mesh.indices.is_empty());
        let before = moved.rust_literal_npt_runtime_evidence().unwrap();
        let mfg = moved
            .manufacturing_mesh_with_receipt(Default::default())
            .unwrap();
        assert!(mfg.edge_pairs > 0 && mfg.vertex_links > 0);
        assert!(
            moved
                .rust_literal_npt_runtime_evidence()
                .unwrap()
                .mesh_builds
                > before.mesh_builds
        );
        assert!((moved.volume_mm3().unwrap() - solid.volume_mm3().unwrap()).abs() < 1e-6);
        assert_eq!(
            solid
                .rust_literal_npt_occurrence()
                .map(|p| p.transforms().len()),
            None
        );
    }
}
