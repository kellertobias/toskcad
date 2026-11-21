use super::*;
fn spec(quarter: bool) -> multicad_rust_kernel::axial_pipe_thread::authored::AuthoredSpec {
    use multicad_rust_kernel::axial_pipe_thread::authored::{AuthoredSection, AuthoredSpec};
    let (tpi, e1, l1, l2, od, max) = if quarter {
        (18., 0.49163, 0.2278, 0.4018, 0.540, 0.088)
    } else {
        (27., 0.37360, 0.1615, 0.2639, 0.405, 0.096)
    };
    let p = 25.4 / tpi;
    let k = 1. / 32.;
    let m = 3f64.sqrt();
    let s = 1. / (m - k) + 1. / (m + k);
    let d = p / s - (0.033 + max) * p;
    AuthoredSpec {
        major_radius_mm: (e1 * 25.4 - 2. * k * l1 * 25.4 + d) / 2.,
        pitch_mm: p,
        length_mm: l2 * 25.4,
        radial_taper: k,
        section: AuthoredSection::NptFlat {
            depth_mm: d,
            flat_width_mm: (p - d * s) / 2.,
        },
        crest_radius_limit_mm: Some(od * 25.4 / 2.),
    }
}
#[test]
fn actual_canonical_clipped_npt_native_source_has_finite_topology_properties_and_mfg() {
    for quarter in [false, true] {
        let spec = spec(quarter);
        let feature = EntityId::new();
        let policy = GeometryPolicy::default();
        let solid = Solid::rust_authored_axial_pipe_thread_with_interrupt(
            spec,
            feature,
            PlaneFrame::default(),
            policy,
            None,
        )
        .unwrap();
        assert_eq!(
            solid.rust_clipped_npt_source().unwrap().source().spec(),
            spec
        );
        assert!(solid.rust_authored_axial_pipe_thread_source().is_none());
        assert!(solid.rust_radial_boundary().is_none());
        assert_eq!(solid.faces().unwrap().len(), 7);
        let top = solid.topology_snapshot(feature, 1).unwrap();
        for (kind, count) in [
            (TopologyKind::Face, 7),
            (TopologyKind::Edge, 17),
            (TopologyKind::Vertex, 12),
        ] {
            assert_eq!(
                top.elements.iter().filter(|e| e.key.kind == kind).count(),
                count
            );
        }
        assert!(solid.volume_mm3().unwrap() > 0.);
        assert!(solid.supports_step_export());
        let file = tempfile::NamedTempFile::new().unwrap();
        solid.write_step(file.path()).unwrap();
        let text = std::fs::read_to_string(file.path()).unwrap();
        {
            use multicad_rust_kernel::step_read::{
                npt_stock::Selected, owned_insert::CompoundDocument,
            };
            let document = CompoundDocument::parse(&text, Default::default(), None).unwrap();
            assert_eq!(document.represented_members().len(), 1);
            let selected = Selected::parse(
                &document,
                document.represented_members()[0].item_entity(),
                Default::default(),
                None,
            )
            .unwrap();
            let material = selected.retain(Default::default(), None).unwrap();
            let volume = material.local_volume_mm3();
            let source = solid.rust_clipped_npt_volume_bounds_mm3().unwrap();
            assert!((volume[0].midpoint(volume[1]) - source[0].midpoint(source[1])).abs() < 1e-5);
            assert_eq!(material.faces().len(), 10);
            assert_eq!(material.edges().len(), 20);
        }
        // Repeated capability queries and exports retain the same immutable
        // authored source and deterministic independently qualified geometry.
        assert!(solid.supports_step_export());
        solid.write_step(file.path()).unwrap();
        assert_eq!(std::fs::read_to_string(file.path()).unwrap(), text);
        assert_eq!(
            solid.rust_clipped_npt_source().unwrap().source().spec(),
            spec
        );
        let mesh = solid
            .manufacturing_mesh_with_receipt(MeshOptions {
                linear_deflection_mm: 0.1,
                angular_deflection_rad: 0.3,
                maximum_triangles: 1_000_000,
            })
            .unwrap();
        assert!(mesh.edge_pairs > 0 && mesh.vertex_links > 0);
        assert_eq!(mesh.maximum_displacement_mm, 0.);
        assert!(solid.rust_clipped_npt_contains([0., 0., 1.]).unwrap());
        assert!(!solid.rust_clipped_npt_contains([10., 0., 1.]).unwrap());
        assert!(solid.rust_clipped_npt_contains([0., 0., 0.]).is_err());
        let share = solid.share().unwrap();
        assert!(std::ptr::eq(
            solid.rust_clipped_npt_source().unwrap(),
            share.rust_clipped_npt_source().unwrap()
        ));
    }
}

#[test]
fn actual_required_radial_allowance_edit_keeps_native_source_properties() {
    for quarter in [false, true] {
        let mut s = spec(quarter);
        s.major_radius_mm -= 0.001;
        let result = Solid::rust_authored_axial_pipe_thread_with_interrupt(
            s,
            EntityId::new(),
            PlaneFrame::default(),
            GeometryPolicy::default(),
            None,
        );
        assert!(
            result.is_ok(),
            "actual edit {quarter}: {:?}",
            result.as_ref().err()
        );
        result
            .unwrap()
            .mesh(MeshOptions {
                linear_deflection_mm: 0.1,
                angular_deflection_rad: 15f64.to_radians(),
                maximum_triangles: 2_000_000,
            })
            .unwrap();
    }
}

#[test]
fn actual_npt_export_receipts_keep_member_mass_and_reject_a_different_source() {
    let policy = GeometryPolicy::default();
    let parts = [false, true].map(|quarter| {
        Solid::rust_authored_axial_pipe_thread_with_interrupt(
            spec(quarter),
            EntityId::new(),
            PlaneFrame::default(),
            policy,
            None,
        )
        .unwrap()
    });
    let ordinary = Solid::box_solid([1., 2., 3.], policy)
        .unwrap()
        .transformed(RigidTransform {
            translation_mm: [30., 0., 0.],
            ..Default::default()
        })
        .unwrap();
    let compound =
        Solid::compound_named(&[&parts[0], &ordinary, &parts[1]], EntityId::new(), policy).unwrap();
    assert!(compound.has_partitioned_step_source());
    let file = tempfile::NamedTempFile::new().unwrap();
    compound.write_step(file.path()).unwrap();
    let text = std::fs::read_to_string(file.path()).unwrap();
    let proof = compound
        .verify_partitioned_step_source(&text, 1e-7, None)
        .unwrap()
        .unwrap();
    assert_eq!(proof.members.solid_count, 3);
    assert_eq!(proof.members.retained_members.len(), 2);
    assert_eq!(proof.members.decoded_scalar_members, [(1, 6.)]);
    assert_eq!(proof.insert_correspondence.len(), 2);
    assert_eq!(proof.insert_correspondence[0].0, 0);
    assert_eq!(proof.insert_correspondence[1].0, 2);
    for (_, correspondence) in &proof.insert_correspondence {
        assert!(correspondence.maximum_position_error_mm <= 1e-7);
        assert!(correspondence.maximum_row_derivative_error_mm <= 1e-6);
        assert!(correspondence.maximum_ruling_derivative_error_mm <= 1e-6);
    }
    assert!(matches!(
        compound.verify_partitioned_step_source(&text, 1e-7, Some(&|| true)),
        Some(Err(GeometryError {
            kind: GeometryErrorKind::Cancelled,
            ..
        }))
    ));
    parts[1].write_step(file.path()).unwrap();
    let different = std::fs::read_to_string(file.path()).unwrap();
    assert!(
        parts[0]
            .verify_partitioned_step_source(&different, 1e-7, None)
            .unwrap()
            .is_err()
    );
    let placed = parts[1]
        .transformed(RigidTransform {
            translation_mm: [1., 0., 0.],
            ..Default::default()
        })
        .unwrap();
    assert!(matches!(
        placed.verify_partitioned_step_source(&different, 1e-7, None),
        Some(Err(GeometryError {
            kind: GeometryErrorKind::Unsupported,
            ..
        }))
    ));
    assert!(placed.write_step(file.path()).is_err());
    assert_eq!(std::fs::read_to_string(file.path()).unwrap(), different);
}
#[test]
fn native_source_scalars_do_not_require_unrequested_fixed_absolute_precision() {
    let mut s = spec(false);
    let factor = 100.;
    s.major_radius_mm *= factor;
    s.pitch_mm *= factor;
    s.length_mm *= factor;
    s.crest_radius_limit_mm = s.crest_radius_limit_mm.map(|r| r * factor);
    let multicad_rust_kernel::axial_pipe_thread::authored::AuthoredSection::NptFlat {
        depth_mm,
        flat_width_mm,
    } = s.section
    else {
        unreachable!()
    };
    s.section = multicad_rust_kernel::axial_pipe_thread::authored::AuthoredSection::NptFlat {
        depth_mm: depth_mm * factor,
        flat_width_mm: flat_width_mm * factor,
    };
    let solid = Solid::rust_authored_axial_pipe_thread_with_interrupt(
        s,
        EntityId::new(),
        PlaneFrame::default(),
        GeometryPolicy::default(),
        None,
    )
    .unwrap();
    assert_eq!(solid.rust_clipped_npt_source().unwrap().source().spec(), s);
    assert!(solid.volume_mm3().unwrap() > 0.);
    assert!(solid.faces().unwrap().iter().all(|f| f.area_mm2 > 0.));
}
#[test]
fn bounded_native_source_vertices_precede_expensive_property_queries() {
    use multicad_rust_kernel::axial_pipe_thread::authored::{AuthoredSection, AuthoredSpec};
    let s = AuthoredSpec {
        major_radius_mm: 300_000.,
        pitch_mm: 0.001,
        length_mm: 0.126,
        radial_taper: 1. / 32.,
        section: AuthoredSection::NptFlat {
            depth_mm: 0.0007,
            flat_width_mm: 0.00005,
        },
        crest_radius_limit_mm: Some(300_000. + 0.0035875),
    };
    let e = Solid::rust_authored_axial_pipe_thread_with_interrupt(
        s,
        EntityId::new(),
        PlaneFrame::default(),
        GeometryPolicy::default(),
        None,
    )
    .err()
    .unwrap();
    assert_eq!(e.kind, GeometryErrorKind::Unsupported);
    assert_eq!(e.operation, "clipped authored NPT");
}
#[test]
fn actual_source_references_edges_cap_wires_and_scalar_enclosures_stay_owned() {
    let policy = GeometryPolicy::default();
    let feature = EntityId::new();
    let solid = Solid::rust_authored_axial_pipe_thread_with_interrupt(
        spec(false),
        feature,
        PlaneFrame::default(),
        policy,
        None,
    )
    .unwrap();
    let snapshot = solid.topology_snapshot(feature, 9).unwrap();
    let edges = solid.display_edges(0.02).unwrap();
    assert_eq!(edges.len(), 17);
    let source = solid.rust_clipped_npt_source().unwrap();
    for edge in edges {
        let key = source.physical_curves()[edge.index as usize];
        for (i, v) in key.endpoints().into_iter().enumerate() {
            let p = source
                .vertex_bounds(v, Default::default(), None)
                .unwrap()
                .position_mm()
                .map(|x| x[0].midpoint(x[1]));
            assert_eq!(
                if i == 0 {
                    edge.points_mm[0]
                } else {
                    *edge.points_mm.last().unwrap()
                },
                p
            );
        }
    }
    let body = solid.clipped_npt_body().unwrap();
    for (index, face) in source.physical_faces().into_iter().enumerate() {
        let element = snapshot
            .elements
            .iter()
            .find(|e| {
                e.key
                    == TopologyKey {
                        kind: TopologyKind::Face,
                        ordinal: index as u32,
                    }
            })
            .unwrap();
        if matches!(face, StockClipFace::StartCap | StockClipFace::EndCap) {
            let wire = body.planar_loops(index as u32);
            assert_eq!(wire.len(), 1);
            assert_eq!(wire[0].len(), 4);
            assert_eq!(element.frames.len(), 1);
            let normal = element.frames[0].frame.normal();
            assert_eq!(
                normal,
                [
                    0.,
                    0.,
                    if face == StockClipFace::StartCap {
                        -1.
                    } else {
                        1.
                    }
                ]
            );
        } else {
            assert!(body.planar_loops(index as u32).is_empty());
            assert!(element.frames.is_empty());
            assert!(element.geometry.is_none());
        }
    }
    assert!(body.planar_loops(u32::MAX).is_empty());
    let v = solid.rust_clipped_npt_volume_bounds_mm3().unwrap();
    let nominal = solid.volume_mm3().unwrap();
    assert!(v[0] <= nominal && nominal <= v[1]);
    let measures = solid.rust_clipped_npt_face_measures().unwrap();
    let faces = solid.faces().unwrap();
    assert_eq!(measures.len(), faces.len());
    for (m, f) in measures.iter().zip(faces) {
        assert!(m.area_mm2()[0] <= f.area_mm2 && f.area_mm2 <= m.area_mm2()[1]);
        for i in 0..3 {
            assert!(
                m.centroid_mm()[i][0] <= f.center_mm[i] && f.center_mm[i] <= m.centroid_mm()[i][1]
            );
        }
    }
    let file = tempfile::NamedTempFile::new().unwrap();
    std::fs::write(file.path(), b"preserved").unwrap();
    solid.write_step(file.path()).unwrap();
    let exported = std::fs::read_to_string(file.path()).unwrap();
    let proof = solid
        .verify_partitioned_step_source(&exported, policy.linear_tolerance_mm, None)
        .unwrap()
        .unwrap();
    assert_eq!(proof.members.retained_members.len(), 1);
    assert_eq!(proof.insert_correspondence.len(), 1);
    assert_eq!(solid.rust_clipped_npt_volume_bounds_mm3().unwrap(), v);
    let foreign =
        Solid::rust_sphere_named(1., PlaneFrame::default(), EntityId::new(), policy).unwrap();
    assert!(foreign.rust_clipped_npt_source().is_none());
    assert!(foreign.rust_clipped_npt_face_measures().is_none());
    assert!(foreign.rust_clipped_npt_volume_bounds_mm3().is_none());
    assert_eq!(
        foreign.rust_clipped_npt_contains([0.; 3]).unwrap_err().kind,
        GeometryErrorKind::Unsupported
    );
}
#[test]
fn native_source_cancel_pose_input_and_mesh_limits_fail_before_publication() {
    let s = spec(false);
    let policy = GeometryPolicy::default();
    let feature = EntityId::new();
    let run = |frame, policy, cancel: Option<&dyn Fn() -> bool>| {
        Solid::rust_authored_axial_pipe_thread_with_interrupt(s, feature, frame, policy, cancel)
    };
    let solid = run(PlaneFrame::default(), policy, None).unwrap();
    let original = solid.topology_snapshot(feature, 0).unwrap();
    let authored_pose = run(
        PlaneFrame {
            origin: [1., 0., 0.],
            ..Default::default()
        },
        policy,
        None,
    )
    .unwrap();
    assert_eq!(
        authored_pose
            .rust_clipped_npt_placement()
            .unwrap()
            .translation_mm,
        [1., 0., 0.]
    );
    let moved = solid
        .transformed(RigidTransform {
            translation_mm: [1., 2., 3.],
            ..Default::default()
        })
        .unwrap();
    assert_eq!(
        moved.rust_clipped_npt_placement().unwrap().translation_mm,
        [1., 2., 3.]
    );
    assert_eq!(solid.topology_snapshot(feature, 0).unwrap(), original);
    assert_eq!(
        solid
            .transformed(RigidTransform::default())
            .unwrap()
            .volume_mm3()
            .unwrap(),
        solid.volume_mm3().unwrap()
    );
    assert_eq!(
        run(
            PlaneFrame::default(),
            GeometryPolicy {
                maximum_coordinate_mm: 3.,
                ..policy
            },
            None
        )
        .err()
        .unwrap()
        .kind,
        GeometryErrorKind::Input
    );
    assert_eq!(
        run(PlaneFrame::default(), policy, Some(&|| true))
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Cancelled
    );
    let polls = std::cell::Cell::new(0usize);
    let observe = || {
        polls.set(polls.get() + 1);
        false
    };
    run(PlaneFrame::default(), policy, Some(&observe)).unwrap();
    for limit in [polls.get() / 2, polls.get() - 1] {
        let n = std::cell::Cell::new(0usize);
        let cancel = || {
            let at = n.get();
            n.set(at + 1);
            at >= limit
        };
        assert_eq!(
            run(PlaneFrame::default(), policy, Some(&cancel))
                .err()
                .unwrap()
                .kind,
            GeometryErrorKind::Cancelled
        );
    }
    assert_eq!(
        solid
            .rust_clipped_npt_contains([f64::NAN, 0., 0.])
            .unwrap_err()
            .kind,
        GeometryErrorKind::Input
    );
    // User decision 2026-11-21: the dedicated law mesher meets requests the
    // former cell-plan mesher refused. A one-triangle frontier and the
    // finest admitted quality under a tighter triangle budget still refuse.
    for q in [
        MeshOptions {
            linear_deflection_mm: 0.1,
            angular_deflection_rad: 0.3,
            maximum_triangles: 1,
        },
        MeshOptions {
            linear_deflection_mm: 0.001,
            angular_deflection_rad: 0.01,
            maximum_triangles: 20_000,
        },
    ] {
        assert_eq!(
            solid.manufacturing_mesh(q).unwrap_err().kind,
            GeometryErrorKind::Limit
        );
    }
    assert_print_quality(
        &solid,
        MeshOptions {
            linear_deflection_mm: 0.001,
            angular_deflection_rad: 0.01,
            maximum_triangles: 1_000_000,
        },
    );
    assert_eq!(
        solid
            .mesh(MeshOptions {
                linear_deflection_mm: f64::NAN,
                angular_deflection_rad: 0.3,
                maximum_triangles: 1_000_000
            })
            .unwrap_err()
            .kind,
        GeometryErrorKind::Input
    );
    assert_eq!(solid.topology_snapshot(feature, 0).unwrap(), original);
}

#[test]
fn source_edge_sampler_enforces_its_actual_numeric_and_storage_contract() {
    let solid = Solid::rust_authored_axial_pipe_thread_with_interrupt(
        spec(false),
        EntityId::new(),
        PlaneFrame::default(),
        GeometryPolicy::default(),
        None,
    )
    .unwrap();
    let bounds = solid.bounds().unwrap();
    assert!(bounds.max[2] > bounds.min[2]);
    let q = Tessellation {
        deflection: 0.1,
        angle: 0.3,
        triangle_limit: 1000,
    };
    // The internal Native mesh entry also accepts raw numeric quality. Test
    // that contract with an actual admitted source, independent of the public
    // Solid option validator; no invalid source or receipt is constructed.
    assert_eq!(
        solid
            .clipped_npt_body()
            .unwrap()
            .mesh(Tessellation {
                deflection: f64::NAN,
                ..q
            })
            .unwrap_err()
            .kind,
        GeometryErrorKind::Input
    );
    assert!(solid.sample_non_analytic_edge(u32::MAX, q).is_none());
    assert_eq!(
        solid
            .sample_non_analytic_edge(
                0,
                Tessellation {
                    deflection: f64::NAN,
                    ..q
                }
            )
            .unwrap()
            .unwrap_err(),
        MeshError::InvalidOptions
    );
    assert_eq!(
        solid
            .sample_non_analytic_edge(
                0,
                Tessellation {
                    triangle_limit: 1,
                    ..q
                }
            )
            .unwrap()
            .unwrap_err(),
        MeshError::InvalidOptions
    );
    assert_eq!(
        solid
            .sample_non_analytic_edge(
                0,
                Tessellation {
                    triangle_limit: 2,
                    ..q
                }
            )
            .unwrap()
            .unwrap_err(),
        MeshError::BudgetExceeded
    );
    assert_eq!(
        solid
            .sample_non_analytic_edge(
                0,
                Tessellation {
                    deflection: f64::MIN_POSITIVE,
                    ..q
                }
            )
            .unwrap()
            .unwrap_err(),
        MeshError::UnresolvedBoundary
    );
    assert!(solid.sample_non_analytic_edge(0, q).unwrap().is_ok());
    // User decision 2026-11-21: this request is now met, not refused.
    assert_print_quality(
        &solid,
        MeshOptions {
            linear_deflection_mm: 0.003,
            angular_deflection_rad: 0.1,
            maximum_triangles: 1_000_000,
        },
    );
}

/// The manufacturing mesh is one closed shell, the law tessellation behind
/// it proves a deviation within the request, and every emitted lateral
/// vertex lies within that deviation of the source's own radial enclosure.
fn assert_print_quality(solid: &Solid, q: MeshOptions) {
    use multicad_rust_kernel::axial_pipe_thread::authored::LawClosure;
    let mesh = solid.manufacturing_mesh_with_receipt(q).unwrap();
    let triangles = mesh.mesh.indices.len() / 3;
    assert!(triangles > 0 && triangles <= q.maximum_triangles as usize);
    assert_eq!(mesh.edge_pairs as usize, 3 * triangles / 2);
    let source = solid.rust_clipped_npt_source().unwrap().source();
    let law = source
        .tessellate_law(
            Tessellation {
                deflection: q.linear_deflection_mm,
                angle: q.angular_deflection_rad,
                triangle_limit: q.maximum_triangles as usize,
            },
            LawClosure::Solid,
            None,
        )
        .unwrap();
    assert!(law.deviation_mm <= q.linear_deflection_mm);
    assert!(law.normal_angle_rad <= q.angular_deflection_rad);
    assert_eq!(law.mesh.face_indices.len(), triangles);
    let spec = source.spec();
    let omega = std::f64::consts::TAU / spec.pitch_mm;
    for v in mesh.mesh.vertices.as_chunks::<3>().0 {
        let r = v[0].hypot(v[1]);
        if r < 1e-9 {
            continue; // Cap fan centre.
        }
        let theta = v[1].atan2(v[0]);
        let phase = (v[2] - theta / omega).rem_euclid(spec.pitch_mm);
        let [lo, hi] = source
            .radial_bounds([phase, phase], [v[2].clamp(0., spec.length_mm); 2])
            .unwrap();
        assert!(
            r >= lo - law.deviation_mm && r <= hi + law.deviation_mm,
            "{r} not in [{lo}, {hi}]"
        );
    }
}

#[test]
fn actual_clipped_source_poses_preserve_material_topology_edges_and_closed_manufacturing() {
    let policy = GeometryPolicy::default();
    for quarter in [false, true] {
        let feature = EntityId::new();
        let spec = spec(quarter);
        let original = Solid::rust_authored_axial_pipe_thread_with_interrupt(
            spec,
            feature,
            PlaneFrame::default(),
            policy,
            None,
        )
        .unwrap();
        let topology = original.topology_snapshot(feature, 0).unwrap();
        let volume = original.volume_mm3().unwrap();
        for rotation in [
            RigidTransform::default().rotation,
            [[0., 0., 1.], [0., 1., 0.], [-1., 0., 0.]],
            [[0.8, -0.36, 0.48], [0.6, 0.48, -0.64], [0., 0.8, 0.6]],
        ] {
            let pose = RigidTransform {
                rotation,
                translation_mm: [17., -8., 13.],
            };
            let moved = original.transformed(pose).unwrap();
            assert_eq!(
                moved.rust_clipped_npt_source().unwrap().source().spec(),
                spec
            );
            assert_eq!(moved.rust_clipped_npt_placement(), Some(pose));
            assert!((moved.volume_mm3().unwrap() - volume).abs() < 1e-6);
            let snapshot = moved.topology_snapshot(feature, 1).unwrap();
            assert_eq!(snapshot.elements.len(), topology.elements.len());
            for (a, b) in snapshot.elements.iter().zip(&topology.elements) {
                assert_eq!(a.key, b.key);
                assert_eq!(a.origins, b.origins);
                assert_eq!(a.boundary, b.boundary);
            }
            let point = |p: [f64; 3]| {
                let q = pose.vector(p);
                std::array::from_fn(|i| q[i] + pose.translation_mm[i])
            };
            assert!(
                moved
                    .rust_clipped_npt_contains(point([0., 0., spec.length_mm / 2.]))
                    .unwrap()
            );
            assert!(
                !moved
                    .rust_clipped_npt_contains(point([20., 0., spec.length_mm / 2.]))
                    .unwrap()
            );
            let bounds = moved.bounds().unwrap();
            for face in moved.faces().unwrap() {
                for i in 0..3 {
                    assert!(
                        bounds.min[i] <= face.center_mm[i] && face.center_mm[i] <= bounds.max[i]
                    );
                }
            }
            let edges = moved.display_edges(0.1).unwrap();
            assert!(!edges.is_empty());
            for edge in &edges {
                for p in &edge.points_mm {
                    for (i, coordinate) in p.iter().enumerate() {
                        assert!(bounds.min[i] <= *coordinate && *coordinate <= bounds.max[i]);
                    }
                }
            }
            let mesh = moved
                .manufacturing_mesh_with_receipt(MeshOptions {
                    linear_deflection_mm: 0.1,
                    angular_deflection_rad: 0.3,
                    maximum_triangles: 1_000_000,
                })
                .unwrap();
            assert!(mesh.edge_pairs > 0 && mesh.vertex_links > 0);
            assert!(!moved.supports_step_export());
            assert_eq!(original.topology_snapshot(feature, 0).unwrap(), topology);
        }
    }
}
