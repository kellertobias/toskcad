use super::*;
fn profile() -> HelicalThread {
    HelicalThread {
        major_radius_mm: 4.08,
        pitch_mm: 1.25,
        length_mm: 4.1,
        depth_mm: 0.65,
        crest_width_mm: 0.15,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        radial_taper: 0.,
        starts: 1,
        left_hand: false,
    }
}
fn relief() -> Relief {
    Relief {
        start_mm: 4.099,
        radius_mm: 4.1,
        tool_length_mm: 3.22,
    }
}
fn stock() -> Solid {
    let ring = [
        [-6.9, 0.],
        [6.9, 0.],
        [6.5, 3.],
        [5.4, 5.],
        [3.3, 7.3],
        [-3.3, 7.3],
        [-5.4, 5.],
        [-6.5, 3.],
    ];
    let outer = (0..8)
        .map(|i| crate::CurvedEdge2d::Line {
            a: ring[i],
            b: ring[(i + 1) % 8],
        })
        .collect();
    Solid::extrude_curved_profile(
        &crate::CurvedProfile {
            frame: PlaneFrame {
                origin: [-7.5, 0., 0.],
                u: [0., 1., 0.],
                v: [0., 0., 1.],
            },
            outer,
            holes: vec![],
        },
        23.,
        GeometryPolicy::default(),
    )
    .unwrap()
}
#[test]
fn real_native_retained_nut_has_source_measures_topology_and_independent_manufacturing_embedding() {
    let policy = GeometryPolicy::default();
    let stock = stock();
    let original = stock
        .linear_body()
        .unwrap()
        .0
        .straight_source()
        .unwrap()
        .clone();
    let solid = stock
        .rust_cut_line8_nut(
            Line8NutCut {
                profile: profile(),
                frame: PlaneFrame::default(),
                bore_from_start_mm: 7.5,
                relief: relief(),
                minimum_wall_mm: 0.5,
            },
            EntityId::new(),
            policy,
            None,
        )
        .unwrap();
    let source = solid.rust_line8_nut_source().unwrap();
    assert_eq!(source.stock().base(), original.base());
    assert_eq!(source.stock().height_mm(), original.height_mm());
    assert_eq!(source.thread().spec().length_mm, 4.1);
    assert_eq!(source.relief().start_mm, 4.099);
    assert!(solid.axial_cylinder().is_none());
    assert!(solid.linear_body().is_none());
    assert!(!solid.supports_step_export());
    assert_eq!(solid.faces().unwrap().len(), 17);
    let topology = solid.topology_snapshot(EntityId::new(), 1).unwrap();
    assert_eq!(topology.elements.len(), 95);
    assert!(solid.volume_mm3().unwrap() < stock.volume_mm3().unwrap());
    assert!(solid.rust_line8_nut_contains([6., 0., 2.]).unwrap());
    assert!(!solid.rust_line8_nut_contains([0., 0., 2.]).unwrap());
    assert!(!solid.rust_line8_nut_contains([0., 3., 7.2]).unwrap());
    let q = MeshOptions {
        linear_deflection_mm: 0.02,
        angular_deflection_rad: 0.2,
        maximum_triangles: 200000,
    };
    let t = std::time::Instant::now();
    let mesh = solid.manufacturing_mesh_with_receipt(q).unwrap();
    eprintln!(
        "Native Line8 independent manufacturing {}T {:?}",
        mesh.mesh.face_indices.len(),
        t.elapsed()
    );
    assert_eq!(mesh.vertex_links as usize, mesh.mesh.vertices.len() / 3);
    assert_eq!(mesh.maximum_displacement_mm, 0.);
    for i in 0..46 {
        let q = Solid::rust_mesh_options(q);
        let points = solid.sample_non_analytic_edge(i, q).unwrap().unwrap();
        assert!(points.len() >= 2);
    }
    let body = solid.line8_nut_body().unwrap();
    for i in 0..17 {
        assert_eq!(
            body.planar_loops(i).is_empty(),
            [11, 12, 13, 14, 16].contains(&i)
        );
    }
    let moved = solid
        .transformed(RigidTransform {
            rotation: [[0., 0., 1.], [0., 1., 0.], [-1., 0., 0.]],
            translation_mm: [3., -2., 7.],
        })
        .unwrap();
    assert!(moved.rust_line8_nut_contains([5., -2., 1.]).unwrap());
    moved.manufacturing_mesh(q).unwrap();
}
#[test]
fn actual_nut_stock_phase_cancel_and_wrong_target_refuse() {
    let p = GeometryPolicy::default();
    let stock = stock();
    assert_eq!(
        stock
            .rust_cut_line8_nut(
                Line8NutCut {
                    profile: profile(),
                    frame: PlaneFrame::default(),
                    bore_from_start_mm: 7.5,
                    relief: relief(),
                    minimum_wall_mm: 0.5
                },
                EntityId::new(),
                p,
                Some(&|| true)
            )
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Cancelled
    );
    let mut phase = PlaneFrame::default();
    phase.origin[0] = 1.;
    assert!(
        stock
            .rust_cut_line8_nut(
                Line8NutCut {
                    profile: profile(),
                    frame: phase,
                    bore_from_start_mm: 7.5,
                    relief: relief(),
                    minimum_wall_mm: 0.5
                },
                EntityId::new(),
                p,
                None
            )
            .is_err()
    );
    let wrong = Solid::cylinder(PlaneFrame::default(), 5., 8., p).unwrap();
    assert!(
        wrong
            .rust_cut_line8_nut(
                Line8NutCut {
                    profile: profile(),
                    frame: PlaneFrame::default(),
                    bore_from_start_mm: 7.5,
                    relief: relief(),
                    minimum_wall_mm: 0.5
                },
                EntityId::new(),
                p,
                None
            )
            .is_err()
    );
}
#[test]
fn actual_resolved_m8_nut_allowances_phase_and_minimum_wall_are_retained() {
    use multicad_core::{
        document::Document,
        model::{ThreadMode, ThreadParameters, length},
        parameters::evaluate,
    };
    let parameters = evaluate(&Document::new("actual Line8 M8 nut").unwrap()).unwrap();
    let request = ThreadParameters {
        length: length(4.1),
        radial_allowance: length(0.08),
        axial_allowance: length(0.04),
        minimum_wall: length(0.5),
        ..Default::default()
    };
    let resolved = request
        .resolve(&parameters, ThreadMode::CutInternal)
        .unwrap();
    let profile = HelicalThread {
        major_radius_mm: resolved.major_radius_mm,
        pitch_mm: resolved.pitch_mm,
        length_mm: resolved.length_mm,
        depth_mm: resolved.profile.depth_mm,
        crest_width_mm: resolved.profile.crest_width_mm,
        lower_flank_deg: resolved.profile.lower_flank_deg,
        upper_flank_deg: resolved.profile.upper_flank_deg,
        radial_taper: resolved.profile.radial_taper,
        starts: resolved.starts,
        left_hand: false,
    };
    let (s, c) = 0.317_f64.portable_sin_cos();
    let cut = Line8NutCut {
        profile,
        frame: PlaneFrame {
            u: [c, s, 0.],
            v: [-s, c, 0.],
            ..Default::default()
        },
        bore_from_start_mm: 7.5,
        relief: relief(),
        minimum_wall_mm: resolved.minimum_wall_mm,
    };
    let blank = stock();
    let p = GeometryPolicy::default();
    let solid = blank
        .rust_cut_line8_nut(cut, EntityId::new(), p, None)
        .unwrap();
    assert_eq!(solid.rust_line8_nut_requested_wall_mm(), Some(0.5));
    let source = solid.rust_line8_nut_source().unwrap();
    assert_eq!(
        source.thread().spec().major_radius_mm,
        profile.major_radius_mm
    );
    assert_eq!(source.thread().spec().length_mm, 4.1);
    assert!(source.minimum_thread_wall_lower_bound_mm() > 0.5);
    let q = MeshOptions {
        linear_deflection_mm: 0.02,
        angular_deflection_rad: 0.2,
        maximum_triangles: 200000,
    };
    solid.manufacturing_mesh(q).unwrap();
    assert!(
        blank
            .rust_cut_line8_nut(
                Line8NutCut {
                    minimum_wall_mm: 4.,
                    ..cut
                },
                EntityId::new(),
                p,
                None
            )
            .is_err()
    );
}

#[test]
fn actual_builder_m8_screw_mates_original_line8_stock_and_relief_without_generic_aliases() {
    use multicad_rust_kernel::{head_profile::Head, profiled_screw::Spec as ScrewSpec};
    let policy = GeometryPolicy::default();
    let pitch = 1.25_f64;
    let depth = 5. * 3_f64.sqrt() / 16. * pitch;
    let head_height = 4.4;
    let underhead = 14. - head_height;
    let runout = 1.25;
    let neck = -5. + head_height;
    let external_start = neck + runout;
    let seat = 4.5;
    let station: f64 = seat - external_start;
    let angle = std::f64::consts::TAU * station.rem_euclid(pitch) / pitch;
    let (s, c) = angle.portable_sin_cos();
    let raw = s.portable_atan2(c);
    let roll = (raw + std::f64::consts::TAU / 2.).rem_euclid(std::f64::consts::TAU)
        - std::f64::consts::TAU / 2.;
    let (s, c) = roll.portable_sin_cos();
    let external = HelicalThread {
        major_radius_mm: 4. - 0.08,
        pitch_mm: pitch,
        length_mm: underhead - runout,
        depth_mm: depth,
        crest_width_mm: pitch / 8. - 0.04,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        radial_taper: 0.,
        starts: 1,
        left_hand: false,
    };
    let screw = Solid::rust_profiled_machine_screw_with_shank(
        external,
        None,
        4.,
        ScrewSpec {
            head: Head::Countersunk {
                radius_mm: 8.,
                height_mm: head_height,
                included_angle_deg: 90.,
                edge_chamfer_mm: 0.05,
                neck_radius_mm: 0.4,
            },
            shank_length_mm: underhead - external.length_mm,
            socket_across_flats_mm: 5.08,
            socket_depth_mm: 3.2,
            minimum_wall_mm: 1e-5,
        },
        EntityId::new(),
        PlaneFrame {
            origin: [0., 0., neck],
            ..PlaneFrame::default()
        },
        policy,
    )
    .unwrap();
    let stock = stock()
        .transformed(RigidTransform {
            translation_mm: [0., 0., seat],
            ..RigidTransform::default()
        })
        .unwrap();
    let nut = stock
        .rust_cut_line8_nut(
            Line8NutCut {
                profile: HelicalThread {
                    major_radius_mm: 4. + 0.08,
                    length_mm: 4.1,
                    crest_width_mm: pitch / 8. + 0.04,
                    ..external
                },
                frame: PlaneFrame {
                    origin: [0., 0., seat],
                    u: [c, s, 0.],
                    v: [-s, c, 0.],
                },
                bore_from_start_mm: 7.5,
                relief: Relief {
                    start_mm: 4.1 - 0.001,
                    radius_mm: 4.10,
                    tool_length_mm: 7.3 - 4.1 + 0.02,
                },
                minimum_wall_mm: 0.5,
            },
            EntityId::new(),
            policy,
            None,
        )
        .unwrap();
    let proof = screw.rust_line8_nut_mating(&nut, None).unwrap();
    assert!(std::ptr::eq(
        proof.nut().source,
        nut.rust_line8_nut_source().unwrap()
    ));
    assert!(proof.complete_turns_lower > 3.);
    assert!(proof.relief_clearance_lower_mm > 0.1);
    assert!(nut.axial_cylinder().is_none());
    assert!(matches!(
        screw
            .rust_line8_nut_mating(&nut, Some(&|| true))
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Cancelled
    ));
    assert!(matches!(
        stock.rust_line8_nut_mating(&nut, None).err().unwrap().kind,
        GeometryErrorKind::Unsupported
    ));
    let offset = screw
        .transformed(RigidTransform {
            translation_mm: [0.1, 0., 0.],
            ..RigidTransform::default()
        })
        .unwrap();
    assert!(matches!(
        offset.rust_line8_nut_mating(&nut, None).err().unwrap().kind,
        GeometryErrorKind::Unsupported
    ));
}
