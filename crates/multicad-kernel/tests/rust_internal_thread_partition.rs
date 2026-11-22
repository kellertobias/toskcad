//! Normal source factories and literal manufacturing/embedding, never injected meshes.
use multicad_core::{
    commands::Session,
    document::{Document, EntityId},
    model::{ThreadFit, ThreadHand, ThreadMode, ThreadParameters, length},
};
use multicad_kernel::{
    GeometryErrorKind, GeometryPolicy, HelicalThread, InternalThreadTreatmentCut, MeshOptions,
    PlaneFrame, RigidTransform, Solid,
    rust_kernel::{
        Circle, Rational,
        authored_internal_treatment::{Intent, Terminal},
        cylinder::AxialCylinder,
    },
};
fn profile(fit: ThreadFit, left: bool, starts: u32) -> HelicalThread {
    let session = Session::new(Document::new("Exact mesh regression").unwrap()).unwrap();
    let p = ThreadParameters {
        major_diameter: length(3.),
        pitch: length(0.5),
        length: length(4.),
        fit,
        hand: if left {
            ThreadHand::Left
        } else {
            ThreadHand::Right
        },
        starts,
        ..Default::default()
    };
    let r = p
        .resolve(session.parameters(), ThreadMode::CutInternal)
        .unwrap();
    HelicalThread {
        major_radius_mm: r.major_radius_mm,
        pitch_mm: r.pitch_mm,
        length_mm: r.length_mm,
        depth_mm: r.profile.depth_mm,
        crest_width_mm: r.profile.crest_width_mm,
        lower_flank_deg: r.profile.lower_flank_deg,
        upper_flank_deg: r.profile.upper_flank_deg,
        radial_taper: r.profile.radial_taper,
        starts: r.starts,
        left_hand: left,
    }
}
fn body(profile: HelicalThread, entry: bool, tail: bool) -> Solid {
    body_radius(profile, entry, tail, Rational::new(5, 2).unwrap())
}
fn body_radius(profile: HelicalThread, entry: bool, tail: bool, radius: Rational) -> Solid {
    let zero = Rational::new(0, 1).unwrap();
    let stock = Solid::from_axial_cylinder(
        AxialCylinder::new(
            7,
            Circle::new([zero; 2], radius).unwrap(),
            zero,
            Rational::new(4, 1).unwrap(),
        )
        .unwrap(),
        PlaneFrame::default(),
        GeometryPolicy::default(),
    )
    .unwrap();
    stock
        .rust_cut_authored_internal_thread_treatments_with_interrupt(
            InternalThreadTreatmentCut {
                profile,
                intent: Intent {
                    entry_depth_mm: entry.then_some(0.5),
                    terminal: tail.then_some(Terminal {
                        width_mm: 0.5,
                        radial_depth_mm: 0.125,
                    }),
                },
                minimum_wall_mm: 0.5,
                frame: PlaneFrame::default(),
            },
            EntityId::new(),
            GeometryPolicy::default(),
            None,
        )
        .unwrap()
}
fn mesh_hash(mesh: &multicad_kernel::Mesh) -> u64 {
    let mut h = 0xcbf29ce484222325_u64;
    for b in mesh
        .vertices
        .iter()
        .flat_map(|v| v.to_bits().to_le_bytes())
        .chain(mesh.indices.iter().flat_map(|v| v.to_le_bytes()))
        .chain(mesh.face_indices.iter().flat_map(|v| v.to_le_bytes()))
    {
        h = (h ^ u64::from(b)).wrapping_mul(0x100000001b3);
    }
    h
}
fn unique(mesh: &multicad_kernel::Mesh) {
    let mut seen = std::collections::BTreeSet::new();
    for p in mesh.vertices.as_chunks::<3>().0 {
        assert!(
            seen.insert(p.map(|v| if v == 0. { 0 } else { v.to_bits() })),
            "duplicate literal emitted point {p:?}"
        );
    }
}

#[test]
fn nearby_positive_terminal_bands_keep_old_exact_binary_domains_and_fresh_world_embedding() {
    let options = MeshOptions {
        linear_deflection_mm: 0.02,
        angular_deflection_rad: 0.3,
        maximum_triangles: 500_000,
    };
    let mut failures = Vec::new();
    for steps in [-4_i32, -2, 0, 1, 2, 4] {
        let mut length = 3.3_f64;
        for _ in 0..steps.unsigned_abs() {
            length = if steps < 0 {
                length.next_down()
            } else {
                length.next_up()
            };
        }
        // Every fixture is in[2,4), hence its exact binary significand has
        // denominator2^51. This is actual stock data, not a rounded decimal.
        assert!((2. ..4.).contains(&length));
        let end = Rational::new(
            ((length.to_bits() & ((1_u64 << 52) - 1)) | (1_u64 << 52)) as i64,
            1_i64 << 51,
        )
        .unwrap();
        assert_eq!(end.as_f64(), length);
        for left in [false, true] {
            let p = HelicalThread {
                length_mm: length,
                ..profile(ThreadFit::Metric6HInternal, left, 1)
            };
            let zero = Rational::new(0, 1).unwrap();
            let stock = Solid::from_axial_cylinder(
                AxialCylinder::new(
                    7,
                    Circle::new([zero; 2], Rational::new(7, 2).unwrap()).unwrap(),
                    zero,
                    end,
                )
                .unwrap(),
                PlaneFrame::default(),
                GeometryPolicy::default(),
            )
            .unwrap();
            let solid = stock
                .rust_cut_authored_internal_thread_treatments_with_interrupt(
                    InternalThreadTreatmentCut {
                        profile: p,
                        intent: Intent {
                            entry_depth_mm: Some(0.4),
                            terminal: Some(Terminal {
                                width_mm: 0.3,
                                radial_depth_mm: 0.05,
                            }),
                        },
                        minimum_wall_mm: 0.5,
                        frame: PlaneFrame::default(),
                    },
                    EntityId::new(),
                    GeometryPolicy::default(),
                    None,
                )
                .unwrap();
            let source = solid.rust_internal_thread_treatment_data().unwrap();
            assert_eq!(source.material.cavity().exact_receiving_extent(), None);
            assert_eq!(source.material.cavity().full_length_mm(), length);
            assert_eq!(
                source.material.cavity().terminal_station_mm(),
                Some(length - 0.3)
            );
            assert_eq!(source.material.stock().height(), end);
            let moved = solid
                .transformed(RigidTransform {
                    rotation: [[0.6, 0., 0.8], [0., 1., 0.], [-0.8, 0., 0.6]],
                    translation_mm: [11., -7., 5.],
                })
                .unwrap();
            for (posed, body) in [(false, &solid), (true, &moved)] {
                match body.manufacturing_mesh_with_receipt(options) {
                    Ok(receipt) => {
                        unique(&receipt.mesh);
                        assert!(receipt.edge_pairs > 0 && receipt.vertex_links > 0);
                        eprintln!(
                            "terminalNeighbor ulps={steps} left={left} posed={posed} facets={}",
                            receipt.mesh.face_indices.len()
                        );
                    }
                    Err(e) => {
                        failures.push(format!("ulps={steps} left={left} posed={posed}: {e:?}"))
                    }
                }
            }
            assert_eq!(solid.rust_internal_thread_treatment_data().unwrap(), source);
        }
    }
    assert!(
        failures.is_empty(),
        "Actual neighboring terminal sources: {failures:#?}"
    );
}
#[test]
fn default_native_nominal_and_six_h_meshes_keep_manufacturing_embedding_and_original_successes() {
    for fit in [ThreadFit::NominalBasic, ThreadFit::Metric6HInternal] {
        for (left, starts, entry, tail, posed) in [
            (false, 1, true, true, false),
            (true, 2, true, true, true),
            (false, 1, true, false, false),
            (true, 1, false, true, false),
        ] {
            // Sourced6H admits one start; nominal actual source also tests two.
            let starts = if matches!(fit, ThreadFit::Metric6HInternal) {
                1
            } else {
                starts
            };
            let original = body(profile(fit, left, starts), entry, tail);
            let solid = if posed {
                original
                    .transformed(RigidTransform {
                        rotation: [[0.8, -0.36, 0.48], [0.6, 0.48, -0.64], [0., 0.8, 0.6]],
                        translation_mm: [10., 20., 30.],
                    })
                    .unwrap()
            } else {
                original
            };
            let data = solid.rust_internal_thread_treatment_data().unwrap();
            let display = solid.mesh(MeshOptions::default()).unwrap();
            unique(&display);
            let manufacturing = solid
                .manufacturing_mesh_with_receipt(MeshOptions::default())
                .unwrap();
            unique(&manufacturing.mesh);
            assert!(manufacturing.edge_pairs > 0 && manufacturing.vertex_links > 0);
            assert_eq!(solid.rust_internal_thread_treatment_data().unwrap(), data);
            if matches!(fit, ThreadFit::Metric6HInternal) {
                // Captured from the original successful64-row implementation,
                // including literal vertices, indices and face identities;
                // re-recorded under portable maths (identical topology,
                // vertices within 4e-15 mm of the host-libm capture).
                let expected = match (entry, tail, posed) {
                    (true, true, false) => 0xcce4d796eeb63a86,
                    (true, true, true) => 0x00e6800d93532b9b,
                    (true, false, false) => 0x4f52b26f299de29f,
                    (false, true, false) => 0x6c7cff25535a9fef,
                    _ => unreachable!(),
                };
                assert_eq!(mesh_hash(&display), expected);
            }
        }
    }
    let p = HelicalThread {
        major_radius_mm: 1.5,
        pitch_mm: 0.3,
        length_mm: 4.,
        depth_mm: 0.15,
        crest_width_mm: 0.08,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        radial_taper: 0.,
        left_hand: false,
        starts: 1,
    };
    let solid = body(p, true, true);
    let display = solid.mesh(MeshOptions::default()).unwrap();
    unique(&display);
    assert_eq!(mesh_hash(&display), 0x021a74c4470c2556);
    let manufacturing = solid
        .manufacturing_mesh_with_receipt(MeshOptions::default())
        .unwrap();
    unique(&manufacturing.mesh);
}
#[test]
fn nonbinary_source_manufactures_and_real_budget_refusal_remains_atomic() {
    let p = HelicalThread {
        major_radius_mm: 1.5,
        pitch_mm: 0.35,
        length_mm: 4.,
        depth_mm: 0.15,
        crest_width_mm: 0.08,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        radial_taper: 0.,
        left_hand: true,
        starts: 3,
    };
    let solid = body(p, true, true);
    let original = solid.rust_internal_thread_treatment_data().unwrap();
    let volume = solid
        .rust_internal_thread_treatment_volume_bounds()
        .unwrap()
        .unwrap();
    let display = solid.mesh(MeshOptions::default()).unwrap();
    unique(&display);
    let manufacturing = solid
        .manufacturing_mesh_with_receipt(MeshOptions::default())
        .unwrap();
    unique(&manufacturing.mesh);
    assert!(!manufacturing.mesh.indices.is_empty());
    let limited = MeshOptions {
        maximum_triangles: 32,
        ..MeshOptions::default()
    };
    assert_eq!(
        solid.mesh(limited).unwrap_err().kind,
        GeometryErrorKind::Limit
    );
    assert_eq!(
        solid
            .manufacturing_mesh_with_receipt(limited)
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Limit
    );
    assert_eq!(
        solid.rust_internal_thread_treatment_data().unwrap(),
        original
    );
    assert_eq!(
        solid
            .rust_internal_thread_treatment_volume_bounds()
            .unwrap()
            .unwrap(),
        volume
    );
}

#[test]
fn actual_wider_stock_keeps_inner_source_and_passes_literal_manufacturing_embedding() {
    for (left, starts, posed) in [(false, 1, false), (true, 2, true)] {
        let p = HelicalThread {
            depth_mm: 0.25,
            crest_width_mm: 0.125,
            ..profile(ThreadFit::NominalBasic, left, starts)
        };
        let mut solid = body_radius(p, true, true, Rational::new(8, 1).unwrap());
        if posed {
            solid = solid
                .transformed(RigidTransform {
                    rotation: [[0.8, -0.36, 0.48], [0.6, 0.48, -0.64], [0., 0.8, 0.6]],
                    translation_mm: [10., 20., 30.],
                })
                .unwrap();
        }
        let original = solid.rust_internal_thread_treatment_data().unwrap();
        let display = solid.mesh(MeshOptions::default()).unwrap();
        unique(&display);
        let manufacturing = solid
            .manufacturing_mesh_with_receipt(MeshOptions::default())
            .unwrap();
        assert!(!manufacturing.mesh.indices.is_empty());
        assert_eq!(
            solid.rust_internal_thread_treatment_data().unwrap(),
            original
        );
    }
}

#[test]
fn deep_terminal_cap_candidate_is_a_closed_literally_embedded_cavity() {
    use multicad_kernel::rust_kernel::{
        authored_flat_thread::{AuthoredThread, SourceOptions},
        authored_internal_treatment::{Cavity, tessellate_candidate},
        display::Tessellation,
        helical_thread::Spec,
    };
    let thread = AuthoredThread::new(
        Spec {
            major_radius_mm: 1.5,
            pitch_mm: 0.5,
            length_mm: 4.,
            depth_mm: 0.25,
            crest_width_mm: 0.125,
            lower_flank_deg: 30.,
            upper_flank_deg: 30.,
            starts: 1,
            left_hand: false,
            root_radius_mm: None,
        },
        SourceOptions {
            position_limit_mm: 1e-7,
            normal_angle_limit_rad: 1e-7,
        },
    )
    .unwrap();
    let cavity = Cavity::new(
        thread,
        Intent {
            entry_depth_mm: Some(0.5),
            terminal: Some(Terminal {
                width_mm: 0.5,
                radial_depth_mm: 100.,
            }),
        },
    )
    .unwrap();
    let options = Tessellation {
        deflection: 0.1,
        angle: 0.3,
        triangle_limit: 1_000_000,
    };
    let candidate = tessellate_candidate(&cavity, options, None).unwrap();
    assert!(candidate.work_used <= 16_000_000);
    let mesh = multicad_kernel::Mesh {
        vertices: candidate.mesh.vertices,
        indices: candidate.mesh.indices,
        face_indices: candidate.mesh.face_indices,
    };
    unique(&mesh);
    // Closed and consistently oriented: every directed edge occurs once and
    // its reverse occurs once.
    let mut edges = std::collections::HashMap::new();
    for t in mesh.indices.as_chunks::<3>().0 {
        for k in 0..3 {
            *edges.entry((t[k], t[(k + 1) % 3])).or_insert(0usize) += 1;
        }
    }
    for (&(a, b), &n) in &edges {
        assert_eq!(n, 1);
        assert_eq!(edges.get(&(b, a)), Some(&1));
    }
    // The enclosed volume is the terminal relief disc plus the threaded bore,
    // bounded by the two full cylinders of the major and terminal radii.
    let p = |i: u32| {
        let i = 3 * i as usize;
        [mesh.vertices[i], mesh.vertices[i + 1], mesh.vertices[i + 2]]
    };
    let volume = mesh
        .indices
        .as_chunks::<3>()
        .0
        .iter()
        .map(|t| {
            let [a, b, c] = [p(t[0]), p(t[1]), p(t[2])];
            a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
                + a[2] * (b[0] * c[1] - b[1] * c[0])
        })
        .sum::<f64>()
        / 6.;
    let r = cavity.terminal_radius_mm().unwrap();
    let pi = std::f64::consts::PI;
    let disc = pi * r * r * 0.5;
    assert!(volume > disc && volume < disc + pi * 1.5 * 1.5 * cavity.full_length_mm());
    let proof = mesh.prove_literal_embedding(64_000_000, None).unwrap();
    assert!(proof.pairs() > 0);
}
