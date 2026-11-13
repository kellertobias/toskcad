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
    let zero = Rational::new(0, 1).unwrap();
    let stock = Solid::from_axial_cylinder(
        AxialCylinder::new(
            7,
            Circle::new([zero; 2], Rational::new(5, 2).unwrap()).unwrap(),
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
                // including literal vertices, indices and face identities.
                let expected = match (entry, tail, posed) {
                    (true, true, false) => 0xc8d6834bec289e5c,
                    (true, true, true) => 0x091fdf1ec200c558,
                    (true, false, false) => 0xfd492d49b3f3bef3,
                    (false, true, false) => 0x7574f20e3ec4298e,
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
    assert_eq!(mesh_hash(&display), 0x09b1cf46b5c1cc7b);
    let manufacturing = solid
        .manufacturing_mesh_with_receipt(MeshOptions::default())
        .unwrap();
    unique(&manufacturing.mesh);
}
#[test]
fn unresolved_nonbinary_limit_remains_atomic_and_is_not_a_geometry_waiver() {
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
    assert_eq!(
        solid.mesh(MeshOptions::default()).unwrap_err().kind,
        GeometryErrorKind::Limit
    );
    assert_eq!(
        solid
            .manufacturing_mesh_with_receipt(MeshOptions::default())
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
