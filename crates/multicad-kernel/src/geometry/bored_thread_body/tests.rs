use super::*;
use multicad_rust_kernel::{Circle, Rational, cylinder::AxialCylinder, manufacturing::certify};

fn drilled(policy: GeometryPolicy) -> Solid {
    let n = |value: i64| Rational::new(value, 1).unwrap();
    let mut stock = Solid::rust_box(
        EntityId::new().uuid_bits(),
        [130, 20, 4].map(n),
        PlaneFrame::default(),
        policy,
    )
    .unwrap();
    for x in [20, 52] {
        let tool = EntityId::new();
        let drill = Solid::from_axial_cylinder_sweep(
            AxialCylinder::new(
                tool.uuid_bits(),
                Circle::new([n(x), n(10)], n(3)).unwrap(),
                n(-1),
                n(6),
            )
            .unwrap(),
            PlaneFrame {
                origin: [x as f64, 10., -1.],
                ..Default::default()
            },
            [0.; 2],
            false,
            Origin {
                source: tool,
                role: "primitive.circle".into(),
            },
            policy,
        )
        .unwrap();
        stock = stock
            .boolean(&drill, BooleanOperation::Subtract, policy)
            .unwrap()
            .unwrap();
    }
    stock
}

fn m8(length_mm: f64) -> HelicalThread {
    HelicalThread {
        major_radius_mm: 4.,
        pitch_mm: 1.25,
        length_mm,
        depth_mm: 0.6,
        crest_width_mm: 0.15,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        radial_taper: 0.,
        starts: 1,
        left_hand: false,
    }
}

fn at(origin: [f64; 3], up: bool) -> PlaneFrame {
    PlaneFrame {
        origin,
        u: [0.6, 0.8, 0.],
        v: if up { [-0.8, 0.6, 0.] } else { [0.8, -0.6, 0.] },
    }
}

#[test]
fn buried_thread_splits_only_its_coaxial_bore_into_collars_steps_and_bands() {
    let policy = GeometryPolicy::default();
    let stock = drilled(policy);
    let feature = EntityId::new();
    for (frame, low) in [
        (at([20., 10., 0.5], true), 0.5),
        (at([20., 10., 3.5], false), 0.5),
    ] {
        let solid = stock
            .rust_cut_internal_thread_in_bore(m8(3.), None, 0.5, feature, frame, policy)
            .unwrap();
        let body = solid.bored_thread_body().unwrap();
        let Wall::Thread(thread) = body.wall else {
            panic!("modeled band")
        };
        assert!(thread.root_radius_mm() > 3.);

        // Exact material: drilled stock minus the tool outside the bore.
        let volume = solid.volume_mm3().unwrap();
        let expected = stock.volume_mm3().unwrap() - (thread.volume_mm3() - PI * 9. * 3.);
        assert!(
            (volume - expected).abs() <= 1e-9 * expected,
            "{volume} {expected}"
        );

        // One closed, certified shell using every published face.
        let options = multicad_rust_kernel::display::Tessellation {
            deflection: 0.01,
            angle: 0.5,
            triangle_limit: 2_000_000,
        };
        let mesh = body.mesh(options).unwrap();
        let local = multicad_rust_kernel::display::DisplayMesh {
            vertices: mesh.vertices.clone(),
            indices: mesh.indices.clone(),
            face_indices: mesh.face_indices.clone(),
        };
        let proof = certify(&local).unwrap();
        assert!(
            (proof.signed_volume_mm3 - volume).abs() < 0.05 * volume.portable_cbrt(),
            "{} {volume}",
            proof.signed_volume_mm3
        );
        let faces = solid.faces().unwrap();
        let used: BTreeSet<u32> = mesh.face_indices.iter().copied().collect();
        assert_eq!(used, (0..faces.len() as u32).collect());
        let area: f64 = faces.iter().map(|f| f.area_mm2).sum();
        assert!(area > 0. && faces.iter().all(|f| f.area_mm2 > 0.));

        // Topology: the original stock references stay; the bore wall splits.
        let topology = solid.topology_snapshot(feature, 1).unwrap();
        let original = stock.topology_snapshot(feature, 1).unwrap();
        let roles = |t: &TopologySnapshot, kind| -> BTreeSet<String> {
            t.elements
                .iter()
                .filter(|e| e.key.kind == kind)
                .flat_map(|e| e.origins.iter().map(|o| o.role.clone()))
                .collect()
        };
        let face_roles = roles(&topology, TopologyKind::Face);
        for role in roles(&original, TopologyKind::Face) {
            assert!(
                face_roles.contains(&role),
                "{role} {:?} {:?}",
                roles(&original, TopologyKind::Face),
                face_roles
            );
        }
        for role in [
            "thread.bore/collar.mouth",
            "thread.bore/collar.exit",
            "thread.bore/step.start",
            "thread.bore/step.end",
        ] {
            assert!(face_roles.contains(role), "{role}");
        }
        assert!(face_roles.iter().any(|r| r.starts_with("thread/")));
        let ordinals: Vec<u32> = topology
            .elements
            .iter()
            .filter(|e| e.key.kind == TopologyKind::Face)
            .map(|e| e.key.ordinal)
            .collect();
        assert_eq!(
            ordinals.iter().copied().collect::<BTreeSet<_>>(),
            (0..faces.len() as u32).collect()
        );
        assert_eq!(ordinals.len(), faces.len());
        // Both untouched mouths of the threaded bore and its collar rims.
        let circles: Vec<[f64; 4]> = topology
            .elements
            .iter()
            .filter_map(|e| match e.geometry {
                Some(AnalyticGeometry::CircularArc {
                    frame, radius_mm, ..
                }) if e.key.kind == TopologyKind::Edge => Some([
                    frame.origin_mm[0],
                    frame.origin_mm[1],
                    frame.origin_mm[2],
                    radius_mm,
                ]),
                _ => None,
            })
            .collect();
        for z in [0., 4., low, low + 3.] {
            assert!(
                circles.iter().any(|c| (c[0] - 20.).abs() < 1e-9
                    && (c[1] - 10.).abs() < 1e-9
                    && (c[2] - z).abs() < 1e-9
                    && (c[3] - 3.).abs() < 1e-9),
                "{z} {circles:?}"
            );
        }
        let bounds = solid.bounds().unwrap();
        assert_eq!(bounds.min, [0.; 3]);
        assert_eq!(bounds.max, [130., 20., 4.]);
        let moved = solid
            .transformed(RigidTransform {
                rotation: [[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]],
                translation_mm: [3., -4., 5.],
            })
            .unwrap();
        assert!((moved.volume_mm3().unwrap() - volume).abs() <= 1e-9 * volume);
        assert_eq!(moved.faces().unwrap().len(), faces.len());
        moved.topology_snapshot(feature, 2).unwrap();
    }
}

#[test]
fn buried_thread_refuses_open_mouths_offsets_wide_bores_and_thin_walls() {
    let policy = GeometryPolicy::default();
    let stock = drilled(policy);
    let cut = |spec, wall, frame| {
        stock
            .rust_cut_internal_thread_in_bore(spec, None, wall, EntityId::new(), frame, policy)
            .err()
            .unwrap()
            .kind
    };
    // Reaching a bore mouth, off the bore axis, or no bore at all.
    for frame in [
        at([20., 10., 0.], true),
        at([20., 10., 1.5], true),
        at([20.5, 10., 0.5], true),
        at([90., 10., 0.5], true),
    ] {
        assert_eq!(cut(m8(3.), 0.5, frame), GeometryErrorKind::Unsupported);
    }
    // The existing bore does not lie inside the root of a smaller thread.
    let small = HelicalThread {
        major_radius_mm: 3.5,
        ..m8(3.)
    };
    assert_eq!(
        cut(small, 0.5, at([20., 10., 0.5], true)),
        GeometryErrorKind::Unsupported
    );
    // The protected collar would leave the stock or reach the other bore.
    for wall in [6.5, 25.] {
        assert_eq!(
            cut(m8(3.), wall, at([20., 10., 0.5], true)),
            GeometryErrorKind::Input
        );
    }
}

#[test]
fn cosmetic_buried_band_is_a_smooth_crest_bore_without_manufacturing() {
    let policy = GeometryPolicy::default();
    let stock = drilled(policy);
    let feature = EntityId::new();
    let frame = at([20., 10., 0.5], true);
    let solid = stock
        .rust_cut_cosmetic_internal_thread_in_bore(3.4, 4., 3., 0.5, feature, frame, policy)
        .unwrap();
    let body = solid.bored_thread_body().unwrap();
    assert!(body.is_cosmetic());
    let volume = solid.volume_mm3().unwrap();
    let expected = stock.volume_mm3().unwrap() - PI * (3.4 * 3.4 - 9.) * 3.;
    assert!(
        (volume - expected).abs() <= 1e-9 * expected,
        "{volume} {expected}"
    );
    let options = multicad_rust_kernel::display::Tessellation {
        deflection: 0.01,
        angle: 0.5,
        triangle_limit: 2_000_000,
    };
    let mesh = body.mesh(options).unwrap();
    let local = multicad_rust_kernel::display::DisplayMesh {
        vertices: mesh.vertices.clone(),
        indices: mesh.indices.clone(),
        face_indices: mesh.face_indices.clone(),
    };
    certify(&local).unwrap();
    let faces = solid.faces().unwrap();
    let used: BTreeSet<u32> = mesh.face_indices.iter().copied().collect();
    assert_eq!(used, (0..faces.len() as u32).collect());
    let topology = solid.topology_snapshot(feature, 1).unwrap();
    assert!(topology.elements.iter().any(|e| {
        e.origins
            .iter()
            .any(|o| o.role == "thread.cosmetic.envelope")
    }));
    assert_eq!(
        solid
            .manufacturing_mesh(Default::default())
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Unsupported
    );
    // The crest must clear the bore and stay inside the major radius.
    for (crest, major) in [(3., 4.), (4., 4.)] {
        assert!(
            stock
                .rust_cut_cosmetic_internal_thread_in_bore(
                    crest, major, 3., 0.5, feature, frame, policy
                )
                .is_err()
        );
    }
}
