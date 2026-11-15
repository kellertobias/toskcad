use super::*;
fn stock() -> Sweep {
    let ring = vec![
        [-6.9, 0.],
        [6.9, 0.],
        [6.5, 3.],
        [5.4, 5.],
        [3.3, 7.3],
        [-3.3, 7.3],
        [-5.4, 5.],
        [-6.5, 3.],
    ];
    Sweep::new(vec![ring.clone()], vec![ring], 23.).unwrap()
}
fn thread() -> AuthoredThread {
    AuthoredThread::new(
        crate::helical_thread::Spec {
            major_radius_mm: 4.08,
            pitch_mm: 1.25,
            length_mm: 4.1,
            depth_mm: 0.65,
            crest_width_mm: 0.15,
            lower_flank_deg: 30.,
            upper_flank_deg: 30.,
            starts: 1,
            left_hand: false,
            root_radius_mm: None,
        },
        crate::authored_flat_thread::SourceOptions {
            position_limit_mm: 1e-7,
            normal_angle_limit_rad: 1e-7,
        },
    )
    .unwrap()
}
fn nut() -> Nut {
    Nut::new(
        stock(),
        7.5,
        thread(),
        [[1., 0.], [0., 1.]],
        Relief {
            start_mm: 4.099,
            radius_mm: 4.1,
            tool_length_mm: 3.22,
        },
        None,
    )
    .unwrap()
}
#[test]
fn actual_stock_tap_and_open_crown_material_volume_and_whole_queries() {
    let nut = nut();
    assert_eq!(nut.stock().base(), stock().base());
    assert_eq!(nut.stock().height_mm(), 23.);
    assert_eq!(nut.thread(), &thread());
    let v = nut.volume_bounds(1e-7, 2048, None).unwrap();
    assert!(v.bounds_mm3()[1] - v.bounds_mm3()[0] < 1e-7);
    let r = 4.1_f64;
    let n = 65536;
    let dz = (7.3 - 4.099) / n as f64;
    let mut removed = 0.;
    for i in 0..n {
        let z = 4.099 + (i as f64 + 0.5) * dz;
        let w = if z < 5. {
            6.5 + (5.4 - 6.5) * (z - 3.) / 2.
        } else {
            5.4 + (3.3 - 5.4) * (z - 5.) / 2.3
        };
        removed += if w >= r {
            std::f64::consts::PI * r * r
        } else {
            2. * (w * (r * r - w * w).sqrt() + r * r * (w / r).asin())
        };
    }
    removed *= dz;
    let avg = thread()
        .bands()
        .into_iter()
        .map(|b| {
            let [a, c] = b.radius_mm;
            (b.phase_mm[1] - b.phase_mm[0]) * (a * a + a * c + c * c) / 3.
        })
        .sum::<f64>()
        / 1.25;
    let oracle = 84.01 * 23. - std::f64::consts::PI * avg * 4.099 - removed;
    let bounds = v.bounds_mm3();
    assert!((bounds[0].midpoint(bounds[1]) - oracle).abs() < 1e-6);
    let naive = 84.01 * 23.
        - std::f64::consts::PI * avg * 4.099
        - std::f64::consts::PI * r * r * (7.3 - 4.099);
    assert!(
        bounds[0] > naive + 0.1,
        "crown opening must retain the actual outer section"
    );
    for (point, expected) in [
        ([6., 0., 2.], Classification::Inside),
        ([0., 0., 2.], Classification::Outside),
        ([0., 0., 6.], Classification::Outside),
        ([0., 3., 7.2], Classification::Outside),
        ([6., 0., 7.2], Classification::Inside),
        ([0., 8., 2.], Classification::Outside),
        ([-8., 0., 2.], Classification::Outside),
    ] {
        assert_eq!(
            nut.classify_box(point.map(|v| [v - 0.001, v + 0.001]), None)
                .unwrap(),
            expected,
            "{point:?}"
        );
    }
    assert_eq!(
        nut.classify_box([[3.5, 4.1], [0., 0.], [2., 2.]], None)
            .unwrap(),
        Classification::Undetermined
    );
}
#[test]
fn retained_nut_source_refuses_invalid_overlap_budget_and_cancel() {
    let n = nut();
    assert_eq!(n.volume_bounds(1e-7, 2047, None).err(), Some(Error::Limit));
    assert_eq!(
        n.volume_bounds(1e-7, 2048, Some(&|| true)).err(),
        Some(Error::Cancelled)
    );
    assert_eq!(
        n.classify_box([[0., 0.]; 3], Some(&|| true)),
        Err(Error::Cancelled)
    );
    assert_eq!(n.classify_box([[1., 0.]; 3], None), Err(Error::Invalid));
    let mut r = n.relief();
    r.radius_mm = 3.9;
    assert!(Nut::new(stock(), 7.5, thread(), [[1., 0.], [0., 1.]], r, None).is_err());
}
#[test]
fn actual_open_crown_nut_candidate_is_one_connected_seventeen_face_shell() {
    for rolled in [false, true] {
        let mut n = nut();
        if rolled {
            let (s, c) = 0.317_f64.sin_cos();
            n = Nut::new(stock(), 7.5, thread(), [[c, -s], [s, c]], n.relief(), None).unwrap();
        }
        let q = crate::display::Tessellation {
            deflection: 0.02,
            angle: 0.2,
            triangle_limit: 200000,
        };
        let start = std::time::Instant::now();
        let mesh = n.mesh_candidate(q, 64_000_000, None).unwrap();
        eprintln!(
            "line8 rolled{rolled}: {}V/{}T {}work {}bytes {:?}",
            mesh.mesh().vertices.len() / 3,
            mesh.mesh().face_indices.len(),
            mesh.work_used(),
            mesh.retained_bytes(),
            start.elapsed()
        );
        assert_eq!(mesh.source().thread(), n.thread());
        assert!(mesh.source_distance_mm() <= q.deflection);
        assert!(mesh.normal_chord() <= q.angle);
        let faces = mesh
            .mesh()
            .face_indices
            .iter()
            .copied()
            .collect::<std::collections::BTreeSet<_>>();
        assert_eq!(faces, (0..17).collect());
        assert_eq!(
            mesh.connectivity().vertex_links as usize,
            mesh.mesh().vertices.len() / 3
        );
        let source = n.volume_bounds(1e-7, 2048, None).unwrap().bounds_mm3();
        assert!((mesh.connectivity().signed_volume_mm3 - source[0].midpoint(source[1])).abs() < 2.);
        assert!(n.mesh_candidate(q, mesh.work_used() - 1, None).is_err());
        assert_eq!(
            n.mesh_candidate(q, 64_000_000, Some(&|| true)).err(),
            Some(Error::Cancelled)
        );
    }
}
#[test]
fn all_seventeen_physical_face_measures_are_source_integrals_including_crown_notches() {
    for rolled in [false, true] {
        let mut n = nut();
        if rolled {
            let (s, c) = 0.317_f64.sin_cos();
            n = Nut::new(stock(), 7.5, thread(), [[c, -s], [s, c]], n.relief(), None).unwrap();
        }
        for face in 0..17 {
            let m = n
                .face_measure(face, Default::default(), None)
                .unwrap_or_else(|e| panic!("face{face}/rolled{rolled}: {e:?}"));
            assert!(m.area_mm2()[0] > 0.);
            assert!(m.area_mm2()[1] - m.area_mm2()[0] < 1e-7);
            assert!(m.first_moments_mm3().iter().all(|v| v[1] - v[0] < 1e-7));
            assert!(m.centroid_mm().iter().all(|v| v[1] - v[0] < 1e-7));
            if [4, 10].contains(&face) {
                assert!(m.centroid_mm()[2][0] <= 7.3 && m.centroid_mm()[2][1] >= 7.3)
            }
            if face == 15 {
                assert!(m.centroid_mm()[2][0] <= 4.099 && m.centroid_mm()[2][1] >= 4.099)
            }
        }
        let left = n.face_measure(4, Default::default(), None).unwrap();
        let right = n.face_measure(10, Default::default(), None).unwrap();
        assert!(right.area_mm2()[0] > left.area_mm2()[1]);
        let shoulder = n.face_measure(15, Default::default(), None).unwrap();
        assert!(shoulder.area_mm2()[0] > 0.1);
        assert!(n.face_measure(17, Default::default(), None).is_err());
        assert_eq!(
            n.face_measure(0, Default::default(), Some(&|| true)).err(),
            Some(Error::Cancelled)
        );
    }
}
#[test]
fn original_finite_topology_and_every_coedge_keep_exact_owners_for_both_hands() {
    for left in [false, true] {
        let mut spec = thread().spec();
        spec.left_hand = left;
        let thread = AuthoredThread::new(
            spec,
            crate::authored_flat_thread::SourceOptions {
                position_limit_mm: 1e-7,
                normal_angle_limit_rad: 1e-7,
            },
        )
        .unwrap();
        let n = Nut::new(
            stock(),
            7.5,
            thread,
            [[1., 0.], [0., 1.]],
            nut().relief(),
            None,
        )
        .unwrap();
        let boundary = n.boundary(131072, None).unwrap();
        assert_eq!(boundary.faces().len(), 17);
        assert_eq!(boundary.edges().len(), 46);
        assert_eq!(boundary.vertices().len(), 32);
        assert!(boundary.retained_bytes() < 65536);
        for &edge in boundary.edges() {
            let whole = n
                .curve_bounds(edge, [0., 1.], 2048, None)
                .unwrap()
                .position_mm();
            assert!(
                n.curve_second_derivative_bound_mm(edge)
                    .unwrap()
                    .is_finite()
            );
            for (i, vertex) in edge.endpoints().unwrap().into_iter().enumerate() {
                let t = i as f64;
                let endpoint = n
                    .curve_bounds(edge, [t, t], 2048, None)
                    .unwrap()
                    .position_mm();
                let owner = n.vertex_bounds(vertex, 2048, None).unwrap().position_mm();
                assert_eq!(endpoint, owner);
                for axis in 0..3 {
                    assert!(
                        whole[axis][0] <= owner[axis][1] && whole[axis][1] >= owner[axis][0],
                        "{edge:?}/axis{axis}: {:?}/{:?}",
                        whole[axis],
                        owner[axis]
                    );
                }
            }
        }
        assert_eq!(n.boundary(131071, None).err(), Some(Error::Limit));
        assert_eq!(
            n.boundary(131072, Some(&|| true)).err(),
            Some(Error::Cancelled)
        );
        assert_eq!(
            n.vertex_bounds(Vertex::Crown(4), 2048, None).err(),
            Some(Error::Invalid)
        );
        assert_eq!(
            n.curve_bounds(Edge::Rail(4), [0., 1.], 2048, None).err(),
            Some(Error::Invalid)
        );
        assert_eq!(
            n.curve_bounds(Edge::Crown(0), [1., 0.], 2048, None).err(),
            Some(Error::Invalid)
        );
    }
}
#[test]
fn actual_full_prefix_stock_plane_wall_bound_is_positive_and_conservative() {
    let n = nut();
    let wall = n.minimum_thread_wall_lower_bound_mm();
    assert!(wall > 1.5 && wall < 1.7);
    let r = thread().maximum_radius_mm();
    for z in [0., 1., 3., 4.099] {
        for y in [-r, r] {
            for side in [1usize, 2, 3, 5, 6, 7] {
                let a = n.stock().base()[0][side];
                let b = n.stock().base()[0][(side + 1) % 8];
                let dy = b[0] - a[0];
                let dz = b[1] - a[1];
                let distance = (dy * (z - a[1]) - dz * (y - a[0])) / dy.hypot(dz);
                assert!(distance >= wall);
            }
        }
    }
}
#[test]
fn actual_resolved_metric_m8_with_axial_and_radial_clearance_keeps_rolled_shell() {
    let mut spec = thread().spec();
    spec.depth_mm = 0.6765823467065926;
    spec.crest_width_mm = 0.19625;
    let source = AuthoredThread::new(
        spec,
        crate::authored_flat_thread::SourceOptions {
            position_limit_mm: 1e-7,
            normal_angle_limit_rad: 1e-7,
        },
    )
    .unwrap();
    let n = Nut::new(
        stock(),
        7.5,
        source,
        [
            [0.9501748447925625, -0.31171744308496674],
            [0.31171744308496674, 0.9501748447925625],
        ],
        nut().relief(),
        None,
    )
    .unwrap();
    n.mesh_candidate(
        crate::display::Tessellation {
            deflection: 0.02,
            angle: 0.2,
            triangle_limit: 200000,
        },
        64_000_000,
        None,
    )
    .unwrap();
}
