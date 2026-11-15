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
        assert_eq!(n.phase_rotation(), [[1., 0.], [0., 1.]]);
        assert_eq!(n.bore_from_start_mm(), 7.5);
        let boundary = n.boundary(131072, None).unwrap();
        assert_eq!(boundary.faces().len(), 17);
        assert_eq!(boundary.edges().len(), 46);
        assert_eq!(boundary.vertices().len(), 32);
        assert!(boundary.retained_bytes() < 65536);
        assert_eq!(boundary.work_used(), 131072);
        let mut uses = std::collections::BTreeMap::<Edge, Vec<bool>>::new();
        let mut face_euler = 0isize;
        for (index, face) in boundary.faces().iter().enumerate() {
            assert_eq!(face.index(), index as u32);
            face_euler += 2 - face.loops().len() as isize;
            let measure = n
                .face_measure(face.index(), Default::default(), None)
                .unwrap();
            assert_eq!(measure.face(), face.index());
            for ring in face.loops() {
                assert!(!ring.is_empty());
                for (i, coedge) in ring.iter().enumerate() {
                    let mut end = coedge.edge.endpoints().unwrap();
                    if !coedge.forward {
                        end.reverse();
                    }
                    let next = ring[(i + 1) % ring.len()];
                    let mut start = next.edge.endpoints().unwrap();
                    if !next.forward {
                        start.reverse();
                    }
                    assert_eq!(
                        end[1], start[0],
                        "face {index} has a disconnected finite wire"
                    );
                    uses.entry(coedge.edge).or_default().push(coedge.forward);
                }
            }
        }
        assert_eq!(uses.len(), boundary.edges().len());
        for orientations in uses.values_mut() {
            orientations.sort_unstable();
            assert_eq!(
                orientations,
                &[false, true],
                "closed shell needs two opposite coedges"
            );
        }
        assert_eq!(
            boundary.vertices().len() as isize - boundary.edges().len() as isize + face_euler,
            0
        );
        for &edge in boundary.edges() {
            let curve = n.curve_bounds(edge, [0., 1.], 2048, None).unwrap();
            assert_eq!(curve.work_used(), 2048);
            let whole = curve.position_mm();
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

#[test]
fn public_source_rejects_real_wrong_stock_phase_and_cut_geometry_atomically() {
    let original = nut();
    let build = |stock, offset, thread, phase, relief, interrupt| {
        Nut::new(stock, offset, thread, phase, relief, interrupt).err()
    };
    let identity = [[1., 0.], [0., 1.]];
    for offset in [f64::NAN, f64::INFINITY] {
        assert_eq!(
            build(stock(), offset, thread(), identity, original.relief(), None),
            Some(Error::Invalid)
        );
    }
    for phase in [
        [[1., f64::NAN], [0., 1.]],
        [[1., 0.], [0., -1.]],
        [[1.0001, 0.], [0., 1.]],
    ] {
        let expected = if phase[0][1].is_nan() {
            Error::Invalid
        } else {
            Error::Unsupported
        };
        assert_eq!(
            build(stock(), 7.5, thread(), phase, original.relief(), None),
            Some(expected)
        );
    }
    let mut spec = thread().spec();
    spec.starts = 2;
    let multiple = AuthoredThread::new(
        spec,
        crate::authored_flat_thread::SourceOptions {
            position_limit_mm: 1e-7,
            normal_angle_limit_rad: 1e-7,
        },
    )
    .unwrap();
    assert_eq!(
        build(stock(), 7.5, multiple, identity, original.relief(), None),
        Some(Error::Unsupported)
    );
    let square = vec![[-7., 0.], [7., 0.], [7., 8.], [-7., 8.]];
    let rectangular = Sweep::new(vec![square.clone()], vec![square], 23.).unwrap();
    assert_eq!(
        build(
            rectangular,
            7.5,
            thread(),
            identity,
            original.relief(),
            None
        ),
        Some(Error::Unsupported)
    );
    let mut top = stock().top().to_vec();
    for point in &mut top[0] {
        point[0] += 0.1;
    }
    let tapered = Sweep::new(stock().base().to_vec(), top, 23.).unwrap();
    assert_eq!(
        build(tapered, 7.5, thread(), identity, original.relief(), None),
        Some(Error::Unsupported)
    );
    let mut concave = stock().base().to_vec();
    concave[0][3][0] = 6.4;
    concave[0][6][0] = -6.4;
    let concave = Sweep::new(concave.clone(), concave, 23.).unwrap();
    assert_eq!(
        build(concave, 7.5, thread(), identity, original.relief(), None),
        Some(Error::Unsupported)
    );
    for relief in [
        Relief {
            start_mm: 0.,
            ..original.relief()
        },
        Relief {
            radius_mm: f64::INFINITY,
            ..original.relief()
        },
        Relief {
            tool_length_mm: -1.,
            ..original.relief()
        },
    ] {
        assert_eq!(
            build(stock(), 7.5, thread(), identity, relief, None),
            Some(Error::Invalid)
        );
    }
    let mut too_large = original.relief();
    too_large.tool_length_mm = 1e8;
    assert_eq!(
        build(stock(), 7.5, thread(), identity, too_large, None),
        Some(Error::Limit)
    );
    for offset in [4., 19.] {
        assert_eq!(
            build(stock(), offset, thread(), identity, original.relief(), None),
            Some(Error::Unsupported)
        );
    }
    let polls = std::cell::Cell::new(0);
    let late = || {
        polls.set(polls.get() + 1);
        polls.get() == 2
    };
    assert_eq!(
        build(
            stock(),
            7.5,
            thread(),
            identity,
            original.relief(),
            Some(&late)
        ),
        Some(Error::Cancelled)
    );
    assert_eq!(polls.get(), 2);
    assert_eq!(nut().bounds_mm(), original.bounds_mm());
    assert_eq!(
        nut().volume_bounds(1e-7, 2048, None).unwrap().bounds_mm3(),
        original
            .volume_bounds(1e-7, 2048, None)
            .unwrap()
            .bounds_mm3()
    );
}

#[test]
fn public_material_queries_reject_nonfinite_precision_and_overflowing_boxes() {
    let source = nut();
    let before = source.volume_bounds(1e-7, 2048, None).unwrap().bounds_mm3();
    for width in [0., -1., f64::NAN, f64::INFINITY] {
        assert_eq!(
            source.volume_bounds(width, 2048, None).err(),
            Some(Error::Invalid)
        );
    }
    assert_eq!(
        source.classify_box([[f64::NEG_INFINITY, f64::INFINITY]; 3], None),
        Err(Error::Invalid)
    );
    assert_eq!(
        source.classify_box([[-f64::MAX, f64::MAX]; 3], None),
        Err(Error::Limit)
    );
    assert_eq!(
        source.volume_bounds(1e-30, 2048, None).err(),
        Some(Error::Unresolved)
    );
    let polls = std::cell::Cell::new(0);
    let late = || {
        polls.set(polls.get() + 1);
        polls.get() == 2
    };
    assert_eq!(
        source.volume_bounds(1e-7, 2048, Some(&late)).err(),
        Some(Error::Cancelled)
    );
    assert_eq!(polls.get(), 2);
    let receipt = source.volume_bounds(1e-7, 2048, None).unwrap();
    assert_eq!(receipt.work_used(), 2048);
    assert_eq!(receipt.bounds_mm3(), before);
}

#[test]
fn finite_owner_queries_enforce_ranges_budgets_and_late_cancel_without_publication() {
    let source = nut();
    let edge = Edge::Crown(0);
    let vertex = Vertex::Crown(0);
    assert_eq!(
        source.vertex_bounds(vertex, 2047, None).err(),
        Some(Error::Limit)
    );
    assert_eq!(
        source.curve_bounds(edge, [0., 1.], 2047, None).err(),
        Some(Error::Limit)
    );
    for range in [[f64::NAN, 1.], [0., f64::INFINITY], [-0.1, 1.], [0., 1.1]] {
        assert_eq!(
            source.curve_bounds(edge, range, 2048, None).err(),
            Some(Error::Invalid)
        );
    }
    for invalid in [
        Edge::StockEnd { end: 2, side: 0 },
        Edge::StockEnd { end: 0, side: 8 },
        Edge::SplitRail { corner: 3, end: 0 },
        Edge::SplitRail { corner: 4, end: 2 },
        Edge::ReliefStart(4),
        Edge::Helix(4),
        Edge::ThreadRing { end: 2, role: 0 },
        Edge::ThreadRing { end: 0, role: 4 },
    ] {
        assert_eq!(invalid.endpoints(), Err(Error::Invalid));
        assert_eq!(
            source.curve_bounds(invalid, [0., 1.], 2048, None).err(),
            Some(Error::Invalid)
        );
    }
    for invalid in [
        Vertex::Stock { end: 2, corner: 0 },
        Vertex::Stock { end: 0, corner: 8 },
        Vertex::ReliefStart(4),
        Vertex::Thread {
            end: 2,
            junction: 0,
        },
        Vertex::Thread {
            end: 0,
            junction: 4,
        },
    ] {
        assert_eq!(
            source.vertex_bounds(invalid, 2048, None).err(),
            Some(Error::Invalid)
        );
    }
    for curve in [false, true] {
        let calls = std::cell::Cell::new(0);
        let late = || {
            calls.set(calls.get() + 1);
            calls.get() == 2
        };
        let result = if curve {
            source.curve_bounds(edge, [0., 1.], 2048, Some(&late))
        } else {
            source.vertex_bounds(vertex, 2048, Some(&late))
        };
        assert_eq!(result.err(), Some(Error::Cancelled));
        assert_eq!(calls.get(), 2);
    }
    let a = source.vertex_bounds(vertex, 2048, None).unwrap();
    let b = source.curve_bounds(edge, [0., 0.], 2048, None).unwrap();
    assert_eq!(a.position_mm(), b.position_mm());
    assert_eq!(a.work_used(), b.work_used());
}
