use super::*;
use crate::{
    display::{DisplayMesh, Tessellation},
    orthogonal::{self, Box3, Operation},
};
use std::collections::{BTreeMap, BTreeSet};
fn fixture() -> (Rc<OrthogonalSolid>, Rc<[OrthogonalThroughBore]>, Vec<Seat>) {
    let stock = Rc::new(
        orthogonal::boolean(
            &[
                Box3 {
                    min: [0.; 3],
                    max: [40., 40., 5.],
                },
                Box3 {
                    min: [0.; 3],
                    max: [40., 5., 40.],
                },
            ],
            &[],
            Operation::Union,
            1e-7,
        )
        .unwrap()
        .unwrap(),
    );
    let boundary = orthogonal_boundary::build(stock.quads()).unwrap();
    let bores = [
        [[20., 20., -0.1], [20., 20., 5.1]],
        [[20., -0.1, 20.], [20., 5.1, 20.]],
    ]
    .map(|[a, b]| {
        orthogonal_pockets::identify_through_bore(&stock, &boundary, a, b, 4.5, 1e-7, 1e-5).unwrap()
    });
    (
        stock,
        Rc::from(bores),
        vec![
            Seat {
                bore: 0,
                end: false,
                radius_mm: 8.05,
                depth_mm: 3.55,
            },
            Seat {
                bore: 1,
                end: false,
                radius_mm: 8.05,
                depth_mm: 3.55,
            },
        ],
    )
}
fn options() -> Options {
    Options {
        clearance_mm: 1e-7,
        minimum_feature_mm: 1e-5,
        maximum_work: 16_000_000,
        maximum_bytes: 64 * 1024 * 1024,
    }
}
fn tessellation(deflection: f64) -> Tessellation {
    Tessellation {
        deflection,
        angle: 0.3,
        triangle_limit: 200_000,
    }
}
fn closed(mesh: &DisplayMesh) {
    let mut edges = BTreeMap::<[u32; 2], Vec<(usize, bool)>>::new();
    let mut links = BTreeMap::<u32, Vec<[u32; 2]>>::new();
    for (i, t) in mesh.indices.as_chunks::<3>().0.iter().enumerate() {
        for j in 0..3 {
            let [a, b] = [t[j], t[(j + 1) % 3]];
            edges
                .entry([a.min(b), a.max(b)])
                .or_default()
                .push((i, a < b));
            links
                .entry(t[j])
                .or_default()
                .push([t[(j + 1) % 3], t[(j + 2) % 3]]);
        }
    }
    let mut adjacency = vec![Vec::new(); mesh.face_indices.len()];
    for v in edges.values() {
        assert_eq!(v.len(), 2);
        assert_ne!(v[0].1, v[1].1);
        adjacency[v[0].0].push(v[1].0);
        adjacency[v[1].0].push(v[0].0);
    }
    let mut stack = vec![0];
    let mut reached = BTreeSet::new();
    while let Some(i) = stack.pop() {
        if reached.insert(i) {
            stack.extend(&adjacency[i]);
        }
    }
    assert_eq!(reached.len(), mesh.face_indices.len());
    let vertices = mesh.vertices.len() / 3;
    assert_eq!(
        vertices as isize - edges.len() as isize + mesh.face_indices.len() as isize,
        -2
    );
    for link in links.values() {
        let mut graph = BTreeMap::<u32, Vec<u32>>::new();
        for &[a, b] in link {
            graph.entry(a).or_default().push(b);
            graph.entry(b).or_default().push(a);
        }
        assert!(graph.values().all(|v| v.len() == 2));
        let mut stack = vec![*graph.keys().next().unwrap()];
        let mut reached = BTreeSet::new();
        while let Some(i) = stack.pop() {
            if reached.insert(i) {
                stack.extend(&graph[&i]);
            }
        }
        assert_eq!(reached.len(), graph.len());
    }
}
fn mesh_volume(mesh: &DisplayMesh) -> f64 {
    let p = |i: u32| {
        let i = i as usize * 3;
        [mesh.vertices[i], mesh.vertices[i + 1], mesh.vertices[i + 2]]
    };
    mesh.indices
        .as_chunks::<3>()
        .0
        .iter()
        .map(|t| {
            let [a, b, c] = [p(t[0]), p(t[1]), p(t[2])];
            (a[0] * (b[1] * c[2] - b[2] * c[1])
                + a[1] * (b[2] * c[0] - b[0] * c[2])
                + a[2] * (b[0] * c[1] - b[1] * c[0]))
                / 6.
        })
        .sum()
}
#[test]
fn actual_connected_angle_stock_retains_two_perpendicular_finite_countersinks() {
    let (stock, bores, seats) = fixture();
    let source = Seats::new(
        Rc::clone(&stock),
        Rc::clone(&bores),
        &seats,
        options(),
        &|| false,
    )
    .unwrap();
    assert!(std::ptr::eq(source.stock(), stock.as_ref()));
    assert_eq!(source.bores(), bores.as_ref());
    assert_eq!(source.seats(), seats);
    assert_eq!(source.stock().solid_count(), 1);
    assert_eq!(source.boundary().shell_count, 1);
    let expected =
        15000. - 2. * PI * 4.5_f64.powi(2) * 5. - 2. * PI * 3.55 * (8.05 - 4.5) * (8.05 + 9.) / 3.;
    let [lo, hi] = source.volume_bounds_mm3();
    assert!(lo < expected && expected < hi);
    assert!(hi - lo < 1e-8);
    let coarse = source.mesh(tessellation(0.1), &|| false).unwrap();
    let fine = source.mesh(tessellation(0.01), &|| false).unwrap();
    closed(&coarse);
    closed(&fine);
    assert!(fine.indices.len() > coarse.indices.len());
    let coarse_error = (mesh_volume(&coarse) - expected).abs();
    let fine_error = (mesh_volume(&fine) - expected).abs();
    eprintln!(
        "actual Angle stock volume={expected} triangles={}/{} MFG volume residual={coarse_error}/{fine_error}",
        coarse.face_indices.len(),
        fine.face_indices.len()
    );
    assert!(fine_error < coarse_error);
    assert!(fine_error < 5.);
    for owner in 0..source.boundary.faces.len() + bores.len() + seats.len() {
        assert!(fine.face_indices.contains(&(owner as u32)));
    }
    assert!(matches!(
        source.mesh(tessellation(0.1), &|| true),
        Err(Error::Cancelled)
    ));
}
#[test]
fn actual_original_stock_and_finite_seat_guards_refuse_atomically() {
    let (stock, bores, seats) = fixture();
    let original = stock.as_ref().clone();
    let qualify =
        |s: &[Seat], o: Options| Seats::new(Rc::clone(&stock), Rc::clone(&bores), s, o, &|| false);
    let good = qualify(&seats, options()).unwrap();
    assert!(matches!(
        qualify(
            &seats,
            Options {
                maximum_work: good.work_used - 1,
                ..options()
            }
        ),
        Err(Error::Limit)
    ));
    assert!(
        qualify(
            &seats,
            Options {
                maximum_work: good.work_used,
                ..options()
            }
        )
        .is_ok()
    );
    assert!(matches!(
        qualify(
            &seats,
            Options {
                maximum_bytes: 1,
                ..options()
            }
        ),
        Err(Error::Limit)
    ));
    for index in 0..seats.len() {
        let mut changed = seats.clone();
        changed[index].depth_mm = 5.;
        assert!(matches!(
            qualify(&changed, options()),
            Err(Error::Unsupported)
        ));
        changed[index].depth_mm = 3.55;
        changed[index].radius_mm = 15.;
        assert!(matches!(
            qualify(&changed, options()),
            Err(Error::Unsupported)
        ));
        changed[index].radius_mm = f64::NAN;
        assert!(matches!(qualify(&changed, options()), Err(Error::Invalid)));
    }
    let mut duplicate = seats.clone();
    duplicate.push(seats[0]);
    assert!(matches!(
        qualify(&duplicate, options()),
        Err(Error::Invalid)
    ));
    let mut wrong = seats.clone();
    wrong[0].bore = bores.len();
    assert!(matches!(qualify(&wrong, options()), Err(Error::Invalid)));
    let polls = std::cell::Cell::new(0);
    assert!(matches!(
        Seats::new(
            Rc::clone(&stock),
            Rc::clone(&bores),
            &seats,
            options(),
            &|| {
                let n = polls.get() + 1;
                polls.set(n);
                n > 4
            }
        ),
        Err(Error::Cancelled)
    ));
    assert!(polls.get() > 4);
    assert_eq!(stock.as_ref(), &original);
    assert_eq!(good.bores(), bores.as_ref());
}
#[test]
fn actual_seat_source_contracts_and_exact_retained_memory_frontier() {
    let (stock, bores, seats) = fixture();
    let qualify =
        |s: &[Seat], o: Options| Seats::new(Rc::clone(&stock), Rc::clone(&bores), s, o, &|| false);
    let good = qualify(&seats, options()).unwrap();
    for o in [
        Options {
            clearance_mm: f64::NAN,
            ..options()
        },
        Options {
            clearance_mm: -1.,
            ..options()
        },
        Options {
            minimum_feature_mm: f64::INFINITY,
            ..options()
        },
        Options {
            minimum_feature_mm: 0.,
            ..options()
        },
        Options {
            maximum_work: 0,
            ..options()
        },
        Options {
            maximum_work: 16_000_001,
            ..options()
        },
        Options {
            maximum_bytes: 0,
            ..options()
        },
        Options {
            maximum_bytes: 64 * 1024 * 1024 + 1,
            ..options()
        },
    ] {
        assert!(matches!(qualify(&seats, o), Err(Error::Invalid)));
    }
    assert!(matches!(qualify(&[], options()), Err(Error::Invalid)));
    assert!(matches!(
        Seats::new(Rc::clone(&stock), Rc::from([]), &seats, options(), &|| {
            false
        }),
        Err(Error::Invalid)
    ));
    assert!(matches!(
        qualify(&[seats[0]; 5], options()),
        Err(Error::Limit)
    ));
    let repeated: Rc<[OrthogonalThroughBore]> = Rc::from(vec![bores[0]; 17]);
    assert!(matches!(
        Seats::new(Rc::clone(&stock), repeated, &seats, options(), &|| false),
        Err(Error::Limit)
    ));
    let bytes = good.accounted_peak_bytes;
    assert!(matches!(
        qualify(
            &seats,
            Options {
                maximum_bytes: bytes - 1,
                ..options()
            }
        ),
        Err(Error::Limit)
    ));
    let exact = qualify(
        &seats,
        Options {
            maximum_bytes: bytes,
            ..options()
        },
    )
    .unwrap();
    assert_eq!(exact.volume_bounds_mm3(), good.volume_bounds_mm3());
    assert!(exact.retained_bytes() < exact.accounted_peak_bytes);
    for depth in [f64::NAN, -1., 0.] {
        let mut changed = seats.clone();
        changed[0].depth_mm = depth;
        assert!(matches!(qualify(&changed, options()), Err(Error::Invalid)));
    }
    let disconnected = Rc::new(
        orthogonal::boolean(
            &[
                Box3 {
                    min: [0.; 3],
                    max: [40.; 3],
                },
                Box3 {
                    min: [100.; 3],
                    max: [101.; 3],
                },
            ],
            &[],
            Operation::Union,
            1e-7,
        )
        .unwrap()
        .unwrap(),
    );
    assert_eq!(disconnected.solid_count(), 2);
    assert!(matches!(
        Seats::new(disconnected, bores, &seats, options(), &|| false),
        Err(Error::Unsupported)
    ));
    assert_eq!(stock.as_ref(), good.stock());
}

#[test]
fn actual_opposite_end_cones_and_reversed_bore_requests_keep_closed_material() {
    let (stock, original_bores, _) = fixture();
    let boundary = orthogonal_boundary::build(stock.quads()).unwrap();
    let reversed = [
        [[20., 20., 5.1], [20., 20., -0.1]],
        [[20., 5.1, 20.], [20., -0.1, 20.]],
    ]
    .map(|[a, b]| {
        orthogonal_pockets::identify_through_bore(&stock, &boundary, a, b, 4.5, 1e-7, 1e-5).unwrap()
    });
    let mut mass = None;
    for bores in [original_bores, Rc::from(reversed)] {
        let seats = (0..2)
            .flat_map(|bore| {
                [false, true].map(move |end| Seat {
                    bore,
                    end,
                    radius_mm: 6.,
                    depth_mm: 0.75,
                })
            })
            .collect::<Vec<_>>();
        let source = Seats::new(Rc::clone(&stock), bores, &seats, options(), &|| false).unwrap();
        let bounds = source.volume_bounds_mm3();
        if let Some(previous) = mass {
            assert_eq!(bounds, previous);
        } else {
            mass = Some(bounds);
        }
        let mesh = source.mesh(tessellation(0.02), &|| false).unwrap();
        closed(&mesh);
        for owner in 0..source.boundary().faces.len() + 2 + 4 {
            assert!(mesh.face_indices.contains(&(owner as u32)));
        }
        assert!((mesh_volume(&mesh) - bounds[0].midpoint(bounds[1])).abs() < 2.);
        for bad in [
            Tessellation {
                deflection: 0.,
                ..tessellation(0.1)
            },
            Tessellation {
                angle: f64::NAN,
                ..tessellation(0.1)
            },
            Tessellation {
                triangle_limit: 0,
                ..tessellation(0.1)
            },
        ] {
            assert!(matches!(source.mesh(bad, &|| false), Err(Error::Invalid)));
        }
        assert!(matches!(
            source.mesh(
                Tessellation {
                    triangle_limit: 1,
                    ..tessellation(0.1)
                },
                &|| false
            ),
            Err(Error::Limit)
        ));
        assert_eq!(source.volume_bounds_mm3(), bounds);
    }
}

#[test]
fn actual_source_and_mesh_cancellation_leave_retained_material_reusable() {
    use std::cell::Cell;
    let (stock, bores, seats) = fixture();
    let polls = Cell::new(0usize);
    let count = || {
        polls.set(polls.get() + 1);
        false
    };
    let source = Seats::new(
        Rc::clone(&stock),
        Rc::clone(&bores),
        &seats,
        options(),
        &count,
    )
    .unwrap();
    let total = polls.get();
    assert!(total > 4);
    let volume = source.volume_bounds_mm3();
    for stop in 1..=total {
        polls.set(0);
        let interrupt = || {
            polls.set(polls.get() + 1);
            polls.get() >= stop
        };
        assert!(
            matches!(
                Seats::new(
                    Rc::clone(&stock),
                    Rc::clone(&bores),
                    &seats,
                    options(),
                    &interrupt
                ),
                Err(Error::Cancelled)
            ),
            "source checkpoint {stop}/{total}"
        );
        assert_eq!(source.bores(), bores.as_ref());
        assert_eq!(source.volume_bounds_mm3(), volume);
    }
    polls.set(0);
    let before = source.mesh(tessellation(0.1), &count).unwrap();
    let total = polls.get();
    assert!(total > before.face_indices.len());
    for stop in [1, 2, 3, total / 2, total - 1, total] {
        polls.set(0);
        let interrupt = || {
            polls.set(polls.get() + 1);
            polls.get() >= stop
        };
        assert!(
            matches!(
                source.mesh(tessellation(0.1), &interrupt),
                Err(Error::Cancelled)
            ),
            "mesh checkpoint {stop}/{total}"
        );
    }
    let after = source.mesh(tessellation(0.1), &|| false).unwrap();
    assert_eq!(before.vertices, after.vertices);
    assert_eq!(before.indices, after.indices);
    assert_eq!(before.face_indices, after.face_indices);
    assert_eq!(source.volume_bounds_mm3(), volume);
}

#[test]
fn actual_unvalidated_bore_bindings_and_far_coordinate_mesh_precision_refuse() {
    let (stock, bores, seats) = fixture();
    for mutate in 0..2 {
        let mut bad = bores.to_vec();
        if mutate == 0 {
            bad[0].faces[0] = usize::MAX;
        } else {
            bad[0].radius_mm = -1.;
        }
        assert!(matches!(
            Seats::new(Rc::clone(&stock), Rc::from(bad), &seats, options(), &|| {
                false
            }),
            Err(Error::Invalid)
        ));
    }
    let shift = 1e9;
    let stock = Rc::new(
        orthogonal::boolean(
            &[
                Box3 {
                    min: [shift, 0., 0.],
                    max: [shift + 40., 40., 5.],
                },
                Box3 {
                    min: [shift, 0., 0.],
                    max: [shift + 40., 5., 40.],
                },
            ],
            &[],
            Operation::Union,
            1e-7,
        )
        .unwrap()
        .unwrap(),
    );
    let boundary = orthogonal_boundary::build(stock.quads()).unwrap();
    let bores = [
        [[shift + 20., 20., -0.1], [shift + 20., 20., 5.1]],
        [[shift + 20., -0.1, 20.], [shift + 20., 5.1, 20.]],
    ]
    .map(|[a, b]| {
        orthogonal_pockets::identify_through_bore(&stock, &boundary, a, b, 4.5, 1e-7, 1e-5).unwrap()
    });
    // This is an actual Pure source; no Native coordinate-policy admission is
    // inferred. Its larger rounding reserve must constrain mesh publication.
    let source = Seats::new(stock, Rc::from(bores), &seats, options(), &|| false).unwrap();
    let volume = source.volume_bounds_mm3();
    assert!(volume[0] > 0.);
    assert!(matches!(
        source.mesh(tessellation(0.001), &|| false),
        Err(Error::Unsupported)
    ));
    let coarse = source.mesh(tessellation(0.1), &|| false).unwrap();
    closed(&coarse);
    assert_eq!(source.volume_bounds_mm3(), volume);
    assert!(matches!(
        source.mesh(
            Tessellation {
                triangle_limit: coarse.face_indices.len() - 1,
                ..tessellation(0.1)
            },
            &|| false
        ),
        Err(Error::Limit)
    ));
}

#[test]
fn seated_stock_receiving_proof_retains_every_original_cell_and_both_members() {
    use crate::{
        normal_sweep::{Edge, Sweep},
        wood_thread_contact as contact,
    };
    let (stock, bores, seats) = fixture();
    let source = Seats::new(stock, bores, &seats, options(), &|| false).unwrap();
    let ring = [[0., 0.], [40., 0.], [40., 40.], [0., 40.]];
    let receiver = Sweep::new(
        vec![
            (0..4)
                .map(|i| Edge::Line {
                    a: ring[i],
                    b: ring[(i + 1) % 4],
                })
                .collect(),
        ],
        100.,
        1e-7,
        1e-5,
    )
    .unwrap();
    let identity = contact::Frame {
        origin_mm: [0.; 3],
        axes: [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
    };
    let stock = contact::SeatedStockPlacement {
        source: &source,
        frame: contact::Frame {
            origin_mm: [40., 40., 20.],
            ..identity
        },
    };
    let receiver = contact::ReceivingSweep {
        source: &receiver,
        frame: identity,
        occurrence: identity,
        reverse: false,
    };
    let options = contact::ReceiverClearanceOptions {
        maximum_overlap_volume_mm3: 0.001,
        maximum_work: 16_000_000,
    };
    for occurrence in [
        identity,
        contact::Frame {
            origin_mm: [40., 0., 0.],
            ..identity
        },
    ] {
        let actual = contact::ReceivingSweep {
            occurrence,
            ..receiver
        };
        let proof =
            contact::certify_seated_stock_receiver_clearance(stock, actual, options, &|| false)
                .unwrap();
        assert!(std::ptr::eq(proof.stock().source, &source));
        assert!(std::ptr::eq(proof.receiver().source, receiver.source));
        assert_eq!(proof.overlap_volume_upper_mm3, 0.);
        assert_eq!(
            proof.source_cell_upper_mm3.len(),
            source.stock().cells().len()
        );
        assert!(proof.source_cell_upper_mm3.iter().all(|v| *v == 0.));
        assert!(matches!(
            contact::certify_seated_stock_receiver_clearance(
                stock,
                actual,
                contact::ReceiverClearanceOptions {
                    maximum_work: proof.work_used - 1,
                    ..options
                },
                &|| false
            ),
            Err(contact::Error::Limit)
        ));
    }
    let penetrating = contact::SeatedStockPlacement {
        frame: contact::Frame {
            origin_mm: [39., 39., 20.],
            ..identity
        },
        ..stock
    };
    assert!(matches!(
        contact::certify_seated_stock_receiver_clearance(penetrating, receiver, options, &|| false),
        Err(contact::Error::Unsupported)
    ));
    let polls = std::cell::Cell::new(0);
    assert!(matches!(
        contact::certify_seated_stock_receiver_clearance(stock, receiver, options, &|| {
            let n = polls.get() + 1;
            polls.set(n);
            n > 12
        }),
        Err(contact::Error::Cancelled)
    ));
    assert!(polls.get() > 12);
}
