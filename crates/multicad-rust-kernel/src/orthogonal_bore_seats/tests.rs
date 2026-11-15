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
