use super::*;
use multicad_precision::PortableMath as _;
type RawGraph = (Vec<[I; 3]>, Vec<[usize; 2]>, Vec<Vec<EdgeUse>>);
fn graph() -> RawGraph {
    let (s, d) = super::super::super::tests::fixture(1);
    let c = Cutter::new(s, d).unwrap();
    let body = c.boundary_topology().unwrap();
    (
        body.vertices.iter().map(|v| v.point).collect(),
        body.edges.iter().map(|e| e.ends).collect(),
        body.faces.iter().map(|f| f.uses.clone()).collect(),
    )
}
#[test]
fn malformed_raw_incidence_and_zero_length_edges_are_rejected() {
    let (p, e, f) = graph();
    assert!(certify(&[], &e, &f).is_err());
    let mut bad = e.clone();
    bad[0] = [0, p.len()];
    assert!(certify(&p, &bad, &f).is_err());
    bad[0] = [0, 0];
    assert!(certify(&p, &bad, &f).is_err());
    let mut points = p.clone();
    points[e[0][0]] = points[e[0][1]];
    assert!(certify(&points, &e, &f).is_err());
    let mut faces = f.clone();
    faces[0].truncate(2);
    assert!(certify(&p, &e, &faces).is_err());
    faces = f.clone();
    faces[0][0].edge = e.len();
    assert!(certify(&p, &e, &faces).is_err());
    faces = f.clone();
    faces[0][1].edge = e.len();
    assert!(certify(&p, &e, &faces).is_err());
    faces = f.clone();
    faces[0][0].forward = !faces[0][0].forward;
    assert!(certify(&p, &e, &faces).is_err());
    faces = f.clone();
    faces.push(faces[0].clone());
    assert!(certify(&p, &e, &faces).is_err());
    faces = f.clone();
    faces.remove(1);
    assert!(certify(&p, &e, &faces).is_err());
    // A closed walk visiting a vertex twice is not a simple face boundary.
    faces = f.clone();
    faces[0].extend(f[0].clone());
    assert!(certify(&p, &e, &faces).is_err());
}
#[test]
fn vertex_links_must_be_connected_cycles_and_component_must_be_single() {
    assert!(!connected(&[]));
    assert!(!connected(&[vec![1], vec![0], vec![3], vec![2]]));
    assert!(connected(&[vec![1], vec![0, 2], vec![1]]));
    assert!(!vertex_link(&[0], &[]));
    assert!(!vertex_link(&[0, 1], &[(2, 1)]));
    assert!(!vertex_link(&[0, 1], &[(0, 2)]));
    assert!(!vertex_link(&[0, 1], &[(0, 0)]));
    assert!(!vertex_link(&[0, 1], &[(0, 1)]));
    // Link of a pinched union of two spheres: each incident edge has degree2,
    // but the two separate cycles do not form one manifold vertex neighborhood.
    assert!(!vertex_link(
        &[0, 1, 2, 3, 4, 5],
        &[(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3)]
    ));
    assert!(vertex_link(&[0, 1, 2], &[(0, 1), (1, 2), (2, 0)]));
}
fn raw_graph(points: Vec<[f64; 3]>, cycles: Vec<Vec<usize>>) -> RawGraph {
    let mut edges = Vec::new();
    let mut faces = Vec::new();
    for cycle in cycles {
        let mut uses = Vec::new();
        for n in 0..cycle.len() {
            let pair = [cycle[n], cycle[(n + 1) % cycle.len()]];
            let found = edges
                .iter()
                .position(|e: &[usize; 2]| *e == pair || *e == [pair[1], pair[0]]);
            let edge = found.unwrap_or_else(|| {
                edges.push(pair);
                edges.len() - 1
            });
            uses.push(EdgeUse {
                edge,
                forward: edges[edge] == pair,
            })
        }
        faces.push(uses)
    }
    (
        points.into_iter().map(|p| p.map(I::point)).collect(),
        edges,
        faces,
    )
}
#[test]
fn actual_disconnected_pinched_and_genus_one_closed_graphs_are_not_the_canonical_body() {
    let tetra = [vec![0, 2, 1], vec![0, 1, 3], vec![0, 3, 2], vec![1, 2, 3]];
    let first = vec![[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
    let (p, e, f) = raw_graph(first.clone(), tetra.to_vec());
    assert_eq!(certify(&p, &e, &f).unwrap().euler_characteristic, 2);
    let mut points = first.clone();
    points.extend(first.iter().map(|[x, y, z]| [x + 3., *y, *z]));
    let mut cycles = tetra.to_vec();
    cycles.extend(tetra.iter().map(|c| c.iter().map(|v| v + 4).collect()));
    let (p, e, f) = raw_graph(points, cycles);
    assert!(certify(&p, &e, &f).is_err());
    let mut points = first;
    points.extend([[-1., 0., 0.], [0., -1., 0.], [0., 0., -1.]]);
    let mut cycles = tetra.to_vec();
    cycles.extend(tetra.iter().map(|c| {
        c.iter()
            .rev()
            .map(|v| if *v == 0 { 0 } else { v + 3 })
            .collect()
    }));
    let (p, e, f) = raw_graph(points, cycles);
    assert!(certify(&p, &e, &f).is_err());
    let mut points = Vec::new();
    for u in 0..3 {
        for v in 0..3 {
            let a = u as f64 * 2. * PI / 3.;
            let b = v as f64 * 2. * PI / 3.;
            points.push([
                (3. + b.portable_cos()) * a.portable_cos(),
                (3. + b.portable_cos()) * a.portable_sin(),
                b.portable_sin(),
            ])
        }
    }
    let mut cycles = Vec::new();
    for u in 0..3 {
        for v in 0..3 {
            cycles.push(vec![
                u * 3 + v,
                ((u + 1) % 3) * 3 + v,
                ((u + 1) % 3) * 3 + (v + 1) % 3,
                u * 3 + (v + 1) % 3,
            ])
        }
    }
    let (p, e, f) = raw_graph(points, cycles);
    assert_eq!(p.len() as i64 - e.len() as i64 + f.len() as i64, 0);
    assert!(certify(&p, &e, &f).is_err());
}
