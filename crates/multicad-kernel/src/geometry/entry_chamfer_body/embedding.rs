//! Literal rounded triangles admit only their indexed common faces.
//!
//! Candidate separating axes are heuristic. Only outward arithmetic signs
//! decide separation; unresolved pairs refuse publication.
use super::placement::{cross_error, product, up};
use super::*;
use crate::geometry::arithmetic::{difference_with_error, dot_error};
#[derive(Clone, Copy, Debug)]
struct Vector {
    value: Point,
    error: Point,
}
fn difference(a: Point, b: Point) -> Vector {
    let (value, error) = difference_with_error(a, b);
    Vector { value, error }
}
fn vector(v: Point) -> Vector {
    Vector {
        value: v,
        error: [0.; 3],
    }
}
fn cross(a: Vector, b: Vector) -> Vector {
    let (value, error) = cross_error(a.value, a.error, b.value, b.error);
    Vector { value, error }
}
fn dot(a: Vector, b: Vector) -> [f64; 2] {
    let (v, mut e) = dot_error(a.value, b.value);
    for i in 0..3 {
        e = up(e, product(a.value[i].abs(), b.error[i]));
        e = up(e, product(b.value[i].abs(), a.error[i]));
        e = up(e, product(a.error[i], b.error[i]));
    }
    [(v - e).next_down(), (v + e).next_up()]
}
fn normalize(v: Point) -> Option<Point> {
    let scale = v.into_iter().map(f64::abs).fold(0., f64::max);
    (scale.is_finite() && scale > 0.).then(|| v.map(|x| x / scale))
}
fn raw_cross(a: Point, b: Point) -> Point {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn sub(a: Point, b: Point) -> Point {
    std::array::from_fn(|i| a[i] - b[i])
}
fn rays_separate(v: [Vector; 4]) -> bool {
    let Some(p) = v
        .iter()
        .map(|v| normalize(v.value))
        .collect::<Option<Vec<_>>>()
    else {
        return false;
    };
    let accepts =
        |k: Point| normalize(k).is_some_and(|k| v.iter().all(|v| dot(vector(k), *v)[0] > 0.));
    if accepts(std::array::from_fn(|i| p.iter().map(|p| p[i]).sum())) {
        return true;
    }
    for i in 0..4 {
        if accepts(p[i]) {
            return true;
        }
        for j in i + 1..4 {
            let d = sub(p[j], p[i]);
            let length = d.iter().map(|x| x * x).sum::<f64>();
            if length > 0. {
                let t =
                    (-p[i].iter().zip(d).map(|(a, b)| a * b).sum::<f64>() / length).clamp(0., 1.);
                if accepts(std::array::from_fn(|k| p[i][k] + t * d[k])) {
                    return true;
                }
            }
            for k in j + 1..4 {
                let n = raw_cross(sub(p[j], p[i]), sub(p[k], p[i]));
                if accepts(n) || accepts(n.map(|x| -x)) {
                    return true;
                }
            }
        }
    }
    false
}
fn separated(a: [Point; 3], b: [Point; 3]) -> bool {
    let ea = [sub(a[1], a[0]), sub(a[2], a[1]), sub(a[0], a[2])];
    let eb = [sub(b[1], b[0]), sub(b[2], b[1]), sub(b[0], b[2])];
    let na = raw_cross(ea[0], ea[1]);
    let nb = raw_cross(eb[0], eb[1]);
    let mut axes = vec![na, nb];
    axes.extend(ea.map(|e| raw_cross(na, e)));
    axes.extend(eb.map(|e| raw_cross(nb, e)));
    for a in ea {
        for b in eb {
            axes.push(raw_cross(a, b));
        }
    }
    axes.into_iter().filter_map(normalize).any(|axis| {
        let bounds = |p: [Point; 3]| {
            let values = p.map(|p| dot(vector(axis), difference(p, a[0])));
            [
                values.iter().map(|p| p[0]).fold(f64::INFINITY, f64::min),
                values
                    .iter()
                    .map(|p| p[1])
                    .fold(f64::NEG_INFINITY, f64::max),
            ]
        };
        let [a, b] = [bounds(a), bounds(b)];
        a[1] < b[0] || b[1] < a[0]
    })
}
fn pair(points: &[[f64; 3]], a: [u32; 3], b: [u32; 3]) -> bool {
    let common = a.into_iter().filter(|v| b.contains(v)).collect::<Vec<_>>();
    let p = |v: u32| points[v as usize];
    match common.as_slice() {
        [] => separated(a.map(p), b.map(p)),
        [at] => {
            let signed = a
                .into_iter()
                .filter(|v| v != at)
                .map(|v| difference(p(v), p(*at)))
                .chain(b.into_iter().filter(|v| v != at).map(|v| {
                    let d = difference(p(v), p(*at));
                    Vector {
                        value: d.value.map(|x| -x),
                        error: d.error,
                    }
                }))
                .collect::<Vec<_>>();
            rays_separate(
                signed
                    .try_into()
                    .expect("two nonshared vertices per triangle"),
            )
        }
        [at, end] => {
            let along = difference(p(*end), p(*at));
            let u = difference(p(*a.iter().find(|v| !common.contains(v)).unwrap()), p(*at));
            let v = difference(p(*b.iter().find(|v| !common.contains(v)).unwrap()), p(*at));
            let one = cross(along, u);
            let two = cross(along, v);
            let det = dot(one, v);
            det[0] > 0. || det[1] < 0. || dot(one, two)[1] < 0.
        }
        _ => false,
    }
}
// A literal projected fan has no nonadjacent wedge intersections when its
// directed rim makes exactly one turn and each triangle advances by < pi.
// This proves actual coordinates, without trusting a semantic face label.
fn fans(points: &[Point], triangles: &[[u32; 3]], faces: &[u32]) -> Vec<Option<usize>> {
    use std::collections::{BTreeMap, BTreeSet};
    let mut groups = BTreeMap::<u32, Vec<usize>>::new();
    for (i, &face) in faces.iter().enumerate() {
        groups.entry(face).or_default().push(i);
    }
    let mut result = vec![None; triangles.len()];
    for (group, ids) in groups.values().enumerate() {
        if ids.len() < 3 {
            continue;
        }
        let apex = triangles[ids[0]]
            .into_iter()
            .find(|v| ids.iter().all(|&i| triangles[i].contains(v)));
        let Some(apex) = apex else {
            continue;
        };
        let mut rim = BTreeMap::new();
        let mut unique = true;
        for &id in ids {
            let t = triangles[id];
            let at = t.iter().position(|&v| v == apex).unwrap();
            if rim.insert(t[(at + 1) % 3], t[(at + 2) % 3]).is_some() {
                unique = false;
                break;
            }
        }
        if !unique {
            continue;
        }
        let start = *rim.keys().next().unwrap();
        let mut current = start;
        let mut seen = BTreeSet::new();
        while seen.insert(current) {
            let Some(&next) = rim.get(&current) else {
                break;
            };
            current = next;
        }
        if current != start || seen.len() != ids.len() {
            continue;
        }
        let center = points[apex as usize];
        let proved = [(0, 1), (1, 2), (2, 0)].into_iter().any(|(x, y)| {
            let half = |id: u32| {
                let p = points[id as usize];
                p[y] > center[y] || (p[y] == center[y] && p[x] > center[x])
            };
            let mut positive = 0;
            let mut negative = 0;
            let mut crossings = 0;
            for (&a, &b) in &rim {
                let a_vector = difference(points[a as usize], center);
                let b_vector = difference(points[b as usize], center);
                let (n, e) = cross_error(
                    a_vector.value,
                    a_vector.error,
                    b_vector.value,
                    b_vector.error,
                );
                let z = 3 - x - y;
                positive += usize::from(n[z] > e[z]);
                negative += usize::from(n[z] < -e[z]);
                crossings += usize::from(!half(a) && half(b));
            }
            (positive == ids.len() || negative == ids.len()) && crossings == 1
        });
        if proved {
            for &id in ids {
                result[id] = Some(group);
            }
        }
    }
    result
}
#[derive(Clone, Copy)]
struct Node {
    bounds: [Point; 2],
    children: Option<[usize; 2]>,
    triangle: usize,
    fan: Option<usize>,
}
fn disjoint(a: [Point; 2], b: [Point; 2]) -> bool {
    (0..3).any(|i| a[1][i] < b[0][i] || b[1][i] < a[0][i])
}
struct Work<'a> {
    remaining: usize,
    interrupt: Option<&'a dyn Fn() -> bool>,
}
impl Work<'_> {
    fn step(&mut self, amount: usize) -> Result<(), GeometryError> {
        if self.interrupt.is_some_and(|f| f()) {
            return Err(failure(source::Error::Cancelled));
        }
        self.remaining = self.remaining.checked_sub(amount).ok_or_else(|| {
            error(
                "Entry mesh embedding",
                GeometryErrorKind::Limit,
                "Literal triangle separation exhausted its bounded work",
            )
        })?;
        Ok(())
    }
}
fn build(
    ids: &mut [usize],
    boxes: &[[Point; 2]],
    fans: &[Option<usize>],
    nodes: &mut Vec<Node>,
    work: &mut Work<'_>,
) -> Result<usize, GeometryError> {
    work.step(1 + ids.len())?;
    let bounds = std::array::from_fn(|side| {
        std::array::from_fn(|axis| {
            ids.iter().map(|&i| boxes[i][side][axis]).fold(
                if side == 0 {
                    f64::INFINITY
                } else {
                    f64::NEG_INFINITY
                },
                if side == 0 { f64::min } else { f64::max },
            )
        })
    });
    let children = if ids.len() == 1 {
        None
    } else {
        let axis = (0..3)
            .max_by(|&a, &b| {
                (bounds[1][a] - bounds[0][a]).total_cmp(&(bounds[1][b] - bounds[0][b]))
            })
            .unwrap();
        ids.sort_unstable_by(|&a, &b| {
            (boxes[a][0][axis] + boxes[a][1][axis])
                .total_cmp(&(boxes[b][0][axis] + boxes[b][1][axis]))
        });
        let middle = ids.len() / 2;
        let (a, b) = ids.split_at_mut(middle);
        Some([
            build(a, boxes, fans, nodes, work)?,
            build(b, boxes, fans, nodes, work)?,
        ])
    };
    let fan = match children {
        None => fans[ids[0]],
        Some([a, b]) => nodes[a].fan.filter(|g| Some(*g) == nodes[b].fan),
    };
    let id = nodes.len();
    nodes.push(Node {
        bounds,
        children,
        triangle: ids[0],
        fan,
    });
    Ok(id)
}
pub(super) fn certify(mesh: &Mesh) -> Result<usize, GeometryError> {
    controlled(mesh, None, 64_000_000)
}
fn controlled(
    mesh: &Mesh,
    interrupt: Option<&dyn Fn() -> bool>,
    budget: usize,
) -> Result<usize, GeometryError> {
    #[cfg(test)]
    eprintln!("literal embedding facets {}", mesh.indices.len() / 3);
    let mut work = Work {
        remaining: budget,
        interrupt,
    };
    work.step(1)?;
    let (points, tail) = mesh.vertices.as_chunks::<3>();
    let (triangles, index_tail) = mesh.indices.as_chunks::<3>();
    let refused = || {
        error(
            "Entry mesh embedding",
            GeometryErrorKind::Unsupported,
            "Literal mesh triangles are intersecting or their separation is unresolved",
        )
    };
    if !tail.is_empty()
        || !index_tail.is_empty()
        || triangles.is_empty()
        || mesh.face_indices.len() != triangles.len()
        || points.iter().flatten().any(|x| !x.is_finite())
        || mesh.indices.iter().any(|&i| i as usize >= points.len())
    {
        return Err(refused());
    }
    if triangles.len() > 1_000_000 {
        return Err(error(
            "Entry mesh embedding",
            GeometryErrorKind::Limit,
            "Literal triangle separation exceeds its memory bound",
        ));
    }
    for t in triangles {
        work.step(1)?;
        let [a, b, c] = t.map(|i| points[i as usize]);
        let n = cross(difference(b, a), difference(c, a));
        if !(0..3).any(|i| n.value[i].abs() > n.error[i]) {
            return Err(refused());
        }
    }
    work.step(triangles.len().saturating_mul(32))?;
    let fans = fans(points, triangles, &mesh.face_indices);
    work.step(1)?;
    let boxes = triangles
        .iter()
        .map(|t| {
            std::array::from_fn(|side| {
                std::array::from_fn(|axis| {
                    t.iter().map(|&i| points[i as usize][axis]).fold(
                        if side == 0 {
                            f64::INFINITY
                        } else {
                            f64::NEG_INFINITY
                        },
                        if side == 0 { f64::min } else { f64::max },
                    )
                })
            })
        })
        .collect::<Vec<[Point; 2]>>();
    let mut nodes = Vec::with_capacity(2 * triangles.len());
    let mut ids = (0..triangles.len()).collect::<Vec<_>>();
    let root = build(&mut ids, &boxes, &fans, &mut nodes, &mut work)?;
    let mut pending = vec![(root, root)];
    let mut pairs = 0;
    while let Some((a, b)) = pending.pop() {
        work.step(1)?;
        let one = nodes[a];
        let two = nodes[b];
        if one.fan.is_some() && one.fan == two.fan {
            continue;
        }
        if a == b {
            if let Some([l, r]) = one.children {
                pending.extend([(l, l), (r, r), (l, r)]);
            }
            continue;
        }
        if disjoint(one.bounds, two.bounds) {
            continue;
        }
        match (one.children, two.children) {
            (None, None) => {
                if fans[one.triangle].is_some() && fans[one.triangle] == fans[two.triangle] {
                    continue;
                }
                pairs += 1;
                if !pair(points, triangles[one.triangle], triangles[two.triangle]) {
                    return Err(refused());
                }
            }
            (Some([l, r]), _) => pending.extend([(l, b), (r, b)]),
            (None, Some([l, r])) => pending.extend([(a, l), (a, r)]),
        }
    }
    Ok(pairs)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn mesh(p: &[Point], triangles: &[[u32; 3]]) -> Mesh {
        Mesh {
            vertices: p.iter().flatten().copied().collect(),
            indices: triangles.iter().flatten().copied().collect(),
            face_indices: vec![0; triangles.len()],
        }
    }
    #[test]
    fn projected_actual_fans_require_one_cycle_one_winding_and_distinct_rays() {
        let mut p = vec![[1e6, -2e6, 3e6]];
        for i in 0..5 {
            let t = std::f64::consts::TAU * i as f64 / 5.;
            p.push([1e6 + t.cos(), -2e6 + t.sin(), 3e6]);
        }
        let cycle = |order: &[u32]| {
            order
                .iter()
                .copied()
                .zip(order.iter().copied().cycle().skip(1))
                .take(order.len())
                .map(|(a, b)| [0, a, b])
                .collect::<Vec<_>>()
        };
        let good = cycle(&[1, 2, 3, 4, 5]);
        assert!(fans(&p, &good, &[0; 5]).iter().all(Option::is_some));
        assert_eq!(certify(&mesh(&p, &good)).unwrap(), 0);
        let reversed = cycle(&[5, 4, 3, 2, 1]);
        assert!(fans(&p, &reversed, &[0; 5]).iter().all(Option::is_some));
        for order in [[1, 3, 5, 2, 4], [1, 2, 4, 3, 5]] {
            let bad = cycle(&order);
            assert!(fans(&p, &bad, &[0; 5]).iter().all(Option::is_none));
            assert!(certify(&mesh(&p, &bad)).is_err());
        }
        p[2] = p[1];
        assert!(fans(&p, &good, &[0; 5]).iter().all(Option::is_none));
        assert!(certify(&mesh(&p, &good)).is_err());
    }
    #[test]
    fn literal_candidates_preserve_only_common_indexed_faces_and_refuse_crossings() {
        let points = [
            [0., 0., 0.],
            [2., 0., 0.],
            [0., 2., 0.],
            [0., 0., 2.],
            [0., -2., 0.],
            [-2., 0., 0.],
            [0.2, 0.2, 0.],
            [1.5, 0.2, 0.],
            [0.2, 1.5, 0.],
        ];
        let a = [0, 1, 2];
        // Shared edge, planar opposite sides, and noncoplanar faces.
        assert!(pair(&points, a, [0, 4, 1]));
        assert!(pair(&points, a, [0, 1, 3]));
        // A shared point alone cannot excuse an overlap of the two cones.
        assert!(pair(&points, a, [0, 3, 5]));
        assert!(!pair(&points, a, [0, 6, 7]));
        assert!(!pair(&points, a, [0, 1, 6]));
        assert!(!pair(&points, a, a));
        assert!(!pair(&points, a, [6, 7, 8]));
        let folded = mesh(&points, &[a, [6, 7, 8]]);
        assert!(certify(&folded).is_err());
        let crossing = mesh(
            &[
                [0., 0., 0.],
                [2., 0., 0.],
                [0., 2., 0.],
                [0.5, 0.5, -1.],
                [0.5, 0.5, 1.],
                [2., 2., 0.],
            ],
            &[[0, 1, 2], [3, 4, 5]],
        );
        assert!(certify(&crossing).is_err());
        let separate = mesh(
            &[
                [0., 0., 0.],
                [1., 0., 0.],
                [0., 1., 0.],
                [2., 2., 0.],
                [3., 2., 0.],
                [2., 3., 0.],
            ],
            &[[0, 1, 2], [3, 4, 5]],
        );
        assert!(certify(&separate).is_ok());
    }
    #[test]
    fn actual_tetrahedra_survive_far_proper_lift_and_cancellation_and_limits_refuse() {
        let p = [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
        let triangles = [[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]];
        for frame in [
            PlaneFrame::default(),
            PlaneFrame {
                origin: [1e6, -2e6, 3e6],
                u: [0.8, 0.6, 0.],
                v: [-0.36, 0.48, 0.8],
            },
        ] {
            let m = mesh(&p.map(|p| position(frame, p)), &triangles);
            assert_eq!(certify(&m).unwrap(), 6);
            assert!(controlled(&m, None, 1).is_err());
            let count = std::cell::Cell::new(0);
            let stop = || {
                count.set(count.get() + 1);
                count.get() >= 15
            };
            assert_eq!(
                controlled(&m, Some(&stop), 16_000_000).err().unwrap().kind,
                GeometryErrorKind::Cancelled
            );
        }
        let mut invalid = mesh(&p, &triangles);
        invalid.indices[0] = 100;
        assert!(certify(&invalid).is_err());
        let degenerate = mesh(&[[0.; 3]; 3], &[[0, 1, 2]]);
        assert!(certify(&degenerate).is_err());
    }
    #[test]
    fn authored_internal_treatment_candidates_require_actual_literal_embedding() {
        use multicad_rust_kernel::{
            authored_flat_thread::{AuthoredThread, SourceOptions},
            authored_internal_treatment::{Cavity, Intent, Terminal, tessellate_candidate},
            display::Tessellation,
            helical_thread::Spec,
        };
        for (left_hand, starts) in [(false, 1), (true, 2)] {
            let pitch_mid = (2.675 + 2.775) / 2.;
            let minor_mid = (2.459 + 2.599) / 2.;
            let thread = AuthoredThread::new(
                Spec {
                    major_radius_mm: 1.5,
                    pitch_mm: 0.5,
                    length_mm: 3.3,
                    depth_mm: (3. - minor_mid) / 2.,
                    crest_width_mm: 0.25 - (3. - pitch_mid) / 3_f64.sqrt(),
                    lower_flank_deg: 30.,
                    upper_flank_deg: 30.,
                    starts,
                    left_hand,
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
                    entry_depth_mm: Some(0.4),
                    terminal: Some(Terminal {
                        width_mm: 0.3,
                        radial_depth_mm: 0.05,
                    }),
                },
            )
            .unwrap();
            let candidate = tessellate_candidate(
                &cavity,
                Tessellation {
                    deflection: 0.1,
                    angle: 0.3,
                    triangle_limit: 1_000_000,
                },
                None,
            )
            .unwrap();
            let mesh = Mesh {
                vertices: candidate.mesh.vertices,
                indices: candidate.mesh.indices,
                face_indices: candidate.mesh.face_indices,
            };
            let timer = std::time::Instant::now();
            let pairs = certify(&mesh).unwrap();
            eprintln!(
                "actual internal candidate literal embedding {left_hand}/{starts}: {pairs} pairs {:?}",
                timer.elapsed()
            );
            assert!(pairs > 0);
            let mut overlap = mesh.clone();
            let first: [u32; 3] = overlap.indices[..3].try_into().unwrap();
            overlap.indices[3..6].copy_from_slice(&first);
            assert!(certify(&overlap).is_err());
            assert_eq!(
                controlled(&mesh, Some(&|| true), 64_000_000)
                    .unwrap_err()
                    .kind,
                GeometryErrorKind::Cancelled
            );
            assert!(controlled(&mesh, None, 1).is_err());
        }
    }
}
