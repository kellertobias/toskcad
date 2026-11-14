//! Source-owned parameter cycles and independently checked literal radial graphs.
use super::*;
pub(super) enum Wall<'a> {
    Wood(&'a ThreadShell),
    Metric(&'a MetricShell),
}
impl Wall<'_> {
    fn mesh(&self) -> &DisplayMesh {
        match self {
            Self::Wood(s) => s.mesh(),
            Self::Metric(s) => s.mesh(),
        }
    }
    fn count(&self, end: ShellEnd) -> usize {
        match self {
            Self::Wood(s) => s.rim_count(end),
            Self::Metric(s) => s.rim_count(end),
        }
    }
    fn angle(&self, id: usize, w: &mut Work<'_>) -> R<num_rational::BigRational> {
        let (q, cost) = match self {
            Self::Wood(s) => s.exact_rim_turns(id, w.remaining, w.interrupt)?,
            Self::Metric(s) => s.exact_rim_turns(id, w.remaining, w.interrupt)?,
        };
        w.charge(cost)?;
        Ok(q)
    }
    fn edge(&self, end: ShellEnd, i: usize) -> R<([usize; 2], f64, f64)> {
        let (ids, r, m, dd, error, lead) = match self {
            Self::Wood(s) => {
                let e = s.rim_edge(end, i).ok_or(E::Unresolved)?;
                let r = I::new(
                    s.source().minimum_radius_mm(),
                    s.source().maximum_radius_mm(),
                );
                (
                    e.vertex_indices,
                    r,
                    I::new(e.phase_slope[0], e.phase_slope[1]),
                    I::point(0.),
                    e.point_error_mm,
                    I::point(s.source().spec().pitch_mm)
                        * I::point(f64::from(s.source().spec().starts)),
                )
            }
            Self::Metric(s) => {
                let e = s.rim_edge(end, i).ok_or(E::Unresolved)?;
                let b = s.source().profile_bounds(e.source_band, [0., 1.])?;
                (
                    e.vertex_indices,
                    I::new(b.radius_mm[0], b.radius_mm[1]),
                    I::new(b.slope[0], b.slope[1]),
                    I::new(b.second_derivative_per_mm[0], b.second_derivative_per_mm[1]),
                    e.point_error_mm,
                    I::point(s.source().spec().pitch_mm)
                        * I::point(f64::from(s.source().spec().starts)),
                )
            }
        };
        let rate = lead / tau();
        // Orthogonal radial/tangential components of the ACTUAL C''.
        let curvature =
            ((dd * rate.square() - r).square() + I::point(4.) * (m * rate).square()).sqrt();
        Ok((ids, upper(curvature)?, error[0].max(error[1])))
    }
    fn minimum(&self) -> f64 {
        match self {
            Self::Wood(s) => s.source().minimum_radius_mm(),
            Self::Metric(s) => s.source().minimum_radius_mm()[0],
        }
    }
}
pub(super) struct Angle {
    token: Box<[u8]>,
    pub bounds: I,
}
impl Angle {
    pub fn new(q: &Q, m: &mut Memory, w: &mut Work<'_>) -> R<Self> {
        let token = Q::encode(&[q], w)?;
        m.charge(token.len())?;
        Ok(Self {
            token,
            bounds: q.enclosure(w)?.1,
        })
    }
    pub fn exact(&self, w: &mut Work<'_>) -> R<Q> {
        Ok(Q::decode(&self.token, w)?)
    }
    pub fn cmp(&self, other: &Self, w: &mut Work<'_>) -> R<Ordering> {
        w.step()?;
        if self.bounds.hi < other.bounds.lo {
            return Ok(Ordering::Less);
        }
        if self.bounds.lo > other.bounds.hi {
            return Ok(Ordering::Greater);
        }
        Ok(self.exact(w)?.cmp(&other.exact(w)?, w)?)
    }
}
pub(super) struct Vertex {
    pub id: usize,
    pub angle: Angle,
}
pub(super) struct Ring {
    pub vertices: Vec<Vertex>,
    pub curve_error: f64,
    pub source_min: f64,
}
#[derive(Clone, Copy)]
struct Link {
    other: usize,
    curvature: f64,
    error: f64,
}
struct Node {
    angle: Angle,
    links: [Link; 2],
    degree: usize,
}
pub(super) fn read(
    wall: Wall<'_>,
    end: ShellEnd,
    offset: usize,
    m: &mut Memory,
    w: &mut Work<'_>,
) -> R<Ring> {
    let count = wall.mesh().vertices.len() / 3;
    let mut nodes: Vec<Option<Node>> = m.allocate(count, w)?;
    w.charge(count)?;
    nodes.resize_with(count, || None);
    for edge in 0..wall.count(end) {
        w.step()?;
        let ([a, b], curvature, error) = wall.edge(end, edge)?;
        insert(&mut nodes, (a, b, curvature, error), &wall, m, w)?;
    }
    order(nodes, offset, wall.minimum(), m, w)
}
fn order(
    mut nodes: Vec<Option<Node>>,
    offset: usize,
    source_min: f64,
    m: &mut Memory,
    w: &mut Work<'_>,
) -> R<Ring> {
    let unique = nodes.iter().filter(|n| n.is_some()).count();
    if unique < 3 {
        return Err(E::Unresolved);
    }
    let mut start: Option<usize> = None;
    for (i, node) in nodes.iter().enumerate() {
        w.step()?;
        if let Some(n) = node {
            if n.degree != 2
                || n.links[0].other == n.links[1].other
                || n.links
                    .iter()
                    .any(|l| l.other == i || l.other >= nodes.len() || nodes[l.other].is_none())
            {
                return Err(E::Unresolved);
            }
            if start.is_none()
                || n.angle
                    .cmp(&nodes[start.unwrap()].as_ref().unwrap().angle, w)?
                    == Ordering::Less
            {
                start = Some(i);
            }
        }
    }
    let start = start.unwrap();
    let mut current = start;
    let mut previous = usize::MAX;
    let mut result: Vec<Vertex> = m.allocate(unique, w)?;
    let mut curve_error = 0_f64;
    loop {
        w.step()?;
        let node = nodes[current].take().ok_or(E::Unresolved)?;
        let next = if previous == usize::MAX {
            let [a, b] = node.links.map(|l| l.other);
            if nodes[a]
                .as_ref()
                .unwrap()
                .angle
                .cmp(&nodes[b].as_ref().unwrap().angle, w)?
                == Ordering::Less
            {
                a
            } else {
                b
            }
        } else {
            node.links
                .iter()
                .find(|l| l.other != previous)
                .ok_or(E::Unresolved)?
                .other
        };
        let edge = node.links.iter().find(|l| l.other == next).unwrap();
        let a = node.angle.exact(w)?;
        let b = if next == start {
            result
                .first()
                .ok_or(E::Unresolved)?
                .angle
                .exact(w)?
                .add(&Q::integer(1, w)?, w)?
        } else {
            let next = nodes[next].as_ref().ok_or(E::Unresolved)?;
            if node.angle.cmp(&next.angle, w)? != Ordering::Less {
                return Err(E::Unresolved);
            }
            next.angle.exact(w)?
        };
        let width = b.sub(&a, w)?;
        if width.sign() != Ordering::Greater || width.cmp(&Q::float(0.5, w)?, w)? != Ordering::Less
        {
            return Err(E::Unresolved);
        }
        let theta = width.enclosure(w)?.1 * tau();
        curve_error = curve_error.max(
            (I::point(edge.curvature) * theta.square() / I::point(8.) + I::point(edge.error)).hi,
        );
        result.push(Vertex {
            id: current + offset,
            angle: node.angle,
        });
        if next == start {
            break;
        }
        // Each iteration consumes one live node. Reaching any earlier node
        // other than start fails through take/ok_or, so no extra cycle bound
        // or unbounded retry is needed.
        previous = current;
        current = next;
    }
    if result.len() != unique {
        return Err(E::Unresolved);
    }
    Ok(Ring {
        vertices: result,
        curve_error,
        source_min,
    })
}
pub(super) struct Graph {
    pub minimum: f64,
    pub maximum: f64,
    pub radial_error: f64,
}
/// A literal positive one-winding polygon is independently a radial graph.
pub(super) fn graph(ring: &Ring, mesh: &DisplayMesh, w: &mut Work<'_>) -> R<Graph> {
    if ring
        .vertices
        .iter()
        .any(|v| v.id >= mesh.vertices.len() / 3)
    {
        return Err(E::Unresolved);
    }
    let mut winding = 0;
    let mut minimum = f64::INFINITY;
    let mut maximum = 0_f64;
    let mut lipschitz = 0_f64;
    for i in 0..ring.vertices.len() {
        w.charge(48)?;
        let a = point(mesh, ring.vertices[i].id);
        let b = point(mesh, ring.vertices[(i + 1) % ring.vertices.len()].id);
        let determinant = det([a[0], a[1]], [b[0], b[1]], w)?;
        if determinant.lo <= 0. {
            return Err(E::Unresolved);
        }
        if a[1] <= 0. && b[1] > 0. {
            winding += 1;
        }
        let a = a.map(I::point);
        let b = b.map(I::point);
        let length = ((b[0] - a[0]).square() + (b[1] - a[1]).square()).sqrt();
        let radius = (a[0].square() + a[1].square())
            .sqrt()
            .hi
            .max((b[0].square() + b[1].square()).sqrt().hi);
        let low = (determinant / I::point(length.hi)).lo;
        let lip = (I::point(radius).square() * I::point(length.hi) / I::point(determinant.lo)).hi;
        minimum = minimum.min(low);
        maximum = maximum.max(radius);
        lipschitz = lipschitz.max(lip);
    }
    if winding != 1
        || minimum <= 0.
        || !lipschitz.is_finite()
        || ring.curve_error >= ring.source_min
    {
        return Err(E::Unresolved);
    }
    // Endpoint angle rounding may cross authored events. The GLOBAL polygon
    // Lipschitz bound repairs angular displacement; it assumes no band equality.
    let delta = I::point(ring.curve_error);
    let radius = I::point(ring.source_min);
    let error = delta + I::point(lipschitz) * delta / (radius.square() - delta.square()).sqrt();
    let radial_error = upper(error)?;
    Ok(Graph {
        minimum,
        maximum,
        radial_error,
    })
}
fn det(a: [f64; 2], b: [f64; 2], w: &mut Work<'_>) -> R<I> {
    let fast = I::point(a[0]) * I::point(b[1]) - I::point(a[1]) * I::point(b[0]);
    if fast.lo > 0. || fast.hi < 0. {
        return Ok(fast);
    }
    let a = a.map(|x| Q::float(x, w));
    let b = b.map(|x| Q::float(x, w));
    let [ax, ay] = a;
    let [bx, by] = b;
    let exact = ax?.mul(&by?, w)?.sub(&ay?.mul(&bx?, w)?, w)?;
    Ok(exact.enclosure(w)?.1)
}
pub(super) fn area(points: [[f64; 3]; 3], w: &mut Work<'_>) -> R<Ordering> {
    let a = points[0];
    let b = points[1];
    let c = points[2];
    let fast = (I::point(b[0]) - I::point(a[0])) * (I::point(c[1]) - I::point(a[1]))
        - (I::point(b[1]) - I::point(a[1])) * (I::point(c[0]) - I::point(a[0]));
    if fast.lo > 0. {
        return Ok(Ordering::Greater);
    }
    if fast.hi < 0. {
        return Ok(Ordering::Less);
    }
    let dx = Q::float(b[0], w)?.sub(&Q::float(a[0], w)?, w)?;
    let dy = Q::float(b[1], w)?.sub(&Q::float(a[1], w)?, w)?;
    let ex = Q::float(c[0], w)?.sub(&Q::float(a[0], w)?, w)?;
    let ey = Q::float(c[1], w)?.sub(&Q::float(a[1], w)?, w)?;
    Ok(dx.mul(&ey, w)?.sub(&dy.mul(&ex, w)?, w)?.sign())
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn literal_radial_graph_rejects_folded_collapsed_double_winding_and_uncertain_rings() {
        let (_, mesh, mut ring, _) = super::super::caps::tests::circles();
        let mut w = Work {
            remaining: 16_000_000,
            used: 0,
            interrupt: None,
        };
        let g = graph(&ring, &mesh, &mut w).unwrap();
        assert!(g.minimum > 1.4 && g.maximum < 1.6 && g.radial_error < 0.1);
        let mut folded = mesh.clone();
        let id = ring.vertices[1].id;
        folded.vertices[id * 3] = -folded.vertices[id * 3];
        assert!(graph(&ring, &folded, &mut w).is_err());
        let mut collapsed = mesh.clone();
        let a = ring.vertices[0].id;
        let b = ring.vertices[1].id;
        let p = point(&collapsed, a);
        collapsed.vertices[b * 3..b * 3 + 3].copy_from_slice(&p);
        assert!(graph(&ring, &collapsed, &mut w).is_err());
        let saved = ring.curve_error;
        ring.curve_error = ring.source_min;
        assert!(graph(&ring, &mesh, &mut w).is_err());
        ring.curve_error = saved;
        let mut m = Memory { used: SCRATCH };
        let extra = ring
            .vertices
            .iter()
            .map(|v| {
                let angle = Angle::new(&v.angle.exact(&mut w).unwrap(), &mut m, &mut w).unwrap();
                Vertex { id: v.id, angle }
            })
            .collect::<Vec<_>>();
        ring.vertices.extend(extra);
        assert!(graph(&ring, &mesh, &mut w).is_err());
        assert_eq!(
            area([[0., 0., 0.], [1., 1., 0.], [2., 2., 0.]], &mut w).unwrap(),
            Ordering::Equal
        );
        let delta = f64::EPSILON;
        assert_eq!(
            area(
                [[1., 1., 0.], [2., 2., 0.], [2., 2. + 2. * delta, 0.]],
                &mut w
            )
            .unwrap(),
            Ordering::Greater
        );
        assert!(det([1., 1.], [1., 1.], &mut w).unwrap().lo == 0.);
        let near = det([1., 1.], [1., 1_f64.next_up()], &mut w).unwrap();
        assert!(near.lo > 0.);
        let last = ring.vertices[0].id;
        ring.vertices[0].id = usize::MAX;
        assert!(graph(&ring, &mesh, &mut w).is_err());
        ring.vertices[0].id = last;
    }
}

fn endpoints(a: usize, b: usize, count: usize) -> R<()> {
    if a == b || a >= count || b >= count {
        Err(E::Unresolved)
    } else {
        Ok(())
    }
}
fn link(node: &mut Node, other: usize, curvature: f64, error: f64) -> R<()> {
    if node.degree == 2 || node.links.iter().any(|l| l.other == other) {
        return Err(E::Unresolved);
    }
    node.links[node.degree] = Link {
        other,
        curvature,
        error,
    };
    node.degree += 1;
    Ok(())
}
fn insert(
    nodes: &mut [Option<Node>],
    edge: (usize, usize, f64, f64),
    wall: &Wall<'_>,
    m: &mut Memory,
    w: &mut Work<'_>,
) -> R<()> {
    let (a, b, curvature, error) = edge;
    endpoints(a, b, nodes.len())?;
    for (v, other) in [(a, b), (b, a)] {
        if nodes[v].is_none() {
            let q = Q::from_rational(wall.angle(v, w)?, w)?;
            let angle = Angle::new(&q, m, w)?;
            nodes[v] = Some(Node {
                angle,
                links: [Link {
                    other: usize::MAX,
                    curvature: 0.,
                    error: 0.,
                }; 2],
                degree: 0,
            });
        }
        link(nodes[v].as_mut().unwrap(), other, curvature, error)?;
    }
    Ok(())
}
#[cfg(test)]
mod cycle_tests {
    use super::*;
    fn nodes(angles: &[f64], m: &mut Memory, w: &mut Work<'_>) -> Vec<Option<Node>> {
        (0..angles.len())
            .map(|i| {
                Some(Node {
                    angle: Angle::new(&Q::float(angles[i], w).unwrap(), m, w).unwrap(),
                    links: [
                        Link {
                            other: (i + angles.len() - 1) % angles.len(),
                            curvature: 1.,
                            error: 0.,
                        },
                        Link {
                            other: (i + 1) % angles.len(),
                            curvature: 1.,
                            error: 0.,
                        },
                    ],
                    degree: 2,
                })
            })
            .collect()
    }
    #[test]
    fn exact_source_cycle_contract_refuses_missing_duplicate_crossed_and_wide_edges() {
        let mut w = Work {
            remaining: 16_000_000,
            used: 0,
            interrupt: None,
        };
        let mut m = Memory { used: SCRATCH };
        assert!(
            order(
                nodes(&[0., 0.25, 0.5, 0.75], &mut m, &mut w),
                0,
                1.,
                &mut m,
                &mut w
            )
            .is_ok()
        );
        assert!(order(nodes(&[0., 0.1], &mut m, &mut w), 0, 1., &mut m, &mut w).is_err());
        let mut bad = nodes(&[0., 0.25, 0.5, 0.75], &mut m, &mut w);
        bad[1].as_mut().unwrap().degree = 1;
        assert!(order(bad, 0, 1., &mut m, &mut w).is_err());
        assert!(
            order(
                nodes(&[0., 0.5, 0.25, 0.75], &mut m, &mut w),
                0,
                1.,
                &mut m,
                &mut w
            )
            .is_err()
        );
        assert!(
            order(
                nodes(&[0., 0.1, 0.2], &mut m, &mut w),
                0,
                1.,
                &mut m,
                &mut w
            )
            .is_err()
        );
        let mut bad = nodes(&[0., 0.2, 0.4, 0.6, 0.8, 0.9], &mut m, &mut w);
        for group in [[0, 1, 2], [3, 4, 5]] {
            for i in 0..3 {
                let n = bad[group[i]].as_mut().unwrap();
                n.links[0].other = group[(i + 2) % 3];
                n.links[1].other = group[(i + 1) % 3];
            }
        }
        assert!(order(bad, 0, 1., &mut m, &mut w).is_err());
        // A strictly short first edge followed by an actual backwards station.
        assert!(
            order(
                nodes(&[0., 0.25, 0.6, 0.5, 0.9], &mut m, &mut w),
                0,
                1.,
                &mut m,
                &mut w
            )
            .is_err()
        );
        // Both cycles have valid sub-half-turn edges; only the first is reached.
        let mut split = nodes(&[0., 0.3, 0.6, 0.7, 0.8, 0.9], &mut m, &mut w);
        for group in [[0, 1, 2], [3, 4, 5]] {
            for i in 0..3 {
                let n = split[group[i]].as_mut().unwrap();
                n.links[0].other = group[(i + 2) % 3];
                n.links[1].other = group[(i + 1) % 3];
            }
        }
        assert!(order(split, 0, 1., &mut m, &mut w).is_err());
        assert!(endpoints(0, 0, 4).is_err());
        assert!(endpoints(0, 4, 4).is_err());
        let mut data = nodes(&[0., 0.25, 0.5, 0.75], &mut m, &mut w);
        let n = data[0].as_mut().unwrap();
        assert!(link(n, 1, 1., 0.).is_err());
        n.degree = 1;
        n.links[0].other = 1;
        assert!(link(n, 1, 1., 0.).is_err());
        let a = Angle::new(&Q::float(0.25, &mut w).unwrap(), &mut m, &mut w).unwrap();
        let b = Angle::new(&Q::float(0.5, &mut w).unwrap(), &mut m, &mut w).unwrap();
        assert_eq!(b.cmp(&a, &mut w).unwrap(), Ordering::Greater);
        assert_eq!(a.cmp(&a, &mut w).unwrap(), Ordering::Equal);
    }
}
