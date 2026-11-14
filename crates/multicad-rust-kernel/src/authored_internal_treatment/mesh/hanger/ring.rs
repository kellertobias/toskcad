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
        if !curvature.valid() {
            return Err(E::Unresolved);
        }
        Ok((ids, curvature.hi, error[0].max(error[1])))
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
        if a == b || a >= count || b >= count {
            return Err(E::Unresolved);
        }
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
            let node = nodes[v].as_mut().unwrap();
            if node.degree == 2 || node.links.iter().any(|l| l.other == other) {
                return Err(E::Unresolved);
            }
            node.links[node.degree] = Link {
                other,
                curvature,
                error,
            };
            node.degree += 1;
        }
    }
    let unique = nodes.iter().filter(|n| n.is_some()).count();
    if unique < 3 {
        return Err(E::Unresolved);
    }
    let mut start: Option<usize> = None;
    for (i, node) in nodes.iter().enumerate() {
        w.step()?;
        if let Some(n) = node {
            if n.degree != 2 {
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
        if result.len() >= unique {
            return Err(E::Unresolved);
        }
        previous = current;
        current = next;
    }
    if result.len() != unique {
        return Err(E::Unresolved);
    }
    Ok(Ring {
        vertices: result,
        curve_error,
        source_min: wall.minimum(),
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
    if !error.valid() {
        return Err(E::Unresolved);
    }
    Ok(Graph {
        minimum,
        maximum,
        radial_error: error.hi,
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
