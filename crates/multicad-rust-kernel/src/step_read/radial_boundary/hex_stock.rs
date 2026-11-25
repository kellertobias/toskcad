//! Hex prism stock joined to one retained radial boundary: a pipe nut (stock
//! around a through cavity, both caps become annuli) or a pipe hex head
//! (stock on the thread's large-end cap, which becomes the shoulder annulus).
//! The radial tables are decoded exactly as their own faces; the stock is
//! decoded from straight edges only. This is a candidate, not an admission.
use super::*;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HexStockKind {
    /// Hex prism around a through cavity; the cavity caps become annuli.
    Nut,
    /// Hex prism on the thread's cap at `span_mm[0]`, which becomes the shoulder.
    Head,
}
/// Finite regular hex prism: axis point at local height zero, unit axis, unit
/// reference through the first corner, across-flats width and axial span.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HexStock {
    pub kind: HexStockKind,
    pub origin: Point,
    pub normal: Point,
    pub reference: Point,
    pub across_flats_mm: f64,
    pub span_mm: [f64; 2],
}
fn dot(a: Point, b: Point) -> f64 {
    (0..3).map(|i| a[i] * b[i]).sum()
}
fn cross(a: Point, b: Point) -> Point {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
impl HexStock {
    /// Lower then upper corner rings, each counterclockwise about `normal`,
    /// starting at the corner on `reference`. `None` for an invalid stock.
    pub fn corners(&self) -> Option<[[Point; 6]; 2]> {
        let Self {
            origin,
            normal,
            reference,
            across_flats_mm: af,
            span_mm: [low, high],
            ..
        } = *self;
        let finite = origin
            .into_iter()
            .chain(normal)
            .chain(reference)
            .chain([af, low, high])
            .all(f64::is_finite);
        if !finite
            || af <= 0.
            || low >= high
            || (dot(normal, normal) - 1.).abs() > 1e-12
            || (dot(reference, reference) - 1.).abs() > 1e-12
            || dot(normal, reference).abs() > 1e-12
        {
            return None;
        }
        let binormal = cross(normal, reference);
        let half = af / 2.;
        let h = half / 3f64.sqrt();
        let polygon = [
            [2. * h, 0.],
            [h, half],
            [-h, half],
            [-2. * h, 0.],
            [-h, -half],
            [h, -half],
        ];
        let ring = |z: f64| {
            polygon.map(|[x, y]| {
                std::array::from_fn(|k| {
                    origin[k] + x * reference[k] + y * binormal[k] + z * normal[k]
                })
            })
        };
        let rings = [ring(low), ring(high)];
        rings
            .iter()
            .flatten()
            .flatten()
            .all(|v| v.is_finite())
            .then_some(rings)
    }
}
/// Decoded candidate only: the radial ownership has no material certificate.
#[derive(Clone)]
pub struct ImportedPipeHexHardware {
    pub boundary: radial::BoundaryOwnership,
    pub kind: HexStockKind,
    /// The twelve distinct stock corners, in first-use order.
    pub corners: Vec<Point>,
}
/// Decode one candidate. This is not an occupied-solid admission.
pub fn read_pipe_hex_hardware(input: &str, tolerance: f64) -> Result<ImportedPipeHexHardware> {
    let mut members = compound::read_members_by(input, tolerance, true, true, |shell, id, _| {
        recognize(shell, id)?.ok_or(ReadError::Unsupported)
    })?;
    members.pop().ok_or(ReadError::Invalid)
}
fn check(condition: bool) -> Result<()> {
    if condition {
        Ok(())
    } else {
        Err(ReadError::Unsupported)
    }
}
/// Whether every edge of this oriented ring lies on a LINE.
fn straight(shell: &Shell<'_, '_>, edges: &[(u32, bool)]) -> Result<bool> {
    for &(edge, _) in edges {
        let [Value::String(_), _, _, curve, _] = shell.args(edge, "EDGE_CURVE")? else {
            return Err(ReadError::Invalid);
        };
        if shell.simple(reference(curve)?)?.name != "LINE" {
            return Ok(false);
        }
    }
    Ok(true)
}
/// Start vertices of an oriented ring, after checking it is a closed chain.
fn ring_vertices(shell: &Shell<'_, '_>, edges: &[(u32, bool)]) -> Result<Vec<u32>> {
    let mut ends = Vec::with_capacity(edges.len());
    for &(edge, forward) in edges {
        let [Value::String(_), a, b, _, _] = shell.args(edge, "EDGE_CURVE")? else {
            return Err(ReadError::Invalid);
        };
        let [a, b] = [reference(a)?, reference(b)?];
        ends.push(if forward { [a, b] } else { [b, a] });
    }
    for i in 0..ends.len() {
        check(ends[i][1] == ends[(i + 1) % ends.len()][0])?;
    }
    Ok(ends.iter().map(|e| e[0]).collect())
}
struct Planar {
    origin: Point,
    outward: Point,
}
fn planar(shell: &Shell<'_, '_>, surface: u32, forward: bool) -> Result<Planar> {
    check(shell.simple(surface)?.name == "PLANE")?;
    let [Value::String(_), placement] = shell.args(surface, "PLANE")? else {
        return Err(ReadError::Invalid);
    };
    let a = axis(shell, reference(placement)?)?;
    check((dot(a.normal, a.normal) - 1.).abs() <= 1e-12)?;
    Ok(Planar {
        origin: a.origin,
        outward: if forward {
            a.normal
        } else {
            a.normal.map(|v| -v)
        },
    })
}
enum Role {
    /// Planar single-loop stock face bounded by straight edges.
    Stock(Vec<u32>),
    /// Planar face with a straight outer ring and one radial hole.
    Annulus(Vec<u32>, Vec<(u32, bool)>),
    /// Every other face belongs to the radial boundary.
    Radial(Vec<(u32, bool)>),
}
/// `None` unless the member is one radial boundary plus a straight hex prism
/// in the nut or head arrangement; other members keep their own readers.
pub(in crate::step_read) fn recognize(
    shell: &mut Shell<'_, '_>,
    id: u32,
) -> Result<Option<ImportedPipeHexHardware>> {
    if shell.simple(id)?.name != "MANIFOLD_SOLID_BREP" {
        return Ok(None);
    }
    let faces = shell.member_faces(id)?;
    if faces.len() > 64 {
        return Ok(None);
    }
    let mut roles = Vec::with_capacity(faces.len());
    let (mut stock, mut annuli) = (0usize, 0usize);
    for &face_id in &faces {
        let (surface, forward, bounds) = face(shell, face_id)?;
        let role = match bounds[..] {
            [single] => {
                let edges = ring(shell, single, true)?;
                if shell.simple(surface)?.name == "PLANE" && straight(shell, &edges)? {
                    stock += 1;
                    Role::Stock(ring_vertices(shell, &edges)?)
                } else {
                    Role::Radial(edges)
                }
            }
            [a, b] => {
                let (outer_first, _, _) = shell.bound(a, false)?;
                let (outer, hole) = if outer_first { (a, b) } else { (b, a) };
                let outer = ring(shell, outer, true)?;
                let hole = ring(shell, hole, false)?;
                if shell.simple(surface)?.name != "PLANE"
                    || !straight(shell, &outer)?
                    || straight(shell, &hole)?
                {
                    return Ok(None);
                }
                annuli += 1;
                Role::Annulus(ring_vertices(shell, &outer)?, hole)
            }
            _ => return Ok(None),
        };
        roles.push((surface, forward, role));
    }
    let kind = match (stock, annuli) {
        (6, 2) => HexStockKind::Nut,
        (7, 1) => HexStockKind::Head,
        _ => return Ok(None),
    };
    let tolerance = shell.tolerance;
    let mut corners = Vec::<(u32, Point)>::new();
    let mut rings = Vec::new();
    let mut decoder = Decoder::default();
    for (surface, forward, role) in roles {
        let stock_ring = match &role {
            Role::Stock(v) | Role::Annulus(v, _) => Some(v.clone()),
            Role::Radial(_) => None,
        };
        if let Some(vertices) = stock_ring {
            let plane = planar(shell, surface, forward)?;
            let mut points = Vec::with_capacity(vertices.len());
            for &v in &vertices {
                let p = shell.vertex(v)?;
                let offset = std::array::from_fn(|k| p[k] - plane.origin[k]);
                check(dot(offset, plane.outward).abs() <= tolerance)?;
                if !corners.iter().any(|(id, _)| *id == v) {
                    corners.push((v, p));
                }
                points.push(p);
            }
            // Every stock ring is counterclockwise about its outward normal.
            let n = points.len();
            let mut newell = [0.; 3];
            for i in 0..n {
                let c = cross(points[i], points[(i + 1) % n]);
                for k in 0..3 {
                    newell[k] += c[k];
                }
            }
            check(dot(newell, plane.outward) > 0.)?;
            rings.push(vertices);
        }
        match role {
            Role::Stock(_) => {}
            Role::Annulus(_, hole) => match kind {
                // Inverse of the writer: a nut annulus keeps the cavity cap's
                // sense and reverses its wire; a head shoulder flips the cap's
                // sense and keeps the wire (material is on the other side).
                HexStockKind::Nut => {
                    let restored = hole.iter().rev().map(|&(e, d)| (e, !d)).collect::<Vec<_>>();
                    decoder.push(shell, surface, forward, &restored)?;
                }
                HexStockKind::Head => decoder.push(shell, surface, !forward, &hole)?,
            },
            Role::Radial(edges) => match kind {
                HexStockKind::Nut => {
                    let restored = edges
                        .iter()
                        .rev()
                        .map(|&(e, d)| (e, !d))
                        .collect::<Vec<_>>();
                    decoder.push(shell, surface, !forward, &restored)?;
                }
                HexStockKind::Head => decoder.push(shell, surface, forward, &edges)?,
            },
        }
    }
    check(corners.len() == 12)?;
    // Two hexagonal rings (annuli or the head top) and six quadrilateral
    // walls, each wall joining one edge of each ring.
    let (hexes, walls): (Vec<_>, Vec<_>) = rings.into_iter().partition(|r| r.len() == 6);
    check(hexes.len() == 2 && walls.len() == 6 && walls.iter().all(|w| w.len() == 4))?;
    check(hexes[0].iter().all(|v| !hexes[1].contains(v)))?;
    for wall in &walls {
        let on = |ring: &Vec<u32>| wall.iter().filter(|v| ring.contains(v)).count();
        check(on(&hexes[0]) == 2 && on(&hexes[1]) == 2)?;
    }
    let boundary =
        radial::BoundaryOwnership::new(decoder.input, tolerance, None).map_err(failure)?;
    Ok(Some(ImportedPipeHexHardware {
        boundary,
        kind,
        corners: corners.into_iter().map(|(_, p)| p).collect(),
    }))
}
