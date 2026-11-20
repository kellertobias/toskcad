//! Replay a smooth convex loft from its planar caps, C¹ cubic rails and the
//! degree (1, 3) walls that share the rail poles.
use super::*;
use crate::planar_boolean::SmoothLoft;

#[derive(Debug)]
pub struct ImportedSmoothLoft {
    pub geometry: SmoothLoft,
    pub origin_mm: Point,
    pub u: Point,
    pub v: Point,
}

struct Use {
    edge: u32,
    start: u32,
    end: u32,
}

struct Edge {
    vertices: [u32; 2],
    rail: Option<SplineEdge>,
    uses: Vec<(u32, bool)>,
}

struct Reader<'s, 'a, 'input> {
    shell: &'s mut Shell<'a, 'input>,
    edges: BTreeMap<u32, Edge>,
}

impl Reader<'_, '_, '_> {
    /// The oriented loop of one single-bound face, with every edge decoded.
    fn face(&mut self, face: u32) -> Result<Option<(u32, bool, Vec<Use>)>> {
        let [Value::String(_), bounds, surface, forward] =
            self.shell.args(face, "ADVANCED_FACE")?
        else {
            return Err(ReadError::Invalid);
        };
        let [bound] = refs(bounds)?
            .try_into()
            .map_err(|_| ReadError::Unsupported)?;
        let (surface, forward) = (reference(surface)?, sense(forward)?);
        let (_, wire, loop_forward) = self.shell.bound(bound, true)?;
        let [Value::String(_), oriented] = self.shell.args(wire, "EDGE_LOOP")? else {
            return Err(ReadError::Invalid);
        };
        let mut oriented = refs(oriented)?;
        if oriented.len() > 66 {
            return Err(ReadError::Limit);
        }
        if !loop_forward {
            oriented.reverse();
        }
        let mut uses = Vec::with_capacity(oriented.len());
        for id in oriented {
            let [
                Value::String(_),
                Value::Derived,
                Value::Derived,
                edge,
                orientation,
            ] = self.shell.args(id, "ORIENTED_EDGE")?
            else {
                return Err(ReadError::Invalid);
            };
            let edge = reference(edge)?;
            let along = sense(orientation)? == loop_forward;
            let vertices = self.edge(edge)?;
            let [start, end] = if along {
                vertices
            } else {
                [vertices[1], vertices[0]]
            };
            if uses.last().is_some_and(|last: &Use| last.end != start) {
                return Err(ReadError::Invalid);
            }
            let record = self.edges.get_mut(&edge).ok_or(ReadError::Invalid)?;
            if record.uses.len() == 2 {
                return Err(ReadError::Invalid);
            }
            record.uses.push((face, along));
            uses.push(Use { edge, start, end });
        }
        if uses.len() < 3 || uses.first().map(|u| u.start) != uses.last().map(|u| u.end) {
            return Err(ReadError::Invalid);
        }
        Ok(Some((surface, forward, uses)))
    }

    fn edge(&mut self, id: u32) -> Result<[u32; 2]> {
        if let Some(edge) = self.edges.get(&id) {
            return Ok(edge.vertices);
        }
        let [Value::String(_), start, end, curve, forward] = self.shell.args(id, "EDGE_CURVE")?
        else {
            return Err(ReadError::Invalid);
        };
        let vertices = [reference(start)?, reference(end)?];
        let forward = sense(forward)?;
        let curve = self.shell.edge_geometry(reference(curve)?, vertices)?;
        let rail = match self.shell.simple(curve)?.name {
            "LINE" => {
                self.shell.line_curve(curve, vertices, forward)?;
                None
            }
            "B_SPLINE_CURVE_WITH_KNOTS" if forward => {
                let rail = self.shell.spline_curve(curve)?;
                let [a, b] = [
                    self.shell.vertex(vertices[0])?,
                    self.shell.vertex(vertices[1])?,
                ];
                let poles = rail.poles_mm();
                if length(sub(a, poles[0])) > self.shell.tolerance
                    || length(sub(b, poles[poles.len() - 1])) > self.shell.tolerance
                {
                    return Err(ReadError::Invalid);
                }
                Some(rail)
            }
            _ => return Err(ReadError::Unsupported),
        };
        self.edges.insert(
            id,
            Edge {
                vertices,
                rail,
                uses: Vec::with_capacity(2),
            },
        );
        Ok(vertices)
    }
}

/// The two pole rows of a degree (1, 3) wall surface and its loft knots.
type WallRows = ([Vec<Point>; 2], Vec<f64>);

fn wall_rows(shell: &Shell<'_, '_>, id: u32) -> Result<Option<WallRows>> {
    use Value::{Enum, List, Number, String};
    let entity = shell.simple(id)?;
    if entity.name != "B_SPLINE_SURFACE_WITH_KNOTS"
        || !matches!(entity.args.get(1..3), Some([Number(1.), Number(3.)]))
    {
        return Ok(None);
    }
    let [
        String(_),
        Number(1.),
        Number(3.),
        List(net),
        Enum(_),
        Enum("F"),
        Enum("F"),
        Enum("F"),
        List(u_mult),
        List(v_mult),
        List(u_knots),
        List(v_knots),
        Enum(_),
    ] = entity.args.as_slice()
    else {
        return Err(ReadError::Invalid);
    };
    let [Number(u0), Number(u1)] = u_knots.as_slice() else {
        return Err(ReadError::Invalid);
    };
    if u_mult.as_slice() != [Number(2.), Number(2.)] || u0 >= u1 {
        return Err(ReadError::Invalid);
    }
    let knots = v_knots
        .iter()
        .map(|value| match value {
            Number(knot) if knot.is_finite() => Ok(*knot),
            _ => Err(ReadError::Invalid),
        })
        .collect::<Result<Vec<_>>>()?;
    let stations = knots.len();
    if !(2..=16).contains(&stations)
        || v_mult.len() != stations
        || v_mult.iter().enumerate().any(|(i, m)| {
            let expected = if i == 0 || i + 1 == stations { 4. } else { 2. };
            !matches!(m, Number(m) if *m == expected)
        })
    {
        return Err(ReadError::Invalid);
    }
    let [List(first), List(second)] = net.as_slice() else {
        return Err(ReadError::Invalid);
    };
    let row = |row: &[Value<'_>]| -> Result<Vec<Point>> {
        if row.len() != 2 * stations {
            return Err(ReadError::Invalid);
        }
        row.iter()
            .map(|pole| shell.point(reference(pole)?))
            .collect()
    };
    Ok(Some(([row(first)?, row(second)?], knots)))
}

pub(super) fn read_member(
    shell: &mut Shell<'_, '_>,
    faces: &[u32],
) -> Result<Option<(Readback, ImportedSmoothLoft)>> {
    if !(5..=66).contains(&faces.len()) {
        return Ok(None);
    }
    // Recognize the shape by its surfaces before decoding any edge.
    let mut caps = Vec::with_capacity(2);
    let mut walls = Vec::with_capacity(faces.len() - 2);
    for &face in faces {
        let [Value::String(_), _, surface, _] = shell.args(face, "ADVANCED_FACE")? else {
            return Err(ReadError::Invalid);
        };
        let surface = reference(surface)?;
        match shell.simple(surface)?.name {
            "PLANE" => caps.push(face),
            "B_SPLINE_SURFACE_WITH_KNOTS" if wall_rows(shell, surface)?.is_some() => {
                walls.push(face)
            }
            _ => return Ok(None),
        }
    }
    let [first_cap, last_cap] = caps.as_slice() else {
        return Ok(None);
    };
    let count = walls.len();
    if faces.iter().collect::<BTreeSet<_>>().len() != faces.len() {
        return Err(ReadError::Invalid);
    }
    let mut reader = Reader {
        shell,
        edges: BTreeMap::new(),
    };
    let mut cap = |face: u32| -> Result<(Plane, Vec<u32>)> {
        let (surface, forward, uses) = reader.face(face)?.ok_or(ReadError::Invalid)?;
        let [Value::String(_), placement] = reader.shell.args(surface, "PLANE")? else {
            return Err(ReadError::Invalid);
        };
        let mut plane = reader.shell.placement(reference(placement)?)?;
        if !forward {
            plane.normal = plane.normal.map(|x| -x);
        }
        if uses.len() != count || uses.iter().any(|u| reader.edges[&u.edge].rail.is_some()) {
            return Err(ReadError::Invalid);
        }
        let vertices: Vec<_> = uses.iter().map(|u| u.start).collect();
        let points = vertices
            .iter()
            .map(|&id| reader.shell.vertex(id))
            .collect::<Result<Vec<_>>>()?;
        let anchor = points[0];
        let area: f64 = points
            .windows(2)
            .map(|pair| {
                dot(
                    plane.normal,
                    cross(sub(pair[0], anchor), sub(pair[1], anchor)),
                )
            })
            .sum();
        if area <= 0.
            || points
                .iter()
                .any(|&p| dot(sub(p, plane.origin), plane.normal).abs() > reader.shell.tolerance)
        {
            return Err(ReadError::Invalid);
        }
        Ok((plane, vertices))
    };
    let (first, first_ring) = cap(*first_cap)?;
    let (last, last_ring) = cap(*last_cap)?;
    let first_ids: BTreeSet<_> = first_ring.iter().copied().collect();
    let last_ids: BTreeSet<_> = last_ring.iter().copied().collect();
    if first_ids.len() != count || last_ids.len() != count || !first_ids.is_disjoint(&last_ids) {
        return Err(ReadError::Invalid);
    }
    // Map each lower vertex to its rail, and keep each wall's surface rows.
    let mut rails = BTreeMap::new();
    let mut surfaces = Vec::with_capacity(count);
    for &face in &walls {
        let (surface, _, uses) = reader.face(face)?.ok_or(ReadError::Invalid)?;
        if uses.len() != 4 {
            return Err(ReadError::Invalid);
        }
        let mut lower = 0;
        let mut upper = 0;
        for u in &uses {
            let edge = &reader.edges[&u.edge];
            match &edge.rail {
                Some(_) => {
                    let [a, b] = edge.vertices;
                    if !first_ids.contains(&a) || !last_ids.contains(&b) {
                        return Err(ReadError::Invalid);
                    }
                    if let Some(previous) = rails.insert(a, (b, u.edge))
                        && previous != (b, u.edge)
                    {
                        return Err(ReadError::Invalid);
                    }
                }
                None => {
                    let ends = [u.start, u.end];
                    if ends.iter().all(|id| first_ids.contains(id)) {
                        lower += 1;
                    } else if ends.iter().all(|id| last_ids.contains(id)) {
                        upper += 1;
                    } else {
                        return Err(ReadError::Invalid);
                    }
                }
            }
        }
        if lower != 1 || upper != 1 {
            return Err(ReadError::Invalid);
        }
        let rows = wall_rows(reader.shell, surface)?.ok_or(ReadError::Invalid)?;
        let rails_here: Vec<_> = uses
            .iter()
            .filter(|u| reader.edges[&u.edge].rail.is_some())
            .map(|u| u.edge)
            .collect();
        surfaces.push((rows, rails_here));
    }
    if rails.len() != count
        || reader.edges.len() != 3 * count
        || reader
            .edges
            .values()
            .any(|edge| !matches!(edge.uses.as_slice(), [a, b] if a.1 != b.1 && a.0 != b.0))
    {
        return Err(ReadError::Invalid);
    }
    let Reader { shell, edges } = reader;
    // The first outward cap runs clockwise along the loft; reverse it.
    let start: Vec<_> = first_ring.iter().rev().copied().collect();
    let end: Vec<_> = start.iter().map(|id| rails[id].0).collect();
    let same_cycle = last_ring
        .iter()
        .position(|id| *id == end[0])
        .is_some_and(|offset| {
            end.iter()
                .enumerate()
                .all(|(i, id)| last_ring[(i + offset) % count] == *id)
        });
    if !same_cycle {
        return Err(ReadError::Invalid);
    }
    let origin_mm = shell.points[&start[0]];
    let normal = first.normal.map(|x| -x);
    let u = first.u;
    let v = cross(normal, u);
    let height = dot(sub(shell.points[&end[0]], origin_mm), normal);
    if dot(normal, last.normal) < 0.
        || length(cross(normal, last.normal)) > 1e-10
        || height <= shell.tolerance
    {
        return Err(ReadError::Invalid);
    }
    let curves: Vec<&SplineEdge> = start
        .iter()
        .map(|id| edges[&rails[id].1].rail.as_ref().expect("rail edge"))
        .collect();
    let knots = curves[0].knots().to_vec();
    if curves.iter().any(|curve| curve.knots() != knots.as_slice()) {
        return Err(ReadError::Invalid);
    }
    // Rail heights are an affine image of their knots.
    let span = knots[knots.len() - 1] - knots[0];
    let heights: Vec<_> = knots
        .iter()
        .map(|knot| (knot - knots[0]) / span * height)
        .collect();
    let mut expected = Vec::with_capacity(2 * heights.len());
    expected.push(0.);
    for pair in heights.windows(2) {
        let step = pair[1] - pair[0];
        expected.extend([pair[0] + step / 3., pair[1] - step / 3.]);
    }
    expected.push(height);
    let mut local = Vec::with_capacity(count);
    for curve in &curves {
        let poles = curve.poles_mm();
        if poles.len() != expected.len() {
            return Err(ReadError::Invalid);
        }
        let mut rail = Vec::with_capacity(poles.len());
        for (&pole, &z) in poles.iter().zip(&expected) {
            let delta = sub(pole, origin_mm);
            if (dot(delta, normal) - z).abs() > shell.tolerance {
                return Err(ReadError::Invalid);
            }
            rail.push([dot(delta, u), dot(delta, v)]);
        }
        local.push(rail);
    }
    // Every wall's rows are the poles of the two rails it is bounded by.
    let near = |a: &[Point], b: &[Point]| {
        a.len() == b.len()
            && a.iter()
                .zip(b)
                .all(|(a, b)| length(sub(*a, *b)) <= shell.tolerance)
    };
    for (([first_row, second_row], wall_knots), wall_rails) in &surfaces {
        let [a, b] = wall_rails.as_slice() else {
            return Err(ReadError::Invalid);
        };
        let [a, b] = [a, b].map(|edge| edges[edge].rail.as_ref().expect("rail edge").poles_mm());
        if wall_knots != &knots
            || !((near(first_row, a) && near(second_row, b))
                || (near(first_row, b) && near(second_row, a)))
        {
            return Err(ReadError::Invalid);
        }
    }
    let geometry =
        SmoothLoft::from_rails(local, heights, shell.tolerance).map_err(|_| ReadError::Invalid)?;
    // Rail poles bound the walls by the convex hull property.
    let mut bounds_mm = [origin_mm; 2];
    for &point in shell
        .points
        .values()
        .chain(curves.iter().flat_map(|curve| curve.poles_mm()))
    {
        include_point(&mut bounds_mm, point);
    }
    let proof = Readback {
        bounds_mm,
        volume_mm3: geometry.volume_mm3(),
        face_count: count + 2,
        edges_mm: edges
            .values()
            .filter(|edge| edge.rail.is_none())
            .map(|edge| edge.vertices.map(|id| shell.points[&id]))
            .collect(),
        circles: Vec::new(),
        semicircles: Vec::new(),
        arcs: Vec::new(),
        ellipses: Vec::new(),
        // Rails in loft corner order, each from the lower cap to the upper.
        splines: curves.into_iter().cloned().collect(),
        solid_count: 1,
    };
    Ok(Some((
        proof,
        ImportedSmoothLoft {
            geometry,
            origin_mm,
            u,
            v,
        },
    )))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::step::write_smooth_loft;

    const TOL: f64 = 1e-7;

    fn rect(w: f64, h: f64, [x, y]: [f64; 2]) -> Vec<[f64; 2]> {
        vec![[x, y], [x + w, y], [x + w, y + h], [x, y + h]]
    }

    fn fixture(clockwise: bool) -> (SmoothLoft, String) {
        let mut sections = vec![
            rect(10., 6., [0., 0.]),
            rect(14., 8., [-2., -1.]),
            rect(12., 7., [1., 0.5]),
            rect(10., 6., [0., 0.]),
        ];
        if clockwise {
            sections.iter_mut().for_each(|ring| ring.reverse());
        }
        let loft = SmoothLoft::interpolate(&sections, &[0., 10., 22., 30.], TOL).unwrap();
        let step =
            write_smooth_loft(&loft, [11., -4., 9.], [0.8, 0.6, 0.], [-0.36, 0.48, 0.8]).unwrap();
        (loft, step)
    }

    #[test]
    fn replays_rails_walls_and_volume_from_exact_step_poles() {
        for clockwise in [false, true] {
            let (loft, step) = fixture(clockwise);
            let proof = read_analytic_solid(&step, TOL).unwrap();
            let (imported_proof, imported) = read_import_solid(&step, TOL).unwrap();
            assert_eq!(proof, imported_proof);
            assert_eq!(proof.face_count, 6);
            assert_eq!(proof.edges_mm.len(), 8);
            assert_eq!(proof.splines.len(), 4);
            assert!(proof.splines.iter().all(|rail| rail.poles_mm().len() == 8));
            assert_eq!(proof.solid_count, 1);
            assert!((proof.volume_mm3 - loft.volume_mm3()).abs() < 1e-9);
            let ImportProfile::SmoothLoft { geometry } = imported.profile else {
                panic!("smooth loft profile");
            };
            assert!((geometry.volume_mm3() - loft.volume_mm3()).abs() < 1e-9);
            assert_eq!(geometry.heights().len(), 4);
            assert!((geometry.height() - 30.).abs() < TOL);
            assert!(!geometry.clockwise());
            // The bounds enclose the widest intermediate section.
            let normal = cross(imported.u, imported.v);
            let world = |[x, y, z]: Point| -> Point {
                std::array::from_fn(|axis| {
                    imported.origin_mm[axis]
                        + imported.u[axis] * x
                        + imported.v[axis] * y
                        + normal[axis] * z
                })
            };
            for (segment, t) in [(0, 0.5), (1, 0.25), (2, 0.9)] {
                for corner in 0..4 {
                    let p = world(geometry.point(corner, 0., segment, t));
                    assert!((0..3).all(|axis| {
                        p[axis] >= proof.bounds_mm[0][axis] - TOL
                            && p[axis] <= proof.bounds_mm[1][axis] + TOL
                    }));
                }
            }
        }
    }

    #[test]
    fn rejects_moved_wall_poles_and_rails_off_their_heights() {
        let (_, step) = fixture(false);
        // Move one interior rail pole: the wall rows no longer match it.
        let line = step
            .lines()
            .find(|line| line.contains("=B_SPLINE_CURVE_WITH_KNOTS('',3,"))
            .unwrap();
        let pole = line
            .split("(#")
            .nth(1)
            .and_then(|rest| rest.split(',').nth(1))
            .unwrap()
            .trim_start_matches('#')
            .to_owned();
        let point = step
            .lines()
            .find(|l| l.starts_with(&format!("#{pole}=CARTESIAN_POINT(")))
            .unwrap();
        let moved = step.replacen(
            point,
            &format!("#{pole}=CARTESIAN_POINT('',(100.,100.,100.));"),
            1,
        );
        assert!(read_import_solid(&moved, TOL).is_err());
        // Reverse a wall surface's knots out of order.
        let wall = step
            .lines()
            .find(|line| line.contains("=B_SPLINE_SURFACE_WITH_KNOTS('',1,3,"))
            .unwrap();
        let bent = step.replacen(wall, &wall.replacen("(4,2,2,4)", "(4,1,3,4)", 1), 1);
        assert!(read_import_solid(&bent, TOL).is_err());
    }
}
