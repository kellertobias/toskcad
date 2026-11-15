//! Whole finite normalized-field bounds, retaining original record owners.
use super::*;
use crate::tapered_thread::interval::I;
pub(super) struct View<'a> {
    pub(super) source: &'a crate::orthogonal_bore_seats::Seats,
    pub(super) placement: crate::axial_pipe_thread::authored::StockClipPlacement,
    pub(super) raw_faces: &'a [usize],
    pub(super) raw_vertices: &'a BTreeMap<[u64; 3], u32>,
    pub(super) raw_chains: &'a [usize],
    pub(super) reversed: &'a [bool],
}
impl<'a> View<'a> {
    pub(super) fn derived(geometry: &'a Geometry) -> Self {
        Self {
            source: &geometry.source,
            placement: geometry.stock.placement,
            raw_faces: &geometry.stock.raw_faces,
            raw_vertices: &geometry.stock.raw_vertices,
            raw_chains: &geometry.raw_chains,
            reversed: &[],
        }
    }
    fn reversed(&self, bore: usize) -> bool {
        self.reversed.get(bore).copied().unwrap_or(false)
    }
}
type V = [I; 3];
#[derive(Clone, Copy)]
struct Axis {
    origin: V,
    n: V,
    u: V,
    v: V,
}
#[derive(Clone, Copy)]
struct Ring {
    center: V,
    a: V,
    b: V,
}
#[derive(Default)]
pub(super) struct Report {
    pub(super) position: f64,
    pub(super) derivative: f64,
    pub(super) sewing: f64,
    pub(super) faces: usize,
    pub(super) edges: usize,
    pub(super) vertices: usize,
}
fn add(a: V, b: V) -> V {
    std::array::from_fn(|i| a[i] + b[i])
}
fn sub_v(a: V, b: V) -> V {
    std::array::from_fn(|i| a[i] - b[i])
}
fn scale(a: V, s: I) -> V {
    a.map(|a| a * s)
}
fn dot_v(a: V, b: V) -> I {
    (0..3).fold(I::point(0.), |s, i| s + a[i] * b[i])
}
fn cross_v(a: V, b: V) -> V {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn norm(a: V) -> f64 {
    a.into_iter()
        .fold(I::point(0.), |s, a| s + I::point(a.magnitude()))
        .hi
}
fn unit_v(a: V) -> std::result::Result<V, Error> {
    let length = dot_v(a, a).sqrt();
    if !length.valid() || length.lo <= 0. {
        return Err(Error::Unsupported);
    }
    let result = a.map(|x| x / length);
    if result.iter().any(|i| !i.valid()) {
        return Err(Error::Unsupported);
    }
    Ok(result)
}
fn literal(
    boundary: &Boundary<'_, '_>,
    id: u32,
    name: &str,
    budget: &mut Budget<'_>,
) -> std::result::Result<Vec<I>, Error> {
    budget.charge(128)?;
    let entity = boundary.shell.simple(id)?;
    if entity.name != name {
        return Err(Error::Invalid);
    }
    let source = entity.source.ok_or(Error::Invalid)?;
    budget.charge(source.len())?;
    let spelling = part21::numeric_spellings(entity).ok_or(Error::Invalid)?;
    let mut values = Vec::new();
    fn collect(values: &[Value<'_>], out: &mut Vec<f64>) -> std::result::Result<(), Error> {
        for value in values {
            match value {
                Value::Number(v) => out.push(*v),
                Value::List(v) => collect(v, out)?,
                Value::Typed(_) => return Err(Error::Unsupported),
                _ => {}
            }
        }
        Ok(())
    }
    collect(&entity.args, &mut values)?;
    if values.len() != spelling.len() || values.len() > 4 {
        return Err(Error::Unsupported);
    }
    let values = values
        .into_iter()
        .zip(spelling)
        .map(|(v, s)| scalar(v, s))
        .collect::<Vec<_>>();
    if values.iter().any(|i| !i.valid()) {
        return Err(Error::Unsupported);
    }
    Ok(values)
}
fn point(
    boundary: &Boundary<'_, '_>,
    id: u32,
    budget: &mut Budget<'_>,
) -> std::result::Result<V, Error> {
    let values = literal(boundary, id, "CARTESIAN_POINT", budget)?;
    values.try_into().map_err(|_| Error::Invalid)
}
fn direction(
    boundary: &Boundary<'_, '_>,
    id: u32,
    budget: &mut Budget<'_>,
) -> std::result::Result<V, Error> {
    let values = literal(boundary, id, "DIRECTION", budget)?;
    unit_v(values.try_into().map_err(|_| Error::Invalid)?)
}
fn axis(
    boundary: &Boundary<'_, '_>,
    id: u32,
    budget: &mut Budget<'_>,
) -> std::result::Result<Axis, Error> {
    budget.charge(2048)?;
    let [_, o, n, r] = boundary.shell.args(id, "AXIS2_PLACEMENT_3D")? else {
        return Err(Error::Invalid);
    };
    let origin = point(boundary, reference(o)?, budget)?;
    let n = direction(boundary, reference(n)?, budget)?;
    let r = direction(boundary, reference(r)?, budget)?;
    let v = unit_v(cross_v(n, r))?;
    let u = cross_v(v, n);
    Ok(Axis { origin, n, u, v })
}
fn edge_geometry(boundary: &Boundary<'_, '_>, id: u32) -> std::result::Result<u32, Error> {
    let [_, _, _, geometry, _] = boundary.shell.args(id, "EDGE_CURVE")? else {
        return Err(Error::Invalid);
    };
    Ok(reference(geometry)?)
}
fn raw_ring(
    boundary: &Boundary<'_, '_>,
    edge: u32,
    budget: &mut Budget<'_>,
) -> std::result::Result<(Ring, Axis), Error> {
    budget.charge(2048)?;
    let geometry = edge_geometry(boundary, edge)?;
    let [_, placement, _] = boundary.shell.args(geometry, "CIRCLE")? else {
        return Err(Error::Invalid);
    };
    let mut axis = axis(boundary, reference(placement)?, budget)?;
    let values = literal(boundary, geometry, "CIRCLE", budget)?;
    let [radius] = values.as_slice() else {
        return Err(Error::Invalid);
    };
    let radius = *radius;
    let vertex = raw_point(boundary, boundary.shell.edges[&edge].vertices[0], budget)?;
    let q = sub_v(vertex, axis.origin);
    let ray = unit_v(add(
        scale(axis.u, dot_v(q, axis.u)),
        scale(axis.v, dot_v(q, axis.v)),
    ))?;
    let [_, _, _, _, forward] = boundary.shell.args(edge, "EDGE_CURVE")? else {
        return Err(Error::Invalid);
    };
    if !sense(forward)? {
        axis.n = scale(axis.n, I::point(-1.));
    }
    Ok((
        Ring {
            center: axis.origin,
            a: scale(ray, radius),
            b: scale(cross_v(axis.n, ray), radius),
        },
        axis,
    ))
}
fn matrix_vector(matrix: [[f64; 3]; 3], v: V) -> V {
    std::array::from_fn(|i| (0..3).fold(I::point(0.), |sum, j| sum + I::point(matrix[i][j]) * v[j]))
}
fn bounds(a: Ring, b: Ring) -> (f64, f64) {
    let center = norm(sub_v(a.center, b.center));
    let radial = I::point(norm(sub_v(a.a, b.a))) + I::point(norm(sub_v(a.b, b.b)));
    let position = (I::point(center) + radial).hi;
    let pi = I::new(
        std::f64::consts::PI.next_down(),
        std::f64::consts::PI.next_up(),
    );
    (position, (I::point(2.) * pi * radial).hi)
}
fn update(report: &mut Report, position: f64, derivative: f64) -> std::result::Result<(), Error> {
    if !position.is_finite() || !derivative.is_finite() {
        return Err(Error::Unsupported);
    }
    report.position = report.position.max(position);
    report.derivative = report.derivative.max(derivative);
    Ok(())
}
fn layers(
    source: &crate::orthogonal_bore_seats::Seats,
    index: usize,
) -> std::result::Result<Vec<(V, I)>, Error> {
    let bore = &source.bores()[index];
    let [a, b] = bore.centers_mm;
    let axis = (0..3).find(|i| a[*i] != b[*i]).ok_or(Error::Unsupported)?;
    let sign = (b[axis] - a[axis]).signum();
    let mut rows = vec![(a.map(I::point), I::point(bore.radius_mm))];
    if let Some(seat) = source.seats().iter().find(|s| s.bore == index && !s.end) {
        rows[0].1 = I::point(seat.radius_mm);
        let mut p = a.map(I::point);
        p[axis] = p[axis] + I::point(sign) * I::point(seat.depth_mm);
        rows.push((p, I::point(bore.radius_mm)));
    }
    if let Some(seat) = source.seats().iter().find(|s| s.bore == index && s.end) {
        let mut p = b.map(I::point);
        p[axis] = p[axis] - I::point(sign) * I::point(seat.depth_mm);
        rows.push((p, I::point(bore.radius_mm)));
        rows.push((b.map(I::point), I::point(seat.radius_mm)));
    } else {
        rows.push((b.map(I::point), I::point(bore.radius_mm)));
    }
    Ok(rows)
}
fn placed(placement: crate::axial_pipe_thread::authored::StockClipPlacement, p: V) -> V {
    add(
        matrix_vector(placement.matrix(), p),
        placement.translation_mm().map(I::point),
    )
}
fn orientation(normal: V, mapped_normal: V) -> std::result::Result<f64, Error> {
    let dot = dot_v(normal, mapped_normal);
    if dot.lo > 0. {
        Ok(1.)
    } else if dot.hi < 0. {
        Ok(-1.)
    } else {
        Err(Error::Unsupported)
    }
}
fn phase(
    geometry: &View<'_>,
    chain: &BoreChain,
    bore: usize,
    boundary: &Boundary<'_, '_>,
    budget: &mut Budget<'_>,
) -> std::result::Result<(V, V, V), Error> {
    let source = geometry.source;
    let center = source.bores()[bore].centers_mm[usize::from(geometry.reversed(bore))];
    let axis_index = (0..3)
        .find(|i| source.bores()[bore].centers_mm[0][*i] != source.bores()[bore].centers_mm[1][*i])
        .ok_or(Error::Unsupported)?;
    let vertex = boundary.shell.edges[&chain.rim_edge_entities()[0]].vertices[0];
    let local = super::local(geometry.placement, boundary.shell.points[&vertex])?;
    let mut radial = std::array::from_fn(|i| I::point(local[i] - center[i]));
    radial[axis_index] = I::point(0.);
    let a = unit_v(radial)?;
    let (_, raw_axis) = raw_ring(boundary, chain.rim_edge_entities()[0], budget)?;
    let raw_n = raw_axis.n;
    let matrix = geometry.placement.matrix();
    let i = (axis_index + 1) % 3;
    let j = (axis_index + 2) % 3;
    let column = |j| std::array::from_fn(|i| I::point(matrix[i][j]));
    let normal = cross_v(column(i), column(j));
    let sign = orientation(raw_n, normal)?;
    let mut k = [I::point(0.); 3];
    k[axis_index] = I::point(sign);
    let b = cross_v(k, a);
    Ok((
        matrix_vector(matrix, a),
        matrix_vector(matrix, b),
        scale(normal, I::point(sign)),
    ))
}
fn raw_wall(
    boundary: &Boundary<'_, '_>,
    face: u32,
    rims: [u32; 2],
    budget: &mut Budget<'_>,
) -> std::result::Result<([Ring; 2], Axis), Error> {
    budget.charge(4096)?;
    let [_, _, surface, _] = boundary.shell.args(face, "ADVANCED_FACE")? else {
        return Err(Error::Invalid);
    };
    let surface = reference(surface)?;
    let entity = boundary.shell.simple(surface)?;
    let (_, placement) = match entity.args.as_slice() {
        [_, placement, ..] => (entity.name, placement),
        _ => return Err(Error::Invalid),
    };
    let axis = axis(boundary, reference(placement)?, budget)?;
    let values = literal(boundary, surface, entity.name, budget)?;
    let (radius, slope) = match (entity.name, values.as_slice()) {
        ("CYLINDRICAL_SURFACE", [radius]) => (*radius, I::point(0.)),
        ("CONICAL_SURFACE", [radius, angle]) => {
            let cosine = angle.cos();
            if cosine.lo <= 0. {
                return Err(Error::Unsupported);
            }
            (*radius, angle.sin() / cosine)
        }
        _ => return Err(Error::Unsupported),
    };
    let (first, _) = raw_ring(boundary, rims[0], budget)?;
    let vertex = raw_point(boundary, boundary.shell.edges[&rims[0]].vertices[0], budget)?;
    let z = dot_v(sub_v(first.center, axis.origin), axis.n);
    let center = add(axis.origin, scale(axis.n, z));
    let q = sub_v(vertex, center);
    let ray = unit_v(sub_v(q, scale(axis.n, dot_v(q, axis.n))))?;
    let side = cross_v(axis.n, ray);
    let mut rows = [Ring {
        center,
        a: ray,
        b: side,
    }; 2];
    for (row, rim) in rows.iter_mut().zip(rims) {
        let (circle, _) = raw_ring(boundary, rim, budget)?;
        let z = dot_v(sub_v(circle.center, axis.origin), axis.n);
        let r = radius + slope * z;
        *row = Ring {
            center: add(axis.origin, scale(axis.n, z)),
            a: scale(ray, r),
            b: scale(side, r),
        };
    }
    Ok((rows, axis))
}
fn plane(
    geometry: &View<'_>,
    boundary: &Boundary<'_, '_>,
    source_face: usize,
    raw_face: usize,
    vertices: &BTreeMap<u32, V>,
    budget: &mut Budget<'_>,
) -> std::result::Result<(f64, f64), Error> {
    budget.charge(4096)?;
    let source = &geometry.source.boundary().faces[source_face];
    let id = boundary.face_entities[raw_face];
    let [_, _, surface, forward] = boundary.shell.args(id, "ADVANCED_FACE")? else {
        return Err(Error::Invalid);
    };
    let [_, placement] = boundary.shell.args(reference(surface)?, "PLANE")? else {
        return Err(Error::Invalid);
    };
    let raw = axis(boundary, reference(placement)?, budget)?;
    let major = source
        .normal
        .iter()
        .position(|n| n.abs() == 1.)
        .ok_or(Error::Unsupported)?;
    let matrix = geometry.placement.matrix();
    let column = |j| std::array::from_fn(|i| I::point(matrix[i][j]));
    let cofactor = scale(
        cross_v(column((major + 1) % 3), column((major + 2) % 3)),
        I::point(source.normal[major]),
    );
    let physical_normal = scale(raw.n, I::point(if sense(forward)? { 1. } else { -1. }));
    if dot_v(cofactor, physical_normal).lo <= 0. {
        return Err(Error::Unsupported);
    }
    let mut error = 0f64;
    for p in source.rings.iter().flatten() {
        let vertex = geometry.raw_vertices[&crate::orthogonal_boundary::point_key(*p)];
        let actual = raw_point(boundary, vertex, budget)?;
        let projected = sub_v(
            actual,
            scale(raw.n, dot_v(sub_v(actual, raw.origin), raw.n)),
        );
        let expected = *vertices.get(&vertex).ok_or(Error::Unsupported)?;
        error = error.max(norm(sub_v(expected, projected)));
    }
    // Use the SAME actual owner-ID triangulation for every source view: the
    // original outward wire, cyclically anchored at its smallest vertex ID.
    // Its signed fan covers the entire outer polygon; holes only remove points.
    // Barycentric interpolation bounds position by max vertex error and each
    // normalized derivative by twice it. This fixed raw field (rather than a
    // source-dependent projection) makes original/raw/healed C1 composition valid.
    Ok((error, (I::point(2.) * I::point(error)).hi))
}
fn line(
    boundary: &Boundary<'_, '_>,
    edge: u32,
    vertices: &BTreeMap<u32, V>,
    budget: &mut Budget<'_>,
) -> std::result::Result<(f64, f64, f64), Error> {
    budget.charge(4096)?;
    let id = edge_geometry(boundary, edge)?;
    let [_, origin, vector] = boundary.shell.args(id, "LINE")? else {
        return Err(Error::Invalid);
    };
    let origin = point(boundary, reference(origin)?, budget)?;
    let vector = reference(vector)?;
    let [_, dir, _] = boundary.shell.args(vector, "VECTOR")? else {
        return Err(Error::Invalid);
    };
    let dir = direction(boundary, reference(dir)?, budget)?;
    let values = literal(boundary, vector, "VECTOR", budget)?;
    let [magnitude] = values.as_slice() else {
        return Err(Error::Invalid);
    };
    if magnitude.lo <= 0. {
        return Err(Error::Unsupported);
    }
    let mut error = 0f64;
    let mut sewing = 0f64;
    for vertex in boundary.shell.edges[&edge].vertices {
        let raw = raw_point(boundary, vertex, budget)?;
        let projected = add(origin, scale(dir, dot_v(sub_v(raw, origin), dir)));
        let expected = *vertices.get(&vertex).ok_or(Error::Unsupported)?;
        error = error.max(norm(sub_v(expected, projected)));
        sewing = sewing.max(norm(sub_v(raw, projected)));
    }
    Ok((error, (I::point(2.) * I::point(error)).hi, sewing))
}
pub(super) fn correspondence(
    geometry: &View<'_>,
    chains: &Chains<'_, '_, '_>,
    budget: &mut Budget<'_>,
) -> std::result::Result<Report, Error> {
    let boundary = chains.boundary();
    let mut report = Report::default();
    let mut vertices = BTreeMap::new();
    let mut circles = BTreeSet::new();
    let mut faces = BTreeSet::new();
    for (key, raw) in geometry.raw_vertices {
        let p = key.map(f64::from_bits);
        let p = world(geometry.placement, p)?;
        if vertices.insert(*raw, p).is_some() {
            return Err(Error::Unsupported);
        }
    }
    for (source_face, &raw_face) in geometry.raw_faces.iter().enumerate() {
        let (position, derivative) =
            plane(geometry, boundary, source_face, raw_face, &vertices, budget)?;
        update(&mut report, position, derivative)?;
        if !faces.insert(raw_face) {
            return Err(Error::Unsupported);
        }
    }
    for (bore, &raw_chain) in geometry.raw_chains.iter().enumerate() {
        let chain = &chains.bores()[raw_chain];
        let mut rows = layers(geometry.source, bore)?;
        if geometry.reversed(bore) {
            rows.reverse();
        }
        if rows.len() != chain.rim_edge_entities().len() {
            return Err(Error::Unsupported);
        }
        let (a, b, normal) = phase(geometry, chain, bore, boundary, budget)?;
        let mut source = Vec::new();
        for ((p, r), &rim) in rows.iter().zip(chain.rim_edge_entities()) {
            budget.charge(4096)?;
            let (raw, axis) = raw_ring(boundary, rim, budget)?;
            let sign = orientation(axis.n, normal)?;
            let ring = Ring {
                center: placed(geometry.placement, *p),
                a: scale(a, *r),
                b: scale(b, *r * I::point(sign)),
            };
            let (position, derivative) = bounds(ring, raw);
            update(&mut report, position, derivative)?;
            if !circles.insert(rim) {
                return Err(Error::Unsupported);
            }
            let vertex = boundary.shell.edges[&rim].vertices[0];
            let point = add(ring.center, ring.a);
            if vertices.insert(vertex, point).is_some() {
                return Err(Error::Unsupported);
            }
            source.push(ring);
        }
        for (wall, &raw_id) in chain.wall_face_entities().iter().enumerate() {
            budget.charge(4096)?;
            let rims = [
                chain.rim_edge_entities()[wall],
                chain.rim_edge_entities()[wall + 1],
            ];
            let (raw, axis) = raw_wall(boundary, raw_id, rims, budget)?;
            let sign = orientation(axis.n, normal)?;
            let raw_travel = dot_v(sub_v(raw[1].center, raw[0].center), axis.n);
            let source_travel = dot_v(sub_v(source[wall + 1].center, source[wall].center), axis.n);
            if !((raw_travel.lo > 0. && source_travel.lo > 0.)
                || (raw_travel.hi < 0. && source_travel.hi < 0.))
            {
                return Err(Error::Unsupported);
            }
            let mut error = 0f64;
            let mut radial = 0f64;
            for (i, raw) in raw.into_iter().enumerate() {
                let ring = Ring {
                    center: source[wall + i].center,
                    a: scale(a, rows[wall + i].1),
                    b: scale(b, rows[wall + i].1 * I::point(sign)),
                };
                let (position, derivative) = bounds(ring, raw);
                error = error.max(position);
                radial = radial.max(derivative);
            }
            update(
                &mut report,
                error,
                radial.max((I::point(2.) * I::point(error)).hi),
            )?;
            let index = boundary
                .face_entities
                .iter()
                .position(|id| *id == raw_id)
                .ok_or(Error::Invalid)?;
            if !faces.insert(index) {
                return Err(Error::Unsupported);
            }
        }
    }
    if faces.len() != boundary.face_entities.len() || vertices.len() != boundary.shell.points.len()
    {
        return Err(Error::Unsupported);
    }
    for (&vertex, &expected) in &vertices {
        budget.charge(2048)?;
        let actual = raw_point(boundary, vertex, budget)?;
        let error = norm(sub_v(expected, actual));
        report.sewing = report.sewing.max(error);
        update(&mut report, error, 0.)?;
    }
    for (&id, edge) in &boundary.shell.edges {
        match edge.curve {
            Curve::Line => {
                let (p, d, s) = line(boundary, id, &vertices, budget)?;
                update(&mut report, p, d)?;
                report.sewing = report.sewing.max(s);
            }
            Curve::Circle(_) => {
                if !circles.contains(&id) {
                    return Err(Error::Unsupported);
                }
            }
            _ => return Err(Error::Unsupported),
        }
    }
    // Every exact derived boundary is shared. Each of its raw incident
    // fields is within position of that common boundary (with the retained
    // orientation/rephasing above), so their whole finite sewing gap is <=2E.
    report.sewing = report
        .sewing
        .max((I::point(2.) * I::point(report.position)).hi);
    report.faces = faces.len();
    report.edges = boundary.shell.edges.len();
    report.vertices = vertices.len();
    Ok(report)
}

#[test]
fn all_original_finite_supports_have_whole_position_and_derivative_bounds() {
    use crate::axial_pipe_thread::authored::StockClipPlacement;
    for placement in [
        StockClipPlacement::new([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]], [0.; 3]).unwrap(),
        StockClipPlacement::new(
            [[0.8, -0.6, 0.], [0.6, 0.8, 0.], [0., 0., 1.]],
            [11., -7., 5.],
        )
        .unwrap(),
    ] {
        let text = super::super::tests::text_placed(placement);
        let document =
            CompoundDocument::parse(&text, owned_insert::DocumentOptions::default(), None).unwrap();
        let boundary = decode(
            &document,
            document.represented_members()[0].item_entity(),
            Options::default(),
            &|| false,
        )
        .unwrap();
        let chains = boundary
            .through_chains(Options::default(), &|| false)
            .unwrap();
        let mut budget = Budget {
            work: 0,
            maximum: 16_000_000,
            interrupt: &|| false,
        };
        let geometry = super::geometry(&chains, 1e-7 / 16., &mut budget).unwrap();
        let view = View::derived(&geometry);
        let report = correspondence(&view, &chains, &mut budget).unwrap();
        eprintln!(
            "actual finite bound C0={} D1={} sewing={} work={}",
            report.position, report.derivative, report.sewing, budget.work
        );
        assert_eq!(report.faces, 12);
        assert_eq!(report.edges, boundary.edge_count());
        assert_eq!(report.vertices, boundary.vertex_count());
        assert!(report.position < 1e-7 / 16.);
        assert!(report.derivative < 1e-7 / 16.);
        assert!(report.sewing < 1e-7 / 8.);
    }
}

pub(super) struct Owners {
    pub(super) faces: Vec<usize>,
    pub(super) vertices: BTreeMap<[u64; 3], u32>,
    pub(super) chains: Vec<usize>,
    pub(super) reversed: Vec<bool>,
}
pub(super) fn original_owners(
    source: &crate::orthogonal_bore_seats::Seats,
    placement: crate::axial_pipe_thread::authored::StockClipPlacement,
    expected_faces: &[u32],
    derived: &Geometry,
    chains: &Chains<'_, '_, '_>,
    position_limit: f64,
    budget: &mut Budget<'_>,
) -> std::result::Result<Owners, Error> {
    let boundary = chains.boundary();
    let stock = &source.boundary();
    if stock.faces.len() != derived.stock.raw_faces.len()
        || expected_faces.len() != stock.faces.len() + source.bores().len() + source.seats().len()
        || expected_faces.len() != boundary.face_entities.len()
    {
        return Err(Error::Unsupported);
    }
    let faces = expected_faces[..stock.faces.len()]
        .iter()
        .map(|id| {
            boundary
                .face_entities
                .iter()
                .position(|raw| raw == id)
                .ok_or(Error::Unsupported)
        })
        .collect::<std::result::Result<Vec<_>, _>>()?;
    let mut vertices = BTreeMap::new();
    let mut used = BTreeSet::new();
    for p in stock.faces.iter().flat_map(|f| &f.rings).flatten() {
        let key = crate::orthogonal_boundary::point_key(*p);
        if vertices.contains_key(&key) {
            continue;
        }
        let expected = world(placement, *p)?;
        let mut found = None;
        for &vertex in derived.stock.raw_vertices.values() {
            budget.charge(256)?;
            let actual = raw_point(boundary, vertex, budget)?;
            let error = norm(sub_v(expected, actual));
            if error <= position_limit && found.replace(vertex).is_some() {
                return Err(Error::Unsupported);
            }
        }
        let original = found.ok_or(Error::Unsupported)?;
        if !used.insert(original) {
            return Err(Error::Unsupported);
        }
        vertices.insert(key, original);
    }
    if used.len() != derived.stock.raw_vertices.len() {
        return Err(Error::Unsupported);
    }
    for (source_face, &raw_face) in stock.faces.iter().zip(&faces) {
        let mut successors = BTreeMap::new();
        for edge in boundary.shell.edges.values() {
            budget.charge(64)?;
            if !matches!(edge.curve, Curve::Line) {
                continue;
            }
            if let Some(u) = edge.uses.iter().find(|u| u.face == raw_face) {
                let [a, b] = if u.forward {
                    edge.vertices
                } else {
                    [edge.vertices[1], edge.vertices[0]]
                };
                if successors.insert(a, b).is_some() {
                    return Err(Error::Unsupported);
                }
            }
        }
        let count = source_face.rings.iter().map(Vec::len).sum::<usize>();
        if successors.len() != count {
            return Err(Error::Unsupported);
        }
        for ring in &source_face.rings {
            for (a, b) in ring.iter().zip(ring.iter().cycle().skip(1)) {
                if successors.get(&vertices[&crate::orthogonal_boundary::point_key(*a)])
                    != Some(&vertices[&crate::orthogonal_boundary::point_key(*b)])
                {
                    return Err(Error::Unsupported);
                }
            }
        }
    }
    let mut matched = Vec::new();
    let mut reversed = Vec::new();
    let mut used = BTreeSet::new();
    for (bore_index, bore) in source.bores().iter().enumerate() {
        let caps = bore.faces.map(|i| expected_faces[i]);
        let mut found = None;
        for (index, chain) in chains.bores().iter().enumerate() {
            budget.charge(256)?;
            let raw = chain.cap_face_entities();
            let reverse = raw == [caps[1], caps[0]];
            let cylinder = chain.wall_face_entities()[chain.cylinder_index()];
            if (raw == caps || reverse)
                && cylinder == expected_faces[stock.faces.len() + bore_index]
                && found.replace((index, reverse)).is_some()
            {
                return Err(Error::Unsupported);
            }
        }
        let (index, reverse) = found.ok_or(Error::Unsupported)?;
        let chain = &chains.bores()[index];
        for (ordinal, seat) in source
            .seats()
            .iter()
            .enumerate()
            .filter(|(_, s)| s.bore == bore_index)
        {
            let end = seat.end ^ reverse;
            let wall = if end {
                chain.wall_face_entities().len() - 1
            } else {
                0
            };
            if chain.wall_face_entities()[wall]
                != expected_faces[stock.faces.len() + source.bores().len() + ordinal]
            {
                return Err(Error::Unsupported);
            }
        }
        if !used.insert(index) {
            return Err(Error::Unsupported);
        }
        matched.push(index);
        reversed.push(reverse);
    }
    if used.len() != chains.bores().len() {
        return Err(Error::Unsupported);
    }
    Ok(Owners {
        faces,
        vertices,
        chains: matched,
        reversed,
    })
}
#[test]
fn original_literal_matrix_is_checked_separately_from_healed_material() {
    use crate::{axial_pipe_thread::authored::StockClipPlacement, step::orthogonal_seats};
    let original = super::super::tests::source();
    for (scale, admitted) in [
        (1., true),
        (1. + 2f64.powi(-42), true),
        (1. + 2f64.powi(-24), false),
    ] {
        let placement = StockClipPlacement::new(
            [[scale, 0., 0.], [0., 1., 0.], [0., 0., 1.]],
            [11., -7., 5.],
        )
        .unwrap();
        let candidate = orthogonal_seats::candidate(
            &original,
            placement,
            orthogonal_seats::Options::default(),
            &|| false,
        )
        .unwrap();
        let document = CompoundDocument::parse(
            candidate.as_str(),
            owned_insert::DocumentOptions::default(),
            None,
        )
        .unwrap();
        let boundary = decode(
            &document,
            document.represented_members()[0].item_entity(),
            Options::default(),
            &|| false,
        )
        .unwrap();
        let chains = boundary
            .through_chains(Options::default(), &|| false)
            .unwrap();
        let mut budget = Budget {
            work: 0,
            maximum: 16_000_000,
            interrupt: &|| false,
        };
        let derived = super::geometry(&chains, 1e-7 / 16., &mut budget).unwrap();
        let data = original_owners(
            &original,
            placement,
            candidate.face_entities(),
            &derived,
            &chains,
            1e-7 / 16.,
            &mut budget,
        )
        .unwrap();
        let view = View {
            source: &original,
            placement,
            raw_faces: &data.faces,
            raw_vertices: &data.vertices,
            raw_chains: &data.chains,
            reversed: &data.reversed,
        };
        let report = correspondence(&view, &chains, &mut budget).unwrap();
        eprintln!(
            "actual original scale={scale} C0={} D1={} work={}",
            report.position, report.derivative, budget.work
        );
        assert_eq!(report.position < 1e-7 && report.derivative < 1e-7, admitted);
        let healed = correspondence(&View::derived(&derived), &chains, &mut budget).unwrap();
        assert!(healed.position < 1e-7 && healed.derivative < 1e-7);
    }
}
