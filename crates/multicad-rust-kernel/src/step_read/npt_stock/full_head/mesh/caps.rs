//! Fresh whole ACTUAL curve/span correspondence to SAME emitted cap chords.
use super::*;
fn evaluate(mut p: [[I; 3]; 12], t: f64) -> [I; 3] {
    let t = I::point(t);
    let one = I::point(1.) - t;
    for n in (1..12).rev() {
        for i in 0..n {
            p[i] = std::array::from_fn(|k| one * p[i][k] + t * p[i + 1][k]);
        }
    }
    p[0]
}
pub(super) fn chord(
    source: &OwnedMaterial,
    sewing: AxialSewingTransfer,
    edge: usize,
    span: usize,
    t: [f64; 2],
    literal: [[f64; 3]; 2],
    w: &mut Work<'_>,
) -> Result<f64> {
    w.step(4096)?;
    if t.iter().any(|v| !v.is_finite() || !(0.0..=1.0).contains(v)) || t[0] == t[1] {
        return Err(Error::Unsupported);
    }
    let curve = source
        .curves()
        .get(
            source
                .edges()
                .get(edge)
                .ok_or(Error::Invalid)?
                .curve_index(),
        )
        .ok_or(Error::Invalid)?;
    if span >= curve.spans() {
        return Err(Error::Invalid);
    }
    let indices = curve
        .point_indices()
        .get(11 * span..11 * span + 12)
        .ok_or(Error::Invalid)?;
    let mut poles = [[I::point(0.); 3]; 12];
    for (p, &id) in poles.iter_mut().zip(indices) {
        let raw = source
            .point(id)
            .ok_or(Error::Invalid)?
            .coordinates()
            .map(|s| s.bounds());
        let actual = source.map_box_mm(raw)?;
        *p = sewing.map_raw_box_mm(actual)?.map(|b| I::new(b[0], b[1]));
    }
    let mut second = [0f64; 3];
    for a in poles.windows(3) {
        for k in 0..3 {
            let d = (a[2][k] - I::point(2.) * a[1][k] + a[0][k]) * I::point(110.);
            if !d.valid() {
                return Err(Error::Unsupported);
            }
            second[k] = second[k].max(d.magnitude());
        }
    }
    let hessian = second
        .into_iter()
        .fold(I::point(0.), |s, x| s + I::point(x));
    let mut endpoints = 0f64;
    for i in 0..2 {
        let exact = evaluate(poles, t[i]);
        let error = (0..3).fold(I::point(0.), |s, k| {
            s + I::point((exact[k] - I::point(literal[i][k])).magnitude())
        });
        if !error.valid() {
            return Err(Error::Unsupported);
        }
        endpoints = endpoints.max(error.hi);
    }
    // Coordinatewise linear-interpolation remainder over the WHOLE span:
    // sup|C-linear(C(a),C(b))| <= (b-a)^2 sup|C''|/8.
    // Convex endpoint perturbation then pairs the complete original arc and
    // SAME literal segment with a two-sided L1 (hence Euclidean) C0 bound.
    let delta = I::point(t[1]) - I::point(t[0]);
    let error = I::point(endpoints) + hessian * delta.square() / I::point(8.);
    if !error.valid() {
        return Err(Error::Unsupported);
    }
    Ok(error.hi)
}
pub(super) fn annulus(
    mesh: &mut DisplayMesh,
    outer: &[u32],
    inner: &[u32],
    face: u32,
    source_normal: [f64; 3],
    angle: f64,
    w: &mut Work<'_>,
) -> Result<f64> {
    let n = outer.len().checked_add(inner.len()).ok_or(Error::Limit)?;
    w.step(n.checked_mul(512).ok_or(Error::Limit)?)?;
    if outer.len() != 6 || inner.len() < 3 {
        return Err(Error::Unsupported);
    }
    let point = |id: u32| -> Result<[f64; 3]> {
        mesh.vertices
            .get(3 * id as usize..3 * id as usize + 3)
            .and_then(|p| p.try_into().ok())
            .ok_or(Error::Invalid)
    };
    let mut winding = I::point(0.);
    for i in 0..inner.len() {
        let a = point(inner[i])?;
        let b = point(inner[(i + 1) % inner.len()])?;
        winding = winding
            + crate::polar_loop11::positive_chord_angle(
                [I::point(a[0]), I::point(a[1])],
                [I::point(b[0]), I::point(b[1])],
                true,
            )
            .ok_or(Error::Unsupported)?;
        for j in 0..outer.len() {
            let p = point(outer[j])?;
            let q = point(outer[(j + 1) % outer.len()])?;
            let side = (I::point(q[0]) - I::point(p[0])) * (I::point(a[1]) - I::point(p[1]))
                - (I::point(q[1]) - I::point(p[1])) * (I::point(a[0]) - I::point(p[0]));
            if !side.valid() || side.lo <= 0. {
                return Err(Error::Unsupported);
            }
        }
    }
    if !winding.valid() || winding.lo <= 4. || winding.hi >= 8. {
        return Err(Error::Unsupported);
    }
    // Pinned Spade2.15.1 bound used by the existing stock-join producer: only
    // actual original vertices, <=3n edges/2n faces. Includes input/copy/flood
    // vectors and the returned DCEL/triangles; no Steiner geometry is added.
    w.storage(
        n.checked_add(8)
            .and_then(|n| n.checked_mul(8192))
            .and_then(|n| n.checked_add(1024 * 1024))
            .ok_or(Error::Limit)?,
    )?;
    let mut vertices = Vec::new();
    vertices.try_reserve_exact(n).map_err(|_| Error::Limit)?;
    vertices.extend(outer);
    vertices.extend(inner);
    let mut points = Vec::new();
    points.try_reserve_exact(n).map_err(|_| Error::Limit)?;
    for &id in &vertices {
        let p = point(id)?;
        points.push(spade::Point2::new(p[0], p[1]));
    }
    let mut edges = Vec::new();
    edges.try_reserve_exact(n).map_err(|_| Error::Limit)?;
    edges.extend((0..outer.len()).map(|i| [i, (i + 1) % outer.len()]));
    edges.extend((0..inner.len()).map(|i| [outer.len() + (i + 1) % inner.len(), outer.len() + i]));
    let region =
        crate::display::triangulate_region(points, &edges, 1).map_err(|_| Error::Unsupported)?;
    let mut bound = 0f64;
    for triangle in region.triangles {
        let mut ids = triangle.map(|i| vertices[i]);
        ids.swap(1, 2); // actual annulus faces down
        let p = ids.map(|id| std::array::from_fn(|k| mesh.vertices[3 * id as usize + k]));
        let error = normal::planar(source_normal, p, w)?;
        if error > angle {
            return Err(Error::Unsupported);
        }
        bound = bound.max(error);
        mesh.indices.extend(ids);
        mesh.face_indices.push(face);
    }
    Ok(bound)
}
