use super::*;
pub(super) struct Geometry {
    pub polygon: [[f64; 2]; 6],
    pub stations: [f64; 2],
    pub planar_sewing: f64,
    pub join_sewing: f64,
    pub minimum_wall: f64,
}
fn ring(head: &HeadInput, face: &HeadFace) -> Vec<usize> {
    face.outer
        .iter()
        .map(|&(edge, d)| head.edges[edge][usize::from(!d)])
        .collect()
}
fn cross(a: [I; 3], b: [I; 3]) -> [I; 3] {
    std::array::from_fn(|k| a[(k + 1) % 3] * b[(k + 2) % 3] - a[(k + 2) % 3] * b[(k + 1) % 3])
}
fn dot(a: [I; 3], b: [I; 3]) -> I {
    (0..3).fold(I::point(0.), |s, k| s + a[k] * b[k])
}
fn norm(a: [I; 3]) -> I {
    dot(a, a).sqrt()
}
fn side(a: [f64; 2], b: [f64; 2], point: [I; 2]) -> (I, I) {
    let d = std::array::from_fn::<_, 2, _>(|k| I::point(b[k]) - I::point(a[k]));
    let p = std::array::from_fn::<_, 2, _>(|k| point[k] - I::point(a[k]));
    (
        d[0] * p[1] - d[1] * p[0],
        (d[0].square() + d[1].square()).sqrt(),
    )
}
pub(super) fn qualify(
    parsed: &Parsed<'_, '_>,
    thread: &Occupied<'_, '_, '_>,
    head: &HeadInput,
    tolerance: f64,
    minimum: f64,
    work: &mut Work<'_>,
) -> Result<Geometry> {
    let annulus = head
        .faces
        .iter()
        .find(|f| !f.inner.is_empty())
        .ok_or(Error::Invalid)?;
    let low = annulus.plane.origin[2];
    require(
        annulus.plane.origin == [0., 0., low]
            && annulus.plane.normal == [0., 0., 1.]
            && annulus.plane.reference == [1., 0., 0.]
            && !annulus.forward
            && annulus.outer.len() == 6,
    )?;
    let candidates = head
        .faces
        .iter()
        .filter(|f| {
            f.inner.is_empty()
                && f.outer.len() == 6
                && f.forward
                && f.plane.normal == [0., 0., 1.]
                && f.plane.origin[2] > low
        })
        .collect::<Vec<_>>();
    let [top] = candidates.as_slice() else {
        return Err(Error::Unsupported);
    };
    let high = top.plane.origin[2];
    require(top.plane.origin == [0., 0., high] && top.plane.reference == [1., 0., 0.])?;
    let lower = ring(head, annulus);
    let upper = ring(head, top);
    require(
        lower.iter().all(|&i| head.vertices[i][2] == low)
            && upper.iter().all(|&i| head.vertices[i][2] == high),
    )?;
    let polygon = std::array::from_fn(|i| {
        let p = head.vertices[upper[i]];
        [p[0], p[1]]
    });
    let mut lower_for_upper = [0; 6];
    for i in 0..6 {
        work.charge(12)?;
        let matched = lower
            .iter()
            .copied()
            .filter(|&k| head.vertices[k][..2] == polygon[i])
            .collect::<Vec<_>>();
        let [id] = matched.as_slice() else {
            return Err(Error::Unsupported);
        };
        lower_for_upper[i] = *id;
    }
    // Every nonadjacent vertex is strictly inside EVERY oriented edge. Local
    // turns alone would also accept self-intersecting stars and are insufficient.
    for i in 0..6 {
        for j in 0..6 {
            work.charge(16)?;
            if j != i && j != (i + 1) % 6 {
                let (s, length) = side(polygon[i], polygon[(i + 1) % 6], polygon[j].map(I::point));
                require(s.valid() && s.lo > 0. && length.valid() && length.lo > 0.)?;
            }
        }
    }
    let side_faces = head
        .faces
        .iter()
        .filter(|f| f.inner.is_empty() && f.outer.len() == 4)
        .collect::<Vec<_>>();
    require(side_faces.len() == 6)?;
    for i in 0..6 {
        work.charge(32)?;
        let j = (i + 1) % 6;
        let expected = [lower_for_upper[i], lower_for_upper[j], upper[j], upper[i]];
        let matching = side_faces
            .iter()
            .filter(|f| {
                let actual = ring(head, f);
                (0..4).any(|shift| (0..4).all(|k| actual[k] == expected[(k + shift) % 4]))
            })
            .count();
        require(matching == 1)?;
    }
    let mut planar = 0_f64;
    for face in &head.faces {
        work.charge(32)?;
        let plane = face.plane;
        require(
            plane
                .origin
                .iter()
                .chain(&plane.normal)
                .chain(&plane.reference)
                .all(|x| x.is_finite() && x.abs() <= 1e7),
        )?;
        let n = plane.normal.map(I::point);
        let length = norm(n);
        let transverse = norm(cross(n, plane.reference.map(I::point)));
        require(length.valid() && length.lo > 0. && transverse.valid() && transverse.lo > 0.)?;
        let vertices = ring(head, face);
        let mut area = [I::point(0.); 3];
        for k in 0..vertices.len() {
            work.charge(48)?;
            let point = head.vertices[vertices[k]].map(I::point);
            let next = head.vertices[vertices[(k + 1) % vertices.len()]].map(I::point);
            let difference = std::array::from_fn(|j| point[j] - I::point(plane.origin[j]));
            let residual = dot(n, difference) / length;
            require(residual.valid())?;
            planar = planar.max(residual.magnitude());
            let c = cross(point, next);
            for j in 0..3 {
                area[j] = area[j] + c[j];
            }
        }
        let outward = if face.forward { n } else { n.map(|x| -x) };
        require(dot(area, outward).lo > 0.)?;
    }
    require(planar <= tolerance / 16.)?;
    let caps = thread.projection().cap_height_bounds_mm();
    require(low > caps[0][1])?;
    let contact = interval(caps[1]) - I::point(low);
    let sewing = I::point(contact.magnitude()) + I::point(planar);
    require(sewing.valid() && sewing.hi <= tolerance / 8.)?;
    let extent = I::point(high) - I::point(low);
    require(extent.valid() && extent.lo > minimum.max(32. * tolerance) + 2. * sewing.hi)?;
    // Projection already proves every actual tensor lies between the original
    // exact cap scalars and winds once. The monotone axial affine map fixing
    // the lower cap and taking that SAME upper cap to `low` moves every point
    // by at most this collar. It does not replace any retained raw coefficient.
    // That declared sewn reference has thread z<=low and head z>=low, with
    // disjoint open interiors and the SAME original rim's boundary ownership.
    let wall_floor = I::point(minimum) + I::point(16. * tolerance) + sewing;
    let mut wall = f64::INFINITY;
    // These are the literal scalar intervals, not their rounded value caches.
    // Every actual spline/tensor image belongs to its Bernstein control hull;
    // checking the complete original point table is conservative for all rims.
    for point in &parsed.points {
        let xy = [
            interval(point.coordinates[0].bounds()),
            interval(point.coordinates[1].bounds()),
        ];
        for i in 0..6 {
            work.charge(24)?;
            let (s, length) = side(polygon[i], polygon[(i + 1) % 6], xy);
            let distance = s / length;
            require(distance.valid() && distance.lo > wall_floor.hi)?;
            wall = wall.min(distance.lo);
        }
    }
    Ok(Geometry {
        polygon,
        stations: [low, high],
        planar_sewing: planar,
        join_sewing: sewing.hi,
        minimum_wall: wall,
    })
}
