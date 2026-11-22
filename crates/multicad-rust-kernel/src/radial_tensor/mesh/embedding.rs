//! Literal embedding of the emitted f64 product grid, after all rounding.
use super::*;
use crate::display::DisplayMesh;
pub(super) enum Verdict {
    Certified,
    HeightUncertain,
    AngularUncertain(usize),
    RadialUncertain(usize),
}
fn vertex(mesh: &DisplayMesh, index: usize) -> [I; 3] {
    std::array::from_fn(|k| I::point(mesh.vertices[3 * index + k]))
}
pub(super) fn certify(
    mesh: &DisplayMesh,
    ring: usize,
    rows: usize,
    work: &mut Work<'_>,
) -> Result<Verdict, Error> {
    let centers = [vertex(mesh, ring * rows), vertex(mesh, ring * rows + 1)];
    let delta = std::array::from_fn(|k| centers[1][k] - centers[0][k]);
    let length = dot(delta, delta).sqrt().positive()?;
    let axis = delta.map(|v| v / length);
    let reference = (0..3)
        .min_by(|&a, &b| axis[a].magnitude().total_cmp(&axis[b].magnitude()))
        .expect("three coordinate axes");
    let reference = std::array::from_fn(|k| I::point(if k == reference { 1. } else { 0. }));
    let transverse = cross(axis, reference);
    let scale = dot(transverse, transverse).sqrt().positive()?;
    let v = transverse.map(|p| p / scale);
    let u = cross(v, axis);
    let mut bottom_max = f64::NEG_INFINITY;
    let mut top_min = f64::INFINITY;
    for end in 0..2 {
        let mut sum = I::point(0.);
        for i in 0..ring {
            work.charge(16)?;
            let at = |i: usize| {
                let p = vertex(mesh, (if end == 0 { 0 } else { rows - 1 }) * ring + i);
                let d = std::array::from_fn(|k| p[k] - centers[0][k]);
                ([dot(d, u), dot(d, v), I::point(0.)], dot(d, axis))
            };
            let (a, z) = at(i);
            let (b, _) = at((i + 1) % ring);
            if end == 0 {
                if z.lo <= 0. {
                    return Ok(Verdict::HeightUncertain);
                }
                bottom_max = bottom_max.max(z.hi);
            } else {
                if z.hi >= length.lo {
                    return Ok(Verdict::HeightUncertain);
                }
                top_min = top_min.min(z.lo);
            }
            let Some(angle) = caps::phase(a, b) else {
                return Ok(Verdict::AngularUncertain(i));
            };
            sum = sum + angle;
        }
        let pi = I::new(
            std::f64::consts::PI.next_down(),
            std::f64::consts::PI.next_up(),
        );
        if sum.lo <= pi.hi || sum.hi >= (I::point(3.) * pi).lo {
            return Err(Error::Unresolved);
        }
    }
    if bottom_max >= top_min {
        return Ok(Verdict::HeightUncertain);
    }
    let side_triangles = 2 * (rows - 1) * ring;
    for triangle in 0..side_triangles {
        work.charge(32)?;
        let indices = std::array::from_fn::<_, 3, _>(|k| mesh.indices[3 * triangle + k] as usize);
        let points = indices.map(|i| vertex(mesh, i));
        let edges = [1, 2].map(|i| std::array::from_fn(|k| points[i][k] - points[0][k]));
        let normal = cross(edges[0], edges[1]);
        for point in points {
            let d = std::array::from_fn(|k| point[k] - centers[0][k]);
            let z = dot(d, axis);
            let radial = std::array::from_fn(|k| d[k] - axis[k] * z);
            if dot(normal, radial).checked()?.lo <= 0. {
                return Ok(Verdict::RadialUncertain(triangle));
            }
        }
    }
    // The exact product-grid annulus and its two cap fans are connected
    // oriented manifolds. Positive radial projection over every complete
    // triangle (the determinant is affine, hence vertex tests suffice)
    // excludes the axis and gives positive local degree. The two strictly
    // ordered, single-wrap boundary graphs fix the global degree to one,
    // so the side is a one-sheet radial graph. Cap rays are unique and lie
    // strictly below/above its boundary; the center-height tests separate
    // the caps as well. Connectivity or signed volume alone is not used
    // as an embedding certificate.
    Ok(Verdict::Certified)
}
fn cross(a: [I; 3], b: [I; 3]) -> [I; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
