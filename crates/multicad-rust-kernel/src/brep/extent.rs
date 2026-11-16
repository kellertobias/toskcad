//! Axis-aligned bounds of the trimmed boundary.
//!
//! Edge curves supply their own extrema. A curved face can also reach
//! further than its edges at an interior point where its normal is parallel
//! to an axis; those points are found on the exact surface and kept when
//! they lie inside the trim region.
use super::domain::Domain;
use super::polygon;
use super::{Boundary, Point, Result, Surface};

fn include(bounds: &mut [Point; 2], point: Point) {
    for axis in 0..3 {
        bounds[0][axis] = bounds[0][axis].min(point[axis]);
        bounds[1][axis] = bounds[1][axis].max(point[axis]);
    }
}

/// A stationary point of one coordinate near `seed`, by Newton's method on
/// the coordinate's gradient.
fn stationary(surface: &Surface, axis: usize, seed: [f64; 2], span: [f64; 2]) -> Option<[f64; 2]> {
    let gradient = |uv: [f64; 2]| {
        let (_, [su, sv]) = surface.evaluate(uv);
        [su[axis], sv[axis]]
    };
    let mut uv = seed;
    for _ in 0..40 {
        let g = gradient(uv);
        let h = [1e-6 * span[0], 1e-6 * span[1]];
        let column = |k: usize| {
            let (mut a, mut b) = (uv, uv);
            a[k] -= h[k];
            b[k] += h[k];
            let (low, high) = (gradient(a), gradient(b));
            [
                (high[0] - low[0]) / (2. * h[k]),
                (high[1] - low[1]) / (2. * h[k]),
            ]
        };
        let (cu, cv) = (column(0), column(1));
        let determinant = cu[0] * cv[1] - cv[0] * cu[1];
        let size = (cu[0].abs() + cv[0].abs()) * (cu[1].abs() + cv[1].abs());
        if !determinant.is_finite() {
            return None;
        }
        let step = if determinant.abs() > 1e-9 * size {
            [
                (g[0] * cv[1] - cv[0] * g[1]) / determinant,
                (cu[0] * g[1] - g[0] * cu[1]) / determinant,
            ]
        } else {
            // On a surface of revolution a coordinate along the axis is
            // stationary along a whole parallel, such as a torus's top
            // circle: step along the parameter whose rate still changes.
            let (du, dv) = (cu[0].abs(), cv[1].abs());
            if du.max(dv) == 0. || cu[1].abs() + cv[0].abs() > 1e-9 * du.max(dv) {
                return None;
            }
            if du >= dv {
                [g[0] / cu[0], 0.]
            } else {
                [0., g[1] / cv[1]]
            }
        };
        uv = [uv[0] - step[0], uv[1] - step[1]];
        if (uv[0] - seed[0]).abs() > 2. * span[0] || (uv[1] - seed[1]).abs() > 2. * span[1] {
            return None;
        }
        if step[0].abs() <= 1e-13 * span[0] && step[1].abs() <= 1e-13 * span[1] {
            // Both rates vanish at a stationary point.
            let (_, [su, sv]) = surface.evaluate(uv);
            let g = gradient(uv);
            let scale = super::norm(su) + super::norm(sv);
            return (g[0].abs() + g[1].abs() <= 1e-9 * scale).then_some(uv);
        }
    }
    None
}

fn interior(
    surface: &Surface,
    domain: &Domain,
    face: usize,
    ends: &[[usize; 2]],
    bounds: &mut [Point; 2],
) -> Result<()> {
    let triangles = polygon::triangulate(domain, surface, face, ends)?;
    let mut low = [f64::INFINITY; 2];
    let mut high = [f64::NEG_INFINITY; 2];
    for corner in triangles.iter().flatten() {
        for axis in 0..2 {
            low[axis] = low[axis].min(corner.uv[axis]);
            high[axis] = high[axis].max(corner.uv[axis]);
        }
    }
    let span = [high[0] - low[0], high[1] - low[1]];
    for axis in 0..3 {
        for i in 0..=6 {
            for j in 0..=6 {
                let seed = [
                    low[0] + span[0] * f64::from(i) / 6.,
                    low[1] + span[1] * f64::from(j) / 6.,
                ];
                let found = stationary(surface, axis, seed, span)
                    .filter(|&uv| polygon::contains(&triangles, uv, surface.periods()));
                if let Some(uv) = found {
                    include(bounds, surface.evaluate(uv).0);
                }
            }
        }
    }
    Ok(())
}

pub(crate) fn bounds(boundary: &Boundary, domains: &[Domain]) -> Result<[Point; 2]> {
    let mut bounds = [[f64::INFINITY; 3], [f64::NEG_INFINITY; 3]];
    let ends: Vec<[usize; 2]> = boundary.edges.iter().map(|edge| edge.vertices).collect();
    for edge in &boundary.edges {
        let mut low = [0.; 3];
        let mut high = [0.; 3];
        for axis in 0..3 {
            [low[axis], high[axis]] = edge.curve.extent(edge.range[0], edge.range[1], axis);
        }
        include(&mut bounds, low);
        include(&mut bounds, high);
    }
    for (index, (face, domain)) in boundary.faces.iter().zip(domains).enumerate() {
        if let Domain::Trimmed {
            pole: Some((v, _)), ..
        } = domain
        {
            include(&mut bounds, face.surface.evaluate([0., *v]).0);
        }
        let ruled = matches!(
            face.surface,
            Surface::Plane(_)
                | Surface::Cylinder { .. }
                | Surface::Cone { .. }
                | Surface::Extrusion { .. }
        );
        if !ruled {
            interior(&face.surface, domain, index, &ends, &mut bounds)?;
        }
    }
    Ok(bounds)
}
