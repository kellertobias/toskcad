//! Original C11 positions and C21 normals. Gaussian remainder derivatives
//! enclose each complete cell; no sampled derivative is a certificate.
use super::*;
use crate::tapered_thread::{Scalar as _, Series};
type V = [I; 3];
type S = Series<9>;
const ZERO: I = I { lo: 0., hi: 0. };
const VZERO: V = [ZERO; 3];
#[derive(Clone)]
struct Cell {
    position: [V; 144],
    normal: [V; 484],
    weight: f64,
    scale: [f64; 2],
    depth: usize,
}
fn basis<const N: usize>(x: I, work: &mut Work<'_>) -> Result<[I; N]> {
    let b = binomials::<N>(work)?;
    let mut a = [p(1.); N];
    let mut c = a;
    for i in 1..N {
        work.step(3)?;
        a[i] = a[i - 1] * x;
        c[i] = c[i - 1] * (p(1.) - x);
    }
    work.step(2 * N)?;
    Ok(std::array::from_fn(|i| p(b[i]) * a[i] * c[N - 1 - i]))
}
// Reuse the first tensor contraction across all four opposite-axis nodes.
// This evaluates the same original polynomial; no quadrature node is dropped.
fn rows<const N: usize>(net: &[V], u: [I; N], work: &mut Work<'_>) -> Result<[V; N]> {
    let mut rows = [VZERO; N];
    for i in 0..N {
        for j in 0..N {
            work.step(6)?;
            for a in 0..3 {
                rows[j][a] = rows[j][a] + u[i] * net[i * N + j][a];
            }
        }
    }
    Ok(rows)
}
fn point<const N: usize>(rows: &[V; N], v: [I; N], work: &mut Work<'_>) -> Result<V> {
    let mut out = VZERO;
    for j in 0..N {
        work.step(6)?;
        for a in 0..3 {
            out[a] = out[a] + rows[j][a] * v[j];
        }
    }
    Ok(out)
}
fn jet(net: &[V], degree: usize, axis: usize, work: &mut Work<'_>) -> Result<[S; 3]> {
    let n = degree + 1;
    let mut current = [VZERO; 484];
    current[..net.len()].copy_from_slice(net);
    let mut sizes = [n, n];
    let mut result = [S::interval(ZERO); 3];
    let mut factor = 1.;
    for k in 0..9 {
        work.step(6 * sizes[0] * sizes[1] + 8)?;
        let hull: V = std::array::from_fn(|a| {
            I::new(
                current[..sizes[0] * sizes[1]]
                    .iter()
                    .map(|v| v[a].lo)
                    .fold(f64::INFINITY, f64::min),
                current[..sizes[0] * sizes[1]]
                    .iter()
                    .map(|v| v[a].hi)
                    .fold(f64::NEG_INFINITY, f64::max),
            )
        });
        for a in 0..3 {
            result[a].c[k] = hull[a] * p(factor);
        }
        if k == 8 {
            break;
        }
        let [nu, nv] = sizes;
        if axis == 0 {
            for i in 0..nu - 1 {
                for j in 0..nv {
                    work.step(3)?;
                    current[i * nv + j] = std::array::from_fn(|a| {
                        current[(i + 1) * nv + j][a] - current[i * nv + j][a]
                    });
                }
            }
        } else {
            for i in 0..nu {
                for j in 0..nv - 1 {
                    work.step(3)?;
                    current[i * (nv - 1) + j] = std::array::from_fn(|a| {
                        current[i * nv + j + 1][a] - current[i * nv + j][a]
                    });
                }
            }
        }
        sizes[axis] -= 1;
        // binomial(degree,k+1), an exactly representable small integer.
        factor = factor * (degree - k) as f64 / (k + 1) as f64;
    }
    Ok(result)
}
fn density(position: [S; 3], normal: [S; 3]) -> Result<[S; 4]> {
    let mut square = normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2];
    square.c[0] = normal[0].c[0].square() + normal[1].c[0].square() + normal[2].c[0].square();
    let area = square.sqrt().checked().map_err(|_| Error::Unsupported)?;
    Ok([
        area,
        position[0] * area,
        position[1] * area,
        position[2] * area,
    ])
}
struct Gaussian {
    weights: [I; 4],
    p: [[I; 12]; 4],
    n: [[I; 22]; 4],
}
impl Gaussian {
    fn new(work: &mut Work<'_>) -> Result<Self> {
        work.step(64)?;
        let root = (p(6.) / p(5.)).sqrt();
        let outer = ((p(3.) + p(2.) * root) / p(7.)).sqrt() / p(2.);
        let inner = ((p(3.) - p(2.) * root) / p(7.)).sqrt() / p(2.);
        let nodes = [
            p(0.5) - outer,
            p(0.5) - inner,
            p(0.5) + inner,
            p(0.5) + outer,
        ];
        let w0 = (p(18.) - p(30.).sqrt()) / p(72.);
        let w1 = (p(18.) + p(30.).sqrt()) / p(72.);
        let mut points = [[ZERO; 12]; 4];
        let mut normals = [[ZERO; 22]; 4];
        for i in 0..4 {
            points[i] = basis(nodes[i], work)?;
            normals[i] = basis(nodes[i], work)?;
        }
        Ok(Self {
            weights: [w0, w1, w1, w0],
            p: points,
            n: normals,
        })
    }
    fn enclosure(
        &self,
        cell: &Cell,
        original: &[[[S; 3]; 2]; 2],
        work: &mut Work<'_>,
    ) -> Result<[I; 4]> {
        let mut remainder = [ZERO; 4];
        for (axis, original_fields) in original.iter().enumerate() {
            // Differentiate the ORIGINAL coefficients once. Restriction can
            // widen coefficient intervals; repeatedly differencing those
            // intervals would create artificial high-order error floors.
            // Chain factors scale the original complete-root derivative hulls
            // into this child's normalized coordinate. Zeroth ranges use the
            // tighter actual restricted controls.
            let mut fields = *original_fields;
            for (field, net) in fields
                .iter_mut()
                .zip([&cell.position[..], &cell.normal[..]])
            {
                work.step(6 * net.len() + 36)?;
                let mut factor = p(1.);
                for k in 0..9 {
                    if k > 0 {
                        factor = factor * p(cell.scale[axis]);
                    }
                    for (a, value) in field.iter_mut().enumerate() {
                        if k == 0 {
                            value.c[0] = I::new(
                                net.iter().map(|x| x[a].lo).fold(f64::INFINITY, f64::min),
                                net.iter()
                                    .map(|x| x[a].hi)
                                    .fold(f64::NEG_INFINITY, f64::max),
                            );
                        } else {
                            value.c[k] = value.c[k] * factor;
                        }
                    }
                }
            }
            let values = density(fields[0], fields[1])?;
            work.step(4096)?; // bounded order-eight products and square-root recurrence
            for i in 0..4 {
                remainder[i] = remainder[i] + p(values[i].c[8].magnitude()) / p(44100.);
            }
        }
        let mut values = [ZERO; 4];
        for i in 0..4 {
            let position_rows = rows(&cell.position, self.p[i], work)?;
            let normal_rows = rows(&cell.normal, self.n[i], work)?;
            for j in 0..4 {
                let pos = point(&position_rows, self.p[j], work)?;
                let n = point(&normal_rows, self.n[j], work)?;
                work.step(64)?;
                let area = (n[0].square() + n[1].square() + n[2].square()).sqrt();
                let f = [area, pos[0] * area, pos[1] * area, pos[2] * area];
                for k in 0..4 {
                    values[k] = values[k] + self.weights[i] * self.weights[j] * f[k];
                }
            }
        }
        for i in 0..4 {
            values[i] = (values[i] + I::new(-remainder[i].hi, remainder[i].hi)) * p(cell.weight);
        }
        if values
            .iter()
            .any(|v| !v.valid() || !v.lo.is_finite() || !v.hi.is_finite())
        {
            return Err(Error::Unsupported);
        }
        Ok(values)
    }
}
fn split<const K: usize>(
    net: &[V; K],
    n: usize,
    axis: usize,
    work: &mut Work<'_>,
) -> Result<[[V; K]; 2]> {
    let mut out = [[VZERO; K]; 2];
    for orthogonal in 0..n {
        let mut row = [VZERO; 22];
        let index = |along: usize| {
            if axis == 0 {
                along * n + orthogonal
            } else {
                orthogonal * n + along
            }
        };
        for i in 0..n {
            row[i] = net[index(i)];
        }
        out[0][index(0)] = row[0];
        out[1][index(n - 1)] = row[n - 1];
        for level in 1..n {
            for i in 0..n - level {
                work.step(6)?;
                row[i] = std::array::from_fn(|a| (row[i][a] + row[i + 1][a]) * p(0.5));
            }
            out[0][index(level)] = row[0];
            out[1][index(n - level - 1)] = row[n - level - 1];
        }
    }
    Ok(out)
}
/// Integrate one original unit-root cell. Normals already include the actual
/// STEP cofactor; positions are already mapped and centered by the caller.
pub(super) fn integrate(
    position: [V; 144],
    normal: &[V],
    area_width: f64,
    centroid_width: f64,
    extent: f64,
    work: &mut Work<'_>,
) -> Result<[I; 4]> {
    let normal: [V; 484] = normal.try_into().map_err(|_| Error::Invalid)?;
    let original = [
        [jet(&position, 11, 0, work)?, jet(&normal, 21, 0, work)?],
        [jet(&position, 11, 1, work)?, jet(&normal, 21, 1, work)?],
    ];
    // At most33 pending DFS cells, one active cell, two children, fixed jets
    // and tables fit2MiB. The immutable source/cache is charged by the caller.
    let mut pending = Vec::new();
    pending.try_reserve_exact(33).map_err(|_| Error::Limit)?;
    pending.push(Cell {
        position,
        normal,
        weight: 1.,
        scale: [1.; 2],
        depth: 0,
    });
    let g = Gaussian::new(work)?;
    let mut sum = [ZERO; 4];
    while let Some(cell) = pending.pop() {
        work.step(1)?;
        let value = g.enclosure(&cell, &original, work);
        let accepted = value.as_ref().is_ok_and(|f| {
            f[0].lo > 0.
                && f[0].width() <= area_width * cell.weight
                && f[0].width() <= centroid_width / (4. * extent.max(1.)) * f[0].lo
                && f[1..]
                    .iter()
                    .all(|m| m.width() <= centroid_width / 4. * f[0].lo)
        });
        if accepted {
            for (s, v) in sum.iter_mut().zip(value?) {
                *s = *s + v;
            }
            continue;
        }
        if let Err(error) = value
            && error != Error::Unsupported
        {
            return Err(error);
        }
        if cell.depth == 32 {
            return Err(Error::Unsupported);
        }
        let axis = cell.depth % 2;
        let [pa, pb] = split(&cell.position, 12, axis, work)?;
        let [na, nb] = split(&cell.normal, 22, axis, work)?;
        for (position, normal) in [(pa, na), (pb, nb)] {
            let mut scale = cell.scale;
            scale[axis] *= 0.5;
            pending.push(Cell {
                position,
                normal,
                weight: cell.weight * 0.5,
                scale,
                depth: cell.depth + 1,
            });
        }
    }
    Ok(sum)
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn original_polynomial_plane_integrals_keep_area_and_centered_moments() {
        let position = std::array::from_fn(|i| {
            [
                p(2.) * p((i / 12) as f64) / p(11.),
                p(3.) * p((i % 12) as f64) / p(11.),
                p(0.5),
            ]
        });
        let normal = [[p(0.), p(0.), p(6.)]; 484];
        let mut work = Work {
            left: 20_000_000,
            used: 0,
            bytes: 2 * 1024 * 1024,
            interrupt: None,
        };
        let values = integrate(position, &normal, 1e-7, 1e-7, 3., &mut work).unwrap();
        for (v, expected) in values.into_iter().zip([6., 6., 9., 3.]) {
            assert!(v.lo <= expected && v.hi >= expected, "{v:?} {expected}");
            assert!(v.width() < 1e-7);
        }
    }
    #[test]
    fn original_parabolic_cylinder_integrals_match_independent_closed_forms() {
        let position = std::array::from_fn(|i| {
            let u = i / 12;
            [
                p(u as f64) / p(11.),
                p((i % 12) as f64) / p(11.),
                p((u * u.saturating_sub(1)) as f64) / p(220.),
            ]
        });
        let normal: [V; 484] =
            std::array::from_fn(|i| [-p((i / 22) as f64) / p(21.), p(0.), p(1.)]);
        let mut work = Work {
            left: 20_000_000,
            used: 0,
            bytes: 2 * 1024 * 1024,
            interrupt: None,
        };
        let values = integrate(position, &normal, 1e-7, 1e-7, 1., &mut work).unwrap();
        // P=(u,v,u²/2), |Pu×Pv|=sqrt(1+u²). Elementary antiderivatives:
        // A=(sqrt2+asinh1)/2, Mx=(2sqrt2-1)/3, My=A/2,
        // Mz=(3sqrt2-asinh1)/16. Decimal constants computed independently.
        for (v, expected) in values.into_iter().zip([
            1.147793574696319_f64,
            0.609_475_708_248_73,
            0.5738967873481595,
            0.21007919375623388,
        ]) {
            assert!(
                v.lo < expected.next_down() && v.hi > expected.next_up(),
                "{v:?} expected={expected}"
            );
            assert!(v.width() < 1e-7);
        }
    }
    #[test]
    fn degenerate_original_normal_owner_and_real_resource_limits_refuse_publication() {
        let position = [VZERO; 144];
        let normal = [VZERO; 484];
        let mut work = Work {
            left: 20_000_000,
            used: 0,
            bytes: 2 * 1024 * 1024,
            interrupt: None,
        };
        assert!(matches!(
            integrate(position, &normal[..483], 1e-7, 1e-7, 1., &mut work),
            Err(Error::Invalid)
        ));
        assert!(matches!(
            integrate(position, &normal, 1e-7, 1e-7, 1., &mut work),
            Err(Error::Unsupported)
        ));
        let mut limited = Work {
            left: 1,
            used: 0,
            bytes: 2 * 1024 * 1024,
            interrupt: None,
        };
        assert!(matches!(
            integrate(position, &normal, 1e-7, 1e-7, 1., &mut limited),
            Err(Error::Limit)
        ));
    }
}
