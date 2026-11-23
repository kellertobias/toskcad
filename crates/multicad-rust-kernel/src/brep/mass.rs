//! Area, centroid and enclosed-volume integrals over a trimmed face.
//!
//! With `F` the integrand on the parameter plane and `G(u, v)` its integral
//! along u from a fixed base, Green's theorem turns the area integral into
//! the boundary integral of `G dv`. Each loop that winds about a periodic
//! surface contributes the full-period integral up to its start, and a face
//! that closes at a pole contributes it up to the pole.
use super::domain::{Domain, UvLoop};
use super::surface::PreparedInverse;
use super::{
    Boundary, Error, Face, GAUSS, NurbsSurface, Point, Result, Surface, cross, distance, dot,
    ensure, norm, scale, sub,
};
use std::cell::RefCell;
use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};

/// Volume flux, area and the three first moments of area.
pub(crate) type Moments = [f64; 5];

fn accumulate(into: &mut Moments, value: Moments, weight: f64) {
    for (slot, value) in into.iter_mut().zip(value) {
        *slot += value * weight;
    }
}

fn integrand(surface: &Surface, uv: [f64; 2], anchor: Point) -> Moments {
    moments(surface.evaluate(uv), anchor)
}

fn moments((s, [su, sv]): (Point, [Point; 2]), anchor: Point) -> Moments {
    let normal = cross(su, sv);
    let area = norm(normal);
    [
        dot(sub(s, anchor), normal) / 3.,
        area,
        s[0] * area,
        s[1] * area,
        s[2] * area,
    ]
}

/// Gauss-Legendre integral of `f` over panels, each cut `refine` times.
fn quadrature(
    panels: &[f64],
    refine: usize,
    mut f: impl FnMut(f64) -> Result<Moments>,
) -> Result<Moments> {
    let mut total = [0.; 5];
    for pair in panels.windows(2) {
        let width = (pair[1] - pair[0]) / refine as f64;
        for piece in 0..refine {
            let middle = pair[0] + width * (piece as f64 + 0.5);
            for (node, weight) in GAUSS {
                accumulate(
                    &mut total,
                    f(middle + 0.5 * width * node)?,
                    0.5 * width * weight,
                );
            }
        }
    }
    Ok(total)
}

/// Bounds the optional memo, not the integral's total workspace.
const MEMO_BYTES: usize = 4 * 1024 * 1024;

/// Multiplicative hash of a parameter's bits; keys are not adversarial.
#[derive(Default)]
struct NodeHasher(u64);

impl Hasher for NodeHasher {
    fn finish(&self) -> u64 {
        self.0
    }

    fn write(&mut self, bytes: &[u8]) {
        for &byte in bytes {
            self.write_u64(self.0 ^ u64::from(byte));
        }
    }

    fn write_u64(&mut self, bits: u64) {
        let mixed = (self.0 ^ bits).wrapping_mul(0x9e37_79b9_7f4a_7c15);
        self.0 = mixed ^ (mixed >> 32);
    }
}

/// The u-axis spline bases already evaluated, keyed by the exact bits of u.
///
/// Every strip starts at the same base and cuts its panels at the same knot
/// breaks with the same arithmetic, so Gauss nodes recur bit for bit across
/// strips. The basis is a pure function of u, so a memoized entry is the
/// value a fresh evaluation returns.
#[derive(Default)]
struct UBasis {
    index: HashMap<u64, usize, BuildHasherDefault<NodeHasher>>,
    spans: Vec<usize>,
    /// Per entry, the `degree + 1` values then the `degree + 1` derivatives.
    values: Vec<f64>,
}

struct Integral<'a> {
    surface: &'a Surface,
    inverse: PreparedInverse<'a>,
    anchor: Point,
    base: [f64; 2],
    refine: usize,
    memo: RefCell<UBasis>,
}

impl Integral<'_> {
    /// `G(u, v)`: the integrand accumulated along u from the base.
    fn strip(&self, [u, v]: [f64; 2]) -> Moments {
        let (low, high, sign) = if u < self.base[0] {
            (u, self.base[0], -1.)
        } else {
            (self.base[0], u, 1.)
        };
        let panels = self.surface.breaks(0, low, high);
        let total = if let Surface::Nurbs { surface, closed } = self.surface {
            // `Surface::evaluate` confines each axis independently, then
            // combines the two bases; v and its basis are shared here.
            let v = Surface::confine_axis(surface, *closed, 1, v);
            let (sv, nv, dv) = surface.basis(1, v);
            // The section of control rows at this v, formed on first use.
            let mut section = vec![None; surface.counts()[0]];
            quadrature(&panels, self.refine, |s| {
                let u = Surface::confine_axis(surface, *closed, 0, s);
                Ok(self.with_basis(surface, u, |along| {
                    let row = |i: usize| {
                        *section[i].get_or_insert_with(|| surface.section_row(i, (sv, &nv, &dv)))
                    };
                    moments(surface.combine_section(along, row), self.anchor)
                }))
            })
        } else {
            quadrature(&panels, self.refine, |s| {
                Ok(integrand(self.surface, [s, v], self.anchor))
            })
        };
        total.unwrap_or_default().map(|g| g * sign)
    }

    /// Call `f` with the u-axis basis at `u`, memoized while room remains.
    fn with_basis<R>(
        &self,
        surface: &NurbsSurface,
        u: f64,
        f: impl FnOnce((usize, &[f64], &[f64])) -> R,
    ) -> R {
        let count = surface.degrees()[0] + 1;
        let mut memo = self.memo.borrow_mut();
        if let Some(&entry) = memo.index.get(&u.to_bits()) {
            let values = &memo.values[2 * count * entry..2 * count * (entry + 1)];
            let (value, derivative) = values.split_at(count);
            return f((memo.spans[entry], value, derivative));
        }
        let (span, value, derivative) = surface.basis(0, u);
        // Index, span and both value runs per entry, with map overhead.
        let entry_bytes = 3 * size_of::<u64>() + 2 * count * size_of::<f64>();
        if (memo.spans.len() + 1) * entry_bytes <= MEMO_BYTES {
            let entry = memo.spans.len();
            memo.index.insert(u.to_bits(), entry);
            memo.spans.push(span);
            memo.values.extend_from_slice(&value[..count]);
            memo.values.extend_from_slice(&derivative[..count]);
        }
        f((span, &value, &derivative))
    }

    /// `H(v)`: the full-period strip accumulated along v from the base.
    fn band(&self, v: f64) -> Moments {
        let period = self.surface.periods()[0].unwrap_or(0.);
        let (low, high, sign) = if v < self.base[1] {
            (v, self.base[1], -1.)
        } else {
            (self.base[1], v, 1.)
        };
        let panels = self.surface.breaks(1, low, high);
        quadrature(&panels, self.refine, |w| {
            Ok(self.strip([self.base[0] + period, w]))
        })
        .unwrap_or_default()
        .map(|h| h * sign)
    }

    fn boundary(&self, boundary: &Boundary, trim: &UvLoop, tolerance: f64) -> Result<Moments> {
        let periods = self.surface.periods();
        let mut total = [0.; 5];
        for (i, segment) in trim.segments.iter().enumerate() {
            let Some(segment) = segment else { continue };
            let curve = &boundary.edges[segment.edge].curve;
            let ends = [trim.at(i, periods), trim.at(i + 1, periods)];
            let piece = quadrature(&segment.t, self.refine, |t| {
                let (point, tangent) = curve.evaluate(t);
                let along = (t - segment.t[0]) / (segment.t[1] - segment.t[0]);
                let seed =
                    [0, 1].map(|axis| ends[0][axis] + (ends[1][axis] - ends[0][axis]) * along);
                let uv = self.inverse.invert(point, Some(seed));
                let (s, [su, sv]) = self.surface.evaluate(uv);
                ensure(distance(s, point) <= tolerance, Error::Invalid)?;
                let (a, b, c) = (dot(su, su), dot(su, sv), dot(sv, sv));
                let (ru, rv) = (dot(tangent, su), dot(tangent, sv));
                let determinant = a * c - b * b;
                let rate = if determinant > 1e-12 * a * c {
                    (a * rv - b * ru) / determinant
                } else {
                    // Where the u derivative stalls, as at a spline end whose
                    // parametrization has zero speed, any motion along u
                    // moves the point by less than the tolerance the volume
                    // error already counts, and a tangent along the v
                    // derivative moves at the rate T.Sv/|Sv|^2. Parallel
                    // derivatives that do not stall leave it undefined.
                    let rate = rv / c;
                    let residual = norm(sub(tangent, scale(sv, rate)));
                    ensure(
                        a.sqrt() <= tolerance && c > 0. && residual <= 1e-9 * norm(tangent),
                        Error::Unsupported,
                    )?;
                    rate
                };
                Ok(self.strip(uv).map(|g| g * rate))
            })?;
            accumulate(&mut total, piece, 1.);
        }
        if trim.winding[0] != 0 {
            accumulate(
                &mut total,
                self.band(trim.points[0].uv[1]),
                -f64::from(trim.winding[0]),
            );
        }
        Ok(total)
    }
}

/// Integrals over the face region, taken counter-clockwise about the
/// surface normal: `[flux, area, moment x, moment y, moment z]`.
pub(crate) fn integrate(
    boundary: &Boundary,
    face: &Face,
    domain: &Domain,
    anchor: Point,
    refine: usize,
    tolerance: f64,
) -> Result<Moments> {
    let surface = &face.surface;
    let mut integral = Integral {
        surface,
        inverse: PreparedInverse::new(surface),
        anchor,
        base: [0.; 2],
        refine,
        memo: RefCell::default(),
    };
    match domain {
        Domain::Complete([low, high]) => {
            integral.base = *low;
            Ok(integral.band(high[1]))
        }
        Domain::Trimmed { loops, pole } => {
            let start = loops[0].points[0].uv;
            let low = loops
                .iter()
                .flat_map(|trim| &trim.points)
                .fold(f64::INFINITY, |low, point| low.min(point.uv[0]));
            // Winding loops share the first loop's period window.
            integral.base = [
                if loops[0].winding[0] != 0 {
                    start[0]
                } else {
                    low
                },
                start[1],
            ];
            let mut total = [0.; 5];
            for trim in loops {
                accumulate(
                    &mut total,
                    integral.boundary(boundary, trim, tolerance)?,
                    1.,
                );
            }
            if let Some((v, _)) = pole {
                accumulate(
                    &mut total,
                    integral.band(*v),
                    f64::from(loops[0].winding[0]),
                );
            }
            Ok(total)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use multicad_precision::PortableMath as _;

    fn patch() -> NurbsSurface {
        let counts = [7, 5];
        let points = (0..counts[0])
            .flat_map(|i| {
                (0..counts[1]).map(move |j| {
                    let (i, j) = (i as f64, j as f64);
                    [
                        i + 0.1 * j * j,
                        j - 0.05 * i,
                        0.3 * (i * 0.7).portable_sin() + 0.2 * j,
                    ]
                })
            })
            .collect::<Vec<_>>();
        let weights = (0..points.len())
            .map(|k| 1. + 0.25 * ((k * 7 % 5) as f64))
            .collect();
        NurbsSurface::new(
            [3, 2],
            counts,
            points,
            Some(weights),
            [vec![0., 0.2, 0.45, 0.7, 1.], vec![0., 0.3, 0.6, 1.]],
            [vec![4, 1, 1, 1, 4], vec![3, 1, 1, 3]],
        )
        .unwrap()
    }

    /// Higham's `gamma_n = n u / (1 - n u)` with the unit roundoff `u`.
    fn gamma(n: usize) -> f64 {
        let nu = n as f64 * f64::EPSILON / 2.;
        nu / (1. - nu)
    }

    /// A first-order bound of one node's `|fresh - memoized|` integrand and
    /// the node's integrand magnitude, derived from the summations alone.
    ///
    /// Both paths form the homogeneous sums `sum_ij a_i b_j C_ij` from
    /// `m = (p + 1)(q + 1)` two-factor products: the fresh path in one chain
    /// (within `gamma_{m+1}` of exact), the section path as inner v sums then
    /// outer u sums (within `gamma_{p+q+2}`, no larger), each relative to
    /// `sum_ij |a_i||b_j||C_ij| <= max|C| sum|a| sum|b|`. So they differ by at
    /// most twice that. The difference is carried to first order through the
    /// rational projection and the moments, and both paths' own roundings of
    /// each later chain of `k` operations add `2 gamma_k` times the
    /// magnitudes that chain combines. Second-order terms are `O(gamma^2)`.
    fn node_rounding(integral: &Integral, [s, v]: [f64; 2]) -> (Moments, Moments) {
        let Surface::Nurbs { surface, closed } = integral.surface else {
            unreachable!("the memo only serves splines")
        };
        let [p, q] = surface.degrees();
        let eta = 2. * gamma((p + 1) * (q + 1) + 1);
        let mut largest = [0_f64; 4];
        for (k, point) in surface.points().iter().enumerate() {
            let w = surface.weights().map_or(1., |weights| weights[k]);
            for c in 0..3 {
                largest[c] = largest[c].max((point[c] * w).abs());
            }
            largest[3] = largest[3].max(w.abs());
        }
        let absolute = |values: &[f64], degree: usize| -> f64 {
            values[..=degree].iter().map(|x| x.abs()).sum()
        };
        let (sv, nv, dv) = surface.basis(1, Surface::confine_axis(surface, *closed, 1, v));
        let (su, nu, du) = surface.basis(0, Surface::confine_axis(surface, *closed, 0, s));
        let (a0, a1, b0, b1) = (
            absolute(&nu, p),
            absolute(&du, p),
            absolute(&nv, q),
            absolute(&dv, q),
        );
        let error = [a0 * b0, a1 * b0, a0 * b1].map(|m| largest.map(|c| eta * c * m));
        let mut sums = [[0.; 4]; 3];
        for i in 0..=p {
            let [row, along_v] = surface.section_row(su + i - p, (sv, &nv, &dv));
            for c in 0..4 {
                sums[0][c] += nu[i] * row[c];
                sums[1][c] += du[i] * row[c];
                sums[2][c] += nu[i] * along_v[c];
            }
        }
        let w = sums[0][3];
        let point: Point = std::array::from_fn(|c| sums[0][c] / w);
        // `S = P / w`: one division per path.
        let d_point: Point = std::array::from_fn(|c| {
            (error[0][c] + point[c].abs() * error[0][3]) / w + 2. * gamma(1) * point[c].abs()
        });
        // `(d - S d_w) (1 / w)`: five operations on `|d| + |S||d_w|` per path.
        let derivative = |k: usize| -> (Point, Point) {
            let (d, e) = (sums[k], error[k]);
            let value = std::array::from_fn(|c| (d[c] - point[c] * d[3]) / w);
            let bound = std::array::from_fn(|c| {
                let size = (d[c].abs() + point[c].abs() * d[3].abs()) / w;
                (e[c] + point[c].abs() * e[3] + d[3].abs() * d_point[c]) / w
                    + size * error[0][3] / w
                    + 2. * gamma(5) * size
            });
            (value, bound)
        };
        let ((along_u, d_along_u), (along_v, d_along_v)) = (derivative(1), derivative(2));
        // Each normal component is a difference of two products.
        let pair = |a: usize, b: usize| -> (f64, f64) {
            let value = along_u[a] * along_v[b] - along_u[b] * along_v[a];
            let size = (along_u[a] * along_v[b]).abs() + (along_u[b] * along_v[a]).abs();
            let bound = d_along_u[a] * along_v[b].abs()
                + along_u[a].abs() * d_along_v[b]
                + d_along_u[b] * along_v[a].abs()
                + along_u[b].abs() * d_along_v[a]
                + 2. * gamma(2) * size;
            (value, bound)
        };
        let [x, y, z] = [pair(1, 2), pair(2, 0), pair(0, 1)];
        let (normal, d_normal) = ([x.0, y.0, z.0], [x.1, y.1, z.1]);
        // The norm: three squares, two additions and a square root.
        let area = norm(normal);
        let d_area = d_normal.iter().sum::<f64>() + 2. * gamma(4) * area;
        // The flux: a subtraction, a three-term dot product and a division.
        let offset = sub(point, integral.anchor);
        let d_flux = ((0..3)
            .map(|c| d_point[c] * normal[c].abs() + offset[c].abs() * d_normal[c])
            .sum::<f64>()
            + 2. * gamma(5) * (0..3).map(|c| (offset[c] * normal[c]).abs()).sum::<f64>())
            / 3.;
        let first = |c: usize| {
            d_point[c] * area + point[c].abs() * d_area + 2. * gamma(1) * point[c].abs() * area
        };
        let value = moments((point, [along_u, along_v]), integral.anchor);
        (
            [d_flux, d_area, first(0), first(1), first(2)],
            value.map(f64::abs),
        )
    }

    /// The bound of `|fresh - memoized|` for the strip to `[u, v]`: the
    /// weighted node bounds plus both paths' Gauss sums over `N` nodes,
    /// each within `gamma_{N+1} sum |weight g|` of its exact sum.
    fn strip_rounding_bound(integral: &Integral, [u, v]: [f64; 2]) -> Moments {
        let base = integral.base[0];
        let (low, high) = if u < base { (u, base) } else { (base, u) };
        let panels = integral.surface.breaks(0, low, high);
        let nodes = (panels.len() - 1) * integral.refine * GAUSS.len();
        // Gauss weights and panel widths are positive.
        let weighted = |part: usize| {
            quadrature(&panels, integral.refine, |s| {
                let (bound, size) = node_rounding(integral, [s, v]);
                Ok([bound, size][part])
            })
            .unwrap()
        };
        let (bound, size) = (weighted(0), weighted(1));
        std::array::from_fn(|k| bound[k] + 2. * gamma(nodes + 1) * size[k])
    }

    /// Asserts the memoized strip agrees with a fresh evaluation within the
    /// rounding bound.
    fn assert_within_rounding(fresh: Moments, memoized: Moments, bound: Moments) {
        for k in 0..5 {
            let difference = (fresh[k] - memoized[k]).abs();
            assert!(
                difference <= bound[k],
                "moment {k}: {} vs {} differ by {difference:e}, bound {:e}",
                fresh[k],
                memoized[k],
                bound[k]
            );
        }
    }

    #[test]
    fn memoized_strips_match_fresh_surface_evaluation_within_rounding() {
        for closed in [[false, false], [true, false], [true, true]] {
            let surface = Surface::Nurbs {
                surface: patch(),
                closed,
            };
            for refine in [1, 2] {
                for base in [0., 0.31] {
                    let integral = Integral {
                        surface: &surface,
                        inverse: PreparedInverse::new(&surface),
                        anchor: [0.4, -0.2, 0.1],
                        base: [base, 0.],
                        refine,
                        memo: RefCell::default(),
                    };
                    let mut nodes = 0;
                    // Repeated and wrapped u, u below the base and v outside
                    // the open range all agree with the uncached evaluation.
                    for &u in &[0.9, 0.05, 0.9, 0.62, 1.3, -0.2, 0.62, 0.31] {
                        for &v in &[0.0, 0.17, 0.5, 0.83, 1.0, 1.4, -0.1] {
                            let (low, high, sign) = if u < base {
                                (u, base, -1.)
                            } else {
                                (base, u, 1.)
                            };
                            let panels = surface.breaks(0, low, high);
                            nodes += (panels.len() - 1) * refine * GAUSS.len();
                            let fresh = quadrature(&panels, refine, |s| {
                                Ok(integrand(&surface, [s, v], integral.anchor))
                            })
                            .unwrap()
                            .map(|g| g * sign);
                            let memoized = integral.strip([u, v]);
                            let bound = strip_rounding_bound(&integral, [u, v]);
                            assert_within_rounding(fresh, memoized, bound);
                        }
                    }
                    let entries = integral.memo.borrow().spans.len();
                    // Most nodes recur across strips and are served by the memo.
                    assert!(entries > 0 && 2 * entries < nodes, "{entries} of {nodes}");
                }
            }
        }
    }

    #[test]
    fn a_full_memo_keeps_evaluating_without_storing() {
        let surface = Surface::Nurbs {
            surface: patch(),
            closed: [false; 2],
        };
        let integral = Integral {
            surface: &surface,
            inverse: PreparedInverse::new(&surface),
            anchor: [0.; 3],
            base: [0.; 2],
            refine: 2,
            memo: RefCell::default(),
        };
        let first = integral.strip([0.8, 0.4]);
        let stored = integral.memo.borrow().spans.len();
        // Pretend the budget is spent: further nodes are evaluated, not kept.
        {
            let mut memo = integral.memo.borrow_mut();
            let count = 4;
            let entry_bytes = 3 * size_of::<u64>() + 2 * count * size_of::<f64>();
            memo.spans.resize(MEMO_BYTES / entry_bytes, 0);
        }
        let full = integral.memo.borrow().spans.len();
        let again = integral.strip([0.9, 0.4]);
        assert!(stored > 0);
        assert_eq!(integral.memo.borrow().spans.len(), full);
        assert_eq!(integral.memo.borrow().index.len(), stored);
        assert_ne!(first, again);
        let fresh = quadrature(&surface.breaks(0, 0., 0.9), 2, |s| {
            Ok(integrand(&surface, [s, 0.4], [0.; 3]))
        })
        .unwrap();
        assert_within_rounding(fresh, again, strip_rounding_bound(&integral, [0.9, 0.4]));
    }
}
