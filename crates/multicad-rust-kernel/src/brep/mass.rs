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
use std::collections::{HashMap, HashSet};
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

/// Add the Gauss-Legendre integral of `f` over one panel cut `refine` times
/// to `total`, node by node exactly as [`quadrature`] does.
fn panel(total: &mut Moments, pair: [f64; 2], refine: usize, f: &mut impl FnMut(f64) -> Moments) {
    let width = (pair[1] - pair[0]) / refine as f64;
    for piece in 0..refine {
        let middle = pair[0] + width * (piece as f64 + 0.5);
        for (node, weight) in GAUSS {
            accumulate(total, f(middle + 0.5 * width * node), 0.5 * width * weight);
        }
    }
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
    /// The v of every strip taken so far, by its exact bits.
    seen: HashSet<u64, BuildHasherDefault<NodeHasher>>,
    /// Running strip totals at v values that recur, as along a boundary
    /// edge on one v row. See [`Prefix`].
    prefixes: HashMap<u64, Prefix, BuildHasherDefault<NodeHasher>>,
    /// Bytes held by `seen` and `prefixes`, bounded like the bases.
    prefix_bytes: usize,
}

/// The running total of strips from the base at one v, after each of its
/// leading whole panels.
///
/// A strip from the base sums its panels in order, node by node, and every
/// panel before its last one ends at a knot break that does not depend on
/// where the strip ends. Two strips at the same v therefore agree bit for
/// bit up to the shorter one's last whole panel, and the longer one may
/// continue from the stored total instead of summing that prefix again.
#[derive(Default)]
struct Prefix {
    /// The end of each stored whole panel.
    ends: Vec<f64>,
    /// The running total after each stored whole panel.
    totals: Vec<Moments>,
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
            let confined = Surface::confine_axis(surface, *closed, 1, v);
            let (sv, nv, dv) = surface.basis(1, confined);
            self.sum(&panels, sign, v, |s| {
                let u = Surface::confine_axis(surface, *closed, 0, s);
                self.with_basis(surface, u, |along| {
                    moments(surface.combine(along, (sv, &nv, &dv)), self.anchor)
                })
            })
        } else {
            self.sum(&panels, sign, v, |s| {
                integrand(self.surface, [s, v], self.anchor)
            })
        };
        total.map(|g| g * sign)
    }

    /// The quadrature of `f` over the strip `panels` at `v`, continuing a
    /// recurring v's stored running total where its panels agree. Strips
    /// below the base start with their partial panel and are summed afresh.
    fn sum(&self, panels: &[f64], sign: f64, v: f64, mut f: impl FnMut(f64) -> Moments) -> Moments {
        let key = v.to_bits();
        let recurring = sign > 0. && {
            let mut memo = self.memo.borrow_mut();
            let memo = &mut *memo;
            if memo.seen.contains(&key) {
                true
            } else {
                if memo.prefix_bytes + size_of::<u64>() <= MEMO_BYTES {
                    memo.seen.insert(key);
                    memo.prefix_bytes += size_of::<u64>();
                }
                false
            }
        };
        let mut total = [0.; 5];
        if !recurring {
            for pair in panels.windows(2) {
                panel(&mut total, [pair[0], pair[1]], self.refine, &mut f);
            }
            return total;
        }
        // Taken out while summing: `f` borrows the memo for the bases.
        let mut prefix = self
            .memo
            .borrow_mut()
            .prefixes
            .remove(&key)
            .unwrap_or_default();
        let whole = panels.len() - 2;
        let agreed = prefix
            .ends
            .iter()
            .zip(&panels[1..=whole])
            .take_while(|(stored, end)| stored.to_bits() == end.to_bits())
            .count();
        if agreed > 0 {
            total = prefix.totals[agreed - 1];
        }
        // Only a prefix that every stored panel agrees with is extended.
        let mut extend = agreed == prefix.ends.len();
        let entry_bytes = size_of::<f64>() + size_of::<Moments>();
        for k in agreed..whole {
            panel(&mut total, [panels[k], panels[k + 1]], self.refine, &mut f);
            if extend {
                let mut memo = self.memo.borrow_mut();
                if memo.prefix_bytes + entry_bytes <= MEMO_BYTES {
                    memo.prefix_bytes += entry_bytes;
                    prefix.ends.push(panels[k + 1]);
                    prefix.totals.push(total);
                } else {
                    extend = false;
                }
            }
        }
        panel(
            &mut total,
            [panels[whole], panels[whole + 1]],
            self.refine,
            &mut f,
        );
        self.memo.borrow_mut().prefixes.insert(key, prefix);
        total
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
                let (uv, (s, [su, sv])) = self.inverse.invert_evaluated(point, Some(seed));
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

    #[test]
    fn memoized_strips_are_bit_identical_to_fresh_surface_evaluation() {
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
                    // the open range all match the uncached evaluation.
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
                            assert_eq!(fresh.map(f64::to_bits), memoized.map(f64::to_bits));
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
        assert_eq!(fresh.map(f64::to_bits), again.map(f64::to_bits));
    }

    #[test]
    fn strips_at_a_recurring_v_continue_its_prefix_bit_for_bit() {
        for closed in [[false, false], [true, false]] {
            let surface = Surface::Nurbs {
                surface: patch(),
                closed,
            };
            for refine in [1, 2] {
                let integral = Integral {
                    surface: &surface,
                    inverse: PreparedInverse::new(&surface),
                    anchor: [0.4, -0.2, 0.1],
                    base: [0.05, 0.],
                    refine,
                    memo: RefCell::default(),
                };
                // Longer, shorter, equal, on-break and wrapped strips at two
                // recurring rows, interleaved, and one strip below the base.
                let us = [0.3, 0.95, 0.2, 0.95, 0.45, 0.7, 1.0, 0.02, 1.35, 0.6, 0.3];
                for &u in &us {
                    for v in [0., 0.6] {
                        let (low, high, sign) = if u < 0.05 {
                            (u, 0.05, -1.)
                        } else {
                            (0.05, u, 1.)
                        };
                        let fresh = quadrature(&surface.breaks(0, low, high), refine, |s| {
                            Ok(integrand(&surface, [s, v], integral.anchor))
                        })
                        .unwrap()
                        .map(|g| g * sign);
                        let continued = integral.strip([u, v]);
                        assert_eq!(fresh.map(f64::to_bits), continued.map(f64::to_bits));
                    }
                }
                let memo = integral.memo.borrow();
                assert_eq!(memo.seen.len(), 2);
                let stored: usize = memo.prefixes.values().map(|p| p.ends.len()).sum();
                assert!(stored > 4, "{stored}");
                assert!(memo.prefix_bytes <= MEMO_BYTES);
            }
        }
    }

    #[test]
    fn a_full_prefix_budget_keeps_summing_without_storing() {
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
        integral.strip([0.3, 0.4]);
        integral.memo.borrow_mut().prefix_bytes = MEMO_BYTES;
        for u in [0.5, 0.9, 0.7] {
            let fresh = quadrature(&surface.breaks(0, 0., u), 2, |s| {
                Ok(integrand(&surface, [s, 0.4], [0.; 3]))
            })
            .unwrap();
            assert_eq!(
                fresh.map(f64::to_bits),
                integral.strip([u, 0.4]).map(f64::to_bits)
            );
        }
        let memo = integral.memo.borrow();
        assert!(memo.prefixes.values().all(|p| p.ends.is_empty()));
        // A v first taken with the budget spent is never marked recurring.
        drop(memo);
        integral.strip([0.5, 0.8]);
        assert!(!integral.memo.borrow().seen.contains(&0.8_f64.to_bits()));
    }
}
