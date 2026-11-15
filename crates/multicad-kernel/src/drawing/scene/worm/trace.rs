//! Zero sets of a function over a rectangular patch domain as chains of
//! cell segments. Every segment is a graph over one parameter inside one
//! grid cell with a sign change across the cell, so any of its points is
//! recovered from the exact function by bracketed root finding.
use super::super::*;

/// The root of a function with a sign change between `a` and `b`, by
/// Illinois regula falsi safeguarded with bisection, to the last few bits.
pub(super) fn root(f: impl Fn(f64) -> f64, mut a: f64, mut b: f64) -> f64 {
    let (mut fa, mut fb) = (f(a), f(b));
    if fa == 0. {
        return a;
    }
    if fb == 0. {
        return b;
    }
    let scale = a.abs().max(b.abs()).max(1e-300);
    let mut side = 0;
    for step in 0..200 {
        if (b - a).abs() <= 2. * f64::EPSILON * scale {
            break;
        }
        let mut x = if step % 3 == 2 {
            0.5 * (a + b)
        } else {
            (a * fb - b * fa) / (fb - fa)
        };
        if !(x > a.min(b) && x < a.max(b)) {
            x = 0.5 * (a + b);
            if x == a || x == b {
                break;
            }
        }
        let fx = f(x);
        if fx == 0. {
            return x;
        }
        if (fx < 0.) == (fa < 0.) {
            a = x;
            fa = fx;
            if side == 1 {
                fb *= 0.5;
            }
            side = 1;
        } else {
            b = x;
            fb = fx;
            if side == 2 {
                fa *= 0.5;
            }
            side = 2;
        }
    }
    if fa.abs() < fb.abs() { a } else { b }
}

/// One piece of a zero curve inside the grid cell `cell = [[u0, u1], [z0, z1]]`
/// from `a` to `b`, a graph over `u` when `along_u` and over `z` otherwise.
#[derive(Clone, Copy, Debug)]
pub(super) struct Seg {
    pub a: [f64; 2],
    pub b: [f64; 2],
    pub cell: [[f64; 2]; 2],
    /// The whole patch domain.
    pub domain: [[f64; 2]; 2],
    pub along_u: bool,
}

impl Seg {
    /// The parameter point at fraction `s` of the segment.
    pub(super) fn at(&self, f: &dyn Fn(f64, f64) -> f64, s: f64) -> [f64; 2] {
        self.solve(f, s).unwrap_or_else(|| {
            // The curve leaves the cell through a corner within rounding.
            std::array::from_fn(|k| self.a[k] + (self.b[k] - self.a[k]) * s)
        })
    }

    /// The curve point over the free parameter at fraction `s`: the root in
    /// the cell, or within one cell beyond it where the curve bulges across
    /// a grid line between two nodes.
    fn solve(&self, f: &dyn Fn(f64, f64) -> f64, s: f64) -> Option<[f64; 2]> {
        if s <= 0. {
            return Some(self.a);
        }
        if s >= 1. {
            return Some(self.b);
        }
        let k = usize::from(!self.along_u);
        let free = self.a[k] + (self.b[k] - self.a[k]) * s;
        let [lo, hi] = self.cell[1 - k];
        let value = |other: f64| {
            if self.along_u {
                f(free, other)
            } else {
                f(other, free)
            }
        };
        let reach = hi - lo;
        let [low, high] = self.domain[1 - k];
        let other = [
            [lo, hi],
            [(lo - reach).max(low), lo],
            [hi, (hi + reach).min(high)],
        ]
        .into_iter()
        .find(|[a, b]| b > a && (value(*a) < 0.) != (value(*b) < 0.))
        .map(|[a, b]| root(value, a, b))?;
        Some(if self.along_u {
            [free, other]
        } else {
            [other, free]
        })
    }
}

/// Chains of the zero set of `f` over `[0, 1] × z`, where `f` already
/// carries the bias that makes a vanishing value positive.
pub(super) fn trace(
    f: &dyn Fn(f64, f64) -> f64,
    z: [f64; 2],
    grid: [usize; 2],
) -> Result<Vec<Vec<Seg>>, GeometryError> {
    let mut grid = grid;
    for _ in 0..4 {
        if let Some(chains) = attempt(f, z, grid) {
            return Ok(chains);
        }
        grid = grid.map(|n| 2 * n);
    }
    Err(err(
        GeometryErrorKind::Limit,
        "Rust worm drawing could not resolve a surface curve within its grid budget",
    ))
}

fn attempt(
    f: &dyn Fn(f64, f64) -> f64,
    z: [f64; 2],
    [nu, nz]: [usize; 2],
) -> Option<Vec<Vec<Seg>>> {
    let us: Vec<f64> = (0..=nu).map(|i| i as f64 / nu as f64).collect();
    let zs: Vec<f64> = (0..=nz)
        .map(|j| {
            if j == nz {
                z[1]
            } else {
                z[0] + (z[1] - z[0]) * j as f64 / nz as f64
            }
        })
        .collect();
    let sign: Vec<Vec<bool>> = us
        .iter()
        .map(|&u| zs.iter().map(|&z| f(u, z) < 0.).collect())
        .collect();
    let any = |wanted: bool| sign.iter().flatten().any(|&s| s == wanted);
    if !any(true) || !any(false) {
        return Some(Vec::new());
    }
    // Edge ids: horizontal edge (i, j) joins nodes (i, j) and (i + 1, j);
    // vertical edge (i, j) joins nodes (i, j) and (i, j + 1).
    let horizontal = |i: usize, j: usize| 2 * (j * (nu + 1) + i);
    let vertical = |i: usize, j: usize| 2 * (j * (nu + 1) + i) + 1;
    let mut crossings: std::collections::HashMap<usize, [f64; 2]> = Default::default();
    let mut cross = |id: usize, along_u: bool, i: usize, j: usize| -> [f64; 2] {
        *crossings.entry(id).or_insert_with(|| {
            if along_u {
                [root(|u| f(u, zs[j]), us[i], us[i + 1]), zs[j]]
            } else {
                [us[i], root(|z| f(us[i], z), zs[j], zs[j + 1])]
            }
        })
    };
    // Segments with the ids of the edges they join.
    let mut segs: Vec<(Seg, [usize; 2])> = Vec::new();
    for i in 0..nu {
        for j in 0..nz {
            let [bl, br, tl, tr] = [
                sign[i][j],
                sign[i + 1][j],
                sign[i][j + 1],
                sign[i + 1][j + 1],
            ];
            let mut found: Vec<(usize, [f64; 2])> = Vec::new();
            // Bottom, right, top, left.
            if bl != br {
                let id = horizontal(i, j);
                found.push((id, cross(id, true, i, j)));
            }
            if br != tr {
                let id = vertical(i + 1, j);
                found.push((id, cross(id, false, i + 1, j)));
            }
            if tl != tr {
                let id = horizontal(i, j + 1);
                found.push((id, cross(id, true, i, j + 1)));
            }
            if bl != tl {
                let id = vertical(i, j);
                found.push((id, cross(id, false, i, j)));
            }
            let cell = [[us[i], us[i + 1]], [zs[j], zs[j + 1]]];
            let pairs: Vec<[usize; 2]> = match found.len() {
                0 => vec![],
                2 => vec![[0, 1]],
                _ => {
                    // A saddle: the centre decides which corners connect.
                    let centre = f(0.5 * (us[i] + us[i + 1]), 0.5 * (zs[j] + zs[j + 1])) < 0.;
                    if centre == bl {
                        vec![[0, 1], [2, 3]]
                    } else {
                        vec![[3, 0], [1, 2]]
                    }
                }
            };
            for [p, q] in pairs {
                let (a, b) = (found[p].1, found[q].1);
                let extent = [
                    (b[0] - a[0]).abs() / (cell[0][1] - cell[0][0]),
                    (b[1] - a[1]).abs() / (cell[1][1] - cell[1][0]),
                ];
                let seg = Seg {
                    a,
                    b,
                    cell,
                    domain: [[0., 1.], z],
                    along_u: extent[0] >= extent[1],
                };
                // The curve must be bracketed across the whole segment.
                let k = usize::from(!seg.along_u);
                if extent[k] > 1e-9 && (1..16).any(|i| seg.solve(f, f64::from(i) / 16.).is_none()) {
                    return None;
                }
                segs.push((seg, [found[p].0, found[q].0]));
            }
        }
    }
    // Link segments that share an edge crossing into chains.
    let mut at: std::collections::HashMap<usize, Vec<usize>> = Default::default();
    for (index, (_, ids)) in segs.iter().enumerate() {
        for id in ids {
            at.entry(*id).or_default().push(index);
        }
    }
    let mut used = vec![false; segs.len()];
    let mut chains = Vec::new();
    let walk = |start: usize, entry: usize, used: &mut Vec<bool>| {
        let mut chain = Vec::new();
        let (mut index, mut entry) = (start, entry);
        loop {
            used[index] = true;
            let (mut seg, ids) = segs[index];
            let exit = if ids[0] == entry {
                ids[1]
            } else {
                std::mem::swap(&mut seg.a, &mut seg.b);
                ids[0]
            };
            chain.push(seg);
            match at[&exit].iter().find(|&&next| !used[next]) {
                Some(&next) => {
                    index = next;
                    entry = exit;
                }
                None => return chain,
            }
        }
    };
    // Open chains start on an edge that only one segment reaches.
    for index in 0..segs.len() {
        for end in 0..2 {
            if !used[index] && at[&segs[index].1[end]].len() == 1 {
                chains.push(walk(index, segs[index].1[end], &mut used));
            }
        }
    }
    for index in 0..segs.len() {
        if !used[index] {
            chains.push(walk(index, segs[index].1[0], &mut used));
        }
    }
    Some(chains)
}
