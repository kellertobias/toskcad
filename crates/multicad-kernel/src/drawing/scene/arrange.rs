//! Clipping of exact paths to the retained region and their arrangement in
//! the paper plane. Chords only locate events; every event parameter is then
//! refined on the exact path.
use super::*;
use multicad_rust_kernel::curved_projection::Ellipse;
use std::collections::{HashMap, HashSet};

/// A clipped range of a candidate path with its certified chords.
pub(in crate::drawing) struct Strand {
    pub member: usize,
    pub path: Arc<dyn Path>,
    pub ts: Vec<f64>,
    pub pts: Vec<P3>,
    /// Index of the plane limit the whole strand lies on.
    pub plane: Option<usize>,
}

/// One limit of the retained region: retained where `value <= 0`.
pub(in crate::drawing) struct Limit {
    pub value: Box<dyn Fn(P3) -> f64>,
    /// Unit normal into the retained side for a plane limit.
    pub inward: Option<P3>,
    /// Center and radius of a circular paper limit. A chord can pass through
    /// such a convex region between two outside ends.
    pub disc: Option<(P2, f64)>,
}

fn paper(p: P3) -> P2 {
    [p[0], p[1]]
}

fn sub(a: P2, b: P2) -> P2 {
    [a[0] - b[0], a[1] - b[1]]
}

fn cross(a: P2, b: P2) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}

fn dot(a: P2, b: P2) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}

pub(in crate::drawing) fn distance(a: P2, b: P2) -> f64 {
    (a[0] - b[0]).hypot(a[1] - b[1])
}

/// Whether two published ellipses are the same paper conic.
pub(in crate::drawing) fn same_conic(a: Ellipse, b: Ellipse, eps: f64) -> bool {
    distance(a.center, b.center) <= eps
        && (a.major - b.major).abs() <= eps
        && (a.minor - b.minor).abs() <= eps
        && ((a.major - a.minor).abs() <= eps || cross(a.u, b.u).abs() * a.major <= eps)
}

fn segment_distance(p: P2, a: P2, b: P2) -> f64 {
    let d = sub(b, a);
    let along = (dot(sub(p, a), d) / dot(d, d).max(1e-300)).clamp(0., 1.);
    distance(p, [a[0] + along * d[0], a[1] + along * d[1]])
}

/// Bisect the exact path for the zero of a limit between two parameters.
fn root(path: &dyn Path, value: &dyn Fn(P3) -> f64, mut a: f64, mut b: f64) -> f64 {
    let inside = value(path.point(a)) <= 0.;
    for _ in 0..80 {
        let mid = 0.5 * (a + b);
        if mid == a || mid == b {
            break;
        }
        if (value(path.point(mid)) <= 0.) == inside {
            a = mid;
        } else {
            b = mid;
        }
    }
    0.5 * (a + b)
}

/// The parts of a path inside every limit, as runs of parameters and points.
pub(in crate::drawing) fn clip(
    path: &Arc<dyn Path>,
    limits: &[&Limit],
    tolerance: f64,
    eps: f64,
) -> Vec<(Vec<f64>, Vec<P3>, Option<usize>)> {
    let ts = path.flatten(tolerance);
    let pts: Vec<P3> = ts.iter().map(|&t| path.point(t)).collect();
    let mut runs = vec![(ts, pts)];
    for limit in limits {
        let mut next = Vec::new();
        for (mut ts, mut pts) in runs {
            if let Some((center, radius)) = limit.disc {
                let mut i = 0;
                while i + 1 < ts.len() {
                    let [a, b] = [paper(pts[i]), paper(pts[i + 1])];
                    let direction = sub(b, a);
                    let length = dot(direction, direction);
                    let station = dot(sub(center, a), direction) / length;
                    if distance(a, center) > radius
                        && distance(b, center) > radius
                        && station > 0.
                        && station < 1.
                        && cross(sub(center, a), direction).abs() / length.sqrt() < radius - eps
                    {
                        let t = ts[i] + (ts[i + 1] - ts[i]) * station;
                        ts.insert(i + 1, t);
                        pts.insert(i + 1, path.point(t));
                        i += 1;
                    }
                    i += 1;
                }
            }
            // A curve can leave and re-enter a limit between two certified
            // stops. The exact path strays from its chord by at most the
            // flattening tolerance, so only chords ending that close to the
            // limit can hide such an excursion; its turning point, located on
            // the exact path, then separates the two crossings.
            let mut i = 0;
            while i + 1 < ts.len() {
                let [g0, g1] = [(limit.value)(pts[i]), (limit.value)(pts[i + 1])];
                if (g0 <= eps) == (g1 <= eps) && g0.abs().min(g1.abs()) <= tolerance + eps {
                    let value = |t: f64| (limit.value)(path.point(t));
                    let [a, b] = [ts[i], ts[i + 1]];
                    // Golden-section search for the extremum away from the ends' side.
                    let sign = if g0 <= eps { -1. } else { 1. };
                    let ratio = 0.618_033_988_749_894_9;
                    let [mut lo, mut hi] = [a, b];
                    let mut x1 = hi - ratio * (hi - lo);
                    let mut x2 = lo + ratio * (hi - lo);
                    let [mut f1, mut f2] = [sign * value(x1), sign * value(x2)];
                    for _ in 0..48 {
                        if f1 < f2 {
                            hi = x2;
                            x2 = x1;
                            f2 = f1;
                            x1 = hi - ratio * (hi - lo);
                            f1 = sign * value(x1);
                        } else {
                            lo = x1;
                            x1 = x2;
                            f1 = f2;
                            x2 = lo + ratio * (hi - lo);
                            f2 = sign * value(x2);
                        }
                    }
                    let t = 0.5 * (lo + hi);
                    let point = path.point(t);
                    if ((limit.value)(point) <= eps) != (g0 <= eps) && t > a && t < b {
                        ts.insert(i + 1, t);
                        pts.insert(i + 1, point);
                        i += 1;
                    }
                }
                i += 1;
            }
            let values: Vec<f64> = pts.iter().map(|&p| (limit.value)(p)).collect();
            let mut run: (Vec<f64>, Vec<P3>) = Default::default();
            for i in 0..ts.len() {
                let inside = values[i] <= eps;
                if inside {
                    run.0.push(ts[i]);
                    run.1.push(pts[i]);
                }
                if i + 1 == ts.len() {
                    break;
                }
                let ahead = values[i + 1] <= eps;
                if inside != ahead {
                    // A vertex within the band is itself the crossing.
                    let (near, far) = if inside { (i, i + 1) } else { (i + 1, i) };
                    if values[near] < 0. {
                        let t = root(path.as_ref(), &limit.value, ts[near], ts[far]);
                        run.0.push(t);
                        run.1.push(path.point(t));
                    }
                    if inside {
                        next.push(std::mem::take(&mut run));
                    }
                }
            }
            next.push(run);
        }
        next.retain(|run: &(Vec<f64>, Vec<P3>)| run.0.len() >= 2);
        runs = next;
    }
    runs.into_iter()
        .filter(|(_, pts)| {
            pts.windows(2)
                .map(|pair| {
                    paths::dot3(
                        paths::add(pair[1], paths::scaled(pair[0], -1.)),
                        paths::add(pair[1], paths::scaled(pair[0], -1.)),
                    )
                    .sqrt()
                })
                .sum::<f64>()
                > eps
        })
        .map(|(ts, pts)| {
            let plane = limits.iter().position(|limit| {
                limit.inward.is_some()
                    && pts.iter().all(|&p| (limit.value)(p).abs() <= eps)
                    && ts.windows(2).all(|pair| {
                        (limit.value)(path.point(0.5 * (pair[0] + pair[1]))).abs() <= eps
                    })
            });
            (ts, pts, plane)
        })
        .collect()
}

/// The parameter within chord `index` whose exact paper point has the same
/// position along the chord as `target`.
pub(in crate::drawing) fn locate(strand: &Strand, index: usize, target: P2, eps: f64) -> f64 {
    let [a, b] = [paper(strand.pts[index]), paper(strand.pts[index + 1])];
    let [mut ta, mut tb] = [strand.ts[index], strand.ts[index + 1]];
    let direction = sub(b, a);
    let length = dot(direction, direction);
    if length == 0. {
        return ta;
    }
    let offset = |t: f64| dot(sub(paper(strand.path.point(t)), target), direction) / length.sqrt();
    let guess = ta + (tb - ta) * (dot(sub(target, a), direction) / length).clamp(0., 1.);
    if offset(guess).abs() <= 1e-6 * eps {
        return guess;
    }
    let [mut fa, mut fb] = [offset(ta), offset(tb)];
    if fa * fb > 0. {
        return guess;
    }
    let mut t = guess;
    for _ in 0..100 {
        let f = offset(t);
        if f.abs() <= 1e-6 * eps {
            break;
        }
        if (f < 0.) == (fa < 0.) {
            ta = t;
            fa = f;
            fb *= 0.5;
        } else {
            tb = t;
            fb = f;
            fa *= 0.5;
        }
        let next = ta + (tb - ta) * fa / (fa - fb);
        t = if next.is_finite() && (next - ta) * (next - tb) <= 0. {
            next
        } else {
            0.5 * (ta + tb)
        };
    }
    t
}

/// Newton refinement of a transversal crossing of two exact paths.
fn crossing(a: (&Strand, usize), b: (&Strand, usize), start: [f64; 2], eps: f64) -> [f64; 2] {
    let ranges = [a, b].map(|(strand, i)| [strand.ts[i], strand.ts[i + 1]]);
    let at = |which: usize, t: f64| paper([a.0, b.0][which].path.point(t));
    let mut t = start;
    for _ in 0..8 {
        let residual = sub(at(0, t[0]), at(1, t[1]));
        if residual[0].hypot(residual[1]) <= 1e-6 * eps {
            break;
        }
        let slope = |which: usize| {
            let [lo, hi] = ranges[which];
            let h = 1e-4 * (hi - lo);
            let [p, q] = [at(which, t[which] - h), at(which, t[which] + h)];
            [(q[0] - p[0]) / (2. * h), (q[1] - p[1]) / (2. * h)]
        };
        let [da, db] = [slope(0), slope(1)];
        let det = cross(da, db);
        if det == 0. || !det.is_finite() {
            break;
        }
        // residual + da * dt - db * du = 0
        let dt = -cross(residual, db) / det;
        let du = -cross(residual, da) / det;
        let next = [t[0] + dt, t[1] + du];
        if (0..2).any(|k| !(ranges[k][0]..=ranges[k][1]).contains(&next[k])) {
            break;
        }
        t = next;
    }
    t
}

/// Parameters at which every strand must be split so that two pieces either
/// coincide in the paper or meet only at their ends. Each list is ascending
/// and includes both strand ends, with the exact point of every stop.
pub(in crate::drawing) fn arrange(strands: &[Strand], eps: f64) -> Vec<Vec<(f64, P3)>> {
    let mut chords: Vec<(usize, usize)> = Vec::new();
    let mut lo = [f64::INFINITY; 2];
    let mut hi = [f64::NEG_INFINITY; 2];
    for (s, strand) in strands.iter().enumerate() {
        for i in 0..strand.pts.len() - 1 {
            chords.push((s, i));
        }
        for p in &strand.pts {
            for k in 0..2 {
                lo[k] = lo[k].min(p[k]);
                hi[k] = hi[k].max(p[k]);
            }
        }
    }
    let mut splits: Vec<Vec<f64>> = strands
        .iter()
        .map(|s| vec![s.ts[0], *s.ts.last().unwrap()])
        .collect();
    // Strands on one paper conic coincide along whole arcs, whatever their
    // parametrisations. Their chords would cross at spurious points, so they
    // are split at each other's ends instead and compared as arcs.
    let conics: Vec<Option<Ellipse>> = strands
        .iter()
        .map(|s| match s.path.paper([s.ts[0], *s.ts.last().unwrap()]) {
            Some(Curve::Ellipse(ellipse)) => Some(ellipse),
            _ => None,
        })
        .collect();
    let mut twins: HashSet<[usize; 2]> = HashSet::new();
    for i in 0..strands.len() {
        for j in i + 1..strands.len() {
            if let (Some(a), Some(b)) = (conics[i], conics[j])
                && same_conic(a, b, eps)
            {
                twins.insert([i, j]);
                for [on, other] in [[i, j], [j, i]] {
                    for end in [strands[other].pts[0], *strands[other].pts.last().unwrap()] {
                        let target = paper(end);
                        let nearest = (0..strands[on].pts.len() - 1).min_by(|&x, &y| {
                            let reach = |k: usize| {
                                segment_distance(
                                    target,
                                    paper(strands[on].pts[k]),
                                    paper(strands[on].pts[k + 1]),
                                )
                            };
                            reach(x).total_cmp(&reach(y))
                        });
                        if let Some(chord) = nearest {
                            let t = locate(&strands[on], chord, target, eps);
                            if distance(paper(strands[on].path.point(t)), target) <= eps {
                                splits[on].push(t);
                            }
                        }
                    }
                }
            }
        }
    }
    if !chords.is_empty() {
        let side = ((chords.len() as f64).sqrt().ceil() as usize).clamp(1, 512);
        let size = [0, 1].map(|k| ((hi[k] - lo[k]) / side as f64).max(eps));
        let cell = |p: P2| [0, 1].map(|k| (((p[k] - lo[k]) / size[k]) as usize).min(side - 1));
        let ends =
            |&(s, i): &(usize, usize)| [paper(strands[s].pts[i]), paper(strands[s].pts[i + 1])];
        let mut grid: HashMap<[usize; 2], Vec<usize>> = HashMap::new();
        let boxes: Vec<[[usize; 2]; 2]> = chords
            .iter()
            .map(|chord| {
                let [a, b] = ends(chord);
                let low =
                    cell([a[0].min(b[0]) - eps, a[1].min(b[1]) - eps].map(|v| v.max(f64::MIN)));
                let high = cell([a[0].max(b[0]) + eps, a[1].max(b[1]) + eps]);
                [[low[0].min(high[0]), low[1].min(high[1])], high]
            })
            .collect();
        for (index, [low, high]) in boxes.iter().enumerate() {
            for x in low[0]..=high[0] {
                for y in low[1]..=high[1] {
                    grid.entry([x, y]).or_default().push(index);
                }
            }
        }
        let mut cells: Vec<_> = grid.into_iter().collect();
        cells.sort_by_key(|(key, _)| *key);
        for (key, members) in cells {
            for (n, &first) in members.iter().enumerate() {
                for &second in &members[n + 1..] {
                    let (ca, cb) = (chords[first], chords[second]);
                    if ca.0 == cb.0 || twins.contains(&[ca.0.min(cb.0), ca.0.max(cb.0)]) {
                        continue;
                    }
                    // Test each pair once, in the first cell both cover.
                    let shared = [0, 1].map(|k| boxes[first][0][k].max(boxes[second][0][k]));
                    if shared != key {
                        continue;
                    }
                    meet(strands, ca, cb, eps, &mut splits);
                }
            }
        }
    }
    strands
        .iter()
        .zip(splits)
        .map(|(strand, mut stops)| {
            stops.sort_by(f64::total_cmp);
            let mut output: Vec<(f64, P3)> = Vec::with_capacity(stops.len());
            for t in stops {
                let point = strand.path.point(t);
                // Stops closer than the tolerance are one event.
                if output.last().is_none_or(|(_, previous)| {
                    distance(paper(*previous), paper(point)).max((point[2] - previous[2]).abs())
                        > eps
                }) {
                    output.push((t, point));
                }
            }
            let end = (*strand.ts.last().unwrap(), *strand.pts.last().unwrap());
            if output.len() > 1 {
                output.pop();
            }
            output.push(end);
            output
        })
        .collect()
}

fn meet(
    strands: &[Strand],
    ca: (usize, usize),
    cb: (usize, usize),
    eps: f64,
    splits: &mut [Vec<f64>],
) {
    let ends = |(s, i): (usize, usize)| [paper(strands[s].pts[i]), paper(strands[s].pts[i + 1])];
    let [p, q] = [ends(ca), ends(cb)];
    for k in 0..2 {
        if p[0][k].max(p[1][k]) + eps < q[0][k].min(q[1][k])
            || q[0][k].max(q[1][k]) + eps < p[0][k].min(p[1][k])
        {
            return;
        }
    }
    let [dp, dq] = [sub(p[1], p[0]), sub(q[1], q[0])];
    let [lp, lq] = [dp[0].hypot(dp[1]), dq[0].hypot(dq[1])];
    if lp <= eps || lq <= eps {
        return;
    }
    // Signed distances of each chord's ends from the other chord's line.
    let from_p = q.map(|point| cross(dp, sub(point, p[0])) / lp);
    let from_q = p.map(|point| cross(dq, sub(point, q[0])) / lq);
    let mut split = |chord: (usize, usize), target: P2| {
        let t = locate(&strands[chord.0], chord.1, target, eps);
        splits[chord.0].push(t);
    };
    // An end of one chord strictly inside the other chord.
    let mut touched = false;
    for (chord, own, other, offsets, direction, length) in
        [(ca, p, q, from_p, dp, lp), (cb, q, p, from_q, dq, lq)]
    {
        for (end, offset) in other.into_iter().zip(offsets) {
            let station = dot(sub(end, own[0]), direction) / length;
            if offset.abs() <= eps && station > eps && station < length - eps {
                split(chord, end);
                touched = true;
            }
        }
    }
    if touched || from_p.iter().all(|d| d.abs() <= eps) || from_q.iter().all(|d| d.abs() <= eps) {
        return;
    }
    if from_p[0] * from_p[1] < 0. && from_q[0] * from_q[1] < 0. {
        let s = from_q[0] / (from_q[0] - from_q[1]);
        let u = from_p[0] / (from_p[0] - from_p[1]);
        let (sa, sb) = (&strands[ca.0], &strands[cb.0]);
        let start = [
            sa.ts[ca.1] + (sa.ts[ca.1 + 1] - sa.ts[ca.1]) * s,
            sb.ts[cb.1] + (sb.ts[cb.1 + 1] - sb.ts[cb.1]) * u,
        ];
        let [ta, tb] = crossing((sa, ca.1), (sb, cb.1), start, eps);
        splits[ca.0].push(ta);
        splits[cb.0].push(tb);
    }
}
