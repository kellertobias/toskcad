//! Exact plane sections of a bevel gear.
//!
//! Each face is a strip of rulings between two boundary curves: the side
//! faces run from the toe rim to the heel rim along apex rays, a cap from
//! its bottom (the bore circle or the back-cone apex) to its rim along
//! back-cone generators, and the bore from its toe circle to its heel
//! circle. A ruling meets the plane where its two ends lie on opposite
//! sides, so the cut of a strip changes only where the plane crosses one of
//! the two boundary curves. Those crossings are the vertices of the section;
//! every vertex ends exactly one cut on each of the two strips that share
//! its edge.
//!
//! A cut is evaluated in one of two charts of its strip: by the ruling
//! parameter, with the distance along the ruling from the plane equation, or
//! by that distance, with the ruling solved from the plane equation. The
//! second chart keeps full precision where the plane runs along the rulings,
//! for instance when it contains the gear axis.
use super::roots::{self, positive};
use super::*;
use std::sync::Mutex;

#[derive(Clone, Copy, Debug)]
enum Sheet {
    Side { sector: u32, piece: usize },
    Cap(End),
    Bore,
}

#[derive(Clone, Copy, Debug)]
enum Chart {
    /// Ruling parameters; the distance follows from the plane.
    Sweep([f64; 2]),
    /// Distances along the ruling; the ruling is solved within a bracket in
    /// which the plane equation is monotone.
    Along { q: [f64; 2], s: [f64; 2] },
}

/// The plane `normal · p = level` in the gear frame.
#[derive(Clone, Copy, Debug)]
struct Plane {
    normal: P3,
    level: f64,
}

impl Core {
    /// Base point and direction of a ruling with their derivatives.
    fn ruling(&self, sheet: Sheet, s: f64) -> [P3; 4] {
        match sheet {
            Sheet::Side { sector, piece } => {
                let (u, du) = self.direction(sector, piece, s);
                [[0.; 3], u, [0.; 3], du]
            }
            Sheet::Cap(end) => {
                let (sin, cos) = s.sin_cos();
                let [sd, cd] = self.pitch;
                [
                    [0., 0., self.bevel.back_cone_apex_mm(end)],
                    [cd * cos, cd * sin, -sd],
                    [0.; 3],
                    [-cd * sin, cd * cos, 0.],
                ]
            }
            Sheet::Bore => {
                let bore = self.bore.unwrap_or(0.);
                let (sin, cos) = s.sin_cos();
                [
                    [bore * cos, bore * sin, 0.],
                    [0., 0., 1.],
                    [-bore * sin, bore * cos, 0.],
                    [0.; 3],
                ]
            }
        }
    }

    /// `(α, β, α', β')` of the plane equation `q α(s) = β(s)` on a strip.
    fn terms(&self, sheet: Sheet, plane: Plane, s: f64) -> [f64; 4] {
        let [base, direction, dbase, ddirection] = self.ruling(sheet, s);
        [
            dot3(plane.normal, direction),
            plane.level - dot3(plane.normal, base),
            dot3(plane.normal, ddirection),
            -dot3(plane.normal, dbase),
        ]
    }

    /// The plane equation `q α(s) = β(s)` of a strip has a constant `β`
    /// (rulings through one apex) or a constant `α` (parallel rulings), so
    /// for a fixed distance it is monotone in `s` where the other term is.
    /// Returns that term's derivative at `s` and a bound of its second
    /// derivative over the strip.
    fn slope(&self, sheet: Sheet, plane: Plane, s: f64) -> [f64; 2] {
        let [_, _, dalpha, dbeta] = self.terms(sheet, plane, s);
        let across = plane.normal[0].hypot(plane.normal[1]);
        match sheet {
            Sheet::Side { piece, .. } => [dalpha, self.shapes[piece].curvature],
            Sheet::Cap(_) => [dalpha, self.pitch[1] * across],
            Sheet::Bore => [dbeta, self.bore.unwrap_or(0.) * across],
        }
    }
}

struct Cut {
    core: Arc<Core>,
    sheet: Sheet,
    plane: Plane,
    chart: Chart,
    flip: bool,
    chords: Mutex<Option<(u64, Vec<f64>)>>,
}

impl Cut {
    fn local(&self, t: f64) -> P3 {
        let t = if self.flip { 1. - t } else { t };
        let (s, q) = match self.chart {
            Chart::Sweep([a, b]) => {
                let s = a + (b - a) * t;
                let [alpha, beta, ..] = self.core.terms(self.sheet, self.plane, s);
                (s, beta / alpha)
            }
            Chart::Along { q: [a, b], s } => {
                let q = a + (b - a) * t;
                (self.solve(q, s), q)
            }
        };
        let [base, direction, ..] = self.core.ruling(self.sheet, s);
        add(base, scaled(direction, q))
    }

    /// The ruling within the bracket whose point at distance `q` lies in
    /// the plane: a safeguarded Newton iteration on a monotone function.
    fn solve(&self, q: f64, [mut lo, mut hi]: [f64; 2]) -> f64 {
        if hi <= lo {
            return lo;
        }
        let eval = |s: f64| {
            let [alpha, beta, dalpha, dbeta] = self.core.terms(self.sheet, self.plane, s);
            (q * alpha - beta, q * dalpha - dbeta)
        };
        let from = positive(eval(lo).0);
        if positive(eval(hi).0) == from {
            // The requested distance lies at an end of the bracket.
            let [a, b] = [eval(lo).0.abs(), eval(hi).0.abs()];
            return if a <= b { lo } else { hi };
        }
        let mut s = 0.5 * (lo + hi);
        for _ in 0..200 {
            let (value, rate) = eval(s);
            if positive(value) == from {
                lo = s;
            } else {
                hi = s;
            }
            let newton = s - value / rate;
            let next = if newton > lo && newton < hi {
                newton
            } else {
                0.5 * (lo + hi)
            };
            if next == s || !(next > lo && next < hi) {
                break;
            }
            s = next;
        }
        s
    }
}

impl Path for Cut {
    fn point(&self, t: f64) -> P3 {
        self.core.to_view(self.local(t))
    }

    fn flatten(&self, tolerance_mm: f64) -> Vec<f64> {
        let mut cache = self
            .chords
            .lock()
            .unwrap_or_else(|poison| poison.into_inner());
        if let Some((key, stops)) = cache.as_ref()
            && *key == tolerance_mm.to_bits()
        {
            return stops.clone();
        }
        let mut stops = vec![0.];
        self.refine([0., 1.], tolerance_mm, 0, &mut stops);
        *cache = Some((tolerance_mm.to_bits(), stops.clone()));
        stops
    }

    fn paper(&self, range: [f64; 2]) -> Option<Curve> {
        let [p, q] = range.map(|t| self.point(t));
        let line = Some(Curve::Line([[p[0], p[1]], [q[0], q[1]]]));
        if matches!(self.chart, Chart::Along { s: [a, b], .. } if a == b) {
            return line;
        }
        let core = &self.core;
        // The quadric (x − a)ᵀ M (x − a) = k of the face in view space.
        let (apex, axis, cosine2, radius2) = match self.sheet {
            Sheet::Bore => ([0.; 3], [0., 0., 1.], None, core.bore?.powi(2)),
            Sheet::Cap(end) => (
                [0., 0., core.bevel.back_cone_apex_mm(end)],
                [0., 0., 1.],
                Some(core.pitch[0].powi(2)),
                0.,
            ),
            Sheet::Side { sector, piece } => match core.pieces[piece] {
                Piece::RootLeft | Piece::RootRight => {
                    ([0.; 3], [0., 0., 1.], Some(core.root.cos().powi(2)), 0.)
                }
                Piece::Tip => ([0.; 3], [0., 0., 1.], Some(core.tip.cos().powi(2)), 0.),
                Piece::MeridianLeft | Piece::MeridianRight => return line,
                Piece::FilletLeft | Piece::FilletRight => {
                    let shape = core.shapes[piece];
                    (
                        [0.; 3],
                        core.turn(sector, shape.axis),
                        Some(shape.cosine.powi(2)),
                        0.,
                    )
                }
                Piece::InvoluteLeft | Piece::InvoluteRight => return None,
            },
        };
        let apex = core.to_view(apex);
        let axis = core.vector(axis);
        let form = |i: usize, j: usize| {
            let unit = f64::from(u8::from(i == j));
            match cosine2 {
                Some(cosine2) => axis[i] * axis[j] - cosine2 * unit,
                None => unit - axis[i] * axis[j],
            }
        };
        let [a, b, c] = [form(0, 0), form(0, 1), form(1, 1)];
        let m = [form(0, 2), form(1, 2)];
        let dz = p[2] - apex[2];
        let det = a * c - b * b;
        // Both eigenvalues of the paper form must share the sign of the
        // constant for an ellipse.
        let mean = 0.5 * (a + c);
        let spread = (0.25 * (a - c) * (a - c) + b * b).sqrt();
        let values = [mean + spread, mean - spread];
        if det.abs() <= 1e-12 || values[0] * values[1] <= 0. {
            return None;
        }
        let solved = [(c * m[0] - b * m[1]) / det, (a * m[1] - b * m[0]) / det];
        let center = [apex[0] - dz * solved[0], apex[1] - dz * solved[1]];
        let constant = radius2 - dz * dz * (form(2, 2) - m[0] * solved[0] - m[1] * solved[1]);
        let turn = 0.5 * (2. * b).atan2(a - c);
        let (sin, cos) = turn.sin_cos();
        // `values[0]` belongs to (cos, sin), `values[1]` to (−sin, cos).
        let radii = [constant / values[0], constant / values[1]];
        if !(radii[0] > 0. && radii[1] > 0.) {
            return None;
        }
        let radii = radii.map(f64::sqrt);
        let directions = [[cos, sin], [-sin, cos]];
        let angle = |point: P3| {
            let offset = [point[0] - center[0], point[1] - center[1]];
            let [x, y] = [0, 1]
                .map(|k| (offset[0] * directions[k][0] + offset[1] * directions[k][1]) / radii[k]);
            (y.atan2(x), x.hypot(y))
        };
        let start = angle(p);
        let mut sweep = 0.;
        let mut previous = start.0;
        let scale = radii[0].min(radii[1]);
        if scale.is_nan() || scale <= 1e-6 * core.scale {
            return None;
        }
        for step in 1..=16 {
            let t = range[0] + (range[1] - range[0]) * f64::from(step) / 16.;
            let (next, radius) = angle(self.point(t));
            // Fall back to chords unless the cut lies on this ellipse.
            if radius.is_nan() || (radius - 1.).abs() * scale > 1e-9 * core.scale {
                return None;
            }
            sweep += (next - previous + PI).rem_euclid(TAU) - PI;
            previous = next;
        }
        if sweep.abs() < 1e-12 {
            return line;
        }
        super::super::paths::conic(
            center,
            directions[0].map(|v| v * radii[0]),
            directions[1].map(|v| v * radii[1]),
            start.0,
            sweep,
        )
    }
}

/// Bounds of the derivatives of a strip's plane terms and rulings.
struct Rates {
    /// |α'|, |α''|, |β'|, |β''|.
    alpha: [f64; 2],
    beta: [f64; 2],
    /// |dir'|, |dir''| and |base'|, |base''|.
    direction: [f64; 2],
    base: [f64; 2],
}

impl Core {
    fn rates(&self, sheet: Sheet, plane: Plane) -> Rates {
        let across = plane.normal[0].hypot(plane.normal[1]);
        match sheet {
            Sheet::Side { piece, .. } => {
                let Shape {
                    speed, curvature, ..
                } = self.shapes[piece];
                Rates {
                    alpha: [speed, curvature],
                    beta: [0.; 2],
                    direction: [speed, curvature],
                    base: [0.; 2],
                }
            }
            Sheet::Cap(_) => Rates {
                alpha: [self.pitch[1] * across; 2],
                beta: [0.; 2],
                direction: [self.pitch[1]; 2],
                base: [0.; 2],
            },
            Sheet::Bore => {
                let bore = self.bore.unwrap_or(0.);
                Rates {
                    alpha: [0.; 2],
                    beta: [bore * across; 2],
                    direction: [0.; 2],
                    base: [bore; 2],
                }
            }
        }
    }
}

impl Cut {
    /// Certified bound of the distance between the cut over a range of its
    /// chart and the chord of that range: an eighth of the squared width
    /// times a bound of the second derivative, or half a bound of the arc
    /// length, whichever is smaller.
    fn sag(&self, [t0, t1]: [f64; 2]) -> f64 {
        let core = &self.core;
        let rates = core.rates(self.sheet, self.plane);
        let [a1, a2] = rates.alpha;
        let [b1, b2] = rates.beta;
        let [d1, d2] = rates.direction;
        let [e1, e2] = rates.base;
        match self.chart {
            Chart::Sweep([a, b]) => {
                // p = base + q dir with q = β / α:
                // q' = (β' − q α') / α, q'' = (β'' − 2 q' α' − q α'') / α,
                // p'' = base'' + q'' dir + 2 q' dir' + q dir''.
                let [lo, hi] = [t0, t1].map(|t| a + (b - a) * t);
                let width = hi - lo;
                let [alpha, beta, ..] = core.terms(self.sheet, self.plane, lo + width / 2.);
                let least = alpha.abs() - a1 * width / 2.;
                if least.is_nan() || least <= 0. {
                    return f64::INFINITY;
                }
                let q = (beta.abs() + b1 * width / 2.) / least;
                let q1 = (b1 + q * a1) / least;
                let q2 = (b2 + 2. * q1 * a1 + q * a2) / least;
                let first = e1 + q1 + q * d1;
                let second = e2 + q2 + 2. * q1 * d1 + q * d2;
                (width * width * second / 8.).min(width * first / 2.)
            }
            Chart::Along { q: [a, b], s } => {
                // p = base(s) + q dir(s) with s(q) from g = q α' − β':
                // s' = −α / g,
                // s'' = −(α' s' g − α (α' + (q α'' − β'') s')) / g²,
                // p'' = (base'' + q dir'') s'² + 2 dir' s'
                //       + (base' + q dir') s''.
                let ends = [t0, t1].map(|t| a + (b - a) * t);
                let width = (ends[1] - ends[0]).abs();
                let most = ends[0].abs().max(ends[1].abs());
                let least = if ends[0] * ends[1] > 0. {
                    ends[0].abs().min(ends[1].abs())
                } else {
                    0.
                };
                let across = s[1] - s[0];
                let extent = width + (most * d1 + e1) * across;
                let [alpha, beta, dalpha, dbeta] =
                    core.terms(self.sheet, self.plane, s[0] + across / 2.);
                let rate = dalpha.abs() + a2 * across / 2.;
                let (slope, reach) = match self.sheet {
                    Sheet::Bore => (dbeta.abs() - b2 * across / 2., alpha.abs()),
                    // On the cut α = β / q with a constant β.
                    _ => (
                        least * (dalpha.abs() - a2 * across / 2.),
                        (alpha.abs() + rate * across / 2.).min(if beta == 0. {
                            0.
                        } else {
                            beta.abs() / least
                        }),
                    ),
                };
                if slope.is_nan() || slope <= 0. || !reach.is_finite() {
                    return extent;
                }
                let s1 = reach / slope;
                let most_slope = most * rate + dbeta.abs() + b2 * across / 2.;
                let s2 = (rate * s1 * most_slope + reach * (rate + (most * a2 + b2) * s1))
                    / (slope * slope);
                let second = (e2 + most * d2) * s1 * s1 + 2. * d1 * s1 + (e1 + most * d1) * s2;
                (width * width * second / 8.).min(extent)
            }
        }
    }

    /// Append the stops of `(a, b]`, halving until each chord is certified.
    fn refine(&self, [a, b]: [f64; 2], tolerance: f64, depth: u32, stops: &mut Vec<f64>) {
        let chart = if self.flip { [1. - b, 1. - a] } else { [a, b] };
        // The stop budget only ends a chart that no bound certifies.
        if depth >= 40 || stops.len() >= 1 << 16 || self.sag(chart) <= tolerance {
            stops.push(b);
            return;
        }
        let mid = 0.5 * (a + b);
        self.refine([a, mid], tolerance, depth + 1, stops);
        self.refine([mid, b], tolerance, depth + 1, stops);
    }
}

/// A crossing of the plane with a boundary curve of a strip.
#[derive(Clone, Copy, Debug)]
struct Stop {
    /// Cell and parameter on the rim loop, or zero and the azimuth.
    cell: usize,
    s: f64,
    /// Azimuth on the sweep from the start of sector 0.
    phi: f64,
    point: P3,
}

enum Arcs {
    Nothing,
    Whole,
    /// Start and end stop of each cut; the flag marks a cut through the
    /// start of the loop.
    Spans(Vec<(usize, usize, bool)>),
}

/// Cuts of a strip from the sorted crossings of its two boundary curves and
/// whether the plane separates the ruling at the start of the loop.
fn sweep(count: usize, initial: bool) -> Result<Arcs, GeometryError> {
    if count == 0 {
        return Ok(if initial { Arcs::Whole } else { Arcs::Nothing });
    }
    if count % 2 == 1 {
        return Err(err(
            GeometryErrorKind::InvalidResult,
            "Bevel section crossings do not pair",
        ));
    }
    let mut state = initial;
    let mut start = None;
    let mut first_end = None;
    let mut spans = Vec::new();
    for index in 0..count {
        state = !state;
        if state {
            start = Some(index);
        } else {
            match start.take() {
                Some(start) => spans.push((start, index, false)),
                None => first_end = Some(index),
            }
        }
    }
    if let (Some(start), Some(end)) = (start, first_end) {
        spans.push((start, end, true));
    }
    Ok(Arcs::Spans(spans))
}

struct Builder<'a> {
    core: &'a Arc<Core>,
    plane: Plane,
    cuts: Vec<Cut>,
}

impl Builder<'_> {
    fn value(&self, point: P3) -> f64 {
        dot3(self.plane.normal, point) - self.plane.level
    }

    fn push(&mut self, sheet: Sheet, chart: Chart) {
        self.cuts.push(Cut {
            core: self.core.clone(),
            sheet,
            plane: self.plane,
            chart,
            flip: false,
            chords: Mutex::new(None),
        });
    }

    /// Distance along the ruling `s` of a strip to the plane.
    fn distance(&self, sheet: Sheet, s: f64) -> f64 {
        let [alpha, beta, ..] = self.core.terms(sheet, self.plane, s);
        beta / alpha
    }

    /// Distance along its ruling of a point of a strip.
    fn along(&self, sheet: Sheet, s: f64, point: P3) -> f64 {
        let [base, direction, ..] = self.core.ruling(sheet, s);
        dot3(sub(point, base), direction)
    }

    /// One cut over `[sa, sb]` with the given end distances, split into
    /// charts that each evaluate it stably.
    fn span(&mut self, sheet: Sheet, [sa, sb]: [f64; 2], [qa, qb]: [f64; 2], depth: u32) {
        let core = self.core;
        let width = sb - sa;
        let mid = sa + width / 2.;
        // |dq/ds| = |q α' − β'| / |α| beyond this loses the point to the
        // rounding of the ruling parameter.
        let limit = 1e2 * core.scale;
        let steep = |s: f64, q: f64| {
            let [alpha, _, dalpha, dbeta] = core.terms(sheet, self.plane, s);
            !q.is_finite() || (q * dalpha - dbeta).abs() > limit * alpha.abs()
        };
        let inner = [0.25, 0.5, 0.75].map(|f| sa + width * f);
        let calm = !steep(sa, qa)
            && !steep(sb, qb)
            && inner.iter().all(|&s| !steep(s, self.distance(sheet, s)));
        if calm {
            return self.push(sheet, Chart::Sweep([sa, sb]));
        }
        let along = Chart::Along {
            q: [qa, qb],
            s: [sa, sb],
        };
        let [rate, bend] = core.slope(sheet, self.plane, mid);
        if rate.abs() > bend * width / 2. {
            return self.push(sheet, along);
        }
        let qm = self.distance(sheet, mid);
        if depth >= 48 || !(mid > sa && mid < sb) || !qm.is_finite() {
            let chart = if steep(sa, qa) || steep(sb, qb) {
                along
            } else {
                Chart::Sweep([sa, sb])
            };
            return self.push(sheet, chart);
        }
        self.span(sheet, [sa, mid], [qa, qm], depth + 1);
        self.span(sheet, [mid, sb], [qm, qb], depth + 1);
    }

    /// Crossings of the plane with one rim, in loop order, and the value at
    /// the start of the loop.
    fn rim_stops(&self, end: End) -> (Vec<Stop>, f64) {
        let core = self.core;
        let count = core.pieces.len();
        let mut stops = Vec::new();
        let eval = |sector: u32, piece: usize, s: f64| {
            let (point, rate) = core.rim(end, sector, piece, s);
            (self.value(point), dot3(self.plane.normal, rate))
        };
        for sector in 0..core.teeth {
            for piece in 0..count {
                let mut found = Vec::new();
                roots::crossings(
                    &|s| eval(sector, piece, s),
                    core.rim_bound(end, piece),
                    [0., 1.],
                    [eval(sector, piece, 0.).0, eval(sector, piece, 1.).0],
                    &mut found,
                );
                for s in found {
                    let point = core.rim(end, sector, piece, s).0;
                    stops.push(Stop {
                        cell: sector as usize * count + piece,
                        s,
                        phi: core.spacing() * f64::from(sector)
                            + core.bevel.sample(core.pieces[piece], s).phi,
                        point,
                    });
                }
            }
        }
        (stops, eval(0, 0, 0.).0)
    }

    /// Crossings of the plane with one bore circle from the sweep start.
    fn bore_stops(&self, end: End, bore: f64, start: f64) -> (Vec<Stop>, f64) {
        let core = self.core;
        let eval = |phi: f64| {
            let [x, y, _] = self.plane.normal;
            let (sin, cos) = phi.sin_cos();
            (
                self.value(core.bottom(end, phi)),
                bore * (y * cos - x * sin),
            )
        };
        let first = eval(start).0;
        let mut found = Vec::new();
        roots::crossings(
            &eval,
            bore * self.plane.normal[0].hypot(self.plane.normal[1]),
            [start, start + TAU],
            [first, first],
            &mut found,
        );
        let stops = found
            .into_iter()
            .map(|phi| Stop {
                cell: 0,
                s: phi,
                phi,
                point: core.bottom(end, phi),
            })
            .collect();
        (stops, first)
    }

    fn side(&mut self, rims: &[(Vec<Stop>, f64); 2]) -> Result<(), GeometryError> {
        let core = self.core;
        let count = core.pieces.len();
        let total = core.teeth as usize * count;
        let mut stops: Vec<Stop> = rims.iter().flat_map(|(stops, _)| stops.clone()).collect();
        stops.sort_by(|a, b| a.cell.cmp(&b.cell).then(a.s.total_cmp(&b.s)));
        let sheet = |cell: usize| Sheet::Side {
            sector: (cell / count) as u32,
            piece: cell % count,
        };
        let initial = positive(rims[0].1) != positive(rims[1].1);
        match sweep(stops.len(), initial)? {
            Arcs::Nothing => {}
            Arcs::Whole => {
                for cell in 0..total {
                    let q = [0., 1.].map(|s| self.distance(sheet(cell), s));
                    self.span(sheet(cell), [0., 1.], q, 0);
                }
            }
            Arcs::Spans(spans) => {
                for (a, b, wraps) in spans {
                    let [a, b] = [stops[a], stops[b]];
                    let cells = if wraps {
                        b.cell + total - a.cell
                    } else {
                        b.cell - a.cell
                    };
                    for step in 0..=cells {
                        let sheet = sheet((a.cell + step) % total);
                        let (lo, q_lo) = if step == 0 {
                            (a.s, norm(a.point))
                        } else {
                            (0., self.distance(sheet, 0.))
                        };
                        let (hi, q_hi) = if step == cells {
                            (b.s, norm(b.point))
                        } else {
                            (1., self.distance(sheet, 1.))
                        };
                        self.span(sheet, [lo, hi], [q_lo, q_hi], 0);
                    }
                }
            }
        }
        Ok(())
    }

    /// Cuts of a strip swept by azimuth from `start`.
    fn round(
        &mut self,
        sheet: Sheet,
        mut stops: Vec<Stop>,
        initial: bool,
        start: f64,
    ) -> Result<(), GeometryError> {
        stops.sort_by(|a, b| a.phi.total_cmp(&b.phi));
        match sweep(stops.len(), initial)? {
            Arcs::Nothing => {}
            Arcs::Whole => {
                for half in [0., PI] {
                    let s = [start + half, start + half + PI];
                    let q = s.map(|s| self.distance(sheet, s));
                    self.span(sheet, s, q, 0);
                }
            }
            Arcs::Spans(spans) => {
                for (a, b, wraps) in spans {
                    let [a, b] = [stops[a], stops[b]];
                    let s = [a.phi, if wraps { b.phi + TAU } else { b.phi }];
                    let q = [
                        self.along(sheet, a.phi, a.point),
                        self.along(sheet, b.phi, b.point),
                    ];
                    self.span(sheet, s, q, 0);
                }
            }
        }
        Ok(())
    }
}

fn open() -> GeometryError {
    err(
        GeometryErrorKind::InvalidResult,
        "Bevel section boundary does not close",
    )
}

pub(super) fn section(
    core: &Arc<Core>,
    depth: f64,
    _tolerance_mm: f64,
) -> Result<Vec<Loop>, GeometryError> {
    let mut plane = Plane {
        normal: core.sight,
        level: depth - core.origin[2],
    };
    // A plane containing a meridian face has a plane function that vanishes
    // identically on the face, its rims and the cap generators there. The
    // section is the region with material immediately behind the plane, so
    // such a plane is displaced by `band` behind it: the face then lies
    // strictly in front of the cut and the cut of the material behind the
    // face converges to the face, which bounds it exactly by its rays and
    // rims. The displacement is far below every drawing tolerance.
    let band = 1e-12 * core.scale;
    if core.meridians().iter().any(|&(_, _, radial, _)| {
        let normal = [-radial[1], radial[0], 0.];
        (1. - dot3(normal, plane.normal).abs()) <= 1e-12 && plane.level.abs() <= band
    }) {
        plane.level -= band;
    }
    let [lo, hi] = core.bounds;
    if depth <= lo[2] || depth >= hi[2] {
        return Ok(Vec::new());
    }
    let mut builder = Builder {
        core,
        plane,
        cuts: Vec::new(),
    };
    let start = -0.5 * core.spacing();
    let rims = End::BOTH.map(|end| builder.rim_stops(end));
    builder.side(&rims)?;
    let circles = core
        .bore
        .map(|bore| End::BOTH.map(|end| builder.bore_stops(end, bore, start)));
    for (index, end) in End::BOTH.into_iter().enumerate() {
        let (mut stops, first) = rims[index].clone();
        let below = match &circles {
            Some(circles) => {
                stops.extend(circles[index].0.iter().copied());
                circles[index].1
            }
            None => builder.value(core.bottom(end, start)),
        };
        builder.round(
            Sheet::Cap(end),
            stops,
            positive(first) != positive(below),
            start,
        )?;
    }
    if let Some(circles) = &circles {
        let stops = circles
            .iter()
            .flat_map(|(stops, _)| stops.clone())
            .collect();
        builder.round(
            Sheet::Bore,
            stops,
            positive(circles[0].1) != positive(circles[1].1),
            start,
        )?;
    }

    // Drop cuts that stay at one point, such as the turn of a cap cut
    // through the apex of its back cone.
    let weld = 1e-9 * core.scale;
    let mut cuts: Vec<Cut> = builder
        .cuts
        .into_iter()
        .filter(|cut| {
            let first = cut.local(0.);
            [0.25, 0.5, 0.75, 1.]
                .into_iter()
                .any(|t| norm(sub(cut.local(t), first)) > weld)
        })
        .collect();
    // Chain the cuts through their shared vertices.
    let reach = 1e-7 * core.scale;
    let ends: Vec<[P3; 2]> = cuts
        .iter()
        .map(|cut| [cut.local(0.), cut.local(1.)])
        .collect();
    let mut used = vec![false; cuts.len()];
    let mut chains: Vec<Vec<usize>> = Vec::new();
    for seed in 0..cuts.len() {
        if used[seed] {
            continue;
        }
        used[seed] = true;
        let mut chain = vec![seed];
        let origin = ends[seed][0];
        let mut tail = ends[seed][1];
        loop {
            let next = (0..cuts.len())
                .filter(|&other| !used[other])
                .flat_map(|other| [(other, 0), (other, 1)])
                .map(|(other, side)| (norm(sub(ends[other][side], tail)), other, side))
                .min_by(|a, b| a.0.total_cmp(&b.0));
            let closing = norm(sub(origin, tail));
            match next {
                Some((gap, other, side)) if gap <= reach && gap <= closing => {
                    used[other] = true;
                    if side == 1 {
                        cuts[other].flip = true;
                    }
                    tail = ends[other][1 - side];
                    chain.push(other);
                }
                _ if closing <= reach => break,
                _ => return Err(open()),
            }
        }
        chains.push(chain);
    }
    // Material lies on the left of each loop as the viewer sees it.
    let step = 1e-6 * core.scale;
    let mut loops = Vec::new();
    let mut slots: Vec<Option<Cut>> = cuts.into_iter().map(Some).collect();
    for mut chain in chains {
        let mut votes = 0_i32;
        for &index in &chain {
            let cut = slots[index].as_ref().ok_or_else(open)?;
            for t in [0.5, 0.25, 0.75] {
                let tangent = sub(cut.local(t + 1e-3), cut.local(t - 1e-3));
                let left = cross(plane.normal, tangent);
                let length = norm(left);
                if length == 0. {
                    continue;
                }
                let at =
                    |sign: f64| core.inside(add(cut.local(t), scaled(left, sign * step / length)));
                match (at(1.), at(-1.)) {
                    (true, false) => votes += 1,
                    (false, true) => votes -= 1,
                    _ => continue,
                }
                break;
            }
        }
        if votes < 0 {
            chain.reverse();
        }
        let mut paths: Loop = Vec::new();
        for index in chain {
            let mut cut = slots[index].take().ok_or_else(open)?;
            if votes < 0 {
                cut.flip = !cut.flip;
            }
            paths.push(Arc::new(cut));
        }
        loops.push(paths);
    }
    Ok(loops)
}
