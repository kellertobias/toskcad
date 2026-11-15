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
        let first = match self.chart {
            Chart::Sweep(_) => 8,
            Chart::Along { .. } => 1,
        };
        let mut stops = vec![0.];
        for i in 0..first {
            let [a, b] = [i, i + 1].map(|k| k as f64 / first as f64);
            self.refine(
                [a, b],
                [self.local(a), self.local(b)],
                tolerance_mm,
                0,
                &mut stops,
            );
        }
        *cache = Some((tolerance_mm.to_bits(), stops.clone()));
        stops
    }

    fn paper(&self, range: [f64; 2]) -> Option<Curve> {
        match self.chart {
            Chart::Along { s: [a, b], .. } if a == b => {
                let [p, q] = range.map(|t| self.point(t));
                Some(Curve::Line([[p[0], p[1]], [q[0], q[1]]]))
            }
            _ => None,
        }
    }
}

impl Cut {
    /// Append the stops of `(a, b]`, halving until the curve stays within
    /// half the tolerance of each chord at its quarter points.
    fn refine(
        &self,
        [a, b]: [f64; 2],
        [p, q]: [P3; 2],
        tolerance: f64,
        depth: u32,
        stops: &mut Vec<f64>,
    ) {
        let chord = sub(q, p);
        let length = dot3(chord, chord);
        let mid = 0.5 * (a + b);
        let middle = self.local(mid);
        let off = |point: P3| {
            let offset = sub(point, p);
            let along = if length > 0. {
                (dot3(offset, chord) / length).clamp(0., 1.)
            } else {
                0.
            };
            norm(sub(offset, scaled(chord, along)))
        };
        let worst = [0.25, 0.75]
            .map(|f| off(self.local(a + (b - a) * f)))
            .into_iter()
            .fold(off(middle), f64::max);
        if worst <= 0.5 * tolerance || depth >= 24 {
            stops.push(b);
            return;
        }
        self.refine([a, mid], [p, middle], tolerance, depth + 1, stops);
        self.refine([mid, b], [middle, q], tolerance, depth + 1, stops);
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
        let limit = 1e4 * core.scale;
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
    let plane = Plane {
        normal: core.sight,
        level: depth - core.origin[2],
    };
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
