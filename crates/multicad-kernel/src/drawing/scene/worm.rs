//! ZA worms and generated worm wheels as exact scene members.
//!
//! Both solids are slabs about their axis with a star-shaped section at
//! every height (`worm/shell.rs`). Everything is evaluated on the source
//! surfaces of the provider: the helicoidal bands of the worm thread and the
//! generating patches of the wheel.
//!
//! - Edges are the patch junctions and the section chains of both end
//!   planes. Outlines are the zero sets of `n · d` on every patch, traced in
//!   the patch domain (`worm/trace.rs`) and evaluated by bracketed root
//!   finding on the exact normal.
//! - A line of sight is intersected with the lateral boundary through the
//!   radial margin `G(t) = R(θ(t), z(t)) - r(t)`. A cell of the line is
//!   discarded when `G` provably keeps its sign on it, and a sign change in
//!   a cell narrower than the resolution is bisected on `G`. Every interval
//!   between consecutive events is decided by the membership predicate at
//!   an interior point off the boundary: `G` keeps its sign there but may
//!   touch zero where the line is tangent to the boundary inside the
//!   material. Sign certificates, each from a bound of its function's
//!   rate along the line over the cell's width:
//!   - `|G| > K w / 2` for the global Lipschitz bound `K`;
//!   - in the flank band of a wheel, where flanks are nearly radial and
//!     `K` is large, the angle to every flank at the point's radius and
//!     height;
//!   - below the cell, `min(R, c) - r` for a radius `c` above it, whose rate
//!     only involves the boundary slopes below `c`.
//!
//!   A discarded cell never holds a sign change, so the cells that reach the
//!   resolution with one, and with them the events, are those of the global
//!   bound alone.
//! - A section plane perpendicular to the axis is the exact section chain;
//!   any other plane is traced on every patch as the zero set of the signed
//!   distance, joined by the chords of both end planes and the bore.
use super::paths::{Segment, add, conic, dot3, scaled};
use super::*;
use shell::{Joint, Shell};
use std::f64::consts::{PI, TAU};
#[cfg(test)]
use trace::root;
use trace::{Seg, root_from, trace};

mod shell;
mod trace;

type V3 = [f64; 3];

fn cross(a: V3, b: V3) -> V3 {
    std::array::from_fn(|i| a[(i + 1) % 3] * b[(i + 2) % 3] - a[(i + 2) % 3] * b[(i + 1) % 3])
}

fn norm(a: V3) -> f64 {
    dot3(a, a).sqrt()
}

fn unit(a: V3) -> V3 {
    scaled(a, 1. / norm(a))
}

fn sub(a: V3, b: V3) -> V3 {
    add(a, scaled(b, -1.))
}

fn limit(message: &'static str) -> GeometryError {
    err(GeometryErrorKind::Limit, message)
}

/// The exact member of a Rust worm or generated wheel, when the solid is one.
pub(in crate::drawing) fn member(
    solid: &Solid,
    view: &View,
) -> Option<Result<Box<dyn Member>, GeometryError>> {
    Some(body(solid, view)?.map(|body| Box::new(body) as Box<dyn Member>))
}

fn body(solid: &Solid, view: &View) -> Option<Result<Body, GeometryError>> {
    let (shell, frame): (Box<dyn Shell>, PlaneFrame) = match solid.rust_worm() {
        Some((worm, frame)) => (Box::new(shell::Screw::new(worm)), frame),
        None => {
            let (wheel, frame) = solid.rust_worm_wheel()?;
            (Box::new(shell::Teeth::new(wheel)), frame)
        }
    };
    let axes = [frame.u, frame.v, frame.normal()].map(|axis| view.vector(axis));
    let [inner, outer] = shell.annulus();
    let scale = outer.max(shell.half());
    let valid =
        inner.is_finite() && scale.is_finite() && shell.lipschitz().iter().all(|k| k.is_finite());
    Some(if valid {
        Ok(Body(Arc::new(Inner {
            shell,
            origin: view.point(frame.origin),
            axes,
            toward: [axes[0][2], axes[1][2], axes[2][2]],
            scale,
            tolerance: solid.geometry_policy().linear_tolerance_mm,
        })))
    } else {
        Err(err(
            GeometryErrorKind::InvalidResult,
            "Rust worm gearing surface bounds are not finite",
        ))
    })
}

struct Inner {
    shell: Box<dyn Shell>,
    /// The local origin and axes in view space.
    origin: P3,
    axes: [P3; 3],
    /// The direction toward the viewer in local coordinates.
    toward: V3,
    scale: f64,
    tolerance: f64,
}

/// The function whose zero set is traced on a patch.
#[derive(Clone, Copy, Debug)]
enum Level {
    /// `n · d` of the unit outward normal.
    Outline,
    /// Signed height above the plane at this view depth less the origin's.
    Plane(f64),
}

impl Inner {
    fn view(&self, p: V3) -> P3 {
        std::array::from_fn(|i| {
            self.origin[i]
                + p[0] * self.axes[0][i]
                + p[1] * self.axes[1][i]
                + p[2] * self.axes[2][i]
        })
    }

    fn local(&self, q: P3) -> V3 {
        let d = sub(q, self.origin);
        self.axes.map(|axis| dot3(axis, d))
    }

    /// The level function with the bias that makes a vanishing value
    /// positive, so that a zero set along a grid line is traced once.
    fn level(&self, patch: usize, level: Level, u: f64, z: f64) -> f64 {
        let jet = self.shell.eval(patch, u, z);
        match level {
            Level::Outline => dot3(unit(jet.normal()), self.toward) + 1e-12,
            Level::Plane(c) => dot3(jet.p, self.toward) - c + 2e-13 * self.scale,
        }
    }

    /// Whether a local point is strictly inside the solid; `slab` also
    /// requires it between the end planes.
    fn inside(&self, [x, y, z]: V3, slab: bool) -> bool {
        let half = self.shell.half();
        let r = x.hypot(y);
        (!slab || z.abs() < half)
            && self.shell.bore().is_none_or(|bore| r > bore)
            && r < self.shell.radius(y.atan2(x), z.clamp(-half, half))
    }

    /// Parameter intervals of the local line `a + t d` strictly inside the
    /// solid, or inside its section at the line's height when `slab` is off.
    fn intervals(&self, a: V3, d: V3, slab: bool) -> Result<Vec<[f64; 2]>, GeometryError> {
        let shell = &self.shell;
        let half = shell.half();
        let [core, outer] = shell.annulus();
        let tau = 1e-10 * self.scale;
        let planar = d[0] * d[0] + d[1] * d[1];
        let moment = a[0] * d[1] - a[1] * d[0];
        let mut range = [f64::NEG_INFINITY, f64::INFINITY];
        if slab {
            if d[2].abs() > 1e-14 {
                let [p, q] = [-half, half].map(|z| (z - a[2]) / d[2]);
                range = [p.min(q), p.max(q)];
            } else if a[2].abs() >= half {
                return Ok(Vec::new());
            }
        }
        // Parameters between which the line is within a radius of the axis.
        let within = |radius: f64| -> Option<[f64; 2]> {
            if planar <= 1e-24 {
                return (a[0].hypot(a[1]) < radius).then_some([f64::NEG_INFINITY, f64::INFINITY]);
            }
            let disc = radius * radius * planar - moment * moment;
            (disc > 0.).then(|| {
                let mid = -(a[0] * d[0] + a[1] * d[1]) / planar;
                let reach = disc.sqrt() / planar;
                [mid - reach, mid + reach]
            })
        };
        // A line within the resolution of the outer cylinder only grazes.
        let (Some(_), Some(reach)) = (within(outer - tau), within(outer)) else {
            return Ok(Vec::new());
        };
        range = [range[0].max(reach[0]), range[1].min(reach[1])];
        if !(range[1] - range[0] > 0. && range[0].is_finite() && range[1].is_finite()) {
            return Ok(Vec::new());
        }
        let mut events = vec![range[0], range[1]];
        // An axial line within the resolution of the core cylinder lies on
        // it: its radius is taken as the core radius exactly, so the root
        // land, whose boundary radius is that same value, is a grazing
        // contact and the thread around it is material.
        let axial = planar <= 1e-24;
        let snapped = axial && (a[0].hypot(a[1]) - core).abs() <= tau;
        let mut zones = vec![range];
        if !snapped && let Some([p, q]) = within(core) {
            zones = vec![[range[0], p.min(range[1])], [q.max(range[0]), range[1]]];
            events.extend([p, q]);
        }
        if let Some(bore) = shell.bore()
            && let Some(pair) = within(bore)
        {
            events.extend(pair);
        }
        let point = |t: f64| add(a, scaled(d, t));
        if snapped && let Some(heights) = shell.core_line(a[1].atan2(a[0]), [-half, half]) {
            // `d` is axial and unit here: the parameter is the height
            // offset along it.
            let output = heights
                .into_iter()
                .map(|[z0, z1]| {
                    let [t0, t1] = [z0, z1].map(|z| (z - a[2]) / d[2]);
                    [t0.min(t1).max(range[0]), t0.max(t1).min(range[1])]
                })
                .filter(|[t0, t1]| t1 - t0 > 1e-12 * self.scale)
                .collect::<Vec<_>>();
            let mut output = output;
            output.sort_by(|x, y| x[0].total_cmp(&y[0]));
            return Ok(output);
        }
        let radius = |x: f64, y: f64| if snapped { core } else { x.hypot(y) };
        let margin = |t: f64| {
            let [x, y, z] = point(t);
            shell.radius(y.atan2(x), z.clamp(-half, half)) - radius(x, y)
        };
        let [turn, rise] = shell.lipschitz();
        let reach = planar.sqrt();
        let bound = (turn * reach / core + rise * d[2].abs() + reach) * (1. + 1e-6);
        let resolution = 2.5e-7 * self.scale;
        // In the flank band the angle to the nearest flank at the point's
        // radius and height changes along the line at a rate bounded without
        // the steep radial slope of nearly radial flanks; a cell whose middle
        // is farther from every flank than that rate over half the cell
        // keeps one side of the boundary. `None` when the cell leaves the
        // band and the radial margin decides.
        let band = shell.flank_band().filter(|_| !snapped && planar > 1e-24);
        let squared = |t: f64| {
            let [x, y, _] = point(t);
            x * x + y * y
        };
        let foot = -(a[0] * d[0] + a[1] * d[1]) / planar;
        // The radial speed along the line is at most its value at either
        // end of a cell over the cell's least radius.
        let speed = |lo: f64, hi: f64, near: f64| {
            let b = a[0] * d[0] + a[1] * d[1];
            (planar * lo + b).abs().max((planar * hi + b).abs()) / near
        };
        let separated = |lo: f64, hi: f64, mid: f64| -> Option<bool> {
            let floor = band?;
            let near = squared(foot.clamp(lo, hi)).sqrt();
            let far = squared(lo).max(squared(hi)).sqrt();
            if !(near > floor && far < outer) {
                return None;
            }
            let [slope, lean] = shell.flank_slopes(near);
            let rate =
                (moment.abs() / (near * near) + slope * speed(lo, hi, near) + lean * d[2].abs())
                    * (1. + 1e-6);
            let [x, y, z] = point(mid);
            Some(shell.separated(
                y.atan2(x),
                x.hypot(y),
                z.clamp(-half, half),
                rate * (hi - lo) / 2.,
            ))
        };
        // Below a radius `c` above the cell, `min(R, c) - r` has the sign of
        // the margin and changes only at the slopes of the boundary below
        // `c`; along the line the angle changes at `|moment| / r²` and the
        // radius at most at the radial speed at either end of the cell.
        let truncated = |lo: f64, hi: f64, mid: f64, at_mid: f64| -> bool {
            if snapped || planar <= 1e-24 {
                return false;
            }
            let near = squared(foot.clamp(lo, hi)).sqrt();
            let far = squared(lo).max(squared(hi)).sqrt();
            let speed = speed(lo, hi, near);
            let r_mid = squared(mid).sqrt();
            let boundary = at_mid + r_mid;
            let lift = far + 1e-9 * self.scale;
            [lift, boundary]
                .into_iter()
                .filter(|c| *c >= lift)
                .any(|c| {
                    let [turn, rise] = shell.slopes_below(c);
                    let rate = (turn * moment.abs() / (near * near) + rise * d[2].abs() + speed)
                        * (1. + 1e-6);
                    (boundary.min(c) - r_mid).abs() > rate * (hi - lo) / 2.
                })
        };
        let mut budget = 4_000_000_usize;
        for [from, to] in zones {
            if to <= from || to.is_nan() || from.is_nan() {
                continue;
            }
            let mut last: Option<(f64, f64)> = None;
            let mut stack: Vec<(f64, f64, Option<f64>, Option<f64>)> = vec![(from, to, None, None)];
            while let Some((lo, hi, at_lo, at_hi)) = stack.pop() {
                let width = hi - lo;
                if width <= resolution {
                    // Neighbouring cells of the resolution share an end.
                    let [at_lo, at_hi] = [(at_lo, lo), (at_hi, hi)].map(|(at, t)| {
                        at.unwrap_or_else(|| match last {
                            Some((s, value)) if s == t => value,
                            _ => margin(t),
                        })
                    });
                    last = Some((hi, at_hi));
                    if (at_lo < 0.) != (at_hi < 0.) {
                        events.push(root_from(margin, [lo, hi], [at_lo, at_hi]));
                    }
                    continue;
                }
                let mid = 0.5 * (lo + hi);
                budget = budget.checked_sub(1).ok_or_else(|| {
                    limit("Rust worm drawing exceeded its line-of-sight evaluation budget")
                })?;
                let at_mid = match separated(lo, hi, mid) {
                    Some(true) => continue,
                    Some(false) => None,
                    None => {
                        let at_mid = margin(mid);
                        // The margin cannot vanish within reach of a larger value.
                        if at_mid.abs() > bound * width / 2. || truncated(lo, hi, mid, at_mid) {
                            continue;
                        }
                        Some(at_mid)
                    }
                };
                stack.push((mid, hi, at_mid, at_hi));
                stack.push((lo, mid, at_lo, at_mid));
            }
        }
        events.retain(|t| *t >= range[0] && *t <= range[1]);
        events.sort_by(f64::total_cmp);
        let mut output: Vec<[f64; 2]> = Vec::new();
        for pair in events.windows(2) {
            // The margin keeps its sign between events, but may touch zero
            // without changing it: a line tangent to the boundary from inside
            // the material, such as an axial line on the waist of a throated
            // root. The interval is decided at the first of a few interior
            // points that is off the boundary.
            let sample = [0.5, 0.25, 0.75, 0.125, 0.875]
                .map(|f| pair[0] + (pair[1] - pair[0]) * f)
                .into_iter()
                .find(|&t| margin(t).abs() > tau)
                .unwrap_or(0.5 * (pair[0] + pair[1]));
            if pair[1] - pair[0] <= 1e-12 * self.scale
                || !self.inside(point(sample), slab)
                || (snapped && margin(sample) <= 0.)
            {
                continue;
            }
            match output.last_mut() {
                Some(last) if last[1] == pair[0] => last[1] = pair[1],
                _ => output.push([pair[0], pair[1]]),
            }
        }
        Ok(output)
    }
}

/// Uniform initial stops refined until the path stays within the tolerance
/// of its chords at the quarter points of every chord.
fn refine(point: &dyn Fn(f64) -> P3, stops: &[f64], tolerance: f64) -> Vec<f64> {
    fn split(
        point: &dyn Fn(f64) -> P3,
        [a, b]: [f64; 2],
        [pa, pb]: [P3; 2],
        tolerance: f64,
        depth: u32,
        output: &mut Vec<f64>,
    ) {
        let mid = 0.5 * (a + b);
        let pm = point(mid);
        let off = |s: f64| {
            let p = point(a + (b - a) * s);
            norm(sub(p, add(scaled(pa, 1. - s), scaled(pb, s))))
        };
        let flat = norm(sub(pm, scaled(add(pa, pb), 0.5))) <= 0.75 * tolerance
            && off(0.25) <= 0.75 * tolerance
            && off(0.75) <= 0.75 * tolerance;
        if !flat && depth < 24 {
            split(point, [a, mid], [pa, pm], tolerance, depth + 1, output);
            split(point, [mid, b], [pm, pb], tolerance, depth + 1, output);
        } else {
            output.push(b);
        }
    }
    let mut output = vec![stops[0]];
    for pair in stops.windows(2) {
        let ends = [point(pair[0]), point(pair[1])];
        split(point, [pair[0], pair[1]], ends, tolerance, 0, &mut output);
    }
    output
}

fn uniform(count: usize) -> Vec<f64> {
    (0..=count).map(|i| i as f64 / count as f64).collect()
}

/// A patch curve at a fixed `u` or a fixed height.
struct Iso {
    inner: Arc<Inner>,
    patch: usize,
    /// The fixed parameter is the height.
    level: bool,
    fixed: f64,
    range: [f64; 2],
    cells: usize,
}

impl Path for Iso {
    fn point(&self, t: f64) -> P3 {
        let free = self.range[0] + (self.range[1] - self.range[0]) * t;
        let (u, z) = if self.level {
            (free, self.fixed)
        } else {
            (self.fixed, free)
        };
        self.inner.view(self.inner.shell.eval(self.patch, u, z).p)
    }

    fn flatten(&self, tolerance_mm: f64) -> Vec<f64> {
        refine(&|t| self.point(t), &uniform(self.cells), tolerance_mm)
    }

    /// A helical junction seen along the axis is an arc of its circle, as
    /// every turn of it is: coincident turns are then compared as arcs of
    /// one conic instead of crossing each other's chords.
    fn paper(&self, [from, to]: [f64; 2]) -> Option<Curve> {
        let inner = &self.inner;
        if self.level || self.fixed != 0. || inner.toward[0].hypot(inner.toward[1]) > 1e-12 {
            return None;
        }
        let [radius, start, rate] = inner.shell.rail(self.patch)?;
        let height = |t: f64| self.range[0] + (self.range[1] - self.range[0]) * t;
        let flat = |p: P3| [p[0], p[1]];
        conic(
            flat(inner.origin),
            flat(scaled(inner.axes[0], radius)),
            flat(scaled(inner.axes[1], radius)),
            start + rate * height(from),
            rate * (height(to) - height(from)),
        )
    }
}

/// `center + cos(a) cosine + sin(a) sine` over the angles from `start`.
struct Round {
    center: P3,
    cosine: P3,
    sine: P3,
    start: f64,
    sweep: f64,
}

impl Round {
    fn angle(&self, t: f64) -> f64 {
        self.start + self.sweep * t
    }
}

impl Path for Round {
    fn point(&self, t: f64) -> P3 {
        let (sin, cos) = self.angle(t).sin_cos();
        add(
            self.center,
            add(scaled(self.cosine, cos), scaled(self.sine, sin)),
        )
    }

    fn flatten(&self, tolerance_mm: f64) -> Vec<f64> {
        // The sagitta of a chord of an ellipse is at most that of its
        // major circle.
        let radius = (dot3(self.cosine, self.cosine) + dot3(self.sine, self.sine)).sqrt();
        let step = if tolerance_mm < radius {
            2. * (1. - tolerance_mm / radius).acos()
        } else {
            PI
        };
        let count = (self.sweep.abs() / step.min(PI / 3.)).ceil().clamp(1., 1e6);
        uniform(count as usize)
    }

    fn paper(&self, [from, to]: [f64; 2]) -> Option<Curve> {
        let flat = |p: P3| [p[0], p[1]];
        conic(
            flat(self.center),
            flat(self.cosine),
            flat(self.sine),
            self.angle(from),
            (to - from) * self.sweep,
        )
        .or_else(|| Some(Curve::Line([from, to].map(|t| flat(self.point(t))))))
    }
}

/// A chain of cell segments of a level set on one patch.
struct Traced {
    inner: Arc<Inner>,
    patch: usize,
    level: Level,
    segs: Vec<Seg>,
    /// Ascending stops of the segment ends, from 0 to 1.
    stops: Vec<f64>,
    reversed: bool,
}

impl Traced {
    fn new(inner: &Arc<Inner>, patch: usize, level: Level, segs: Vec<Seg>) -> Self {
        let at = |uz: [f64; 2]| inner.shell.eval(patch, uz[0], uz[1]).p;
        let mut stops = vec![0.];
        for seg in &segs {
            // A vanishing segment keeps a strictly ascending stop.
            let length = norm(sub(at(seg.b), at(seg.a))).max(1e-12 * inner.scale);
            stops.push(stops[stops.len() - 1] + length);
        }
        let total = stops[stops.len() - 1];
        stops.iter_mut().for_each(|s| *s /= total);
        *stops.last_mut().unwrap() = 1.;
        Self {
            inner: inner.clone(),
            patch,
            level,
            segs,
            stops,
            reversed: false,
        }
    }

    /// The patch parameters at `t`.
    fn parameters(&self, t: f64) -> [f64; 2] {
        let t = if self.reversed { 1. - t } else { t }.clamp(0., 1.);
        let index = self
            .stops
            .partition_point(|stop| *stop <= t)
            .clamp(1, self.segs.len())
            - 1;
        let [a, b] = [self.stops[index], self.stops[index + 1]];
        self.segs[index].at(
            &|u, z| self.inner.level(self.patch, self.level, u, z),
            (t - a) / (b - a),
        )
    }

    fn jet(&self, t: f64) -> shell::Jet3 {
        let [u, z] = self.parameters(t);
        self.inner.shell.eval(self.patch, u, z)
    }
}

impl Path for Traced {
    fn point(&self, t: f64) -> P3 {
        self.inner.view(self.jet(t).p)
    }

    fn flatten(&self, tolerance_mm: f64) -> Vec<f64> {
        let stops: Vec<f64> = if self.reversed {
            self.stops.iter().rev().map(|s| 1. - s).collect()
        } else {
            self.stops.clone()
        };
        refine(&|t| self.point(t), &stops, tolerance_mm)
    }

    fn paper(&self, _: [f64; 2]) -> Option<Curve> {
        None
    }
}

/// A range of another path.
struct Part {
    path: Arc<dyn Path>,
    range: [f64; 2],
}

impl Part {
    fn inner(&self, t: f64) -> f64 {
        self.range[0] + (self.range[1] - self.range[0]) * t
    }
}

impl Path for Part {
    fn point(&self, t: f64) -> P3 {
        self.path.point(self.inner(t))
    }

    fn flatten(&self, tolerance_mm: f64) -> Vec<f64> {
        let [lo, hi] = self.range;
        let mut inner: Vec<f64> = self
            .path
            .flatten(tolerance_mm)
            .into_iter()
            .map(|t| (t - lo) / (hi - lo))
            .filter(|t| *t > 1e-9 && *t < 1. - 1e-9)
            .collect();
        if hi < lo {
            inner.reverse();
        }
        let mut stops = vec![0.];
        stops.extend(inner);
        stops.push(1.);
        stops
    }

    fn paper(&self, [from, to]: [f64; 2]) -> Option<Curve> {
        self.path.paper([self.inner(from), self.inner(to)])
    }
}

struct Body(Arc<Inner>);

impl Body {
    /// Split a path where its paper image turns back along its longer paper
    /// extent, so that every part is injective in the paper.
    fn monotone(&self, path: Arc<dyn Path>, tolerance: f64, output: &mut Vec<Arc<dyn Path>>) {
        let ts = path.flatten(tolerance);
        let points: Vec<P3> = ts.iter().map(|&t| path.point(t)).collect();
        let extent = |k: usize| {
            let (lo, hi) = points
                .iter()
                .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), p| {
                    (lo.min(p[k]), hi.max(p[k]))
                });
            hi - lo
        };
        let axis = usize::from(extent(1) > extent(0));
        let noise = 1e-11 * self.0.scale;
        let mut splits = vec![0.];
        let mut heading = 0.;
        let mut anchor = 0;
        for i in 1..points.len() {
            let step = points[i][axis] - points[i - 1][axis];
            if step.abs() <= noise {
                continue;
            }
            if heading * step < 0. {
                // The extreme lies between the last two significant steps.
                let (mut lo, mut hi) = (ts[anchor.max(1) - 1], ts[i]);
                let sign = if heading > 0. { 1. } else { -1. };
                let value = |t: f64| sign * path.point(t)[axis];
                for _ in 0..200 {
                    if hi - lo <= 1e-13 {
                        break;
                    }
                    let [p, q] = [lo + 0.382 * (hi - lo), hi - 0.382 * (hi - lo)];
                    if value(p) < value(q) {
                        lo = p;
                    } else {
                        hi = q;
                    }
                }
                let at = 0.5 * (lo + hi);
                if at > splits[splits.len() - 1] + 1e-9 && at < 1. - 1e-9 {
                    splits.push(at);
                }
            }
            heading = step;
            anchor = i;
        }
        splits.push(1.);
        if splits.len() == 2 {
            output.push(path);
            return;
        }
        for pair in splits.windows(2) {
            output.push(Arc::new(Part {
                path: path.clone(),
                range: [pair[0], pair[1]],
            }));
        }
    }

    /// The section chain of one patch at a height, counterclockwise about
    /// the local axis.
    fn rim(&self, patch: usize, z: f64) -> Arc<dyn Path> {
        self.rim_range(patch, z, [0., 1.])
    }

    /// A range of patch parameters of the section chain of one patch.
    fn rim_range(&self, patch: usize, z: f64, [from, to]: [f64; 2]) -> Arc<dyn Path> {
        let inner = &self.0;
        match inner.shell.arc(patch, z) {
            // Arcs are traversed uniformly in the patch parameter.
            Some([radius, start, sweep]) => Arc::new(Round {
                center: inner.view([0., 0., z]),
                cosine: scaled(inner.axes[0], radius),
                sine: scaled(inner.axes[1], radius),
                start: start + sweep * from,
                sweep: sweep * (to - from),
            }),
            None => Arc::new(Iso {
                inner: inner.clone(),
                patch,
                level: true,
                fixed: z,
                range: [from, to],
                cells: 2,
            }),
        }
    }

    fn circle(&self, radius: f64, z: f64, clockwise: bool) -> Arc<dyn Path> {
        Arc::new(Round {
            center: self.0.view([0., 0., z]),
            cosine: scaled(self.0.axes[0], radius),
            sine: scaled(self.0.axes[1], radius),
            start: 0.,
            sweep: if clockwise { -TAU } else { TAU },
        })
    }

    fn traced(&self, patch: usize, level: Level) -> Result<Vec<Traced>, GeometryError> {
        let inner = &self.0;
        let half = inner.shell.half();
        Ok(trace(
            &|u, z| inner.level(patch, level, u, z),
            [-half, half],
            inner.shell.grid(patch),
        )?
        .into_iter()
        .map(|segs| Traced::new(inner, patch, level, segs))
        .collect())
    }

    /// Boundary pieces of the section by the plane at a view depth, each
    /// with the material on its left as seen by the viewer.
    fn cut(&self, depth: f64) -> Result<Vec<Arc<dyn Path>>, GeometryError> {
        let inner = &self.0;
        let shell = &inner.shell;
        let half = shell.half();
        let toward = inner.toward;
        let level = depth - inner.origin[2];
        let lateral = toward[0].hypot(toward[1]);
        let mut pieces: Vec<Arc<dyn Path>> = Vec::new();
        let reverse = |path: Arc<dyn Path>| -> Arc<dyn Path> {
            Arc::new(Part {
                path,
                range: [1., 0.],
            })
        };
        if lateral < 1e-9 {
            // Material immediately behind the plane decides at an end plane.
            let z = level / toward[2];
            let behind = z - toward[2].signum() * 2. * inner.tolerance;
            if behind.abs() >= half {
                return Ok(pieces);
            }
            let z = z.clamp(-half, half);
            let up = toward[2] > 0.;
            for patch in 0..shell.patches() {
                let rim = self.rim(patch, z);
                pieces.push(if up { rim } else { reverse(rim) });
            }
            if let Some(bore) = shell.bore() {
                pieces.push(self.circle(bore, z, up));
            }
            return Ok(pieces);
        }
        for patch in 0..shell.patches() {
            for mut path in self.traced(patch, Level::Plane(level))? {
                // Travel along `d × n` keeps the material on the left.
                let h = 1e-4;
                let tangent = sub(path.jet(0.5 + h).p, path.jet(0.5 - h).p);
                let wanted = cross(toward, path.jet(0.5).normal());
                path.reversed = dot3(tangent, wanted) < 0.;
                pieces.push(Arc::new(path));
            }
        }
        // Chords of both end planes.
        let along = [-toward[1] / lateral, toward[0] / lateral, 0.];
        for side in [-1., 1.] {
            let z = side * half;
            let reach = (level - toward[2] * z) / (lateral * lateral);
            let base = [toward[0] * reach, toward[1] * reach, z];
            for [from, to] in inner.intervals(base, along, false)? {
                let [p, q] = [from, to].map(|t| inner.view(add(base, scaled(along, t))));
                pieces.push(Arc::new(if side > 0. {
                    Segment(q, p)
                } else {
                    Segment(p, q)
                }));
            }
        }
        if let Some(bore) = shell.bore() {
            let heading = toward[1].atan2(toward[0]);
            if toward[2].abs() < 1e-12 {
                let reach = level / (lateral * bore);
                if reach.abs() < 1. {
                    for sign in [-1., 1.] {
                        let angle = heading + sign * reach.acos();
                        let (sin, cos) = angle.sin_cos();
                        // The bore wall faces the axis: travel along
                        // `d × (-radial)`.
                        let rise = -(toward[0] * sin - toward[1] * cos);
                        let [p, q] = [-half, half]
                            .map(|z| inner.view([bore * cos, bore * sin, z * rise.signum()]));
                        pieces.push(Arc::new(Segment(p, q)));
                    }
                }
            } else {
                // The plane meets the bore in the ellipse
                // `z(θ) = (level - lateral bore cos(θ - heading)) / d_z`.
                let bounds = [-1., 1.].map(|side: f64| {
                    ((level - side * toward[2] * half) / (lateral * bore)).clamp(-1., 1.)
                });
                let [low, high] = [bounds[0].min(bounds[1]), bounds[0].max(bounds[1])];
                if high > low {
                    let center = inner.view([0., 0., level / toward[2]]);
                    let slope = -lateral * bore / toward[2];
                    let (sin, cos) = heading.sin_cos();
                    let cosine = add(
                        scaled(inner.axes[0], bore * cos),
                        add(
                            scaled(inner.axes[1], bore * sin),
                            scaled(inner.axes[2], slope),
                        ),
                    );
                    let sine = add(
                        scaled(inner.axes[0], -bore * sin),
                        scaled(inner.axes[1], bore * cos),
                    );
                    // Seen by the viewer the hole runs clockwise, which is
                    // a decreasing angle when the axis faces the viewer.
                    let turn = -toward[2].signum();
                    for [from, to] in [[high.acos(), low.acos()], [-low.acos(), -high.acos()]] {
                        let (start, sweep) = if turn > 0. {
                            (from, to - from)
                        } else {
                            (to, from - to)
                        };
                        pieces.push(Arc::new(Round {
                            center,
                            cosine,
                            sine,
                            start,
                            sweep,
                        }));
                    }
                }
            }
        }
        Ok(pieces)
    }
}

impl Member for Body {
    fn bounds(&self) -> [P3; 2] {
        let inner = &self.0;
        let reach =
            |k: usize, sign: f64| inner.shell.support(inner.axes.map(|axis| sign * axis[k]));
        [
            std::array::from_fn(|k| inner.origin[k] - reach(k, -1.)),
            std::array::from_fn(|k| inner.origin[k] + reach(k, 1.)),
        ]
    }

    fn candidates(&self, tolerance_mm: f64) -> Result<Vec<Candidate>, GeometryError> {
        let inner = &self.0;
        let shell = &inner.shell;
        let half = shell.half();
        let mut found: Vec<(Kind, Arc<dyn Path>)> = Vec::new();
        let patches = shell.patches();
        // An outline that ends on an edge is where the visibility of that
        // edge may change, so the edge is split exactly there: junctions at
        // the height of the end, end-plane rims at its patch parameter.
        let mut junction_stops: Vec<Vec<f64>> = vec![Vec::new(); patches];
        let mut rim_stops: Vec<[Vec<f64>; 2]> = vec![[Vec::new(), Vec::new()]; patches];
        let near = 1e-9;
        for patch in 0..patches {
            for path in self.traced(patch, Level::Outline)? {
                for [u, z] in [path.parameters(0.), path.parameters(1.)] {
                    if u <= near {
                        junction_stops[patch].push(z);
                    } else if u >= 1. - near {
                        junction_stops[(patch + 1) % patches].push(z);
                    }
                    if (z + half).abs() <= near * half {
                        rim_stops[patch][0].push(u);
                    } else if (z - half).abs() <= near * half {
                        rim_stops[patch][1].push(u);
                    }
                }
                found.push((Kind::Outline, Arc::new(path)));
            }
        }
        // Ascending interior stops between the two ends of a range.
        let ranges = |mut stops: Vec<f64>, [from, to]: [f64; 2]| {
            let margin = near * (to - from);
            stops.retain(|s| *s > from + margin && *s < to - margin);
            stops.sort_by(f64::total_cmp);
            stops.dedup_by(|a, b| (*a - *b).abs() <= margin);
            let mut edges = vec![from];
            edges.extend(stops);
            edges.push(to);
            edges
                .windows(2)
                .map(|pair| [pair[0], pair[1]])
                .collect::<Vec<_>>()
        };
        for patch in 0..patches {
            let kind = match shell.joint(patch) {
                Joint::Edge => Some(Kind::Edge),
                Joint::Smooth => Some(Kind::Smooth),
                Joint::Seam => None,
            };
            if let Some(kind) = kind {
                for range in ranges(std::mem::take(&mut junction_stops[patch]), [-half, half]) {
                    found.push((
                        kind,
                        Arc::new(Iso {
                            inner: inner.clone(),
                            patch,
                            level: false,
                            fixed: 0.,
                            range,
                            cells: (shell.grid(patch)[1] / 32).max(2),
                        }),
                    ));
                }
            }
            for (end, z) in [-half, half].into_iter().enumerate() {
                let stops = std::mem::take(&mut rim_stops[patch][end]);
                for range in ranges(stops, [0., 1.]) {
                    found.push((Kind::Edge, self.rim_range(patch, z, range)));
                }
            }
        }
        if let Some(bore) = shell.bore() {
            for z in [-half, half] {
                found.push((Kind::Edge, self.circle(bore, z, false)));
            }
            let toward = inner.toward;
            if toward[0].hypot(toward[1]) > 1e-9 {
                let heading = toward[1].atan2(toward[0]);
                for side in [-1., 1.] {
                    let (sin, cos) = (heading + side * PI / 2.).sin_cos();
                    let [p, q] = [-half, half].map(|z| inner.view([bore * cos, bore * sin, z]));
                    found.push((Kind::Outline, Arc::new(Segment(p, q))));
                }
            }
        }
        let mut output = Vec::new();
        for (kind, path) in found {
            let mut parts = Vec::new();
            self.monotone(path, tolerance_mm, &mut parts);
            output.extend(parts.into_iter().map(|path| Candidate { kind, path }));
        }
        Ok(output)
    }

    fn sight(&self, paper: P2) -> Result<Vec<[f64; 2]>, GeometryError> {
        let inner = &self.0;
        inner.intervals(inner.local([paper[0], paper[1], 0.]), inner.toward, true)
    }

    /// The sight-line machinery on an arbitrary line: the same isolated
    /// roots of the radial margin, each interval decided by the provider's
    /// predicate, along the unit direction and rescaled to `direction`.
    fn line(&self, start: P3, direction: P3) -> Option<Result<Vec<[f64; 2]>, GeometryError>> {
        let inner = &self.0;
        let length = norm(direction);
        if !(length.is_finite() && length > 0.) || start.iter().any(|v| !v.is_finite()) {
            return Some(Err(err(
                GeometryErrorKind::Input,
                "Rust worm line query needs a finite start and a nonzero direction",
            )));
        }
        let unit = inner.axes.map(|axis| dot3(axis, direction) / length);
        Some(
            inner
                .intervals(inner.local(start), unit, true)
                .map(|spans| {
                    spans
                        .into_iter()
                        .map(|[a, b]| [a / length, b / length])
                        .collect()
                }),
        )
    }

    fn section(&self, depth: f64, tolerance_mm: f64) -> Result<Vec<Loop>, GeometryError> {
        let mut output: Loop = Vec::new();
        for path in self.cut(depth)? {
            self.monotone(path, tolerance_mm, &mut output);
        }
        Ok(if output.is_empty() {
            vec![]
        } else {
            vec![output]
        })
    }
}

#[cfg(test)]
mod tests;
