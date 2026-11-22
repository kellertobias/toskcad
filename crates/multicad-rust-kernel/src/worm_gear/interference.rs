//! Common material volume of a placed worm and generated wheel with
//! perpendicular axes.
//!
//! The wheel is a stack of exact axial sections `S(v)`, and every section
//! plane contains the worm axis direction. In plane coordinates `(s, t)`,
//! with `t` along the worm axis, the worm material `W` is the intersection
//! of the infinite thread `T`, the end slab `C` and the outside of the bore
//! `B`. Along a line of constant `s` the thread phase is affine in `t`, so
//! `∂T` consists of explicit flank graphs `t = rate·θ(s) + φ(r(s)) + n·p`
//! between the root and crest radii and of segments of constant `s` on the
//! crest, root land and bore circles. `∂C` is the pair of graphs `t = ±L/2`.
//!
//! The section area of `S ∩ T ∩ C ∩ B` is the sum, over the boundary of
//! every member, of `∮ -t ds` along the parts that lie inside all other
//! members (Green's theorem with every boundary oriented counterclockwise).
//! Segments of constant `s` contribute nothing, so only the exact wheel
//! section pieces, the flank graphs and the end graphs are integrated. Each
//! of them is split where it crosses another member's boundary:
//!
//! - wheel pieces where the exact worm margin changes sign, found on a
//!   uniform sample and refined by bisection;
//! - flank graphs at their own end-plane crossings and at the wheel
//!   crossings found from the wheel side, which are points on the graph;
//! - end graphs at the flank crossings and at the crest, root and bore
//!   radii.
//!
//! Between splits each part is smooth and lies wholly inside or outside the
//! other members, decided at its midpoint, and is integrated with a
//! sixteen-point Gauss rule. The volume is the integral of the section area
//! over the common axial span by adaptive Gauss-Kronrod quadrature; the
//! reported error is the summed Kronrod-Gauss difference. Both are
//! estimates of a quadrature over exact sections, not interval bounds.
use super::clearance::Placement;
use super::jet::Jet;
use super::wheel::{PIECES, Piece, Section, Wheel, turn};
use super::worm::Worm;
use super::*;
use std::f64::consts::{FRAC_PI_2, PI};

/// Uniform samples per boundary piece for crossing detection.
const SAMPLES: usize = 64;
/// Distance at which a wheel crossing point is taken to lie on a graph.
const ON_GRAPH_MM: f64 = 1e-9;
const MAX_INTERVALS: usize = 4096;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Error {
    /// Invalid placements, axes that are not perpendicular, or a worm
    /// profile outside the straight-flank recipe.
    Invalid,
    Limit,
    Cancelled,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Volume {
    pub volume_mm3: f64,
    /// Summed Kronrod-Gauss difference of the axial quadrature.
    pub error_mm3: f64,
}

fn gauss() -> &'static [(f64, f64); 16] {
    static NODES: std::sync::OnceLock<[(f64, f64); 16]> = std::sync::OnceLock::new();
    NODES.get_or_init(|| {
        const N: usize = 16;
        let mut nodes = [(0., 0.); N];
        for (i, node) in nodes.iter_mut().enumerate() {
            let mut x = (PI * (i as f64 + 0.75) / (N as f64 + 0.5)).portable_cos();
            let mut derivative = 1.;
            for _ in 0..32 {
                let (mut p0, mut p1) = (1., x);
                for k in 2..=N {
                    let k = k as f64;
                    (p0, p1) = (p1, ((2. * k - 1.) * x * p1 - (k - 1.) * p0) / k);
                }
                derivative = N as f64 * (x * p1 - p0) / (x * x - 1.);
                x -= p1 / derivative;
            }
            *node = (x, 2. / ((1. - x * x) * derivative * derivative));
        }
        nodes
    })
}

fn gauss_integral(f: impl Fn(f64) -> f64, a: f64, b: f64) -> f64 {
    let (mid, half) = ((a + b) / 2., (b - a) / 2.);
    half * gauss()
        .iter()
        .map(|&(x, w)| w * f(mid + half * x))
        .sum::<f64>()
}

/// Bisect a sign change of `f` between `lo` (where `inside(lo)`) and `hi`.
fn bisect(f: impl Fn(f64) -> bool, mut lo: f64, mut hi: f64) -> f64 {
    let start = f(lo);
    for _ in 0..60 {
        let mid = (lo + hi) / 2.;
        if f(mid) == start {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    (lo + hi) / 2.
}

/// Parameters in `[a, b]` where `inside` changes, from a uniform sample.
fn changes(inside: impl Fn(f64) -> bool, a: f64, b: f64, out: &mut Vec<f64>) {
    let at = |i: usize| a + (b - a) * i as f64 / SAMPLES as f64;
    let mut previous = inside(a);
    for i in 1..=SAMPLES {
        let current = inside(at(i));
        if current != previous {
            out.push(bisect(&inside, at(i - 1), at(i)));
        }
        previous = current;
    }
}

/// A worm boundary graph `t = g(s)` over `[s0, s1]`.
#[derive(Clone, Copy, Debug, PartialEq)]
enum Graph {
    /// Lower flank of turn `n`: material above.
    Rise(i64),
    /// Upper flank of turn `n`: material below.
    Fall(i64),
    /// End plane `t = sign · L/2`.
    End(f64),
}

impl Graph {
    /// `+1` where the material lies above the graph.
    fn below(self) -> bool {
        matches!(self, Self::Rise(_)) || self == Self::End(-1.)
    }
}

/// One wheel boundary piece of a sector.
#[derive(Clone, Copy, Debug)]
enum Boundary {
    Arc { radius: f64, from: f64, to: f64 },
    Generated(Piece),
}

struct Pair<'a> {
    worm: &'a Worm,
    wheel: &'a Wheel,
    worm_at: Placement,
    wheel_at: Placement,
    m: f64,
    /// Unit wheel axis in worm local coordinates; it has no axial part.
    normal: [f64; 2],
    /// `normal · p` of the wheel origin in worm local coordinates.
    offset: f64,
    rate: f64,
    pitch: f64,
    half: f64,
    major: f64,
    root: f64,
    bore: f64,
    /// Phases `φ_fall(r) = fall[0] + (r - major) fall[1]` of the falling
    /// flank and `φ_rise` of the rising flank.
    fall: [f64; 2],
    rise: [f64; 2],
}

/// The section plane at one wheel height in worm local coordinates.
struct Plane {
    /// Unit in-plane normal toward the plane, and the plane distance `c ≥ 0`.
    normal: [f64; 2],
    c: f64,
    /// In-plane direction of increasing `s`.
    e: Point,
    /// Wheel local mm `(x, y)` to `(s, t)`: offset and Jacobian rows.
    origin: [f64; 2],
    jacobian: [[f64; 2]; 2],
    /// Orientation of that map.
    sign: f64,
}

impl<'a> Pair<'a> {
    fn new(
        worm: &'a Worm,
        worm_at: Placement,
        wheel: &'a Wheel,
        wheel_at: Placement,
    ) -> Result<Self, Error> {
        if !worm_at.valid() || !wheel_at.valid() {
            return Err(Error::Invalid);
        }
        let axis = worm_at.axes.map(|a| dot(a, wheel_at.axes[2]));
        let thread = worm.thread();
        if axis[2].abs() > 1e-9 || thread.root_circle().is_some() {
            return Err(Error::Invalid);
        }
        let length = axis[0].portable_hypot(axis[1]);
        let normal = [axis[0] / length, axis[1] / length];
        let origin = worm_at.local(wheel_at.origin_mm);
        let spec = thread.spec();
        let [_, crest_end, _, _, crest_start, _] = thread.profile_breaks();
        let [fall, rise] = [1, 3].map(|band| thread.slope_at_phase(0., band));
        Ok(Self {
            worm,
            wheel,
            worm_at,
            wheel_at,
            m: wheel.module_mm(),
            normal,
            offset: normal[0] * origin[0] + normal[1] * origin[1],
            rate: thread.helical_rate(),
            pitch: spec.pitch_mm,
            half: worm.length_mm() / 2.,
            major: spec.major_radius_mm,
            root: thread.root_radius_mm(),
            bore: worm.bore_radius_mm().unwrap_or(0.),
            fall: [crest_end, 1. / fall],
            rise: [crest_start, 1. / rise],
        })
    }

    /// Plane coordinates of the wheel section at wheel height `z_mm`, or
    /// `None` when it misses the worm crest cylinder.
    fn plane(&self, z_mm: f64) -> Option<Plane> {
        let c = self.offset + z_mm;
        if c.abs() >= self.major {
            return None;
        }
        let flip = if c < 0. { -1. } else { 1. };
        let normal = self.normal.map(|x| x * flip);
        let e = [-normal[1], normal[0], 0.];
        let local = |x: f64, y: f64| {
            let p = self.worm_at.local(self.wheel_at.world([x, y, z_mm]));
            [dot(p, e), p[2]]
        };
        let o = local(0., 0.);
        let [ux, uy] = [local(1., 0.), local(0., 1.)];
        let jacobian = [[ux[0] - o[0], uy[0] - o[0]], [ux[1] - o[1], uy[1] - o[1]]];
        let det = jacobian[0][0] * jacobian[1][1] - jacobian[0][1] * jacobian[1][0];
        Some(Plane {
            normal,
            c: c.abs(),
            e,
            origin: o,
            jacobian,
            sign: det.signum(),
        })
    }

    fn worm_point(&self, plane: &Plane, [s, t]: [f64; 2]) -> Point {
        let [nx, ny] = plane.normal;
        [
            plane.c * nx + s * plane.e[0],
            plane.c * ny + s * plane.e[1],
            t,
        ]
    }

    /// Wheel section coordinates in module units.
    fn wheel_point(&self, plane: &Plane, st: [f64; 2]) -> [f64; 2] {
        let q = self
            .wheel_at
            .local(self.worm_at.world(self.worm_point(plane, st)));
        [q[0] / self.m, q[1] / self.m]
    }

    fn in_wheel(&self, plane: &Plane, section: &Section, st: [f64; 2]) -> bool {
        let [x, y] = self.wheel_point(plane, st);
        let radius = x.portable_hypot(y);
        let bore = self.wheel.bore_radius_mm().unwrap_or(0.) / self.m;
        radius > bore
            && radius < self.wheel.blank()
            && self.wheel.section_space_margin(section, x, y) > 0.
    }

    /// Thread angle along the plane line through `s`.
    fn theta(&self, plane: &Plane, s: f64) -> f64 {
        plane.normal[1].portable_atan2(plane.normal[0]) - self.worm.phase_rad()
            + s.portable_atan2(plane.c)
    }

    fn graph(&self, plane: &Plane, graph: Graph, s: f64) -> f64 {
        let r = s.portable_hypot(plane.c);
        let lift = self.rate * self.theta(plane, s) - self.half;
        match graph {
            Graph::Rise(n) => {
                lift + self.rise[0] + (r - self.major) * self.rise[1] + (n - 1) as f64 * self.pitch
            }
            Graph::Fall(n) => {
                lift + self.fall[0] + (r - self.major) * self.fall[1] + n as f64 * self.pitch
            }
            Graph::End(sign) => sign * self.half,
        }
    }

    /// Worm margins of the end slab, the thread and the bore at `(s, t)`.
    fn worm_margins(&self, plane: &Plane, st: [f64; 2]) -> [f64; 3] {
        let p = self.worm_point(plane, st);
        let r = p[0].portable_hypot(p[1]);
        let q = self.worm.to_thread(p);
        let height = q[2].clamp(0., 2. * self.half);
        [
            self.half - st[1].abs(),
            self.worm.thread().radius(q[1].portable_atan2(q[0]), height) - r,
            r - self.bore,
        ]
    }

    fn in_worm(&self, plane: &Plane, st: [f64; 2]) -> bool {
        self.worm_margins(plane, st).iter().all(|&m| m > 0.)
    }

    fn to_plane(&self, plane: &Plane, [x, y]: [f64; 2]) -> [f64; 2] {
        let j = plane.jacobian;
        [
            plane.origin[0] + j[0][0] * x + j[0][1] * y,
            plane.origin[1] + j[1][0] * x + j[1][1] * y,
        ]
    }

    /// Point and derivative of a sector boundary piece in wheel mm.
    fn boundary(
        &self,
        section: &Section,
        sector: f64,
        boundary: Boundary,
        u: f64,
    ) -> ([f64; 2], [f64; 2]) {
        let [x, y] = match boundary {
            Boundary::Arc { radius, from, to } => {
                let (sin, cos) = (from + (to - from) * u).portable_sin_cos();
                let sweep = radius * (to - from);
                [
                    Jet {
                        value: radius * cos,
                        d: [-sweep * sin, 0.],
                    },
                    Jet {
                        value: radius * sin,
                        d: [sweep * cos, 0.],
                    },
                ]
            }
            Boundary::Generated(piece) => {
                self.wheel.piece_point(section, piece, Jet::variable(u, 0))
            }
        };
        let [x, y] = turn([x, y], <Jet as super::jet::Real>::constant(sector));
        let m = self.m;
        ([x.value * m, y.value * m], [x.d[0] * m, y.d[0] * m])
    }

    /// Sector indices whose wedge may meet the worm box in this plane.
    fn sectors(&self, plane: &Plane) -> (Vec<u32>, bool) {
        let teeth = self.wheel.teeth();
        let span = self.wheel.span();
        let reach = (self.major * self.major - plane.c * plane.c).sqrt();
        let corners = [
            [-reach, -self.half],
            [reach, -self.half],
            [reach, self.half],
            [-reach, self.half],
        ]
        .map(|st| self.wheel_point(plane, st));
        // The wheel axis meets the plane at `origin`; the nearest box point
        // to it bounds the smallest section radius met.
        let origin = plane.origin;
        let nearest = [
            origin[0].clamp(-reach, reach),
            origin[1].clamp(-self.half, self.half),
        ];
        let closest = (nearest[0] - origin[0]).portable_hypot(nearest[1] - origin[1]) / self.m;
        let bore = self.wheel.bore_radius_mm().unwrap_or(0.) / self.m;
        let with_bore = bore > 0. && closest < bore;
        if closest >= self.wheel.blank() {
            return (vec![], false);
        }
        if closest == 0. {
            return ((0..teeth).collect(), with_bore);
        }
        let [cx, cy] = corners
            .iter()
            .fold([0., 0.], |acc, c| [acc[0] + c[0] / 4., acc[1] + c[1] / 4.]);
        let center = cy.portable_atan2(cx);
        let wrap = |a: f64| (a + PI).rem_euclid(2. * PI) - PI;
        let offsets = corners.map(|[x, y]| wrap(y.portable_atan2(x) - center));
        let low = offsets.iter().copied().fold(f64::INFINITY, f64::min) - span;
        let high = offsets.iter().copied().fold(f64::NEG_INFINITY, f64::max) + span;
        let chosen = (0..teeth)
            .filter(|&k| {
                let offset = wrap(self.wheel.sector_angle(k) - center);
                offset >= low - 1e-12 && offset <= high + 1e-12
            })
            .collect();
        (chosen, with_bore)
    }

    /// Exact section area of the common material at wheel height `z_mm`.
    fn area(&self, z_mm: f64) -> f64 {
        let Some(plane) = self.plane(z_mm) else {
            return 0.;
        };
        let section = self.wheel.section(z_mm / self.m);
        let (sectors, with_bore) = self.sectors(&plane);
        let span = self.wheel.span();
        let blank = self.wheel.blank();
        let [low, high] = section.exit;
        let mut area = 0.;
        let mut hits: Vec<[f64; 2]> = vec![];
        let mut cuts = vec![];
        for k in sectors {
            let angle = self.wheel.sector_angle(k);
            let mut pieces = vec![Boundary::Arc {
                radius: blank,
                from: -span,
                to: low,
            }];
            pieces.extend(PIECES.map(Boundary::Generated));
            pieces.push(Boundary::Arc {
                radius: blank,
                from: high,
                to: span,
            });
            if with_bore {
                pieces.push(Boundary::Arc {
                    radius: self.wheel.bore_radius_mm().unwrap_or(0.) / self.m,
                    from: span,
                    to: -span,
                });
            }
            for piece in pieces {
                let at = |u: f64| self.to_plane(&plane, self.boundary(&section, angle, piece, u).0);
                let inside = |u: f64| self.in_worm(&plane, at(u));
                cuts.clear();
                cuts.push(0.);
                changes(inside, 0., 1., &mut cuts);
                cuts.push(1.);
                for &u in &cuts[1..cuts.len() - 1] {
                    hits.push(at(u));
                }
                for w in cuts.windows(2) {
                    let [a, b] = [w[0], w[1]];
                    if b <= a || !inside((a + b) / 2.) {
                        continue;
                    }
                    area += plane.sign
                        * gauss_integral(
                            |u| {
                                let (p, d) = self.boundary(&section, angle, piece, u);
                                let [_, t] = self.to_plane(&plane, p);
                                let j = plane.jacobian;
                                -t * (j[0][0] * d[0] + j[0][1] * d[1])
                            },
                            a,
                            b,
                        );
                }
            }
        }
        area + self.graph_area(&plane, &section, &hits)
    }

    /// Contribution of the worm flank and end graphs inside the wheel.
    fn graph_area(&self, plane: &Plane, section: &Section, hits: &[[f64; 2]]) -> f64 {
        let reach = (self.major * self.major - plane.c * plane.c).sqrt();
        let inner =
            (plane.c < self.root).then(|| (self.root * self.root - plane.c * plane.c).sqrt());
        let flank_domains: Vec<[f64; 2]> = match inner {
            Some(inner) => vec![[-reach, -inner], [inner, reach]],
            None => vec![[-reach, reach]],
        };
        let spread = self.rate.abs() * FRAC_PI_2 + self.pitch;
        let base = self.rate * self.theta(plane, 0.);
        let first = ((-base - spread - self.pitch) / self.pitch).floor() as i64 - 1;
        let last = ((2. * self.half - base + spread + self.pitch) / self.pitch).ceil() as i64 + 1;
        let mut graphs: Vec<(Graph, [f64; 2])> = vec![];
        for n in first..=last {
            for &domain in &flank_domains {
                graphs.push((Graph::Rise(n), domain));
                graphs.push((Graph::Fall(n), domain));
            }
        }
        // Flank crossings of the end planes split both graphs.
        let mut end_cuts: Vec<f64> = [self.root, self.bore]
            .iter()
            .filter(|&&r| r > plane.c)
            .flat_map(|&r| {
                let s = (r * r - plane.c * plane.c).sqrt();
                [-s, s]
            })
            .collect();
        let mut flank_cuts: Vec<Vec<f64>> = vec![vec![]; graphs.len()];
        for (i, &(graph, [a, b])) in graphs.iter().enumerate() {
            for sign in [-1., 1.] {
                let level = sign * self.half;
                let above = |s: f64| self.graph(plane, graph, s) > level;
                let before = flank_cuts[i].len();
                changes(above, a, b, &mut flank_cuts[i]);
                end_cuts.extend_from_slice(&flank_cuts[i][before..]);
            }
        }
        for sign in [-1., 1.] {
            graphs.push((Graph::End(sign), [-reach, reach]));
            flank_cuts.push(end_cuts.clone());
        }
        let mut area = 0.;
        for ((graph, [a, b]), mut cuts) in graphs.into_iter().zip(flank_cuts) {
            cuts.extend(hits.iter().filter_map(|&[s, t]| {
                ((a..=b).contains(&s) && (self.graph(plane, graph, s) - t).abs() <= ON_GRAPH_MM)
                    .then_some(s)
            }));
            cuts.retain(|&s| s > a && s < b);
            cuts.push(a);
            cuts.push(b);
            cuts.sort_by(f64::total_cmp);
            for w in cuts.windows(2) {
                let [s0, s1] = [w[0], w[1]];
                if s1 <= s0 {
                    continue;
                }
                let s = (s0 + s1) / 2.;
                let st = [s, self.graph(plane, graph, s)];
                let [end, thread, bore] = self.worm_margins(plane, st);
                let others = match graph {
                    Graph::End(_) => thread > 0. && bore > 0.,
                    _ => end > 0. && bore > 0.,
                };
                if !others || !self.in_wheel(plane, section, st) {
                    continue;
                }
                let integral = gauss_integral(|s| self.graph(plane, graph, s), s0, s1);
                area += if graph.below() { -integral } else { integral };
            }
        }
        area
    }
}

/// Volume of the common material of a placed worm and generated wheel.
/// `tolerance_mm3` bounds the estimated axial quadrature error.
pub fn interference_volume(
    worm: &Worm,
    worm_at: Placement,
    wheel: &Wheel,
    wheel_at: Placement,
    tolerance_mm3: f64,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<Volume, Error> {
    if tolerance_mm3.is_nan() || tolerance_mm3 <= 0. {
        return Err(Error::Invalid);
    }
    let pair = Pair::new(worm, worm_at, wheel, wheel_at)?;
    let width = wheel.half_width() * pair.m;
    let low = (-pair.major - pair.offset).max(-width);
    let high = (pair.major - pair.offset).min(width);
    if high <= low {
        return Ok(Volume {
            volume_mm3: 0.,
            error_mm3: 0.,
        });
    }
    // Topology changes where the plane meets the bore, root and axis.
    let mut splits = vec![low, high];
    for r in [0., pair.root, pair.bore] {
        for z in [r - pair.offset, -r - pair.offset] {
            if z > low && z < high {
                splits.push(z);
            }
        }
    }
    splits.sort_by(f64::total_cmp);
    splits.dedup();
    let cancelled = || interrupt.is_some_and(|f| f());
    let mut intervals = vec![];
    for w in splits.windows(2) {
        if cancelled() {
            return Err(Error::Cancelled);
        }
        intervals.push(kronrod(&|z| pair.area(z), w[0], w[1]));
    }
    loop {
        let error: f64 = intervals.iter().map(|i| i.error).sum();
        if error <= tolerance_mm3 {
            let volume = intervals.iter().map(|i| i.value).sum();
            return Ok(Volume {
                volume_mm3: volume,
                error_mm3: error,
            });
        }
        if intervals.len() >= MAX_INTERVALS {
            return Err(Error::Limit);
        }
        if cancelled() {
            return Err(Error::Cancelled);
        }
        let (worst, _) = intervals
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.error.total_cmp(&b.1.error))
            .expect("nonempty");
        let Interval { a, b, .. } = intervals.swap_remove(worst);
        let mid = (a + b) / 2.;
        intervals.push(kronrod(&|z| pair.area(z), a, mid));
        intervals.push(kronrod(&|z| pair.area(z), mid, b));
    }
}

struct Interval {
    a: f64,
    b: f64,
    value: f64,
    error: f64,
}

/// Seven-point Gauss and fifteen-point Kronrod rules on `[a, b]`.
fn kronrod(f: &dyn Fn(f64) -> f64, a: f64, b: f64) -> Interval {
    const X: [f64; 8] = [
        0.991_455_371_120_812_6,
        0.949_107_912_342_758_5,
        0.864_864_423_359_769_1,
        0.741_531_185_599_394_4,
        0.586_087_235_467_691_1,
        0.405_845_151_377_397_2,
        0.207_784_955_007_898_5,
        0.,
    ];
    const K: [f64; 8] = [
        0.022_935_322_010_529_22,
        0.063_092_092_629_978_55,
        0.104_790_010_322_250_18,
        0.140_653_259_715_525_92,
        0.169_004_726_639_267_9,
        0.190_350_578_064_785_4,
        0.204_432_940_075_298_9,
        0.209_482_141_084_727_83,
    ];
    const G: [f64; 4] = [
        0.129_484_966_168_869_7,
        0.279_705_391_489_276_7,
        0.381_830_050_505_118_9,
        0.417_959_183_673_469_4,
    ];
    let (mid, half) = ((a + b) / 2., (b - a) / 2.);
    let center = f(mid);
    let mut kronrod = K[7] * center;
    let mut gauss = G[3] * center;
    for i in 0..7 {
        let pair = f(mid - half * X[i]) + f(mid + half * X[i]);
        kronrod += K[i] * pair;
        if i % 2 == 1 {
            gauss += G[i / 2] * pair;
        }
    }
    Interval {
        a,
        b,
        value: kronrod * half,
        error: ((kronrod - gauss) * half).abs(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::worm_gear::{wheel::tests::spec as wheel_spec, worm::tests::spec as worm_spec};

    /// Wheel axis +Z at the origin, worm axis along wheel +Y at center
    /// distance `(teeth/2 + 5) m`, the wheel turned by `extra` from mesh.
    fn pose(teeth: u32, drive: f64, extra: f64) -> [Placement; 2] {
        let center = f64::from(teeth) / 2. + 5.;
        let (sw, cw) = drive.portable_sin_cos();
        let turn = -drive / f64::from(teeth) + extra;
        let (s, c) = turn.portable_sin_cos();
        [
            Placement {
                origin_mm: [center, 0., 0.],
                axes: [[cw, 0., -sw], [-sw, 0., -cw], [0., 1., 0.]],
            },
            Placement {
                origin_mm: [0.; 3],
                axes: [[c, s, 0.], [-s, c, 0.], [0., 0., 1.]],
            },
        ]
    }

    fn members() -> (Worm, Wheel) {
        let worm = Worm::new(worm_spec(1., 16., 1., 0.05), 1e-3).unwrap();
        let wheel = Wheel::new(wheel_spec(1., 32, 2., 2., 0.05), 1e-3, 1e7).unwrap();
        (worm, wheel)
    }

    fn moved(p: Placement, by: &dyn Fn(Point) -> Point, origin: Point) -> Placement {
        let base = by([0.; 3]);
        Placement {
            origin_mm: by(origin),
            axes: p.axes.map(|a| sub(by(a), base)),
        }
    }

    #[test]
    fn misphased_section_areas_match_a_dense_membership_grid() {
        let (worm, wheel) = members();
        let [a, b] = pose(32, 0.4, std::f64::consts::PI / 32.);
        let pair = Pair::new(&worm, a, &wheel, b).unwrap();
        for z in [-0.83, -0.2, 0.35, 0.9] {
            let plane = pair.plane(z).unwrap();
            let section = wheel.section(z);
            let reach = (pair.major * pair.major - plane.c * plane.c).sqrt();
            let n = 700;
            let [ds, dt] = [2. * reach / n as f64, 2. * pair.half / n as f64];
            let mut grid = 0.;
            for i in 0..n {
                for j in 0..n {
                    let st = [
                        -reach + (i as f64 + 0.5) * ds,
                        -pair.half + (j as f64 + 0.5) * dt,
                    ];
                    if pair.in_worm(&plane, st) && pair.in_wheel(&plane, &section, st) {
                        grid += ds * dt;
                    }
                }
            }
            let exact = pair.area(z);
            assert!(exact > 0.05, "section {z} area {exact}");
            assert!(
                (exact - grid).abs() <= 2e-3 * exact,
                "section {z}: {exact} vs grid {grid}"
            );
        }
    }

    #[test]
    fn mating_pose_has_no_volume_and_misphase_is_rigidly_invariant() {
        let (worm, wheel) = members();
        let [a, b] = pose(32, 0.4, 0.);
        let clear = interference_volume(&worm, a, &wheel, b, 1e-7, None).unwrap();
        assert_eq!(clear.volume_mm3, 0.);
        let [a, b] = pose(32, 0.4, std::f64::consts::PI / 32.);
        let found = interference_volume(&worm, a, &wheel, b, 1e-7, None).unwrap();
        assert!(
            found.volume_mm3 > 0.01 && found.error_mm3 <= 1e-7,
            "{found:?}"
        );
        let rigid = |p: Point| {
            let q = rotate_z(p, 0.7);
            let (s, c) = 0.3_f64.portable_sin_cos();
            [
                q[0] + 3.,
                c * q[1] - s * q[2] - 2.,
                s * q[1] + c * q[2] + 5.,
            ]
        };
        let again = interference_volume(
            &worm,
            moved(a, &rigid, a.origin_mm),
            &wheel,
            moved(b, &rigid, b.origin_mm),
            1e-7,
            None,
        )
        .unwrap();
        assert!(
            (again.volume_mm3 - found.volume_mm3).abs() <= 1e-6,
            "{again:?} vs {found:?}"
        );
        // Interruption and invalid poses stay typed.
        let stop = || true;
        assert_eq!(
            interference_volume(&worm, a, &wheel, b, 1e-7, Some(&stop)),
            Err(Error::Cancelled)
        );
        let tilted = Placement {
            axes: [[1., 0., 0.], [0., 0.8, 0.6], [0., -0.6, 0.8]],
            ..b
        };
        assert_eq!(
            interference_volume(&worm, a, &wheel, tilted, 1e-7, None),
            Err(Error::Invalid)
        );
        assert_eq!(
            interference_volume(&worm, a, &wheel, b, 0., None),
            Err(Error::Invalid)
        );
    }
}
