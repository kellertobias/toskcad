//! Regularized subtraction of one finite oblique cylinder from a box.
//!
//! The result keeps its actual supports: six stock planes, the two tool cap
//! planes and the tool cylinder. Every edge is a retained intersection of
//! two supports (a line, a circle arc or an ellipse arc) and every vertex is
//! one retained intersection of three supports with a single numerical
//! source. No mesh, polygonal cylinder or fitted curve defines the boundary.
//!
//! Certification is by margin. Each classification that admits or rejects a
//! boundary piece must clear `RESOLUTION` times the largest source magnitude;
//! each transversal root must clear the same bound from tangency. Floating
//! evaluation errors are bounded by the conditioning that `MINIMUM_SLOPE`
//! admits and stay well below that margin. Anything closer is refused as
//! unresolved, so a small positive penetration never becomes contact. Exact
//! face-plane separation, including true contact, is decided in rational
//! arithmetic on the literal binary operands.
use crate::display::{DisplayMesh, MeshError, Tessellation, triangulate_region};
use num_rational::BigRational;
use spade::Point2;
use std::{
    collections::{BTreeMap, BTreeSet},
    f64::consts::TAU,
};

/// Relative margin that every admitted classification and root must clear.
pub const RESOLUTION: f64 = 1. / (1_u64 << 30) as f64;
/// Smallest admitted nonzero axis component along a stock normal.
pub const MINIMUM_SLOPE: f64 = 1e-3;
/// Wall vertices closer than this angle share one tessellation station.
const STATION_RAD: f64 = 1e-9;
const MAXIMUM_MM: f64 = 1e7;
const PLANES: usize = 6;
const CAP_START: usize = 6;
const CAP_END: usize = 7;
const WALL: usize = 8;

type V = [f64; 3];

fn dot(a: V, b: V) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: V, b: V) -> V {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// `a + scale * b`.
fn mix(a: V, scale: f64, b: V) -> V {
    std::array::from_fn(|i| a[i] + scale * b[i])
}

fn sub(a: V, b: V) -> V {
    mix(a, -1., b)
}

fn unit(axis: usize) -> V {
    std::array::from_fn(|i| f64::from(u8::from(i == axis)))
}

fn root(parents: &mut [usize], mut node: usize) -> usize {
    while parents[node] != node {
        node = parents[node];
    }
    node
}

fn join(parents: &mut [usize], a: usize, b: usize) {
    let (a, b) = (root(parents, a), root(parents, b));
    parents[a.max(b)] = a.min(b);
}

/// Finite tool in the stock chart. `axes` are the two radial directions and
/// the travel axis, a right-handed orthonormal frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Tool {
    pub base_mm: V,
    pub axes: [V; 3],
    pub radius_mm: f64,
    pub height_mm: f64,
    /// The operands are the exact source values, not a rounded chart of
    /// them. Only then is separation with true contact decided exactly;
    /// otherwise every decision must clear the certified margin.
    pub literal: bool,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Error {
    Invalid,
    Unsupported,
    Unresolved,
    Cancelled,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Support {
    Stock { axis: usize, high: bool },
    ToolStart,
    ToolEnd,
    ToolSide,
}

fn support(surface: usize) -> Support {
    match surface {
        CAP_START => Support::ToolStart,
        CAP_END => Support::ToolEnd,
        WALL => Support::ToolSide,
        _ => Support::Stock {
            axis: surface / 2,
            high: surface % 2 == 1,
        },
    }
}

/// `center + cos(angle) * cos_mm + sin(angle) * sin_mm` over the tool angle.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Arc {
    pub center_mm: V,
    pub cos_mm: V,
    pub sin_mm: V,
    pub start_rad: f64,
    pub sweep_rad: f64,
}

/// Principal form of an arc: orthonormal axes, semiaxes and ellipse angles.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Principal {
    pub center_mm: V,
    pub u: V,
    pub v: V,
    pub semiaxes_mm: [f64; 2],
    pub start_rad: f64,
    pub sweep_rad: f64,
}

impl Arc {
    pub fn point(&self, angle: f64) -> V {
        let (sin, cos) = angle.sin_cos();
        mix(mix(self.center_mm, cos, self.cos_mm), sin, self.sin_mm)
    }

    fn tangent(&self, angle: f64) -> V {
        let (sin, cos) = angle.sin_cos();
        mix(mix([0.; 3], -sin, self.cos_mm), cos, self.sin_mm)
    }

    pub fn principal(&self) -> Principal {
        let (a, b) = (self.cos_mm, self.sin_mm);
        let turn = 0.5 * (2. * dot(a, b)).atan2(dot(a, a) - dot(b, b));
        let (sin, cos) = turn.sin_cos();
        let major = mix(mix([0.; 3], cos, a), sin, b);
        let minor = mix(mix([0.; 3], -sin, a), cos, b);
        let semiaxes_mm = [major, minor].map(|axis| dot(axis, axis).sqrt());
        Principal {
            center_mm: self.center_mm,
            u: major.map(|x| x / semiaxes_mm[0]),
            v: minor.map(|x| x / semiaxes_mm[1]),
            semiaxes_mm,
            start_rad: self.start_rad - turn,
            sweep_rad: self.sweep_rad,
        }
    }

    /// Integrals of the travel coordinate and its square over an angle range.
    fn travel_moments(&self, chart: &Chart, [from, to]: [f64; 2]) -> [f64; 2] {
        let w = chart.axes[2];
        let p = dot(sub(self.center_mm, chart.base), w);
        let (q, s) = (dot(self.cos_mm, w), dot(self.sin_mm, w));
        let sweep = to - from;
        let (sin, cos) = (to.sin() - from.sin(), to.cos() - from.cos());
        let sin2 = (2. * to).sin() - (2. * from).sin();
        let cos2 = (2. * to).cos() - (2. * from).cos();
        let first = p * sweep + q * sin - s * cos;
        let second = p * p * sweep
            + 2. * p * (q * sin - s * cos)
            + (q * q + s * s) * sweep / 2.
            + (q * q - s * s) * sin2 / 4.
            - q * s * cos2 / 2.;
        [first, second]
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Edge {
    pub ends: [usize; 2],
    pub supports: [Support; 2],
    /// Incident face per support, in the same order.
    pub faces: [usize; 2],
    /// A straight segment between its ends when absent.
    pub curve: Option<Arc>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Face {
    pub support: Support,
    /// Edge uses with material on the left seen from outside the solid.
    pub loops: Vec<Vec<(usize, bool)>>,
    pub area_mm2: f64,
    /// Planar centroid, or the axis point at the mean travel of a wall face.
    pub center_mm: V,
}

#[derive(Clone, Copy, Debug)]
struct Chart {
    size: V,
    base: V,
    axes: [V; 3],
    radius: f64,
    height: f64,
    resolution: f64,
}

impl Chart {
    fn local(&self, point: V) -> V {
        let offset = sub(point, self.base);
        self.axes.map(|axis| dot(offset, axis))
    }

    fn level(&self, surface: usize) -> f64 {
        if surface.is_multiple_of(2) {
            0.
        } else {
            self.size[surface / 2]
        }
    }

    fn station(&self, surface: usize) -> f64 {
        if surface == CAP_END { self.height } else { 0. }
    }

    /// Positive on the side of a support that bounds its closed source.
    fn value(&self, surface: usize, point: V) -> f64 {
        let [x, y, travel] = self.local(point);
        match surface {
            CAP_START => travel,
            CAP_END => self.height - travel,
            WALL => self.radius - x.hypot(y),
            _ if surface.is_multiple_of(2) => point[surface / 2],
            _ => self.size[surface / 2] - point[surface / 2],
        }
    }

    /// Outward normal of the result on a support.
    fn normal(&self, surface: usize, point: V) -> V {
        let [x, y, _] = self.local(point);
        let [e1, e2, w] = self.axes;
        match surface {
            CAP_START => w,
            CAP_END => w.map(|c| -c),
            WALL => mix(mix([0.; 3], -x / self.radius, e1), -y / self.radius, e2),
            _ if surface.is_multiple_of(2) => unit(surface / 2).map(|c| -c),
            _ => unit(surface / 2),
        }
    }

    /// Complete intersection of the tool cylinder with a cap or sloped plane.
    fn arc(&self, surface: usize) -> Arc {
        let [e1, e2, w] = self.axes;
        let (travel, k1, k2) = if surface >= PLANES {
            (self.station(surface), 0., 0.)
        } else {
            let k = surface / 2;
            let travel = (self.level(surface) - self.base[k]) / w[k];
            (travel, e1[k] / w[k], e2[k] / w[k])
        };
        Arc {
            center_mm: mix(self.base, travel, w),
            cos_mm: mix(e1, -k1, w).map(|c| c * self.radius),
            sin_mm: mix(e2, -k2, w).map(|c| c * self.radius),
            start_rad: 0.,
            sweep_rad: TAU,
        }
    }

    /// Signed clearance of a candidate edge point from every other support.
    fn margin(&self, point: V, [s1, s2]: [usize; 2]) -> f64 {
        let least = |range: std::ops::Range<usize>| {
            range
                .filter(|&s| s != s1 && s != s2)
                .map(|s| self.value(s, point))
                .fold(f64::INFINITY, f64::min)
        };
        let inside = least(PLANES..WALL + 1);
        least(0..PLANES).min(if s2 >= PLANES { inside } else { -inside })
    }

    /// Least value over the supports that bound a common point of three
    /// others: the stock planes and, with a tool support, the tool itself.
    /// A clearly negative value places the point off the finite boundary.
    fn clearance(&self, s: [usize; 3], point: V) -> f64 {
        let tool = s.iter().any(|&surface| surface >= PLANES);
        (0..=WALL)
            .filter(|other| !s.contains(other) && (*other < PLANES || tool))
            .map(|other| self.value(other, point))
            .fold(f64::INFINITY, f64::min)
    }

    fn basis(&self, surface: usize) -> [V; 2] {
        let [e1, e2, _] = self.axes;
        let k = surface / 2;
        match surface {
            CAP_START => [e1, e2],
            CAP_END => [e2, e1],
            _ if surface.is_multiple_of(2) => [unit((k + 2) % 3), unit((k + 1) % 3)],
            _ => [unit((k + 1) % 3), unit((k + 2) % 3)],
        }
    }
}

fn exact(value: f64) -> BigRational {
    BigRational::from_float(value).expect("validated finite operand")
}

/// Whether one stock plane separates the closed tool from the open box. The
/// test is exact for the literal operands, so true contact is separated.
fn separated(chart: &Chart) -> bool {
    let [e1, e2, w] = chart.axes;
    let zero = exact(0.);
    (0..3).any(|k| {
        let run = exact(chart.height) * exact(w[k]);
        let (low, high) = if run < zero {
            (run, zero.clone())
        } else {
            (zero.clone(), run)
        };
        let reach = exact(chart.radius);
        let reach = &reach * &reach * (exact(e1[k]) * exact(e1[k]) + exact(e2[k]) * exact(e2[k]));
        let above = exact(chart.base[k]) + low - exact(chart.size[k]);
        let below = -(exact(chart.base[k]) + high);
        [above, below]
            .into_iter()
            .any(|gap| gap >= zero && &gap * &gap >= reach)
    })
}

struct Draft {
    supports: [usize; 2],
    ends: [usize; 2],
    curve: Option<Arc>,
    middle: V,
}

struct Arrangement {
    chart: Chart,
    /// Candidate triple points with their tool angle.
    points: Vec<(V, f64)>,
    triples: BTreeMap<[usize; 3], Vec<usize>>,
    resolved: bool,
}

fn turn_between(a: f64, b: f64) -> f64 {
    ((a - b + TAU / 2.).rem_euclid(TAU) - TAU / 2.).abs()
}

impl Arrangement {
    /// Tool angles at which the wall meets a stock plane on one cap.
    fn wall_angles(&mut self, plane: usize, cap: usize) -> Vec<f64> {
        let c = self.chart;
        let [e1, e2, w] = c.axes;
        let k = plane / 2;
        let center = mix(c.base, c.station(cap), w);
        let offset = c.level(plane) - center[k];
        let amplitude = c.radius * e1[k].hypot(e2[k]);
        let side = if offset < 0. { -1. } else { 1. };
        let phase = (side * e2[k]).atan2(side * e1[k]);
        // A grazing rim is unresolved unless its nearest point is clearly
        // off the finite boundary.
        let (sin, cos) = phase.sin_cos();
        let nearest = mix(mix(center, c.radius * cos, e1), c.radius * sin, e2);
        let grazing = (offset.abs() - amplitude).abs() <= c.resolution;
        self.resolved &= !grazing || c.clearance([plane, cap, WALL], nearest) < -c.resolution;
        if offset.abs() >= amplitude {
            return Vec::new();
        }
        let half = (offset.abs() / amplitude).acos();
        vec![
            (phase - half).rem_euclid(TAU),
            (phase + half).rem_euclid(TAU),
        ]
    }

    /// Every common point of three supports that lies on the finite part of
    /// each of them, computed once per triple.
    fn triple(&mut self, mut s: [usize; 3]) -> Vec<usize> {
        s.sort_unstable();
        if let Some(known) = self.triples.get(&s) {
            return known.clone();
        }
        let c = self.chart;
        let [e1, e2, w] = c.axes;
        let planes = s.iter().filter(|&&surface| surface < PLANES).count();
        let distinct = s[..planes]
            .windows(2)
            .all(|pair| pair[0] / 2 != pair[1] / 2);
        let mut found: Vec<(V, Option<f64>)> = Vec::new();
        if planes == 3 && distinct {
            found.push((s.map(|surface| c.level(surface)), None));
        } else if planes == 2 && distinct {
            let [i, j] = [s[0] / 2, s[1] / 2];
            let m = 3 - i - j;
            let mut start = [0.; 3];
            start[i] = c.level(s[0]);
            start[j] = c.level(s[1]);
            let [x, y, travel] = c.local(start);
            if s[2] == WALL {
                // A wall point on a stock plane parallel to the axis lies on
                // one of that plane's generators and takes its exact angle.
                let generators: Vec<f64> = [s[0], s[1]]
                    .into_iter()
                    .filter(|&surface| w[surface / 2] == 0.)
                    .flat_map(|surface| self.wall_angles(surface, CAP_START))
                    .collect();
                let lead = e1[m] * e1[m] + e2[m] * e2[m];
                let near = -(x * e1[m] + y * e2[m]) / lead;
                let reach = (x + near * e1[m]).hypot(y + near * e2[m]);
                // A grazing edge line is unresolved unless its nearest
                // point is clearly off the finite boundary.
                let grazing = (reach - c.radius).abs() <= c.resolution;
                let nearest = c.clearance(s, mix(start, near, unit(m)));
                self.resolved &= !grazing || nearest < -c.resolution;
                if reach < c.radius {
                    let half = ((c.radius * c.radius - reach * reach) / lead).sqrt();
                    for along in [near - half, near + half] {
                        let point = mix(start, along, unit(m));
                        let [x, y, _] = c.local(point);
                        let angle = y.atan2(x).rem_euclid(TAU);
                        let snapped = generators.iter().copied().min_by(|p, q| {
                            turn_between(*p, angle).total_cmp(&turn_between(*q, angle))
                        });
                        found.push((point, Some(snapped.unwrap_or(angle))));
                    }
                }
            } else if w[m] != 0. {
                let along = (c.station(s[2]) - travel) / w[m];
                found.push((mix(start, along, unit(m)), None));
            } else {
                // The cap plane is parallel to this stock edge line.
                self.resolved &= (travel - c.station(s[2])).abs() > c.resolution;
            }
        } else if planes == 1 && s[2] == WALL {
            for angle in self.wall_angles(s[0], s[1]) {
                let (sin, cos) = angle.sin_cos();
                let center = mix(c.base, c.station(s[1]), w);
                let point = mix(mix(center, c.radius * cos, e1), c.radius * sin, e2);
                found.push((point, Some(angle)));
            }
        }
        let mut ids = Vec::with_capacity(found.len());
        for (point, angle) in found {
            let [x, y, _] = c.local(point);
            let residual = s.map(|surface| c.value(surface, point).abs());
            self.resolved &= residual.into_iter().all(|v| v <= c.resolution / 64.);
            // A common point of the unbounded supports is a vertex candidate
            // only inside the stock and, with a tool support, on the tool.
            let clearance = c.clearance(s, point);
            self.resolved &= clearance.abs() > c.resolution;
            if clearance > 0. {
                ids.push(self.points.len());
                self.points
                    .push((point, angle.unwrap_or_else(|| y.atan2(x).rem_euclid(TAU))));
            }
        }
        self.triples.insert(s, ids.clone());
        ids
    }

    /// Split every support intersection curve at its triple points and keep
    /// the pieces whose middle clears every other support on the solid side.
    fn drafts(&mut self, check: &dyn Fn() -> Result<(), Error>) -> Result<Vec<Draft>, Error> {
        let c = self.chart;
        let w = c.axes[2];
        let mut drafts = Vec::new();
        for s1 in 0..WALL {
            for s2 in s1 + 1..=WALL {
                check()?;
                let ids: Vec<usize> = (0..=WALL)
                    .filter(|&s| s != s1 && s != s2)
                    .flat_map(|s| self.triple([s1, s2, s]))
                    .collect();
                let closed = s2 == WALL && (s1 >= PLANES || w[s1 / 2] != 0.);
                let arc = closed.then(|| c.arc(s1));
                let direction = if s2 < PLANES {
                    unit((6 - s1 / 2 - s2 / 2) % 3)
                } else if s2 < WALL {
                    cross(w, unit(s1 / 2 % 3))
                } else {
                    w
                };
                let mut groups: BTreeMap<u64, Vec<(f64, usize)>> = BTreeMap::new();
                for id in ids {
                    let (point, angle) = self.points[id];
                    let entry = match (closed, s2 == WALL) {
                        (true, _) => (0, angle),
                        (false, true) => (angle.to_bits(), dot(point, direction)),
                        (false, false) => (0, dot(point, direction)),
                    };
                    groups.entry(entry.0).or_default().push((entry.1, id));
                }
                if let Some(arc) = arc
                    && groups.is_empty()
                {
                    groups.insert(0, vec![(0., self.points.len())]);
                    self.points.push((arc.point(0.), 0.));
                }
                for mut group in groups.into_values() {
                    group.sort_by(|a, b| a.0.total_cmp(&b.0));
                    let count = group.len();
                    for i in 0..count - usize::from(!closed) {
                        let [(from, a), (to, b)] = [group[i], group[(i + 1) % count]];
                        let ends = [self.points[a].0, self.points[b].0];
                        let (middle, curve) = if let Some(arc) = arc {
                            let sweep_rad = if count == 1 {
                                TAU
                            } else {
                                (to - from).rem_euclid(TAU)
                            };
                            let piece = Arc {
                                start_rad: from,
                                sweep_rad,
                                ..arc
                            };
                            (arc.point(from + sweep_rad / 2.), Some(piece))
                        } else {
                            (mix(ends[0], 0.5, sub(ends[1], ends[0])), None)
                        };
                        let margin = c.margin(middle, [s1, s2]);
                        self.resolved &= margin.abs() > c.resolution;
                        if margin > 0. {
                            drafts.push(Draft {
                                supports: [s1, s2],
                                ends: [a, b],
                                curve,
                                middle,
                            });
                        }
                    }
                }
            }
        }
        Ok(drafts)
    }
}

#[derive(Clone, Copy, Debug)]
struct Span {
    start: f64,
    sweep: f64,
    /// Lower and upper wall arcs over this angle range.
    arcs: Option<[usize; 2]>,
}

/// A box with one oblique cylindrical cut and its retained boundary.
#[derive(Clone, Debug)]
pub struct ObliqueBore {
    chart: Chart,
    tool: Tool,
    vertices: Vec<V>,
    edges: Vec<Edge>,
    faces: Vec<Face>,
    spans: Vec<Span>,
    /// First station and number of spans of every edge; zero for lines.
    reach: Vec<[usize; 2]>,
    stretch: f64,
    volume: f64,
    shells: usize,
}

#[derive(Clone, Debug)]
pub enum Cut {
    /// Disjoint or zero-volume contact: the regularized stock is unchanged.
    Unchanged,
    /// The tool contains the whole stock.
    Empty,
    Body(Box<ObliqueBore>),
}

/// Subtract a finite cylinder from the box `[0, size]` in the stock chart.
pub fn subtract(
    size_mm: V,
    tool: Tool,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<Cut, Error> {
    let check = || {
        if interrupt.is_some_and(|stop| stop()) {
            return Err(Error::Cancelled);
        }
        Ok(())
    };
    check()?;
    let [e1, e2, w] = tool.axes;
    let lengths = [size_mm[0], size_mm[1], size_mm[2]]
        .into_iter()
        .chain([tool.radius_mm, tool.height_mm]);
    let frame = [
        dot(e1, e1) - 1.,
        dot(e2, e2) - 1.,
        dot(e1, e2),
        dot(sub(cross(e1, e2), w), sub(cross(e1, e2), w)).sqrt(),
    ];
    let scale = lengths
        .clone()
        .chain(tool.base_mm.map(f64::abs))
        .fold(0., f64::max);
    // NaN operands fail these comparisons and are refused with the rest.
    if !(lengths.clone().all(|v| v > 0.)
        && tool.base_mm.into_iter().all(|v| v.abs() <= MAXIMUM_MM)
        && scale <= MAXIMUM_MM
        && lengths.fold(f64::INFINITY, f64::min) >= scale * RESOLUTION * 1024.
        && frame.into_iter().all(|v| v.abs() <= 1e-12))
    {
        return Err(Error::Invalid);
    }
    if w.into_iter().filter(|&c| c == 0.).count() > 1 {
        return Err(Error::Unsupported);
    }
    if w.into_iter().any(|c| c != 0. && c.abs() < MINIMUM_SLOPE) {
        return Err(Error::Unresolved);
    }
    let chart = Chart {
        size: size_mm,
        base: tool.base_mm,
        axes: tool.axes,
        radius: tool.radius_mm,
        height: tool.height_mm,
        resolution: scale * RESOLUTION,
    };
    if tool.literal && separated(&chart) {
        return Ok(Cut::Unchanged);
    }
    let mut arrangement = Arrangement {
        chart,
        points: Vec::new(),
        triples: BTreeMap::new(),
        resolved: true,
    };
    let drafts = arrangement.drafts(&check)?;
    let tool_edges = drafts.iter().filter(|d| d.supports[1] >= PLANES).count();
    let body = assemble(&mut arrangement, tool, drafts);
    check()?;
    if !arrangement.resolved {
        return Err(Error::Unresolved);
    }
    if tool_edges == 0 {
        return Ok(if body.edges.is_empty() {
            Cut::Empty
        } else {
            Cut::Unchanged
        });
    }
    if body.shells != 1 {
        return Err(Error::Unsupported);
    }
    Ok(Cut::Body(Box::new(body)))
}

/// Link the admitted edges into oriented loops, faces and wall spans. Every
/// consistency failure clears `resolved`; nothing here panics on a refused
/// arrangement.
fn assemble(arrangement: &mut Arrangement, tool: Tool, drafts: Vec<Draft>) -> ObliqueBore {
    let chart = arrangement.chart;
    let [_, _, w] = chart.axes;
    let mut resolved = arrangement.resolved;

    // Publish only the triple points that bound an admitted edge.
    let mut ids = BTreeMap::new();
    let mut vertices = Vec::new();
    let mut angles = Vec::new();
    let mut edges = Vec::with_capacity(drafts.len());
    let mut uses: [Vec<(usize, bool)>; WALL + 1] = Default::default();
    for (index, draft) in drafts.iter().enumerate() {
        let ends = draft.ends.map(|candidate| {
            *ids.entry(candidate).or_insert_with(|| {
                let (point, angle) = arrangement.points[candidate];
                vertices.push(point);
                angles.push(angle);
                vertices.len() - 1
            })
        });
        let [s1, s2] = draft.supports;
        let tangent = match draft.curve {
            Some(arc) => arc.tangent(arc.start_rad + arc.sweep_rad / 2.),
            None => sub(vertices[ends[1]], vertices[ends[0]]),
        };
        // Only the cap and wall of the removed tool meet in a reflex edge.
        let turn = if s1 >= PLANES { 1. } else { -1. };
        let [n1, n2] = draft.supports.map(|s| chart.normal(s, draft.middle));
        let forward = turn * dot(n1, cross(tangent, n2)) > 0.;
        uses[s1].push((index, forward));
        uses[s2].push((index, !forward));
        edges.push(Edge {
            ends,
            supports: draft.supports.map(support),
            faces: [0; 2],
            curve: draft.curve,
        });
    }
    for (a, p) in vertices.iter().enumerate() {
        for q in &vertices[..a] {
            let apart = sub(*p, *q).map(f64::abs).into_iter().fold(0., f64::max);
            resolved &= apart > chart.resolution;
        }
    }

    // Stations cluster the wall vertices around the tool axis.
    let on_wall: BTreeSet<usize> = uses[WALL].iter().flat_map(|u| edges[u.0].ends).collect();
    let mut ring: Vec<(f64, usize)> = on_wall.into_iter().map(|v| (angles[v], v)).collect();
    ring.sort_by(|a, b| a.0.total_cmp(&b.0));
    let count = ring.len();
    let gap = |i: usize| (ring[(i + 1) % count].0 - ring[i].0).rem_euclid(TAU);
    let first = (0..count)
        .max_by(|&a, &b| gap(a).total_cmp(&gap(b)))
        .map_or(0, |i| (i + 1) % count);
    let mut stations: Vec<f64> = Vec::new();
    let mut station_of = BTreeMap::new();
    for k in 0..count {
        let i = (first + k) % count;
        if k == 0 || gap((i + count - 1) % count) > STATION_RAD {
            stations.push(ring[i].0);
        }
        station_of.insert(ring[i].1, stations.len() - 1);
    }
    let count = stations.len();
    let mut cover = vec![Vec::new(); count];
    let mut reach = vec![[0; 2]; edges.len()];
    let mut stretch: f64 = 0.;
    for &(index, _) in &uses[WALL] {
        let edge = &edges[index];
        if let Some(arc) = edge.curve {
            let [from, to] = edge.ends.map(|v| station_of[&v]);
            let length = if edge.ends[0] == edge.ends[1] {
                count
            } else {
                (to + count - from) % count
            };
            resolved &= length > 0;
            reach[index] = [from, length];
            for i in 0..length {
                cover[(from + i) % count].push(index);
            }
            stretch = stretch.max(
                dot(arc.cos_mm, arc.cos_mm)
                    .sqrt()
                    .hypot(dot(arc.sin_mm, arc.sin_mm).sqrt()),
            );
        }
    }
    let travel = |index: usize, angle: f64| {
        let arc = edges[index].curve.expect("wall spans cover arcs");
        dot(sub(arc.point(angle), chart.base), w)
    };
    let spans: Vec<Span> = (0..count)
        .map(|i| {
            let start = stations[i];
            let sweep = if count == 1 {
                TAU
            } else {
                (stations[(i + 1) % count] - start).rem_euclid(TAU)
            };
            resolved &= matches!(cover[i].len(), 0 | 2);
            let arcs = if let [a, b] = cover[i][..] {
                let middle = start + sweep / 2.;
                let [low, high] = [travel(a, middle), travel(b, middle)];
                resolved &= (high - low).abs() > chart.resolution;
                Some(if low < high { [a, b] } else { [b, a] })
            } else {
                None
            };
            Span { start, sweep, arcs }
        })
        .collect();

    // Loops, grouped into faces per support.
    let mut faces: Vec<Face> = Vec::new();
    let mut volume = 0.;
    for (surface, uses) in uses.iter().enumerate() {
        let start = |&(e, forward): &(usize, bool)| edges[e].ends[usize::from(!forward)];
        let end = |&(e, forward): &(usize, bool)| edges[e].ends[usize::from(forward)];
        let next: BTreeMap<usize, usize> = uses
            .iter()
            .enumerate()
            .map(|(i, u)| (start(u), i))
            .collect();
        resolved &= next.len() == uses.len();
        let mut seen = vec![false; uses.len()];
        let mut loops: Vec<Vec<(usize, bool)>> = Vec::new();
        for first in 0..uses.len() {
            let mut ring = Vec::new();
            let mut at = first;
            while !seen[at] {
                seen[at] = true;
                ring.push(uses[at]);
                resolved &= next.contains_key(&end(&uses[at]));
                at = next.get(&end(&uses[at])).copied().unwrap_or(first);
            }
            resolved &= ring.is_empty() || at == first;
            if !ring.is_empty() {
                loops.push(ring);
            }
        }
        // Signed area and first moment of each planar loop about one point.
        let normal = chart.normal(surface, chart.base);
        let reference = uses.first().map_or([0.; 3], |u| vertices[start(u)]);
        let measures: Vec<(f64, V)> = loops
            .iter()
            .map(|ring| {
                let mut area = 0.;
                let mut moment = [0.; 3];
                for u in ring {
                    let [p, q] = [vertices[start(u)], vertices[end(u)]];
                    let piece = dot(normal, cross(sub(p, reference), sub(q, reference))) / 2.;
                    area += piece;
                    moment = mix(moment, piece / 3., mix(mix(reference, 1., p), 1., q));
                    if let Some(arc) = edges[u.0].curve {
                        let sweep = arc.sweep_rad;
                        let lens = (sweep - sweep.sin()) / 2.;
                        let sign = if u.1 { 1. } else { -1. };
                        let piece = sign * lens * dot(normal, cross(arc.cos_mm, arc.sin_mm));
                        let lever = 2. * (sweep / 2.).sin().powi(3) / (3. * lens);
                        let middle = arc.point(arc.start_rad + sweep / 2.);
                        let centroid = mix(arc.center_mm, lever, sub(middle, arc.center_mm));
                        area += piece;
                        moment = mix(moment, piece, centroid);
                    }
                }
                (area, moment)
            })
            .collect();
        let mut groups: Vec<usize> = (0..loops.len()).collect();
        if surface == WALL {
            let owner: BTreeMap<usize, usize> = loops
                .iter()
                .enumerate()
                .flat_map(|(i, ring)| ring.iter().map(move |u| (u.0, i)))
                .collect();
            for [low, high] in spans.iter().filter_map(|span| span.arcs) {
                join(&mut groups, owner[&low], owner[&high]);
            }
        } else {
            let outer = measures.iter().position(|m| m.0 > 0.).unwrap_or(0);
            let outers = measures.iter().filter(|m| m.0 > 0.).count();
            for (i, measure) in measures.iter().enumerate() {
                if measure.0 <= 0. {
                    resolved &= outers == 1;
                    groups[i] = outer;
                }
            }
        }
        let first_face = faces.len();
        let mut face_of = BTreeMap::new();
        for i in 0..loops.len() {
            let group = root(&mut groups, i);
            let face = *face_of.entry(group).or_insert_with(|| {
                faces.push(Face {
                    support: support(surface),
                    loops: Vec::new(),
                    area_mm2: 0.,
                    center_mm: [0.; 3],
                });
                faces.len() - 1
            });
            for u in &loops[i] {
                let side = usize::from(edges[u.0].supports[0] != support(surface));
                edges[u.0].faces[side] = face;
            }
            faces[face].loops.push(loops[i].clone());
            faces[face].area_mm2 += measures[i].0;
            faces[face].center_mm = mix(faces[face].center_mm, 1., measures[i].1);
        }
        if surface == WALL {
            // Wall faces are measured exactly from their angle spans.
            let mut moments = vec![[0.; 2]; faces.len() - first_face];
            for span in &spans {
                if let Some([low, high]) = span.arcs {
                    let range = [span.start, span.start + span.sweep];
                    let arc = |index: usize| edges[index].curve.expect("wall spans cover arcs");
                    let [l, h] = [low, high].map(|i| arc(i).travel_moments(&chart, range));
                    let slot = &mut moments[edges[low].faces[1] - first_face];
                    *slot = [slot[0] + h[0] - l[0], slot[1] + (h[1] - l[1]) / 2.];
                }
            }
            for (face, [area, moment]) in faces[first_face..].iter_mut().zip(moments) {
                face.area_mm2 = chart.radius * area;
                face.center_mm = mix(chart.base, moment / area, w);
                volume -= chart.radius * face.area_mm2;
            }
        } else {
            let distance = dot(sub(reference, chart.base), normal);
            for face in &mut faces[first_face..] {
                face.center_mm = face.center_mm.map(|c| c / face.area_mm2);
                volume += distance * face.area_mm2;
            }
        }
        resolved &= faces[first_face..].iter().all(|face| face.area_mm2 > 0.);
    }
    volume /= 3.;
    let scale = chart.resolution / RESOLUTION;
    let stock = chart.size[0] * chart.size[1] * chart.size[2];
    let slack = chart.resolution * scale * scale;
    resolved &= edges.is_empty() || (volume > 0. && volume < stock + slack);

    let mut shells: Vec<usize> = (0..faces.len()).collect();
    for edge in &edges {
        join(&mut shells, edge.faces[0], edge.faces[1]);
    }
    let shells = (0..faces.len())
        .filter(|&face| root(&mut shells, face) == face)
        .count();
    arrangement.resolved = resolved;
    ObliqueBore {
        chart,
        tool,
        vertices,
        edges,
        faces,
        spans,
        reach,
        stretch,
        volume,
        shells,
    }
}

impl ObliqueBore {
    pub fn size_mm(&self) -> V {
        self.chart.size
    }

    pub fn tool(&self) -> Tool {
        self.tool
    }

    pub fn vertices(&self) -> &[V] {
        &self.vertices
    }

    pub fn edges(&self) -> &[Edge] {
        &self.edges
    }

    pub fn faces(&self) -> &[Face] {
        &self.faces
    }

    /// Divergence integral of the retained boundary about the tool base.
    pub fn volume_mm3(&self) -> f64 {
        self.volume
    }

    /// Floating evaluation bound of `volume_mm3` for the admitted conditioning.
    pub fn volume_error_mm3(&self) -> f64 {
        let scale = self.chart.resolution / RESOLUTION;
        self.chart.resolution * scale * scale
    }

    /// Absolute clearance that every admitted classification cleared.
    pub fn resolution_mm(&self) -> f64 {
        self.chart.resolution
    }

    /// Vertices minus edges plus faces, less the extra loops of ringed faces.
    pub fn euler_characteristic(&self) -> i64 {
        let rings: usize = self.faces.iter().map(|face| face.loops.len()).sum();
        (self.vertices.len() + 2 * self.faces.len()) as i64 - (self.edges.len() + rings) as i64
    }

    /// A point and the outward normal of a planar support.
    pub fn plane(&self, face: Support) -> Option<(V, V)> {
        let c = &self.chart;
        let surface = match face {
            Support::Stock { axis, high } => 2 * axis + usize::from(high),
            Support::ToolStart => CAP_START,
            Support::ToolEnd => CAP_END,
            Support::ToolSide => return None,
        };
        let point = if surface < PLANES {
            unit(surface / 2).map(|x| x * c.level(surface))
        } else {
            mix(c.base, c.station(surface), c.axes[2])
        };
        Some((point, c.normal(surface, point)))
    }

    /// Exact range of the solid along a direction in the stock chart.
    pub fn extent(&self, direction: V) -> [f64; 2] {
        let mut values: Vec<f64> = self.vertices.iter().map(|p| dot(*p, direction)).collect();
        for arc in self.edges.iter().filter_map(|edge| edge.curve) {
            let [q, s] = [dot(arc.cos_mm, direction), dot(arc.sin_mm, direction)];
            for sign in [-1., 1.] {
                let peak = (sign * s).atan2(sign * q);
                if (peak - arc.start_rad).rem_euclid(TAU) < arc.sweep_rad {
                    values.push(dot(arc.center_mm, direction) + sign * q.hypot(s));
                }
            }
        }
        [f64::INFINITY, f64::NEG_INFINITY].map(|seed: f64| {
            let pick: fn(f64, f64) -> f64 = if seed > 0. { f64::min } else { f64::max };
            values.iter().copied().fold(seed, pick)
        })
    }

    /// Tessellate every face from one shared sampling of each retained edge.
    pub fn mesh(&self, options: Tessellation) -> Result<DisplayMesh, MeshError> {
        options.check()?;
        let limit = options.triangle_limit;
        // The second derivative of every wall arc is bounded by `stretch`.
        // Two chords per span keep a face between an arc and its own chord
        // from collapsing below the requested deflection.
        let step = options
            .angle
            .min((8. * options.deflection / self.stretch).sqrt());
        let counts: Vec<f64> = self
            .spans
            .iter()
            .map(|span| span.arcs.map_or(0., |_| (span.sweep / step).ceil().max(2.)))
            .collect();
        if counts.iter().sum::<f64>() > limit as f64 {
            return Err(MeshError::BudgetExceeded);
        }
        let stations = self.spans.len();
        let mut mesh = DisplayMesh {
            vertices: self.vertices.iter().flatten().copied().collect(),
            indices: Vec::new(),
            face_indices: Vec::new(),
        };
        let mut samples: Vec<Vec<usize>> = Vec::with_capacity(self.edges.len());
        for (edge, [from, length]) in self.edges.iter().zip(&self.reach) {
            let mut line = vec![edge.ends[0]];
            if let Some(arc) = edge.curve {
                let mut angles: Vec<f64> = (0..*length)
                    .flat_map(|i| {
                        let j = (from + i) % stations;
                        let span = self.spans[j];
                        let count = counts[j];
                        (1..=count as usize)
                            .map(move |m| span.start + span.sweep * m as f64 / count)
                    })
                    .collect();
                angles.pop();
                for angle in angles {
                    line.push(mesh.vertices.len() / 3);
                    mesh.vertices.extend(arc.point(angle));
                }
            }
            line.push(edge.ends[1]);
            samples.push(line);
        }
        for (j, span) in self.spans.iter().enumerate() {
            let Some([low, high]) = span.arcs else {
                continue;
            };
            let count = counts[j] as usize;
            let strip = |edge: usize| {
                let from = self.reach[edge][0];
                let skip: f64 = (0..(j + stations - from) % stations)
                    .map(|i| counts[(from + i) % stations])
                    .sum();
                &samples[edge][skip as usize..=skip as usize + count]
            };
            let (l, h) = (strip(low), strip(high));
            let face = self.edges[low].faces[1] as u32;
            for m in 0..count {
                for triangle in [[l[m], h[m], h[m + 1]], [l[m], h[m + 1], l[m + 1]]] {
                    if triangle[0] != triangle[1] && triangle[1] != triangle[2] {
                        mesh.triangle(triangle, face);
                    }
                }
            }
        }
        for (index, face) in self.faces.iter().enumerate() {
            let surface = match face.support {
                Support::Stock { axis, high } => 2 * axis + usize::from(high),
                Support::ToolStart => CAP_START,
                Support::ToolEnd => CAP_END,
                Support::ToolSide => continue,
            };
            let [u, v] = self.chart.basis(surface);
            let mut ids: Vec<usize> = Vec::new();
            let mut rings: Vec<[usize; 2]> = Vec::new();
            for ring in &face.loops {
                let first = ids.len();
                for &(edge, forward) in ring {
                    let line = &samples[edge];
                    let last = line.len() - 1;
                    ids.extend((0..last).map(|i| line[if forward { i } else { last - i }]));
                }
                let length = ids.len() - first;
                rings.extend((0..length).map(|i| [first + i, first + (i + 1) % length]));
            }
            let holes = face.loops.len() - 1;
            if mesh.face_indices.len() + ids.len() + 2 * holes - 2 > limit {
                return Err(MeshError::BudgetExceeded);
            }
            let points = ids
                .iter()
                .map(|&id| {
                    let p = [0, 1, 2].map(|i| mesh.vertices[3 * id + i]);
                    Point2::new(dot(p, u), dot(p, v))
                })
                .collect();
            for triangle in triangulate_region(points, &rings, holes)?.triangles {
                mesh.triangle(triangle.map(|i| ids[i]), index as u32);
            }
        }
        Ok(mesh)
    }
}

#[cfg(test)]
mod tests;
