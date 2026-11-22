//! Box stock with one oblique cylindrical cut (CAD-354). The retained edges
//! are the provider's exact lines and cylinder-plane arcs; occlusion and
//! sections evaluate the box planes and the tool cylinder directly.
use super::paths::{Lift, ProfilePath, Segment, add, dot3, scaled};
use super::*;
use multicad_precision::PortableMath as _;
use multicad_rust_kernel::{normal_sweep, oblique_bore::Support, spur_gear::Edge};
use std::f64::consts::{PI, TAU};

pub(in crate::drawing) struct Bore {
    origin: P3,
    axes: [P3; 3],
    size: [f64; 3],
    base: P3,
    /// Two radial directions and the travel axis of the tool.
    tool: [P3; 3],
    radius: f64,
    height: f64,
    edges: Vec<Candidate>,
    tolerance: f64,
}

fn unit_circle() -> Edge {
    Edge::Planar(normal_sweep::Edge::Circle {
        center: [0.; 2],
        radius_mm: 1.,
        clockwise: false,
    })
}

fn sub(a: P3, b: P3) -> P3 {
    add(a, scaled(b, -1.))
}

/// An arc `origin + cos(a) * x + sin(a) * y` over angles, split where its
/// paper image turns back when the arc is seen edge-on.
fn arcs(lift: Lift, [from, to]: [f64; 2]) -> Vec<Arc<dyn Path>> {
    let [x, y] = [[lift.x[0], lift.x[1]], [lift.y[0], lift.y[1]]];
    let scale = (x[0].portable_hypot(x[1])).max(y[0].portable_hypot(y[1]));
    let straight = (x[0] * y[1] - x[1] * y[0]).abs() <= 1e-12 * scale * scale;
    let mut stops = vec![from];
    if straight {
        let d = if x[0].portable_hypot(x[1]) >= y[0].portable_hypot(y[1]) {
            x
        } else {
            y
        };
        let turn = (y[0] * d[0] + y[1] * d[1]).portable_atan2(x[0] * d[0] + x[1] * d[1]);
        let [lo, hi] = [from.min(to), from.max(to)];
        let mut k = ((lo - turn) / PI).ceil();
        while turn + k * PI < hi {
            if turn + k * PI > lo {
                stops.push(turn + k * PI);
            }
            k += 1.;
        }
        if to < from {
            stops[1..].reverse();
        }
    }
    stops.push(to);
    stops
        .windows(2)
        .filter(|pair| (pair[1] - pair[0]).abs() > 1e-12)
        .map(|pair| {
            Arc::new(ProfilePath {
                edge: unit_circle(),
                range: [pair[0] / TAU, pair[1] / TAU],
                lift,
                straight,
            }) as Arc<dyn Path>
        })
        .collect()
}

impl Bore {
    pub(in crate::drawing) fn new(solid: &Solid, view: &View) -> Option<Self> {
        let (body, frame) = solid.oblique_bore_body()?;
        let placed = body.drawing_box_and_tool(frame);
        let vector = |v: [f64; 3]| -> P3 {
            view.vector(std::array::from_fn(|i| {
                (0..3).map(|k| placed.axes[k][i] * v[k]).sum()
            }))
        };
        let point = |p: [f64; 3]| -> P3 { add(view.point(placed.origin), vector(p)) };
        let geometry = body.geometry();
        let chart = geometry.tool();
        let mut edges = Vec::new();
        for edge in geometry.edges() {
            let kind =
                if edge.supports[0] == edge.supports[1] && edge.supports[0] == Support::ToolSide {
                    Kind::Smooth
                } else {
                    Kind::Edge
                };
            match edge.curve {
                None => {
                    let [a, b] = edge.ends.map(|index| point(geometry.vertices()[index]));
                    edges.push(Candidate {
                        kind,
                        path: Arc::new(Segment(a, b)),
                    });
                }
                Some(arc) => {
                    let lift = Lift {
                        origin: point(arc.center_mm),
                        x: vector(arc.cos_mm),
                        y: vector(arc.sin_mm),
                    };
                    for path in arcs(lift, [arc.start_rad, arc.start_rad + arc.sweep_rad]) {
                        edges.push(Candidate { kind, path });
                    }
                }
            }
        }
        Some(Self {
            origin: view.point(placed.origin),
            axes: placed.axes.map(|axis| view.vector(axis)),
            size: placed.size,
            base: view.point(placed.base),
            tool: chart.axes.map(vector),
            radius: placed.radius,
            height: placed.height,
            edges,
            tolerance: solid.geometry_policy().linear_tolerance_mm,
        })
    }

    /// Parameters of the line `start + t * direction` inside the box.
    fn in_box(&self, start: P3, direction: P3) -> Option<[f64; 2]> {
        let mut span = [f64::NEG_INFINITY, f64::INFINITY];
        for (axis, size) in self.axes.iter().zip(self.size) {
            let at = dot3(sub(start, self.origin), *axis);
            let rate = dot3(direction, *axis);
            if rate.abs() <= 1e-14 {
                if at <= self.tolerance || at >= size - self.tolerance {
                    return None;
                }
                continue;
            }
            let ends = [-at / rate, (size - at) / rate];
            span = [
                span[0].max(ends[0].min(ends[1])),
                span[1].min(ends[0].max(ends[1])),
            ];
        }
        (span[1] > span[0]).then_some(span)
    }

    /// Parameters of the line strictly inside the finite tool cylinder.
    fn in_tool(&self, start: P3, direction: P3) -> Option<[f64; 2]> {
        let w = self.tool[2];
        let offset = sub(start, self.base);
        let [s0, rate] = [dot3(offset, w), dot3(direction, w)];
        let radial = sub(offset, scaled(w, s0));
        let across = sub(direction, scaled(w, rate));
        let a = dot3(across, across);
        let b = dot3(radial, across);
        let c = dot3(radial, radial) - self.radius * self.radius;
        let mut span = if a <= 1e-24 {
            if c >= 0. {
                return None;
            }
            [f64::NEG_INFINITY, f64::INFINITY]
        } else {
            let discriminant = b * b - a * c;
            if discriminant <= 0. {
                return None;
            }
            let root = discriminant.sqrt();
            [(-b - root) / a, (-b + root) / a]
        };
        if rate.abs() <= 1e-14 {
            if s0 <= 0. || s0 >= self.height {
                return None;
            }
        } else {
            let ends = [-s0 / rate, (self.height - s0) / rate];
            span = [
                span[0].max(ends[0].min(ends[1])),
                span[1].min(ends[0].max(ends[1])),
            ];
        }
        (span[1] > span[0]).then_some(span)
    }
}

impl Member for Bore {
    fn bounds(&self) -> [P3; 2] {
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for corner in 0..8 {
            let mut p = self.origin;
            for (bit, (axis, size)) in self.axes.iter().zip(self.size).enumerate() {
                if corner >> bit & 1 == 1 {
                    p = add(p, scaled(*axis, size));
                }
            }
            for k in 0..3 {
                lo[k] = lo[k].min(p[k]);
                hi[k] = hi[k].max(p[k]);
            }
        }
        [lo, hi]
    }

    fn candidates(&self, _: f64) -> Result<Vec<Candidate>, GeometryError> {
        let mut output = self.edges.clone();
        // The wall's outline generators, where its normal is perpendicular
        // to the sight, over the part of each generator inside the stock.
        let w = self.tool[2];
        let side = [w[1], -w[0], 0.];
        let length = side[0].portable_hypot(side[1]);
        if length > 1e-9 {
            for sign in [-1., 1.] {
                let start = add(self.base, scaled(side, sign * self.radius / length));
                if let Some([a, b]) = self.in_box(start, w) {
                    let [a, b] = [a.max(0.), b.min(self.height)];
                    if b > a {
                        output.push(Candidate {
                            kind: Kind::Outline,
                            path: Arc::new(Segment(
                                add(start, scaled(w, a)),
                                add(start, scaled(w, b)),
                            )),
                        });
                    }
                }
            }
        }
        Ok(output)
    }

    fn sight(&self, [x, y]: P2) -> Result<Vec<[f64; 2]>, GeometryError> {
        let [start, direction] = [[x, y, 0.], [0., 0., 1.]];
        let Some([a, b]) = self.in_box(start, direction) else {
            return Ok(Vec::new());
        };
        let Some([c, d]) = self.in_tool(start, direction) else {
            return Ok(vec![[a, b]]);
        };
        Ok([[a, b.min(c)], [a.max(d), b]]
            .into_iter()
            .filter(|[lo, hi]| hi > lo)
            .collect())
    }

    fn line(&self, start: P3, direction: P3) -> Option<Result<Vec<[f64; 2]>, GeometryError>> {
        let Some([a, b]) = self.in_box(start, direction) else {
            return Some(Ok(Vec::new()));
        };
        let Some([c, d]) = self.in_tool(start, direction) else {
            return Some(Ok(vec![[a, b]]));
        };
        Some(Ok([[a, b.min(c)], [a.max(d), b]]
            .into_iter()
            .filter(|[lo, hi]| hi > lo)
            .collect()))
    }

    fn section(&self, depth: f64, _: f64) -> Result<Vec<Loop>, GeometryError> {
        let mut output: Loop = Vec::new();
        let mut segment = |a: P3, b: P3| {
            if dot3(sub(b, a), sub(b, a)).sqrt() > self.tolerance {
                output.push(Arc::new(Segment(a, b)) as Arc<dyn Path>);
            }
        };
        let normal = [0., 0., 1.];
        let cross = |a: P3, b: P3| -> P3 {
            std::array::from_fn(|i| {
                a[(i + 1) % 3] * b[(i + 2) % 3] - a[(i + 2) % 3] * b[(i + 1) % 3]
            })
        };
        // A point of the line where the plane meets another plane.
        let meet = |point: P3, plane: P3| -> Option<(P3, P3)> {
            let direction = cross(normal, plane);
            let length = dot3(direction, direction);
            if length <= 1e-18 {
                return None;
            }
            // Solve z = depth and plane . (p - point) = 0 in the span of both normals.
            let c = dot3(plane, normal);
            let [d1, d2] = [depth, dot3(plane, point)];
            let k1 = (d1 - d2 * c) / (1. - c * c);
            let k2 = (d2 - d1 * c) / (1. - c * c);
            Some((
                add(scaled(normal, k1), scaled(plane, k2)),
                scaled(direction, 1. / length.sqrt()),
            ))
        };
        // Stock faces outside the tool.
        for (index, (axis, size)) in self.axes.iter().zip(self.size).enumerate() {
            for level in [0., size] {
                let Some((start, direction)) = meet(add(self.origin, scaled(*axis, level)), *axis)
                else {
                    continue;
                };
                let mut span = [f64::NEG_INFINITY, f64::INFINITY];
                for other in (0..3).filter(|other| *other != index) {
                    let at = dot3(sub(start, self.origin), self.axes[other]);
                    let rate = dot3(direction, self.axes[other]);
                    if rate.abs() <= 1e-14 {
                        if at < 0. || at > self.size[other] {
                            span = [0., -1.];
                        }
                        continue;
                    }
                    let ends = [-at / rate, (self.size[other] - at) / rate];
                    span = [
                        span[0].max(ends[0].min(ends[1])),
                        span[1].min(ends[0].max(ends[1])),
                    ];
                }
                if span[1] <= span[0] {
                    continue;
                }
                let at = |t: f64| add(start, scaled(direction, t));
                match self.in_tool(start, direction) {
                    Some([c, d]) if c < span[1] && d > span[0] => {
                        if c > span[0] {
                            segment(at(span[0]), at(c));
                        }
                        if d < span[1] {
                            segment(at(d), at(span[1]));
                        }
                    }
                    _ => segment(at(span[0]), at(span[1])),
                }
            }
        }
        // Tool caps inside the stock.
        let [e1, e2, w] = self.tool;
        for station in [0., self.height] {
            let center = add(self.base, scaled(w, station));
            let Some((start, direction)) = meet(center, w) else {
                continue;
            };
            let along = dot3(sub(center, start), direction);
            let foot = add(start, scaled(direction, along));
            let reach = self.radius * self.radius - dot3(sub(foot, center), sub(foot, center));
            if reach <= 0. {
                continue;
            }
            let half = reach.sqrt();
            if let Some([a, b]) = self.in_box(foot, direction) {
                let [a, b] = [a.max(-half), b.min(half)];
                if b > a {
                    segment(
                        add(foot, scaled(direction, a)),
                        add(foot, scaled(direction, b)),
                    );
                }
            }
        }
        // The wall between its caps and inside the stock.
        if w[2].abs() <= 1e-12 {
            let [p, q] = [self.radius * e1[2], self.radius * e2[2]];
            let reach = (depth - self.base[2]) / p.portable_hypot(q);
            if reach.abs() < 1. {
                for sign in [-1., 1.] {
                    let (sin, cos) =
                        (q.portable_atan2(p) + sign * reach.portable_acos()).portable_sin_cos();
                    let start = add(
                        self.base,
                        add(scaled(e1, self.radius * cos), scaled(e2, self.radius * sin)),
                    );
                    if let Some([a, b]) = self.in_box(start, w) {
                        let [a, b] = [a.max(0.), b.min(self.height)];
                        if b > a {
                            segment(add(start, scaled(w, a)), add(start, scaled(w, b)));
                        }
                    }
                }
            }
            return Ok(if output.is_empty() {
                vec![]
            } else {
                vec![output]
            });
        }
        let lift = Lift {
            origin: add(self.base, scaled(w, (depth - self.base[2]) / w[2])),
            x: scaled(sub(e1, scaled(w, e1[2] / w[2])), self.radius),
            y: scaled(sub(e2, scaled(w, e2[2] / w[2])), self.radius),
        };
        let at = |angle: f64| {
            let (sin, cos) = angle.portable_sin_cos();
            add(lift.origin, add(scaled(lift.x, cos), scaled(lift.y, sin)))
        };
        // Each limit is `a + b cos + c sin <= 0` along the wall curve.
        let mut limits: Vec<[f64; 3]> = Vec::new();
        let mut plane = |point: P3, outward: P3| {
            limits.push([
                dot3(sub(lift.origin, point), outward),
                dot3(lift.x, outward),
                dot3(lift.y, outward),
            ]);
        };
        plane(self.base, scaled(w, -1.));
        plane(add(self.base, scaled(w, self.height)), w);
        for (axis, size) in self.axes.iter().zip(self.size) {
            plane(self.origin, scaled(*axis, -1.));
            plane(add(self.origin, scaled(*axis, size)), *axis);
        }
        let mut stops = vec![0., TAU];
        for [a, b, c] in &limits {
            let reach = -a / b.portable_hypot(*c);
            if reach.abs() < 1. {
                for sign in [-1., 1.] {
                    stops.push(
                        (c.portable_atan2(*b) + sign * reach.portable_acos()).rem_euclid(TAU),
                    );
                }
            }
        }
        stops.sort_by(f64::total_cmp);
        let mut kept: Vec<[f64; 2]> = Vec::new();
        for pair in stops.windows(2) {
            let mid = 0.5 * (pair[0] + pair[1]);
            let (sin, cos) = mid.portable_sin_cos();
            if pair[1] - pair[0] > 1e-12
                && limits.iter().all(|[a, b, c]| a + b * cos + c * sin < 0.)
            {
                match kept.last_mut() {
                    Some(last) if last[1] == pair[0] => last[1] = pair[1],
                    _ => kept.push([pair[0], pair[1]]),
                }
            }
        }
        for [from, to] in kept {
            // Halves keep a full ellipse from closing on itself.
            for range in [[from, 0.5 * (from + to)], [0.5 * (from + to), to]] {
                if dot3(
                    sub(at(range[0]), at(range[1])),
                    sub(at(range[0]), at(range[1])),
                )
                .sqrt()
                    > self.tolerance
                {
                    output.extend(arcs(lift, range));
                }
            }
        }
        Ok(if output.is_empty() {
            vec![]
        } else {
            vec![output]
        })
    }
}
