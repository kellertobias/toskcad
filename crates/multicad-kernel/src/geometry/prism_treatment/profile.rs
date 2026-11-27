//! Planar line-and-arc profile geometry for prism edge treatments: inward
//! offsets, the corners where consecutive offsets meet, and vertical
//! corner rounding.
//!
//! A profile ring runs counterclockwise, so the material lies to the left
//! of every edge. The inset of a point from an edge is its signed distance
//! from the edge's carrier line or circle, positive into the material. The
//! inward offset of an edge by `delta` is the carrier curve where the inset
//! equals `delta`, and the gradient of the inset is the unit inward normal.
use multicad_precision::PortableMath as _;
use std::f64::consts::{PI, TAU};

pub(super) type P2 = [f64; 2];

/// Turning below this many radians counts as a tangent continuation.
const SMOOTH: f64 = 1e-9;

pub(super) fn sub(a: P2, b: P2) -> P2 {
    [a[0] - b[0], a[1] - b[1]]
}
pub(super) fn add(a: P2, b: P2) -> P2 {
    [a[0] + b[0], a[1] + b[1]]
}
pub(super) fn scale(a: P2, s: f64) -> P2 {
    [a[0] * s, a[1] * s]
}
pub(super) fn dot(a: P2, b: P2) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}
pub(super) fn cross(a: P2, b: P2) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}
pub(super) fn norm(a: P2) -> f64 {
    a[0].portable_hypot(a[1])
}
pub(super) fn unit(a: P2) -> P2 {
    scale(a, 1. / norm(a))
}
fn left(a: P2) -> P2 {
    [-a[1], a[0]]
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) enum Seg {
    Line {
        a: P2,
        b: P2,
    },
    /// A clockwise arc turns away from the material, which lies outside it.
    Arc {
        c: P2,
        r: f64,
        a: P2,
        b: P2,
        cw: bool,
    },
}

/// An inward offset curve.
#[derive(Clone, Copy, Debug)]
pub(super) enum Carrier {
    Line { point: P2, direction: P2 },
    Circle { center: P2, radius: f64 },
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Corner {
    /// The edges continue tangentially.
    Smooth,
    /// The material angle is below 180 degrees.
    Convex,
    Concave,
}

impl Seg {
    pub(super) fn start(self) -> P2 {
        match self {
            Self::Line { a, .. } | Self::Arc { a, .. } => a,
        }
    }
    pub(super) fn end(self) -> P2 {
        match self {
            Self::Line { b, .. } | Self::Arc { b, .. } => b,
        }
    }
    fn tangent(self, p: P2) -> P2 {
        match self {
            Self::Line { a, b } => unit(sub(b, a)),
            Self::Arc { c, cw, .. } => {
                let t = unit(left(sub(p, c)));
                if cw { scale(t, -1.) } else { t }
            }
        }
    }
    pub(super) fn start_tangent(self) -> P2 {
        self.tangent(self.start())
    }
    pub(super) fn end_tangent(self) -> P2 {
        self.tangent(self.end())
    }
    /// The unit inward normal of the carrier at `p`: the inset gradient.
    pub(super) fn inward(self, p: P2) -> P2 {
        match self {
            Self::Line { a, b } => left(unit(sub(b, a))),
            Self::Arc { c, cw, .. } => {
                let out = unit(sub(p, c));
                if cw { out } else { scale(out, -1.) }
            }
        }
    }
    /// Signed distance of `p` from the carrier, positive into the material.
    pub(super) fn inset(self, p: P2) -> f64 {
        match self {
            Self::Line { a, .. } => dot(self.inward(a), sub(p, a)),
            Self::Arc { c, r, cw, .. } => {
                let d = norm(sub(p, c));
                if cw { d - r } else { r - d }
            }
        }
    }
    pub(super) fn carrier(self, delta: f64) -> Carrier {
        match self {
            Self::Line { a, b } => Carrier::Line {
                point: add(a, scale(self.inward(a), delta)),
                direction: unit(sub(b, a)),
            },
            Self::Arc { c, r, cw, .. } => Carrier::Circle {
                center: c,
                radius: if cw { r + delta } else { r - delta },
            },
        }
    }
    /// The closest carrier point.
    pub(super) fn foot(self, p: P2) -> P2 {
        match self {
            Self::Line { .. } => sub(p, scale(self.inward(p), self.inset(p))),
            Self::Arc { c, r, .. } => add(c, scale(unit(sub(p, c)), r)),
        }
    }
    /// The turn of the arc from its start to `p`, in its own direction and
    /// within `(-pi, 2 pi - pi)` of the start, or the distance along a line.
    pub(super) fn position(self, p: P2) -> f64 {
        match self {
            Self::Line { a, b } => dot(sub(p, a), unit(sub(b, a))),
            Self::Arc { c, a, cw, .. } => {
                let (from, to) = (sub(a, c), sub(p, c));
                let angle = cross(from, to).portable_atan2(dot(from, to));
                let angle = if cw { -angle } else { angle };
                // Points past the end wrap to negative turns only within
                // half a turn of the start.
                let sweep = self.sweep();
                let turned = angle.rem_euclid(TAU);
                if turned > sweep + (TAU - sweep) / 2. {
                    turned - TAU
                } else {
                    turned
                }
            }
        }
    }
    /// The turn of an arc, in `(0, 2 pi)`.
    pub(super) fn sweep(self) -> f64 {
        match self {
            Self::Line { .. } => 0.,
            Self::Arc { c, a, b, cw, .. } => {
                let (from, to) = (sub(a, c), sub(b, c));
                let angle = cross(from, to).portable_atan2(dot(from, to));
                let angle = if cw { -angle } else { angle };
                let turned = angle.rem_euclid(TAU);
                if turned == 0. { TAU } else { turned }
            }
        }
    }
    /// Arc length along the edge from its start to a carrier point.
    pub(super) fn along(self, p: P2) -> f64 {
        match self {
            Self::Line { .. } => self.position(p),
            Self::Arc { r, .. } => r * self.position(p),
        }
    }
    pub(super) fn length(self) -> f64 {
        match self {
            Self::Line { a, b } => norm(sub(b, a)),
            Self::Arc { r, .. } => r * self.sweep(),
        }
    }
    pub(super) fn with_ends(self, a: P2, b: P2) -> Self {
        match self {
            Self::Line { .. } => Self::Line { a, b },
            Self::Arc { c, r, cw, .. } => Self::Arc { c, r, a, b, cw },
        }
    }
    /// The midpoint direction of an arc from its center.
    pub(super) fn bisector(self) -> P2 {
        match self {
            Self::Line { a, b } => left(unit(sub(b, a))),
            Self::Arc { c, a, cw, .. } => {
                let half = self.sweep() / 2. * if cw { -1. } else { 1. };
                let (s, k) = half.portable_sin_cos();
                let d = unit(sub(a, c));
                [k * d[0] - s * d[1], s * d[0] + k * d[1]]
            }
        }
    }
}

/// The intersection of two carriers nearest `near`.
pub(super) fn meet(p: Carrier, q: Carrier, near: P2) -> Option<P2> {
    let candidates: Vec<P2> = match (p, q) {
        (
            Carrier::Line {
                point: a,
                direction: d,
            },
            Carrier::Line {
                point: b,
                direction: e,
            },
        ) => {
            let det = cross(d, e);
            if det.abs() < SMOOTH {
                return None;
            }
            vec![add(a, scale(d, cross(sub(b, a), e) / det))]
        }
        (
            Carrier::Line {
                point: a,
                direction: d,
            },
            Carrier::Circle {
                center: c,
                radius: r,
            },
        )
        | (
            Carrier::Circle {
                center: c,
                radius: r,
            },
            Carrier::Line {
                point: a,
                direction: d,
            },
        ) => {
            let w = sub(a, c);
            let along = dot(w, d);
            let disc = along * along - (dot(w, w) - r * r);
            if disc < 0. {
                return None;
            }
            let root = disc.sqrt();
            vec![
                add(a, scale(d, -along - root)),
                add(a, scale(d, -along + root)),
            ]
        }
        (
            Carrier::Circle {
                center: c0,
                radius: r0,
            },
            Carrier::Circle {
                center: c1,
                radius: r1,
            },
        ) => {
            let d = sub(c1, c0);
            let l = norm(d);
            if l == 0. {
                return None;
            }
            let a = (r0 * r0 - r1 * r1 + l * l) / (2. * l);
            let h2 = r0 * r0 - a * a;
            if h2 < 0. {
                return None;
            }
            let (m, across) = (add(c0, scale(d, a / l)), scale(left(d), h2.sqrt() / l));
            vec![sub(m, across), add(m, across)]
        }
    };
    candidates
        .into_iter()
        .min_by(|a, b| norm(sub(*a, near)).total_cmp(&norm(sub(*b, near))))
}

/// How the ring turns from `a` into `b` at their shared vertex; `None` for
/// a cusp.
pub(super) fn corner(a: Seg, b: Seg) -> Option<Corner> {
    let (ta, tb) = (a.end_tangent(), b.start_tangent());
    let turn = cross(ta, tb);
    if turn.abs() <= SMOOTH {
        (dot(ta, tb) > 0.).then_some(Corner::Smooth)
    } else if turn > 0. {
        Some(Corner::Convex)
    } else {
        Some(Corner::Concave)
    }
}

/// Where the offset of `a` by `da` meets the offset of `b` by `db` beside
/// their shared vertex.
pub(super) fn inset_corner(a: Seg, b: Seg, kind: Corner, [da, db]: [f64; 2]) -> Option<P2> {
    let v = b.start();
    if da == 0. && db == 0. {
        return Some(v);
    }
    if kind == Corner::Smooth {
        return (da == db).then(|| add(v, scale(b.inward(v), da)));
    }
    meet(a.carrier(da), b.carrier(db), v)
}

/// The rate of an inset corner at `p` when the insets change at `rates`.
pub(super) fn inset_corner_rate(a: Seg, b: Seg, p: P2, [ra, rb]: [f64; 2]) -> Option<P2> {
    let (ga, gb) = (a.inward(p), b.inward(p));
    let det = cross(ga, gb);
    (det.abs() > SMOOTH).then(|| {
        [
            (ra * gb[1] - ga[1] * rb) / det,
            (ga[0] * rb - gb[0] * ra) / det,
        ]
    })
}

/// Where the offsets meet at the vertex between `a` and `b`: the end of
/// the offset of `a` and the start of the offset of `b`. They differ only
/// at a tangent vertex where one edge is offset and the other is not; a
/// straight step on the common normal then joins them.
pub(super) fn joint_ends(a: Seg, b: Seg, kind: Corner, [da, db]: [f64; 2]) -> Option<[P2; 2]> {
    if kind == Corner::Smooth {
        let (v, m) = (b.start(), b.inward(b.start()));
        return Some([add(v, scale(m, da)), add(v, scale(m, db))]);
    }
    inset_corner(a, b, kind, [da, db]).map(|p| [p, p])
}

/// Twice the area enclosed by the rings of offset edges. Edge `k` runs
/// along its carrier from `starts[k]` to `ends[k]`, and a straight chord
/// joins `ends[k]` to `starts[next[k]]`.
pub(super) fn twice_area(
    ring: &[Seg],
    deltas: &[f64],
    starts: &[P2],
    ends: &[P2],
    next: &[usize],
) -> f64 {
    (0..ring.len())
        .map(|k| {
            let [p, q] = [starts[k], ends[k]];
            let chord = cross(q, starts[next[k]]);
            chord
                + match ring[k].carrier(deltas[k]) {
                    Carrier::Circle { center: c, radius } => {
                        let Seg::Arc { cw, .. } = ring[k] else {
                            unreachable!("circle carrier")
                        };
                        let turn = if radius == 0. {
                            0.
                        } else {
                            let (from, to) = (sub(p, c), sub(q, c));
                            let angle = cross(from, to).portable_atan2(dot(from, to));
                            if cw {
                                -(-angle).rem_euclid(TAU)
                            } else {
                                angle.rem_euclid(TAU)
                            }
                        };
                        radius * radius * turn + c[0] * (q[1] - p[1]) - c[1] * (q[0] - p[0])
                    }
                    Carrier::Line { .. } => cross(p, q),
                }
        })
        .sum()
}

/// A vertical corner treatment of the prism.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) enum CornerCut {
    Round(f64),
    Bevel(f64),
}

/// Where a segment of the cut ring comes from.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Source {
    /// Original edge `i`, possibly shortened.
    Edge(usize),
    /// The new edge at original vertex `k`.
    Corner(usize),
}

/// Round or bevel the vertices `selected` (vertex `k` starts edge `k`).
pub(super) fn cut_corners(
    ring: &[Seg],
    selected: &[bool],
    cut: CornerCut,
    minimum: f64,
) -> Result<Vec<(Seg, Source)>, &'static str> {
    let n = ring.len();
    let mut inserted = vec![None; n];
    for k in (0..n).filter(|&k| selected[k]) {
        let (a, b) = (ring[(k + n - 1) % n], ring[k]);
        let v = b.start();
        let kind = corner(a, b).ok_or("the selected side edge meets in a cusp")?;
        if kind == Corner::Smooth {
            return Err(
                "the selected side edge joins two tangent faces and has no corner to treat",
            );
        }
        let concave = kind == Corner::Concave;
        let seg = match cut {
            CornerCut::Round(radius) => {
                let offset = if concave { -radius } else { radius };
                let center = meet(a.carrier(offset), b.carrier(offset), v)
                    .ok_or("the side fillet radius does not fit the corner")?;
                Seg::Arc {
                    c: center,
                    r: radius,
                    a: a.foot(center),
                    b: b.foot(center),
                    cw: concave,
                }
            }
            CornerCut::Bevel(distance) => {
                let (Seg::Line { .. }, Seg::Line { .. }) = (a, b) else {
                    return Err("side chamfers need two flat faces at the corner");
                };
                Seg::Line {
                    a: sub(v, scale(a.end_tangent(), distance)),
                    b: add(v, scale(b.start_tangent(), distance)),
                }
            }
        };
        inserted[k] = Some(seg);
    }
    let mut cut_ring = Vec::with_capacity(n + selected.iter().filter(|&&s| s).count());
    for (i, &seg) in ring.iter().enumerate() {
        if let Some(corner) = inserted[i] {
            cut_ring.push((corner, Source::Corner(i)));
        }
        let start = inserted[i].map_or(seg.start(), |corner| corner.end());
        let end = inserted[(i + 1) % n].map_or(seg.end(), |corner| corner.start());
        let (from, to) = (seg.along(start), seg.along(end));
        let length = seg.length();
        let slack = 1e-9 * length.max(1.);
        if from < -slack || to > length + slack || to - from < minimum {
            return Err("the side edge treatment is too large for the neighbouring faces");
        }
        cut_ring.push((seg.with_ends(start, end), Source::Edge(i)));
    }
    Ok(cut_ring)
}

/// Gauss-Legendre nodes and weights on `[-1, 1]`.
pub(super) fn gauss_legendre(count: usize) -> Vec<(f64, f64)> {
    (1..=count)
        .map(|i| {
            let mut x = (PI * (i as f64 - 0.25) / (count as f64 + 0.5)).portable_cos();
            for _ in 0..100 {
                let (mut p0, mut p1) = (1., x);
                for k in 2..=count {
                    let k = k as f64;
                    (p0, p1) = (p1, ((2. * k - 1.) * x * p1 - (k - 1.) * p0) / k);
                }
                let slope = count as f64 * (x * p1 - p0) / (x * x - 1.);
                let step = p1 / slope;
                x -= step;
                if step.abs() < 1e-15 {
                    return (x, 2. / ((1. - x * x) * slope * slope));
                }
            }
            unreachable!("Legendre roots converge")
        })
        .collect()
}
