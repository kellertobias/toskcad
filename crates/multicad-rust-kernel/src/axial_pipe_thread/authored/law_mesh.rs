//! Dedicated tessellation of the authored tapered pipe law.
//!
//! The law is a radial graph `r = A(q, z)` of the helical phase
//! `q = z - theta / omega` with `A = a(q) + k z`, optionally capped by a stock
//! crest limit. The mesh lives in the `(q, theta)` chart: helical lines at
//! every profile breakpoint, stock crease and sub-band, sampled on one common
//! `theta` grid and joined by meridian edges at fixed `theta`. Every cell then
//! lies in one smooth piece; the flat caps `z = 0, L` and the stock onset
//! circle are straight lines in this chart and clip cells exactly. Per
//! triangle a second-order interpolation bound gives the chord deviation and
//! an explicit normal bound gives the facet angle; vertex radii are checked
//! against the law's own interval enclosure. The grid is refined until both
//! bounds meet the request or the triangle frontier refuses.
use super::*;
use crate::display::{DisplayMesh, Tessellation};
use crate::manufacturing::{self, MeshProof};

/// Material inside the law surface (a pipe end) or outside it (a bore).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LawClosure {
    /// Lateral surface and two planar cap fans, outward normals.
    Solid,
    /// Lateral surface only, oriented for material outside; rims returned.
    Bore,
}
/// Face role of the emitted triangles. `face_indices` hold
/// [`LawMesh::role_index`] values.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LawRole {
    /// Authored band, numbered as in the raw retained boundary.
    Band(usize),
    /// NPT flat root (a cone); rounded sections have a root band instead.
    Root,
    /// Stock crest cylinder of an active crest limit.
    Clip,
    StartCap,
    EndCap,
}
pub struct LawMesh {
    pub mesh: DisplayMesh,
    band_count: usize,
    /// Rim vertex ids at `z = 0` and `z = L`, counterclockwise about `+z`.
    pub rims: [Vec<u32>; 2],
    /// Proven distance bound between every emitted facet and the authored
    /// law surface it replaces, including vertex evaluation error.
    pub deviation_mm: f64,
    /// Proven bound of the angle between each facet normal and every exact
    /// law normal over that facet.
    pub normal_angle_rad: f64,
    /// Connectivity proof; present for a closed solid only.
    pub proof: Option<MeshProof>,
    pub rows: usize,
    pub columns: usize,
}
impl LawMesh {
    pub fn role_index(&self, role: LawRole) -> u32 {
        let n = self.band_count as u32;
        match role {
            LawRole::Band(b) => b as u32,
            LawRole::Root => n,
            LawRole::Clip => n + 1,
            LawRole::StartCap => n + 2,
            LawRole::EndCap => n + 3,
        }
    }
    pub fn role(&self, index: u32) -> Option<LawRole> {
        let n = self.band_count as u32;
        Some(match index {
            i if i < n => LawRole::Band(i as usize),
            i if i == n => LawRole::Root,
            i if i == n + 1 => LawRole::Clip,
            i if i == n + 2 => LawRole::StartCap,
            i if i == n + 3 => LawRole::EndCap,
            _ => return None,
        })
    }
}

#[derive(Clone, Copy)]
enum Profile {
    Constant(f64),
    Line {
        q0: f64,
        a0: f64,
        slope: f64,
    },
    /// `a = a0 + sign sqrt(r^2 - x^2) - k x`, `x = q - qc`.
    Arc {
        qc: f64,
        a0: f64,
        radius: f64,
        sign: f64,
    },
}
/// A column boundary: fixed phase, or the stock crease `base +- c+(z)/slope`.
#[derive(Clone, Copy)]
enum Track {
    Fixed(f64),
    Crease { base: f64, sign: f64, slope: f64 },
}
#[derive(Clone, Copy)]
struct Piece {
    profile: Profile,
    start: Track,
    end: Track,
    role: u32,
    /// Crest piece whose material above the onset is the stock cylinder.
    clipped: bool,
    /// |a'| and |a''| bounds over the piece.
    d1: f64,
    d2: f64,
    /// Bound of |a''| / (1 + (a' + k)^2) over the piece.
    turn: f64,
    r_max: f64,
    /// Rounded band: sub-lines are uniform in the arc angle
    /// `phi`, `x = q - qc = radius sin(phi)`, `phi` in `[phi0, phi1]`.
    arc: Option<(f64, f64, f64, f64)>,
}
struct Law {
    k: f64,
    omega: f64,
    pitch: f64,
    length: f64,
    limit: Option<f64>,
    major: f64,
    onset: Option<f64>,
    pieces: Vec<Piece>,
    minor_min: f64,
}
impl Law {
    fn cylinder(&self, piece: &Piece, z: f64) -> bool {
        piece.clipped && self.onset.is_some_and(|o| z >= o)
    }
    fn excess(&self, z: f64) -> f64 {
        self.limit
            .map_or(0., |l| (self.major + self.k * z - l).max(0.))
    }
    fn track(&self, t: Track, z: f64) -> f64 {
        match t {
            Track::Fixed(q) => q,
            Track::Crease { base, sign, slope } => base + sign * self.excess(z) / slope,
        }
    }
    /// Phase of a track on the helical chart line `theta`: the crease
    /// `q = base + sign c+(z)/slope` with `z = q + theta/omega` is linear in
    /// `theta` above the onset and constant below it.
    fn track_theta(&self, t: Track, theta: f64) -> f64 {
        match t {
            Track::Fixed(q) => q,
            Track::Crease { base, sign, slope } => {
                let Some(limit) = self.limit else {
                    return base;
                };
                let z = base + theta / self.omega;
                if self.major + self.k * z - limit <= 0. {
                    return base;
                }
                (base + sign * (self.major - limit + self.k * theta / self.omega) / slope)
                    / (1. - sign * self.k / slope)
            }
        }
    }
    /// Radius and dA/dq, dA/dz of the piece at `(q, z)`.
    fn radius(&self, piece: &Piece, q: f64, z: f64) -> (f64, f64, f64) {
        if self.cylinder(piece, z) {
            return (self.limit.expect("clipped piece has a limit"), 0., 0.);
        }
        let (a, da) = match piece.profile {
            Profile::Constant(a) => (a, 0.),
            Profile::Line { q0, a0, slope } => (a0 + slope * (q - q0), slope),
            Profile::Arc {
                qc,
                a0,
                radius,
                sign,
            } => {
                let x = q - qc;
                let s = (radius * radius - x * x).max(0.).sqrt();
                (a0 + sign * s - self.k * x, -sign * x / s - self.k)
            }
        };
        (a + self.k * z, da, self.k)
    }
}

fn build_law(thread: &AuthoredThread, closure: LawClosure) -> Result<Law, Error> {
    let spec = thread.spec;
    let p = spec.pitch_mm;
    let r = spec.major_radius_mm;
    let k = spec.radial_taper;
    let l = spec.length_mm;
    let d = thread.depth();
    let omega = std::f64::consts::TAU / p;
    let limit = match closure {
        LawClosure::Solid => spec.crest_radius_limit_mm,
        // A bore's crest limit would remove cavity, not stock: not a bore law.
        LawClosure::Bore if spec.crest_radius_limit_mm.is_some() => {
            return Err(Error::Unsupported);
        }
        LawClosure::Bore => None,
    };
    let major_max = r + (k * l).max(0.);
    let minor_min = r - d + (k * l).min(0.);
    let onset = limit.and_then(|lim| {
        // Clip is active where r + k z > limit; only increasing tapers clip.
        (k > 0. && major_max > lim).then(|| ((lim - r) / k).max(0.))
    });
    let capped = |v: f64| limit.map_or(v, |lim| v.min(lim));
    let mut pieces = Vec::new();
    match spec.section {
        AuthoredSection::NptFlat {
            flat_width_mm: w, ..
        } => {
            let m = 3f64.sqrt();
            let h = w / 2.;
            let fall = m + k;
            let rise = m - k;
            let flank = |role, profile, start, end, slope: f64| Piece {
                profile,
                start,
                end,
                role,
                clipped: false,
                d1: slope.abs(),
                d2: 0.,
                turn: 0.,
                r_max: capped(major_max),
                arc: None,
            };
            pieces.push(Piece {
                profile: Profile::Constant(r),
                start: Track::Crease {
                    base: -h,
                    sign: -1.,
                    slope: rise,
                },
                end: Track::Crease {
                    base: h,
                    sign: 1.,
                    slope: fall,
                },
                role: 1,
                clipped: onset.is_some(),
                d1: 0.,
                d2: 0.,
                turn: 0.,
                r_max: capped(major_max),
                arc: None,
            });
            pieces.push(flank(
                2,
                Profile::Line {
                    q0: h,
                    a0: r,
                    slope: -fall,
                },
                Track::Crease {
                    base: h,
                    sign: 1.,
                    slope: fall,
                },
                Track::Fixed(h + d / fall),
                fall,
            ));
            pieces.push(Piece {
                r_max: major_max - d,
                ..flank(
                    3,
                    Profile::Constant(r - d),
                    Track::Fixed(h + d / fall),
                    Track::Fixed(p - h - d / rise),
                    0.,
                )
            });
            pieces.push(flank(
                0,
                Profile::Line {
                    q0: p - h,
                    a0: r,
                    slope: rise,
                },
                Track::Fixed(p - h - d / rise),
                Track::Crease {
                    base: p - h,
                    sign: -1.,
                    slope: rise,
                },
                rise,
            ));
        }
        AuthoredSection::BsptRounded { .. } => {
            if limit.is_some() {
                return Err(Error::Unsupported);
            }
            let b = thread.rounded()?;
            let mid = |v: I| v.lo + (v.hi - v.lo) * 0.5;
            let (ra, normal, sin, m, valley) = (
                mid(b.radius),
                mid(b.normal),
                mid(b.sin),
                mid(b.m),
                mid(b.valley),
            );
            let h = mid(b.radius * b.cos);
            let centre = r - ra * normal;
            let tangent = centre + ra * sin;
            let s = (ra * ra - h * h).sqrt();
            let arc_d1 = h / s + k.abs();
            let arc_d2 = ra * ra / (s * s * s);
            let arc = |role, profile: Profile, a: f64, b: f64| {
                let Profile::Arc { qc, .. } = profile else {
                    unreachable!("arc piece")
                };
                let phi = |q: f64| ((q - qc) / ra).clamp(-1., 1.).portable_asin();
                Piece {
                    profile,
                    start: Track::Fixed(a),
                    end: Track::Fixed(b),
                    role,
                    clipped: false,
                    d1: arc_d1,
                    d2: arc_d2,
                    // |a''|/(1+(a'+k)^2) = 1/sqrt(ra^2 - x^2) on these arcs.
                    turn: 1. / s,
                    r_max: major_max,
                    arc: Some((qc, ra, phi(a), phi(b))),
                }
            };
            let line = |role, profile, a, b, slope: f64| Piece {
                profile,
                start: Track::Fixed(a),
                end: Track::Fixed(b),
                role,
                clipped: false,
                d1: slope.abs(),
                d2: 0.,
                turn: 0.,
                r_max: major_max,
                arc: None,
            };
            pieces.push(arc(
                0,
                Profile::Arc {
                    qc: 0.,
                    a0: centre,
                    radius: ra,
                    sign: 1.,
                },
                0.,
                h,
            ));
            pieces.push(line(
                1,
                Profile::Line {
                    q0: h,
                    a0: tangent - k * h,
                    slope: -(m + k),
                },
                h,
                valley - h,
                m + k,
            ));
            pieces.push(arc(
                2,
                Profile::Arc {
                    qc: valley,
                    a0: r - d + ra * normal,
                    radius: ra,
                    sign: -1.,
                },
                valley - h,
                valley + h,
            ));
            pieces.push(line(
                3,
                Profile::Line {
                    q0: p - h,
                    a0: tangent + k * h,
                    slope: m - k,
                },
                valley + h,
                p - h,
                m - k,
            ));
            pieces.push(arc(
                4,
                Profile::Arc {
                    qc: p,
                    a0: centre,
                    radius: ra,
                    sign: 1.,
                },
                p - h,
                p,
            ));
        }
    }
    if minor_min.is_nan()
        || minor_min <= 0.
        || pieces
            .iter()
            .any(|p| !p.d1.is_finite() || !p.d2.is_finite())
    {
        return Err(Error::Unresolved);
    }
    Ok(Law {
        k,
        omega,
        pitch: p,
        length: l,
        limit,
        major: r,
        onset,
        pieces,
        minor_min,
    })
}

type V3 = [f64; 3];
fn sub(a: V3, b: V3) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn cross(a: V3, b: V3) -> V3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn dot(a: V3, b: V3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn norm(a: V3) -> f64 {
    dot(a, a).sqrt()
}
fn angle(a: V3, b: V3) -> f64 {
    // atan2 of |a x b| and a.b is accurate for small angles.
    norm(cross(a, b)).portable_atan2(dot(a, b))
}

/// Bound of |dA/dq| = |a' + k| over a piece.
fn twist(law: &Law, piece: &Piece) -> f64 {
    match piece.profile {
        Profile::Constant(_) => law.k.abs(),
        Profile::Line { slope, .. } => (slope + law.k).abs(),
        Profile::Arc { .. } => piece.d1 - law.k.abs(),
    }
}

/// Helical chart lines: column `l` is sub-line `index[l]` of `piece[l]`.
struct Lines {
    piece: Vec<usize>,
    index: Vec<usize>,
    counts: Vec<usize>,
}
impl Lines {
    fn len(&self) -> usize {
        self.piece.len()
    }
    /// Phase of line `l` (0..=len; `len` is line 0 one pitch on) at `theta`.
    fn q(&self, law: &Law, l: usize, theta: f64) -> f64 {
        if l == self.len() {
            // The same physical line one turn on: (q + p, theta) is
            // (q, theta + 2 pi), and creases drift in q along a turn.
            return self.q(law, 0, theta + std::f64::consts::TAU) + law.pitch;
        }
        let piece = &law.pieces[self.piece[l]];
        let f = self.index[l] as f64 / self.counts[self.piece[l]] as f64;
        if let Some((qc, ra, phi0, phi1)) = piece.arc {
            if self.index[l] == 0 {
                return law.track_theta(piece.start, theta);
            }
            return qc + ra * (phi0 + (phi1 - phi0) * f).portable_sin();
        }
        let a = law.track_theta(piece.start, theta);
        let b = law.track_theta(piece.end, theta);
        a + (b - a) * f
    }
}
/// A `2 pi` periodic theta grid: `per` nodes in `[base, base + 2 pi)`.
struct Thetas {
    nodes: Vec<f64>,
}
impl Thetas {
    fn at(&self, k: i64) -> f64 {
        let n = self.nodes.len() as i64;
        let turn = k.div_euclid(n);
        self.nodes[k.rem_euclid(n) as usize] + std::f64::consts::TAU * turn as f64
    }
    fn per(&self) -> i64 {
        self.nodes.len() as i64
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
enum Key {
    Grid(usize, i64),
    /// Cut `c` of the meridian edge of strip `s` at theta index `k`.
    Meridian(usize, i64, u8),
    /// Cut `c` of the helical edge of line `l` from `k` to `k + 1`.
    Helical(usize, i64, u8),
}
#[derive(Clone, Copy)]
struct Corner {
    key: Key,
    q: f64,
    theta: f64,
}
struct Builder<'a> {
    thread: &'a AuthoredThread,
    law: &'a Law,
    lines: &'a Lines,
    thetas: &'a Thetas,
    ids: std::collections::HashMap<Key, u32>,
    vertices: Vec<f64>,
    /// Vertex z and phase, for cap rims and checks.
    heights: Vec<f64>,
    vertex_error: f64,
    snap: f64,
}
impl Builder<'_> {
    fn canonical(&self, key: Key) -> Key {
        let n = self.lines.len();
        let m = self.thetas.per();
        match key {
            // Line `n` at theta k is line 0 at theta k + one turn.
            Key::Grid(l, k) if l == n => Key::Grid(0, k + m),
            Key::Helical(l, k, c) if l == n => Key::Helical(0, k + m, c),
            other => other,
        }
    }
    fn vertex(&mut self, corner: Corner) -> Result<u32, Error> {
        let key = self.canonical(corner.key);
        if let Some(&id) = self.ids.get(&key) {
            return Ok(id);
        }
        let law = self.law;
        let mut z = corner.q + corner.theta / law.omega;
        // Cut and snapped vertices lie on their cap plane exactly.
        for cut in [0., law.length] {
            if (z - cut).abs() <= 2. * self.snap {
                z = cut;
            }
        }
        // Any adjacent piece gives the same radius on a shared line.
        let piece = self.piece_at(corner.q, z);
        let (radius, _, _) = law.radius(piece, corner.q, z);
        let p = law.pitch;
        let phase = corner.q.rem_euclid(p).clamp(0., p);
        let enclosure = self
            .thread
            .radial_bounds([phase, phase], [z.clamp(0., law.length); 2])?;
        let error = (radius - enclosure[0])
            .abs()
            .max((enclosure[1] - radius).abs());
        if !error.is_finite() {
            return Err(Error::Unresolved);
        }
        self.vertex_error = self.vertex_error.max(error);
        let (s, c) = corner.theta.portable_sin_cos();
        let id = (self.vertices.len() / 3) as u32;
        self.vertices.extend([radius * c, radius * s, z]);
        self.heights.push(z);
        self.ids.insert(key, id);
        Ok(id)
    }
    /// The piece containing phase `q` on its own chart line.
    fn piece_at(&self, q: f64, z: f64) -> &Piece {
        let law = self.law;
        let theta_of = |q: f64| law.omega * (z - q);
        for piece in &law.pieces {
            for shift in [-law.pitch, 0., law.pitch] {
                let qq = q + shift;
                let t = theta_of(qq);
                let a = law.track_theta(piece.start, t);
                let b = law.track_theta(piece.end, t);
                if qq >= a - 1e-12 && qq <= b + 1e-12 {
                    return piece;
                }
            }
        }
        &law.pieces[0]
    }
    fn grid(&self, l: usize, k: i64) -> Corner {
        let theta = self.thetas.at(k);
        let mut q = self.lines.q(self.law, l, theta);
        // Snap a grid vertex lying within `snap` of a cap onto it along its
        // meridian: still an exact law point, never a sliver cell.
        let z = q + theta / self.law.omega;
        for cut in [0., self.law.length] {
            if (z - cut).abs() <= self.snap {
                q = cut - theta / self.law.omega;
            }
        }
        Corner {
            key: Key::Grid(l, k),
            q,
            theta,
        }
    }
}

/// Clip a convex param polygon to `z >= cut` (`keep_above`) or `z <= cut`,
/// creating keyed cut vertices on its edges.
fn clip(
    law: &Law,
    polygon: &[(Corner, Option<Key>)],
    cut: f64,
    keep_above: bool,
    tag: u8,
    on: f64,
) -> Vec<(Corner, Option<Key>)> {
    let z = |c: &Corner| c.q + c.theta / law.omega;
    // A corner within `on` of the cut lies on it: kept on both sides and
    // never duplicated by a cut vertex beside it.
    let inside = |c: &Corner| {
        if keep_above {
            z(c) >= cut - on
        } else {
            z(c) <= cut + on
        }
    };
    let strict = |c: &Corner| {
        if keep_above {
            z(c) > cut + on
        } else {
            z(c) < cut - on
        }
    };
    let n = polygon.len();
    let mut out = Vec::with_capacity(n + 2);
    for i in 0..n {
        let (a, edge) = polygon[i];
        let (b, _) = polygon[(i + 1) % n];
        let crossing = (strict(&a) && !inside(&b)) || (!inside(&a) && strict(&b));
        if inside(&a) {
            out.push(polygon[i]);
        }
        if crossing {
            let (za, zb) = (z(&a), z(&b));
            let t = (cut - za) / (zb - za);
            let key = match edge.expect("cell edges carry keys") {
                Key::Meridian(s, k, _) => Key::Meridian(s, k, tag),
                Key::Helical(l, k, _) => Key::Helical(l, k, tag),
                Key::Grid(..) => unreachable!("edge keys are not grid keys"),
            };
            let (q, theta) = if a.theta == b.theta {
                // Meridian: exact phase of the cut at this theta.
                (cut - a.theta / law.omega, a.theta)
            } else {
                (a.q + (b.q - a.q) * t, a.theta + (b.theta - a.theta) * t)
            };
            out.push((Corner { key, q, theta }, Some(key)));
        }
    }
    // Each output edge keeps the key of the source edge it lies on; a cut
    // chord between two cut vertices needs none.
    let m = out.len();
    let mut keyed = Vec::with_capacity(m);
    for i in 0..m {
        let (c, _) = out[i];
        let (next, _) = out[(i + 1) % m];
        let edge = polygon
            .iter()
            .enumerate()
            .find_map(|(j, (pc, key))| {
                let (qc, _) = polygon[(j + 1) % n];
                let on = |x: &Corner| {
                    x.key == pc.key
                        || x.key == qc.key
                        || matches!((x.key, key), (Key::Meridian(s, k, _), Some(Key::Meridian(s2, k2, _))) if s == *s2 && k == *k2)
                        || matches!((x.key, key), (Key::Helical(l, k, _), Some(Key::Helical(l2, k2, _))) if l == *l2 && k == *k2)
                };
                (on(&c) && on(&next)).then_some(*key)
            })
            .flatten();
        keyed.push((c, edge));
    }
    keyed
}

struct Emitted {
    mesh: DisplayMesh,
    rims: [Vec<u32>; 2],
    deviation: f64,
    angle: f64,
}

#[allow(clippy::too_many_arguments)]
fn emit(
    thread: &AuthoredThread,
    law: &Law,
    alpha: f64,
    counts: &[usize],
    closure: LawClosure,
    triangle_limit: usize,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<Emitted, Error> {
    let w = law.omega;
    let p = law.pitch;
    let l = law.length;
    // Lines.
    let mut lines = Lines {
        piece: vec![],
        index: vec![],
        counts: counts.to_vec(),
    };
    for (piece, &count) in counts.iter().enumerate() {
        for i in 0..count {
            lines.piece.push(piece);
            lines.index.push(i);
        }
    }
    let n = lines.len();
    // Theta grid with the crease onset kinks as nodes.
    let per = (std::f64::consts::TAU / alpha).ceil().max(16.);
    if per > 1e6 {
        return Err(Error::Limit);
    }
    let per = per as usize;
    let mut nodes = (0..per)
        .map(|j| std::f64::consts::TAU * j as f64 / per as f64)
        .collect::<Vec<_>>();
    if let Some(onset) = law.onset {
        for piece in &law.pieces {
            for t in [piece.start, piece.end] {
                if let Track::Crease { base, .. } = t {
                    let kink = (w * (onset - base)).rem_euclid(std::f64::consts::TAU);
                    nodes.push(kink);
                }
            }
        }
    }
    nodes.sort_by(f64::total_cmp);
    nodes.dedup_by(|a, b| (*a - *b).abs() < 1e-12);
    let thetas = Thetas { nodes };
    let m = thetas.per();
    // Phase range of all lines, generously.
    let q_lo = law
        .pieces
        .iter()
        .flat_map(|pc| [pc.start, pc.end])
        .map(|t| law.track(t, 0.).min(law.track(t, l)))
        .fold(f64::INFINITY, f64::min)
        - 1e-9;
    let q_hi = q_lo + 2. * p;
    let k_lo = ((-w * q_hi) / std::f64::consts::TAU * m as f64).floor() as i64 - 2 * m;
    let k_hi = ((w * (l - q_lo)) / std::f64::consts::TAU * m as f64).ceil() as i64 + 2 * m;
    let estimate = 2 * n * ((k_hi - k_lo) as usize) / 2;
    if estimate > 4 * triangle_limit.max(1) {
        return Err(Error::Limit);
    }
    let mut builder = Builder {
        thread,
        law,
        lines: &lines,
        thetas: &thetas,
        ids: std::collections::HashMap::new(),
        vertices: Vec::new(),
        heights: Vec::new(),
        vertex_error: 0.,
        snap: 1e-9 * (l + p),
    };
    let mut indices = Vec::new();
    let mut faces = Vec::new();
    let mut deviation = 0f64;
    let mut worst = 0f64;
    let minor = law.minor_min;
    let cancelled = || interrupt.is_some_and(|f| f());
    for s in 0..n {
        if cancelled() {
            return Err(Error::Cancelled);
        }
        let piece = &law.pieces[lines.piece[s]];
        let b = twist(law, piece);
        for k in k_lo..k_hi {
            let corners = [
                builder.grid(s, k),
                builder.grid(s + 1, k),
                builder.grid(s + 1, k + 1),
                builder.grid(s, k + 1),
            ];
            let zs = corners.map(|c| c.q + c.theta / w);
            if zs.iter().all(|&z| z <= 0.) || zs.iter().all(|&z| z >= l) {
                continue;
            }
            let polygon = vec![
                (corners[0], Some(Key::Meridian(s, k, 0))),
                (corners[1], Some(Key::Helical(s + 1, k, 0))),
                (corners[2], Some(Key::Meridian(s, k + 1, 0))),
                (corners[3], Some(Key::Helical(s, k, 0))),
            ];
            let on = builder.snap;
            let mut polygon = clip(law, &polygon, 0., true, 1, on);
            if polygon.len() >= 3 {
                polygon = clip(law, &polygon, l, false, 2, on);
            }
            if polygon.len() < 3 {
                continue;
            }
            // The stock onset circle splits a clipped crest cell exactly.
            let parts = match law.onset {
                Some(onset) if piece.clipped && onset > 0. && onset < l => [
                    clip(law, &polygon, onset, false, 3, on),
                    clip(law, &polygon, onset, true, 3, on),
                ],
                _ => [polygon, vec![]],
            };
            for part in parts {
                if part.len() < 3 {
                    continue;
                }
                let zmid =
                    part.iter().map(|(c, _)| c.q + c.theta / w).sum::<f64>() / part.len() as f64;
                let cylinder = law.cylinder(piece, zmid);
                let role = if cylinder {
                    thread.band_count() as u32 + 1
                } else {
                    piece.role
                };
                let ids = part
                    .iter()
                    .map(|(c, _)| builder.vertex(*c))
                    .collect::<Result<Vec<_>, _>>()?;
                // Param extents of the cell part bound every fan triangle.
                let (mut dq, mut dt) = (0f64, 0f64);
                for a in &part {
                    for c in &part {
                        dq = dq.max((a.0.q - c.0.q).abs());
                        dt = dt.max((a.0.theta - c.0.theta).abs());
                    }
                }
                let (r_max, d1, d2, bb, turn_k, kk) = if cylinder {
                    (law.limit.expect("clip"), 0., 0., 0., 0., 0.)
                } else if let Some((qc, ra, _, _)) = piece.arc {
                    // Local arc bounds at the strip's largest |x|.
                    let x = part
                        .iter()
                        .map(|(c, _)| (c.q - qc).abs())
                        .fold(0., f64::max)
                        .min(ra * (1. - 1e-12));
                    let sx = (ra * ra - x * x).sqrt();
                    let slope = x / sx;
                    (
                        piece.r_max,
                        slope + law.k.abs(),
                        ra * ra / (sx * sx * sx),
                        slope,
                        1. / sx,
                        law.k.abs(),
                    )
                } else {
                    (piece.r_max, piece.d1, piece.d2, b, piece.turn, law.k.abs())
                };
                // |D2 S[d,d]| in (q, theta): |a''|dq^2 + 2|a'+k|dq dt + (A + 2|k|/w)dt^2.
                let chord =
                    0.5 * (d2 * dq * dq + 2. * bb * dq * dt + (r_max + 2. * kk / w) * dt * dt);
                deviation = deviation.max(chord);
                // Normal turn over the whole part: frame dt, profile
                // K dq, and the slow u = a'/(omega A) drift.
                let radial = bb * dq + kk * dt / w;
                let turn = dt
                    + turn_k * (1. + (d1 / (w * minor)).powi(2)) * dq
                    + d2 * dq / (w * minor)
                    + d1 * radial / (w * minor * minor);
                // Best fan of the convex part: maximise the worst triangle.
                let point = |id: u32| {
                    let i = 3 * id as usize;
                    [
                        builder.vertices[i],
                        builder.vertices[i + 1],
                        builder.vertices[i + 2],
                    ]
                };
                let exact = |c: &Corner| -> V3 {
                    let z = c.q + c.theta / w;
                    let (radius, rq, _) = law.radius(piece, c.q, z);
                    // At fixed (theta, z): R_theta = -a'/w, R_z = a' + k;
                    // the stock cylinder has neither.
                    let (radius, rt, rz) = if cylinder {
                        (law.limit.expect("clip"), 0., 0.)
                    } else {
                        (radius, -rq / w, rq + law.k)
                    };
                    let (sn, co) = c.theta.portable_sin_cos();
                    let e = [co, sn, 0.];
                    let t = [-sn, co, 0.];
                    let st: V3 = std::array::from_fn(|i| rt * e[i] + radius * t[i]);
                    let sz: V3 = std::array::from_fn(|i| rz * e[i] + [0., 0., 1.][i]);
                    let out = cross(st, sz);
                    let length = norm(out);
                    let out = out.map(|v| v / length);
                    if closure == LawClosure::Bore {
                        out.map(|v| -v)
                    } else {
                        out
                    }
                };
                let normals = part.iter().map(|(c, _)| exact(c)).collect::<Vec<_>>();
                let count = part.len();
                let mut best: Option<(f64, usize)> = None;
                for root in 0..count {
                    let mut fan_worst = 0f64;
                    let mut ok = true;
                    for i in 1..count - 1 {
                        let tri = [root, (root + i) % count, (root + i + 1) % count];
                        let pts = tri.map(|t| point(ids[t]));
                        let facet = cross(sub(pts[1], pts[0]), sub(pts[2], pts[0]));
                        if norm(facet).is_nan() || norm(facet) <= 0. {
                            ok = false;
                            break;
                        }
                        let mut a = 0f64;
                        // Exact normals at the vertices and edge midpoints:
                        // every facet point is within a third of the
                        // extents of one of these samples.
                        let mut samples = tri.map(|t| normals[t]).to_vec();
                        for e in 0..3 {
                            let (u, v) = (&part[tri[e]].0, &part[tri[(e + 1) % 3]].0);
                            samples.push(exact(&Corner {
                                key: u.key,
                                q: 0.5 * (u.q + v.q),
                                theta: 0.5 * (u.theta + v.theta),
                            }));
                        }
                        for n in samples {
                            // Orientation is fixed below; compare unsigned.
                            let g = angle(facet, n);
                            a = a.max(g.min(std::f64::consts::PI - g));
                        }
                        fan_worst = fan_worst.max(a);
                    }
                    if ok && best.is_none_or(|(b, _)| fan_worst < b) {
                        best = Some((fan_worst, root));
                    }
                }
                let Some((facet_angle, root)) = best else {
                    return Err(Error::Unresolved);
                };
                worst = worst.max(facet_angle + turn / 3.);
                for i in 1..count - 1 {
                    let mut tri = [root, (root + i) % count, (root + i + 1) % count];
                    let pts = tri.map(|t| point(ids[t]));
                    let facet = cross(sub(pts[1], pts[0]), sub(pts[2], pts[0]));
                    if dot(facet, normals[tri[0]]) < 0. {
                        tri.swap(1, 2);
                    }
                    indices.extend(tri.map(|t| ids[t]));
                    faces.push(role);
                }
            }
            if faces.len() > triangle_limit {
                return Err(Error::Limit);
            }
        }
    }
    let mut vertices = builder.vertices;
    let heights = builder.heights;
    let vertex_error = builder.vertex_error
        + 8. * f64::EPSILON * law.pieces.iter().map(|p| p.r_max).fold(0., f64::max);
    // Cap rims: every vertex on z = 0 or z = L, counterclockwise.
    let mut rims = [Vec::new(), Vec::new()];
    for (id, &z) in heights.iter().enumerate() {
        for (end, cut) in [0., l].into_iter().enumerate() {
            if z == cut {
                rims[end].push(id as u32);
            }
        }
    }
    for rim in &mut rims {
        let angle_of = |id: &u32| {
            let i = 3 * *id as usize;
            vertices[i + 1].portable_atan2(vertices[i])
        };
        rim.sort_by(|a, b| angle_of(a).total_cmp(&angle_of(b)));
        if rim.len() < 3 {
            return Err(Error::Unresolved);
        }
    }
    if closure == LawClosure::Solid {
        let band_count = thread.band_count() as u32;
        let centre = (vertices.len() / 3) as u32;
        vertices.extend([0., 0., 0., 0., 0., l]);
        for (end, rim) in rims.iter().enumerate() {
            let c = centre + end as u32;
            for i in 0..rim.len() {
                let a = rim[i];
                let b = rim[(i + 1) % rim.len()];
                if end == 0 {
                    indices.extend([c, b, a]);
                } else {
                    indices.extend([c, a, b]);
                }
                faces.push(band_count + 2 + end as u32);
            }
        }
        if faces.len() > triangle_limit {
            return Err(Error::Limit);
        }
    }
    // Exact capacities: callers account retained bytes by capacity.
    vertices.shrink_to_fit();
    indices.shrink_to_fit();
    faces.shrink_to_fit();
    Ok(Emitted {
        mesh: DisplayMesh {
            vertices,
            indices,
            face_indices: faces,
        },
        rims,
        deviation: (deviation + vertex_error).next_up(),
        angle: worst,
    })
}

/// Sub-line counts per piece for a theta spacing `alpha`, from the chord and
/// normal budgets; `None` when `alpha` alone exhausts a budget.
fn plan(law: &Law, alpha: f64, chord: f64, angle: f64, turn: f64) -> Option<Vec<usize>> {
    let w = law.omega;
    let minor = law.minor_min;
    let mut counts = Vec::with_capacity(law.pieces.len());
    for piece in &law.pieces {
        let b = twist(law, piece);
        let a = piece.r_max + 2. * law.k.abs() / w;
        let rest = 2. * chord - a * alpha * alpha;
        if rest <= 0. {
            return None;
        }
        // d2 x^2 + 2 b alpha x <= rest.
        let x = if piece.d2 > 0. {
            let bb = b * alpha;
            (-bb + (bb * bb + piece.d2 * rest).sqrt()) / piece.d2
        } else if b > 0. {
            rest / (2. * b * alpha)
        } else {
            f64::INFINITY
        };
        let turn_rate = piece.turn * (1. + (piece.d1 / (w * minor)).powi(2))
            + piece.d2 / (w * minor)
            + piece.d1 * b / (w * minor * minor);
        // Facet angle adds about half a sub-band turn to a third of it;
        // `turn` is the share of the angle a sub-band may use.
        let budget = turn * angle - alpha;
        if budget <= 0. {
            return None;
        }
        let y = if turn_rate > 0. {
            budget / turn_rate
        } else {
            f64::INFINITY
        };
        let mut width = 0f64;
        for z in [0., law.onset.unwrap_or(0.), law.length] {
            width = width.max(law.track(piece.end, z) - law.track(piece.start, z));
        }
        let count = if let Some((_, ra, phi0, phi1)) = piece.arc {
            // Uniform in phi: normal turn dphi per sub-band, chord
            // ra dphi^2 / cos(phi) + 2 ra dphi alpha + A alpha^2.
            let cos = (phi0.abs().max(phi1.abs())).portable_cos().max(1e-6);
            let c2 = ra / cos;
            let c1 = 2. * ra * alpha;
            let by_chord = (-c1 + (c1 * c1 + 4. * c2 * rest).sqrt()) / (2. * c2);
            let dphi = by_chord.min(budget);
            ((phi1 - phi0).abs() / dphi).ceil().max(1.)
        } else {
            (width / x.min(y)).ceil().max(1.)
        };
        if !count.is_finite() || count > 1e5 {
            return None;
        }
        counts.push(count as usize);
    }
    Some(counts)
}

impl AuthoredThread {
    /// Tessellate the authored law within `tessellation.deflection` and
    /// `tessellation.angle`. The crease limit is honoured for a solid.
    pub fn tessellate_law(
        &self,
        tessellation: Tessellation,
        closure: LawClosure,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<LawMesh, Error> {
        tessellation.check().map_err(|_| Error::Invalid)?;
        let law = build_law(self, closure)?;
        let mut chord = 0.92 * tessellation.deflection;
        let mut turn = 1.;
        let angle = tessellation.angle;
        let r_max = law.pieces.iter().map(|p| p.r_max).fold(0., f64::max);
        // Theta spacing ceiling: the chord of the widest circle and the
        // frame turn; lowered when a facet angle bound fails, because the
        // tilt of cap-cut facets scales with alpha^2.
        let mut alpha_limit = (2. * chord / r_max).sqrt().min(0.15 * angle);
        for _ in 0..10 {
            let mut best: Option<(f64, f64, Vec<usize>)> = None;
            let mut alpha = alpha_limit * 0.98;
            for _ in 0..24 {
                if let Some(counts) = plan(&law, alpha, chord, angle, turn) {
                    let cost = counts.iter().sum::<usize>() as f64 / alpha;
                    if best.as_ref().is_none_or(|(c, _, _)| cost < *c) {
                        best = Some((cost, alpha, counts));
                    }
                }
                alpha *= 0.85;
            }
            let (_, alpha, counts) = best.ok_or(Error::Limit)?;
            let emitted = emit(
                self,
                &law,
                alpha,
                &counts,
                closure,
                tessellation.triangle_limit,
                interrupt,
            )?;
            if emitted.deviation <= tessellation.deflection && emitted.angle <= angle {
                let proof = if closure == LawClosure::Solid {
                    Some(manufacturing::certify(&emitted.mesh).map_err(|_| Error::Unresolved)?)
                } else {
                    None
                };
                return Ok(LawMesh {
                    rows: (std::f64::consts::TAU / alpha).ceil() as usize,
                    columns: counts.iter().sum(),
                    mesh: emitted.mesh,
                    band_count: self.band_count(),
                    rims: emitted.rims,
                    deviation_mm: emitted.deviation,
                    normal_angle_rad: emitted.angle,
                    proof,
                });
            }
            if emitted.deviation > tessellation.deflection {
                chord *= 0.8;
            }
            if emitted.angle > angle {
                // Both the frame step and the sub-band turn shrink.
                alpha_limit = alpha * 0.8;
                turn *= 0.85;
            }
        }
        Err(Error::Limit)
    }
}
#[cfg(test)]
mod tests;
