//! Spherical bevel gears as exact scene members.
//!
//! In the gear's frame the solid is the set of points `r u` whose direction
//! `u` lies in the spherical tooth domain and whose distance lies between
//! the toe and heel back cones, less the bore. Every face is ruled:
//!
//! * side faces are cones through the pitch apex over the pieces of the
//!   domain boundary: circular cones of the root and tip, circular cones of
//!   the root fillets, meridian planes and spherical-involute cones;
//! * both caps are circular back cones about the axis;
//! * the bore is a cylinder.
//!
//! The surface normal is constant along each ruling, so outlines are whole
//! rulings, and every line or plane meets a ruled face where a function of
//! the ruling parameter alone changes sign. Those functions are evaluated
//! from the provider's exact boundary law and their sign changes isolated
//! with certified derivative bounds; no mesh or polygon stands in.
use super::paths::{Lift, ProfilePath, Segment, add, dot3, scaled};
use super::*;
use multicad_rust_kernel::{
    bevel_gear::{Bevel, Containment, End, Piece},
    normal_sweep,
    spur_gear::Edge,
};
use std::f64::consts::{PI, TAU};

mod roots;
mod section;
#[cfg(test)]
mod tests;

fn cross(a: P3, b: P3) -> P3 {
    std::array::from_fn(|i| a[(i + 1) % 3] * b[(i + 2) % 3] - a[(i + 2) % 3] * b[(i + 1) % 3])
}

fn norm(a: P3) -> f64 {
    dot3(a, a).sqrt()
}

fn sub(a: P3, b: P3) -> P3 {
    std::array::from_fn(|i| a[i] - b[i])
}

fn polar(gamma: f64, phi: f64) -> P3 {
    let (sg, cg) = gamma.sin_cos();
    let (sp, cp) = phi.sin_cos();
    [sg * cp, sg * sp, cg]
}

fn end_index(end: End) -> usize {
    usize::from(end == End::Toe)
}

fn flank(piece: Piece) -> bool {
    !matches!(piece, Piece::RootLeft | Piece::RootRight | Piece::Tip)
}

fn involute(piece: Piece) -> bool {
    matches!(piece, Piece::InvoluteLeft | Piece::InvoluteRight)
}

fn meridian(piece: Piece) -> bool {
    matches!(piece, Piece::MeridianLeft | Piece::MeridianRight)
}

/// Shape of one boundary piece on the unit sphere, in its own parameter.
#[derive(Clone, Copy, Debug)]
struct Shape {
    /// Bound of `|u'|`.
    speed: f64,
    /// Bound of `|u''|`.
    curvature: f64,
    /// Turning rate: the angle swept about a circle's axis, or the involute
    /// roll interval.
    turn: f64,
    /// Axis of a circular piece in sector 0 and the cosine of its radius.
    axis: P3,
    cosine: f64,
}

/// A placed bevel gear in one view.
pub(super) struct Core {
    bevel: Bevel,
    origin: P3,
    axes: [P3; 3],
    /// Unit direction toward the viewer in the gear frame.
    sight: P3,
    teeth: u32,
    pieces: Vec<Piece>,
    /// Directions of the piece joints of sector 0.
    joints: Vec<P3>,
    shapes: Vec<Shape>,
    /// Sine and cosine of the pitch cone angle.
    pitch: [f64; 2],
    root: f64,
    tip: f64,
    /// Heel and toe distances along the pitch generator.
    ends: [f64; 2],
    bore: Option<f64>,
    /// Half azimuth of a tooth at the root cone and at the tip cone.
    foot: f64,
    crest: f64,
    /// Largest distance of the solid from the pitch apex.
    scale: f64,
    bounds: [P3; 2],
}

impl Core {
    fn new(bevel: Bevel, frame: PlaneFrame, view: &View) -> Self {
        let axes = [frame.u, frame.v, frame.normal()].map(|axis| view.vector(axis));
        let sight = [axes[0][2], axes[1][2], axes[2][2]];
        let sight = scaled(sight, 1. / norm(sight));
        let spec = bevel.spec();
        let pieces = bevel.pieces().to_vec();
        let joints = pieces
            .iter()
            .map(|&piece| bevel.sample(piece, 0.).direction)
            .collect();
        let base_cos = spec.base_cone_angle_rad.cos();
        let shapes = pieces
            .iter()
            .map(|&piece| {
                if involute(piece) {
                    // u' = c sin(b t) Δt T with |T| = 1 and c = cos β,
                    // b = sin β; the roll t is largest at the tip, where
                    // cos(b t) = cos(tip) / c.
                    let tip_end = if piece == Piece::InvoluteRight {
                        0.
                    } else {
                        1.
                    };
                    let speed = norm(bevel.sample(piece, tip_end).derivative);
                    let rise = (1. - (spec.tip_cone_angle_rad.cos() / base_cos).powi(2)).sqrt();
                    let turn = speed / (base_cos * rise);
                    Shape {
                        speed,
                        // |d²u/dt²| <= 1 before the reparametrisation.
                        curvature: turn * turn,
                        turn,
                        axis: [0., 0., 1.],
                        cosine: 0.,
                    }
                } else {
                    // A circle on the sphere: its tangent turns about the
                    // axis by the swept angle.
                    let [a, b] = [0., 0.5].map(|s| bevel.sample(piece, s));
                    let speed = norm(a.derivative);
                    let normal = cross(a.derivative, b.derivative);
                    let turn = 2. * norm(normal).atan2(dot3(a.derivative, b.derivative));
                    let mut axis = scaled(normal, 1. / norm(normal));
                    if dot3(axis, a.direction) < 0. {
                        axis = scaled(axis, -1.);
                    }
                    Shape {
                        speed,
                        curvature: speed * turn,
                        turn,
                        axis,
                        cosine: dot3(axis, a.direction),
                    }
                }
            })
            .collect();
        let origin = view.point(frame.origin);
        let mut bounds = [origin; 2];
        for i in 0..3 {
            let direction = [axes[0][i], axes[1][i], axes[2][i]];
            bounds[1][i] += bevel.support(direction);
            bounds[0][i] -= bevel.support(scaled(direction, -1.));
        }
        let pitch = spec.pitch_cone_angle_rad;
        Self {
            origin,
            axes,
            sight,
            teeth: spec.teeth,
            joints,
            shapes,
            pitch: [pitch.sin(), pitch.cos()],
            root: spec.root_cone_angle_rad,
            tip: spec.tip_cone_angle_rad,
            ends: End::BOTH.map(|end| bevel.end_distance_mm(end)),
            bore: bevel.bore_radius_mm(),
            foot: bevel.half_width(spec.root_cone_angle_rad),
            crest: bevel.sample(Piece::Tip, 1.).phi,
            scale: spec.outer_cone_distance_mm / (spec.tip_cone_angle_rad - pitch).cos(),
            bounds,
            pieces,
            bevel,
        }
    }

    fn to_view(&self, p: P3) -> P3 {
        std::array::from_fn(|i| {
            self.origin[i]
                + p[0] * self.axes[0][i]
                + p[1] * self.axes[1][i]
                + p[2] * self.axes[2][i]
        })
    }

    fn vector(&self, v: P3) -> P3 {
        std::array::from_fn(|i| {
            v[0] * self.axes[0][i] + v[1] * self.axes[1][i] + v[2] * self.axes[2][i]
        })
    }

    fn to_local(&self, q: P3) -> P3 {
        let offset = sub(q, self.origin);
        self.axes.map(|axis| dot3(axis, offset))
    }

    fn spacing(&self) -> f64 {
        TAU / f64::from(self.teeth)
    }

    /// Rotate a sector-0 vector into a sector.
    fn turn(&self, sector: u32, p: P3) -> P3 {
        let (s, c) = (self.spacing() * f64::from(sector)).sin_cos();
        [c * p[0] - s * p[1], s * p[0] + c * p[1], p[2]]
    }

    /// The one direction of a joint between two pieces; the last joint of a
    /// sector is the first of the next.
    fn joint(&self, sector: u32, joint: usize) -> P3 {
        if joint == self.pieces.len() {
            self.turn((sector + 1) % self.teeth, self.joints[0])
        } else {
            self.turn(sector, self.joints[joint])
        }
    }

    /// Direction of a boundary piece and its parameter derivative. Both ends
    /// evaluate to the shared joint directions.
    fn direction(&self, sector: u32, piece: usize, s: f64) -> (P3, P3) {
        let sample = self.bevel.sample(self.pieces[piece], s);
        let u = if s == 0. {
            self.joint(sector, piece)
        } else if s == 1. {
            self.joint(sector, piece + 1)
        } else {
            self.turn(sector, sample.direction)
        };
        (u, self.turn(sector, sample.derivative))
    }

    /// Distance from the apex to a back cone along a unit direction.
    fn reach(&self, end: End, u: P3) -> f64 {
        self.ends[end_index(end)] / (self.pitch[0] * u[0].hypot(u[1]) + self.pitch[1] * u[2])
    }

    /// Point of a cap edge and its parameter derivative.
    fn rim(&self, end: End, sector: u32, piece: usize, s: f64) -> (P3, P3) {
        let (u, du) = self.direction(sector, piece, s);
        let radial = u[0].hypot(u[1]);
        let cosine = self.pitch[0] * radial + self.pitch[1] * u[2];
        let rate = self.pitch[0] * (u[0] * du[0] + u[1] * du[1]) / radial + self.pitch[1] * du[2];
        let r = self.ends[end_index(end)] / cosine;
        (
            scaled(u, r),
            add(scaled(du, r), scaled(u, -r * rate / cosine)),
        )
    }

    /// Bound of the second parameter derivative of a cap edge.
    ///
    /// With `p = r(γ) u`, `r = E / cos(γ − δ)`, `|u'| <= S`, `|u''| <= C`
    /// and `|γ'| <= |u'|`, differentiating `cos γ = u_z` twice gives
    /// `|γ''| <= (C + S²) / sin γ`; then
    /// `p'' = (r''γ'² + r'γ'')u + 2r'γ'u' + r u''` with `sin γ >= sin root`.
    fn rim_bound(&self, end: End, piece: usize) -> f64 {
        let delta = self.pitch[0].atan2(self.pitch[1]);
        let m = (delta - self.root).abs().max((self.tip - delta).abs());
        let (sm, cm) = m.sin_cos();
        let k = self.ends[end_index(end)];
        let r0 = k / cm;
        let r1 = k * sm / (cm * cm);
        let r2 = k * (1. + sm * sm) / cm.powi(3);
        let Shape {
            speed, curvature, ..
        } = self.shapes[piece];
        let gamma2 = (curvature + speed * speed) / self.root.sin();
        r2 * speed * speed + r1 * gamma2 + 2. * r1 * speed * speed + r0 * curvature
    }

    /// Azimuth of the tooth center nearest to an azimuth, and the offset.
    fn tooth(&self, phi: f64) -> (f64, f64) {
        let center = (phi / self.spacing()).round() * self.spacing();
        (center, phi - center)
    }

    /// Largest polar angle of the tooth domain at an azimuth offset from a
    /// tooth center. The half tooth width never increases with γ.
    fn top(&self, offset: f64) -> f64 {
        let offset = offset.abs();
        if offset >= self.foot {
            return self.root;
        }
        if offset <= self.crest {
            return self.tip;
        }
        let (mut lo, mut hi) = (self.root, self.tip);
        for _ in 0..64 {
            let mid = 0.5 * (lo + hi);
            if self.bevel.half_width(mid) >= offset {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    }

    /// Bottom of a cap generator: on the bore or at the back-cone apex.
    fn bottom(&self, end: End, phi: f64) -> P3 {
        match self.bore {
            Some(bore) => [
                bore * phi.cos(),
                bore * phi.sin(),
                self.bevel.back_cone_height_mm(end, bore),
            ],
            None => [0., 0., self.bevel.back_cone_apex_mm(end)],
        }
    }

    fn inside(&self, local: P3) -> bool {
        self.bevel.classify(local, 0.) == Containment::Inside
    }

    /// Meridian flank planes of the gear: sector, piece, unit radial
    /// direction of the half-plane and the polar band of the face, ascending.
    fn meridians(&self) -> Vec<(u32, usize, P3, [f64; 2])> {
        let mut output = Vec::new();
        for sector in 0..self.teeth {
            for (piece, _) in self
                .pieces
                .iter()
                .enumerate()
                .filter(|(_, p)| meridian(**p))
            {
                let [a, b] = [0., 1.].map(|s| self.direction(sector, piece, s).0);
                let length = a[0].hypot(a[1]);
                let radial = [a[0] / length, a[1] / length, 0.];
                let [ga, gb] = [a, b].map(|u| u[2].clamp(-1., 1.).acos());
                output.push((sector, piece, radial, [ga.min(gb), ga.max(gb)]));
            }
        }
        output
    }

    /// Whether a point lies on a meridian face whose plane contains the
    /// line of sight: such a point is on the boundary, never strictly
    /// inside, whatever the rounding of its azimuth.
    fn on_face(planes: &[(u32, usize, P3, [f64; 2])], local: P3) -> bool {
        planes.iter().any(|&(_, _, radial, [lo, hi])| {
            let length = norm(local);
            let gamma = (local[2] / length).clamp(-1., 1.).acos();
            dot3(local, radial) > 0. && lo <= gamma && gamma <= hi
        })
    }

    /// Outline function `det[d, u, T]` of a fillet or involute cone with its
    /// derivative, `T` the unit tangent of the piece.
    fn limb(&self, sector: u32, piece: usize, s: f64) -> (f64, f64) {
        let d = self.sight;
        let (u, du) = self.direction(sector, piece, s);
        let shape = self.shapes[piece];
        let tangent = scaled(du, 1. / norm(du));
        let value = dot3(d, cross(u, tangent));
        let rate = if involute(self.pieces[piece]) {
            // dT/ds = Δt (z × T): the tangent of the base cone's rolling
            // great circle turns about the axis.
            shape.turn * dot3(d, cross(u, cross([0., 0., 1.], tangent)))
        } else {
            // u'' = −ω² (u − cos ρ c) on a circle about c.
            let axis = self.turn(sector, shape.axis);
            shape.turn * shape.turn * shape.cosine * dot3(d, cross(u, axis)) / shape.speed
        };
        (value, rate)
    }

    /// The value of the outline function at a joint of the smooth flank
    /// chain, shared by both neighbors. A meridian plane has one normal.
    fn limb_joint(&self, sector: u32, joint: usize) -> f64 {
        let [before, after] = [joint - 1, joint];
        for piece in [before, after] {
            if meridian(self.pieces[piece]) {
                return self.limb(sector, piece, 0.5).0;
            }
        }
        for (piece, s) in [(before, 1.), (after, 0.)] {
            if flank(self.pieces[piece]) && !involute(self.pieces[piece]) {
                return self.limb(sector, piece, s).0;
            }
        }
        if involute(self.pieces[before]) {
            self.limb(sector, before, 1.).0
        } else {
            self.limb(sector, after, 0.).0
        }
    }
}

fn unit_circle() -> Edge {
    Edge::Planar(normal_sweep::Edge::Circle {
        center: [0.; 2],
        radius_mm: 1.,
        clockwise: false,
    })
}

/// An arc `origin + cos(a) x + sin(a) y` over angles, split where its paper
/// image turns back when the arc is seen edge-on, and into at most half
/// turns so that no path closes on itself.
fn arcs(lift: Lift, [from, to]: [f64; 2]) -> Vec<Arc<dyn Path>> {
    let [x, y] = [[lift.x[0], lift.x[1]], [lift.y[0], lift.y[1]]];
    let scale = (x[0].hypot(x[1])).max(y[0].hypot(y[1]));
    let straight = (x[0] * y[1] - x[1] * y[0]).abs() <= 1e-12 * scale * scale;
    let mut stops = vec![from];
    let first = if straight {
        let d = if x[0].hypot(x[1]) >= y[0].hypot(y[1]) {
            x
        } else {
            y
        };
        (y[0] * d[0] + y[1] * d[1]).atan2(x[0] * d[0] + x[1] * d[1])
    } else {
        from
    };
    let mut k = ((from - first) / PI).ceil();
    while first + k * PI < to {
        if first + k * PI > from {
            stops.push(first + k * PI);
        }
        k += 1.;
    }
    stops.push(to);
    stops
        .windows(2)
        .filter(|pair| pair[1] - pair[0] > 1e-12)
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

/// A fillet or involute cap edge.
struct Rim {
    core: Arc<Core>,
    end: End,
    sector: u32,
    piece: usize,
}

impl Path for Rim {
    fn point(&self, t: f64) -> P3 {
        self.core
            .to_view(self.core.rim(self.end, self.sector, self.piece, t).0)
    }

    fn flatten(&self, tolerance_mm: f64) -> Vec<f64> {
        // A chord of parameter width h stays within h² |p''| / 8.
        let second = self.core.rim_bound(self.end, self.piece);
        let count = (second / (8. * tolerance_mm)).sqrt().ceil().clamp(1., 1e6) as usize;
        (0..=count).map(|i| i as f64 / count as f64).collect()
    }

    fn paper(&self, _: [f64; 2]) -> Option<Curve> {
        None
    }
}

pub(super) struct Gear {
    core: Arc<Core>,
}

pub(in crate::drawing) fn member(
    solid: &Solid,
    view: &View,
) -> Option<Result<Box<dyn Member>, GeometryError>> {
    let (body, frame) = solid.bevel_body()?;
    Some(Ok(Box::new(Gear {
        core: Arc::new(Core::new(body.geometry.clone(), frame, view)),
    })))
}

impl Gear {
    /// Azimuths where `radius cos(φ − heading) = level`.
    fn azimuths(radius: f64, heading: f64, level: f64) -> Vec<f64> {
        if radius <= level.abs() || radius <= 1e-12 {
            return Vec::new();
        }
        let half = (level / radius).acos();
        vec![heading - half, heading + half]
    }
}

impl Member for Gear {
    fn bounds(&self) -> [P3; 2] {
        self.core.bounds
    }

    fn candidates(&self, _: f64) -> Result<Vec<Candidate>, GeometryError> {
        let core = &self.core;
        let count = core.pieces.len();
        let mut output = Vec::new();
        let mut push = |kind: Kind, path: Arc<dyn Path>| output.push(Candidate { kind, path });
        let line =
            |a: P3, b: P3| -> Arc<dyn Path> { Arc::new(Segment(core.to_view(a), core.to_view(b))) };
        // A circle about the axis through the back-cone point at polar γ.
        let ring = |end: End, gamma: f64, range: [f64; 2]| {
            let r = core.reach(end, polar(gamma, 0.));
            let radius = r * gamma.sin();
            arcs(
                Lift {
                    origin: core.to_view([0., 0., r * gamma.cos()]),
                    x: core.vector([radius, 0., 0.]),
                    y: core.vector([0., radius, 0.]),
                },
                range,
            )
        };
        let ruling = |u: P3| {
            line(
                scaled(u, core.reach(End::Toe, u)),
                scaled(u, core.reach(End::Heel, u)),
            )
        };
        let d = core.sight;
        let heading = d[1].atan2(d[0]);
        let across = d[0].hypot(d[1]);
        for sector in 0..core.teeth {
            let center = core.spacing() * f64::from(sector);
            for end in End::BOTH {
                for (index, &piece) in core.pieces.iter().enumerate() {
                    match piece {
                        // Both root lands between two teeth are one arc.
                        Piece::RootLeft => {}
                        Piece::RootRight => {
                            let range = [center + core.foot, center + core.spacing() - core.foot];
                            for path in ring(end, core.root, range) {
                                push(Kind::Edge, path);
                            }
                        }
                        Piece::Tip => {
                            for path in
                                ring(end, core.tip, [center - core.crest, center + core.crest])
                            {
                                push(Kind::Edge, path);
                            }
                        }
                        Piece::MeridianLeft | Piece::MeridianRight => {
                            let [a, b] = [0., 1.].map(|s| core.rim(end, sector, index, s).0);
                            push(Kind::Edge, line(a, b));
                        }
                        _ => push(
                            Kind::Edge,
                            Arc::new(Rim {
                                core: core.clone(),
                                end,
                                sector,
                                piece: index,
                            }),
                        ),
                    }
                }
            }
            // Apex rays between side faces: sharp only at the tip.
            for joint in 1..count {
                let sharp =
                    core.pieces[joint - 1] == Piece::Tip || core.pieces[joint] == Piece::Tip;
                push(
                    if sharp { Kind::Edge } else { Kind::Smooth },
                    ruling(core.joint(sector, joint)),
                );
            }
            // Outline rulings of the fillet and involute cones.
            for index in (0..count)
                .filter(|&index| flank(core.pieces[index]) && !meridian(core.pieces[index]))
            {
                let shape = core.shapes[index];
                let ends = [index, index + 1].map(|joint| core.limb_joint(sector, joint));
                let mut found = Vec::new();
                roots::crossings(
                    &|s| core.limb(sector, index, s),
                    2. * shape.turn * shape.turn,
                    [0., 1.],
                    ends,
                    &mut found,
                );
                for s in found {
                    push(Kind::Outline, ruling(core.direction(sector, index, s).0));
                }
            }
        }
        // Outline generators of the root and tip cones: the normal
        // (cos γ cos φ, cos γ sin φ, −sin γ) is perpendicular to the sight.
        for (gamma, tip) in [(core.root, false), (core.tip, true)] {
            for phi in Self::azimuths(across, heading, gamma.tan() * d[2]) {
                let offset = core.tooth(phi).1.abs();
                if (tip && offset <= core.crest) || (!tip && offset >= core.foot) {
                    push(Kind::Outline, ruling(polar(gamma, phi)));
                }
            }
        }
        // Back cones: the normal is the pitch direction at the azimuth.
        let [sin, cos] = core.pitch;
        for end in End::BOTH {
            for phi in Self::azimuths(across * sin, heading, -cos * d[2]) {
                let gamma = core.top(core.tooth(phi).1);
                let u = polar(gamma, phi);
                push(
                    Kind::Outline,
                    line(core.bottom(end, phi), scaled(u, core.reach(end, u))),
                );
            }
        }
        if let Some(bore) = core.bore {
            for end in End::BOTH {
                let z = core.bevel.back_cone_height_mm(end, bore);
                let lift = Lift {
                    origin: core.to_view([0., 0., z]),
                    x: core.vector([bore, 0., 0.]),
                    y: core.vector([0., bore, 0.]),
                };
                for path in arcs(lift, [0., TAU]) {
                    push(Kind::Edge, path);
                }
            }
            for phi in Self::azimuths(across, heading, 0.) {
                push(
                    Kind::Outline,
                    line(core.bottom(End::Toe, phi), core.bottom(End::Heel, phi)),
                );
            }
        }
        Ok(output)
    }

    fn sight(&self, [x, y]: P2) -> Result<Vec<[f64; 2]>, GeometryError> {
        let core = &self.core;
        let [lo, hi] = core.bounds;
        if x < lo[0] || x > hi[0] || y < lo[1] || y > hi[1] {
            return Ok(Vec::new());
        }
        let start = core.to_local([x, y, 0.]);
        let d = core.sight;
        let mut stops = Vec::new();
        // Root and tip cones: p_z² = cos²γ |p|².
        for gamma in [core.root, core.tip] {
            let c = gamma.cos().powi(2);
            roots::quadratic(
                d[2] * d[2] - c,
                start[2] * d[2] - c * dot3(start, d),
                start[2] * start[2] - c * dot3(start, start),
                &mut stops,
            );
        }
        // Back cones: sin²δ ρ² = (E − cos δ z)².
        let [sin, cos] = core.pitch;
        let flat = [start[0] * d[0] + start[1] * d[1], d[0] * d[0] + d[1] * d[1]];
        let radial = start[0] * start[0] + start[1] * start[1];
        for end in core.ends {
            let k = end - cos * start[2];
            roots::quadratic(
                sin * sin * flat[1] - cos * cos * d[2] * d[2],
                sin * sin * flat[0] + k * cos * d[2],
                sin * sin * radial - k * k,
                &mut stops,
            );
        }
        if let Some(bore) = core.bore {
            roots::quadratic(flat[1], flat[0], radial - bore * bore, &mut stops);
        }
        // A line lying in a meridian plane meets the meridian face along a
        // whole stretch: its plane function vanishes identically, so the
        // face is decided by exact case analysis instead of root finding.
        // The stretch ends where the line crosses the cones of the face's
        // polar band.
        let band = 1e-12 * core.scale;
        let planes: Vec<_> = core
            .meridians()
            .into_iter()
            .filter(|&(_, _, radial, _)| {
                let normal = [-radial[1], radial[0], 0.];
                dot3(normal, d).abs() <= 1e-12 && dot3(normal, start).abs() <= band
            })
            .collect();
        for &(_, _, _, gammas) in &planes {
            for gamma in gammas {
                let c = gamma.cos().powi(2);
                roots::quadratic(
                    d[2] * d[2] - c,
                    start[2] * d[2] - c * dot3(start, d),
                    start[2] * start[2] - c * dot3(start, start),
                    &mut stops,
                );
            }
        }
        // Flank cones: the ruling u lies in the plane through the apex and
        // the line where m · u = 0 with m = start × d.
        let moment = cross(start, d);
        let length = norm(moment);
        if length > 1e-12 * core.scale {
            let moment = scaled(moment, 1. / length);
            let middle = 0.5 * (core.root + core.tip);
            // Every flank direction of a tooth lies within this angle of
            // its center direction at the middle polar angle.
            let cap = 0.5 * (core.tip - core.root) + core.foot * core.tip.sin() + 1e-9;
            for sector in 0..core.teeth {
                let center = polar(middle, core.spacing() * f64::from(sector));
                if cap < 0.5 * PI && dot3(moment, center).abs() > cap.sin() {
                    continue;
                }
                for (index, _) in core
                    .pieces
                    .iter()
                    .enumerate()
                    .filter(|(_, piece)| flank(**piece))
                    // The fillet and involute next to a meridian are tangent
                    // to its plane along the joint ruling and lie on one
                    // side of it, so they meet the plane only there: at the
                    // band cones already among the stops.
                    .filter(|(index, _)| {
                        !planes
                            .iter()
                            .any(|&(k, piece, ..)| k == sector && index.abs_diff(piece) <= 1)
                    })
                {
                    let eval = |s: f64| {
                        let (u, du) = core.direction(sector, index, s);
                        (dot3(moment, u), dot3(moment, du))
                    };
                    let mut found = Vec::new();
                    roots::crossings(
                        &eval,
                        core.shapes[index].curvature,
                        [0., 1.],
                        [eval(0.).0, eval(1.).0],
                        &mut found,
                    );
                    for s in found {
                        let u = core.direction(sector, index, s).0;
                        let turn = cross(d, u);
                        let size = dot3(turn, turn);
                        if size > 1e-24 {
                            stops.push(-dot3(cross(start, u), turn) / size);
                        }
                    }
                }
            }
        }
        stops.retain(|z| z.is_finite());
        stops.sort_by(f64::total_cmp);
        stops.dedup_by(|b, a| *b - *a <= 1e-11 * core.scale);
        let mut output: Vec<[f64; 2]> = Vec::new();
        for pair in stops.windows(2) {
            // A probe can land exactly on a point where the line only
            // grazes a surface inside the interval; two of three decide.
            let inside = [0.25, 0.5, 0.75]
                .into_iter()
                .filter(|f| {
                    let local = add(start, scaled(d, pair[0] + (pair[1] - pair[0]) * f));
                    !Core::on_face(&planes, local) && core.inside(local)
                })
                .count();
            if inside >= 2 {
                match output.last_mut() {
                    Some(last) if last[1] == pair[0] => last[1] = pair[1],
                    _ => output.push([pair[0], pair[1]]),
                }
            }
        }
        Ok(output)
    }

    fn section(&self, depth: f64, tolerance_mm: f64) -> Result<Vec<Loop>, GeometryError> {
        section::section(&self.core, depth, tolerance_mm)
    }
}
