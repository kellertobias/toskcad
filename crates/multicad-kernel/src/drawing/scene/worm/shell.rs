//! The lateral boundary of a ZA worm or a generated worm wheel as a closed
//! chain of exact parametric patches.
//!
//! Both solids are slabs `|z| < half` about the local axis whose section at
//! every height is star shaped about that axis: `bore < r < R(θ, z)`. Each
//! patch covers `(u, z)` in `[0, 1] × [-half, half]`, keeps the height as its
//! second parameter and runs counterclockwise in `u`, so the chain of patches
//! at one height is the section boundary and `∂S/∂u × ∂S/∂z` points out of
//! the material.
use super::V3;
use multicad_rust_kernel::worm_gear::{
    Wheel, Worm,
    wheel::{PIECES, Piece},
};
use std::f64::consts::{PI, TAU};
use std::sync::Mutex;

/// A patch point with its partial derivatives, in local millimetres.
#[derive(Clone, Copy, Debug)]
pub(super) struct Jet3 {
    pub p: V3,
    pub du: V3,
    pub dz: V3,
}

impl Jet3 {
    /// Outward normal, not normalized.
    pub(super) fn normal(&self) -> V3 {
        super::cross(self.du, self.dz)
    }
}

/// How a patch meets the previous patch of the chain along `u = 0`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Joint {
    Edge,
    Smooth,
    /// Both patches lie on one surface.
    Seam,
}

pub(super) trait Shell: Send + Sync {
    fn half(&self) -> f64;
    fn bore(&self) -> Option<f64>;
    fn patches(&self) -> usize;
    fn eval(&self, patch: usize, u: f64, z: f64) -> Jet3;
    fn joint(&self, patch: usize) -> Joint;
    /// Cells along `u` and `z` that resolve implicit curves on a patch.
    fn grid(&self, patch: usize) -> [usize; 2];
    /// Radius, start angle and sweep when the patch section at height `z` is
    /// a circular arc about the axis traversed uniformly in `u`.
    fn arc(&self, patch: usize, z: f64) -> Option<[f64; 3]>;
    /// The section boundary radius `R(θ, z)`.
    fn radius(&self, theta: f64, z: f64) -> f64;
    /// Radii below which every section is material down to the bore, and
    /// above which none is.
    fn annulus(&self) -> [f64; 2];
    /// Bounds of `|∂R/∂θ|` and `|∂R/∂z|`.
    fn lipschitz(&self) -> [f64; 2];
    /// Maximum of `direction · p` over the solid.
    fn support(&self, direction: V3) -> f64;

    /// Heights in `[from, to]` strictly inside the solid along the axial
    /// line at angle `theta` on the core cylinder, when the boundary meets
    /// that line along whole stretches (the worm's root land). `None` when
    /// it only touches the line at isolated heights.
    fn core_line(&self, _theta: f64, _range: [f64; 2]) -> Option<Vec<[f64; 2]>> {
        None
    }
}

/// ZA worm: four helicoidal bands per turn of the single-start thread.
pub(super) struct Screw {
    worm: Worm,
    /// Rail phases in counterclockwise order, closing one pitch.
    phases: [f64; 5],
    radii: [f64; 5],
    rate: f64,
    turn: f64,
    half: f64,
}

impl Screw {
    pub(super) fn new(worm: &Worm) -> Self {
        let thread = worm.thread();
        let [_, crest_end, upper_end, lower_end, crest_start, pitch] = thread.profile_breaks();
        let phases = [
            crest_start,
            lower_end,
            upper_end,
            crest_end,
            crest_start - pitch,
        ];
        Self {
            worm: worm.clone(),
            phases,
            radii: phases.map(|phase| thread.radius_at_phase(phase)),
            rate: thread.helical_rate(),
            turn: worm.phase_rad(),
            half: worm.length_mm() / 2.,
        }
    }
}

impl Shell for Screw {
    fn half(&self) -> f64 {
        self.half
    }

    fn bore(&self) -> Option<f64> {
        self.worm.bore_radius_mm()
    }

    fn patches(&self) -> usize {
        4
    }

    fn eval(&self, patch: usize, u: f64, z: f64) -> Jet3 {
        let [s0, s1] = [self.phases[patch], self.phases[patch + 1]];
        let [r0, r1] = [self.radii[patch], self.radii[patch + 1]];
        let s = s0 + (s1 - s0) * u;
        // The profile is affine between its breaks.
        let r = r0 + (r1 - r0) * u;
        let theta = (z + self.half - s) / self.rate + self.turn;
        let (sin, cos) = theta.sin_cos();
        let spin = -(s1 - s0) / self.rate;
        let dr = r1 - r0;
        Jet3 {
            p: [r * cos, r * sin, z],
            du: [dr * cos - r * sin * spin, dr * sin + r * cos * spin, 0.],
            dz: [-r * sin / self.rate, r * cos / self.rate, 1.],
        }
    }

    fn joint(&self, _: usize) -> Joint {
        Joint::Edge
    }

    fn grid(&self, patch: usize) -> [usize; 2] {
        let turns = 2. * self.half / (self.rate.abs() * TAU);
        [
            if self.radii[patch] == self.radii[patch + 1] {
                4
            } else {
                8
            },
            // A height cell spans 1/256 of a turn along the rails.
            (256. * turns).ceil() as usize,
        ]
    }

    fn arc(&self, patch: usize, z: f64) -> Option<[f64; 3]> {
        let [s0, s1] = [self.phases[patch], self.phases[patch + 1]];
        (self.radii[patch] == self.radii[patch + 1]).then(|| {
            [
                self.radii[patch],
                (z + self.half - s0) / self.rate + self.turn,
                (s0 - s1) / self.rate,
            ]
        })
    }

    fn radius(&self, theta: f64, z: f64) -> f64 {
        self.worm.thread().radius(theta - self.turn, z + self.half)
    }

    fn annulus(&self) -> [f64; 2] {
        [self.worm.root_radius_mm(), self.worm.major_radius_mm()]
    }

    fn lipschitz(&self) -> [f64; 2] {
        let slope = self.worm.thread().maximum_profile_slope();
        [self.rate.abs() * slope, slope]
    }

    fn support(&self, direction: V3) -> f64 {
        self.worm.support(direction)
    }

    /// Along the line the profile phase is `s = z + half - rate (θ - turn)`;
    /// the boundary radius equals the core radius exactly on the root band
    /// `[upper_end, lower_end]` of each pitch and exceeds it elsewhere.
    fn core_line(&self, theta: f64, [from, to]: [f64; 2]) -> Option<Vec<[f64; 2]>> {
        let thread = self.worm.thread();
        let [_, _, upper_end, lower_end, _, pitch] = thread.profile_breaks();
        let shift = self.half - self.rate * (theta - self.turn);
        // Root bands in height: `[upper_end - shift, lower_end - shift] + k pitch`.
        let first = ((from + shift - lower_end) / pitch).floor() as i64;
        let mut output = Vec::new();
        let mut start = from;
        let mut k = first;
        loop {
            let [low, high] = [upper_end, lower_end].map(|b| b - shift + k as f64 * pitch);
            if low >= to {
                break;
            }
            if low > start {
                output.push([start, low]);
            }
            start = start.max(high);
            k += 1;
        }
        if to > start {
            output.push([start, to]);
        }
        Some(output)
    }
}

/// Exit of one flank from the blank at a height: the cutter radius, its
/// height derivative, and the section angle with its height derivative.
#[derive(Clone, Copy, Debug)]
struct Exit {
    radius: f64,
    slope: f64,
    angle: f64,
    turn: f64,
}

/// One section of sector 0 in module units.
#[derive(Clone, Copy, Debug)]
struct Station {
    v: f64,
    exits: [Exit; 2],
    /// Section angles of the six piece ends from the minus to the plus exit.
    joints: [f64; 6],
}

/// Generated worm wheel: per sector a tip land, the five generated pieces
/// of the tooth space and the second tip land.
pub(super) struct Teeth {
    wheel: Wheel,
    m: f64,
    half: f64,
    span: f64,
    blank: f64,
    tip: f64,
    root: f64,
    bounds: [f64; 2],
    joints: [Joint; 7],
    last: Mutex<Option<Station>>,
}

const CHAIN: usize = 7;

fn cross2(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}

fn dot2(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}

impl Teeth {
    pub(super) fn new(wheel: &Wheel) -> Self {
        let m = wheel.module_mm();
        let [root, tip] = wheel.flank_radii();
        let mut teeth = Self {
            wheel: wheel.clone(),
            m,
            half: wheel.half_width_mm() / m,
            span: PI / f64::from(wheel.teeth()),
            blank: wheel.blank_radius(),
            tip,
            root,
            bounds: [0.; 2],
            joints: [Joint::Edge; CHAIN],
            last: Mutex::new(None),
        };
        teeth.joints[0] = Joint::Seam;
        // Slopes of the polar section graph, sampled over every piece.
        let mut bounds = [0_f64; 2];
        for piece in 1..6 {
            for i in 0..=32 {
                for j in 0..=16 {
                    let u = f64::from(i) / 32.;
                    let v = teeth.half * (f64::from(j) / 8. - 1.);
                    let [p, du, dv] = teeth.piece(&teeth.station(v), piece, u);
                    let turn = cross2(p, du);
                    let r = p[0].hypot(p[1]);
                    bounds[0] = bounds[0].max((r * dot2(p, du) / turn).abs());
                    bounds[1] = bounds[1]
                        .max(((dot2(p, dv) - dot2(p, du) * cross2(p, dv) / turn) / r).abs());
                }
            }
        }
        teeth.bounds = [2. * bounds[0] * m, 2. * bounds[1]];
        // A junction is smooth where the unit normals of both pieces agree
        // across the whole width.
        for patch in 2..6 {
            let mut worst = 0_f64;
            for j in 0..=16 {
                let z = teeth.half * m * (f64::from(j) / 8. - 1.);
                let [a, b] = [(patch - 1, 1.), (patch, 0.)]
                    .map(|(patch, u)| super::unit(teeth.eval(patch, u, z).normal()));
                worst = worst.max(super::norm(super::cross(a, b)));
            }
            if worst < 1e-7 {
                teeth.joints[patch] = Joint::Smooth;
            }
        }
        teeth
    }

    /// The exit of a flank: Newton on the squared radius, which decreases
    /// strictly along the cutter radius, safeguarded by bisection.
    fn exit(&self, v: f64, piece: Piece) -> Exit {
        let target = self.blank * self.blank;
        let value = |a: f64| {
            let [p, da, _] = self.wheel.surface_frame(piece, a, v);
            (dot2(p, p) - target, 2. * dot2(p, da))
        };
        let (mut lo, mut hi) = (self.root, self.tip);
        let mut a = 0.5 * (lo + hi);
        for _ in 0..80 {
            let (f, slope) = value(a);
            if f > 0. {
                lo = a;
            } else {
                hi = a;
            }
            let next = a - f / slope;
            let next = if next > lo && next < hi {
                next
            } else {
                0.5 * (lo + hi)
            };
            if (next - a).abs() <= 4. * f64::EPSILON * a {
                a = next;
                break;
            }
            a = next;
        }
        let [p, da, dv] = self.wheel.surface_frame(piece, a, v);
        let slope = -dot2(p, dv) / dot2(p, da);
        let rise = [dv[0] + da[0] * slope, dv[1] + da[1] * slope];
        Exit {
            radius: a,
            slope,
            angle: p[1].atan2(p[0]),
            turn: cross2(p, rise) / dot2(p, p),
        }
    }

    fn station(&self, v: f64) -> Station {
        let mut last = self.last.lock().unwrap_or_else(|e| e.into_inner());
        if let Some(station) = *last
            && station.v == v
        {
            return station;
        }
        let exits = [Piece::FlankMinus, Piece::FlankPlus].map(|piece| self.exit(v, piece));
        let mut station = Station {
            v,
            exits,
            joints: [0.; 6],
        };
        station.joints[0] = exits[0].angle;
        station.joints[5] = exits[1].angle;
        for joint in 1..5 {
            let [p, ..] = self.piece(&station, joint + 1, 0.);
            station.joints[joint] = p[1].atan2(p[0]);
        }
        *last = Some(station);
        station
    }

    /// Chain piece `1..=5` of sector 0 over `(u, v)` in module units.
    fn piece(&self, station: &Station, patch: usize, u: f64) -> [[f64; 2]; 3] {
        let piece = PIECES[patch - 1];
        let v = station.v;
        let side = match piece {
            Piece::FlankMinus => 0,
            Piece::FlankPlus => 1,
            _ => return self.wheel.surface_frame(piece, u, v),
        };
        let exit = station.exits[side];
        // The retained flank runs linearly in cutter radius between its
        // exit and the cutter tip.
        let weight = if side == 0 { 1. - u } else { u };
        let a = self.tip - (self.tip - exit.radius) * weight;
        let [p, da, dv] = self.wheel.surface_frame(piece, a, v);
        let along = (self.tip - exit.radius) * if side == 0 { 1. } else { -1. };
        let lift = exit.slope * weight;
        [
            p,
            [da[0] * along, da[1] * along],
            [dv[0] + da[0] * lift, dv[1] + da[1] * lift],
        ]
    }

    /// Chain patch of sector 0 in module units.
    fn chain(&self, station: &Station, patch: usize, u: f64) -> [[f64; 2]; 3] {
        let [low, high] = station.exits;
        let (from, to, from_turn, to_turn) = match patch {
            0 => (-self.span, low.angle, 0., low.turn),
            6 => (high.angle, self.span, high.turn, 0.),
            _ => return self.piece(station, patch, u),
        };
        let angle = from + (to - from) * u;
        let (sin, cos) = angle.sin_cos();
        let tangent = [-self.blank * sin, self.blank * cos];
        let turn = from_turn + (to_turn - from_turn) * u;
        [
            [self.blank * cos, self.blank * sin],
            [tangent[0] * (to - from), tangent[1] * (to - from)],
            [tangent[0] * turn, tangent[1] * turn],
        ]
    }

    fn sector(&self, patch: usize) -> f64 {
        2. * self.span * (patch / CHAIN) as f64
    }

    /// Radius of the tooth space boundary of sector 0 at a section angle
    /// strictly between the exits, in module units.
    fn space(&self, station: &Station, angle: f64) -> f64 {
        let patch = 1
            + (1..5)
                .rev()
                .find(|&joint| station.joints[joint] <= angle)
                .unwrap_or(0);
        // The section angle increases strictly along every piece.
        let (mut lo, mut hi) = (0_f64, 1_f64);
        let mut u = (angle - station.joints[patch - 1])
            / (station.joints[patch] - station.joints[patch - 1]);
        u = u.clamp(0., 1.);
        let mut radius = self.blank;
        for _ in 0..80 {
            let [p, du, _] = self.piece(station, patch, u);
            radius = p[0].hypot(p[1]);
            let f = p[1].atan2(p[0]) - angle;
            if f < 0. {
                lo = u;
            } else {
                hi = u;
            }
            let next = u - f * radius * radius / cross2(p, du);
            let next = if next > lo && next < hi {
                next
            } else {
                0.5 * (lo + hi)
            };
            if (next - u).abs() <= 4. * f64::EPSILON || hi - lo <= 4. * f64::EPSILON {
                let [p, ..] = self.piece(station, patch, next);
                return p[0].hypot(p[1]);
            }
            u = next;
        }
        radius
    }
}

impl Shell for Teeth {
    fn half(&self) -> f64 {
        self.half * self.m
    }

    fn bore(&self) -> Option<f64> {
        self.wheel.bore_radius_mm()
    }

    fn patches(&self) -> usize {
        CHAIN * self.wheel.teeth() as usize
    }

    fn eval(&self, patch: usize, u: f64, z: f64) -> Jet3 {
        let station = self.station(z / self.m);
        let [p, du, dv] = self.chain(&station, patch % CHAIN, u);
        let (sin, cos) = self.sector(patch).sin_cos();
        let turn = |[x, y]: [f64; 2], scale: f64, z: f64| {
            [scale * (cos * x - sin * y), scale * (sin * x + cos * y), z]
        };
        Jet3 {
            p: turn(p, self.m, z),
            du: turn(du, self.m, 0.),
            dz: turn(dv, 1., 1.),
        }
    }

    fn joint(&self, patch: usize) -> Joint {
        self.joints[patch % CHAIN]
    }

    fn grid(&self, patch: usize) -> [usize; 2] {
        match patch % CHAIN {
            0 | 6 => [3, 6],
            _ => [8, 8],
        }
    }

    fn arc(&self, patch: usize, z: f64) -> Option<[f64; 3]> {
        let chain = patch % CHAIN;
        if ![0, 3, 6].contains(&chain) {
            return None;
        }
        let station = self.station(z / self.m);
        let [a, b] = [0., 1.].map(|u| self.chain(&station, chain, u)[0]);
        let from = a[1].atan2(a[0]);
        Some([
            self.m * a[0].hypot(a[1]),
            from + self.sector(patch),
            b[1].atan2(b[0]) - from,
        ])
    }

    fn radius(&self, theta: f64, z: f64) -> f64 {
        let pitch = 2. * self.span;
        let local = theta - (theta / pitch).round() * pitch;
        let station = self.station(z / self.m);
        let [low, high] = station.exits;
        self.m
            * if local <= low.angle || local >= high.angle {
                self.blank
            } else {
                self.space(&station, local)
            }
    }

    fn annulus(&self) -> [f64; 2] {
        // The throat of the root lies at the mid plane.
        let [p, ..] = self.wheel.surface_frame(Piece::Torus, 0.5, 0.);
        [self.m * p[0].hypot(p[1]), self.m * self.blank]
    }

    fn lipschitz(&self) -> [f64; 2] {
        self.bounds
    }

    fn support(&self, direction: V3) -> f64 {
        self.wheel.support(direction)
    }
}
