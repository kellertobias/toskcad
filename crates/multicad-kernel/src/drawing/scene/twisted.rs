//! Helical and continuous herringbone gears as exact scene members.
//!
//! The solid is the transverse profile rotated by a phase that is linear in
//! the axial coordinate within each band. Every answer is evaluated on that
//! law: a flank is the surface `S(s, z) = (Rot(phase(z)) c(s), z)` of a
//! profile edge `c`, with the outward normal
//! `(Rot(phase) nu(s), twist * c(s) . tau(s))` for the edge's unit tangent
//! `tau` and outward normal `nu`.
use super::paths::{dot2, dot3};
use super::*;
use multicad_rust_kernel::spur_gear::Edge;
use std::f64::consts::TAU;

mod cut;
mod flank;
mod outline;
mod sight;
mod solve;
mod trace;

#[cfg(test)]
mod tests;

use flank::{Flank, rotate, turn};
use trace::Step;

/// Below this transverse or axial part of the sight direction the view is
/// treated as axial or edge-on, as for the prism members.
const AXIAL: f64 = 1e-9;

/// An axial range with a linear phase law.
#[derive(Clone, Copy, Debug)]
pub(super) struct Band {
    pub z: [f64; 2],
    /// Phase at `z[0]`.
    pub phase: f64,
    pub twist: f64,
}

impl Band {
    pub(super) fn phase_at(&self, z: f64) -> f64 {
        self.phase + self.twist * (z - self.z[0])
    }

    fn height(&self) -> f64 {
        self.z[1] - self.z[0]
    }
}

/// Rigid placement of the gear frame in view space.
#[derive(Clone, Copy, Debug)]
pub(super) struct Frame {
    pub origin: P3,
    pub axes: [P3; 3],
}

impl Frame {
    pub(super) fn place(&self, local: P3) -> P3 {
        std::array::from_fn(|i| {
            self.origin[i]
                + local[0] * self.axes[0][i]
                + local[1] * self.axes[1][i]
                + local[2] * self.axes[2][i]
        })
    }

    pub(super) fn vector(&self, local: P3) -> P3 {
        std::array::from_fn(|i| {
            local[0] * self.axes[0][i] + local[1] * self.axes[1][i] + local[2] * self.axes[2][i]
        })
    }

    pub(super) fn local(&self, view: P3) -> P3 {
        let offset: P3 = std::array::from_fn(|i| view[i] - self.origin[i]);
        self.axes.map(|axis| dot3(axis, offset))
    }

    /// The transverse plane at `z`, turned by `phase`.
    pub(super) fn level(&self, z: f64, phase: f64) -> paths::Lift {
        let (sin, cos) = phase.sin_cos();
        paths::Lift {
            origin: self.place([0., 0., z]),
            x: std::array::from_fn(|i| cos * self.axes[0][i] + sin * self.axes[1][i]),
            y: std::array::from_fn(|i| cos * self.axes[1][i] - sin * self.axes[0][i]),
        }
    }
}

/// Which zero set of a flank is traced.
#[derive(Clone, Copy, Debug)]
pub(super) enum Level {
    /// The plane `toward . p = offset` in gear coordinates.
    Cut(f64),
    /// The outline: the normal is perpendicular to the sight direction.
    Fold,
}

/// One flank surface: a profile edge swept through one band.
#[derive(Clone, Copy, Debug)]
pub(super) struct Patch {
    pub flank: Flank,
    pub band: Band,
    pub frame: Frame,
    /// Direction toward the viewer in gear coordinates.
    pub toward: P3,
}

impl Patch {
    fn z(&self, zeta: f64) -> f64 {
        self.band.z[0] + zeta * self.band.height()
    }

    pub(super) fn local(&self, s: f64, zeta: f64) -> P3 {
        let z = self.z(zeta);
        let [x, y] = rotate(self.flank.point(s), self.band.phase_at(z));
        [x, y, z]
    }

    pub(super) fn surface(&self, s: f64, zeta: f64) -> P3 {
        self.frame.place(self.local(s, zeta))
    }

    /// Outward normal in view space, not normalized.
    pub(super) fn normal(&self, s: f64, zeta: f64) -> P3 {
        let tangent = self.flank.tangent(s);
        let [x, y] = rotate([tangent[1], -tangent[0]], self.band.phase_at(self.z(zeta)));
        self.frame
            .vector([x, y, self.band.twist * dot2(self.flank.point(s), tangent)])
    }

    /// Value and gradient of a level function over the unit square of the
    /// edge fraction and the band fraction.
    pub(super) fn field(&self, level: Level, s: f64, zeta: f64) -> (f64, [f64; 2]) {
        let height = self.band.height();
        let z = self.z(zeta);
        let twist = self.band.twist;
        // The transverse sight direction in profile coordinates.
        let e = rotate([self.toward[0], self.toward[1]], -self.band.phase_at(z));
        let [c, c1, c2, _] = self.flank.at(s);
        match level {
            Level::Cut(offset) => (
                dot2(e, c) + self.toward[2] * z - offset,
                [
                    dot2(e, c1),
                    height * (twist * dot2(e, turn(c)) + self.toward[2]),
                ],
            ),
            Level::Fold => {
                let lead = twist * self.toward[2];
                (
                    -dot2(e, turn(c1)) + lead * dot2(c, c1),
                    [
                        -dot2(e, turn(c2)) + lead * (dot2(c1, c1) + dot2(c, c2)),
                        height * twist * dot2(e, c1),
                    ],
                )
            }
        }
    }

    pub(super) fn bounds(&self, level: Level) -> trace::Bounds {
        let height = self.band.height();
        let twist = self.band.twist.abs();
        let slope = self.toward[0].hypot(self.toward[1]);
        let [c0, c1, c2, c3] = self.flank.size;
        match level {
            Level::Cut(_) => trace::Bounds {
                first: [
                    slope * c1,
                    height * (twist * slope * c0 + self.toward[2].abs()),
                ],
                second: [
                    slope * c2,
                    height * twist * slope * c1,
                    height * height * twist * twist * slope * c0,
                ],
            },
            Level::Fold => {
                let lead = twist * self.toward[2].abs();
                trace::Bounds {
                    first: [
                        slope * c2 + lead * (c1 * c1 + c0 * c2),
                        height * twist * slope * c1,
                    ],
                    second: [
                        slope * c3 + lead * (3. * c1 * c2 + c0 * c3),
                        height * twist * slope * c2,
                        height * height * twist * twist * slope * c1,
                    ],
                }
            }
        }
    }
}

/// A traced curve of a flank, evaluated on the flank itself.
pub(super) struct Traced {
    pub patch: Patch,
    pub level: Level,
    pub steps: Arc<Vec<Step>>,
    /// Range of the chain parameter, possibly descending.
    pub range: [f64; 2],
}

impl Traced {
    /// Edge fraction and band fraction at a path parameter.
    pub(super) fn parameters(&self, t: f64) -> [f64; 2] {
        let count = self.steps.len() as f64;
        let at = (self.range[0] + (self.range[1] - self.range[0]) * t).clamp(0., 1.) * count;
        let index = (at.floor() as usize).min(self.steps.len() - 1);
        let patch = self.patch;
        let level = self.level;
        self.steps[index].point(&|s, zeta| patch.field(level, s, zeta), at - index as f64)
    }
}

fn distance(a: P3, b: P3) -> f64 {
    let d: P3 = std::array::from_fn(|i| a[i] - b[i]);
    dot3(d, d).sqrt()
}

/// Ascending parameters whose chords stay within the tolerance at the
/// midpoint and both quarter points of every chord.
pub(super) fn chords(point: &dyn Fn(f64) -> P3, seeds: Vec<f64>, tolerance: f64) -> Vec<f64> {
    fn refine(
        point: &dyn Fn(f64) -> P3,
        [a, b]: [f64; 2],
        [pa, pb]: [P3; 2],
        tolerance: f64,
        depth: u32,
        output: &mut Vec<f64>,
    ) {
        let mid = 0.5 * (a + b);
        let pm = point(mid);
        let deviation = [0.25, 0.5, 0.75]
            .into_iter()
            .map(|f| {
                let on_chord: P3 = std::array::from_fn(|i| pa[i] + (pb[i] - pa[i]) * f);
                let on_curve = if f == 0.5 { pm } else { point(a + (b - a) * f) };
                distance(on_chord, on_curve)
            })
            .fold(0., f64::max);
        if deviation <= tolerance || depth >= 20 {
            output.push(b);
        } else {
            refine(point, [a, mid], [pa, pm], tolerance, depth + 1, output);
            refine(point, [mid, b], [pm, pb], tolerance, depth + 1, output);
        }
    }
    let mut output = vec![seeds[0]];
    for pair in seeds.windows(2) {
        refine(
            point,
            [pair[0], pair[1]],
            [point(pair[0]), point(pair[1])],
            tolerance,
            0,
            &mut output,
        );
    }
    output
}

impl Path for Traced {
    fn point(&self, t: f64) -> P3 {
        let [s, zeta] = self.parameters(t);
        self.patch.surface(s, zeta)
    }

    fn flatten(&self, tolerance_mm: f64) -> Vec<f64> {
        // Seed with the chart ends inside the range.
        let count = self.steps.len() as f64;
        let [a, b] = self.range;
        let mut seeds = vec![0.];
        let mut inner: Vec<f64> = (1..self.steps.len())
            .map(|i| (i as f64 / count - a) / (b - a))
            .filter(|t| *t > 1e-12 && *t < 1. - 1e-12)
            .collect();
        inner.sort_by(f64::total_cmp);
        seeds.extend(inner);
        seeds.push(1.);
        chords(&|t| self.point(t), seeds, tolerance_mm)
    }

    fn paper(&self, _: [f64; 2]) -> Option<Curve> {
        None
    }
}

/// The helix traced by one profile point through a range of one band.
pub(super) struct Helix {
    pub frame: Frame,
    pub band: Band,
    pub point: P2,
    pub z: [f64; 2],
    /// The helix is seen along its axis, as a circular arc.
    pub axial: bool,
}

impl Path for Helix {
    fn point(&self, t: f64) -> P3 {
        let z = self.z[0] + (self.z[1] - self.z[0]) * t;
        let [x, y] = rotate(self.point, self.band.phase_at(z));
        self.frame.place([x, y, z])
    }

    fn flatten(&self, tolerance_mm: f64) -> Vec<f64> {
        let radius = self.point[0].hypot(self.point[1]);
        let step = if tolerance_mm < radius {
            2. * (1. - tolerance_mm / radius).acos()
        } else {
            TAU / 6.
        };
        let turn = (self.band.twist * (self.z[1] - self.z[0])).abs();
        let count = ((turn / step.min(TAU / 6.)).ceil() as usize).clamp(1, 1 << 20);
        (0..=count).map(|i| i as f64 / count as f64).collect()
    }

    fn paper(&self, [from, to]: [f64; 2]) -> Option<Curve> {
        if !self.axial {
            return None;
        }
        let z = [from, to].map(|t| self.z[0] + (self.z[1] - self.z[0]) * t);
        let radius = self.point[0].hypot(self.point[1]);
        let start = self.point[1].atan2(self.point[0]) + self.band.phase_at(z[0]);
        let [x, y] = [self.frame.axes[0], self.frame.axes[1]];
        paths::conic(
            [self.frame.origin[0], self.frame.origin[1]],
            [radius * x[0], radius * x[1]],
            [radius * y[0], radius * y[1]],
            start,
            self.band.twist * (z[1] - z[0]),
        )
    }
}

/// A coaxial cylinder of the hub: an extension beyond the tooth band.
#[derive(Clone, Copy, Debug)]
pub(super) struct Extension {
    pub radius: f64,
    pub z: [f64; 2],
    /// Which end of `z` carries the free end face.
    pub outer: usize,
}

pub(super) struct Gear {
    pub frame: Frame,
    pub toward: P3,
    /// Length and azimuth of the transverse part of `toward`.
    pub slope: f64,
    pub azimuth: f64,
    pub flanks: Vec<Flank>,
    pub bands: Vec<Band>,
    pub height: f64,
    /// Zero without a bore.
    pub bore: f64,
    pub extensions: Vec<Extension>,
    /// Full axial range of the material.
    pub axis: [f64; 2],
    /// Least and greatest radius of the tooth outline.
    pub radius: [f64; 2],
    /// Angular tooth pitch.
    pub pitch: f64,
    pub tolerance: f64,
    bounds: [P3; 2],
}

fn invalid(message: &'static str) -> GeometryError {
    err(GeometryErrorKind::InvalidResult, message)
}

impl Gear {
    fn new(solid: &Solid, view: &View) -> Option<Result<Self, GeometryError>> {
        let (body, placement) = solid.herringbone_body()?;
        let profile = body.profile();
        let frame = Frame {
            origin: view.point(placement.origin),
            axes: [placement.u, placement.v, placement.normal()].map(|axis| view.vector(axis)),
        };
        let toward = frame.axes.map(|axis| axis[2]);
        let height = profile.height_mm();
        let bands = if body.is_helical() {
            vec![Band {
                z: [0., height],
                phase: 0.,
                twist: body.phase_at(1.),
            }]
        } else {
            let middle = 0.5 * height;
            let twist = body.phase_at(middle) / middle;
            vec![
                Band {
                    z: [0., middle],
                    phase: 0.,
                    twist,
                },
                Band {
                    z: [middle, height],
                    phase: twist * middle,
                    twist: -twist,
                },
            ]
        };
        let flanks: Vec<Flank> = profile
            .outline()
            .iter()
            .map(|edge| Flank::new(edge.edge))
            .collect();
        let radius = flanks.iter().fold([f64::INFINITY, 0_f64], |[lo, hi], f| {
            [lo.min(f.radius[0]), hi.max(f.radius[1])]
        });
        let extensions = profile
            .hub_extensions()
            .into_iter()
            .flatten()
            .map(|extension| {
                let [a, b] = [extension.gear_z_mm, extension.end_z_mm];
                Extension {
                    radius: extension.radius_mm,
                    z: [a.min(b), a.max(b)],
                    outer: usize::from(b > a),
                }
            })
            .collect();
        let bounds = [-1., 1.].map(|sign: f64| {
            std::array::from_fn(|k| {
                let direction = frame.axes.map(|axis| sign * axis[k]);
                frame.origin[k] + sign * body.support(direction)
            })
        });
        let gear = Self {
            frame,
            toward,
            slope: toward[0].hypot(toward[1]),
            azimuth: toward[1].atan2(toward[0]),
            flanks,
            bands,
            height,
            bore: profile
                .bore()
                .and_then(multicad_rust_kernel::normal_sweep::Edge::radius_mm)
                .unwrap_or(0.),
            extensions,
            axis: profile.axial_bounds(),
            radius,
            pitch: TAU / f64::from(profile.outline().last()?.tooth + 1),
            tolerance: solid.geometry_policy().linear_tolerance_mm,
            bounds,
        };
        Some(
            if gear.flanks.is_empty()
                || !bounds.iter().flatten().all(|v| v.is_finite())
                || gear
                    .bands
                    .iter()
                    .any(|band| !(band.twist.is_finite() && band.twist != 0.))
            {
                Err(invalid("Invalid twisted gear for a drawing scene"))
            } else {
                Ok(gear)
            },
        )
    }

    pub(super) fn axial(&self) -> bool {
        self.slope <= AXIAL
    }

    pub(super) fn edge_on(&self) -> bool {
        self.toward[2].abs() <= AXIAL
    }

    pub(super) fn patch(&self, band: usize, flank: usize) -> Patch {
        Patch {
            flank: self.flanks[flank],
            band: self.bands[band],
            frame: self.frame,
            toward: self.toward,
        }
    }

    /// Phase of the tooth profile at an axial coordinate of the tooth band.
    pub(super) fn phase_at(&self, z: f64) -> f64 {
        let band = self
            .bands
            .iter()
            .find(|band| z <= band.z[1])
            .unwrap_or(&self.bands[self.bands.len() - 1]);
        band.phase_at(z)
    }

    pub(super) fn circle(radius: f64, clockwise: bool) -> Edge {
        Edge::Planar(multicad_rust_kernel::normal_sweep::Edge::Circle {
            center: [0.; 2],
            radius_mm: radius,
            clockwise,
        })
    }
}

impl Member for Gear {
    fn bounds(&self) -> [P3; 2] {
        self.bounds
    }

    fn candidates(&self, _: f64) -> Result<Vec<Candidate>, GeometryError> {
        outline::candidates(self)
    }

    fn sight(&self, paper: P2) -> Result<Vec<[f64; 2]>, GeometryError> {
        let base = self.frame.local([paper[0], paper[1], 0.]);
        self.spans(base, self.toward, &|p| self.contains(p))
    }

    fn section(&self, depth: f64, _: f64) -> Result<Vec<Loop>, GeometryError> {
        cut::section(self, depth)
    }
}

/// The scene member of a helical or continuous herringbone gear.
pub(in crate::drawing) fn member(
    solid: &Solid,
    view: &View,
) -> Option<Result<Box<dyn Member>, GeometryError>> {
    Some(Gear::new(solid, view)?.map(|gear| Box::new(gear) as Box<dyn Member>))
}
