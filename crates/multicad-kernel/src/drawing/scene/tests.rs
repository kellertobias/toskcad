//! Engine regressions that need members built directly in view space.
use super::paths::{Lift, ProfilePath};
use super::prism::{Band, Prism};
use super::*;
use multicad_rust_kernel::{normal_sweep, spur_gear::Edge};
use std::f64::consts::{PI, TAU};

/// A member that publishes extra candidates on top of another member's.
struct Extra {
    inner: Prism,
    extra: Vec<Candidate>,
}

impl Member for Extra {
    fn bounds(&self) -> [P3; 2] {
        self.inner.bounds()
    }
    fn candidates(&self, tolerance_mm: f64) -> Result<Vec<Candidate>, GeometryError> {
        let mut candidates = self.inner.candidates(tolerance_mm)?;
        candidates.extend(self.extra.iter().cloned());
        Ok(candidates)
    }
    fn sight(&self, paper: P2) -> Result<Vec<[f64; 2]>, GeometryError> {
        self.inner.sight(paper)
    }
    fn section(&self, depth: f64, tolerance_mm: f64) -> Result<Vec<Loop>, GeometryError> {
        self.inner.section(depth, tolerance_mm)
    }
}

const RADIUS: f64 = 5.;

fn circle() -> Vec<Edge> {
    vec![Edge::Planar(normal_sweep::Edge::Circle {
        center: [0.; 2],
        radius_mm: RADIUS,
        clockwise: false,
    })]
}

fn cylinder(origin: P3, axes: [P3; 3], axis: [f64; 2]) -> Prism {
    Prism::new(
        origin,
        axes,
        vec![Band {
            rings: vec![circle()],
            axis,
        }],
        1e-7,
    )
    .unwrap()
}

fn polar(angle: f64) -> P2 {
    [RADIUS * angle.cos(), RADIUS * angle.sin()]
}

/// An arc of the rim from `from` sweeping by `sweep`, as a planar arc edge.
fn arc(lift: Lift, from: f64, sweep: f64) -> Candidate {
    Candidate {
        kind: Kind::Edge,
        path: Arc::new(ProfilePath {
            edge: Edge::Planar(normal_sweep::Edge::Arc {
                center: [0.; 2],
                a: polar(from),
                b: polar(from + sweep),
                clockwise: sweep < 0.,
            }),
            range: [0., 1.],
            lift,
            straight: false,
        }),
    }
}

fn draw(members: Vec<Box<dyn Member>>, hidden_lines: bool) -> Vec<ProjectedVector> {
    let members: Vec<SceneMember> = members
        .into_iter()
        .map(|member| SceneMember {
            source: 0,
            member,
            hatch: None,
        })
        .collect();
    project(
        &members,
        DrawingProjection {
            hidden_lines,
            ..Default::default()
        },
        GeometryPolicy::default(),
    )
    .unwrap()
}

type Conic = ([f64; 2], [f64; 2], [f64; 2], f64, f64, [f64; 2]);

fn conics(vectors: &[ProjectedVector], role: VectorRole) -> Vec<Conic> {
    vectors
        .iter()
        .filter(|vector| vector.role == role)
        .filter_map(|vector| match vector.curve {
            VectorCurve::Ellipse {
                center,
                u,
                v,
                major_mm,
                minor_mm,
                start_rad,
                end_rad,
            } => Some((center, u, v, major_mm, minor_mm, [start_rad, end_rad])),
            _ => None,
        })
        .collect()
}

/// Whether a paper point lies strictly inside the trim of a published arc.
fn within((center, u, v, major, minor, [start, end]): Conic, p: P2) -> bool {
    let d = [p[0] - center[0], p[1] - center[1]];
    let angle = ((d[0] * v[0] + d[1] * v[1]) / minor).atan2((d[0] * u[0] + d[1] * u[1]) / major);
    let travel = (angle - start).rem_euclid(TAU);
    travel > 1e-9 && travel < end - start - 1e-9
}

fn at((center, u, v, major, minor, _): Conic, angle: f64) -> P2 {
    let (sin, cos) = angle.sin_cos();
    [
        center[0] + major * cos * u[0] + minor * sin * v[0],
        center[1] + major * cos * u[1] + minor * sin * v[1],
    ]
}

/// No two arcs of one role overlap, and the arcs of each conic of this
/// center and size cover the given total sweep.
fn assert_once(arcs: &[Conic], center: P2, major: f64, minor: f64, sweep: f64, label: &str) {
    let same: Vec<Conic> = arcs
        .iter()
        .copied()
        .filter(|c| {
            (c.0[0] - center[0]).hypot(c.0[1] - center[1]) < 1e-9
                && (c.3 - major).abs() < 1e-9
                && (c.4 - minor).abs() < 1e-9
        })
        .collect();
    for (i, a) in same.iter().enumerate() {
        for (j, b) in same.iter().enumerate() {
            if i != j {
                for q in [0.25, 0.5, 0.75] {
                    let p = at(*a, a.5[0] + q * (a.5[1] - a.5[0]));
                    assert!(!within(*b, p), "{label}: arcs {i} and {j} overlap");
                }
            }
        }
    }
    let total: f64 = same.iter().map(|c| c.5[1] - c.5[0]).sum();
    assert!((total - sweep).abs() < 1e-9, "{label}: {total} of {sweep}");
}

/// Seen along the axis, both rims and two further copies of the top rim,
/// split at other angles and one of them traversed clockwise, coincide: the
/// rim publishes once, as a full visible circle, and nothing hidden.
#[test]
fn one_circle_in_several_parametrisations_publishes_once() {
    let axes = [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
    let inner = cylinder([0.; 3], axes, [0., 1.]);
    let top = Lift {
        origin: [0., 0., 1.],
        x: axes[0],
        y: axes[1],
    };
    let extra = vec![
        arc(top, 1., PI),
        arc(top, 1. + PI, PI),
        arc(top, 2.5, -2.),
        arc(top, 0.5, -(TAU - 2.)),
    ];
    let vectors = draw(vec![Box::new(Extra { inner, extra })], true);
    assert!(conics(&vectors, VectorRole::Hidden).is_empty());
    let visible = conics(&vectors, VectorRole::Visible);
    assert_once(&visible, [0.; 2], RADIUS, RADIUS, TAU, "axial");
}

/// Two coaxial cylinders of one source stacked with their profile charts
/// turned against each other, seen obliquely: their rims at the shared level
/// project onto one ellipse with different parametrisations. Each arc of it
/// publishes once, the front visible and the back hidden.
#[test]
fn one_ellipse_from_two_members_publishes_once() {
    let tilt = 0.6_f64;
    let (s, c) = tilt.sin_cos();
    let axes = [[1., 0., 0.], [0., c, s], [0., -s, c]];
    let (st, ct) = 0.7_f64.sin_cos();
    let turned = [
        std::array::from_fn(|i| ct * axes[0][i] + st * axes[1][i]),
        std::array::from_fn(|i| -st * axes[0][i] + ct * axes[1][i]),
        axes[2],
    ];
    let lower = cylinder([0.; 3], axes, [0., 1.]);
    let upper = cylinder([0.; 3], turned, [1., 2.]);
    let vectors = draw(vec![Box::new(lower), Box::new(upper)], true);
    // The rim at w = 1 is centered at the paper image of the axis point.
    let center = [0., -s];
    let [major, minor] = [RADIUS, RADIUS * c];
    let visible = conics(&vectors, VectorRole::Visible);
    let hidden = conics(&vectors, VectorRole::Hidden);
    assert_once(&visible, center, major, minor, PI, "visible");
    assert_once(&hidden, center, major, minor, PI, "hidden");
}
