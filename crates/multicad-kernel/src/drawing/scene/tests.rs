//! Engine regressions that need members built directly in view space.
use super::paths::{Lift, ProfilePath, dot3};
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

/// Material along lines through a thread member agrees with dense marching of
/// the source solid's own containment test, for both thread families, both
/// hands and lines in every direction, including lines along the axis and in
/// a cap plane.
#[test]
fn thread_material_along_lines_agrees_with_source_containment() {
    use crate::{HelicalThread, Solid};
    use multicad_core::document::EntityId;
    use multicad_rust_kernel::rounded_thread::Spec;
    let policy = GeometryPolicy::default();
    let frame = PlaneFrame {
        origin: [1.5, -2., 0.5],
        u: [0.8, 0.6, 0.],
        v: [-0.36, 0.48, 0.8],
    };
    let mut solids = Vec::new();
    for left_hand in [false, true] {
        solids.push(two_start(left_hand, frame));
        solids.push(rounded_bore(left_hand, frame));
        solids.push(chamfered(left_hand, left_hand, frame));
        solids.push(screw(left_hand, left_hand, frame));
        solids.push(treated(left_hand, 1 + u32::from(left_hand), frame));
        solids.push(grooved(!left_hand, 1 + u32::from(left_hand), frame));
    }
    for left_hand in [false, true] {
        let root = 1. / (4. * 3_f64.sqrt());
        solids.push(
            Solid::rust_helical_thread_end_relief(
                HelicalThread {
                    major_radius_mm: 3.,
                    pitch_mm: 1.,
                    length_mm: 6.6,
                    depth_mm: (1. - 0.125) / (2. * 30_f64.to_radians().tan()) - root,
                    crest_width_mm: 0.125,
                    lower_flank_deg: 30.,
                    upper_flank_deg: 30.,
                    radial_taper: 0.,
                    starts: 1,
                    left_hand,
                },
                Some(root),
                0.6,
                0.1,
                EntityId::new(),
                frame,
                policy,
            )
            .unwrap(),
        );
        let p = 1.27;
        let a = 27.5_f64.to_radians();
        let r = p * a.cos() / (12. * (1. - a.sin()));
        solids.push(
            Solid::rust_rounded_thread(
                Spec {
                    major_radius_mm: 3.175,
                    pitch_mm: p,
                    length_mm: 3.7 * p,
                    depth_mm: p / (3. * a.tan()),
                    crest_radius_mm: r,
                    root_radius_mm: r,
                    flank_angle_deg: 27.5,
                    starts: 1,
                    left_hand,
                },
                EntityId::new(),
                frame,
                policy,
            )
            .unwrap(),
        );
    }
    let view = View::new(PlaneFrame::default());
    let axis = frame.normal();
    let world = |l: [f64; 3]| -> P3 {
        std::array::from_fn(|i| {
            frame.origin[i] + l[0] * frame.u[i] + l[1] * frame.v[i] + l[2] * axis[i]
        })
    };
    let lines: Vec<(P3, P3)> = {
        let mut lines = Vec::new();
        // Generic directions through and beside the thread.
        for (k, d) in [
            [0.3, -0.5, 0.81],
            [1., 0.2, -0.1],
            [0.05, 0.99, 0.1],
            [-0.7, 0.7, 0.14],
        ]
        .into_iter()
        .enumerate()
        {
            let length: f64 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
            let d = d.map(|v: f64| v / length.sqrt());
            for offset in [0., 1.3, 2.4, 2.9] {
                let start = world([
                    offset * (k as f64 * 0.7).cos(),
                    offset * (k as f64 * 0.7).sin(),
                    2.,
                ]);
                lines.push((
                    start,
                    std::array::from_fn(|i| d[0] * frame.u[i] + d[1] * frame.v[i] + d[2] * axis[i]),
                ));
            }
        }
        // Along the axis at several radii, and in the entry cap plane.
        for r in [0.5, 2.5, 2.8, 3.05, 3.3, 3.6] {
            lines.push((world([r, 0.3, -1.]), axis));
        }
        lines.push((world([-4., 0.37, 1e-3]), frame.u));
        // In a treated bore's groove shoulder plane, and just beside it.
        for z in [5.3, 5.3 - 2e-3, 5.3 + 2e-3] {
            lines.push((world([-4.6, 0.41, z]), frame.u));
        }
        lines
    };
    let mut covered = 0;
    for solid in &solids {
        let member = families::member(solid, &view).unwrap().unwrap();
        for &(start, direction) in &lines {
            let spans = member.line(start, direction).unwrap().unwrap();
            // March from t = -12 to 12 in steps of 2 um.
            let steps = 12_000;
            let inside =
                |t: f64| contains(solid, std::array::from_fn(|i| start[i] + t * direction[i]));
            let mut mismatches = 0;
            for i in 0..=steps {
                let t = -12. + 24. * i as f64 / steps as f64;
                let near = spans.iter().flatten().any(|end| (end - t).abs() < 1e-4);
                if near {
                    continue;
                }
                let claimed = spans.iter().any(|[a, b]| *a < t && t < *b);
                if claimed != inside(t) {
                    mismatches += 1;
                }
            }
            assert_eq!(mismatches, 0, "{start:?} {direction:?}: {spans:?}");
            covered += usize::from(!spans.is_empty());
        }
    }
    // Most lines meet material; the march is not vacuous.
    assert!(covered >= 3 * solids.len() * lines.len() / 4, "{covered}");
}

/// Every outline candidate of a rounded thread lies on the source surface
/// with its normal perpendicular to the sight, in side, oblique and
/// near-axial views, for both hands.
#[test]
fn rounded_thread_outlines_satisfy_the_contour_condition() {
    use crate::Solid;
    use multicad_core::document::EntityId;
    use multicad_rust_kernel::rounded_thread::Spec;
    let p = 1.27;
    let a = 27.5_f64.to_radians();
    let r = p * a.cos() / (12. * (1. - a.sin()));
    for left_hand in [false, true] {
        let solid = Solid::rust_rounded_thread(
            Spec {
                major_radius_mm: 3.175,
                pitch_mm: p,
                length_mm: 3.7 * p,
                depth_mm: p / (3. * a.tan()),
                crest_radius_mm: r,
                root_radius_mm: r,
                flank_angle_deg: 27.5,
                starts: 1,
                left_hand,
            },
            EntityId::new(),
            PlaneFrame::default(),
            GeometryPolicy::default(),
        )
        .unwrap();
        let thread = solid.rust_rounded_thread_data().unwrap().thread;
        let s = 0.5_f64.sqrt();
        let t = (1_f64 / 6.).sqrt();
        for [u, v] in [
            [[1., 0., 0.], [0., 0., 1.]],
            [[s, -s, 0.], [t, t, -2. * t]],
            [[1., 0., 0.], [0., 0.9995_f64.sqrt(), (0.0005_f64).sqrt()]],
        ] {
            let frame = PlaneFrame {
                origin: [0.2, 0.1, 2.35],
                u,
                v,
            };
            let view = View::new(frame);
            let toward = view.axes[2];
            let member = families::member(&solid, &view).unwrap().unwrap();
            let mut checked = 0;
            for candidate in member.candidates(0.0025).unwrap() {
                if candidate.kind != Kind::Outline {
                    continue;
                }
                for k in 0..=8 {
                    let q = candidate.path.point(k as f64 / 8.);
                    // Back to world (the default thread frame is the world).
                    let w: [f64; 3] = std::array::from_fn(|i| {
                        frame.origin[i] + (0..3).map(|j| view.axes[j][i] * q[j]).sum::<f64>()
                    });
                    let theta = w[1].atan2(w[0]);
                    let radius = thread.radius(theta, w[2]);
                    assert!((w[0].hypot(w[1]) - radius).abs() < 1e-9, "{w:?}");
                    let normal = thread.normal(theta, w[2]);
                    let dot: f64 = (0..3).map(|i| normal[i] * toward[i]).sum();
                    assert!(dot.abs() < 1e-7, "{w:?}: {dot}");
                    checked += 1;
                }
            }
            assert!(checked > 100, "{checked}");
        }
    }
}

/// Visibility of thread drawings against marching the source solid's own
/// containment test: every sampled point of every candidate curve that is
/// visible by marching lies on a published visible vector, and every hidden
/// one on a hidden or a coincident visible vector.
#[test]
fn thread_visibility_agrees_with_source_containment_marching() {
    use crate::{DrawingSolid, HelicalThread, Solid};
    use multicad_core::document::EntityId;
    use multicad_rust_kernel::rounded_thread::Spec;
    let policy = GeometryPolicy::default();
    let p = 1.27;
    let a = 27.5_f64.to_radians();
    let r = p * a.cos() / (12. * (1. - a.sin()));
    let rounded = Solid::rust_rounded_thread(
        Spec {
            major_radius_mm: 3.175,
            pitch_mm: p,
            length_mm: 3.7 * p,
            depth_mm: p / (3. * a.tan()),
            crest_radius_mm: r,
            root_radius_mm: r,
            flank_angle_deg: 27.5,
            starts: 1,
            left_hand: false,
        },
        EntityId::new(),
        PlaneFrame::default(),
        policy,
    )
    .unwrap();
    let root = 1. / (4. * 3_f64.sqrt());
    let relief = Solid::rust_helical_thread_end_relief(
        HelicalThread {
            major_radius_mm: 3.,
            pitch_mm: 1.,
            length_mm: 6.6,
            depth_mm: (1. - 0.125) / (2. * 30_f64.to_radians().tan()) - root,
            crest_width_mm: 0.125,
            lower_flank_deg: 30.,
            upper_flank_deg: 30.,
            radial_taper: 0.,
            starts: 1,
            left_hand: true,
        },
        Some(root),
        0.6,
        0.1,
        EntityId::new(),
        PlaneFrame::default(),
        policy,
    )
    .unwrap();
    let s = 0.5_f64.sqrt();
    let t = (1_f64 / 6.).sqrt();
    let double = two_start(false, PlaneFrame::default());
    let bore = rounded_bore(true, PlaneFrame::default());
    let treatment = treated(false, 1, PlaneFrame::default());
    let groove = grooved(true, 2, PlaneFrame::default());
    for (solid, center) in [
        (&rounded, 2.35),
        (&relief, 3.3),
        (&double, 3.3),
        (&bore, 2.35),
        (&treatment, 3.25),
        (&groove, 4.6),
    ] {
        for [u, v] in [
            [[1., 0., 0.], [0., 0., 1.]],
            [[s, -s, 0.], [t, t, -2. * t]],
            [[1., 0., 0.], [0., 1., 0.]],
        ] {
            let frame = PlaneFrame {
                origin: [0.2, 0.1, center],
                u,
                v,
            };
            let view = View::new(frame);
            let options = DrawingProjection {
                frame,
                hidden_lines: true,
                ..Default::default()
            };
            let vectors = project_solids(&[DrawingSolid { solid, hatch: None }], options)
                .unwrap()
                .unwrap();
            let world = |q: P3| -> [f64; 3] {
                std::array::from_fn(|i| {
                    frame.origin[i] + (0..3).map(|j| view.axes[j][i] * q[j]).sum::<f64>()
                })
            };
            // Material strictly in front of a view point, by marching.
            let covered = |q: P3| {
                (1..4000).any(|i| {
                    let ahead = [q[0], q[1], q[2] + 1e-4 + 2e-3 * f64::from(i)];
                    contains(solid, world(ahead))
                        && contains(solid, world([q[0], q[1], ahead[2] - 5e-5]))
                })
            };
            let near = |role: VectorRole, point: P2| {
                vectors.iter().filter(|v| v.role == role).any(|v| {
                    v.preview_mm.windows(2).any(|w| {
                        let [[a0, a1], [b0, b1]] = [w[0], w[1]];
                        let [dx, dy] = [b0 - a0, b1 - a1];
                        let along = (((point[0] - a0) * dx + (point[1] - a1) * dy)
                            / (dx * dx + dy * dy).max(1e-300))
                        .clamp(0., 1.);
                        (point[0] - a0 - along * dx).hypot(point[1] - a1 - along * dy) < 0.02
                    })
                })
            };
            let member = families::member(solid, &view).unwrap().unwrap();
            let (mut seen, mut wrong) = (0, Vec::new());
            for candidate in member.candidates(0.0025).unwrap() {
                if candidate.kind == Kind::Smooth {
                    continue;
                }
                for k in 1..12 {
                    let tt = k as f64 / 12.;
                    let q = candidate.path.point(tt);
                    let truth = covered(q);
                    // Skip points near a change of visibility along the curve.
                    let neighbours = [tt - 0.02, tt + 0.02]
                        .map(|t| covered(candidate.path.point(t.clamp(0., 1.))));
                    if neighbours.iter().any(|n| *n != truth) {
                        continue;
                    }
                    seen += 1;
                    let point = [q[0], q[1]];
                    let ok = if truth {
                        near(VectorRole::Hidden, point) || near(VectorRole::Visible, point)
                    } else {
                        near(VectorRole::Visible, point)
                    };
                    if !ok {
                        wrong.push((
                            candidate.kind,
                            q,
                            truth,
                            near(VectorRole::Hidden, point),
                            u,
                            v,
                        ));
                    }
                }
            }
            assert!(seen > 100, "{seen}");
            assert!(
                wrong.is_empty(),
                "{} of {seen}: {:?}",
                wrong.len(),
                &wrong[..wrong.len().min(6)]
            );
        }
    }
}

/// Every published section boundary point of both thread families lies in
/// the cut plane and on the source boundary: the helicoid, a cap, the neck,
/// its far disk or the shoulder. Axial sections enclose the exact area.
#[test]
fn thread_sections_lie_on_the_source_boundary() {
    use crate::{DrawingHatch, DrawingSolid, HelicalThread, Solid};
    use multicad_core::{document::EntityId, drawings::VectorCurve, materials::HatchPattern};
    use multicad_rust_kernel::rounded_thread::Spec;
    let policy = GeometryPolicy::default();
    let p = 1.27;
    let a = 27.5_f64.to_radians();
    let r = p * a.cos() / (12. * (1. - a.sin()));
    let rounded = Solid::rust_rounded_thread(
        Spec {
            major_radius_mm: 3.175,
            pitch_mm: p,
            length_mm: 3.7 * p,
            depth_mm: p / (3. * a.tan()),
            crest_radius_mm: r,
            root_radius_mm: r,
            flank_angle_deg: 27.5,
            starts: 1,
            left_hand: true,
        },
        EntityId::new(),
        PlaneFrame::default(),
        policy,
    )
    .unwrap();
    let root = 1. / (4. * 3_f64.sqrt());
    let relief = Solid::rust_helical_thread_end_relief(
        HelicalThread {
            major_radius_mm: 3.,
            pitch_mm: 1.,
            length_mm: 6.6,
            depth_mm: (1. - 0.125) / (2. * 30_f64.to_radians().tan()) - root,
            crest_width_mm: 0.125,
            lower_flank_deg: 30.,
            upper_flank_deg: 30.,
            radial_taper: 0.,
            starts: 1,
            left_hand: false,
        },
        Some(root),
        0.6,
        0.1,
        EntityId::new(),
        PlaneFrame::default(),
        policy,
    )
    .unwrap();
    let rounded_thread = rounded.rust_rounded_thread_data().unwrap().thread;
    let (relief_data, _) = relief.rust_thread_end_relief_data().unwrap();
    let double_data = two_start(true, PlaneFrame::default())
        .rust_thread_end_relief_data()
        .unwrap()
        .0;
    let s = 0.5_f64.sqrt();
    let t = (1_f64 / 6.).sqrt();
    let frames = |center: f64| {
        [
            ([0.2, 0.1, center], [[1., 0., 0.], [0., 0., 1.]]),
            ([0., 0., center], [[s, s, 0.], [0., 0., 1.]]),
            ([0.2, 0.1, center], [[s, -s, 0.], [t, t, -2. * t]]),
            ([0.3, -0.2, center + 0.37], [[1., 0., 0.], [0., 1., 0.]]),
        ]
    };
    let double = two_start(true, PlaneFrame::default());
    let bore = rounded_bore(false, PlaneFrame::default());
    let bore_thread = *bore.rust_rounded_thread_bore_data().unwrap().bore.thread();
    let chamfer = chamfered(true, true, PlaneFrame::default());
    let chamfer_data = chamfer.rust_thread_entry_chamfer_data().unwrap().source;
    let socket_screw = screw(true, false, PlaneFrame::default());
    let hex_screw = screw(false, true, PlaneFrame::default());
    let treatment = treated(true, 2, PlaneFrame::default());
    let treatment_data = treatment.rust_internal_thread_treatment_data().unwrap();
    let groove = grooved(false, 1, PlaneFrame::default());
    let groove_data = groove.rust_internal_thread_treatment_data().unwrap();
    for (label, solid, center) in [
        ("groove", &groove, 3.25),
        ("groove tail", &groove, 5.5),
        ("rounded", &rounded, 2.35),
        ("relief", &relief, 3.3),
        ("two-start", &double, 3.3),
        ("bore", &bore, 2.35),
        ("chamfer", &chamfer, 0.8),
        ("socket screw", &socket_screw, 0.5),
        ("hex screw", &hex_screw, 1.5),
        ("treatment", &treatment, 3.25),
    ] {
        for (origin, [u, v]) in frames(center) {
            let frame = PlaneFrame { origin, u, v };
            let view = View::new(frame);
            let options = DrawingProjection {
                frame,
                section_offset_mm: Some(0.),
                section_only: true,
                ..Default::default()
            };
            let vectors = project_solids(
                &[DrawingSolid {
                    solid,
                    hatch: Some(DrawingHatch {
                        angle_rad: 0.,
                        spacing_mm: 1.,
                        line_width_mm: 0.18,
                        pattern: HatchPattern::Solid,
                    }),
                }],
                options,
            )
            .unwrap()
            .unwrap_or_else(|error| panic!("{label} {origin:?} {u:?}: {error:?}"));
            let mut checked = 0;
            for vector in vectors.iter().filter(|v| v.role == VectorRole::Section) {
                for q in &vector.preview_mm {
                    let w: [f64; 3] = std::array::from_fn(|i| {
                        frame.origin[i] + view.axes[0][i] * q[0] + view.axes[1][i] * q[1]
                    });
                    let [x, y, z] = w;
                    let radius = x.hypot(y);
                    let theta = y.atan2(x);
                    let on = |value: f64| value.abs() < 1e-6;
                    let ok = if label.ends_with("screw") {
                        let solid = if label == "socket screw" {
                            &socket_screw
                        } else {
                            &hex_screw
                        };
                        let (source, placed) = solid.rust_machine_screw_geometry().unwrap();
                        // The screw's thread frame is offset from its mounting frame.
                        let normal = placed.normal();
                        let offset: [f64; 3] = std::array::from_fn(|i| w[i] - placed.origin[i]);
                        let dot =
                            |a: [f64; 3]| a.iter().zip(offset).map(|(p, q)| p * q).sum::<f64>();
                        let [x, y, z] = [dot(placed.u), dot(placed.v), dot(normal)];
                        let radius = x.hypot(y);
                        let theta = y.atan2(x);
                        let thread = source.thread();
                        let length = thread.spec().length_mm;
                        let (top, seat) = (source.head_top_z(), source.head_seat_z());
                        let spec = source.spec();
                        let head = source.head_radius(theta);
                        let socket = source.socket_radius(theta);
                        let depth = spec.head.socket().map_or(0., |(_, depth)| depth);
                        let bound = thread.radius(theta, z.clamp(0., length));
                        let shank = spec.shank_radius_mm;
                        (on(radius - bound) && z > -1e-9 && z < length + 1e-9)
                            || (on(z - length) && radius <= bound + 1e-6)
                            || (on(z)
                                && radius >= bound - 1e-6
                                && radius
                                    <= if spec.shank_length_mm > 0. {
                                        shank
                                    } else {
                                        head
                                    } + 1e-6)
                            || (spec.shank_length_mm > 0.
                                && ((on(radius - shank) && z > seat - 1e-9 && z < 1e-9)
                                    || (on(z - seat)
                                        && radius >= shank - 1e-6
                                        && radius <= head + 1e-6)))
                            || (on(radius - head) && z > top - 1e-9 && z < seat + 1e-9)
                            || (on(z - top)
                                && radius <= head + 1e-6
                                && socket.is_none_or(|inner| radius >= inner - 1e-6))
                            || socket.is_some_and(|inner| {
                                (on(radius - inner) && z > top - 1e-9 && z < top + depth + 1e-9)
                                    || (on(z - top - depth) && radius <= inner + 1e-6)
                            })
                    } else if label == "chamfer" {
                        let material = chamfer_data.material_length_mm();
                        let top = chamfer_data.full_length_mm();
                        let neck = chamfer_data.terminal_relief().unwrap().neck_radius_mm();
                        let bound = chamfer_data.radius(theta, z.clamp(0., top));
                        (on(radius - bound) && z > -1e-9 && z < top + 1e-9)
                            || (on(z) && radius <= bound + 1e-6)
                            || (on(z - material)
                                && radius >= neck - 1e-6
                                && radius <= chamfer_data.radius(theta, material) + 1e-6)
                            || (on(z - top) && radius <= neck + 1e-6)
                    } else if label.starts_with("groove") {
                        let cavity = groove_data.material.cavity();
                        let length = cavity.full_length_mm();
                        let stock = groove_data.material.stock().circle().radius().as_f64();
                        let station = cavity.terminal_station_mm().unwrap();
                        let tail = cavity.terminal_radius_mm().unwrap();
                        let bound = cavity.radius_mm(theta, z.clamp(0., length));
                        // The thread's void radius just below the shoulder.
                        let below = cavity.radius_mm(theta, station - 1e-9);
                        (on(radius - bound) && z > -1e-9 && z < length + 1e-9)
                            || (on(radius - stock) && z > -1e-9 && z < length + 1e-9)
                            || ((on(z) || on(z - length))
                                && radius >= bound - 1e-6
                                && radius <= stock + 1e-6)
                            || (on(z - station) && radius >= below - 1e-6 && radius <= tail + 1e-6)
                    } else if label == "treatment" {
                        let cavity = treatment_data.material.cavity();
                        let length = cavity.full_length_mm();
                        let stock = treatment_data.material.stock().circle().radius().as_f64();
                        let bound = cavity.radius_mm(theta, z.clamp(0., length));
                        (on(radius - bound) && z > -1e-9 && z < length + 1e-9)
                            || (on(radius - stock) && z > -1e-9 && z < length + 1e-9)
                            || ((on(z) || on(z - length))
                                && radius >= bound - 1e-6
                                && radius <= stock + 1e-6)
                    } else if label == "bore" {
                        let bound = bore_thread.radius(theta, z.clamp(0., 3.7 * p));
                        let length = 3.7 * p;
                        (on(radius - bound) && z > -1e-9 && z < length + 1e-9)
                            || (on(radius - 5.) && z > -1e-9 && z < length + 1e-9)
                            || ((on(z) || on(z - length))
                                && radius >= bound - 1e-6
                                && radius <= 5. + 1e-6)
                    } else if label == "rounded" {
                        let bound = rounded_thread.radius(theta, z.clamp(0., 3.7 * p));
                        let length = 3.7 * p;
                        (on(radius - bound) && z > -1e-9 && z < length + 1e-9)
                            || ((on(z) || on(z - length)) && radius <= bound + 1e-6)
                    } else {
                        let relief_data = if label == "relief" {
                            relief_data
                        } else {
                            double_data
                        };
                        let thread = relief_data.material_thread();
                        let [junction, top] = [
                            relief_data.material_length_mm(),
                            relief_data.full_length_mm(),
                        ];
                        let neck = relief_data.neck_radius_mm();
                        let bound = thread.radius(theta, z.clamp(0., junction));
                        (on(radius - bound) && z > -1e-9 && z < junction + 1e-9)
                            || (on(z) && radius <= bound + 1e-6)
                            || (on(z - junction) && radius >= neck - 1e-6 && radius <= bound + 1e-6)
                            || (on(radius - neck) && z > junction - 1e-9 && z < top + 1e-9)
                            || (on(z - top) && radius <= neck + 1e-6)
                    };
                    assert!(ok, "{label} {origin:?} {u:?}: {w:?} r {radius}");
                    checked += 1;
                }
            }
            assert!(checked > 20, "{label}: {checked}");
            let area: f64 = vectors
                .iter()
                .filter_map(|v| match v.curve {
                    VectorCurve::HatchTriangle { points: [a, b, c] } => Some(
                        ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])).abs() / 2.,
                    ),
                    _ => None,
                })
                .sum();
            // Independent area of the cut by containment on a fine grid.
            let step = 0.01;
            let mut inside = 0usize;
            // The grid covers the published section with a margin.
            let points = vectors.iter().flat_map(|v| v.preview_mm.iter());
            let lo = points
                .clone()
                .fold([f64::INFINITY; 2], |m, p| [m[0].min(p[0]), m[1].min(p[1])]);
            let hi = points.fold([f64::NEG_INFINITY; 2], |m, p| {
                [m[0].max(p[0]), m[1].max(p[1])]
            });
            let first = [0, 1].map(|k| ((lo[k] - 0.5) / step).floor() as i64);
            let last = [0, 1].map(|k| ((hi[k] + 0.5) / step).ceil() as i64);
            for i in first[0]..last[0] {
                for j in first[1]..last[1] {
                    let q = [(i as f64 + 0.5) * step, (j as f64 + 0.5) * step];
                    let w: [f64; 3] = std::array::from_fn(|k| {
                        frame.origin[k] + view.axes[0][k] * q[0] + view.axes[1][k] * q[1]
                    });
                    inside += usize::from(contains(solid, w));
                }
            }
            let grid = inside as f64 * step * step;
            assert!(
                (area - grid).abs() < 0.01 * grid,
                "{label} {origin:?} {u:?}: {area} {grid}"
            );
            if u == [1., 0., 0.] && v == [0., 1., 0.] && !label.ends_with("screw") {
                // Across the axis: the constant cross-section area.
                let area: f64 = vectors
                    .iter()
                    .filter_map(|v| match v.curve {
                        VectorCurve::HatchTriangle { points: [a, b, c] } => Some(
                            ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])).abs()
                                / 2.,
                        ),
                        _ => None,
                    })
                    .sum();
                let exact = if label == "chamfer" {
                    // The axial cut at 0.37 mm above the centre lies past the
                    // cone, where the area is the thread's.
                    let thread = chamfer_data.source_thread();
                    thread.volume_mm3() / thread.spec().length_mm
                } else if label == "bore" {
                    PI * 25. - bore_thread.volume_mm3() / (3.7 * p)
                } else if label == "groove tail" {
                    // Past the shoulder: the annulus outside the tail.
                    let cavity = groove_data.material.cavity();
                    let tail = cavity.terminal_radius_mm().unwrap();
                    let stock = groove_data.material.stock().circle().radius().as_f64();
                    PI * (stock * stock - tail * tail)
                } else if label == "treatment" || label == "groove" {
                    // Past the entry cone: the stock disc less the cavity's
                    // section, pi times the mean of R^2 over one pitch.
                    let data = if label == "groove" {
                        &groove_data
                    } else {
                        &treatment_data
                    };
                    let cavity = data.material.cavity();
                    let thread = cavity.source_thread();
                    let mean: f64 = thread
                        .bands()
                        .iter()
                        .map(|band| {
                            let [a, b] = band.radius_mm;
                            (band.phase_mm[1] - band.phase_mm[0]) * (a * a + a * b + b * b) / 3.
                        })
                        .sum::<f64>()
                        / thread.spec().pitch_mm;
                    let stock = data.material.stock().circle().radius().as_f64();
                    PI * (stock * stock - mean)
                } else if label == "rounded" {
                    rounded_thread.volume_mm3() / (3.7 * p)
                } else {
                    let data = if label == "relief" {
                        relief_data
                    } else {
                        double_data
                    };
                    let thread = data.material_thread();
                    thread.volume_mm3() / thread.spec().length_mm
                };
                assert!(
                    (area - exact).abs() < 2e-3 * exact,
                    "{label}: {area} {exact}"
                );
            }
        }
    }
}

/// A two-start ISO-like thread with a terminal relief.
fn two_start(left_hand: bool, frame: PlaneFrame) -> crate::Solid {
    let root = 0.5 / (4. * 3_f64.sqrt());
    crate::Solid::rust_helical_thread_end_relief(
        crate::HelicalThread {
            major_radius_mm: 3.,
            pitch_mm: 0.5,
            length_mm: 6.6,
            depth_mm: (0.5 - 0.0625) / (2. * 30_f64.to_radians().tan()) - root,
            crest_width_mm: 0.0625,
            lower_flank_deg: 30.,
            upper_flank_deg: 30.,
            radial_taper: 0.,
            starts: 2,
            left_hand,
        },
        Some(root),
        0.6,
        0.1,
        multicad_core::document::EntityId::new(),
        frame,
        GeometryPolicy::default(),
    )
    .unwrap()
}

/// A rounded thread bore through a radius-5 stock.
fn rounded_bore(left_hand: bool, frame: PlaneFrame) -> crate::Solid {
    let p = 1.27;
    let a = 27.5_f64.to_radians();
    let r = p * a.cos() / (12. * (1. - a.sin()));
    crate::Solid::rust_rounded_thread_bore(
        multicad_rust_kernel::rounded_thread::Spec {
            major_radius_mm: 3.175,
            pitch_mm: p,
            length_mm: 3.7 * p,
            depth_mm: p / (3. * a.tan()),
            crest_radius_mm: r,
            root_radius_mm: r,
            flank_angle_deg: 27.5,
            starts: 1,
            left_hand,
        },
        5.,
        multicad_core::document::EntityId::new(),
        frame,
        GeometryPolicy::default(),
    )
    .unwrap()
}

/// An authored internal thread treatment: a radius-4.5 stock less an ISO-like
/// internal thread cavity with an entry cone, placed by `frame`.
fn treated(left_hand: bool, starts: u32, frame: PlaneFrame) -> crate::Solid {
    treated_with(left_hand, starts, frame, None)
}

/// As `treated`, with a terminal groove 1.2 mm wide and 0.4 mm deep: past
/// its station at 5.3 mm the void is the radius-3.4 tail cylinder.
fn grooved(left_hand: bool, starts: u32, frame: PlaneFrame) -> crate::Solid {
    treated_with(
        left_hand,
        starts,
        frame,
        Some(
            multicad_rust_kernel::authored_internal_treatment::Terminal {
                width_mm: 1.2,
                radial_depth_mm: 0.4,
            },
        ),
    )
}

fn treated_with(
    left_hand: bool,
    starts: u32,
    frame: PlaneFrame,
    terminal: Option<multicad_rust_kernel::authored_internal_treatment::Terminal>,
) -> crate::Solid {
    use multicad_rust_kernel::{Circle, Rational, cylinder::AxialCylinder};
    let policy = GeometryPolicy::default();
    let zero = Rational::new(0, 1).unwrap();
    let stock = AxialCylinder::new(
        7,
        Circle::new([zero; 2], Rational::new(9, 2).unwrap()).unwrap(),
        zero,
        Rational::new(13, 2).unwrap(),
    )
    .unwrap();
    let pitch = 0.5 * f64::from(starts);
    let treated = crate::Solid::from_axial_cylinder(stock, PlaneFrame::default(), policy)
        .unwrap()
        .rust_cut_authored_internal_thread_treatments_with_interrupt(
            crate::InternalThreadTreatmentCut {
                profile: crate::HelicalThread {
                    major_radius_mm: 3.,
                    pitch_mm: pitch,
                    length_mm: 6.5,
                    depth_mm: 5. * 3_f64.sqrt() / 16. * pitch,
                    crest_width_mm: pitch / 8.,
                    lower_flank_deg: 30.,
                    upper_flank_deg: 30.,
                    radial_taper: 0.,
                    starts,
                    left_hand,
                },
                intent: multicad_rust_kernel::authored_internal_treatment::Intent {
                    entry_depth_mm: Some(0.8),
                    terminal,
                },
                minimum_wall_mm: 0.5,
                frame: PlaneFrame::default(),
            },
            multicad_core::document::EntityId::new(),
            policy,
            None,
        )
        .unwrap();
    let n = frame.normal();
    treated
        .transformed(crate::RigidTransform {
            rotation: std::array::from_fn(|i| [frame.u[i], frame.v[i], n[i]]),
            translation_mm: frame.origin,
        })
        .unwrap()
}

/// Source containment: the thread solids' own test, or for a treated bore
/// the actual stock less the authored cavity's nominal evaluation.
fn contains(solid: &crate::Solid, w: P3) -> bool {
    match solid.rust_internal_thread_treatment_data() {
        Some(data) => {
            let frame = data.frame;
            let offset: P3 = std::array::from_fn(|i| w[i] - frame.origin[i]);
            let local = [frame.u, frame.v, frame.normal()].map(|axis| dot3(axis, offset));
            let stock = data.material.stock().circle().radius().as_f64();
            let length = data.material.cavity().full_length_mm();
            local[0].hypot(local[1]) < stock
                && local[2] > 0.
                && local[2] < length
                && !data.material.cavity().contains_void(local, 0.)
        }
        None => solid.rust_thread_contains(w).unwrap(),
    }
}

/// An independent depth-map oracle for thread drawings: the first material
/// along each pixel's sight line is found by marching the source solid's
/// containment test alone. Wherever that depth jumps between neighbouring
/// pixels, or a pixel's sight line starts or stops meeting material, an
/// occluding contour or edge is seen there, and a published visible vector
/// must pass within two and a half pixels.
#[test]
fn thread_occluding_contours_agree_with_a_depth_map_oracle() {
    use crate::DrawingSolid;
    let relief = two_start(false, PlaneFrame::default());
    let bore = rounded_bore(true, PlaneFrame::default());
    let p = 1.27;
    let a = 27.5_f64.to_radians();
    let r = p * a.cos() / (12. * (1. - a.sin()));
    let rounded = crate::Solid::rust_rounded_thread(
        multicad_rust_kernel::rounded_thread::Spec {
            major_radius_mm: 3.175,
            pitch_mm: p,
            length_mm: 3.7 * p,
            depth_mm: p / (3. * a.tan()),
            crest_radius_mm: r,
            root_radius_mm: r,
            flank_angle_deg: 27.5,
            starts: 1,
            left_hand: false,
        },
        multicad_core::document::EntityId::new(),
        PlaneFrame::default(),
        GeometryPolicy::default(),
    )
    .unwrap();
    let s = 0.5_f64.sqrt();
    let t = (1_f64 / 6.).sqrt();
    let treatment = treated(false, 2, PlaneFrame::default());
    let groove = grooved(false, 1, PlaneFrame::default());
    for (label, solid, center) in [
        ("rounded", &rounded, 2.35),
        ("relief", &relief, 3.3),
        ("bore", &bore, 2.35),
        ("treatment", &treatment, 3.25),
        ("groove", &groove, 4.6),
    ] {
        for [u, v] in [[[1., 0., 0.], [0., 0., 1.]], [[s, -s, 0.], [t, t, -2. * t]]] {
            let frame = PlaneFrame {
                origin: [0.1, -0.05, center],
                u,
                v,
            };
            let view = View::new(frame);
            let options = DrawingProjection {
                frame,
                ..Default::default()
            };
            let vectors = project_solids(&[DrawingSolid { solid, hatch: None }], options)
                .unwrap()
                .unwrap();
            let world = |q: P3| -> [f64; 3] {
                std::array::from_fn(|i| {
                    frame.origin[i] + (0..3).map(|j| view.axes[j][i] * q[j]).sum::<f64>()
                })
            };
            let inside = |q: P3| contains(solid, world(q));
            // First material from the viewer, refined by bisection.
            let depth = |x: f64, y: f64| -> Option<f64> {
                let (top, bottom, step) = (8., -8., 0.01);
                let mut z = top;
                while z > bottom {
                    let next = z - step;
                    if inside([x, y, next]) {
                        let [mut hi, mut lo] = [z, next];
                        for _ in 0..30 {
                            let mid = 0.5 * (hi + lo);
                            if inside([x, y, mid]) {
                                lo = mid
                            } else {
                                hi = mid
                            }
                        }
                        return Some(lo);
                    }
                    z = next;
                }
                None
            };
            let pixel = 0.05;
            let count = 200;
            let grid: Vec<Vec<Option<f64>>> = (0..count)
                .map(|i| {
                    (0..count)
                        .map(|j| {
                            depth(
                                -5. + pixel * (i as f64 + 0.5),
                                -5. + pixel * (j as f64 + 0.5),
                            )
                        })
                        .collect()
                })
                .collect();
            let near_visible = |x: f64, y: f64| {
                vectors
                    .iter()
                    .filter(|v| v.role == VectorRole::Visible)
                    .any(|v| {
                        v.preview_mm.windows(2).any(|w| {
                            let [[a0, a1], [b0, b1]] = [w[0], w[1]];
                            let [dx, dy] = [b0 - a0, b1 - a1];
                            let along = (((x - a0) * dx + (y - a1) * dy)
                                / (dx * dx + dy * dy).max(1e-300))
                            .clamp(0., 1.);
                            (x - a0 - along * dx).hypot(y - a1 - along * dy) < 2.5 * pixel
                        })
                    })
            };
            let (mut seen, mut missing) = (0, Vec::new());
            for i in 0..count - 1 {
                for j in 0..count - 1 {
                    for (k, l) in [(i + 1, j), (i, j + 1)] {
                        let jump = match (grid[i][j], grid[k][l]) {
                            (Some(a), Some(b)) => (a - b).abs() > 0.3,
                            (None, None) => false,
                            _ => true,
                        };
                        if !jump {
                            continue;
                        }
                        seen += 1;
                        let x = -5. + pixel * (0.5 * (i + k) as f64 + 0.5);
                        let y = -5. + pixel * (0.5 * (j + l) as f64 + 0.5);
                        if !near_visible(x, y) {
                            missing.push([x, y]);
                        }
                    }
                }
            }
            assert!(seen > 100, "{label}: {seen}");
            assert!(
                missing.is_empty(),
                "{label} {u:?}: {} of {seen}: {:?}",
                missing.len(),
                &missing[..missing.len().min(6)]
            );
        }
    }
}

/// An ISO-like thread with an entry chamfer and, optionally, a terminal relief.
fn chamfered(left_hand: bool, terminal: bool, frame: PlaneFrame) -> crate::Solid {
    let root = 1. / (4. * 3_f64.sqrt());
    crate::Solid::rust_helical_thread_entry_chamfer(
        crate::HelicalThread {
            major_radius_mm: 3.,
            pitch_mm: 1.,
            length_mm: 6.6,
            depth_mm: (1. - 0.125) / (2. * 30_f64.to_radians().tan()) - root,
            crest_width_mm: 0.125,
            lower_flank_deg: 30.,
            upper_flank_deg: 30.,
            radial_taper: 0.,
            starts: 1,
            left_hand,
        },
        Some(root),
        0.3,
        terminal.then_some(multicad_rust_kernel::thread_end_relief::Intent {
            width_mm: 0.6,
            radial_depth_mm: 0.1,
        }),
        multicad_core::document::EntityId::new(),
        frame,
        GeometryPolicy::default(),
    )
    .unwrap()
}

/// An ISO-like machine screw: a plain shank under a socket cap head, or a
/// hex head seated directly on the thread.
fn screw(socket: bool, left_hand: bool, frame: PlaneFrame) -> crate::Solid {
    use multicad_rust_kernel::headed_screw::{ScrewHead, ScrewSpec};
    let root = 1. / (4. * 3_f64.sqrt());
    crate::Solid::rust_machine_screw(
        crate::HelicalThread {
            major_radius_mm: 3.,
            pitch_mm: 1.,
            length_mm: 6.6,
            depth_mm: (1. - 0.125) / (2. * 30_f64.to_radians().tan()) - root,
            crest_width_mm: 0.125,
            lower_flank_deg: 30.,
            upper_flank_deg: 30.,
            radial_taper: 0.,
            starts: 1,
            left_hand,
        },
        Some(root),
        if socket {
            ScrewSpec {
                head: ScrewHead::SocketCap {
                    radius_mm: 4.5,
                    height_mm: 3.,
                    socket_across_flats_mm: 3.,
                    socket_depth_mm: 1.5,
                },
                shank_length_mm: 2.,
                shank_radius_mm: 3.,
            }
        } else {
            ScrewSpec {
                head: ScrewHead::ExternalHex {
                    across_flats_mm: 10.,
                    height_mm: 4.,
                },
                shank_length_mm: 0.,
                shank_radius_mm: 3.,
            }
        },
        multicad_core::document::EntityId::new(),
        frame,
        GeometryPolicy::default(),
    )
    .unwrap()
}

/// Side cuts just off the axis of each section-capable thread family close
/// with a valid hatch of the area found by containment on a grid. There the
/// plane meets the meridian lines at a grazing angle, so a profile point
/// moves by an angle's error times the radius squared over the offset, and
/// a pair of arc touches nearly coincide.
#[test]
fn thread_near_axis_side_cuts_close_with_the_grid_area() {
    use crate::{DrawingHatch, DrawingSolid};
    use multicad_core::{drawings::VectorCurve, materials::HatchPattern};
    for (label, solid, center) in [
        ("two-start", two_start(true, PlaneFrame::default()), 3.3),
        ("bore", rounded_bore(false, PlaneFrame::default()), 2.35),
        ("chamfer", chamfered(true, true, PlaneFrame::default()), 0.8),
        (
            "socket screw",
            screw(true, false, PlaneFrame::default()),
            0.5,
        ),
        ("treatment", treated(true, 2, PlaneFrame::default()), 3.25),
        ("groove", grooved(false, 1, PlaneFrame::default()), 3.25),
    ] {
        for offset in [1e-6_f64, 1e-4, 1e-3, 1e-2] {
            for angle in [0., 1., 2.6] {
                let (sin, cos) = f64::sin_cos(angle);
                let frame = PlaneFrame {
                    origin: [-offset * sin, offset * cos, center],
                    u: [cos, sin, 0.],
                    v: [0., 0., 1.],
                };
                let view = View::new(frame);
                let options = DrawingProjection {
                    frame,
                    section_offset_mm: Some(0.),
                    section_only: true,
                    ..Default::default()
                };
                let vectors = project_solids(
                    &[DrawingSolid {
                        solid: &solid,
                        hatch: Some(DrawingHatch {
                            angle_rad: 0.,
                            spacing_mm: 1.,
                            line_width_mm: 0.18,
                            pattern: HatchPattern::Solid,
                        }),
                    }],
                    options,
                )
                .unwrap()
                .unwrap_or_else(|error| panic!("{label} {offset} {angle}: {error:?}"));
                let area: f64 = vectors
                    .iter()
                    .filter_map(|v| match v.curve {
                        VectorCurve::HatchTriangle { points: [a, b, c] } => Some(
                            ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])).abs()
                                / 2.,
                        ),
                        _ => None,
                    })
                    .sum();
                let step = 0.02;
                let points = vectors.iter().flat_map(|v| v.preview_mm.iter());
                let lo = points
                    .clone()
                    .fold([f64::INFINITY; 2], |m, p| [m[0].min(p[0]), m[1].min(p[1])]);
                let hi = points.fold([f64::NEG_INFINITY; 2], |m, p| {
                    [m[0].max(p[0]), m[1].max(p[1])]
                });
                let first = [0, 1].map(|k| ((lo[k] - 0.5) / step).floor() as i64);
                let last = [0, 1].map(|k| ((hi[k] + 0.5) / step).ceil() as i64);
                let mut inside = 0usize;
                for i in first[0]..last[0] {
                    for j in first[1]..last[1] {
                        let q = [(i as f64 + 0.5) * step, (j as f64 + 0.5) * step];
                        let w: [f64; 3] = std::array::from_fn(|k| {
                            frame.origin[k] + view.axes[0][k] * q[0] + view.axes[1][k] * q[1]
                        });
                        inside += usize::from(contains(&solid, w));
                    }
                }
                let grid = inside as f64 * step * step;
                assert!(
                    (area - grid).abs() < 0.01 * grid,
                    "{label} {offset} {angle}: {area} {grid}"
                );
            }
        }
    }
}
