//! Certification of the pipe thread member against the source's own
//! certified radial bounds: material along lines, visibility by marching,
//! contour conditions, section boundaries and areas, and a depth-map oracle.
use super::paths::dot3;
use super::*;
use crate::{DrawingHatch, DrawingSolid, RigidTransform, Solid};
use multicad_core::{
    document::EntityId,
    drawings::{VectorCurve, VectorRole},
    materials::HatchPattern,
};
use multicad_rust_kernel::axial_pipe_thread::authored::{AuthoredSection, AuthoredSpec};
use std::f64::consts::TAU;

/// A BSPT-like rounded pipe thread (the inventory's 1/4" row).
fn bspt(frame: PlaneFrame) -> Solid {
    let p = 25.4 / 28.;
    build(
        AuthoredSpec {
            major_radius_mm: 9.728 / 2. - 4. / 32.,
            pitch_mm: p,
            length_mm: 6.5,
            radial_taper: 1. / 32.,
            section: AuthoredSection::BsptRounded {
                depth_mm: p / (3. * 27.5_f64.to_radians().tan()),
            },
            crest_radius_limit_mm: None,
        },
        frame,
    )
}

/// An NPT-like flat-crested pipe thread, tapering inward along its axis.
fn npt(frame: PlaneFrame, taper: f64) -> Solid {
    let p = 25.4 / 18.;
    build(
        AuthoredSpec {
            major_radius_mm: 6.8,
            pitch_mm: p,
            length_mm: 5.5,
            radial_taper: taper,
            section: AuthoredSection::NptFlat {
                depth_mm: 0.8 * p,
                flat_width_mm: 0.038 * p,
            },
            crest_radius_limit_mm: None,
        },
        frame,
    )
}

fn build(spec: AuthoredSpec, frame: PlaneFrame) -> Solid {
    Solid::rust_authored_axial_pipe_thread_with_interrupt(
        spec,
        EntityId::new(),
        frame,
        GeometryPolicy::default(),
        None,
    )
    .unwrap()
}

fn posed() -> PlaneFrame {
    PlaneFrame {
        origin: [1.5, -2., 0.5],
        u: [0.8, 0.6, 0.],
        v: [-0.36, 0.48, 0.8],
    }
}

/// The source's construction-chart point of a world point.
fn chart(solid: &Solid, w: P3) -> P3 {
    let (source, occurrence) = solid.rust_authored_axial_pipe_thread_source().unwrap();
    let invert = |frame: PlaneFrame, p: P3| -> P3 {
        let d: P3 = std::array::from_fn(|i| p[i] - frame.origin[i]);
        [frame.u, frame.v, frame.normal()].map(|axis| dot3(axis, d))
    };
    invert(source.construction_frame, invert(occurrence, w))
}

/// Containment by the source's certified radial bounds alone: `None` where
/// the point is within the bounds' width of the boundary.
fn contains(solid: &Solid, w: P3) -> Option<bool> {
    let (source, _) = solid.rust_authored_axial_pipe_thread_source().unwrap();
    let spec = source.transport.spec();
    let [x, y, z] = chart(solid, w);
    if z <= 0. || z >= spec.length_mm {
        return Some(false);
    }
    let k = spec.pitch_mm / TAU;
    let q = (z - k * y.atan2(x)).rem_euclid(spec.pitch_mm);
    let [lo, hi] = source.transport.radial_bounds([q, q], [z, z]).ok()?;
    let r = x.hypot(y);
    if r < lo {
        Some(true)
    } else if r > hi {
        Some(false)
    } else {
        None
    }
}

fn inside(solid: &Solid, w: P3) -> bool {
    contains(solid, w).unwrap_or(false)
}

fn solids() -> Vec<(&'static str, Solid)> {
    let moved = bspt(PlaneFrame::default())
        .transformed(RigidTransform {
            rotation: [[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]],
            translation_mm: [0.3, -0.2, 0.1],
        })
        .unwrap();
    vec![
        ("bspt", bspt(PlaneFrame::default())),
        ("npt", npt(PlaneFrame::default(), 1. / 32.)),
        ("npt inward", npt(PlaneFrame::default(), -1. / 32.)),
        ("bspt posed", bspt(posed())),
        ("bspt moved", moved),
    ]
}

/// Material along lines agrees with marching the source's radial bounds,
/// in generic directions, along the axis and in a cap plane.
#[test]
fn pipe_material_along_lines_agrees_with_source_bounds() {
    let view = View::new(PlaneFrame::default());
    for (label, solid) in solids() {
        let (source, occurrence) = solid.rust_authored_axial_pipe_thread_source().unwrap();
        let construction = source.construction_frame;
        let place = |frame: PlaneFrame, p: P3, point: bool| -> P3 {
            let n = frame.normal();
            std::array::from_fn(|i| {
                (if point { frame.origin[i] } else { 0. })
                    + p[0] * frame.u[i]
                    + p[1] * frame.v[i]
                    + p[2] * n[i]
            })
        };
        let world = |l: P3| place(occurrence, place(construction, l, true), true);
        let direction = |l: P3| place(occurrence, place(construction, l, false), false);
        let mut lines = Vec::new();
        for (k, d) in [
            [0.3, -0.5, 0.81],
            [1., 0.2, -0.1],
            [0.05, 0.99, 0.1],
            [-0.7, 0.7, 0.14],
        ]
        .into_iter()
        .enumerate()
        {
            let length = dot3(d, d).sqrt();
            let d = d.map(|v| v / length);
            for offset in [0., 3.9, 4.6, 5.9, 6.7] {
                let angle = k as f64 * 0.7;
                lines.push((
                    world([offset * angle.cos(), offset * angle.sin(), 3.]),
                    direction(d),
                ));
            }
        }
        for r in [1., 4.2, 4.5, 5.0, 6.0, 6.5, 6.9] {
            lines.push((world([r, 0.3, -1.]), direction([0., 0., 1.])));
        }
        lines.push((world([-8., 0.37, 1e-3]), direction([1., 0., 0.])));
        lines.push((world([-8., -0.41, 5.5 - 1e-3]), direction([1., 0., 0.])));
        let member = families::member(&solid, &view).unwrap().unwrap();
        let mut covered = 0;
        for &(start, d) in &lines {
            let spans = member.line(start, d).unwrap().unwrap();
            let steps = 16_000;
            let mut mismatches = Vec::new();
            for i in 0..=steps {
                let t = -16. + 32. * i as f64 / steps as f64;
                if spans.iter().flatten().any(|end| (end - t).abs() < 1e-4) {
                    continue;
                }
                let Some(truth) = contains(&solid, std::array::from_fn(|i| start[i] + t * d[i]))
                else {
                    continue;
                };
                let claimed = spans.iter().any(|[a, b]| *a < t && t < *b);
                if claimed != truth {
                    mismatches.push(t);
                }
            }
            assert!(
                mismatches.is_empty(),
                "{label} {start:?} {d:?}: {:?} {spans:?}",
                &mismatches[..mismatches.len().min(4)]
            );
            covered += usize::from(!spans.is_empty());
        }
        assert!(covered >= 3 * lines.len() / 4, "{label}: {covered}");
    }
}

fn frames(center: f64) -> Vec<PlaneFrame> {
    let s = 0.5_f64.sqrt();
    let t = (1_f64 / 6.).sqrt();
    vec![
        PlaneFrame {
            origin: [0.2, 0.1, center],
            u: [1., 0., 0.],
            v: [0., 0., 1.],
        },
        PlaneFrame {
            origin: [0., 0., center],
            u: [s, s, 0.],
            v: [0., 0., 1.],
        },
        PlaneFrame {
            origin: [0.2, 0.1, center],
            u: [s, -s, 0.],
            v: [t, t, -2. * t],
        },
        PlaneFrame {
            origin: [0.3, -0.2, center + 0.37],
            u: [1., 0., 0.],
            v: [0., 1., 0.],
        },
    ]
}

fn world(view: &View, frame: PlaneFrame, q: P3) -> P3 {
    std::array::from_fn(|i| frame.origin[i] + (0..3).map(|j| view.axes[j][i] * q[j]).sum::<f64>())
}

fn near(vectors: &[ProjectedVector], role: Option<VectorRole>, point: P2, reach: f64) -> bool {
    vectors
        .iter()
        .filter(|v| role.is_none_or(|role| v.role == role))
        .any(|v| {
            v.preview_mm.windows(2).any(|w| {
                let [[a0, a1], [b0, b1]] = [w[0], w[1]];
                let [dx, dy] = [b0 - a0, b1 - a1];
                let along = (((point[0] - a0) * dx + (point[1] - a1) * dy)
                    / (dx * dx + dy * dy).max(1e-300))
                .clamp(0., 1.);
                (point[0] - a0 - along * dx).hypot(point[1] - a1 - along * dy) < reach
            })
        })
}

/// Every outline candidate lies on the thread surface with the surface
/// normal, by central differences of the source's radial bounds,
/// perpendicular to the sight.
#[test]
fn pipe_outlines_satisfy_the_contour_condition() {
    for (label, solid) in solids().into_iter().take(3) {
        let (source, _) = solid.rust_authored_axial_pipe_thread_source().unwrap();
        let spec = source.transport.spec();
        let k = spec.pitch_mm / TAU;
        let radius = |theta: f64, z: f64| {
            let q = (z - k * theta).rem_euclid(spec.pitch_mm);
            let [lo, hi] = source.transport.radial_bounds([q, q], [z, z]).unwrap();
            0.5 * (lo + hi)
        };
        for frame in frames(3.).into_iter().take(3) {
            let view = View::new(frame);
            let toward = view.axes[2];
            let member = families::member(&solid, &view).unwrap().unwrap();
            let mut checked = 0;
            for candidate in member.candidates(0.0025).unwrap() {
                if candidate.kind != Kind::Outline {
                    continue;
                }
                for i in 1..8 {
                    let w = world(&view, frame, candidate.path.point(i as f64 / 8.));
                    let [x, y, z] = w;
                    let theta = y.atan2(x);
                    let r = x.hypot(y);
                    let bound = radius(theta, z);
                    assert!((r - bound).abs() < 1e-8, "{label}: {w:?} {r} {bound}");
                    // G = r - F(theta, z): grad in Cartesian coordinates.
                    let h = 1e-6;
                    let df_theta = (radius(theta + h, z) - radius(theta - h, z)) / (2. * h);
                    let df_z = (radius(theta, z + h) - radius(theta, z - h)) / (2. * h);
                    let (sin, cos) = theta.sin_cos();
                    let normal = [cos + sin * df_theta / r, sin - cos * df_theta / r, -df_z];
                    let length = dot3(normal, normal).sqrt();
                    let dot = dot3(normal, toward) / length;
                    assert!(dot.abs() < 1e-5, "{label} {:?}: {w:?} {dot}", frame.u);
                    checked += 1;
                }
            }
            assert!(checked > 100, "{label}: {checked}");
        }
    }
}

/// Visibility by marching the source's radial bounds: every sampled point
/// of every non-smooth candidate that is visible lies on a visible vector,
/// every hidden one on a hidden or coincident visible vector.
#[test]
fn pipe_visibility_agrees_with_source_bounds_marching() {
    for (label, solid) in solids().into_iter().take(3) {
        let all = frames(3.);
        for frame in [all[0], all[2], all[3]] {
            let view = View::new(frame);
            let options = DrawingProjection {
                frame,
                hidden_lines: true,
                ..Default::default()
            };
            let vectors = project_solids(
                &[DrawingSolid {
                    solid: &solid,
                    hatch: None,
                }],
                options,
            )
            .unwrap()
            .unwrap();
            let covered = |q: P3| {
                (1..6000).any(|i| {
                    let ahead = [q[0], q[1], q[2] + 1e-4 + 2e-3 * f64::from(i)];
                    inside(&solid, world(&view, frame, ahead))
                        && inside(&solid, world(&view, frame, [q[0], q[1], ahead[2] - 5e-5]))
                })
            };
            let member = families::member(&solid, &view).unwrap().unwrap();
            let (mut seen, mut wrong) = (0, Vec::new());
            for candidate in member.candidates(0.0025).unwrap() {
                if candidate.kind == Kind::Smooth {
                    continue;
                }
                for i in 1..12 {
                    let t = i as f64 / 12.;
                    let q = candidate.path.point(t);
                    let truth = covered(q);
                    let neighbours = [t - 0.02, t + 0.02]
                        .map(|t| covered(candidate.path.point(t.clamp(0., 1.))));
                    if neighbours.iter().any(|n| *n != truth) {
                        continue;
                    }
                    seen += 1;
                    let point = [q[0], q[1]];
                    let ok = if truth {
                        near(&vectors, None, point, 0.02)
                    } else {
                        near(&vectors, Some(VectorRole::Visible), point, 0.02)
                    };
                    if !ok {
                        wrong.push((candidate.kind, q, truth));
                    }
                }
            }
            assert!(seen > 100, "{label}: {seen}");
            assert!(
                wrong.is_empty(),
                "{label} {:?}: {} of {seen}: {:?}",
                frame.u,
                wrong.len(),
                &wrong[..wrong.len().min(6)]
            );
        }
    }
}

/// Every published section boundary point lies in the cut plane and on the
/// source boundary (the thread surface within its bounds' width, or a cap
/// within the surface), and the hatch area matches a containment grid.
#[test]
fn pipe_sections_lie_on_the_source_boundary() {
    for (label, solid) in solids() {
        let (source, _) = solid.rust_authored_axial_pipe_thread_source().unwrap();
        let spec = source.transport.spec();
        let k = spec.pitch_mm / TAU;
        let bounds = |theta: f64, z: f64| {
            let z = z.clamp(0., spec.length_mm);
            let q = (z - k * theta).rem_euclid(spec.pitch_mm);
            source.transport.radial_bounds([q, q], [z, z]).unwrap()
        };
        for frame in frames(3.) {
            // The test frames are given in the construction chart.
            let (_, occurrence) = solid.rust_authored_axial_pipe_thread_source().unwrap();
            let construction = source.construction_frame;
            let place = |f: PlaneFrame, p: P3, point: bool| -> P3 {
                let n = f.normal();
                std::array::from_fn(|i| {
                    (if point { f.origin[i] } else { 0. })
                        + p[0] * f.u[i]
                        + p[1] * f.v[i]
                        + p[2] * n[i]
                })
            };
            let to_world =
                |p: P3, point: bool| place(occurrence, place(construction, p, point), point);
            let frame = PlaneFrame {
                origin: to_world(frame.origin, true),
                u: to_world(frame.u, false),
                v: to_world(frame.v, false),
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
            .unwrap_or_else(|error| panic!("{label} {:?}: {error:?}", frame.u));
            let mut checked = 0;
            for vector in vectors.iter().filter(|v| v.role == VectorRole::Section) {
                for q in &vector.preview_mm {
                    let w = world(&view, frame, [q[0], q[1], 0.]);
                    let [x, y, z] = chart(&solid, w);
                    let radius = x.hypot(y);
                    let theta = y.atan2(x);
                    let [lo, hi] = bounds(theta, z);
                    let on = |value: f64| value.abs() < 1e-6;
                    let ok = (radius > lo - 1e-6
                        && radius < hi + 1e-6
                        && z > -1e-9
                        && z < spec.length_mm + 1e-9)
                        || ((on(z) || on(z - spec.length_mm)) && radius <= hi + 1e-6);
                    assert!(
                        ok,
                        "{label} {:?}: {w:?} r {radius} [{lo}, {hi}] z {z}",
                        frame.u
                    );
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
            let step = 0.01;
            let points = vectors.iter().flat_map(|v| v.preview_mm.iter());
            let lo = points
                .clone()
                .fold([f64::INFINITY; 2], |m, p| [m[0].min(p[0]), m[1].min(p[1])]);
            let hi = points.fold([f64::NEG_INFINITY; 2], |m, p| {
                [m[0].max(p[0]), m[1].max(p[1])]
            });
            let first = [0, 1].map(|i| ((lo[i] - 0.5) / step).floor() as i64);
            let last = [0, 1].map(|i| ((hi[i] + 0.5) / step).ceil() as i64);
            let mut count = 0usize;
            for i in first[0]..last[0] {
                for j in first[1]..last[1] {
                    let q = [(i as f64 + 0.5) * step, (j as f64 + 0.5) * step, 0.];
                    count += usize::from(inside(&solid, world(&view, frame, q)));
                }
            }
            let grid = count as f64 * step * step;
            assert!(
                (area - grid).abs() < 0.01 * grid,
                "{label} {:?}: {area} {grid}",
                frame.u
            );
        }
    }
}

/// An independent depth-map oracle: wherever the first material along
/// neighbouring pixels' sight lines jumps, or starts or stops, a published
/// visible vector passes within two and a half pixels.
#[test]
fn pipe_occluding_contours_agree_with_a_depth_map_oracle() {
    for (label, solid) in solids().into_iter().take(2) {
        let all = frames(3.);
        for frame in [all[0], all[2]] {
            let view = View::new(frame);
            let options = DrawingProjection {
                frame,
                ..Default::default()
            };
            let vectors = project_solids(
                &[DrawingSolid {
                    solid: &solid,
                    hatch: None,
                }],
                options,
            )
            .unwrap()
            .unwrap();
            let solid_at = |q: P3| inside(&solid, world(&view, frame, q));
            let depth = |x: f64, y: f64| -> Option<f64> {
                let (top, bottom, step) = (9., -9., 0.01);
                let mut z = top;
                while z > bottom {
                    let next = z - step;
                    if solid_at([x, y, next]) {
                        let [mut hi, mut lo] = [z, next];
                        for _ in 0..30 {
                            let mid = 0.5 * (hi + lo);
                            if solid_at([x, y, mid]) {
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
            let count = if label == "bspt" { 200 } else { 280 };
            let half = 0.5 * pixel * count as f64;
            let grid: Vec<Vec<Option<f64>>> = (0..count)
                .map(|i| {
                    (0..count)
                        .map(|j| {
                            depth(
                                -half + pixel * (i as f64 + 0.5),
                                -half + pixel * (j as f64 + 0.5),
                            )
                        })
                        .collect()
                })
                .collect();
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
                        // A steep visible flank also steps by more than the
                        // threshold between pixels, but its step shrinks with
                        // the spacing; an occluding contour's does not. Keep
                        // only steps that survive a 64-fold refinement.
                        let at = |i: usize, f: f64| -half + pixel * (i as f64 + 0.5) + f;
                        let samples: Vec<Option<f64>> = (0..=64)
                            .map(|m| {
                                let f = m as f64 / 64.;
                                depth(
                                    at(i, f * pixel * (k - i) as f64),
                                    at(j, f * pixel * (l - j) as f64),
                                )
                            })
                            .collect();
                        let step = samples.windows(2).any(|w| match (w[0], w[1]) {
                            (Some(a), Some(b)) => (a - b).abs() > 0.15,
                            (None, None) => false,
                            _ => true,
                        });
                        if !step {
                            continue;
                        }
                        seen += 1;
                        let x = -half + pixel * (0.5 * (i + k) as f64 + 0.5);
                        let y = -half + pixel * (0.5 * (j + l) as f64 + 0.5);
                        if !near(&vectors, Some(VectorRole::Visible), [x, y], 2.5 * pixel) {
                            missing.push([x, y]);
                        }
                    }
                }
            }
            assert!(seen > 100, "{label}: {seen}");
            assert!(
                missing.is_empty(),
                "{label} {:?}: {} of {seen}: {:?}",
                frame.u,
                missing.len(),
                &missing[..missing.len().min(6)]
            );
        }
    }
}

/// The hatch area of the side cut through the plane at `offset` from the
/// axis whose horizontal direction is at `angle`, and the area of the same
/// cut by containment on a grid of `step`.
fn side_cut(solid: &Solid, offset: f64, angle: f64, step: f64) -> (f64, f64) {
    let (sin, cos) = angle.sin_cos();
    let frame = PlaneFrame {
        origin: [-offset * sin, offset * cos, 3.],
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
    .unwrap_or_else(|error| panic!("{offset} {angle}: {error:?}"));
    let area: f64 = vectors
        .iter()
        .filter_map(|v| match v.curve {
            VectorCurve::HatchTriangle { points: [a, b, c] } => {
                Some(((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])).abs() / 2.)
            }
            _ => None,
        })
        .sum();
    let points = vectors.iter().flat_map(|v| v.preview_mm.iter());
    let lo = points
        .clone()
        .fold([f64::INFINITY; 2], |m, p| [m[0].min(p[0]), m[1].min(p[1])]);
    let hi = points.fold([f64::NEG_INFINITY; 2], |m, p| {
        [m[0].max(p[0]), m[1].max(p[1])]
    });
    let first = [0, 1].map(|i| ((lo[i] - 0.5) / step).floor() as i64);
    let last = [0, 1].map(|i| ((hi[i] + 0.5) / step).ceil() as i64);
    let mut count = 0usize;
    for i in first[0]..last[0] {
        for j in first[1]..last[1] {
            let q = [(i as f64 + 0.5) * step, (j as f64 + 0.5) * step, 0.];
            count += usize::from(inside(solid, world(&view, frame, q)));
        }
    }
    (area, count as f64 * step * step)
}

/// Side cuts just off the axis close with a valid hatch of the grid's area.
/// There the plane meets the meridian lines at a grazing angle, so a profile
/// point moves by an angle's error times the radius squared over the offset.
/// The first cut was refused: its crest arc touches the meridian line
/// 3e-5 mm above the cap, where the branch's two roots meet; turning the cut
/// moves that touch up and down the crest relative to the cap.
#[test]
fn pipe_near_axis_side_cuts_close_with_the_grid_area() {
    let (area, grid) = side_cut(&bspt(PlaneFrame::default()), 1e-3, 0., 0.01);
    assert!((area - grid).abs() < 0.01 * grid, "reported: {area} {grid}");
    for (label, solid) in solids().into_iter().take(3) {
        for offset in [1e-6, 1e-5, 1e-4, 1e-3, 1e-2] {
            for angle in [0., 0.07, 1., 2.513, 2.583, 4.4] {
                let (area, grid) = side_cut(&solid, offset, angle, 0.02);
                assert!(
                    (area - grid).abs() < 0.01 * grid,
                    "{label} {offset} {angle}: {area} {grid}"
                );
            }
        }
    }
}
