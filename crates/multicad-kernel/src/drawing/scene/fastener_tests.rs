//! Certification of the plain helical thread, internal thread and hex nut
//! members against the kernel solids' own containment and exact source
//! descriptors, in the manner of the other thread members (`tests.rs`).
use super::paths::dot3;
use super::*;
use crate::{DrawingHatch, DrawingSolid, HelicalThread, Solid};
use multicad_core::{document::EntityId, drawings::VectorCurve, materials::HatchPattern};
use multicad_precision::PortableMath as _;
use std::f64::consts::{PI, TAU};

fn iso(length_mm: f64, left_hand: bool, starts: u32, pitch: f64) -> HelicalThread {
    HelicalThread {
        major_radius_mm: 3.,
        pitch_mm: pitch,
        length_mm,
        depth_mm: 5. * 3_f64.sqrt() / 16. * pitch,
        crest_width_mm: pitch / 8.,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        radial_taper: 0.,
        starts,
        left_hand,
    }
}

/// A plain ISO-like thread with a rounded root.
fn plain(left_hand: bool, starts: u32, frame: PlaneFrame) -> Solid {
    let pitch = 0.5 * f64::from(starts);
    let root = pitch / (4. * 3_f64.sqrt());
    Solid::rust_helical_thread(
        HelicalThread {
            depth_mm: (pitch - pitch / 8.) / (2. * 30_f64.to_radians().portable_tan()) - root,
            ..iso(6.6, left_hand, starts, pitch)
        },
        Some(root),
        EntityId::new(),
        frame,
        GeometryPolicy::default(),
    )
    .unwrap()
}

/// A radius-5, 4 mm stock cylinder with an internal thread of `length`,
/// entered at its base or, with `top`, at its far cap.
fn internal(length: f64, left_hand: bool, top: bool, frame: PlaneFrame) -> Solid {
    let policy = GeometryPolicy::default();
    let n = frame.normal();
    let entry = if top {
        PlaneFrame {
            origin: std::array::from_fn(|i| frame.origin[i] + 4. * n[i]),
            u: frame.u,
            v: frame.v.map(|x| -x),
        }
    } else {
        frame
    };
    Solid::cylinder(frame, 5., 4., policy)
        .unwrap()
        .rust_cut_internal_thread(
            iso(length, left_hand, 1, 1.),
            0.2,
            EntityId::new(),
            entry,
            policy,
        )
        .unwrap()
}

/// The inventory's 5.5 mm across-flats hex nut.
fn nut(left_hand: bool, frame: PlaneFrame) -> Solid {
    Solid::rust_hex_nut(
        HelicalThread {
            major_radius_mm: 1.55,
            pitch_mm: 0.5,
            length_mm: 2.4,
            depth_mm: 0.28,
            crest_width_mm: 0.08,
            ..iso(2.4, left_hand, 1, 1.)
        },
        5.5,
        0.2,
        EntityId::new(),
        frame,
        GeometryPolicy::default(),
    )
    .unwrap()
}

/// Every fixture with the view centres where its sections are certified.
fn fixtures(frame: PlaneFrame) -> Vec<(&'static str, Solid, Vec<f64>)> {
    vec![
        ("plain", plain(false, 1, frame), vec![3.3]),
        ("plain two-start left", plain(true, 2, frame), vec![2.1]),
        ("overrun", internal(5.2, false, false, frame), vec![2.]),
        ("exact left", internal(4., true, false, frame), vec![1.6]),
        ("floor", internal(2.5, false, false, frame), vec![1.2, 3.2]),
        (
            "floor top left",
            internal(2.5, true, true, frame),
            vec![0.9, 2.7],
        ),
        ("overrun top", internal(5.2, false, true, frame), vec![2.4]),
        ("nut", nut(false, frame), vec![1.2]),
        ("nut left", nut(true, frame), vec![0.9]),
    ]
}

fn contains(solid: &Solid, w: P3) -> bool {
    solid.rust_thread_contains(w).unwrap()
}

fn local(frame: PlaneFrame, w: P3) -> P3 {
    let offset: P3 = std::array::from_fn(|i| w[i] - frame.origin[i]);
    [frame.u, frame.v, frame.normal()].map(|axis| dot3(axis, offset))
}

fn hex_value(across_flats: f64, [x, y]: [f64; 2]) -> f64 {
    (0..6)
        .map(|k| {
            let (s, c) = (f64::from(k) * PI / 3.).portable_sin_cos();
            x * c + y * s
        })
        .fold(f64::NEG_INFINITY, f64::max)
        - 0.5 * across_flats
}

/// Whether a world point lies on the source boundary, read from the exact
/// descriptor: the helicoid, the caps, the stock wall, a partial-depth
/// void's floor, or the nut's flats.
fn on_boundary(solid: &Solid, w: P3) -> bool {
    let on = |value: f64| value.abs() < 1e-6;
    let within = |z: f64, a: f64, b: f64| z > a - 1e-9 && z < b + 1e-9;
    if let Some((thread, frame)) = solid.rust_helical_thread_data() {
        let [x, y, z] = local(frame, w);
        let (r, theta) = (x.portable_hypot(y), y.portable_atan2(x));
        let length = thread.spec().length_mm;
        let bound = thread.radius(theta, z.clamp(0., length));
        return (on(r - bound) && within(z, 0., length))
            || ((on(z) || on(z - length)) && r <= bound + 1e-6);
    }
    if let Some(data) = solid.rust_internal_thread_data() {
        let [x, y, z] = local(data.frame, w);
        let (r, theta) = (x.portable_hypot(y), y.portable_atan2(x));
        let length = data.thread.spec().length_mm;
        let (stock, top) = (data.stock_radius_mm, data.stock_length_mm);
        let floor = top > length;
        let bound = data.thread.radius(theta, z.clamp(0., length));
        return (on(r - bound) && within(z, 0., length))
            || (on(r - stock) && within(z, 0., top))
            || (on(z) && r >= bound - 1e-6 && r <= stock + 1e-6)
            || (on(z - top) && r <= stock + 1e-6 && (floor || r >= bound - 1e-6))
            || (floor && on(z - length) && r <= bound + 1e-6);
    }
    let (nut, frame) = solid.rust_hex_nut_geometry().unwrap();
    let [x, y, z] = local(frame, w);
    let (r, theta) = (x.portable_hypot(y), y.portable_atan2(x));
    let length = nut.thread().spec().length_mm;
    let bound = nut.thread().radius(theta, z.clamp(0., length));
    let flat = hex_value(nut.across_flats_mm(), [x, y]);
    (on(r - bound) && within(z, 0., length))
        || (on(flat) && within(z, 0., length))
        || ((on(z) || on(z - length)) && r >= bound - 1e-6 && flat <= 1e-6)
}

/// The exact area of a cut across the axis at local level `z`.
fn cross_section(solid: &Solid, z: f64) -> f64 {
    if let Some((thread, _)) = solid.rust_helical_thread_data() {
        return thread.volume_mm3() / thread.spec().length_mm;
    }
    if let Some(data) = solid.rust_internal_thread_data() {
        let length = data.thread.spec().length_mm;
        let disc = PI * data.stock_radius_mm * data.stock_radius_mm;
        return if z < length {
            disc - data.thread.volume_mm3() / length
        } else {
            disc
        };
    }
    let (nut, _) = solid.rust_hex_nut_geometry().unwrap();
    let thread = nut.thread();
    nut.outer_area_mm2() - thread.volume_mm3() / thread.spec().length_mm
}

fn frame_of(solid: &Solid) -> PlaneFrame {
    solid
        .rust_helical_thread_data()
        .map(|(_, frame)| frame)
        .or_else(|| solid.rust_internal_thread_data().map(|data| data.frame))
        .or_else(|| solid.rust_hex_nut_geometry().map(|(_, frame)| frame))
        .unwrap()
}

#[test]
fn fastener_material_along_lines_agrees_with_source_containment() {
    let frame = PlaneFrame {
        origin: [1.5, -2., 0.5],
        u: [0.8, 0.6, 0.],
        v: [-0.36, 0.48, 0.8],
    };
    let view = View::new(PlaneFrame::default());
    for (label, solid, _) in fixtures(frame) {
        let thread = frame_of(&solid);
        let axis = thread.normal();
        let world = |l: [f64; 3]| -> P3 {
            std::array::from_fn(|i| {
                thread.origin[i] + l[0] * thread.u[i] + l[1] * thread.v[i] + l[2] * axis[i]
            })
        };
        let mut lines: Vec<(P3, P3)> = Vec::new();
        for (k, d) in [
            [0.3, -0.5, 0.81],
            [1., 0.2, -0.1],
            [0.05, 0.99, 0.1],
            [-0.7, 0.7, 0.14],
        ]
        .into_iter()
        .enumerate()
        {
            let length: f64 = d.iter().map(|x: &f64| x * x).sum::<f64>().sqrt();
            let d = d.map(|v: f64| v / length);
            for offset in [0., 1.3, 2.4, 2.9, 4.2] {
                let start = world([
                    offset * (k as f64 * 0.7).portable_cos(),
                    offset * (k as f64 * 0.7).portable_sin(),
                    1.7,
                ]);
                lines.push((
                    start,
                    std::array::from_fn(|i| {
                        d[0] * thread.u[i] + d[1] * thread.v[i] + d[2] * axis[i]
                    }),
                ));
            }
        }
        // Along the axis at several radii, across a flat, near the entry
        // cap and on both sides of a partial-depth void's floor.
        for r in [0.5, 1.45, 1.6, 2.5, 2.8, 3.05, 3.3, 4.9] {
            lines.push((world([r, 0.3, -1.]), axis));
        }
        lines.push((world([2.7, -0.2, -1.]), axis));
        for z in [1e-3, 2.5 - 2e-3, 2.5 + 2e-3, 4. - 1e-3] {
            lines.push((world([-6., 0.41, z]), thread.u));
            lines.push((world([0.37, -6., z]), thread.v));
        }
        let member = families::member(&solid, &view).unwrap().unwrap();
        let mut covered = 0;
        for &(start, direction) in &lines {
            let spans = member.line(start, direction).unwrap().unwrap();
            let steps = 12_000;
            let mut mismatches = 0;
            for i in 0..=steps {
                let t = -12. + 24. * i as f64 / steps as f64;
                if spans.iter().flatten().any(|end| (end - t).abs() < 1e-4) {
                    continue;
                }
                let claimed = spans.iter().any(|[a, b]| *a < t && t < *b);
                let truth = contains(&solid, std::array::from_fn(|k| start[k] + t * direction[k]));
                mismatches += usize::from(claimed != truth);
            }
            assert_eq!(mismatches, 0, "{label} {start:?} {direction:?}: {spans:?}");
            covered += usize::from(!spans.is_empty());
        }
        assert!(covered >= lines.len() / 2, "{label}: {covered}");
    }
}

/// Every outline candidate lies on the helicoid `r = R(z - k theta)` with
/// the exact surface normal perpendicular to the sight.
#[test]
fn fastener_outlines_satisfy_the_contour_condition() {
    let s = 0.5_f64.sqrt();
    let t = (1_f64 / 6.).sqrt();
    for (label, solid, centers) in fixtures(PlaneFrame::default()) {
        let thread = solid
            .rust_helical_thread_data()
            .map(|(thread, _)| thread)
            .or_else(|| solid.rust_internal_thread_data().map(|data| data.thread))
            .or_else(|| solid.rust_hex_nut_geometry().map(|(nut, _)| *nut.thread()))
            .unwrap();
        let length = thread.spec().length_mm;
        let k = thread.helical_rate();
        for [u, v] in [
            [[1., 0., 0.], [0., 0., 1.]],
            [[s, -s, 0.], [t, t, -2. * t]],
            [[1., 0., 0.], [0., 0.9995_f64.sqrt(), (0.0005_f64).sqrt()]],
        ] {
            let frame = PlaneFrame {
                origin: [0.2, 0.1, centers[0]],
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
                for step in 0..=8 {
                    let q = candidate.path.point(step as f64 / 8.);
                    let w: P3 = std::array::from_fn(|i| {
                        frame.origin[i] + (0..3).map(|j| view.axes[j][i] * q[j]).sum::<f64>()
                    });
                    // Top-entry fixtures use a flipped thread frame.
                    let [x, y, z] = local(frame_of(&solid), w);
                    let (r, theta) = (x.portable_hypot(y), y.portable_atan2(x));
                    if let Some(data) = solid.rust_internal_thread_data()
                        && (r - data.stock_radius_mm).abs() < 1e-9
                    {
                        // The stock cylinder's silhouette generator.
                        let f = frame_of(&solid);
                        let radial: P3 = std::array::from_fn(|i| {
                            f.u[i] * theta.portable_cos() + f.v[i] * theta.portable_sin()
                        });
                        assert!(
                            z > -1e-9 && z < data.stock_length_mm + 1e-9,
                            "{label} {w:?}"
                        );
                        assert!(dot3(radial, toward).abs() < 1e-7, "{label} {w:?}");
                        checked += 1;
                        continue;
                    }
                    assert!(z > -1e-9 && z < length + 1e-9, "{label} {w:?}");
                    assert!((r - thread.radius(theta, z)).abs() < 1e-9, "{label} {w:?}");
                    let phase = z.rem_euclid(thread.spec().pitch_mm) - k * theta.rem_euclid(TAU);
                    let (c, sn) = (theta.portable_cos(), theta.portable_sin());
                    let f = frame_of(&solid);
                    let axes = [f.u, f.v, f.normal()];
                    // At a profile break (an outline ending on a crest or
                    // root edge) either one-sided normal applies.
                    let dot = [-1e-9, 0., 1e-9]
                        .map(|shift| {
                            let [_, slope] = thread.radius_and_slope_at_phase(phase + shift);
                            let local_normal =
                                [c - slope * k / r * sn, sn + slope * k / r * c, -slope];
                            let normal: P3 = std::array::from_fn(|i| {
                                (0..3).map(|j| axes[j][i] * local_normal[j]).sum()
                            });
                            dot3(normal, toward) / dot3(normal, normal).sqrt()
                        })
                        .into_iter()
                        .fold(f64::INFINITY, |m, d| if d.abs() < m.abs() { d } else { m });
                    assert!(dot.abs() < 1e-7, "{label} {w:?}: {dot}");
                    checked += 1;
                }
            }
            assert!(checked > 50, "{label}: {checked}");
        }
    }
}

#[test]
fn fastener_visibility_agrees_with_source_containment_marching() {
    let s = 0.5_f64.sqrt();
    let t = (1_f64 / 6.).sqrt();
    for (label, solid, centers) in fixtures(PlaneFrame::default()) {
        for [u, v] in [
            [[1., 0., 0.], [0., 0., 1.]],
            [[s, -s, 0.], [t, t, -2. * t]],
            [[1., 0., 0.], [0., 1., 0.]],
        ] {
            let frame = PlaneFrame {
                origin: [0.2, 0.1, centers[0]],
                u,
                v,
            };
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
            let world = |q: P3| -> P3 {
                std::array::from_fn(|i| {
                    frame.origin[i] + (0..3).map(|j| view.axes[j][i] * q[j]).sum::<f64>()
                })
            };
            let covered = |q: P3| {
                (1..6000).any(|i| {
                    let ahead = [q[0], q[1], q[2] + 1e-4 + 2e-3 * f64::from(i)];
                    contains(&solid, world(ahead))
                        && contains(&solid, world([q[0], q[1], ahead[2] - 5e-5]))
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
                        (point[0] - a0 - along * dx).portable_hypot(point[1] - a1 - along * dy)
                            < 0.02
                    })
                })
            };
            let member = families::member(&solid, &view).unwrap().unwrap();
            let (mut seen, mut wrong) = (0, Vec::new());
            for candidate in member.candidates(0.0025).unwrap() {
                if candidate.kind == Kind::Smooth {
                    continue;
                }
                for k in 1..12 {
                    let tt = k as f64 / 12.;
                    let q = candidate.path.point(tt);
                    let truth = covered(q);
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
                        wrong.push((candidate.kind, q, truth));
                    }
                }
            }
            assert!(seen > 100, "{label}: {seen}");
            assert!(
                wrong.is_empty(),
                "{label} {u:?}: {} of {seen}: {:?}",
                wrong.len(),
                &wrong[..wrong.len().min(6)]
            );
        }
    }
}

fn hatched() -> Option<DrawingHatch> {
    Some(DrawingHatch {
        angle_rad: 0.,
        spacing_mm: 1.,
        line_width_mm: 0.18,
        pattern: HatchPattern::Solid,
    })
}

fn hatch_area(vectors: &[ProjectedVector]) -> f64 {
    vectors
        .iter()
        .filter_map(|v| match v.curve {
            VectorCurve::HatchTriangle { points: [a, b, c] } => {
                Some(((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])).abs() / 2.)
            }
            _ => None,
        })
        .sum()
}

/// The area of the cut by containment on a grid covering the published
/// section with a margin.
fn grid_area(
    solid: &Solid,
    frame: PlaneFrame,
    view: &View,
    vectors: &[ProjectedVector],
    step: f64,
) -> f64 {
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
            let w: P3 = std::array::from_fn(|k| {
                frame.origin[k] + view.axes[0][k] * q[0] + view.axes[1][k] * q[1]
            });
            inside += usize::from(contains(solid, w));
        }
    }
    inside as f64 * step * step
}

/// Every published section boundary point lies in the cut plane on the
/// source boundary; the hatch has the grid area of the cut, and across the
/// axis the exact cross-section area.
#[test]
fn fastener_sections_lie_on_the_source_boundary() {
    let s = 0.5_f64.sqrt();
    let t = (1_f64 / 6.).sqrt();
    for (label, solid, centers) in fixtures(PlaneFrame::default()) {
        for center in centers {
            for (origin, [u, v]) in [
                ([0.2, 0.1, center], [[1., 0., 0.], [0., 0., 1.]]),
                ([0., 0., center], [[s, s, 0.], [0., 0., 1.]]),
                ([0.2, 0.1, center], [[s, -s, 0.], [t, t, -2. * t]]),
                ([0.3, -0.2, center + 0.37], [[1., 0., 0.], [0., 1., 0.]]),
            ] {
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
                        solid: &solid,
                        hatch: hatched(),
                    }],
                    options,
                )
                .unwrap()
                .unwrap_or_else(|error| panic!("{label} {origin:?} {u:?}: {error:?}"));
                let mut checked = 0;
                for vector in vectors.iter().filter(|v| v.role == VectorRole::Section) {
                    for q in &vector.preview_mm {
                        let w: P3 = std::array::from_fn(|i| {
                            frame.origin[i] + view.axes[0][i] * q[0] + view.axes[1][i] * q[1]
                        });
                        assert!(on_boundary(&solid, w), "{label} {origin:?} {u:?}: {w:?}");
                        checked += 1;
                    }
                }
                assert!(checked > 20, "{label}: {checked}");
                let area = hatch_area(&vectors);
                let grid = grid_area(&solid, frame, &view, &vectors, 0.01);
                assert!(
                    (area - grid).abs() < 0.01 * grid,
                    "{label} {origin:?} {u:?}: {area} {grid}"
                );
                if v == [0., 1., 0.] {
                    let z = local(frame_of(&solid), origin)[2];
                    let exact = cross_section(&solid, z);
                    assert!(
                        (area - exact).abs() < 2e-3 * exact,
                        "{label} {z}: {area} {exact}"
                    );
                }
            }
        }
    }
}

/// Side cuts just off the axis close with the grid area.
#[test]
fn fastener_near_axis_side_cuts_close_with_the_grid_area() {
    for (label, solid, centers) in fixtures(PlaneFrame::default()) {
        for offset in [1e-6_f64, 1e-4, 1e-3, 1e-2] {
            for angle in [0., 1., 2.6] {
                let (sin, cos) = f64::portable_sin_cos(angle);
                let frame = PlaneFrame {
                    origin: [-offset * sin, offset * cos, centers[0]],
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
                        hatch: hatched(),
                    }],
                    options,
                )
                .unwrap()
                .unwrap_or_else(|error| panic!("{label} {offset} {angle}: {error:?}"));
                let area = hatch_area(&vectors);
                let grid = grid_area(&solid, frame, &view, &vectors, 0.02);
                assert!(
                    (area - grid).abs() < 0.01 * grid,
                    "{label} {offset} {angle}: {area} {grid}"
                );
            }
        }
    }
}

/// Wherever the first-material depth found by marching the source solid's
/// containment jumps between neighbouring pixels, a published visible
/// vector passes within two and a half pixels.
#[test]
fn fastener_occluding_contours_agree_with_a_depth_map_oracle() {
    let s = 0.5_f64.sqrt();
    let t = (1_f64 / 6.).sqrt();
    for (label, solid, centers) in
        fixtures(PlaneFrame::default())
            .into_iter()
            .filter(|(label, ..)| {
                ["plain", "overrun", "floor", "floor top left", "nut"].contains(label)
            })
    {
        for [u, v] in [[[1., 0., 0.], [0., 0., 1.]], [[s, -s, 0.], [t, t, -2. * t]]] {
            let frame = PlaneFrame {
                origin: [0.1, -0.05, centers[0]],
                u,
                v,
            };
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
            let world = |q: P3| -> P3 {
                std::array::from_fn(|i| {
                    frame.origin[i] + (0..3).map(|j| view.axes[j][i] * q[j]).sum::<f64>()
                })
            };
            let inside = |q: P3| contains(&solid, world(q));
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
            let (pixel, count, edge) = (0.05, 240, -6.);
            let grid: Vec<Vec<Option<f64>>> = (0..count)
                .map(|i| {
                    (0..count)
                        .map(|j| {
                            depth(
                                edge + pixel * (i as f64 + 0.5),
                                edge + pixel * (j as f64 + 0.5),
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
                            (x - a0 - along * dx).portable_hypot(y - a1 - along * dy) < 2.5 * pixel
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
                        let x = edge + pixel * (0.5 * (i + k) as f64 + 0.5);
                        let y = edge + pixel * (0.5 * (j + l) as f64 + 0.5);
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

/// Every curved path of the three families flattens into chords whose proven
/// bound meets half the tolerance and holds the sampled deviation.
#[test]
fn fastener_chord_bounds_hold_and_meet_the_tolerance() {
    let s = 0.5_f64.sqrt();
    let t = (1_f64 / 6.).sqrt();
    let tolerance = 0.0025;
    for (label, solid, centers) in fixtures(PlaneFrame::default()) {
        let mut total = helix::Tally::default();
        for [u, v] in [
            [[1., 0., 0.], [0., 0., 1.]],
            [[s, s, 0.], [0., 0., 1.]],
            [[s, -s, 0.], [t, t, -2. * t]],
            [[1., 0., 0.], [0., 1., 0.]],
            [[1., 0., 0.], [0., 0.9995_f64.sqrt(), 0.0005_f64.sqrt()]],
        ] {
            for center in &centers {
                let view = View::new(PlaneFrame {
                    origin: [0.2, 0.1, *center],
                    u,
                    v,
                });
                let member = families::member(&solid, &view).unwrap().unwrap();
                let mut paths: Vec<_> = member
                    .candidates(tolerance)
                    .unwrap()
                    .into_iter()
                    .map(|candidate| candidate.path)
                    .collect();
                paths.extend(member.section(0., tolerance).unwrap().into_iter().flatten());
                for path in paths {
                    helix::tally(&*path, tolerance, &mut total, label);
                }
            }
        }
        assert!(total.chords > 100, "{label}: {total:?}");
        assert_eq!(total.unproven, 0, "{label}: {total:?}");
    }
}
