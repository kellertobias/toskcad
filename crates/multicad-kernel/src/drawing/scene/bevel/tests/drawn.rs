//! Drawing-level checks of the bevel views and sections that OpenCascade
//! 7.9.3 refuses for its own bevel ("Projected curve family is
//! unsupported"), from sources independent of the scene member: the
//! provider's membership predicate marched along sight lines, the provider's
//! own topology edges, and scan-line quadrature of the cut region.
use super::*;
use crate::{DrawingHatch, DrawingProjection, DrawingSolid, ProjectedVector, project_drawing};
use multicad_core::drawings::{VectorCurve, VectorRole};
use multicad_rust_kernel::bevel_gear::{Curve as TopologyCurve, Piece as TopologyPiece};

/// Strict membership of one world point in one source.
type Oracle = Box<dyn Fn([f64; 3]) -> bool>;

fn oracle(solid: &Solid) -> Oracle {
    if let Some((body, frame)) = solid.bevel_body() {
        let bevel = body.geometry.clone();
        let axes = [frame.u, frame.v, frame.normal()];
        let origin = frame.origin;
        return Box::new(move |p| {
            let offset = sub(p, origin);
            bevel.classify(axes.map(|axis| dot3(axis, offset)), 0.) == Containment::Inside
        });
    }
    let b = solid.bounds().unwrap();
    Box::new(move |p| (0..3).all(|k| b.min[k] < p[k] && p[k] < b.max[k]))
}

/// Sharp edges of one source in world coordinates, as dense samples.
fn edges(solid: &Solid) -> Vec<Vec<[f64; 3]>> {
    let samples = 10;
    if let Some((body, frame)) = solid.bevel_body() {
        let bevel = &body.geometry;
        let topology = bevel.topology();
        let pieces = bevel.pieces();
        let world = |p: [f64; 3]| -> [f64; 3] {
            let n = frame.normal();
            std::array::from_fn(|i| {
                frame.origin[i] + p[0] * frame.u[i] + p[1] * frame.v[i] + p[2] * n[i]
            })
        };
        return topology
            .edges
            .iter()
            .filter(|edge| match edge.curve {
                TopologyCurve::Profile { .. } | TopologyCurve::BoreArc { .. } => true,
                TopologyCurve::Ray { joint, .. } => {
                    joint > 0
                        && (pieces[joint - 1] == TopologyPiece::Tip
                            || pieces[joint] == TopologyPiece::Tip)
                }
                _ => false,
            })
            .map(|edge| {
                (1..samples)
                    .map(|k| {
                        let s = f64::from(k) / f64::from(samples);
                        world(bevel.curve_point(&topology, edge, s))
                    })
                    .collect()
            })
            .collect();
    }
    let b = solid.bounds().unwrap();
    let corner = |i: usize| -> [f64; 3] {
        std::array::from_fn(|k| if i >> k & 1 == 1 { b.max[k] } else { b.min[k] })
    };
    let mut output = Vec::new();
    for a in 0..8 {
        for bit in 0..3 {
            let c = a | 1 << bit;
            if c != a {
                let [p, q] = [corner(a), corner(c)];
                output.push(
                    (1..4 * samples)
                        .map(|k| {
                            let f = f64::from(k) / f64::from(4 * samples);
                            std::array::from_fn(|i| p[i] + (q[i] - p[i]) * f)
                        })
                        .collect(),
                );
            }
        }
    }
    output
}

fn distance_to(point: P2, vectors: &[&ProjectedVector]) -> f64 {
    let mut best = f64::INFINITY;
    for vector in vectors {
        for pair in vector.preview_mm.windows(2) {
            let [a, b] = [pair[0], pair[1]];
            let [u, v] = [b[0] - a[0], b[1] - a[1]];
            let along = (((point[0] - a[0]) * u + (point[1] - a[1]) * v)
                / (u * u + v * v).max(1e-300))
            .clamp(0., 1.);
            best =
                best.min((point[0] - a[0] - along * u).portable_hypot(point[1] - a[1] - along * v));
        }
    }
    best
}

/// Distance to the nearest chord of `vectors` and the sine of the angle
/// between that chord and `direction`.
fn nearest(point: P2, direction: P2, vectors: &[&ProjectedVector]) -> (f64, f64) {
    let mut best = (f64::INFINITY, 1.);
    for vector in vectors {
        for pair in vector.preview_mm.windows(2) {
            let [a, b] = [pair[0], pair[1]];
            let [u, v] = [b[0] - a[0], b[1] - a[1]];
            let length = u.portable_hypot(v);
            if length == 0. {
                continue;
            }
            let along =
                (((point[0] - a[0]) * u + (point[1] - a[1]) * v) / (length * length)).clamp(0., 1.);
            let gap = (point[0] - a[0] - along * u).portable_hypot(point[1] - a[1] - along * v);
            if gap < best.0 {
                best = (gap, ((u * direction[1] - v * direction[0]) / length).abs());
            }
        }
    }
    best
}

#[derive(Default, Debug)]
struct Report {
    edge_samples: usize,
    ambiguous: usize,
    visible_missing: usize,
    hidden_missing: usize,
    hidden_drawn_visible: usize,
    /// Hidden samples under a visible curve explained neither by a visible
    /// edge through the same spot nor by the silhouette of the material in
    /// front passing within the drawing tolerance.
    hidden_unexplained: usize,
    section_points: usize,
    section_off_boundary: usize,
    boundary_points: usize,
    boundary_missing: usize,
    worst_boundary_gap: f64,
    area: [f64; 2],
}

/// Check one drawing of `solids` against the independent sources.
fn check(label: &str, solids: &[&Solid], options: DrawingProjection) -> Report {
    let hatch = DrawingHatch {
        pattern: multicad_core::materials::HatchPattern::Solid,
        angle_rad: 0.,
        spacing_mm: 1.,
        line_width_mm: 0.1,
    };
    let inputs: Vec<_> = solids
        .iter()
        .map(|&solid| DrawingSolid {
            solid,
            hatch: Some(hatch),
        })
        .collect();
    let vectors = project_drawing(&inputs, options).unwrap_or_else(|e| panic!("{label}: {e:?}"));
    let view = View::new(options.frame);
    let toward = view.axes[2];
    let world = |p: P3| -> [f64; 3] {
        std::array::from_fn(|i| {
            view.origin[i] + (0..3).map(|k| p[k] * view.axes[k][i]).sum::<f64>()
        })
    };
    let top = options.section_offset_mm.unwrap();
    let far = options.section_depth_mm.map(|d| top - d);
    let oracles: Vec<Oracle> = solids.iter().map(|s| oracle(s)).collect();
    let in_window = |p: P2, margin: f64| {
        [options.section_u_mm, options.section_v_mm]
            .iter()
            .enumerate()
            .all(|(k, w)| w.is_none_or(|[a, b]| a + margin <= p[k] && p[k] <= b - margin))
            && options
                .detail
                .is_none_or(|(c, r)| (p[0] - c[0]).portable_hypot(p[1] - c[1]) <= r - margin)
    };
    let tolerance = 3. * options.chord_tolerance_mm;
    let mut report = Report::default();
    let of = |source: usize, roles: &[VectorRole]| -> Vec<&ProjectedVector> {
        vectors
            .iter()
            .filter(|v| v.source_index == source && roles.contains(&v.role))
            .collect()
    };
    // Every solid's view-space depth extent, for marching.
    let reach = solids
        .iter()
        .map(|s| {
            let b = s.bounds().unwrap();
            (0..8)
                .map(|i| {
                    let corner: [f64; 3] =
                        std::array::from_fn(|k| if i >> k & 1 == 1 { b.max[k] } else { b.min[k] });
                    view.point(corner)[2]
                })
                .fold(f64::NEG_INFINITY, f64::max)
        })
        .fold(f64::NEG_INFINITY, f64::max);

    // Samples of material in front of a view-space point within the
    // retained depth, marched in 0.01 mm steps and capped at three.
    let front = |p: P3| {
        let step = 0.01;
        let mut count = 0;
        let mut z = p[2] + step;
        while z < top.min(reach) && count < 3 {
            let q = world([p[0], p[1], z]);
            count += oracles.iter().filter(|inside| inside(q)).count().min(1);
            z += step;
        }
        count
    };
    // Edges: visibility by marching the membership of every source in front
    // of each sample within the retained depth.
    if !options.section_only {
        for (source, solid) in solids.iter().enumerate() {
            let visible = of(source, &[VectorRole::Visible]);
            let drawn = of(source, &[VectorRole::Visible, VectorRole::Hidden]);
            let mut seen: Vec<[P2; 2]> = Vec::new();
            let mut covered: Vec<P3> = Vec::new();
            for edge in edges(solid) {
                let mut previous: Option<P2> = None;
                let projected: Vec<P3> = edge.iter().map(|&point| view.point(point)).collect();
                for (index, &p) in projected.iter().enumerate() {
                    // Paper direction of the edge at the sample.
                    let [a, b] = [
                        index.saturating_sub(1),
                        (index + 1).min(projected.len() - 1),
                    ]
                    .map(|k| projected[k]);
                    let span = (b[0] - a[0]).portable_hypot(b[1] - a[1]);
                    let direction = if span > 0. {
                        [(b[0] - a[0]) / span, (b[1] - a[1]) / span]
                    } else {
                        [0., 0.]
                    };
                    let eta = 0.02;
                    if p[2] > top - eta
                        || far.is_some_and(|far| p[2] < far + eta)
                        || !in_window([p[0], p[1]], eta)
                    {
                        continue;
                    }
                    report.edge_samples += 1;
                    let count = front(p);
                    let paper = [p[0], p[1]];
                    let last = previous.take();
                    match count {
                        0 => {
                            if let Some(last) = last {
                                seen.push([last, paper]);
                            }
                            previous = Some(paper);
                            if distance_to(paper, &visible) > tolerance {
                                report.visible_missing += 1;
                            }
                        }
                        3 => {
                            if options.hidden_lines && distance_to(paper, &drawn) > tolerance {
                                report.hidden_missing += 1;
                            }
                            // Under a visible curve along it, not merely
                            // crossing it.
                            let (gap, sine) = nearest(paper, direction, &visible);
                            if gap <= tolerance
                                && span > 0.
                                && sine <= 10_f64.to_radians().portable_sin()
                            {
                                report.hidden_drawn_visible += 1;
                                covered.push(p);
                            }
                        }
                        _ => report.ambiguous += 1,
                    }
                }
            }
            // A rear edge legitimately lies under a visible curve where a
            // front edge projects onto it, such as a bore circle seen
            // edge-on. Visible edge samples are 1/10 of an edge apart, so
            // look along the curve within a sample spacing.
            for point in covered {
                let paper = [point[0], point[1]];
                let on_visible_edge = seen.iter().any(|&[a, b]| {
                    let [u, v] = [b[0] - a[0], b[1] - a[1]];
                    let along = (((paper[0] - a[0]) * u + (paper[1] - a[1]) * v)
                        / (u * u + v * v).max(1e-300))
                    .clamp(0., 1.);
                    (paper[0] - a[0] - along * u).portable_hypot(paper[1] - a[1] - along * v)
                        <= tolerance
                });
                let silhouette = (0..16).any(|k| {
                    let (sin, cos) = (f64::from(k) * TAU / 16.).portable_sin_cos();
                    front([
                        point[0] + tolerance * cos,
                        point[1] + tolerance * sin,
                        point[2],
                    ]) == 0
                });
                let explained = on_visible_edge || silhouette;
                if !explained {
                    report.hidden_unexplained += 1;
                    println!("SUSPECT {label} source {source} at {point:?}");
                }
            }
        }
    }

    // Sections: the cut region of each source in the plane, with material
    // immediately behind it, inside the window and detail.
    for (source, inside) in oracles.iter().enumerate() {
        let region = |p: P2| in_window(p, 0.) && inside(world([p[0], p[1], top - 1e-7]));
        let sections = of(source, &[VectorRole::Section]);
        for vector in &sections {
            for &q in &vector.preview_mm {
                report.section_points += 1;
                // On the boundary: material and void both within 0.004 mm.
                // A cut vertex can open a narrow wedge, so probe densely in
                // angle at several radii.
                let mut states = [false; 2];
                for radius in [4e-3, 1e-3, 1e-4, 1e-5] {
                    for k in 0..64 {
                        let (sin, cos) = (f64::from(k) * TAU / 64.).portable_sin_cos();
                        states[usize::from(region([q[0] + radius * cos, q[1] + radius * sin]))] =
                            true;
                    }
                }
                if !(states[0] && states[1]) {
                    report.section_off_boundary += 1;
                    println!("OFF {label} source {source} at {q:?}");
                }
            }
        }
        let triangles: f64 = vectors
            .iter()
            .filter(|v| v.source_index == source)
            .filter_map(|v| match v.curve {
                VectorCurve::HatchTriangle { points: [a, b, c] } => Some(
                    0.5 * ((b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1])).abs(),
                ),
                _ => None,
            })
            .sum();
        // Scan lines across the source's paper bounds.
        let b = solids[source].bounds().unwrap();
        let mut lo = [f64::INFINITY; 2];
        let mut hi = [f64::NEG_INFINITY; 2];
        for i in 0..8 {
            let corner: [f64; 3] =
                std::array::from_fn(|k| if i >> k & 1 == 1 { b.max[k] } else { b.min[k] });
            let p = view.point(corner);
            for k in 0..2 {
                lo[k] = lo[k].min(p[k]);
                hi[k] = hi[k].max(p[k]);
            }
        }
        // Rows tile the window's v-range exactly: a horizontal border then
        // never cuts through a row of the midpoint rule.
        if let Some([a, b]) = options.section_v_mm {
            lo[1] = lo[1].max(a);
            hi[1] = hi[1].min(b);
        }
        let (rows, columns) = (300, 600);
        let dy = (hi[1] - lo[1]) / f64::from(rows);
        let dx = (hi[0] - lo[0] + 2e-3) / f64::from(columns);
        let mut area = 0.;
        let drawn_section = !sections.is_empty();
        for row in 0..rows {
            let y = lo[1] + dy * (f64::from(row) + 0.5);
            let mut previous = (lo[0] - 1e-3, false);
            for column in 1..=columns {
                let x = lo[0] - 1e-3 + dx * f64::from(column);
                let state = region([x, y]);
                if state != previous.1 {
                    let (mut a, mut c) = (previous.0, x);
                    for _ in 0..40 {
                        let mid = 0.5 * (a + c);
                        if region([mid, y]) == previous.1 {
                            a = mid;
                        } else {
                            c = mid;
                        }
                    }
                    area += if state { -a } else { a } * dy;
                    if drawn_section {
                        report.boundary_points += 1;
                        let gap = distance_to([a, y], &sections);
                        report.worst_boundary_gap = report.worst_boundary_gap.max(gap);
                        if gap > tolerance {
                            report.boundary_missing += 1;
                        }
                    }
                }
                previous = (x, state);
            }
        }
        if drawn_section {
            report.area[0] += triangles;
            report.area[1] += area;
        }
    }
    let _ = toward;
    println!("{{\"case\":\"{label}\",\"report\":\"{report:?}\"}}");
    report
}

const SCENE: RigidTransform = RigidTransform {
    rotation: [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
    translation_mm: [2., 0., 20.],
};

/// The cases of `rust_scene_bevel.rs` that OpenCascade refuses.
fn refused(transform: RigidTransform, origin: [f64; 3]) -> Vec<(String, DrawingProjection)> {
    let column = |k: usize| -> [f64; 3] { std::array::from_fn(|i| transform.rotation[i][k]) };
    let [x, y, z] = [column(0), column(1), column(2)];
    let mix = |a: f64, b: f64, c: f64| -> [f64; 3] {
        std::array::from_fn(|i| a * x[i] + b * y[i] + c * z[i])
    };
    let s = 0.5_f64.sqrt();
    let t = (1_f64 / 6.).sqrt();
    let side = PlaneFrame { origin, u: x, v: z };
    let oblique = PlaneFrame {
        origin,
        u: mix(s, -s, 0.),
        v: mix(t, t, -2. * t),
    };
    let section = |frame| DrawingProjection {
        frame,
        section_offset_mm: Some(0.),
        ..Default::default()
    };
    let mut cases = vec![(
        "side section hidden".to_string(),
        DrawingProjection {
            hidden_lines: true,
            section_offset_mm: Some(0.3),
            ..section(side)
        },
    )];
    let base = section(oblique);
    cases.extend([
        ("oblique section".to_string(), base),
        (
            "oblique section hidden".to_string(),
            DrawingProjection {
                hidden_lines: true,
                section_offset_mm: Some(0.3),
                ..base
            },
        ),
        (
            "oblique section-only".to_string(),
            DrawingProjection {
                section_only: true,
                ..base
            },
        ),
        (
            "oblique finite depth".to_string(),
            DrawingProjection {
                section_depth_mm: Some(2.),
                ..base
            },
        ),
        (
            "oblique window".to_string(),
            DrawingProjection {
                section_u_mm: Some([-14., 12.5]),
                section_v_mm: Some([-3.5, 9.]),
                section_depth_mm: Some(4.5),
                ..base
            },
        ),
        (
            "oblique section detail".to_string(),
            DrawingProjection {
                detail: Some(([9., 1.5], 7.)),
                ..base
            },
        ),
    ]);
    cases
}

fn center(solid: &Solid) -> [f64; 3] {
    let b = solid.bounds().unwrap();
    std::array::from_fn(|i| (b.min[i] + b.max[i]) / 2.)
}

#[test]
fn views_and_sections_refused_by_open_cascade_match_marched_membership() {
    let make = |spec, transform: RigidTransform| {
        Solid::rust_spherical_bevel_gear_named(spec, EntityId::new(), GeometryPolicy::default())
            .unwrap()
            .transformed(transform)
            .unwrap()
    };
    let bored = make(specification(2., 20, 30, 8., 3.), RigidTransform::default());
    let placed = make(specification(2., 16, 16, 6., 0.), PLACED);
    let stock = Solid::rust_box(
        EntityId::new().uuid_bits(),
        [20, 12, 8].map(|v| multicad_rust_kernel::Rational::new(v, 1).unwrap()),
        PlaneFrame::default(),
        GeometryPolicy::default(),
    )
    .unwrap();
    let scene = make(specification(2., 20, 30, 8., 3.), SCENE);
    let mut totals = Report::default();
    let mut cases = 0;
    for (label, solids, transform, origin) in [
        (
            "bored",
            vec![&bored],
            RigidTransform::default(),
            center(&bored),
        ),
        ("solid placed", vec![&placed], PLACED, center(&placed)),
        (
            "box+bevel",
            vec![&stock, &scene],
            RigidTransform::default(),
            [2., 0., 45.],
        ),
    ] {
        for (name, options) in refused(transform, origin) {
            let report = check(&format!("{label} {name}"), &solids, options);
            cases += 1;
            totals.edge_samples += report.edge_samples;
            totals.ambiguous += report.ambiguous;
            totals.visible_missing += report.visible_missing;
            totals.hidden_missing += report.hidden_missing;
            totals.hidden_drawn_visible += report.hidden_drawn_visible;
            totals.hidden_unexplained += report.hidden_unexplained;
            totals.section_points += report.section_points;
            totals.section_off_boundary += report.section_off_boundary;
            totals.boundary_points += report.boundary_points;
            totals.boundary_missing += report.boundary_missing;
            totals.worst_boundary_gap = totals.worst_boundary_gap.max(report.worst_boundary_gap);
            let [drawn, scanned] = report.area;
            assert!(
                (drawn - scanned).abs() <= 2e-3 * scanned + 2e-2,
                "{label} {name}: hatch {drawn} against scanned {scanned}"
            );
            assert!(
                report.section_points > 0 && report.boundary_points > 0,
                "{label} {name}"
            );
        }
    }
    println!("{{\"cases\":{cases},\"totals\":\"{totals:?}\"}}");
    assert_eq!(cases, 21);
    assert_eq!(totals.visible_missing, 0);
    assert_eq!(totals.hidden_missing, 0);
    assert_eq!(totals.section_off_boundary, 0);
    assert_eq!(totals.hidden_unexplained, 0);
    assert_eq!(totals.boundary_missing, 0);
}
