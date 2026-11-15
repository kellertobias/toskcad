//! Helical and continuous herringbone gears drawn through the exact-curve
//! scene engine (CAD-346, findings F10, F11 and F17): non-axial, cropped and
//! sectioned views and overlapping scenes. Expectations are derived here from
//! the gear specification, independently of the drawing code.
use multicad_core::{document::EntityId, drawings::VectorRole};
use multicad_kernel::{
    DrawingProjection, DrawingSolid, GearHub, GearSweep, GeometryPolicy, InvoluteGear, PlaneFrame,
    ProjectedVector, Solid, project_drawing,
    rust_kernel::{
        Rational,
        spur_gear::{Edge, Profile, Spec},
    },
};
use std::f64::consts::PI;

const WIDTH: f64 = 4.;
const BORE: f64 = 3.;
const PITCH: f64 = 12.;
const TEETH: u32 = 24;

fn spec(sweep: GearSweep, twist: f64, hub: Option<GearHub>) -> InvoluteGear {
    let a = 20_f64.to_radians();
    InvoluteGear {
        hub,
        teeth: TEETH,
        pitch_radius_mm: PITCH,
        base_radius_mm: PITCH * a.cos(),
        tip_radius_mm: PITCH + 1.,
        root_radius_mm: PITCH - 1.25,
        root_fillet_radius_mm: 0.2,
        half_base_angle_rad: (PI / 2. - 0.05) / (2. * PITCH) + a.tan() - a,
        face_width_mm: WIDTH,
        bore_radius_mm: BORE,
        twist_rad_per_mm: twist,
        sweep,
    }
}

fn gear(spec: InvoluteGear) -> Solid {
    let build = match spec.sweep {
        GearSweep::Helical => Solid::rust_helical_gear_named,
        GearSweep::GappedChevron { .. } => Solid::rust_double_helical_named,
        _ => Solid::rust_herringbone_named,
    };
    build(spec, EntityId::new(), GeometryPolicy::default()).unwrap()
}

fn hub() -> Option<GearHub> {
    Some(GearHub {
        radius_mm: 6.,
        axial_start_mm: -2.,
        axial_end_mm: 5.5,
    })
}

fn gears() -> Vec<(&'static str, InvoluteGear)> {
    vec![
        ("helical", spec(GearSweep::Helical, 0.05, None)),
        ("helical left hub", spec(GearSweep::Helical, -0.04, hub())),
        ("herringbone", spec(GearSweep::Chevron, 0.05, None)),
        (
            "herringbone left hub",
            spec(GearSweep::Chevron, -0.04, hub()),
        ),
        ("double-helical", spec(GAPPED, 0.05, None)),
        ("double-helical left hub", spec(GAPPED, -0.04, hub())),
    ]
}

const GAPPED: GearSweep = GearSweep::GappedChevron { gap_mm: 1. };

/// Axial length of each tooth band of a chevron.
fn band(spec: &InvoluteGear) -> f64 {
    match spec.sweep {
        GearSweep::GappedChevron { gap_mm } => 0.5 * (WIDTH - gap_mm),
        _ => 0.5 * WIDTH,
    }
}

fn frames(origin: [f64; 3]) -> [(&'static str, PlaneFrame); 3] {
    let s = 0.5_f64.sqrt();
    let t = (1_f64 / 6.).sqrt();
    [
        ("axial", [[1., 0., 0.], [0., 1., 0.]]),
        ("side", [[1., 0., 0.], [0., 0., 1.]]),
        ("oblique", [[s, -s, 0.], [t, t, -2. * t]]),
    ]
    .map(|(name, [u, v])| (name, PlaneFrame { origin, u, v }))
}

/// Plain and hidden-line views, then sections through `origin`: whole,
/// section-only, finite depth, lateral window and detail.
fn cases(origin: [f64; 3]) -> Vec<(String, DrawingProjection)> {
    let mut cases = Vec::new();
    for (name, frame) in frames(origin) {
        let view = DrawingProjection {
            frame,
            ..Default::default()
        };
        let section = DrawingProjection {
            section_offset_mm: Some(0.),
            ..view
        };
        cases.extend([
            (name.to_string(), view),
            (
                format!("{name} hidden"),
                DrawingProjection {
                    hidden_lines: true,
                    ..view
                },
            ),
            (
                format!("{name} detail"),
                DrawingProjection {
                    detail: Some(([0.5, 0.25], 3.)),
                    hidden_lines: true,
                    ..view
                },
            ),
            (format!("{name} section"), section),
            (
                format!("{name} section hidden"),
                DrawingProjection {
                    hidden_lines: true,
                    section_offset_mm: Some(0.3),
                    ..section
                },
            ),
            (
                format!("{name} section-only"),
                DrawingProjection {
                    section_only: true,
                    ..section
                },
            ),
            (
                format!("{name} finite depth"),
                DrawingProjection {
                    section_depth_mm: Some(1.),
                    ..section
                },
            ),
            (
                format!("{name} window"),
                DrawingProjection {
                    section_u_mm: Some([-4., 4.5]),
                    section_v_mm: Some([-1.5, 3.]),
                    section_depth_mm: Some(2.5),
                    ..section
                },
            ),
            (
                format!("{name} section detail"),
                DrawingProjection {
                    detail: Some(([1., 0.5], 4.)),
                    ..section
                },
            ),
        ]);
    }
    cases
}

fn draw(solids: &[&Solid], options: DrawingProjection) -> Vec<ProjectedVector> {
    let inputs: Vec<_> = solids
        .iter()
        .map(|&solid| DrawingSolid { solid, hatch: None })
        .collect();
    project_drawing(&inputs, options).unwrap_or_else(|error| panic!("{error:?}"))
}

fn center(solid: &Solid) -> [f64; 3] {
    let b = solid.bounds().unwrap();
    std::array::from_fn(|i| (b.min[i] + b.max[i]) / 2.)
}

fn count(vectors: &[ProjectedVector], role: VectorRole) -> usize {
    vectors.iter().filter(|vector| vector.role == role).count()
}

/// Every class of view that refused before now publishes curves of the
/// expected roles, and nothing leaves the gear's projected extent.
#[test]
fn every_view_class_of_a_twisted_gear_draws() {
    let only = std::env::var("SCENE_CASE").ok();
    for (name, spec) in gears() {
        let solid = gear(spec);
        let origin = center(&solid);
        for (case, options) in cases(origin) {
            let label = format!("{name} {case}");
            if only
                .as_ref()
                .is_some_and(|only| !label.contains(only.as_str()))
            {
                continue;
            }
            let started = std::time::Instant::now();
            let vectors = draw(&[&solid], options);
            println!(
                "{{\"case\":\"{label}\",\"visible\":{},\"hidden\":{},\"section\":{},\"seconds\":{:.2}}}",
                count(&vectors, VectorRole::Visible),
                count(&vectors, VectorRole::Hidden),
                count(&vectors, VectorRole::Section),
                started.elapsed().as_secs_f64()
            );
            assert!(!vectors.is_empty(), "{label}");
            let reach = (PITCH + 1.).hypot(8.) + 1e-6;
            for vector in &vectors {
                assert!(vector.preview_mm.len() >= 2, "{label}");
                for p in &vector.preview_mm {
                    assert!(p[0].hypot(p[1]) <= reach, "{label}: {p:?}");
                }
            }
            if options.section_offset_mm.is_some() {
                assert!(count(&vectors, VectorRole::Section) > 0, "{label}");
            }
            if options.section_only {
                assert_eq!(
                    count(&vectors, VectorRole::Section),
                    vectors.len(),
                    "{label}"
                );
            }
            if !options.hidden_lines {
                assert_eq!(count(&vectors, VectorRole::Hidden), 0, "{label}");
            }
        }
    }
}

/// Exact membership in an unplaced gear, independent of the drawing code:
/// the winding number of the exact profile curves about the point turned
/// back by the phase, with every chord refined until it provably subtends
/// the same angle as its curve.
struct Oracle {
    edges: Vec<Edge>,
    spec: InvoluteGear,
}

fn swept(edge: Edge, w: [f64; 2], [a, b]: [f64; 2], [pa, pb]: [[f64; 2]; 2], depth: u32) -> f64 {
    let chord = (pb[0] - pa[0]).hypot(pb[1] - pa[1]);
    let mid = [0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1])];
    if (w[0] - mid[0]).hypot(w[1] - mid[1]) > chord || depth > 60 {
        let [p, q] = [[pa[0] - w[0], pa[1] - w[1]], [pb[0] - w[0], pb[1] - w[1]]];
        return (p[0] * q[1] - p[1] * q[0]).atan2(p[0] * q[0] + p[1] * q[1]);
    }
    let m = 0.5 * (a + b);
    let pm = edge.point(m);
    swept(edge, w, [a, m], [pa, pm], depth + 1) + swept(edge, w, [m, b], [pm, pb], depth + 1)
}

impl Oracle {
    fn new(spec: InvoluteGear) -> Self {
        let profile = Profile::new(
            Spec {
                teeth: spec.teeth,
                pitch_radius_mm: spec.pitch_radius_mm,
                base_radius_mm: spec.base_radius_mm,
                tip_radius_mm: spec.tip_radius_mm,
                root_radius_mm: spec.root_radius_mm,
                root_fillet_radius_mm: spec.root_fillet_radius_mm,
                half_base_angle_rad: spec.half_base_angle_rad,
                face_width_mm: spec.face_width_mm,
                bore_radius_mm: spec.bore_radius_mm,
            },
            1e-6,
            1e6,
        )
        .unwrap();
        Self {
            edges: profile.outline().iter().map(|edge| edge.edge).collect(),
            spec,
        }
    }

    fn inside(&self, [x, y, z]: [f64; 3]) -> bool {
        let r = x.hypot(y);
        if r <= self.spec.bore_radius_mm {
            return false;
        }
        if !(z > 0. && z < WIDTH) {
            return self.spec.hub.is_some_and(|hub| {
                r < hub.radius_mm && z > hub.axial_start_mm && z < hub.axial_end_mm
            });
        }
        if r < self.spec.root_radius_mm - 1e-3 {
            return true;
        }
        // The central gap of a double-helical gear keeps the root core.
        if z > band(&self.spec) && z < WIDTH - band(&self.spec) {
            return r < self.spec.root_radius_mm;
        }
        if r > self.spec.tip_radius_mm + 1e-3 {
            return false;
        }
        let (sin, cos) = (-phase(&self.spec, z)).sin_cos();
        let w = [x * cos - y * sin, x * sin + y * cos];
        let mut total = 0.;
        for &edge in &self.edges {
            let mut previous = edge.point(0.);
            for i in 1..=4 {
                let [a, b] = [f64::from(i - 1) / 4., f64::from(i) / 4.];
                let next = edge.point(b);
                total += swept(edge, w, [a, b], [previous, next], 0);
                previous = next;
            }
        }
        total.abs() > PI
    }
}

/// The world point of a paper point at a depth toward the viewer.
fn world(frame: PlaneFrame, [x, y]: [f64; 2], depth: f64) -> [f64; 3] {
    let n = frame.normal();
    std::array::from_fn(|i| frame.origin[i] + x * frame.u[i] + y * frame.v[i] + depth * n[i])
}

/// Distance from a point to the nearest preview chord of the vectors.
fn distance_to(vectors: &[&ProjectedVector], [x, y]: [f64; 2]) -> f64 {
    vectors
        .iter()
        .flat_map(|vector| vector.preview_mm.windows(2))
        .map(|w| {
            let [[a, b], [c, d]] = [w[0], w[1]];
            let [u, v] = [c - a, d - b];
            let along = (((x - a) * u + (y - b) * v) / (u * u + v * v).max(1e-300)).clamp(0., 1.);
            (x - a - along * u).hypot(y - b - along * v)
        })
        .fold(f64::INFINITY, f64::min)
}

/// Looks along +Y: paper u is world X, paper v is world Z.
fn side(origin: [f64; 3]) -> PlaneFrame {
    PlaneFrame {
        origin,
        u: [1., 0., 0.],
        v: [0., 0., 1.],
    }
}

fn phase(spec: &InvoluteGear, z: f64) -> f64 {
    spec.twist_rad_per_mm
        * if spec.sweep == GearSweep::Helical {
            z
        } else {
            z.min(WIDTH - z).min(band(spec))
        }
}

/// A section plane through the axis cuts each involute flank where the
/// phase law turns the flank into the plane: at radius `r` on the side with
/// world azimuth `side`, tooth `k` meets the plane where
/// `phase(z) = side - k pitch -+ (half - inv(r))`. Every published vertex of
/// the section boundary satisfies this law or lies on an end face, the bore,
/// the tip or the root region, and every point of the law is drawn.
#[test]
fn a_section_through_the_axis_follows_the_exact_flank_law() {
    for (name, spec) in [
        ("helical", spec(GearSweep::Helical, 0.05, None)),
        ("herringbone", spec(GearSweep::Chevron, -0.05, None)),
        ("double-helical", spec(GAPPED, 0.05, None)),
    ] {
        let solid = gear(spec);
        let vectors = draw(
            &[&solid],
            DrawingProjection {
                frame: side([0.; 3]),
                section_offset_mm: Some(0.),
                section_only: true,
                ..Default::default()
            },
        );
        assert!(vectors.iter().all(|v| v.role == VectorRole::Section));
        let base = spec.base_radius_mm;
        let pitch = 2. * PI / f64::from(TEETH);
        let involute = |r: f64| {
            let a = (base / r).acos();
            a.tan() - a
        };
        let wrap = |a: f64| (a + PI).rem_euclid(2. * PI) - PI;
        let mut on_flank = 0;
        let mut worst: f64 = 0.;
        for vector in &vectors {
            for &[x, z] in &vector.preview_mm {
                let r = x.abs();
                assert!((-1e-9..=WIDTH + 1e-9).contains(&z), "{name}: {x} {z}");
                let end_face = z.abs() <= 1e-9 || (z - WIDTH).abs() <= 1e-9;
                let coaxial = [BORE, spec.tip_radius_mm]
                    .iter()
                    .any(|radius| (r - radius).abs() <= 1e-9);
                if end_face || coaxial {
                    continue;
                }
                assert!(
                    r > spec.root_radius_mm - 1e-9 && r < spec.tip_radius_mm + 1e-9,
                    "{name}: section vertex [{x}, {z}] is inside the body"
                );
                if r < base + 1e-3 {
                    // Root arc, fillet and radial flank below the base circle.
                    continue;
                }
                let azimuth = if x > 0. { 0. } else { PI };
                let local = azimuth - phase(&spec, z);
                let tooth = (local / pitch).round();
                let offset = wrap(local - tooth * pitch).abs();
                let residual = (offset - (spec.half_base_angle_rad - involute(r))).abs() * r;
                worst = worst.max(residual);
                on_flank += 1;
                assert!(
                    residual <= 1e-9,
                    "{name}: section vertex [{x}, {z}] is {residual} mm off the flank law"
                );
            }
        }
        assert!(on_flank > 50, "{name}: only {on_flank} flank vertices");
        // Completeness: points of the law are on the drawing.
        let all: Vec<&ProjectedVector> = vectors.iter().collect();
        let mut expected = 0;
        for (azimuth, sign) in [(0., 1.), (PI, -1.)] {
            for tooth in 0..TEETH {
                for flank in [-1., 1.] {
                    for step in 1..20 {
                        let r = base + (spec.tip_radius_mm - base) * f64::from(step) / 20.;
                        let target = wrap(
                            azimuth
                                - f64::from(tooth) * pitch
                                - flank * (spec.half_base_angle_rad - involute(r)),
                        );
                        // Solve phase(z) = target (mod a turn) in each band.
                        let k = spec.twist_rad_per_mm;
                        let bands: Vec<(f64, f64, [f64; 2])> = if spec.sweep == GearSweep::Helical {
                            vec![(0., k, [0., WIDTH])]
                        } else {
                            let length = band(&spec);
                            vec![
                                (0., k, [0., length]),
                                (k * WIDTH, -k, [WIDTH - length, WIDTH]),
                            ]
                        };
                        for (offset, slope, range) in bands {
                            let z = (target - offset) / slope;
                            if z > range[0] + 1e-3 && z < range[1] - 1e-3 {
                                expected += 1;
                                let gap = distance_to(&all, [sign * r, z]);
                                assert!(
                                    gap <= 0.01,
                                    "{name}: the flank law point [{}, {z}] is {gap} mm off the section",
                                    sign * r
                                );
                            }
                        }
                    }
                }
            }
        }
        assert!(
            expected > 20,
            "{name}: only {expected} law points fall in the band"
        );
        println!(
            "{{\"case\":\"{name} axis section\",\"flank_vertices\":{on_flank},\"worst_law_residual_mm\":{worst:e},\"law_points\":{expected}}}"
        );
    }
}

/// A box between the viewer and half of the gear hides exactly the curves
/// behind it: there every gear curve of the scene is hidden, elsewhere the
/// gear draws as it does alone, and no gear curve is lost or invented.
#[test]
fn a_box_in_front_hides_exactly_the_covered_half_of_a_twisted_gear() {
    for (name, spec) in [
        ("helical", spec(GearSweep::Helical, 0.05, None)),
        ("herringbone hub", spec(GearSweep::Chevron, -0.04, hub())),
        ("double-helical hub", spec(GAPPED, 0.05, hub())),
    ] {
        let solid = gear(spec);
        // The viewer of the side view is at negative Y.
        let cover = stock([20, 5, 30], [-20., -30., -10.]);
        let options = DrawingProjection {
            frame: side([0.; 3]),
            hidden_lines: true,
            ..Default::default()
        };
        let alone = draw(&[&solid], options);
        let scene = draw(&[&cover, &solid], options);
        let of_gear: Vec<&ProjectedVector> = scene.iter().filter(|v| v.source_index == 1).collect();
        let role = |vectors: &[ProjectedVector], role: VectorRole| -> Vec<ProjectedVector> {
            vectors.iter().filter(|v| v.role == role).cloned().collect()
        };
        let margin = 0.02;
        let mut worst: f64 = 0.;
        let mut covered = 0;
        for vector in &of_gear {
            for pair in vector.preview_mm.windows(2) {
                let mid = [
                    0.5 * (pair[0][0] + pair[1][0]),
                    0.5 * (pair[0][1] + pair[1][1]),
                ];
                if mid[0] < -margin {
                    covered += 1;
                    assert_eq!(
                        vector.role,
                        VectorRole::Hidden,
                        "{name}: {mid:?} is covered"
                    );
                    let all: Vec<&ProjectedVector> = alone.iter().collect();
                    worst = worst.max(distance_to(&all, mid));
                } else if mid[0] > margin {
                    let same = role(&alone, vector.role);
                    let same: Vec<&ProjectedVector> = same.iter().collect();
                    worst = worst.max(distance_to(&same, mid));
                }
            }
        }
        assert!(covered > 100, "{name}: only {covered} covered chords");
        for vector in &alone {
            for pair in vector.preview_mm.windows(2) {
                let mid = [
                    0.5 * (pair[0][0] + pair[1][0]),
                    0.5 * (pair[0][1] + pair[1][1]),
                ];
                if mid[0].abs() <= margin {
                    continue;
                }
                let found: Vec<&ProjectedVector> = of_gear
                    .iter()
                    .copied()
                    .filter(|v| mid[0] < 0. || v.role == vector.role)
                    .collect();
                worst = worst.max(distance_to(&found, mid));
            }
        }
        println!(
            "{{\"case\":\"{name} behind a box\",\"covered_chords\":{covered},\"worst_mm\":{worst:e}}}"
        );
        // Chords of either drawing stay within a quarter of the chord
        // tolerance of the same exact curves, split at different events.
        assert!(
            worst <= 0.5 * options.chord_tolerance_mm,
            "{name}: the scene differs from the gear alone by {worst} mm"
        );
    }
}

/// The solid hatch of a section is the cut region itself: its triangles lie
/// in the material of the exact solid and their total area matches a
/// scan-line quadrature of the oracle on the plane, for oblique and
/// lengthwise sections of every family, windows included.
#[test]
fn section_hatches_cover_exactly_the_cut_material() {
    use multicad_core::{drawings::VectorCurve, materials::HatchPattern};
    use multicad_kernel::DrawingHatch;
    let s = 0.5_f64.sqrt();
    let t = (1_f64 / 6.).sqrt();
    let oblique = PlaneFrame {
        origin: [0.3, -0.2, 2.2],
        u: [s, -s, 0.],
        v: [t, t, -2. * t],
    };
    let mut worst: f64 = 0.;
    for (name, spec) in gears() {
        let solid = gear(spec);
        let oracle = Oracle::new(spec);
        for (view, frame, window) in [
            ("side", side([0.4, 1.1, 2.]), None),
            ("oblique", oblique, None),
            ("oblique window", oblique, Some(([-6., 5.], [-4., 7.]))),
        ] {
            let vectors = project_drawing(
                &[DrawingSolid {
                    solid: &solid,
                    hatch: Some(DrawingHatch {
                        angle_rad: 0.,
                        spacing_mm: 1.,
                        line_width_mm: 0.18,
                        pattern: HatchPattern::Solid,
                    }),
                }],
                DrawingProjection {
                    frame,
                    section_offset_mm: Some(0.),
                    section_only: true,
                    section_u_mm: window.map(|w| w.0),
                    section_v_mm: window.map(|w| w.1),
                    ..Default::default()
                },
            )
            .unwrap();
            let triangles: Vec<[[f64; 2]; 3]> = vectors
                .iter()
                .filter_map(|vector| match vector.curve {
                    VectorCurve::HatchTriangle { points } => Some(points),
                    _ => None,
                })
                .collect();
            let area: f64 = triangles
                .iter()
                .map(|[a, b, c]| {
                    ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])).abs() / 2.
                })
                .sum();
            let inside = |p: [f64; 2]| {
                window.is_none_or(|([u0, u1], [v0, v1])| {
                    p[0] > u0 && p[0] < u1 && p[1] > v0 && p[1] < v1
                }) && oracle.inside(world(frame, p, 0.))
            };
            // Every triangle of real size lies in the material.
            let mut outside = 0;
            for [a, b, c] in &triangles {
                let size = ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])).abs();
                let centroid = [(a[0] + b[0] + c[0]) / 3., (a[1] + b[1] + c[1]) / 3.];
                if size > 1e-6 && !inside(centroid) {
                    outside += 1;
                }
            }
            // Scan lines across the plane, crossings bisected on the oracle;
            // the lines span the region's exact extent across them, found by
            // bisection, so that no line straddles its first or last edge.
            let reach = PITCH + 2.;
            let width = |y: f64| {
                let mut previous = (-reach, inside([-reach, y]));
                let mut entered = -reach;
                let mut total = 0.;
                for i in 1..=800 {
                    let x = -reach + 2. * reach * f64::from(i) / 800.;
                    let state = inside([x, y]);
                    if state != previous.1 {
                        let [mut a, mut b] = [previous.0, x];
                        for _ in 0..40 {
                            let mid = 0.5 * (a + b);
                            if inside([mid, y]) == previous.1 {
                                a = mid;
                            } else {
                                b = mid;
                            }
                        }
                        if state {
                            entered = 0.5 * (a + b);
                        } else {
                            total += 0.5 * (a + b) - entered;
                        }
                    }
                    previous = (x, state);
                }
                total
            };
            let coarse: Vec<f64> = (0..=200)
                .map(|k| -reach + 2. * reach * f64::from(k) / 200.)
                .collect();
            let filled: Vec<bool> = coarse.iter().map(|&y| width(y) > 0.).collect();
            let first = filled.iter().position(|f| *f).unwrap();
            let last = filled.iter().rposition(|f| *f).unwrap();
            let edge = |[mut outside, mut within]: [f64; 2]| {
                for _ in 0..40 {
                    let mid = 0.5 * (outside + within);
                    if width(mid) > 0. {
                        within = mid;
                    } else {
                        outside = mid;
                    }
                }
                0.5 * (outside + within)
            };
            let low = edge([coarse[first - 1], coarse[first]]);
            let high = edge([coarse[last + 1], coarse[last]]);
            const LINES: usize = 400;
            let dy = (high - low) / LINES as f64;
            let expected: f64 = (0..LINES)
                .map(|line| width(low + (line as f64 + 0.5) * dy) * dy)
                .sum();
            let error = (area - expected).abs() / expected;
            worst = worst.max(error);
            println!(
                "{{\"case\":\"{name} {view} hatch\",\"triangles\":{},\"area\":{area},\"quadrature\":{expected},\"outside\":{outside}}}",
                triangles.len()
            );
            assert!(expected > 1., "{name} {view}");
            assert!(
                outside == 0,
                "{name} {view}: {outside} triangles outside the material"
            );
            assert!(
                error <= 3e-3,
                "{name} {view}: hatch {area} against {expected}"
            );
        }
    }
    println!("{{\"worst_relative_hatch_area_error\":{worst:e}}}");
}

/// Depth of the nearest material along the sight line through a paper
/// point, from the oracle alone.
fn front(
    oracle: &Oracle,
    frame: PlaneFrame,
    paper: [f64; 2],
    [far, near]: [f64; 2],
) -> Option<f64> {
    const STEP: f64 = 0.04;
    let at = |depth: f64| oracle.inside(world(frame, paper, depth));
    let count = ((near - far) / STEP).ceil() as usize;
    let mut previous = near;
    for i in 1..=count {
        let depth = near - (near - far) * i as f64 / count as f64;
        if at(depth) {
            let [mut a, mut b] = [previous, depth];
            for _ in 0..40 {
                let mid = 0.5 * (a + b);
                if at(mid) {
                    b = mid;
                } else {
                    a = mid;
                }
            }
            return Some(0.5 * (a + b));
        }
        previous = depth;
    }
    None
}

/// Wherever the depth of the nearest material jumps between two sight
/// lines 1 µm apart, the drawing has a visible curve: every occluding
/// contour of the exact solid is drawn. The depth map comes from the oracle
/// alone, on a grid refined by bisection toward each jump.
#[test]
fn every_occluding_contour_of_the_exact_solid_is_a_visible_curve() {
    let s = 0.5_f64.sqrt();
    let t = (1_f64 / 6.).sqrt();
    for (name, spec, frame, grid) in [
        (
            "helical side",
            spec(GearSweep::Helical, 0.05, None),
            side([0., 0., 2.]),
            0.31,
        ),
        (
            "herringbone hub oblique",
            spec(GearSweep::Chevron, -0.04, hub()),
            PlaneFrame {
                origin: [0., 0., 2.],
                u: [s, -s, 0.],
                v: [t, t, -2. * t],
            },
            0.47,
        ),
        (
            "double-helical oblique",
            spec(GAPPED, 0.05, None),
            PlaneFrame {
                origin: [0., 0., 2.],
                u: [s, -s, 0.],
                v: [t, t, -2. * t],
            },
            0.43,
        ),
    ] {
        let solid = gear(spec);
        let oracle = Oracle::new(spec);
        let vectors = draw(
            &[&solid],
            DrawingProjection {
                frame,
                ..Default::default()
            },
        );
        let visible: Vec<&ProjectedVector> = vectors.iter().collect();
        assert!(visible.iter().all(|v| v.role == VectorRole::Visible));
        let reach = (PITCH + 1.).hypot(8.);
        let depths = [-reach, reach];
        let depth = |paper: [f64; 2]| front(&oracle, frame, paper, depths);
        let jump = |a: Option<f64>, b: Option<f64>| match (a, b) {
            (None, None) => false,
            (Some(a), Some(b)) => (a - b).abs() > 0.25,
            _ => true,
        };
        let count = (2. * reach / grid).ceil() as i64;
        let node = |i: i64, j: i64| [i as f64 * grid + 0.013, j as f64 * grid + 0.007];
        let map: Vec<Vec<Option<f64>>> = (-count..=count)
            .map(|i| (-count..=count).map(|j| depth(node(i, j))).collect())
            .collect();
        let get = |i: i64, j: i64| map[(i + count) as usize][(j + count) as usize];
        let mut contours = 0;
        let mut worst: f64 = 0.;
        for i in -count..count {
            for j in -count..count {
                for (di, dj) in [(1, 0), (0, 1)] {
                    let [mut a, mut b] = [node(i, j), node(i + di, j + dj)];
                    let [mut da, mut db] = [get(i, j), get(i + di, j + dj)];
                    if !jump(da, db) {
                        continue;
                    }
                    // Follow the jump down to a 1 µm pair of sight lines; a
                    // steep but continuous surface loses it on the way.
                    let mut real = true;
                    while (b[0] - a[0]).hypot(b[1] - a[1]) > 1e-3 {
                        let mid = [0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1])];
                        let dm = depth(mid);
                        if jump(da, dm) {
                            b = mid;
                            db = dm;
                        } else if jump(dm, db) {
                            a = mid;
                            da = dm;
                        } else {
                            real = false;
                            break;
                        }
                    }
                    if !real {
                        continue;
                    }
                    contours += 1;
                    let gap = distance_to(&visible, a);
                    worst = worst.max(gap);
                    assert!(
                        gap <= 0.03,
                        "{name}: the depth jumps from {da:?} to {db:?} at {a:?}, {gap} mm from every visible curve"
                    );
                }
            }
        }
        println!(
            "{{\"case\":\"{name}\",\"occluding_contour_crossings\":{contours},\"worst_distance_mm\":{worst:e}}}"
        );
        assert!(contours > 300, "{name}: only {contours} contour crossings");
    }
}

fn stock(size: [i64; 3], origin: [f64; 3]) -> Solid {
    Solid::rust_box(
        EntityId::new().uuid_bits(),
        size.map(|n| Rational::new(n, 1).unwrap()),
        PlaneFrame {
            origin,
            ..Default::default()
        },
        GeometryPolicy::default(),
    )
    .unwrap()
}

/// Dump drawings as polylines for visual inspection:
/// `SCENE_DUMP=<dir> SCENE_CASE=<filter> cargo test ... dump_twisted_drawings -- --ignored`.
#[test]
#[ignore = "visual inspection aid"]
fn dump_twisted_drawings() {
    let Ok(directory) = std::env::var("SCENE_DUMP") else {
        return;
    };
    let only = std::env::var("SCENE_CASE").ok();
    for (name, spec) in gears() {
        let solid = gear(spec);
        for (case, options) in cases(center(&solid)) {
            let label = format!("{name} {case}");
            if only
                .as_ref()
                .is_some_and(|only| !label.contains(only.as_str()))
            {
                continue;
            }
            let mut text = String::new();
            for vector in draw(&[&solid], options) {
                text.push_str(&format!("{:?}", vector.role));
                for p in &vector.preview_mm {
                    text.push_str(&format!(" {} {}", p[0], p[1]));
                }
                text.push('\n');
            }
            std::fs::write(
                std::path::Path::new(&directory).join(format!("{}.txt", label.replace(' ', "_"))),
                text,
            )
            .unwrap();
        }
    }
}

/// The largest distance from a preview point or chord midpoint of `from` to
/// the previews of `to` with the same source and, if required, role.
#[cfg(feature = "occt")]
fn deviation(
    from: &[ProjectedVector],
    to: &[ProjectedVector],
    same_role: bool,
) -> (f64, [f64; 2], VectorRole, usize) {
    let mut worst = (0., [0.; 2], VectorRole::Visible, 0);
    for vector in from {
        let spans: Vec<_> = to
            .iter()
            .filter(|other| {
                other.source_index == vector.source_index
                    && (!same_role || other.role == vector.role)
            })
            .flat_map(|other| other.preview_mm.windows(2).map(|w| [w[0], w[1]]))
            .collect();
        for pair in vector.preview_mm.windows(2) {
            let mid = [
                pair[0][0].midpoint(pair[1][0]),
                pair[0][1].midpoint(pair[1][1]),
            ];
            for [x, y] in [pair[0], mid, pair[1]] {
                let distance = spans
                    .iter()
                    .map(|&[[a, b], [c, d]]| {
                        let [u, v] = [c - a, d - b];
                        let along = (((x - a) * u + (y - b) * v) / (u * u + v * v).max(1e-300))
                            .clamp(0., 1.);
                        (x - a - along * u).hypot(y - b - along * v)
                    })
                    .fold(f64::INFINITY, f64::min);
                if distance > worst.0 {
                    worst = (distance, [x, y], vector.role, vector.source_index);
                }
            }
        }
    }
    worst
}

/// Measure every case against OpenCascade's drawing of the independently
/// built native gear: `SCENE_CASE=<filter> cargo test --features occt ...
/// report_native_deviations -- --ignored --nocapture`.
#[cfg(feature = "occt")]
#[test]
#[ignore = "measurement report"]
fn report_native_deviations() {
    let only = std::env::var("SCENE_CASE").ok();
    let dump = std::env::var("SCENE_DUMP").ok();
    for (name, spec) in gears() {
        let solid = gear(spec);
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("gear.step");
        solid.write_step(&path).unwrap();
        let native = if std::env::var("SCENE_BUILT").is_ok() {
            Solid::involute_gear_named(spec, EntityId::new(), GeometryPolicy::default()).unwrap()
        } else {
            Solid::read_step(&path, solid.geometry_policy()).unwrap()
        };
        for (case, options) in cases(center(&solid)) {
            let label = format!("{name} {case}");
            if only
                .as_ref()
                .is_some_and(|only| !label.contains(only.as_str()))
            {
                continue;
            }
            let rust = draw(&[&solid], options);
            let complete = draw(
                &[&solid],
                DrawingProjection {
                    tangent_edges: true,
                    hidden_lines: true,
                    ..options
                },
            );
            let reference = match project_drawing(
                &[DrawingSolid {
                    solid: &native,
                    hatch: None,
                }],
                options,
            ) {
                Ok(reference) => reference,
                Err(error) => {
                    println!(
                        "{{\"case\":\"{label}\",\"native_error\":\"{}\"}}",
                        error.message
                    );
                    continue;
                }
            };
            let forward = deviation(&rust, &reference, true);
            let any = deviation(&rust, &reference, false);
            let backward = deviation(&reference, &complete, false);
            let backward_role = deviation(&reference, &rust, true);
            println!(
                "{{\"case\":\"{label}\",\"rust\":{},\"native\":{},\"rust_to_native_role\":{:?},\"rust_to_native_any\":{:.4},\"native_to_rust_any\":{:?},\"native_to_rust_role\":{:.4}}}",
                rust.len(),
                reference.len(),
                forward,
                any.0,
                backward,
                backward_role.0
            );
            if let Some(directory) = &dump {
                for (kind, vectors) in [("rust", &rust), ("native", &reference)] {
                    let mut text = String::new();
                    for vector in vectors {
                        text.push_str(&format!("{:?}", vector.role));
                        for p in &vector.preview_mm {
                            text.push_str(&format!(" {} {}", p[0], p[1]));
                        }
                        text.push('\n');
                    }
                    std::fs::write(
                        std::path::Path::new(directory)
                            .join(format!("{}.{kind}.txt", label.replace(' ', "_"))),
                        text,
                    )
                    .unwrap();
                }
            }
        }
    }
}

/// Depth of the nearest retained material along the sight line through a
/// paper point, by an independent membership predicate of world points.
#[cfg(feature = "occt")]
fn front_of(
    inside: &dyn Fn([f64; 3]) -> bool,
    frame: PlaneFrame,
    paper: [f64; 2],
    [far, near]: [f64; 2],
) -> Option<f64> {
    // Fine enough to find a thin wedge of retained material at a cut.
    let step = (0.02_f64).min((near - far) / 2000.);
    let at = |depth: f64| inside(world(frame, paper, depth));
    let count = ((near - far) / step).ceil() as usize;
    let mut previous = near;
    for i in 0..=count {
        let depth = near - (near - far) * i as f64 / count as f64;
        let depth = if i == count { far + 1e-9 } else { depth };
        if at(depth) {
            if i == 0 {
                return Some(depth);
            }
            let [mut a, mut b] = [previous, depth];
            for _ in 0..40 {
                let mid = 0.5 * (a + b);
                if at(mid) {
                    b = mid;
                } else {
                    a = mid;
                }
            }
            return Some(0.5 * (a + b));
        }
        previous = depth;
    }
    None
}

/// Whether a visible chord through `mid` along `along` is an edge of the
/// retained solid as the viewer sees it, by the predicate alone: across the
/// chord the depth of the nearest retained material jumps, appears or
/// disappears (an occluding contour or a rim) or changes slope (a crease).
#[cfg(feature = "occt")]
fn seen_edge(
    inside: &dyn Fn([f64; 3]) -> bool,
    options: DrawingProjection,
    mid: [f64; 2],
    along: [f64; 2],
) -> bool {
    let reach = 40.;
    let near = options.section_offset_mm.unwrap_or(reach);
    let far = options
        .section_offset_mm
        .zip(options.section_depth_mm)
        .map_or(-reach, |(top, depth)| top - depth);
    let length = along[0].hypot(along[1]);
    let normal = [-along[1] / length, along[0] / length];
    let delta = 0.01;
    let depth = |k: f64| {
        front_of(
            inside,
            options.frame,
            [
                mid[0] + k * delta * normal[0],
                mid[1] + k * delta * normal[1],
            ],
            [far, near],
        )
    };
    match [-2., -1., 1., 2.].map(depth) {
        [_, None, None, _] => false,
        [_, Some(_), None, _] | [_, None, Some(_), _] => true,
        [m2, Some(m1), Some(p1), p2] => {
            if (m1 - p1).abs() > 0.05 {
                return true;
            }
            match (m2, p2) {
                (Some(m2), Some(p2)) => ((p2 - p1) - (m1 - m2)).abs() / delta > 0.1,
                _ => true,
            }
        }
    }
}

/// Compare one drawing of Rust solids with OpenCascade's drawing of their
/// STEP imports by the rules of `twisted_gears_match_native_open_cascade`.
/// A Rust chord without any native curve of its source is accepted only
/// where `certify` holds; returns how many were.
#[cfg(feature = "occt")]
fn compare_native(
    label: &str,
    (solids, natives): (&[&Solid], &[&Solid]),
    options: DrawingProjection,
    certify: &dyn Fn(VectorRole, [f64; 2], [f64; 2]) -> bool,
    (failures, worst): (&mut Vec<String>, &mut [f64; 4]),
) -> usize {
    let rust = draw(solids, options);
    let complete = draw(
        solids,
        DrawingProjection {
            tangent_edges: true,
            hidden_lines: true,
            ..options
        },
    );
    let reference = draw(natives, options);
    let bound = 3. * options.chord_tolerance_mm;
    let of = |vectors: &[ProjectedVector], role: VectorRole| -> Vec<ProjectedVector> {
        vectors.iter().filter(|v| v.role == role).cloned().collect()
    };
    let mut certified = 0;
    let mut unmatched = Vec::new();
    let matched: Vec<ProjectedVector> = rust
        .iter()
        .filter(|vector| {
            let found = deviation(std::slice::from_ref(vector), &reference, false).0 <= bound;
            if !found {
                unmatched.push((*vector).clone());
            }
            found
        })
        .cloned()
        .collect();
    for vector in &unmatched {
        let all: Vec<&ProjectedVector> = reference
            .iter()
            .filter(|other| other.source_index == vector.source_index)
            .collect();
        for pair in vector.preview_mm.windows(2) {
            let mid = [
                0.5 * (pair[0][0] + pair[1][0]),
                0.5 * (pair[0][1] + pair[1][1]),
            ];
            if distance_to(&all, mid) <= bound {
                continue;
            }
            let along = [pair[1][0] - pair[0][0], pair[1][1] - pair[0][1]];
            if certify(vector.role, mid, along) {
                certified += 1;
            } else {
                failures.push(format!(
                    "{label}: Rust {:?} of source {} at {mid:?} has no native curve",
                    vector.role, vector.source_index
                ));
            }
        }
    }
    // Visible Rust chords off every native visible curve: OpenCascade drew
    // them hidden or not at all; each must be certified.
    let mut visible = of(&matched, VectorRole::Visible);
    let native_visible = of(&reference, VectorRole::Visible);
    visible.retain(|vector| {
        let near: Vec<&ProjectedVector> = native_visible
            .iter()
            .filter(|other| other.source_index == vector.source_index)
            .collect();
        vector.preview_mm.windows(2).all(|pair| {
            let mid = [
                0.5 * (pair[0][0] + pair[1][0]),
                0.5 * (pair[0][1] + pair[1][1]),
            ];
            let along = [pair[1][0] - pair[0][0], pair[1][1] - pair[0][1]];
            if [pair[0], mid, pair[1]]
                .into_iter()
                .all(|point| distance_to(&near, point) <= bound)
            {
                return true;
            }
            if certify(VectorRole::Visible, mid, along) {
                certified += 1;
            } else {
                failures.push(format!(
                    "{label}: Rust visible chord of source {} at {mid:?} is not a seen edge",
                    vector.source_index
                ));
            }
            false
        })
    });
    let checks = [
        deviation(&matched, &reference, false),
        deviation(&visible, &reference, true),
        deviation(&of(&rust, VectorRole::Section), &reference, true),
        deviation(&reference, &complete, false),
        deviation(&of(&reference, VectorRole::Section), &rust, true),
    ];
    println!(
        "{{\"case\":\"{label}\",\"rust_to_native_mm\":{},\"rust_visible_to_native_visible_mm\":{},\"rust_section_to_native_section_mm\":{},\"native_to_rust_mm\":{},\"native_section_to_rust_section_mm\":{},\"uncounted_rust_curves\":{}}}",
        checks[0].0,
        checks[1].0,
        checks[2].0,
        checks[3].0,
        checks[4].0,
        unmatched.len()
    );
    for (index, check) in checks.iter().enumerate() {
        worst[index.min(3)] = worst[index.min(3)].max(check.0);
        if check.0 > bound {
            failures.push(format!("{label}: check {index} {check:?}"));
        }
    }
    certified
}

/// Overlapping scenes of a box and a twisted gear against OpenCascade: the
/// box passes through the tooth band, so every view mixes occlusion by both
/// solids and sections cut both.
#[cfg(feature = "occt")]
#[test]
fn overlapping_box_and_twisted_gear_scenes_match_native_open_cascade() {
    let only = std::env::var("SCENE_CASE").ok();
    let mut failures = Vec::new();
    let mut worst = [0_f64; 4];
    let mut compared = 0;
    let mut certified = 0;
    let mut unavailable = 0;
    for (name, spec) in [
        ("herringbone", spec(GearSweep::Chevron, 0.05, None)),
        ("double-helical", spec(GAPPED, -0.04, hub())),
    ] {
        let solid = gear(spec);
        let oracle = Oracle::new(spec);
        let cover = stock([6, 30, 3], [-2., -15., 1.5]);
        let directory = tempfile::tempdir().unwrap();
        let natives: Vec<Solid> = [&cover, &solid]
            .iter()
            .enumerate()
            .map(|(index, solid)| {
                let path = directory.path().join(format!("solid-{index}.step"));
                solid.write_step(&path).unwrap();
                Solid::read_step(&path, solid.geometry_policy()).unwrap()
            })
            .collect();
        for (case, options) in cases([0., 0., 2.]) {
            let label = format!("box + {name} {case}");
            if only
                .as_ref()
                .is_some_and(|only| !label.contains(only.as_str()))
            {
                continue;
            }
            let native: Vec<DrawingSolid> = natives
                .iter()
                .map(|solid| DrawingSolid { solid, hatch: None })
                .collect();
            if let Err(error) = project_drawing(&native, options) {
                println!(
                    "{{\"case\":\"{label}\",\"native_error\":\"{}\"}}",
                    error.message
                );
                unavailable += 1;
                continue;
            }
            compared += 1;
            let inside = |p: [f64; 3]| {
                (p[0] > -2. && p[0] < 4. && p[1] > -15. && p[1] < 15. && p[2] > 1.5 && p[2] < 4.5)
                    || oracle.inside(p)
            };
            certified += compare_native(
                &label,
                (&[&cover, &solid], &[&natives[0], &natives[1]]),
                options,
                &|role, point, along| {
                    role == VectorRole::Visible && seen_edge(&inside, options, point, along)
                },
                (&mut failures, &mut worst),
            );
        }
    }
    println!(
        "{{\"cases\":{compared},\"native_unavailable\":{unavailable},\"worst_mm\":{worst:?},\"oracle_certified_chords\":{certified}}}"
    );
    assert!(failures.is_empty(), "{failures:#?}");
    assert!(
        only.is_some() || compared >= 30,
        "{compared} cases compared"
    );
}

/// Every non-axial, cropped and sectioned view of the twisted gears against
/// OpenCascade 7.9.3 hidden-line removal and sections of the same solids,
/// read from the Rust STEP export, within three chord tolerances:
///
/// - every Rust curve lies on a native curve, visible and section curves on
///   native curves of the same role;
/// - every native curve lies on a Rust curve of the drawing with hidden
///   lines and tangent-continuous seams, and native section curves on Rust
///   section curves.
///
/// Measured deviations of the native reference, none of which widens the
/// bound:
///
/// - the imported flank seams are not marked smooth, so OpenCascade draws
///   them as sharp edges; they are compared with the Rust seams;
/// - OpenCascade publishes rear flank edges, seams and the bore outline of a
///   detail view as visible although material lies in front (the member's
///   sight intervals are certified against exact membership in the unit
///   tests); a hidden Rust curve may lie on a native curve of either role;
/// - in the finite-depth oblique sections of the hubbed gears OpenCascade
///   omits a stretch of the far cut's rim on one tooth (up to 0.34 mm from
///   its other curves). Such a Rust curve is accepted only where the oracle
///   certifies it as a rim of the far cut seen against empty space.
///
/// The uncropped axial views stay with the axial provider and its own
/// regressions; OpenCascade fails on some axial herringbone sections with
/// "Drawing curve has no finite trim domain".
#[cfg(feature = "occt")]
#[test]
fn twisted_gears_match_native_open_cascade() {
    let only = std::env::var("SCENE_CASE").ok();
    let mut failures = Vec::new();
    let mut worst = [0_f64; 4];
    let mut certified = 0;
    let mut compared = 0;
    for (name, spec) in gears() {
        let solid = gear(spec);
        let oracle = Oracle::new(spec);
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("gear.step");
        solid.write_step(&path).unwrap();
        let native = Solid::read_step(&path, solid.geometry_policy()).unwrap();
        for (case, options) in cases(center(&solid)) {
            let label = format!("{name} {case}");
            let cropped = options.detail.is_some() || options.section_u_mm.is_some();
            if only
                .as_ref()
                .is_some_and(|only| !label.contains(only.as_str()))
                || (case.starts_with("axial")
                    && !cropped
                    && !matches!(spec.sweep, GearSweep::GappedChevron { .. }))
            {
                continue;
            }
            compared += 1;
            certified += compare_native(
                &label,
                (&[&solid], &[&native]),
                options,
                &|role, point, along| {
                    role == VectorRole::Visible
                        && seen_edge(&|p| oracle.inside(p), options, point, along)
                },
                (&mut failures, &mut worst),
            );
        }
    }
    println!(
        "{{\"cases\":{compared},\"worst_mm\":{worst:?},\"oracle_certified_chords\":{certified}}}"
    );
    assert!(failures.is_empty(), "{failures:#?}");
    assert!(
        only.is_some() || compared >= 80,
        "{compared} cases compared"
    );
}
