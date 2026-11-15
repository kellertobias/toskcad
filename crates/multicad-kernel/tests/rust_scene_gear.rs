//! Analytic regressions for spur gears and internal rings drawn through the
//! exact-curve scene engine (CAD-346). Expectations are derived here from the
//! exact tooth profile, independently of the engine: a dense sampling of the
//! profile decides occlusion along each sight line.
use multicad_core::{
    document::EntityId,
    drawings::{VectorCurve, VectorRole},
    materials::HatchPattern,
};
use multicad_kernel::{
    DrawingHatch, DrawingProjection, DrawingSolid, GearHub, GearSweep, GeometryErrorKind,
    GeometryPolicy, InvoluteGear, PlaneFrame, ProjectedVector, RigidTransform, Solid,
    project_drawing,
    rust_kernel::spur_gear::{Edge, Profile, Spec},
};
use std::f64::consts::PI;

const TEETH: u32 = 24;
const WIDTH: f64 = 4.;
const BORE: f64 = 3.;

fn spec() -> InvoluteGear {
    let a = 20_f64.to_radians();
    let pitch = 12.;
    InvoluteGear {
        hub: None,
        teeth: TEETH,
        pitch_radius_mm: pitch,
        base_radius_mm: pitch * a.cos(),
        tip_radius_mm: pitch + 1.,
        root_radius_mm: pitch - 1.25,
        root_fillet_radius_mm: 0.2,
        half_base_angle_rad: (PI / 2. - 0.05) / (2. * pitch) + a.tan() - a,
        face_width_mm: WIDTH,
        bore_radius_mm: BORE,
        twist_rad_per_mm: 0.,
        sweep: GearSweep::Spur,
    }
}

fn gear(hub: Option<GearHub>) -> Solid {
    Solid::rust_spur_gear_named(
        InvoluteGear { hub, ..spec() },
        EntityId::new(),
        GeometryPolicy::default(),
    )
    .unwrap()
}

/// The exact profile, built directly from the pure kernel.
fn profile() -> Profile {
    let s = spec();
    Profile::new(
        Spec {
            teeth: s.teeth,
            pitch_radius_mm: s.pitch_radius_mm,
            base_radius_mm: s.base_radius_mm,
            tip_radius_mm: s.tip_radius_mm,
            root_radius_mm: s.root_radius_mm,
            root_fillet_radius_mm: s.root_fillet_radius_mm,
            half_base_angle_rad: s.half_base_angle_rad,
            face_width_mm: s.face_width_mm,
            bore_radius_mm: s.bore_radius_mm,
        },
        1e-5,
        1e7,
    )
    .unwrap()
}

const SAMPLES: usize = 4000;

/// Dense closed polylines of the outline and the bore.
fn rings(profile: &Profile) -> Vec<Vec<[f64; 2]>> {
    let mut outline = Vec::new();
    for edge in profile.outline() {
        outline.extend((0..SAMPLES).map(|i| edge.edge.point(i as f64 / SAMPLES as f64)));
    }
    let bore = Edge::Planar(profile.bore().unwrap());
    vec![
        outline,
        (0..SAMPLES)
            .map(|i| bore.point(i as f64 / SAMPLES as f64))
            .collect(),
    ]
}

/// Material intervals of the line `x = at` through the sampled profile.
fn column(rings: &[Vec<[f64; 2]>], at: f64) -> Vec<[f64; 2]> {
    let mut stops = Vec::new();
    for ring in rings {
        for (a, b) in ring.iter().zip(ring.iter().cycle().skip(1)) {
            if (a[0] < at) != (b[0] < at) {
                stops.push(a[1] + (b[1] - a[1]) * (at - a[0]) / (b[0] - a[0]));
            }
        }
    }
    stops.sort_by(f64::total_cmp);
    stops.chunks(2).map(|pair| [pair[0], pair[1]]).collect()
}

fn draw(solids: &[&Solid], options: DrawingProjection) -> Vec<ProjectedVector> {
    let inputs: Vec<_> = solids
        .iter()
        .map(|&solid| DrawingSolid {
            solid,
            hatch: Some(DrawingHatch {
                angle_rad: 0.37,
                spacing_mm: 0.25,
                line_width_mm: 0.18,
                pattern: HatchPattern::Lines,
            }),
        })
        .collect();
    project_drawing(&inputs, options).unwrap()
}

/// Looks along +Y: paper u is world X, paper v is world Z, the viewer is at
/// negative Y.
fn side() -> PlaneFrame {
    PlaneFrame {
        origin: [0.; 3],
        u: [1., 0., 0.],
        v: [0., 0., 1.],
    }
}

fn line(vector: &ProjectedVector) -> Option<[[f64; 2]; 2]> {
    match vector.curve {
        VectorCurve::Line { a, b } => Some([a, b]),
        _ => None,
    }
}

/// Stations of the full-height generators with a role, ascending.
fn generators(vectors: &[ProjectedVector], role: VectorRole) -> Vec<f64> {
    let mut stations: Vec<f64> = vectors
        .iter()
        .filter(|vector| vector.role == role)
        .filter_map(line)
        .filter(|[a, b]| (a[0] - b[0]).abs() < 1e-9)
        .map(|[a, b]| {
            assert!((a[1].min(b[1])).abs() < 1e-9 && (a[1].max(b[1]) - WIDTH).abs() < 1e-9);
            a[0]
        })
        .collect();
    stations.sort_by(f64::total_cmp);
    stations
}

fn same(found: &[f64], expected: &[f64], label: &str) {
    assert_eq!(
        found.len(),
        expected.len(),
        "{label}: {found:?} {expected:?}"
    );
    for (a, b) in found.iter().zip(expected) {
        assert!((a - b).abs() < 2e-6, "{label}: {a} {b}");
    }
}

/// Every sharp profile vertex and every point where the profile's tangent is
/// along the sight is a generator; it is visible exactly when no material
/// lies between it and the viewer. Coincident stations publish once, visible
/// winning. The bore's two outlines are hidden inside the body.
#[test]
fn side_view_generators_follow_the_exact_profile_and_its_occlusion() {
    let profile = profile();
    let rings = rings(&profile);
    let mut candidates: Vec<[f64; 2]> = Vec::new();
    let edges: Vec<Edge> = profile.outline().iter().map(|edge| edge.edge).collect();
    for (index, edge) in edges.iter().enumerate() {
        let next = edges[(index + 1) % edges.len()];
        let [a, b] = [edge.tangent(1.), next.tangent(0.)];
        if (a[0] * b[1] - a[1] * b[0]).abs() > 1e-9 {
            candidates.push(edge.point(1.));
        }
        // Tangent along Y: the X coordinate turns. Bisect the exact tangent.
        let steps = 61;
        for step in 0..steps {
            let [mut lo, mut hi] = [step as f64 / steps as f64, (step + 1) as f64 / steps as f64];
            let sign = |t: f64| edge.tangent(t)[0] > 0.;
            if sign(lo + 1e-12) == sign(hi - 1e-12) {
                continue;
            }
            for _ in 0..60 {
                let mid = 0.5 * (lo + hi);
                if sign(mid) == sign(lo + 1e-12) {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            if lo > 1e-9 && hi < 1. - 1e-9 {
                candidates.push(edge.point(lo));
            }
        }
    }
    candidates.extend([[-BORE, 0.], [BORE, 0.]]);
    let mut visible: Vec<f64> = Vec::new();
    let mut hidden: Vec<f64> = Vec::new();
    for [x, y] in candidates {
        // The viewer is at negative Y. Probe both sides of a tangent contact.
        let covered = [-1e-7, 1e-7].into_iter().all(|shift| {
            column(&rings, x + shift)
                .into_iter()
                .any(|[a, b]| b.min(y) - a > 1e-4)
        });
        if covered { &mut hidden } else { &mut visible }.push(x);
    }
    let unique = |mut list: Vec<f64>| {
        list.sort_by(f64::total_cmp);
        list.dedup_by(|a, b| (*a - *b).abs() < 1e-6);
        list
    };
    let visible = unique(visible);
    let hidden: Vec<f64> = unique(hidden)
        .into_iter()
        .filter(|x| !visible.iter().any(|v| (v - x).abs() < 1e-6))
        .collect();
    assert!(visible.len() > 20 && hidden.len() > 10);

    let solid = gear(None);
    let plain = draw(
        &[&solid],
        DrawingProjection {
            frame: side(),
            ..Default::default()
        },
    );
    assert!(
        plain
            .iter()
            .all(|vector| vector.role == VectorRole::Visible)
    );
    same(
        &generators(&plain, VectorRole::Visible),
        &visible,
        "visible",
    );
    let full = draw(
        &[&solid],
        DrawingProjection {
            frame: side(),
            hidden_lines: true,
            ..Default::default()
        },
    );
    same(
        &generators(&full, VectorRole::Visible),
        &visible,
        "visible with hidden lines",
    );
    same(&generators(&full, VectorRole::Hidden), &hidden, "hidden");
    // The drawing is the tip-to-tip rectangle of the tooth band.
    let reach = rings[0].iter().map(|p| p[0]).fold(0., f64::max);
    for (k, expected) in [[-reach, reach], [0., WIDTH]].into_iter().enumerate() {
        let values = plain
            .iter()
            .flat_map(|vector| vector.preview_mm.iter().map(|p| p[k]));
        let lo = values.clone().fold(f64::INFINITY, f64::min);
        let hi = values.fold(f64::NEG_INFINITY, f64::max);
        assert!((lo - expected[0]).abs() < 1e-5 && (hi - expected[1]).abs() < 1e-5);
    }
}

/// A section plane parallel to the axis cuts rectangles whose widths are the
/// exact chords of the profile at that offset: tooth cut widths 1 mm inside
/// the tip circle, and both walls through the bore.
#[test]
fn side_sections_are_the_exact_profile_chords_and_hatch_only_material() {
    let profile = profile();
    let rings = rings(&profile);
    let solid = gear(None);
    for offset in [12., 0.] {
        // The plane y = offset; turning the profile a quarter maps it to x.
        let turned: Vec<Vec<[f64; 2]>> = rings
            .iter()
            .map(|ring| ring.iter().map(|[x, y]| [*y, *x]).collect())
            .collect();
        let chords = column(&turned, offset);
        assert_eq!(chords.len(), if offset == 0. { 2 } else { 3 });
        let vectors = draw(
            &[&solid],
            DrawingProjection {
                frame: PlaneFrame {
                    origin: [0., offset, 0.],
                    ..side()
                },
                section_offset_mm: Some(0.),
                section_only: true,
                ..Default::default()
            },
        );
        let mut walls: Vec<f64> = vectors
            .iter()
            .filter(|vector| vector.role == VectorRole::Section)
            .filter_map(line)
            .filter(|[a, b]| (a[0] - b[0]).abs() < 1e-9)
            .map(|[a, _]| a[0])
            .collect();
        walls.sort_by(f64::total_cmp);
        let expected: Vec<f64> = chords.iter().flatten().copied().collect();
        same(&walls, &expected, "section walls");
        let section = vectors
            .iter()
            .filter(|vector| vector.role == VectorRole::Section)
            .count();
        assert_eq!(section, 4 * chords.len());
        let strokes: Vec<_> = vectors
            .iter()
            .filter(|vector| vector.role == VectorRole::Hatch)
            .filter_map(line)
            .collect();
        assert!(!strokes.is_empty());
        for [a, b] in strokes {
            for q in [0.01, 0.5, 0.99] {
                let p = [a[0] + (b[0] - a[0]) * q, a[1] + (b[1] - a[1]) * q];
                assert!(p[1] > -1e-9 && p[1] < WIDTH + 1e-9);
                assert!(
                    chords
                        .iter()
                        .any(|[lo, hi]| p[0] > lo - 1e-6 && p[0] < hi + 1e-6),
                    "hatch leaves the material at {p:?}"
                );
            }
        }
    }
}

/// A plane tilted by `tilt` about X through the band cuts the affine image of
/// the whole profile: its area is the profile area over cos(tilt), measured
/// here by the solid hatch triangulation of the published region.
#[test]
fn oblique_section_area_is_the_profile_area_over_the_cosine() {
    let tilt = 0.12_f64;
    let (sin, cos) = tilt.sin_cos();
    let solid = gear(None);
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
            frame: PlaneFrame {
                origin: [0., 0., WIDTH / 2.],
                u: [1., 0., 0.],
                v: [0., cos, sin],
            },
            section_offset_mm: Some(0.),
            section_only: true,
            ..Default::default()
        },
    )
    .unwrap();
    let area: f64 = vectors
        .iter()
        .filter_map(|vector| match vector.curve {
            VectorCurve::HatchTriangle { points: [a, b, c] } => {
                Some(((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])).abs() / 2.)
            }
            _ => None,
        })
        .sum();
    let exact = profile().area_mm2() / cos;
    assert!((area - exact).abs() < 2e-4 * exact, "{area} {exact}");
    // The bore becomes an ellipse with semi-axes bore and bore / cos.
    let rims: Vec<_> = vectors
        .iter()
        .filter_map(|vector| match vector.curve {
            VectorCurve::Ellipse {
                major_mm, minor_mm, ..
            } if major_mm > 2. && major_mm < 4. => Some([major_mm, minor_mm]),
            _ => None,
        })
        .collect();
    assert!(!rims.is_empty());
    for [major, minor] in rims {
        assert!((major - BORE / cos).abs() < 1e-9 && (minor - BORE).abs() < 1e-9);
    }
}

/// A plane across the upper hub cuts the hub ring; the tooth band behind it
/// stays visible and the bore, continuous through hub and band, has no rim
/// at their junction.
#[test]
fn hub_sections_cut_the_hub_ring_in_front_of_the_tooth_band() {
    let solid = gear(Some(GearHub {
        radius_mm: 6.,
        axial_start_mm: -1.5,
        axial_end_mm: 6.,
    }));
    let vectors = draw(
        &[&solid],
        DrawingProjection {
            frame: PlaneFrame {
                origin: [0., 0., 5.],
                ..Default::default()
            },
            section_offset_mm: Some(0.),
            hidden_lines: true,
            ..Default::default()
        },
    );
    let circles = |role: VectorRole| {
        let mut radii: Vec<f64> = vectors
            .iter()
            .filter(|vector| vector.role == role)
            .filter_map(|vector| match vector.curve {
                VectorCurve::Ellipse {
                    major_mm, minor_mm, ..
                } => {
                    assert!((major_mm - minor_mm).abs() < 1e-12);
                    Some(major_mm)
                }
                _ => None,
            })
            .collect();
        radii.sort_by(f64::total_cmp);
        radii.dedup_by(|a, b| (*a - *b).abs() < 1e-9);
        radii
    };
    assert_eq!(circles(VectorRole::Section), [3., 6.]);
    // Visible: the cut rims, and the tip and root arcs of the band behind.
    let visible = circles(VectorRole::Visible);
    for radius in [3., 6., 10.75, 13.] {
        assert!(
            visible.iter().any(|r| (r - radius).abs() < 1e-9),
            "{visible:?}"
        );
    }
    // The lower hub's rims coincide with the visible cut rim in this view and
    // publish once; the bore has no hidden rim at either junction.
    assert!(circles(VectorRole::Hidden).is_empty());
    for vector in vectors.iter().filter(|v| v.role == VectorRole::Hatch) {
        let [a, b] = line(vector).unwrap();
        let r = (0.5 * (a[0] + b[0])).hypot(0.5 * (a[1] + b[1]));
        assert!(r > 3. - 1e-6 && r < 6. + 1e-6);
    }
}

/// A box in front of the gear hides every generator exactly up to its top
/// face, and the detail circle and the lateral window bound the drawing.
#[test]
fn a_box_in_front_hides_generators_up_to_its_top_and_crops_bound_the_drawing() {
    let top = 1.25;
    let cover = Solid::box_solid([40., 10., 6.], GeometryPolicy::default())
        .unwrap()
        .transformed(RigidTransform {
            rotation: [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
            translation_mm: [-20., -30., top - 6.],
        })
        .unwrap();
    let solid = gear(None);
    let vectors = draw(
        &[&cover, &solid],
        DrawingProjection {
            frame: side(),
            hidden_lines: true,
            ..Default::default()
        },
    );
    let mut seen = 0;
    for [a, b] in vectors
        .iter()
        .filter(|v| v.source_index == 1 && v.role == VectorRole::Visible)
        .filter_map(line)
        .filter(|[a, b]| (a[0] - b[0]).abs() < 1e-9)
    {
        assert!((a[1].min(b[1]) - top).abs() < 1e-9 && (a[1].max(b[1]) - WIDTH).abs() < 1e-9);
        seen += 1;
    }
    assert!(seen > TEETH);
    // The same stations continue hidden behind the box down to the cap.
    assert!(
        vectors
            .iter()
            .filter(|v| v.source_index == 1 && v.role == VectorRole::Hidden)
            .filter_map(line)
            .any(|[a, b]| (a[0] - b[0]).abs() < 1e-9
                && a[1].min(b[1]).abs() < 1e-9
                && (a[1].max(b[1]) - top).abs() < 1e-9)
    );

    let detail = ([9., 2.], 2.5);
    let cropped = draw(
        &[&solid],
        DrawingProjection {
            frame: side(),
            detail: Some(detail),
            ..Default::default()
        },
    );
    assert!(!cropped.is_empty());
    for p in cropped.iter().flat_map(|vector| &vector.preview_mm) {
        assert!((p[0] - detail.0[0]).hypot(p[1] - detail.0[1]) <= detail.1 + 1e-9);
    }
    let window = draw(
        &[&solid],
        DrawingProjection {
            frame: PlaneFrame {
                origin: [0., 0., 2.],
                ..side()
            },
            section_offset_mm: Some(0.),
            section_u_mm: Some([4., 12.5]),
            section_v_mm: Some([-1., 1.5]),
            section_depth_mm: Some(6.),
            ..Default::default()
        },
    );
    for role in [VectorRole::Visible, VectorRole::Section, VectorRole::Hatch] {
        assert!(window.iter().any(|vector| vector.role == role));
    }
    for p in window.iter().flat_map(|vector| &vector.preview_mm) {
        assert!(
            p[0] >= 4. - 1e-9 && p[0] <= 12.5 + 1e-9 && p[1] >= -1. - 1e-9 && p[1] <= 1.5 + 1e-9
        );
    }
    // The window's lateral faces are drawn where material lies behind them:
    // the left border along the whole window height.
    assert!(
        window
            .iter()
            .filter(|v| v.role == VectorRole::Visible)
            .filter_map(line)
            .any(|[a, b]| (a[0] - 4.).abs() < 1e-9
                && (b[0] - 4.).abs() < 1e-9
                && (a[1] - b[1]).abs() > 2.5 - 1e-9)
    );
}

/// A plane exactly on the junction of the tooth band and the upper hub: the
/// cut face is the hub ring, and the junction circle keeps separating it from
/// the coplanar cap of the band, as in the native common with a half-space.
/// The whole cap is the section face; the bore stays one continuous hole.
#[test]
fn a_plane_on_the_hub_junction_keeps_the_junction_circle() {
    let solid = gear(Some(GearHub {
        radius_mm: 6.,
        axial_start_mm: -1.5,
        axial_end_mm: 6.,
    }));
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
            frame: PlaneFrame {
                origin: [0., 0., WIDTH],
                ..Default::default()
            },
            section_offset_mm: Some(0.),
            ..Default::default()
        },
    )
    .unwrap();
    for role in [VectorRole::Section, VectorRole::Visible] {
        let mut sweeps = [0.; 2];
        for vector in vectors.iter().filter(|vector| vector.role == role) {
            if let VectorCurve::Ellipse {
                major_mm,
                start_rad,
                end_rad,
                ..
            } = vector.curve
            {
                for (sweep, radius) in sweeps.iter_mut().zip([BORE, 6.]) {
                    if (major_mm - radius).abs() < 1e-9 {
                        *sweep += end_rad - start_rad;
                    }
                }
            }
        }
        // Each circle once: the bore rim and the junction circle.
        for sweep in sweeps {
            assert!((sweep - 2. * PI).abs() < 1e-9, "{role:?}: {sweeps:?}");
        }
    }
    let area: f64 = vectors
        .iter()
        .filter_map(|vector| match vector.curve {
            VectorCurve::HatchTriangle { points: [a, b, c] } => {
                Some(((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])).abs() / 2.)
            }
            _ => None,
        })
        .sum();
    let exact = profile().area_mm2();
    assert!((area - exact).abs() < 2e-4 * exact, "{area} {exact}");
}

/// A section-only plane lying in either cap of the tooth band meets the gear
/// in that cap, as the native common of the solid and the plane does.
#[test]
fn a_section_only_plane_in_a_cap_cuts_the_cap_profile() {
    let solid = gear(None);
    for level in [0., WIDTH] {
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
                frame: PlaneFrame {
                    origin: [0., 0., level],
                    ..Default::default()
                },
                section_offset_mm: Some(0.),
                section_only: true,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(
            vectors
                .iter()
                .all(|vector| matches!(vector.role, VectorRole::Section | VectorRole::Hatch))
        );
        let area: f64 = vectors
            .iter()
            .filter_map(|vector| match vector.curve {
                VectorCurve::HatchTriangle { points: [a, b, c] } => {
                    Some(((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])).abs() / 2.)
                }
                _ => None,
            })
            .sum();
        let exact = profile().area_mm2();
        assert!(
            (area - exact).abs() < 2e-4 * exact,
            "{level}: {area} {exact}"
        );
    }
}

/// The rim pieces of each edge-on cap publish as one straight line, and a
/// rigidly placed gear seen through the equally placed frame gives the same
/// drawing as the gear at the origin.
#[test]
fn edge_on_caps_are_single_lines_and_tilted_placements_draw_identically() {
    let reach = rings(&profile())[0].iter().map(|p| p[0]).fold(0., f64::max);
    let solid = gear(None);
    let options = DrawingProjection {
        frame: side(),
        hidden_lines: true,
        ..Default::default()
    };
    let plain = draw(&[&solid], options);
    for level in [0., WIDTH] {
        let caps: Vec<_> = plain
            .iter()
            .filter(|vector| vector.role == VectorRole::Visible)
            .filter_map(line)
            .filter(|[a, b]| (a[1] - level).abs() < 1e-9 && (b[1] - level).abs() < 1e-9)
            .collect();
        assert_eq!(caps.len(), 1, "{caps:?}");
        let [a, b] = caps[0];
        assert!((a[0].min(b[0]) + reach).abs() < 1e-5 && (a[0].max(b[0]) - reach).abs() < 1e-5);
    }
    // Hidden cap stretches are covered by the visible line.
    assert!(
        !plain
            .iter()
            .filter(|vector| vector.role == VectorRole::Hidden)
            .filter_map(line)
            .any(|[a, b]| (a[1] - b[1]).abs() < 1e-9)
    );

    // A proper rotation (about (1, 2, 3) by 0.7 rad) and a translation.
    let n = 14_f64.sqrt();
    let k = [1. / n, 2. / n, 3. / n];
    let (sin, cos) = 0.7_f64.sin_cos();
    let rotation: [[f64; 3]; 3] = std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            let cross = [[0., -k[2], k[1]], [k[2], 0., -k[0]], [-k[1], k[0], 0.]][i][j];
            f64::from(u8::from(i == j)) * cos + sin * cross + (1. - cos) * k[i] * k[j]
        })
    });
    let translation_mm = [7., -3., 5.];
    let turn = |v: [f64; 3]| -> [f64; 3] {
        std::array::from_fn(|i| (0..3).map(|j| rotation[i][j] * v[j]).sum())
    };
    let placed = solid
        .transformed(RigidTransform {
            rotation,
            translation_mm,
        })
        .unwrap();
    let s = 0.5_f64.sqrt();
    let t = (1_f64 / 6.).sqrt();
    for (label, frame) in [
        ("side", side()),
        (
            "oblique",
            PlaneFrame {
                origin: [0., 0., 2.],
                u: [s, -s, 0.],
                v: [t, t, -2. * t],
            },
        ),
    ] {
        for section in [None, Some(0.4)] {
            let options = DrawingProjection {
                frame,
                hidden_lines: true,
                section_offset_mm: section,
                ..Default::default()
            };
            let moved = PlaneFrame {
                origin: std::array::from_fn(|i| turn(frame.origin)[i] + translation_mm[i]),
                u: turn(frame.u),
                v: turn(frame.v),
            };
            let expected = draw(&[&solid], options);
            let found = draw(
                &[&placed],
                DrawingProjection {
                    frame: moved,
                    ..options
                },
            );
            // Every preview point of one drawing lies on a curve of the
            // other with the same role, far inside the chord tolerance.
            for (from, to) in [(&found, &expected), (&expected, &found)] {
                for vector in from.iter() {
                    let spans: Vec<[[f64; 2]; 2]> = to
                        .iter()
                        .filter(|other| other.role == vector.role)
                        .flat_map(|other| other.preview_mm.windows(2).map(|w| [w[0], w[1]]))
                        .collect();
                    for p in &vector.preview_mm {
                        let distance = spans
                            .iter()
                            .map(|[a, b]| {
                                let d = [b[0] - a[0], b[1] - a[1]];
                                let along = (((p[0] - a[0]) * d[0] + (p[1] - a[1]) * d[1])
                                    / (d[0] * d[0] + d[1] * d[1]).max(1e-300))
                                .clamp(0., 1.);
                                (p[0] - a[0] - along * d[0]).hypot(p[1] - a[1] - along * d[1])
                            })
                            .fold(f64::INFINITY, f64::min);
                        assert!(
                            distance < 1e-6,
                            "{label} {section:?} {:?}: {distance}",
                            vector.role
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn invalid_requests_and_exhausted_budgets_fail_without_partial_output() {
    let solid = gear(None);
    let request = |options| {
        project_drawing(
            &[DrawingSolid {
                solid: &solid,
                hatch: None,
            }],
            options,
        )
    };
    let view = DrawingProjection {
        frame: side(),
        ..Default::default()
    };
    for invalid in [
        DrawingProjection {
            chord_tolerance_mm: 0.,
            ..view
        },
        DrawingProjection {
            section_depth_mm: Some(1.),
            ..view
        },
        DrawingProjection {
            detail: Some(([0., 0.], -1.)),
            ..view
        },
    ] {
        assert_eq!(request(invalid).unwrap_err().kind, GeometryErrorKind::Input);
    }
    for limited in [
        DrawingProjection {
            maximum_curves: 20,
            ..view
        },
        DrawingProjection {
            maximum_points: 40,
            ..view
        },
    ] {
        assert_eq!(request(limited).unwrap_err().kind, GeometryErrorKind::Limit);
    }
}

/// Cold release-build timings of representative scene drawings.
#[test]
#[ignore = "timing report; run with --release --ignored --nocapture"]
fn report_scene_drawing_timings() {
    let solid = gear(None);
    let s = 0.5_f64.sqrt();
    let t = (1_f64 / 6.).sqrt();
    for (label, options) in [
        (
            "side view",
            DrawingProjection {
                frame: side(),
                ..Default::default()
            },
        ),
        (
            "oblique view with hidden lines",
            DrawingProjection {
                frame: PlaneFrame {
                    origin: [0.; 3],
                    u: [s, -s, 0.],
                    v: [t, t, -2. * t],
                },
                hidden_lines: true,
                ..Default::default()
            },
        ),
        (
            "windowed side section",
            DrawingProjection {
                frame: PlaneFrame {
                    origin: [0., 0., 2.],
                    ..side()
                },
                section_offset_mm: Some(0.),
                section_u_mm: Some([4., 12.5]),
                section_v_mm: Some([-1., 1.5]),
                ..Default::default()
            },
        ),
    ] {
        let start = std::time::Instant::now();
        let vectors = draw(&[&solid], options);
        println!(
            "{{\"case\":\"24-tooth spur {label}\",\"milliseconds\":{:.3},\"vectors\":{}}}",
            start.elapsed().as_secs_f64() * 1e3,
            vectors.len()
        );
    }
}
