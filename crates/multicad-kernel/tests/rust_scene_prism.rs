//! Analytic regressions for cylinders, annuli and rounded slots drawn through
//! the exact-curve scene engine (CAD-346): the retained oblique sections and
//! windowed side sections that the family providers refuse.
use multicad_core::{
    document::EntityId,
    drawings::{VectorCurve, VectorRole},
    materials::HatchPattern,
    references::Origin,
};
use multicad_kernel::{
    DrawingHatch, DrawingProjection, DrawingSolid, GeometryPolicy, PlaneFrame, ProjectedVector,
    Solid, project_drawing,
    rust_kernel::{Circle, Rational, cylinder::AxialCylinder},
};
use std::f64::consts::PI;

const RADIUS: f64 = 5.;
const HEIGHT: f64 = 8.;

fn cylinder() -> Solid {
    let rational = |value: i64| Rational::new(value, 1).unwrap();
    Solid::from_axial_cylinder(
        AxialCylinder::new(
            EntityId::new().uuid_bits(),
            Circle::new([rational(0); 2], rational(5)).unwrap(),
            rational(0),
            rational(8),
        )
        .unwrap(),
        PlaneFrame::default(),
        GeometryPolicy::default(),
    )
    .unwrap()
}

fn draw(solid: &Solid, options: DrawingProjection, pattern: HatchPattern) -> Vec<ProjectedVector> {
    project_drawing(
        &[DrawingSolid {
            solid,
            hatch: Some(DrawingHatch {
                angle_rad: 0.6,
                spacing_mm: 0.5,
                line_width_mm: 0.18,
                pattern,
            }),
        }],
        options,
    )
    .unwrap()
}

/// A frame through `origin` whose normal is tilted by `tilt` from the
/// cylinder axis toward -Y; paper u is world X.
fn tilted(origin: [f64; 3], tilt: f64) -> PlaneFrame {
    let (sin, cos) = tilt.sin_cos();
    PlaneFrame {
        origin,
        u: [1., 0., 0.],
        v: [0., cos, sin],
    }
}

fn triangles(vectors: &[ProjectedVector]) -> f64 {
    vectors
        .iter()
        .filter_map(|vector| match vector.curve {
            VectorCurve::HatchTriangle { points: [a, b, c] } => {
                Some(((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])).abs() / 2.)
            }
            _ => None,
        })
        .sum()
}

/// A retained cut by a plane tilted from the axis: the cut face is the wall
/// ellipse with semi-axes radius and radius / cos(tilt); the rim behind the
/// plane is the cap circle foreshortened to radius * cos(tilt), visible
/// outside the cut face and hidden behind it; the outlines are the two
/// generators at one radius from the axis.
#[test]
fn a_retained_oblique_cylinder_section_has_exact_conics_and_generators() {
    let tilt = 0.35_f64;
    let (sin, cos) = tilt.sin_cos();
    let solid = cylinder();
    let frame = tilted([0., 0., HEIGHT / 2.], tilt);
    let options = DrawingProjection {
        frame,
        section_offset_mm: Some(0.),
        hidden_lines: true,
        ..Default::default()
    };
    let vectors = draw(&solid, options, HatchPattern::Solid);
    let conics = |role: VectorRole| -> Vec<([f64; 2], f64, f64)> {
        vectors
            .iter()
            .filter(|vector| vector.role == role)
            .filter_map(|vector| match vector.curve {
                VectorCurve::Ellipse {
                    center,
                    major_mm,
                    minor_mm,
                    ..
                } => Some((center, major_mm, minor_mm)),
                _ => None,
            })
            .collect()
    };
    let face = conics(VectorRole::Section);
    assert!(!face.is_empty());
    for (center, major, minor) in face {
        assert!(center[0].hypot(center[1]) < 1e-9);
        assert!((major - RADIUS / cos).abs() < 1e-9 && (minor - RADIUS).abs() < 1e-9);
    }
    // The whole wall ellipse lies between the caps here, so the cut face is
    // the full ellipse and its area is exact.
    assert!(RADIUS * sin / cos < HEIGHT / 2.);
    let exact = PI * RADIUS * RADIUS / cos;
    let area = triangles(&vectors);
    assert!((area - exact).abs() < 2e-4 * exact, "{area} {exact}");
    // The retained cap is the one behind the plane, seen foreshortened.
    let rear: Vec<_> = conics(VectorRole::Visible)
        .into_iter()
        .chain(conics(VectorRole::Hidden))
        .filter(|(center, ..)| center[1].abs() > 1e-6)
        .collect();
    assert!(!rear.is_empty());
    for (center, major, minor) in rear {
        assert!(center[0].abs() < 1e-9 && (center[1].abs() - HEIGHT / 2. * sin).abs() < 1e-9);
        assert!((major - RADIUS).abs() < 1e-9 && (minor - RADIUS * cos).abs() < 1e-9);
    }
    let mut outlines: Vec<[[f64; 2]; 2]> = vectors
        .iter()
        .filter(|vector| vector.role == VectorRole::Visible)
        .filter_map(|vector| match vector.curve {
            VectorCurve::Line { a, b } => Some([a, b]),
            _ => None,
        })
        .collect();
    outlines.sort_by(|a, b| a[0][0].total_cmp(&b[0][0]));
    assert_eq!(outlines.len(), 2);
    for ([a, b], side) in outlines.into_iter().zip([-RADIUS, RADIUS]) {
        assert!((a[0] - side).abs() < 1e-9 && (b[0] - side).abs() < 1e-9);
        // From the cut plane, where the outline meets the face, to the cap.
        assert!((a[1].abs().max(b[1].abs()) - HEIGHT / 2. * sin).abs() < 1e-9);
        assert!(a[1].abs().min(b[1].abs()) < 1e-9);
    }
}

/// Section only of a rounded slot by an oblique plane: the round ends cut in
/// elliptical arcs with semi-axes radius and radius / cos(tilt).
#[test]
fn an_oblique_section_of_a_rounded_slot_keeps_exact_elliptical_ends() {
    let owner = EntityId::new();
    let slot = Solid::rust_extrude_rounded_slot(
        PlaneFrame::default(),
        [[0., 0.], [6., 0.]],
        2.,
        5.,
        owner,
        std::array::from_fn(|edge| Origin {
            source: owner,
            role: format!("profile.{edge}"),
        }),
        GeometryPolicy::default(),
    )
    .unwrap();
    let tilt = 0.3_f64;
    let vectors = draw(
        &slot,
        DrawingProjection {
            frame: tilted([3., 0., 2.5], tilt),
            section_offset_mm: Some(0.),
            section_only: true,
            ..Default::default()
        },
        HatchPattern::Solid,
    );
    assert!(
        vectors
            .iter()
            .all(|vector| matches!(vector.role, VectorRole::Section | VectorRole::Hatch))
    );
    let mut ends = 0;
    for vector in &vectors {
        if let VectorCurve::Ellipse {
            major_mm, minor_mm, ..
        } = vector.curve
        {
            assert!((major_mm - 2. / tilt.cos()).abs() < 1e-9 && (minor_mm - 2.).abs() < 1e-9);
            ends += 1;
        }
    }
    assert!(ends >= 2);
    // Stadium area: a 6 x 4 rectangle and a circle of radius 2, over cos.
    let exact = (24. + 4. * PI) / tilt.cos();
    let area = triangles(&vectors);
    assert!((area - exact).abs() < 2e-4 * exact, "{area} {exact}");
}

/// A windowed lengthwise section of an annulus: the cut walls inside the
/// window, and the window's lateral faces drawn where material lies behind.
#[test]
fn a_windowed_annulus_section_draws_cut_walls_and_window_faces() {
    let ring = Solid::rust_annulus(
        3.,
        6.,
        2.5,
        EntityId::new(),
        PlaneFrame::default(),
        GeometryPolicy::default(),
    )
    .unwrap();
    let window = [[-4.5, 5.], [-1., 2.]];
    let vectors = draw(
        &ring,
        DrawingProjection {
            frame: PlaneFrame {
                origin: [0., 0., 1.25],
                u: [1., 0., 0.],
                v: [0., 0., 1.],
            },
            section_offset_mm: Some(0.),
            section_u_mm: Some(window[0]),
            section_v_mm: Some(window[1]),
            ..Default::default()
        },
        HatchPattern::Solid,
    );
    for p in vectors.iter().flat_map(|vector| &vector.preview_mm) {
        for k in 0..2 {
            assert!(p[k] >= window[k][0] - 1e-9 && p[k] <= window[k][1] + 1e-9);
        }
    }
    // Cut walls: x in [-4.5, -3] and [3, 5], z from the window's lower edge
    // to the top cap at v = 1.25.
    let exact = (1.5 + 2.) * 2.25;
    let area = triangles(&vectors);
    assert!((area - exact).abs() < 1e-9, "{area} {exact}");
    // The right window border lies in the wall: visible over the wall height.
    let lines: Vec<[[f64; 2]; 2]> = vectors
        .iter()
        .filter(|vector| vector.role == VectorRole::Visible)
        .filter_map(|vector| match vector.curve {
            VectorCurve::Line { a, b } => Some([a, b]),
            _ => None,
        })
        .collect();
    assert!(lines.iter().any(|[a, b]| (a[0] - 5.).abs() < 1e-9
        && (b[0] - 5.).abs() < 1e-9
        && (a[1].min(b[1]) + 1.).abs() < 1e-9
        && (a[1].max(b[1]) - 1.25).abs() < 1e-9));
    // Through the bore, the far wall's lower window face shows as the bottom
    // border between the two cut walls.
    assert!(lines.iter().any(|[a, b]| (a[1] + 1.).abs() < 1e-9
        && (b[1] + 1.).abs() < 1e-9
        && a[0].min(b[0]) <= -3. + 1e-9
        && a[0].max(b[0]) >= 3. - 1e-9));
}

fn arcs(vectors: &[ProjectedVector], role: VectorRole, source: usize) -> Vec<(f64, f64)> {
    vectors
        .iter()
        .filter(|vector| vector.role == role && vector.source_index == source)
        .filter_map(|vector| match vector.curve {
            VectorCurve::Ellipse {
                major_mm,
                start_rad,
                end_rad,
                ..
            } => Some((major_mm, end_rad - start_rad)),
            _ => None,
        })
        .collect()
}

/// A section-only plane lying in a cap meets the solid in that cap from
/// either side, as the native common of the solid and the plane does. A
/// retained section whose plane lies in the cap toward the viewer keeps the
/// whole solid uncut, and one in the far cap retains nothing.
#[test]
fn planes_coincident_with_a_cap_cut_the_cap_or_nothing() {
    let solid = cylinder();
    for level in [0., HEIGHT] {
        let frame = PlaneFrame {
            origin: [0., 0., level],
            ..Default::default()
        };
        let only = draw(
            &solid,
            DrawingProjection {
                frame,
                section_offset_mm: Some(0.),
                section_only: true,
                ..Default::default()
            },
            HatchPattern::Solid,
        );
        let rim = arcs(&only, VectorRole::Section, 0);
        assert!(
            rim.iter()
                .all(|(radius, _)| (radius - RADIUS).abs() < 1e-12)
        );
        let sweep: f64 = rim.iter().map(|(_, sweep)| sweep).sum();
        assert!((sweep - 2. * PI).abs() < 1e-9, "{level}: {sweep}");
        // The hatch region is an inscribed polygon at the chord tolerance.
        let area = triangles(&only);
        let exact = PI * RADIUS * RADIUS;
        assert!(area <= exact + 1e-9 && area > exact * (1. - 2e-3), "{area}");
        let retained = draw(
            &solid,
            DrawingProjection {
                frame,
                section_offset_mm: Some(0.),
                ..Default::default()
            },
            HatchPattern::Solid,
        );
        if level == 0. {
            assert!(retained.is_empty());
        } else {
            // Nothing is cut off: the cap is an ordinary visible face.
            assert!(
                retained
                    .iter()
                    .all(|vector| vector.role == VectorRole::Visible)
            );
            let sweep: f64 = arcs(&retained, VectorRole::Visible, 0)
                .iter()
                .map(|(_, sweep)| sweep)
                .sum();
            assert!((sweep - 2. * PI).abs() < 1e-9);
        }
    }
}

/// A window border 0.1 um inside the rim: the rim leaves and re-enters the
/// window within one certified chord. Both crossings are found on the exact
/// circle, so the cut face keeps the exact arc and gains the exact border
/// chord.
#[test]
fn a_rim_grazing_the_window_is_clipped_at_both_exact_crossings() {
    let border = RADIUS - 1e-4;
    let vectors = draw(
        &cylinder(),
        DrawingProjection {
            frame: PlaneFrame {
                origin: [0., 0., HEIGHT / 2.],
                ..Default::default()
            },
            section_offset_mm: Some(0.),
            section_only: true,
            section_u_mm: Some([-10., border]),
            ..Default::default()
        },
        HatchPattern::None,
    );
    let half = (border / RADIUS).acos();
    assert!(half < 0.01);
    let sweep: f64 = arcs(&vectors, VectorRole::Section, 0)
        .iter()
        .map(|(_, sweep)| sweep)
        .sum();
    assert!((sweep - (2. * PI - 2. * half)).abs() < 1e-9, "{sweep}");
    let chords: Vec<[[f64; 2]; 2]> = vectors
        .iter()
        .filter_map(|vector| match vector.curve {
            VectorCurve::Line { a, b } => Some([a, b]),
            _ => None,
        })
        .collect();
    assert_eq!(chords.len(), 1);
    let [a, b] = chords[0];
    let reach = (RADIUS * RADIUS - border * border).sqrt();
    assert!((a[0] - border).abs() < 1e-12 && (b[0] - border).abs() < 1e-12);
    assert!((a[1].abs() - reach).abs() < 1e-9 && (a[1] + b[1]).abs() < 1e-9);
}

/// An annulus above a cylinder whose rim is internally tangent to its bore:
/// the rim shows only outside the annulus, over the exact arc, and the single
/// point seen through the bore publishes no degenerate vector.
#[test]
fn a_rim_tangent_to_an_occluding_bore_publishes_no_degenerate_arc() {
    let ring = Solid::rust_annulus(
        3.,
        6.,
        2.5,
        EntityId::new(),
        PlaneFrame {
            origin: [2., 0., 20.],
            ..Default::default()
        },
        GeometryPolicy::default(),
    )
    .unwrap();
    let solid = cylinder();
    for hidden_lines in [false, true] {
        let vectors = project_drawing(
            &[&solid, &ring].map(|solid| DrawingSolid { solid, hatch: None }),
            DrawingProjection {
                hidden_lines,
                ..Default::default()
            },
        )
        .unwrap();
        for vector in &vectors {
            if let VectorCurve::Ellipse {
                major_mm,
                start_rad,
                end_rad,
                ..
            } = vector.curve
            {
                assert!(major_mm * (end_rad - start_rad) > 1e-5, "{vector:?}");
            }
        }
        // |p - (2, 0)| > 6 on the rim: cos(angle) < -7 / 20.
        let visible = arcs(&vectors, VectorRole::Visible, 0);
        assert_eq!(visible.len(), 1);
        let exact = 2. * PI - 2. * (-0.35_f64).acos();
        assert!((visible[0].1 - exact).abs() < 1e-9, "{visible:?}");
    }
}

/// A box edge piercing a spur gear's body: the edge turns hidden exactly
/// where it enters the root cylinder and visible again where it leaves it,
/// although nothing crosses it in the paper at those points.
#[test]
fn an_edge_piercing_another_body_changes_visibility_at_the_exact_entry() {
    use multicad_kernel::{GearSweep, InvoluteGear, RigidTransform};
    let a = 20_f64.to_radians();
    let pitch = 12.;
    let gear = Solid::rust_spur_gear_named(
        InvoluteGear {
            hub: None,
            teeth: 24,
            pitch_radius_mm: pitch,
            base_radius_mm: pitch * a.cos(),
            tip_radius_mm: pitch + 1.,
            root_radius_mm: pitch - 1.25,
            root_fillet_radius_mm: 0.2,
            half_base_angle_rad: (PI / 2. - 0.05) / (2. * pitch) + a.tan() - a,
            face_width_mm: 4.,
            bore_radius_mm: 3.,
            twist_rad_per_mm: 0.,
            sweep: GearSweep::Spur,
        },
        EntityId::new(),
        GeometryPolicy::default(),
    )
    .unwrap();
    // x in [-4, 4], y in [-6, 6], z in [1.5, 9.5]: the lower front edge runs
    // along x at y = -6, z = 1.5, inside the root cylinder (radius 10.75) and
    // outside the bore.
    let rational = |value: i64| Rational::new(value, 1).unwrap();
    let block = Solid::rust_box(
        EntityId::new().uuid_bits(),
        [8, 12, 8].map(rational),
        PlaneFrame::default(),
        GeometryPolicy::default(),
    )
    .unwrap()
    .transformed(RigidTransform {
        rotation: [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
        translation_mm: [-4., -6., 1.5],
    })
    .unwrap();
    // Seen from (1, 1, 1): the vertical block edges at (4, -6) and (-4, 6)
    // are silhouettes of the block. Each runs from z = 1.5 inside the gear
    // body (radius 7.2 between bore and root) through the top cap at z = 4,
    // far from any rim in the paper, and is visible above it.
    let s = 0.5_f64.sqrt();
    let t = (1_f64 / 6.).sqrt();
    let [u, v] = [[s, -s, 0.], [t, t, -2. * t]];
    let vectors = project_drawing(
        &[&gear, &block].map(|solid| DrawingSolid { solid, hatch: None }),
        DrawingProjection {
            frame: PlaneFrame {
                origin: [0.; 3],
                u,
                v,
            },
            hidden_lines: true,
            ..Default::default()
        },
    )
    .unwrap();
    let paper = |p: [f64; 3]| -> [f64; 2] {
        [u, v].map(|axis| axis[0] * p[0] + axis[1] * p[1] + axis[2] * p[2])
    };
    for [x, y] in [[4., -6.], [-4., 6.]] {
        for (role, [low, high]) in [
            (VectorRole::Hidden, [1.5, 4.]),
            (VectorRole::Visible, [4., 9.5]),
        ] {
            let [from, to] = [paper([x, y, low]), paper([x, y, high])];
            let near = |p: [f64; 2], q: [f64; 2]| (p[0] - q[0]).hypot(p[1] - q[1]) < 1e-9;
            let found = vectors
                .iter()
                .filter(|vector| vector.source_index == 1 && vector.role == role)
                .any(|vector| match vector.curve {
                    VectorCurve::Line { a, b } => {
                        (near(a, from) && near(b, to)) || (near(a, to) && near(b, from))
                    }
                    _ => false,
                });
            assert!(
                found,
                "{role:?} edge at ({x}, {y}) from z = {low} to {high}"
            );
        }
    }
}
