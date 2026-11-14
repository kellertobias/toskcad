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
