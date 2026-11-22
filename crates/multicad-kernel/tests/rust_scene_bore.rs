//! Analytic regressions for plain views and general retained sections of a
//! box with an oblique bore (CAD-354) drawn through the scene engine
//! (CAD-346): a 10 x 8 x 6 mm stock bored with radius 1 mm through its
//! center along (1, 0.5, 1).
use multicad_core::{
    document::EntityId,
    drawings::{VectorCurve, VectorRole},
    materials::HatchPattern,
    references::Origin,
};
use multicad_kernel::{
    BooleanOperation, DrawingHatch, DrawingProjection, DrawingSolid, GeometryPolicy, PlaneFrame,
    ProjectedVector, Solid, project_drawing,
    rust_kernel::{Circle, Rational, cylinder::AxialCylinder},
};
use std::f64::consts::PI;

type V = [f64; 3];
const CENTER: V = [5., 4., 3.];

fn cross(a: V, b: V) -> V {
    std::array::from_fn(|i| a[(i + 1) % 3] * b[(i + 2) % 3] - a[(i + 2) % 3] * b[(i + 1) % 3])
}

fn normalized(v: V) -> V {
    let length = v.iter().map(|x| x * x).sum::<f64>().sqrt();
    v.map(|x| x / length)
}

fn axis() -> V {
    normalized([1., 0.5, 1.])
}

fn frame(origin: V, normal: V) -> PlaneFrame {
    let w = normalized(normal);
    let seed = if w[2].abs() < 0.9 {
        [0., 0., 1.]
    } else {
        [1., 0., 0.]
    };
    let u = normalized(cross(seed, w));
    PlaneFrame {
        origin,
        u,
        v: cross(w, u),
    }
}

fn bored() -> Solid {
    let policy = GeometryPolicy::default();
    let rational = |value: i64| Rational::new(value, 1).unwrap();
    let stock = Solid::rust_box(
        EntityId::new().uuid_bits(),
        [10, 8, 6].map(rational),
        PlaneFrame::default(),
        policy,
    )
    .unwrap();
    let id = EntityId::new();
    let base = std::array::from_fn(|i| CENTER[i] - 10. * axis()[i]);
    let tool = Solid::from_axial_cylinder_sweep(
        AxialCylinder::new(
            id.uuid_bits(),
            Circle::new([rational(0); 2], rational(1)).unwrap(),
            rational(0),
            rational(20),
        )
        .unwrap(),
        frame(base, axis()),
        [0.; 2],
        false,
        Origin {
            source: id,
            role: "primitive.circle".into(),
        },
        policy,
    )
    .unwrap();
    stock
        .boolean(&tool, BooleanOperation::Subtract, policy)
        .unwrap()
        .unwrap()
}

fn draw(options: DrawingProjection, pattern: HatchPattern) -> Vec<ProjectedVector> {
    project_drawing(
        &[DrawingSolid {
            solid: &bored(),
            hatch: Some(DrawingHatch {
                angle_rad: 0.3,
                spacing_mm: 0.4,
                line_width_mm: 0.18,
                pattern,
            }),
        }],
        options,
    )
    .unwrap()
}

fn ellipses(vectors: &[ProjectedVector], role: VectorRole) -> Vec<([f64; 2], f64, f64)> {
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
}

/// Seen along the bore, both mouths project onto the bore circle and the
/// wall has no outline.
#[test]
fn the_view_along_the_bore_shows_one_exact_circle() {
    let vectors = draw(
        DrawingProjection {
            frame: frame(CENTER, axis()),
            hidden_lines: true,
            ..Default::default()
        },
        HatchPattern::None,
    );
    let rims = ellipses(&vectors, VectorRole::Visible);
    assert!(!rims.is_empty());
    for (center, major, minor) in rims.iter().chain(&ellipses(&vectors, VectorRole::Hidden)) {
        assert!(center[0].hypot(center[1]) < 1e-9, "{center:?}");
        assert!((major - 1.).abs() < 1e-9 && (minor - 1.).abs() < 1e-9);
    }
    // The visible arcs close the circle.
    let sweep: f64 = vectors
        .iter()
        .filter(|vector| vector.role == VectorRole::Visible)
        .filter_map(|vector| match vector.curve {
            VectorCurve::Ellipse {
                start_rad, end_rad, ..
            } => Some(end_rad - start_rad),
            _ => None,
        })
        .sum();
    assert!((sweep - 2. * PI).abs() < 1e-9, "{sweep}");
}

/// From above, the mouth on the top face is the ellipse with semi-axes
/// radius and radius / cos(axis, Z); the wall's two outline generators run
/// hidden inside the stock, parallel to the projected axis at one radius.
#[test]
fn the_top_view_shows_the_mouth_ellipse_and_hidden_wall_outlines() {
    let vectors = draw(
        DrawingProjection {
            hidden_lines: true,
            ..Default::default()
        },
        HatchPattern::None,
    );
    let major = 1. / axis()[2];
    let visible = ellipses(&vectors, VectorRole::Visible);
    assert!(!visible.is_empty());
    for (_, a, b) in visible
        .iter()
        .chain(&ellipses(&vectors, VectorRole::Hidden))
    {
        assert!((a - major).abs() < 1e-9 && (b - 1.).abs() < 1e-9, "{a} {b}");
    }
    let along = {
        let length = axis()[0].hypot(axis()[1]);
        [axis()[0] / length, axis()[1] / length]
    };
    let mut offsets: Vec<f64> = vectors
        .iter()
        .filter(|vector| vector.role == VectorRole::Hidden)
        .filter_map(|vector| match vector.curve {
            VectorCurve::Line { a, b } => Some([a, b]),
            _ => None,
        })
        .filter(|[a, b]| {
            let d = [b[0] - a[0], b[1] - a[1]];
            (d[0] * along[1] - d[1] * along[0]).abs() < 1e-9 * d[0].hypot(d[1])
        })
        .map(|[a, _]| (a[0] - CENTER[0]) * along[1] - (a[1] - CENTER[1]) * along[0])
        .collect();
    offsets.sort_by(f64::total_cmp);
    assert_eq!(offsets.len(), 2, "{offsets:?}");
    assert!((offsets[0] + 1.).abs() < 1e-9 && (offsets[1] - 1.).abs() < 1e-9);
    // Without hidden lines neither the outlines nor any hidden rim appear.
    let plain = draw(DrawingProjection::default(), HatchPattern::None);
    assert!(
        plain
            .iter()
            .all(|vector| vector.role == VectorRole::Visible)
    );
}

/// A retained horizontal cut is not along the bore: its face is the stock
/// rectangle minus the bore ellipse, hatched only in material.
#[test]
fn a_retained_cut_across_the_bore_publishes_the_exact_face_and_hatch() {
    let options = DrawingProjection {
        frame: PlaneFrame {
            origin: [0., 0., 3.],
            ..Default::default()
        },
        section_offset_mm: Some(0.),
        hidden_lines: true,
        ..Default::default()
    };
    let major = 1. / axis()[2];
    let vectors = draw(options, HatchPattern::Solid);
    let cut = ellipses(&vectors, VectorRole::Section);
    assert!(!cut.is_empty());
    for (center, a, b) in cut {
        assert!((center[0] - CENTER[0]).hypot(center[1] - CENTER[1]) < 1e-9);
        assert!((a - major).abs() < 1e-9 && (b - 1.).abs() < 1e-9);
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
    let exact = 80. - PI * major;
    assert!((area - exact).abs() < 1e-3, "{area} {exact}");
    // The cut face boundary is also the visible outline of the retained part,
    // and the lower mouth behind it is hidden.
    assert!(!ellipses(&vectors, VectorRole::Visible).is_empty());
    assert!(!ellipses(&vectors, VectorRole::Hidden).is_empty());
}
