//! Oblique cylindrical cuts in box stock exchange through Rust STEP.
use multicad_core::document::EntityId;
use multicad_kernel::{BooleanOperation, GeometryPolicy, PlaneFrame, Solid};

type V = [f64; 3];

fn policy() -> GeometryPolicy {
    GeometryPolicy::default()
}

fn cross(a: V, b: V) -> V {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn normalized(v: V) -> V {
    let length = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    v.map(|x| x / length)
}

fn frame(origin: V, travel: V) -> PlaneFrame {
    let w = normalized(travel);
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

fn cut(stock: &Solid, tool: PlaneFrame, radius: f64, height: f64) -> Solid {
    let tool = Solid::cylinder(tool, radius, height, policy()).unwrap();
    stock
        .boolean(&tool, BooleanOperation::Subtract, policy())
        .unwrap()
        .unwrap()
}

fn round_trip(name: &str, solid: &Solid) {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("cut.step");
    assert!(solid.supports_step_export(), "{name}");
    solid
        .write_step(&path)
        .unwrap_or_else(|e| panic!("{name}: {e}"));
    let text = std::fs::read_to_string(&path).unwrap();
    let back = Solid::import_rust_step(&text, EntityId::new(), policy())
        .unwrap_or_else(|e| panic!("{name}: {e}"));
    let [a, b] = [solid, &back].map(|s| s.volume_mm3().unwrap());
    // One closed exterior shell on the retained analytic supports.
    assert_eq!(text.matches("=MANIFOLD_SOLID_BREP(").count(), 1, "{name}");
    assert_eq!(text.matches("=CLOSED_SHELL(").count(), 1, "{name}");
    assert!(!text.contains("B_SPLINE"), "{name}");
    assert_eq!(
        text.matches("=CYLINDRICAL_SURFACE(").count(),
        solid
            .faces()
            .unwrap()
            .iter()
            .filter(|f| f.plane.is_none())
            .count(),
        "{name}"
    );
    assert!(text.contains("=ELLIPSE("), "{name}");
    assert_eq!(back.solid_count().unwrap(), 1, "{name}");
    assert_eq!(
        back.faces().unwrap().len(),
        solid.faces().unwrap().len(),
        "{name}"
    );
    assert!((a - b).abs() <= 1e-9 * a, "{name}: {a} {b}");
    let [p, q] = [solid, &back].map(|s| s.bounds().unwrap());
    for i in 0..3 {
        assert!(
            (p.min[i] - q.min[i]).abs() < 1e-9 && (p.max[i] - q.max[i]).abs() < 1e-9,
            "{name}"
        );
    }
}

#[test]
fn oblique_cuts_round_trip_with_exact_volume_and_bounds() {
    let stock = Solid::box_solid([10., 8., 6.], policy()).unwrap();
    let mcp = PlaneFrame {
        origin: [4., 4., 3.],
        u: [0.8, 0., -0.6],
        ..Default::default()
    };
    round_trip("mcp blind", &cut(&stock, mcp, 2., 5.));
    round_trip(
        "through",
        &cut(&stock, frame([3., 4., -2.], [0.2, 0.1, 1.]), 1., 20.),
    );
    round_trip(
        "blind",
        &cut(&stock, frame([5., 4., 8.], [0.2, -0.1, -1.]), 1., 5.),
    );
    round_trip(
        "slot",
        &cut(&stock, frame([-2., -1., 5.5], [1., 1., 0.]), 1., 30.),
    );
    round_trip(
        "notch",
        &cut(&stock, frame([-2., 6., 3.], [1., 1., 0.]), 1., 30.),
    );
    round_trip(
        "corner",
        &cut(&stock, frame([9.2, 4., 3.], [0.3, 0., 1.]), 1., 6.),
    );
    round_trip(
        "edge",
        &cut(&stock, frame([9.5, 7.5, 3.], [0.3, 0.2, 1.]), 1., 6.),
    );
    round_trip(
        "through corner",
        &cut(&stock, frame([9., 1., -1.], [0.2, -0.2, 1.]), 1.5, 9.),
    );
}
