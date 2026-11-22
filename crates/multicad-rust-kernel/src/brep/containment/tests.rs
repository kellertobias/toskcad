use super::*;
use crate::step_read::brep::read_solids;

const TOLERANCE: f64 = 1e-7;

macro_rules! fixture {
    ($name:literal) => {
        include_str!(concat!(
            "../../../../multicad-app/tests/fixtures/exchange/",
            $name
        ))
    };
}

fn solid(text: &str) -> Solid {
    let mut all = read_solids(text, TOLERANCE, None).unwrap();
    assert_eq!(all.len(), 1);
    all.remove(0)
}

fn cube(center: Point, side: f64) -> [Point; 2] {
    [center.map(|v| v - side / 2.), center.map(|v| v + side / 2.)]
}

/// Clearance, then parity at a deflection a quarter of the clearance.
fn side(solid: &Solid, center: Point, size: f64) -> BoxSide {
    let bounds = cube(center, size);
    let Some(clearance) = solid.box_clearance(bounds, TOLERANCE, None).unwrap() else {
        return BoxSide::Unresolved;
    };
    let deflection = (clearance / 4.).min(0.5);
    let tessellation = mesh::tessellate(
        solid,
        mesh::Options {
            deflection_mm: deflection,
            angle_rad: 0.5,
            triangle_limit: 4_000_000,
        },
        None,
    )
    .unwrap();
    solid.box_parity(bounds, clearance, &tessellation, deflection)
}

#[test]
fn analytic_supports_settle_boxes_by_their_exact_distance() {
    let cylinder = solid(fixture!("corpus/occt-cylinder.step"));
    let [lo, hi] = cylinder.bounds_mm();
    let z = 0.5 * (lo[2] + hi[2]);
    let [x, y] = [0.5 * (lo[0] + hi[0]), 0.5 * (lo[1] + hi[1])];
    assert!((hi[2] - lo[2] - 15.).abs() < 1e-9);
    assert_eq!(side(&cylinder, [x, y, z], 1.), BoxSide::Inside);
    assert_eq!(side(&cylinder, [x + 5.9, y, z], 0.01), BoxSide::Inside);
    assert_eq!(side(&cylinder, [x + 6.1, y, z], 0.01), BoxSide::Outside);
    assert_eq!(side(&cylinder, [x, y, hi[2] + 0.2], 0.01), BoxSide::Outside);
    // On the wall's infinite support but beyond the trimmed side face: the
    // face's own axial slab settles it.
    assert_eq!(
        side(&cylinder, [x + 6., y, hi[2] + 0.5], 0.01),
        BoxSide::Outside
    );
    // A box across the wall or a cap is never settled.
    assert_eq!(side(&cylinder, [x + 6., y, z], 0.01), BoxSide::Unresolved);
    assert_eq!(side(&cylinder, [x, y, hi[2]], 0.01), BoxSide::Unresolved);
    // The lower bound is the true distance up to the box reach.
    let gap = cylinder
        .box_clearance(cube([x + 6.5, y, z], 0.), TOLERANCE, None)
        .unwrap()
        .unwrap();
    assert!((gap - 0.5).abs() < 1e-9, "{gap}");
}

#[test]
fn spline_extrusions_settle_by_halved_bezier_hulls() {
    // The cubic from (20, 0) through controls (20, 12) and (0, 12) to the
    // origin, closed by the X axis and swept 6 mm along Z; its top at
    // t = 1/2 is (10, 9).
    let prism = solid(fixture!("corpus/occt-spline-prism.step"));
    assert_eq!(side(&prism, [10., 4.5, 3.], 1.), BoxSide::Inside);
    assert_eq!(side(&prism, [10., 8.95, 3.], 0.01), BoxSide::Inside);
    assert_eq!(side(&prism, [10., 9.05, 3.], 0.01), BoxSide::Outside);
    assert_eq!(side(&prism, [10., 9., 3.], 0.01), BoxSide::Unresolved);
    assert_eq!(side(&prism, [10., -0.05, 3.], 0.01), BoxSide::Outside);
}

#[test]
fn voids_and_invalid_boxes() {
    let hollow = solid(fixture!("corpus/occt-box-void.step"));
    let boundary = hollow.boundary();
    let void_vertices = boundary.shells[1]
        .faces
        .iter()
        .flat_map(|&f| &boundary.faces[f].loops)
        .flat_map(|l| &l.uses)
        .flat_map(|&(e, _)| boundary.edges[e].vertices)
        .map(|v| boundary.vertices[v]);
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for p in void_vertices {
        for i in 0..3 {
            lo[i] = lo[i].min(p[i]);
            hi[i] = hi[i].max(p[i]);
        }
    }
    let middle = std::array::from_fn(|i| 0.5 * (lo[i] + hi[i]));
    assert_eq!(side(&hollow, middle, 0.1), BoxSide::Outside);
    // Just below the void floor lies in the wall material.
    let wall = [middle[0], middle[1], lo[2] - 0.05];
    assert_eq!(side(&hollow, wall, 0.01), BoxSide::Inside);
    assert_eq!(
        hollow.box_clearance([[1., 0., 0.], [0., 0., 0.]], TOLERANCE, None),
        Err(super::super::Error::Invalid)
    );
    let stop = || true;
    assert_eq!(
        hollow.box_clearance(cube(middle, 0.1), TOLERANCE, Some(&stop)),
        Err(super::super::Error::Cancelled)
    );
}

#[test]
fn a_box_near_the_boundary_escapes_within_clearance_balls_to_coarse_parity() {
    let prism = solid(fixture!("corpus/occt-spline-prism.step"));
    let coarse = 0.25;
    let tessellation = mesh::tessellate(
        &prism,
        mesh::Options {
            deflection_mm: coarse,
            angle_rad: 0.5,
            triangle_limit: 100_000,
        },
        None,
    )
    .unwrap();
    for (center, expected) in [
        ([10., 8.95, 3.], BoxSide::Inside),
        ([10., 9.05, 3.], BoxSide::Outside),
        ([10., 0.04, 3.], BoxSide::Inside),
    ] {
        let bounds = cube(center, 0.01);
        let clearance = prism
            .box_clearance(bounds, TOLERANCE, None)
            .unwrap()
            .unwrap();
        // Too close for the coarse tessellation at the box itself.
        assert!(clearance < 4. * coarse);
        assert_eq!(
            prism.box_parity(bounds, clearance, &tessellation, coarse),
            BoxSide::Unresolved
        );
        let (point, free) = prism
            .escape(center, clearance, 4. * coarse, TOLERANCE, None)
            .unwrap();
        assert!(free >= 4. * coarse);
        assert_eq!(
            prism.box_parity([point, point], free, &tessellation, coarse),
            expected,
            "{center:?} -> {point:?}"
        );
    }
    // A point boxed in by a thin wall never reaches the wanted clearance.
    let bounds = cube([10., 4.5, 3.], 0.01);
    let clearance = prism
        .box_clearance(bounds, TOLERANCE, None)
        .unwrap()
        .unwrap();
    let (_, free) = prism
        .escape([10., 4.5, 3.], clearance, 100., TOLERANCE, None)
        .unwrap();
    assert!(free >= clearance && free < 4.5, "{free}");
}
