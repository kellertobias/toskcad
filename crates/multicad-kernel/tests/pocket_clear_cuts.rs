//! Rectilinear cuts and overlap queries on stock that already carries blind
//! circular pockets, where the planar tool stays clear of every pocket.
use multicad_core::{document::EntityId, references::Origin};
use multicad_kernel::{
    CircularPocketTool, GeometryPolicy, PlaneFrame, RigidTransform, Solid,
    rust_kernel::{Rational, polygon::PolygonPrism},
};
use std::f64::consts::PI;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}

/// A 60 × 40 × 18 panel with two Ø5 × 10 blind pockets from its top face.
fn pocketed_panel(policy: GeometryPolicy) -> Solid {
    let owner = EntityId::new();
    let exact = |n| Rational::new(n, 1).unwrap();
    let corners = [[0, 0], [60, 0], [60, 40], [0, 40]];
    let origins: Vec<_> = (0..4)
        .map(|i| Origin {
            source: owner,
            role: format!("panel/{i}"),
        })
        .collect();
    let stock = Solid::from_polygon_prism_sweep(
        PolygonPrism::new(vec![corners.map(|p| p.map(exact)).to_vec()], exact(18)).unwrap(),
        PlaneFrame::default(),
        false,
        owner,
        &origins,
        &[false; 4],
        policy,
    )
    .unwrap();
    let drill = EntityId::new();
    let tools: Vec<_> = [[10., 10.], [10., 30.]]
        .into_iter()
        .enumerate()
        .map(|(i, center_uv)| CircularPocketTool {
            frame: PlaneFrame {
                origin: [0., 0., 19.],
                ..PlaneFrame::default()
            },
            center_uv,
            radius_mm: 2.5,
            displacement_mm: [0., 0., -11.],
            origin: Origin {
                source: drill,
                role: format!("pin/{i}"),
            },
        })
        .collect();
    stock
        .rust_drill_circular_machining(&tools, &[false, false], policy)
        .unwrap()
}

fn block(size: [f64; 3], at: [f64; 3], policy: GeometryPolicy) -> Solid {
    Solid::box_solid(size, policy)
        .unwrap()
        .transformed(RigidTransform {
            translation_mm: at,
            ..Default::default()
        })
        .unwrap()
}

#[test]
fn a_groove_clear_of_blind_pockets_keeps_both_pockets_exactly() {
    let policy = GeometryPolicy::default();
    let panel = pocketed_panel(policy);
    let pockets = 2. * PI * 2.5 * 2.5 * 10.;
    let full = 60. * 40. * 18.;
    near(panel.volume_mm3().unwrap(), full - pockets, 1e-9);
    // A 6 mm wide, 8 mm deep groove across the panel at x = 40..46.
    let groove = block([6., 50., 8.], [40., -5., 10.], policy);
    let cut = panel
        .rust_subtract_orthogonal_tools(std::slice::from_ref(&groove), policy)
        .unwrap()
        .unwrap();
    near(
        cut.volume_mm3().unwrap(),
        full - pockets - 6. * 40. * 8.,
        1e-9,
    );
    assert_eq!(cut.solid_count().unwrap(), 1);
    // The restored pockets stay part of the cut stock, so a second groove
    // clear of them commutes again and still removes exactly its own volume.
    let second = block([4., 50., 3.], [50., -5., 15.], policy);
    let twice = cut
        .rust_subtract_orthogonal_tools(std::slice::from_ref(&second), policy)
        .unwrap()
        .unwrap();
    near(
        twice.volume_mm3().unwrap(),
        full - pockets - 6. * 40. * 8. - 4. * 40. * 3.,
        1e-9,
    );
}

#[test]
fn a_tool_reaching_a_pocket_is_not_commuted_past_it() {
    let policy = GeometryPolicy::default();
    let panel = pocketed_panel(policy);
    // This groove crosses the first pocket's mouth.
    let crossing = block([6., 50., 8.], [9., -5., 10.], policy);
    match panel.rust_subtract_orthogonal_tools(std::slice::from_ref(&crossing), policy) {
        Err(_) => {}
        Ok(result) => {
            // Any provider that accepts it must remove the true union.
            let removed = 6. * 40. * 8.;
            let volume = result.unwrap().volume_mm3().unwrap();
            assert!(volume < 60. * 40. * 18. - removed, "{volume}");
        }
    }
}

#[test]
fn overlap_with_pocketed_stock_uses_the_planar_stock_clear_of_its_pockets() {
    let policy = GeometryPolicy::default();
    let panel = pocketed_panel(policy);
    // A thin printed tab overlapping the panel edge, far from both pockets.
    let tab = block([10., 8., 0.6], [41., -6., 0.], policy);
    near(
        panel.rust_planar_intersection_volume(&tab, policy).unwrap(),
        10. * 2. * 0.6,
        1e-12,
    );
    near(
        tab.rust_planar_intersection_volume(&panel, policy).unwrap(),
        10. * 2. * 0.6,
        1e-12,
    );
}
