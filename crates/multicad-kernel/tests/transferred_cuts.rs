use multicad_core::{
    document::EntityId,
    references::{AnalyticGeometry, Origin, TopologyKind},
};

use multicad_kernel::*;
fn p() -> GeometryPolicy {
    Default::default()
}
fn origin(id: EntityId, role: &str) -> Origin {
    Origin {
        source: id,
        role: role.into(),
    }
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-5, "{a} != {b}");
}
fn rect(frame: PlaneFrame, a: [f64; 2], b: [f64; 2]) -> CurvedProfile {
    let points = [a, [b[0], a[1]], b, [a[0], b[1]]];
    CurvedProfile {
        frame,
        outer: (0..4)
            .map(|i| CurvedEdge2d::Line {
                a: points[i],
                b: points[(i + 1) % 4],
            })
            .collect(),
        holes: vec![],
    }
}
fn circle(frame: PlaneFrame, center: [f64; 2], r: f64) -> CurvedProfile {
    CurvedProfile {
        frame,
        outer: vec![CurvedEdge2d::Circle {
            center,
            radius_mm: r,
        }],
        holes: vec![],
    }
}
fn names(n: usize) -> Vec<Origin> {
    (0..n)
        .map(|_| origin(EntityId::new(), "source-curve"))
        .collect()
}
fn xy(z: f64) -> PlaneFrame {
    PlaneFrame {
        origin: [0., 0., z],
        ..Default::default()
    }
}

#[test]
fn rust_transfer_preserves_disjoint_through_bores_and_rejects_a_cut_into_their_wall() {
    use multicad_kernel::rust_kernel::Rational;

    let stock_id = EntityId::new();
    let stock = Solid::rust_box(
        stock_id.uuid_bits(),
        [12, 8, 3].map(|n| Rational::new(n, 1).unwrap()),
        PlaneFrame::default(),
        p(),
    )
    .unwrap();
    let round_id = EntityId::new();
    let entry = xy(3.);
    let through = stock.transfer_through_depth(entry, true, p()).unwrap();
    let bored = stock
        .rust_cut_transfer_profile(
            &circle(entry, [3., 4.], 1.),
            through,
            true,
            round_id,
            &[origin(round_id, "mounting-bore")],
            p(),
        )
        .unwrap();
    let before = 288. - std::f64::consts::PI * 3.;
    near(bored.volume_mm3().unwrap(), before);

    let slot_id = EntityId::new();
    let edges = names(4);
    let result = bored
        .rust_cut_transfer_profile(
            &rect(entry, [7., 3.], [9., 5.]),
            through,
            true,
            slot_id,
            &edges,
            p(),
        )
        .unwrap();
    near(result.volume_mm3().unwrap(), before - 12.);
    let key_origin = |solid: &Solid| {
        solid
            .topology_snapshot(stock_id, 1)
            .unwrap()
            .elements
            .into_iter()
            .filter(|element| element.key.kind == TopologyKind::Edge)
            .filter(|element| {
                element
                    .origins
                    .iter()
                    .any(|origin| origin.source == round_id)
            })
            .filter_map(|element| element.geometry)
            .filter(|geometry| matches!(geometry, AnalyticGeometry::CircularArc { .. }))
            .count()
    };
    assert_eq!(key_origin(&bored), 2);
    assert_eq!(key_origin(&result), 2);

    assert!(
        bored
            .rust_cut_transfer_profile(
                &rect(entry, [3., 3.], [5., 5.]),
                through,
                true,
                slot_id,
                &edges,
                p(),
            )
            .is_err()
    );
    near(bored.volume_mm3().unwrap(), before);
    near(stock.volume_mm3().unwrap(), 288.);
}
