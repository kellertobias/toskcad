use super::*;
use std::f64::consts::PI;

fn policy() -> GeometryPolicy {
    GeometryPolicy::default()
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-7 * expected.abs().max(1.),
        "{actual} != {expected}"
    );
}

fn block(size: [f64; 3], at: [f64; 3]) -> Solid {
    Solid::box_solid(size, policy())
        .unwrap()
        .transformed(RigidTransform {
            translation_mm: at,
            ..Default::default()
        })
        .unwrap()
}

fn rod(center: [f64; 2], radius: f64) -> Solid {
    Solid::cylinder(
        PlaneFrame {
            origin: [center[0], center[1], -1.],
            ..Default::default()
        },
        radius,
        20.,
        policy(),
    )
    .unwrap()
}

fn cut(a: &Solid, b: &Solid) -> Solid {
    a.boolean(b, BooleanOperation::Subtract, policy())
        .unwrap()
        .unwrap()
}

const PLATE: f64 = 120. * 90. * 18.;
const BORES: f64 = 3. * PI * 6.1 * 6.1 * 18.;
/// A 7 mm wide, 40 mm long slot with two half-round ends.
const SLOT: f64 = (7. * 40. + PI * 3.5 * 3.5) * 18.;

/// The drill guide body: three bushing bores, then two round-ended slots,
/// each cut as a middle block and two end cylinders.
#[test]
fn bores_and_round_ended_through_slots_give_one_exact_profile_sweep() {
    let mut body = [30., 60., 90.]
        .into_iter()
        .fold(block([120., 90., 18.], [0.; 3]), |body, x| {
            cut(&body, &rod([x, 65.], 6.1))
        });
    for (x, left) in [(12., 8.5), (108., 104.5)] {
        body = cut(&body, &block([7., 40., 20.], [left, 20., -1.]));
        assert!(body.normal_sweep_body().is_some());
        for y in [20., 60.] {
            body = cut(&body, &rod([x, y], 3.5));
        }
    }
    near(body.volume_mm3().unwrap(), PLATE - BORES - 2. * SLOT);
    let bounds = body.bounds().unwrap();
    near(bounds.min[2], 0.);
    near(bounds.max[2], 18.);
    // Each slot end is one half-round wall: four in all, beside the three
    // bore walls.
    let curved = body
        .faces()
        .unwrap()
        .into_iter()
        .filter(|face| face.plane.is_none())
        .count();
    assert_eq!(curved, 7);
    let (_, _, _, rings) = body.normal_sweep_section().unwrap();
    assert_eq!(rings.len(), 6);
}

#[test]
fn bored_box_stock_takes_a_clear_through_window_as_a_hole() {
    let plate = block([60., 40., 10.], [0.; 3]);
    let window = block([10., 10., 14.], [20., 15., -2.]);
    let cut = cut(&cut(&plate, &rod([45., 20.], 4.)), &window);
    near(cut.volume_mm3().unwrap(), (2400. - 16. * PI - 100.) * 10.);
    assert!(cut.normal_sweep_body().is_some());
}

#[test]
fn overlaps_tangencies_and_blind_tools_are_not_profile_cuts() {
    let policy = policy();
    let slotted = cut(
        &cut(&block([120., 90., 18.], [0.; 3]), &rod([60., 65.], 6.1)),
        &block([7., 40., 20.], [8.5, 20., -1.]),
    );
    let declined = |tool: &Solid| {
        slotted
            .through_profile_cut(tool, BooleanOperation::Subtract, policy)
            .is_none()
    };
    // Wider or off-centre end cylinders do not round the slot exactly.
    assert!(declined(&rod([12., 20.], 3.6)));
    assert!(declined(&rod([12.5, 20.], 3.5)));
    // A window crossing the bore, or touching the outer edge.
    assert!(declined(&block([10., 10., 20.], [55., 60., -1.])));
    assert!(declined(&block([10., 10., 20.], [110., 40., -1.])));
    // A blind tool stops inside the stock.
    assert!(declined(&block([10., 10., 10.], [80., 30., 9.])));
    let tool = rod([12., 20.], 3.5);
    for operation in [BooleanOperation::Union, BooleanOperation::Intersect] {
        assert!(
            slotted
                .through_profile_cut(&tool, operation, policy)
                .is_none()
        );
    }
}

#[test]
fn a_second_through_hole_in_a_named_sweep_keeps_consistent_metadata() {
    let profile = CurvedProfile {
        frame: PlaneFrame::default(),
        outer: [[0., 0.], [60., 0.], [60., 40.], [0., 40.]]
            .iter()
            .zip([[60., 0.], [60., 40.], [0., 40.], [0., 0.]])
            .map(|(&a, b)| CurvedEdge2d::Line { a, b })
            .collect(),
        holes: vec![vec![CurvedEdge2d::Circle {
            center: [10., 10.],
            radius_mm: 2.,
        }]],
    };
    let plate = Solid::extrude_curved_profile_vector(&profile, [0., 0., 8.], policy()).unwrap();
    let once = cut(&plate, &rod([30., 20.], 3.));
    let twice = cut(&once, &rod([50., 20.], 3.));
    near(
        twice.volume_mm3().unwrap(),
        (2400. - 4. * PI - 18. * PI) * 8.,
    );
    assert!(twice.topology_snapshot(EntityId::new(), 0).is_ok());
}
