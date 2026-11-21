//! Exact minimum distance between parallel circular bodies, and between
//! parallel planar faces bounded by circular arcs.
use multicad_core::{document::EntityId, references::Origin};
use multicad_kernel::{
    rust_kernel::{Circle, Rational, polygon::PolygonPrism},
    *,
};

fn p() -> GeometryPolicy {
    GeometryPolicy::default()
}
fn at(origin: [f64; 3]) -> PlaneFrame {
    PlaneFrame {
        origin,
        ..Default::default()
    }
}
fn length(a: [f64; 3], b: [f64; 3]) -> f64 {
    (0..3).map(|i| (a[i] - b[i]).powi(2)).sum::<f64>().sqrt()
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-9, "{a} != {b}");
}
fn whole(a: &Solid, b: &Solid) -> MinimumDistance {
    a.minimum_distance(BRepSelection::Whole, b, BRepSelection::Whole)
        .unwrap()
}
fn check(result: &MinimumDistance) {
    for w in &result.witnesses {
        near(length(w.a_mm, w.b_mm), result.distance_mm);
    }
}

#[test]
fn parallel_rings_split_into_radial_and_axial_gaps_with_material_witnesses() {
    let tube = Solid::rust_annulus(6., 10., 10., EntityId::new(), at([0.; 3]), p()).unwrap();
    let pin = |origin, radius, height| Solid::cylinder(at(origin), radius, height, p()).unwrap();
    // Off-center pin inside the bore: the radial gap along the offset.
    let inside = whole(&tube, &pin([1., 0., 2.], 2., 6.));
    near(inside.distance_mm, 3.);
    check(&inside);
    let w = inside.witnesses[0];
    near(w.a_mm[0].hypot(w.a_mm[1]), 6.);
    near((w.b_mm[0] - 1.).hypot(w.b_mm[1]), 2.);
    // An axial overlap leaves a continuum of closest pairs.
    assert_eq!(inside.total_solutions, None);
    assert!(!inside.witnesses_complete && !inside.inner_solution);
    // Radially overlapping sections above the tube: the axial gap alone.
    let above = whole(&tube, &pin([8., 0., 12.], 1., 3.));
    near(above.distance_mm, 2.);
    check(&above);
    near(above.witnesses[0].a_mm[2], 10.);
    near(above.witnesses[0].b_mm[2], 12.);
    // Outside, in either order.
    let outside = pin([20., 0., 0.], 2., 10.);
    for (result, swapped) in [
        (whole(&tube, &outside), false),
        (whole(&outside, &tube), true),
    ] {
        near(result.distance_mm, 8.);
        check(&result);
        let w = result.witnesses[0];
        let (ring, rod) = if swapped {
            (w.b_mm, w.a_mm)
        } else {
            (w.a_mm, w.b_mm)
        };
        near(ring[0].hypot(ring[1]), 10.);
        near((rod[0] - 20.).hypot(rod[1]), 2.);
    }
    // Both gaps positive off center: one unique closest pair.
    let corner = whole(&tube, &pin([1., 0., 12.], 2., 3.));
    near(corner.distance_mm, 13_f64.sqrt());
    assert_eq!(corner.total_solutions, Some(1));
    assert!(corner.witnesses_complete);
    let w = corner.witnesses[0];
    for (actual, expected) in w.a_mm.into_iter().zip([6., 0., 10.]) {
        near(actual, expected);
    }
    for (actual, expected) in w.b_mm.into_iter().zip([3., 0., 12.]) {
        near(actual, expected);
    }
    // Overlapping material.
    let crossing = whole(&tube, &pin([8., 0., 2.], 1., 3.));
    near(crossing.distance_mm, 0.);
    assert!(crossing.inner_solution);
    // An opposite axial sense: the ring occupies z in [0, 10] downward.
    let flipped = Solid::rust_annulus(
        6.,
        10.,
        10.,
        EntityId::new(),
        PlaneFrame {
            origin: [0., 0., 10.],
            u: [1., 0., 0.],
            v: [0., -1., 0.],
        },
        p(),
    )
    .unwrap();
    let centered = whole(&flipped, &pin([0., 0., 12.], 1., 3.));
    near(centered.distance_mm, 29_f64.sqrt());
    check(&centered);
    assert_eq!(centered.total_solutions, None);
}

fn plate() -> Solid {
    let q = |n| Rational::new(n, 1).unwrap();
    let prism = PolygonPrism::with_circular_holes(
        vec![
            vec![[0, 0], [20, 0], [20, 16], [0, 16]]
                .into_iter()
                .map(|p| p.map(q))
                .collect(),
        ],
        vec![Circle::new([q(10), q(8)], q(1)).unwrap()],
        q(4),
    )
    .unwrap();
    let f = EntityId::new();
    let origins: Vec<_> = (0..5)
        .map(|i| Origin {
            source: f,
            role: format!("edge.{i}"),
        })
        .collect();
    Solid::from_polygon_prism_sweep(prism, at([0.; 3]), false, f, &origins, &[false; 5], p())
        .unwrap()
}
fn block(size: [i64; 3], origin: [f64; 3]) -> Solid {
    Solid::rust_box(
        EntityId::new().uuid_bits(),
        size.map(|n| Rational::new(n, 1).unwrap()),
        at(origin),
        p(),
    )
    .unwrap()
}
fn face(solid: &Solid, normal: [f64; 3]) -> BRepSelection {
    let face = solid
        .faces()
        .unwrap()
        .into_iter()
        .find(|f| f.plane.is_some_and(|(_, n)| n == normal))
        .unwrap();
    BRepSelection::Face(face.index)
}

#[test]
fn parallel_faces_with_hole_arcs_measure_their_separation_or_refuse() {
    let plate = plate();
    let top = face(&plate, [0., 0., 1.]);
    let measure = |other: &Solid| plate.minimum_distance(top, other, face(other, [0., 0., -1.]));
    let over = block([4, 4, 2], [2., 2., 10.]);
    let result = measure(&over).unwrap();
    near(result.distance_mm, 6.);
    check(&result);
    let w = result.witnesses[0];
    near(w.a_mm[2], 4.);
    near(w.b_mm[2], 10.);
    // The witness lies on both faces: inside the plate top, off the hole,
    // and under the block.
    assert!((0. ..=20.).contains(&w.a_mm[0]) && (0. ..=16.).contains(&w.a_mm[1]));
    assert!((w.a_mm[0] - 10.).hypot(w.a_mm[1] - 8.) > 1.);
    assert!((2. ..=6.).contains(&w.b_mm[0]) && (2. ..=6.).contains(&w.b_mm[1]));
    assert_eq!(result.total_solutions, None);
    // A block covering the hole still meets the plate material.
    let covering = block([6, 6, 2], [7., 5., 10.]);
    near(measure(&covering).unwrap().distance_mm, 6.);
    // Projections that meet only over the hole, or not at all, are refused:
    // the separation is not then attained.
    for other in [
        block([1, 1, 1], [9.5, 7.5, 10.]),
        block([4, 4, 2], [40., 0., 10.]),
    ] {
        let refusal = measure(&other).unwrap_err();
        assert_eq!(refusal.kind, GeometryErrorKind::Unsupported);
        assert!(
            refusal.message.contains("overlapping normal projections"),
            "{}",
            refusal.message
        );
    }
}

#[test]
fn mating_worm_and_generated_wheel_measure_a_certified_flank_gap() {
    use multicad_kernel::{GeneratedZaWormWheel, RigidTransform, ZaWorm};
    let worm = Solid::rust_za_worm_named(
        ZaWorm {
            axial_module_mm: 1.,
            threaded_length_mm: 16.,
            bore_radius_mm: 1.,
            axial_thinning_mm: 0.05,
        },
        EntityId::new(),
        p(),
    )
    .unwrap();
    let wheel = Solid::rust_generated_za_worm_wheel_named(
        GeneratedZaWormWheel {
            teeth: 32,
            axial_module_mm: 1.,
            face_width_mm: 2.,
            bore_radius_mm: 2.,
            pitch_thinning_mm: 0.05,
        },
        EntityId::new(),
        p(),
    )
    .unwrap();
    let place = |angle: f64| {
        let (s, c) = angle.sin_cos();
        wheel
            .transformed(RigidTransform {
                rotation: [[c, -s, 0.], [0., 0., 1.], [-s, -c, 0.]],
                translation_mm: [-21., 0., 0.],
            })
            .unwrap()
    };
    let placed = place(0.);
    let (_, witness) = worm
        .certify_rust_worm_pair(&placed, p(), 0.2, None)
        .unwrap();
    let found = worm
        .minimum_distance(BRepSelection::Whole, &placed, BRepSelection::Whole)
        .unwrap();
    check(&found);
    // The clearance witness is an actual pair, so it is never closer; the
    // gap lies between the total thinning projected onto the flank normal
    // and the total thinning.
    assert!(found.distance_mm <= witness.unwrap().distance_mm + 1e-12);
    assert!(found.distance_mm > 0.04 && found.distance_mm < 0.1);
    // Swapped order keeps the witness owners.
    let swapped = placed
        .minimum_distance(BRepSelection::Whole, &worm, BRepSelection::Whole)
        .unwrap();
    near(swapped.distance_mm, found.distance_mm);
    near(
        length(swapped.witnesses[0].a_mm, found.witnesses[0].b_mm),
        0.,
    );
    // Half a tooth out of phase the solids overlap: no distance is claimed.
    let wrong = place(std::f64::consts::PI / 32.);
    let refusal = worm
        .minimum_distance(BRepSelection::Whole, &wrong, BRepSelection::Whole)
        .unwrap_err();
    assert_eq!(refusal.kind, GeometryErrorKind::Unsupported);
}
