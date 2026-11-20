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

fn moved(solid: Solid, translation_mm: [f64; 3]) -> Solid {
    solid
        .transformed(RigidTransform {
            translation_mm,
            ..Default::default()
        })
        .unwrap()
}

fn block(size: [f64; 3], at: [f64; 3]) -> Solid {
    moved(Solid::box_solid(size, policy()).unwrap(), at)
}

fn common(a: &Solid, b: &Solid) -> Option<Solid> {
    a.boolean(b, BooleanOperation::Intersect, policy()).unwrap()
}

/// Area of the disk of radius `r` on the far side of the chord at `d`.
fn segment(r: f64, d: f64) -> f64 {
    r * r * (d / r).acos() - d * (r * r - d * d).sqrt()
}

#[test]
fn a_box_face_through_a_cylinder_leaves_an_exact_circular_segment() {
    let cylinder = Solid::cylinder(Default::default(), 5., 3., policy()).unwrap();
    let clip = block([5., 12., 5.], [2., -6., -1.]);
    let clipped = common(&cylinder, &clip).unwrap();
    near(clipped.volume_mm3().unwrap(), 3. * segment(5., 2.));
    let bounds = clipped.bounds().unwrap();
    near(bounds.min[0], 2.);
    near(bounds.max[0], 5.);
    near(bounds.max[2] - bounds.min[2], 3.);
    // The section is one chord and one arc of the source circle.
    let arcs = clipped
        .faces()
        .unwrap()
        .into_iter()
        .filter(|face| face.plane.is_none())
        .count();
    assert_eq!(arcs, 1);
    // The operation is symmetric in its operands.
    let swapped = common(&clip, &cylinder).unwrap();
    near(swapped.volume_mm3().unwrap(), clipped.volume_mm3().unwrap());
}

#[test]
fn a_corner_in_the_disk_gives_two_chords_and_an_arc() {
    let cylinder = Solid::cylinder(Default::default(), 5., 4., policy()).unwrap();
    let clip = block([10., 10., 2.], [1., 1., 1.]);
    let clipped = common(&cylinder, &clip).unwrap();
    // Region x > 1, y > 1 of the disk, swept over z in [1, 3]:
    // the integral of sqrt(25 - x^2) - 1 for x from 1 to sqrt(24).
    let x = 24f64.sqrt();
    let primitive = |x: f64| 0.5 * x * (25. - x * x).sqrt() + 12.5 * (x / 5.).asin();
    let area = primitive(x) - primitive(1.) - (x - 1.);
    near(clipped.volume_mm3().unwrap(), 2. * area);
    let bounds = clipped.bounds().unwrap();
    near(bounds.min[2], 1.);
    near(bounds.max[2], 3.);
}

#[test]
fn a_disk_inside_or_outside_the_section_trims_or_vanishes() {
    // A cylinder along Y through a plate profile in the XZ plane.
    let cylinder = Solid::cylinder(
        PlaneFrame {
            origin: [240., -150., 780.],
            u: [0., 0., 1.],
            v: [1., 0., 0.],
        },
        8.,
        500.,
        policy(),
    )
    .unwrap();
    let plate = block([240., 40., 140.], [180., -40., 760.]);
    let trimmed = common(&cylinder, &plate).unwrap();
    near(trimmed.volume_mm3().unwrap(), 64. * PI * 40.);
    let bounds = trimmed.bounds().unwrap();
    near(bounds.min[1], -40.);
    near(bounds.max[1], 0.);
    assert!(trimmed.axial_cylinder_frame().is_some());
    let apart = block([240., 40., 140.], [300., -40., 800.]);
    assert!(common(&cylinder, &apart).is_none());
    // Beyond the cylinder's axial range there is nothing in common either.
    let beyond = block([240., 40., 140.], [180., 400., 760.]);
    assert!(common(&cylinder, &beyond).is_none());
}

#[test]
fn tangent_contacts_and_skew_axes_are_not_claimed() {
    let cylinder = Solid::cylinder(Default::default(), 5., 3., policy()).unwrap();
    for clip in [
        // The box face is tangent to the cylinder.
        block([5., 12., 5.], [5., -6., -1.]),
        // A box corner lies on the circle.
        block([5., 5., 5.], [3., 4., -1.]),
        // The box end touches the cylinder cap within tolerance.
        block([5., 12., 3. + 1e-9], [2., -6., 0.]),
    ] {
        assert!(
            cylinder
                .parallel_extrusion_intersection(&clip, BooleanOperation::Intersect, policy())
                .is_none()
        );
    }
    let tilted = Solid::cylinder(
        PlaneFrame {
            u: [1., 0., 0.],
            v: [0., 0.6, 0.8],
            ..Default::default()
        },
        5.,
        3.,
        policy(),
    )
    .unwrap();
    let clip = block([5., 12., 5.], [2., -6., -1.]);
    assert!(
        tilted
            .parallel_extrusion_intersection(&clip, BooleanOperation::Intersect, policy())
            .is_none()
    );
    for operation in [BooleanOperation::Union, BooleanOperation::Subtract] {
        assert!(
            cylinder
                .parallel_extrusion_intersection(&clip, operation, policy())
                .is_none()
        );
    }
}

#[test]
fn polygon_prisms_with_holes_classify_a_clear_disk() {
    let r = |n| Rational::new(n, 1).unwrap();
    let feature = EntityId::new();
    let origins: Vec<_> = (0..9)
        .map(|i| Origin {
            source: feature,
            role: format!("edge/{i}"),
        })
        .collect();
    let outer = [[0, 0], [40, 0], [40, 40], [0, 40]]
        .map(|p| p.map(r))
        .to_vec();
    let window = [[20, 20], [20, 30], [30, 30], [30, 20]]
        .map(|p| p.map(r))
        .to_vec();
    let prism = PolygonPrism::with_circular_holes(
        vec![outer, window],
        vec![Circle::new([r(8), r(30)], r(3)).unwrap()],
        r(6),
    )
    .unwrap();
    let plate = Solid::from_polygon_prism_sweep(
        prism,
        PlaneFrame::default(),
        false,
        feature,
        &origins,
        &[false; 9],
        policy(),
    )
    .unwrap();
    let rod = |x: f64, y: f64| {
        Solid::cylinder(
            PlaneFrame {
                origin: [x, y, -5.],
                ..Default::default()
            },
            2.,
            20.,
            policy(),
        )
        .unwrap()
    };
    near(
        common(&rod(10., 10.), &plate)
            .unwrap()
            .volume_mm3()
            .unwrap(),
        4. * PI * 6.,
    );
    // Inside the square window and inside the round hole: no material.
    assert!(common(&rod(25., 25.), &plate).is_none());
    assert!(common(&rod(8., 30.), &plate).is_none());
    // Across the window edge of a holed profile the provider declines.
    assert!(
        rod(20., 25.)
            .parallel_extrusion_intersection(&plate, BooleanOperation::Intersect, policy())
            .is_none()
    );
}

#[test]
fn rectilinear_cells_holding_the_whole_disk_trim_the_cylinder() {
    // An L-shaped beam: a long bar with a raised block on one end.
    let bar = block([150., 40., 40.], [225., 230., 760.]);
    let riser = block([30., 40., 30.], [345., 230., 800.]);
    let beam = bar
        .boolean(&riser, BooleanOperation::Union, policy())
        .unwrap()
        .unwrap();
    assert!(beam.orthogonal_body().is_some());
    let along_y = |x: f64, z: f64| {
        Solid::cylinder(
            PlaneFrame {
                origin: [x, -150., z],
                u: [0., 0., 1.],
                v: [1., 0., 0.],
            },
            8.,
            500.,
            policy(),
        )
        .unwrap()
    };
    let through = common(&along_y(240., 780.), &beam).unwrap();
    near(through.volume_mm3().unwrap(), 64. * PI * 40.);
    let bounds = through.bounds().unwrap();
    near(bounds.min[1], 230.);
    near(bounds.max[1], 270.);
    // Above the bar and beside the riser nothing is in common.
    assert!(common(&along_y(240., 815.), &beam).is_none());
    // A disk across the riser's side is only partly covered.
    assert!(
        along_y(345., 815.)
            .parallel_extrusion_intersection(&beam, BooleanOperation::Intersect, policy())
            .is_none()
    );
    // A vertical cylinder sees two stacked cells as one interval.
    let vertical = Solid::cylinder(
        PlaneFrame {
            origin: [360., 250., 700.],
            ..Default::default()
        },
        5.,
        200.,
        policy(),
    )
    .unwrap();
    let stacked = common(&vertical, &beam).unwrap();
    near(stacked.volume_mm3().unwrap(), 25. * PI * 70.);

    // A bolt hole far from a cylinder leaves its cells authoritative; one
    // inside the cylinder's region withholds them.
    let bolt = CircularPocketTool {
        frame: PlaneFrame {
            origin: [300., 250., 750.],
            ..Default::default()
        },
        center_uv: [0.; 2],
        radius_mm: 4.,
        displacement_mm: [0., 0., 60.],
        origin: Origin {
            source: EntityId::new(),
            role: "bolt".into(),
        },
    };
    let drilled = beam
        .rust_drill_circular_machining(&[bolt], &[true], policy())
        .unwrap();
    let through = common(&along_y(240., 780.), &drilled).unwrap();
    near(through.volume_mm3().unwrap(), 64. * PI * 40.);
    assert!(
        along_y(300., 780.)
            .parallel_extrusion_intersection(&drilled, BooleanOperation::Intersect, policy())
            .is_none()
    );

    // A guide rod through a coaxial clearance bore meets no material; a
    // bore that only grazes the rod's disk is not resolved.
    let clearance = |radius_mm| CircularPocketTool {
        frame: PlaneFrame {
            origin: [240., 220., 780.],
            u: [0., 0., 1.],
            v: [1., 0., 0.],
        },
        center_uv: [0.; 2],
        radius_mm,
        displacement_mm: [0., 60., 0.],
        origin: Origin {
            source: EntityId::new(),
            role: "clearance".into(),
        },
    };
    let bored = beam
        .rust_drill_circular_machining(&[clearance(8.5)], &[true], policy())
        .unwrap();
    assert!(common(&along_y(240., 780.), &bored).is_none());
    let tight = beam
        .rust_drill_circular_machining(&[clearance(6.)], &[true], policy())
        .unwrap();
    assert!(
        along_y(240., 780.)
            .parallel_extrusion_intersection(&tight, BooleanOperation::Intersect, policy())
            .is_none()
    );
}
