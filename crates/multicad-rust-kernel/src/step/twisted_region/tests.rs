use super::*;

fn spiral(angle_rad: [f64; 2], radius_mm: [f64; 2]) -> Spiral {
    Spiral {
        angle_rad,
        radius_mm,
    }
}

/// A curved triangle in polar coordinates: two affine spirals and one arc.
fn triangle(offset: f64) -> Vec<Spiral> {
    vec![
        spiral([offset, offset + 1.], [1., 1.5]),
        spiral([offset + 1., offset + 0.5], [1.5, 2.]),
        spiral([offset + 0.5, offset], [2., 1.]),
    ]
}

#[test]
fn twisted_spiral_loops_read_back_as_closed_solids_of_area_times_height() {
    // A twisted extrusion keeps every section congruent, so its volume is
    // the exact section area times the rise, independent of the twist.
    // This loop runs clockwise; the encoder orients it itself.
    let area: f64 = -triangle(0.)
        .iter()
        .map(|piece| piece.signed_area())
        .sum::<f64>();
    assert!(area > 0.1);
    let tolerance = 1e-7;
    for (twist, origin) in [
        (0.37, [0.; 3]),
        (-2.9, [40., -15., 7.]),
        (12.5, [0., 0., -3.]),
    ] {
        let loops = [triangle(0.3), triangle(0.3 + std::f64::consts::PI)];
        let heights = [1.25, 3.75];
        let text = write_twisted_spirals(
            &loops,
            heights,
            twist,
            origin,
            [1., 0., 0.],
            [0., 1., 0.],
            tolerance,
        )
        .unwrap();
        let solids = crate::step_read::brep::read_solids(&text, tolerance, None).unwrap();
        assert_eq!(solids.len(), 2);
        for solid in &solids {
            let expected = area * (heights[1] - heights[0]);
            assert!(
                (solid.volume_mm3() - expected).abs() < 1e-6,
                "{} vs {expected}",
                solid.volume_mm3()
            );
            let [low, high] = solid.bounds_mm();
            assert!((low[2] - origin[2] - heights[0]).abs() < 1e-9);
            assert!((high[2] - origin[2] - heights[1]).abs() < 1e-9);
            assert!(high[0] - origin[0] <= 2. + 1e-6 && low[0] - origin[0] >= -2. - 1e-6);
        }
    }
    // A reversed loop describes the same material.
    let reversed: Vec<_> = triangle(0.).iter().rev().map(|p| p.reversed()).collect();
    let text = write_twisted_spirals(
        &[reversed],
        [0., 1.],
        0.5,
        [0.; 3],
        [1., 0., 0.],
        [0., 1., 0.],
        tolerance,
    )
    .unwrap();
    let solids = crate::step_read::brep::read_solids(&text, tolerance, None).unwrap();
    assert!((solids[0].volume_mm3() - area).abs() < 1e-6);
}

#[test]
fn open_degenerate_and_unbounded_loops_are_refused() {
    let frame = ([0.; 3], [1., 0., 0.], [0., 1., 0.]);
    let write = |loops: &[Vec<Spiral>], heights, twist| {
        write_twisted_spirals(loops, heights, twist, frame.0, frame.1, frame.2, 1e-7)
    };
    let mut open = triangle(0.);
    open[2].radius_mm[1] = 1.1;
    assert_eq!(
        write(&[open], [0., 1.], 1.),
        Err(StepError::InvalidGeometry)
    );
    let mut radial = triangle(0.);
    radial[0].angle_rad[1] = radial[0].angle_rad[0];
    assert_eq!(
        write(&[radial], [0., 1.], 1.),
        Err(StepError::InvalidGeometry)
    );
    assert_eq!(
        write(&[triangle(0.)], [1., 1.], 1.),
        Err(StepError::InvalidGeometry)
    );
    assert_eq!(
        write(&[triangle(0.)], [0., 1.], 0.),
        Err(StepError::InvalidGeometry)
    );
    assert_eq!(write(&[], [0., 1.], 1.), Err(StepError::InvalidGeometry));
    let mut negative = triangle(0.);
    negative[1].radius_mm[0] = -1.;
    assert_eq!(
        write(&[negative], [0., 1.], 1.),
        Err(StepError::InvalidGeometry)
    );
}
