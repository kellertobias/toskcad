use super::*;
use std::collections::BTreeMap;

fn shaft(steps: &[(f64, [f64; 2])]) -> Revolution {
    let mut ring = vec![[0., steps[0].1[0]]];
    for &(radius, [low, high]) in steps {
        ring.push([radius, low]);
        ring.push([radius, high]);
    }
    ring.push([0., steps[steps.len() - 1].1[1]]);
    Revolution::new(vec![ring], TAU).unwrap()
}

fn tool(base: V, axis: V, radius: f64, height: f64) -> Tool {
    Tool {
        base_mm: base,
        axis,
        radius_mm: radius,
        height_mm: height,
        literal: true,
    }
}

fn body(cut: Result<Cut, Error>) -> CrossBore {
    match cut {
        Ok(Cut::Body(body)) => *body,
        other => panic!("expected a cut, got {other:?}"),
    }
}

/// Removed volume by slices across the tool travel: each slice of the tool
/// section is a disc clipped by two chords of the shaft, whose area is
/// elementary. Independent of the angle integrals of the provider.
fn removed_by_slices(radius: f64, offset: f64, rho: f64, range: [f64; 2]) -> f64 {
    let segment = |d: f64| {
        if d >= rho {
            0.
        } else if d <= -rho {
            PI * rho * rho
        } else {
            rho * rho * (d / rho).portable_acos() - d * (rho * rho - d * d).sqrt()
        }
    };
    let n = 400_000;
    let [a, b] = [range[0].max(-radius), range[1].min(radius)];
    let h = (b - a) / f64::from(n);
    (0..n)
        .map(|k| {
            let x = a + (f64::from(k) + 0.5) * h;
            let chord = (radius * radius - x * x).max(0.).sqrt();
            PI * rho * rho - segment(chord - offset) - segment(chord + offset)
        })
        .sum::<f64>()
        * h
}

fn closed(mesh: &DisplayMesh) -> bool {
    let mut directed: BTreeMap<(u32, u32), i32> = BTreeMap::new();
    for t in mesh.indices.as_chunks::<3>().0 {
        for i in 0..3 {
            let (a, b) = (t[i], t[(i + 1) % 3]);
            *directed.entry((a, b)).or_default() += 1;
        }
    }
    directed
        .iter()
        .all(|(&(a, b), &count)| count == 1 && directed.get(&(b, a)) == Some(&1))
}

fn signed_volume(mesh: &DisplayMesh) -> f64 {
    let p = |i: u32| -> V { std::array::from_fn(|k| mesh.vertices[3 * i as usize + k]) };
    mesh.indices
        .as_chunks::<3>()
        .0
        .iter()
        .map(|t| dot(p(t[0]), cross(p(t[1]), p(t[2]))) / 6.)
        .sum()
}

fn options() -> Tessellation {
    Tessellation {
        deflection: 0.01,
        angle: 0.2,
        triangle_limit: 2_000_000,
    }
}

#[test]
fn radial_through_bore_matches_crossing_cylinders() {
    let hub = shaft(&[(9., [-15., 15.])]);
    let cut = body(subtract(
        &hub,
        tool([-50., 0., 2.], [1., 0., 0.], 2.75, 100.),
        None,
    ));
    assert_eq!(cut.ends(), [End::Through; 2]);
    // The crossing-cylinder volume in the form used by the workbench oracle.
    let crossing = {
        let n = 4096;
        let step = PI / f64::from(n);
        (0..n)
            .map(|k| {
                let t = -PI / 2. + (f64::from(k) + 0.5) * step;
                t.portable_cos().powi(2) * (81. - 2.75f64.powi(2) * t.portable_sin().powi(2)).sqrt()
            })
            .sum::<f64>()
            * step
            * 4.
            * 2.75
            * 2.75
    };
    assert!((cut.removed_mm3() - crossing).abs() < 1e-10 * crossing);
    assert!(cut.volume_error_mm3() < 1e-9);
    let expected = PI * 81. * 30. - crossing;
    assert!((cut.volume_mm3() - expected).abs() <= cut.volume_error_mm3() + 1e-10 * expected);
    let slices = removed_by_slices(9., 0., 2.75, [-50., 50.]);
    assert!((cut.removed_mm3() - slices).abs() < 1e-6 * slices);
    // Both holes are symmetric about the axis and at the tool level.
    let wall = cut.bored_wall();
    let holes = TAU * 9. * 30. - wall.area_mm2;
    assert!(holes > 2. * PI * 2.75 * 2.75);
    assert!(wall.centroid_mm.iter().take(2).all(|c| c.abs() < 1e-9));
    let tool_wall = cut.tool_wall();
    assert!((tool_wall.centroid_mm[2] - 2.).abs() < 1e-12);
    assert!(tool_wall.centroid_mm[0].abs() < 1e-9);
}

#[test]
fn offset_through_bore_in_a_stepped_shaft_hub() {
    let screw = shaft(&[(9., [0., 30.]), (5., [30., 540.])]);
    let cut = body(subtract(
        &screw,
        tool([0., -40., 15.], [0., 1., 0.], 2.75, 80.),
        None,
    ));
    // Offset of the tool axis is along +x for travel +y.
    let shifted = body(subtract(
        &screw,
        tool([3., -40., 15.], [0., 1., 0.], 2.75, 80.),
        None,
    ));
    assert_eq!(cut.wall().1, 9.);
    let slices = removed_by_slices(9., 3., 2.75, [-40., 40.]);
    assert!((shifted.removed_mm3() - slices).abs() < 1e-6 * slices);
    assert!(shifted.removed_mm3() < cut.removed_mm3());
    let full = PI * (81. * 30. + 25. * 510.);
    assert!((shifted.volume_mm3() - (full - shifted.removed_mm3())).abs() < 1e-9);
    // The far hole is the mirror of the near one across the offset plane.
    let [near, far] = shifted.holes.map(Option::unwrap);
    assert!((near.area_mm2 - far.area_mm2).abs() < 1e-12 * near.area_mm2);
    assert!((near.centroid_mm[1] + far.centroid_mm[1]).abs() < 1e-9);
    let mesh = shifted.mesh(options()).unwrap();
    assert!(closed(&mesh));
    let error = (signed_volume(&mesh) - shifted.volume_mm3()).abs();
    assert!(error < 2e-3 * shifted.volume_mm3(), "{error}");
}

#[test]
fn blind_bore_keeps_a_planar_bottom() {
    let hub = shaft(&[(9., [-15., 15.])]);
    // Enters from -x and stops 1 mm past the axis.
    let cut = body(subtract(
        &hub,
        tool([-20., 0., 0.], [1., 0., 0.], 2.75, 21.),
        None,
    ));
    assert_eq!(cut.ends(), [End::Through, End::Cap]);
    let slices = removed_by_slices(9., 0., 2.75, [-20., 1.]);
    assert!((cut.removed_mm3() - slices).abs() < 1e-6 * slices);
    let (center, normal, radius) = cut.cap(1).unwrap();
    assert_eq!(
        (center, normal, radius),
        ([1., 0., 0.], [-1., 0., 0.], 2.75)
    );
    assert!(cut.cap(0).is_none() && cut.holes[1].is_none());
    // Reversed travel: the same bore entered from +x.
    let reversed = body(subtract(
        &hub,
        tool([20., 0., 0.], [-1., 0., 0.], 2.75, 21.),
        None,
    ));
    assert_eq!(reversed.ends(), [End::Through, End::Cap]);
    assert!((reversed.removed_mm3() - cut.removed_mm3()).abs() < 1e-10);
    let mesh = cut.mesh(options()).unwrap();
    assert!(closed(&mesh));
    assert_eq!(
        *mesh.face_indices.iter().max().unwrap(),
        cut.shaft_faces() + 1
    );
    assert!((signed_volume(&mesh) - cut.volume_mm3()).abs() < 2e-3 * cut.volume_mm3());
    // Material behind the bottom stays, the bore itself is empty.
    assert!(cut.contains([1.5, 0., 0.], 1e-9));
    assert!(!cut.contains([0.5, 0., 0.], 1e-9));
    assert!(!cut.contains([-8., 0., 0.], 1e-9));
    assert!(cut.contains([-8., 0., 2.9], 1e-9));
    assert!(cut.contains([1., 0., 0.], 1e-6));
}

#[test]
fn misses_and_contacts_leave_the_shaft_unchanged() {
    let hub = shaft(&[(9., [-15., 15.])]);
    let unchanged = |t: Tool| matches!(subtract(&hub, t, None), Ok(Cut::Unchanged));
    // Radially apart, beyond an end, before or after the shaft along travel.
    assert!(unchanged(tool([-50., 20., 0.], [1., 0., 0.], 2.75, 100.)));
    assert!(unchanged(tool([-50., 0., 30.], [1., 0., 0.], 2.75, 100.)));
    assert!(unchanged(tool([-50., 0., 0.], [1., 0., 0.], 2.75, 30.)));
    assert!(unchanged(tool([12., 0., 0.], [1., 0., 0.], 2.75, 30.)));
    // Exact contacts: tangent from outside, flush with an end plane, a cap
    // on the outer tangent plane.
    assert!(unchanged(tool([-50., 11.75, 0.], [1., 0., 0.], 2.75, 100.)));
    assert!(unchanged(tool([-50., 0., 17.75], [1., 0., 0.], 2.75, 100.)));
    assert!(unchanged(tool([-50., 0., 0.], [1., 0., 0.], 2.75, 41.)));
    assert!(unchanged(tool([9., 0., 0.], [1., 0., 0.], 2.75, 41.)));
}

#[test]
fn tangency_and_unsupported_configurations_are_refused() {
    let hub = shaft(&[(9., [-15., 15.])]);
    let run = |t: Tool| subtract(&hub, t, None).map(|_| ()).unwrap_err();
    // Internal tangency of the tool wall to the shaft wall.
    assert_eq!(
        run(tool([-50., 6.25, 0.], [1., 0., 0.], 2.75, 100.)),
        Error::Unsupported
    );
    // A chordal cut, and a section crossing an end of the bored wall.
    assert_eq!(
        run(tool([-50., 8., 0.], [1., 0., 0.], 2.75, 100.)),
        Error::Unsupported
    );
    assert_eq!(
        run(tool([-50., 0., 13.], [1., 0., 0.], 2.75, 100.)),
        Error::Unsupported
    );
    // A blind bottom that crosses the shaft wall, and an enclosed void.
    assert_eq!(
        run(tool([-50., 0., 0.], [1., 0., 0.], 2.75, 58.8)),
        Error::Unsupported
    );
    assert_eq!(
        run(tool([-2., 0., 0.], [1., 0., 0.], 2.75, 4.)),
        Error::Unsupported
    );
    // Within the margin of tangency and not literal: unresolved.
    let mut near = tool([-50., 6.25 + 1e-12, 0.], [1., 0., 0.], 2.75, 100.);
    near.literal = false;
    assert_eq!(run(near), Error::Unresolved);
    // A tool that is not perpendicular to the shaft is invalid here.
    assert_eq!(
        run(tool([-50., 0., 0.], [0.6, 0., 0.8], 2.75, 100.)),
        Error::Invalid
    );
    // A hollow or partial revolution is not a stepped solid shaft.
    let tube = Revolution::new(vec![vec![[4., 0.], [9., 0.], [9., 10.], [4., 10.]]], TAU).unwrap();
    assert!(stepped_walls(&tube).is_none());
}

#[test]
fn shoulder_band_crossing_is_refused_but_a_far_tool_misses() {
    let screw = shaft(&[(9., [0., 30.]), (5., [30., 540.])]);
    let run = |t: Tool| subtract(&screw, t, None);
    assert_eq!(
        run(tool([0., -40., 29.], [0., 1., 0.], 2.75, 80.)).unwrap_err(),
        Error::Unsupported
    );
    // A bore in the narrow shank is cut in that wall.
    let shank = body(run(tool([0., -40., 100.], [0., 1., 0.], 1.5, 80.)));
    assert_eq!(shank.wall().1, 5.);
    assert!(matches!(
        run(tool([0., -40., 30.], [0., 1., 0.], 2.75, 30.)),
        Ok(Cut::Unchanged)
    ));
}

#[test]
fn interrupted_cut_is_cancelled() {
    let hub = shaft(&[(9., [-15., 15.])]);
    let stop = || true;
    assert_eq!(
        subtract(
            &hub,
            tool([-50., 0., 0.], [1., 0., 0.], 2.75, 100.),
            Some(&stop)
        )
        .unwrap_err(),
        Error::Cancelled
    );
}

#[test]
fn rims_lie_on_both_supports() {
    let hub = shaft(&[(9., [-15., 15.])]);
    let cut = body(subtract(
        &hub,
        tool([-50., 2., 1.], [1., 0., 0.], 2.75, 100.),
        None,
    ));
    for side in 0..2 {
        for point in cut.sample_rim(side, options()).unwrap() {
            assert!((point[0].portable_hypot(point[1]) - 9.).abs() < 1e-12);
            assert!(((point[1] - 2.).portable_hypot(point[2] - 1.) - 2.75).abs() < 1e-12);
        }
    }
    let [low, high] = cut.extent([1., 0., 0.]);
    assert_eq!([low, high], [-9., 9.]);
}
