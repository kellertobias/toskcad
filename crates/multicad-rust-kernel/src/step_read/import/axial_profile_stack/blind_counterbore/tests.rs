use super::*;
#[test]
fn blind_counterbore_recognition_rejects_stale_or_unearned_literal_geometry() {
    let stack = ProfileStack::blind_counterbore(
        RectangularBlindCounterbore {
            min_mm: [-10.; 2],
            max_mm: [10.; 2],
            height_mm: 20.,
            bore_radius_mm: 2.,
            bore_depth_mm: 16.,
            seat_radius_mm: 6.,
            seat_depth_mm: 1.,
            minimum_wall_mm: 1.,
            minimum_floor_mm: 1.,
        },
        1e-7,
        1e-5,
    )
    .unwrap();
    let text = crate::step::write_axial_profile_stack(&stack, [0.; 3], [1., 0., 0.], [0., 1., 0.])
        .unwrap();
    let entities = parse_entities(&text, 1e-7).unwrap();
    let mut shell = Shell::new(&entities, 1e-7, true);
    let faces = shell.validated_faces().unwrap();
    let mut read = shell.readback(&faces).unwrap();
    let frame = Frame {
        origin: [0.; 3],
        axes: [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
    };
    assert!(shell.rectangular_blind_counterbore(&faces, &read).is_some());
    read.volume_mm3 += 1.;
    assert!(shell.rectangular_blind_counterbore(&faces, &read).is_none());
    read.volume_mm3 -= 1.;
    let arcs = std::mem::take(&mut read.arcs);
    assert!(shell.rectangular_blind_counterbore(&faces, &read).is_none());
    assert!(
        shell
            .match_blind_counterbore(&faces, &read, &frame)
            .is_none()
    );
    read.arcs = arcs;
    let point_id = *shell.points.keys().next().unwrap();
    let original = shell.points[&point_id];
    shell.points.get_mut(&point_id).unwrap()[2] = 5.;
    assert!(
        shell
            .match_blind_counterbore(&faces, &read, &frame)
            .is_none()
    );
    shell.points.insert(point_id, original);
    let arc = read.arcs[0];
    for a in &mut read.arcs {
        a.radius_mm = 100.;
    }
    assert!(
        shell
            .match_blind_counterbore(&faces, &read, &frame)
            .is_none()
    );
    read.arcs[0] = arc;
    read = shell.readback(&faces).unwrap();
    shell.points.remove(&point_id);
    assert!(shell.rectangular_blind_counterbore(&faces, &read).is_none());
    shell.points.insert(point_id, original);
    // A permuted/stale face table must not be accepted from bounds/mass.
    let mut permuted = faces;
    permuted.rotate_left(9);
    assert!(
        shell
            .rectangular_blind_counterbore(&permuted, &read)
            .is_none()
    );
    let mut caps = permuted;
    for face in &mut caps {
        if let Face::Planar { normal, .. } = face {
            *normal = [0., 0., 1.];
        }
    }
    assert!(shell.rectangular_blind_counterbore(&caps, &read).is_none());
}
