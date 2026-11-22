use super::*;
fn tool(left: bool) -> crate::helical_thread::Thread {
    crate::helical_thread::Thread::new(crate::helical_thread::Spec {
        major_radius_mm: 3.,
        pitch_mm: 1.,
        length_mm: 6.6,
        depth_mm: 0.5,
        crest_width_mm: 0.125,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts: 1,
        left_hand: left,
        root_radius_mm: None,
    })
    .unwrap()
}
#[test]
fn steep_actual_cone_operator_uses_the_earned_full_height_domain_and_whole_readback() {
    for left in [false, true] {
        for length in [0.02, 0.005] {
            let body = Chamfer::from_material(tool(left), 2.6, length, None).unwrap();
            assert!(
                (body.source_thread().spec().major_radius_mm - body.entry_radius_mm()) / length
                    > 10.
            );
            let placement =
                Placement::new([17., -8., 13.], [0.8, 0.6, 0.], [-0.36, 0.48, 0.8]).unwrap();
            let staged = prepare(&body, &placement, 1e-7, None).unwrap();
            assert!(staged.groups.iter().any(|g| g.active));
            for group in &staged.groups {
                for trim in &group.trims {
                    assert!(
                        trim.height_poles
                            .into_iter()
                            .chain(trim.poles.map(|p| p[2]))
                            .all(|z| z >= 0. && z <= length)
                    );
                }
            }
            assert!(staged.points.values().all(|p| p[2] >= 0. && p[2] <= length));
            let output = crate::step::write_thread_entry_chamfer(
                &body,
                [17., -8., 13.],
                [0.8, 0.6, 0.],
                [-0.36, 0.48, 0.8],
                1e-7,
            )
            .unwrap();
            crate::step_read::verify_thread_entry_chamfer(
                &output,
                &body,
                [17., -8., 13.],
                [0.8, 0.6, 0.],
                [-0.36, 0.48, 0.8],
                1e-7,
            )
            .unwrap();
        }
    }
}
#[test]
fn near_right_angle_without_a_positive_whole_cosine_bound_refuses_before_emission() {
    let body = Chamfer::from_material(tool(false), 2.6, 1e-20, None).unwrap();
    assert_eq!(
        ((3. - 2.6) / body.entry_length_mm()).portable_atan(),
        std::f64::consts::FRAC_PI_2
    );
    assert!(matches!(cone_angle(&body, 1e-7), Err(Error::Unsupported)));
    assert!(matches!(
        crate::step::write_thread_entry_chamfer(&body, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7),
        Err(Error::Unsupported)
    ));
}
