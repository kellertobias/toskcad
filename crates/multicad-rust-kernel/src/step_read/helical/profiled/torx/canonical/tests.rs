use super::*;
use crate::{
    head_profile::Head,
    helical_thread::{Spec, Thread},
    profiled_screw::{Screw, Spec as ScrewSpec},
};
fn text(left: bool, plain: f64, pose: [Point; 3]) -> String {
    let thread = Thread::new(Spec {
        major_radius_mm: 1.5,
        pitch_mm: 0.5,
        length_mm: 3.,
        depth_mm: 0.2,
        crest_width_mm: 0.0625,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts: if left { 2 } else { 1 },
        left_hand: left,
        root_radius_mm: None,
    })
    .unwrap();
    let screw = Screw::new_with_shank_thread_union(
        thread,
        1.5,
        ScrewSpec {
            head: Head::Crowned {
                radius_mm: 2.8,
                height_mm: 2.4,
                crown_radius_mm: 5.,
                top_land_radius_mm: 2.05,
                neck_radius_mm: 0.125,
            },
            shank_length_mm: plain,
            socket_across_flats_mm: 2.,
            socket_depth_mm: 1.,
            minimum_wall_mm: 0.1,
        },
    )
    .unwrap();
    crate::step::write_profiled_screw(&screw, pose[0], pose[1], pose[2], 1e-7).unwrap()
}
fn shell<T>(text: &str, f: impl FnOnce(&mut Shell<'_, '_>, Vec<u32>) -> T) -> T {
    let entities = crate::step_read::parse_entities(text, 1e-7).unwrap();
    let faces = entities
        .values()
        .flatten()
        .find(|e| e.name == "CLOSED_SHELL")
        .map(|e| refs(&e.args[1]).unwrap())
        .unwrap_or_default();
    f(&mut Shell::new(&entities, 1e-7, true), faces)
}
fn axis(shell: &Shell<'_, '_>, faces: &[u32], count: usize) -> u32 {
    let (s, _, _) = face(shell, faces[count]).unwrap();
    reference(&shell.args(s, "PLANE").unwrap()[1]).unwrap()
}
#[test]
fn actual_full_support_datums_and_inverse_transports_follow_both_hands_and_proper_poses() {
    for left in [false, true] {
        for plain in [0., 1.] {
            for pose in [
                [[0.; 3], [1., 0., 0.], [0., 1., 0.]],
                [[30., -20., 10.], [0.8, 0.6, 0.], [-0.36, 0.48, 0.8]],
                [
                    [10_000., -20_000., 30_000.],
                    [0.8, 0.6, 0.],
                    [-0.36, 0.48, 0.8],
                ],
            ] {
                let text = text(left, plain, pose);
                shell(&text, |s, faces| {
                    let frame = Frame::new(pose).unwrap();
                    let count = if left { 8 } else { 4 };
                    let proof = Canonical::new(s, &faces, frame, count).unwrap();
                    let id = axis(s, &faces, count);
                    proof.parallel_axis(s, id).unwrap();
                    proof.coplanar_origins(s, id, id).unwrap();
                    proof.same_axis_line(s, id, id).unwrap();
                    proof.coaxial_with_thread(s, id).unwrap();
                    let local = proof.local_point(pose[0], true).unwrap();
                    assert!(local.into_iter().all(|i| i.lo <= 0. && i.hi >= 0.));
                    let v = proof.local_vector(pose[1], true).unwrap();
                    assert!(v[0].lo <= 1. && v[0].hi >= 1.);
                    let vector = proof.project_vector(pose[2].map(I::point), false).unwrap();
                    assert!(vector.into_iter().all(I::valid));
                    let error = proof.transport(s, &faces[count..]).unwrap();
                    // Independent composed finite support points exercise the
                    // joint origin AND normalized-basis transport reserve.
                    let axis_float = s.placement(id).unwrap();
                    for conversion in [false, true] {
                        let axis_bounds = Axis::read(s, id, conversion).unwrap();
                        for parameters in [[0.2, 0.4, 1.], [1., -0.5, -0.75]] {
                            let world = std::array::from_fn(|k| {
                                axis_bounds.origin[k]
                                    + (0..3).fold(I::point(0.), |sum, j| {
                                        sum + axis_bounds.axes[j][k] * I::point(parameters[j])
                                    })
                            });
                            let expanded = proof.project_point(world, conversion).unwrap();
                            let local_axis = proof.local_axis(s, id, conversion).unwrap();
                            let bounds = std::array::from_fn(|k| {
                                local_axis.origin[k]
                                    + (0..3).fold(I::point(0.), |sum, j| {
                                        sum + local_axis.axes[j][k] * I::point(parameters[j])
                                    })
                            });
                            // Both enclose the same point; expanding in world
                            // coordinates adds dependency/rounding slack that
                            // is not itself a physical displacement.
                            for k in 0..3 {
                                assert!(
                                    expanded[k].lo <= bounds[k].hi
                                        && bounds[k].lo <= expanded[k].hi
                                );
                            }
                            let axes = [axis_float.u, axis_float.v, axis_float.normal]
                                .map(|v| frame.rows.map(|row| dot(row, v) / frame.determinant));
                            let origin = frame.local(axis_float.origin);
                            let actual = std::array::from_fn(|k| {
                                origin[k] + (0..3).map(|j| axes[j][k] * parameters[j]).sum::<f64>()
                            });
                            let bound = if conversion {
                                error.numeric_conversion_mm
                            } else {
                                error.arithmetic_mm
                            };
                            assert!(
                                distance(actual, bounds) <= bound,
                                "joint support displacement must fit: delta={} bound={} conversion={} pose={:?} parameters={:?}",
                                distance(actual, bounds),
                                bound,
                                conversion,
                                pose,
                                parameters
                            );
                        }
                    }

                    assert!(
                        error.total_mm().is_finite()
                            && error.total_mm() > 0.
                            && error.total_mm() < 1e-6,
                        "{}",
                        error.total_mm()
                    );
                });
            }
        }
    }
}
#[test]
fn raw_support_offsets_and_tilts_cannot_borrow_caller_tolerance() {
    let pose = [[0.; 3], [1., 0., 0.], [0., 1., 0.]];
    let original = text(false, 0., pose);
    let extra = "#900000=CARTESIAN_POINT('',(0.,0.,0.000001));#900001=DIRECTION('',(0.,0.,1.));#900002=DIRECTION('',(1.,0.,0.));#900003=AXIS2_PLACEMENT_3D('',#900000,#900001,#900002);#900004=CARTESIAN_POINT('',(0.000001,0.,0.));#900005=AXIS2_PLACEMENT_3D('',#900004,#900001,#900002);#900006=DIRECTION('',(0.000001,0.,1.));#900007=AXIS2_PLACEMENT_3D('',#900004,#900006,#900002);";
    let mut text = original;
    let end = text.rfind("ENDSEC;").unwrap();
    text.insert_str(end, extra);
    shell(&text, |s, faces| {
        s.tolerance = 1.;
        assert!(Axis::read(s, 900003, true).is_ok());
        assert!(Axis::read(s, 900007, true).is_ok());
        let proof = Canonical::new(s, &faces, Frame::new(pose).unwrap(), 4).unwrap();
        let id = axis(s, &faces, 4);
        assert!(proof.coplanar_origins(s, 900003, id).is_err());
        assert!(proof.same_axis_line(s, 900005, id).is_err());
        assert!(proof.coaxial_with_thread(s, 900005).is_err());
        assert!(proof.parallel_axis(s, 900007).is_err());
        assert!(Canonical::new(s, &[], Frame::new(pose).unwrap(), 4).is_err());
        assert!(Canonical::new(s, &faces, Frame::new(pose).unwrap(), usize::MAX).is_err());
    });
}
#[test]
fn directed_axis_normalization_and_malformed_records_fail_without_a_geometry_tolerance() {
    let input = "ISO-10303-21;HEADER;FILE_DESCRIPTION(('test'),'2;1');FILE_NAME('','',(),(),'','','');FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));ENDSEC;DATA;#1=CARTESIAN_POINT('',(0.,0.,0.));#2=DIRECTION('',(0.,0.,1.E-320));#3=DIRECTION('',(1.E-320,0.,0.));#4=AXIS2_PLACEMENT_3D('',#1,#2,#3);#5=DIRECTION('',(0.,0.,0.));#6=AXIS2_PLACEMENT_3D('',#1,#5,#3);#7=AXIS2_PLACEMENT_3D('',#1,#2,#2);#8=DIRECTION('',(1.,0.));#9=AXIS2_PLACEMENT_3D('');#10=DIRECTION('');#11=PLANE('');#12=ADVANCED_FACE('',(),#11,.T.);#13=B_SPLINE_SURFACE_WITH_KNOTS('',9.,9.,.UNSPECIFIED.);#15=ADVANCED_FACE('',(),#13,.T.);#17=PLANE('',#4);#18=ADVANCED_FACE('',(),#17,.T.);#19=AXIS2_PLACEMENT_3D('',#1,#2,#5);ENDSEC;END-ISO-10303-21;";
    shell(input, |s, _| {
        let a = Axis::read(s, 4, false).unwrap();
        assert!(a.axes[2][2].lo <= 1. && a.axes[2][2].hi >= 1.);
        assert!(Axis::read(s, 6, false).is_err());
        assert!(Axis::read(s, 7, false).is_err());
        assert!(Axis::read(s, 9, false).is_err());
        assert!(Axis::read(s, 19, false).is_err());
        assert!(raw_direction(s, 8, false).is_err());
        assert!(raw_direction(s, 10, false).is_err());
        let frame = Frame::new([[0.; 3], [1., 0., 0.], [0., 1., 0.]]).unwrap();
        // Malformed serialized records fail before any source geometry exists.
        assert!(Canonical::new(s, &[15, 12], frame, 1).is_err());
        assert!(Canonical::new(s, &[15, 18], frame, 1).is_err());
        assert!(unit([I::point(f64::INFINITY); 3]).is_err());
        assert!(unit([I::new(-1., 1.); 3]).is_err());
    });
}
