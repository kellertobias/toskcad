use super::*;
use crate::{
    head_profile::{Head, Profile, Segment},
    profiled_screw::{SourceSpec, Spec},
};

pub(super) fn thread(starts: u32) -> Thread {
    Thread::new(crate::helical_thread::Spec {
        major_radius_mm: 1.5,
        pitch_mm: 0.5,
        length_mm: 0.5,
        depth_mm: 0.2,
        crest_width_mm: 0.0625,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts,
        left_hand: false,
        root_radius_mm: None,
    })
    .unwrap()
}
pub(super) fn spec(chamfer: f64) -> Spec {
    Spec {
        head: Head::Countersunk {
            radius_mm: 3.,
            height_mm: 1.65,
            included_angle_deg: 90.,
            edge_chamfer_mm: chamfer,
            neck_radius_mm: 0.1,
        },
        shank_length_mm: 1.,
        socket_across_flats_mm: 2.,
        socket_depth_mm: 0.6,
        minimum_wall_mm: 0.1,
    }
}
fn corners(file: &mut StepFile, thread: &Thread, p: &Placement) -> Vec<[Vertex; 2]> {
    let [_, first, _, _, last, _] = thread.profile_breaks();
    (0..thread.spec().starts)
        .map(|k| {
            [first, last].map(|phase| {
                let angle =
                    -(phase + f64::from(k) * thread.spec().pitch_mm) / thread.helical_rate();
                let (s, c) = angle.sin_cos();
                file.vertex(p.point_f64(
                    [
                        thread.spec().major_radius_mm * c,
                        thread.spec().major_radius_mm * s,
                    ],
                    0.,
                ))
            })
        })
        .collect()
}

#[test]
fn unresolved_positive_head_dimensions_and_distorted_frames_refuse() {
    let small = Screw::new(thread(1), spec(1e-9)).unwrap();
    assert_eq!(admission(&small, [0.; 3], 1e-7), Err(StepError::Limit));
    // Positive endpoint coupling residuals cannot make a sub-resolution
    // literal support radius pass merely because its trim spans round up.
    let r = 1e-6;
    let e = 2e-14;
    let c = 1.5 + r + e;
    let profile = Profile::from_retained_segments(&[Segment {
        start: [c, -r - e],
        end: [1.5, 0.],
        curve: Curve::Arc {
            center: [c, 0.],
            radius_mm: r,
            sweep: -std::f64::consts::FRAC_PI_2,
        },
        role: Role::Neck,
    }])
    .unwrap();
    let small = Screw::from_source(
        thread(1),
        profile,
        SourceSpec {
            head: None,
            shank_length_mm: 0.,
            minimum_wall_mm: 0.1,
        },
        crate::screw_recess::Recess::canonical_flat_hex(1., 2e-7).unwrap(),
    )
    .unwrap();
    assert_eq!(admission(&small, [0.; 3], 1e-7), Err(StepError::Limit));
    let body = Screw::new(thread(1), spec(0.1)).unwrap();
    let frame = Placement::new([0.; 3], [1., 0., 0.], [0., 1. + 1e-12, 0.]).unwrap();
    assert_eq!(
        placement(&body, &frame, 1e-7),
        Err(StepError::InvalidPlacement)
    );
    let close = Screw::new_with_shank_radius(thread(1), 1.5 + 1e-8, spec(0.1)).unwrap();
    assert_eq!(admission(&close, [0.; 3], 1e-7), Err(StepError::Limit));
    let profile = Profile::from_retained_segments(&[Segment {
        start: [2., -2.],
        end: [2., 0.],
        curve: Curve::Line,
        role: Role::Shank,
    }])
    .unwrap();
    let cylinder = Screw::from_source(
        thread(1),
        profile,
        SourceSpec {
            head: None,
            shank_length_mm: 0.,
            minimum_wall_mm: 0.1,
        },
        crate::screw_recess::Recess::canonical_flat_hex(1., 0.6).unwrap(),
    )
    .unwrap();
    assert_eq!(
        admission(&cylinder, [0.; 3], 1e-7),
        Err(StepError::Unsupported)
    );
}

#[test]
fn physical_junction_and_drive_encoders_refuse_collapsed_geometry() {
    let p = Placement::new([0.; 3], [1., 0., 0.], [0., 1., 0.]).unwrap();
    for second in [false, true] {
        let thread = thread(2);
        let mut file = StepFile::new();
        let mut rim = corners(&mut file, &thread, &p);
        let index = usize::from(!second);
        let mut displaced = rim[0][index].point;
        displaced[2] += 1.;
        rim[0][index] = file.vertex(displaced);
        assert!(
            junction(
                &mut file,
                &p,
                &thread,
                thread.spec().major_radius_mm,
                0.,
                &rim
            )
            .is_err()
        );
    }
    let thread = thread(1);
    let profile = Profile::from_retained_segments(&[Segment {
        start: [3., -1.5],
        end: [1.5, 0.],
        curve: Curve::Line,
        role: Role::Cone,
    }])
    .unwrap();
    let source = Screw::from_source(
        thread,
        profile,
        SourceSpec {
            head: None,
            shank_length_mm: 0.,
            minimum_wall_mm: 0.1,
        },
        crate::screw_recess::Recess::canonical_flat_hex(1.5, 0.6).unwrap(),
    )
    .unwrap();
    let mut file = StepFile::new();
    let rim = corners(&mut file, &thread, &p);
    let join = junction(
        &mut file,
        &p,
        &thread,
        thread.spec().major_radius_mm,
        0.,
        &rim,
    )
    .unwrap();
    let collapsed = Placement {
        normal: [0.; 3],
        ..p
    };
    assert!(head(&mut file, &collapsed, (&source).into(), join, 1e-7, None).is_err());
    let mut file = StepFile::new();
    let rim = corners(&mut file, &thread, &p);
    assert!(junction(&mut file, &p, &thread, f64::INFINITY, 0., &rim).is_err());
}

#[test]
fn far_world_collapsed_crest_refuses_public_export_and_writer_retries_cleanly() {
    let tiny = Thread::new(crate::helical_thread::Spec {
        crest_width_mm: 1e-14,
        ..thread(1).spec()
    })
    .unwrap();
    let body = Screw::new(tiny, spec(0.1)).unwrap();
    let mut writer = crate::step::Writer::new();
    let failed =
        writer.append_profiled_screw(&body, [100000.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7);
    assert_eq!(failed, Err(StepError::InvalidGeometry));
    let qualified = Screw::new(thread(1), spec(0.1)).unwrap();
    writer
        .append_profiled_screw(&qualified, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7)
        .unwrap();
    let text = writer.finish().unwrap();
    assert_eq!(text.matches("MANIFOLD_SOLID_BREP(").count(), 1);
    crate::step_read::read_import_solid(&text, 1e-7).unwrap();
}

#[test]
#[cfg(all(target_os = "macos", target_arch = "aarch64"))]
fn shared_step_outer_preserves_frozen_existing_head_file_bytes() {
    use crate::{
        hex_drive::HexDrive,
        hexalobular_drive::{Cone, Drive, TangentSpec},
        screw_recess::Recess,
    };
    let mut fingerprints = Vec::new();
    for family in 0..3 {
        for union in [false, true] {
            let mut ts = thread(1).spec();
            ts.major_radius_mm = 1.47;
            ts.left_hand = union;
            let t = Thread::new(ts).unwrap();
            let h = if family == 0 {
                spec(0.1).head
            } else {
                Head::Crowned {
                    radius_mm: 2.8,
                    height_mm: 2.4,
                    crown_radius_mm: 3.,
                    top_land_radius_mm: 2.1,
                    neck_radius_mm: 0.1,
                }
            };
            let drive = match family {
                0 => Recess::canonical_flat_hex(2., 0.6).unwrap(),
                1 => {
                    let angle = 0.193_f64;
                    let canonical = HexDrive::canonical(2., 1e-12).unwrap();
                    let corners = canonical.corners_mm().map(|[x, y]| {
                        [
                            x * angle.cos() - y * angle.sin(),
                            x * angle.sin() + y * angle.cos(),
                        ]
                    });
                    Recess::flat_hex(
                        HexDrive::from_retained_corners(corners, 2., 1e-12).unwrap(),
                        0.6,
                    )
                    .unwrap()
                }
                _ => Recess::conical_hexalobular(
                    Drive::from_authored_tangent(
                        TangentSpec {
                            major_diameter_mm: 2.8,
                            minor_diameter_mm: 2.05,
                            outer_lobe_radius_mm: 0.28,
                        },
                        Cone {
                            tip_depth_mm: 1.9,
                            depth_per_radius: 1. / 3f64.sqrt(),
                        },
                    )
                    .unwrap(),
                    1e-7,
                    100_000,
                )
                .unwrap(),
            };
            let s = Screw::new_with_recess(t, 1.5, h, if union { 0. } else { 1. }, 0.01, drive)
                .unwrap();
            let (o, u, v) = if union {
                ([105., -20., 60.], [0.8, 0.6, 0.], [0.36, -0.48, -0.8])
            } else {
                ([0.; 3], [1., 0., 0.], [0., 1., 0.])
            };
            let text = crate::step::write_profiled_screw(&s, o, u, v, 1e-7).unwrap();
            let mut hash = 0xcbf29ce484222325u64;
            for byte in text.bytes() {
                hash = (hash ^ u64::from(byte)).wrapping_mul(0x100000001b3);
            }
            fingerprints.push((text.len(), hash));
        }
    }
    // Captured from the unrefactored source at291514ca with the pinned
    // arm64 macOS toolchain. Cross-platform behavior has independent geometry
    // readback controls; platform transcendental rounding need not be equal.
    assert_eq!(
        fingerprints,
        vec![
            (212947, 12384164169844171353),
            (268920, 10569167352034862003),
            (213716, 17712331813652097306),
            (269663, 17769957927807670706),
            (477654, 5975792124410495725),
            (528452, 12947489326318806268),
        ]
    );
}
