pub(in crate::step_read::phillips_cutter::physical::body::derived_thread) fn actual_source()
-> crate::phillips_drive::PhillipsScrew {
    use crate::{
        head_profile::{Head, Profile},
        helical_thread::{Spec as ThreadSpec, Thread},
        phillips_drive::{Cavity, Cutter, Domain, Inner, PhillipsScrew, Spec},
        profiled_screw::SourceSpec,
    };
    // Literal packet captured from the actual ISO7045 M2 normal factory. It
    // creates the file fixture; none of these values enters the decoder proof.
    let head = Head::Crowned {
        radius_mm: 2.,
        height_mm: 1.6,
        crown_radius_mm: 3.2,
        top_land_radius_mm: 1.235645633015642,
        neck_radius_mm: 0.1,
    };
    let cavity = Cavity::new(
        Cutter::new(
            Spec {
                mouth_diameter_mm: 1.9,
                mouth_radius_mm: 0.3,
                b_width_mm: 0.61,
                g_diameter_mm: 0.81,
                tilt_angle_deg: 7.,
                outer_flank_angle_deg: 26.5,
                floor_angle_deg: 28.,
                outer_corner_angle_deg: 92.,
                inner: Inner::Rounded { radius_mm: 0.25 },
            },
            Domain {
                corner_start_depth_mm: -3.8,
                corner_end_depth_mm: 3.208443210824065,
                corner_reach_mm: 16.,
            },
        )
        .unwrap(),
        1e-5,
        100_000,
    )
    .unwrap();
    PhillipsScrew::new(
        Thread::new(ThreadSpec {
            major_radius_mm: 1.,
            pitch_mm: 0.4,
            length_mm: 11.5,
            depth_mm: 0.21650635094610965,
            crest_width_mm: 0.05,
            lower_flank_deg: 30.,
            upper_flank_deg: 30.,
            starts: 1,
            left_hand: false,
            root_radius_mm: None,
        })
        .unwrap(),
        Profile::new_for_thread_union(head, 1., 0.5).unwrap(),
        SourceSpec {
            head: Some(head),
            shank_length_mm: 0.5,
            minimum_wall_mm: 1e-5,
        },
        cavity,
        100_000,
        &|| false,
    )
    .unwrap()
}
pub(super) fn actual_text() -> String {
    let source = actual_source();
    crate::step::write_phillips_screw(&source, [0., 0., 0.5], [1., 0., 0.], [0., 1., 0.], 1e-7)
        .unwrap()
}
