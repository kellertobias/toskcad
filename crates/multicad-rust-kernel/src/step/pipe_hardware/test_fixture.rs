//! Real authored/owned fixture shared by public-file and independent helper controls.
use crate::axial_pipe_thread::authored::{
    AuthoredSection, AuthoredSpec, AuthoredThread,
    hardware::shell::{self, occupied},
};
pub(super) fn engineering() -> occupied::OccupiedHead {
    engineering_height(3.)
}
pub(super) fn engineering_height(height: f64) -> occupied::OccupiedHead {
    engineering_dimensions(6., height)
}
pub(super) fn engineering_dimensions(length: f64, height: f64) -> occupied::OccupiedHead {
    // Engineering NPT-flat law, NOT a standard/Core preset acceptance tuple.
    // Canonical Core external NPT forwards an active crest limit and remains
    // outside this deliberately unclipped, fixed-frame first source family.
    let thread = AuthoredThread::new(
        AuthoredSpec {
            major_radius_mm: 5.,
            pitch_mm: 1.,
            length_mm: length,
            radial_taper: 1. / 32.,
            section: AuthoredSection::NptFlat {
                depth_mm: 0.6,
                flat_width_mm: 0.2,
            },
            crest_radius_limit_mm: None,
        },
        1e-7,
        1e-6,
    )
    .unwrap();
    let boundary = thread
        .raw_boundary_ownership(1e-7, None)
        .unwrap()
        .ownership
        .bind_surfaces(1e-7, None)
        .unwrap()
        .certify_axial(1e-7, 1e-6, None)
        .unwrap()
        .certify_caps(None)
        .unwrap()
        .certify_ordering(None)
        .unwrap()
        .certify_material(1e-7, None)
        .unwrap();
    occupied::certify(
        shell::external_npt_head(&thread, &boundary, 16., height, None).unwrap(),
        1e-7,
        1e-6,
        None,
    )
    .unwrap()
}
