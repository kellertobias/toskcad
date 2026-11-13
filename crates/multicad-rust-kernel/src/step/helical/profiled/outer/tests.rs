use super::*;
use crate::{
    head_profile::Head,
    phillips_drive::{Cavity, Cutter},
    profiled_screw::{SourceSpec, Spec},
};
#[test]
fn borrowed_outer_sources_keep_actual_profile_original_thread_and_visible_phase() {
    let mut spec = super::super::tests::spec(0.1);
    spec.head = Head::Crowned {
        radius_mm: 2.8,
        height_mm: 2.4,
        crown_radius_mm: 3.,
        top_land_radius_mm: 2.1,
        neck_radius_mm: 0.1,
    };
    let old = Screw::new(super::super::tests::thread(2), spec).unwrap();
    let (s, d) = crate::phillips_drive::test_fixture(1);
    let phi = PhillipsScrew::new(
        *old.thread(),
        *old.profile(),
        SourceSpec::from(spec),
        Cavity::new(Cutter::new(s, d).unwrap(), 1e-5, 100_000).unwrap(),
        100_000,
        &|| false,
    )
    .unwrap();
    let data = Data::from(&phi);
    assert!(std::ptr::eq(data.profile, phi.profile()));
    assert!(std::ptr::eq(data.thread, phi.thread()));
    assert_eq!(data.plain, phi.spec().shank_length_mm);
    assert_eq!(data.visible, phi.visible_thread_interval_mm());
    assert_eq!(data.top_z(), phi.head_top_z());
    let legacy = Data::from(&old);
    assert_eq!(legacy.profile, data.profile);
    assert_eq!(legacy.thread, data.thread);
    let p = Placement::new([105., -20., 60.], [0.8, 0.6, 0.], [0.36, -0.48, -0.8]).unwrap();
    admission(data, p.origin, 1e-7).unwrap();
    placement(data, &p, 1e-7).unwrap();
    // Head authoring hints are absent from this geometry seam altogether.
    let _: Spec = spec;
}
