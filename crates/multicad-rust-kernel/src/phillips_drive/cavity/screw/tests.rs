use super::*;
use crate::{head_profile::Head, helical_thread};
pub(super) fn fixture(index: usize, plain: f64, left: bool) -> PhillipsScrew {
    let cavity = super::super::tests::cavity(index);
    let (rs, pitch, head) = match index {
        0 => (
            1.,
            0.4,
            Head::Crowned {
                radius_mm: 2.,
                height_mm: 1.6,
                crown_radius_mm: 3.2,
                top_land_radius_mm: 1.25,
                neck_radius_mm: 0.1,
            },
        ),
        1 => (
            1.5,
            0.5,
            Head::Crowned {
                radius_mm: 2.8,
                height_mm: 2.4,
                crown_radius_mm: 3.,
                top_land_radius_mm: 2.1,
                neck_radius_mm: 0.1,
            },
        ),
        _ => (
            1.5,
            0.5,
            Head::Countersunk {
                radius_mm: 2.75,
                height_mm: 1.65,
                included_angle_deg: 90.,
                edge_chamfer_mm: 0.05,
                neck_radius_mm: 0.1,
            },
        ),
    };
    let thread = Thread::new(helical_thread::Spec {
        major_radius_mm: rs - 0.03,
        pitch_mm: pitch,
        length_mm: 8.,
        depth_mm: 0.6 * pitch,
        crest_width_mm: 0.125 * pitch,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts: 1,
        left_hand: left,
        root_radius_mm: None,
    })
    .unwrap();
    let profile = Profile::new_for_thread_union(head, rs, plain).unwrap();
    PhillipsScrew::new(
        thread,
        profile,
        SourceSpec {
            head: Some(head),
            shank_length_mm: plain,
            minimum_wall_mm: 1e-5,
        },
        cavity,
        100_000,
        &|| false,
    )
    .unwrap()
}
#[test]
fn real_cavities_keep_original_thread_profile_and_exact_side_face_ownership() {
    for index in 0..3 {
        for plain in [0., 0.2, 1.] {
            for left in [false, true] {
                let s = fixture(index, plain, left);
                let t = *s.thread();
                let p = *s.profile();
                let spec = s.spec();
                assert_eq!(t.spec().left_hand, left);
                assert_eq!(spec.shank_length_mm, plain);
                assert_eq!(s.head_top_z(), p.top()[1] - plain);
                assert_eq!(
                    s.visible_thread_interval_mm(),
                    [p.end()[1] - plain, t.spec().length_mm]
                );
                assert!(s.nominal_volume_mm3() > 0.);
                let faces = s.faces();
                assert!(!faces.contains(&ScrewFace::Cavity(FaceKey::Mouth)));
                let side = s.cavity().side_face_keys();
                assert_eq!(side.len(), 28);
                for key in side {
                    assert_eq!(
                        faces
                            .iter()
                            .filter(|&&f| f == ScrewFace::Cavity(*key))
                            .count(),
                        1
                    );
                }
                assert!(s.host_envelope().minimum_radial_clearance_mm() > spec.minimum_wall_mm);
                let clone = s.clone();
                assert_eq!(clone.thread(), &t);
                assert_eq!(clone.profile(), &p);
                assert_eq!(clone.spec(), spec);
                assert_eq!(clone.faces(), faces);
                // The source cavity still owns real H0 inner cylinders or V planes.
                assert_eq!(clone.cavity().source(), s.cavity().source());
            }
        }
    }
}
#[test]
fn wrong_union_source_and_host_breach_remain_closed() {
    let s = fixture(2, 0., false);
    let mut spec = s.spec();
    spec.shank_length_mm = 1.;
    assert!(matches!(
        PhillipsScrew::new(
            *s.thread(),
            *s.profile(),
            spec,
            s.cavity().clone(),
            100_000,
            &|| false
        ),
        Err(Error::Invalid)
    ));
    spec = s.spec();
    spec.minimum_wall_mm = 1.;
    assert!(matches!(
        PhillipsScrew::new(
            *s.thread(),
            *s.profile(),
            spec,
            s.cavity().clone(),
            100_000,
            &|| false
        ),
        Err(Error::Unresolved)
    ));
    assert!(matches!(
        PhillipsScrew::new(
            *s.thread(),
            *s.profile(),
            s.spec(),
            s.cavity().clone(),
            100_000,
            &|| true
        ),
        Err(Error::Cancelled)
    ));
}
