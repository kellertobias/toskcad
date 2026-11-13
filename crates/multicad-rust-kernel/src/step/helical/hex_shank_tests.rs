use super::*;
use crate::{headed_screw::ScrewSpec, helical_thread::Spec, step_read};
fn source(left_hand: bool, rounded: bool, plain: f64) -> Screw {
    let root = rounded.then_some(1. / (4. * 3_f64.sqrt()));
    Screw::new(
        Thread::new(Spec {
            major_radius_mm: 3.,
            pitch_mm: 1.,
            length_mm: 3.7,
            depth_mm: root.map_or(5. * 3_f64.sqrt() / 16., |r| 0.875 * 3_f64.sqrt() / 2. - r),
            crest_width_mm: 0.125,
            lower_flank_deg: 30.,
            upper_flank_deg: 30.,
            starts: 1,
            left_hand,
            root_radius_mm: root,
        })
        .unwrap(),
        ScrewSpec {
            head: ScrewHead::ExternalHex {
                across_flats_mm: 10.,
                height_mm: 4.,
            },
            shank_length_mm: plain,
            shank_radius_mm: 3.,
        },
    )
    .unwrap()
}
#[test]
fn hex_plain_shank_step_retains_real_cylinder_seat_head_and_independent_physical_mass() {
    for left in [false, true] {
        for rounded in [false, true] {
            let screw = source(left, rounded, 2.3);
            for [o, u, v] in [
                [[0.; 3], [1., 0., 0.], [0., 1., 0.]],
                [[1000.3, -2000.7, 0.1], [0., 1., 0.], [0., 0., 1.]],
                [[2., -3., 4.], [0.8, 0.6, 0.], [-0.36, 0.48, 0.8]],
            ] {
                let text = write_helical_screw(&screw, o, u, v, 1e-7).unwrap();
                let volume = step_read::verify_helical_screw(&text, &screw, o, u, v, 1e-7).unwrap();
                assert!((volume - screw.volume_mm3()).abs() < 1e-6);
                let (physical, imported) = step_read::read_import_solid(&text, 1e-7).unwrap();
                assert_eq!(physical.solid_count, 1);
                assert!((physical.volume_mm3 - screw.volume_mm3()).abs() < 1e-6);
                assert!(step_read::helical_screw_import_matches(
                    &imported,
                    &screw,
                    [o, u, v],
                    1e-7
                ));
                let step_read::ImportProfile::HelicalScrew { screw: actual, .. } = imported.profile
                else {
                    panic!("must retain the actual head and plain shaft")
                };
                assert!((actual.shank_length_mm - 2.3).abs() < 1e-9);
                assert_eq!(actual.shank_radius_mm, 3.);
                assert!(
                    matches!(actual.head,ScrewHead::ExternalHex {across_flats_mm,height_mm} if (across_flats_mm-10.).abs()<1e-9 && (height_mm-4.).abs()<1e-9)
                );
                let wrong = source(left, rounded, 2.4);
                assert!(step_read::verify_helical_screw(&text, &wrong, o, u, v, 1e-7).is_err());
                assert!(!step_read::helical_screw_import_matches(
                    &imported,
                    &wrong,
                    [o, u, v],
                    1e-7
                ));
            }
        }
    }
}

#[test]
fn hex_plain_shank_source_limits_and_writer_rollback_keep_previous_members() {
    let screw = source(false, false, 2.3);
    let unequal = Screw::new(
        *screw.thread(),
        ScrewSpec {
            shank_radius_mm: 3.1,
            ..screw.spec()
        },
    )
    .unwrap();
    let tiny = source(false, false, 1e-9);
    let thin = Screw::new(
        *screw.thread(),
        ScrewSpec {
            head: ScrewHead::ExternalHex {
                across_flats_mm: 6.00000001,
                height_mm: 4.,
            },
            ..screw.spec()
        },
    )
    .unwrap();
    let low = Screw::new(
        *screw.thread(),
        ScrewSpec {
            head: ScrewHead::ExternalHex {
                across_flats_mm: 10.,
                height_mm: 1e-9,
            },
            ..screw.spec()
        },
    )
    .unwrap();
    let mut writer = Writer::new();
    let [u, v] = [[1., 0., 0.], [0., 1., 0.]];
    writer
        .append_helical_screw(&screw, [0.; 3], u, v, 1e-7)
        .unwrap();
    assert_eq!(
        writer.append_helical_screw(&unequal, [0.; 3], u, v, 1e-7),
        Err(StepError::Unsupported)
    );
    for invalid in [tiny, thin, low] {
        assert!(!supports_helical_screw(&invalid, [0.; 3], 1e-7));
        assert_eq!(
            writer.append_helical_screw(&invalid, [0.; 3], u, v, 1e-7),
            Err(StepError::Limit)
        );
    }
    writer
        .append_helical_screw(&screw, [20., 0., 0.], u, v, 1e-7)
        .unwrap();
    let text = writer.finish().unwrap();
    assert_eq!(text.matches("MANIFOLD_SOLID_BREP").count(), 2);
    let solids = step_read::read_import_compound(&text, 1e-7).unwrap();
    assert_eq!(solids.len(), 2);
}

#[test]
fn malformed_finite_head_chart_refuses_before_thread_records() {
    for length in [0., 2.3] {
        let screw = source(false, false, length);
        let error = encode(
            &mut StepFile::new(),
            screw.thread(),
            Closure::screw(&screw).unwrap(),
            Placement {
                origin: [0.; 3],
                u: [0.; 3],
                v: [0.; 3],
                normal: [0., 0., 1.],
            },
            1e-7,
        )
        .err()
        .unwrap();
        assert_eq!(error, StepError::InvalidGeometry);
    }
}
