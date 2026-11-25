use super::*;
use crate::axial_pipe_thread::authored::{AuthoredSection, AuthoredSpec, AuthoredThread};
fn clip(quarter: bool) -> StockClipPlan {
    let (tpi, e1, l1, l2, od, max) = if quarter {
        (18., 0.49163, 0.2278, 0.4018, 0.540, 0.088)
    } else {
        (27., 0.37360, 0.1615, 0.2639, 0.405, 0.096)
    };
    let pitch = 25.4 / tpi;
    let k = 1. / 32.;
    let m = 3f64.sqrt();
    let slopes = 1. / (m - k) + 1. / (m + k);
    let depth = pitch / slopes - (0.033 + max) * pitch;
    AuthoredThread::new(
        AuthoredSpec {
            major_radius_mm: (e1 * 25.4 - 2. * k * l1 * 25.4 + depth) / 2.,
            pitch_mm: pitch,
            length_mm: l2 * 25.4,
            radial_taper: k,
            section: AuthoredSection::NptFlat {
                depth_mm: depth,
                flat_width_mm: (pitch - depth * slopes) / 2.,
            },
            crest_radius_limit_mm: Some(od * 25.4 / 2.),
        },
        1e-7,
        1e-6,
    )
    .unwrap()
    .stock_clip_plan(Default::default(), None)
    .unwrap()
}
fn factors(length: f64, placed: bool) -> Vec<OccurrenceFactor> {
    // The native hex-head chain: neck to origin, thread along -z, then frame.
    let mut f = vec![(
        [[1., 0., 0.], [0., -1., 0.], [0., 0., -1.]],
        [0., 0., length],
    )];
    if placed {
        f.push((
            [[0.8, -0.36, 0.48], [0.6, 0.48, -0.64], [0., 0.8, 0.6]],
            [17., -8., 13.],
        ));
    }
    f
}
const HEAD: HexHeadIntent = HexHeadIntent {
    across_flats_mm: 18.,
    height_mm: 3.,
};
#[test]
fn authored_clipped_hex_head_export_qualifies_its_own_mapped_bytes() {
    for quarter in [false, true] {
        let clip = clip(quarter);
        let length = clip.source().spec().length_mm;
        for placed in [false, true] {
            let f = factors(length, placed);
            let export =
                export_authored_hex_head(&clip, HEAD, &f, Default::default(), None).unwrap();
            assert!(export.as_str().contains("MAPPED_ITEM"));
            let r = export.receipt();
            assert!(r.aggregate_work <= 272_000_000 && r.accounted_peak_bytes <= 64 * 1024 * 1024);
            assert!(r.source.world_sewn_reference_position_error_mm <= 1e-7);
            let v = export.world_sewn_volume_mm3();
            assert!(0. < v[0] && v[0] <= v[1]);
            let again = qualify_authored_hex_head_file(
                export.as_str(),
                &clip,
                HEAD,
                &f,
                Default::default(),
                None,
            )
            .unwrap();
            assert_eq!(again.world_sewn_volume_mm3, v);
            assert!(again.source.world_sewn_reference_position_error_mm <= 1e-7);
            // Wrong head intent and a 1 um moved chain are different sources.
            let wrong = HexHeadIntent {
                across_flats_mm: 17.,
                ..HEAD
            };
            assert!(
                qualify_authored_hex_head_file(
                    export.as_str(),
                    &clip,
                    wrong,
                    &f,
                    Default::default(),
                    None
                )
                .is_err()
            );
            let mut moved = f.clone();
            moved.last_mut().unwrap().1[0] += 1e-3;
            assert!(
                qualify_authored_hex_head_file(
                    export.as_str(),
                    &clip,
                    HEAD,
                    &moved,
                    Default::default(),
                    None
                )
                .is_err()
            );
        }
    }
}
#[test]
fn authored_clipped_hex_head_refuses_empty_chains_bad_options_and_cancellation() {
    let clip = clip(false);
    assert!(matches!(
        export_authored_hex_head(&clip, HEAD, &[], Default::default(), None),
        Err(Error::Invalid)
    ));
    let bad = HexHeadIntent {
        across_flats_mm: f64::NAN,
        ..HEAD
    };
    let f = factors(clip.source().spec().length_mm, false);
    assert!(export_authored_hex_head(&clip, bad, &f, Default::default(), None).is_err());
    assert!(matches!(
        export_authored_hex_head(&clip, HEAD, &f, Default::default(), Some(&|| true)),
        Err(Error::Cancelled)
    ));
}
