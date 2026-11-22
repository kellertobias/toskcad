use super::*;
use crate::{
    Rational, manufacturing,
    rounded_thread::{Spec, Thread},
};
fn source(left_hand: bool) -> BlindBore {
    let pitch = 1.27;
    let angle = 27.5_f64.to_radians();
    let rounding = pitch * angle.cos() / (12. * (1. - angle.sin()));
    BlindBore::new(
        Thread::new(Spec {
            major_radius_mm: 3.175,
            pitch_mm: pitch,
            length_mm: 3. * pitch,
            depth_mm: pitch / (3. * angle.tan()),
            crest_radius_mm: rounding,
            root_radius_mm: rounding,
            flank_angle_deg: 27.5,
            starts: 1,
            left_hand,
        })
        .unwrap(),
        Rational::new(127, 20).unwrap(),
        Rational::new(127, 25).unwrap(),
    )
    .unwrap()
}
fn quality() -> Tessellation {
    Tessellation {
        deflection: 0.05,
        angle: 0.3,
        triangle_limit: 200_000,
    }
}
#[test]
fn actual_blind_floor_closes_the_original_rounded_rim_for_both_hands() {
    for hand in [false, true] {
        let bore = source(hand);
        let mesh = tessellate_blind_bore_with_interrupt(&bore, quality(), None).unwrap();
        let proof = manufacturing::certify(&mesh).unwrap();
        assert!(proof.signed_volume_mm3 > 0.);
        for face in 0..9 {
            assert!(
                mesh.face_indices.contains(&face),
                "missing actual face {face}"
            );
        }
        let height = bore.stock_height_mm().as_f64();
        let tip = bore.cutter().spec().length_mm;
        let mut floor_radii = Vec::new();
        for (tri, &face) in mesh
            .indices
            .as_chunks::<3>()
            .0
            .iter()
            .zip(&mesh.face_indices)
        {
            if face == 1 || face == 8 {
                let points = tri.map(|id| {
                    let k = 3 * id as usize;
                    [mesh.vertices[k], mesh.vertices[k + 1], mesh.vertices[k + 2]]
                });
                for p in points {
                    assert_eq!(p[2], if face == 1 { tip } else { height });
                }
                let area_z = (points[1][0] - points[0][0]) * (points[2][1] - points[0][1])
                    - (points[1][1] - points[0][1]) * (points[2][0] - points[0][0]);
                assert!(if face == 1 { area_z < 0. } else { area_z > 0. });
                if face == 1 {
                    for p in points {
                        if p[0] != 0. || p[1] != 0. {
                            floor_radii.push(p[0].hypot(p[1]));
                        }
                    }
                }
            }
        }
        let min = floor_radii.iter().copied().fold(f64::INFINITY, f64::min);
        let max = floor_radii.iter().copied().fold(0., f64::max);
        assert!(
            max - min > 0.1,
            "floor must use the noncircular actual tip rim"
        );
        let analytic = std::f64::consts::PI * bore.stock_radius_mm().as_f64().powi(2) * height
            - bore.cutter().volume_mm3();
        // Independent source-volume comparison supplements the actual whole
        // side/planar certificates; it cannot authorize a wrong boundary.
        assert!((proof.signed_volume_mm3 - analytic).abs() < analytic * 0.02);
    }
}

#[test]
fn blind_mesh_limits_and_late_cancel_leave_actual_source_retryable() {
    let bore = source(false);
    let q = quality();
    assert!(
        tessellate_blind_bore_with_interrupt(
            &bore,
            Tessellation {
                triangle_limit: 20,
                ..q
            },
            None
        )
        .is_err()
    );
    let calls = std::cell::Cell::new(0usize);
    let count = || {
        calls.set(calls.get() + 1);
        false
    };
    let original = tessellate_blind_bore_with_interrupt(&bore, q, Some(&count)).unwrap();
    let total = calls.get();
    assert!(total > 3);
    for stop in [0, total / 2, total - 1] {
        calls.set(0);
        let cancel = || {
            let n = calls.get();
            calls.set(n + 1);
            n >= stop
        };
        assert_eq!(
            tessellate_blind_bore_with_interrupt(&bore, q, Some(&cancel)),
            Err(Error::Cancelled)
        );
    }
    let retry = tessellate_blind_bore_with_interrupt(&bore, q, None).unwrap();
    assert_eq!(original, retry);
    assert_eq!(bore.stock_height_mm(), Rational::new(127, 25).unwrap());
}
