use super::*;
use crate::{
    display::{Tessellation, counterbored_tapped_orthogonal},
    helical_thread::Spec,
    manufacturing::certify,
    orthogonal::{Box3, Operation, boolean},
    orthogonal_boundary::build,
    orthogonal_pockets::identify_pocket,
};
fn thread(left_hand: bool, starts: u32) -> Thread {
    Thread::new(Spec {
        major_radius_mm: 1.5,
        pitch_mm: 0.5,
        length_mm: 3.5,
        depth_mm: 0.27,
        crest_width_mm: 0.08,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts,
        left_hand,
        root_radius_mm: None,
    })
    .unwrap()
}
#[test]
fn finite_counterbore_keeps_original_thread_and_actual_closed_stock_mesh() {
    for left in [false, true] {
        for starts in [1, 2] {
            let original = thread(left, starts);
            let seat = Counterbore::new(original, 3.5, 0.8).unwrap();
            assert_eq!(seat.thread(), &original);
            assert_eq!(seat.visible_thread_interval_mm(), [0.8, 3.5]);
            assert_eq!(seat.depth_mm(), 0.8);
            assert!(seat.contains_void([2., 0., 0.4], 0.));
            assert!(!seat.contains_void([2., 0., 1.], 0.));
            let stock = boolean(
                &[Box3 {
                    min: [0.; 3],
                    max: [20., 16., 10.],
                }],
                &[],
                Operation::Union,
                1e-7,
            )
            .unwrap()
            .unwrap();
            let boundary = build(stock.quads()).unwrap();
            let pocket = identify_pocket(
                &stock,
                &boundary,
                [6., 7., 11.],
                [6., 7., 6.5],
                3.5,
                1e-7,
                1e-5,
            )
            .unwrap();
            let options = Tessellation {
                deflection: 0.01,
                angle: 0.2,
                triangle_limit: 200_000,
            };
            let mesh = counterbored_tapped_orthogonal(
                &boundary,
                &seat,
                pocket,
                [[0.8, 0.6, 0.], [0.6, -0.8, 0.]],
                options,
            )
            .unwrap();
            let proof = certify(&mesh).unwrap();
            let expected = stock.volume_mm3() - seat.removed_volume_mm3();
            assert!(
                (proof.signed_volume_mm3 - expected).abs() < 0.15,
                "{} vs {expected}",
                proof.signed_volume_mm3
            );
            let bound = crate::helical_thread_mesh::counterbore_volume_error_bound(&seat, options);
            assert!(bound.is_finite() && bound > 0.);
            assert!(
                (proof.signed_volume_mm3 - expected).abs() <= bound,
                "actual volume error {} exceeds {bound}",
                (proof.signed_volume_mm3 - expected).abs()
            );
            assert!(mesh.face_indices.contains(&seat.cylinder_face()));
            assert!(mesh.face_indices.contains(&seat.shoulder_face()));
            assert!(mesh.face_indices.contains(&0));
            for triangle in mesh
                .indices
                .as_chunks::<3>()
                .0
                .iter()
                .zip(&mesh.face_indices)
                .filter(|(_, face)| **face == seat.shoulder_face())
            {
                for &index in triangle.0 {
                    assert!(
                        (mesh.vertices.as_chunks::<3>().0[index as usize][2] - 9.2).abs() < 1e-12
                    );
                }
            }
        }
    }
}
#[test]
fn invalid_or_unresolved_seats_cannot_replace_the_thread() {
    let thread = thread(false, 1);
    for (r, d) in [
        (1.5, 0.5),
        (1.4, 0.5),
        (3., 0.),
        (3., 3.5),
        (3., 4.),
        (f64::NAN, 0.5),
        (3., f64::INFINITY),
    ] {
        assert_eq!(
            Counterbore::new(thread, r, d).err(),
            Some(ThreadError::Invalid)
        );
    }
    assert_eq!(
        Counterbore::new(thread, 1e8, 0.5).err(),
        Some(ThreadError::Limit)
    );
    let seat = Counterbore::new(thread, 3.5, 0.8).unwrap();
    assert!(!seat.contains_void([0., 0., -0.1], 0.));
    assert!(!seat.contains_void([0., 0., 3.6], 0.));
    assert!(!seat.contains_void([f64::NAN, 0., 0.], 0.));
    assert!(!seat.contains_void([0.; 3], -1.));
    assert!(!seat.contains_void([0.; 3], f64::NAN));
    let options = Tessellation {
        deflection: 0.01,
        angle: 0.2,
        triangle_limit: 1,
    };
    assert!(crate::helical_thread_mesh::open_counterbored_cavity(&seat, options).is_err());
    assert_eq!(seat.thread(), &thread);
}
