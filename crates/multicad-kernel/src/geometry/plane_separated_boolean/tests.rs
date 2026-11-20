use super::*;
use crate::{CurvedEdge2d, CurvedProfile};

fn policy() -> GeometryPolicy {
    GeometryPolicy::default()
}

fn wedge(corners: &[[f64; 2]], shift: f64) -> Solid {
    let n = corners.len();
    Solid::extrude_curved_profile_vector(
        &CurvedProfile {
            frame: PlaneFrame {
                origin: [shift, 0., 0.],
                ..Default::default()
            },
            outer: (0..n)
                .map(|i| CurvedEdge2d::Line {
                    a: corners[i],
                    b: corners[(i + 1) % n],
                })
                .collect(),
            holes: vec![],
        },
        [0., 0., 10.],
        policy(),
    )
    .unwrap()
}

const LOWER: [[f64; 2]; 3] = [[0., 0.], [10., 0.], [0., 10.]];
const UPPER: [[f64; 2]; 3] = [[10., 0.], [10., 10.], [0., 10.]];

#[test]
fn halves_meeting_on_a_diagonal_plane_share_no_interior() {
    let (a, b) = (wedge(&LOWER, 0.), wedge(&UPPER, 0.));
    assert!(
        a.boolean(&b, BooleanOperation::Intersect, policy())
            .unwrap()
            .is_none()
    );
    let kept = a
        .boolean(&b, BooleanOperation::Subtract, policy())
        .unwrap()
        .unwrap();
    assert!((kept.volume_mm3().unwrap() - 500.).abs() < 1e-9);
    // A separating plane of the second operand also certifies the pair.
    assert!(
        b.plane_separated_boolean(&a, BooleanOperation::Intersect, policy())
            .is_some()
    );
}

#[test]
fn overlap_unions_and_curved_operands_are_not_certified() {
    let a = wedge(&LOWER, 0.);
    // Shifted 1 mm into the lower half, the halves overlap.
    let overlapping = wedge(&UPPER, -1.);
    assert!(
        a.plane_separated_boolean(&overlapping, BooleanOperation::Intersect, policy())
            .is_none()
    );
    let b = wedge(&UPPER, 0.);
    assert!(
        a.plane_separated_boolean(&b, BooleanOperation::Union, policy())
            .is_none()
    );
    let rod = Solid::cylinder(
        PlaneFrame {
            origin: [20., 20., 0.],
            ..Default::default()
        },
        2.,
        5.,
        policy(),
    )
    .unwrap();
    assert!(
        a.plane_separated_boolean(&rod, BooleanOperation::Intersect, policy())
            .is_none()
    );
}
