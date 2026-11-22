use super::*;
use crate::brep::{Frame, NurbsSurface};

fn actual_box() -> Boundary {
    let text = include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../multicad-app/tests/fixtures/exchange/corpus/occt-placed-box.step"
    ));
    crate::step_read::brep::read_solids(text, 1e-7, None)
        .unwrap()
        .remove(0)
        .boundary()
        .clone()
}
fn nurbs_plane(frame: Frame, offset: f64) -> Surface {
    let points = [-128., 128.]
        .into_iter()
        .flat_map(|u| {
            [-128., 128.].into_iter().map(move |v| {
                std::array::from_fn(|i| {
                    frame.origin[i] + u * frame.x[i] + v * frame.y[i] + offset * frame.z[i]
                })
            })
        })
        .collect();
    let surface = NurbsSurface::new(
        [1, 1],
        [2, 2],
        points,
        None,
        [vec![0., 1.], vec![0., 1.]],
        [vec![2, 2], vec![2, 2]],
    )
    .unwrap();
    Surface::Nurbs {
        surface,
        closed: [false, false],
    }
}
fn replace_plane(boundary: &mut Boundary, offset: f64) {
    let Surface::Plane(frame) = boundary.faces[0].surface else {
        panic!()
    };
    boundary.faces[0].surface = nurbs_plane(frame, offset);
}

#[test]
fn unresolved_nurbs_orientation_never_claims_normalized_topology() {
    let mut boundary = actual_box();
    replace_plane(&mut boundary, 0.);
    for trim in &mut boundary.faces[0].loops {
        trim.uses.reverse();
        for used in &mut trim.uses {
            used.1 = !used.1;
        }
    }
    // These stored uses need the original full orientation normalization.
    // Checking their parity prematurely would falsely reject this geometry.
    assert_eq!(boundary.validate(), Err(Error::Invalid));
    assert!(validate_incidence(boundary.clone(), 1e-7, None).is_ok());
    // Only the Full path normalizes and then earns an actual solid.
    let full = Solid::new(boundary, 1e-7, None).unwrap();
    assert!(full.volume_mm3() > 0.);
}
#[test]
fn actual_parallel_nurbs_support_disagreement_still_refuses() {
    let mut boundary = actual_box();
    replace_plane(&mut boundary, 0.125);
    // The actual positive bilinear net lies in a parallel plane displaced
    // by 0.125mm, so its original box edges cannot lie on this support.
    assert_eq!(
        validate_incidence(boundary, 1e-7, None),
        Err(Error::Invalid)
    );
}
#[test]
fn refusal_only_actual_boundary_keeps_cancellation() {
    let mut boundary = actual_box();
    replace_plane(&mut boundary, 0.);
    let stop = || true;
    assert_eq!(
        validate_incidence(boundary, 1e-7, Some(&stop)),
        Err(Error::Cancelled)
    );
}
