//! Literal inverse-coordinate enclosure for retained authored material.
use super::*;
use crate::geometry::internal_thread_treatment_body::interval::Interval as B;
use multicad_rust_kernel::helical_thread::BoxClassification;
fn dot(a: [B; 3], b: [B; 3]) -> B {
    (0..3).fold(B::ZERO, |sum, i| sum + a[i] * b[i])
}
fn cross(a: [B; 3], b: [B; 3]) -> [B; 3] {
    std::array::from_fn(|i| {
        let j = (i + 1) % 3;
        let k = (i + 2) % 3;
        a[j] * b[k] - a[k] * b[j]
    })
}
fn unresolved() -> GeometryError {
    error(
        OPERATION,
        GeometryErrorKind::Unsupported,
        "Phillips point lies in an unresolved source or inverse-coordinate boundary enclosure",
    )
}
// This numerical primitive refuses singular literal frames independently of
// the caller's CAD frame policy. It never substitutes a transpose for A^-1.
fn local_box(frame: PlaneFrame, world: Point) -> Result<[B; 3], GeometryError> {
    let axes = [frame.u, frame.v, frame.normal()].map(|a| a.map(B::point));
    let cofactors = [
        cross(axes[1], axes[2]),
        cross(axes[2], axes[0]),
        cross(axes[0], axes[1]),
    ];
    let determinant = dot(axes[0], cofactors[0]);
    if !determinant.valid() || determinant.lo <= 0. {
        return Err(unresolved());
    }
    let offset = std::array::from_fn(|i| B::point(world[i]) - B::point(frame.origin[i]));
    let local = cofactors.map(|c| dot(c, offset) / determinant);
    Ok(local)
}
pub(in crate::geometry::helical_body) fn contains(
    source: &PhillipsScrew,
    frame: PlaneFrame,
    world: Point,
    policy: GeometryPolicy,
) -> Result<bool, GeometryError> {
    policy.validate()?;
    policy.point(world)?;
    frame.checked(policy)?;
    let local = local_box(frame, world)?;
    if local
        .into_iter()
        .any(|v| !v.valid() || v.width() > policy.linear_tolerance_mm / 8.)
    {
        return Err(unresolved());
    }
    match source
        .classify_box(local.map(|v| [v.lo, v.hi]))
        .map_err(failure)?
    {
        BoxClassification::Inside => Ok(true),
        BoxClassification::Outside => Ok(false),
        BoxClassification::UnresolvedBoundary => Err(unresolved()),
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn actual_raw_inverse_retains_material_void_and_uncertain_world_boundary() {
        let policy = GeometryPolicy::default();
        for index in 0..6 {
            for left in [false, true] {
                let source = super::super::tests::fixture(index, 0.2, left);
                for frame in [
                    PlaneFrame::default(),
                    PlaneFrame {
                        origin: [13., -9., 5.],
                        u: [0.8, 0.6, 0.],
                        v: [-0.6, 0.8, 0.],
                    },
                    PlaneFrame {
                        origin: [-7., 4., 2.],
                        u: [1. + 1e-12, 0., 0.],
                        v: [0., 1., 0.],
                    },
                ] {
                    for (p, expected) in [
                        ([0., 0., 4.], true),
                        ([0., 0., source.head_top_z() + 0.1], false),
                        ([20., 0., 0.], false),
                    ] {
                        assert_eq!(
                            contains(&source, frame, position(frame, p), policy).unwrap(),
                            expected
                        );
                    }
                    let world = position(frame, [0., 0., source.head_top_z()]);
                    assert_eq!(
                        contains(&source, frame, world, policy).unwrap_err().kind,
                        GeometryErrorKind::Unsupported
                    );
                }
                assert!(
                    contains(&source, PlaneFrame::default(), [f64::NAN, 0., 0.], policy).is_err()
                );
                assert!(
                    contains(
                        &source,
                        PlaneFrame {
                            u: [0.; 3],
                            ..PlaneFrame::default()
                        },
                        [0.; 3],
                        policy
                    )
                    .is_err()
                );
                let far = PlaneFrame {
                    origin: [1e6; 3],
                    ..PlaneFrame::default()
                };
                assert!(contains(&source, far, position(far, [0., 0., 4.]), policy).unwrap());
                assert_eq!(
                    contains(
                        &source,
                        far,
                        position(far, [0., 0., 4.]),
                        GeometryPolicy {
                            linear_tolerance_mm: 1e-16,
                            ..policy
                        }
                    )
                    .unwrap_err()
                    .kind,
                    GeometryErrorKind::Input
                );
                // Valid world coordinates and frame, but the directed inverse
                // cannot meet the policy's coordinate-width allowance.
                let opposite = PlaneFrame {
                    origin: [-1e7; 3],
                    ..PlaneFrame::default()
                };
                assert_eq!(
                    contains(&source, opposite, [1e7; 3], policy)
                        .unwrap_err()
                        .kind,
                    GeometryErrorKind::Unsupported
                );
            }
        }
    }
}

#[cfg(test)]
mod inverse_tests {
    use super::*;
    #[test]
    fn singular_literal_coordinate_maps_have_no_certified_inverse() {
        for frame in [
            PlaneFrame {
                u: [0.; 3],
                ..PlaneFrame::default()
            },
            PlaneFrame {
                v: [f64::NAN, 0., 0.],
                ..PlaneFrame::default()
            },
        ] {
            assert_eq!(
                local_box(frame, [0.; 3]).unwrap_err().kind,
                GeometryErrorKind::Unsupported
            );
        }
        let bounds = local_box(PlaneFrame::default(), [1., 2., 3.]).unwrap();
        for (bound, expected) in bounds.into_iter().zip([1., 2., 3.]) {
            assert!(bound.lo <= expected && expected <= bound.hi);
        }
    }
}
