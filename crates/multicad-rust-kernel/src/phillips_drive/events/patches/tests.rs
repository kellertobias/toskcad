use super::*;
use std::cell::Cell;
fn zero(x: I) {
    assert!(x.valid() && x.lo <= 0. && x.hi >= 0., "{x:?}")
}
fn support(surface: Surface, p: [I; 3]) {
    match surface {
        Surface::Plane {
            normal_bounds,
            offset_mm_bounds,
        } => zero(
            normal_bounds
                .map(b)
                .into_iter()
                .zip(p)
                .fold(I::point(0.), |s, (n, x)| s + n * x)
                - b(offset_mm_bounds),
        ),
        Surface::Cylinder {
            origin_mm_bounds,
            axis_bounds,
            radius_mm,
            ..
        } => {
            let o = origin_mm_bounds.map(b);
            let a = axis_bounds.map(b);
            let d = std::array::from_fn::<_, 3, _>(|i| p[i] - o[i]);
            let q = d
                .into_iter()
                .zip(a)
                .fold(I::point(0.), |s, (d, a)| s + d * a);
            zero(
                d.into_iter().fold(I::point(0.), |s, v| s + v.square())
                    - q.square()
                    - I::point(radius_mm).square(),
            )
        }
        Surface::Torus {
            major_radius_mm_bounds,
            minor_radius_mm,
            center_depth_mm,
            ..
        } => {
            let r = (p[0].square() + p[1].square()).sqrt();
            zero(
                (r - b(major_radius_mm_bounds)).square()
                    + (p[2] - I::point(center_depth_mm)).square()
                    - I::point(minor_radius_mm).square(),
            )
        }
        Surface::Cone {
            radius_at_zero_mm_bounds,
            radius_decrease_per_depth_bounds,
        } => zero(
            (p[0].square() + p[1].square()).sqrt() + p[2] * b(radius_decrease_per_depth_bounds)
                - b(radius_at_zero_mm_bounds),
        ),
    }
}
fn triangle_quad(p: [[f64; 3]; 4], u: f64, v: f64) -> [f64; 3] {
    if v <= u {
        std::array::from_fn(|i| (1. - u) * p[0][i] + (u - v) * p[1][i] + v * p[2][i])
    } else {
        std::array::from_fn(|i| (1. - v) * p[0][i] + u * p[2][i] + (v - u) * p[3][i])
    }
}
#[test]
fn every_metric_patch_retains_actual_support_normal_and_full_triangle_error() {
    for fixture in 0..=5 {
        let (s, d) = super::super::super::tests::fixture(fixture);
        let source = Cutter::new(s, d).unwrap();
        let atlas = source.surface_patches().unwrap();
        assert_eq!(atlas.patches().len(), 88);
        assert!(atlas.radial_margin_mm_bounds()[0] > 0.);
        for (index, patch) in atlas.patches().iter().enumerate() {
            let _key = patch.source_identity();
            let surface = atlas.topology().face_surface(patch.face_index()).unwrap();
            let center = atlas.patch_bounds(index, [[0.5; 2]; 2]).unwrap();
            let du = center.du_mm.map(b);
            let dv = center.dv_mm.map(b);
            let n = center.normal.map(b);
            let cross = std::array::from_fn::<_, 3, _>(|i| {
                du[(i + 1) % 3] * dv[(i + 2) % 3] - du[(i + 2) % 3] * dv[(i + 1) % 3]
            });
            let oriented = cross
                .into_iter()
                .zip(n)
                .fold(I::point(0.), |s, (v, n)| s + v * n)
                * I::point(atlas.patch_winding(index).unwrap() as f64);
            assert!(oriented.lo > 0.);
            for [u, v] in [[0., 0.], [0., 1.], [1., 0.], [1., 1.], [0.3, 0.7]] {
                let sample = atlas.patch_bounds(index, [[u, u], [v, v]]).unwrap();
                support(surface, sample.point_mm.map(b));
                zero(
                    sample
                        .normal
                        .map(b)
                        .into_iter()
                        .fold(I::point(0.), |s, n| s + n.square())
                        - I::point(1.),
                );
            }
            let uv = [[0.2, 0.3], [0.4, 0.5]];
            let coords = [[0.2, 0.4], [0.3, 0.4], [0.3, 0.5], [0.2, 0.5]];
            let corners = coords.map(|[u, v]| {
                atlas
                    .patch_bounds(index, [[u, u], [v, v]])
                    .unwrap()
                    .point_mm
                    .map(|i| b(i).midpoint())
            });
            let error = atlas.certify_quad(index, uv, corners).unwrap();
            assert!(error.positional_error_mm() > 0.);
            for [u, v] in [[0.2, 0.7], [0.5, 0.5], [0.8, 0.3]] {
                let p = atlas
                    .patch_bounds(index, [[0.2 + u * 0.1; 2], [0.4 + v * 0.1; 2]])
                    .unwrap()
                    .point_mm
                    .map(b);
                let q = triangle_quad(corners, u, v);
                let distance = norm(std::array::from_fn(|i| I::point(q[i]) - p[i]));
                assert!(
                    distance.hi <= error.positional_error_mm() + 1e-11,
                    "{fixture} {index} {distance:?} {error:?}"
                );
            }
        }
    }
}
#[test]
fn bad_patch_requests_and_cancel_fail_without_changing_authored_source() {
    let (s, d) = super::super::super::tests::fixture(1);
    let c = Cutter::new(s, d).unwrap();
    let before = c.clone();
    let n = Cell::new(0);
    let atlas = c
        .surface_patches_with_cancel(&|| {
            n.set(n.get() + 1);
            false
        })
        .unwrap();
    for stop in 1..=n.get() {
        let call = Cell::new(0);
        assert!(matches!(
            c.surface_patches_with_cancel(&|| {
                call.set(call.get() + 1);
                call.get() == stop
            }),
            Err(Error::Cancelled)
        ))
    }
    for (index, uv) in [
        (88, [[0., 1.]; 2]),
        (0, [[-1., 0.]; 2]),
        (0, [[0., 2.]; 2]),
        (0, [[1., 0.]; 2]),
        (0, [[f64::NAN, 1.]; 2]),
    ] {
        assert!(matches!(atlas.patch_bounds(index, uv), Err(Error::Invalid)))
    }
    let mut corners = [[0.; 3]; 4];
    corners[0][0] = f64::NAN;
    assert!(matches!(
        atlas.certify_quad(0, [[0., 1.]; 2], corners),
        Err(Error::Invalid)
    ));
    corners[0][0] = 1e10;
    assert!(matches!(
        atlas.certify_quad(0, [[0., 1.]; 2], corners),
        Err(Error::Limit)
    ));
    assert!(
        face_owner(
            atlas.topology().faces(),
            FaceKey::Envelope {
                role: EnvelopeRole::Floor,
                quarter: 4
            }
        )
        .is_err()
    );
    assert_eq!(c, before);
}

#[test]
fn saddle_requires_full_mixed_triangle_discrepancy_and_overflow_refuses() {
    // P(u,v)=[u,v,u*v] has zero pure second derivatives and an exact bilinear
    // patch. Replacing it by two literal triangles still errs by1/4 at center.
    let corners = [[0., 0., 0.], [1., 0., 0.], [1., 1., 1.], [0., 1., 0.]];
    let n = [I::point(0.), I::point(0.), I::point(1.)];
    let bound = quad_error(I::point(0.), 0., corners, n, n).unwrap();
    assert!(bound.bilinear_to_triangles_mm >= 0.25);
    assert!(bound.bilinear_to_triangles_mm < 0.25 + 1e-14);
    let p = triangle_quad(corners, 0.5, 0.5);
    assert_eq!(p[2] - 0.25, 0.25);
    assert_eq!(
        quad_error(I::new(0., f64::INFINITY), 0., corners, n, n),
        Err(Error::Unresolved)
    );
    let (s, d) = super::super::super::tests::fixture(0);
    let c = Cutter::new(s, d).unwrap();
    let atlas = c.surface_patches().unwrap();
    assert!(matches!(atlas.patch_winding(88), Err(Error::Invalid)));
}

#[test]
fn owner_schema_and_actual_cylinder_normal_chart_refuse_unsupported_data() {
    let (s, d) = super::super::super::tests::fixture(0);
    let c = Cutter::new(s, d).unwrap();
    let atlas = c.surface_patches().unwrap();
    assert!(FaceOwners::read(&atlas.topology().faces()[1..]).is_err());
    let mut duplicate = atlas.topology().faces().to_vec();
    duplicate.push(duplicate[0].clone());
    assert!(FaceOwners::read(&duplicate).is_err());
    let mut no_envelope = atlas.topology().faces().to_vec();
    no_envelope.retain(|f| {
        f.source_identity()
            != FaceKey::Envelope {
                role: EnvelopeRole::Floor,
                quarter: 1,
            }
    });
    assert!(FaceOwners::read(&no_envelope).is_err());
    let mut no_wall = atlas.topology().faces().to_vec();
    no_wall.retain(|f| {
        f.source_identity()
            != FaceKey::Wall {
                kind: WallKind::Inner,
                quarter: 2,
                lower: true,
            }
    });
    assert!(FaceOwners::read(&no_wall).is_err());
    let k = &c.corners[0];
    assert!(cylinder_normal(k, I::new(0., 0.25), 0.25).is_err());
    let normal = cylinder_normal(k, k.join_y, 0.25).unwrap();
    zero(normal.into_iter().fold(I::point(0.), |s, n| s + n.square()) - I::point(1.));
}
