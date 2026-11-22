use super::*;

fn maps() -> (Affine, Affine) {
    let t = 2f64.powi(-24);
    let s = -2f64.powi(-25);
    (
        Affine::new([[0.25, -0.5, 0.75], [1., 0., 0.], [t, 1., 0.], [0., 0., 1.]]).unwrap(),
        Affine::new([[-0.5, 0.25, 0.125], [1., s, 0.], [0., 1., 0.], [0., 0., 1.]]).unwrap(),
    )
}

#[test]
fn actual_composed_shears_enclose_independent_exact_dyadic_cofactors() {
    let (construction, occurrence) = maps();
    let normals = PreparedNormals::new(construction, occurrence);
    let t = 2f64.powi(-24);
    let s = -2f64.powi(-25);
    // All operations in this independent closed-form result are exactly
    // representable dyadics: t*s=-2^-49 and 1+t*s fits in53 significand bits.
    let exact = [[1. + t * s, -t, 0.], [-s, 1., 0.], [0., 0., 1.]];
    for (column, expected) in normals.columns.into_iter().zip(exact) {
        for (bound, value) in column.into_iter().zip(expected) {
            assert!(bound.valid() && bound.lo <= value && value <= bound.hi);
        }
    }
}

#[test]
fn actual_emitted_facets_keep_orientation_and_refuse_flips_collapse_and_angular_error() {
    let (construction, occurrence) = maps();
    let normals = PreparedNormals::new(construction, occurrence);
    for local in [
        [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]],
        [[-0.5, 0.25, 0.75], [0.125, 0.5, 1.], [-0.25, 1., 0.5]],
    ] {
        let world = local.map(|p| occurrence.literal(construction.literal(p)));
        normals.require_normal(local, world, 1e-8).unwrap();
        let mut flipped = world;
        flipped.swap(1, 2);
        assert!(normals.require_normal(local, flipped, 1e-8).is_err());
        assert!(normals.require_normal(local, [world[0]; 3], 1e-8).is_err());
        let mut tilted = world;
        tilted[1][2] += 0.125;
        assert!(normals.require_normal(local, tilted, 1e-8).is_err());
    }
}
