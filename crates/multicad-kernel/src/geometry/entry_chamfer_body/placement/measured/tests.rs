use super::*;
use multicad_rust_kernel::{Rational, display, prism::Prism};
fn options() -> Tessellation {
    Tessellation {
        deflection: 0.1,
        angle: 0.3,
        triangle_limit: 1000,
    }
}
fn frame() -> PlaneFrame {
    PlaneFrame {
        origin: [10., 20., 30.],
        u: [0.8, 0.6, 0.],
        v: [-0.36, 0.48, 0.8],
    }
}
fn linear_source(thin: bool) -> Prism {
    Prism::rectangular(
        7,
        [
            Rational::new(1, if thin { 1_i64 << 40 } else { 8 }).unwrap(),
            Rational::new(1, 4).unwrap(),
            Rational::new(1, 2).unwrap(),
        ],
    )
    .unwrap()
}
fn build(source: &Prism, quality: Tessellation) -> Result<SourceMesh, GeometryError> {
    Ok(SourceMesh {
        mesh: crate::rust_display::mesh(display::rust_prism(source, quality))?,
        chord_mm: 0.,
        normal_rad: 0.,
    })
}
#[test]
fn measured_actual_linear_source_keeps_arrays_and_qualifies_raw_placement() {
    let source = linear_source(false);
    let local = qualify(PlaneFrame::default(), options(), &mut |q| build(&source, q)).unwrap();
    let world = qualify(frame(), options(), &mut |q| build(&source, q)).unwrap();
    assert_eq!(local.indices, world.indices);
    assert_eq!(local.face_indices, world.face_indices);
    assert_ne!(local.vertices, world.vertices);
}
#[test]
fn earned_residuals_must_be_finite_positive_and_fit_the_actual_world_budget() {
    let source = linear_source(false);
    for (chord, normal) in [
        (f64::NAN, 0.),
        (-1., 0.),
        (0.11, 0.),
        (0., f64::INFINITY),
        (0., -1.),
        (0., 0.31),
    ] {
        assert!(
            qualify(PlaneFrame::default(), options(), &mut |q| {
                let mut c = build(&source, q)?;
                c.chord_mm = chord;
                c.normal_rad = normal;
                Ok(c)
            })
            .is_err()
        );
    }
    // Even valid conservative source bounds need reserve for the actual affine
    // distortion and emitted point/facet rounding; no bound is silently reduced.
    for (chord, normal) in [(0.1, 0.), (0., 0.3)] {
        assert!(
            qualify(frame(), options(), &mut |q| {
                let mut c = build(&source, q)?;
                c.chord_mm = chord;
                c.normal_rad = normal;
                Ok(c)
            })
            .is_err()
        );
    }
    assert!(
        qualify(frame(), options(), &mut |_| Err(error(
            "source refusal",
            GeometryErrorKind::Limit,
            "real source budget refusal"
        )))
        .is_err()
    );
}
#[test]
fn actual_frame_and_rounding_refusals_do_not_normalize_or_repair_source() {
    let source = linear_source(false);
    let called = std::cell::Cell::new(false);
    assert!(
        qualify(
            PlaneFrame {
                u: [1.01, 0., 0.],
                ..frame()
            },
            options(),
            &mut |q| {
                called.set(true);
                build(&source, q)
            }
        )
        .is_err()
    );
    assert!(!called.get());
    assert!(
        qualify(
            PlaneFrame {
                origin: [f64::INFINITY, 0., 0.],
                ..frame()
            },
            options(),
            &mut |q| build(&source, q)
        )
        .is_err()
    );
    assert!(
        qualify(
            frame(),
            Tessellation {
                deflection: 1e-16,
                ..options()
            },
            &mut |q| build(&source, q)
        )
        .is_err()
    );
    let thin = linear_source(true);
    assert!(
        qualify(
            PlaneFrame {
                origin: [1e6, -2e6, 3e6],
                ..frame()
            },
            options(),
            &mut |q| build(&thin, q)
        )
        .is_err()
    );
}
