//! Transfer a certified local tessellation without spending its world error budget twice.
use super::*;
pub(in crate::geometry) mod measured;
use crate::geometry::arithmetic::{difference_with_error, dot_error, sum_error};
use multicad_rust_kernel::display::Tessellation;

pub(super) fn up(a: f64, b: f64) -> f64 {
    (a + b).next_up()
}
pub(super) fn product(a: f64, b: f64) -> f64 {
    (a * b).next_up()
}
fn refused() -> GeometryError {
    error(
        "Entry chamfer mesh placement",
        GeometryErrorKind::Unsupported,
        "Requested world mesh accuracy cannot be certified",
    )
}
pub(super) fn cross_error(a: Point, ae: Point, b: Point, be: Point) -> (Point, Point) {
    let pairs = [(1, 2), (2, 0), (0, 1)];
    let mut value = [0.; 3];
    let mut error = [0.; 3];
    for (i, (j, k)) in pairs.into_iter().enumerate() {
        let (v, mut e) = dot_error([a[j], a[k], 0.], [b[k], -b[j], 0.]);
        for (j, k) in [(j, k), (k, j)] {
            e = up(e, product(a[j].abs(), be[k]));
            e = up(e, product(b[k].abs(), ae[j]));
            e = up(e, product(ae[j], be[k]));
        }
        value[i] = v;
        error[i] = e;
    }
    (value, error)
}
fn chart(axes: [Point; 3]) -> Result<(), GeometryError> {
    // ||AᵀA-I||₂ <= its entrywise L1 norm. This proves ||A||<2,
    // sigma_min²>1/2, ||A||<1.00001 and condition(A)<1.00001, without normalizing A.
    let mut gram = 0.;
    for (i, a) in axes.iter().enumerate() {
        for (j, b) in axes.iter().enumerate() {
            let (v, e) = dot_error(*a, *b);
            let (v, d) = sum_error(v, if i == j { -1. } else { 0. });
            gram = up(gram, up(v.abs(), up(e, d)));
        }
    }
    let (normal, ne) = cross_error(axes[0], [0.; 3], axes[1], [0.; 3]);
    let (det, mut de) = dot_error(normal, axes[2]);
    for i in 0..3 {
        de = up(de, product(ne[i], axes[2][i].abs()));
    }
    if !gram.is_finite() || gram > 1e-6 || (det - de).next_down() <= 0. {
        return Err(refused());
    }
    Ok(())
}
fn placed(frame: PlaneFrame, axes: [Point; 3], p: Point) -> (Point, f64) {
    let point = position(frame, p);
    let mut error = 0.;
    for i in 0..3 {
        let (d, de) = dot_error(axes.map(|a| a[i]), p);
        let (v, ve) = sum_error(frame.origin[i], d);
        let (delta, residual) = sum_error(point[i], -v);
        error = up(error, up(delta.abs(), up(residual, up(de, ve))));
    }
    (point, error)
}
// Compare exact affine and actual emitted facet normals directly. This keeps
// the shared translation out of the error bound for small, distant facets.
fn normal_transfer(
    axes: [Point; 3],
    a: Point,
    ae: Point,
    b: Point,
    be: Point,
    world: [Point; 3],
    angle: f64,
) -> bool {
    let linear = |v: Point, e: Point| {
        let mut value = [0.; 3];
        let mut error = [0.; 3];
        for i in 0..3 {
            let (d, mut de) = dot_error(axes.map(|a| a[i]), v);
            for j in 0..3 {
                de = up(de, product(axes[j][i].abs(), e[j]));
            }
            value[i] = d;
            error[i] = de;
        }
        (value, error)
    };
    let (a, ae) = linear(a, ae);
    let (b, be) = linear(b, be);
    let (exact, ee) = cross_error(a, ae, b, be);
    let (a, ae) = difference_with_error(world[1], world[0]);
    let (b, be) = difference_with_error(world[2], world[0]);
    let (actual, me) = cross_error(a, ae, b, be);
    let (d, mut de) = dot_error(exact, actual);
    for i in 0..3 {
        de = up(de, product(exact[i].abs(), me[i]));
        de = up(de, product(actual[i].abs(), ee[i]));
        de = up(de, product(ee[i], me[i]));
    }
    let lower = (d - de).next_down();
    let (c, ce) = cross_error(exact, ee, actual, me);
    let upper = (0..3).fold(0., |sum, i| up(sum, up(c[i].abs(), ce[i])));
    lower > 0. && (upper / lower).next_up() <= angle
}
pub(super) fn mesh(
    source: &Chamfer,
    frame: PlaneFrame,
    options: Tessellation,
) -> Result<Mesh, GeometryError> {
    qualified_source_mesh(frame, options, |quality| {
        crate::rust_display::mesh(source::tessellate(source, quality))
    })
}
/// Transfer a source-certified local candidate, retaining all original actual
/// point, normal and literal-intersection gates for either owned source.
pub(in crate::geometry) fn qualified_source_mesh(
    frame: PlaneFrame,
    options: Tessellation,
    mut build: impl FnMut(Tessellation) -> Result<Mesh, GeometryError>,
) -> Result<Mesh, GeometryError> {
    source_mesh(frame, options, &mut build, AngleBudget::Standard)
}
/// Same complete transfer proof with a smaller accepted literal normal error.
/// This does not change caller quality or any source/chord/embedding condition.
pub(in crate::geometry) fn qualified_source_mesh_with_small_normal_reserve(
    frame: PlaneFrame,
    options: Tessellation,
    mut build: impl FnMut(Tessellation) -> Result<Mesh, GeometryError>,
) -> Result<Mesh, GeometryError> {
    source_mesh(frame, options, &mut build, AngleBudget::SmallNormalReserve)
}
#[derive(Clone, Copy)]
enum AngleBudget {
    Standard,
    SmallNormalReserve,
}
impl AngleBudget {
    fn source(self, angle: f64) -> f64 {
        match self {
            Self::Standard => angle * 0.75,
            Self::SmallNormalReserve => angle * 0.99,
        }
    }
    fn transfer(self, angle: f64) -> f64 {
        match self {
            Self::Standard => angle / 5.,
            Self::SmallNormalReserve => angle / 200.,
        }
    }
}
fn source_mesh(
    frame: PlaneFrame,
    options: Tessellation,
    build: &mut dyn FnMut(Tessellation) -> Result<Mesh, GeometryError>,
    angle_budget: AngleBudget,
) -> Result<Mesh, GeometryError> {
    if frame == PlaneFrame::default() {
        let mesh = build(options)?;
        super::embedding::certify(&mesh)?;
        return Ok(mesh);
    }
    let axes = [frame.u, frame.v, frame.normal()];
    chart(axes)?;
    let local = Tessellation {
        deflection: options.deflection / 2.,
        angle: angle_budget.source(options.angle),
        ..options
    };
    if local.angle < 0.01 {
        return Err(refused());
    }
    let mut mesh = build(local)?;
    let points = mesh.vertices.as_chunks::<3>().0;
    let world = points
        .iter()
        .map(|p| placed(frame, axes, *p))
        .collect::<Vec<_>>();
    if world.iter().any(|(p, e)| {
        !e.is_finite() || p.iter().any(|v| !v.is_finite()) || *e > options.deflection / 4.
    }) {
        return Err(refused());
    }
    // The exact affine image has normal A^-T n. Its source/facet angle is
    // <= condition(A)*local.angle. Half the chord budget plus point rounding
    // uses less than 0.751 of the request. The standard angular split uses
    // <0.951 of the request. The small-normal-error split uses
    // 0.99*1.00001+1/200 <0.996; it requires the correspondingly stricter
    // actual literal-facet transfer below. Neither policy loosens source quality.
    // Certify the additional emitted-point error
    // triangle by triangle, relative to a lower bound on the exact local area.
    for ids in mesh.indices.as_chunks::<3>().0 {
        let [i, j, k] = ids.map(|i| i as usize);
        let (a, ae) = difference_with_error(points[j], points[i]);
        let (b, be) = difference_with_error(points[k], points[i]);
        let (n, ne) = cross_error(a, ae, b, be);
        let lower = (0..3)
            .map(|i| (n[i].abs() - ne[i]).next_down().max(0.))
            .fold(0., f64::max);
        let length = |v: Point, e: Point| (0..3).fold(0., |sum, i| up(sum, up(v[i].abs(), e[i])));
        let e = world[i].1.max(world[j].1).max(world[k].1);
        let perturbation = up(
            product(product(4., e), up(length(a, ae), length(b, be))),
            product(4., product(e, e)),
        );
        // asin(ratio) <= 2 ratio for ratio <= 1/4. Each named policy's
        // actual transfer reserve plus its conditioned source angle is <request.
        if !rounding_angle_fits(perturbation, lower, angle_budget.transfer(options.angle))
            && !normal_transfer(
                axes,
                a,
                ae,
                b,
                be,
                [world[i].0, world[j].0, world[k].0],
                angle_budget.transfer(options.angle),
            )
        {
            return Err(refused());
        }
    }
    mesh.vertices = world.into_iter().flat_map(|(p, _)| p).collect();
    super::embedding::certify(&mesh)?;
    Ok(mesh)
}

// A zero/underflowed certified area must use the direct normal proof. A
// negative denominator can otherwise make a negative ratio look harmless.
fn rounding_angle_fits(perturbation: f64, normal_lower: f64, angle: f64) -> bool {
    let denominator = (0.5 * normal_lower).next_down();
    if !normal_lower.is_finite()
        || normal_lower <= 0.
        || !denominator.is_finite()
        || denominator <= 0.
        || !perturbation.is_finite()
        || perturbation < 0.
    {
        return false;
    }
    let ratio = (perturbation / denominator).next_up();
    ratio.is_finite() && (0. ..=0.25).contains(&ratio) && product(2., ratio) <= angle
}

#[cfg(test)]
mod tests {
    #[test]
    fn underflowed_area_cannot_pass_a_negative_rounding_ratio() {
        let tiny = f64::from_bits(1);
        let old_ratio = (tiny / (0.5 * 0_f64).next_down()).next_up();
        assert!(old_ratio.is_finite() && old_ratio < 0.);
        for lower in [0., tiny, 2. * tiny, f64::INFINITY] {
            assert!(!rounding_angle_fits(tiny, lower, 0.3), "lower {lower:e}");
        }
        assert!(!rounding_angle_fits(f64::INFINITY, 1., 0.3));
        assert!(!rounding_angle_fits(-tiny, 1., 0.3));
        assert!(!rounding_angle_fits(1., 1., 0.3));
        assert!(rounding_angle_fits(1e-8, 1e-3, 0.3));
    }
    #[test]
    fn original_entry_source_uses_the_unchanged_standard_transfer() {
        let root = 1. / (4. * 3_f64.sqrt());
        let body = Solid::rust_helical_thread_entry_chamfer(
            HelicalThread {
                major_radius_mm: 3.,
                pitch_mm: 1.,
                length_mm: 6.6,
                depth_mm: 0.875 * 3_f64.sqrt() / 2. - root,
                crest_width_mm: 0.125,
                lower_flank_deg: 30.,
                upper_flank_deg: 30.,
                radial_taper: 0.,
                starts: 1,
                left_hand: false,
            },
            Some(root),
            0.3,
            None,
            EntityId::new(),
            PlaneFrame::default(),
            GeometryPolicy::default(),
        )
        .unwrap();
        let mesh = body.manufacturing_mesh(MeshOptions::default()).unwrap();
        assert!(!mesh.indices.is_empty());
        assert!(
            mesh.indices
                .iter()
                .all(|&i| (i as usize) < mesh.vertices.len())
        );
    }
    use super::*;
    #[test]
    fn proper_chart_and_actual_normals_are_required_without_normalizing_away_errors() {
        let axes = [[0., 1., 0.], [-1., 0., 0.], [0., 0., 1.]];
        chart(axes).unwrap();
        assert!(chart([[1., 0., 0.], [0., 1., 0.], [0., 0., -1.]]).is_err());
        assert!(chart([[1.01, 0., 0.], [0., 1., 0.], [0., 0., 1.]]).is_err());
        let a = [1., 0., 0.];
        let b = [0., 1., 0.];
        assert!(normal_transfer(
            axes,
            a,
            [0.; 3],
            b,
            [0.; 3],
            [[1e6, 1e6, 1e6], [1e6, 1e6 + 1., 1e6], [1e6 - 1., 1e6, 1e6]],
            0.01
        ));
        assert!(!normal_transfer(
            axes,
            a,
            [0.; 3],
            b,
            [0.; 3],
            [[0.; 3], [0., 1., 0.], [1., 0., 0.]],
            0.01
        ));
        assert!(!normal_transfer(
            axes,
            a,
            [0.; 3],
            b,
            [0.; 3],
            [[0.; 3]; 3],
            0.01
        ));
    }
    #[test]
    fn small_normal_reserve_preserves_actual_source_quality_and_refuses_literal_rounding() {
        use multicad_rust_kernel::{Rational, display, prism::Prism};
        let source = Prism::rectangular(
            7,
            [
                Rational::new(1, 8).unwrap(),
                Rational::new(1, 4).unwrap(),
                Rational::new(1, 2).unwrap(),
            ],
        )
        .unwrap();
        let frame = PlaneFrame {
            origin: [10., 20., 30.],
            u: [0.8, 0.6, 0.],
            v: [-0.36, 0.48, 0.8],
        };
        let options = Tessellation {
            deflection: 0.1,
            angle: 0.3,
            triangle_limit: 1000,
        };
        let build = |quality| crate::rust_display::mesh(display::rust_prism(&source, quality));
        let old = qualified_source_mesh(frame, options, build).unwrap();
        let new = qualified_source_mesh_with_small_normal_reserve(frame, options, build).unwrap();
        assert_eq!(old.vertices, new.vertices);
        assert_eq!(old.indices, new.indices);
        assert_eq!(old.face_indices, new.face_indices);
        assert!(
            qualified_source_mesh_with_small_normal_reserve(
                frame,
                Tessellation {
                    angle: 0.01,
                    ..options
                },
                build
            )
            .is_err()
        );
        let thin = Prism::rectangular(
            7,
            [
                Rational::new(1, 1_000_000_000_000).unwrap(),
                Rational::new(1, 8).unwrap(),
                Rational::new(1, 4).unwrap(),
            ],
        )
        .unwrap();
        assert!(
            qualified_source_mesh_with_small_normal_reserve(
                PlaneFrame {
                    origin: [1e6, -2e6, 3e6],
                    ..frame
                },
                options,
                |quality| crate::rust_display::mesh(display::rust_prism(&thin, quality))
            )
            .is_err()
        );
    }
    #[test]
    fn actual_dyadic_source_transfer_refuses_tight_angles_point_rounding_and_collapsed_thin_facets()
    {
        use multicad_rust_kernel::{Rational, display, prism::Prism};
        let frame = PlaneFrame {
            origin: [1e6, -2e6, 3e6],
            u: [0.8, 0.6, 0.],
            v: [-0.36, 0.48, 0.8],
        };
        let options = Tessellation {
            deflection: 0.1,
            angle: 0.3,
            triangle_limit: 1000,
        };
        let source = Prism::rectangular(
            7,
            [
                Rational::new(1, 8).unwrap(),
                Rational::new(1, 4).unwrap(),
                Rational::new(1, 2).unwrap(),
            ],
        )
        .unwrap();
        let build = |quality| crate::rust_display::mesh(display::rust_prism(&source, quality));
        assert!(qualified_source_mesh(frame, options, build).is_ok());
        assert!(
            qualified_source_mesh(
                frame,
                Tessellation {
                    angle: 0.01,
                    ..options
                },
                build
            )
            .is_err()
        );
        assert!(
            qualified_source_mesh(
                frame,
                Tessellation {
                    deflection: 1e-12,
                    ..options
                },
                build
            )
            .is_err()
        );
        let tiny = Prism::rectangular(
            7,
            [
                Rational::new(1, 1_i64 << 40).unwrap(),
                Rational::new(1, 1).unwrap(),
                Rational::new(1, 1).unwrap(),
            ],
        )
        .unwrap();
        // This is real exact source material. At this distance its tiny side
        // collapses in emitted binary64 coordinates; no normal/mesh is repaired.
        let refused = qualified_source_mesh(frame, options, |quality| {
            crate::rust_display::mesh(display::rust_prism(&tiny, quality))
        })
        .unwrap_err();
        assert_eq!(refused.kind, GeometryErrorKind::Unsupported);
    }
}
