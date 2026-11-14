//! Transport a matched integer-grid file after its exact coefficient unit has
//! been earned by coupled source correspondence. The actual decoded operator
//! and source raw frame remain separate authorities.
use super::step::{fail, read_error};
use super::*;
use multicad_rust_kernel::{step_read::owned_insert as read, wood_insert::Side};
fn interval([lo, hi]: [f64; 2]) -> B {
    B::new(lo, hi)
}
fn checked(v: B) -> Result<B, GeometryError> {
    v.checked().ok_or_else(|| {
        fail(
            GeometryErrorKind::Unsupported,
            "Placed STEP enclosure overflow",
        )
    })
}
// Each term is a real square, so the exact sum is nonnegative even if
// outward addition rounds its lower endpoint one subnormal below zero.
fn root_sum(s: B) -> Result<f64, GeometryError> {
    let s = checked(s)?;
    Ok(checked(B::new(s.lo.max(0.), s.hi).sqrt())?.hi)
}
fn norm(v: [B; 3]) -> Result<f64, GeometryError> {
    root_sum(v.into_iter().fold(B::ZERO, |s, x| s + x.square()))
}
pub(super) fn transport(
    source: &source::Insert,
    frame: PlaneFrame,
    boundary: &source::Boundary,
    mapping: &read::PhysicalMapping<'_, '_>,
    mut proof: WoodInsertStepVerification,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<WoodInsertStepVerification, GeometryError> {
    if mapping.handedness().map_err(read_error)? != 1 {
        return Err(fail(
            GeometryErrorKind::InvalidResult,
            "Insert mapping reverses material orientation",
        ));
    }
    // The caller earned ACTUAL scale == owned coefficient unit exactly.
    // Work in physical-local mm: normalized EXPRESS columns are orthogonal,
    // so their exact Euclidean operator norm is one (not the tiny grid scale).
    let actual = [frame.u, frame.v, frame.normal()];
    let file = mapping.axes().map(|a| a.map(interval));
    let delta: [[B; 3]; 3] =
        std::array::from_fn(|j| std::array::from_fn(|k| file[j][k] - B::point(actual[j][k])));
    let bounds = source.bounds_mm().map(interval);
    let origin = mapping.origin_mm().map(interval);
    let placement: [B; 3] = std::array::from_fn(|k| {
        (0..3).fold(origin[k] - B::point(frame.origin[k]), |s, j| {
            s + delta[j][k] * bounds[j]
        })
    });
    let position =
        checked(B::point(norm(placement)?) + B::point(proof.maximum_position_error_mm))?.hi;
    if position > proof.verified_tolerance_mm {
        return Err(fail(
            GeometryErrorKind::InvalidResult,
            "Actual STEP placement exceeds source position tolerance",
        ));
    }
    // Frobenius norm bounds the actual linear-map difference in every direction.
    let matrix = root_sum(
        delta
            .into_iter()
            .flatten()
            .fold(B::ZERO, |s, x| s + x.square()),
    )?;
    let mut derivative = 0_f64;
    let mut work = 0usize;
    for side in [Side::Outer, Side::Inner] {
        for (_, curve) in boundary.source_boundary(side).source_curves() {
            if interrupt.is_some_and(|f| f()) {
                return Err(fail(
                    GeometryErrorKind::Cancelled,
                    "Coupled source derivative verification superseded",
                ));
            }
            work = work.saturating_add(128);
            if work > 16_000_000 {
                return Err(fail(
                    GeometryErrorKind::Limit,
                    "Coupled source derivative work limit",
                ));
            }
            // This is the ORIGINAL immutable curve's affine r/z/turn law,
            // including hidden cell seams. No rounded short trim is divided.
            // |P'|² = dr² + (tau*r*dv)² + dz², with true tau < 7.
            let [base, delta] = curve.affine_rzv_bounds().map(|row| row.map(B::from));
            let radius = base[0] + delta[0] * B::new(0., 1.);
            let angular = B::point(7.) * radius * delta[2];
            let squared = delta[0].square() + angular.square() + delta[1].square();
            derivative = derivative.max(root_sum(squared)?);
        }
    }
    let row_error = checked(
        B::point(proof.maximum_row_derivative_error_mm) + B::point(matrix) * B::point(derivative),
    )?
    .hi;
    let ruling_bound = norm(bounds.map(|b| B::point(b.hi) - B::point(b.lo)))?;
    let ruling_error = checked(
        B::point(proof.maximum_ruling_derivative_error_mm)
            + B::point(matrix) * B::point(ruling_bound),
    )?
    .hi;
    // Source frame has literal u/v/cross columns. Its actual determinant is
    // independent of STEP's normalized orthogonal operator and uniform scale.
    let [u, v, n] = actual.map(|a| a.map(B::point));
    let cross = [
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    ];
    let det = checked((0..3).fold(B::ZERO, |s, k| s + cross[k] * n[k]))?;
    if det.lo <= 0. {
        return Err(fail(
            GeometryErrorKind::Unsupported,
            "Source frame determinant is unresolved",
        ));
    }
    let decoded = checked(
        interval(proof.decoded_volume_bounds_mm3)
            * interval(mapping.volume_scale_bounds().map_err(read_error)?),
    )?;
    let source_volume = checked(interval(proof.source_volume_bounds_mm3) * det)?;
    let local = [proof.decoded_bounds_mm[0], proof.decoded_bounds_mm[1]];
    let mapped = mapping
        .map_box_mm(std::array::from_fn(|k| [local[0][k], local[1][k]]))
        .map_err(read_error)?;
    proof.decoded_bounds_mm = [mapped.map(|b| b[0]), mapped.map(|b| b[1])];
    proof.decoded_volume_bounds_mm3 = [decoded.lo, decoded.hi];
    proof.source_volume_bounds_mm3 = [source_volume.lo, source_volume.hi];
    proof.maximum_position_error_mm = position;
    proof.maximum_row_derivative_error_mm = row_error;
    proof.maximum_ruling_derivative_error_mm = ruling_error;
    if interrupt.is_some_and(|f| f()) {
        return Err(fail(
            GeometryErrorKind::Cancelled,
            "Placed STEP verification superseded",
        ));
    }
    Ok(proof)
}
