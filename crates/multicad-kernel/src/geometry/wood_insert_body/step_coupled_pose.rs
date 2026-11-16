//! Transport a matched integer-grid file after its exact coefficient unit has
//! been earned by coupled source correspondence. The actual decoded operator
//! and source raw frame remain separate authorities.
use super::step::{fail, read_error};
use super::*;
use multicad_rust_kernel::wood_insert::Side;
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
// Transfer a nonnegative source-error bound through a measured operator
// difference. Overflow is a numerical refusal, never a cleared uncertainty.
fn propagated_error(previous: f64, operator: f64, bound: f64) -> Result<f64, GeometryError> {
    // Exact unit transport needs no multiplication reserve; this preserves
    // the original C0 addition order/bounds as well as its physical units.
    let transferred = if operator == 1. {
        B::point(bound)
    } else {
        B::point(operator) * B::point(bound)
    };
    let result = B::point(previous) + transferred;
    Ok(checked(result)?.hi)
}
// These numeric stages are also defined for finite transport inputs beyond
// admitted frame magnitudes: unrepresentable output must refuse, not saturate.
fn matrix_norm(delta: [[B; 3]; 3]) -> Result<f64, GeometryError> {
    root_sum(
        delta
            .into_iter()
            .flatten()
            .fold(B::ZERO, |s, x| s + x.square()),
    )
}
fn derivative_errors(
    previous: [f64; 2],
    matrix: f64,
    derivative: f64,
    ruling: f64,
) -> Result<[f64; 2], GeometryError> {
    let row = propagated_error(previous[0], matrix, derivative)?;
    propagated_error(previous[1], matrix, ruling).map(|ruling| [row, ruling])
}
fn scaled_volume(volume: B, scale: B) -> Result<B, GeometryError> {
    checked(volume * scale)
}
fn source_determinant(frame: PlaneFrame) -> Result<B, GeometryError> {
    let [u, v, n] = [frame.u, frame.v, frame.normal()].map(|a| a.map(B::point));
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
    Ok(det)
}
pub(super) fn transport(
    frame: PlaneFrame,
    local: super::step_coupled::QualifiedLocal<'_, '_, '_, '_>,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<WoodInsertStepVerification, GeometryError> {
    // The private local receipt earned SAME-mapping handedness==1 in
    // owned_insert/correspondence/coupled.rs before publication. It also
    // carries the source/boundary and raw, not already world, verification.
    let source = local.source();
    let boundary = local.boundary();
    let mapping = local.mapping();
    let mut proof = local.proof();
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
    let placement_error = norm(placement)?;
    let position = propagated_error(proof.maximum_position_error_mm, 1., placement_error)?;
    if position > proof.verified_tolerance_mm {
        return Err(fail(
            GeometryErrorKind::InvalidResult,
            "Actual STEP placement exceeds source position tolerance",
        ));
    }
    // Frobenius norm bounds the actual linear-map difference in every direction.
    let matrix = matrix_norm(delta)?;
    let mut derivative = 0_f64;
    // BoundaryOptions::validate (authored_internal_treatment/mesh/boundary/api.rs)
    // caps each sealed component at 8192 source curves. Even counting both
    // independently gives 2*8192*128=2,097,152 logical units, below 16M.
    // The insert constructor also subtracts outer curves from inner's budget.
    for side in [Side::Outer, Side::Inner] {
        for (_, curve) in boundary.source_boundary(side).source_curves() {
            if interrupt.is_some_and(|f| f()) {
                return Err(fail(
                    GeometryErrorKind::Cancelled,
                    "Coupled source derivative verification superseded",
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
    let ruling_bound = norm(bounds.map(|b| B::point(b.hi) - B::point(b.lo)))?;
    let [row_error, ruling_error] = derivative_errors(
        [
            proof.maximum_row_derivative_error_mm,
            proof.maximum_ruling_derivative_error_mm,
        ],
        matrix,
        derivative,
        ruling_bound,
    )?;
    // Source frame has literal u/v/cross columns. Its actual determinant is
    // independent of STEP's normalized orthogonal operator and uniform scale.
    let det = source_determinant(frame)?;
    let decoded_scale = interval(mapping.volume_scale_bounds().map_err(read_error)?);
    let decoded = scaled_volume(interval(proof.decoded_volume_bounds_mm3), decoded_scale)?;
    let source_volume = scaled_volume(interval(proof.source_volume_bounds_mm3), det)?;
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

#[cfg(test)]
mod controls {
    use super::*;
    #[test]
    fn actual_numeric_metric_and_error_transport_refuse_unrepresentable_bounds() {
        let rounded_zero = root_sum(B::new(-f64::from_bits(1), 0.)).unwrap();
        assert!((0. ..=f64::from_bits(1)).contains(&rounded_zero));
        assert!(norm([B::point(3.), B::point(4.), B::point(0.)]).unwrap() >= 5.);
        assert!(norm([B::point(f64::MAX); 3]).is_err());
        assert!(checked(B::new(f64::NAN, 1.)).is_err());
        assert!(propagated_error(f64::MAX, 2., f64::MAX).is_err());
        assert!(propagated_error(0.1, 1., 0.2).unwrap() >= 0.1 + 0.2);
        assert!(propagated_error(0.1, 0.5, 0.4).unwrap() >= 0.1 + 0.2);
        assert!(matrix_norm([[B::point(f64::MAX); 3]; 3]).is_err());
        assert!(matrix_norm([[B::ZERO; 3]; 3]).unwrap() < 1e-100);
        assert!(derivative_errors([0., 0.], 2., f64::MAX, 1.).is_err());
        assert!(derivative_errors([0., 0.], 2., 1., f64::MAX).is_err());
        let errors = derivative_errors([0.1, 0.2], 0.5, 0.4, 0.6).unwrap();
        assert!(errors[0] >= 0.3 && errors[1] >= 0.5);
        assert!(scaled_volume(B::point(f64::MAX), B::point(2.)).is_err());
        let scaled = scaled_volume(B::point(2.), B::point(3.)).unwrap();
        assert!(scaled.lo <= 6. && scaled.hi >= 6.);
        let det = source_determinant(PlaneFrame {
            u: [1. + 1e-10, 0., 0.],
            ..Default::default()
        })
        .unwrap();
        assert!(det.lo > 1.);
        assert!(
            source_determinant(PlaneFrame {
                u: [0.; 3],
                ..Default::default()
            })
            .is_err()
        );
    }
}
