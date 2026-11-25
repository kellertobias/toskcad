//! Placed Hanger readback: transport the local qualified proofs through the
//! decoded file operator (docs/architecture/hanger-placed-step-transport.md).
use super::*;
use crate::geometry::internal_thread_treatment_body::interval::Interval as B;
fn interval([lo, hi]: [f64; 2]) -> B {
    B::new(lo, hi)
}
fn upper(v: B) -> Result<f64, GeometryError> {
    v.checked()
        .map(|v| v.hi)
        .filter(|h| h.is_finite())
        .ok_or_else(|| {
            failure(
                GeometryErrorKind::Unsupported,
                "Placed hanger STEP bound overflow",
            )
        })
}
fn root(s: B) -> Result<f64, GeometryError> {
    let s = s.checked().ok_or_else(|| {
        failure(
            GeometryErrorKind::Unsupported,
            "Placed hanger STEP bound overflow",
        )
    })?;
    upper(B::new(s.lo.max(0.), s.hi).sqrt())
}
#[allow(clippy::too_many_arguments)]
pub(super) fn placed(
    qualified: &writer::QualifiedHangerStep,
    actual: &read::CompoundDocument<'_>,
    member: read::Member,
    frame: PlaneFrame,
    policy: GeometryPolicy,
    mut proof: HealedHangerStepVerification,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<HealedHangerStepVerification, GeometryError> {
    // The Hanger's own decoder reads this mapped item's operator by the same
    // rules as a placed singleton; scale is exactly one there. Decoded decimal
    // components are widened to their neighbouring floats.
    let (_, placement) = multicad_rust_kernel::step_read::owned_hanger::decode_mapped_item(
        actual.literal_input(),
        member.item_entity(),
        Default::default(),
        interrupt,
    )
    .map_err(|e| {
        failure(
            match e {
                multicad_rust_kernel::step_read::owned_hanger::Error::Limit => {
                    GeometryErrorKind::Limit
                }
                multicad_rust_kernel::step_read::owned_hanger::Error::Cancelled => {
                    GeometryErrorKind::Cancelled
                }
                multicad_rust_kernel::step_read::owned_hanger::Error::Unsupported => {
                    GeometryErrorKind::Unsupported
                }
                _ => GeometryErrorKind::InvalidResult,
            },
            format!("Actual hanger operator decoding failed: {e:?}"),
        )
    })?;
    let widen = |x: f64| B::new(x.next_down(), x.next_up());
    let a: [[B; 3]; 3] = placement.axes().map(|c| c.map(widen));
    let origin = placement.origin_mm().map(widen);
    let r = [frame.u, frame.v, frame.normal()].map(|c| c.map(B::point));
    // 1. Gershgorin bound on the Gram matrix: |A|_2^2 <= max_i sum_j |G_ij|.
    let mut gram = 0f64;
    for i in 0..3 {
        let mut row = 0f64;
        for j in 0..3 {
            let g = (0..3).fold(B::ZERO, |s, k| s + a[i][k] * a[j][k]);
            let g = g.checked().ok_or_else(|| {
                failure(
                    GeometryErrorKind::Unsupported,
                    "Placed hanger STEP bound overflow",
                )
            })?;
            row = upper(B::point(row) + B::point(g.lo.abs().max(g.hi.abs())))?;
        }
        gram = gram.max(row);
    }
    let norm_a = root(B::point(gram))?;
    // 2. Frobenius bound of D = A - R.
    let d: [[B; 3]; 3] = std::array::from_fn(|j| std::array::from_fn(|k| a[j][k] - r[j][k]));
    let norm_d = root(d.iter().flatten().fold(B::ZERO, |s, x| s + x.square()))?;
    // 3. Position: |A|e0 + |o - t + D c| + |D|_F r over the local source box.
    let local = qualified.source().bounds_mm();
    let centre = local.map(|v| (interval([v[0], v[0]]) + interval([v[1], v[1]])) / B::point(2.));
    let half = root(local.iter().fold(B::ZERO, |s, v| {
        s + ((B::point(v[1]) - B::point(v[0])) / B::point(2.)).square()
    }))?;
    let offset: [B; 3] = std::array::from_fn(|k| {
        (0..3).fold(origin[k] - B::point(frame.origin[k]), |s, j| {
            s + d[j][k] * centre[j]
        })
    });
    let offset = root(offset.into_iter().fold(B::ZERO, |s, x| s + x.square()))?;
    let position = upper(
        B::point(norm_a) * B::point(proof.source_to_healed_position_error_mm)
            + B::point(offset)
            + B::point(norm_d) * B::point(half),
    )?;
    if position > policy.linear_tolerance_mm {
        return Err(failure(
            GeometryErrorKind::InvalidResult,
            "Actual STEP placement exceeds the source position tolerance",
        ));
    }
    // 4. Derivatives: |A|e + |D|_F sup|dS| from the decoded support hull.
    let hull = qualified.derivative_hull(interrupt).map_err(writer_error)?;
    let transport =
        |e: f64, sup: f64| upper(B::point(norm_a) * B::point(e) + B::point(norm_d) * B::point(sup));
    let [eu, ev] = proof.source_support_derivative_error_mm;
    proof.source_support_derivative_error_mm = [
        transport(eu, hull.tensor_u_mm.max(hull.plane_mm))?,
        transport(ev, hull.tensor_v_mm.max(hull.plane_mm))?,
    ];
    proof.original_coedge_derivative_error_mm =
        transport(proof.original_coedge_derivative_error_mm, hull.curve_mm)?;
    let raw = upper(
        B::point(norm_a) * B::point(proof.source_to_raw_position_error_mm)
            + B::point(offset)
            + B::point(norm_d) * B::point(half),
    )?;
    proof.source_to_raw_position_error_mm = raw;
    proof.source_to_healed_position_error_mm = position;
    // 5. Volume scales by det A, which must be positive.
    let det = (0..3).fold(B::ZERO, |s, k| {
        let (i, j) = ((k + 1) % 3, (k + 2) % 3);
        s + a[0][k] * (a[1][i] * a[2][j] - a[1][j] * a[2][i])
    });
    let det = det.checked().filter(|v| v.lo > 0.).ok_or_else(|| {
        failure(
            GeometryErrorKind::InvalidResult,
            "Hanger mapping is singular or reflecting",
        )
    })?;
    let healed = interval(proof.healed_volume_bounds_mm3) * det;
    proof.healed_volume_bounds_mm3 = [healed.lo, healed.hi];
    // World enclosure: the authored frame image of the local box, widened.
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for corner in 0..8 {
        let p: [f64; 3] = std::array::from_fn(|j| local[j][(corner >> j) & 1]);
        for k in 0..3 {
            let w = (0..3).fold(B::point(frame.origin[k]), |s, j| {
                s + r[j][k] * B::point(p[j])
            });
            lo[k] = lo[k].min((w.lo - position).next_down());
            hi[k] = hi[k].max((w.hi + position).next_up());
        }
    }
    proof.correspondence_enclosure_mm = [lo, hi];
    Ok(proof)
}
