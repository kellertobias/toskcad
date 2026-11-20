//! Conservative whole-prism relations against a retained planar-and-boss union.
//!
//! Without a quarter-round band or axial cuts, the retained shell is exactly
//! the union of the chart's convex planar pieces and its attached
//! cylindrical bosses. A probe inside one of them, or outside all of them by
//! more than the construction tolerance, shares one unchanged operand. No
//! cut surface is computed.
use super::*;
use crate::geometry::helical_body::contained_boolean::{
    chart_bounds, probe_bounds, whole_operand_result,
};
use crate::geometry::internal_thread_treatment_body::interval::Interval as B;
use multicad_rust_kernel::helical_thread::BoxClassification;

/// `Some(true)` inside the finite boss cylinder, `Some(false)` outside it.
fn boss_relation(boss: &AuthoredBoss, local: [B; 3], margin: f64) -> Option<bool> {
    let [p0, p1] = boss.ends;
    let axis: [f64; 3] = std::array::from_fn(|i| p1[i] - p0[i]);
    let length = dot(axis, axis).sqrt();
    if !length.is_finite() || length <= 0. {
        return None;
    }
    // An orthonormal chart about the axis keeps each coordinate a linear
    // form of the box, avoiding the |v|^2 - t^2 dependency blow-up. Its
    // few-ulp non-orthogonality is covered by the slack below.
    let unit = axis.map(|a| a / length);
    let helper = if unit[0].abs() < 0.9 {
        [1., 0., 0.]
    } else {
        [0., 1., 0.]
    };
    let first = {
        let c = cross(unit, helper);
        let n = dot(c, c).sqrt();
        c.map(|x| x / n)
    };
    let second = cross(unit, first);
    let v: [B; 3] = std::array::from_fn(|i| local[i] - B::point(p0[i]));
    let along =
        |direction: [f64; 3]| (0..3).fold(B::ZERO, |sum, i| sum + v[i] * B::point(direction[i]));
    let t = along(unit);
    let radial = along(first).square() + along(second).square();
    if !t.valid() || !radial.valid() {
        return None;
    }
    let magnitude = p0
        .into_iter()
        .chain(p1)
        .chain(local.iter().flat_map(|b| [b.lo, b.hi]))
        .map(f64::abs)
        .fold(boss.radius.max(1.), f64::max);
    let slack = margin + 256. * f64::EPSILON * magnitude;
    let outer = (boss.radius + slack) * (boss.radius + slack);
    let inner = (boss.radius - slack).max(0.).powi(2);
    if t.hi < -slack || t.lo > length + slack || radial.lo > outer {
        Some(false)
    } else if t.lo > slack && t.hi < length - slack && radial.hi < inner {
        Some(true)
    } else {
        None
    }
}

impl Solid {
    pub(in crate::geometry) fn contained_warped_union_boolean(
        &self,
        other: &Self,
        operation: BooleanOperation,
        policy: GeometryPolicy,
    ) -> Option<Result<Option<Self>, GeometryError>> {
        if policy.boolean_fuzzy_mm != 0. || self.policy != policy || other.policy != policy {
            return None;
        }
        let (probe, source_is_left, (body, frame)) = match self.warped_union_body() {
            Some(found) => (other, true, found),
            None => (self, false, other.warped_union_body()?),
        };
        let chart = body.chart.as_ref()?;
        if chart.band.is_some() || !chart.axial.is_empty() {
            return None;
        }
        let local = chart_bounds(
            probe_bounds(probe)?,
            frame.origin,
            [frame.u, frame.v, frame.normal()],
        )
        .map(|[lo, hi]| B::new(lo, hi));
        let bounds = local.map(|b| [b.lo, b.hi]);
        let margin = policy.linear_tolerance_mm;
        // Inside any source settles it; otherwise every source must be
        // proved outside, and any unresolved source leaves the probe alone.
        for piece in &chart.pieces {
            match piece.classify_box(bounds, margin)? {
                true => {
                    return whole_operand_result(
                        self,
                        other,
                        source_is_left,
                        BoxClassification::Inside,
                        operation,
                    );
                }
                false => {}
            }
        }
        for boss in &chart.bosses {
            match boss_relation(boss, local, margin)? {
                true => {
                    return whole_operand_result(
                        self,
                        other,
                        source_is_left,
                        BoxClassification::Inside,
                        operation,
                    );
                }
                false => {}
            }
        }
        whole_operand_result(
            self,
            other,
            source_is_left,
            BoxClassification::Outside,
            operation,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn head_and_stock_probes_resolve_against_their_sources() {
        let policy = GeometryPolicy::default();
        let owner = EntityId::new();
        let head = Solid::rust_box(
            owner.uuid_bits(),
            [10, 10, 4].map(|n| Rational::new(n, 1).unwrap()),
            PlaneFrame::default(),
            policy,
        )
        .unwrap()
        .transformed(RigidTransform {
            translation_mm: [-5., -5., -4.],
            ..Default::default()
        })
        .unwrap();
        let shaft =
            crate::geometry::construction::cylinder(PlaneFrame::default(), 2., 10., owner, policy)
                .unwrap();
        let union =
            Solid::rust_join_planar_head_and_round_stock(&head, &shaft, owner, policy).unwrap();
        let probe = |c: [f64; 3]| {
            Solid::box_solid([2e-3; 3], policy)
                .unwrap()
                .transformed(RigidTransform {
                    translation_mm: c.map(|v| v - 1e-3),
                    ..Default::default()
                })
                .unwrap()
        };
        for (center, inside) in [
            ([4.99, 0., -2.], true),
            ([5.01, 0., -2.], false),
            ([1.99, 0., 5.], true),
            ([2.01, 0., 5.], false),
            ([0., 0., 9.99], true),
            ([0., 0., 10.01], false),
            ([3., 3., 1.], false),
        ] {
            let hit = union
                .boolean(&probe(center), BooleanOperation::Intersect, policy)
                .unwrap();
            assert_eq!(hit.is_some(), inside, "{center:?}");
            let kept = probe(center)
                .boolean(&union, BooleanOperation::Subtract, policy)
                .unwrap();
            assert_eq!(kept.is_none(), inside, "{center:?}");
        }
        // On the stock wall nothing is claimed.
        assert!(
            union
                .contained_warped_union_boolean(
                    &probe([2., 0., 5.]),
                    BooleanOperation::Intersect,
                    policy
                )
                .is_none()
        );
    }
}
