//! Booleans of polygon prisms that share one exact cross-section.
//!
//! Two polygon prisms with identical rings and circular holes on the same
//! plane frame, offset only along that frame's normal, occupy the same
//! profile over two exact height intervals. Their Boolean is the profile
//! over the interval result: a union of overlapping or touching intervals, a
//! subtraction that removes one end, or the common interval. Results that
//! would split into two pieces, or operands with pockets or bore treatments,
//! are left to the other providers.
use super::*;

type Interval = (Rational, Rational);

fn exact(value: f64) -> Option<Rational> {
    let rational = rational_coordinate(value, 0.)?;
    (rational.as_f64() == value).then_some(rational)
}

fn plain(placed: &PlacedPolygon) -> bool {
    placed.pockets.is_empty() && placed.bore_treatments.is_none()
}

/// The occupied interval of `placed` along `base`'s normal, measured from
/// `base`'s origin, when the two prisms share one exact profile.
fn interval(base: &PlacedPolygon, placed: &PlacedPolygon) -> Option<Interval> {
    if base.frame.u != placed.frame.u
        || base.frame.v != placed.frame.v
        || base.prism.rings() != placed.prism.rings()
        || base.prism.circular_holes() != placed.prism.circular_holes()
    {
        return None;
    }
    let normal = base.frame.normal();
    let delta: [f64; 3] = std::array::from_fn(|i| placed.frame.origin[i] - base.frame.origin[i]);
    let offset = dot(delta, normal);
    if (0..3).any(|i| delta[i] != offset * normal[i]) {
        return None;
    }
    let offset = exact(offset)?;
    let height = placed.prism.height();
    Some(if placed.reverse_sweep {
        (offset.checked_sub(height).ok()?, offset)
    } else {
        (offset, offset.checked_add(height).ok()?)
    })
}

/// The exact interval Boolean, `Some(None)` for an empty result and `None`
/// when the result is not one interval.
fn combine(a: Interval, b: Interval, operation: BooleanOperation) -> Option<Option<Interval>> {
    let overlap = (a.0.max(b.0), a.1.min(b.1));
    match operation {
        BooleanOperation::Union => {
            (b.0 <= a.1 && a.0 <= b.1).then_some(Some((a.0.min(b.0), a.1.max(b.1))))
        }
        BooleanOperation::Intersect => (overlap.0 < overlap.1).then_some(Some(overlap)),
        BooleanOperation::Subtract => {
            if overlap.0 >= overlap.1 {
                return None;
            }
            match (b.0 <= a.0, a.1 <= b.1) {
                (true, true) => Some(None),
                (true, false) => Some(Some((b.1, a.1))),
                (false, true) => Some(Some((a.0, b.0))),
                (false, false) => None,
            }
        }
    }
}

impl Solid {
    pub(super) fn stacked_prism_boolean(
        &self,
        other: &Self,
        operation: BooleanOperation,
        policy: GeometryPolicy,
    ) -> Option<Result<Option<Self>, GeometryError>> {
        if policy.boolean_fuzzy_mm != 0. || policy != self.policy || policy != other.policy {
            return None;
        }
        let (a, b) = (self.polygon_prism()?, other.polygon_prism()?);
        if !plain(a) || !plain(b) {
            return None;
        }
        let zero = Rational::new(0, 1).expect("zero is exact");
        let own = if a.reverse_sweep {
            (zero.checked_sub(a.prism.height()).ok()?, zero)
        } else {
            (zero, a.prism.height())
        };
        let Some((low, high)) = combine(own, interval(a, b)?, operation)? else {
            return Some(Ok(None));
        };
        let height = high.checked_sub(low).ok()?;
        let prism = PolygonPrism::with_circular_holes(
            a.prism.rings().to_vec(),
            a.prism.circular_holes().to_vec(),
            height,
        )
        .ok()?;
        let base = if a.reverse_sweep { high } else { low }.as_f64();
        let normal = a.frame.normal();
        let frame = PlaneFrame {
            origin: std::array::from_fn(|i| a.frame.origin[i] + base * normal[i]),
            ..a.frame
        };
        let (origins, reversed): (Vec<_>, Vec<_>) = a.edges.iter().cloned().unzip();
        Some(
            Self::from_polygon_prism_sweep(
                prism,
                frame,
                a.reverse_sweep,
                a.feature,
                &origins,
                &reversed,
                policy,
            )
            .map(Some),
        )
    }
}

#[cfg(test)]
mod tests;
