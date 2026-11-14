//! Whole-source containment/disjointness only. No cut surface is approximated
//! or manufactured by this path: a proved result shares an unchanged operand.
use super::*;
use crate::geometry::internal_thread_treatment_body::interval::Interval as B;
use multicad_rust_kernel::helical_thread::BoxClassification;

fn exact(value: Rational) -> B {
    // Rational::as_f64 rounds once. Both neighbours enclose the original
    // rational, including decimals whose binary representation is inexact.
    let value = value.as_f64();
    B::new(value.next_down(), value.next_up())
}
fn world_bounds(local: [B; 3], frame: PlaneFrame) -> [B; 3] {
    let axes = [frame.u, frame.v, frame.normal()];
    std::array::from_fn(|i| {
        (0..3).fold(B::point(frame.origin[i]), |sum, j| {
            sum + B::point(axes[j][i]) * local[j]
        })
    })
}
fn probe_bounds(probe: &Solid) -> Option<[B; 3]> {
    if let Some((source, frame)) = probe.rust_prism_frame() {
        if source.rounded_end().is_some() {
            return None;
        }
        // Holes remove material and cannot enlarge this actual outer source.
        return Some(world_bounds(
            source.dimensions().map(|d| B::new(0., exact(d).hi)),
            frame,
        ));
    }
    let source = probe.polygon_prism()?;
    let mut local = [B::new(f64::INFINITY, f64::NEG_INFINITY); 3];
    for point in &source.prism.rings()[0] {
        for i in 0..2 {
            let coordinate = exact(point[i]);
            local[i].lo = local[i].lo.min(coordinate.lo);
            local[i].hi = local[i].hi.max(coordinate.hi);
        }
    }
    // Linear polygon edges, circular holes and subtracted pockets are inside
    // the retained outer polygon's box. Reverse sweeps retain their real sign.
    let height = exact(source.prism.height());
    local[2] = if source.reverse_sweep {
        B::new(-height.hi, 0.)
    } else {
        B::new(0., height.hi)
    };
    Some(world_bounds(local, source.frame))
}
fn cross(a: [B; 3], b: [B; 3]) -> [B; 3] {
    std::array::from_fn(|i| a[(i + 1) % 3] * b[(i + 2) % 3] - a[(i + 2) % 3] * b[(i + 1) % 3])
}
fn thread_bounds(world: [B; 3], frame: PlaneFrame) -> [[f64; 2]; 3] {
    let [u, v, n] = [frame.u, frame.v, frame.normal()].map(|a| a.map(B::point));
    let rows = [cross(v, n), cross(n, u), cross(u, v)];
    let determinant = (0..3).fold(B::ZERO, |sum, i| sum + u[i] * rows[0][i]);
    let offset: [B; 3] = std::array::from_fn(|i| world[i] - B::point(frame.origin[i]));
    let local =
        rows.map(|row| (0..3).fold(B::ZERO, |sum, i| sum + row[i] * offset[i]) / determinant);
    // The source predicate below validates these enclosures itself. Invalid
    // arithmetic therefore remains unresolved, without a second admission path.
    local.map(|b| [b.lo, b.hi])
}
impl Solid {
    fn disjoint_thread_mating(
        &self,
        other: &Self,
        operation: BooleanOperation,
    ) -> Option<Result<Option<Self>, GeometryError>> {
        let retain_left = match operation {
            BooleanOperation::Intersect => false,
            BooleanOperation::Subtract => true,
            BooleanOperation::Union => return None,
        };
        let (a, af) = self.helical_body()?;
        let (b, bf) = other.helical_body()?;
        if af != bf {
            return None;
        }
        // The profiles use the same actual raw chart, not merely nearly
        // coaxial metadata. A blind bore needs the entire shaft in its void's
        // finite axial interval; a through nut has no stock beyond its thread.
        let (shaft, cavity) = match (&a.form, &b.form) {
            (HelicalForm::Shaft, HelicalForm::Bore(_))
                if a.geometry.spec().length_mm <= b.geometry.spec().length_mm =>
            {
                (&a.geometry, &b.geometry)
            }
            (HelicalForm::Bore(_), HelicalForm::Shaft)
                if b.geometry.spec().length_mm <= a.geometry.spec().length_mm =>
            {
                (&b.geometry, &a.geometry)
            }
            (HelicalForm::Shaft, HelicalForm::Nut(_)) => (&a.geometry, &b.geometry),
            (HelicalForm::Nut(_), HelicalForm::Shaft) => (&b.geometry, &a.geometry),
            _ => return None,
        };
        shaft
            .prove_profile_containment(cavity, Default::default(), None)
            .ok()?;
        // A nonnegative whole-profile gap proves no common INTERIOR. Exact
        // coincident boundary contact cannot produce a three-dimensional solid.
        Some(if retain_left {
            self.share().map(Some)
        } else {
            Ok(None)
        })
    }
    pub(in crate::geometry) fn contained_thread_boolean(
        &self,
        other: &Self,
        operation: BooleanOperation,
        policy: GeometryPolicy,
    ) -> Option<Result<Option<Self>, GeometryError>> {
        if policy.boolean_fuzzy_mm != 0. || self.policy != policy || other.policy != policy {
            return None;
        }
        if let Some(result) = self.disjoint_thread_mating(other, operation) {
            return Some(result);
        }
        let (thread, probe, thread_is_left, body, frame) =
            if let Some((body, frame)) = self.helical_body() {
                (self, other, true, body, frame)
            } else {
                let (body, frame) = other.helical_body()?;
                (other, self, false, body, frame)
            };
        if !matches!(body.form, HelicalForm::Shaft) {
            return None;
        }
        let bounds = thread_bounds(probe_bounds(probe)?, frame);
        let relation = body.geometry.classify_box(bounds).ok()?;
        match (relation, operation, thread_is_left) {
            (BoxClassification::Inside, BooleanOperation::Intersect, _) => {
                Some(probe.share().map(Some))
            }
            (BoxClassification::Inside, BooleanOperation::Union, _) => {
                Some(thread.share().map(Some))
            }
            (BoxClassification::Inside, BooleanOperation::Subtract, false)
            | (BoxClassification::Outside, BooleanOperation::Intersect, _) => Some(Ok(None)),
            (BoxClassification::Outside, BooleanOperation::Subtract, _) => {
                Some(self.share().map(Some))
            }
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests;
