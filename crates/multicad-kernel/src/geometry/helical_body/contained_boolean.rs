//! Whole-source containment/disjointness only. No cut surface is approximated
//! or manufactured by this path: a proved result shares an unchanged operand.
use super::*;
use crate::geometry::coaxial::Coaxial;
use crate::geometry::internal_thread_treatment_body::interval::Interval as B;
use multicad_rust_kernel::helical_thread::{BoxClassification, ProfileContainmentError};
mod common;
mod headed;

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
pub(in crate::geometry) fn probe_bounds(probe: &Solid) -> Option<[B; 3]> {
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
    chart_bounds(world, frame.origin, [frame.u, frame.v, frame.normal()])
}
pub(in crate::geometry) fn chart_bounds(
    world: [B; 3],
    origin: Point,
    axes: [Point; 3],
) -> [[f64; 2]; 3] {
    let [u, v, n] = axes.map(|a| a.map(B::point));
    let rows = [cross(v, n), cross(n, u), cross(u, v)];
    let determinant = (0..3).fold(B::ZERO, |sum, i| sum + u[i] * rows[0][i]);
    let offset: [B; 3] = std::array::from_fn(|i| world[i] - B::point(origin[i]));
    let local =
        rows.map(|row| (0..3).fold(B::ZERO, |sum, i| sum + row[i] * offset[i]) / determinant);
    // The source predicate below validates these enclosures itself. Invalid
    // arithmetic therefore remains unresolved, without a second admission path.
    local.map(|b| [b.lo, b.hi])
}
fn bore_relation(
    body: &HelicalBody,
    stock: &BoreStock,
    bounds: [[f64; 2]; 3],
) -> Option<BoxClassification> {
    use BoxClassification::*;
    // The cavity is the actual finite source (including a retained prefix).
    // A whole box inside it contains no occupied bore material.
    let cavity = body.geometry.classify_box(bounds).ok()?;
    if cavity == Inside {
        return Some(Outside);
    }
    let (origin, axes) = stock.frame_in_thread.source_chart();
    let [x, y, z] = chart_bounds(bounds.map(|[lo, hi]| B::new(lo, hi)), origin, axes)
        .map(|[lo, hi]| B::new(lo, hi));
    let radius_squared = B::point(stock.radius()).square();
    let rho_squared = x.square() + y.square();
    if !rho_squared.valid() || !z.valid() {
        return None;
    }
    if rho_squared.lo > radius_squared.hi || z.hi < 0. || z.lo > stock.length() {
        Some(Outside)
    } else if cavity == Outside
        && rho_squared.hi < radius_squared.lo
        && z.lo > 0.
        && z.hi < stock.length()
    {
        Some(Inside)
    } else {
        Some(UnresolvedBoundary)
    }
}
/// Which axial interval must contain the other for a profile proof to
/// cover every pair of material points.
#[derive(Clone, Copy)]
enum Span {
    ShaftInCavity,
    CavityInShaft,
    Unbounded,
}
impl Span {
    /// `offset` is the cavity start along the shaft axis, within `error`.
    fn admits(self, offset: f64, error: f64, shaft: f64, cavity: f64) -> bool {
        match self {
            Span::ShaftInCavity => offset <= error && offset + cavity >= shaft - error,
            Span::CavityInShaft => offset >= -error && offset + cavity <= shaft + error,
            Span::Unbounded => true,
        }
    }
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
        // The external thread in its own chart, then the cavity, and which
        // finite axial interval must contain the other. A bore needs the
        // entire shaft in its void's interval; a through nut has no stock
        // beyond its thread; a nut on a headed screw must sit wholly on the
        // threaded length, away from the shank and head.
        let (shaft, frame, cavity, cavity_frame, rule, nut) = match (&a.form, &b.form) {
            (HelicalForm::Shaft, HelicalForm::Bore(_)) => {
                (&a.geometry, af, &b.geometry, bf, Span::ShaftInCavity, None)
            }
            (HelicalForm::Bore(_), HelicalForm::Shaft) => {
                (&b.geometry, bf, &a.geometry, af, Span::ShaftInCavity, None)
            }
            (HelicalForm::Shaft, HelicalForm::Nut(nut)) => (
                &a.geometry,
                af,
                &b.geometry,
                bf,
                Span::Unbounded,
                Some(nut.across_flats_mm()),
            ),
            (HelicalForm::Nut(nut), HelicalForm::Shaft) => (
                &b.geometry,
                bf,
                &a.geometry,
                af,
                Span::Unbounded,
                Some(nut.across_flats_mm()),
            ),
            (HelicalForm::Screw(screw), HelicalForm::Nut(nut)) => (
                screw.thread(),
                af,
                nut.thread(),
                bf,
                Span::CavityInShaft,
                Some(nut.across_flats_mm()),
            ),
            (HelicalForm::Nut(nut), HelicalForm::Screw(screw)) => (
                screw.thread(),
                bf,
                nut.thread(),
                af,
                Span::CavityInShaft,
                Some(nut.across_flats_mm()),
            ),
            _ => return None,
        };
        let [shaft_length, cavity_length] = [shaft, cavity].map(|thread| thread.spec().length_mm);
        let (proof, shift, axial) = if frame == cavity_frame {
            if !rule.admits(0., 0., shaft_length, cavity_length) {
                return None;
            }
            let proof = shaft
                .prove_profile_containment(cavity, Default::default(), None)
                .map(|_| ());
            (proof, [0.; 2], 0.)
        } else {
            // A cavity rotated about or moved along the common axis is the
            // same profile law at a certified phase shift.
            let placement = Coaxial::of(frame, cavity_frame, shaft_length + cavity_length)?;
            if !rule.admits(
                placement.axial_mm,
                placement.error_mm,
                shaft_length,
                cavity_length,
            ) {
                return None;
            }
            let shift = placement.phase_shift(cavity.helical_rate(), cavity.root_radius_mm())?;
            let proof = shaft
                .prove_shifted_profile_containment(cavity, shift, Default::default(), None)
                .ok()
                .filter(|proof| proof.minimum_gap_mm() >= placement.radial_mm)
                .map(|_| ())
                .ok_or(ProfileContainmentError::Unresolved);
            (proof, shift, placement.axial_mm)
        };
        if proof.is_err() {
            // Overlapping laws: only an intersection with a through nut,
            // whose material is its hex less the threaded void over its
            // own length, has a constructed common material. The external
            // crest must stay inside the hex's inscribed circle.
            let across_flats = nut?;
            if retain_left || shaft.spec().major_radius_mm >= across_flats / 2. {
                return None;
            }
            let span = [axial.max(0.), (axial + cavity_length).min(shaft_length)];
            return Self::coaxial_thread_common(
                shaft,
                frame,
                cavity,
                shift[0].midpoint(shift[1]),
                span,
                self.policy,
            );
        }
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
        let (probe, thread_is_left, body, frame) = if let Some((body, frame)) = self.helical_body()
        {
            (other, true, body, frame)
        } else {
            let (body, frame) = other.helical_body()?;
            (self, false, body, frame)
        };
        let bounds = thread_bounds(probe_bounds(probe)?, frame);
        let relation = match &body.form {
            HelicalForm::Shaft => body.geometry.classify_box(bounds).ok()?,
            HelicalForm::Bore(stock) => bore_relation(body, stock, bounds)?,
            HelicalForm::Screw(source) => headed::screw_relation(source, bounds),
            HelicalForm::Nut(source) => headed::nut_relation(source, bounds),
            _ => return None,
        };
        whole_operand_result(self, other, thread_is_left, relation, operation)
    }
}

/// Only an earned whole-probe relation may reuse either original operand.
pub(in crate::geometry) fn whole_operand_result(
    lhs: &Solid,
    rhs: &Solid,
    source_is_left: bool,
    relation: BoxClassification,
    operation: BooleanOperation,
) -> Option<Result<Option<Solid>, GeometryError>> {
    let (source, probe) = if source_is_left {
        (lhs, rhs)
    } else {
        (rhs, lhs)
    };
    match (relation, operation, source_is_left) {
        (BoxClassification::Inside, BooleanOperation::Intersect, _) => {
            Some(probe.share().map(Some))
        }
        (BoxClassification::Inside, BooleanOperation::Union, _) => Some(source.share().map(Some)),
        (BoxClassification::Inside, BooleanOperation::Subtract, false)
        | (BoxClassification::Outside, BooleanOperation::Intersect, _) => Some(Ok(None)),
        (BoxClassification::Outside, BooleanOperation::Subtract, _) => Some(lhs.share().map(Some)),
        _ => None,
    }
}

#[cfg(test)]
mod tests;
