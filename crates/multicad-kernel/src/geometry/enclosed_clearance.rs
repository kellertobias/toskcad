//! Intersections and cuts of axial bodies proved empty by coaxial cylinders.
//!
//! A placed fastener or a polygonal revolution lies inside a stack of
//! closed cylindrical slabs about its own axis: between two consecutive
//! meridian heights, every boundary point lies on a non-horizontal surface
//! or meridian edge crossing that open interval, whose endpoints bound its
//! radius, so the bounded cross-section lies in the disk of the largest such
//! endpoint radius. A shoulder at a slab height bounds only a plane of zero
//! volume, so the closed slabs cover all material up to sets that cannot
//! carry a common interior. When the other operand
//! shares no interior with any slab, by the ordinary exact Boolean of that
//! operand with an exact cylinder, it shares none with the axial body
//! either. The intersection is then empty and a cut leaves the target
//! unchanged. Any overlapping or refused slab returns `None`, leaving the
//! pair to the remaining providers and their refusal.
use super::*;
use multicad_rust_kernel::cosmetic_fastener::Surface;

/// One closed slab `[low, high]` of radius `radius` about the axis.
#[derive(Clone, Copy, Debug, PartialEq)]
struct Slab {
    low: f64,
    high: f64,
    radius: f64,
}

/// Slabs covering every edge `(r₀, z₀)–(r₁, z₁)` of a closed meridian.
fn slabs(edges: &[([f64; 2], [f64; 2])]) -> Option<Vec<Slab>> {
    let mut heights: Vec<f64> = edges.iter().flat_map(|(a, b)| [a[1], b[1]]).collect();
    if heights.iter().any(|z| !z.is_finite())
        || edges
            .iter()
            .any(|(a, b)| !(a[0].is_finite() && b[0].is_finite() && a[0] >= 0. && b[0] >= 0.))
    {
        return None;
    }
    heights.sort_by(f64::total_cmp);
    heights.dedup();
    let slabs: Vec<Slab> = heights
        .windows(2)
        .map(|pair| {
            let [low, high] = [pair[0], pair[1]];
            let radius = edges
                .iter()
                .filter(|(a, b)| a[1].min(b[1]) < high && a[1].max(b[1]) > low)
                .map(|(a, b)| a[0].max(b[0]))
                .fold(0., f64::max);
            Slab { low, high, radius }
        })
        .collect();
    (!slabs.is_empty()).then_some(slabs)
}

fn fastener_edges(surfaces: impl Iterator<Item = Surface>) -> Vec<([f64; 2], [f64; 2])> {
    surfaces
        .map(|surface| {
            let [a, b] = surface.rings;
            (
                [a.section.outer_radius(), a.z],
                [b.section.outer_radius(), b.z],
            )
        })
        .collect()
}

impl Solid {
    /// The axis frame and covering slabs of an axial body.
    fn enclosing_slabs(&self) -> Option<(PlaneFrame, Vec<Slab>)> {
        match &self.body {
            SolidBody::CosmeticFastener(body, frame) => {
                Some((*frame, slabs(&fastener_edges(body.geometry.surfaces()))?))
            }
            SolidBody::Revolved(body, frame) => {
                // Holes lie inside the outer meridian, which alone bounds
                // the material.
                let outer = body.geometry.rings().first()?;
                let edges: Vec<_> = outer
                    .iter()
                    .zip(outer.iter().cycle().skip(1))
                    .map(|(&a, &b)| (a, b))
                    .collect();
                Some((*frame, slabs(&edges)?))
            }
            _ => None,
        }
    }

    /// Exact cylinders that together enclose this axial body, each slightly
    /// widened and lengthened by a rounding margin.
    fn enclosing_cylinders(&self, policy: GeometryPolicy) -> Option<Vec<Self>> {
        let (frame, slabs) = self.enclosing_slabs()?;
        let normal = frame.normal();
        let scale = frame
            .origin
            .into_iter()
            .chain(
                slabs
                    .iter()
                    .flat_map(|slab| [slab.low, slab.high, slab.radius]),
            )
            .map(f64::abs)
            .fold(1., f64::max);
        let margin = 64. * f64::EPSILON * scale;
        slabs
            .into_iter()
            .filter(|slab| slab.radius > 0.)
            .map(|slab| {
                let low = slab.low - margin;
                let origin = std::array::from_fn(|i| frame.origin[i] + normal[i] * low);
                construction::cylinder(
                    PlaneFrame { origin, ..frame },
                    slab.radius + margin,
                    slab.high + margin - low,
                    EntityId::new(),
                    policy,
                )
                .ok()
            })
            .collect()
    }

    pub(in crate::geometry) fn enclosed_clearance_boolean(
        &self,
        other: &Self,
        operation: BooleanOperation,
        policy: GeometryPolicy,
    ) -> Option<Result<Option<Self>, GeometryError>> {
        if operation == BooleanOperation::Union
            || policy.boolean_fuzzy_mm != 0.
            || self.policy != policy
            || other.policy != policy
        {
            return None;
        }
        let (axial, probe) = match (self.enclosing_slabs(), other.enclosing_slabs()) {
            (Some(_), _) => (self, other),
            (None, Some(_)) => (other, self),
            (None, None) => return None,
        };
        for cylinder in axial.enclosing_cylinders(policy)? {
            match probe.boolean_geometry(&cylinder, BooleanOperation::Intersect, policy) {
                Ok(None) => {}
                Ok(Some(_)) | Err(_) => return None,
            }
        }
        Some(match operation {
            BooleanOperation::Intersect => Ok(None),
            _ => self.share().map(Some),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn slabs_take_the_largest_endpoint_radius_of_every_edge_they_meet() {
        // A stepped meridian: radius 3 up to z = 2, then radius 1 to z = 5.
        let profile = [[0., 0.], [3., 0.], [3., 2.], [1., 2.], [1., 5.], [0., 5.]];
        let edges: Vec<_> = profile
            .iter()
            .zip(profile.iter().cycle().skip(1))
            .map(|(&a, &b)| (a, b))
            .collect();
        assert_eq!(
            slabs(&edges).unwrap(),
            [
                Slab {
                    low: 0.,
                    high: 2.,
                    radius: 3.
                },
                // The shoulder at z = 2 bounds no volume above it.
                Slab {
                    low: 2.,
                    high: 5.,
                    radius: 1.
                },
            ]
        );
        assert!(slabs(&[([f64::NAN, 0.], [1., 1.])]).is_none());
        assert!(slabs(&[([1., 0.], [1., 0.])]).is_none());
    }
}
