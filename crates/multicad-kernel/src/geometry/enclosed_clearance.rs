//! Intersections and cuts of axial bodies proved empty by coaxial cylinders.
//!
//! A placed fastener, a modeled machine screw or a polygonal revolution
//! lies inside a stack of
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

/// A modeled machine screw in its thread chart: the thread within its
/// crest cylinder over `[0, L]`, the plain shank over `[-s, 0]` and the head
/// within its round rim or hex corners over `[-s - h, -s]`. A socket only
/// removes head material.
fn screw_slabs(screw: &multicad_rust_kernel::headed_screw::Screw) -> Option<Vec<Slab>> {
    use multicad_rust_kernel::headed_screw::ScrewHead;
    let thread = screw.thread().spec();
    let spec = screw.spec();
    let head = match spec.head {
        ScrewHead::SocketCap { radius_mm, .. } => radius_mm,
        ScrewHead::ExternalHex {
            across_flats_mm, ..
        } => across_flats_mm / 3_f64.sqrt(),
    };
    let seat = screw.head_seat_z();
    let mut slabs = vec![Slab {
        low: screw.head_top_z(),
        high: seat,
        radius: head,
    }];
    if spec.shank_length_mm > 0. {
        slabs.push(Slab {
            low: seat,
            high: 0.,
            radius: spec.shank_radius_mm,
        });
    }
    slabs.push(Slab {
        low: 0.,
        high: thread.length_mm,
        radius: thread.major_radius_mm,
    });
    slabs
        .iter()
        .all(|slab| {
            [slab.low, slab.high, slab.radius]
                .iter()
                .all(|v| v.is_finite())
                && slab.low < slab.high
                && slab.radius > 0.
        })
        .then_some(slabs)
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
            _ => {
                let (screw, frame) = self.rust_machine_screw_geometry()?;
                Some((frame, screw_slabs(screw)?))
            }
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
    fn machine_screws_lie_in_their_thread_shank_and_head_cylinders() {
        use multicad_rust_kernel::headed_screw::{Screw, ScrewHead, ScrewSpec};
        use multicad_rust_kernel::helical_thread::{Spec, Thread};
        let thread = Thread::new(Spec {
            major_radius_mm: 1.5,
            pitch_mm: 0.5,
            length_mm: 8.,
            depth_mm: 0.27,
            crest_width_mm: 0.06,
            lower_flank_deg: 30.,
            upper_flank_deg: 30.,
            starts: 1,
            left_hand: false,
            root_radius_mm: None,
        })
        .unwrap();
        let socket = ScrewHead::SocketCap {
            radius_mm: 2.75,
            height_mm: 3.,
            socket_across_flats_mm: 2.5,
            socket_depth_mm: 1.3,
        };
        let hex = ScrewHead::ExternalHex {
            across_flats_mm: 5.5,
            height_mm: 2.,
        };
        let cap = Screw::new(
            thread,
            ScrewSpec {
                head: socket,
                shank_length_mm: 2.,
                shank_radius_mm: 1.5,
            },
        )
        .unwrap();
        assert_eq!(
            screw_slabs(&cap).unwrap(),
            [
                Slab {
                    low: -5.,
                    high: -2.,
                    radius: 2.75
                },
                Slab {
                    low: -2.,
                    high: 0.,
                    radius: 1.5
                },
                Slab {
                    low: 0.,
                    high: 8.,
                    radius: 1.5
                },
            ]
        );
        // Without a shank the head sits on the thread; hex corners bound it.
        let flush = Screw::new(
            thread,
            ScrewSpec {
                head: hex,
                shank_length_mm: 0.,
                shank_radius_mm: 1.5,
            },
        )
        .unwrap();
        let slabs = screw_slabs(&flush).unwrap();
        assert_eq!(slabs.len(), 2);
        assert_eq!([slabs[0].low, slabs[0].high], [-2., 0.]);
        assert!((slabs[0].radius - 5.5 / 3_f64.sqrt()).abs() < 1e-15);
        assert_eq!(slabs[1].radius, 1.5);
    }

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
