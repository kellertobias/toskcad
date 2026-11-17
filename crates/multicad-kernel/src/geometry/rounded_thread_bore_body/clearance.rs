//! Certified clearance between a rounded external thread and a coaxial
//! rounded through receiver.
use super::*;
use crate::geometry::coaxial::Coaxial;

impl Solid {
    /// The through receiver has material only within its own axial interval
    /// and outside its cavity profile there; an external thread whose whole
    /// meridian clears that profile at the actual phase shares no interior
    /// with it, wherever along the axis the two overlap.
    pub(in crate::geometry) fn rounded_thread_clearance_boolean(
        &self,
        other: &Self,
        operation: BooleanOperation,
        policy: GeometryPolicy,
    ) -> Option<Result<Option<Self>, GeometryError>> {
        if policy.boolean_fuzzy_mm != 0.
            || self.policy != policy
            || other.policy != policy
            || operation == BooleanOperation::Union
        {
            return None;
        }
        let ((thread, frame), (bore, bore_frame)) =
            match (self.rounded_thread_body(), other.rounded_thread_bore_body()) {
                (Some(thread), Some(bore)) => (thread, bore),
                _ => (
                    other.rounded_thread_body()?,
                    self.rounded_thread_bore_body()?,
                ),
            };
        let external = &thread.geometry;
        let cavity = bore.geometry.thread();
        let span = external.spec().length_mm + cavity.spec().length_mm;
        let placement = Coaxial::of(frame, bore_frame, span)?;
        let [lo, hi] = placement.phase_shift(cavity.helical_rate(), cavity.root_radius_mm())?;
        // The clearance must hold across the whole shift enclosure.
        let slack = 0.5 * (hi - lo);
        let clearance = external
            .prove_clearance(cavity, 0.5 * (lo + hi))
            .ok()
            .filter(|gap| *gap >= placement.radial_mm + slack * cavity.meridian_slope_limit())?;
        debug_assert!(clearance >= 0.);
        Some(match operation {
            BooleanOperation::Subtract => self.share().map(Some),
            _ => Ok(None),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use multicad_rust_kernel::rounded_thread::Spec;

    fn spec(major: f64) -> Spec {
        let angle = 27.5_f64.to_radians();
        Spec {
            major_radius_mm: major,
            pitch_mm: 1.,
            length_mm: 4.,
            depth_mm: 1. / (2. * angle.tan()) - 0.274 * (1. / angle.sin() - 1.),
            crest_radius_mm: 0.137,
            root_radius_mm: 0.137,
            flank_angle_deg: 27.5,
            starts: 1,
            left_hand: false,
        }
    }

    fn p() -> GeometryPolicy {
        GeometryPolicy::default()
    }

    #[test]
    fn coaxial_clear_pairs_have_no_common_interior_in_either_order() {
        let frame = PlaneFrame {
            origin: [3., -2., 7.],
            u: [0., 1., 0.],
            v: [0., 0., 1.],
        };
        let shaft = Solid::rust_rounded_thread(spec(5.), EntityId::new(), frame, p()).unwrap();
        // The same chart, then one advanced by a whole pitch along the axis.
        for origin in [frame.origin, [4., -2., 7.]] {
            let nut = Solid::rust_rounded_thread_bore(
                spec(5.05),
                8.,
                EntityId::new(),
                PlaneFrame { origin, ..frame },
                p(),
            )
            .unwrap();
            for (a, b) in [(&shaft, &nut), (&nut, &shaft)] {
                assert!(
                    a.boolean(b, BooleanOperation::Intersect, p())
                        .unwrap()
                        .is_none()
                );
                let kept = a
                    .boolean(b, BooleanOperation::Subtract, p())
                    .unwrap()
                    .unwrap();
                assert!((kept.volume_mm3().unwrap() - a.volume_mm3().unwrap()).abs() < 1e-9);
                assert!(
                    a.rounded_thread_clearance_boolean(b, BooleanOperation::Union, p())
                        .is_none()
                );
            }
        }
    }

    #[test]
    fn misphased_tight_or_unrelated_pairs_are_not_claimed() {
        let frame = PlaneFrame::default();
        let shaft = Solid::rust_rounded_thread(spec(5.), EntityId::new(), frame, p()).unwrap();
        let bore = |major, origin: [f64; 3]| {
            Solid::rust_rounded_thread_bore(
                spec(major),
                8.,
                EntityId::new(),
                PlaneFrame { origin, ..frame },
                p(),
            )
            .unwrap()
        };
        // Half a pitch along the axis puts crest against crest.
        let misphased = bore(5.05, [0., 0., 0.5]);
        assert!(
            shaft
                .rounded_thread_clearance_boolean(&misphased, BooleanOperation::Intersect, p())
                .is_none()
        );
        // A cavity smaller than the external thread overlaps it.
        let tight = bore(4.95, [0., 0., 0.]);
        assert!(
            shaft
                .rounded_thread_clearance_boolean(&tight, BooleanOperation::Intersect, p())
                .is_none()
        );
        // Off-axis, fuzzy, or a pair without a receiver.
        let offset = bore(5.05, [0.01, 0., 0.]);
        assert!(
            shaft
                .rounded_thread_clearance_boolean(&offset, BooleanOperation::Intersect, p())
                .is_none()
        );
        let fuzzy = GeometryPolicy {
            boolean_fuzzy_mm: 1e-6,
            ..p()
        };
        let clear = bore(5.05, [0., 0., 0.]);
        assert!(
            shaft
                .rounded_thread_clearance_boolean(&clear, BooleanOperation::Intersect, fuzzy)
                .is_none()
        );
        assert!(
            shaft
                .rounded_thread_clearance_boolean(&shaft, BooleanOperation::Intersect, p())
                .is_none()
        );
    }
}
