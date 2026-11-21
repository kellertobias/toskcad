//! The actual common material of an external thread and a coaxial threaded
//! cavity whose profiles overlap.
//!
//! Both sources are invariant under the same screw motion, so their common
//! material over the shared axial span is the screw sweep of the meridian
//! lenses where the cavity law lies inside the external law. Every section
//! of that sweep is one fixed section turned about the axis at the helical
//! rate, so each lens is the twisted extrusion of its exact section, whose
//! pieces are affine-radius spirals. The result is retained through its
//! certified exchange boundary and read back with the generic reader.
use super::*;
use multicad_rust_kernel::step::{Spiral, write_twisted_spirals};

impl Solid {
    /// `None` leaves the pair to the remaining providers when the lenses
    /// are unresolved, empty or outside the encoder's budget.
    pub(super) fn coaxial_thread_common(
        shaft: &Thread,
        frame: PlaneFrame,
        cavity: &Thread,
        shift_mm: f64,
        span_mm: [f64; 2],
        policy: GeometryPolicy,
    ) -> Option<Result<Option<Self>, GeometryError>> {
        let tolerance = policy.linear_tolerance_mm;
        let [z0, z1] = span_mm;
        if z1 - z0 <= policy.minimum_feature_mm {
            return None;
        }
        let lenses = shaft.common_lenses(cavity, shift_mm, tolerance).ok()?;
        if lenses.is_empty() {
            return None;
        }
        let rate = shaft.helical_rate();
        // On the section at z0 a phase lies at θ = (z0 - phase) / rate.
        let loops: Vec<Vec<Spiral>> = lenses
            .iter()
            .map(|lens| {
                lens.pieces
                    .iter()
                    .map(|piece| Spiral {
                        angle_rad: piece.phase_mm.map(|phase| (z0 - phase) / rate),
                        radius_mm: piece.radius_mm,
                    })
                    .collect()
            })
            .collect();
        let text = write_twisted_spirals(
            &loops,
            span_mm,
            rate.recip(),
            frame.origin,
            frame.u,
            frame.v,
            tolerance,
        )
        .ok()?;
        Some(Self::read_thread_common(&text, policy))
    }

    fn read_thread_common(
        text: &str,
        policy: GeometryPolicy,
    ) -> Result<Option<Self>, GeometryError> {
        let invalid = |message: &str| {
            error(
                "Boolean",
                GeometryErrorKind::InvalidResult,
                format!("Rust thread interference {message}"),
            )
        };
        let solids = multicad_rust_kernel::step_read::brep::read_solids(
            text,
            policy.linear_tolerance_mm,
            None,
        )
        .map_err(|cause| invalid(&format!("boundary did not read back: {cause:?}")))?;
        let feature = EntityId::new();
        let members = solids
            .into_iter()
            .map(|solid| Self::from_assembled_brep(solid, feature, policy))
            .collect::<Result<Vec<_>, _>>()?;
        match members.len() {
            0 => Err(invalid("produced no solid")),
            1 => Ok(members.into_iter().next()),
            _ => {
                let members: Vec<&Self> = members.iter().collect();
                RustCompound::assemble(&members, feature, policy, None).map(Some)
            }
        }
    }
}
