//! Actual profiled-screw receiver clearance, separate from nut/mating proofs.
use super::*;
impl Solid {
    pub fn rust_profiled_screw_receiver_clearance<'a>(
        &'a self,
        receiver: &'a Self,
        occurrence: RigidTransform,
        maximum_overlap_volume_mm3: f64,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<source::ProfiledScrewReceiverClearance<'a>, GeometryError> {
        self.policy.validate()?;
        occurrence.validate(self.policy)?;
        if self.policy != receiver.policy || self.policy.boolean_fuzzy_mm != 0. {
            return Err(error(
                "Screw receiving stock",
                GeometryErrorKind::Input,
                "Receiving clearance requires equal actual policies without Boolean fuzz",
            ));
        }
        let unsupported = || {
            error(
                "Screw receiving stock",
                GeometryErrorKind::Unsupported,
                "Receiving clearance requires the actual profiled screw and an original unpierced normal sweep",
            )
        };
        let (screw, frame) = self
            .rust_profiled_screw_geometry()
            .ok_or_else(unsupported)?;
        let (member, placed) = receiver.normal_sweep_body().ok_or_else(unsupported)?;
        let (sweep, reverse) = member.line8_receiving_source().ok_or_else(unsupported)?;
        source::certify_profiled_screw_receiver_clearance(
            source::ProfiledScrewPlacement {
                source: screw,
                frame: source::Frame {
                    origin_mm: frame.origin,
                    axes: [frame.u, frame.v, frame.normal()],
                },
            },
            source::ReceivingSweep {
                source: sweep,
                frame: source::Frame {
                    origin_mm: placed.origin,
                    axes: [placed.u, placed.v, placed.normal()],
                },
                occurrence: source::Frame {
                    origin_mm: occurrence.translation_mm,
                    axes: std::array::from_fn(|j| {
                        std::array::from_fn(|i| occurrence.rotation[i][j])
                    }),
                },
                reverse,
            },
            source::ReceiverClearanceOptions {
                maximum_overlap_volume_mm3,
                maximum_work: 16_000_000,
            },
            &|| interrupt.is_some_and(|f| f()),
        )
        .map_err(|cause| {
            error(
                "Screw receiving stock",
                match cause {
                    source::Error::Invalid => GeometryErrorKind::Input,
                    source::Error::Unsupported => GeometryErrorKind::Unsupported,
                    source::Error::Limit => GeometryErrorKind::Limit,
                    source::Error::Cancelled => GeometryErrorKind::Cancelled,
                },
                format!("Original screw/receiving section proof failed: {cause:?}"),
            )
        })
    }
}
