//! Bounds original seated-stock cells against each actual receiving member.
use super::*;
use multicad_rust_kernel::wood_thread_contact as contact;
impl Solid {
    pub fn rust_seated_stock_receiver_clearance<'a>(
        &'a self,
        receiver: &'a Self,
        occurrence: RigidTransform,
        maximum_overlap_volume_mm3: f64,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<contact::SeatedStockReceiverClearance<'a>, GeometryError> {
        self.policy.validate()?;
        occurrence.validate(self.policy)?;
        if self.policy != receiver.policy || self.policy.boolean_fuzzy_mm != 0. {
            return Err(input(
                "Receiving clearance requires the original exact shared policies",
            ));
        }
        let unsupported =
            || unsupported("Expected actual seated stock and original unpierced receiving sweep");
        let (body, frame) = self.orthogonal_seat_body().ok_or_else(unsupported)?;
        let (member, placed) = receiver.normal_sweep_body().ok_or_else(unsupported)?;
        let (sweep, reverse) = member.line8_receiving_source().ok_or_else(unsupported)?;
        contact::certify_seated_stock_receiver_clearance(
            contact::SeatedStockPlacement {
                source: &body.source,
                frame: contact::Frame {
                    origin_mm: frame.origin,
                    axes: [frame.u, frame.v, frame.normal()],
                },
            },
            contact::ReceivingSweep {
                source: sweep,
                frame: contact::Frame {
                    origin_mm: placed.origin,
                    axes: [placed.u, placed.v, placed.normal()],
                },
                occurrence: contact::Frame {
                    origin_mm: occurrence.translation_mm,
                    axes: std::array::from_fn(|j| {
                        std::array::from_fn(|i| occurrence.rotation[i][j])
                    }),
                },
                reverse,
            },
            contact::ReceiverClearanceOptions {
                maximum_overlap_volume_mm3,
                maximum_work: 16_000_000,
            },
            &|| interrupt.is_some_and(|f| f()),
        )
        .map_err(|cause| {
            error(
                "Seated-stock receiving clearance",
                match cause {
                    contact::Error::Invalid => GeometryErrorKind::Input,
                    contact::Error::Unsupported => GeometryErrorKind::Unsupported,
                    contact::Error::Limit => GeometryErrorKind::Limit,
                    contact::Error::Cancelled => GeometryErrorKind::Cancelled,
                },
                format!("Original stock cell/receiver support proof failed: {cause:?}"),
            )
        })
    }
}
