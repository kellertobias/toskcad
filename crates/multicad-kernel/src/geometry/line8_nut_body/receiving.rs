//! Receiving stock clearance is independent of screw/nut mating.
use super::*;
mod round;
mod screw;
use multicad_rust_kernel::wood_thread_contact as source;
impl Solid {
    /// The original unpierced sweep a receiving member keeps, its placement
    /// and the transverse tool that cut it, if any. A cut member's material
    /// is a subset of its original sweep, so original-sweep clearance bounds
    /// remain bounds for it.
    pub(crate) fn receiving_sweep(
        &self,
    ) -> Option<(
        &super::super::normal_sweep_body::NormalSweepBody,
        PlaneFrame,
        Option<&super::super::imported_brep_body::TransverseRecord>,
    )> {
        if let Some((body, frame)) = self.normal_sweep_body() {
            return Some((body, frame, None));
        }
        let record = self.transverse_record()?;
        let (body, frame) = record.source.normal_sweep_body()?;
        Some((body, frame, Some(record)))
    }

    /// Bound occupied overlap with ONE actual receiving member in its original
    /// occurrence pose. Call for every receiving member before admitting hardware.
    pub fn rust_line8_receiver_clearance<'a>(
        &'a self,
        receiver: &'a Self,
        occurrence: RigidTransform,
        maximum_overlap_volume_mm3: f64,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<source::Line8ReceiverClearance<'a>, GeometryError> {
        self.policy.validate()?;
        occurrence.validate(self.policy)?;
        if self.policy != receiver.policy || self.policy.boolean_fuzzy_mm != 0. {
            return Err(error(
                "Line8 receiving stock",
                GeometryErrorKind::Input,
                "Receiving clearance requires equal actual source policies without Boolean fuzz",
            ));
        }
        let unsupported = || {
            error(
                "Line8 receiving stock",
                GeometryErrorKind::Unsupported,
                "Receiving clearance requires actual retained Line8 stock and an original unpierced normal sweep",
            )
        };
        let SolidBody::Line8Nut(nut) = &self.body else {
            return Err(unsupported());
        };
        let (member, placed, _) = receiver.receiving_sweep().ok_or_else(unsupported)?;
        let (sweep, reverse) = member.line8_receiving_source().ok_or_else(unsupported)?;
        let literal = |matrix: [[f64; 3]; 3], origin_mm: Point| source::Frame {
            origin_mm,
            axes: std::array::from_fn(|j| std::array::from_fn(|i| matrix[i][j])),
        };
        source::certify_line8_receiver_clearance(
            source::ClockedLine8Placement {
                source: &nut.source,
                stock_frame: literal(nut.stock.matrix(), nut.stock.translation_mm()),
                outer_frame: literal(nut.outer.rotation, nut.outer.translation_mm),
                minimum_wall_mm: nut.requested_wall_mm,
            },
            source::ReceivingSweep {
                source: sweep,
                frame: source::Frame {
                    origin_mm: placed.origin,
                    axes: [placed.u, placed.v, placed.normal()],
                },
                occurrence: literal(occurrence.rotation, occurrence.translation_mm),
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
                "Line8 receiving stock",
                match cause {
                    source::Error::Invalid => GeometryErrorKind::Input,
                    source::Error::Cancelled => GeometryErrorKind::Cancelled,
                    source::Error::Limit => GeometryErrorKind::Limit,
                    source::Error::Unsupported => GeometryErrorKind::Unsupported,
                },
                format!("Original receiving section/occurrence proof failed: {cause:?}"),
            )
        })
    }
}
#[cfg(test)]
mod tests;
