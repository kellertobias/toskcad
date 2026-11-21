//! Round hardware against a receiving member, including one cut by a
//! transverse stepped bore: the tool bands removed that member's material.
use super::*;
impl Solid {
    /// Bound occupied overlap of this round hardware with ONE receiving
    /// member in its original occurrence pose. The hardware is an actual
    /// profiled screw or a cross-bored round barrel, taken as its original
    /// uncut cylinder; the receiver is an original unpierced sweep or one
    /// cut by a transverse stepped bore.
    pub fn rust_round_hardware_receiver_clearance(
        &self,
        receiver: &Self,
        occurrence: RigidTransform,
        maximum_overlap_volume_mm3: f64,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<source::BoredReceiverClearance, GeometryError> {
        const OPERATION: &str = "Round receiving stock";
        self.policy.validate()?;
        occurrence.validate(self.policy)?;
        if self.policy != receiver.policy || self.policy.boolean_fuzzy_mm != 0. {
            return Err(error(
                OPERATION,
                GeometryErrorKind::Input,
                "Receiving clearance requires equal actual policies without Boolean fuzz",
            ));
        }
        let unsupported = || {
            error(
                OPERATION,
                GeometryErrorKind::Unsupported,
                "Round receiving clearance requires an actual profiled screw or cross-bored barrel \
                 and an original or transversely bored unpierced normal sweep",
            )
        };
        let failure = |cause: source::Error| {
            error(
                OPERATION,
                match cause {
                    source::Error::Invalid => GeometryErrorKind::Input,
                    source::Error::Unsupported => GeometryErrorKind::Unsupported,
                    source::Error::Limit => GeometryErrorKind::Limit,
                    source::Error::Cancelled => GeometryErrorKind::Cancelled,
                },
                format!("Round hardware/receiving section proof failed: {cause:?}"),
            )
        };
        let (pieces, frame) = if let Some((screw, frame)) = self.rust_profiled_screw_geometry() {
            (
                source::profiled_screw_pieces(screw).map_err(failure)?,
                source::Frame {
                    origin_mm: frame.origin,
                    axes: [frame.u, frame.v, frame.normal()],
                },
            )
        } else if let Some(record) = self.transverse_record()
            && let Some((_, frame)) = record.source.axial_cylinder_frame()
            && let Some((origin, _, radius, height)) = record.source.rust_cylinder_axis()
        {
            (
                vec![source::HardwarePiece {
                    radius_mm: radius,
                    span_mm: [0., height],
                }],
                source::Frame {
                    origin_mm: origin,
                    axes: [frame.u, frame.v, frame.normal()],
                },
            )
        } else {
            return Err(unsupported());
        };
        let (member, placed, record) = receiver.receiving_sweep().ok_or_else(unsupported)?;
        let (sweep, reverse) = member.line8_receiving_source().ok_or_else(unsupported)?;
        let bands: Vec<source::ToolBand> = record
            .map(|record| {
                record
                    .bands
                    .iter()
                    .map(|&[start_mm, end_mm, radius_mm]| source::ToolBand {
                        start_mm,
                        end_mm,
                        radius_mm,
                    })
                    .collect()
            })
            .unwrap_or_default();
        let bore = record.map(|record| source::ReceivingBore {
            origin_mm: record.origin_mm,
            axis: record.axis,
            bands: &bands,
        });
        source::certify_round_pieces_clearance(
            &pieces,
            frame,
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
            bore,
            source::ReceiverClearanceOptions {
                maximum_overlap_volume_mm3,
                maximum_work: 16_000_000,
            },
            &|| interrupt.is_some_and(|f| f()),
        )
        .map_err(failure)
    }
}
