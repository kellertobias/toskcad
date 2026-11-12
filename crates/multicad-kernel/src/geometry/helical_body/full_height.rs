//! Shared exact full-height receiver qualification and source-cap charts.
use super::*;
pub(in crate::geometry) struct FullHeightStock {
    caps: [TopologyElement; 2],
    length: f64,
}
pub(in crate::geometry) fn cancelled() -> GeometryError {
    error(
        "Rust receiving thread treatment",
        GeometryErrorKind::Cancelled,
        "Receiving thread evaluation superseded",
    )
}
fn unsupported(message: &str) -> GeometryError {
    error(
        "Rust receiving thread treatment",
        GeometryErrorKind::Unsupported,
        message,
    )
}
impl FullHeightStock {
    #[allow(clippy::too_many_arguments)]
    pub(in crate::geometry) fn qualify(
        stock: &Solid,
        spec: HelicalThread,
        root_radius_mm: Option<f64>,
        feature: EntityId,
        frame: PlaneFrame,
        policy: GeometryPolicy,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<(Self, PlaneFrame), GeometryError> {
        if interrupt.is_some_and(|f| f()) {
            return Err(cancelled());
        }
        // This existing provider qualifies the literal full tool, source stock,
        // original policy and complete proper-frame map before any replacement.
        let qualified =
            stock.rust_cut_external_thread(spec, root_radius_mm, feature, frame, policy)?;
        let intervals = qualified
            .rust_receiving_thread_intervals()
            .expect("qualified external cylinder cut interval");
        if intervals.requested_tool_length_mm != intervals.stock_length_mm
            || intervals.physical_thread_length_mm != intervals.stock_length_mm
        {
            return Err(unsupported(
                "Receiving thread treatments require the complete literal stock height; partial tails and overruns are not certified",
            ));
        }
        let placed = intervals.frame;
        let snapshot = qualified.topology_snapshot(feature, 0)?;
        let anchors = stock_frame::store_anchors(placed, snapshot.clone(), policy)?;
        let caps = std::array::from_fn(|cap| {
            let mut element = snapshot
                .elements
                .iter()
                .find(|e| {
                    e.key
                        == TopologyKey {
                            kind: TopologyKind::Face,
                            ordinal: cap as u32,
                        }
                })
                .expect("the qualified complete external cut has both finite caps")
                .clone();
            // Only actual stock ancestry transfers. The new material owns its
            // changed finite wires; removed circles/walls acquire no aliases.
            element.origins.retain(|o| o.source != feature);
            element.frames.retain(|a| a.origin.source != feature);
            for anchor in &mut element.frames {
                anchor.frame = anchors[&anchor.origin];
            }
            element
        });
        Ok((
            Self {
                caps,
                length: intervals.stock_length_mm,
            },
            placed,
        ))
    }
    pub(in crate::geometry) fn extend_caps(
        &self,
        frame: PlaneFrame,
        targets: [u32; 2],
        elements: &mut [TopologyElement],
    ) {
        for (source, ordinal) in self.caps.iter().zip(targets) {
            let target = elements
                .iter_mut()
                .find(|e| {
                    e.key
                        == TopologyKey {
                            kind: TopologyKind::Face,
                            ordinal,
                        }
                })
                .expect("published actual receiving cap");
            target.origins.extend(source.origins.iter().cloned());
            target
                .frames
                .extend(source.frames.iter().map(|a| AnchoredFrame {
                    origin: a.origin.clone(),
                    frame: boundary::placed_frame(frame, a.frame.origin_mm, a.frame.u, a.frame.v),
                }));
        }
    }
    pub(in crate::geometry) fn intervals(
        &self,
        frame: PlaneFrame,
        physical_thread_length_mm: f64,
    ) -> ReceivingThreadIntervals {
        ReceivingThreadIntervals {
            frame,
            requested_tool_length_mm: self.length,
            physical_thread_length_mm,
            stock_length_mm: self.length,
        }
    }
}
