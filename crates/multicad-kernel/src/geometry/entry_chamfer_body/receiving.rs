//! A full-height stock cut reuses its qualified chart and actual entry material.
use super::super::helical_body::{ReceivingThreadIntervals, stock_frame};
use super::*;

pub(super) struct ReceivingStock {
    caps: [TopologyElement; 2],
    length: f64,
}
fn unsupported(message: &str) -> GeometryError {
    error(
        "Rust external receiving entry",
        GeometryErrorKind::Unsupported,
        message,
    )
}
impl Solid {
    /// Cut entry-chamfered material over the complete height of a coaxial
    /// cylinder. Partial cuts and overruns need distinct finite-tail proofs.
    #[allow(clippy::too_many_arguments)]
    pub fn rust_cut_external_thread_entry_chamfer_with_interrupt(
        &self,
        spec: HelicalThread,
        root_radius_mm: Option<f64>,
        radial_depth_mm: f64,
        terminal: Option<TerminalIntent>,
        feature: EntityId,
        frame: PlaneFrame,
        policy: GeometryPolicy,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, GeometryError> {
        if interrupt.is_some_and(|f| f()) {
            return Err(failure(source::Error::Cancelled));
        }
        // This existing provider qualifies the literal full tool, source stock,
        // original policy and complete proper-frame map before any replacement.
        let qualified =
            self.rust_cut_external_thread(spec, root_radius_mm, feature, frame, policy)?;
        let intervals = qualified
            .rust_receiving_thread_intervals()
            .expect("qualified external cylinder cut interval");
        if intervals.requested_tool_length_mm != intervals.stock_length_mm
            || intervals.physical_thread_length_mm != intervals.stock_length_mm
        {
            return Err(unsupported(
                "Entry-chamfered receiving cuts require the complete literal stock height; partial tails and overruns are not certified",
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
        Self::make_entry_chamfer(
            spec,
            root_radius_mm,
            radial_depth_mm,
            terminal,
            feature,
            placed,
            policy,
            interrupt,
            Some(ReceivingStock {
                caps,
                length: intervals.stock_length_mm,
            }),
        )
    }
}
impl ReceivingStock {
    pub(super) fn extend_topology(
        &self,
        body: &Chamfer,
        frame: PlaneFrame,
        elements: &mut [TopologyElement],
    ) {
        use source::Face;
        for (index, source) in self.caps.iter().enumerate() {
            let role = if index == 0 {
                Face::EntryCap
            } else if body.terminal_relief().is_some() {
                Face::FarCap
            } else {
                Face::EndCap
            };
            let ordinal = (0..body.face_count())
                .find(|&i| body.face(i) == Some(role))
                .expect("retained source cap");
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
                    frame: ReferenceFrame {
                        origin_mm: position(frame, a.frame.origin_mm),
                        u: placed_direction(frame, a.frame.u),
                        v: placed_direction(frame, a.frame.v),
                    },
                }));
        }
    }
}
impl EntryChamferBody {
    pub(in crate::geometry) fn receiving_intervals(
        &self,
        frame: PlaneFrame,
    ) -> Option<ReceivingThreadIntervals> {
        self.receiving
            .as_ref()
            .map(|stock| ReceivingThreadIntervals {
                frame,
                requested_tool_length_mm: stock.length,
                physical_thread_length_mm: self.geometry.material_length_mm(),
                stock_length_mm: stock.length,
            })
    }
}

fn placed_direction(frame: PlaneFrame, p: Point) -> Point {
    let n = frame.normal();
    std::array::from_fn(|i| frame.u[i] * p[0] + frame.v[i] * p[1] + n[i] * p[2])
}
