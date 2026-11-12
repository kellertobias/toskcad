//! A full-height stock cut reuses its qualified chart and actual entry material.
use super::super::helical_body::ReceivingThreadIntervals;
use super::*;

pub(super) type ReceivingStock = super::super::helical_body::full_height::FullHeightStock;
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
        let (stock, placed) = ReceivingStock::qualify(
            self,
            spec,
            root_radius_mm,
            feature,
            frame,
            policy,
            interrupt,
        )?;
        Self::make_entry_chamfer(
            spec,
            root_radius_mm,
            radial_depth_mm,
            terminal,
            feature,
            placed,
            policy,
            interrupt,
            Some(stock),
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
        let far = if body.terminal_relief().is_some() {
            source::Face::FarCap
        } else {
            source::Face::EndCap
        };
        let ordinal = (0..body.face_count())
            .find(|&i| body.face(i) == Some(far))
            .expect("retained source cap");
        self.extend_caps(frame, [0, ordinal], elements);
    }
}
impl EntryChamferBody {
    pub(in crate::geometry) fn receiving_intervals(
        &self,
        frame: PlaneFrame,
    ) -> Option<ReceivingThreadIntervals> {
        self.receiving
            .as_ref()
            .map(|stock| stock.intervals(frame, self.geometry.material_length_mm()))
    }
}
