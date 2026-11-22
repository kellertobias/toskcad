//! Supplied author/thread correspondence, distinct from inferred repair and
//! independently occupied material. No head/cavity/mass admission is inferred.
use super::*;
use crate::phillips_drive::PhillipsScrew;
/// The COMPLETE original entry curve in each retained U-band. World C0/Du
/// use its actual angle knots and actual end-B coefficient, never another V.
#[derive(Clone, Copy, Debug)]
pub struct SuppliedThreadRimTransfer {
    pub face_entity: u32,
    pub surface_entity: u32,
    pub original_rail_entities: [u32; 2],
    pub original_v_start_cell: usize,
    pub original_u_cells: usize,
    pub distance_mm: f64,
    /// Derivative in each original normalized U-cell of this whole rim.
    pub derivative_mm: f64,
}
pub struct SuppliedThreadCorrespondence<'source, 'a, 'input> {
    source: &'source PhillipsScrew,
    body: &'a Body<'input>,
    frame: [Point; 3],
    transfers: Vec<AnalyticThreadTransfer>,
    rims: Vec<SuppliedThreadRimTransfer>,
    new_work: usize,
    inherited_work: usize,
    retained: usize,
}
impl<'source, 'a, 'input> SuppliedThreadCorrespondence<'source, 'a, 'input> {
    pub fn source(&self) -> &'source PhillipsScrew {
        self.source
    }
    pub fn body(&self) -> &'a Body<'input> {
        self.body
    }
    pub fn frame(&self) -> [Point; 3] {
        self.frame
    }
    /// Directed determinant of the actual supplied columns, including its
    /// rounded cross-product normal. No unit/rigid normalization is assumed.
    pub fn frame_determinant_bounds(&self) -> [f64; 2] {
        let value = frame_determinant(self.frame);
        [value.lo, value.hi]
    }
    pub fn transfers(&self) -> &[AnalyticThreadTransfer] {
        &self.transfers
    }
    pub fn entry_rim_transfers(&self) -> &[SuppliedThreadRimTransfer] {
        &self.rims
    }
    pub fn new_work_used(&self) -> usize {
        self.new_work
    }
    pub fn inherited_head_work_used(&self) -> usize {
        self.inherited_work
    }
    /// Decoded head prerequisite plus this fresh correspondence only. Caller
    /// author construction precedes this request and is not claimed here.
    pub fn total_work_used(&self) -> usize {
        self.inherited_work + self.new_work
    }
    /// Conservative complete author +decoded graph +receipt capacity envelope.
    pub fn retained_budget_bytes(&self) -> usize {
        self.retained
    }
}
#[derive(Clone, Copy)]
pub(super) struct SuppliedLaw<'s> {
    pub source: &'s PhillipsScrew,
    pub frame: [Point; 3],
}
pub(super) fn resident_bytes(source: &PhillipsScrew) -> usize {
    // Profile, Thread and HostEnvelopeProof are inline. Cavity owns exactly a
    // Cutter Arc and a side-key slice Arc; both actual reachable payloads and
    // counter/alignment envelopes are included even if another owner shares it.
    std::mem::size_of_val(source)
        + std::mem::size_of_val(source.cavity().source())
        + std::mem::size_of_val(source.cavity().side_face_keys())
        + 128
}
pub(super) fn valid_frame(frame: [Point; 3]) -> bool {
    if frame
        .iter()
        .flatten()
        .any(|x| !x.is_finite() || x.abs() > 1e7)
    {
        return false;
    }
    let determinant = frame_determinant(frame);
    determinant.lo > 0. && determinant.hi.is_finite()
}
fn frame_determinant(frame: [Point; 3]) -> I {
    let target = Target::raw(frame);
    let n = target.columns[3];
    let [_, u, v] = frame;
    let cross = std::array::from_fn::<_, 3, _>(|k| {
        I::point(u[(k + 1) % 3]) * I::point(v[(k + 2) % 3])
            - I::point(u[(k + 2) % 3]) * I::point(v[(k + 1) % 3])
    });
    (0..3).fold(I::point(0.), |s, k| s + I::point(n[k]) * cross[k])
}
impl<'h, 'q, 'r, 'a, 'input> ClosedRepairedHead<'h, 'q, 'r, 'a, 'input> {
    /// Prove the actual supplied cached-rate thread and world frame against
    /// every original raw coefficient cell. The decoded Body is never changed.
    /// Actual author head/cavity correspondence and material closure/mass are
    /// separate trust contracts; this cannot construct or admit native geometry.
    pub fn supplied_thread_correspondence<'source>(
        &self,
        source: &'source PhillipsScrew,
        frame: [Point; 3],
        options: AnalyticThreadOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> std::result::Result<SuppliedThreadCorrespondence<'source, 'a, 'input>, PhillipsReadError>
    {
        let value = self.build_thread_repair(
            options,
            interrupt,
            false,
            None,
            Some(SuppliedLaw { source, frame }),
        )?;
        Ok(SuppliedThreadCorrespondence {
            source,
            body: value.body,
            frame: value.frame,
            transfers: value.transfers,
            rims: value.rims,
            new_work: value.new_work,
            inherited_work: value.inherited_work,
            retained: value.retained,
        })
    }
}
#[cfg(test)]
mod tests;
