//! Actual decoded degree11 nets, with shared rail/end coefficients and complete
//! original-cell C0/C1 transfer bounds. No material, mass or export admission.
use super::*;
use crate::tapered_thread::RetainedTensor11;
mod occupied;
pub use occupied::{Degree11OccupiedOptions, OccupiedDegree11Thread};
mod material;
pub use material::{
    Degree11ClosedMaterial, Degree11JoinCollar, Degree11MaterialFace, Degree11MaterialOptions,
    Degree11MaterialVolume, Degree11MaterialVolumeOptions, Degree11OwnedMaterialOptions,
    Degree11RimUse, OwnedDegree11Material, SuppliedDegree11MaterialMatch,
    SuppliedDegree11MaterialOptions,
};
/// A separately bounded transfer stage after the head has already qualified.
/// `maximum_work` bounds NEW work only; inherited head work is reported
/// separately. This does not change DerivedThreadOptions' whole-request budget.
#[derive(Clone, Copy)]
pub struct Degree11ThreadOptions {
    pub maximum_distance_mm: f64,
    pub maximum_cell_derivative_mm: f64,
    pub maximum_controls: usize,
    pub maximum_work: usize,
    pub maximum_bytes: usize,
}
impl Default for Degree11ThreadOptions {
    fn default() -> Self {
        Self {
            maximum_distance_mm: 1e-7,
            maximum_cell_derivative_mm: 1e-7,
            maximum_controls: MAX_POINTS,
            maximum_work: WORK,
            maximum_bytes: BYTES,
        }
    }
}
/// Immutable defining net; source IDs retain their original topology identity.
pub struct Degree11ThreadFace {
    source: RetainedTensor11,
    original_controls: Vec<u32>,
    forward: bool,
    transfer: ThreadCarrierTransfer,
}
impl Degree11ThreadFace {
    pub fn source(&self) -> &RetainedTensor11 {
        &self.source
    }
    /// Row-major IDs in the original visible control grid, before measured
    /// rail/end corrections. Coincident coordinates do not merge owners.
    pub fn original_control_entities(&self) -> &[u32] {
        &self.original_controls
    }
    pub fn forward(&self) -> bool {
        self.forward
    }
    pub fn transfer(&self) -> ThreadCarrierTransfer {
        self.transfer
    }
}
/// This packet retains the actual file Body. It cannot enter the distinct D9
/// occupied providers and does not certify finite cap, winding or head contact.
pub struct CoherentDegree11ThreadCarrier<'head, 'a, 'input> {
    body: &'a Body<'input>,
    origin: &'head ClosedRepairedHead<'head, 'head, 'head, 'a, 'input>,
    chart: CommonMouthChart,
    faces: Vec<Degree11ThreadFace>,
    stations: [f64; 2],
    new_work: usize,
    inherited_head_work: usize,
    retained: usize,
}
impl<'head, 'a, 'input> CoherentDegree11ThreadCarrier<'head, 'a, 'input> {
    pub fn body(&self) -> &'a Body<'input> {
        self.body
    }
    /// Immutable provenance, not a profile reconstructed from thread headers.
    /// closed_repaired_meridian computes fixed line/arc parts from this exact
    /// Profile; its options only refuse, never choose different coefficients.
    /// Holds the complete originating receipt and its immutable common/cavity
    /// prerequisite graph alive. Pointer identity cannot be recycled while
    /// this borrow exists; equivalent separately rebuilt receipts are distinct.
    pub fn originating_head(&self) -> &'head ClosedRepairedHead<'head, 'head, 'head, 'a, 'input> {
        self.origin
    }
    pub fn original_head_profile(&self) -> &'a crate::head_profile::Profile {
        self.origin.source().profile_proposal()
    }
    pub fn chart(&self) -> CommonMouthChart {
        self.chart
    }
    pub fn faces(&self) -> &[Degree11ThreadFace] {
        &self.faces
    }
    pub fn cap_stations_mm(&self) -> [f64; 2] {
        self.stations
    }
    pub fn new_work_used(&self) -> usize {
        self.new_work
    }
    pub fn inherited_head_work_used(&self) -> usize {
        self.inherited_head_work
    }
    /// Both already completed stages; not a claim of a sixteen-million import.
    pub fn total_work_used(&self) -> usize {
        self.new_work + self.inherited_head_work
    }
    /// Inherited head plus retained candidate. Caller sharing the Body must
    /// count it once across sibling packets; temporary scratch is not retained.
    pub fn retained_budget_bytes(&self) -> usize {
        self.retained
    }
}
struct Candidate {
    poles: Vec<Point>,
    original_controls: Vec<u32>,
    dimensions: [usize; 2],
}
fn valid(o: Degree11ThreadOptions) -> bool {
    valid_options(DerivedThreadOptions {
        maximum_distance_mm: o.maximum_distance_mm,
        maximum_cell_derivative_mm: o.maximum_cell_derivative_mm,
        maximum_controls: o.maximum_controls,
        maximum_work: o.maximum_work,
        maximum_bytes: o.maximum_bytes,
    })
}
impl<'h, 'q, 'r, 'a, 'input> ClosedRepairedHead<'h, 'q, 'r, 'a, 'input> {
    pub fn degree11_thread_carrier<'head>(
        &'head self,
        options: Degree11ThreadOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> std::result::Result<CoherentDegree11ThreadCarrier<'head, 'a, 'input>, PhillipsReadError>
    {
        let mut work = Work {
            remaining: options.maximum_work,
            interrupt,
            cancelled: Cell::new(false),
        };
        let result = (|| {
            work.step(1)?;
            if !valid(options) {
                return Err(ReadError::Invalid);
            }
            let base = self
                .retained_budget_bytes()
                .saturating_add(std::mem::size_of::<CoherentDegree11ThreadCarrier<'_, '_, '_>>());
            if base.saturating_add(SCRATCH) > options.maximum_bytes {
                return Err(ReadError::Limit);
            }
            let mapping = chart::Chart::new(self);
            let mut shell = Shell::new(&self.body().entities, self.body().tolerance, true);
            let raw = raw(
                self,
                &mut shell,
                &mut work,
                options.maximum_bytes - base - SCRATCH,
            )?;
            let count: usize = raw
                .iter()
                .map(|r| r.ids.len() * (r.ids[0].len() - r.entry))
                .sum();
            if count > options.maximum_controls {
                return Err(ReadError::Limit);
            }
            let raw_count: usize = raw.iter().map(|r| r.ids.len() * r.ids[0].len()).sum();
            // Bounds simultaneously live raw IDs/headers, Shell's point map,
            // candidate point+ID vectors, difference coefficients, Arc creation
            // copies and output headers. No unbounded per-cell subdivision.
            let peak = base
                .saturating_add(SCRATCH)
                .saturating_add(raw_count.saturating_mul(128))
                .saturating_add(count.saturating_mul(160));
            if peak > options.maximum_bytes {
                return Err(ReadError::Limit);
            }
            let stations = [
                self.virtual_bottom_station()[1],
                midpoint(mapping.local(shell.point(*raw[0].ids[0].last().unwrap())?)?[2]),
            ];
            require(stations[1] > stations[0])?;
            let mut candidates = Vec::with_capacity(raw.len());
            for r in &raw {
                let dimensions = [r.ids.len(), r.ids[0].len() - r.entry];
                let n = dimensions[0] * dimensions[1];
                let mut poles = Vec::with_capacity(n);
                let mut original_controls = Vec::with_capacity(n);
                for row in &r.ids {
                    for &id in &row[r.entry..] {
                        work.step(100)?;
                        poles.push(mapping.local(shell.point(id)?)?.map(midpoint));
                        original_controls.push(id);
                    }
                }
                for row in poles.chunks_exact_mut(dimensions[1]) {
                    row[0][2] = stations[0];
                    row[dimensions[1] - 1][2] = stations[1];
                }
                candidates.push(Candidate {
                    poles,
                    original_controls,
                    dimensions,
                });
            }
            for i in 0..candidates.len() {
                let next = (i + 1) % candidates.len();
                // verify_prefix checks the entire original rail-ID ring. Its
                // tensor_surface/Proof::poles checks every row's exact count
                // against ONE axial basis, not just geometric proximity or
                // compatible knots. Outer->Common preserves ordered receipts;
                // raw additionally verifies each retained entry column equals
                // that common basis.visible_start. These are inherited source
                // invariants; neither mismatch is a new user-facing refusal.
                assert_eq!(
                    raw[i].rails[1], raw[next].rails[0],
                    "qualified rail owner ring"
                );
                assert_eq!(
                    candidates[i].dimensions[1], candidates[next].dimensions[1],
                    "qualified common axial grid"
                );
                let height = candidates[i].dimensions[1];
                work.step(height)?;
                for j in 0..height {
                    let point = candidates[next].poles[j];
                    let last = candidates[i].poles.len() - height + j;
                    candidates[i].poles[last] = point;
                }
            }
            let mut faces = Vec::with_capacity(raw.len());
            for (i, (r, c)) in raw.iter().zip(candidates).enumerate() {
                let [width, height] = c.dimensions;
                let mut difference = Vec::with_capacity(c.poles.len());
                let mut distance = 0_f64;
                for (&id, &p) in c.original_controls.iter().zip(&c.poles) {
                    work.step(200)?;
                    let original = shell.point(id)?.map(|x| I::new(x.next_down(), x.next_up()));
                    let actual = mapping.world(p.map(I::point));
                    let d = std::array::from_fn(|a| original[a] - actual[a]);
                    distance = distance.max(reduce::l1(d));
                    difference.push(d);
                }
                let mut derivative = [0_f64; 2];
                for u in 0..width {
                    for v in 0..height {
                        work.step(30)?;
                        let d = difference[u * height + v];
                        for (axis, next) in [
                            (0, (u + 1 < width).then_some((u + 1) * height + v)),
                            (1, (v + 1 < height).then_some(u * height + v + 1)),
                        ] {
                            if let Some(next) = next {
                                derivative[axis] =
                                    derivative[axis].max(reduce::l1(std::array::from_fn(|a| {
                                        I::point(11.) * (difference[next][a] - d[a])
                                    })));
                            }
                        }
                    }
                }
                // Bernstein convex hull and derivative-control differences
                // cover each COMPLETE ORIGINAL normalized degree11 cell.
                if distance > options.maximum_distance_mm
                    || derivative
                        .into_iter()
                        .any(|d| d > options.maximum_cell_derivative_mm)
                {
                    return Err(ReadError::Unsupported);
                }
                let spans = [(width - 1) / 11, (height - 1) / 11];
                // The sealed raw/common-head receipts retain complete degree11
                // rows and an entry knot before the final visible cell. Filling
                // above preserves their exact (11*spans+1) dimensions. Aggregate
                // count<=262144 also implies the unchanged per-net524288 cap.
                // Chart::local checked finite abs<1e7 before midpoint selection;
                // halves cannot overflow and their sum lies in the same finite
                // interval. Rail copies preserve those values. The first endZ
                // is the actual retained meridian's final station, not recipe
                // metadata: Profile::from_transported_segments already checked
                // all physical stations finite with abs<=1e7. The other endZ
                // comes through that same checked local chart. Revalidation
                // therefore cannot return a user/source error here.
                let source = RetainedTensor11::new(spans, c.poles)
                    .expect("validated degree-eleven carrier grid");
                faces.push(Degree11ThreadFace {
                    source,
                    original_controls: c.original_controls,
                    forward: r.forward,
                    transfer: ThreadCarrierTransfer {
                        face_entity: r.face,
                        surface_entity: r.surface,
                        source_face_index: 2 + i,
                        original_v_start_cell: r.entry / 11,
                        original_rail_entities: r.rails,
                        distance_mm: distance,
                        derivative_mm: derivative,
                        original_cells: spans,
                        source_cells: spans,
                    },
                });
            }
            let retained = base
                + faces.capacity() * std::mem::size_of::<Degree11ThreadFace>()
                + faces
                    .iter()
                    .map(|f| {
                        std::mem::size_of_val(f.source.poles())
                            + f.original_controls.capacity() * std::mem::size_of::<u32>()
                            + 2 * std::mem::size_of::<usize>()
                    })
                    .sum::<usize>();
            Ok(CoherentDegree11ThreadCarrier {
                body: self.body(),
                origin: self,
                chart: self.chart(),
                faces,
                stations,
                new_work: options.maximum_work - work.remaining,
                inherited_head_work: self.work_used(),
                retained,
            })
        })();
        result.map_err(|e| bound_error(e, &work))
    }
}

#[cfg(test)]
pub(super) mod tests;
