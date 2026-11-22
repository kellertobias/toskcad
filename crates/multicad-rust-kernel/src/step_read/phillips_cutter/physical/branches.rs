//! Independent declared-intersection section proofs. An edge tube alone is
//! not a closed finite-loop, occupied material, flux or import certificate.
use super::*;
type BranchResult<T> = std::result::Result<T, PhillipsReadError>;
mod functional;
mod lifts;
mod matrix;
pub use lifts::{FiniteLift, FiniteLifts, LiftKind, LiftOptions};
mod derived;
mod section;
mod sewn;
pub use derived::{
    ClosedDerivedCavityCertificate, DerivedCavityOptions, DerivedCavityRepresentative,
    DerivedEdgeCorrespondence, DerivedFaceCorrespondence, DerivedVertexCorrespondence,
    RawCavitySupportCorrespondence, RawCavitySupportFace, RawCavitySupportOptions,
};
pub use sewn::{SewnEdgeSource, SewnEdgeSources, SewnEvaluation, SewnOptions, SewnParameter};
use supports::chart::Chart;
const WORK: usize = 16_000_000;
const BYTES: usize = 64 * 1024 * 1024;
const CELLS: usize = 65536;
const SCRATCH: usize = 64 * 1024;
#[derive(Clone, Copy)]
pub struct Options {
    pub maximum_position_error_mm: f64,
    pub maximum_cells: usize,
    /// Includes the already earned C0 receipt work; no budget is refunded.
    pub maximum_work: usize,
    pub maximum_bytes: usize,
}
impl Default for Options {
    fn default() -> Self {
        Self {
            maximum_position_error_mm: 1e-7,
            maximum_cells: CELLS,
            maximum_work: WORK,
            maximum_bytes: BYTES,
        }
    }
}
struct Budget<'a> {
    remaining: usize,
    interrupt: Option<&'a dyn Fn() -> bool>,
}
impl Budget<'_> {
    fn charge(&mut self, n: usize) -> BranchResult<()> {
        if self.interrupt.is_some_and(|f| f()) {
            return Err(PhillipsReadError::Cancelled);
        }
        self.remaining = self
            .remaining
            .checked_sub(n)
            .ok_or(PhillipsReadError::Limit)?;
        Ok(())
    }
}
#[derive(Clone, Copy)]
struct Cell {
    poles: [[I; 3]; 4],
    parameter: [f64; 2],
    depth: u8,
}
#[derive(Clone, Copy, Debug)]
pub struct EdgeSection {
    span: usize,
    proof: section::Section,
}
pub struct EdgeTube {
    edge: usize,
    sections: std::ops::Range<usize>,
    functional: functional::Functional,
}
/// The cubic/support association has a unique nearby implicit branch on each
/// recorded section. Endpoints of different edges are NOT coupled here.
/// This is neither a closed finite-face filling nor a material/mass receipt.
pub struct Tubes<'a> {
    source: CubicSupportBounds<'a>,
    edges: Vec<EdgeTube>,
    sections: Vec<EdgeSection>,
    work_used: usize,
}
fn add(a: usize, b: usize) -> usize {
    a.saturating_add(b)
}
fn allocation<T>(v: &Vec<T>) -> usize {
    v.capacity().saturating_mul(std::mem::size_of::<T>())
}
fn arena(b: &Boundary) -> usize {
    let mut size = std::mem::size_of::<Boundary>();
    size = add(size, b.exact_support_keys_capacity_bytes());
    size = add(size, allocation(&b.vertices));
    size = add(size, allocation(&b.edges));
    size = add(size, allocation(&b.faces));
    size = add(size, allocation(&b.mouth.material_uses));
    size = add(size, allocation(&b.mouth.material_loop.uses));
    for face in &b.faces {
        size = add(size, allocation(&face.uses));
        if let Some(loop_) = &face.material_loop {
            size = add(size, allocation(&loop_.uses));
        }
    }
    for edge in &b.edges {
        if let Curve::CubicIntersection {
            pole_entities,
            spans,
            ..
        } = &edge.curve
        {
            size = add(size, allocation(pole_entities));
            size = add(size, allocation(spans));
        }
    }
    size
}
fn memory(
    base: usize,
    capacities: [usize; 3],
    limit: usize,
    budget: &mut Budget<'_>,
) -> BranchResult<()> {
    budget.charge(1)?;
    let bytes = add(
        base,
        add(
            capacities[0].saturating_mul(std::mem::size_of::<EdgeTube>()),
            add(
                capacities[1].saturating_mul(std::mem::size_of::<EdgeSection>()),
                capacities[2].saturating_mul(std::mem::size_of::<Chart>()),
            ),
        ),
    );
    if bytes > limit {
        Err(PhillipsReadError::Limit)
    } else {
        Ok(())
    }
}
fn scratch(bytes: usize, budget: &mut Budget<'_>) -> BranchResult<()> {
    budget.charge(1)?;
    if bytes > SCRATCH / 2 {
        Err(PhillipsReadError::Limit)
    } else {
        Ok(())
    }
}
impl<'a> CubicSupportBounds<'a> {
    pub fn intersection_tubes(
        self,
        options: Options,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> BranchResult<Tubes<'a>> {
        if !options.maximum_position_error_mm.is_finite()
            || options.maximum_position_error_mm <= 0.
            || options.maximum_cells == 0
            || options.maximum_cells > CELLS
            || options.maximum_work == 0
            || options.maximum_work > WORK
            || options.maximum_bytes == 0
            || options.maximum_bytes > BYTES
        {
            return Err(PhillipsReadError::Invalid);
        }
        if interrupt.is_some_and(|f| f()) {
            return Err(PhillipsReadError::Cancelled);
        }
        let mut budget = Budget {
            remaining: options
                .maximum_work
                .checked_sub(self.work_used())
                .ok_or(PhillipsReadError::Limit)?,
            interrupt,
        };
        let boundary = self.boundary();
        let base = add(
            add(
                add(arena(boundary), self.retained_capacity_bytes()),
                std::mem::size_of::<Tubes<'_>>(),
            ),
            SCRATCH,
        );
        let planned = [
            self.edges().len(),
            options.maximum_cells.min(budget.remaining / 10_000),
            boundary.faces().len(),
        ];
        memory(base, planned, options.maximum_bytes, &mut budget)?;
        let mut edges = Vec::new();
        edges
            .try_reserve_exact(planned[0])
            .map_err(|_| PhillipsReadError::Limit)?;
        let mut sections = Vec::new();
        sections
            .try_reserve_exact(planned[1])
            .map_err(|_| PhillipsReadError::Limit)?;
        let mut charts = Vec::new();
        charts
            .try_reserve_exact(planned[2])
            .map_err(|_| PhillipsReadError::Limit)?;
        memory(
            base,
            [edges.capacity(), sections.capacity(), charts.capacity()],
            options.maximum_bytes,
            &mut budget,
        )?;
        for face in boundary.faces() {
            budget.charge(256)?;
            charts.push(Chart::new(&face.support).map_err(PhillipsReadError::from)?);
        }
        // Depth-first dyadic subdivision retains at most25 cells. All bounded
        // polynomial temporaries and this stack fit the fixed64KiB scratch.
        let mut pending = Vec::new();
        pending
            .try_reserve_exact(25)
            .map_err(|_| PhillipsReadError::Limit)?;
        scratch(allocation(&pending), &mut budget)?;
        let mut cells = 0usize;
        for (edge_index, edge) in boundary.edges().iter().enumerate() {
            let Curve::CubicIntersection { spans, .. } = &edge.curve else {
                continue;
            };
            budget.charge(1 + spans.len() * 16)?;
            let functional = functional::Functional::new(spans)?;
            functional.monotone(spans)?;
            let surfaces = edge.owners.map(|i| &boundary.faces()[i].support);
            let owners = edge.owners.map(|i| &charts[i]);
            let start = sections.len();
            for (span_index, span) in spans.iter().enumerate() {
                pending.push(Cell {
                    poles: span.poles_mm.map(|p| p.map(support::adjacent)),
                    parameter: [0., 1.],
                    depth: 0,
                });
                while let Some(cell) = pending.pop() {
                    // Two degree<=12 implicit polynomials, complete derivative
                    // boxes and a 2x2 proof use a fixed bounded operation count.
                    budget.charge(10_000)?;
                    cells += 1;
                    if cells > options.maximum_cells {
                        return Err(PhillipsReadError::Limit);
                    }
                    if let Some(proof) = section::certify(
                        surfaces,
                        owners,
                        functional,
                        cell.poles,
                        cell.parameter,
                        options.maximum_position_error_mm,
                    ) {
                        sections.push(EdgeSection {
                            span: span_index,
                            proof,
                        });
                    } else {
                        if cell.depth == 24 {
                            return Err(PhillipsReadError::Unsupported);
                        }
                        let mid = cell.parameter[0].midpoint(cell.parameter[1]);
                        let [left, right] = section::split(cell.poles);
                        pending.push(Cell {
                            poles: right,
                            parameter: [mid, cell.parameter[1]],
                            depth: cell.depth + 1,
                        });
                        pending.push(Cell {
                            poles: left,
                            parameter: [cell.parameter[0], mid],
                            depth: cell.depth + 1,
                        });
                    }
                }
            }
            edges.push(EdgeTube {
                edge: edge_index,
                sections: start..sections.len(),
                functional,
            });
        }
        budget.charge(1)?;
        Ok(Tubes {
            source: self,
            edges,
            sections,
            work_used: options.maximum_work - budget.remaining,
        })
    }
}
impl EdgeTube {
    pub fn edge_index(&self) -> usize {
        self.edge
    }
    pub fn section_range(&self) -> std::ops::Range<usize> {
        self.sections.clone()
    }
    /// Fixed exact chosen linear functional; this is not a unit direction.
    pub fn section_direction(&self) -> [f64; 3] {
        self.functional.direction
    }
    pub fn pivot_coordinate(&self) -> usize {
        self.functional.pivot
    }
}
impl EdgeSection {
    pub fn span_index(&self) -> usize {
        self.span
    }
    /// Exact dyadic cell within the original Bezier span's local [0,1].
    pub fn parameter_range(&self) -> [f64; 2] {
        self.proof.parameter
    }
    pub fn branch_coordinate_bounds_mm(&self) -> [[f64; 2]; 3] {
        self.proof.source_box
    }
    /// Componentwise error relative to the actual literal cubic at the SAME
    /// fixed section. This is derived from residual and conditioning.
    pub fn position_error_bounds_mm(&self) -> [f64; 3] {
        self.proof.correction_mm
    }
    /// Componentwise derivative error per unit of the ORIGINAL local Bezier
    /// span parameter, including the actual implicit branch derivative.
    pub fn derivative_error_bounds_mm(&self) -> [f64; 3] {
        self.proof.derivative_correction_mm
    }
    pub fn contraction_bound(&self) -> f64 {
        self.proof.contraction
    }
}
impl<'a> Tubes<'a> {
    pub fn boundary(&self) -> &'a Boundary {
        self.source.boundary()
    }
    pub fn edges(&self) -> &[EdgeTube] {
        &self.edges
    }
    pub fn sections(&self) -> &[EdgeSection] {
        &self.sections
    }
    pub fn work_used(&self) -> usize {
        self.work_used
    }
    /// Actual retained receipt capacities, excluding its borrowed parsed arena.
    pub fn retained_capacity_bytes(&self) -> usize {
        add(
            std::mem::size_of::<Self>(),
            add(
                allocation(&self.edges),
                add(
                    allocation(&self.sections),
                    self.source
                        .retained_capacity_bytes()
                        .saturating_sub(std::mem::size_of::<CubicSupportBounds<'_>>()),
                ),
            ),
        )
    }
}
#[cfg(test)]
mod tests;

#[cfg(test)]
pub(super) mod representative;
