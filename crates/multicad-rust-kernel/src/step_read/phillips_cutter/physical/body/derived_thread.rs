//! A distinct repaired polynomial source, measured against actual decoded nets.
//! Header laws are not source authority. Raw finite-curve matching is inherited;
//! this receipt alone is neither the complete head material nor decoded mass.
use super::*;
use crate::{
    radial_tensor as radial,
    tapered_thread::{RetainedTensor, interval::I},
};
mod analytic;
mod chart;
mod degree11;
mod graph;
mod join;
mod reduce;
pub use analytic::{
    AnalyticThreadOptions, AnalyticThreadRepair, AnalyticThreadTransfer, AuthoredPeriodicOptions,
    CommonChartPeriodicThreadRepair, OwnedPeriodicMaterial, PeriodicAuthoredApproximation,
    PeriodicClosedMaterial, PeriodicCurveOwner, PeriodicFaceProperties, PeriodicFiniteCurve,
    PeriodicFiniteCurves, PeriodicFiniteKind, PeriodicFiniteOptions, PeriodicLiteralCertificate,
    PeriodicMaterialFace, PeriodicMaterialOptions, PeriodicMaterialProvenance,
    PeriodicMaterialVolume, PeriodicOwnedMaterialQuery, PeriodicOwnedOptions,
    PeriodicShoulderOptions, PeriodicShoulderRepair, PeriodicThreadBounds, PeriodicThreadRepair,
    PeriodicThreadSource, PeriodicVolumeOptions,
};
pub use degree11::{CoherentDegree11ThreadCarrier, Degree11ThreadFace, Degree11ThreadOptions};
pub use join::{
    ContactThreadCarrier, Degree11ContactOptions, Degree11DomainEnd, Degree11RimContact,
    Degree11RimContactId, Degree11RimKnot, Degree11SignDomain, Degree11ThreadContacts,
    DyadicParameter, RimContact, RimContactId, RimDomainEnd, RimDomainRole, RimKnot, RimSignDomain,
    ThreadContactOptions,
};
const WORK: usize = 16_000_000;
const BYTES: usize = 64 * 1024 * 1024;
const SCRATCH: usize = 1024 * 1024;
const MAX_POINTS: usize = 262_144;
#[derive(Clone, Copy)]
pub struct DerivedThreadOptions {
    pub maximum_distance_mm: f64,
    /// World derivative discrepancy per ORIGINAL normalized Bernstein cell.
    pub maximum_cell_derivative_mm: f64,
    /// Candidate control quota; may only lower the hard limit of 262,144.
    pub maximum_controls: usize,
    pub maximum_work: usize,
    pub maximum_bytes: usize,
}
impl Default for DerivedThreadOptions {
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
#[derive(Clone, Copy, Debug, Default)]
pub struct ThreadCarrierTransfer {
    pub face_entity: u32,
    pub surface_entity: u32,
    pub source_face_index: usize,
    /// First original axial Bernstein cell retained after the actual trim.
    pub original_v_start_cell: usize,
    pub original_rail_entities: [u32; 2],
    /// Whole-cell world C0 L1 bound, also bounding Euclidean displacement.
    pub distance_mm: f64,
    /// World derivative L1 bounds for each ORIGINAL normalized Bernstein cell,
    /// ordered u,v; these are not derivatives in the global STEP knot units.
    pub derivative_mm: [f64; 2],
    pub original_cells: [usize; 2],
    pub source_cells: [usize; 2],
}
/// A distinct local-chart polynomial boundary. Both end disks are explicitly
/// derived closures, not a claim that the original file was exactly sewn.
/// Transfer records cover the visible tensor cells; finite raw curve matching
/// remains the separately identified prerequisite.
pub struct DerivedThreadCarrier<'a, 'input> {
    body: &'a Body<'input>,
    chart: CommonMouthChart,
    input: radial::BoundaryInput,
    bounds: Vec<ThreadCarrierTransfer>,
    stations: [f64; 2],
    entry_radius_upper: f64,
    local_enclosure: [[f64; 2]; 3],
    work: usize,
    retained: usize,
}
impl<'a, 'input> DerivedThreadCarrier<'a, 'input> {
    pub fn body(&self) -> &'a Body<'input> {
        self.body
    }
    pub fn chart(&self) -> CommonMouthChart {
        self.chart
    }
    pub fn source(&self) -> &radial::BoundaryInput {
        &self.input
    }
    pub fn transfers(&self) -> &[ThreadCarrierTransfer] {
        &self.bounds
    }
    pub fn cap_stations_mm(&self) -> [f64; 2] {
        self.stations
    }
    /// Conservative radius in the common local chart, not an affine-world wall.
    pub fn entry_radius_upper_mm(&self) -> f64 {
        self.entry_radius_upper
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_budget_bytes(&self) -> usize {
        self.retained
    }
    /// Separate bounded operation: the existing ownership/trim/axial/cap/order/
    /// occupied-source providers retain their own hard work/size limits. This
    /// consumes the candidate without duplicating its entire input. No final
    /// head/cavity join, raw-file closure or manufacturing publication follows.
    pub fn certify_occupied(
        self,
        tolerance: f64,
        minimum: f64,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> std::result::Result<OccupiedDerivedThread<'a, 'input>, crate::tapered_thread::Error> {
        let boundary = radial::BoundaryOwnership::new(self.input, tolerance, interrupt)?
            .bind_surfaces(tolerance, interrupt)?
            .certify_axial(tolerance, minimum, interrupt)?
            .certify_caps(interrupt)?
            .certify_ordering(interrupt)?
            .certify_material(tolerance, interrupt)?;
        Ok(OccupiedDerivedThread {
            body: self.body,
            chart: self.chart,
            boundary,
            bounds: self.bounds,
            stations: self.stations,
            entry_radius_upper: self.entry_radius_upper,
            local_enclosure: self.local_enclosure,
            construction_work: self.work,
        })
    }
}
pub struct OccupiedDerivedThread<'a, 'input> {
    body: &'a Body<'input>,
    chart: CommonMouthChart,
    boundary: radial::RadialBoundary,
    bounds: Vec<ThreadCarrierTransfer>,
    stations: [f64; 2],
    entry_radius_upper: f64,
    local_enclosure: [[f64; 2]; 3],
    construction_work: usize,
}
impl<'a, 'input> OccupiedDerivedThread<'a, 'input> {
    pub fn body(&self) -> &'a Body<'input> {
        self.body
    }
    pub fn chart(&self) -> CommonMouthChart {
        self.chart
    }
    pub fn boundary(&self) -> &radial::RadialBoundary {
        &self.boundary
    }
    pub fn transfers(&self) -> &[ThreadCarrierTransfer] {
        &self.bounds
    }
    pub fn cap_stations_mm(&self) -> [f64; 2] {
        self.stations
    }
    pub fn entry_radius_upper_mm(&self) -> f64 {
        self.entry_radius_upper
    }
    pub(crate) fn local_enclosure_bounds(&self) -> [[f64; 2]; 3] {
        self.local_enclosure
    }
    /// Carrier construction/transfer only. The separate occupied provider keeps
    /// its independently declared hard work budget; this is not its work count.
    pub fn construction_work_used(&self) -> usize {
        self.construction_work
    }
    /// Original Body once plus the actual occupied source and transfer payload.
    /// Shared-Body aggregate callers use additional_budget_bytes instead.
    pub fn retained_budget_bytes(&self) -> usize {
        self.body
            .retained_budget_bytes()
            .saturating_add(self.additional_budget_bytes())
    }
    pub(crate) fn additional_budget_bytes(&self) -> usize {
        (std::mem::size_of::<Self>() - std::mem::size_of::<radial::RadialBoundary>())
            .saturating_add(self.boundary.retained_budget_bytes())
            .saturating_add(
                self.bounds
                    .capacity()
                    .saturating_mul(std::mem::size_of::<ThreadCarrierTransfer>()),
            )
    }
}
struct Grid {
    face: u32,
    surface: u32,
    forward: bool,
    ids: Vec<Vec<u32>>,
}
struct Raw {
    face: u32,
    surface: u32,
    forward: bool,
    ids: Vec<Vec<u32>>,
    rails: [u32; 2],
    entry: usize,
}
#[derive(Clone, Copy)]
struct Proposal {
    limits: DerivedThreadOptions,
    u_parts: usize,
}
struct Net {
    forward: bool,
    spans: [usize; 2],
    poles: Vec<Point>,
}
fn grids(
    body: &Body<'_>,
    shell: &mut Shell<'_, '_>,
    work: &mut Work<'_>,
    maximum_bytes: usize,
) -> Result<Vec<Grid>> {
    let headers = body
        .tensor_face_count
        .saturating_mul(2 * std::mem::size_of::<Grid>());
    if headers > maximum_bytes {
        return Err(ReadError::Limit);
    }
    let mut result = Vec::with_capacity(body.tensor_face_count);
    let mut held = result.capacity() * std::mem::size_of::<Grid>();
    for &face in body.tensor_face_entities() {
        work.step(128)?;
        let (surface, forward, _) = helical::face(shell, face)?;
        let args = shell.args(surface, "B_SPLINE_SURFACE_WITH_KNOTS")?;
        let rs = list(args.get(3).ok_or(ReadError::Invalid)?)?;
        // The original parser caps the entire input at64MiB, hence even
        // one-byte numeric/reference tokens cannot overflow usize32 here.
        let mut count = 0usize;
        for row in rs {
            count += list(row)?.len();
        }
        let next = held
            .saturating_add(count.saturating_mul(2 * std::mem::size_of::<u32>()))
            .saturating_add(rs.len().saturating_mul(2 * std::mem::size_of::<Vec<u32>>()));
        if next > maximum_bytes {
            return Err(ReadError::Limit);
        }
        work.step(count)?;
        let mut ids = Vec::with_capacity(rs.len());
        for row in rs {
            let values = list(row)?;
            let mut points = Vec::with_capacity(values.len());
            for value in values {
                points.push(reference(value)?);
            }
            ids.push(points);
        }
        held = held
            .saturating_add(ids.capacity() * std::mem::size_of::<Vec<u32>>())
            .saturating_add(
                ids.iter()
                    .map(|r| r.capacity() * std::mem::size_of::<u32>())
                    .sum::<usize>(),
            );
        // Explicit reservations never grow. The preceding2x request bound
        // covers these actual capacities; held remains the charge for later rows.
        require(ids.len() > 11 && (ids.len() - 1).is_multiple_of(11))?;
        let height = ids[0].len();
        require(
            height > 11 && (height - 1).is_multiple_of(11) && ids.iter().all(|r| r.len() == height),
        )?;
        result.push(Grid {
            face,
            surface,
            forward,
            ids,
        });
    }
    Ok(result)
}
fn entry_index(
    shell: &Shell<'_, '_>,
    row: &[u32],
    endpoint: Point,
    work: &mut Work<'_>,
) -> Result<usize> {
    work.step(row.len() * 64)?;
    let mut entry = None;
    // With endpoint multiplicity12/interior11, only these controls are
    // actual curve knot points. Interior Bezier poles are not trim stations.
    for index in (0..row.len()).step_by(11) {
        if shell.point(row[index])? == endpoint && entry.replace(index).is_some() {
            return Err(ReadError::Unsupported);
        }
    }
    entry.ok_or(ReadError::Unsupported)
}
fn raw(
    head: &ClosedRepairedHead<'_, '_, '_, '_, '_>,
    shell: &mut Shell<'_, '_>,
    work: &mut Work<'_>,
    maximum_bytes: usize,
) -> Result<Vec<Raw>> {
    let body = head.body();
    let grids = grids(body, shell, work, maximum_bytes)?;
    let mut result = Vec::with_capacity(grids.len());
    for (grid, t) in grids
        .into_iter()
        .zip(head.common_mouth().tensor_transfers())
    {
        // Same immutable Body, exact source index and the successful prefix
        // receipt provide four coedges, same hand, common axial basis and rail
        // ownership. Do not reinterpret/reject alternate POINT entity aliases.
        let point = shell.args(t.entry_vertex_entity, "VERTEX_POINT")?;
        let endpoint = shell.point(reference(&point[1])?)?;
        let entry = entry_index(shell, &grid.ids[0], endpoint, work)?;
        require(entry == t.entry_control_column)?;
        result.push(Raw {
            face: grid.face,
            surface: grid.surface,
            forward: grid.forward,
            ids: grid.ids,
            rails: t.long_rails,
            entry,
        });
    }
    Ok(result)
}
// Only called with Chart::local's finite abs<1e7 enclosure; halving
// before adding cannot overflow. Any chosen-point error is measured below.
fn midpoint(x: I) -> f64 {
    x.lo / 2. + x.hi / 2.
}
fn bound_error(e: ReadError, w: &Work<'_>) -> PhillipsReadError {
    if w.cancelled.get() {
        PhillipsReadError::Cancelled
    } else {
        e.into()
    }
}
fn valid_options(options: DerivedThreadOptions) -> bool {
    options.maximum_distance_mm.is_finite()
        && options.maximum_distance_mm > 0.
        && options.maximum_distance_mm <= 1e-7
        && options.maximum_cell_derivative_mm.is_finite()
        && options.maximum_cell_derivative_mm > 0.
        && options.maximum_cell_derivative_mm <= 1e-7
        && options.maximum_controls <= MAX_POINTS
        && options.maximum_work <= WORK
        && options.maximum_bytes <= BYTES
}
impl<'h, 'q, 'r, 'a, 'input> ClosedRepairedHead<'h, 'q, 'r, 'a, 'input> {
    pub fn derived_thread_carrier(
        &self,
        options: DerivedThreadOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> std::result::Result<DerivedThreadCarrier<'a, 'input>, PhillipsReadError> {
        self.with_candidate(options, interrupt, 1, |_, head_radius, entry_radius, _| {
            if (I::point(head_radius) - entry_radius).lo <= options.maximum_distance_mm {
                return Err(ReadError::Unsupported);
            }
            Ok(())
        })
        .map(|(source, ())| source)
    }
    fn with_candidate<C>(
        &self,
        options: DerivedThreadOptions,
        interrupt: Option<&dyn Fn() -> bool>,
        // Sealed callers use1 (original public path) or4 (nominal prototype).
        u_parts: usize,
        certify: impl FnOnce(&[Net], f64, I, &mut Work<'_>) -> Result<C>,
    ) -> std::result::Result<(DerivedThreadCarrier<'a, 'input>, C), PhillipsReadError> {
        let mut work = Work {
            remaining: options.maximum_work,
            interrupt,
            cancelled: Cell::new(false),
        };
        let result = (|| {
            work.step(1)?;
            if !valid_options(options) {
                return Err(ReadError::Invalid);
            }
            work.step(self.work_used())?;
            let base = self
                .retained_budget_bytes()
                .saturating_add(std::mem::size_of::<DerivedThreadCarrier<'_, '_>>());
            if base.saturating_add(SCRATCH) > options.maximum_bytes {
                return Err(ReadError::Limit);
            }
            let body = self.body();
            let mapping = chart::Chart::new(self);
            let elevation = reduce::Elevation::new();
            let mut shell = Shell::new(&body.entities, body.tolerance, true);
            let raw = raw(
                self,
                &mut shell,
                &mut work,
                options
                    .maximum_bytes
                    .saturating_sub(base)
                    .saturating_sub(SCRATCH),
            )?;
            let raw_points = raw
                .iter()
                .map(|r| r.ids.len() * r.ids[0].len())
                .sum::<usize>();
            // At most four exact subdivisions in each axis: candidate size
            // <=16*rawControls. The64MiB ID frontier bounds raw controls
            // <=16,777,216, so this fits usize even on32bit targets.
            let candidate_points = raw
                .iter()
                .map(|r| {
                    (u_parts * 9 * ((r.ids.len() - 1) / 11) + 1)
                        * (36 * ((r.ids[0].len() - 1 - r.entry) / 11) + 1)
                })
                .sum::<usize>();
            if candidate_points > options.maximum_controls {
                return Err(ReadError::Limit);
            }
            // Per point reserve covers simultaneous raw IDs, interval/local
            // grids, candidate nets, immutable tensor/curve storage and Vec
            // headers; bounded fixed scratch covers <=12x12 reduction arrays.
            let peak = base
                .saturating_add(raw_points.saturating_mul(128))
                .saturating_add(candidate_points.saturating_mul(128))
                .saturating_add(SCRATCH);
            if peak > options.maximum_bytes {
                return Err(ReadError::Limit);
            }
            let entry = self.virtual_bottom_station()[1];
            let first = &raw[0];
            let tip = midpoint(mapping.local(shell.point(*first.ids[0].last().unwrap())?)?[2]);
            require(tip > entry)?;
            let stations = [entry, tip];
            let mut nets = Vec::with_capacity(raw.len());
            let mut bounds = Vec::with_capacity(raw.len());
            for r in &raw {
                let (net, mut bound) = construct(
                    r,
                    &shell,
                    &mapping,
                    &elevation,
                    stations,
                    Proposal {
                        limits: options,
                        u_parts,
                    },
                    &mut work,
                )?;
                bound.source_face_index = 2 + nets.len();
                nets.push(net);
                bounds.push(bound);
            }
            let mut entry_radius = I::point(0.);
            for net in &nets {
                let height = 9 * net.spans[1] + 1;
                for row in net.poles.chunks_exact(height) {
                    work.step(1)?;
                    let p = row[0];
                    let radius = (I::point(p[0]).square() + I::point(p[1]).square()).sqrt();
                    require(radius.valid() && radius.lo > 0.)?;
                    entry_radius.hi = entry_radius.hi.max(radius.hi);
                }
            }
            let contact = certify(
                &nets,
                self.virtual_bottom_station()[0],
                entry_radius,
                &mut work,
            )?;
            // Private Net is the defining inventory: graph::build emits ONLY
            // these Tensor cells, their identical Spline rails/rims, and two
            // finite planar cap fans. Include both actual axis cap centers.
            // No generic analytic Cone/Circle is filtered or approximated.
            let mut local_enclosure = [
                [0., 0.],
                [0., 0.],
                [stations[0].min(stations[1]), stations[0].max(stations[1])],
            ];
            for net in &nets {
                work.step(net.poles.len())?;
                for p in &net.poles {
                    for axis in 0..3 {
                        local_enclosure[axis][0] = local_enclosure[axis][0].min(p[axis]);
                        local_enclosure[axis][1] = local_enclosure[axis][1].max(p[axis]);
                    }
                }
            }
            let input = graph::build(nets, stations, &mut work)?;
            let payload = graph::bytes(&input)
                .saturating_add(bounds.capacity() * std::mem::size_of::<ThreadCarrierTransfer>());
            let retained = base.saturating_add(payload);
            // The earlier128B/raw reserve covers96B of interval grids plus
            // <8B IDs/row headers and <1B descriptors. Candidate Vec+temporary
            // Arc<=48B/control; rails+two rims<3.7B/control and all descriptors
            // <4B/control (each face>=10x37). Actual capacities are charged
            // above; this stronger live preflight also bounds final retention.
            work.step(1)?;
            Ok((
                DerivedThreadCarrier {
                    body,
                    chart: self.chart(),
                    input,
                    bounds,
                    stations,
                    entry_radius_upper: entry_radius.hi,
                    local_enclosure,
                    work: options.maximum_work - work.remaining,
                    retained,
                },
                contact,
            ))
        })();
        result.map_err(|e| bound_error(e, &work))
    }
}
fn restrict_v(p: [[[I; 3]; 12]; 12]) -> [[[[I; 3]; 12]; 12]; 4] {
    let first = p.map(reduce::split);
    let halves = [first.map(|p| p.0), first.map(|p| p.1)];
    let a = halves[0].map(reduce::split);
    let b = halves[1].map(reduce::split);
    [
        a.map(|p| p.0),
        a.map(|p| p.1),
        b.map(|p| p.0),
        b.map(|p| p.1),
    ]
}
// Unqualified proposal coefficients must fit the retained geometry domain
// before transport or publication. This is not an occupied-source assertion.
fn proposal_domain(q: &[[Point; 10]; 10]) -> Result<()> {
    require(
        q.iter()
            .flatten()
            .flatten()
            .all(|x| x.is_finite() && x.abs() < 1e7),
    )
}
fn construct(
    r: &Raw,
    shell: &Shell<'_, '_>,
    mapping: &chart::Chart,
    elevation: &reduce::Elevation,
    stations: [f64; 2],
    proposal: Proposal,
    work: &mut Work<'_>,
) -> Result<(Net, ThreadCarrierTransfer)> {
    let Proposal {
        limits: options,
        u_parts,
    } = proposal;
    let [width, height] = [r.ids.len(), r.ids[0].len()];
    let spans = [
        u_parts * ((width - 1) / 11),
        4 * ((height - 1 - r.entry) / 11),
    ];
    let h = 9 * spans[1] + 1;
    let w = 9 * spans[0] + 1;
    let mut output = vec![[f64::NAN; 3]; w * h];
    let mut original = Vec::with_capacity(width);
    let mut local = Vec::with_capacity(width);
    for row in &r.ids {
        let mut a = Vec::with_capacity(height);
        let mut b = Vec::with_capacity(height);
        for &id in row {
            // Fixed 3x3 interval inverse product, actual point leaf lookup and
            // validity guards; <=100 scalar operations per original control.
            work.step(100)?;
            let p = shell.point(id)?;
            a.push(p.map(|x| I::new(x.next_down(), x.next_up())));
            b.push(mapping.local(p)?);
        }
        original.push(a);
        local.push(b);
    }
    let mut bound = ThreadCarrierTransfer {
        face_entity: r.face,
        surface_entity: r.surface,
        original_v_start_cell: r.entry / 11,
        original_rail_entities: r.rails,
        original_cells: [spans[0] / u_parts, spans[1] / 4],
        source_cells: spans,
        ..Default::default()
    };
    for u in (0..width - 1).step_by(11) {
        for v in (r.entry..height - 1).step_by(11) {
            // Three binary restrictions for each of two 12x12 interval grids.
            // Each has 12*66*3 weighted coordinate averages;20k includes copies.
            work.step(20_000)?;
            let world = std::array::from_fn(|i| std::array::from_fn(|j| original[u + i][v + j]));
            let local = std::array::from_fn(|i| std::array::from_fn(|j| local[u + i][v + j]));
            // Original source unchanged: split the actual interval polynomial
            // by exact1/2 de Casteljau, then measure each new proposal back
            // in the ORIGINAL cell parameter, not its quarter's parameter.
            let us = if u_parts == 1 {
                vec![world]
            } else {
                restrict_v(std::array::from_fn(|j| {
                    std::array::from_fn(|i| world[i][j])
                }))
                .map(|p| std::array::from_fn(|i| std::array::from_fn(|j| p[j][i])))
                .to_vec()
            };
            let ls = if u_parts == 1 {
                vec![local]
            } else {
                restrict_v(std::array::from_fn(|j| {
                    std::array::from_fn(|i| local[i][j])
                }))
                .map(|p| std::array::from_fn(|i| std::array::from_fn(|j| p[j][i])))
                .to_vec()
            };
            for (uchild, (world, local)) in us.into_iter().zip(ls).enumerate() {
                // The outer20k charge covers original V splitting, or the new
                // U split. Each U-quarter's two V restrictions add20k.
                if u_parts != 1 {
                    work.step(20_000)?;
                }
                let worlds = restrict_v(world);
                let locals = restrict_v(local);
                for (child, (actual, local)) in worlds.into_iter().zip(locals).enumerate() {
                    // Work units bound fixed three-coordinate coefficient operations,
                    // not CPU FLOPs. With precomputed weights: <3828 proposal,
                    //1584 elevation,700 transport,144 differences,1344 hull and
                    //400 storage/check operations;10k covers this complete cell.
                    work.step(10_000)?;
                    let nominal = local.map(|row| row.map(|p| p.map(|x| x.lo / 2. + x.hi / 2.)));
                    let columns = std::array::from_fn::<_, 12, _>(|j| {
                        reduce::reduce(std::array::from_fn(|i| nominal[i][j]))
                    });
                    let mut q = std::array::from_fn::<_, 10, _>(|i| {
                        reduce::reduce(std::array::from_fn(|j| columns[j][i]))
                    });
                    for row in &mut q {
                        if v == r.entry && child == 0 {
                            row[0][2] = stations[0];
                        }
                        if v + 11 == height - 1 && child == 3 {
                            row[9][2] = stations[1];
                        }
                    }
                    proposal_domain(&q)?;
                    let raised_v =
                        q.map(|row| elevation.apply(row.map(|p| mapping.world(p.map(I::point)))));
                    let raised = std::array::from_fn::<_, 12, _>(|j| {
                        elevation.apply(std::array::from_fn(|i| raised_v[i][j]))
                    });
                    let difference = std::array::from_fn::<_, 12, _>(|i| {
                        std::array::from_fn::<_, 12, _>(|j| {
                            std::array::from_fn::<_, 3, _>(|a| actual[i][j][a] - raised[j][i][a])
                        })
                    });
                    for i in 0..12 {
                        for j in 0..12 {
                            let e = reduce::l1(difference[i][j]);
                            require(e.is_finite())?;
                            bound.distance_mm = bound.distance_mm.max(e);
                            if i < 11 {
                                let e = reduce::l1(std::array::from_fn(|a| {
                                    I::point(11. * u_parts as f64)
                                        * (difference[i + 1][j][a] - difference[i][j][a])
                                }));
                                require(e.is_finite())?;
                                bound.derivative_mm[0] = bound.derivative_mm[0].max(e);
                            }
                            if j < 11 {
                                let e = reduce::l1(std::array::from_fn(|a| {
                                    I::point(44.) * (difference[i][j + 1][a] - difference[i][j][a])
                                }));
                                require(e.is_finite())?;
                                bound.derivative_mm[1] = bound.derivative_mm[1].max(e);
                            }
                        }
                    }
                    if bound.distance_mm > options.maximum_distance_mm
                        || bound
                            .derivative_mm
                            .iter()
                            .any(|&x| x > options.maximum_cell_derivative_mm)
                    {
                        return Err(ReadError::Unsupported);
                    }
                    for (i, row) in q.into_iter().enumerate() {
                        for (j, p) in row.into_iter().enumerate() {
                            let index = ((u / 11 * u_parts + uchild) * 9 + i) * h
                                + (((v - r.entry) / 11 * 4 + child) * 9 + j);
                            if output[index][0].is_nan() {
                                output[index] = p
                            } else {
                                require(output[index] == p)?;
                            }
                        }
                    }
                }
            }
        }
    }
    require(output.iter().flatten().all(|x| x.is_finite()))?;
    Ok((
        Net {
            forward: r.forward,
            spans,
            poles: output,
        },
        bound,
    ))
}

#[cfg(test)]
mod tests;
