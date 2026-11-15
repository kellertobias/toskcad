//! Original finite coedge/curve ownership and bounded transfer to the DISTINCT
//! common-chart periodic source. No exact raw-file incidence or material flag.
use super::*;
mod circle;
type FResult<T> = std::result::Result<T, PhillipsReadError>;
const FINITE_TOTAL: usize = 48_000_000;
const FINITE_SCRATCH: usize = 8 * 1024 * 1024;
#[derive(Clone, Copy)]
pub struct PeriodicFiniteOptions {
    pub maximum_distance_mm: f64,
    pub maximum_cell_derivative_mm: f64,
    pub maximum_work: usize,
    pub maximum_total_work: usize,
    pub maximum_bytes: usize,
}
impl Default for PeriodicFiniteOptions {
    fn default() -> Self {
        Self {
            maximum_distance_mm: 1e-7,
            maximum_cell_derivative_mm: 1e-7,
            maximum_work: WORK,
            maximum_total_work: FINITE_TOTAL,
            maximum_bytes: BYTES,
        }
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PeriodicFiniteKind {
    Polynomial,
    Circle,
}
#[derive(Debug)]
pub struct PeriodicFiniteCurve {
    pub face_entity: u32,
    pub edge_entity: u32,
    pub curve_entity: u32,
    pub role: usize,
    pub forward_in_face: bool,
    pub curve_same_sense: bool,
    pub vertex_entities: [u32; 2],
    pub vertex_point_entities: [u32; 2],
    pub kind: PeriodicFiniteKind,
    pub original_curve_controls: Vec<u32>,
    pub original_surface_controls: Vec<u32>,
    pub trim_control_indices: [usize; 2],
    pub distance_mm: f64,
    /// Polynomial: each original normalized Bernstein cell. Circle: the
    /// complete oriented finite arc's normalized parameter.
    pub cell_derivative_mm: f64,
    /// Polynomial vertex/control mismatch, or explicit projection of an
    /// original circle endpoint vertex onto its actual circle support.
    pub endpoint_incidence_reserve_mm: f64,
}
pub struct PeriodicFiniteCurves<'s, 'a, 'input> {
    thread: &'s CommonChartPeriodicThreadRepair<'a, 'input>,
    curves: Vec<PeriodicFiniteCurve>,
    stations: [f64; 2],
    new_work: usize,
    retained: usize,
}
impl<'s, 'a, 'input> PeriodicFiniteCurves<'s, 'a, 'input> {
    pub fn thread(&self) -> &'s CommonChartPeriodicThreadRepair<'a, 'input> {
        self.thread
    }
    pub fn curves(&self) -> &[PeriodicFiniteCurve] {
        &self.curves
    }
    pub fn cap_stations_mm(&self) -> [f64; 2] {
        self.stations
    }
    pub fn new_work_used(&self) -> usize {
        self.new_work
    }
    pub fn inherited_work_used(&self) -> usize {
        self.thread.total_work_used()
    }
    pub fn total_work_used(&self) -> usize {
        self.new_work + self.inherited_work_used()
    }
    pub fn retained_budget_bytes(&self) -> usize {
        self.retained
    }
}
fn valid(o: PeriodicFiniteOptions) -> bool {
    o.maximum_distance_mm.is_finite()
        && o.maximum_distance_mm > 0.
        && o.maximum_distance_mm <= 1e-7
        && o.maximum_cell_derivative_mm.is_finite()
        && o.maximum_cell_derivative_mm > 0.
        && o.maximum_cell_derivative_mm <= 1e-7
        && o.maximum_work <= WORK
        && o.maximum_total_work <= FINITE_TOTAL
        && o.maximum_bytes <= BYTES
}
fn point(body: &Body<'_>, id: u32, work: &mut Work<'_>) -> FResult<[I; 3]> {
    let p = body.exact_point_leaf(id, work.remaining, work.interrupt)?;
    work.step(p.work_used())?;
    Ok(p.coordinate_bounds_mm().map(|x| I::new(x[0], x[1])))
}
fn numeric(body: &Body<'_>, id: u32, work: &mut Work<'_>) -> FResult<Vec<I>> {
    let p = body.original_numeric_leaves(id, work.remaining, work.interrupt)?;
    work.step(p.work_used())?;
    Ok(p.leaves()
        .iter()
        .map(|p| {
            let b = p.bounds();
            I::new(b[0], b[1])
        })
        .collect())
}
fn polynomial(
    body: &Body<'_>,
    shell: &Shell<'_, '_>,
    record: &mut PeriodicFiniteCurve,
    field: &AnalyticThreadTransfer,
    work: &mut Work<'_>,
) -> FResult<()> {
    let args = shell.args(record.curve_entity, "B_SPLINE_CURVE_WITH_KNOTS")?;
    let ids = refs(&args[2])?;
    assert_eq!(ids.len(), record.original_surface_controls.len());
    // Body::outer_carrier_bounds called outer::integers::validate on EVERY
    // actual side-face coedge before this head receipt. It proved the original
    // INTEGER degree/multiplicities; helical::trimmed_spline_edge proved the exact
    // list counts and increasing knots. This measurement uses those same entity
    // IDs and original normalized cells, never rounded integer metadata.
    assert!(record.curve_same_sense);
    let mut differences = Vec::with_capacity(ids.len());
    let mut distance = 0_f64;
    for (&id, &surface) in ids.iter().zip(&record.original_surface_controls) {
        work.step(200)?;
        let a = point(body, id, work)?;
        let b = point(body, surface, work)?;
        let d = std::array::from_fn(|k| a[k] - b[k]);
        distance = distance.max(reduce::l1(d));
        differences.push(d);
    }
    let derivative = differences
        .windows(2)
        .map(|p| reduce::l1(std::array::from_fn(|k| I::point(11.) * (p[1][k] - p[0][k]))))
        .fold(0_f64, f64::max);
    let mut endpoint = 0_f64;
    for (k, &index) in record.trim_control_indices.iter().enumerate() {
        let a = point(body, record.vertex_point_entities[k], work)?;
        let b = point(body, ids[index], work)?;
        endpoint = endpoint.max(reduce::l1(std::array::from_fn(|k| a[k] - b[k])));
    }
    record.kind = PeriodicFiniteKind::Polynomial;
    record.original_curve_controls = ids;
    record.distance_mm = (I::point(field.distance_mm) + I::point(distance) + I::point(endpoint)).hi;
    record.cell_derivative_mm =
        (I::point(field.derivative_mm[usize::from(!record.role.is_multiple_of(2))])
            + I::point(derivative))
        .hi;
    record.endpoint_incidence_reserve_mm = endpoint;
    Ok(())
}
impl<'h, 'q, 'r, 'a, 'input> ClosedRepairedHead<'h, 'q, 'r, 'a, 'input> {
    pub fn periodic_thread_finite_curves<'s>(
        &self,
        thread: &'s CommonChartPeriodicThreadRepair<'a, 'input>,
        options: PeriodicFiniteOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> FResult<PeriodicFiniteCurves<'s, 'a, 'input>> {
        let mut work = Work {
            remaining: options.maximum_work,
            interrupt,
            cancelled: Cell::new(false),
        };
        let result = (|| {
            work.step(1)?;
            if !valid(options) {
                return Err(PhillipsReadError::Invalid);
            }
            if !std::ptr::eq(self.body(), thread.body()) || self.chart() != thread.chart() {
                return Err(PhillipsReadError::Invalid);
            }
            let allowance = options.maximum_work.min(
                options
                    .maximum_total_work
                    .checked_sub(thread.total_work_used())
                    .ok_or(PhillipsReadError::Limit)?,
            );
            work.remaining = allowance.checked_sub(1).ok_or(PhillipsReadError::Limit)?;
            let base = thread.retained_budget_bytes()
                + std::mem::size_of::<PeriodicFiniteCurves<'_, '_, '_>>();
            // Exact dimensions are already earned by the same immutable field receipt.
            // Preflight all sixteen records and up to2x curve-ID capacity plus surface
            // IDs before rereading raw grids. The fixed8MiB scratch covers raw-ID grids
            // (262144 controls, >=12 rows/columns), numeric-leaf collectors and one
            // polynomial difference array together. No input-dependent late allocation.
            let n = thread.transfers().len() * 4;
            let payload = n * std::mem::size_of::<PeriodicFiniteCurve>()
                + thread
                    .transfers()
                    .iter()
                    .map(|t| {
                        12 * (2 * (11 * t.original_cells[0] + 1)
                            + 2 * (11 * t.original_cells[1] + 1))
                    })
                    .sum::<usize>();
            let retained = base.saturating_add(payload);
            if retained.saturating_add(FINITE_SCRATCH) > options.maximum_bytes {
                return Err(PhillipsReadError::Limit);
            }
            let body = self.body();
            let mut shell = Shell::new(&body.entities, body.tolerance, true);
            let raw = raw(self, &mut shell, &mut work, FINITE_SCRATCH)?;
            let count: usize = raw.iter().map(|r| r.ids.len() * r.ids[0].len()).sum();
            // The supplied opaque common-chart receipt already counted/capped
            // these SAME immutable Body grids before publication.
            assert!(count <= MAX_POINTS);
            let factor_count: usize =
                raw.iter().map(|r| r.ids.len()).sum::<usize>() + raw[0].ids[0].len();
            work.step(factor_count * 1024)?;
            let splines = crate::step::helical::prepare_retained_reference(
                &body.thread,
                body.visible[0],
                body.tolerance,
            )
            .map_err(helical::step_error)?;
            let entry = splines.heights[raw[0].entry];
            let exit = *splines
                .heights
                .last()
                .expect("validated nonempty axial basis");
            let mut curves = Vec::with_capacity(raw.len() * 4);
            for (index, r) in raw.iter().enumerate() {
                assert_eq!(splines.heights[r.entry], entry);
                let (_, forward, bounds) = helical::face(&shell, r.face)?;
                let uses = helical::orient(&helical::ring(&shell, bounds[0], true)?, forward);
                // The exact four-use order is a sealed inherited prefix invariant.
                assert_eq!(uses.len(), 4);
                assert_eq!([uses[3].0, uses[1].0], r.rails);
                for (role, (edge, forward_in_face)) in uses.into_iter().enumerate() {
                    work.step(256)?;
                    let a = shell.args(edge, "EDGE_CURVE")?;
                    let vertices = [reference(&a[1])?, reference(&a[2])?];
                    let geometry = reference(&a[3])?;
                    let same = sense(&a[4])?;
                    let vertex_points = vertices
                        .map(|v| shell.args(v, "VERTEX_POINT").and_then(|a| reference(&a[1])))
                        .into_iter()
                        .collect::<Result<Vec<_>>>()?;
                    let surface = match role {
                        0 => r.ids.iter().map(|row| row[r.entry]).collect(),
                        1 => r.ids.last().expect("validated radial rows").clone(),
                        2 => r
                            .ids
                            .iter()
                            .map(|row| *row.last().expect("validated axial columns"))
                            .collect(),
                        _ => r.ids[0].clone(),
                    };
                    let trim = if !role.is_multiple_of(2) {
                        [r.entry, surface.len() - 1]
                    } else {
                        [0, surface.len() - 1]
                    };
                    let mut record = PeriodicFiniteCurve {
                        face_entity: r.face,
                        edge_entity: edge,
                        curve_entity: geometry,
                        role,
                        forward_in_face,
                        curve_same_sense: same,
                        vertex_entities: vertices,
                        vertex_point_entities: [vertex_points[0], vertex_points[1]],
                        kind: PeriodicFiniteKind::Polynomial,
                        original_curve_controls: Vec::new(),
                        original_surface_controls: surface,
                        trim_control_indices: trim,
                        distance_mm: 0.,
                        cell_derivative_mm: 0.,
                        endpoint_incidence_reserve_mm: 0.,
                    };
                    // outer::Body::outer_carrier_bounds verified every exact
                    // face/bound/coedge ID reread above: verify_prefix requires
                    // trimmed_spline_edge for all rails and every exit rim
                    // (exit_crests=false), and spline_edge for noncrest entry
                    // rims. The sole skipped entry crest goes through
                    // profiled::junction -> tangent_shank::arc, requiring this
                    // exact EDGE_CURVE geometry to be CIRCLE. No other retained
                    // kind can exist in this sealed receipt; source-facing
                    // schemas/INTEGER/ownership checks remain in those paths.
                    let kind = shell.simple(geometry)?.name;
                    if kind == "CIRCLE" {
                        circle::qualify(
                            body,
                            &shell,
                            thread,
                            index,
                            [entry, exit],
                            &mut record,
                            &mut work,
                        )?;
                    } else {
                        assert_eq!(kind, "B_SPLINE_CURVE_WITH_KNOTS");
                        polynomial(
                            body,
                            &shell,
                            &mut record,
                            &thread.transfers()[index],
                            &mut work,
                        )?;
                    }

                    if record.distance_mm > options.maximum_distance_mm
                        || record.cell_derivative_mm > options.maximum_cell_derivative_mm
                    {
                        return Err(PhillipsReadError::Unsupported);
                    }
                    curves.push(record);
                }
            }
            work.step(1)?;
            Ok(PeriodicFiniteCurves {
                thread,
                curves,
                stations: [entry, exit],
                new_work: allowance - work.remaining,
                retained,
            })
        })();
        if work.cancelled.get() {
            Err(PhillipsReadError::Cancelled)
        } else {
            result
        }
    }
}
#[cfg(test)]
mod tests;
