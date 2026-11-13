//! Indexed, source-certified tessellation of retained radial material.
use super::*;
use crate::display::{DisplayMesh, Tessellation};
use std::collections::{BTreeMap, BTreeSet};
mod cap;
mod certificate;
mod cone;
mod edge;
mod embedding;
mod normal_jet;
pub(super) mod parallel;
mod preflight;
mod spline;
mod strip;
mod tensor;
use cap::Cap;
use cone::Cone;
use spline::Spline;
use tensor::Tensor;
const MAX_TRIANGLES: usize = 4_000_000;
#[derive(Clone, Copy)]
struct Limits {
    work: usize,
    passes: usize,
    points: usize,
    certificates: usize,
    strips: usize,
    normal_strips: bool,
    compensate_v: bool,
    preflight: bool,
    probe_cells: usize,
    compensation_cells: usize,
    parallel_jobs: usize,
    parallel_sources: usize,
    parallel_bytes: usize,
    parallel_stack: usize,
}
const DEFAULT_LIMITS: Limits = Limits {
    work: 512_000_000,
    passes: 64,
    points: 131072,
    certificates: 65536,
    strips: 128,
    normal_strips: true,
    compensate_v: true,
    preflight: true,
    probe_cells: 512,
    compensation_cells: 65536,
    parallel_jobs: 4096,
    parallel_sources: 65536,
    parallel_bytes: 32 * 1024 * 1024,
    parallel_stack: 1024 * 1024,
};
#[derive(Clone, Copy, PartialEq, Eq)]
enum CellCheck {
    Match,
    Chord,
    Normal,
}
struct CandidateCell {
    side: usize,
    range: [[f64; 2]; 2],
    vertices: [Point; 4],
    actual: [[usize; 3]; 2],
}
enum Side {
    Tensor(Tensor),
    Cone(Cone),
}
struct Prepared<'a> {
    body: &'a RadialBoundary,
    sides: Vec<Side>,
    rails: Vec<Spline>,
    caps: [Cap; 2],
    options: Tessellation,
    points: BTreeMap<[u64; 3], Point>,
    certificates: BTreeSet<[u64; 18]>,
    limits: Limits,
    // Coupled-grid heuristics are qualified for closed tensor-side cycles.
    // Mixed analytic-cone cycles retain their prior grid steering.
    tensor_cycle: bool,
}
struct Grid {
    u: Vec<Vec<f64>>,
    v: Vec<f64>,
}
struct Ring {
    entries: Vec<(usize, usize)>,
    offsets: Vec<usize>,
}
struct Refinement {
    u: Vec<BTreeSet<u64>>,
    v: BTreeSet<u64>,
}
impl RadialBoundary {
    /// Sample an actual finite source edge with whole-curve chord bounds to
    /// the emitted world-space segments, including source-owned endpoints.
    pub fn sample_edge_with_interrupt(
        &self,
        edge: usize,
        deflection: f64,
        point_limit: usize,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Vec<Point>, Error> {
        edge::sample(self, edge, deflection, point_limit, interrupt)
    }
    pub fn tessellate_with_interrupt(
        &self,
        options: Tessellation,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<DisplayMesh, Error> {
        build(self, options, interrupt)
    }
}
fn build(
    body: &RadialBoundary,
    options: Tessellation,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<DisplayMesh, Error> {
    build_with_limits(body, options, interrupt, DEFAULT_LIMITS)
}
fn build_with_limits(
    body: &RadialBoundary,
    options: Tessellation,
    interrupt: Option<&dyn Fn() -> bool>,
    limits: Limits,
) -> Result<DisplayMesh, Error> {
    if !options.deflection.is_finite()
        || options.deflection <= 0.
        || !options.angle.is_finite()
        || options.angle <= 0.
        || options.angle > std::f64::consts::PI
        || options.triangle_limit == 0
    {
        return Err(Error::Invalid);
    }
    if options.triangle_limit > MAX_TRIANGLES {
        return Err(Error::Limit);
    }
    if (I::point(4.) * I::point(body.error_mm)).hi >= options.deflection {
        return Err(Error::Unresolved);
    }
    let mut work = Work {
        remaining: limits.work,
        interrupt,
    };
    work.step()?;
    let mut grid = Grid {
        u: body
            .ordered
            .caps
            .kinds
            .iter()
            .map(|kind| match kind {
                caps::SideKind::Tensor(_) => vec![0., 1.],
                caps::SideKind::Cone => vec![0., 0.5, 1.],
            })
            .collect(),
        v: vec![0., 1.],
    };
    counts(&grid, options.triangle_limit)?;
    #[cfg(test)]
    let probe = parallel_probe::start(parallel_probe::Stage::Prepare, &work);
    let mut prepared = Prepared::new(body, options, limits, &mut work)?;
    #[cfg(test)]
    parallel_probe::end(probe, &work);
    for side in &prepared.sides {
        if let Side::Cone(cone) = side {
            if cone.closure_error >= options.deflection {
                return Err(Error::Unresolved);
            }
            cone.normal_budget(options.angle)?;
        }
    }
    let radius = (I::point(body.ordered.caps.axial.chart.radius_range.lo)
        - I::point(body.error_mm))
    .positive()?;
    let shift = (I::point(options.deflection) / I::point(32.))
        .lo
        .min((radius * I::point(options.angle) / I::point(64.)).lo);
    let shift = I::point(shift).positive()?.lo;
    for pass in 0..limits.passes {
        counts(&grid, options.triangle_limit)?;
        let ring = Ring::new(&grid);
        let mut restart = false;
        for attempt in 0..3 {
            work.step()?;
            let mut changes = Refinement::new(grid.u.len());
            #[cfg(test)]
            let probe = parallel_probe::start(parallel_probe::Stage::Emit, &work);
            let mesh = prepared.emit(
                &grid,
                &ring,
                shift * (1_u32 << attempt) as f64,
                options.triangle_limit,
                &mut work,
            )?;
            #[cfg(test)]
            parallel_probe::end(probe, &work);
            if prepared.limits.preflight
                && prepared.tensor_cycle
                && mesh.indices.len() / 3 >= 32768
                && {
                    #[cfg(test)]
                    let probe = parallel_probe::start(parallel_probe::Stage::Preflight, &work);
                    let result = prepared.normal_preflight(
                        &mesh,
                        &grid,
                        &ring,
                        pass,
                        &mut changes,
                        &mut work,
                    )?;
                    #[cfg(test)]
                    parallel_probe::end(probe, &work);
                    result
                }
            {
                changes.apply(&mut grid, &mut work)?;
                restart = true;
                break;
            }
            #[cfg(test)]
            let probe = parallel_probe::start(parallel_probe::Stage::Check, &work);
            prepared.check(&mesh, &grid, &ring, pass, &mut changes, &mut work)?;
            #[cfg(test)]
            parallel_probe::end(probe, &work);
            if changes.empty() && {
                #[cfg(test)]
                let probe = parallel_probe::start(parallel_probe::Stage::Embedding, &work);
                let result = changes.embedding(&prepared, &grid, &ring, &mesh, &mut work)?;
                #[cfg(test)]
                parallel_probe::end(probe, &work);
                result
            } {
                work.step()?;
                #[cfg(test)]
                parallel_probe::finished(limits.work - work.remaining);
                return Ok(mesh);
            }
            if changes.empty() {
                continue; // Retry a larger, still certified cap-center shift.
            }
            changes.apply(&mut grid, &mut work)?;
            restart = true;
            break;
        }
        restart.then_some(()).ok_or(Error::Unresolved)?;
    }
    Err(Error::Limit)
}
fn counts(grid: &Grid, limit: usize) -> Result<(usize, usize), Error> {
    let ring = grid
        .u
        .iter()
        .try_fold(0_usize, |sum, nodes| sum.checked_add(nodes.len() - 1))
        .ok_or(Error::Limit)?;
    let vertices = ring.checked_mul(grid.v.len()).ok_or(Error::Limit)?;
    let triangles = vertices.checked_mul(2).ok_or(Error::Limit)?;
    if triangles > limit || triangles > MAX_TRIANGLES {
        return Err(Error::Limit);
    }
    Ok((vertices.checked_add(2).ok_or(Error::Limit)?, triangles))
}
impl Ring {
    fn new(grid: &Grid) -> Self {
        let mut entries = Vec::new();
        let mut offsets = Vec::new();
        for (side, nodes) in grid.u.iter().enumerate() {
            offsets.push(entries.len());
            entries.extend((0..nodes.len() - 1).map(|index| (side, index)));
        }
        Self { entries, offsets }
    }
    fn index(&self, side: usize, u: usize, row: usize) -> usize {
        row * self.entries.len() + (self.offsets[side] + u) % self.entries.len()
    }
}
impl<'a> Prepared<'a> {
    fn new(
        body: &'a RadialBoundary,
        options: Tessellation,
        limits: Limits,
        work: &mut Work<'_>,
    ) -> Result<Self, Error> {
        let rails = body
            .ordered
            .rails
            .iter()
            .map(|p| Spline::new(p.clone()))
            .collect::<Vec<_>>();
        let mut sides = Vec::with_capacity(body.side_faces().len());
        for (side, kind) in body.ordered.caps.kinds.iter().enumerate() {
            work.step()?;
            sides.push(match kind {
                caps::SideKind::Tensor(index) => Side::Tensor(Tensor::new(
                    &body.ordered.caps.axial.tensors[*index],
                    *index,
                    &body.ordered.caps.axial.chart,
                    work,
                )?),
                caps::SideKind::Cone => Side::Cone(Cone::new(body, side, &rails, work)?),
            });
        }
        let caps = [Cap::new(body, 0, work)?, Cap::new(body, 1, work)?];
        let tensor_cycle = sides.iter().all(|s| matches!(s, Side::Tensor(_)));
        Ok(Self {
            body,
            sides,
            tensor_cycle,
            rails,
            caps,
            options,
            points: BTreeMap::new(),
            certificates: BTreeSet::new(),
            limits,
        })
    }
    fn point(
        &self,
        side: usize,
        u: f64,
        v: f64,
        row: Option<&cone::Row>,
        work: &mut Work<'_>,
    ) -> Result<Point, Error> {
        let p = if v == 0. || v == 1. {
            self.caps[usize::from(v == 1.)].point(self.body, usize::from(v == 1.), side, u, work)?
        } else if u == 0. {
            let p = self.rails[side].bounds([v, v], work)?[0];
            let chart = &self.body.ordered.caps.axial.chart;
            let offset = std::array::from_fn(|k| {
                I::point(chart.bottom_point[k]) - I::point(chart.axis_origin[k])
            });
            let station = dot(offset, chart.axes[2]);
            std::array::from_fn(|k| {
                I::point(chart.axis_origin[k])
                    + chart.axes[0][k] * p[0]
                    + chart.axes[1][k] * p[1]
                    + chart.axes[2][k] * (station + I::point(chart.nominal_height) * I::point(v))
            })
        } else {
            match &self.sides[side] {
                Side::Tensor(t) => t.at_rectified(self.body, [u, v], work)?,
                Side::Cone(c) => {
                    c.jet_at_row(
                        self.body,
                        I::point(u),
                        u >= 0.5,
                        row.expect("interior cone sampling owns its row"),
                        work,
                    )?
                    .position
                }
            }
        };
        for value in p {
            value.checked()?;
        }
        Ok(p.map(|v| v.midpoint()))
    }
    fn emit(
        &mut self,
        grid: &Grid,
        ring: &Ring,
        shift: f64,
        limit: usize,
        work: &mut Work<'_>,
    ) -> Result<DisplayMesh, Error> {
        let (vertices, triangles) = counts(grid, limit)?;
        let mut mesh = DisplayMesh {
            vertices: Vec::with_capacity(3 * vertices),
            indices: Vec::with_capacity(3 * triangles),
            face_indices: Vec::with_capacity(triangles),
        };
        for &v in &grid.v {
            let mut rows = vec![None; self.sides.len()];
            for &(side, index) in &ring.entries {
                work.charge(64)?;
                let u = grid.u[side][index];
                let key = [side as u64, u.to_bits(), v.to_bits()];
                let p = if let Some(&p) = self.points.get(&key) {
                    p
                } else {
                    if v > 0.
                        && v < 1.
                        && u > 0.
                        && rows[side].is_none()
                        && let Side::Cone(c) = &self.sides[side]
                    {
                        rows[side] = Some(c.row(I::point(v), work)?);
                    }
                    let p = self.point(side, u, v, rows[side].as_ref(), work)?;
                    // Optional, request-local acceleration. Reaching a
                    // cache cap never replaces or skips the actual proof.
                    if self.points.len() < self.limits.points {
                        self.points.insert(key, p);
                    }
                    p
                };
                mesh.vertices.extend(p);
            }
        }
        for cap in &self.caps {
            work.step()?;
            mesh.vertices.extend(std::array::from_fn::<_, 3, _>(|k| {
                (cap.center[k] + cap.normal[k] * I::point(shift)).midpoint()
            }));
        }
        for row in 0..grid.v.len() - 1 {
            for (i, &(side, _)) in ring.entries.iter().enumerate() {
                work.step()?;
                let n = ring.entries.len();
                let corners = [
                    row * n + i,
                    row * n + (i + 1) % n,
                    (row + 1) * n + (i + 1) % n,
                    (row + 1) * n + i,
                ];
                let face = match &self.sides[side] {
                    Side::Tensor(t) => t.face,
                    Side::Cone(c) => c.face,
                } as u32;
                for triangle in quad_triangles(corners.map(|i| point(&mesh, i))) {
                    mesh.triangle(triangle.map(|i| corners[i]), face);
                }
            }
        }
        let center = ring.entries.len() * grid.v.len();
        for i in 0..ring.entries.len() {
            work.step()?;
            let next = (i + 1) % ring.entries.len();
            mesh.triangle([center, next, i], self.caps[0].face as u32);
            let top = (grid.v.len() - 1) * ring.entries.len();
            mesh.triangle([center + 1, top + i, top + next], self.caps[1].face as u32);
        }
        Ok(mesh)
    }
    fn check(
        &mut self,
        mesh: &DisplayMesh,
        grid: &Grid,
        ring: &Ring,
        pass: usize,
        changes: &mut Refinement,
        work: &mut Work<'_>,
    ) -> Result<(), Error> {
        let options = self.options;
        let compensate = self.limits.compensate_v && self.tensor_cycle;
        let mut accepted = Vec::new();
        #[cfg(test)]
        let probe = parallel_probe::start(parallel_probe::Stage::Cells, work);
        let parallel = parallel::check(self, mesh, grid, ring, pass, changes, &mut accepted, work)?;
        if !parallel {
            let mut strips = self
                .sides
                .iter()
                .map(|side| match side {
                    Side::Tensor(t) => Some(strip::Cache::new(t, self.limits.strips)),
                    Side::Cone(_) => None,
                })
                .collect::<Vec<_>>();
            for row in 0..grid.v.len() - 1 {
                for (side, nodes) in grid.u.iter().enumerate() {
                    if let Some(cache) = &mut strips[side] {
                        cache.advance(grid.v[row], work)?;
                    }
                    // These immutable rail enclosures depend only on this
                    // side and longitudinal row. Reuse them across U cells;
                    // every U blend and actual-facet proof remains separate.
                    let cone_rows = if let Side::Cone(c) = &self.sides[side] {
                        match c.rows([grid.v[row], grid.v[row + 1]], work) {
                            Ok(rows) => Some(rows),
                            Err(Error::Unresolved) => {
                                changes.mark(self, grid, [side, 0, row], 1, work)?;
                                continue;
                            }
                            Err(error) => return Err(error),
                        }
                    } else {
                        None
                    };
                    for u in 0..nodes.len() - 1 {
                        work.step()?;
                        let range = [[nodes[u], nodes[u + 1]], [grid.v[row], grid.v[row + 1]]];
                        let corner_indices = [
                            ring.index(side, u, row),
                            ring.index(side, u + 1, row),
                            ring.index(side, u + 1, row + 1),
                            ring.index(side, u, row + 1),
                        ];
                        let corners = corner_indices.map(|i| point(mesh, i));
                        let cell = row * ring.entries.len() + ring.offsets[side] + u;
                        let actual = std::array::from_fn::<_, 2, _>(|triangle| {
                            std::array::from_fn::<_, 3, _>(|k| {
                                let index = mesh.indices[6 * cell + 3 * triangle + k] as usize;
                                corner_indices
                                    .iter()
                                    .position(|&corner| corner == index)
                                    .expect("constructed quad uses its own corners")
                            })
                        });
                        let topology = actual
                            .into_iter()
                            .flatten()
                            .enumerate()
                            .fold(0_u64, |key, (i, v)| key | ((v as u64) << (2 * i)));
                        let mut key = [0_u64; 18];
                        key[..6].copy_from_slice(&[
                            side as u64,
                            range[0][0].to_bits(),
                            range[0][1].to_bits(),
                            range[1][0].to_bits(),
                            range[1][1].to_bits(),
                            topology,
                        ]);
                        for (slot, value) in key[6..].iter_mut().zip(corners.into_iter().flatten())
                        {
                            *slot = value.to_bits();
                        }
                        // Only successful full certificates are retained. The
                        // immutable owner/options and exact UV/topology/point
                        // key prevent reuse across a different actual facet.
                        work.charge(128)?;
                        if self.certificates.contains(&key) {
                            if compensate && accepted.len() < self.limits.compensation_cells {
                                work.charge(4)?;
                                accepted.push([side, u, row]);
                            }
                            continue;
                        }
                        let verdict = match self.cell_check(
                            CandidateCell {
                                side,
                                range,
                                vertices: corners,
                                actual,
                            },
                            cone_rows.as_ref(),
                            strips[side].as_mut(),
                            work,
                        ) {
                            Ok(v) => Some(v),
                            Err(Error::Unresolved) => None,
                            Err(e) => return Err(e),
                        };
                        if compensate
                            && verdict == Some(CellCheck::Normal)
                            && accepted.len() < self.limits.compensation_cells
                        {
                            work.charge(4)?;
                            accepted.push([side, u, row]);
                        }
                        if verdict == Some(CellCheck::Match) {
                            if compensate && accepted.len() < self.limits.compensation_cells {
                                work.charge(4)?;
                                accepted.push([side, u, row]);
                            }
                            work.charge(128)?;
                            if self.certificates.len() < self.limits.certificates {
                                self.certificates.insert(key);
                            }
                        } else {
                            let axis = self.hint(side, range, pass, work)?;
                            changes.mark(self, grid, [side, u, row], axis, work)?;
                        }
                    }
                }
            }
        }
        #[cfg(test)]
        parallel_probe::end(probe, work);
        if compensate {
            changes.compensate(self, grid, mesh, ring, &accepted, work)?;
        }
        let centers = [
            point(mesh, ring.entries.len() * grid.v.len()),
            point(mesh, ring.entries.len() * grid.v.len() + 1),
        ];
        for (end, &center) in centers.iter().enumerate() {
            let center_error = (0..3)
                .fold(I::point(0.), |sum, k| {
                    sum + I::point((self.caps[end].center[k] - I::point(center[k])).magnitude())
                })
                .checked()?
                .hi;
            if (I::point(center_error) + I::point(self.body.error_mm)).hi > options.deflection {
                return Err(Error::Unresolved);
            }
            for (side, nodes) in grid.u.iter().enumerate() {
                for u in 0..nodes.len() - 1 {
                    work.step()?;
                    let row = if end == 0 { 0 } else { grid.v.len() - 1 };
                    let edge = [
                        point(mesh, ring.index(side, u, row)),
                        point(mesh, ring.index(side, u + 1, row)),
                    ];
                    let triangle = if end == 0 {
                        [center, edge[1], edge[0]]
                    } else {
                        [center, edge[0], edge[1]]
                    };
                    let error =
                        self.caps[end].curve_error(side, [nodes[u], nodes[u + 1]], edge, work)?;
                    let chord =
                        (I::point(error.max(center_error)) + I::point(self.body.error_mm)).hi;
                    let normal = uncertain(certificate::normal_matches(
                        vec![vec![self.caps[end].normal]],
                        [[0., 1.]; 2],
                        [[0., 1.]; 2],
                        certificate::facet_normal(triangle),
                        options.angle,
                        work,
                    ))?;
                    if chord > options.deflection || !normal {
                        changes.mark(self, grid, [side, u, 0], 0, work)?;
                    }
                }
            }
        }
        Ok(())
    }
    #[allow(clippy::too_many_arguments)]
    #[cfg(test)]
    fn cell_matches(
        &self,
        side: usize,
        range: [[f64; 2]; 2],
        vertices: [Point; 4],
        actual: [[usize; 3]; 2],
        cone_rows: Option<&cone::Rows>,
        strips: Option<&mut strip::Cache<'_>>,
        work: &mut Work<'_>,
    ) -> Result<bool, Error> {
        self.cell_check(
            CandidateCell {
                side,
                range,
                vertices,
                actual,
            },
            cone_rows,
            strips,
            work,
        )
        .map(|v| v == CellCheck::Match)
    }
    fn cell_check(
        &self,
        cell: CandidateCell,
        cone_rows: Option<&cone::Rows>,
        mut strips: Option<&mut strip::Cache<'_>>,
        work: &mut Work<'_>,
    ) -> Result<CellCheck, Error> {
        let CandidateCell {
            side,
            range,
            vertices,
            actual,
        } = cell;
        let options = self.options;
        let uv = [
            [range[0][0], range[1][0]],
            [range[0][1], range[1][0]],
            [range[0][1], range[1][1]],
            [range[0][0], range[1][1]],
        ];
        match &self.sides[side] {
            Side::Tensor(t) => {
                // Both literal triangles cover this same immutable source
                // cell. Reuse only source-dependent enclosures; each facet
                // still owns its affine chord and normal certificates.
                let mut patches = t
                    .patches(range, work)?
                    .into_iter()
                    .map(|patch| (patch, None, None, None))
                    .collect::<Vec<_>>();
                for indices in actual {
                    let triangle = indices.map(|i| vertices[i]);
                    let affine = certificate::Affine::new(indices.map(|i| uv[i]), triangle)?;
                    let normal = certificate::facet_normal(triangle);
                    for (patch, source_normal, source_cell, normal_coefficients) in &mut patches {
                        let cell = if let Some(cell) = source_cell {
                            cell
                        } else {
                            source_cell.insert(if let Some(cache) = strips.as_deref_mut() {
                                cache.cell(patch, range, work)?
                            } else {
                                certificate::SourceCell::new(
                                    &patch.controls,
                                    patch.domain,
                                    range,
                                    work,
                                )?
                            })
                        };
                        if (I::point(cell.chord_error(&affine, work)?)
                            + I::point(self.body.error_mm))
                        .hi > options.deflection
                        {
                            return Ok(CellCheck::Chord);
                        }
                        if uncertain(patch.taylor_normal_box(range, work).and_then(|value| {
                            certificate::normal_matches(
                                vec![vec![value]],
                                [[0., 1.]; 2],
                                [[0., 1.]; 2],
                                normal,
                                options.angle,
                                work,
                            )
                        }))? {
                            continue;
                        }
                        let source_normal = match *source_normal {
                            Some(value) => value,
                            None => {
                                let value = patch.normal_box(range, work);
                                *source_normal = Some(value);
                                value
                            }
                        };
                        let fast = uncertain(source_normal.and_then(|value| {
                            certificate::normal_matches(
                                vec![vec![value]],
                                [[0., 1.]; 2],
                                [[0., 1.]; 2],
                                normal,
                                options.angle,
                                work,
                            )
                        }))?;
                        if !fast {
                            let matches = if self.limits.normal_strips {
                                let coefficients = if let Some(coefficients) = normal_coefficients {
                                    coefficients
                                } else {
                                    normal_coefficients.insert(
                                        if let Some(cache) = strips.as_deref_mut() {
                                            cache.normal_coefficients(patch, range, work)?
                                        } else {
                                            certificate::SourceStrip::new(
                                                &patch.normal,
                                                patch.domain,
                                                range[0],
                                                work,
                                            )?
                                            .normal_coefficients(range[1], work)?
                                        },
                                    )
                                };
                                certificate::normal_coefficients_match(
                                    coefficients,
                                    normal,
                                    options.angle,
                                    work,
                                )?
                            } else {
                                certificate::normal_matches(
                                    patch.normal.clone(),
                                    patch.domain,
                                    range,
                                    normal,
                                    options.angle,
                                    work,
                                )?
                            };
                            if !matches {
                                return Ok(CellCheck::Normal);
                            }
                        }
                    }
                }
            }
            Side::Cone(c) => {
                let (error, normal) = c.cell_error(
                    self.body,
                    range[0],
                    cone_rows.expect("cone rows bound to this side"),
                    vertices,
                    range[0][0] >= 0.5,
                    work,
                )?;
                if error > options.deflection {
                    return Ok(CellCheck::Chord);
                }
                for indices in actual {
                    if !certificate::normal_matches(
                        vec![vec![normal]],
                        [[0., 1.]; 2],
                        [[0., 1.]; 2],
                        certificate::facet_normal(indices.map(|i| vertices[i])),
                        c.normal_budget(options.angle)?,
                        work,
                    )? {
                        return Ok(CellCheck::Normal);
                    }
                }
            }
        }
        Ok(CellCheck::Match)
    }
    fn hint(
        &self,
        side: usize,
        range: [[f64; 2]; 2],
        pass: usize,
        work: &mut Work<'_>,
    ) -> Result<usize, Error> {
        let options = self.options;
        let score = match &self.sides[side] {
            Side::Tensor(t) => t.hint(range, options.deflection, options.angle, work)?,
            Side::Cone(c) => match c.jet(
                self.body,
                range.map(|d| I::new(d[0], d[1])),
                range[0][0] >= 0.5,
                work,
            ) {
                Ok(jet) => std::array::from_fn(|i| {
                    jet.angular_rate[i].magnitude() * (range[i][1] - range[i][0])
                }),
                Err(Error::Unresolved) => return Ok(1),
                Err(error) => return Err(error),
            },
        };
        Ok(refinement_axis(score, pass))
    }
}
fn refinement_axis(score: [f64; 2], pass: usize) -> usize {
    if score[0] == score[1] {
        pass % 2
    } else {
        usize::from(score[1] > score[0])
    }
}
impl Refinement {
    // A global V split can spoil a previously sound, highly skew quad on
    // another band. Physical edge ratios only propose compensating U knots;
    // every resulting actual facet still needs its complete source proof.
    fn compensate(
        &mut self,
        source: &Prepared<'_>,
        grid: &Grid,
        mesh: &DisplayMesh,
        ring: &Ring,
        cells: &[[usize; 3]],
        work: &mut Work<'_>,
    ) -> Result<(), Error> {
        if self.v.is_empty() {
            return Ok(());
        }
        for &[side, u, row] in cells {
            work.charge(64)?;
            let vparts = 1 + self
                .v
                .range(grid.v[row].to_bits()..grid.v[row + 1].to_bits())
                .count();
            if vparts == 1 {
                continue;
            }
            let uparts = 1 + self.u[side]
                .range(grid.u[side][u].to_bits()..grid.u[side][u + 1].to_bits())
                .count();
            let p = [
                ring.index(side, u, row),
                ring.index(side, u + 1, row),
                ring.index(side, u + 1, row + 1),
                ring.index(side, u, row + 1),
            ]
            .map(|i| point(mesh, i));
            let length = |a: Point, b: Point| {
                a.into_iter()
                    .zip(b)
                    .map(|(x, y)| (x - y).powi(2))
                    .sum::<f64>()
                    .sqrt()
            };
            let uedge = length(p[0], p[1]).max(length(p[3], p[2])) / uparts as f64;
            let vedge = length(p[0], p[3]).min(length(p[1], p[2])) / vparts as f64;
            if uedge.is_finite() && vedge.is_finite() && vedge > 0. && uedge > 2. * vedge {
                self.mark(source, grid, [side, u, row], 0, work)?;
            }
        }
        Ok(())
    }

    /// Plan recovery from the literal candidate proof. This does not
    /// certify source approximation: the owner must finish its complete
    /// source checks before it can publish any accepted candidate.
    fn embedding(
        &mut self,
        source: &Prepared<'_>,
        grid: &Grid,
        ring: &Ring,
        mesh: &DisplayMesh,
        work: &mut Work<'_>,
    ) -> Result<bool, Error> {
        match embedding::certify(mesh, ring.entries.len(), grid.v.len(), work)? {
            embedding::Verdict::Certified => return Ok(true),
            embedding::Verdict::HeightUncertain => {}
            embedding::Verdict::AngularUncertain(index) => {
                let (side, u) = ring.entries[index];
                self.mark(source, grid, [side, u, 0], 0, work)?;
            }
            embedding::Verdict::RadialUncertain(triangle) => {
                let cell = triangle / 2;
                let row = cell / ring.entries.len();
                let (side, u) = ring.entries[cell % ring.entries.len()];
                self.mark(source, grid, [side, u, row], 0, work)?;
                self.mark(source, grid, [side, u, row], 1, work)?;
            }
        }
        Ok(false)
    }
    fn new(sides: usize) -> Self {
        Self {
            u: vec![BTreeSet::new(); sides],
            v: BTreeSet::new(),
        }
    }
    fn empty(&self) -> bool {
        self.v.is_empty() && self.u.iter().all(BTreeSet::is_empty)
    }
    fn mark(
        &mut self,
        source: &Prepared<'_>,
        grid: &Grid,
        cell: [usize; 3],
        axis: usize,
        work: &mut Work<'_>,
    ) -> Result<(), Error> {
        let [side, u, v] = cell;
        work.charge(64)?;
        let range = if axis == 0 {
            [grid.u[side][u], grid.u[side][u + 1]]
        } else {
            [grid.v[v], grid.v[v + 1]]
        };
        let split = match &source.sides[side] {
            Side::Tensor(t) => choose_split(&t.knots[axis], range),
            Side::Cone(c) => {
                if axis == 0 {
                    range[0] + (range[1] - range[0]) * 0.5
                } else {
                    c.split_at(range)
                }
            }
        };
        if !split.is_finite() || split <= range[0] || split >= range[1] {
            return Err(Error::Unresolved);
        }
        if axis == 0 {
            self.u[side].insert(split.to_bits());
        } else {
            self.v.insert(split.to_bits());
        }
        Ok(())
    }
    fn apply(self, grid: &mut Grid, work: &mut Work<'_>) -> Result<(), Error> {
        for (nodes, added) in grid.u.iter_mut().zip(self.u) {
            *nodes = merge(nodes, added, work)?;
        }
        grid.v = merge(&grid.v, self.v, work)?;
        Ok(())
    }
}
fn choose_split(knots: &[f64], range: [f64; 2]) -> f64 {
    let middle = range[0] + (range[1] - range[0]) * 0.5;
    let at = knots.partition_point(|&k| k < middle);
    [at.checked_sub(1), Some(at)]
        .into_iter()
        .flatten()
        .filter_map(|i| knots.get(i).copied())
        .filter(|&x| x > range[0] && x < range[1])
        .min_by(|a, b| (a - middle).abs().total_cmp(&(b - middle).abs()))
        .unwrap_or(middle)
}
fn merge(old: &[f64], added: BTreeSet<u64>, work: &mut Work<'_>) -> Result<Vec<f64>, Error> {
    let capacity = old.len().checked_add(added.len()).ok_or(Error::Limit)?;
    let mut result = Vec::with_capacity(capacity);
    let mut a = old.iter().copied().peekable();
    let mut b = added.into_iter().map(f64::from_bits).peekable();
    while a.peek().is_some() || b.peek().is_some() {
        work.step()?;
        let value = match (a.peek(), b.peek()) {
            (Some(&x), Some(&y)) if x < y => a.next().expect("peeked original node"),
            (Some(&x), Some(&y)) if x == y => {
                b.next();
                a.next().expect("peeked equal node")
            }
            (_, Some(_)) => b.next().expect("peeked refinement node"),
            _ => a.next().expect("at least one remaining node"),
        };
        result.push(value);
    }
    Ok(result)
}
fn point(mesh: &DisplayMesh, index: usize) -> Point {
    std::array::from_fn(|k| mesh.vertices[3 * index + k])
}
fn uncertain(result: Result<bool, Error>) -> Result<bool, Error> {
    match result {
        Err(Error::Unresolved) => Ok(false),
        other => other,
    }
}

/// A geometric sampling heuristic only. The actual selected triangles
/// still undergo all whole-source and literal embedding certificates.
fn quad_triangles(points: [Point; 4]) -> [[usize; 3]; 2] {
    let length = |a: usize, b: usize| {
        (0..3)
            .map(|k| (points[a][k] - points[b][k]).powi(2))
            .sum::<f64>()
    };
    if length(1, 3) < length(0, 2) {
        [[1, 3, 0], [1, 2, 3]]
    } else {
        [[0, 1, 2], [0, 2, 3]]
    }
}

#[cfg(test)]
mod tests;

#[cfg(test)]
pub(super) mod parallel_probe;
