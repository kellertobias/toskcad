//! Bounded approximate normal stencils steer only unaccepted candidates.
//! A passing stencil always falls through to complete source certificates.
//! Uniform V refinement chooses a grid; it asserts no source periodicity.
use super::*;
#[cfg(test)]
mod tests;
impl Prepared<'_> {
    pub(super) fn normal_preflight(
        &mut self,
        mesh: Option<&DisplayMesh>,
        grid: &Grid,
        ring: &Ring,
        pass: usize,
        changes: &mut Refinement,
        work: &mut Work<'_>,
    ) -> Result<bool, Error> {
        let mut points = std::mem::take(&mut self.points);
        let result = (|| {
            let last = grid.v.len() - 2;
            let mut rows = BTreeSet::from([0, last / 2, last]);
            work.charge(4 * grid.v.len())?;
            // Refinement makes the longitudinal grid nonuniform. Index-based
            // probes alone can all land in short rows and miss its widest cell.
            let widest = grid
                .v
                .windows(2)
                .enumerate()
                .max_by(|(_, a), (_, b)| (a[1] - a[0]).total_cmp(&(b[1] - b[0])))
                .expect("constructed longitudinal grid has cells")
                .0;
            rows.insert(widest);
            'knots: for side in &self.sides {
                let Side::Tensor(tensor) = side else {
                    continue;
                };
                for &knot in &tensor.knots[1] {
                    if rows.len() >= self.limits.probe_cells.max(4) {
                        break 'knots;
                    }
                    work.charge(64 + grid.v.len().ilog2() as usize)?;
                    let at = grid.v.binary_search_by(|v| v.total_cmp(&knot));
                    let row = match at {
                        Ok(i) => i.min(last),
                        Err(i) => i.saturating_sub(1).min(last),
                    };
                    rows.insert(row);
                }
            }
            let mut representative = self.certificate_stencil(
                mesh,
                grid,
                ring,
                (&rows, pass),
                changes,
                &mut points,
                work,
            )?;
            if changes.u.iter().any(|u| !u.is_empty()) {
                // A source certificate has requested transverse recovery. Resolve
                // that proposal before a sampled V failure can double the entire
                // longitudinal grid and hide the same U defect on the next pass.
                changes.v.clear();
                return Ok(true);
            }
            let mut probes = 0;
            'probes: for &row in &rows {
                for (side, nodes) in grid.u.iter().enumerate() {
                    if !matches!(self.sides[side], Side::Tensor(_)) {
                        continue;
                    }
                    let count = nodes.len() - 1;
                    let indices = if row == 0 {
                        (0..count)
                            .step_by(count.div_ceil(8))
                            .collect::<BTreeSet<_>>()
                    } else {
                        BTreeSet::from([3 * count / 8, 5 * count / 8])
                    };
                    for u in indices {
                        if probes == self.limits.probe_cells {
                            break 'probes;
                        }
                        work.step()?;
                        probes += 1;
                        if let Some(axis) = self.probe_normal(
                            mesh,
                            grid,
                            ring,
                            [side, u, row],
                            pass,
                            &mut points,
                            work,
                        )? {
                            changes.mark(self, grid, [side, u, row], axis, work)?;
                            if axis == 1 && representative.is_none() {
                                representative = Some([side, u, row]);
                            }
                        }
                    }
                }
            }
            if changes.u.iter().any(|u| !u.is_empty()) {
                changes.v.clear();
                return Ok(true);
            }
            if let Some([side, u, failed]) = representative {
                let failed_range = [grid.v[failed], grid.v[failed + 1]];
                for row in 0..=last {
                    if same_or_coarser([grid.v[row], grid.v[row + 1]], failed_range) {
                        changes.mark(self, grid, [side, u, row], 1, work)?;
                    }
                }
            }
            if changes.empty() {
                return Ok(false);
            }
            // Unaccepted stencil cells do not establish a need for transverse
            // compensation. The complete check still inspects every literal facet
            // and can propose recovery from actual certified/failed neighbors.
            Ok(true)
        })();
        self.points = points;
        result
    }
    /// Full source certificates steer a bounded stencil only. A passing
    /// stencil never certifies the remaining grid or skips final validation.
    #[allow(clippy::too_many_arguments)]
    fn certificate_stencil(
        &mut self,
        mesh: Option<&DisplayMesh>,
        grid: &Grid,
        ring: &Ring,
        selection: (&BTreeSet<usize>, usize),
        changes: &mut Refinement,
        points: &mut BTreeMap<[u64; 3], Point>,
        work: &mut Work<'_>,
    ) -> Result<Option<[usize; 3]>, Error> {
        let (rows, pass) = selection;
        let mut representative: Option<[usize; 3]> = None;
        let mut selected = BTreeSet::new();
        let last = grid.v.len() - 2;
        'base: for row in BTreeSet::from([0, last / 2, last]) {
            for (side, nodes) in grid.u.iter().enumerate() {
                if !matches!(self.sides[side], Side::Tensor(_)) {
                    continue;
                }
                for u in 0..nodes.len() - 1 {
                    if selected.len() == self.limits.probe_cells {
                        break 'base;
                    }
                    work.step()?;
                    selected.insert([row, side, u]);
                }
            }
        }
        'extra: for round in 0..2 {
            for &row in rows {
                for (side, nodes) in grid.u.iter().enumerate() {
                    if !matches!(self.sides[side], Side::Tensor(_)) {
                        continue;
                    }
                    if selected.len() == self.limits.probe_cells {
                        break 'extra;
                    }
                    work.step()?;
                    selected.insert([row, side, if round == 0 { 0 } else { nodes.len() - 2 }]);
                }
            }
        }
        let mut caches = self
            .sides
            .iter()
            .map(|s| match s {
                Side::Tensor(t) => Some(strip::Cache::new(t, self.limits.strips)),
                Side::Cone(_) => None,
            })
            .collect::<Vec<_>>();
        for [row, side, u] in selected {
            let nodes = &grid.u[side];
            let count = nodes.len() - 1;
            let cache = caches[side]
                .as_mut()
                .expect("selected tensor owns source cache");
            cache.advance(grid.v[row], work)?;
            let candidate = self.stencil_cell(mesh, grid, ring, [side, u, row], points, work)?;
            let CandidateCell { range, actual, .. } = candidate;
            let topology = actual
                .into_iter()
                .flatten()
                .enumerate()
                .fold(0u64, |key, (i, v)| key | ((v as u64) << (2 * i)));
            let mut key = [0; 18];
            key[..6].copy_from_slice(&[
                side as u64,
                range[0][0].to_bits(),
                range[0][1].to_bits(),
                range[1][0].to_bits(),
                range[1][1].to_bits(),
                topology,
            ]);
            for (slot, value) in key[6..]
                .iter_mut()
                .zip(candidate.vertices.into_iter().flatten())
            {
                *slot = value.to_bits();
            }
            work.charge(128)?;
            if self.certificates.contains(&key) {
                continue;
            }
            match self.cell_check(candidate, None, Some(cache), work) {
                Ok(CellCheck::Match) => {
                    work.charge(128)?;
                    if self.certificates.len() < self.limits.certificates {
                        self.certificates.insert(key);
                    }
                }
                verdict => {
                    if let Err(error) = verdict
                        && error != Error::Unresolved
                    {
                        return Err(error);
                    }
                    let axis = self.hint(side, range, pass, work)?;
                    if axis == 0 {
                        let failed_range = [nodes[u], nodes[u + 1]];
                        for column in 0..count {
                            if same_or_coarser([nodes[column], nodes[column + 1]], failed_range) {
                                changes.mark(self, grid, [side, column, row], 0, work)?;
                            }
                        }
                    } else {
                        changes.mark(self, grid, [side, u, row], axis, work)?;
                        if representative.is_none_or(|[_, _, old]| {
                            grid.v[row + 1] - grid.v[row] > grid.v[old + 1] - grid.v[old]
                        }) {
                            representative = Some([side, u, row]);
                        }
                    }
                }
            }
        }
        Ok(representative)
    }
    /// Match the literal quad emitted for this grid. The final upper-U
    /// corner belongs to the next side's canonical lower seam, exactly as
    /// Ring/emit omit duplicated upper seam vertices.
    #[allow(clippy::too_many_arguments)]
    fn stencil_cell(
        &self,
        mesh: Option<&DisplayMesh>,
        grid: &Grid,
        ring: &Ring,
        cell: [usize; 3],
        points: &mut BTreeMap<[u64; 3], Point>,
        work: &mut Work<'_>,
    ) -> Result<CandidateCell, Error> {
        let [side, u, row] = cell;
        let range = [
            [grid.u[side][u], grid.u[side][u + 1]],
            [grid.v[row], grid.v[row + 1]],
        ];
        let (vertices, actual) = if let Some(mesh) = mesh {
            let corners = [
                ring.index(side, u, row),
                ring.index(side, u + 1, row),
                ring.index(side, u + 1, row + 1),
                ring.index(side, u, row + 1),
            ];
            let index = row * ring.entries.len() + ring.offsets[side] + u;
            let actual = std::array::from_fn(|tri| {
                std::array::from_fn(|k| {
                    corners
                        .iter()
                        .position(|&i| i == mesh.indices[6 * index + 3 * tri + k] as usize)
                        .expect("constructed quad corners")
                })
            });
            (corners.map(|i| point(mesh, i)), actual)
        } else {
            let uv = [
                [range[0][0], range[1][0]],
                [range[0][1], range[1][0]],
                [range[0][1], range[1][1]],
                [range[0][0], range[1][1]],
            ];
            let mut vertices = [[0.; 3]; 4];
            for (vertex, [u, v]) in vertices.iter_mut().zip(uv) {
                work.charge(64)?;
                let (owner, u) = if u == 1. {
                    ((side + 1) % self.sides.len(), 0.)
                } else {
                    (side, u)
                };
                let key = [owner as u64, u.to_bits(), v.to_bits()];
                *vertex = if let Some(&point) = points.get(&key) {
                    point
                } else {
                    let point = self.point(owner, u, v, None, work)?;
                    if points.len() < self.limits.points {
                        points.insert(key, point);
                    }
                    point
                };
            }
            (vertices, quad_triangles(vertices))
        };
        Ok(CandidateCell {
            side,
            range,
            vertices,
            actual,
        })
    }
    #[allow(clippy::too_many_arguments)]
    pub(super) fn probe_normal(
        &self,
        mesh: Option<&DisplayMesh>,
        grid: &Grid,
        ring: &Ring,
        cell: [usize; 3],
        pass: usize,
        points: &mut BTreeMap<[u64; 3], Point>,
        work: &mut Work<'_>,
    ) -> Result<Option<usize>, Error> {
        work.step()?;
        let [side, u, row] = cell;
        let Side::Tensor(t) = &self.sides[side] else {
            return Ok(None);
        };
        let range = [
            [grid.u[side][u], grid.u[side][u + 1]],
            [grid.v[row], grid.v[row + 1]],
        ];
        let uv = [
            [range[0][0], range[1][0]],
            [range[0][1], range[1][0]],
            [range[0][1], range[1][1]],
            [range[0][0], range[1][1]],
        ];
        let candidate = self.stencil_cell(mesh, grid, ring, cell, points, work)?;
        let mut reject = false;
        for local in candidate.actual {
            let center = std::array::from_fn(|k| local.iter().map(|&i| uv[i][k]).sum::<f64>() / 3.);
            let Some(normal) = t.normal_hint(center, work)? else {
                continue;
            };
            let p = local.map(|i| candidate.vertices[i]);
            let a = std::array::from_fn::<_, 3, _>(|k| p[1][k] - p[0][k]);
            let b = std::array::from_fn::<_, 3, _>(|k| p[2][k] - p[0][k]);
            let facet = [
                a[1] * b[2] - a[2] * b[1],
                a[2] * b[0] - a[0] * b[2],
                a[0] * b[1] - a[1] * b[0],
            ];
            let norm = |v: [f64; 3]| v.into_iter().map(|x| x * x).sum::<f64>().sqrt();
            let fnorm = norm(facet);
            let nnorm = norm(normal);
            if !fnorm.is_finite() || !nnorm.is_finite() || fnorm == 0. || nnorm == 0. {
                continue;
            }
            let dot = facet
                .into_iter()
                .zip(normal)
                .map(|(a, b)| a * b)
                .sum::<f64>();
            let cross = [
                facet[1] * normal[2] - facet[2] * normal[1],
                facet[2] * normal[0] - facet[0] * normal[2],
                facet[0] * normal[1] - facet[1] * normal[0],
            ];
            let cnorm = norm(cross);
            if !cnorm.is_finite() || !dot.is_finite() {
                continue;
            }
            let angle = cnorm.atan2(dot);
            if angle.is_finite() && angle > self.options.angle * 1.000001 {
                reject = true;
            }
        }
        if reject {
            Ok(Some(self.hint(side, range, pass, work)?))
        } else {
            Ok(None)
        }
    }
}

/// A planning cohort, never a source or facet acceptance bound. Adjacent
/// endpoint intervals group rounded copies of an otherwise equal grid step;
/// strict subtraction-width comparisons can leave whole nearly-equal rows
/// unrefined merely because the endpoint magnitudes differ.
fn same_or_coarser(candidate: [f64; 2], failed: [f64; 2]) -> bool {
    let width = |r: [f64; 2]| {
        I::new(r[1].next_down(), r[1].next_up()) - I::new(r[0].next_down(), r[0].next_up())
    };
    width(candidate).hi >= width(failed).lo
}
