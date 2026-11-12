//! Bounded approximate normal stencils steer only unaccepted candidates.
//! A passing stencil always falls through to complete source certificates.
//! Uniform V refinement chooses a grid; it asserts no source periodicity.
use super::*;
#[cfg(test)]
mod tests;
impl Prepared<'_> {
    pub(super) fn normal_preflight(
        &self,
        mesh: &DisplayMesh,
        grid: &Grid,
        ring: &Ring,
        pass: usize,
        changes: &mut Refinement,
        work: &mut Work<'_>,
    ) -> Result<bool, Error> {
        let last = grid.v.len() - 2;
        let rows = BTreeSet::from([0, last / 2, last]);
        let mut probes = 0;
        let mut representative = None;
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
                    if let Some(axis) =
                        self.probe_normal(mesh, grid, ring, [side, u, row], pass, work)?
                    {
                        changes.mark(self, grid, [side, u, row], axis, work)?;
                        if axis == 1 && representative.is_none() {
                            representative = Some([side, u, row]);
                        }
                    }
                }
            }
        }
        if let Some([side, u, _]) = representative {
            for row in 0..=last {
                changes.mark(self, grid, [side, u, row], 1, work)?;
            }
        }
        if changes.empty() {
            return Ok(false);
        }
        let mut compensation = Vec::new();
        'compensate: for row in rows {
            for (side, nodes) in grid.u.iter().enumerate() {
                for u in 0..nodes.len() - 1 {
                    if compensation.len() == self.limits.compensation_cells {
                        break 'compensate;
                    }
                    work.charge(4)?;
                    compensation.push([side, u, row]);
                }
            }
        }
        changes.compensate(self, grid, mesh, ring, &compensation, work)?;
        Ok(true)
    }
    fn probe_normal(
        &self,
        mesh: &DisplayMesh,
        grid: &Grid,
        ring: &Ring,
        cell: [usize; 3],
        pass: usize,
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
        let corners = [
            ring.index(side, u, row),
            ring.index(side, u + 1, row),
            ring.index(side, u + 1, row + 1),
            ring.index(side, u, row + 1),
        ];
        let cell = row * ring.entries.len() + ring.offsets[side] + u;
        let mut reject = false;
        for tri in 0..2 {
            let ids =
                std::array::from_fn::<_, 3, _>(|k| mesh.indices[6 * cell + 3 * tri + k] as usize);
            let local = ids.map(|i| {
                corners
                    .iter()
                    .position(|&c| c == i)
                    .expect("constructed quad uses its own corners")
            });
            let center = std::array::from_fn(|k| local.iter().map(|&i| uv[i][k]).sum::<f64>() / 3.);
            let Some(normal) = t.normal_hint(center, work)? else {
                continue;
            };
            let p = ids.map(|i| point(mesh, i));
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
