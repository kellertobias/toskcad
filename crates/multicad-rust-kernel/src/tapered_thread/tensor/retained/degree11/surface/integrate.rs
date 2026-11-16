use super::*;
struct Cell {
    geometry: TensorCell11,
    normal: Arc<NormalPolynomial11>,
    value: Values,
    depth: u8,
}
// A degree11 cell owns144 interval triples plus bounded row headers. 16KiB
// covers its actual capacities, including a subdivision's old/new overlap.
const CELL_BYTES: usize = 16 * 1024;
// Two22x22 restriction grids, root copies, derivative grids through order8,
// source position jets, Gaussian arrays and Series temporaries coexist within
// this512KiB envelope. Original normal storage is charged separately.
const SCRATCH: usize = 512 * 1024;
fn reservation(
    source: &RetainedTensor11,
    capacity: usize,
    normals: usize,
    live: usize,
) -> Result<usize, Error> {
    source
        .poles
        .len()
        .checked_mul(std::mem::size_of::<Point>())
        .and_then(|n| n.checked_add(std::mem::size_of::<RetainedTensor11>()))
        // Source Arc counters/alignment and the fixed calculation ledger.
        .and_then(|n| n.checked_add(64))
        .and_then(|n| n.checked_add(capacity.checked_mul(std::mem::size_of::<Cell>())?))
        .and_then(|n| n.checked_add(normals.checked_mul(NORMAL11_BYTES + 64)?))
        .and_then(|n| n.checked_add(live.checked_mul(CELL_BYTES)?))
        .and_then(|n| n.checked_add(SCRATCH))
        .ok_or(Error::Limit)
}
fn output(value: Values, center: Point, options: SurfaceOptions11) -> Option<SurfaceProperties11> {
    let area = value[0];
    if !area.valid() || area.lo <= 0. {
        return None;
    }
    let raw = std::array::from_fn::<_, 3, _>(|k| value[k + 1] + I::point(center[k]) * area);
    let centroid = std::array::from_fn::<_, 3, _>(|k| I::point(center[k]) + value[k + 1] / area);
    if area.width() > options.area_width_mm2
        || raw
            .iter()
            .any(|x| !x.valid() || x.width() > options.first_moment_width_mm3)
        || centroid
            .iter()
            .any(|x| !x.valid() || x.width() > options.centroid_width_mm)
    {
        return None;
    }
    Some(SurfaceProperties11 {
        area_mm2: pair(area),
        first_moment_mm3: raw.map(pair),
        centroid_mm: centroid.map(pair),
        centered_first_moment_mm3: std::array::from_fn(|k| pair(value[k + 1])),
        center_mm: center,
        work_used: 0,
        remaining_work: 0,
        accounted_peak_bytes: 0,
        evaluated_cells: 0,
    })
}
fn score(cell: &Cell, center: Point, options: SurfaceOptions11) -> f64 {
    // Refinement priority only; final acceptance always uses the complete
    // recomposed intervals, including centroid division and world origin.
    let mut score = cell.value[0].width() / options.area_width_mm2;
    for (k, origin) in center.iter().enumerate() {
        score = score.max(
            (cell.value[k + 1].width() + origin.abs() * cell.value[0].width())
                / options.first_moment_width_mm3,
        );
    }
    score
}
pub(super) fn calculate(
    source: &RetainedTensor11,
    options: SurfaceOptions11,
    budget: &mut Budget<'_>,
) -> Result<SurfaceProperties11, Error> {
    let roots = source.spans[0]
        .checked_mul(source.spans[1])
        .ok_or(Error::Limit)?;
    if roots > options.max_cells {
        return Err(Error::Limit);
    }
    budget.bytes(
        reservation(source, options.max_cells, 0, 0)?,
        options.max_bytes,
    )?;
    budget.charge(source.poles.len().checked_mul(6).ok_or(Error::Limit)? + 128)?;
    let mut bounds = [[f64::INFINITY, f64::NEG_INFINITY]; 3];
    for p in source.poles.iter() {
        for k in 0..3 {
            bounds[k][0] = bounds[k][0].min(p[k]);
            bounds[k][1] = bounds[k][1].max(p[k]);
        }
    }
    let center = bounds.map(|b| b[0] + (b[1] - b[0]) * 0.5);
    let mut cells = Vec::with_capacity(options.max_cells);
    for u in 0..source.spans[0] {
        for v in 0..source.spans[1] {
            let n = cells.len() + 1;
            budget.bytes(
                reservation(source, cells.capacity(), n, n + 3)?,
                options.max_bytes,
            )?;
            budget.charge(NORMAL11_WORK)?;
            let geometry = source.cell([u, v])?;
            let normal = Arc::new(NormalPolynomial11::new(&geometry)?);
            let value = quadrature::enclosure(&geometry, &normal, center, budget)?;
            cells.push(Cell {
                geometry,
                normal,
                value,
                depth: 0,
            });
        }
    }
    loop {
        budget.charge(cells.len().checked_mul(128).ok_or(Error::Limit)? + 128)?;
        let mut total = [I::point(0.); 4];
        // Rebuild every iteration: interval subtraction never accumulates a
        // fictitious historical uncertainty floor in the retained live sum.
        for cell in &cells {
            for (sum, value) in total.iter_mut().zip(cell.value) {
                *sum = *sum + value;
            }
        }
        if let Some(mut result) = output(total, center, options) {
            result.work_used = options.max_work - budget.remaining;
            result.remaining_work = budget.remaining;
            result.accounted_peak_bytes = budget.peak;
            result.evaluated_cells = budget.evaluated;
            return Ok(result);
        }
        let index = cells
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| {
                score(a, center, options).total_cmp(&score(b, center, options))
            })
            .map(|(i, _)| i)
            .ok_or(Error::Unresolved)?;
        if cells[index].depth >= options.max_depth {
            return Err(Error::Unresolved);
        }
        if cells.len() == options.max_cells {
            return Err(Error::Limit);
        }
        budget.bytes(
            reservation(source, cells.capacity(), roots, cells.len() + 3)?,
            options.max_bytes,
        )?;
        // Actual degree11 split: at most12 rows, two144-control children,
        // directed12-level de Casteljau arithmetic and transient transpose.
        budget.charge(262_144)?;
        let old = cells.swap_remove(index);
        let children = old.geometry.subdivide(usize::from(old.depth % 2 == 1))?;
        for geometry in children {
            let value = quadrature::enclosure(&geometry, &old.normal, center, budget)?;
            cells.push(Cell {
                geometry,
                normal: old.normal.clone(),
                value,
                depth: old.depth + 1,
            });
        }
    }
}
