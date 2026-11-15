//! Exact input words preserve numerical and local topology-owner identity.
use super::*;
pub(super) fn words(
    source: &Seats,
    placement: StockClipPlacement,
) -> Result<Vec<u64>, QualificationError> {
    let b = source.boundary();
    let mut count = 4
        + 12
        + 2
        + 1
        + source.stock().cells().len() * 6
        + 1
        + source.stock().quads().len() * 13
        + 2
        + 1
        + source.bores().len() * 9
        + 1
        + source.seats().len() * 4;
    for f in &b.faces {
        count += 1 + 3 + 1 + 1 + f.quad_indices.len();
        for r in &f.rings {
            count += 1 + 3 * r.len();
        }
    }
    if count > 16_384 {
        return Err(QualificationError::Limit);
    }
    let mut out = Vec::new();
    out.try_reserve_exact(count)
        .map_err(|_| QualificationError::Limit)?;
    out.extend([
        1,
        source.work_used as u64,
        source.accounted_peak_bytes as u64,
        source.retained_bytes() as u64,
    ]);
    let mut point = |p: Point| out.extend(p.map(f64::to_bits));
    for row in placement.matrix() {
        point(row)
    }
    point(placement.translation_mm());
    let stock = source.stock();
    out.push(stock.solid_count() as u64);
    out.push(stock.volume_mm3().to_bits());
    out.push(stock.cells().len() as u64);
    for cell in stock.cells() {
        out.extend(cell.min.map(f64::to_bits));
        out.extend(cell.max.map(f64::to_bits));
    }
    out.push(stock.quads().len() as u64);
    for q in stock.quads() {
        out.push(q.component as u64);
        for p in q.corners {
            out.extend(p.map(f64::to_bits));
        }
    }
    let b = source.boundary();
    out.push(b.shell_count as u64);
    out.push(b.faces.len() as u64);
    for f in &b.faces {
        out.push(f.component as u64);
        out.extend(f.normal.map(f64::to_bits));
        out.push(f.rings.len() as u64);
        for r in &f.rings {
            out.push(r.len() as u64);
            for p in r {
                out.extend(p.map(f64::to_bits));
            }
        }
        out.push(f.quad_indices.len() as u64);
        out.extend(f.quad_indices.iter().map(|n| *n as u64));
    }
    out.push(source.bores().len() as u64);
    for b in source.bores() {
        out.extend(b.faces.map(|i| i as u64));
        for p in b.centers_mm {
            out.extend(p.map(f64::to_bits));
        }
        out.push(b.radius_mm.to_bits());
    }
    out.push(source.seats().len() as u64);
    for s in source.seats() {
        out.extend([
            s.bore as u64,
            u64::from(s.end),
            s.radius_mm.to_bits(),
            s.depth_mm.to_bits(),
        ]);
    }
    assert_eq!(out.len(), count, "exact source snapshot word count");
    Ok(out)
}
