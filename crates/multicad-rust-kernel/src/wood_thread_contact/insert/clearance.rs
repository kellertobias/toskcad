//! Complete actual screw versus the retained blind pilot, independently of the
//! insert's shorter finite cavity. No stock collision is waived by mating.
use super::*;
use mating::MatingExternal;
#[derive(Clone, Copy)]
pub struct ScrewClearanceRequirements {
    pub minimum_tip_clearance_mm: f64,
    pub linear_tolerance_mm: f64,
    pub minimum_feature_mm: f64,
}
#[derive(Clone, Debug)]
pub struct ScrewClearance {
    pub tip_clearance_mm: [f64; 2],
    pub prefix_entry_clearance_mm: [f64; 2],
    pub minimum_radial_clearance_mm: f64,
}
pub fn certify_insert_screw_clearance(
    stock: Stock<'_>,
    external: MatingExternal<'_>,
    req: ScrewClearanceRequirements,
    interrupt: &dyn Fn() -> bool,
) -> R<ScrewClearance> {
    checkpoint(interrupt)?;
    if [
        req.minimum_tip_clearance_mm,
        req.linear_tolerance_mm,
        req.minimum_feature_mm,
    ]
    .into_iter()
    .any(|x| !x.is_finite() || x < 0. || x > 1e7)
        || !(1e-7..=1e-3).contains(&req.linear_tolerance_mm)
        || req.minimum_feature_mm == 0.
    {
        return Err(Error::Invalid);
    }
    if stock.material.cells().len() > 4096 {
        return Err(Error::Limit);
    }
    if stock.material.solid_count() != 1 {
        return Err(Error::Unsupported);
    }
    let boundary = orthogonal_boundary::build(stock.material.quads()).map_err(boundary_error)?;
    orthogonal_pockets::validate(
        stock.material,
        &boundary,
        &[stock.pilot],
        req.linear_tolerance_mm,
        req.minimum_feature_mm,
    )
    .map_err(pocket_error)?;
    let metric = scale(stock.frame)?;
    scale(external.frame)?;
    let (origin, columns) = relative(stock.frame, external.frame, interrupt)?;
    let face = boundary.faces.get(stock.pilot.face).ok_or(Error::Invalid)?;
    let (axis, sign) = inward_axis(face.normal)?;
    receiving_halfspace(&stock, axis, sign, interrupt)?;
    if columns[2][axis] != sign as i8 {
        return Err(Error::Unsupported);
    }
    let mating::Parts {
        thread,
        visible: [_, end],
        prefix,
    } = external.source.parts()?;
    let station = bounded(q(sign)? * (&origin[axis] - q(stock.pilot.center_mm[axis])?))?;
    let prefix_gap = bounded(-&station - prefix)?;
    if prefix_gap < Q::zero() {
        return Err(Error::Unsupported);
    }
    let tip = bounded(&station + end)?;
    let gap = bounded(q(stock.pilot.depth_mm)? - tip)?;
    // Exact binary-rational datum/end supports establish the real floor gap;
    // the retained stock-frame metric converts its local separation to world.
    if gap < Q::zero() || &gap * q(metric.minimum_scale.lo)? < q(req.minimum_tip_clearance_mm)? {
        return Err(Error::Unsupported);
    }
    let radial = [(axis + 1) % 3, (axis + 2) % 3];
    let offset: Q = radial
        .into_iter()
        .map(|i| (&origin[i] - q(stock.pilot.center_mm[i]).unwrap()).abs())
        .sum();
    let radius = thread
        .radius_at_phase_interval(I::new(0., thread.spec().pitch_mm))
        .ok_or(Error::Unsupported)?;
    let remaining = bounded(q(stock.pilot.radius_mm)? - q(radius.hi)? - offset)?;
    if remaining <= Q::zero()
        || &remaining * q(metric.minimum_scale.lo)? < q(req.minimum_feature_mm)?
    {
        return Err(Error::Unsupported);
    }
    let mut tip_clearance = mating::physical_gap(&gap, &metric)?;
    // Intersect with the independently earned exact threshold inequality.
    tip_clearance[0] = tip_clearance[0].max(req.minimum_tip_clearance_mm);
    checkpoint(interrupt)?;
    Ok(ScrewClearance {
        tip_clearance_mm: tip_clearance,
        prefix_entry_clearance_mm: mating::physical_gap(&prefix_gap, &metric)?,
        minimum_radial_clearance_mm: (enclosure(&remaining)? * metric.minimum_scale)
            .lo
            .max(req.minimum_feature_mm),
    })
}
#[cfg(test)]
mod tests;
