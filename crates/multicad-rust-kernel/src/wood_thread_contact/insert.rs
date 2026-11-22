//! Actual authored insert in retained drilled stock, not a headed-screw alias.
use super::*;
use crate::wood_insert::{Classification, Insert};
mod clamp;
mod clearance;
mod clocked;
pub use clamp::{
    InsertClampContact, InsertClampRequirements, InsertClampSources, certify_insert_clamp,
};
pub use clocked::{
    ClockedInsertMating, ClockedLine8Mating, ClockedLine8Placement, ClockedMatingBounds,
    ClockedMatingOptions, certify_clocked_insert_mating, certify_clocked_insert_screw_clearance,
    certify_clocked_line8_mating,
};
pub(super) mod mating;
pub use clearance::{ScrewClearance, ScrewClearanceRequirements, certify_insert_screw_clearance};
pub use mating::{
    InsertMating, MatingExternal, MatingRequirements, MatingScrew, certify_insert_mating,
};
#[derive(Clone, Debug)]
pub struct InsertContact {
    pub engagement_mm: [f64; 2],
    pub tip_clearance_mm: [f64; 2],
    pub thread_witness_mm: [[f64; 2]; 3],
    pub witness_volume_lower_mm3: f64,
}
/// Actual original insert and a small proven nonnegative placement inset.
pub struct InsetInsertContact {
    pub pilot_radius_mm: f64,
    pub pilot_depth_mm: f64,
    pub pilot_mouth_mm: [[f64; 2]; 3],
    pub contact: InsertContact,
    pub source_chart_inset_mm: [f64; 2],
    pub inset_mm: [f64; 2],
}
/// Flush contact retains the existing zero-inset contract.
pub fn certify_insert(
    stock: Stock<'_>,
    insert: &Insert,
    insert_frame: Frame,
    req: Requirements,
    interrupt: &dyn Fn() -> bool,
) -> R<InsertContact> {
    certify_inset_insert(stock, insert, insert_frame, req, 0., interrupt).map(|proof| proof.contact)
}
/// Prove a real small inset without changing stock, source length or phase.
/// The maximum is a caller-lowerable authoring tolerance, not material fuzz.
pub fn certify_inset_insert(
    stock: Stock<'_>,
    insert: &Insert,
    insert_frame: Frame,
    req: Requirements,
    maximum_inset_mm: f64,
    interrupt: &dyn Fn() -> bool,
) -> R<InsetInsertContact> {
    checkpoint(interrupt)?;
    let values = [
        maximum_inset_mm,
        req.minimum_wall_mm,
        req.minimum_floor_mm,
        req.minimum_engagement_mm,
        req.tip_clearance_mm,
        req.linear_tolerance_mm,
        req.minimum_feature_mm,
    ];
    if values
        .into_iter()
        .any(|x| !x.is_finite() || x < 0. || x > 1e7)
        || req.minimum_wall_mm == 0.
        || req.minimum_floor_mm == 0.
        || req.minimum_engagement_mm == 0.
        || req.minimum_feature_mm == 0.
        || !(1e-7..=1e-3).contains(&req.linear_tolerance_mm)
        || maximum_inset_mm > req.linear_tolerance_mm
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
    let sigma = metric.minimum_scale;
    let det = scale(insert_frame)?.determinant;
    let (origin, columns) = relative(stock.frame, insert_frame, interrupt)?;
    let face = boundary.faces.get(stock.pilot.face).ok_or(Error::Invalid)?;
    let (axis, sign) = inward_axis(face.normal)?;
    receiving_halfspace(&stock, axis, sign, interrupt)?;
    let inset = bounded(q(sign)? * (&origin[axis] - q(stock.pilot.center_mm[axis])?))?;
    if columns[2][axis] != sign as i8
        || inset < Q::zero()
        || &inset * q(metric.maximum_scale.hi)? > q(maximum_inset_mm)?
    {
        return Err(Error::Unsupported);
    }
    let inset_interval = if inset.is_zero() {
        ip(0.)
    } else {
        enclosure(&inset)?
    };
    let physical_inset = if inset.is_zero() {
        ip(0.)
    } else {
        inset_interval * I::new(metric.minimum_scale.lo, metric.maximum_scale.hi)
    };
    let radial = [(axis + 1) % 3, (axis + 2) % 3];
    let offset: Q = radial
        .into_iter()
        .map(|i| (&origin[i] - q(stock.pilot.center_mm[i]).unwrap()).abs())
        .sum();
    let eccentricity = enclosure(&offset)?;
    let outer = insert.outer();
    // The sealed Insert proves inner.max + minimumWall <= outer.min exactly.
    // Requiring pilot >= outer.min + center offset therefore keeps the entire
    // actual cavity inside the pilot; it cannot remove intentional wood bite.
    if q(stock.pilot.radius_mm)? < q(outer.minimum_radius_mm())? + &offset
        || q(stock.pilot.radius_mm)? + &offset >= q(outer.maximum_radius_mm())?
    {
        return Err(Error::Unsupported);
    }
    let n = vector(insert_frame.axes[2])?;
    let length = insert.full_length_mm();
    let length_norm = insert_frame.axes[2]
        .map(ip)
        .into_iter()
        .fold(ip(0.), |s, v| s + v.square())
        .sqrt();
    let mut engagement = ip(length) * length_norm;
    let squared = bounded(q(length)?.pow(2) * dot(&n, &n)?)?;
    let gap = bounded(q(stock.pilot.depth_mm)? - &inset - q(length)?)?;
    if squared < q(req.minimum_engagement_mm)?.pow(2)
        || gap < Q::zero()
        || &gap * q(sigma.lo)? < q(req.tip_clearance_mm)?
    {
        return Err(Error::Unsupported);
    }
    engagement.lo = engagement.lo.max(req.minimum_engagement_mm);
    finite(engagement)?;
    let gap = enclosure(&gap)?;
    // Sealed positive length, bounded frame singular values and the exact
    // gap inequality above prove these finite ordered bounds directly.
    let clearance = I::new(
        (gap * sigma).lo.max(req.tip_clearance_mm),
        (gap * metric.maximum_scale).hi,
    );
    let protected =
        (ip(outer.maximum_radius_mm()) + eccentricity + ip(req.minimum_wall_mm) / sigma).hi;
    let floor = (ip(req.minimum_floor_mm) / sigma).hi;
    let mut end = stock.pilot.center_mm;
    let endpoint = ip(end[axis]) + ip(sign) * (ip(stock.pilot.depth_mm) + ip(floor));
    end[axis] = if sign > 0. { endpoint.hi } else { endpoint.lo };
    let mut approach = stock.pilot.center_mm;
    approach[axis] = (ip(approach[axis])
        - ip(sign) * ip(8. * req.linear_tolerance_mm + req.minimum_feature_mm))
    .midpoint();
    checkpoint(interrupt)?;
    let protected_pocket = orthogonal_pockets::identify_pocket(
        stock.material,
        &boundary,
        approach,
        end,
        protected,
        req.linear_tolerance_mm,
        req.minimum_feature_mm,
    )
    .map_err(pocket_error)?;
    same_mouth(protected_pocket, stock.pilot)?;
    let z = length / 2.;
    let phase = outer.bands()[0].phase_mm[1] / 2.;
    let theta = ((ip(z) - ip(phase)) * outer.angular_gain()).midpoint();
    let radius = stock.pilot.radius_mm.midpoint(outer.maximum_radius_mm());
    let center = [
        radius * theta.portable_cos(),
        radius * theta.portable_sin(),
        z,
    ];
    let mut half = ((outer.maximum_radius_mm() - stock.pilot.radius_mm) / 32.)
        .min(phase / 32.)
        .min(length / 32.)
        .min(1e-3);
    for _ in 0..12 {
        checkpoint(interrupt)?;
        let witness = center.map(|c| [(c - half).next_down(), (c + half).next_up()]);
        let classification = insert
            .classify_box(witness, Some(interrupt))
            .map_err(insert_error)?;
        if classification == Classification::Inside {
            let mut local = [ip(0.); 3];
            for i in 0..3 {
                local[i] = enclosure(&origin[i])?;
                for j in 0..3 {
                    local[i] = local[i]
                        + I::new(witness[j][0], witness[j][1]) * ip(f64::from(columns[j][i]));
                }
            }
            let delta = radial.map(|i| local[i] - ip(stock.pilot.center_mm[i]));
            let axial = (local[axis] - ip(stock.pilot.center_mm[axis])) * ip(sign);
            let volume = ip(2. * half) * ip(2. * half) * ip(2. * half) * det;
            if (delta[0].square() + delta[1].square()).lo > ip(stock.pilot.radius_mm).square().hi
                && axial.lo > 0.
                && axial.hi < stock.pilot.depth_mm
                && volume.valid()
                && volume.lo > 0.
            {
                let pilot_mouth_mm =
                    world_point(stock.frame, vector(stock.pilot.center_mm)?, interrupt)?;
                checkpoint(interrupt)?;
                return Ok(InsetInsertContact {
                    pilot_mouth_mm,
                    pilot_radius_mm: stock.pilot.radius_mm,
                    pilot_depth_mm: stock.pilot.depth_mm,
                    contact: InsertContact {
                        engagement_mm: [engagement.lo, engagement.hi],
                        tip_clearance_mm: [clearance.lo, clearance.hi],
                        thread_witness_mm: witness,
                        witness_volume_lower_mm3: volume.lo,
                    },
                    source_chart_inset_mm: [inset_interval.lo.max(0.), inset_interval.hi],
                    inset_mm: [physical_inset.lo.max(0.), physical_inset.hi],
                });
            }
        }
        half /= 2.;
    }
    Err(Error::Unsupported)
}
fn receiving_halfspace(
    stock: &Stock<'_>,
    axis: usize,
    sign: f64,
    interrupt: &dyn Fn() -> bool,
) -> R<()> {
    // A local pilot alone never excludes overhanging actual stock above entry.
    let mouth = q(stock.pilot.center_mm[axis])?;
    for cell in stock.material.cells() {
        checkpoint(interrupt)?;
        if sign > 0. && q(cell.min[axis])? < mouth || sign < 0. && q(cell.max[axis])? > mouth {
            return Err(Error::Unsupported);
        }
    }
    Ok(())
}
fn insert_error(e: crate::wood_insert::Error) -> Error {
    match e {
        crate::wood_insert::Error::Cancelled => Error::Cancelled,
        crate::wood_insert::Error::Limit => Error::Limit,
        crate::wood_insert::Error::Invalid => Error::Invalid,
        crate::wood_insert::Error::Unsupported | crate::wood_insert::Error::Unresolved => {
            Error::Unsupported
        }
    }
}
#[cfg(test)]
mod tests;
