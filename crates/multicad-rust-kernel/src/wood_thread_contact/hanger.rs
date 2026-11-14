//! Actual authored Hanger wood bite and retained-nut mating. A wood-thread
//! portion's intentional material overlap is distinct from all other hardware
//! and stock collisions; neither a recipe nor an envelope is material proof.
use super::*;
use crate::authored_hanger_bolt::{Classification, Hanger};
mod mating;
pub use mating::{HangerMating, HangerMatingRequirements, certify_hanger_nut};
#[cfg(test)]
mod tests;
#[derive(Clone, Debug)]
pub struct HangerContact {
    pub engagement_mm: [f64; 2],
    /// Distance of the actual parallel wood-tip and pilot-floor planes.
    pub tip_clearance_mm: [f64; 2],
    pub source_axial_offset_mm: [f64; 2],
    pub source_eccentricity_mm: [f64; 2],
    pub thread_witness_mm: [[f64; 2]; 3],
    pub witness_volume_lower_mm3: f64,
}
fn exact_interval(v: &Q) -> R<I> {
    let f = v.to_f64().ok_or(Error::Unsupported)?;
    if q(f)? == *v { Ok(ip(f)) } else { enclosure(v) }
}
fn source_error(e: crate::authored_hanger_bolt::Error) -> Error {
    match e {
        crate::authored_hanger_bolt::Error::Invalid => Error::Invalid,
        crate::authored_hanger_bolt::Error::Limit => Error::Limit,
        crate::authored_hanger_bolt::Error::Cancelled => Error::Cancelled,
        crate::authored_hanger_bolt::Error::Unsupported
        | crate::authored_hanger_bolt::Error::Unresolved => Error::Unsupported,
    }
}
/// One actual blind pilot and the original full wood segment. Alignment
/// tolerance is an intent gate only: exact actual offsets independently enter
/// every physical stock, prefix, engagement, floor and witness proof.
pub fn certify_hanger_wood(
    stock: Stock<'_>,
    source: &Hanger,
    frame: Frame,
    req: Requirements,
    interrupt: &dyn Fn() -> bool,
) -> R<HangerContact> {
    checkpoint(interrupt)?;
    let values = [
        req.minimum_wall_mm,
        req.minimum_floor_mm,
        req.minimum_engagement_mm,
        req.tip_clearance_mm,
        req.linear_tolerance_mm,
        req.minimum_feature_mm,
    ];
    if values
        .into_iter()
        .any(|v| !v.is_finite() || v < 0. || v > 1e7)
        || req.minimum_wall_mm == 0.
        || req.minimum_floor_mm == 0.
        || req.minimum_engagement_mm == 0.
        || req.minimum_feature_mm == 0.
        || !(1e-7..=1e-3).contains(&req.linear_tolerance_mm)
    {
        return Err(Error::Invalid);
    }
    if stock.material.cells().len() > 4096 {
        return Err(Error::Limit);
    }
    if stock.material.solid_count() != 1 {
        return Err(Error::Unsupported);
    }
    let metric = scale(stock.frame)?;
    let placed = scale(frame)?;
    let boundary = orthogonal_boundary::build(stock.material.quads()).map_err(boundary_error)?;
    orthogonal_pockets::validate(
        stock.material,
        &boundary,
        &[stock.pilot],
        req.linear_tolerance_mm,
        req.minimum_feature_mm,
    )
    .map_err(pocket_error)?;
    let (origin, columns) = relative(stock.frame, frame, interrupt)?;
    let face = boundary.faces.get(stock.pilot.face).ok_or(Error::Invalid)?;
    let (axis, sign) = inward_axis(face.normal)?;
    // All stock is on the actual receiving side, not only a local pilot.
    let mouth = q(stock.pilot.center_mm[axis])?;
    for cell in stock.material.cells() {
        checkpoint(interrupt)?;
        if sign > 0. && q(cell.min[axis])? < mouth || sign < 0. && q(cell.max[axis])? > mouth {
            return Err(Error::Unsupported);
        }
    }
    if columns[2][axis] != sign as i8 {
        return Err(Error::Unsupported);
    }
    let axial = bounded(q(sign)? * (&origin[axis] - &mouth))?;
    let radial = [(axis + 1) % 3, (axis + 2) % 3];
    let radial_offset: Q = radial
        .into_iter()
        .map(|i| (&origin[i] - q(stock.pilot.center_mm[i]).unwrap()).abs())
        .sum();
    let eccentric = bounded(radial_offset)?;
    let e = exact_interval(&eccentric)?;
    let a = exact_interval(&axial)?;
    // Exact physical origin-distance intent. This never proves equality or
    // substitutes for the independent physical containment checks below.
    let delta = finite_origin_delta(
        frame.origin_mm,
        stock.frame.origin_mm,
        stock.frame.axes,
        stock.pilot.center_mm,
    );
    if dot(&delta, &delta)? > q(req.linear_tolerance_mm)?.pow(2) {
        return Err(Error::Unsupported);
    }
    let wood = source.wood();
    let length = q(wood.full_length_mm())?;
    // The sealed Hanger guarantees metric.max < plainR < wood.min. Therefore
    // every non-wood portion is bounded by this real plain cylinder. If actual
    // positive axial offset brings any prefix into stock, it fits the pilot.
    if q(source.spec().plain_radius_mm)? + &eccentric > q(stock.pilot.radius_mm)?
        || q(stock.pilot.radius_mm)? + &eccentric >= q(wood.maximum_radius_mm())?
    {
        return Err(Error::Unsupported);
    }
    let end = bounded(&axial + &length)?;
    let begin = axial.clone().max(Q::zero());
    let overlap = bounded(&end - &begin)?;
    let gap = bounded(q(stock.pilot.depth_mm)? - &end)?;
    let n = vector(frame.axes[2])?;
    if overlap <= Q::zero()
        || bounded(overlap.pow(2) * dot(&n, &n)?)? < q(req.minimum_engagement_mm)?.pow(2)
        || gap < Q::zero()
        || bounded(&gap * q(metric.minimum_scale.lo)?)? < q(req.tip_clearance_mm)?
    {
        return Err(Error::Unsupported);
    }
    let norm = frame.axes[2]
        .map(ip)
        .into_iter()
        .fold(ip(0.), |s, x| s + x.square())
        .sqrt();
    let mut engagement = finite(exact_interval(&overlap)? * norm)?;
    engagement.lo = engagement.lo.max(req.minimum_engagement_mm);
    // Singular bounds enclose the actual normal-plane distance as well as the
    // axial-column distance; the receipt promises the former, not a fittedaxis.
    let mut clearance =
        finite(exact_interval(&gap)? * I::new(metric.minimum_scale.lo, metric.maximum_scale.hi))?;
    clearance.lo = clearance.lo.max(req.tip_clearance_mm);
    let protected =
        (ip(wood.maximum_radius_mm()) + e + ip(req.minimum_wall_mm) / metric.minimum_scale).hi;
    let floor = (ip(req.minimum_floor_mm) / metric.minimum_scale).hi;
    let mut endpoint = stock.pilot.center_mm;
    let z = ip(endpoint[axis]) + ip(sign) * (ip(stock.pilot.depth_mm) + ip(floor));
    endpoint[axis] = if sign > 0. { z.hi } else { z.lo };
    let mut approach = stock.pilot.center_mm;
    approach[axis] = (ip(approach[axis])
        - ip(sign) * ip(8. * req.linear_tolerance_mm + req.minimum_feature_mm))
    .midpoint();
    checkpoint(interrupt)?;
    let complete = orthogonal_pockets::identify_pocket(
        stock.material,
        &boundary,
        approach,
        endpoint,
        protected,
        req.linear_tolerance_mm,
        req.minimum_feature_mm,
    )
    .map_err(pocket_error)?;
    same_mouth(complete, stock.pilot)?;
    // A proposal only. Publication needs the actual whole Hanger predicate and
    // actual mapped stock/pilot inequalities for the entire open box.
    let source_begin = (-axial.clone()).max(Q::zero());
    let center_z = (exact_interval(&source_begin)? + ip(wood.full_length_mm())).midpoint() * 0.5;
    let phase = wood.bands()[0].phase_mm[1] * 0.5;
    let theta = ((ip(center_z) - ip(phase)) * wood.angular_gain()).midpoint();
    let radius = (ip(stock.pilot.radius_mm) + e + ip(wood.maximum_radius_mm())).midpoint() * 0.5;
    let center = [radius * theta.cos(), radius * theta.sin(), center_z];
    let mut half = ((ip(wood.maximum_radius_mm()) - ip(stock.pilot.radius_mm) - e).lo / 32.)
        .min(phase / 32.)
        .min(overlap.to_f64().ok_or(Error::Unsupported)? / 32.)
        .min(1e-3);
    half = positive(ip(half))?.lo;
    for _ in 0..12 {
        checkpoint(interrupt)?;
        let witness = center.map(|v| {
            let b = ip(v) + I::new(-half, half);
            [b.lo, b.hi]
        });
        let placement = WitnessPlacement {
            origin: &origin,
            columns: &columns,
            axis,
            sign,
            determinant: placed.determinant,
        };
        if let Some(volume) = material_witness(source, stock.pilot, witness, placement, interrupt)?
        {
            checkpoint(interrupt)?;
            return Ok(HangerContact {
                engagement_mm: [engagement.lo, engagement.hi],
                tip_clearance_mm: [clearance.lo, clearance.hi],
                source_axial_offset_mm: [a.lo, a.hi],
                source_eccentricity_mm: [e.lo, e.hi],
                thread_witness_mm: witness,
                witness_volume_lower_mm3: volume,
            });
        }
        half *= 0.5;
    }
    Err(Error::Unsupported)
}

struct WitnessPlacement<'a> {
    origin: &'a [Q; 3],
    columns: &'a [[i8; 3]; 3],
    axis: usize,
    sign: f64,
    determinant: I,
}
fn material_witness(
    source: &Hanger,
    pilot: OrthogonalPocket,
    witness: [[f64; 2]; 3],
    map: WitnessPlacement<'_>,
    interrupt: &dyn Fn() -> bool,
) -> R<Option<f64>> {
    if source
        .classify_box(witness, Some(interrupt))
        .map_err(source_error)?
        == Classification::Inside
    {
        witness_in_stock(pilot, witness, map)
    } else {
        Ok(None)
    }
}
/// Source membership alone is insufficient: the whole box must lie in actual
/// stock material outside the finite pilot, with positive placed volume.
fn witness_in_stock(
    pilot: OrthogonalPocket,
    witness: [[f64; 2]; 3],
    map: WitnessPlacement<'_>,
) -> R<Option<f64>> {
    let mut delta = [ip(0.); 3];
    for (i, target) in delta.iter_mut().enumerate() {
        *target = exact_interval(&bounded(&map.origin[i] - q(pilot.center_mm[i])?)?)?;
        for (j, b) in witness.iter().enumerate() {
            *target = *target + I::new(b[0], b[1]) * ip(f64::from(map.columns[j][i]));
        }
    }
    let radial = [(map.axis + 1) % 3, (map.axis + 2) % 3];
    let depth = delta[map.axis] * ip(map.sign);
    let volume = witness
        .into_iter()
        .fold(ip(1.), |v, b| v * (ip(b[1]) - ip(b[0])))
        * map.determinant;
    if (delta[radial[0]].square() + delta[radial[1]].square()).lo > ip(pilot.radius_mm).square().hi
        && depth.lo > 0.
        && depth.hi < pilot.depth_mm
        && volume.valid()
        && volume.lo > 0.
    {
        Ok(Some(volume.lo))
    } else {
        Ok(None)
    }
}

/// Three products of finite binary64 pairs and two finite origins are dyadic
/// with denominator<=2^2148 and magnitude<2^2052. Their three squared sums have
/// combined numerator/denominator length<12710bits, below the existing16384
/// guard, even before the stronger public frame/dimension admission bounds.
/// No cofactor-ratio or general rational bound is removed by this specialization.
fn finite_origin_delta(
    point: [f64; 3],
    origin: [f64; 3],
    axes: [[f64; 3]; 3],
    center: [f64; 3],
) -> [Q; 3] {
    // Public frame scale and actual pocket validation already proved every
    // field finite. Revalidating here would add an impossible failure path.
    let mut delta = point.map(|v| Q::from_float(v).unwrap());
    let origin = origin.map(|v| Q::from_float(v).unwrap());
    let center = center.map(|v| Q::from_float(v).unwrap());
    let axes = axes.map(|a| a.map(|v| Q::from_float(v).unwrap()));
    for (i, d) in delta.iter_mut().enumerate() {
        *d = &*d - &origin[i] - (0..3).map(|j| &axes[j][i] * &center[j]).sum::<Q>();
    }
    delta
}
