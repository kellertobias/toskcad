//! Source-owned clearance through an actual stock bore or annular material.
//! A containing prefix cylinder proves empty-space clearance only; it never
//! becomes manufacturing material or replaces either authored thread source.
use super::*;
use crate::{annulus::Annulus, orthogonal_pockets::OrthogonalThroughBore};
#[cfg(test)]
mod tests;
#[derive(Clone, Copy)]
pub struct ThroughStock<'a> {
    pub material: &'a OrthogonalSolid,
    pub bore: OrthogonalThroughBore,
    pub frame: Frame,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PassageRequirements {
    pub minimum_radial_clearance_mm: f64,
    pub linear_tolerance_mm: f64,
    pub minimum_feature_mm: f64,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum WoodSide {
    BeforeStart,
    AfterEnd,
}
#[derive(Clone, Debug)]
pub struct HangerPassage {
    pub radial_clearance_lower_mm: f64,
    /// Enclosure of the source PREFIX ENVELOPE radius in the bore chart,
    /// not a measured occupied-source support or a material classification.
    pub source_prefix_envelope_radius_mm: [f64; 2],
    /// Normal distance of the bounding wood cylinder from the chosen slab
    /// plane, not an upper bound on the actual wood material's clearance.
    pub wood_envelope_plane_gap_mm: [f64; 2],
    pub wood_side: WoodSide,
}
/// Certify a one-sided enclosure by exact bracketing. In particular a real
/// value just below MAX may round to MAX; an unnecessary next_up would turn
/// its valid upper endpoint into infinity. No conversion accuracy is assumed.
fn directed(value: &Q) -> R<I> {
    let f = value.to_f64().ok_or(Error::Unsupported)?;
    if !f.is_finite() {
        return Err(Error::Unsupported);
    }
    let fq = q(f)?;
    let (lo, hi) = match fq.cmp(value) {
        std::cmp::Ordering::Less => (f, f.next_up()),
        std::cmp::Ordering::Greater => (f.next_down(), f),
        std::cmp::Ordering::Equal => (f, f),
    };
    if !lo.is_finite() || !hi.is_finite() || q(lo)? > *value || q(hi)? < *value {
        Err(Error::Unsupported)
    } else {
        Ok(I::new(lo, hi))
    }
}
/// The caller earned 0<=gap-radius*sqrt(square)<=gap by exact support
/// comparison. Intersecting the raw outward result with that independent
/// enclosure recovers finite endpoints even if a redundant upper overflows.
fn gap_interval(gap: &Q, radius: &Q, square: &Q) -> R<I> {
    let bound = directed(gap)?;
    let raw = bound - directed(radius)? * directed(square)?.sqrt();
    Ok(I::new(raw.lo.max(0.), raw.hi.min(bound.hi)))
}
/// The caller proved 0<=center+radius*sqrt(square)<=cap by exact available
/// radius comparison. This is an independent enclosure, not a fitted source.
fn prefix_interval(center: &Q, radius: &Q, square: &Q, cap: &Q) -> R<I> {
    let upper = directed(cap)?;
    let raw = directed(center)? + directed(radius)? * directed(square)?.sqrt();
    Ok(I::new(raw.lo.max(0.), raw.hi.min(upper.hi)))
}
fn requirements(req: PassageRequirements) -> R<()> {
    if [req.minimum_radial_clearance_mm, req.minimum_feature_mm]
        .into_iter()
        .any(|x| !x.is_finite() || x <= 0. || x > 1e7)
        || !(1e-7..=1e-3).contains(&req.linear_tolerance_mm)
    {
        Err(Error::Invalid)
    } else {
        Ok(())
    }
}
fn margin(value: f64, metric: &Metric) -> R<Q> {
    if metric.minimum_scale.lo == 1. && metric.minimum_scale.hi == 1. {
        q(value)
    } else {
        q(finite(ip(value) / ip(metric.minimum_scale.lo))?.hi)
    }
}
struct Domain {
    axis: usize,
    center: [Q; 3],
    span: [Q; 2],
    radius: Q,
    frame: Frame,
}
fn wood_separation(
    source: &Hanger,
    o: &[Q; 3],
    a: &[[Q; 3]; 3],
    domain: &Domain,
    metric: &Metric,
) -> R<(WoodSide, I)> {
    let k = domain.axis;
    let length = q(source.wood().full_length_mm())?;
    let along = bounded(&a[2][k] * length)?;
    let radial_squared = bounded(&a[0][k] * &a[0][k] + &a[1][k] * &a[1][k])?;
    let radius = q(source.wood().maximum_radius_mm())?;
    let radius_squared = bounded(&radius * &radius * &radial_squared)?;
    let before = bounded(&domain.span[0] - &o[k] - along.clone().max(Q::zero()))?;
    let after = bounded(&o[k] + along.min(Q::zero()) - &domain.span[1])?;
    for (side, gap) in [(WoodSide::BeforeStart, before), (WoodSide::AfterEnd, after)] {
        // An exact squared support inequality admits real plane equality,
        // never an interval merely containing zero or a proximity fit.
        if gap < Q::zero() {
            continue;
        }
        let support_squared = bounded(&gap * &gap)?;
        if support_squared >= radius_squared {
            // Equality is proved in exact source arithmetic, including an
            // oblique cylinder tangent. No rounded interval-zero inference.
            if support_squared == radius_squared {
                return Ok((side, ip(0.)));
            }
            let actual = gap_interval(&gap, &radius, &radial_squared)?;
            let world = finite(actual * I::new(metric.minimum_scale.lo, metric.maximum_scale.hi))?;
            return Ok((side, world));
        }
    }
    Err(Error::Unsupported)
}
fn passage(
    source: &Hanger,
    frame: Frame,
    domain: Domain,
    req: PassageRequirements,
    interrupt: &dyn Fn() -> bool,
) -> R<HangerPassage> {
    scale(frame)?;
    let metric = scale(domain.frame)?;
    let (o, a) = affine(domain.frame, frame, interrupt)?;
    let (wood_side, wood_gap) = wood_separation(source, &o, &a, &domain, &metric)?;
    let radial = [(domain.axis + 1) % 3, (domain.axis + 2) % 3];
    let zero = Q::zero();
    let gram = |u: usize, v: usize| {
        bounded(
            radial
                .iter()
                .fold(Q::zero(), |s, &j| s + &a[u][j] * &a[v][j]),
        )
    };
    let g00 = gram(0, 0)?;
    let g11 = gram(1, 1)?;
    let g01 = gram(0, 1)?.abs();
    // Gershgorin bounds the exact two-column radial operator. This is an
    // enclosure of the ACTUAL affine image, without normalizing its axes.
    let lambda_squared = bounded((g00 + &g01).max(g11 + g01))?;
    let z0 = bounded(-q(source.spec().plain_length_mm)? - q(source.metric().full_length_mm())?)?;
    let mut center = Q::zero();
    for j in radial {
        checkpoint(interrupt)?;
        let at_zero = bounded(&o[j] - &domain.center[j])?;
        let at_tip = bounded(&at_zero + &a[2][j] * &z0)?;
        center = bounded(center + at_zero.abs().max(at_tip.abs()))?;
    }
    // Strict Hanger support dominance proves BOTH actual negative sections
    // lie inside this real plain-radius cylinder over their full finite z.
    let radius = q(source.spec().plain_radius_mm)?;
    let available =
        bounded(&domain.radius - &center - margin(req.minimum_radial_clearance_mm, &metric)?)?;
    if available < zero
        || bounded(&available * &available)? < bounded(&radius * &radius * &lambda_squared)?
    {
        return Err(Error::Unsupported);
    }
    let envelope = prefix_interval(&center, &radius, &lambda_squared, &domain.radius)?;
    let chart_gap = finite(directed(&domain.radius)? - envelope)?;
    let physical = finite(chart_gap * I::new(metric.minimum_scale.lo, metric.maximum_scale.hi))?;
    checkpoint(interrupt)?;
    Ok(HangerPassage {
        radial_clearance_lower_mm: physical.lo.max(req.minimum_radial_clearance_mm),
        source_prefix_envelope_radius_mm: [envelope.lo, envelope.hi],
        wood_envelope_plane_gap_mm: [wood_gap.lo, wood_gap.hi],
        wood_side,
    })
}
/// Bind this theorem to the ORIGINAL stock and retained cut domain. The Native
/// caller must preserve that ownership; a recipe radius alone is insufficient.
pub fn certify_hanger_through_stock(
    stock: ThroughStock<'_>,
    source: &Hanger,
    frame: Frame,
    minimum_wall_mm: f64,
    req: PassageRequirements,
    interrupt: &dyn Fn() -> bool,
) -> R<HangerPassage> {
    checkpoint(interrupt)?;
    requirements(req)?;
    if !minimum_wall_mm.is_finite() || minimum_wall_mm <= 0. || minimum_wall_mm > 1e7 {
        return Err(Error::Invalid);
    }
    if stock.material.cells().len() > 4096 {
        return Err(Error::Limit);
    }
    if stock.material.solid_count() != 1 {
        return Err(Error::Unsupported);
    }
    let metric = scale(stock.frame)?;
    let boundary = orthogonal_boundary::build(stock.material.quads()).map_err(boundary_error)?;
    orthogonal_pockets::validate_with_through_bores(
        stock.material,
        &boundary,
        &[],
        &[stock.bore],
        req.linear_tolerance_mm,
        req.minimum_feature_mm,
    )
    .map_err(pocket_error)?;
    let wall = margin(minimum_wall_mm, &metric)?;
    let mut protected = stock.bore;
    protected.radius_mm = directed(&bounded(q(protected.radius_mm)? + wall)?)?.hi;
    orthogonal_pockets::validate_with_through_bores(
        stock.material,
        &boundary,
        &[],
        &[protected],
        req.linear_tolerance_mm,
        req.minimum_feature_mm,
    )
    .map_err(pocket_error)?;
    let face = boundary
        .faces
        .get(stock.bore.faces[0])
        .ok_or(Error::Invalid)?;
    let (axis, _) = inward_axis(face.normal)?;
    let [p0, p1] = stock.bore.centers_mm;
    let a = q(p0[axis])?;
    let b = q(p1[axis])?;
    let span = [a.clone().min(b.clone()), a.max(b)];
    for cell in stock.material.cells() {
        checkpoint(interrupt)?;
        if q(cell.min[axis])? < span[0] || q(cell.max[axis])? > span[1] {
            return Err(Error::Unsupported);
        }
    }
    passage(
        source,
        frame,
        Domain {
            axis,
            center: vector(p0)?,
            span,
            radius: q(stock.bore.radius_mm)?,
            frame: stock.frame,
        },
        req,
        interrupt,
    )
}
/// An actual annular material occupies only its finite axial slab outside its
/// retained inner cylinder. No annulus is inferred from a label or a bbox.
pub fn certify_hanger_annulus(
    ring: &Annulus,
    ring_frame: Frame,
    source: &Hanger,
    frame: Frame,
    req: PassageRequirements,
    interrupt: &dyn Fn() -> bool,
) -> R<HangerPassage> {
    checkpoint(interrupt)?;
    requirements(req)?;
    if bounded(q(ring.outer_radius_mm())? - q(ring.inner_radius_mm())?)?
        < q(req.minimum_feature_mm)?
        || q(ring.height_mm())? < q(req.minimum_feature_mm)?
    {
        return Err(Error::Unsupported);
    }
    passage(
        source,
        frame,
        Domain {
            axis: 2,
            center: [Q::zero(), Q::zero(), Q::zero()],
            span: [Q::zero(), q(ring.height_mm())?],
            radius: q(ring.inner_radius_mm())?,
            frame: ring_frame,
        },
        req,
        interrupt,
    )
}
