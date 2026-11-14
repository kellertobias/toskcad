//! Contact with an actual flat-ended modeled screw and one retained pilot.
//! This certifies geometry and a positive material witness, not installation
//! suitability, preload, pull-out strength or a manufacturer wood profile.
use crate::{
    headed_screw::{Screw, ScrewHead},
    helical_thread::BoxClassification,
    orthogonal::OrthogonalSolid,
    orthogonal_boundary,
    orthogonal_pockets::{self, OrthogonalPocket},
    tapered_thread::interval::I,
};
use num_rational::BigRational as Q;
use num_traits::{Signed, ToPrimitive, Zero};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Error {
    Invalid,
    Unsupported,
    Limit,
    Cancelled,
}
#[derive(Clone, Copy)]
pub struct Frame {
    pub origin_mm: [f64; 3],
    pub axes: [[f64; 3]; 3],
}
#[derive(Clone, Copy)]
pub struct Requirements {
    pub minimum_wall_mm: f64,
    pub minimum_floor_mm: f64,
    pub minimum_engagement_mm: f64,
    pub tip_clearance_mm: f64,
    pub linear_tolerance_mm: f64,
    pub minimum_feature_mm: f64,
}
pub struct Stock<'a> {
    pub material: &'a OrthogonalSolid,
    pub pilot: OrthogonalPocket,
    pub frame: Frame,
}
#[derive(Clone, Debug)]
pub struct Contact {
    pub engagement_mm: [f64; 2],
    /// Normal separation of the actual parallel thread-tip and pilot-floor
    /// supporting planes. Distance to the finite floor is at least this gap.
    pub tip_clearance_mm: [f64; 2],
    /// A complete source-coordinate open box lies in both actual materials.
    pub thread_witness_mm: [[f64; 2]; 3],
    /// This is a contained witness, not total intersected material volume.
    pub witness_volume_lower_mm3: f64,
    pub source_chart_bearing_radial_land_mm: [f64; 2],
}
mod hanger;
mod insert;
pub use hanger::{
    HangerContact, HangerMating, HangerMatingRequirements, HangerPassage, PassageRequirements,
    ThroughStock, WoodSide, certify_hanger_annulus, certify_hanger_nut,
    certify_hanger_through_stock, certify_hanger_wood,
};
mod separation;
pub use insert::{
    InsertContact, InsertMating, MatingExternal, MatingRequirements, MatingScrew, ScrewClearance,
    ScrewClearanceRequirements, certify_insert, certify_insert_mating,
    certify_insert_screw_clearance,
};
pub use separation::{separated_from_rectangular_prism, separated_from_stock};

type R<T> = Result<T, Error>;
fn checkpoint(f: &dyn Fn() -> bool) -> R<()> {
    if f() { Err(Error::Cancelled) } else { Ok(()) }
}
fn q(v: f64) -> R<Q> {
    Q::from_float(v).ok_or(Error::Invalid)
}
fn bounded(v: Q) -> R<Q> {
    if v.numer().bits() + v.denom().bits() > 16384 {
        Err(Error::Limit)
    } else {
        Ok(v)
    }
}
fn dot(a: &[Q; 3], b: &[Q; 3]) -> R<Q> {
    bounded((0..3).map(|i| &a[i] * &b[i]).sum())
}
fn cross(a: &[Q; 3], b: &[Q; 3]) -> R<[Q; 3]> {
    Ok([
        bounded(&a[1] * &b[2] - &a[2] * &b[1])?,
        bounded(&a[2] * &b[0] - &a[0] * &b[2])?,
        bounded(&a[0] * &b[1] - &a[1] * &b[0])?,
    ])
}
fn vector(v: [f64; 3]) -> R<[Q; 3]> {
    Ok([q(v[0])?, q(v[1])?, q(v[2])?])
}
fn enclosure(v: &Q) -> R<I> {
    let x = v.to_f64().ok_or(Error::Unsupported)?;
    let r = I::new(x.next_down(), x.next_up());
    if !r.valid() || q(r.lo)? > *v || q(r.hi)? < *v {
        Err(Error::Unsupported)
    } else {
        Ok(r)
    }
}
fn ip(v: f64) -> I {
    I::point(v)
}
struct Metric {
    minimum_scale: I,
    maximum_scale: I,
    determinant: I,
}
fn exact_orthonormal(frame: Frame) -> R<bool> {
    let a = [
        vector(frame.axes[0])?,
        vector(frame.axes[1])?,
        vector(frame.axes[2])?,
    ];
    for i in 0..3 {
        for j in 0..3 {
            if dot(&a[i], &a[j])? != q(if i == j { 1. } else { 0. })? {
                return Ok(false);
            }
        }
    }
    Ok(true)
}
fn scale(frame: Frame) -> R<Metric> {
    if frame
        .origin_mm
        .into_iter()
        .any(|x| !x.is_finite() || x.abs() > 1e7)
        || frame
            .axes
            .into_iter()
            .flatten()
            .any(|x| !x.is_finite() || x.abs() > 2.)
    {
        return Err(Error::Invalid);
    }
    let mut e = ip(0.);
    for i in 0..3 {
        for j in 0..3 {
            let g = (0..3).fold(ip(0.), |s, k| {
                s + ip(frame.axes[i][k]) * ip(frame.axes[j][k])
            });
            e = e + ip((g - ip(if i == j { 1. } else { 0. })).magnitude());
        }
    }
    let a = frame.axes.map(|a| a.map(ip));
    let d = a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
    if !e.valid() || e.hi > 1e-6 || !d.valid() || d.lo <= 0. {
        return Err(Error::Unsupported);
    }
    // Exact Gram identity plus the proved positive determinant means this
    // literal matrix is an isometry with determinant1. No frame is fitted.
    if exact_orthonormal(frame)? {
        return Ok(Metric {
            minimum_scale: ip(1.),
            maximum_scale: ip(1.),
            determinant: ip(1.),
        });
    }
    Ok(Metric {
        minimum_scale: (ip(1.) - e).sqrt(),
        maximum_scale: (ip(1.) + e).sqrt(),
        determinant: d,
    })
}
/// Exact inverse of the literal stock frame. The returned affine map is never
/// replaced by a fitted/normalized rotation or inferred from nearby points.
fn affine(stock: Frame, screw: Frame, f: &dyn Fn() -> bool) -> R<([Q; 3], [[Q; 3]; 3])> {
    let a = [
        vector(stock.axes[0])?,
        vector(stock.axes[1])?,
        vector(stock.axes[2])?,
    ];
    let c = [
        cross(&a[1], &a[2])?,
        cross(&a[2], &a[0])?,
        cross(&a[0], &a[1])?,
    ];
    let d = dot(&a[0], &c[0])?;
    if d <= Q::zero() {
        return Err(Error::Unsupported);
    }
    let offset: [Q; 3] =
        std::array::from_fn(|i| q(screw.origin_mm[i]).unwrap() - q(stock.origin_mm[i]).unwrap());
    let origin = [
        bounded(dot(&c[0], &offset)? / &d)?,
        bounded(dot(&c[1], &offset)? / &d)?,
        bounded(dot(&c[2], &offset)? / &d)?,
    ];
    let mut columns: [[Q; 3]; 3] = std::array::from_fn(|_| std::array::from_fn(|_| Q::zero()));
    for (j, column) in columns.iter_mut().enumerate() {
        checkpoint(f)?;
        let b = vector(screw.axes[j])?;
        for i in 0..3 {
            column[i] = bounded(dot(&c[i], &b)? / &d)?;
        }
    }
    Ok((origin, columns))
}
fn relative(stock: Frame, screw: Frame, f: &dyn Fn() -> bool) -> R<([Q; 3], [[i8; 3]; 3])> {
    let (origin, actual) = affine(stock, screw, f)?;
    let mut columns = [[0_i8; 3]; 3];
    for j in 0..3 {
        for i in 0..3 {
            columns[j][i] = if actual[j][i] == q(1.)? {
                1
            } else if actual[j][i] == q(-1.)? {
                -1
            } else if actual[j][i].is_zero() {
                0
            } else {
                return Err(Error::Unsupported);
            };
        }
        if columns[j].iter().filter(|&&x| x != 0).count() != 1 {
            return Err(Error::Unsupported);
        }
    }
    for i in 0..3 {
        if columns.iter().filter(|c| c[i] != 0).count() != 1 {
            return Err(Error::Unsupported);
        }
    }
    Ok((origin, columns))
}
fn finite(value: I) -> R<I> {
    if value.valid() {
        Ok(value)
    } else {
        Err(Error::Unsupported)
    }
}
fn positive(value: I) -> R<I> {
    let value = finite(value)?;
    if value.lo > 0. {
        Ok(value)
    } else {
        Err(Error::Unsupported)
    }
}
fn boundary_error(cause: orthogonal_boundary::BoundaryError) -> Error {
    match cause {
        orthogonal_boundary::BoundaryError::Limit => Error::Limit,
        orthogonal_boundary::BoundaryError::Invalid
        | orthogonal_boundary::BoundaryError::Unresolved => Error::Unsupported,
    }
}
fn pocket_error(cause: orthogonal_pockets::PocketError) -> Error {
    match cause {
        orthogonal_pockets::PocketError::BudgetExceeded => Error::Limit,
        orthogonal_pockets::PocketError::InvalidInput => Error::Invalid,
        orthogonal_pockets::PocketError::UnresolvedBoundary => Error::Unsupported,
    }
}
fn boolean_error(cause: crate::orthogonal::BooleanError) -> Error {
    match cause {
        crate::orthogonal::BooleanError::Limit => Error::Limit,
        crate::orthogonal::BooleanError::Invalid | crate::orthogonal::BooleanError::Unresolved => {
            Error::Unsupported
        }
    }
}
fn thread_error(_: crate::helical_thread::ThreadError) -> Error {
    Error::Unsupported
}
fn inward_axis(normal: [f64; 3]) -> R<(usize, f64)> {
    let axis = normal
        .iter()
        .position(|n| n.abs() == 1.)
        .ok_or(Error::Unsupported)?;
    if normal
        .iter()
        .enumerate()
        .any(|(i, &v)| i != axis && v != 0.)
    {
        return Err(Error::Unsupported);
    }
    Ok((axis, -normal[axis]))
}
fn same_mouth(actual: OrthogonalPocket, expected: OrthogonalPocket) -> R<()> {
    if actual.face != expected.face || actual.center_mm != expected.center_mm {
        Err(Error::Unsupported)
    } else {
        Ok(())
    }
}
fn thread_radius(s: &Screw) -> R<f64> {
    let t = s.thread();
    let spec = t.spec();
    let b = t.profile_breaks();
    let mut r = ip(spec.major_radius_mm);
    for (band, k) in [(1, 1), (3, 4)] {
        let mid = b[band].midpoint(b[band + 1]);
        if !(b[band] < mid && mid < b[band + 1]) {
            return Err(Error::Unsupported);
        }
        let slope = t.radius_and_slope_at_phase(mid)[1];
        let v = ip(spec.major_radius_mm) + ip(slope) * (I::new(b[band], b[band + 1]) - ip(b[k]));
        let v = finite(v)?;
        r = I::new(r.lo.min(v.lo), r.hi.max(v.hi));
    }
    let upper = (r.hi + 64. * f64::EPSILON * r.magnitude().max(1.)).next_up();
    Ok(finite(ip(upper))?.hi)
}
fn head_radii(s: &Screw) -> R<(f64, f64)> {
    match s.spec().head {
        ScrewHead::SocketCap { radius_mm, .. } => Ok((radius_mm, radius_mm)),
        ScrewHead::ExternalHex {
            across_flats_mm, ..
        } => {
            let v = crate::headed_screw::hex_vertices(across_flats_mm);
            let mut outer = 0_f64;
            let mut inner = f64::INFINITY;
            for i in 0..6 {
                let a = v[i].map(ip);
                let b = v[(i + 1) % 6].map(ip);
                outer = outer.max((a[0].square() + a[1].square()).sqrt().hi);
                let edge = [b[0] - a[0], b[1] - a[1]];
                let length = (edge[0].square() + edge[1].square()).sqrt();
                let distance = (a[0] * b[1] - a[1] * b[0]) / length;
                inner = inner.min(positive(distance)?.lo);
            }
            Ok((inner, outer))
        }
    }
}
fn gap(interval: [f64; 2], center: I) -> I {
    let lower = (ip(interval[0]) - center)
        .lo
        .max((center - ip(interval[1])).lo)
        .max(0.);
    ip(lower)
}

/// Bounded standalone proof. At most two existing 4M stock-cylinder proofs,
/// 4096 actual cells, fixed-size exact frame arithmetic and 12 box predicates.
pub fn certify(
    stock: Stock<'_>,
    screw: &Screw,
    screw_frame: Frame,
    req: Requirements,
    interrupt: &dyn Fn() -> bool,
) -> R<Contact> {
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
        .any(|x| !x.is_finite() || x < 0. || x > 1e7)
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
    if stock.material.solid_count() != 1 || screw.thread().spec().root_radius_mm.is_some() {
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
    let det = scale(screw_frame)?.determinant;
    let (origin, columns) = relative(stock.frame, screw_frame, interrupt)?;
    let face = boundary.faces.get(stock.pilot.face).ok_or(Error::Invalid)?;
    let (axis, sign) = inward_axis(face.normal)?;
    if columns[2][axis] != sign as i8 {
        return Err(Error::Unsupported);
    }
    let spec = screw.spec();
    let t = screw.thread();
    let plain = q(spec.shank_length_mm)?;
    // A real shared bearing plane, not an interval that happens to contain0.
    if &origin[axis] - q(sign)? * plain != q(stock.pilot.center_mm[axis])? {
        return Err(Error::Unsupported);
    }
    let radial: [usize; 2] = [(axis + 1) % 3, (axis + 2) % 3];
    let offset = radial.map(|i| enclosure(&(&origin[i] - q(stock.pilot.center_mm[i]).unwrap())));
    let offset = [offset[0]?, offset[1]?];
    let eccentricity = offset
        .into_iter()
        .fold(ip(0.), |s, v| s + ip(v.magnitude()));
    // The retained flat branch uses the original R/depth coefficients;
    // also cover the scalar getter's rounded subtraction. Exact dyadic
    // comparisons preserve a genuine shared-radius fit without a waiver.
    let offset_exact: Q = radial
        .into_iter()
        .map(|i| (&origin[i] - q(stock.pilot.center_mm[i]).unwrap()).abs())
        .sum();
    let core = (q(t.spec().major_radius_mm)? - q(t.spec().depth_mm)?).max(q(t.root_radius_mm())?);
    if q(stock.pilot.radius_mm)? < core + &offset_exact
        || spec.shank_length_mm > 0.
            && q(spec.shank_radius_mm)? + &offset_exact > q(stock.pilot.radius_mm)?
    {
        return Err(Error::Unsupported);
    }
    let length = (screw_frame.axes[2]
        .map(ip)
        .into_iter()
        .fold(ip(0.), |s, v| s + v.square()))
    .sqrt();
    let mut engagement = ip(t.spec().length_mm) * length;
    let n = vector(screw_frame.axes[2])?;
    let actual_length_squared = bounded(q(t.spec().length_mm)?.pow(2) * dot(&n, &n)?)?;
    let required_length_squared = q(req.minimum_engagement_mm)?.pow(2);
    let axial_gap =
        bounded(q(stock.pilot.depth_mm)? - q(spec.shank_length_mm)? - q(t.spec().length_mm)?)?;
    if actual_length_squared < required_length_squared
        || axial_gap < Q::zero()
        || &axial_gap * q(sigma.lo)? < q(req.tip_clearance_mm)?
    {
        return Err(Error::Unsupported);
    }
    // Intersect two independently valid bounds: actual norm enclosure and
    // exact squared-length inequality. This preserves true threshold equality.
    engagement.lo = engagement.lo.max(req.minimum_engagement_mm);
    finite(engagement)?;
    let axial_gap = enclosure(&axial_gap)?;
    let clearance_bounds = I::new(
        (axial_gap * sigma).lo.max(req.tip_clearance_mm),
        (axial_gap * metric.maximum_scale).hi,
    );
    let clearance = finite(clearance_bounds)?;
    let radius = thread_radius(screw)?;
    if stock.pilot.radius_mm >= t.spec().major_radius_mm {
        return Err(Error::Unsupported);
    }
    let protected = (ip(radius) + eccentricity + ip(req.minimum_wall_mm) / sigma)
        .hi
        .max(stock.pilot.radius_mm);
    let floor = (ip(req.minimum_floor_mm) / sigma).hi;
    let mut end = stock.pilot.center_mm;
    end[axis] = (ip(end[axis]) + ip(sign) * (ip(stock.pilot.depth_mm) + ip(floor))).lo;
    if sign > 0. {
        end[axis] = (ip(stock.pilot.center_mm[axis]) + ip(stock.pilot.depth_mm) + ip(floor)).hi;
    }
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
    let (head_inner, head_outer) = head_radii(screw)?;
    let head_seat = q(stock.pilot.center_mm[axis])?;
    let head_top = &head_seat - q(sign)? * q(spec.head.height_mm())?;
    for cell in stock.material.cells() {
        checkpoint(interrupt)?;
        let lo = q(cell.min[axis])?;
        let hi = q(cell.max[axis])?;
        if sign > 0. && (lo >= head_seat || hi <= head_top)
            || sign < 0. && (hi <= head_seat || lo >= head_top)
        {
            continue;
        }
        let a = gap(
            [cell.min[radial[0]], cell.max[radial[0]]],
            enclosure(&origin[radial[0]])?,
        );
        let b = gap(
            [cell.min[radial[1]], cell.max[radial[1]]],
            enclosure(&origin[radial[1]])?,
        );
        if (a.square() + b.square()).lo <= ip(head_outer).square().hi {
            return Err(Error::Unsupported);
        }
    }
    let land = [
        (ip(radius.max(spec.shank_radius_mm)) + eccentricity).hi,
        head_inner.min((ip(protected) - eccentricity).lo),
    ];
    if (ip(land[1]) - ip(land[0])).lo * sigma.lo < req.minimum_feature_mm {
        return Err(Error::Unsupported);
    }
    // Proposals retain the original rate and phase. The strict whole-box
    // source predicate below, not this steering arithmetic, admits a witness.
    let z = t.spec().length_mm / 2.;
    let phase = t.spec().crest_width_mm / 4.;
    let theta = ((z - phase) / t.helical_rate()).rem_euclid(std::f64::consts::TAU);
    let r = stock.pilot.radius_mm.midpoint(t.spec().major_radius_mm);
    let center = [r * theta.cos(), r * theta.sin(), z];
    let mut half = ((t.spec().major_radius_mm - stock.pilot.radius_mm) / 32.)
        .min(t.spec().crest_width_mm / 32.)
        .min(t.spec().length_mm / 32.)
        .min(1e-3);
    for _ in 0..12 {
        checkpoint(interrupt)?;
        let witness = center.map(|c| [(c - half).next_down(), (c + half).next_up()]);
        if t.classify_box(witness).map_err(thread_error)? == BoxClassification::Inside {
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
            let v = ip(2. * half) * ip(2. * half) * ip(2. * half) * det;
            if (delta[0].square() + delta[1].square()).lo > ip(stock.pilot.radius_mm).square().hi
                && axial.lo > 0.
                && axial.hi < stock.pilot.depth_mm
                && v.valid()
                && v.lo > 0.
            {
                checkpoint(interrupt)?;
                return Ok(Contact {
                    engagement_mm: [engagement.lo, engagement.hi],
                    tip_clearance_mm: [clearance.lo, clearance.hi],
                    thread_witness_mm: witness,
                    witness_volume_lower_mm3: v.lo,
                    source_chart_bearing_radial_land_mm: land,
                });
            }
        }
        half /= 2.;
    }
    Err(Error::Unsupported)
}
/// An inward cell in an actual exact rectangular prism. Rational supports
/// are not silently rounded into a different occupied stock provider.
pub fn rectangular_obstruction_witness(
    prism: &crate::prism::Prism,
    material_frame: Frame,
    screw: &Screw,
    screw_frame: Frame,
    interrupt: &dyn Fn() -> bool,
) -> R<Option<f64>> {
    checkpoint(interrupt)?;
    if prism.rounded_end().is_some() || prism.inner_loops().len() != 0 {
        return Ok(None);
    }
    let mut upper = [0.; 3];
    for (i, length) in prism.dimensions().into_iter().enumerate() {
        let (n, d) = length.ratio();
        upper[i] = positive(enclosure(&Q::new(n.into(), d.into()))?)?.lo;
    }
    let material = crate::orthogonal::boolean(
        &[crate::orthogonal::Box3 {
            min: [0.; 3],
            max: upper,
        }],
        &[],
        crate::orthogonal::Operation::Union,
        0.,
    )
    .map_err(boolean_error)?
    .ok_or(Error::Unsupported)?;
    obstruction_witness(&material, material_frame, screw, screw_frame, interrupt)
}
/// A positive actual-cell witness in the original thread material. None means
/// unknown; it is NEVER a proof of disjointness and grants no collision waiver.
pub fn obstruction_witness(
    material: &OrthogonalSolid,
    material_frame: Frame,
    screw: &Screw,
    screw_frame: Frame,
    interrupt: &dyn Fn() -> bool,
) -> R<Option<f64>> {
    checkpoint(interrupt)?;
    if material.cells().len() > 4096 {
        return Err(Error::Limit);
    }
    let det = scale(material_frame)?.determinant;
    scale(screw_frame)?;
    let (origin, columns) = affine(screw_frame, material_frame, interrupt)?;
    let origin = [
        enclosure(&origin[0])?,
        enclosure(&origin[1])?,
        enclosure(&origin[2])?,
    ];
    let columns = [
        [
            enclosure(&columns[0][0])?,
            enclosure(&columns[0][1])?,
            enclosure(&columns[0][2])?,
        ],
        [
            enclosure(&columns[1][0])?,
            enclosure(&columns[1][1])?,
            enclosure(&columns[1][2])?,
        ],
        [
            enclosure(&columns[2][0])?,
            enclosure(&columns[2][1])?,
            enclosure(&columns[2][2])?,
        ],
    ];
    for cell in material.cells() {
        checkpoint(interrupt)?;
        let local: [I; 3] = std::array::from_fn(|i| {
            (0..3).fold(origin[i], |s, j| {
                s + columns[j][i] * I::new(cell.min[j], cell.max[j])
            })
        });
        if local.into_iter().all(I::valid)
            && screw
                .thread()
                .classify_box(local.map(|v| [v.lo, v.hi]))
                .map_err(thread_error)?
                == BoxClassification::Inside
        {
            let volume = (0..3).fold(det, |v, i| v * (ip(cell.max[i]) - ip(cell.min[i])));
            let volume = positive(volume)?;
            checkpoint(interrupt)?;
            return Ok(Some(volume.lo));
        }
    }
    checkpoint(interrupt)?;
    Ok(None)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        headed_screw::ScrewSpec,
        helical_thread::{Spec, Thread},
        orthogonal::{self, Box3, Operation},
    };
    fn identity() -> Frame {
        Frame {
            origin_mm: [0.; 3],
            axes: [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
        }
    }
    pub(super) fn fixture(
        left: bool,
        starts: u32,
    ) -> (
        OrthogonalSolid,
        OrthogonalPocket,
        Screw,
        Frame,
        Requirements,
    ) {
        let material = orthogonal::boolean(
            &[Box3 {
                min: [0.; 3],
                max: [20.; 3],
            }],
            &[],
            Operation::Union,
            1e-7,
        )
        .unwrap()
        .unwrap();
        let b = orthogonal_boundary::build(material.quads()).unwrap();
        let pilot = orthogonal_pockets::identify_pocket(
            &material,
            &b,
            [10., 10., 20.001],
            [10., 10., 4.],
            2.,
            1e-7,
            1e-4,
        )
        .unwrap();
        let thread = Thread::new(Spec {
            major_radius_mm: 3.,
            pitch_mm: 2.,
            length_mm: 12.,
            depth_mm: 1.,
            crest_width_mm: 0.4,
            lower_flank_deg: 30.,
            upper_flank_deg: 30.,
            starts,
            left_hand: left,
            root_radius_mm: None,
        })
        .unwrap();
        let screw = Screw::new(
            thread,
            ScrewSpec {
                head: ScrewHead::ExternalHex {
                    across_flats_mm: 10.,
                    height_mm: 4.,
                },
                shank_length_mm: 0.,
                shank_radius_mm: 3.,
            },
        )
        .unwrap();
        (
            material,
            pilot,
            screw,
            Frame {
                origin_mm: [10., 10., 20.],
                axes: [[1., 0., 0.], [0., -1., 0.], [0., 0., -1.]],
            },
            Requirements {
                minimum_wall_mm: 1.,
                minimum_floor_mm: 1.,
                minimum_engagement_mm: 8.,
                tip_clearance_mm: 1.,
                linear_tolerance_mm: 1e-7,
                minimum_feature_mm: 1e-4,
            },
        )
    }
    #[test]
    fn actual_closed_pilot_and_original_flat_threads_prove_positive_bite_and_shared_bearing() {
        for left in [false, true] {
            for starts in [1, 2] {
                let (material, pilot, screw, frame, req) = fixture(left, starts);
                let proof = certify(
                    Stock {
                        material: &material,
                        pilot,
                        frame: identity(),
                    },
                    &screw,
                    frame,
                    req,
                    &|| false,
                )
                .unwrap();
                assert!(proof.engagement_mm[0] <= 12. && proof.engagement_mm[1] >= 12.);
                assert!(proof.witness_volume_lower_mm3 > 0.);
                assert_eq!(
                    screw
                        .thread()
                        .classify_box(proof.thread_witness_mm)
                        .unwrap(),
                    BoxClassification::Inside
                );
                assert!(proof.tip_clearance_mm[0] > 3.99);
                assert!(
                    proof.source_chart_bearing_radial_land_mm[1]
                        > proof.source_chart_bearing_radial_land_mm[0]
                );
            }
        }
    }
    fn pilot(material: &OrthogonalSolid, radius: f64, bottom: bool) -> OrthogonalPocket {
        let boundary = orthogonal_boundary::build(material.quads()).unwrap();
        orthogonal_pockets::identify_pocket(
            material,
            &boundary,
            [10., 10., if bottom { -0.001 } else { 20.001 }],
            [10., 10., if bottom { 16. } else { 4. }],
            radius,
            1e-7,
            1e-4,
        )
        .unwrap()
    }
    fn run(
        material: &OrthogonalSolid,
        pilot: OrthogonalPocket,
        screw: &Screw,
        frame: Frame,
        req: Requirements,
    ) -> R<Contact> {
        certify(
            Stock {
                material,
                pilot,
                frame: identity(),
            },
            screw,
            frame,
            req,
            &|| false,
        )
    }
    #[test]
    fn exact_shared_radius_and_plane_contacts_refuse_real_ulp_penetration_and_invalid_envelopes() {
        let (material, p, screw, frame, req) = fixture(false, 1);
        let original = screw;
        for case in 0..13 {
            let mut p = p;
            let mut f = frame;
            let mut req = req;
            match case {
                0 => p.radius_mm = 2_f64.next_down(),
                1 => p.radius_mm = 3.,
                2 => req.minimum_engagement_mm = 13.,
                3 => req.tip_clearance_mm = 5.,
                4 => req.minimum_wall_mm = 8.,
                5 => req.minimum_floor_mm = 5.,
                6 => f.origin_mm[0] = 10_f64.next_up(),
                7 => f.origin_mm[2] = 20_f64.next_up(),
                8 => f.origin_mm[2] = 20_f64.next_down(),
                9 => f.axes[0][0] = 1. + 1e-12,
                10 => f.axes[2][2] = 1.,
                11 => req.minimum_feature_mm = 1.1,
                _ => p.face = usize::MAX,
            }
            assert!(run(&material, p, &screw, f, req).is_err(), "case{case}");
        }
        assert_eq!(screw, original);
        let mut p = p;
        p.radius_mm = 2.1;
        let mut shifted = frame;
        shifted.origin_mm[0] = 10_f64.next_up();
        assert!(run(&material, p, &screw, shifted, req).is_ok());
        let mut socket = screw.spec();
        socket.head = ScrewHead::SocketCap {
            radius_mm: 5.,
            height_mm: 4.,
            socket_across_flats_mm: 3.,
            socket_depth_mm: 1.,
        };
        assert!(
            run(
                &material,
                p,
                &Screw::new(*screw.thread(), socket).unwrap(),
                frame,
                req
            )
            .is_ok()
        );
        let mut plain = screw.spec();
        plain.shank_length_mm = 1.;
        let mut plain_frame = frame;
        plain_frame.origin_mm[2] = 19.;
        assert!(
            run(
                &material,
                p,
                &Screw::new(*screw.thread(), plain).unwrap(),
                plain_frame,
                req
            )
            .is_err()
        );
        let mut small = screw.spec();
        small.head = ScrewHead::ExternalHex {
            across_flats_mm: 6.0001,
            height_mm: 4.,
        };
        assert!(
            run(
                &material,
                p,
                &Screw::new(*screw.thread(), small).unwrap(),
                frame,
                req
            )
            .is_err()
        );
        let mut short = screw.thread().spec();
        short.length_mm = 2.;
        let mut short_req = req;
        short_req.minimum_engagement_mm = 1.;
        assert!(
            run(
                &material,
                p,
                &Screw::new(Thread::new(short).unwrap(), screw.spec()).unwrap(),
                frame,
                short_req
            )
            .is_ok()
        );
    }
    #[test]
    fn actual_stock_voids_and_head_obstructions_cannot_borrow_the_pilot_exception() {
        let (base, _, screw, frame, req) = fixture(false, 1);
        let cut = orthogonal::boolean(
            base.cells(),
            &[Box3 {
                min: [12.5, 9.8, 8.],
                max: [13., 10.2, 12.],
            }],
            Operation::Subtract,
            1e-7,
        )
        .unwrap()
        .unwrap();
        assert!(run(&cut, pilot(&cut, 2., false), &screw, frame, req).is_err());
        for (near, top) in [(true, 23.), (false, 24.)] {
            let tab = if near {
                Box3 {
                    min: [14.5, 9.5, 20.],
                    max: [15., 10.5, top],
                }
            } else {
                Box3 {
                    min: [0., 0., 20.],
                    max: [2., 2., top],
                }
            };
            let joined = orthogonal::boolean(base.cells(), &[tab], Operation::Union, 1e-7)
                .unwrap()
                .unwrap();
            assert_eq!(
                run(&joined, pilot(&joined, 2., false), &screw, frame, req).is_err(),
                near
            );
        }
        let p = pilot(&base, 2., true);
        let f = Frame {
            origin_mm: [10., 10., 0.],
            ..identity()
        };
        assert!(run(&base, p, &screw, f, req).is_ok());
        let mut rounded = screw.thread().spec();
        rounded.depth_mm = (2. - 0.4) / (2. * 30_f64.to_radians().tan()) - 0.1;
        rounded.root_radius_mm = Some(0.1);
        let rounded = Screw::new(Thread::new(rounded).unwrap(), screw.spec()).unwrap();
        assert!(run(&base, p, &rounded, f, req).is_err());
        let separate = orthogonal::boolean(
            base.cells(),
            &[Box3 {
                min: [30.; 3],
                max: [31.; 3],
            }],
            Operation::Union,
            1e-7,
        )
        .unwrap()
        .unwrap();
        assert!(run(&separate, p, &screw, f, req).is_err());
    }
    #[test]
    fn literal_shared_raw_metric_frames_preserve_geometry_without_normalizing() {
        let (material, mut p, screw, mut frame, req) = fixture(true, 1);
        p.radius_mm = 2.1;
        let stock = Frame {
            origin_mm: [10., 20., 30.],
            axes: [[0.8, 0.6, 0.], [-0.6, 0.8, 0.], [0., 0., 1.]],
        };
        frame.origin_mm = [12., 34., 50.];
        frame.axes = [
            stock.axes[0],
            stock.axes[1].map(|x| -x),
            stock.axes[2].map(|x| -x),
        ];
        let result = certify(
            Stock {
                material: &material,
                pilot: p,
                frame: stock,
            },
            &screw,
            frame,
            req,
            &|| false,
        )
        .unwrap();
        assert!(result.witness_volume_lower_mm3 > 0.);
        let raw = Frame {
            axes: [[1. + 1e-12, 0., 0.], [0., 1., 0.], [0., 0., 1. + 1e-12]],
            ..identity()
        };
        let screw_frame = Frame {
            origin_mm: [10. * (1. + 1e-12), 10., 20. * (1. + 1e-12)],
            axes: [
                raw.axes[0],
                raw.axes[1].map(|x| -x),
                raw.axes[2].map(|x| -x),
            ],
        };
        let active = certify(
            Stock {
                material: &material,
                pilot: p,
                frame: raw,
            },
            &screw,
            screw_frame,
            req,
            &|| false,
        )
        .unwrap();
        assert!(
            active.engagement_mm[0] > 12.,
            "actual raw axis scale remains in the receipt"
        );
        let mut displaced = screw_frame;
        displaced.origin_mm[2] = displaced.origin_mm[2].next_up();
        assert!(
            certify(
                Stock {
                    material: &material,
                    pilot: p,
                    frame: raw
                },
                &screw,
                displaced,
                req,
                &|| false
            )
            .is_err()
        );
        for bad in [
            Frame {
                axes: [[2., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
                ..identity()
            },
            Frame {
                axes: [[-1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
                ..identity()
            },
            Frame {
                origin_mm: [f64::NAN; 3],
                ..identity()
            },
        ] {
            assert!(
                certify(
                    Stock {
                        material: &material,
                        pilot: p,
                        frame: bad
                    },
                    &screw,
                    frame,
                    req,
                    &|| false
                )
                .is_err()
            );
        }
    }
    #[test]
    fn actual_obstruction_witnesses_never_turn_unknown_foreign_material_into_disjointness() {
        let (_, _, screw, frame, _) = fixture(false, 1);
        for (bounds, hit) in [
            (
                Box3 {
                    min: [9., 9., 14.],
                    max: [11., 11., 16.],
                },
                true,
            ),
            (
                Box3 {
                    min: [30.; 3],
                    max: [31.; 3],
                },
                false,
            ),
            (
                Box3 {
                    min: [0.; 3],
                    max: [20.; 3],
                },
                false,
            ),
        ] {
            let material = orthogonal::boolean(&[bounds], &[], Operation::Union, 1e-7)
                .unwrap()
                .unwrap();
            assert_eq!(
                obstruction_witness(&material, identity(), &screw, frame, &|| false)
                    .unwrap()
                    .is_some(),
                hit
            );
        }
        let prism = crate::prism::Prism::rectangular(7, [crate::Rational::new(21, 10).unwrap(); 3])
            .unwrap();
        let placed = Frame {
            origin_mm: [9., 9., 14.],
            ..identity()
        };
        assert!(
            rectangular_obstruction_witness(&prism, placed, &screw, frame, &|| false)
                .unwrap()
                .is_some()
        );
    }
    #[test]
    fn invalid_and_delayed_cancelled_queries_preserve_actual_sources() {
        use std::cell::Cell;
        let (material, p, screw, frame, req) = fixture(false, 1);
        let calls = Cell::new(0);
        certify(
            Stock {
                material: &material,
                pilot: p,
                frame: identity(),
            },
            &screw,
            frame,
            req,
            &|| {
                calls.set(calls.get() + 1);
                false
            },
        )
        .unwrap();
        for target in [1, 4, calls.get()] {
            let n = Cell::new(0);
            assert_eq!(
                certify(
                    Stock {
                        material: &material,
                        pilot: p,
                        frame: identity()
                    },
                    &screw,
                    frame,
                    req,
                    &|| {
                        n.set(n.get() + 1);
                        n.get() >= target
                    }
                )
                .unwrap_err(),
                Error::Cancelled
            );
        }
        for value in [f64::NAN, -1., 0., 1e8] {
            let mut invalid = req;
            invalid.minimum_wall_mm = value;
            assert_eq!(
                run(&material, p, &screw, frame, invalid).unwrap_err(),
                Error::Invalid
            );
        }
        assert!(run(&material, p, &screw, frame, req).is_ok());
    }
    #[test]
    fn real_source_capacity_and_intrinsic_precision_have_bounded_honest_refusals() {
        let (_, p, screw, frame, req) = fixture(false, 1);
        let mut cells = vec![Box3 {
            min: [0.; 3],
            max: [2800., 1., 1.],
        }];
        cells.extend((0..1400).map(|i| Box3 {
            min: [f64::from(i * 2), 1., 0.],
            max: [f64::from(i * 2 + 1), 2., 1.],
        }));
        let large = orthogonal::boolean(&cells, &[], Operation::Union, 1e-7)
            .unwrap()
            .unwrap();
        assert!(large.cells().len() > 4096);
        assert_eq!(
            run(&large, p, &screw, frame, req).unwrap_err(),
            Error::Limit
        );
        assert_eq!(
            obstruction_witness(&large, identity(), &screw, frame, &|| false).unwrap_err(),
            Error::Limit
        );
        let (material, p, screw, frame, req) = fixture(false, 1);
        let mut narrow = screw.thread().spec();
        narrow.crest_width_mm = 1e-15;
        let narrow = Screw::new(Thread::new(narrow).unwrap(), screw.spec()).unwrap();
        assert_eq!(
            run(&material, p, &narrow, frame, req).unwrap_err(),
            Error::Unsupported
        );
        let mut small = screw.thread().spec();
        small.major_radius_mm = 0.001;
        small.depth_mm = 2_f64.powi(-54);
        small.pitch_mm = 0.5;
        small.crest_width_mm = 0.2;
        small.length_mm = 1.;
        let small = Thread::new(small).unwrap();
        let thin = Screw::new(
            small,
            ScrewSpec {
                head: ScrewHead::ExternalHex {
                    across_flats_mm: 0.004,
                    height_mm: 0.001,
                },
                shank_length_mm: 0.,
                shank_radius_mm: 0.001,
            },
        )
        .unwrap();
        assert_eq!(thread_radius(&thin).unwrap_err(), Error::Unsupported);
    }
    #[test]
    fn exact_frame_and_scalar_numerical_primitives_refuse_nonfinite_and_unbounded_requests() {
        assert_eq!(q(f64::NAN).unwrap_err(), Error::Invalid);
        let huge = Q::from_integer(num_bigint::BigInt::from(1_u8) << 20000);
        assert_eq!(bounded(huge).unwrap_err(), Error::Limit);
        let huge = Q::from_integer(num_bigint::BigInt::from(1_u8) << 2000);
        assert_eq!(enclosure(&huge).unwrap_err(), Error::Unsupported);
        assert_eq!(
            enclosure(&q(f64::MAX).unwrap()).unwrap_err(),
            Error::Unsupported
        );
        let third = Q::new(1.into(), 3.into());
        let interval = enclosure(&third).unwrap();
        assert!(q(interval.lo).unwrap() <= third && q(interval.hi).unwrap() >= third);
        let mut singular = identity();
        singular.axes[2] = [0.; 3];
        assert!(relative(singular, identity(), &|| false).is_err());
        let mut duplicate = identity();
        duplicate.axes[1] = duplicate.axes[0];
        assert!(relative(identity(), duplicate, &|| false).is_err());
        let mut zero = identity();
        zero.axes[0] = [0.; 3];
        assert!(relative(identity(), zero, &|| false).is_err());
    }
    #[test]
    fn retained_normal_mouth_and_error_carriers_keep_their_original_failure_contracts() {
        assert!(finite(ip(f64::NAN)).is_err());
        assert!(positive(ip(0.)).is_err());
        assert!(positive(ip(1.)).is_ok());
        assert!(inward_axis([0., 0., 0.5]).is_err());
        assert!(inward_axis([1., 1., 0.]).is_err());
        for (cause, expected) in [
            (
                orthogonal_boundary::BoundaryError::Invalid,
                Error::Unsupported,
            ),
            (
                orthogonal_boundary::BoundaryError::Unresolved,
                Error::Unsupported,
            ),
            (orthogonal_boundary::BoundaryError::Limit, Error::Limit),
        ] {
            assert_eq!(boundary_error(cause), expected)
        }
        for (cause, expected) in [
            (
                orthogonal_pockets::PocketError::InvalidInput,
                Error::Invalid,
            ),
            (
                orthogonal_pockets::PocketError::UnresolvedBoundary,
                Error::Unsupported,
            ),
            (
                orthogonal_pockets::PocketError::BudgetExceeded,
                Error::Limit,
            ),
        ] {
            assert_eq!(pocket_error(cause), expected)
        }
        for (cause, expected) in [
            (orthogonal::BooleanError::Invalid, Error::Unsupported),
            (orthogonal::BooleanError::Unresolved, Error::Unsupported),
            (orthogonal::BooleanError::Limit, Error::Limit),
        ] {
            assert_eq!(boolean_error(cause), expected)
        }
        let (material, p, screw, mut frame, req) = fixture(false, 1);
        let error = screw
            .thread()
            .classify_box([[1., 0.], [0.; 2], [0.; 2]])
            .unwrap_err();
        assert_eq!(thread_error(error), Error::Unsupported);
        let boundary = orthogonal_boundary::build(material.quads()).unwrap();
        let different = orthogonal_pockets::identify_pocket(
            &material,
            &boundary,
            [5., 5., 20.001],
            [5., 5., 4.],
            2.,
            1e-7,
            1e-4,
        )
        .unwrap();
        assert!(same_mouth(different, p).is_err());
        frame.axes = [[0., 1., 0.], [0., 0., -1.], [-1., 0., 0.]];
        assert!(run(&material, p, &screw, frame, req).is_err());
    }
    #[test]
    fn uncertain_witnesses_refuse_and_subnormal_volume_is_never_declared_positive() {
        let (material, mut p, screw, mut frame, req) = fixture(false, 1);
        p.radius_mm = 2.9;
        frame.origin_mm[0] -= 0.5;
        assert!(run(&material, p, &screw, frame, req).is_err());
        // The sealed material constructor itself refuses a cube whose
        // numerical volume would underflow. The scalar proof also refuses
        // the directed underflow enclosure, without forging such a source.
        assert!(
            orthogonal::boolean(
                &[Box3 {
                    min: [0.; 3],
                    max: [1e-110; 3]
                }],
                &[],
                Operation::Union,
                0.
            )
            .is_err()
        );
        assert!(positive(ip(1e-110) * ip(1e-110) * ip(1e-110)).is_err());
        let (_, _, screw, frame, _) = fixture(false, 1);
        let mut prism =
            crate::prism::Prism::rectangular(7, [crate::Rational::new(10, 1).unwrap(); 3]).unwrap();
        prism = prism
            .with_inner_loop(
                9,
                crate::Circle::new(
                    [crate::Rational::new(5, 1).unwrap(); 2],
                    crate::Rational::new(1, 1).unwrap(),
                )
                .unwrap(),
            )
            .unwrap();
        assert!(
            rectangular_obstruction_witness(&prism, identity(), &screw, frame, &|| false)
                .unwrap()
                .is_none()
        );
    }
    #[test]
    fn exact_metric_thresholds_and_raw_clearance_receipts_keep_true_equality() {
        let (material, pilot, screw, frame, mut req) = fixture(false, 1);
        req.minimum_engagement_mm = 12.;
        req.tip_clearance_mm = 4.;
        let proof = run(&material, pilot, &screw, frame, req).unwrap();
        assert!(proof.engagement_mm[0] <= 12. && proof.engagement_mm[1] >= 12.);
        assert!(proof.tip_clearance_mm[0] <= 4. && proof.tip_clearance_mm[1] >= 4.);
        let mut failed = req;
        failed.minimum_engagement_mm = 12_f64.next_up();
        assert_eq!(
            run(&material, pilot, &screw, frame, failed).unwrap_err(),
            Error::Unsupported
        );
        failed = req;
        failed.tip_clearance_mm = 4_f64.next_up();
        assert_eq!(
            run(&material, pilot, &screw, frame, failed).unwrap_err(),
            Error::Unsupported
        );
        let mut end = screw.thread().spec();
        end.length_mm = 16.;
        let end = Screw::new(Thread::new(end).unwrap(), screw.spec()).unwrap();
        req.minimum_engagement_mm = 16.;
        req.tip_clearance_mm = 0.;
        let touching = run(&material, pilot, &end, frame, req).unwrap();
        assert!(touching.tip_clearance_mm[0] <= 0. && touching.tip_clearance_mm[1] >= 0.);
        let mut breakthrough = end.thread().spec();
        breakthrough.length_mm = 16_f64.next_up();
        let breakthrough = Screw::new(Thread::new(breakthrough).unwrap(), end.spec()).unwrap();
        assert_eq!(
            run(&material, pilot, &breakthrough, frame, req).unwrap_err(),
            Error::Unsupported
        );

        let strain = 2_f64.powi(-40);
        let raw = Frame {
            axes: [[1. + strain, 0., 0.], [0., 1., 0.], [0., 0., 1. + strain]],
            ..identity()
        };
        let placed = Frame {
            origin_mm: [10. * (1. + strain), 10., 20. * (1. + strain)],
            axes: [
                raw.axes[0],
                raw.axes[1].map(|x| -x),
                raw.axes[2].map(|x| -x),
            ],
        };
        req.minimum_engagement_mm = 8.;
        req.tip_clearance_mm = 1.;
        let proof = certify(
            Stock {
                material: &material,
                pilot,
                frame: raw,
            },
            &screw,
            placed,
            req,
            &|| false,
        )
        .unwrap();
        let actual = Q::from_integer(4.into()) * q(1. + strain).unwrap();
        assert!(
            q(proof.tip_clearance_mm[0]).unwrap() <= actual
                && q(proof.tip_clearance_mm[1]).unwrap() >= actual
        );
        assert!(
            proof.tip_clearance_mm[1] > 4.,
            "upper bound includes actual stretched cap separation"
        );
    }
}
