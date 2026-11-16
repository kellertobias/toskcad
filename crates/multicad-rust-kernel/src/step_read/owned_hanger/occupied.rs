//! Decoded-support-derived healed material clipped by decoded planes.
//! Original finite coedges remain bounded correspondence records; their raw
//! occupied set and mass are not asserted. No authoring recipe is used.
use super::*;
mod query;
mod volume;
use super::{
    bernstein::{self as poly, Algebra, Polynomial},
    polar::{PolarStrip, PolarSupports},
};
use num_bigint::{BigInt, Sign};
use num_rational::BigRational;
use num_traits::{ToPrimitive, Zero};
pub use query::{PointLocation, PointQuery, QueryOptions};
pub use volume::{HealedVolume, VolumeOptions};
#[derive(Clone, Copy)]
pub struct HealingOptions {
    pub maximum_work: usize,
    pub maximum_bytes: usize,
    pub maximum_sewing_error_mm: f64,
}
impl Default for HealingOptions {
    fn default() -> Self {
        Self {
            maximum_work: 64_000_000,
            maximum_bytes: 64 * 1024 * 1024,
            maximum_sewing_error_mm: 1e-7,
        }
    }
}
pub struct HealedMaterial<'p, 'g, 'd, 'l, 'i> {
    polar: &'p PolarSupports<'g, 'd, 'i>,
    lifts: &'l SurfaceLifts<'d, 'i>,
    norm: I,
    jacobian: f64,
    sewing: f64,
    volume: I,
    radii: [I; 3],
    work: usize,
    bytes: usize,
}
fn p(x: f64) -> I {
    I::point(x)
}
fn physical(n: &BigInt, d: &BigInt, b: &mut Budget<'_>) -> Result<I, Error> {
    let words = (n.bits().max(d.bits()) as usize).div_ceil(64);
    b.charge(128 * (1 + words * words))?;
    let q = BigRational::new(n.clone(), d.clone());
    let v = q.to_f64().ok_or(Error::Unsupported)?;
    let cache = BigRational::from_float(v).ok_or(Error::Unsupported)?;
    Ok(match q.cmp(&cache) {
        std::cmp::Ordering::Equal => p(v),
        std::cmp::Ordering::Less => I::new(v.next_down(), v),
        std::cmp::Ordering::Greater => I::new(v, v.next_up()),
    })
}
fn integer(v: &PoleInteger, b: &mut Budget<'_>) -> Result<BigInt, Error> {
    b.charge(v.numerator_magnitude_le().len() + 1)?;
    Ok(BigInt::from_bytes_le(
        match v.signum() {
            -1 => Sign::Minus,
            0 => Sign::NoSign,
            _ => Sign::Plus,
        },
        v.numerator_magnitude_le(),
    ))
}
fn angular(
    polar: &PolarSupports<'_, '_, '_>,
    algebra: &Algebra,
    b: &mut Budget<'_>,
) -> Result<(I, f64, f64), Error> {
    let mut norm = I::new(f64::INFINITY, 0.);
    let mut minimum = f64::INFINITY;
    let mut derivative = 0_f64;
    for piece in &polar.angular {
        let x = piece.bounds.iter().map(|v| v[0]).collect::<Polynomial>();
        let y = piece.bounds.iter().map(|v| v[1]).collect::<Polynomial>();
        let mut squared = algebra.multiply(&x, &x, b)?;
        let yy = algebra.multiply(&y, &y, b)?;
        poly::add(&mut squared, &yy, b)?;
        let lo = poly::lower(&squared);
        let hi = squared.iter().map(|v| v.hi).fold(0., f64::max);
        if lo <= 0. || !hi.is_finite() {
            return Err(Error::Unsupported);
        }
        norm.lo = norm.lo.min(lo);
        norm.hi = norm.hi.max(hi);
        let dx = poly::derivative(&x, b)?;
        let dy = poly::derivative(&y, b)?;
        let a = algebra.multiply(&x, &dy, b)?;
        let z = algebra.multiply(&y, &dx, b)?;
        let j = poly::subtract(&a, &z, b)?;
        let lo = poly::lower(&j);
        if lo <= 0. {
            return Err(Error::Unsupported);
        }
        minimum = minimum.min(lo);
        let width = piece.knots[1].clone() - piece.knots[0].clone();
        let width = physical(width.numer(), width.denom(), b)?;
        derivative =
            derivative.max(((p(poly::magnitude(&dx)) + p(poly::magnitude(&dy))) / width).hi);
        let quarter = (&piece.knots[0] * BigRational::from_integer(4.into()))
            .to_integer()
            .to_usize()
            .ok_or(Error::Unsupported)?;
        if quarter >= 4 {
            return Err(Error::Unsupported);
        }
        for q in &piece.points {
            if (matches!(quarter, 0 | 3) && q[0] < BigInt::zero())
                || (matches!(quarter, 1 | 2) && q[0] > BigInt::zero())
                || (quarter < 2 && q[1] < BigInt::zero())
                || (quarter >= 2 && q[1] > BigInt::zero())
            {
                return Err(Error::Unsupported);
            }
        }
    }
    if !norm.valid() || !minimum.is_finite() {
        return Err(Error::Unsupported);
    }
    Ok((norm.sqrt(), minimum, derivative))
}
fn strip_speed(strip: &PolarStrip, lead: I, norm: I, derivative: f64) -> Result<f64, Error> {
    let r = strip
        .radius
        .iter()
        .map(|r| r.magnitude())
        .fold(0., f64::max);
    let q = strip.phase.iter().map(|q| q.magnitude()).fold(0., f64::max);
    let minw = strip
        .weights
        .iter()
        .map(|w| w.lo)
        .fold(f64::INFINITY, f64::min);
    let maxw = strip.weights.iter().map(|w| w.hi).fold(0., f64::max);
    if minw <= 0. || lead.lo <= 0. {
        return Err(Error::Unsupported);
    }
    let factor =
        p(4. * (strip.radius.len() - 1) as f64) * (p(maxw) / p(minw)) * (p(maxw) / p(minw));
    let speed = factor * p(r) * p(norm.hi) + p(r) * factor * p(q) / lead * p(derivative);
    if !speed.valid() {
        return Err(Error::Unsupported);
    }
    Ok(speed.hi)
}
fn side_u(d: &Decoded<'_>, pc: &Pcurve) -> Result<Option<usize>, Error> {
    let c = d.curves.get(&pc.curve).ok_or(Error::Invalid)?;
    let mut value = None;
    for &point in &c.points {
        let u = d.points.get(&point).ok_or(Error::Invalid)?.coordinates[0];
        let v = if u.is_zero() {
            0
        } else if u.is_one() {
            1
        } else {
            return Ok(None);
        };
        if value.is_some_and(|x| x != v) {
            return Ok(None);
        }
        value = Some(v);
    }
    Ok(value)
}
fn ordered_u<'i>(
    d: &Decoded<'i>,
    edge: &Edge,
    pc: &Pcurve,
    forward: bool,
    b: &mut Budget<'_>,
) -> Result<[Scalar<'i>; 2], Error> {
    let c = d.curves.get(&pc.curve).ok_or(Error::Invalid)?;
    if c.knots.len() != 2 {
        return Err(Error::Unsupported);
    }
    let mut first = None;
    let mut previous = None;
    let mut direction = None;
    for &id in &c.points {
        let u = d.points.get(&id).ok_or(Error::Invalid)?.coordinates[0];
        numbers::unit_coordinate(u, b)?;
        if let Some(previous) = previous {
            let order = numbers::compare(previous, u, b)?;
            if order == std::cmp::Ordering::Equal || direction.is_some_and(|d| d != order) {
                return Err(Error::Unsupported);
            }
            direction = Some(order);
        }
        if first.is_none() {
            first = Some(u)
        }
        previous = Some(u);
    }
    let ends = [
        first.ok_or(Error::Invalid)?,
        previous.ok_or(Error::Invalid)?,
    ];
    Ok(if forward == edge.same_sense {
        ends
    } else {
        [ends[1], ends[0]]
    })
}
fn thread_domain(
    polar: &PolarSupports<'_, '_, '_>,
    strip: &PolarStrip,
    leads: [I; 2],
    norm: I,
    derivative: f64,
    b: &mut Budget<'_>,
) -> Result<f64, Error> {
    let d = polar.grid.source();
    let tensor = &d.tensors[strip.tensor];
    let face = d
        .faces
        .iter()
        .find(|f| f.surface == tensor.entity)
        .ok_or(Error::Invalid)?;
    if face.loops.len() != 1 {
        return Err(Error::Unsupported);
    }
    let ring = &face.loops[0];
    let mut sides = [None, None];
    for (i, c) in ring.iter().enumerate() {
        b.charge(1)?;
        let edge = &d.edges[&c.edge];
        let pc = edge
            .pcurves
            .iter()
            .find(|pc| pc.surface == tensor.entity)
            .ok_or(Error::Invalid)?;
        if let Some(u) = side_u(d, pc)?
            && sides[u].replace(i).is_some()
        {
            return Err(Error::Unsupported);
        }
    }
    let a = sides[0].ok_or(Error::Unsupported)?;
    let z = sides[1].ok_or(Error::Unsupported)?;
    let lead = leads[usize::from(strip.lower == 2)];
    let speed = strip_speed(strip, lead, norm, derivative)?;
    let mut omission = p(0.);
    for (start, end, increasing) in [(a, z, true), (z, a, false)] {
        let side = &d.edges[&ring[start].edge];
        let vertex = side.vertices[usize::from(ring[start].forward)];
        let height = d.points[&d.vertices[&vertex]].coordinates[2];
        let mut cursor = (start + 1) % ring.len();
        let mut previous = Scalar::new(
            if increasing { "0." } else { "1." },
            if increasing { 0. } else { 1. },
        );
        let mut count = 0;
        while cursor != end {
            b.charge(1)?;
            count += 1;
            if count > ring.len() {
                return Err(Error::Invalid);
            }
            let c = &ring[cursor];
            let edge = &d.edges[&c.edge];
            let pc = edge
                .pcurves
                .iter()
                .find(|pc| pc.surface == tensor.entity)
                .ok_or(Error::Invalid)?;
            let plane = edge
                .pcurves
                .iter()
                .find(|pc| pc.surface != tensor.entity)
                .and_then(|pc| d.planes.iter().find(|p| p.entity == pc.surface))
                .ok_or(Error::Unsupported)?;
            let planez = d.points[&plane.origin].coordinates[2];
            if numbers::compare(height, planez, b)? != std::cmp::Ordering::Equal {
                return Err(Error::Unsupported);
            }
            let u = ordered_u(d, edge, pc, c.forward, b)?;
            let direction = numbers::compare(u[0], u[1], b)?;
            if (increasing && direction != std::cmp::Ordering::Less)
                || (!increasing && direction != std::cmp::Ordering::Greater)
            {
                return Err(Error::Unsupported);
            }
            omission = omission + p((u[0].interval() - previous.interval()).magnitude()) * p(speed);
            previous = u[1];
            cursor = (cursor + 1) % ring.len();
        }
        if count == 0 {
            return Err(Error::Unsupported);
        }
        let target = p(if increasing { 1. } else { 0. });
        omission = omission + p((target - previous.interval()).magnitude()) * p(speed);
    }
    if !omission.valid() {
        return Err(Error::Unsupported);
    }
    Ok(omission.hi)
}
fn plain_domains(polar: &PolarSupports<'_, '_, '_>, b: &mut Budget<'_>) -> Result<(), Error> {
    let d = polar.grid.source();
    for (index, t) in d.tensors.iter().enumerate() {
        if polar.strips.iter().any(|s| s.tensor == index) {
            continue;
        }
        let face = d
            .faces
            .iter()
            .find(|f| f.surface == t.entity)
            .ok_or(Error::Invalid)?;
        let mut caps = [0_usize; 2];
        let mut sides = [0_usize; 2];
        for coedge in face.loops.iter().flatten() {
            b.charge(1)?;
            let edge = &d.edges[&coedge.edge];
            let pc = edge
                .pcurves
                .iter()
                .find(|p| p.surface == t.entity)
                .ok_or(Error::Invalid)?;
            let curve = &d.curves[&pc.curve];
            if curve.degree != 1 || curve.points.len() != 2 {
                return Err(Error::Unsupported);
            }
            let a = d.points[&curve.points[0]].coordinates();
            let z = d.points[&curve.points[1]].coordinates();
            if let Some(u) = side_u(d, pc)? {
                caps[u] += 1;
                let plane = edge
                    .pcurves
                    .iter()
                    .find(|p| p.surface != t.entity)
                    .and_then(|pc| d.planes.iter().find(|p| p.entity == pc.surface))
                    .ok_or(Error::Unsupported)?;
                if plane.entity != polar.planes[1 + u].0
                    || numbers::compare(a[1], t.knots_v[0], b)? != std::cmp::Ordering::Equal
                    || numbers::compare(z[1], *t.knots_v.last().ok_or(Error::Invalid)?, b)?
                        != std::cmp::Ordering::Equal
                {
                    return Err(Error::Unsupported);
                }
            } else {
                if numbers::compare(a[1], z[1], b)? != std::cmp::Ordering::Equal
                    || !(a[0].is_zero() && z[0].is_one())
                {
                    return Err(Error::Unsupported);
                }
                let side = if numbers::compare(a[1], t.knots_v[0], b)? == std::cmp::Ordering::Equal
                {
                    0
                } else if numbers::compare(a[1], *t.knots_v.last().ok_or(Error::Invalid)?, b)?
                    == std::cmp::Ordering::Equal
                {
                    1
                } else {
                    return Err(Error::Unsupported);
                };
                sides[side] += 1;
            }
        }
        if caps != [1, 1] || sides != [1, 1] {
            return Err(Error::Unsupported);
        }
    }
    Ok(())
}
fn weighted_curve(
    d: &Decoded<'_>,
    curve: &Curve<'_>,
    b: &mut Budget<'_>,
) -> Result<[Polynomial; 4], Error> {
    let scale = curve
        .weights
        .iter()
        .map(|w| w.interval().hi)
        .fold(1., f64::max);
    let mut out = std::array::from_fn(|_| Vec::with_capacity(curve.points.len()));
    for (i, &id) in curve.points.iter().enumerate() {
        b.charge(96)?;
        let q = d.points.get(&id).ok_or(Error::Invalid)?.coordinates();
        if q.len() != 3 {
            return Err(Error::Invalid);
        }
        let w = curve.weights.get(i).map_or(p(1.), Scalar::interval) / p(scale);
        for axis in 0..3 {
            out[axis].push(q[axis].interval() * w)
        }
        out[3].push(w);
    }
    Ok(out)
}
fn cap_edge(
    d: &Decoded<'_>,
    edge: &Edge,
    lift: &EdgeLift,
    algebra: &Algebra,
    b: &mut Budget<'_>,
) -> Result<(), Error> {
    let c = &d.curves[&edge.curve];
    let net = weighted_curve(d, c, b)?;
    let domain = lift.parameter_bounds().map(|v| I::new(v[0], v[1]));
    // The earned end residuals enclose the exact corrections from the
    // transmitted rational curve to its two shared literal vertices. Proving
    // all coordinatewise choices in this box also proves the actual correction.
    let errors = lift.vertex_errors_mm();
    let delta = errors.map(|e| [I::new(-e, e); 2]);
    let normals = [
        [1., 0.],
        [0., 1.],
        [-1., 0.],
        [0., -1.],
        [1., 1.],
        [-1., 1.],
        [-1., -1.],
        [1., -1.],
    ];
    let mut halfplanes = [true; 8];
    let mut spans = 0;
    for (span, knots) in c.knots.windows(2).enumerate() {
        let a = numbers::exact_interval(knots[0], b)?;
        let z = numbers::exact_interval(knots[1], b)?;
        if z.hi <= domain[0].lo || a.lo >= domain[1].hi {
            continue;
        }
        let lo = I::new(domain[0].lo.max(a.lo), domain[0].hi.max(a.hi));
        let hi = I::new(domain[1].lo.min(z.lo), domain[1].hi.min(z.hi));
        let u = (lo - a) / (z - a);
        let v = (hi - a) / (z - a);
        let u = I::new(u.lo.max(0.), u.hi.min(1.));
        let v = I::new(v.lo.max(0.), v.hi.min(1.));
        spans += 1;
        let mut restricted = Vec::new();
        for axis in &net {
            restricted.push(poly::restrict(
                &axis[span * c.degree..span * c.degree + c.degree + 1],
                u,
                v,
                b,
            )?);
        }
        let s = [
            (lo - domain[0]) / (domain[1] - domain[0]),
            (hi - domain[0]) / (domain[1] - domain[0]),
        ];
        let mut xy = Vec::new();
        for axis in 0..2 {
            let correction = s.map(|s| delta[0][axis] * (p(1.) - s) + delta[1][axis] * s);
            let add = algebra.multiply(&restricted[3], &correction, b)?;
            let mut base = algebra.multiply(&restricted[axis], &[p(1.), p(1.)], b)?;
            poly::add(&mut base, &add, b)?;
            xy.push(base);
        }
        for (i, n) in normals.iter().enumerate() {
            let mut projection = xy[0].clone();
            poly::scale(&mut projection, p(n[0]), b)?;
            let mut y = xy[1].clone();
            poly::scale(&mut y, p(n[1]), b)?;
            poly::add(&mut projection, &y, b)?;
            halfplanes[i] &= poly::lower(&projection) > 0.;
        }
    }
    // A common open half-plane excludes the origin and gives a homotopy,
    // relative to its exact shared endpoints, from this corrected finite curve
    // to the vertex chord. The complete chord polygon is checked independently;
    // no positive angular derivative of the approximate raw curve is required.
    if !halfplanes.iter().any(|&v| v) || spans == 0 {
        return Err(Error::Unsupported);
    }
    Ok(())
}
fn polygon_winding(
    polar: &PolarSupports<'_, '_, '_>,
    ring: &[Coedge],
    b: &mut Budget<'_>,
) -> Result<i32, Error> {
    let d = polar.grid.source();
    let mut winding = 0;
    for c in ring {
        let edge = &d.edges[&c.edge];
        let ids = if c.forward {
            edge.vertices
        } else {
            [edge.vertices[1], edge.vertices[0]]
        };
        let a = polar
            .grid
            .points()
            .get(&d.vertices[&ids[0]])
            .ok_or(Error::Invalid)?;
        let z = polar
            .grid
            .points()
            .get(&d.vertices[&ids[1]])
            .ok_or(Error::Invalid)?;
        let x = integer(&a[0], b)?;
        let y = integer(&a[1], b)?;
        let xx = integer(&z[0], b)?;
        let yy = integer(&z[1], b)?;
        b.charge(
            32 * (1
                + ((x.bits().max(y.bits()).max(xx.bits()).max(yy.bits()) as usize).div_ceil(64))
                    .pow(2)),
        )?;
        let determinant = &x * &yy - &y * &xx;
        if y <= BigInt::zero() && yy > BigInt::zero() && determinant > BigInt::zero() {
            winding += 1
        }
        if yy <= BigInt::zero() && y > BigInt::zero() && determinant < BigInt::zero() {
            winding -= 1
        }
    }
    Ok(winding)
}
impl<'p, 'g, 'd, 'l, 'i> HealedMaterial<'p, 'g, 'd, 'l, 'i> {
    pub fn prepare(
        polar: &'p PolarSupports<'g, 'd, 'i>,
        lifts: &'l SurfaceLifts<'d, 'i>,
        options: HealingOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, Error> {
        if !std::ptr::eq(polar.grid.source(), lifts.source())
            || options.maximum_work == 0
            || options.maximum_work > 128_000_000
            || options.maximum_bytes == 0
            || options.maximum_bytes > 64 * 1024 * 1024
            || !options.maximum_sewing_error_mm.is_finite()
            || options.maximum_sewing_error_mm <= 0.
            || options.maximum_sewing_error_mm > 1e-7
        {
            return Err(Error::Invalid);
        }
        let d = polar.grid.source();
        let inherited = polar.retained_bytes() + lifts.retained_bytes() - d.retained_budget_bytes();
        let mut b = Budget {
            work: 0,
            remaining: options.maximum_work,
            bytes: inherited,
            limit: options.maximum_bytes,
            interrupt,
        };
        b.charge(1)?;
        b.retain(4 * 1024 * 1024 + size_of::<Self>())?;
        let algebra = Algebra::new(&mut b)?;
        let (norm, jacobian, derivative) = angular(polar, &algebra, &mut b)?;
        let unit = integer(polar.grid.denominator(), &mut b)?;
        let leads = [
            physical(&polar.leads[0], &unit, &mut b)?,
            physical(&polar.leads[1], &unit, &mut b)?,
        ];
        let plain = physical(&polar.plain_radius, &unit, &mut b)?;
        let mut radius = [I::new(f64::INFINITY, 0.), plain, I::new(f64::INFINITY, 0.)];
        let mut omission = 0_f64;
        for strip in &polar.strips {
            b.charge(1)?;
            let group = if strip.lower == 0 { 0 } else { 2 };
            for r in &strip.radius {
                radius[group].lo = radius[group].lo.min(r.lo);
                radius[group].hi = radius[group].hi.max(r.hi);
            }
            omission = omission.max(thread_domain(
                polar, strip, leads, norm, derivative, &mut b,
            )?);
        }
        plain_domains(polar, &mut b)?;
        let radii = radius.map(|r| r * norm);
        if radii.iter().any(|r| !r.valid() || r.lo <= 0.) {
            return Err(Error::Unsupported);
        }
        let mut sewing = p(omission);
        let mut largest = 0_f64;
        for lift in lifts.edges() {
            let edge = &d.edges[&lift.edge_entity()];
            let mut projection = p(1.);
            for pc in &edge.pcurves {
                if let Some(strip) = polar
                    .strips
                    .iter()
                    .find(|s| d.tensors[s.tensor].entity == pc.surface)
                {
                    let rmax = strip.radius.iter().map(|r| r.hi).fold(0., f64::max);
                    let lead = leads[usize::from(strip.lower != 0)];
                    let factor = p(1.) + p(rmax) * p(derivative) / lead;
                    projection.hi = projection.hi.max(factor.hi);
                }
            }
            largest = largest.max(
                (p(lift.field_bounds()[0]) * projection
                    + p(lift.vertex_errors_mm().iter().copied().fold(0., f64::max)))
                .hi,
            );
        }
        sewing = sewing + p(largest);
        if !sewing.valid() || sewing.hi > options.maximum_sewing_error_mm {
            return Err(Error::Unsupported);
        }
        if (p(radii[0].hi) + p(2.) * sewing).hi >= radii[1].lo
            || (p(radii[1].hi) + p(2.) * sewing).hi >= radii[2].lo
        {
            return Err(Error::Unsupported);
        }
        let expected = [vec![0], vec![1, 0], vec![2, 1], vec![2]];
        for (plane_index, (plane, _)) in polar.planes.iter().enumerate() {
            let mut found = d.faces.iter().filter(|f| f.surface == *plane);
            let face = found.next().ok_or(Error::Unsupported)?;
            if found.next().is_some()
                || face.same_sense != (plane_index == 3)
                || face.loops.len() != expected[plane_index].len()
            {
                return Err(Error::Unsupported);
            }
            for (ring_index, ring) in face.loops.iter().enumerate() {
                let sign =
                    if face.same_sense { 1. } else { -1. } * if ring_index == 0 { 1. } else { -1. };
                let winding = polygon_winding(polar, ring, &mut b)?;
                if winding != sign as i32 {
                    return Err(Error::Unsupported);
                }
                for coedge in ring {
                    let edge = &d.edges[&coedge.edge];
                    let plane_surface = d
                        .planes
                        .iter()
                        .find(|p| p.entity == *plane)
                        .ok_or(Error::Invalid)?;
                    let height = d.points[&plane_surface.origin].coordinates[2];
                    for vertex in edge.vertices {
                        if numbers::compare(
                            d.points[&d.vertices[&vertex]].coordinates[2],
                            height,
                            &mut b,
                        )? != std::cmp::Ordering::Equal
                        {
                            return Err(Error::Unsupported);
                        }
                    }
                    let side = edge
                        .pcurves
                        .iter()
                        .find(|p| p.surface != *plane)
                        .ok_or(Error::Invalid)?;
                    let tensor = d
                        .tensors
                        .iter()
                        .position(|t| t.entity == side.surface)
                        .ok_or(Error::Unsupported)?;
                    let group = polar
                        .strips
                        .iter()
                        .find(|s| s.tensor == tensor)
                        .map_or(1, |s| s.lower);
                    if group != expected[plane_index][ring_index] {
                        return Err(Error::Unsupported);
                    }
                    let lift = lifts
                        .edges()
                        .iter()
                        .find(|l| l.edge_entity() == edge.entity)
                        .ok_or(Error::Invalid)?;
                    cap_edge(d, edge, lift, &algebra, &mut b)?;
                }
            }
        }
        let pi = I::new(
            std::f64::consts::PI.next_down(),
            std::f64::consts::PI.next_up(),
        );
        let mut volume = p(0.);
        for (i, radius) in radii.iter().enumerate() {
            let height = physical(
                &(&polar.planes[i + 1].1 - &polar.planes[i].1),
                &unit,
                &mut b,
            )?;
            if height.lo <= 0. {
                return Err(Error::Unsupported);
            }
            volume = volume + pi * *radius * *radius * height;
        }
        if !volume.valid() || volume.lo <= 0. {
            return Err(Error::Unsupported);
        }
        b.bytes -= 4 * 1024 * 1024;
        Ok(Self {
            polar,
            lifts,
            norm,
            jacobian,
            sewing: sewing.hi,
            volume,
            radii,
            work: b.work,
            bytes: b.bytes,
        })
    }
    pub fn polar_supports(&self) -> &PolarSupports<'g, 'd, 'i> {
        self.polar
    }
    pub fn finite_lifts(&self) -> &SurfaceLifts<'d, 'i> {
        self.lifts
    }
    pub fn sewing_error_bound_mm(&self) -> f64 {
        self.sewing
    }
    pub fn healed_volume_bounds_mm3(&self) -> [f64; 2] {
        [self.volume.lo, self.volume.hi]
    }
    pub fn component_radius_bounds_mm(&self) -> [[f64; 2]; 3] {
        self.radii.map(|r| [r.lo, r.hi])
    }
    pub fn angular_norm_bounds(&self) -> [f64; 2] {
        [self.norm.lo, self.norm.hi]
    }
    pub fn minimum_angular_jacobian(&self) -> f64 {
        self.jacobian
    }
    pub fn admission_work_used(&self) -> usize {
        self.work
    }
    pub fn aggregate_known_work_used(&self) -> usize {
        self.lifts.work_used()
            + self.polar.grid.preparation_work_used()
            + self.polar.preparation_work_used()
            + self.work
    }
    pub fn retained_budget_bytes(&self) -> usize {
        self.bytes
    }
}
