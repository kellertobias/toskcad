//! Actual finite counterbored clamp versus original regular screw and insert.
//! Related receiver/receiver material intersection is deliberately not certified.
use super::*;
use crate::{
    Rational, chamfered_bores::ChamferedBores, headed_screw::Screw, polygon::PolygonPrism,
};
pub struct InsertClampSources<'a> {
    pub prism: &'a PolygonPrism,
    pub cut: &'a ChamferedBores,
    pub clamp_frame: Frame,
    pub screw: &'a Screw,
    pub screw_frame: Frame,
    pub insert: &'a Insert,
    pub insert_frame: Frame,
}
#[derive(Clone, Copy)]
pub struct InsertClampRequirements {
    pub minimum_wall_mm: f64,
    pub minimum_bearing_land_mm: f64,
    pub maximum_bearing_gap_mm: f64,
    pub linear_tolerance_mm: f64,
    pub minimum_feature_mm: f64,
}
pub struct InsertClampContact {
    pub bore_radius_mm: f64,
    pub seat_radius_mm: f64,
    pub seat_depth_mm: f64,
    pub source_height_mm: f64,
    pub bore_mouths_mm: [[[f64; 2]; 3]; 2],
    pub source_bearing_radii_mm: [f64; 2],
    pub bearing_land_lower_mm: f64,
    /// Nonnegative corresponding-face gap; never a permitted intrusion.
    pub bearing_face_gap_mm: [f64; 2],
    pub shaft_radial_clearance_lower_mm: f64,
    pub head_radial_clearance_lower_mm: f64,
    pub insert_to_clamp_gap_mm: [f64; 2],
}
fn qr(v: Rational) -> Q {
    let (n, d) = v.ratio();
    Q::new(n.into(), d.into())
}
fn exact(v: &Q) -> R<I> {
    let x = v.to_f64().ok_or(Error::Unsupported)?;
    if q(x)? == *v { Ok(ip(x)) } else { enclosure(v) }
}
struct Map {
    o: [Q; 3],
    a: [[Q; 3]; 3],
    low: I,
    high: I,
    target: Metric,
}
impl Map {
    fn new(from: Frame, to: Frame, interrupt: &dyn Fn() -> bool) -> R<Self> {
        let source = scale(from)?;
        let target = scale(to)?;
        let (o, a) = affine(to, from, interrupt)?;
        let low = positive(source.minimum_scale / target.maximum_scale)?;
        let high = positive(source.maximum_scale / target.minimum_scale)?;
        Ok(Self {
            o,
            a,
            low,
            high,
            target,
        })
    }
    fn center_reach(&self, z: [f64; 2], center: &[Q; 2], interrupt: &dyn Fn() -> bool) -> R<I> {
        let mut reach = Q::zero();
        for (i, coordinate) in center.iter().enumerate() {
            let x = &self.o[i] - coordinate + &self.a[2][i] * q(z[0])?;
            let y = &self.o[i] - coordinate + &self.a[2][i] * q(z[1])?;
            reach += x.abs().max(y.abs());
        }
        let result = exact(&bounded(reach)?)?;
        // Exact rational evaluation may consume large binary exponents.
        // Poll before publishing this completed source-domain query.
        checkpoint(interrupt)?;
        Ok(result)
    }
}
fn rectangle(prism: &PolygonPrism, interrupt: &dyn Fn() -> bool) -> R<[[Q; 2]; 2]> {
    let [ring] = prism.rings() else {
        return Err(Error::Unsupported);
    };
    if ring.len() != 4 {
        return Err(Error::Unsupported);
    }
    let mut low = ring[0].map(qr);
    let mut high = low.clone();
    for p in ring {
        checkpoint(interrupt)?;
        for i in 0..2 {
            let x = qr(p[i]);
            if x.abs() > q(1e7)? {
                return Err(Error::Limit);
            }
            low[i] = low[i].clone().min(x.clone());
            high[i] = high[i].clone().max(x);
        }
    }
    // A genuine four-edge axis-aligned rectangle, not its bounding box.
    // The sealed simple positive ring has four distinct vertices/edges.
    for (p, r) in ring.iter().zip(ring.iter().cycle().skip(1)) {
        if (p[0] == r[0]) == (p[1] == r[1])
            || (0..2).any(|i| qr(p[i]) != low[i] && qr(p[i]) != high[i])
        {
            return Err(Error::Unsupported);
        }
    }
    Ok([low, high])
}
/// Earn complete finite head/shaft clearance and an actual projected bearing
/// annulus. The insert is independently proved behind the ACTUAL rational end
/// plane, not a rounded world endpoint. Geometry is never moved by this query.
pub fn certify_insert_clamp(
    input: InsertClampSources<'_>,
    req: InsertClampRequirements,
    interrupt: &dyn Fn() -> bool,
) -> R<InsertClampContact> {
    let InsertClampSources {
        prism,
        cut,
        clamp_frame,
        screw,
        screw_frame,
        insert,
        insert_frame,
    } = input;
    checkpoint(interrupt)?;
    if [
        req.minimum_wall_mm,
        req.minimum_bearing_land_mm,
        req.maximum_bearing_gap_mm,
        req.linear_tolerance_mm,
        req.minimum_feature_mm,
    ]
    .into_iter()
    .any(|v| !v.is_finite() || !(0. ..=1e7).contains(&v))
        || req.minimum_wall_mm == 0.
        || req.minimum_bearing_land_mm == 0.
        || req.minimum_feature_mm == 0.
        || !(1e-7..=1e-3).contains(&req.linear_tolerance_mm)
        || req.maximum_bearing_gap_mm > req.linear_tolerance_mm
    {
        return Err(Error::Invalid);
    }
    if !cut.belongs_to(prism) || !cut.chamfers().is_empty() || !cut.countersinks().is_empty() {
        return Err(Error::Unsupported);
    }
    let [seat] = cut.counterbores() else {
        return Err(Error::Unsupported);
    };
    let [circle] = prism.circular_holes() else {
        return Err(Error::Unsupported);
    };
    if seat.end {
        return Err(Error::Unsupported);
    }
    let [low, high] = rectangle(prism, interrupt)?;
    let height = qr(prism.height());
    if height > q(1e7)? {
        return Err(Error::Limit);
    }
    let center = circle.center().map(qr);
    // counterbored() owns one stack per circular hole. With exactly one flat
    // start cut its actual layers are [mouth,shoulderOuter,shoulderInner,exit].
    let stack = &cut.stacks()[0];
    let inner = &stack.layers[2];
    // Each original Rational has <=64-bit limbs; its nonzero rounded cache
    // is >=2^-64 and <=1e7+1. The common dyadic denominator is <=2^117;
    // combining two original denominators leaves <1024 total rational bits,
    // far below bounded()'s 16384-bit guard. No query-dependent size grows here.
    let center_error: Q = (0..2)
        .map(|i| (&center[i] - q(stack.center_mm[i]).unwrap()).abs())
        .sum();
    let bore_lower = qr(circle.radius()).min(q(inner.radius_mm)? - &center_error);
    let bore_upper = qr(circle.radius()).max(q(inner.radius_mm)? + &center_error);
    let seat_lower = q(seat.radius_mm)? - &center_error;
    let map = Map::new(screw_frame, clamp_frame, interrupt)?;
    let wall = (ip(req.minimum_wall_mm) / map.target.minimum_scale).hi;
    for i in 0..2 {
        let protected = q(seat.radius_mm)? + &center_error + q(wall)?;
        if &center[i] - &protected < low[i] || &center[i] + protected > high[i] {
            return Err(Error::Unsupported);
        }
    }
    if map.a[2][2] <= Q::zero() {
        return Err(Error::Unsupported);
    }
    let spec = screw.spec();
    let plain = q(spec.shank_length_mm)?;
    let (head_inner, head_outer) = head_radii(screw)?;
    let plane = bounded(&map.o[2] - &map.a[2][2] * &plain)?;
    let tilt = bounded((map.a[0][2].abs() + map.a[1][2].abs()) * q(head_outer)?)?;
    let seat_z = q(seat.depth_mm)?;
    let gap_low = bounded(&seat_z - &plane - &tilt)?;
    let gap_high = bounded(&seat_z - plane + tilt)?;
    // The entire actual head is on the free side of the shoulder; tolerance
    // bounds only a real nonnegative gap, never interpenetration.
    if gap_low < Q::zero() {
        return Err(Error::Unsupported);
    }
    let face_gap = if gap_low.is_zero() && gap_high.is_zero() {
        ip(0.)
    } else {
        I::new(exact(&gap_low)?.lo.max(0.), exact(&gap_high)?.hi)
            * I::new(map.target.minimum_scale.lo, map.target.maximum_scale.hi)
    };
    if face_gap.hi > req.maximum_bearing_gap_mm {
        return Err(Error::Unsupported);
    }
    let head_z_low = finite(-ip(spec.shank_length_mm) - ip(spec.head.height_mm()))?.lo;
    let head_reach = map.high * ip(head_outer)
        + map.center_reach([head_z_low, -spec.shank_length_mm], &center, interrupt)?;
    let head_clearance = (exact(&seat_lower)? - head_reach) * map.target.minimum_scale;
    if head_clearance.lo < req.minimum_feature_mm {
        return Err(Error::Unsupported);
    }
    let r = screw.thread().source_coordinate_bounds_mm()[0][1].max(if spec.shank_length_mm > 0. {
        spec.shank_radius_mm
    } else {
        0.
    });
    let shaft_reach = map.high * ip(r)
        + map.center_reach(
            [-spec.shank_length_mm, screw.thread().spec().length_mm],
            &center,
            interrupt,
        )?;
    let shaft_clearance = (exact(&bore_lower)? - shaft_reach) * map.target.minimum_scale;
    if shaft_clearance.lo < req.minimum_feature_mm {
        return Err(Error::Unsupported);
    }
    let off = map.center_reach([-spec.shank_length_mm; 2], &center, interrupt)?;
    // Removing the Z row from the full map loses at most its squared norm.
    // This earns the actual 2D head-plane projection singular-value bound.
    let zrow = exact(&(&map.a[0][2] * &map.a[0][2] + &map.a[1][2] * &map.a[1][2]))?;
    let xy_low = positive(map.low.square() - zrow)?.sqrt();
    let inner_radius = screw.thread().source_coordinate_bounds_mm()[0][1].max(spec.shank_radius_mm);
    let land_low = (map.high * ip(inner_radius) + off)
        .hi
        .max(exact(&bore_upper)?.hi);
    let land_high = (xy_low * ip(head_inner) - off)
        .lo
        .min(exact(&seat_lower)?.lo);
    let land = finite((ip(land_high) - ip(land_low)) * map.target.minimum_scale)?;
    if land.lo < req.minimum_bearing_land_mm {
        return Err(Error::Unsupported);
    }
    let other = Map::new(insert_frame, clamp_frame, interrupt)?;
    if other.a[2][2] <= Q::zero() {
        return Err(Error::Unsupported);
    }
    // All binary64 frame inputs share a 2^1074 denominator. Cofactor/det
    // ratios have <=3232-bit common numerator/denominator; adding two entries
    // retains the SAME determinant denominator. Multiplying the finite radius
    // adds <=1100/1074 bits, hence <8700 total bits (<16384), without a
    // repeated size guard. The actual constructor scale proof remains above.
    let cap = (other.a[0][2].abs() + other.a[1][2].abs()) * q(insert.outer().maximum_radius_mm())?;
    let near = bounded(&other.o[2] - &cap - &height)?;
    let far = bounded(&other.o[2] + cap - &height)?;
    if near < Q::zero() {
        return Err(Error::Unsupported);
    }
    let insert_gap = if near.is_zero() && far.is_zero() {
        ip(0.)
    } else {
        I::new(exact(&near)?.lo.max(0.), exact(&far)?.hi)
            * I::new(other.target.minimum_scale.lo, other.target.maximum_scale.hi)
    };
    let bore_mouths_mm = [
        world_point(
            clamp_frame,
            [center[0].clone(), center[1].clone(), Q::zero()],
            interrupt,
        )?,
        world_point(
            clamp_frame,
            [center[0].clone(), center[1].clone(), height],
            interrupt,
        )?,
    ];
    checkpoint(interrupt)?;
    Ok(InsertClampContact {
        bore_mouths_mm,
        bore_radius_mm: inner.radius_mm,
        seat_radius_mm: seat.radius_mm,
        seat_depth_mm: seat.depth_mm,
        source_height_mm: prism.height().as_f64(),
        source_bearing_radii_mm: [land_low, land_high],
        bearing_land_lower_mm: land.lo,
        bearing_face_gap_mm: [face_gap.lo.max(0.), face_gap.hi],
        shaft_radial_clearance_lower_mm: shaft_clearance.lo,
        head_radial_clearance_lower_mm: head_clearance.lo,
        insert_to_clamp_gap_mm: [insert_gap.lo.max(0.), insert_gap.hi],
    })
}

#[cfg(test)]
mod tests;
