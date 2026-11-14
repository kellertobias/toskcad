//! Join actual source certificates into a geometric Hanger stack. Material
//! halfspaces and cap lands are independent proofs; tolerance bounds a real
//! nonnegative gap and never grants permission for intrusion.
use super::super::insert::mating::phase_error;
use super::*;
use crate::{annulus::Annulus, hex_nut::Nut};
#[cfg(test)]
mod tests;
#[derive(Clone, Copy)]
pub struct PlacedAnnulus<'a> {
    pub source: &'a Annulus,
    pub frame: Frame,
}
pub struct StackSources<'a> {
    pub wood: Stock<'a>,
    pub clamp: ThroughStock<'a>,
    pub hanger: &'a Hanger,
    pub hanger_frame: Frame,
    /// Ordered from the clamp toward the nut.
    pub washers: &'a [PlacedAnnulus<'a>],
    pub nut: &'a Nut,
    pub nut_frame: Frame,
}
#[derive(Clone, Copy)]
pub struct StackRequirements {
    pub wood: Requirements,
    pub minimum_clamp_wall_mm: f64,
    pub passage: PassageRequirements,
    pub mating: HangerMatingRequirements,
    pub minimum_bearing_land_mm: f64,
    pub maximum_bearing_gap_mm: f64,
}
#[derive(Clone, Debug)]
pub struct BearingContact {
    /// Actual A-face source-chart radii, not assumed world radii.
    pub source_annular_radii_mm: [f64; 2],
    pub radial_land_lower_mm: f64,
    /// Distance between corresponding real face points, along the actual B
    /// axial column. This is an upper bound, not a minimum-distance estimate.
    pub correspondence_gap_upper_mm: f64,
}
#[derive(Clone, Debug)]
pub struct HangerStack {
    pub wood: HangerContact,
    pub nut: HangerMating,
    pub clamp: HangerPassage,
    pub washer_passages: Vec<HangerPassage>,
    /// Wood/clamp, then the ordered washer chain, then the nut.
    pub bearings: Vec<BearingContact>,
}
#[derive(Clone, Copy)]
enum Material<'a> {
    Cells(&'a OrthogonalSolid),
    Ring(&'a Annulus),
    Nut(&'a Nut),
}
struct Placed<'a> {
    material: Material<'a>,
    frame: Frame,
    metric: Metric,
}
impl<'a> Placed<'a> {
    fn new(material: Material<'a>, frame: Frame) -> R<Self> {
        Ok(Self {
            material,
            frame,
            metric: scale(frame)?,
        })
    }
}
#[derive(Clone)]
struct Face {
    frame: Frame,
    axis: usize,
    center: [Q; 3],
    outward: i8,
    inner: f64,
    outer: f64,
    minimum_scale: I,
}
struct Map {
    origin: [Q; 3],
    columns: [[Q; 3]; 3],
}
/// Exact decoding for coordinates already earned finite by their source
/// constructor. Unknown public inputs still pass q()/source validation.
fn sealed_binary(value: f64) -> Q {
    let (mantissa, exponent, sign) = num_traits::Float::integer_decode(value);
    let numerator = num_bigint::BigInt::from(mantissa) * num_bigint::BigInt::from(sign);
    if exponent >= 0 {
        Q::from_integer(numerator << exponent as usize)
    } else {
        Q::new(
            numerator,
            num_bigint::BigInt::from(1) << usize::from(exponent.unsigned_abs()),
        )
    }
}
fn map(to: Frame, from: Frame, interrupt: &dyn Fn() -> bool) -> R<Map> {
    let (origin, columns) = affine(to, from, interrupt)?;
    Ok(Map { origin, columns })
}
fn coordinate(map: &Map, axis: usize, point: [Q; 3]) -> R<Q> {
    bounded((0..3).fold(map.origin[axis].clone(), |s, j| {
        s + &map.columns[j][axis] * &point[j]
    }))
}
/// Complete B material lies in A's outward closed halfspace. The uncut cell
/// union and convex prism/cylinder bounds only prove separation; they are
/// never treated as occupied material when proving the bearing land.
fn outside(a: &Face, b: &Placed<'_>, interrupt: &dyn Fn() -> bool) -> R<Map> {
    let m = map(a.frame, b.frame, interrupt)?;
    let k = a.axis;
    let sign = q(f64::from(a.outward))?;
    match b.material {
        Material::Cells(stock) => {
            for cell in stock.cells() {
                checkpoint(interrupt)?;
                for corner in 0..8 {
                    let p = std::array::from_fn(|j| {
                        if corner & (1 << j) == 0 {
                            cell.min[j]
                        } else {
                            cell.max[j]
                        }
                    });
                    let p = vector(p)?;
                    if bounded(&sign * (coordinate(&m, k, p)? - &a.center[k]))? < Q::zero() {
                        return Err(Error::Unsupported);
                    }
                }
            }
        }
        Material::Nut(nut) => {
            for xy in crate::headed_screw::hex_vertices(nut.across_flats_mm()) {
                checkpoint(interrupt)?;
                for z in [0., nut.thread().spec().length_mm] {
                    if bounded(
                        &sign * (coordinate(&m, k, [q(xy[0])?, q(xy[1])?, q(z)?])? - &a.center[k]),
                    )? < Q::zero()
                    {
                        return Err(Error::Unsupported);
                    }
                }
            }
        }
        Material::Ring(ring) => {
            let base = bounded(&sign * (&m.origin[k] - &a.center[k]))?;
            let along = bounded(&sign * &m.columns[2][k] * q(ring.height_mm())?)?;
            let gap = bounded(base + along.min(Q::zero()))?;
            let square =
                bounded(&m.columns[0][k] * &m.columns[0][k] + &m.columns[1][k] * &m.columns[1][k])?;
            let r = q(ring.outer_radius_mm())?;
            if gap < Q::zero() || bounded(&gap * &gap)? < bounded(&r * &r * square)? {
                return Err(Error::Unsupported);
            }
        }
    }
    Ok(m)
}
fn nut_land(nut: &Nut, interrupt: &dyn Fn() -> bool) -> R<[f64; 2]> {
    let thread = nut.thread();
    let allowance = phase_error::measure(thread, interrupt)?;
    let slope = thread
        .slope_at_phase(0., 1)
        .abs()
        .max(thread.slope_at_phase(0., 3).abs());
    // Sealed Thread dimensions are <=1e7 and flank angles lie in [1,80]
    // degrees. Even the root-circle centre bound is <9e8; phase_error earns
    // allowance<.001 and the retained flank slopes are <58. This positive
    // outward sum is therefore finite independently of a requested fit.
    let inner =
        (ip(thread.source_coordinate_bounds_mm()[0][1]) + ip(slope) * ip(allowance.total_mm)).hi;
    let vertices = crate::headed_screw::hex_vertices(nut.across_flats_mm());
    let mut outer = f64::INFINITY;
    for (p, r) in vertices.iter().zip(vertices.iter().cycle().skip(1)) {
        checkpoint(interrupt)?;
        let [px, py] = [q(p[0])?, q(p[1])?];
        let [rx, ry] = [q(r[0])?, q(r[1])?];
        let det = bounded(&px * &ry - &py * &rx)?;
        let dx = bounded(rx - px)?;
        let dy = bounded(ry - py)?;
        let length = bounded(&dx * &dx + &dy * &dy)?;
        // Nut::new earns x=AF/2>major>0 and caps AF/sqrt(3)<=1e7.
        // hex_vertices uses x and y=x/sqrt(3)>0, then exact sign flips and
        // doubling. All six literal edges have positive determinant and
        // nonzero squared length, including subnormal admissible dimensions.
        let radius = positive(passage::directed(&bounded(&det * &det / length)?)?.sqrt())?;
        outer = outer.min(radius.lo);
    }
    if inner >= outer {
        return Err(Error::Unsupported);
    }
    Ok([inner, outer])
}
fn cap(b: &Placed<'_>, axis: usize, center: [Q; 3], outward: i8, radii: [f64; 2]) -> Face {
    Face {
        frame: b.frame,
        axis,
        center,
        outward,
        inner: radii[0],
        outer: radii[1],
        minimum_scale: b.metric.minimum_scale,
    }
}
fn ring_cap(b: &Placed<'_>, ring: &Annulus, end: bool) -> R<Face> {
    Ok(cap(
        b,
        2,
        [
            Q::zero(),
            Q::zero(),
            q(if end { ring.height_mm() } else { 0. })?,
        ],
        if end { 1 } else { -1 },
        [ring.inner_radius_mm(), ring.outer_radius_mm()],
    ))
}
fn member_near_cap(a: &Face, b: &Placed<'_>, m: &Map, interrupt: &dyn Fn() -> bool) -> R<Face> {
    let derivative = bounded(q(f64::from(a.outward))? * &m.columns[2][a.axis])?;
    if derivative.abs() < q(0.5)? {
        return Err(Error::Unsupported);
    }
    let end = derivative < Q::zero();
    match b.material {
        Material::Ring(ring) => ring_cap(b, ring, end),
        Material::Nut(nut) => Ok(cap(
            b,
            2,
            [
                Q::zero(),
                Q::zero(),
                sealed_binary(if end {
                    nut.thread().spec().length_mm
                } else {
                    0.
                }),
            ],
            if end { 1 } else { -1 },
            nut_land(nut, interrupt)?,
        )),
        Material::Cells(_) => Err(Error::Invalid),
    }
}
fn bearing(
    a: &Face,
    bface: &Face,
    b: &Placed<'_>,
    m: &Map,
    req: StackRequirements,
    interrupt: &dyn Fn() -> bool,
) -> R<BearingContact> {
    let k = a.axis;
    let j = bface.axis;
    let derivative = bounded(&m.columns[j][k] * q(f64::from(a.outward * bface.outward))?)?;
    if derivative > -q(0.5)? {
        return Err(Error::Unsupported);
    }
    let rb = [(j + 1) % 3, (j + 2) % 3];
    // Capacity theorem for these actual maps, not arbitrary supplied Qs:
    // scale() earns |axis|<=2, |origin|<=1e7. With D=2^1074, axes M/D
    // have integer |M|<=2D. det(M) and adj(M)*N each have <3229 bits;
    // adj(M)*(OB-OA) has <3252. All relative columns/origins share det(M).
    // A finite binary64 source coordinate has numerator <2^2098 over D.
    // Thus coordinate minus another source coordinate has <5332 numerator
    // and <4304 denominator bits (<9636 total). Two column squares share
    // det(M)^2: <6459 numerator + <6458 denominator bits (<12917 total).
    // BigRational::new/arithmetic normalize, so cancellation can only reduce
    // these common-denominator bounds; both are below bounded()'s16384 cap.
    // Faces retain original binary64 source centres, never arbitrary Qs.
    let center =
        q(f64::from(a.outward))? * (coordinate(m, k, bface.center.clone())? - &a.center[k]);
    let square =
        &m.columns[rb[0]][k] * &m.columns[rb[0]][k] + &m.columns[rb[1]][k] * &m.columns[rb[1]][k];
    let gap_upper = if center.is_zero() && square.is_zero() {
        0.
    } else {
        finite(
            (passage::directed(&center)? + ip(bface.outer) * passage::directed(&square)?.sqrt())
                * b.metric.maximum_scale
                / positive(passage::directed(&m.columns[j][k].abs())?)?,
        )?
        .hi
    };
    if gap_upper > req.maximum_bearing_gap_mm {
        return Err(Error::Unsupported);
    }
    let inverse = map(bface.frame, a.frame, interrupt)?;
    let ca = coordinate(&inverse, rb[0], a.center.clone())? - &bface.center[rb[0]];
    let cb = coordinate(&inverse, rb[1], a.center.clone())? - &bface.center[rb[1]];
    let eccentric = bounded(ca.abs() + cb.abs())?;
    let ra = [(k + 1) % 3, (k + 2) % 3];
    let gram = |u: usize, v: usize| {
        bounded(rb.iter().fold(Q::zero(), |s, &i| {
            s + &inverse.columns[u][i] * &inverse.columns[v][i]
        }))
    };
    let g00 = gram(ra[0], ra[0])?;
    let g11 = gram(ra[1], ra[1])?;
    let g01 = gram(ra[0], ra[1])?.abs();
    let low = bounded((&g00 - &g01).min(&g11 - &g01))?;
    let high = bounded((g00 + &g01).max(g11 + g01))?;
    let lambda_low = positive(passage::directed(&low)?.sqrt())?.lo;
    let lambda_high = positive(passage::directed(&high)?.sqrt())?.hi;
    let e = passage::directed(&eccentric)?;
    let inner = a
        .inner
        .max(finite((ip(bface.inner) + e) / ip(lambda_low))?.hi);
    let outer = a
        .outer
        .min(finite((ip(bface.outer) - e) / ip(lambda_high))?.lo);
    let width = bounded(q(outer)? - q(inner)?)?;
    if width <= Q::zero()
        || bounded(&width * q(a.minimum_scale.lo)?)? < q(req.minimum_bearing_land_mm)?
    {
        return Err(Error::Unsupported);
    }
    let physical = finite(passage::directed(&width)? * a.minimum_scale)?;
    Ok(BearingContact {
        source_annular_radii_mm: [inner, outer],
        radial_land_lower_mm: physical.lo.max(req.minimum_bearing_land_mm),
        correspondence_gap_upper_mm: gap_upper,
    })
}
/// All related material pairs are proved here. Unrelated scene components still
/// require their own source-owned separation proof in the caller.
pub fn certify_hanger_stack(
    input: StackSources<'_>,
    req: StackRequirements,
    interrupt: &dyn Fn() -> bool,
) -> R<HangerStack> {
    checkpoint(interrupt)?;
    if !req.minimum_bearing_land_mm.is_finite()
        || req.minimum_bearing_land_mm <= 0.
        || req.minimum_bearing_land_mm > 1e7
        || !req.maximum_bearing_gap_mm.is_finite()
        || req.maximum_bearing_gap_mm < 0.
        || req.maximum_bearing_gap_mm > 1e7
    {
        return Err(Error::Invalid);
    }
    if input.washers.len() > 62 {
        return Err(Error::Limit);
    }
    let wood = certify_hanger_wood(
        Stock {
            material: input.wood.material,
            pilot: input.wood.pilot,
            frame: input.wood.frame,
        },
        input.hanger,
        input.hanger_frame,
        req.wood,
        interrupt,
    )?;
    let clamp_passage = certify_hanger_through_stock(
        input.clamp,
        input.hanger,
        input.hanger_frame,
        req.minimum_clamp_wall_mm,
        req.passage,
        interrupt,
    )?;
    let nut = certify_hanger_nut(
        input.hanger,
        input.hanger_frame,
        input.nut,
        input.nut_frame,
        req.mating,
        interrupt,
    )?;
    let wood_body = Placed::new(Material::Cells(input.wood.material), input.wood.frame)?;
    let clamp_body = Placed::new(Material::Cells(input.clamp.material), input.clamp.frame)?;
    let boundary =
        orthogonal_boundary::build(input.wood.material.quads()).map_err(boundary_error)?;
    let mouth = boundary
        .faces
        .get(input.wood.pilot.face)
        .ok_or(Error::Invalid)?;
    let (axis, inward) = inward_axis(mouth.normal)?;
    // The positive original wood-material witness already earns pilotR<woodR.
    let wood_face = cap(
        &wood_body,
        axis,
        vector(input.wood.pilot.center_mm)?,
        -(inward as i8),
        [
            input.wood.pilot.radius_mm,
            input.hanger.wood().maximum_radius_mm(),
        ],
    );
    let bounds =
        orthogonal_boundary::build(input.clamp.material.quads()).map_err(boundary_error)?;
    let face = bounds
        .faces
        .get(input.clamp.bore.faces[0])
        .ok_or(Error::Invalid)?;
    let (clamp_axis, _) = inward_axis(face.normal)?;
    let [p0, p1] = input.clamp.bore.centers_mm;
    let (start, end) = if p0[clamp_axis] < p1[clamp_axis] {
        (p0, p1)
    } else {
        (p1, p0)
    };
    // Consume the actual material band earned by the existing source proof,
    // rather than repeating protected-radius arithmetic and validation.
    let clamp_land = clamp_passage.source_cap_annular_land_mm;
    let wood_end = clamp_passage.wood_side == WoodSide::AfterEnd;
    let clamp_near = cap(
        &clamp_body,
        clamp_axis,
        vector(if wood_end { end } else { start })?,
        if wood_end { 1 } else { -1 },
        clamp_land,
    );
    // Same radii as the checked near cap; do not validate that predicate twice.
    let clamp_far = cap(
        &clamp_body,
        clamp_axis,
        vector(if wood_end { start } else { end })?,
        if wood_end { -1 } else { 1 },
        clamp_land,
    );
    let projection = outside(&wood_face, &clamp_body, interrupt)?;
    let mut bearings = Vec::with_capacity(input.washers.len() + 2);
    bearings.push(bearing(
        &wood_face,
        &clamp_near,
        &clamp_body,
        &projection,
        req,
        interrupt,
    )?);
    let mut planes = Vec::with_capacity(input.washers.len() + 2);
    planes.push(wood_face);
    planes.push(clamp_far);
    let mut washer_passages = Vec::with_capacity(input.washers.len());
    for placed in input.washers {
        checkpoint(interrupt)?;
        washer_passages.push(certify_hanger_annulus(
            placed.source,
            placed.frame,
            input.hanger,
            input.hanger_frame,
            req.passage,
            interrupt,
        )?);
        let b = Placed::new(Material::Ring(placed.source), placed.frame)?;
        for a in &planes[..planes.len() - 1] {
            outside(a, &b, interrupt)?;
        }
        let a = planes.last().unwrap();
        let projection = outside(a, &b, interrupt)?;
        let near = member_near_cap(a, &b, &projection, interrupt)?;
        bearings.push(bearing(a, &near, &b, &projection, req, interrupt)?);
        planes.push(ring_cap(&b, placed.source, near.outward < 0)?);
    }
    let b = Placed::new(Material::Nut(input.nut), input.nut_frame)?;
    for a in &planes[..planes.len() - 1] {
        outside(a, &b, interrupt)?;
    }
    let a = planes.last().unwrap();
    let projection = outside(a, &b, interrupt)?;
    let near = member_near_cap(a, &b, &projection, interrupt)?;
    bearings.push(bearing(a, &near, &b, &projection, req, interrupt)?);
    checkpoint(interrupt)?;
    Ok(HangerStack {
        wood,
        nut,
        clamp: clamp_passage,
        washer_passages,
        bearings,
    })
}
