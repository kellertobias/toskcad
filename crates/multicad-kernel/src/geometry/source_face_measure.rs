//! Actual affine surface measures, with bounded area and centroid reweighting.
use super::*;
use crate::geometry::arithmetic::{dot_error, sum_error};
use crate::geometry::internal_thread_treatment_body::interval::Interval as B;
use multicad_rust_kernel::authored_internal_treatment::FacePropertyOptions;
#[derive(Clone, Copy)]
pub(super) struct Measures {
    pub area_mm2: B,
    pub first_moment_mm3: [B; 3],
    pub centroid_mm: [B; 3],
}
pub(super) struct Receipt<P> {
    pub local: P,
    pub world: Measures,
    pub area_jacobian: B,
    pub local_centroid_reweighting_error_mm: [f64; 3],
}
pub(super) struct Request {
    pub frame: PlaneFrame,
    pub options: FacePropertyOptions,
    pub planar: bool,
    pub domain: [B; 3],
    pub whole_domain: [B; 3],
    pub area_upper_mm2: f64,
}
fn refused() -> GeometryError {
    error(
        "Rust source face measure",
        GeometryErrorKind::Unsupported,
        "Actual world face measure exceeds its precision bound",
    )
}
#[derive(Clone, Copy)]
pub(super) struct Placement {
    pub origin: [f64; 3],
    pub axes: [[f64; 3]; 3],
}
impl From<PlaneFrame> for Placement {
    fn from(frame: PlaneFrame) -> Self {
        Self {
            origin: frame.origin,
            axes: [frame.u, frame.v, frame.normal()],
        }
    }
}
pub(super) fn cardinal(axes: [[f64; 3]; 3]) -> bool {
    let mut occupied = 0;
    for axis in axes {
        let mut row = None;
        for (i, v) in axis.into_iter().enumerate() {
            if v == 1. || v == -1. {
                if row.replace(i).is_some() {
                    return false;
                }
            } else if v != 0. {
                return false;
            }
        }
        let Some(row) = row else { return false };
        if occupied & (1 << row) != 0 {
            return false;
        }
        occupied |= 1 << row;
    }
    true
}
pub(super) fn positive(value: B) -> Result<B, GeometryError> {
    let value = value.checked().ok_or_else(refused)?;
    if value.lo <= 0. {
        Err(refused())
    } else {
        Ok(value)
    }
}
pub(super) fn metric(placement: Placement, planar: bool) -> Result<(B, bool), GeometryError> {
    let axes = placement.axes;
    if axes.into_iter().flatten().any(|v| !v.is_finite()) {
        return Err(input("Face placement axes must be finite"));
    }
    let mut e = 0.;
    for (i, a) in axes.iter().enumerate() {
        for (j, b) in axes.iter().enumerate() {
            let (v, err) = dot_error(*a, *b);
            let (v, res) = sum_error(v, if i == j { -1. } else { 0. });
            e = up(e, up(v.abs(), up(err, res)));
        }
    }
    let pairs = [(1, 2), (2, 0), (0, 1)];
    let values =
        pairs.map(|(j, k)| dot_error([axes[0][j], axes[0][k], 0.], [axes[1][k], -axes[1][j], 0.]));
    let cross = values.map(|(value, _)| value);
    let errors = values.map(|(_, error)| error);
    let (d, mut de) = dot_error(cross, axes[2]);
    for i in 0..3 {
        de = up(de, product(errors[i], axes[2][i].abs()));
    }
    if !e.is_finite() || e > 1e-6 || (d - de).next_down() <= 0. {
        return Err(refused());
    }
    if cardinal(axes) {
        return Ok((B::ONE, true));
    }
    if planar {
        let normal = std::array::from_fn::<_, 3, _>(|i| {
            B::new(
                (cross[i] - errors[i]).next_down(),
                (cross[i] + errors[i]).next_up(),
            )
        });
        // The Gram bound above puts both columns near unit length and their
        // cross product away from zero, also bounding every intermediate.
        let jac = normal
            .into_iter()
            .fold(B::ZERO, |sum, v| sum + v.square())
            .sqrt();
        return Ok((jac, true));
    }
    // If ||A^T A-I||_2<=e, every product of two singular values lies in
    // [1-e,1+e]. Hence |cof(A)n| has this bound for every actual unit normal.
    Ok((B::new((1. - e).next_down(), (1. + e).next_up()), false))
}
pub(super) fn world_point(placement: Placement, p: [B; 3]) -> [B; 3] {
    std::array::from_fn(|i| {
        B::point(placement.origin[i])
            + (0..3).fold(B::ZERO, |sum, k| {
                sum + B::point(placement.axes[k][i]) * p[k]
            })
    })
}
pub(super) fn local_options(
    options: FacePropertyOptions,
    jac: B,
    area_upper: f64,
    row: f64,
    world_extent: f64,
) -> Result<FacePropertyOptions, GeometryError> {
    let a = B::point(options.area_tolerance_mm2) / (B::point(4.) * B::point(jac.hi));
    let translation = B::point(options.first_moment_tolerance_mm3)
        / (B::point(4.) * B::point(jac.hi) * B::point(world_extent + 1.));
    let c = B::point(options.centroid_tolerance_mm) / (B::point(4.) * B::point(row));
    let moment_c = B::point(options.first_moment_tolerance_mm3)
        / (B::point(4.) * B::point(jac.hi * area_upper + 1.) * B::point(row));
    let moment = B::point(options.first_moment_tolerance_mm3)
        / (B::point(4.) * B::point(jac.hi) * B::point(row));
    let result = FacePropertyOptions {
        area_tolerance_mm2: a.lo.min(translation.lo),
        first_moment_tolerance_mm3: moment.lo,
        centroid_tolerance_mm: c.lo.min(moment_c.lo),
        ..options
    };
    if [
        result.area_tolerance_mm2,
        result.first_moment_tolerance_mm3,
        result.centroid_tolerance_mm,
    ]
    .into_iter()
    .any(|v| !v.is_finite() || v <= 0.)
    {
        Err(refused())
    } else {
        Ok(result)
    }
}

pub(super) fn transport<P>(
    request: Request,
    policy: GeometryPolicy,
    query: impl FnOnce(FacePropertyOptions) -> Result<(P, Measures), GeometryError>,
) -> Result<Receipt<P>, GeometryError> {
    let Request {
        frame,
        options,
        planar,
        domain,
        whole_domain,
        area_upper_mm2: upper,
    } = request;
    policy.validate()?;
    if options.maximum_cells == 0
        || options.maximum_cells > 65536
        || [
            options.area_tolerance_mm2,
            options.first_moment_tolerance_mm3,
            options.centroid_tolerance_mm,
        ]
        .into_iter()
        .any(|v| !v.is_finite() || v <= 0.)
    {
        return Err(input("Invalid world face-property accuracy or cell budget"));
    }
    policy.point(frame.origin)?;
    let placement = Placement::from(frame);
    let (jac, constant) = metric(placement, planar)?;
    let world_bounds = world_point(placement, whole_domain);
    for bounds in world_bounds {
        bounds.checked().ok_or_else(refused)?;
    }
    policy.point(world_bounds.map(|b| b.lo))?;
    policy.point(world_bounds.map(|b| b.hi))?;
    let row = (0..3)
        .map(|i| (0..3).fold(0., |sum, k| up(sum, placement.axes[k][i].abs())))
        .fold(1., f64::max);
    let extent = world_bounds
        .into_iter()
        .map(B::magnitude)
        .fold(0., f64::max);
    let (local, values) = query(local_options(options, jac, upper, row, extent)?)?;
    let area = values.area_mm2;
    let world_area = positive(jac * area)?;
    let center = jac.lo + (jac.hi - jac.lo) * 0.5;
    let delta = if constant {
        0.
    } else {
        (B::point(center) - jac).magnitude().next_up()
    };
    let local_c = values.centroid_mm;
    let mut reweight = [0.; 3];
    let mut weighted = local_c;
    for k in 0..3 {
        let distance = (domain[k] - local_c[k]).magnitude();
        reweight[k] = (B::point(delta) / B::point(jac.lo) * B::point(distance)).hi;
        weighted[k] = (local_c[k] + B::new(-reweight[k], reweight[k]))
            .intersect(domain[k])
            .ok_or_else(refused)?;
    }
    // The area factor varies with the actual normal. Its covariance with x-C
    // is bounded by delta/Jlo times the source-coordinate distance, since
    // integral(x-C)dS=0 and integral(J)dS>=Jlo*integral(dS).
    let centroid = world_point(placement, weighted);
    let weighted_moment = std::array::from_fn::<_, 3, _>(|k| {
        if constant {
            values.first_moment_mm3[k] * jac
        } else {
            B::point(center) * values.first_moment_mm3[k]
                + B::new(-1., 1.) * B::point(delta) * area * B::point(domain[k].magnitude())
        }
    });
    let direct = std::array::from_fn::<_, 3, _>(|i| {
        B::point(placement.origin[i]) * world_area
            + (0..3).fold(B::ZERO, |sum, k| {
                sum + B::point(placement.axes[k][i]) * weighted_moment[k]
            })
    });
    let mut moments = [B::ZERO; 3];
    for i in 0..3 {
        moments[i] = (centroid[i] * world_area)
            .intersect(direct[i])
            .ok_or_else(refused)?;
    }
    let values = centroid.into_iter().chain(moments).chain([world_area]);
    if !values.into_iter().all(B::valid)
        || world_area.width() > options.area_tolerance_mm2
        || moments
            .iter()
            .any(|b| b.width() > options.first_moment_tolerance_mm3)
        || centroid
            .iter()
            .any(|b| b.width() > options.centroid_tolerance_mm)
    {
        return Err(refused());
    }
    Ok(Receipt {
        local,
        world: Measures {
            area_mm2: world_area,
            first_moment_mm3: moments,
            centroid_mm: centroid,
        },
        area_jacobian: jac,
        local_centroid_reweighting_error_mm: reweight,
    })
}
fn up(a: f64, b: f64) -> f64 {
    (a + b).next_up()
}
fn product(a: f64, b: f64) -> f64 {
    (a * b).next_up()
}

pub(super) fn midpoint(values: [B; 3]) -> Result<([f64; 3], f64), GeometryError> {
    if values.into_iter().any(|v| !v.valid()) {
        return Err(refused());
    }
    let point = values.map(|v| v.lo.midpoint(v.hi));
    let e = values
        .into_iter()
        .zip(point)
        .fold(0., |s, (b, p)| up(s, (B::point(p) - b).magnitude()));
    Ok((point, e))
}
pub(super) fn center(values: [B; 3], policy: GeometryPolicy) -> Result<[f64; 3], GeometryError> {
    let (p, e) = midpoint(values)?;
    if e > policy.linear_tolerance_mm / 8. {
        return Err(refused());
    }
    Ok(p)
}

pub(super) type PlaneSupport = ([f64; 3], [f64; 3]);
pub(super) fn plane(
    frame: PlaneFrame,
    z: B,
    radius: f64,
    up: bool,
    policy: GeometryPolicy,
) -> Result<PlaneSupport, GeometryError> {
    let u = frame.u.map(B::point);
    let v = frame.v.map(B::point);
    let cross = [(1, 2), (2, 0), (0, 1)].map(|(j, k)| u[j] * v[k] - u[k] * v[j]);
    let length = cross
        .into_iter()
        .fold(B::ZERO, |s, v| s + v.square())
        .sqrt();
    if !length.valid() || length.lo <= 0. {
        return Err(refused());
    }
    let orientation = if up { 1. } else { -1. };
    let unit = cross.map(|v| B::point(orientation) * v / length);
    let (normal, normal_error) = midpoint(unit)?;
    let point =
        std::array::from_fn(|i| B::point(frame.origin[i]) + B::point(frame.normal()[i]) * z);
    let (point, point_error) = midpoint(point)?;
    let extent = B::point(radius)
        * frame
            .u
            .into_iter()
            .chain(frame.v)
            .fold(B::ZERO, |s, v| s + B::point(v.abs()));
    let error =
        B::point(point_error) * (B::ONE + B::point(normal_error)) + B::point(normal_error) * extent;
    if normal_error > 1e-12 || !error.valid() || error.hi > policy.linear_tolerance_mm / 8. {
        return Err(refused());
    }
    policy.point(point)?;
    Ok((point, normal))
}
