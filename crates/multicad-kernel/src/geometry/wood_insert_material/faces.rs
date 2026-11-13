//! Measures and descriptors of the actual raw affine source image.
//! These APIs do not establish finite boundary incidence or Solid admission.
use super::*;
use crate::geometry::arithmetic::{dot_error, sum_error};
use crate::geometry::internal_thread_treatment_body::interval::Interval as B;
use multicad_rust_kernel::authored_internal_treatment::Enclosure;
use multicad_rust_kernel::wood_insert::{FaceKey, FaceProperties, FacePropertyOptions};
#[cfg(test)]
mod tests;
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WoodInsertFaceProperties {
    pub local: FaceProperties,
    pub world: FaceProperties,
    pub area_jacobian: Enclosure,
    pub local_centroid_reweighting_error_mm: [f64; 3],
}
#[derive(Clone, Copy)]
struct Placement {
    origin: [f64; 3],
    axes: [[f64; 3]; 3],
}
impl From<PlaneFrame> for Placement {
    fn from(frame: PlaneFrame) -> Self {
        Self {
            origin: frame.origin,
            axes: [frame.u, frame.v, frame.normal()],
        }
    }
}
fn refused() -> GeometryError {
    error(
        "Rust wood insert face measure",
        GeometryErrorKind::Unsupported,
        "Actual world face measure exceeds its precision bound",
    )
}
impl WoodInsertMaterial {
    /// Bounded world measures; the literal frame is never normalized.
    pub fn world_face_properties(
        &self,
        face: FaceKey,
        options: FacePropertyOptions,
        policy: GeometryPolicy,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<super::WoodInsertFaceProperties, GeometryError> {
        checkpoint(interrupt)?;
        policy.validate()?;
        policy.length(self.source.minimum_wall_mm())?;
        let r = self.source.outer().maximum_radius_mm();
        self.frame
            .checked_metric_extent([r, r, self.source.full_length_mm()], policy)?;
        transport(
            &self.source,
            Placement::from(self.frame),
            face,
            options,
            policy,
            interrupt,
        )
    }
}
fn cardinal(axes: [[f64; 3]; 3]) -> bool {
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
fn positive(value: B) -> Result<B, GeometryError> {
    let value = value.checked().ok_or_else(refused)?;
    if value.lo <= 0. {
        Err(refused())
    } else {
        Ok(value)
    }
}
fn metric(placement: Placement, planar: bool) -> Result<(B, bool), GeometryError> {
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
fn planar(face: FaceKey) -> bool {
    matches!(face, FaceKey::StartCap | FaceKey::EndCap)
}
fn coordinates(material: &Insert, face: FaceKey) -> [B; 3] {
    let r = material.outer().maximum_radius_mm();
    let z = match face {
        FaceKey::StartCap => B::ZERO,
        FaceKey::EndCap => B::point(material.full_length_mm()),
        _ => B::new(0., material.full_length_mm()),
    };
    [B::new(-r, r), B::new(-r, r), z]
}
fn world_point(placement: Placement, p: [B; 3]) -> [B; 3] {
    std::array::from_fn(|i| {
        B::point(placement.origin[i])
            + (0..3).fold(B::ZERO, |sum, k| {
                sum + B::point(placement.axes[k][i]) * p[k]
            })
    })
}
fn area_upper(material: &Insert, face: FaceKey) -> Result<f64, GeometryError> {
    let pi = B::new(
        std::f64::consts::PI.next_down(),
        std::f64::consts::PI.next_up(),
    );
    let thread = match face {
        FaceKey::StartCap | FaceKey::EndCap => {
            return Ok(positive(pi * B::point(material.outer().maximum_radius_mm()).square())?.hi);
        }
        FaceKey::OuterThread { .. } => material.outer(),
        FaceKey::InnerThread { .. } => material.inner(),
    };
    let spec = thread.spec();
    let gain = B::point(2.) * pi / (B::point(spec.pitch_mm) * B::point(f64::from(spec.starts)));
    // Complete actual wall is a conservative upper bound for any one role.
    // This bound supplies numerical budgets, never face ownership or geometry.
    let mut sum = B::ZERO;
    for b in thread.bands() {
        let dq = B::point(b.phase_mm[1]) - B::point(b.phase_mm[0]);
        let m = (B::point(b.radius_mm[1]) - B::point(b.radius_mm[0])) / dq;
        let r = B::point(b.radius_mm[0].max(b.radius_mm[1]));
        sum = sum
            + B::point(spec.length_mm)
                * dq
                * ((r * gain).square() * (B::ONE + m.square()) + m.square()).sqrt();
    }
    Ok(positive(sum * B::point(f64::from(spec.starts)))?.hi)
}
fn local_options(
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
fn transport(
    material: &Insert,
    placement: Placement,
    face: FaceKey,
    options: FacePropertyOptions,
    policy: GeometryPolicy,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<WoodInsertFaceProperties, GeometryError> {
    checkpoint(interrupt)?;
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
    if material.face_index(face).is_none() {
        return Err(input("Face does not belong to the actual material source"));
    }
    policy.point(placement.origin)?;
    let (jac, constant) = metric(placement, planar(face))?;
    let domain = coordinates(material, face);
    let world_bounds = world_point(
        placement,
        [
            B::new(
                -material.outer().maximum_radius_mm(),
                material.outer().maximum_radius_mm(),
            ),
            B::new(
                -material.outer().maximum_radius_mm(),
                material.outer().maximum_radius_mm(),
            ),
            B::new(0., material.full_length_mm()),
        ],
    );
    for bounds in world_bounds {
        bounds.checked().ok_or_else(refused)?;
    }
    policy.point(world_bounds.map(|b| b.lo))?;
    policy.point(world_bounds.map(|b| b.hi))?;
    let row = (0..3)
        .map(|i| (0..3).fold(0., |sum, k| up(sum, placement.axes[k][i].abs())))
        .fold(1., f64::max);
    let upper = area_upper(material, face)?;
    let extent = world_bounds
        .into_iter()
        .map(B::magnitude)
        .fold(0., f64::max);
    let local = material
        .local_face_properties(
            face,
            local_options(options, jac, upper, row, extent)?,
            interrupt,
        )
        .map_err(failure)?;
    let area = B::from(local.area_mm2);
    let world_area = positive(jac * area)?;
    let center = jac.lo + (jac.hi - jac.lo) * 0.5;
    let delta = if constant {
        0.
    } else {
        (B::point(center) - jac).magnitude().next_up()
    };
    let local_c = local.centroid_mm.map(B::from);
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
            B::from(local.first_moment_mm3[k]) * jac
        } else {
            B::point(center) * B::from(local.first_moment_mm3[k])
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
    let result = WoodInsertFaceProperties {
        local,
        world: FaceProperties {
            face,
            area_mm2: world_area.source(),
            first_moment_mm3: moments.map(B::source),
            centroid_mm: centroid.map(B::source),
        },
        area_jacobian: jac.source(),
        local_centroid_reweighting_error_mm: reweight,
    };
    checkpoint(interrupt)?;
    Ok(result)
}

fn up(a: f64, b: f64) -> f64 {
    (a + b).next_up()
}
fn product(a: f64, b: f64) -> f64 {
    (a * b).next_up()
}
fn midpoint(values: [B; 3]) -> Result<([f64; 3], f64), GeometryError> {
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
fn center(values: [B; 3], policy: GeometryPolicy) -> Result<[f64; 3], GeometryError> {
    let (p, e) = midpoint(values)?;
    if e > policy.linear_tolerance_mm / 8. {
        return Err(refused());
    }
    Ok(p)
}
type PlaneSupport = ([f64; 3], [f64; 3]);
fn plane(
    material: &Insert,
    frame: PlaneFrame,
    face: FaceKey,
    policy: GeometryPolicy,
) -> Result<Option<PlaneSupport>, GeometryError> {
    let z = match face {
        FaceKey::StartCap => 0.,
        FaceKey::EndCap => material.full_length_mm(),
        _ => return Ok(None),
    };
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
    let orientation = if face == FaceKey::StartCap { -1. } else { 1. };
    let unit = cross.map(|v| B::point(orientation) * v / length);
    let (normal, normal_error) = midpoint(unit)?;
    let point = std::array::from_fn(|i| {
        B::point(frame.origin[i]) + B::point(frame.normal()[i]) * B::point(z)
    });
    let (point, point_error) = midpoint(point)?;
    let extent = B::point(material.outer().maximum_radius_mm())
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
    Ok(Some((point, normal)))
}
impl WoodInsertMaterial {
    /// Source-role descriptors only; no finite boundary incidence is implied.
    pub fn faces(
        &self,
        policy: GeometryPolicy,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Vec<Face>, GeometryError> {
        checkpoint(interrupt)?;
        policy.validate()?;
        let r = self.source.outer().maximum_radius_mm();
        let l = self.source.full_length_mm();
        policy.length(self.source.minimum_wall_mm())?;
        self.frame.checked_metric_extent([r, r, l], policy)?;
        let scale = r.max(l).max(1.);
        let extent = self
            .frame
            .origin
            .into_iter()
            .map(f64::abs)
            .fold(1., f64::max)
            + 2. * (2. * r + l);
        let area_tolerance_mm2 = policy.linear_tolerance_mm * scale / 8.;
        let options = FacePropertyOptions {
            area_tolerance_mm2,
            centroid_tolerance_mm: policy.linear_tolerance_mm / 8.,
            first_moment_tolerance_mm3: area_tolerance_mm2 * (extent + scale),
            ..Default::default()
        };
        let mut result = Vec::with_capacity(self.source.face_count() as usize);
        for index in 0..self.source.face_count() {
            let key = self
                .source
                .face_key(index)
                .expect("sealed source face ordinal");
            let receipt = self.world_face_properties(key, options, policy, interrupt)?;
            let center_mm = center(receipt.world.centroid_mm.map(B::from), policy)?;
            let face = Face {
                index,
                plane: plane(&self.source, self.frame, key, policy)?,
                center_mm,
                area_mm2: receipt
                    .world
                    .area_mm2
                    .lower
                    .midpoint(receipt.world.area_mm2.upper),
            };
            policy.point(face.center_mm)?;
            result.push(face);
        }
        checkpoint(interrupt)?;
        Ok(result)
    }
}
