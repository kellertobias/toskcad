//! Measures of the actual affine image, not an assumed Euclidean isometry.
use super::interval::Interval as B;
use super::*;
use source::{FaceKey, FaceProperties, FacePropertyOptions};
#[cfg(test)]
mod tests;
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InternalThreadFaceProperties {
    pub local: FaceProperties,
    pub world: FaceProperties,
    pub area_jacobian: source::Enclosure,
    /// Bound on reweighting of each local centroid coordinate. It is zero
    /// only when the actual map has a source-proved constant area Jacobian.
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
    unsupported("Actual world face measure cannot be certified at the requested precision")
}
impl Solid {
    /// Typed source handoff only; this does not enable general face/topology APIs.
    pub fn rust_internal_thread_treatment_world_face_properties_with_interrupt(
        &self,
        face: FaceKey,
        options: FacePropertyOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Option<InternalThreadFaceProperties>, GeometryError> {
        checkpoint(interrupt)?;
        self.internal_thread_treatment_body()
            .map(|(body, frame)| {
                body.world_face_properties(frame, face, options, self.policy, interrupt)
            })
            .transpose()
    }
}
impl InternalThreadTreatmentBody {
    pub(super) fn world_face_properties(
        &self,
        frame: PlaneFrame,
        face: FaceKey,
        options: FacePropertyOptions,
        policy: GeometryPolicy,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<InternalThreadFaceProperties, GeometryError> {
        checkpoint(interrupt)?;
        self.validate_frame(frame, policy)?;
        transport(
            &self.material,
            Placement::from(frame),
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
    matches!(
        face,
        FaceKey::StartCap | FaceKey::EndCap | FaceKey::TerminalShoulder
    )
}
fn coordinates(material: &Material, face: FaceKey) -> [B; 3] {
    let radius = B::new(radius_bounds(material)[0], radius_bounds(material)[1]);
    let c = material.cavity();
    let z = match face {
        FaceKey::StartCap => B::ZERO,
        FaceKey::EndCap => B::point(c.full_length_mm()),
        FaceKey::TerminalShoulder => B::point(c.terminal_station_mm().expect("validated shoulder")),
        FaceKey::TerminalCylinder => B::new(
            c.terminal_station_mm().expect("validated cylinder"),
            c.full_length_mm(),
        ),
        _ => B::new(0., c.full_length_mm()),
    };
    [
        B::new(-radius.hi, radius.hi),
        B::new(-radius.hi, radius.hi),
        z,
    ]
}
fn world_point(placement: Placement, p: [B; 3]) -> [B; 3] {
    std::array::from_fn(|i| {
        B::point(placement.origin[i])
            + (0..3).fold(B::ZERO, |sum, k| {
                sum + B::point(placement.axes[k][i]) * p[k]
            })
    })
}
fn area_upper(material: &Material, face: FaceKey) -> Result<f64, GeometryError> {
    let c = material.cavity();
    let r = B::point(radius_upper(material));
    let l = B::point(c.full_length_mm());
    let pi = B::new(
        std::f64::consts::PI.next_down(),
        std::f64::consts::PI.next_up(),
    );
    let tau = B::point(2.) * pi;
    let area = match face {
        FaceKey::MaterialOuterWall => tau * r * l,
        FaceKey::StartCap | FaceKey::EndCap | FaceKey::TerminalShoulder => pi * r.square(),
        FaceKey::TerminalCylinder => {
            tau * r * (l - B::point(c.terminal_station_mm().ok_or_else(refused)?))
        }
        FaceKey::EntryCone => {
            let a = B::point(c.entry_radius_mm().ok_or_else(refused)?);
            let depth = B::point(c.intent().entry_depth_mm.ok_or_else(refused)?);
            let k = (a - B::point(c.source_thread().minimum_radius_mm())) / depth;
            tau * a * depth * (B::ONE + k.square()).sqrt()
        }
        FaceKey::Thread { .. } => {
            let spec = c.source_thread().spec();
            let lead = B::point(spec.pitch_mm) * B::point(f64::from(spec.starts));
            let gain = tau / lead;
            let mut sum = B::ZERO;
            for period in 0..spec.starts {
                for (band, b) in c.source_thread().bands().into_iter().enumerate() {
                    if c.thread_face_key(i64::from(period), band)
                        .map_err(failure)?
                        != face
                    {
                        continue;
                    }
                    let dq = B::point(b.phase_mm[1]) - B::point(b.phase_mm[0]);
                    let m = (B::point(b.radius_mm[1]) - B::point(b.radius_mm[0])) / dq;
                    let radius = B::point(b.radius_mm[0].max(b.radius_mm[1]));
                    sum = sum
                        + l * dq
                            * ((radius * gain).square() * (B::ONE + m.square()) + m.square())
                                .sqrt();
                }
            }
            sum
        }
    };
    Ok(positive(area)?.hi)
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
    material: &Material,
    placement: Placement,
    face: FaceKey,
    options: FacePropertyOptions,
    policy: GeometryPolicy,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<InternalThreadFaceProperties, GeometryError> {
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
    let world_bounds = world_point(placement, coordinates(material, FaceKey::MaterialOuterWall));
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
        .local_face_properties_with_interrupt(
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
    let result = InternalThreadFaceProperties {
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
