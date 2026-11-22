//! Author actual numerical seating poses, never a geometry acceptance waiver.
//! Final Native material/phase/floor proofs remain mandatory. Directed finite
//! projections measure the computed poses, including their actual trig outputs.
use super::*;
mod insert;
pub(super) use insert::pose as insert_pose;
#[cfg(test)]
mod tests;
#[derive(Clone, Copy)]
struct I {
    lo: f64,
    hi: f64,
}
impl I {
    fn p(x: f64) -> Self {
        Self { lo: x, hi: x }
    }
    fn neg(self) -> Self {
        Self {
            lo: -self.hi,
            hi: -self.lo,
        }
    }
    fn add(self, b: Self) -> Self {
        if self.lo == 0. && self.hi == 0. {
            return b;
        }
        if b.lo == 0. && b.hi == 0. {
            return self;
        }
        Self {
            lo: (self.lo + b.lo).next_down(),
            hi: (self.hi + b.hi).next_up(),
        }
    }
    fn sub(self, b: Self) -> Self {
        if self.lo == self.hi && self.lo == b.lo && b.lo == b.hi {
            Self::p(0.)
        } else {
            self.add(b.neg())
        }
    }
    fn mul(self, b: Self) -> Self {
        if self.lo == 0. && self.hi == 0. || b.lo == 0. && b.hi == 0. {
            return Self::p(0.);
        }
        if self.lo == 1. && self.hi == 1. {
            return b;
        }
        if b.lo == 1. && b.hi == 1. {
            return self;
        }
        if self.lo == -1. && self.hi == -1. {
            return b.neg();
        }
        if b.lo == -1. && b.hi == -1. {
            return self.neg();
        }
        let p = [
            self.lo * b.lo,
            self.lo * b.hi,
            self.hi * b.lo,
            self.hi * b.hi,
        ];
        Self {
            lo: p.into_iter().fold(f64::INFINITY, f64::min).next_down(),
            hi: p.into_iter().fold(f64::NEG_INFINITY, f64::max).next_up(),
        }
    }
    /// Only the squared length of the validated axis column is used here.
    /// RigidPose::validate earns its interval lower bound above .99.
    fn axis_div(self, b: Self) -> Self {
        if b.lo == 1. && b.hi == 1. {
            return self;
        }
        self.mul(Self {
            lo: b.hi.recip().next_down(),
            hi: b.lo.recip().next_up(),
        })
    }
    fn abs(self) -> f64 {
        self.lo.abs().max(self.hi.abs())
    }
    fn midpoint(self) -> f64 {
        self.lo + (self.hi - self.lo) * 0.5
    }
}
fn dot(a: [f64; 3], b: [I; 3]) -> I {
    (0..3).fold(I::p(0.), |s, i| s.add(I::p(a[i]).mul(b[i])))
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
#[derive(Clone, Copy)]
pub(super) struct Options {
    maximum_gap_mm: f64,
    maximum_attempts: usize,
}
impl Options {
    pub(super) fn canonical() -> Self {
        // Minimum supported material-verification precision. This is a real
        // maximum gap contract, never permission to penetrate a component.
        Self {
            maximum_gap_mm: 1e-7,
            maximum_attempts: 4,
        }
    }
}
fn extent(part: &ResolvedMechanicalPart) -> Result<(f64, f64), DocumentError> {
    match part {
        ResolvedMechanicalPart::Washer {
            outside_radius_mm,
            thickness_mm,
            ..
        } => Ok((*outside_radius_mm, *thickness_mm)),
        ResolvedMechanicalPart::HexNut {
            across_flats_mm,
            thread,
        } => Ok((*across_flats_mm, thread.length_mm)),
        ResolvedMechanicalPart::CustomConnection { parameters } => {
            if let model::ResolvedCustomConnectionHardware::HangerBolt {
                wood_thread,
                plain_length_mm,
                plain_radius_mm,
                metric_thread,
                ..
            } = parameters.as_ref()
            {
                Ok((
                    wood_thread
                        .major_radius_mm
                        .max(*plain_radius_mm)
                        .max(metric_thread.major_radius_mm),
                    I::p(wood_thread.length_mm)
                        .add(I::p(*plain_length_mm))
                        .add(I::p(metric_thread.length_mm))
                        .hi,
                ))
            } else {
                Err(invalid("Hanger seating requires an actual hanger source"))
            }
        }
        _ => Err(invalid(
            "Hanger seating requires actual hanger, nut and washer sources",
        )),
    }
}
fn residual(
    axis: RigidPose,
    pose: RigidPose,
    mount: RigidPose,
    target: I,
    radius: f64,
    height: f64,
) -> f64 {
    let direction = axis.vector([0., 0., 1.]);
    let denominator = dot(direction, direction.map(I::p));
    // These are the literal published source-frame operations; no column is
    // normalized. Its third column is the literal cross of published U,V.
    let origin = pose.point(mount.translation_mm);
    let u = pose.vector(mount.vector([1., 0., 0.]));
    let v = pose.vector(mount.vector([0., 1., 0.]));
    let n = cross(u, v);
    let offset = dot(
        direction,
        std::array::from_fn(|i| I::p(origin[i]).sub(I::p(axis.translation_mm[i]))),
    )
    .axis_div(denominator)
    .sub(target)
    .abs();
    let radial = I::p(dot(direction, u.map(I::p)).axis_div(denominator).abs())
        .add(I::p(
            dot(direction, v.map(I::p)).axis_div(denominator).abs(),
        ))
        .mul(I::p(radius));
    let axial = I::p(
        dot(direction, n.map(I::p))
            .axis_div(denominator)
            .sub(I::p(1.))
            .abs(),
    )
    .mul(I::p(height));
    I::p(offset).add(radial).add(axial).hi
}
fn datum(
    rows: &[model::ScrewConnectionHardwareRequest],
    parts: &[ResolvedMechanicalPart],
    station: f64,
) -> Result<(f64, f64, model::ThreadHand), DocumentError> {
    let sources = rows
        .iter()
        .zip(parts)
        .filter(|(r, _)| r.role == Role::HangerBolt)
        .map(|(_, p)| p)
        .collect::<Vec<_>>();
    let [ResolvedMechanicalPart::CustomConnection { parameters }] = sources.as_slice() else {
        return Err(invalid(
            "Hanger seating requires exactly one actual hanger source",
        ));
    };
    let model::ResolvedCustomConnectionHardware::HangerBolt {
        plain_length_mm,
        metric_thread,
        ..
    } = parameters.as_ref()
    else {
        return Err(invalid(
            "Hanger seating requires an actual continuous metric datum",
        ));
    };
    Ok((
        station - plain_length_mm - metric_thread.length_mm,
        metric_thread.lead_mm,
        metric_thread.hand,
    ))
}
/// Only actual source pose authoring occurs here. In particular a positive
/// wood shift is checked later against original pilot depth/floor/engagement;
/// it cannot silently repair or override a failed material certificate.
pub(super) fn pose(
    rows: &[model::ScrewConnectionHardwareRequest],
    parts: &[ResolvedMechanicalPart],
    role: EntityId,
    grip: f64,
    axis: RigidPose,
    p: &EvaluatedParameters,
    options: Options,
) -> Result<(f64, RigidPose), DocumentError> {
    axis.validate()?;
    if parts.len() != rows.len()
        || parts.len() > 64
        || !options.maximum_gap_mm.is_finite()
        || options.maximum_gap_mm <= 0.
        || options.maximum_attempts == 0
        || options.maximum_attempts > 4
    {
        return Err(invalid(
            "Hanger numerical seating inputs exceed their contract",
        ));
    }
    let limits = parts.iter().map(extent).collect::<Result<Vec<_>, _>>()?;
    let radius = limits.iter().map(|x| x.0).fold(0., f64::max);
    let height = limits.iter().map(|x| x.1).fold(0., f64::max);
    let mut tails = vec![I::p(0.); parts.len() + 1];
    let mut count = 0;
    for i in (0..parts.len()).rev() {
        tails[i] = tails[i + 1];
        if rows[i].role == Role::NutWasher {
            tails[i] = tails[i].add(I::p(washer(&parts[i])?));
            count += 1;
        }
    }
    let mut gap = 0.;
    let mut wood_shift = 0.;
    for _ in 0..options.maximum_attempts {
        let hanger_target = I::p(grip).add(I::p(wood_shift));
        let hanger_station = hanger_target.midpoint();
        let reference = datum(rows, parts, hanger_station)?;
        let mut poses = Vec::with_capacity(parts.len());
        let mut left = count;
        let mut error = 0_f64;
        let mut wood_error = 0_f64;
        for (i, (row, part)) in rows.iter().zip(parts).enumerate() {
            let (target, roll) = match row.role {
                Role::HangerBolt => (hanger_target, 0.),
                Role::NutWasher => {
                    let target = tails[i].neg().sub(I::p(gap).mul(I::p(f64::from(left))));
                    left -= 1;
                    (target, 0.)
                }
                Role::Nut => {
                    let target = tails[0]
                        .neg()
                        .sub(I::p(limits[i].1))
                        .sub(I::p(gap).mul(I::p(f64::from(count + 1))));
                    (
                        target,
                        signed_turn(reference.2, target.midpoint() - reference.0, reference.1)?,
                    )
                }
                _ => return Err(invalid("Hanger numerical seating has an incompatible role")),
            };
            let control = model::NumericTransform {
                translation: [length(0.), length(0.), length(target.midpoint())],
                angle: angle(roll.to_degrees()),
                ..Default::default()
            }
            .resolve(p)?;
            let mount = part.mount_pose();
            let pose = mount.inverse()?.then(control)?.then(axis)?;
            let r = residual(axis, pose, mount, target, radius, height);
            if row.role == Role::HangerBolt {
                wood_error = r;
            } else {
                error = error.max(r);
            }
            poses.push((row.role_id, target.midpoint(), pose));
        }
        if I::p(gap).add(I::p(2.).mul(I::p(error))).hi > options.maximum_gap_mm
            || I::p(wood_shift).add(I::p(wood_error)).hi > options.maximum_gap_mm
        {
            return Err(invalid(
                "The available contact tolerance cannot contain actual pose rounding",
            ));
        }
        let needed = I::p(4.).mul(I::p(error)).hi;
        let wood_needed = I::p(4.).mul(I::p(wood_error)).hi;
        if needed <= gap && wood_needed <= wood_shift {
            return poses
                .into_iter()
                .find(|x| x.0 == role)
                .map(|x| (x.1, x.2))
                .ok_or_else(|| invalid("Hanger seating role is absent from its connection"));
        }
        gap = I::p(gap.max(needed)).mul(I::p(1.25)).hi;
        wood_shift = I::p(wood_shift.max(wood_needed)).mul(I::p(1.25)).hi;
    }
    Err(invalid(
        "Hanger numerical seating exhausted its bounded pose attempts",
    ))
}
