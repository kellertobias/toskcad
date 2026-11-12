//! Retained analytic meridians for countersunk and crowned screw heads.
//!
//! Coordinates are `[radius, axial position]`, with the nominal bearing at
//! zero and the head on negative Z. The final station is the plain-shank
//! junction, not an approximation of the helical thread.

use crate::helical_thread::ThreadError;
use std::f64::consts::{FRAC_PI_2, PI, TAU};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Role {
    Chamfer,
    Rim,
    Cone,
    Crown,
    Skirt,
    Bearing,
    Neck,
    Shank,
}

impl Role {
    pub fn source_role(self) -> &'static str {
        match self {
            Self::Chamfer => "top.chamfer",
            Self::Rim => "rim",
            Self::Cone => "cone",
            Self::Crown => "crown",
            Self::Skirt => "skirt",
            Self::Bearing => "bearing",
            Self::Neck => "neck.blend",
            Self::Shank => "shank",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Curve {
    Line,
    /// Literal support radius; trim vertices may differ only within the
    /// profile's exposed arithmetic reserve. The signed angle keeps direction.
    Arc {
        center: [f64; 2],
        radius_mm: f64,
        sweep: f64,
    },
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Segment {
    pub start: [f64; 2],
    pub end: [f64; 2],
    pub curve: Curve,
    pub role: Role,
}

impl Segment {
    pub fn point(self, fraction: f64) -> [f64; 2] {
        // Keep neighboring source faces on the very same endpoint values.
        if fraction == 0. {
            return self.start;
        }
        if fraction == 1. {
            return self.end;
        }
        match self.curve {
            Curve::Line => {
                std::array::from_fn(|i| self.start[i] + fraction * (self.end[i] - self.start[i]))
            }
            Curve::Arc {
                center,
                radius_mm: radius,
                sweep,
            } => {
                let angle = self.angle(center) + sweep * fraction;
                [
                    center[0] + radius * angle.cos(),
                    center[1] + radius * angle.sin(),
                ]
            }
        }
    }

    fn angle(self, center: [f64; 2]) -> f64 {
        (self.start[1] - center[1]).atan2(self.start[0] - center[0])
    }

    /// Analytic disk integral, without tessellation or numerical quadrature.
    pub fn volume_mm3(self) -> f64 {
        match self.curve {
            Curve::Line => {
                let [a, z] = self.start;
                let [b, end] = self.end;
                PI * (end - z) * (a * a + a * b + b * b) / 3.
            }
            Curve::Arc {
                center: [c, z],
                radius_mm: radius,
                sweep,
            } => {
                let start = self.angle([c, z]);
                if c == 0. {
                    let dz = 2. * radius * (start + sweep / 2.).cos() * (sweep / 2.).sin();
                    return PI
                        * dz
                        * ((radius * start.cos()).powi(2)
                            - radius * start.sin() * dz
                            - dz.powi(2) / 3.);
                }
                let integral = |t: f64| {
                    let s = t.sin();
                    radius * c * c * s
                        + c * radius * radius * (t + (2. * t).sin() / 2.)
                        + radius.powi(3) * (s - s.powi(3) / 3.)
                };
                PI * (integral(start + sweep) - integral(start))
            }
        }
    }

    /// Area of this revolved face, including horizontal bearing annuli.
    pub fn area_mm2(self) -> f64 {
        match self.curve {
            Curve::Line => {
                PI * (self.start[0] + self.end[0])
                    * (self.end[0] - self.start[0]).hypot(self.end[1] - self.start[1])
            }
            Curve::Arc {
                center,
                radius_mm: radius,
                sweep,
            } => {
                let angle = self.angle(center);
                if center[0] == 0. {
                    return TAU
                        * radius
                        * (2. * radius * (angle + sweep / 2.).cos() * (sweep / 2.).sin()).abs();
                }
                TAU * radius
                    * (center[0] * sweep + radius * ((angle + sweep).sin() - angle.sin())).abs()
            }
        }
    }

    /// Centroid of the complete revolved face, in the source meridian frame.
    pub fn center_mm(self) -> [f64; 3] {
        let z = match self.curve {
            Curve::Line => {
                let [r, z] = self.start;
                let [end_r, end_z] = self.end;
                z + (end_z - z) * (r + 2. * end_r) / (3. * (r + end_r))
            }
            Curve::Arc {
                center: [c, z],
                radius_mm: radius,
                sweep,
            } => {
                let a = self.angle([c, z]);
                if c == 0. {
                    let dz = 2. * radius * (a + sweep / 2.).cos() * (sweep / 2.).sin();
                    return [0., 0., z + radius * a.sin() + dz / 2.];
                }
                let b = a + sweep;
                let denominator = c * sweep + radius * (b.sin() - a.sin());
                z + (radius * c * (a.cos() - b.cos())
                    + radius * radius * (b.sin().powi(2) - a.sin().powi(2)) / 2.)
                    / denominator
            }
        };
        [0., 0., z]
    }

    /// Support of the revolved face; `radial` must be a nonnegative norm.
    pub fn support(self, radial: f64, axial: f64) -> f64 {
        let dot = |[r, z]: [f64; 2]| radial.mul_add(r, axial * z);
        let mut value = dot(self.start).max(dot(self.end));
        if let Curve::Arc {
            center,
            radius_mm: radius,
            sweep,
        } = self.curve
        {
            let start = self.angle(center);
            for angle in [start, start + sweep] {
                value = value.max(dot([
                    center[0] + radius * angle.cos(),
                    center[1] + radius * angle.sin(),
                ]));
            }
            let candidate = axial.atan2(radial);
            let advance = if sweep > 0. {
                (candidate - start).rem_euclid(TAU)
            } else {
                (start - candidate).rem_euclid(TAU)
            };
            if advance <= sweep.abs() {
                value = value.max(dot(center) + radius * radial.hypot(axial));
            }
        }
        value
    }

    fn radius_at(self, z: f64) -> f64 {
        if self.start[1] != self.end[1] {
            if z == self.start[1] {
                return self.start[0];
            }
            if z == self.end[1] {
                return self.end[0];
            }
        }
        match self.curve {
            Curve::Line => {
                if self.start[1] == self.end[1] {
                    self.start[0].max(self.end[0])
                } else {
                    self.start[0]
                        + (z - self.start[1]) / (self.end[1] - self.start[1])
                            * (self.end[0] - self.start[0])
                }
            }
            Curve::Arc {
                center,
                radius_mm: radius,
                ..
            } => {
                let offset = (z - center[1]).abs();
                let radial = ((radius - offset) * (radius + offset)).max(0.).sqrt();
                if self.role == Role::Crown {
                    center[0] + radial
                } else {
                    center[0] - radial
                }
            }
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Head {
    Countersunk {
        radius_mm: f64,
        height_mm: f64,
        included_angle_deg: f64,
        edge_chamfer_mm: f64,
        neck_radius_mm: f64,
    },
    /// Pan and button heads share this source geometry and differ in dimensions.
    Crowned {
        radius_mm: f64,
        height_mm: f64,
        crown_radius_mm: f64,
        top_land_radius_mm: f64,
        neck_radius_mm: f64,
    },
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Profile {
    segments: [Option<Segment>; 5],
    top: [f64; 2],
    end: [f64; 2],
    rounding_error_mm: f64,
}

impl Profile {
    /// Retain a complete physical meridian without reconstructing an authored
    /// head recipe. Shared stations must already have one source identity.
    /// The two radial caps and the axis close this monotone outer boundary.
    pub fn from_retained_segments(segments: &[Segment]) -> Result<Self, ThreadError> {
        if segments.is_empty() || segments.len() > 5 {
            return Err(ThreadError::Invalid);
        }
        let top = segments[0].start;
        let end = segments.last().expect("nonempty meridian").end;
        let mut scale = 0_f64;
        for segment in segments {
            let points = match segment.curve {
                Curve::Line => [segment.start, segment.end, [0.; 2]],
                Curve::Arc {
                    center, radius_mm, ..
                } => {
                    if !radius_mm.is_finite() || radius_mm <= 0. {
                        return Err(ThreadError::Invalid);
                    }
                    scale = scale.max(radius_mm);
                    [segment.start, segment.end, center]
                }
            };
            for point in points {
                if !point.into_iter().all(f64::is_finite) {
                    return Err(ThreadError::Invalid);
                }
                scale = point.into_iter().map(f64::abs).fold(scale, f64::max);
            }
        }
        if scale > 1e7 {
            return Err(ThreadError::Limit);
        }
        let reserve = 128. * f64::EPSILON * scale;
        if top[0] < 0.
            || end[0] <= reserve
            || top[1] >= 0.
            || end[1] < 0.
            || end[1] - top[1] <= reserve
            || (top[0] > 0. && top[0] <= reserve)
        {
            return Err(ThreadError::Invalid);
        }
        let mut result = Self {
            segments: [None; 5],
            top,
            end,
            rounding_error_mm: reserve,
        };
        for (i, &segment) in segments.iter().enumerate() {
            if (i > 0 && segments[i - 1].end != segment.start)
                || segments[..i].iter().any(|s| s.role == segment.role)
                || segment.start[0].min(segment.end[0]) < 0.
                || segment.end[0] <= reserve
                || !retained_segment_valid(segment, reserve)
            {
                return Err(ThreadError::Invalid);
            }
            let area = segment.area_mm2();
            let volume = segment.volume_mm3();
            let center = segment.center_mm();
            if !area.is_finite()
                || area <= 0.
                || !volume.is_finite()
                || (segment.role != Role::Bearing && volume <= 0.)
                || !center.into_iter().all(f64::is_finite)
                || center[2] < segment.start[1] - reserve
                || center[2] > segment.end[1] + reserve
            {
                return Err(ThreadError::Invalid);
            }
            result.segments[i] = Some(segment);
        }
        Ok(result)
    }

    pub fn new(head: Head, shank: f64, plain: f64) -> Result<Self, ThreadError> {
        Self::construct(head, shank, plain, false)
    }

    /// Retain the complete neck even when it overlaps the authored thread.
    /// The screw union factory must separately prove that hidden thread material
    /// is covered by this actual profile.
    pub fn new_for_thread_union(head: Head, shank: f64, plain: f64) -> Result<Self, ThreadError> {
        Self::construct(head, shank, plain, true)
    }

    fn construct(head: Head, shank: f64, plain: f64, union: bool) -> Result<Self, ThreadError> {
        let (radius, height, neck) = match head {
            Head::Countersunk {
                radius_mm,
                height_mm,
                neck_radius_mm,
                ..
            }
            | Head::Crowned {
                radius_mm,
                height_mm,
                neck_radius_mm,
                ..
            } => (radius_mm, height_mm, neck_radius_mm),
        };
        if ![radius, height, neck, shank, plain]
            .into_iter()
            .all(f64::is_finite)
            || shank <= 0.
            || radius <= shank
            || height <= 0.
            || neck < 0.
            || plain < 0.
        {
            return Err(ThreadError::Invalid);
        }
        if radius.max(height + plain).max(neck) > 1e7 {
            return Err(ThreadError::Limit);
        }
        let top_z = -height;
        // Reserve several dozen rounding steps at the complete source scale.
        // Stable factored square roots below avoid amplifying cancellation.
        // This budget alone is not a manufacturing certificate: the placed
        // facade must check source/trim residuals against its geometry policy.
        let mut reserve = 64. * f64::EPSILON * (radius + height + plain + neck);
        let mut result = Self {
            segments: [None; 5],
            top: [radius, top_z],
            end: [shank, plain],
            rounding_error_mm: reserve,
        };
        let mut count = 0;
        let mut push = |start, end, curve, role| {
            if start != end {
                result.segments[count] = Some(Segment {
                    start,
                    end,
                    curve,
                    role,
                });
                count += 1;
            }
        };
        match head {
            Head::Countersunk {
                included_angle_deg,
                edge_chamfer_mm: ch,
                ..
            } => {
                if !included_angle_deg.is_finite()
                    || !(0. ..180.).contains(&included_angle_deg)
                    || included_angle_deg == 0.
                    || !ch.is_finite()
                    || ch < 0.
                    || ch >= radius
                {
                    return Err(ThreadError::Invalid);
                }
                let theta = included_angle_deg.to_radians() / 2.;
                let rim = height - (radius - shank) / theta.tan();
                let neck_z = neck * (theta / 2.).tan();
                let plain = if union { plain.max(neck_z) } else { plain };
                result.end[1] = plain;
                let center = [shank + neck, neck_z];
                let tangent = [
                    center[0] - neck * theta.cos(),
                    center[1] - neck * theta.sin(),
                ];
                if rim < ch || plain < neck_z || top_z + rim > tangent[1] || radius <= tangent[0] {
                    return Err(ThreadError::Invalid);
                }
                if plain > neck_z && plain - neck_z <= reserve {
                    return Err(ThreadError::Invalid);
                }
                result.top[0] = radius - ch;
                push(result.top, [radius, top_z + ch], Curve::Line, Role::Chamfer);
                push(
                    [radius, top_z + ch],
                    [radius, top_z + rim],
                    Curve::Line,
                    Role::Rim,
                );
                push([radius, top_z + rim], tangent, Curve::Line, Role::Cone);
                if neck > 0. {
                    push(
                        tangent,
                        [shank, neck_z],
                        Curve::Arc {
                            center,
                            radius_mm: neck,
                            sweep: -theta,
                        },
                        Role::Neck,
                    );
                }
                push([shank, neck_z], [shank, plain], Curve::Line, Role::Shank);
            }
            Head::Crowned {
                crown_radius_mm: crown,
                top_land_radius_mm: land,
                ..
            } => {
                let plain = if union { plain.max(neck) } else { plain };
                result.end[1] = plain;
                if !crown.is_finite()
                    || !land.is_finite()
                    || crown < radius
                    || land < 0.
                    || land >= radius
                    || radius <= shank + neck
                    || plain < neck
                {
                    return Err(ThreadError::Invalid);
                }
                if crown > 1e7 {
                    return Err(ThreadError::Limit);
                }
                reserve += 64. * f64::EPSILON * crown;
                if plain > neck && plain - neck <= reserve {
                    return Err(ThreadError::Invalid);
                }
                let center_z = top_z + ((crown - land) * (crown + land)).sqrt();
                let skirt_z = center_z - ((crown - radius) * (crown + radius)).sqrt();
                if skirt_z >= 0. {
                    return Err(ThreadError::Invalid);
                }
                result.top[0] = land;
                let start_angle = (top_z - center_z).atan2(land);
                let end_angle = (skirt_z - center_z).atan2(radius);
                push(
                    result.top,
                    [radius, skirt_z],
                    Curve::Arc {
                        center: [0., center_z],
                        radius_mm: crown,
                        sweep: end_angle - start_angle,
                    },
                    Role::Crown,
                );
                push([radius, skirt_z], [radius, 0.], Curve::Line, Role::Skirt);
                push([radius, 0.], [shank + neck, 0.], Curve::Line, Role::Bearing);
                if neck > 0. {
                    push(
                        [shank + neck, 0.],
                        [shank, neck],
                        Curve::Arc {
                            center: [shank + neck, neck],
                            radius_mm: neck,
                            sweep: -FRAC_PI_2,
                        },
                        Role::Neck,
                    );
                }
                push([shank, neck], [shank, plain], Curve::Line, Role::Shank);
            }
        }
        if height + result.end[1] > 1e7 {
            return Err(ThreadError::Limit);
        }
        result.rounding_error_mm = reserve;
        let extra = match head {
            Head::Countersunk {
                edge_chamfer_mm, ..
            } => edge_chamfer_mm,
            Head::Crowned {
                top_land_radius_mm, ..
            } => top_land_radius_mm,
        };
        if [height, shank, radius - shank]
            .into_iter()
            .any(|v| v <= reserve)
            || [neck, plain, extra]
                .into_iter()
                .any(|v| v > 0. && v <= reserve)
            || result
                .segments()
                .any(|s| s.role != Role::Bearing && s.end[1] - s.start[1] <= reserve)
        {
            return Err(ThreadError::Invalid);
        }
        Ok(result)
    }

    pub fn segments(&self) -> impl Iterator<Item = Segment> + '_ {
        self.segments.iter().flatten().copied()
    }

    pub fn top(&self) -> [f64; 2] {
        self.top
    }
    pub fn end(&self) -> [f64; 2] {
        self.end
    }
    /// Reserve for source arithmetic; policy admission must keep this below
    /// its positional tolerance before publishing manufacturing geometry.
    pub fn rounding_error_mm(&self) -> f64 {
        self.rounding_error_mm
    }

    pub fn volume_mm3(&self) -> f64 {
        self.segments().map(Segment::volume_mm3).sum()
    }

    pub fn support(&self, direction: [f64; 3]) -> f64 {
        let radial = direction[0].hypot(direction[1]);
        self.segments()
            .map(|s| s.support(radial, direction[2]))
            .fold(f64::NEG_INFINITY, f64::max)
    }

    pub fn radius_at(&self, z: f64) -> Option<f64> {
        if !z.is_finite() || z < self.top[1] || z > self.end[1] {
            return None;
        }
        self.segments()
            .filter(|s| z >= s.start[1] && z <= s.end[1])
            .map(|s| s.radius_at(z))
            .reduce(f64::max)
    }

    /// Prove a blind circular envelope fits with the requested radial wall.
    /// A hex socket uses its corner radius, not its across-flats radius.
    pub fn admits_socket(&self, radius: f64, depth: f64, wall: f64) -> bool {
        if ![radius, depth, wall].into_iter().all(f64::is_finite)
            || radius <= 0.
            || depth <= 0.
            || wall <= 0.
            || self.top[1] + depth >= self.end[1]
        {
            return false;
        }
        let floor = self.top[1] + depth;
        let minimum = self
            .segments()
            .filter(|s| s.start[1] <= floor)
            .map(|s| {
                s.radius_at(s.start[1])
                    .min(s.radius_at(s.end[1].min(floor)))
            })
            .fold(f64::INFINITY, f64::min);
        minimum >= radius + wall + 2. * self.rounding_error_mm
    }
}

fn retained_segment_valid(segment: Segment, reserve: f64) -> bool {
    let dr = segment.end[0] - segment.start[0];
    let dz = segment.end[1] - segment.start[1];
    match segment.curve {
        Curve::Line => match segment.role {
            Role::Chamfer => dr > reserve && dz > reserve,
            Role::Cone => dr < -reserve && dz > reserve,
            Role::Rim | Role::Skirt | Role::Shank => dr == 0. && dz > reserve,
            Role::Bearing => dr < -reserve && dz == 0.,
            Role::Crown | Role::Neck => false,
        },
        Curve::Arc {
            center,
            radius_mm: radius,
            sweep,
        } => {
            if !sweep.is_finite() || center[0] < 0. || dz <= reserve {
                return false;
            }
            let start_radius = (segment.start[0] - center[0]).hypot(segment.start[1] - center[1]);
            if radius <= reserve || radius > 1e7 || (start_radius - radius).abs() > reserve {
                return false;
            }
            let a = segment.angle(center);
            let b = a + sweep;
            let angular = 32. * f64::EPSILON;
            let domain = match segment.role {
                Role::Crown => {
                    sweep > 0. && a >= -FRAC_PI_2 - angular && b <= angular && dr > reserve
                }
                Role::Neck => {
                    sweep < 0.
                        && a <= -FRAC_PI_2 + angular
                        && b >= -PI - angular
                        && dr < -reserve
                        && center[0] > radius
                }
                _ => false,
            };
            if !domain {
                return false;
            }
            let target = [center[0] + radius * b.cos(), center[1] + radius * b.sin()];
            (target[0] - segment.end[0]).hypot(target[1] - segment.end[1]) <= reserve
        }
    }
}

#[cfg(test)]
mod tests;
