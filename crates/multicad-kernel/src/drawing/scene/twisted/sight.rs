//! Material along a straight line, from the exact roots of every flank.
use super::super::paths::dot2;
use super::super::profile::merge;
use super::super::*;
use super::flank::{Flank, Shape, rotate};
use super::solve::{Search, isolate_touching};
use super::{Band, Gear, Region};
use std::f64::consts::{PI, TAU};

/// Relative distance from a coaxial circle within which a point is on it.
const ON_CIRCLE: f64 = 1e-12;

/// A line that stays within this fraction of the gear's size of a surface
/// along a stretch only grazes it there.
const GRAZING: f64 = 1e-9;

fn wrap(angle: f64) -> f64 {
    (angle + PI).rem_euclid(TAU) - PI
}

/// Polar angle of the involute frame's x axis and its handedness.
fn involute_frame(x: P2, y: P2) -> (f64, f64) {
    (x[1].atan2(x[0]), (x[0] * y[1] - x[1] * y[0]).signum())
}

/// Angles at which a flank meets the circle of a radius about the axis, and
/// whether counterclockwise travel along the circle enters the material.
fn circle_hits(flank: &Flank, radius: f64, output: &mut Vec<(f64, bool)>) {
    let slack = 1e-12 * radius;
    if radius < flank.radius[0] - slack || radius > flank.radius[1] + slack {
        return;
    }
    match flank.shape {
        Shape::Line { a, b } => {
            let d = [b[0] - a[0], b[1] - a[1]];
            let dd = dot2(d, d);
            let ad = dot2(a, d);
            let disc = ad * ad - dd * (dot2(a, a) - radius * radius);
            if disc < 0. {
                return;
            }
            for sign in [-1., 1.] {
                let s = (-ad + sign * disc.sqrt()) / dd;
                if (-1e-12..=1. + 1e-12).contains(&s) {
                    let p = [a[0] + s * d[0], a[1] + s * d[1]];
                    let outward = dot2(p, d);
                    if outward != 0. {
                        output.push((p[1].atan2(p[0]), outward > 0.));
                    }
                }
            }
        }
        Shape::Arc {
            center,
            radius: own,
            start,
            sweep,
        } => {
            let offset = center[0].hypot(center[1]);
            if offset == 0. {
                return;
            }
            let cosine = (radius * radius + offset * offset - own * own) / (2. * radius * offset);
            if cosine.abs() > 1. {
                return;
            }
            let toward = center[1].atan2(center[0]);
            for sign in [-1., 1.] {
                let angle = toward + sign * cosine.acos();
                let (sin, cos) = angle.sin_cos();
                let p = [radius * cos, radius * sin];
                let spoke = [p[0] - center[0], p[1] - center[1]];
                let travel = ((spoke[1].atan2(spoke[0]) - start) * sweep.signum()).rem_euclid(TAU);
                let slack = 1e-12;
                if travel <= sweep.abs() + slack || TAU - travel <= slack {
                    let outward = sweep.signum() * (p[1] * spoke[0] - p[0] * spoke[1]);
                    if outward != 0. {
                        output.push((angle, outward > 0.));
                    }
                }
            }
        }
        Shape::Involute {
            base, span, x, y, ..
        } => {
            let (zero, hand) = involute_frame(x, y);
            let t = ((radius / base).powi(2) - 1.).max(0.).sqrt();
            output.push((zero + hand * (t - t.atan()), span > 0.));
        }
    }
}

impl Gear {
    /// Whether a point of the profile plane is strictly inside the tooth
    /// outline: at its radius, the next flank counterclockwise is left.
    pub(super) fn contains_profile(&self, point: P2) -> bool {
        let radius = point[0].hypot(point[1]);
        let angle = point[1].atan2(point[0]);
        // On the root circle, within rounding, a point lies on a root arc
        // (the boundary) or under a tooth (inside).
        if (radius - self.radius[0]).abs() <= ON_CIRCLE * self.radius[0] {
            return !self.on_coaxial_arc(self.radius[0], angle);
        }
        if radius <= self.radius[0] {
            return true;
        }
        if radius >= self.radius[1] * (1. - ON_CIRCLE) {
            return false;
        }
        // The next flank is at most one tooth gap ahead; look two pitches
        // ahead first and around the whole outline only if none is there.
        let ahead = 2. * self.pitch;
        for all in [false, true] {
            let mut hits = Vec::new();
            for flank in &self.flanks {
                let [middle, half] = flank.azimuth;
                if all || (middle - angle + half).rem_euclid(TAU) <= ahead + 2. * half {
                    circle_hits(flank, radius, &mut hits);
                }
            }
            if let Some((_, enters)) = hits
                .into_iter()
                .map(|(at, enters)| ((at - angle).rem_euclid(TAU), enters))
                .filter(|(gap, _)| all || *gap <= ahead)
                .min_by(|a, b| a.0.total_cmp(&b.0))
            {
                return !enters;
            }
        }
        false
    }

    /// Whether a profile angle lies on a coaxial arc of the given radius.
    pub(super) fn on_coaxial_arc(&self, radius: f64, angle: f64) -> bool {
        self.flanks.iter().any(|flank| {
            flank
                .coaxial()
                .is_some_and(|own| (own - radius).abs() <= ON_CIRCLE * radius)
                && match flank.shape {
                    Shape::Arc { start, sweep, .. } => {
                        let low = start.min(start + sweep);
                        (angle - low).rem_euclid(TAU) <= sweep.abs()
                    }
                    _ => false,
                }
        })
    }

    /// Whether a point in gear coordinates is strictly inside the solid:
    /// material on both sides of its axial level.
    pub(super) fn contains(&self, [x, y, z]: P3) -> bool {
        let radius = x.hypot(y);
        // Points within rounding of a coaxial cylinder lie on its surface.
        radius > self.bore * (1. + ON_CIRCLE)
            && [-1., 1.]
                .into_iter()
                .all(|side| match self.material(z, side) {
                    Some(Region::Profile(phase)) => self.contains_profile(rotate([x, y], -phase)),
                    Some(Region::Disk(outer)) => radius < outer * (1. - ON_CIRCLE),
                    None => false,
                })
    }

    /// Parameters at which the line `base + t direction` can cross the
    /// boundary: every plane, coaxial cylinder and flank surface, the flanks
    /// unbounded along their profile curves.
    pub(super) fn stops(
        &self,
        base: P3,
        direction: P3,
    ) -> Result<(Vec<f64>, Vec<[f64; 2]>), GeometryError> {
        let a = [base[0], base[1]];
        let b = [direction[0], direction[1]];
        let rise = direction[2];
        let slope2 = dot2(b, b);
        let slope = slope2.sqrt();
        let transverse = slope > 1e-300;
        let nearest = if transverse { -dot2(a, b) / slope2 } else { 0. };
        let closest2 = if transverse {
            (dot2(a, a) - dot2(a, b) * dot2(a, b) / slope2).max(0.)
        } else {
            dot2(a, a)
        };
        let reach = |radius: f64| {
            (transverse && radius * radius > closest2)
                .then(|| (radius * radius - closest2).sqrt() / slope)
        };
        let mut stops = Vec::new();
        if rise != 0. {
            stops.extend(self.levels().into_iter().map(|z| (z - base[2]) / rise));
        }
        if transverse {
            // Where the line passes closest to the axis it can touch a
            // coaxial cylinder without crossing it.
            stops.push(nearest);
        }
        let mut radii = vec![self.bore];
        radii.extend(self.cylinders.iter().map(|cylinder| cylinder.radius));
        radii.extend(self.flanks.iter().filter_map(Flank::coaxial));
        radii.sort_by(f64::total_cmp);
        radii.dedup();
        for radius in radii {
            // A line tangent to a cylinder only touches it, at `nearest`.
            if radius * radius <= closest2 * (1. + 1e-12) {
                continue;
            }
            if let Some(half) = reach(radius) {
                stops.extend([nearest - half, nearest + half]);
            }
        }
        // A line along the axis on a coaxial arc's cylinder leaves the arc
        // where its profile angle passes one of the arc's ends.
        if !transverse && rise != 0. {
            let radius = closest2.sqrt();
            let angle = a[1].atan2(a[0]);
            for band in &self.bands {
                for flank in &self.flanks {
                    let (Some(own), Shape::Arc { start, sweep, .. }) =
                        (flank.coaxial(), flank.shape)
                    else {
                        continue;
                    };
                    if (own - radius).abs() > ON_CIRCLE * own {
                        continue;
                    }
                    for end in [start, start + sweep] {
                        // phase(z) = angle - end (mod a turn), z in the band.
                        let ends = band.z.map(|z| band.phase_at(z));
                        let [lo, hi] = [ends[0].min(ends[1]), ends[0].max(ends[1])];
                        let target = angle - end;
                        let mut k = ((lo - target) / TAU).ceil();
                        while target + k * TAU <= hi {
                            let phase = target + k * TAU;
                            let z = band.z[0] + (phase - band.phase) / band.twist;
                            stops.push((z - base[2]) / rise);
                            k += 1.;
                        }
                    }
                }
            }
        }
        let mut budget = 4_000_000_usize;
        let mut contacts = Vec::new();
        for band in &self.bands {
            let range = if rise != 0. {
                let ends = band.z.map(|z| (z - base[2]) / rise);
                [ends[0].min(ends[1]), ends[0].max(ends[1])]
            } else if base[2] >= band.z[0] && base[2] <= band.z[1] {
                [f64::NEG_INFINITY, f64::INFINITY]
            } else {
                continue;
            };
            for flank in &self.flanks {
                if flank.coaxial().is_some() {
                    continue;
                }
                if let Shape::Involute { .. } = flank.shape {
                    self.involute_roots(
                        flank,
                        band,
                        (base, direction),
                        range,
                        &mut budget,
                        (&mut stops, &mut contacts),
                    )?;
                    continue;
                }
                let pad = 1e-9 * flank.radius[1];
                let rings = [[flank.radius[0] - pad, flank.radius[1] + pad]];
                let mut windows: Vec<([f64; 2], [f64; 2])> = Vec::new();
                for ring in rings {
                    let [inner, outer] = ring;
                    if transverse {
                        let Some(far) = reach(outer) else { continue };
                        match reach(inner) {
                            Some(near) => windows.extend([
                                ([nearest - far, nearest - near], ring),
                                ([nearest + near, nearest + far], ring),
                            ]),
                            None => windows.push(([nearest - far, nearest + far], ring)),
                        }
                    } else if closest2.sqrt() >= inner && closest2.sqrt() <= outer {
                        windows.push((range, ring));
                    }
                }
                for (window, ring) in windows {
                    let window = [window[0].max(range[0]), window[1].min(range[1])];
                    if window[1] <= window[0] || !(window[1] - window[0]).is_finite() {
                        continue;
                    }
                    // The line's azimuth is monotone and its phase linear, so
                    // its image in the profile stays within their combined
                    // swing of its azimuth at the window's start.
                    let at = |t: f64| [a[0] + t * b[0], a[1] + t * b[1]];
                    let [p, q] = [at(window[0]), at(window[1])];
                    let swing = (p[0] * q[1] - p[1] * q[0]).atan2(dot2(p, q)).abs()
                        + (band.twist * rise * (window[1] - window[0])).abs();
                    let start = p[1].atan2(p[0]) - band.phase_at(base[2] + window[0] * rise);
                    if wrap(flank.azimuth[0] - start).abs() > swing + flank.azimuth[1] + 1e-9 {
                        continue;
                    }
                    self.flank_roots(
                        flank,
                        band,
                        (base, direction),
                        window,
                        ring,
                        &mut budget,
                        (&mut stops, &mut contacts),
                    )?;
                }
            }
        }
        // Only the part of the line within the solid's extent can matter.
        let mut extent = [f64::NEG_INFINITY, f64::INFINITY];
        if rise != 0. {
            let ends = self.axis.map(|z| (z - base[2]) / rise);
            extent = [ends[0].min(ends[1]), ends[0].max(ends[1])];
        }
        if transverse {
            let Some(far) = reach(self.radius[1] * (1. + 1e-9)) else {
                return Ok((Vec::new(), Vec::new()));
            };
            extent = [extent[0].max(nearest - far), extent[1].min(nearest + far)];
        }
        let pad = 1e-9 * (1. + extent[0].abs().max(extent[1].abs()));
        stops.retain(|t| *t >= extent[0] - pad && *t <= extent[1] + pad);
        Ok((stops, contacts))
    }

    /// Crossings of the line `base + t direction` with an involute flank
    /// swept through a band, through the flank's rulings. The swept involute
    /// is developable: with the foot angle `a = zero + hand u + phase(z)` of
    /// the involute parameter `u`, the points of constant `a` form the line
    /// `H(a) + v D(a)`, `H(a) = (b e(a), z(a))` on the base helix and
    /// `D(a) = (b k J e(a), 1)` its tangent, `k` the band's twist. The line
    /// meets the ruling of `a` exactly where `(H(a) - base) . (D(a) x d)`
    /// vanishes, a smooth function of `a` with bounded derivatives (unlike
    /// the involute's polar angle at its base circle). A line along a ruling
    /// runs in the surface: a contact.
    fn involute_roots(
        &self,
        flank: &Flank,
        band: &Band,
        (base, direction): (P3, P3),
        range: [f64; 2],
        budget: &mut usize,
        (stops, contacts): (&mut Vec<f64>, &mut Vec<[f64; 2]>),
    ) -> Result<(), GeometryError> {
        let Shape::Involute {
            base: circle,
            start,
            span,
            x,
            y,
        } = flank.shape
        else {
            return Ok(());
        };
        let (zero, hand) = involute_frame(x, y);
        let k = band.twist;
        let [t0, t1] = [start.min(start + span), start.max(start + span)];
        let [z0, z1] = band.z;
        // Foot angle of the ruling through (u, z) and its base-helix height.
        let foot = |u: f64, z: f64| zero + hand * u + band.phase_at(z);
        let height = |a: f64| z0 + (a - zero - band.phase) / k;
        let corners = [foot(t0, z0), foot(t0, z1), foot(t1, z0), foot(t1, z1)];
        let lo = corners.iter().copied().fold(f64::INFINITY, f64::min);
        let hi = corners.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let d = direction;
        let ruling = |a: f64| {
            let (sin, cos) = a.sin_cos();
            (
                [circle * cos, circle * sin, height(a)],
                [-circle * k * sin, circle * k * cos, 1.],
            )
        };
        let cross = |p: P3, q: P3| -> P3 {
            [
                p[1] * q[2] - p[2] * q[1],
                p[2] * q[0] - p[0] * q[2],
                p[0] * q[1] - p[1] * q[0],
            ]
        };
        let dot = |p: P3, q: P3| p[0] * q[0] + p[1] * q[1] + p[2] * q[2];
        let sub = |p: P3, q: P3| -> P3 { [p[0] - q[0], p[1] - q[1], p[2] - q[2]] };
        // The meeting of the line with the ruling of `a`, if on the flank.
        let meet = |a: f64| -> Option<f64> {
            let (h, axis) = ruling(a);
            let w = sub(h, base);
            let normal = cross(d, axis);
            let size = dot(normal, normal);
            if size == 0. {
                return None;
            }
            let t = dot(cross(w, axis), normal) / size;
            let v = dot(cross(w, d), normal) / size;
            let z = height(a) + v;
            let u = -hand * k * v;
            let slack = 1e-12 * (1. + t1.abs());
            (z >= z0 - 1e-12 * (1. + z1.abs())
                && z <= z1 + 1e-12 * (1. + z1.abs())
                && u >= t0 - slack
                && u <= t1 + slack
                && t >= range[0]
                && t <= range[1])
                .then_some(t)
        };
        // A line parallel to a ruling: along it when it passes through it.
        let planar = d[0].hypot(d[1]);
        let pitch = circle * k.abs() * d[2].abs();
        if (planar - pitch).abs() <= 1e-12 * (1. + planar) && d[2] != 0. {
            // J e(a) = d_xy / (b k d_z), so e(a) = -J (d_xy) / (b k d_z).
            let scale = circle * k * d[2];
            let a0 = (-d[0] / scale).atan2(d[1] / scale);
            let turns = ((lo - a0) / std::f64::consts::TAU).ceil();
            let mut a = a0 + turns * std::f64::consts::TAU;
            while a <= hi {
                let (h, axis) = ruling(a);
                let w = sub(h, base);
                let off = cross(w, axis);
                if dot(off, off).sqrt() <= 1e-9 * dot(axis, axis).sqrt() * (1. + dot(w, w).sqrt()) {
                    // The stretch of the ruling on the flank, in line terms.
                    let mut ends = Vec::new();
                    for (u, z) in [(t0, None), (t1, None), (0., Some(z0)), (0., Some(z1))] {
                        let v = match z {
                            Some(z) => z - height(a),
                            None => -u / (hand * k),
                        };
                        let z = height(a) + v;
                        let u = -hand * k * v;
                        if z >= z0 - 1e-12 && z <= z1 + 1e-12 && u >= t0 - 1e-12 && u <= t1 + 1e-12
                        {
                            let point = [h[0] + v * axis[0], h[1] + v * axis[1], z];
                            ends.push(dot(sub(point, base), d) / dot(d, d));
                        }
                    }
                    if let (Some(first), Some(last)) = (
                        ends.iter().copied().reduce(f64::min),
                        ends.iter().copied().reduce(f64::max),
                    ) && last > first
                    {
                        let stretch = [first.max(range[0]), last.min(range[1])];
                        if stretch[1] > stretch[0] {
                            stops.extend(stretch);
                            contacts.push(stretch);
                        }
                    }
                }
                a += std::f64::consts::TAU;
            }
        }
        // Bounds of the coplanarity function and its first two derivatives
        // over the foot angles of the flank.
        let q = planar;
        let w_bound = [lo, hi]
            .map(|a| (height(a) - base[2]).abs())
            .into_iter()
            .fold(0., f64::max);
        let [bx, by] = [base[0].abs(), base[1].abs()];
        let p1 = [circle + bx.max(by), circle, circle];
        let q1 = [
            circle * k.abs() * d[2].abs() + q,
            circle * k.abs() * d[2].abs(),
        ];
        let p3 = [w_bound, 1. / k.abs()];
        let q3 = circle * k.abs() * q;
        let first = 2. * (p1[1] * q1[0] + p1[0] * q1[1]) + p3[1] * q3 + p3[0] * q3;
        let second = 2. * (p1[2] * q1[0] + 2. * p1[1] * q1[1] + p1[0] * q1[1])
            + 2. * p3[1] * q3
            + p3[0] * q3;
        let coplanar = |a: f64| {
            let (h, axis) = ruling(a);
            dot(sub(h, base), cross(axis, d))
        };
        let mut roots = Vec::new();
        isolate_touching(
            &mut |a| coplanar(a),
            &Search {
                slope: first,
                leaf: 1e-13 * (1. + lo.abs().max(hi.abs())),
                // A line in the tangent plane along a ruling touches the
                // surface: the function has a double root there, which
                // rounding must not split into many.
                noise: 1e-13 * (2. * p1[0] * q1[0] + p3[0] * q3),
                second,
                lift: 0.,
                wrapped: false,
            },
            [lo, hi],
            budget,
            &mut roots,
            &mut Vec::new(),
        )?;
        stops.extend(roots.into_iter().filter_map(meet));
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn flank_roots(
        &self,
        flank: &Flank,
        band: &Band,
        (base, direction): (P3, P3),
        window: [f64; 2],
        [_, outer]: [f64; 2],
        budget: &mut usize,
        (stops, found): (&mut Vec<f64>, &mut Vec<[f64; 2]>),
    ) -> Result<(), GeometryError> {
        let mut contacts = Vec::new();
        let slope = direction[0].hypot(direction[1]);
        // Speed of the line's image in profile coordinates.
        let speed = slope + (band.twist * direction[2]).abs() * outer;
        let image = |t: f64| {
            let z = base[2] + t * direction[2];
            rotate(
                [base[0] + t * direction[0], base[1] + t * direction[1]],
                -band.phase_at(z),
            )
        };
        let leaf = 1e-9 * (1. + window[0].abs().max(window[1].abs()));
        // Bound of the image's acceleration: it turns at a uniform rate.
        let rate = (band.twist * direction[2]).abs();
        let bend = 2. * rate * slope + rate * rate * outer;
        let result = match flank.shape {
            Shape::Line { a, b } => {
                let d = [b[0] - a[0], b[1] - a[1]];
                let length = d[0].hypot(d[1]);
                let normal = [d[1] / length, -d[0] / length];
                isolate_touching(
                    &mut |t| {
                        let p = image(t);
                        dot2(normal, [p[0] - a[0], p[1] - a[1]])
                    },
                    &Search {
                        slope: speed,
                        leaf,
                        noise: GRAZING * outer,
                        second: bend,
                        lift: 0.,
                        wrapped: false,
                    },
                    window,
                    budget,
                    stops,
                    &mut contacts,
                )
            }
            Shape::Arc { center, radius, .. } => isolate_touching(
                &mut |t| {
                    let p = image(t);
                    (p[0] - center[0]).hypot(p[1] - center[1]) - radius
                },
                &Search {
                    slope: speed,
                    leaf,
                    noise: GRAZING * outer,
                    // Of the squared distance from the center.
                    second: 2. * speed * speed + 2. * (outer + center[0].hypot(center[1])) * bend,
                    lift: radius,
                    wrapped: false,
                },
                window,
                budget,
                stops,
                &mut contacts,
            ),
            // Involutes are intersected through their rulings.
            Shape::Involute { .. } => Ok(()),
        };
        result?;
        // A stretch along the surface of the unbounded curve counts only
        // where it lies on the flank itself.
        // A transversal root is near zero over a short stretch too; only a
        // stretch long against the gear is a line running in the surface.
        let shortest = 1e-5 * (1. + window[0].abs().max(window[1].abs()));
        found.extend(
            contacts
                .into_iter()
                .filter(|[a, b]| b - a >= shortest && flank.covers(image(0.5 * (a + b)))),
        );
        Ok(())
    }

    pub(super) fn spans(
        &self,
        base: P3,
        direction: P3,
        inside: &dyn Fn(P3) -> bool,
    ) -> Result<Vec<[f64; 2]>, GeometryError> {
        let (mut stops, contacts) = self.stops(base, direction)?;
        stops.extend(contacts.iter().flatten());
        stops.retain(|t| t.is_finite());
        stops.sort_by(f64::total_cmp);
        let scale = stops
            .iter()
            .fold(1_f64, |scale, t| scale.max(t.abs()))
            .max(self.radius[1]);
        // Shorter contacts are grazing.
        let grazing = GRAZING * scale;
        let mut output = Vec::new();
        for pair in stops.windows(2) {
            let middle = 0.5 * (pair[0] + pair[1]);
            // Along a stretch where the line runs in a flank it only touches
            // the boundary.
            if pair[1] - pair[0] <= grazing
                || contacts
                    .iter()
                    .any(|[a, b]| *a - grazing <= middle && middle <= *b + grazing)
            {
                continue;
            }
            // A surface can touch the line at one point between two stops;
            // two of three interior points decide.
            let votes = [0.3, 0.5, 0.7]
                .into_iter()
                .filter(|f| {
                    let t = pair[0] + (pair[1] - pair[0]) * f;
                    inside(std::array::from_fn(|i| base[i] + t * direction[i]))
                })
                .count();
            if votes >= 2 {
                output.push([pair[0], pair[1]]);
            }
        }
        Ok(merge(output, 2. * grazing))
    }
}
