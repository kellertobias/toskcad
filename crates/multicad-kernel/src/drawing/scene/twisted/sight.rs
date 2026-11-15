//! Material along a straight line, from the exact roots of every flank.
use super::super::paths::dot2;
use super::super::profile::merge;
use super::super::*;
use super::flank::{Flank, Shape, rotate};
use super::solve::{Search, isolate};
use super::{Band, Gear};
use std::f64::consts::{PI, TAU};

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
        if radius <= self.radius[0] {
            return true;
        }
        if radius >= self.radius[1] {
            return false;
        }
        let angle = point[1].atan2(point[0]);
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

    /// Whether a point in gear coordinates is strictly inside the solid.
    pub(super) fn contains(&self, [x, y, z]: P3) -> bool {
        let radius = x.hypot(y);
        if radius <= self.bore {
            return false;
        }
        if z > 0. && z < self.height {
            return self.contains_profile(rotate([x, y], -self.phase_at(z)));
        }
        self.extensions.iter().any(|extension| {
            radius < extension.radius
                && ((z > extension.z[0] && z < extension.z[1])
                    || (z == extension.z[1 - extension.outer]))
        })
    }

    /// Parameters at which the line `base + t direction` can cross the
    /// boundary: every plane, coaxial cylinder and flank surface, the flanks
    /// unbounded along their profile curves.
    fn stops(&self, base: P3, direction: P3) -> Result<Vec<f64>, GeometryError> {
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
            let mut levels = vec![self.axis[0], self.axis[1]];
            levels.extend(self.bands.iter().flat_map(|band| band.z));
            stops.extend(levels.into_iter().map(|z| (z - base[2]) / rise));
        }
        if transverse {
            // Where the line passes closest to the axis it can touch a
            // coaxial cylinder without crossing it.
            stops.push(nearest);
        }
        let mut radii = vec![self.bore];
        radii.extend(self.extensions.iter().map(|extension| extension.radius));
        radii.extend(self.flanks.iter().filter_map(Flank::coaxial));
        radii.sort_by(f64::total_cmp);
        radii.dedup();
        for radius in radii {
            if let Some(half) = reach(radius) {
                stops.extend([nearest - half, nearest + half]);
            }
        }
        let mut budget = 4_000_000_usize;
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
                let pad = 1e-9 * flank.radius[1];
                // An involute's angle has unbounded curvature at its base
                // circle only: search a thin ring there apart from the rest.
                let mut rings = vec![[flank.radius[0] - pad, flank.radius[1] + pad]];
                if let Shape::Involute { base: circle, .. } = flank.shape {
                    let split = circle * 0.06_f64.hypot(1.);
                    if split > rings[0][0] && split < rings[0][1] {
                        rings = vec![[rings[0][0], split], [split, rings[0][1]]];
                    }
                }
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
                        &mut stops,
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
                return Ok(Vec::new());
            };
            extent = [extent[0].max(nearest - far), extent[1].min(nearest + far)];
        }
        let pad = 1e-9 * (1. + extent[0].abs().max(extent[1].abs()));
        stops.retain(|t| *t >= extent[0] - pad && *t <= extent[1] + pad);
        Ok(stops)
    }

    #[allow(clippy::too_many_arguments)]
    fn flank_roots(
        &self,
        flank: &Flank,
        band: &Band,
        (base, direction): (P3, P3),
        window: [f64; 2],
        [inner, outer]: [f64; 2],
        budget: &mut usize,
        stops: &mut Vec<f64>,
    ) -> Result<(), GeometryError> {
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
        match flank.shape {
            Shape::Line { a, b } => {
                let d = [b[0] - a[0], b[1] - a[1]];
                let length = d[0].hypot(d[1]);
                let normal = [d[1] / length, -d[0] / length];
                isolate(
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
                )
            }
            Shape::Arc { center, radius, .. } => isolate(
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
            ),
            Shape::Involute {
                base: circle, x, y, ..
            } => {
                let (zero, hand) = involute_frame(x, y);
                let inner = inner.max(circle);
                let roll = ((outer / circle).powi(2) - 1.).max(0.).sqrt();
                let first = ((inner / circle).powi(2) - 1.).max(0.).sqrt();
                isolate(
                    &mut |t| {
                        let p = image(t);
                        let u = ((dot2(p, p) / (circle * circle)) - 1.).max(0.).sqrt();
                        wrap(p[1].atan2(p[0]) - zero - hand * (u - u.atan()))
                    },
                    &Search {
                        slope: (speed + roll * slope) / inner,
                        leaf,
                        noise: GRAZING,
                        // The involute angle's curvature is unbounded at
                        // the base circle; away from it both terms are.
                        second: if first > 0.05 {
                            bend / inner
                                + 2. * (speed / inner).powi(2)
                                + (1. / (circle * circle * first) + roll / (inner * inner))
                                    * slope
                                    * slope
                                + roll * slope * slope / (inner * inner)
                        } else {
                            f64::INFINITY
                        },
                        lift: 0.,
                        wrapped: true,
                    },
                    window,
                    budget,
                    stops,
                )
            }
        }
    }

    /// Ascending disjoint parameter intervals of the line on which the
    /// membership predicate holds. Candidate ends are exact boundary roots;
    /// membership is decided inside every candidate interval.
    pub(super) fn spans(
        &self,
        base: P3,
        direction: P3,
        inside: &dyn Fn(P3) -> bool,
    ) -> Result<Vec<[f64; 2]>, GeometryError> {
        let mut stops = self.stops(base, direction)?;
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
            if pair[1] - pair[0] <= grazing {
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
