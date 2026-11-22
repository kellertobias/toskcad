//! Touching conjugate flanks of a worm and the wheel it generates.
//!
//! In each wheel section `v` the generating cutter acts as a rack. A
//! cutter flank point at cutter radius `r` sits, at rack pose `τ`, at
//! `Φ(r, τ) = R(-τ/ρₚ) (C - q(r), τ - z(r))` in sector-0 coordinates, with
//! `q = sqrt(r² - v²)`, the axial flank position `z(r)` of
//! `Wheel::contact`, center distance `C` and pitch radius `ρₚ`. A worm
//! without thinning, on a wheel without thinning, has exactly these flanks:
//! its lower flank is the minus cutter flank and its upper flank the plus
//! flank, at one pose `τ₀` for sector 0 and `τ₀ + kπ` for sector `k`. The
//! pose is measured from the placed worm and checked on both flanks.
//!
//! At fixed wheel radius `ρ`, let `Θ(ρ, τ)` be the polar angle where the
//! cutter flank at pose `τ` crosses `ρ`. The generated flank is `Θ` at the
//! pose `τ*(ρ)` where `∂Θ/∂τ = 0`, the contact of `Wheel::contact`. When
//! `σ ∂²Θ/∂τ² > 0` on every pose between `τ*` and the worm pose, the worm
//! flank angle differs from the wheel flank angle by `σ ∫∫ ∂²Θ/∂τ² ≥ 0`
//! toward the tooth space (`σ = 1` for the minus flank). Both flanks are
//! radial graphs over the flank range, proved by a negative radius
//! derivative, so a point between them at that radius lies in the closed
//! tooth space. The curvature is enclosed by interval Taylor series in `τ`,
//! with the cutter radius `r(ρ, τ)` solved implicitly from
//! `(C - q)² + (τ - z)² = ρ²`.
use super::clearance::Placement;
use super::interval::{Interval, Series};
use super::jet::Real;
use super::pair_distance::{V3, WormMap, k, mid, square};
use super::wheel::{CUTTER_TIP, WORM_ROOT, Wheel, turn};
use super::worm::Worm;
use super::*;
use std::f64::consts::PI;

type S3 = Series<Interval, 3>;
type S2 = Series<Interval, 2>;

/// Cells of the cutter range searched for an implicit radius.
const RADIUS_CELLS: usize = 32;
/// Worm crest radius in module units.
const WORM_CREST: f64 = 6.;

pub(super) struct Contact<'a> {
    wheel: &'a Wheel,
    pub(super) map: WormMap<'a>,
    /// Rack pose of sector 0.
    pose: f64,
    /// Cutter side of the worm's lower flank.
    lower: f64,
    center: f64,
    pitch: f64,
}

/// Hull of the cutter-range cells that may hold a root, refined by
/// bisection: a discarded cell is proved free of roots.
fn solve(mut keep: impl FnMut(Interval) -> bool) -> Option<Interval> {
    let mut live: Vec<Interval> = super::wheel_bounds::cells(WORM_ROOT, CUTTER_TIP, RADIUS_CELLS)
        .filter(|c| keep(*c))
        .collect();
    for _ in 0..40 {
        if live.is_empty() || live.len() > 64 {
            break;
        }
        let width: f64 = live.iter().map(|c| c.hi - c.lo).sum();
        if width < 1e-12 {
            break;
        }
        live = live
            .into_iter()
            .flat_map(|c| [Interval::new(c.lo, mid(c)), Interval::new(mid(c), c.hi)])
            .filter(|c| keep(*c))
            .collect();
    }
    live.into_iter().reduce(hull)
}

fn hull(a: Interval, b: Interval) -> Interval {
    Interval::new(a.lo.min(b.lo), a.hi.max(b.hi))
}

impl<'a> Contact<'a> {
    /// The rack pose of a worm without thinning on its generated wheel, or
    /// `None` when the placed flanks are not those cutter flanks.
    pub(super) fn new(
        worm: &'a Worm,
        worm_at: Placement,
        wheel: &'a Wheel,
        wheel_at: Placement,
    ) -> Option<Self> {
        if worm.spec().axial_thinning_mm != 0. || wheel.spec().pitch_thinning_mm != 0. {
            return None;
        }
        let m = wheel.module_mm();
        if (worm.spec().axial_module_mm - m).abs() > 1e-12 * m {
            return None;
        }
        let [center, pitch] = wheel.generating_radii();
        let map = WormMap::new(worm, worm_at, wheel, wheel_at);
        let mut contact = Self {
            wheel,
            map,
            pose: 0.,
            lower: 1.,
            center,
            pitch,
        };
        // A worm turned end for end meets the wheel with its lower flank on
        // the plus cutter side.
        for lower in [1., -1.] {
            contact.lower = lower;
            if let Some(pose) = contact.common_pose() {
                contact.pose = pose;
                return Some(contact);
            }
        }
        None
    }

    /// Cutter side of a worm flank band face, or `None` for other faces.
    pub(super) fn side(&self, face: u32) -> Option<f64> {
        match face {
            2 => Some(self.lower),
            4 => Some(-self.lower),
            _ => None,
        }
    }

    /// The common sector-0 rack pose of sampled points on both flanks.
    fn common_pose(&self) -> Option<f64> {
        let wheel = self.wheel;
        let m = wheel.module_mm();
        let bands = super::pair_distance::worm_patches(self.map.worm);
        let mut poses = vec![];
        for (band, side) in [(0, self.lower), (2, -self.lower)] {
            let (patch, [s, z]) = bands[band];
            for i in 1..4 {
                for j in 0..64 {
                    let x = [
                        s.lo + (s.hi - s.lo) * f64::from(i) / 4.,
                        z.lo + (z.hi - z.lo) * f64::from(j) / 64.,
                    ];
                    let (local, _) = self.map.point(patch, x);
                    // Only the near side of the worm faces the wheel.
                    if local[2].abs() >= wheel.half_width()
                        || local[0].portable_hypot(local[1]) >= wheel.blank()
                    {
                        continue;
                    }
                    let thread = self.map.thread_point(patch, x);
                    let r = thread[0].portable_hypot(thread[1]) / m;
                    poses.push(self.measure(local, r, side)?);
                }
            }
        }
        let first = *poses.first()?;
        poses
            .iter()
            .all(|p| (p - first).abs() <= 1e-9)
            .then_some(first)
    }

    fn sector(&self, angle: f64) -> f64 {
        (angle / (2. * self.wheel.span())).round()
    }

    /// Sector-0 rack pose of one worm flank point, checked against `Φ`.
    fn measure(&self, local: Point, r: f64, side: f64) -> Option<f64> {
        let k = self.sector(local[1].portable_atan2(local[0]));
        let q = rotate_z(local, -k * 2. * self.wheel.span());
        let v = q[2];
        let rho = q[0].portable_hypot(q[1]);
        let [qq, z, _] = self.wheel.contact(r, v, side);
        let x = self.center - qq;
        let reach = (rho * rho - x * x).max(0.).sqrt();
        [1., -1.].into_iter().find_map(|sign| {
            let tau = z + sign * reach;
            let [px, py] = turn([x, tau - z], -tau / self.pitch);
            (norm(sub([px, py, v], q)) < 1e-9).then_some(tau - k * PI)
        })
    }

    /// Cutter radii whose generated flank of `side` reaches radius `rho`
    /// at sections `v`, after proving the radius decreases along the flank.
    fn flank_radius(&self, rho: Interval, v: Interval, side: f64) -> Option<Interval> {
        let target = square(rho);
        let mut monotone = true;
        let found = solve(|cell| {
            let [x, y] = self.wheel.flank(S2::variable(cell), S2::lift(v), side);
            let r2 = x * x + y * y;
            monotone &= r2.0[1].hi < 0.;
            r2.0[0].hi >= target.lo && r2.0[0].lo <= target.hi
        });
        found.filter(|_| monotone)
    }

    /// Cutter radii on the rack flank at poses `tau` crossing radius `rho`.
    fn rack_radius(
        &self,
        rho: Interval,
        tau: Interval,
        v: Interval,
        side: f64,
    ) -> Option<Interval> {
        solve(|cell| {
            let g = self.rack_gap(cell, tau, v, rho, side);
            g.lo <= 0. && g.hi >= 0.
        })
    }

    /// `(C - q)² + (τ - z)² - ρ²` over any `Real`.
    fn rack_gap<T: Real>(&self, r: T, tau: T, v: T, rho: T, side: f64) -> T {
        let [q, z, _] = self.wheel.contact(r, v, side);
        let x = T::constant(self.center) - q;
        let y = tau - z;
        x * x + y * y - rho * rho
    }

    /// Whether every point of a worm flank box image is in the closed tooth
    /// space of its sector. `side` is that flank's cutter side.
    pub(super) fn covers(&self, value: V3, side: f64) -> bool {
        self.space(value, side, false)
    }

    /// Whether every point of a worm crest box image is in the closed tooth
    /// space of its sector.
    pub(super) fn covers_crest(&self, value: V3) -> bool {
        self.space(value, 0., true)
    }

    /// The worm box image lies in the closed tooth space: flank points
    /// (`side`) or crest points (`on_crest`).
    fn space(&self, value: V3, side: f64, on_crest: bool) -> bool {
        let span = self.wheel.span();
        let sector = self.sector(mid(value[1]).portable_atan2(mid(value[0])));
        let (s, c) = k(-sector * 2. * span).sin_cos();
        let x = c * value[0] - s * value[1];
        let y = s * value[0] + c * value[1];
        let v = value[2];
        let h = self.wheel.half_width();
        if v.lo > h || v.hi < -h {
            return true;
        }
        if x.lo <= 0. {
            return false;
        }
        // Points beyond the end planes are outside the wheel; the rest have
        // sections in the clamped range.
        let v = Interval::new(v.lo.max(-h), v.hi.min(h));
        let rho = (square(x) + square(y)).sqrt();
        let blank = self.wheel.blank();
        if rho.lo > blank {
            return true;
        }
        let rho = Interval::new(rho.lo, rho.hi.min(blank));
        let theta = (y / x).atan();
        if theta.lo <= -span || theta.hi >= span {
            return false;
        }
        let tau = k(self.pose + sector * PI);
        // Rack coordinates of crest points: the axial position a = τ - Y.
        let crest = on_crest.then(|| {
            let (s, c) = (tau / k(self.pitch)).sin_cos();
            tau - (s * x + c * y)
        });
        [1., -1.].into_iter().all(|flank| {
            let Some(radius) = self.flank_radius(rho, v, flank) else {
                return false;
            };
            // Below the cutter tip the boundary leaves the flank.
            if radius.hi >= CUTTER_TIP {
                return false;
            }
            let [q, z, w] = self.wheel.contact(radius, v, flank);
            let x_f = k(self.center) - q;
            let wheel_angle = (w / x_f).atan() - (z + w) / k(self.pitch);
            // The minus flank needs θ ≥ θ_F, the plus flank θ ≤ θ_F.
            let gap = theta - wheel_angle;
            if flank > 0. && gap.lo > 0. || flank < 0. && gap.hi < 0. {
                return true;
            }
            if crest.is_none() && flank != side {
                return false;
            }
            match crest {
                None => self.reaches(rho, v, tau, z + w, flank),
                // A crest point at radius ρ is the rail point of this flank
                // at pose τ + (z_rail - a), turned toward the space by
                // `(z_rail - a)/ρₚ` (rack points move by -τ/ρₚ).
                Some(a) => {
                    let [_, rail, _] = self.wheel.contact(k(WORM_CREST), v, flank);
                    // Crest points lie between the two rails, so the shift
                    // has the sign of `flank`; meet the enclosure with it.
                    let shift = rail - a;
                    let shift = if flank > 0. {
                        Interval::new(shift.lo.max(0.), shift.hi)
                    } else {
                        Interval::new(shift.lo, shift.hi.min(0.))
                    };
                    shift.lo <= shift.hi && self.reaches(rho, v, tau + shift, z + w, flank)
                }
            }
        })
    }

    /// Enclosures of `∂Θ/∂τ` and of a positive multiple of `∂²Θ/∂τ²` over
    /// radii `rho`, sections `v` and poses `tau`.
    fn slopes(
        &self,
        rho: Interval,
        v: Interval,
        tau: Interval,
        side: f64,
    ) -> Option<[Interval; 2]> {
        let r0 = self.rack_radius(rho, tau, v, side)?;
        let slope = self
            .rack_gap(
                S2::variable(r0),
                S2::lift(tau),
                S2::lift(v),
                S2::lift(rho),
                side,
            )
            .0[1];
        if slope.lo <= 0. && slope.hi >= 0. {
            return None;
        }
        let t = S3::variable(tau);
        let (vs, rs) = (S3::lift(v), S3::lift(rho));
        let mut r = S3::lift(r0);
        for order in 1..3 {
            let g = self.rack_gap(r, t, vs, rs, side);
            r.0[order] = -(g.0[order] / slope);
        }
        let [q, z, _] = self.wheel.contact(r, vs, side);
        let a = (t - z) / (S3::constant(self.center) - q);
        let [a0, a1, a2] = a.0;
        let base = k(1.) + a0 * a0;
        // Θ = atan(a) - τ/ρₚ: its τ coefficient, and the numerator of its
        // τ² coefficient, which has the sign of the curvature.
        Some([a1 / base - k(1. / self.pitch), a2 * base - a0 * a1 * a1])
    }

    /// The worm pose `tau` reaches the wheel angle from the space side:
    /// `σ Θ` falls monotonically from `tau` toward the contact poses `star`,
    /// proved by the sign of `∂Θ/∂τ` away from `star` and by `σ ∂²Θ/∂τ² > 0`
    /// on pose cells meeting `star`, where `∂Θ/∂τ` vanishes.
    fn reaches(
        &self,
        rho: Interval,
        v: Interval,
        tau: Interval,
        star: Interval,
        side: f64,
    ) -> bool {
        let span = hull(tau, star);
        let mut pending: Vec<(Interval, u32)> = super::wheel_bounds::cells(span.lo, span.hi, 4)
            .map(|c| (c, 0))
            .collect();
        while let Some((cell, depth)) = pending.pop() {
            let meets = cell.hi >= star.lo && cell.lo <= star.hi;
            let proved = self.slopes(rho, v, cell, side).is_some_and(|[d1, d2]| {
                if meets {
                    if side > 0. { d2.lo > 0. } else { d2.hi < 0. }
                } else if cell.hi < star.lo {
                    // Before the contact pose σΘ decreases.
                    if side > 0. { d1.hi < 0. } else { d1.lo > 0. }
                } else if side > 0. {
                    d1.lo > 0.
                } else {
                    d1.hi < 0.
                }
            });
            if !proved {
                if depth >= 8 {
                    return false;
                }
                let [a, b] = [
                    Interval::new(cell.lo, mid(cell)),
                    Interval::new(mid(cell), cell.hi),
                ];
                pending.extend([(a, depth + 1), (b, depth + 1)]);
            }
        }
        true
    }

    /// A contact point on a worm flank and the identical wheel flank point,
    /// in world mm: the generated point whose contact pose is the worm pose.
    pub(super) fn witness(&self) -> Option<[Point; 2]> {
        let h = self.wheel.half_width();
        let m = self.wheel.module_mm();
        let teeth = i64::from(self.wheel.teeth());
        let mut best: Option<(f64, Point, f64)> = None;
        for side in [1., -1.] {
            for step in 0..=40 {
                let r = WORM_ROOT + 0.05 + (5.95 - WORM_ROOT) * f64::from(step) / 40.;
                let pose_at = |v: f64| {
                    let [_, z, w] = self.wheel.contact(r, v, side);
                    z + w
                };
                for k in -teeth..=teeth {
                    let tau = self.pose + k as f64 * PI;
                    let (mut lo, mut hi) = (-h, h);
                    let flo = pose_at(lo) - tau;
                    if flo * (pose_at(hi) - tau) > 0. {
                        continue;
                    }
                    for _ in 0..80 {
                        let v = (lo + hi) / 2.;
                        if (pose_at(v) - tau) * flo > 0. {
                            lo = v;
                        } else {
                            hi = v;
                        }
                    }
                    let v = (lo + hi) / 2.;
                    let [x, y] = self.wheel.flank(r, v, side);
                    let point = rotate_z([x, y, v], k as f64 * 2. * self.wheel.span());
                    if best.is_none_or(|(b, ..)| v.abs() < b) {
                        best = Some((v.abs(), point, side));
                    }
                }
            }
        }
        let (_, point, side) = best?;
        let (worm_at, wheel_at) = self.map.to_world;
        let on_wheel = wheel_at.world(point.map(|c| c * m));
        // The worm flank point at this pose is `Φ(r, τ)`, the same point.
        // Take the worm's own flank point at its thread height and phase.
        let worm = self.map.worm;
        let thread = worm.to_thread(worm_at.local(on_wheel));
        let rate = worm.thread().helical_rate();
        let pitch = worm.thread().spec().pitch_mm;
        let phases = worm.rail_phases();
        let band = if side == self.lower { 0 } else { 2 };
        let (lo, hi) = (phases[band + 1], phases[band]);
        let raw = thread[2] - rate * thread[1].portable_atan2(thread[0]);
        let phase = lo + (raw - lo).rem_euclid(pitch);
        if phase > hi {
            return None;
        }
        let on_worm = worm_at.world(worm.to_local(worm.wall_point(phase, thread[2])));
        Some([on_worm, on_wheel])
    }
}

#[cfg(test)]
mod tests {
    use super::super::clearance::{ResultKind, certify};
    use super::super::{wheel::tests::spec as wheel_spec, worm::tests::spec as worm_spec};
    use super::*;

    fn pose(teeth: u32, drive: f64, extra: f64) -> [Placement; 2] {
        let center = f64::from(teeth) / 2. + 5.;
        let (sw, cw) = drive.portable_sin_cos();
        let turn = -drive / f64::from(teeth) + extra;
        let (s, c) = turn.portable_sin_cos();
        [
            Placement {
                origin_mm: [center, 0., 0.],
                axes: [[cw, 0., -sw], [-sw, 0., -cw], [0., 1., 0.]],
            },
            Placement {
                origin_mm: [0.; 3],
                axes: [[c, s, 0.], [-s, c, 0.], [0., 0., 1.]],
            },
        ]
    }

    #[test]
    fn unthinned_conjugate_flanks_touch_and_any_turn_overlaps() {
        let worm = Worm::new(worm_spec(1., 16., 1., 0.), 1e-3).unwrap();
        let wheel = Wheel::new(wheel_spec(1., 32, 2., 2., 0.), 1e-3, 1e7).unwrap();
        for drive in [2.5, 0., 0.7] {
            let [a, b] = pose(32, drive, 0.);
            let found = certify(&worm, a, &wheel, b, 1e-7, 0., None).unwrap();
            assert_eq!(found.result, ResultKind::Contact, "drive {drive}");
            let witness = found.nearest.unwrap();
            assert!(witness.distance_mm < 1e-9, "{witness:?}");
            assert!(worm.margin(a.local(witness.points_mm[0])).abs() < 1e-9);
            assert!(wheel.margin(b.local(witness.points_mm[1])).abs() < 1e-9);
            // Without backlash either turn of the wheel crosses a flank.
            for extra in [2e-3, -2e-3] {
                let [a, b] = pose(32, drive, extra);
                let turned = certify(&worm, a, &wheel, b, 1e-7, 0., None).unwrap();
                assert_eq!(turned.result, ResultKind::CommonInterior, "{drive} {extra}");
            }
        }
        // Thinned members keep positive clearance instead of contact.
        let thin = Worm::new(worm_spec(1., 16., 1., 0.05), 1e-3).unwrap();
        let [a, b] = pose(32, 0., 0.);
        assert!(Contact::new(&thin, a, &wheel, b).is_none());
    }
}
