//! Phase-resolved whole-box classification against the authored pipe law.
//!
//! Material is `rho <= R(q, z)` for `0 <= z <= length`, with the helical
//! phase `q = z - pitch * theta / (2 pi)` taken modulo the pitch. A box's
//! angular range (when it avoids the axis) and axial range give a phase
//! enclosure; `radial_bounds` encloses the law there. Only a whole-box
//! relation is published: boxes reaching the law surface or a cap are
//! `None`. Callers dilate the box by any representation error first.
use super::*;

/// Whole-box relation to the authored material.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AuthoredBoxRelation {
    Inside,
    Outside,
}

/// Conservative angular enclosure of a box in the xy plane avoiding the axis.
fn angle_range(x: [f64; 2], y: [f64; 2]) -> Option<[f64; 2]> {
    if x[0] <= 0. && x[1] >= 0. && y[0] <= 0. && y[1] >= 0. {
        return None;
    }
    let center = (0.5 * (y[0] + y[1])).atan2(0.5 * (x[0] + x[1]));
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    for px in x {
        for py in y {
            // Wrap each corner around the centre direction: a convex box
            // that avoids the origin subtends less than a half turn.
            let mut d = py.atan2(px) - center;
            if d > PI {
                d -= 2. * PI;
            } else if d < -PI {
                d += 2. * PI;
            }
            lo = lo.min(d);
            hi = hi.max(d);
        }
    }
    // atan2 and the subtractions are within a few ulps of a half turn.
    let pad = 64. * f64::EPSILON * PI;
    (hi - lo < PI - pad).then_some([center + lo - pad, center + hi + pad])
}

impl AuthoredThread {
    /// Classify every point of a closed box in the authored construction chart.
    pub fn classify_box(
        &self,
        bounds: [[f64; 2]; 3],
    ) -> Result<Option<AuthoredBoxRelation>, Error> {
        let [x, y, z] = bounds.map(interval).map(|r| r.ok());
        let (Some(xi), Some(yi), Some(zi)) = (x, y, z) else {
            return Err(Error::Invalid);
        };
        let length = self.spec.length_mm;
        let pitch = self.spec.pitch_mm;
        if zi.hi < 0. || zi.lo > length {
            return Ok(Some(AuthoredBoxRelation::Outside));
        }
        let axial = [zi.lo.max(0.), zi.hi.min(length)];
        let full = [[0., pitch]];
        let mut pieces = Vec::with_capacity(2);
        match angle_range(bounds[0], bounds[1]) {
            None => pieces.extend(full),
            Some([t0, t1]) => {
                let rate = I::point(pitch) / (I::point(2.) * pi());
                let q = I::new(zi.lo, zi.hi) - rate * I::new(t0, t1);
                if !q.valid() || q.width() >= pitch {
                    pieces.extend(full);
                } else {
                    let turns = (q.lo / pitch).floor();
                    let shifted = q - I::point(turns) * I::point(pitch);
                    let [lo, hi] = [shifted.lo.max(0.), shifted.hi];
                    if hi <= pitch {
                        pieces.push([lo, hi]);
                    } else {
                        pieces.push([lo.min(pitch), pitch]);
                        pieces.push([0., (hi - pitch).min(pitch)]);
                    }
                }
            }
        }
        let mut radius: Option<I> = None;
        for phase in pieces {
            let r = interval(self.radial_bounds(phase, axial)?)?;
            radius = Some(radius.map_or(r, |old| hull(old, r)));
        }
        let radius = radius.ok_or(Error::Unresolved)?;
        let rho = xi.square() + yi.square();
        if !rho.valid() || !radius.valid() {
            return Err(Error::Unresolved);
        }
        if rho.lo > I::point(radius.hi.max(0.)).square().hi {
            return Ok(Some(AuthoredBoxRelation::Outside));
        }
        if zi.lo > 0.
            && zi.hi < length
            && radius.lo > 0.
            && rho.hi < I::point(radius.lo).square().lo
        {
            return Ok(Some(AuthoredBoxRelation::Inside));
        }
        Ok(None)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn npt() -> AuthoredThread {
        // NPT 1/8-27, sourced truncation 0.0645 pitch, 1/32 radial taper.
        let pitch = 25.4 / 27.;
        let (m, k, truncation) = (3_f64.sqrt(), 1. / 32., 0.0645);
        let fundamental = pitch * (m * m - k * k) / (2. * m);
        let depth = fundamental - 2. * truncation * pitch;
        let width = truncation * pitch * (1. / (m - k) + 1. / (m + k));
        AuthoredThread::new(
            AuthoredSpec {
                major_radius_mm: 5.,
                pitch_mm: pitch,
                length_mm: 6.,
                radial_taper: k,
                section: AuthoredSection::NptFlat {
                    depth_mm: depth,
                    flat_width_mm: width,
                },
                crest_radius_limit_mm: None,
            },
            1e-6,
            1e-4,
        )
        .unwrap()
    }

    fn cube(center: [f64; 3], half: f64) -> [[f64; 2]; 3] {
        center.map(|c| [c - half, c + half])
    }

    #[test]
    fn probes_beside_the_authored_law_are_resolved_per_phase() {
        let t = npt();
        let s = t.spec();
        // The law at the theta = 0 meridian: radius R(q = z, z).
        for z in [1.1, 2.7, 4.3] {
            let q = z % s.pitch_mm;
            let r = t.radial_bounds([q; 2], [z; 2]).unwrap();
            assert!(r[1] - r[0] < 1e-9);
            let inside = t.classify_box(cube([r[0] - 1e-3, 0., z], 1e-4)).unwrap();
            let outside = t.classify_box(cube([r[1] + 1e-3, 0., z], 1e-4)).unwrap();
            assert_eq!(inside, Some(AuthoredBoxRelation::Inside), "z={z}");
            assert_eq!(outside, Some(AuthoredBoxRelation::Outside), "z={z}");
            // On the surface nothing is claimed.
            assert_eq!(t.classify_box(cube([r[0], 0., z], 1e-4)).unwrap(), None);
        }
        // The same law a quarter turn round: q = z - pitch / 4.
        let z = 3.;
        let q = (z - s.pitch_mm / 4.).rem_euclid(s.pitch_mm);
        let r = t.radial_bounds([q; 2], [z; 2]).unwrap()[0];
        assert_eq!(
            t.classify_box(cube([0., r - 1e-3, z], 1e-4)).unwrap(),
            Some(AuthoredBoxRelation::Inside)
        );
        assert_eq!(
            t.classify_box(cube([0., r + 1e-3, z], 1e-4)).unwrap(),
            Some(AuthoredBoxRelation::Outside)
        );
    }

    #[test]
    fn caps_axis_and_far_boxes() {
        let t = npt();
        // Core material on the axis, away from the caps.
        assert_eq!(
            t.classify_box(cube([0., 0., 3.], 0.5)).unwrap(),
            Some(AuthoredBoxRelation::Inside)
        );
        // A box touching a cap plane is not claimed inside.
        assert_eq!(t.classify_box(cube([0., 0., 0.], 0.1)).unwrap(), None);
        // Past either cap, or well beyond the major radius.
        for center in [[0., 0., -1.], [0., 0., 7.], [9., 0., 3.], [-6., -6., 3.]] {
            assert_eq!(
                t.classify_box(cube(center, 0.1)).unwrap(),
                Some(AuthoredBoxRelation::Outside),
                "{center:?}"
            );
        }
        assert_eq!(
            t.classify_box([[1., 0.], [0., 1.], [0., 1.]]),
            Err(Error::Invalid)
        );
    }
}
