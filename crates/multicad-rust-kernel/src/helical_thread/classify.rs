//! Strict membership of an entire coordinate box in the original source.
use super::*;
use crate::tapered_thread::interval::I;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BoxClassification {
    Inside,
    Outside,
    UnresolvedBoundary,
}

impl Thread {
    /// Both one-sided literal support laws at finite phase joins. The caller
    /// supplies a closed subinterval of one original axial period.
    pub(crate) fn radius_at_phase_interval(&self, phase: I) -> Option<I> {
        radius(self, phase)
    }
}

// Keep both one-sided intervals at the angular seam instead of replacing a
// tiny seam-crossing box by an entire turn. The retained rounded helical rate
// is applied to BOTH intervals; it is never replaced by an exact lead.
fn angle(x: I, y: I) -> [Option<I>; 2] {
    let pi = I::new(PI.next_down(), PI.next_up());
    let half = pi / I::point(2.);
    let tau = pi * I::point(2.);
    let a = if x.lo > 0. {
        (y / x).atan()
    } else if x.hi < 0. && y.lo >= 0. {
        (y / x).atan() + pi
    } else if x.hi < 0. && y.hi <= 0. {
        (y / x).atan() - pi
    } else if x.hi < 0. {
        // The negative-X seam is continuous in the wrapped [0,2pi] chart.
        (y / x).atan() + pi
    } else if y.lo > 0. {
        half - (x / y).atan()
    } else if y.hi < 0. {
        -half - (x / y).atan()
    } else {
        return [Some(I::new(0., tau.hi)), None];
    };
    // Enclose the actual atan2 evaluator independently of interval atan's
    // endpoint library allowance, as for retained meridian source proofs.
    let a = a + I::new(-16. * f64::EPSILON, 16. * f64::EPSILON);
    if a.lo >= 0. {
        [Some(a), None]
    } else if a.hi < 0. {
        [Some(a + tau), None]
    } else {
        [Some(I::new(0., a.hi)), Some(I::new(a.lo, 0.) + tau)]
    }
}
fn radius(thread: &Thread, phase: I) -> Option<I> {
    let mut result = I::new(f64::INFINITY, f64::NEG_INFINITY);
    for (band, pair) in thread.breaks.windows(2).enumerate() {
        let p = I::new(phase.lo.max(pair[0]), phase.hi.min(pair[1]));
        if p.lo > p.hi {
            continue;
        }
        let r = match band {
            1 => {
                I::point(thread.spec.major_radius_mm)
                    + I::point(thread.flanks[0]) * (p - I::point(thread.breaks[1]))
            }
            3 => {
                I::point(thread.spec.major_radius_mm)
                    + I::point(thread.flanks[1]) * (p - I::point(thread.breaks[4]))
            }
            2 => match thread.root {
                Some(root) => {
                    I::point(root.center_radius)
                        - (I::point(root.radius).square() - (p - I::point(root.phase)).square())
                            .sqrt()
                }
                None => I::point(thread.spec.major_radius_mm) - I::point(thread.spec.depth_mm),
            },
            _ => I::point(thread.spec.major_radius_mm),
        };
        // Source admission is not a substitute for a directed numeric proof.
        r.valid().then_some(())?;
        result.lo = result.lo.min(r.lo);
        result.hi = result.hi.max(r.hi);
    }
    result.valid().then_some(result)
}
fn phase_radius(thread: &Thread, q: I, quotient: I) -> Option<I> {
    let result = if quotient.width() >= 1. {
        radius(thread, I::new(0., thread.spec.pitch_mm))?
    } else {
        let mut r = I::new(f64::INFINITY, f64::NEG_INFINITY);
        let first = quotient.lo.floor() as i64 - 1;
        let last = quotient.hi.floor() as i64 + 1;
        for k in first..=last {
            let offset = I::point(k as f64) * I::point(thread.spec.pitch_mm);
            // Retain both global source-station evaluation paths. A
            // rounded k*p+break must not silently become an exact shift.
            let mut shift_error = 0_f64;
            for knot in thread.breaks {
                let exact = offset + I::point(knot);
                for literal in [
                    (k as f64).mul_add(thread.spec.pitch_mm, knot),
                    k as f64 * thread.spec.pitch_mm + knot,
                ] {
                    shift_error = shift_error.max((I::point(literal) - exact).magnitude());
                }
            }
            let phase = q - offset + I::new(-shift_error, shift_error);
            let phase = I::new(phase.lo.max(0.), phase.hi.min(thread.spec.pitch_mm));
            if phase.lo <= phase.hi {
                let piece = radius(thread, phase)?;
                r.lo = r.lo.min(piece.lo);
                r.hi = r.hi.max(piece.hi);
            }
        }
        r
    };
    result.valid().then_some(result)
}
impl Thread {
    /// Classify every point of a closed local coordinate box. Surface or
    /// phase uncertainty remains explicit; this is not a distance query.
    pub fn classify_box(&self, bounds: [[f64; 2]; 3]) -> Result<BoxClassification, ThreadError> {
        use BoxClassification::*;
        let p = bounds.map(|[lo, hi]| I::new(lo, hi));
        if p.into_iter().any(|i| !i.valid()) {
            return Err(ThreadError::Invalid);
        }
        let [x, y, z] = p;
        let rho = (x.square() + y.square()).sqrt();
        if z.hi < 0. || z.lo > self.spec.length_mm || rho.lo > self.spec.major_radius_mm {
            return Ok(Outside);
        }
        // Axial rejection above bounds the only boxes for which this is used.
        // Each of at most two angular pieces needs at most one whole period.
        let radii = angle(x, y).into_iter().flatten().try_fold(
            I::new(f64::INFINITY, f64::NEG_INFINITY),
            |radii, angle| {
                let q = z - I::point(self.helical_rate()) * angle;
                let quotient = q / I::point(self.spec.pitch_mm);
                if !q.valid() || !quotient.valid() {
                    return None;
                }
                let piece = phase_radius(self, q, quotient)?;
                Some(I::new(radii.lo.min(piece.lo), radii.hi.max(piece.hi)))
            },
        );
        if !rho.valid() {
            return Ok(UnresolvedBoundary);
        }
        // Failed piece proofs remain unresolved. Every accepted piece is
        // independently checked; angle always supplies at least one interval.
        Ok(radii.map_or(UnresolvedBoundary, |radii| {
            if rho.lo > radii.hi {
                Outside
            } else if z.lo > 0. && z.hi < self.spec.length_mm && rho.hi < radii.lo {
                Inside
            } else {
                UnresolvedBoundary
            }
        }))
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    fn fixture(rounded: bool, left_hand: bool, starts: u32) -> Thread {
        Thread::new(Spec {
            major_radius_mm: 2.,
            pitch_mm: 0.7,
            length_mm: 4.2,
            depth_mm: if rounded {
                let [lower, upper] = [30_f64.to_radians(), 35_f64.to_radians()];
                (0.7 - 0.1) / (lower.tan() + upper.tan())
                    - 0.06 * (((upper - lower) * 0.5).cos() / ((upper + lower) * 0.5).sin() - 1.)
            } else {
                0.35
            },
            crest_width_mm: 0.1,
            lower_flank_deg: 30.,
            upper_flank_deg: 35.,
            starts,
            left_hand,
            root_radius_mm: rounded.then_some(0.06),
        })
        .unwrap()
    }
    #[test]
    fn small_boxes_across_both_angular_seams_keep_the_actual_wrapped_source_law() {
        for rounded in [false, true] {
            for left in [false, true] {
                for starts in [1, 2] {
                    let thread = fixture(rounded, left, starts);
                    for theta in [0., PI] {
                        for z in [0.23, 0.47, 1.13] {
                            let radius = thread.radius(theta, z);
                            for (offset, occupied) in [(-0.001, true), (0.001, false)] {
                                let centre = [
                                    (radius + offset) * theta.cos(),
                                    (radius + offset) * theta.sin(),
                                    z,
                                ];
                                let bounds = centre.map(|v| [v - 1e-6, v + 1e-6]);
                                assert_eq!(
                                    thread.classify_box(bounds),
                                    Ok(if occupied {
                                        BoxClassification::Inside
                                    } else {
                                        BoxClassification::Outside
                                    }),
                                    "rounded={rounded},left={left},starts={starts},theta={theta},z={z}"
                                );
                                for corner in 0..8 {
                                    assert_eq!(
                                        thread.contains(
                                            std::array::from_fn(
                                                |i| bounds[i][usize::from(corner & (1 << i) != 0)]
                                            ),
                                            0.
                                        ),
                                        occupied
                                    );
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    #[test]
    fn original_flat_and_rounded_boxes_preserve_inside_outside_and_boundary_uncertainty() {
        use BoxClassification::*;
        for rounded in [false, true] {
            for left in [false, true] {
                for starts in [1, 2] {
                    let t = fixture(rounded, left, starts);
                    assert_eq!(t.classify_box([[0., 0.], [0., 0.], [1., 2.]]), Ok(Inside));
                    assert_eq!(t.classify_box([[3., 4.], [0., 0.], [1., 2.]]), Ok(Outside));
                    assert_eq!(
                        t.classify_box([[0., 0.], [0., 0.], [-2., -1.]]),
                        Ok(Outside)
                    );
                    assert_eq!(t.classify_box([[0., 0.], [0., 0.], [5., 6.]]), Ok(Outside));
                    assert_eq!(
                        t.classify_box([[0., 0.], [0., 0.], [0., 0.]]),
                        Ok(UnresolvedBoundary)
                    );
                    assert_eq!(
                        t.classify_box([[-3., 3.], [-3., 3.], [-1., 5.]]),
                        Ok(UnresolvedBoundary)
                    );
                    assert_eq!(
                        t.classify_box([[0., 0.], [0., 0.], [-f64::MAX, f64::MAX]]),
                        Ok(UnresolvedBoundary)
                    );
                    assert_eq!(
                        t.classify_box([[-f64::MAX, f64::MAX]; 3]),
                        Ok(UnresolvedBoundary)
                    );
                    assert_eq!(
                        t.classify_box([[-f64::MAX, f64::MAX], [-f64::MAX, f64::MAX], [1., 2.]]),
                        Ok(UnresolvedBoundary)
                    );

                    assert_eq!(
                        t.classify_box([[f64::NAN, 1.], [0., 0.], [1., 1.]]),
                        Err(ThreadError::Invalid)
                    );
                    assert_eq!(
                        t.classify_box([[1., 0.], [0., 0.], [1., 1.]]),
                        Err(ThreadError::Invalid)
                    );
                    // All these separated small boxes now classify strictly,
                    // including both seams. Check independent nominal corners.
                    for n in 0..48 {
                        let theta = n as f64 * TAU / 48.;
                        let z = 0.9 + n as f64 * 0.0123;
                        let r = t.radius(theta, z);
                        for offset in [-0.04, 0.04] {
                            let p = [(r + offset) * theta.cos(), (r + offset) * theta.sin(), z];
                            let b = p.map(|v| [v - 1e-7, v + 1e-7]);
                            let class = t.classify_box(b).unwrap();
                            assert_eq!(class, if offset < 0. { Inside } else { Outside });
                            for bits in 0..8 {
                                let q = std::array::from_fn(|i| b[i][(bits >> i) & 1]);
                                assert_eq!(t.contains(q, 0.), class == Inside);
                            }
                        }
                    }
                }
            }
        }
    }
}
#[cfg(test)]
mod precision_tests {
    use super::*;
    #[test]
    fn genuine_public_tiny_root_never_panics_or_mutates_the_source() {
        let a = 1_f64.to_radians();
        let root = 2_f64.powi(-52).next_up();
        let thread = Thread::new(Spec {
            major_radius_mm: 20.,
            pitch_mm: 0.5,
            length_mm: 4.,
            crest_width_mm: 0.0625,
            lower_flank_deg: 1.,
            upper_flank_deg: 1.,
            depth_mm: (0.5 - 0.0625) / (2. * a.tan()) - root * (1. / a.sin() - 1.),
            root_radius_mm: Some(root),
            starts: 1,
            left_hand: false,
        })
        .unwrap();
        let original = thread;
        assert_eq!(
            thread.classify_box([[0., 0.], [0., 0.], [1., 2.]]),
            Ok(BoxClassification::Inside)
        );
        for z in thread.profile_breaks() {
            assert!(thread.classify_box([[8., 8.], [0., 0.], [z, z]]).is_ok());
        }
        assert_eq!(thread, original);
    }
}
