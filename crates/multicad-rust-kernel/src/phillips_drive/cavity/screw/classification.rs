//! Whole-box material classification; no midpoint or nominal frame admission.
use super::*;
use crate::head_profile::{Curve, Role};
use crate::helical_thread::BoxClassification as C;

fn profile(profile: &Profile, p: [I; 3]) -> C {
    profile_bounds(profile, p).unwrap_or(C::UnresolvedBoundary)
}
// None means the numerical/source enclosure cannot certify a predicate.
// Keep that failure separate from a proved geometric boundary intersection.
fn profile_bounds(profile: &Profile, p: [I; 3]) -> Option<C> {
    let precision = crate::head_profile::source_volume::precision_envelope(profile)?;
    let d = precision.deviation.hi;
    let top = I::point(profile.top()[1]);
    let end = I::point(profile.end()[1]);
    let rho = (p[0].square() + p[1].square()).sqrt();
    if p[2].hi < (top - I::point(d)).lo || p[2].lo > (end + I::point(d)).hi {
        return Some(C::Outside);
    }
    // An exposed meridian tube can move both axially and radially. Query the
    // entire expanded axial interval before shrinking/expanding radial bounds.
    let z = p[2] + I::new(-d, d);
    let mut radial = I::new(f64::INFINITY, f64::NEG_INFINITY);
    for segment in profile.segments() {
        let q = I::new(z.lo.max(segment.start[1]), z.hi.min(segment.end[1]));
        if q.lo > q.hi {
            continue;
        }
        let r = match segment.curve {
            Curve::Line if segment.start[1] == segment.end[1] => I::new(
                segment.start[0].min(segment.end[0]),
                segment.start[0].max(segment.end[0]),
            ),
            Curve::Line => {
                I::point(segment.start[0])
                    + (q - I::point(segment.start[1]))
                        / (I::point(segment.end[1]) - I::point(segment.start[1]))
                        * (I::point(segment.end[0]) - I::point(segment.start[0]))
            }
            Curve::Arc {
                center, radius_mm, ..
            } => {
                let radial =
                    (I::point(radius_mm).square() - (q - I::point(center[1])).square()).sqrt();
                if segment.role == Role::Crown {
                    I::point(center[0]) + radial
                } else {
                    I::point(center[0]) - radial
                }
            }
        };
        r.valid().then_some(())?;
        radial.lo = radial.lo.min(r.lo);
        radial.hi = radial.hi.max(r.hi);
    }
    (radial.valid() && rho.valid()).then_some(())?;
    Some(if rho.lo > (radial + I::point(d)).hi {
        C::Outside
    } else if p[2].lo > (top + I::point(d)).hi
        && p[2].hi < (end - I::point(d)).lo
        && rho.hi < (radial - I::point(d)).lo
    {
        C::Inside
    } else {
        C::UnresolvedBoundary
    })
}
impl PhillipsScrew {
    /// Classify all points of a closed source-coordinate box. This retains
    /// the original outer source and finite cutter, including source reserve.
    pub fn classify_box(&self, bounds: [[f64; 2]; 3]) -> Result<C, Error> {
        let p = bounds.map(|[lo, hi]| I::new(lo, hi));
        if p.into_iter().any(|i| !i.valid()) {
            return Err(Error::Invalid);
        }
        let thread = self
            .thread
            .classify_box(bounds)
            .map_err(|_| Error::Invalid)?;
        let mut head = p;
        head[2] = head[2] + I::point(self.spec.shank_length_mm);
        let head = profile(&self.profile, head);
        let outer = match (thread, head) {
            (C::Inside, _) | (_, C::Inside) => C::Inside,
            (C::Outside, C::Outside) => C::Outside,
            _ => C::UnresolvedBoundary,
        };
        if outer == C::Outside {
            return Ok(outer);
        }
        let mut cutter = p;
        cutter[2] =
            cutter[2] - (I::point(self.profile.top()[1]) - I::point(self.spec.shank_length_mm));
        let cavity = classify::enclosure(self.cavity.source(), cutter)?;
        Ok(match (outer, cavity) {
            (_, Classification::Inside) => C::Outside,
            (C::Inside, Classification::Outside) => C::Inside,
            _ => C::UnresolvedBoundary,
        })
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn actual_metric_void_material_and_whole_box_boundary_remain_distinct() {
        for index in 0..3 {
            for left in [false, true] {
                let screw = super::super::tests::fixture(index, 0.2, left);
                let top = screw.head_top_z();
                for (point, expected) in [
                    ([0., 0., 4.], C::Inside),
                    ([0., 0., top + 0.1], C::Outside),
                    ([20., 0., 0.], C::Outside),
                    ([0., 0., 9.], C::Outside),
                ] {
                    assert_eq!(
                        screw.classify_box(point.map(|v| [v - 1e-8, v + 1e-8])),
                        Ok(expected)
                    );
                }
                assert_eq!(
                    screw.classify_box([[0., 0.], [0., 0.], [top - 1e-7, top + 1e-7]]),
                    Ok(C::UnresolvedBoundary)
                );
                assert_eq!(
                    screw.classify_box([[f64::NAN, 0.], [0., 0.], [0., 0.]]),
                    Err(Error::Invalid)
                );
                let r = screw.profile().top()[0];
                assert_eq!(
                    screw.classify_box([
                        [r - 0.001, r - 0.0009],
                        [0., 0.],
                        [top + 0.0001, top + 0.0002]
                    ]),
                    Ok(C::Inside)
                );
            }
        }
    }
}

#[cfg(test)]
mod numeric_tests {
    use super::*;
    use crate::head_profile::TransportError;
    #[test]
    fn actual_bearing_neck_and_transported_source_uncertainty_remain_explicit() {
        let screw = super::super::tests::fixture(1, 1., false);
        let source = screw.profile();
        for z in [0., 0.05, 0.1, 0.5] {
            assert_eq!(
                profile(source, [I::point(0.), I::point(0.), I::point(z)]),
                C::Inside
            );
        }
        let loose = Profile::from_transported_segments(
            &source.segments().collect::<Vec<_>>(),
            TransportError {
                numeric_conversion_mm: 0.005,
                arithmetic_mm: 0.,
            },
        )
        .unwrap();
        assert_eq!(profile(&loose, [I::point(0.); 3]), C::UnresolvedBoundary);
        assert_eq!(
            profile(source, [I::new(-f64::MAX, f64::MAX); 3]),
            C::UnresolvedBoundary
        );
    }
}
