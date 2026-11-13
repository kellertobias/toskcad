//! Authored/retained source mass bounds, never independent file readback.
use super::*;
use crate::tapered_thread::interval::I;
#[derive(Clone, Copy, Debug)]
pub(crate) struct PrecisionEnvelope {
    pub volume: I,
    pub radial: I,
    pub height: I,
    pub deviation: I,
}
/// The established finite-meridian tube bound, shared with physical head
/// readback. Its positional source reserve does not certify surface area.
pub(crate) fn precision_envelope(profile: &Profile) -> Option<PrecisionEnvelope> {
    let error = profile.rounding_error_mm();
    let mut length =
        I::point(profile.top()[0]) + I::point(profile.end()[0]) + I::point(profile.end()[1])
            - I::point(profile.top()[1]);
    let mut radial = I::point(0.);
    let mut count = 3.;
    for segment in profile.segments() {
        let part = match segment.curve {
            Curve::Line => {
                radial.hi = radial.hi.max(segment.start[0]).max(segment.end[0]);
                ((I::point(segment.end[0]) - I::point(segment.start[0])).square()
                    + (I::point(segment.end[1]) - I::point(segment.start[1])).square())
                .sqrt()
            }
            Curve::Arc {
                radius_mm,
                sweep,
                center,
            } => {
                if radius_mm <= 32. * error {
                    return None;
                }
                radial.hi = radial
                    .hi
                    .max((I::point(center[0].abs()) + I::point(radius_mm)).hi);
                I::point(radius_mm) * I::point(sweep.abs())
            }
        };
        length = length + part;
        count += 1.;
    }
    // r>32e bounds inverse angular conditioning. A complete meridian
    // displacement below16e covers support/trim rounding. A finite planar
    // curve's d-tube area is at most2Ld+pi*d² per segment; revolution adds
    // at most2pi(R+d). This includes concave necks and finite cap seams.
    let d = I::point(16.) * I::point(error);
    let pi = I::new(PI.next_down(), PI.next_up());
    let r = radial + d;
    let volume =
        I::point(2.) * pi * r * (I::point(2.) * length * d + I::point(count) * pi * d.square());
    let height = I::point(profile.end()[1]) - I::point(profile.top()[1]) + I::point(2.) * d;
    Some(PrecisionEnvelope {
        volume,
        radial: r,
        height,
        deviation: d,
    })
}
fn cube(x: I) -> I {
    x * x * x
}
pub(crate) fn angle(segment: Segment, center: [f64; 2]) -> I {
    let x = I::point(segment.start[0]) - I::point(center[0]);
    let z = I::point(segment.start[1]) - I::point(center[1]);
    let pi = I::new(PI.next_down(), PI.next_up());
    let a = if x.lo > 0. {
        (z / x).atan()
    } else if x.hi < 0. {
        (z / x).atan() - pi
    } else {
        // The qualified circle has r>32e. A start direction whose X
        // enclosure meets zero is therefore strictly below the center on
        // these retained Crown/Neck domains; division by -Z is conditioned.
        (-x / z).atan() - pi / I::point(2.)
    };
    // Include the actual floating atan2 evaluation used by Segment::point,
    // independently of the interval atan endpoint library allowance.
    a + I::new(-16. * f64::EPSILON, 16. * f64::EPSILON)
}
/// Encloses the nominal actual supports and their exposed finite source
/// precision. No mesh, authoring hint or caller tolerance changes this bound.
pub(crate) fn volume_bounds(profile: &Profile) -> Option<I> {
    let precision = precision_envelope(profile)?;
    let pi = I::new(PI.next_down(), PI.next_up());
    let mut volume = I::point(0.);
    for segment in profile.segments() {
        let v = match segment.curve {
            Curve::Line => {
                let a = I::point(segment.start[0]);
                let b = I::point(segment.end[0]);
                pi * (I::point(segment.end[1]) - I::point(segment.start[1]))
                    * (a.square() + a * b + b.square())
                    / I::point(3.)
            }
            Curve::Arc {
                center,
                radius_mm,
                sweep,
            } => {
                let r = I::point(radius_mm);
                let c = I::point(center[0]);
                let a = angle(segment, center);
                let sweep = I::point(sweep);
                if center[0] == 0. {
                    let dz = I::point(2.)
                        * r
                        * (a + sweep / I::point(2.)).cos()
                        * (sweep / I::point(2.)).sin();
                    pi * dz
                        * ((r * a.cos()).square() - r * a.sin() * dz - dz.square() / I::point(3.))
                } else {
                    let primitive = |t: I| {
                        let s = t.sin();
                        r * c.square() * s
                            + c * r.square() * (t + (I::point(2.) * t).sin() / I::point(2.))
                            + cube(r) * (s - cube(s) / I::point(3.))
                    };
                    pi * (primitive(a + sweep) - primitive(a))
                }
            }
        };
        volume = volume + v;
    }
    let result = volume + I::new(-precision.volume.hi, precision.volume.hi);
    (result.valid() && result.lo > 0.).then_some(result)
}

#[cfg(test)]
mod tests;
