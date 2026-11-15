//! Original-law phase moments; no derived literal circle is substituted.
use super::*;
#[derive(Clone, Copy, Debug)]
pub(super) struct Moments {
    first: I,
    second: I,
}
fn line(a: (I, I), b: (I, I)) -> Moments {
    let width = b.0 - a.0;
    Moments {
        first: width * (a.1 + b.1) / I::point(2.),
        second: width * (a.1.square() + a.1 * b.1 + b.1.square()) / I::point(3.),
    }
}
fn arc(center: I, sign: f64, radius: I, rounded: Rounded, k: I) -> Moments {
    // q=q0+R*sin(theta), a=A+s*R*cos(theta)-k*R*sin(theta).
    // Integrate a^j dq on [-beta,beta]. Odd terms vanish exactly;
    // all remaining terms are the elementary complete sin/cos integrals.
    // beta=pi/2-alpha, so sin(beta)=cos(alpha), cos(beta)=sin(alpha).
    let sine = rounded.cos;
    let cube = sine * sine.square();
    let cosine_square = rounded.beta + sine * rounded.sin;
    let sine_square_cosine = I::point(2.) * cube / I::point(3.);
    let cosine_cube = I::point(2.) * (sine - cube / I::point(3.));
    Moments {
        first: radius * (I::point(2.) * center * sine + I::point(sign) * radius * cosine_square),
        second: radius
            * (I::point(2.) * center.square() * sine
                + I::point(2. * sign) * center * radius * cosine_square
                + radius.square() * cosine_cube
                + k.square() * radius.square() * sine_square_cosine),
    }
}
impl Moments {
    pub(super) fn new(thread: &AuthoredThread) -> Result<Self, Error> {
        let r = I::point(thread.spec.major_radius_mm);
        let d = I::point(thread.depth());
        let p = I::point(thread.spec.pitch_mm);
        let k = I::point(thread.spec.radial_taper);
        let mut result = Self {
            first: I::point(0.),
            second: I::point(0.),
        };
        let mut add = |value: Self| {
            result.first = result.first + value.first;
            result.second = result.second + value.second;
        };
        match thread.spec.section {
            AuthoredSection::NptFlat { .. } => {
                for band in 0..3 {
                    add(line(
                        thread.endpoint(band, false)?,
                        thread.endpoint(band, true)?,
                    ));
                }
                let gap = p - (thread.endpoint(2, true)?.0 - thread.endpoint(0, false)?.0);
                add(Self {
                    first: gap * (r - d),
                    second: gap * (r - d).square(),
                });
            }
            AuthoredSection::BsptRounded { .. } => {
                let b = thread.rounded()?;
                add(arc(r - b.radius * b.normal, 1., b.radius, b, k));
                add(arc(r - d + b.radius * b.normal, -1., b.radius, b, k));
                for band in [1, 3] {
                    add(line(
                        thread.endpoint(band, false)?,
                        thread.endpoint(band, true)?,
                    ));
                }
            }
        }
        if !result.first.valid()
            || !result.second.valid()
            || result.first.lo <= 0.
            || result.second.lo <= 0.
        {
            return Err(Error::Unresolved);
        }
        Ok(result)
    }
    pub(super) fn area(self, thread: &AuthoredThread, z: f64) -> I {
        let shift = I::point(thread.spec.radial_taper) * I::point(z);
        pi() / I::point(thread.spec.pitch_mm)
            * (self.second
                + I::point(2.) * shift * self.first
                + shift.square() * I::point(thread.spec.pitch_mm))
    }
    pub(super) fn volume(self, thread: &AuthoredThread) -> I {
        let length = I::point(thread.spec.length_mm);
        let shift = I::point(thread.spec.radial_taper) * length;
        pi() * length / I::point(thread.spec.pitch_mm)
            * (self.second
                + shift * self.first
                + shift.square() * I::point(thread.spec.pitch_mm) / I::point(3.))
    }
}
