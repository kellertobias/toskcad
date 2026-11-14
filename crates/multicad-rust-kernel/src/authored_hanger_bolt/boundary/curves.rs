//! Bounds of the original finite source laws; no mesh point or fitted circle
//! supplies a curve. Endpoint identity is the parent graph's symbolic token.
use super::super::classify::metric_error;
use super::*;
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CurveBounds {
    pub point_mm: [[f64; 2]; 3],
    pub first_derivative_mm: [[f64; 2]; 3],
    pub second_derivative_mm: [[f64; 2]; 3],
}
fn tau() -> I {
    I::new(
        std::f64::consts::TAU.next_down(),
        std::f64::consts::TAU.next_up(),
    )
}
fn enclosure(v: [f64; 2]) -> I {
    I::new(v[0], v[1])
}
struct Polar {
    r: I,
    dr: I,
    ddr: I,
    theta: I,
    dt: I,
    z: I,
    dz: I,
}
impl Polar {
    fn finish(self) -> Result<CurveBounds, Error> {
        let Self {
            r,
            dr,
            ddr,
            theta,
            dt,
            z,
            dz,
        } = self;
        let (c, s) = (theta.cos(), theta.sin());
        let a = ddr - r * dt.square();
        let b = I::point(2.) * dr * dt;
        let values = [
            [r * c, r * s, z],
            [dr * c - r * s * dt, dr * s + r * c * dt, dz],
            [a * c - b * s, a * s + b * c, I::point(0.)],
        ];
        for row in &values {
            checked(row)?;
        }
        Ok(CurveBounds {
            point_mm: values[0].map(pair),
            first_derivative_mm: values[1].map(pair),
            second_derivative_mm: values[2].map(pair),
        })
    }
}
impl Boundary {
    fn dimensions(&self, c: Component) -> (I, I, I, bool) {
        match c {
            Component::Metric => {
                let s = self.source.metric().spec();
                (
                    I::point(s.pitch_mm),
                    I::point(s.pitch_mm) * I::point(f64::from(s.starts)),
                    I::point(s.length_mm),
                    s.left_hand,
                )
            }
            Component::Wood => {
                let s = self.source.wood().spec();
                (
                    I::point(s.pitch_mm),
                    I::point(s.pitch_mm) * I::point(f64::from(s.starts)),
                    I::point(s.length_mm),
                    s.left_hand,
                )
            }
        }
    }
    fn z(&self, c: Component, z: I) -> I {
        match c {
            Component::Wood => z,
            Component::Metric => {
                z - I::point(self.source.metric().full_length_mm())
                    - I::point(self.source.spec().plain_length_mm)
            }
        }
    }
    fn knot(&self, c: Component, k: Knot) -> (I, I) {
        let i = match k {
            Knot::CrestEnd => 0,
            Knot::UpperRoot => 1,
            Knot::LowerRoot => 2,
            Knot::CrestStart => 3,
        };
        match c {
            Component::Wood => {
                let b = self.source.wood().bands()[i];
                (I::point(b.phase_mm[1]), I::point(b.radius_mm[1]))
            }
            Component::Metric => {
                use crate::authored_metric_thread::Knot as K;
                let p = self
                    .source
                    .metric()
                    .knot([K::CrestEnd, K::UpperRoot, K::LowerRoot, K::CrestStart][i]);
                (enclosure(p.phase_mm), enclosure(p.radius_mm))
            }
        }
    }
    fn cap_profile(
        &self,
        c: Component,
        band: ThreadBandRole,
        u: I,
    ) -> Result<(I, I, I, I, I), Error> {
        let index = match band {
            ThreadBandRole::Upper => 1,
            ThreadBandRole::Root => 2,
            ThreadBandRole::Lower => 3,
            ThreadBandRole::Crest => 4,
        };
        match c {
            Component::Wood => {
                let bands = self.source.wood().bands();
                let b = bands[index];
                let q0 = I::point(b.phase_mm[0]);
                let dq = if band == ThreadBandRole::Crest {
                    (I::point(self.source.wood().spec().pitch_mm) - q0)
                        + I::point(bands[0].phase_mm[1])
                } else {
                    I::point(b.phase_mm[1]) - q0
                };
                let r0 = I::point(b.radius_mm[0]);
                let dr = I::point(b.radius_mm[1]) - r0;
                Ok((q0, dq, r0 + dr * u, dr, I::point(0.)))
            }
            Component::Metric => {
                let t = self.source.metric();
                if band == ThreadBandRole::Crest {
                    let (q0, r) = self.knot(c, Knot::CrestStart);
                    return Ok((
                        q0,
                        I::point(t.spec().crest_width_mm),
                        r,
                        I::point(0.),
                        I::point(0.),
                    ));
                }
                let p = t
                    .profile_bounds(index, [u.lo, u.hi])
                    .map_err(metric_error)?;
                let (a, _) = self.knot(c, KNOTS[index - 1]);
                let (b, _) = self.knot(c, KNOTS[index]);
                let dq = b - a;
                Ok((
                    a,
                    dq,
                    enclosure(p.radius_mm),
                    enclosure(p.slope) * dq,
                    enclosure(p.second_derivative_per_mm) * dq.square(),
                ))
            }
        }
    }
    /// Enclose actual source position and derivatives with respect to the
    /// finite curve's normalized parameter. These are local chart quantities;
    /// placement/renderer tolerance and literal rounding remain separate.
    pub fn curve_bounds(
        &self,
        key: CurveKey,
        parameter: [f64; 2],
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<CurveBounds, Error> {
        checkpoint(interrupt)?;
        let u = enclosure(parameter);
        if !u.valid() || u.lo < 0. || u.hi > 1. || !self.curves.iter().any(|c| c.key == key) {
            return Err(Error::Invalid);
        }
        let polar = match key {
            CurveKey::PlainCircle { end } => Polar {
                r: I::point(self.source.spec().plain_radius_mm),
                dr: I::point(0.),
                ddr: I::point(0.),
                theta: tau() * u,
                dt: tau(),
                z: I::point(if end == End::Start {
                    -self.source.spec().plain_length_mm
                } else {
                    0.
                }),
                dz: I::point(0.),
            },
            CurveKey::Rail {
                component,
                period,
                knot,
            } => {
                let (pitch, lead, length, left) = self.dimensions(component);
                let (phase, r) = self.knot(component, knot);
                let rate = tau() / lead * I::point(if left { -1. } else { 1. });
                let z = length * u;
                Polar {
                    r,
                    dr: I::point(0.),
                    ddr: I::point(0.),
                    theta: rate * (z - phase - I::point(f64::from(period)) * pitch),
                    dt: rate * length,
                    z: self.z(component, z),
                    dz: length,
                }
            }
            CurveKey::Cap {
                component,
                end,
                period,
                band,
            } => {
                let (pitch, lead, length, left) = self.dimensions(component);
                let (q0, dq, r, dr, ddr) = self.cap_profile(component, band, u)?;
                let z = if end == End::Start {
                    I::point(0.)
                } else {
                    length
                };
                let rate = tau() / lead * I::point(if left { -1. } else { 1. });
                Polar {
                    r,
                    dr,
                    ddr,
                    theta: rate * (z - q0 - dq * u - I::point(f64::from(period)) * pitch),
                    dt: -rate * dq,
                    z: self.z(component, z),
                    dz: I::point(0.),
                }
            }
        };
        let result = polar.finish()?;
        checkpoint(interrupt)?;
        Ok(result)
    }
    pub fn vertex_bounds(
        &self,
        key: VertexKey,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<[[f64; 2]; 3], Error> {
        checkpoint(interrupt)?;
        if !self.vertices.contains(&key) {
            return Err(Error::Invalid);
        }
        let bounds = match key {
            VertexKey::Plain { end } => [
                I::point(self.source.spec().plain_radius_mm),
                I::point(0.),
                I::point(if end == End::Start {
                    -self.source.spec().plain_length_mm
                } else {
                    0.
                }),
            ],
            VertexKey::Thread {
                component,
                end,
                period,
                knot,
            } => {
                let (pitch, lead, length, left) = self.dimensions(component);
                let (q, r) = self.knot(component, knot);
                let z = if end == End::Start {
                    I::point(0.)
                } else {
                    length
                };
                let angle = (z - q - I::point(f64::from(period)) * pitch) * tau() / lead
                    * I::point(if left { -1. } else { 1. });
                [r * angle.cos(), r * angle.sin(), self.z(component, z)]
            }
        };
        checked(&bounds)?;
        checkpoint(interrupt)?;
        Ok(bounds.map(pair))
    }
}

fn checked(values: &[I]) -> Result<(), Error> {
    if values.iter().any(|v| !v.valid()) {
        Err(Error::Unresolved)
    } else {
        Ok(())
    }
}
#[cfg(test)]
mod numeric_tests {
    use super::*;
    #[test]
    fn numeric_enclosure_and_actual_metric_domain_failures_remain_explicit() {
        assert_eq!(checked(&[I::point(f64::NAN)]), Err(Error::Unresolved));
        assert_eq!(checked(&[I::new(2., 1.)]), Err(Error::Unresolved));
        assert!(checked(&[I::point(0.)]).is_ok());
        let source = Hanger::new(super::super::super::tests::fixture(false, false), None).unwrap();
        let b = source.boundary(BoundaryOptions::default(), None).unwrap();
        assert!(matches!(
            b.cap_profile(Component::Metric, ThreadBandRole::Root, I::new(-0.1, 0.)),
            Err(Error::Invalid)
        ));
    }
}
