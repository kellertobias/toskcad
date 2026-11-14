//! Actual canonical rounded metric density in the original one-lead q/z chart.
//! Root coefficients enclose the owned shared Q(sqrt3) expressions. Neither a
//! legacy Thread nor a tolerance-fitted circle supplies any geometry.
use super::{
    series::{S, sin_cos},
    *,
};
use crate::authored_metric_thread::{AuthoredMetricThread, Knot};
#[derive(Clone, Copy)]
enum Kind {
    Wall,
    Disk {
        station: I,
        body_z: I,
        outer: Option<I>,
    },
}
struct Part<'a> {
    source: &'a AuthoredMetricThread,
    band: usize,
    q0: I,
    dq: I,
    r0: I,
    dr: I,
    root_radius: I,
    root_center: I,
    kind: Kind,
    shift: I,
}
fn interval(a: [f64; 2]) -> I {
    I::new(a[0], a[1])
}
fn tau() -> I {
    I::new(
        std::f64::consts::TAU.next_down(),
        std::f64::consts::TAU.next_up(),
    )
}
fn root3() -> I {
    I::point(3.).sqrt()
}
impl integrate::Density for Part<'_> {
    fn field(&self, u: I) -> Result<[S; 4], Error> {
        let t = S::variable(u);
        let mut root = None;
        let r = if self.band == 2 {
            let x = t.sub(S::constant(I::point(0.5))).scale(root3());
            let h = S::constant(I::point(1.)).sub(x.square()).sqrt()?;
            root = Some((x, h));
            S::constant(self.root_center).sub(h.scale(self.root_radius))
        } else {
            S::constant(self.r0).add(t.scale(self.dr))
        };
        let spec = self.source.spec();
        let lead = I::point(spec.pitch_mm) * I::point(f64::from(spec.starts));
        let abs_s = lead / tau();
        let s = if spec.left_hand { -abs_s } else { abs_s };
        let theta = |z: I| sin_cos((z - self.q0) / s, -self.dq / s, u);
        let values = match self.kind {
            Kind::Wall => {
                // On the circle H²+X²=1 and H>0 exactly. This homogeneous
                // identity avoids dependency from 1+(X/H)², without changing J.
                let metric = if let Some((x, h)) = root {
                    r.scale(I::point(1.) / s)
                        .square()
                        .add(x.square())
                        .sqrt()?
                        .mul(h.positive_reciprocal()?)
                } else {
                    let m = if self.band == 1 {
                        -root3()
                    } else if self.band == 3 {
                        root3()
                    } else {
                        I::point(0.)
                    };
                    r.scale(I::point(1.) / s)
                        .square()
                        .scale(I::point(1.) + m.square())
                        .add(S::constant(m.square()))
                        .sqrt()?
                }
                .scale(self.dq);
                let length = I::point(spec.length_mm);
                let (sin0, cos0) = theta(I::point(0.));
                let (sin1, cos1) = theta(length);
                let area = metric.scale(length);
                let xy = metric.mul(r).scale(s);
                [
                    area,
                    xy.mul(sin1.sub(sin0)),
                    xy.mul(cos0.sub(cos1)),
                    area.scale(length / I::point(2.) + self.shift),
                ]
            }
            Kind::Disk {
                station,
                body_z,
                outer,
            } => {
                let (sine, cosine) = theta(station);
                let squared = outer.map_or_else(
                    || r.square(),
                    |v| S::constant(v).sub(r).mul(S::constant(v).add(r)),
                );
                let area = squared.scale(self.dq / (I::point(2.) * abs_s));
                // The complete outer circle has zero XY moments over the full
                // one-lead domain. The actual inner profile moment is retained.
                let radial = r.square().mul(r).scale(
                    self.dq / (I::point(3.) * abs_s)
                        * I::point(if outer.is_some() { -1. } else { 1. }),
                );
                [
                    area,
                    radial.mul(cosine),
                    radial.mul(sine),
                    area.scale(body_z),
                ]
            }
        };
        let mut result = [S::constant(I::point(0.)); 4];
        for (out, value) in result.iter_mut().zip(values) {
            *out = value.checked()?;
        }
        Ok(result)
    }
}
#[derive(Clone, Copy)]
pub(super) struct Chart {
    pub shift: I,
    pub cap_body_z: Option<I>,
    pub outer: Option<I>,
}
pub(super) fn evaluate(
    source: &AuthoredMetricThread,
    face: FaceKey,
    chart: Chart,
    options: FacePropertyOptions,
    work: &mut integrate::Work<'_>,
) -> Result<[I; 4], Error> {
    let spec = source.spec();
    let Chart {
        shift,
        cap_body_z,
        outer,
    } = chart;
    let kind = match face {
        FaceKey::StartCap => Kind::Disk {
            station: I::point(0.),
            body_z: cap_body_z.unwrap_or(shift),
            outer,
        },
        FaceKey::EndCap => Kind::Disk {
            station: I::point(spec.length_mm),
            body_z: cap_body_z.unwrap_or(I::point(spec.length_mm) + shift),
            outer,
        },
        FaceKey::Thread { start, .. } if start < spec.starts && outer.is_none() => Kind::Wall,
        _ => return Err(Error::Invalid),
    };
    let keys = [
        Knot::Begin,
        Knot::CrestEnd,
        Knot::UpperRoot,
        Knot::LowerRoot,
        Knot::CrestStart,
        Knot::End,
    ];
    let knots = keys.map(|k| source.knot(k));
    let radius = interval(source.root_radius_mm());
    // Both intervals enclose their actual owned expressions, rather than
    // selecting independent rounded center/radius inputs.
    let center = interval(source.minimum_radius_mm()) + radius;
    let mut parts = Vec::with_capacity(5 * spec.starts as usize);
    for period in 0..spec.starts {
        for band in 0..5 {
            work.checkpoint()?;
            let group = (period + [0, 0, 0, 1, 1][band]) % spec.starts;
            let start = if spec.left_hand {
                group
            } else {
                (spec.starts - group) % spec.starts
            };
            let role = [
                ThreadBandRole::Crest,
                ThreadBandRole::Upper,
                ThreadBandRole::Root,
                ThreadBandRole::Lower,
                ThreadBandRole::Crest,
            ][band];
            if matches!(kind, Kind::Wall) && face != (FaceKey::Thread { start, band: role }) {
                continue;
            }
            let q0 = interval(knots[band].phase_mm);
            let r0 = interval(knots[band].radius_mm);
            // Equal authored crest radii are one constant expression, not two
            // interval endpoints to subtract and then use as an artificial slope.
            let dr = if band == 0 || band == 4 {
                I::point(0.)
            } else {
                interval(knots[band + 1].radius_mm) - r0
            };
            parts.push(Part {
                source,
                band,
                q0: q0 + I::point(f64::from(period)) * I::point(spec.pitch_mm),
                dq: interval(knots[band + 1].phase_mm) - q0,
                r0,
                dr,
                root_radius: radius,
                root_center: center,
                kind,
                shift,
            });
        }
    }
    let densities: Vec<&dyn integrate::Density> =
        parts.iter().map(|p| p as &dyn integrate::Density).collect();
    let extent = source
        .maximum_radius_mm()
        .max(spec.length_mm + shift.magnitude())
        .max(1.);
    integrate::parts(&densities, extent, options, work)
}
impl AuthoredMetricThread {
    /// Bounded local area and first moments of the actual finite source.
    /// This does not imply a finite graph or area invariance under placement.
    pub fn local_face_properties(
        &self,
        face: FaceKey,
        options: FacePropertyOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<FaceProperties, Error> {
        let mut work = integrate::Work::new(options, interrupt)?;
        let values = evaluate(
            self,
            face,
            Chart {
                shift: I::point(0.),
                cap_body_z: None,
                outer: None,
            },
            options,
            &mut work,
        )?;
        let result = finish(face, values, options)?;
        work.checkpoint()?;
        Ok(result)
    }
}
