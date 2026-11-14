//! The one-lead q/z chart covers the same physical source as q/v modulo
//! (Q,v=1)=(Q+lead,v=0). No finite mesh or legacy Thread supplies a domain.
use super::{
    series::{S, sin_cos},
    *,
};
#[derive(Clone, Copy)]
enum Kind {
    Wall,
    Cone,
    Disk { z: f64, outer: Option<I> },
}
#[derive(Clone, Copy)]
pub(super) struct Part {
    kind: Kind,
    q0: I,
    dq: I,
    r0: I,
    dr: I,
    s: I,
    abs_s: I,
    m: I,
    g0: I,
    gd: I,
    a: I,
    k: I,
    f: I,
}
fn pi() -> I {
    I::new(
        std::f64::consts::PI.next_down(),
        std::f64::consts::PI.next_up(),
    )
}
fn constant(v: f64) -> S {
    S::constant(I::point(v))
}
fn circle(radius: I, outer: Option<I>, z: f64) -> [I; 4] {
    let squared = outer.map_or_else(|| radius.square(), |r| (r - radius) * (r + radius));
    let area = pi() * squared;
    [area, I::point(0.), I::point(0.), area * I::point(z)]
}
fn cylinder(radius: I, z0: f64, z1: f64) -> [I; 4] {
    let area = I::point(2.) * pi() * radius * (I::point(z1) - I::point(z0));
    [
        area,
        I::point(0.),
        I::point(0.),
        area * (I::point(z0) + I::point(z1)) / I::point(2.),
    ]
}
/// Borrowed actual material domain; bare threads have neither treatment.
#[derive(Clone, Copy)]
pub(super) struct Domain<'a> {
    thread: &'a AuthoredThread,
    entry: Option<Entry>,
    tail: Option<Tail>,
}
impl<'a> Domain<'a> {
    pub fn treated(body: &'a Cavity) -> Self {
        Self {
            thread: &body.thread,
            entry: body.entry,
            tail: body.tail,
        }
    }
    pub fn bare(thread: &'a AuthoredThread) -> Self {
        Self {
            thread,
            entry: None,
            tail: None,
        }
    }
}
pub(super) fn evaluate(
    body: Domain<'_>,
    stock: Option<I>,
    face: FaceKey,
    options: FacePropertyOptions,
    work: &mut integrate::Work<'_>,
) -> Result<[I; 4], Error> {
    let length = body.thread.full_length_mm();
    let kind = match face {
        FaceKey::MaterialOuterWall => {
            return Ok(cylinder(stock.ok_or(Error::Invalid)?, 0., length));
        }
        FaceKey::TerminalCylinder => {
            let t = body.tail.ok_or(Error::Invalid)?;
            return Ok(cylinder(I::point(t.radius), t.station, length));
        }
        FaceKey::StartCap => {
            if let Some(e) = body.entry {
                return Ok(circle(I::point(e.radius), stock, 0.));
            } else {
                Kind::Disk {
                    z: 0.,
                    outer: stock,
                }
            }
        }
        FaceKey::EndCap => {
            if let Some(t) = body.tail {
                return Ok(circle(I::point(t.radius), stock, length));
            } else {
                Kind::Disk {
                    z: length,
                    outer: stock,
                }
            }
        }
        FaceKey::Thread { .. } => Kind::Wall,
        FaceKey::EntryCone => {
            body.entry.ok_or(Error::Invalid)?;
            Kind::Cone
        }
        FaceKey::TerminalShoulder => {
            let t = body.tail.ok_or(Error::Invalid)?;
            Kind::Disk {
                z: t.station,
                outer: Some(I::point(t.radius)),
            }
        }
    };
    let spec = body.thread.spec();
    let lead = I::point(spec.pitch_mm) * I::point(f64::from(spec.starts));
    let mut s = lead / (I::point(2.) * pi());
    if spec.left_hand {
        s = -s
    }
    let abs_s = if spec.left_hand { -s } else { s };
    let (a, k) = body.entry.map_or((I::point(0.), I::point(1.)), |e| {
        (I::point(e.radius), I::new(e.slope[0], e.slope[1]))
    });
    let f = I::point(body.tail.map_or(length, |t| t.station));
    let mut parts = Vec::with_capacity(5 * spec.starts as usize);
    for period in 0..spec.starts {
        for (band, b) in body.thread.bands().into_iter().enumerate() {
            work.checkpoint()?;
            if matches!(kind, Kind::Wall)
                && thread_face_key(body.thread, i64::from(period), band)? != face
            {
                continue;
            }
            let dq = I::point(b.phase_mm[1]) - I::point(b.phase_mm[0]);
            let dr = if b.radius_mm[0] == b.radius_mm[1] {
                I::point(0.)
            } else {
                I::point(b.radius_mm[1]) - I::point(b.radius_mm[0])
            };
            let r0 = I::point(b.radius_mm[0]);
            let (g0, gd) = if body.entry.is_some() {
                ((a - r0) / k, -dr / k)
            } else {
                (I::point(0.), I::point(0.))
            };
            parts.push(Part {
                kind,
                q0: I::point(b.phase_mm[0]) + I::point(f64::from(period)) * I::point(spec.pitch_mm),
                dq,
                r0,
                dr,
                s,
                abs_s,
                m: dr / dq,
                g0,
                gd,
                a,
                k,
                f,
            });
        }
    }
    let extent = stock
        .map_or(
            body.thread
                .maximum_radius_mm()
                .max(body.entry.map_or(0., |e| e.radius))
                .max(body.tail.map_or(0., |t| t.radius)),
            |r| r.hi,
        )
        .max(length)
        .max(1.);
    integrate::parts(&parts, extent, options, work)
}
impl Part {
    pub fn field(&self, u: I) -> Result<[S; 4], Error> {
        let variable = S::variable(u);
        let r = variable.scale(self.dr).add(S::constant(self.r0));
        let g = variable.scale(self.gd).add(S::constant(self.g0));
        let f = S::constant(self.f);
        let angle = |z0: I, zd: I| sin_cos((z0 - self.q0) / self.s, (zd - self.dq) / self.s, u);
        let values = match self.kind {
            Kind::Wall => {
                let metric = r
                    .scale(I::point(1.) / self.s)
                    .square()
                    .scale(I::point(1.) + self.m.square())
                    .add(S::constant(self.m.square()))
                    .sqrt()?
                    .scale(self.dq);
                let (sinf, cosf) = angle(self.f, I::point(0.));
                let (sing, cosg) = angle(self.g0, self.gd);
                let area = f.sub(g).mul(metric);
                let xy = metric.mul(r).scale(self.s);
                [
                    area,
                    xy.mul(sinf.sub(sing)),
                    xy.mul(cosg.sub(cosf)),
                    f.square().sub(g.square()).mul(metric).scale(I::point(0.5)),
                ]
            }
            Kind::Cone => {
                let factor = (I::point(1.) + self.k.square()).sqrt() / self.abs_s * self.dq;
                let a = S::constant(self.a);
                let area = g.mul(a.add(r)).scale(factor * I::point(0.5));
                let z = g
                    .square()
                    .mul(a.add(r.scale(I::point(2.))))
                    .scale(factor / I::point(6.));
                let (sing, cosg) = angle(self.g0, self.gd);
                let (sin0, cos0) = angle(I::point(0.), I::point(0.));
                let ks = self.k * self.s;
                let ks2 = I::point(2.) * ks.square();
                let b = r.square().sub(S::constant(ks2));
                let b0 = a.square().sub(S::constant(ks2));
                let x = b.mul(sing).sub(b0.mul(sin0)).scale(self.s).sub(
                    r.mul(cosg)
                        .sub(a.mul(cos0))
                        .scale(I::point(2.) * self.k * self.s.square()),
                );
                let y = b0.mul(cos0).sub(b.mul(cosg)).scale(self.s).sub(
                    r.mul(sing)
                        .sub(a.mul(sin0))
                        .scale(I::point(2.) * self.k * self.s.square()),
                );
                [area, x.scale(factor), y.scale(factor), z]
            }
            Kind::Disk { z, outer } => {
                let (sine, cosine) = angle(I::point(z), I::point(0.));
                let squared = outer.map_or_else(
                    || r.square(),
                    |v| S::constant(v).sub(r).mul(S::constant(v).add(r)),
                );
                let area = squared.scale(self.dq / (I::point(2.) * self.abs_s));
                // The complete outer disk's x/y moments vanish over the full
                // one-lead source domain. Inner profile contributions remain.
                let radial = r.square().mul(r).scale(
                    self.dq / (I::point(3.) * self.abs_s)
                        * I::point(if outer.is_some() { -1. } else { 1. }),
                );
                [
                    area,
                    radial.mul(cosine),
                    radial.mul(sine),
                    area.scale(I::point(z)),
                ]
            }
        };
        let mut checked = [constant(0.); 4];
        for (out, v) in checked.iter_mut().zip(values) {
            *out = v.checked()?
        }
        Ok(checked)
    }
}
