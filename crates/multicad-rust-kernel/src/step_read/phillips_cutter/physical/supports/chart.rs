//! True AXIS2 semantics and directed physical implicit-support evaluation.
use super::polynomial::{product, sub};
use super::*;
#[derive(Clone, Copy)]
pub(super) enum Kind {
    Plane,
    Cylinder(I),
    Cone(I, I),
    Torus(I, I),
}
pub(in super::super) struct Chart {
    pub(super) frame: support::Frame,
    pub(super) kind: Kind,
}
fn interval(x: f64) -> I {
    support::adjacent(x)
}
impl Chart {
    pub(in super::super) fn new(surface: &Surface) -> Result<Self> {
        let (axis, kind) = match surface {
            Surface::Plane { axis } => (*axis, Kind::Plane),
            Surface::Cylinder { axis, radius_mm } => (*axis, Kind::Cylinder(interval(*radius_mm))),
            Surface::Cone {
                axis,
                radius_mm,
                semi_angle_rad,
            } => {
                let a = interval(*semi_angle_rad);
                let cosine = a.cos();
                require(cosine.lo > 0.)?;
                (*axis, Kind::Cone(interval(*radius_mm), a.sin() / cosine))
            }
            Surface::Torus {
                axis,
                major_radius_mm,
                minor_radius_mm,
            } => {
                let (a, b) = (interval(*major_radius_mm), interval(*minor_radius_mm));
                if (a.square() - b.square()).lo <= 0. {
                    return Err(ReadError::Unsupported);
                }
                (*axis, Kind::Torus(a, b))
            }
        };
        let frame = support::Frame::literal(&LiteralAxis {
            origin: axis.origin,
            normal: axis.normal,
            reference: axis.reference,
        })
        .map_err(source_error)?;
        Ok(Self { frame, kind })
    }

    pub(in super::super) fn local(&self, point: Point) -> [I; 3] {
        self.project(point.map(interval))
    }
    pub(in super::super) fn project(&self, point: [I; 3]) -> [I; 3] {
        let delta = std::array::from_fn::<_, 3, _>(|i| point[i] - self.frame.origin[i]);
        self.frame.axes.map(|axis| {
            axis.into_iter()
                .zip(delta)
                .fold(I::point(0.), |s, (a, b)| s + a * b)
        })
    }
    pub(in super::super) fn world(&self, point: [I; 3]) -> [I; 3] {
        std::array::from_fn(|i| {
            self.frame.origin[i]
                + (0..3).fold(I::point(0.), |s, j| s + self.frame.axes[j][i] * point[j])
        })
    }
    pub(in super::super) fn implicit(&self, [x, y, z]: [I; 3]) -> I {
        let s = x.square() + y.square();
        match self.kind {
            Kind::Plane => z,
            Kind::Cylinder(r) => s - r.square(),
            Kind::Cone(r, k) => s - (r + k * z).square(),
            Kind::Torus(a, b) => {
                (s + z.square() + a.square() - b.square()).square() - I::point(4.) * a.square() * s
            }
        }
    }
    /// Local cubic input; all mixed terms retained through the full degree12
    /// Bernstein polynomial. Differentiation uses n*(c[i+1]-c[i]) on [0,1].
    pub(in super::super) fn implicit_curve(&self, p: [[I; 3]; 4]) -> Vec<I> {
        let [x, y, z] = [0, 1, 2].map(|i| p.map(|p| p[i]));
        let radial = add(&square(&x), &square(&y));
        match self.kind {
            Kind::Plane => z.to_vec(),
            Kind::Cylinder(r) => radial.iter().map(|&v| v - r.square()).collect(),
            Kind::Cone(r, k) => sub(&radial, &square(&z.map(|v| r + k * v))),
            Kind::Torus(a, b) => {
                let norm = add(&radial, &square(&z))
                    .into_iter()
                    .map(|v| v + a.square() - b.square())
                    .collect::<Vec<_>>();
                let elevated = product(&radial, &[I::point(1.); 7]);
                sub(
                    &square(&norm),
                    &elevated
                        .iter()
                        .map(|&v| I::point(4.) * a.square() * v)
                        .collect::<Vec<_>>(),
                )
            }
        }
    }
    pub(in super::super) fn gradient(&self, [x, y, z]: [I; 3]) -> [I; 3] {
        let two = I::point(2.);
        let four = I::point(4.);
        match self.kind {
            Kind::Plane => [I::point(0.), I::point(0.), I::point(1.)],
            Kind::Cylinder(_) => [two * x, two * y, I::point(0.)],
            Kind::Cone(r, k) => [two * x, two * y, -two * k * (r + k * z)],
            Kind::Torus(a, b) => {
                let t = x.square() + y.square() + z.square() + a.square() - b.square();
                [
                    four * x * (t - two * a.square()),
                    four * y * (t - two * a.square()),
                    four * z * t,
                ]
            }
        }
    }
    pub(in super::super) fn world_gradient(&self, point: [I; 3]) -> [I; 3] {
        let gradient = self.gradient(point);
        std::array::from_fn(|i| {
            (0..3).fold(I::point(0.), |s, j| s + self.frame.axes[j][i] * gradient[j])
        })
    }
}
