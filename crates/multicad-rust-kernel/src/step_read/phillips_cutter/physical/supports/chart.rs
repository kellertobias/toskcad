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
        self.project_vector(std::array::from_fn(|i| point[i] - self.frame.origin[i]))
    }
    pub(in super::super) fn project_vector(&self, delta: [I; 3]) -> [I; 3] {
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
    pub(in super::super) fn homogeneous_degree(&self) -> usize {
        match self.kind {
            Kind::Plane => 1,
            Kind::Cylinder(_) | Kind::Cone(..) => 2,
            Kind::Torus(..) => 4,
        }
    }
    /// Exact rational quadratic parameterization, not a polynomial fit. The
    /// result is W^degree * F(numerator/W), through Bernstein degree8.
    pub(in super::super) fn implicit_homogeneous(&self, n: [[I; 3]; 3], w: [I; 3]) -> Vec<I> {
        let [x, y, z] = [0, 1, 2].map(|i| n.map(|p| p[i]));
        let radial = add(&square(&x), &square(&y));
        let ww = square(&w);
        match self.kind {
            Kind::Plane => z.to_vec(),
            Kind::Cylinder(r) => sub(
                &radial,
                &ww.iter().map(|&v| v * r.square()).collect::<Vec<_>>(),
            ),
            Kind::Cone(r, k) => sub(
                &radial,
                &square(&std::array::from_fn::<_, 3, _>(|i| r * w[i] + k * z[i])),
            ),
            Kind::Torus(a, b) => {
                let norm = add(
                    &add(&radial, &square(&z)),
                    &ww.iter()
                        .map(|&v| v * (a.square() - b.square()))
                        .collect::<Vec<_>>(),
                );
                sub(
                    &square(&norm),
                    &product(&radial, &ww)
                        .iter()
                        .map(|&v| I::point(4.) * a.square() * v)
                        .collect::<Vec<_>>(),
                )
            }
        }
    }
    /// Uniform whole-support metric; deliberately conservative. Finite chart
    /// filling and branch/winding validity are independent responsibilities.
    pub(in super::super) fn uniform_distance(&self, residual: I, local_box: [I; 3]) -> Result<f64> {
        let denominator = match self.kind {
            Kind::Plane => return finite_bound(residual.magnitude()),
            Kind::Cylinder(r) => r,
            Kind::Cone(r, k) => r + k * local_box[2],
            Kind::Torus(a, b) => b * (a.square() - b.square()),
        };
        if !denominator.valid() || denominator.lo <= 0. {
            return Err(ReadError::Unsupported);
        }
        finite_bound((I::point(residual.magnitude()) / denominator).hi)
    }
    /// Bound one finite affine segment against a positive cone ruling, also
    /// at its apex. The ruling is only a comparison for an error bound; this
    /// method never replaces the transmitted line or its endpoint data.
    pub(in super::super) fn cone_line_distance(&self, p: [[I; 3]; 2]) -> Result<f64> {
        let Kind::Cone(r, k) = self.kind else {
            return Err(ReadError::Unsupported);
        };
        if k.lo <= 0. {
            return Err(ReadError::Unsupported);
        }
        let radial = p.map(|q| (q[0].square() + q[1].square()).sqrt());
        if radial.into_iter().any(|i| !i.valid()) {
            return Err(ReadError::Limit);
        }
        let index = usize::from(radial[1].lo > radial[0].lo);
        if radial[index].lo <= 0. {
            return Err(ReadError::Unsupported);
        }
        let unit = [p[index][0] / radial[index], p[index][1] / radial[index]];
        let mut upper = 0_f64;
        for point in p {
            let q = r + k * point[2];
            if !q.valid() {
                return Err(ReadError::Limit);
            }
            let clipped = I::new(q.lo.max(0.), q.hi.max(0.));
            let target = [clipped * unit[0], clipped * unit[1], (clipped - r) / k];
            let deviation = support::norm(std::array::from_fn(|i| point[i] - target[i]));
            if !deviation.valid() {
                return Err(ReadError::Limit);
            }
            upper = upper.max(deviation.hi);
        }
        finite_bound(upper)
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
