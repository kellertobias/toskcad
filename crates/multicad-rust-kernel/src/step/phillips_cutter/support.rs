//! Directed whole-support source-to-literal envelopes. Requested precision is
//! never an input to the derived numeric error.
use super::*;
pub(crate) type BPoint = [I; 3];
pub(crate) fn b([a, z]: [f64; 2]) -> I {
    I::new(a, z)
}
pub(crate) fn adjacent(x: f64) -> I {
    I::new(x.next_down(), x.next_up())
}
pub(crate) fn norm(p: BPoint) -> I {
    p.into_iter()
        .fold(I::point(0.), |s, x| s + x.square())
        .sqrt()
}
pub(crate) fn cross(a: BPoint, b: BPoint) -> BPoint {
    std::array::from_fn(|i| a[(i + 1) % 3] * b[(i + 2) % 3] - a[(i + 2) % 3] * b[(i + 1) % 3])
}
fn absolute(x: I) -> I {
    I::new(0., x.magnitude())
}
fn difference(a: BPoint, b: BPoint) -> BPoint {
    std::array::from_fn(|i| a[i] - b[i])
}
pub(super) fn unit(p: BPoint) -> Result<BPoint, Error> {
    let scale = p.into_iter().map(I::magnitude).fold(0_f64, f64::max);
    if !scale.is_finite() || scale == 0. {
        return Err(Error::InvalidGeometry);
    }
    let p = p.map(|v| I::new((v.lo / scale).next_down(), (v.hi / scale).next_up()));
    let n = norm(p);
    if !n.valid() || n.lo <= 0. {
        return Err(Error::Unsupported);
    }
    // Scaled components are bounded by one plus rounding. A positive square
    // root is at least sqrt(min_subnormal), so this quotient cannot overflow.
    Ok(p.map(|v| v / n))
}
#[derive(Clone, Copy)]
pub(crate) struct Frame {
    pub origin: BPoint,
    pub axes: [BPoint; 3],
}
impl Frame {
    pub fn new(origin: BPoint, normal: BPoint, reference: BPoint) -> Result<Self, Error> {
        let n = unit(normal)?;
        let v = unit(cross(n, unit(reference)?))?;
        Ok(Self {
            origin,
            axes: [cross(v, n), v, n],
        })
    }
    fn source(origin: BPoint, normal: BPoint) -> Result<Self, Error> {
        let k = (0..3)
            .min_by(|&a, &z| normal[a].magnitude().total_cmp(&normal[z].magnitude()))
            .ok_or(Error::InvalidGeometry)?;
        let seed = std::array::from_fn(|i| I::point(if i == k { 1. } else { 0. }));
        Self::new(origin, normal, cross(seed, normal))
    }
    pub fn literal(raw: &LiteralAxis) -> Result<Self, Error> {
        Self::new(
            raw.origin.map(adjacent),
            raw.normal.map(adjacent),
            raw.reference.map(adjacent),
        )
    }
}
#[derive(Clone, Copy, Debug)]
pub(crate) struct LiteralAxis {
    pub origin: Point,
    pub normal: Point,
    pub reference: Point,
}
impl LiteralAxis {
    pub(super) fn emit(&self, file: &mut StepFile) -> Entity {
        let p = file.point(self.origin);
        file.axis(p, self.normal, self.reference)
    }
}
#[derive(Clone, Copy)]
pub(crate) struct Pose {
    pub origin: Point,
    pub axes: [Point; 3],
}
impl Pose {
    pub fn new(origin: Point, u: Point, v: Point) -> Result<Self, Error> {
        let p = Placement::new(origin, u, v)?;
        Ok(Self {
            origin: p.origin,
            axes: [p.u, p.v, p.normal],
        })
    }
    pub fn vector(&self, [x, y, z]: Point) -> Point {
        std::array::from_fn(|i| self.axes[0][i] * x + self.axes[1][i] * y + self.axes[2][i] * z)
    }
    pub fn point(&self, p: Point) -> Point {
        std::array::from_fn(|i| {
            self.origin[i]
                + self.axes[0][i] * p[0]
                + self.axes[1][i] * p[1]
                + self.axes[2][i] * p[2]
        })
    }
}
pub(crate) fn map_vector(p: &Pose, v: BPoint) -> BPoint {
    std::array::from_fn(|i| {
        I::point(p.axes[0][i]) * v[0]
            + I::point(p.axes[1][i]) * v[1]
            + I::point(p.axes[2][i]) * v[2]
    })
}
pub(crate) fn map_point(p: &Pose, v: BPoint) -> BPoint {
    let v = map_vector(p, v);
    std::array::from_fn(|i| I::point(p.origin[i]) + v[i])
}
pub(crate) fn literal_axis(p: &Pose, frame: Frame) -> LiteralAxis {
    LiteralAxis {
        origin: p.point(frame.origin.map(I::midpoint)),
        normal: p.vector(frame.axes[2].map(I::midpoint)),
        reference: p.vector(frame.axes[0].map(I::midpoint)),
    }
}
#[derive(Clone, Copy, Debug)]
pub(crate) enum Kind {
    Plane,
    Cylinder(f64),
    Torus(f64, f64),
    Cone(f64, f64),
}
#[derive(Clone, Copy, Debug)]
pub(crate) struct LiteralSupport {
    pub axis: LiteralAxis,
    pub kind: Kind,
    pub forward: bool,
}
impl LiteralSupport {
    pub(super) fn emit(&self, file: &mut StepFile) -> Entity {
        let axis = self.axis.emit(file);
        match self.kind {
            Kind::Plane => file.entity(format_args!("PLANE('',{axis})")),
            Kind::Cylinder(r) => {
                file.entity(format_args!("CYLINDRICAL_SURFACE('',{axis},{r:.16E})"))
            }
            Kind::Torus(a, z) => file.entity(format_args!(
                "TOROIDAL_SURFACE('',{axis},{a:.16E},{z:.16E})"
            )),
            Kind::Cone(r, a) => {
                file.entity(format_args!("CONICAL_SURFACE('',{axis},{r:.16E},{a:.16E})"))
            }
        }
    }
}
#[derive(Clone, Copy)]
enum SourceKind {
    Plane,
    Cylinder(I),
    Torus(I, I),
    Cone(I, I),
}
#[derive(Clone, Copy)]
pub(crate) struct Support {
    frame: Frame,
    kind: SourceKind,
    forward: bool,
}
impl Support {
    pub fn new(surface: Surface) -> Result<Self, Error> {
        let zero = I::point(0.);
        let z = [zero, zero, I::point(1.)];
        let x = [I::point(1.), zero, zero];
        let (frame, kind, forward) = match surface {
            Surface::Plane {
                normal_bounds,
                offset_mm_bounds,
            } => {
                let n = normal_bounds.map(b);
                (
                    Frame::source(n.map(|v| v * b(offset_mm_bounds)), n)?,
                    SourceKind::Plane,
                    true,
                )
            }
            Surface::Cylinder {
                origin_mm_bounds,
                axis_bounds,
                radius_mm,
                outward_sign,
            } => (
                Frame::source(origin_mm_bounds.map(b), axis_bounds.map(b))?,
                SourceKind::Cylinder(I::point(radius_mm)),
                outward_sign > 0,
            ),
            Surface::Torus {
                major_radius_mm_bounds,
                minor_radius_mm,
                center_depth_mm,
                outward_sign,
            } => (
                Frame::new([zero, zero, I::point(center_depth_mm)], z, x)?,
                SourceKind::Torus(b(major_radius_mm_bounds), I::point(minor_radius_mm)),
                outward_sign > 0,
            ),
            Surface::Cone {
                radius_at_zero_mm_bounds,
                radius_decrease_per_depth_bounds,
            } => (
                Frame::new([zero; 3], z.map(|v| -v), x)?,
                SourceKind::Cone(
                    b(radius_at_zero_mm_bounds),
                    b(radius_decrease_per_depth_bounds),
                ),
                true,
            ),
        };
        Ok(Self {
            frame,
            kind,
            forward,
        })
    }
    pub fn literal(&self, p: &Pose) -> LiteralSupport {
        let kind = match self.kind {
            SourceKind::Plane => Kind::Plane,
            SourceKind::Cylinder(r) => Kind::Cylinder(r.midpoint()),
            SourceKind::Torus(a, z) => Kind::Torus(a.midpoint(), z.midpoint()),
            SourceKind::Cone(r, k) => Kind::Cone(r.midpoint(), k.atan().midpoint()),
        };
        LiteralSupport {
            axis: literal_axis(p, self.frame),
            kind,
            forward: self.forward,
        }
    }
    /// A parameterwise C0 enclosure on the complete finite support domain. It
    /// includes literal decimal conversion, AXIS2 normalization and frame strain.
    pub fn distance_bound(
        &self,
        p: &Pose,
        raw: &LiteralSupport,
        source_box: [[f64; 3]; 2],
    ) -> Result<f64, Error> {
        if raw.forward != self.forward {
            return Err(Error::Unsupported);
        }
        let literal = Frame::literal(&raw.axis)?;
        let o = map_point(p, self.frame.origin);
        let a = self.frame.axes.map(|v| map_vector(p, v));
        let center = norm(difference(o, literal.origin));
        let length = norm(std::array::from_fn(|i| {
            I::new(source_box[0][i], source_box[1][i]) - self.frame.origin[i]
        }));
        let radius = length;
        let vector_error = |source: BPoint, actual: BPoint| norm(difference(source, actual));
        let scaled = |v: BPoint, s: I| v.map(|x| x * s);
        let delta = match (self.kind, raw.kind) {
            (SourceKind::Plane, Kind::Plane) => {
                let n = literal.axes[2];
                let dot = |v: BPoint| {
                    v.into_iter()
                        .zip(n)
                        .fold(I::point(0.), |s, (a, b)| s + a * b)
                };
                if dot(a[2]).lo <= 0. {
                    return Err(Error::Unsupported);
                }
                absolute(dot(difference(o, literal.origin)))
                    + length * (absolute(dot(a[0])) + absolute(dot(a[1])))
            }
            (SourceKind::Cylinder(r), Kind::Cylinder(actual)) if actual > 0. => {
                let actual = adjacent(actual);
                center
                    + vector_error(scaled(a[0], r), scaled(literal.axes[0], actual))
                    + vector_error(scaled(a[1], r), scaled(literal.axes[1], actual))
                    + length * vector_error(a[2], literal.axes[2])
            }
            (SourceKind::Torus(major, minor), Kind::Torus(actual_major, actual_minor))
                if actual_major > actual_minor && actual_minor > 0. =>
            {
                let major_raw = adjacent(actual_major);
                let minor_raw = adjacent(actual_minor);
                center
                    + vector_error(scaled(a[0], major), scaled(literal.axes[0], major_raw))
                    + vector_error(scaled(a[1], major), scaled(literal.axes[1], major_raw))
                    + vector_error(scaled(a[0], minor), scaled(literal.axes[0], minor_raw))
                    + vector_error(scaled(a[1], minor), scaled(literal.axes[1], minor_raw))
                    + vector_error(scaled(a[2], minor), scaled(literal.axes[2], minor_raw))
            }
            (SourceKind::Cone(reference, k), Kind::Cone(actual_reference, angle))
                if actual_reference >= 0. && angle > 0. && angle < std::f64::consts::FRAC_PI_2 =>
            {
                let angle = adjacent(angle);
                let cot = angle.cos() / angle.sin();
                let actual_radius = adjacent(actual_reference);
                let source_apex = std::array::from_fn(|i| o[i] - a[2][i] * reference / k);
                let raw_apex = std::array::from_fn(|i| {
                    literal.origin[i] - literal.axes[2][i] * actual_radius * cot
                });
                norm(difference(source_apex, raw_apex))
                    + radius
                        * (vector_error(a[0], literal.axes[0])
                            + vector_error(a[1], literal.axes[1])
                            + vector_error(
                                scaled(a[2], I::point(1.) / k),
                                scaled(literal.axes[2], cot),
                            ))
            }
            _ => return Err(Error::Unsupported),
        };
        if !delta.valid() {
            return Err(Error::Unsupported);
        }
        Ok(delta.hi)
    }
}
