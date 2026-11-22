//! Orientation-preserving affine images of a retained boundary.
//!
//! Points, lines and spline control nets map directly; affine maps commute
//! with rational B-spline evaluation. A circle or ellipse maps to the ellipse
//! spanned by its two mapped radius vectors, re-expressed on its principal
//! axes with a parameter phase shift, so every edge range stays exact. A
//! cylinder maps to the linear extrusion of its mapped base ellipse. Cones,
//! spheres, tori and revolutions keep their families only under a
//! similarity; any other image of them is refused rather than approximated.
use super::{
    Boundary, Curve, Edge, Error, Face, Frame, Point, Result, Surface, add, cross, dot, ensure,
    norm, scale,
};

/// `point -> matrix * point + translation` with a positive determinant.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Affine {
    pub matrix: [[f64; 3]; 3],
    pub translation: Point,
}

/// Relative agreement that admits a uniform scale or orthogonality.
const SIMILAR: f64 = 1e-12;

impl Affine {
    fn linear(&self, v: Point) -> Point {
        self.matrix.map(|row| dot(row, v))
    }

    fn point(&self, p: Point) -> Point {
        add(self.linear(p), self.translation)
    }

    fn determinant(&self) -> f64 {
        let [a, b, c] = self.matrix;
        dot(a, cross(b, c))
    }

    /// The uniform factor when the linear part is a rotation times it.
    fn similarity(&self) -> Option<f64> {
        let columns: [Point; 3] =
            std::array::from_fn(|j| self.linear(std::array::from_fn(|i| f64::from(i == j))));
        let s = norm(columns[0]);
        let uniform = columns.iter().all(|c| (norm(*c) - s).abs() <= SIMILAR * s);
        let orthogonal = (0..3).all(|i| {
            (0..3)
                .filter(|&j| j != i)
                .all(|j| dot(columns[i], columns[j]).abs() <= SIMILAR * s * s)
        });
        (uniform && orthogonal).then_some(s)
    }

    fn similar_frame(&self, frame: &Frame, s: f64) -> Frame {
        Frame {
            origin: self.point(frame.origin),
            x: scale(self.linear(frame.x), 1. / s),
            y: scale(self.linear(frame.y), 1. / s),
            z: scale(self.linear(frame.z), 1. / s),
        }
    }

    /// A right-handed orthonormal frame on the mapped plane of `frame`.
    fn plane_frame(&self, frame: &Frame) -> Result<Frame> {
        let x = unit(self.linear(frame.x))?;
        let z = unit(cross(self.linear(frame.x), self.linear(frame.y)))?;
        Ok(Frame {
            origin: self.point(frame.origin),
            x,
            y: cross(z, x),
            z,
        })
    }

    /// The image of `origin + a cos t x + b sin t y` as a conic on its
    /// principal axes, and the phase `phi` with `image(t) = conic(t - phi)`.
    fn conic(&self, frame: &Frame, [a, b]: [f64; 2]) -> Result<(Curve, f64)> {
        let p = scale(self.linear(frame.x), a);
        let q = scale(self.linear(frame.y), b);
        let [pp, pq, qq] = [dot(p, p), dot(p, q), dot(q, q)];
        let magnitude = pp.max(qq);
        if (pp - qq).abs() <= SIMILAR * magnitude && pq.abs() <= SIMILAR * magnitude {
            // Already a circle: the mapped radii are orthogonal and equal.
            let radius = 0.5 * (pp.sqrt() + qq.sqrt());
            let x = unit(p)?;
            let z = unit(cross(p, q))?;
            let frame = Frame {
                origin: self.point(frame.origin),
                x,
                y: cross(z, x),
                z,
            };
            return Ok((Curve::Circle { frame, radius }, 0.));
        }
        let phi = 0.5 * (2. * pq).atan2(pp - qq);
        let (sin, cos) = phi.sin_cos();
        let major = add(scale(p, cos), scale(q, sin));
        let minor = add(scale(p, -sin), scale(q, cos));
        let x = unit(major)?;
        let z = unit(cross(major, minor))?;
        let frame = Frame {
            origin: self.point(frame.origin),
            x,
            y: cross(z, x),
            z,
        };
        let semiaxes = [norm(major), norm(minor)];
        ensure(
            semiaxes.iter().all(|s| s.is_finite() && *s > 0.),
            Error::Invalid,
        )?;
        Ok((Curve::Ellipse { frame, semiaxes }, phi))
    }

    /// The mapped curve and the shift subtracted from its parameters.
    fn curve(&self, curve: &Curve, similar: Option<f64>) -> Result<(Curve, f64)> {
        Ok(match (curve, similar) {
            (Curve::Line { origin, direction }, _) => (
                Curve::Line {
                    origin: self.point(*origin),
                    direction: self.linear(*direction),
                },
                0.,
            ),
            (Curve::Nurbs { curve, closed }, _) => (
                Curve::Nurbs {
                    curve: curve.mapped(&|p| self.point(p)),
                    closed: *closed,
                },
                0.,
            ),
            (Curve::Circle { frame, radius }, Some(s)) => (
                Curve::Circle {
                    frame: self.similar_frame(frame, s),
                    radius: radius * s,
                },
                0.,
            ),
            (Curve::Ellipse { frame, semiaxes }, Some(s)) => (
                Curve::Ellipse {
                    frame: self.similar_frame(frame, s),
                    semiaxes: semiaxes.map(|a| a * s),
                },
                0.,
            ),
            (Curve::Hyperbola { frame, semiaxes }, Some(s)) => (
                Curve::Hyperbola {
                    frame: self.similar_frame(frame, s),
                    semiaxes: semiaxes.map(|a| a * s),
                },
                0.,
            ),
            (Curve::Parabola { frame, focal }, Some(s)) => (
                Curve::Parabola {
                    frame: self.similar_frame(frame, s),
                    focal: focal * s,
                },
                0.,
            ),
            (Curve::Circle { frame, radius }, None) => self.conic(frame, [*radius; 2])?,
            (Curve::Ellipse { frame, semiaxes }, None) => self.conic(frame, *semiaxes)?,
            (Curve::Hyperbola { .. } | Curve::Parabola { .. }, None) => {
                return Err(Error::Unsupported);
            }
        })
    }

    fn surface(&self, surface: &Surface, similar: Option<f64>) -> Result<Surface> {
        Ok(match (surface, similar) {
            (Surface::Plane(frame), _) => Surface::Plane(self.plane_frame(frame)?),
            (Surface::Nurbs { surface, closed }, _) => Surface::Nurbs {
                surface: surface.mapped(&|p| self.point(p)),
                closed: *closed,
            },
            (Surface::Extrusion { curve, direction }, _) => Surface::Extrusion {
                curve: self.curve(curve, similar)?.0,
                direction: self.linear(*direction),
            },
            (Surface::Cylinder { frame, radius }, Some(s)) => Surface::Cylinder {
                frame: self.similar_frame(frame, s),
                radius: radius * s,
            },
            (Surface::Cylinder { frame, radius }, None) => {
                // `origin + r radial(u) + v z` maps to the extrusion of the
                // mapped base circle along the mapped axis. The conic phase
                // only renames u, which the face domain recovers from its
                // trim edges; the parameter orientation is unchanged.
                Surface::Extrusion {
                    curve: self.conic(frame, [*radius; 2])?.0,
                    direction: self.linear(frame.z),
                }
            }
            (
                Surface::Cone {
                    frame,
                    radius,
                    half_angle,
                },
                Some(s),
            ) => Surface::Cone {
                frame: self.similar_frame(frame, s),
                radius: radius * s,
                half_angle: *half_angle,
            },
            (Surface::Sphere { frame, radius }, Some(s)) => Surface::Sphere {
                frame: self.similar_frame(frame, s),
                radius: radius * s,
            },
            (
                Surface::Torus {
                    frame,
                    major,
                    minor,
                },
                Some(s),
            ) => Surface::Torus {
                frame: self.similar_frame(frame, s),
                major: major * s,
                minor: minor * s,
            },
            (
                Surface::Revolution {
                    curve,
                    origin,
                    axis,
                },
                Some(s),
            ) => Surface::Revolution {
                curve: self.curve(curve, Some(s))?.0,
                origin: self.point(*origin),
                axis: scale(self.linear(*axis), 1. / s),
            },
            (
                Surface::Cone { .. }
                | Surface::Sphere { .. }
                | Surface::Torus { .. }
                | Surface::Revolution { .. },
                None,
            ) => return Err(Error::Unsupported),
        })
    }
}

fn unit(v: Point) -> Result<Point> {
    let length = norm(v);
    ensure(length.is_finite() && length > 0., Error::Invalid)?;
    Ok(scale(v, 1. / length))
}

impl Boundary {
    /// The boundary's image under `map`, which must preserve orientation.
    /// Faces keep their loops and senses; validation of the image is left to
    /// [`super::Solid::new`], which recomputes every domain and property.
    pub fn affine_image(&self, map: Affine) -> Result<Self> {
        let determinant = map.determinant();
        ensure(
            map.matrix
                .iter()
                .flatten()
                .chain(&map.translation)
                .all(|v| v.is_finite())
                && determinant.is_finite()
                && determinant > 0.,
            Error::Invalid,
        )?;
        let similar = map.similarity();
        let edges = self
            .edges
            .iter()
            .map(|edge| {
                let (curve, phase) = map.curve(&edge.curve, similar)?;
                Ok(Edge {
                    curve,
                    vertices: edge.vertices,
                    range: edge.range.map(|t| t - phase),
                })
            })
            .collect::<Result<_>>()?;
        let faces = self
            .faces
            .iter()
            .map(|face| {
                Ok(Face {
                    surface: map.surface(&face.surface, similar)?,
                    same_sense: face.same_sense,
                    loops: face.loops.clone(),
                })
            })
            .collect::<Result<_>>()?;
        Ok(Self {
            vertices: self.vertices.iter().map(|&p| map.point(p)).collect(),
            edges,
            faces,
            shells: self.shells.clone(),
        })
    }
}
