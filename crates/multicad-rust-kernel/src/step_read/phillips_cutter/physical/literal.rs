//! Retain literal numeric supports. No gauge or source expression is fitted.
use super::*;
use std::f64::consts::PI;

#[derive(Clone, Copy, Debug)]
pub struct Axis {
    pub origin: Point,
    pub normal: Point,
    pub reference: Point,
}
impl From<LiteralAxis> for Axis {
    fn from(a: LiteralAxis) -> Self {
        Self {
            origin: a.origin,
            normal: a.normal,
            reference: a.reference,
        }
    }
}
#[derive(Clone, Debug)]
pub enum Surface {
    Plane {
        axis: Axis,
    },
    Cylinder {
        axis: Axis,
        radius_mm: f64,
    },
    Cone {
        axis: Axis,
        radius_mm: f64,
        semi_angle_rad: f64,
    },
    Torus {
        axis: Axis,
        major_radius_mm: f64,
        minor_radius_mm: f64,
    },
}
pub(super) fn surface(shell: &Shell<'_, '_>, id: u32, forward: bool) -> Result<Surface> {
    let s = geometry::support(shell, id, forward)?;
    support::Frame::literal(&s.axis).map_err(source_error)?;
    let axis = s.axis.into();
    Ok(match s.kind {
        Kind::Plane => Surface::Plane { axis },
        Kind::Cylinder(r) => {
            require(r.is_finite() && r > 0.)?;
            Surface::Cylinder { axis, radius_mm: r }
        }
        Kind::Cone(r, a) => {
            require(r.is_finite() && r >= 0. && a.is_finite() && a > 0. && a < PI / 2.)?;
            Surface::Cone {
                axis,
                radius_mm: r,
                semi_angle_rad: a,
            }
        }
        Kind::Torus(a, b) => {
            require(a.is_finite() && b.is_finite() && b > 0. && a > b)?;
            Surface::Torus {
                axis,
                major_radius_mm: a,
                minor_radius_mm: b,
            }
        }
    })
}
#[derive(Clone, Debug)]
pub enum Curve {
    Line {
        origin: Point,
        direction: Point,
        magnitude: f64,
        vector_entity: u32,
        direction_entity: u32,
    },
    Circle {
        axis: Axis,
        radius_mm: f64,
    },
    Ellipse {
        axis: Axis,
        semiaxes_mm: [f64; 2],
    },
    /// The file supplies surface associations, not a fitted implicit-curve
    /// recipe. PCURVE associations are not accepted by this narrow parser.
    CubicIntersection {
        spline_entity: u32,
        associated_surface_entities: [u32; 2],
        pole_entities: Vec<u32>,
        spans: Vec<CubicTrimSpan>,
    },
}
pub(super) fn curve(
    shell: &Shell<'_, '_>,
    edge: graph::Edge,
    owners: [u32; 2],
    endpoints: [u32; 2],
    work: &mut Work<'_>,
) -> Result<Curve> {
    if shell.simple(edge.geometry)?.name == "INTERSECTION_CURVE" {
        let ends = if edge.forward {
            endpoints
        } else {
            [endpoints[1], endpoints[0]]
        };
        let spans = geometry::spline(shell, edge.geometry, owners, ends, work)?;
        let args = shell.args(edge.geometry, "INTERSECTION_CURVE")?;
        let spline_entity = reference(&args[1])?;
        let associated = refs(&args[2])?;
        let poles = refs(&shell.args(spline_entity, "B_SPLINE_CURVE_WITH_KNOTS")?[2])?;
        return Ok(Curve::CubicIntersection {
            spline_entity,
            associated_surface_entities: [associated[0], associated[1]],
            pole_entities: poles,
            spans,
        });
    }
    Ok(match geometry::analytic(shell, edge)? {
        geometry::Analytic::Line {
            origin,
            direction,
            vector_entity,
            direction_entity,
            magnitude,
            ..
        } => {
            require(
                origin.into_iter().chain(direction).all(f64::is_finite)
                    && support::norm(direction.map(I::point)).lo > 0.,
            )?;
            Curve::Line {
                origin,
                direction,
                magnitude,
                vector_entity,
                direction_entity,
            }
        }
        geometry::Analytic::Circle { axis, radius, .. } => {
            support::Frame::literal(&axis).map_err(source_error)?;
            require(radius.is_finite() && radius > 0.)?;
            Curve::Circle {
                axis: axis.into(),
                radius_mm: radius,
            }
        }
        geometry::Analytic::Ellipse { axis, semiaxes, .. } => {
            support::Frame::literal(&axis).map_err(source_error)?;
            require(
                semiaxes.into_iter().all(|r| r.is_finite() && r > 0.) && semiaxes[0] >= semiaxes[1],
            )?;
            Curve::Ellipse {
                axis: axis.into(),
                semiaxes_mm: semiaxes,
            }
        }
    })
}
