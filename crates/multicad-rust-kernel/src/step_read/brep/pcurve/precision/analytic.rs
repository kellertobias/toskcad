//! Raw REAL analytic operators versus the retained binary64 reference.
//! Periodic reserves compare complete coefficient fields after an increasing
//! natural-domain reparameterization, not individual evaluation samples.
use super::super::{analytic as retained, uv};
use super::{
    Budget, I, RawFrame, ReadError, Result, Source, V, cached, difference, frame, list, norm,
    point, reference, scalar, translated_from, translation_origin, unwrapped, vector as raw_vector,
};
use crate::brep::{Curve, Edge, Surface};
use num_rational::BigRational;
use retained::Image;

fn add(a: V, b: V) -> V {
    std::array::from_fn(|i| a[i] + b[i])
}
fn scale(a: V, b: I) -> V {
    a.map(|x| x * b)
}
fn dot(a: V, b: V) -> I {
    (0..3).fold(I::point(0.), |sum, i| sum + a[i] * b[i])
}
fn cross(a: V, b: V) -> V {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn finite(a: V) -> bool {
    a.into_iter().all(I::valid)
}
fn lift(f: RawFrame, p: V) -> V {
    add(
        f.origin,
        add(scale(f.axes[0], p[0]), scale(f.axes[1], p[1])),
    )
}
fn lift_vector(f: RawFrame, p: V) -> V {
    add(scale(f.axes[0], p[0]), scale(f.axes[1], p[1]))
}

// A Lipschitz expansion about one finite point avoids dependence on locating
// trigonometric extrema by a floating quotient when a raw angle is large.
fn trig(value: I, cosine: bool) -> Result<I> {
    if !value.valid() {
        return Err(ReadError::Unsupported);
    }
    let center = value.midpoint();
    if !center.is_finite() {
        return Err(ReadError::Unsupported);
    }
    let radius = (value - I::point(center)).magnitude();
    let at = if cosine {
        I::point(center).cos()
    } else {
        I::point(center).sin()
    };
    let result = I::new(
        (at.lo - radius).next_down().max(-1.),
        (at.hi + radius).next_up().min(1.),
    );
    if result.valid() {
        Ok(result)
    } else {
        Err(ReadError::Unsupported)
    }
}
fn radial(f: RawFrame, u: I) -> Result<V> {
    Ok(add(
        scale(f.axes[0], trig(u, true)?),
        scale(f.axes[1], trig(u, false)?),
    ))
}
fn raw_curve(source: &Source<'_>, id: u32, dim: usize, budget: &mut Budget<'_>) -> Result<Image> {
    budget.charge(1024)?;
    let entity = source.simple(id)?;
    match (entity.name, entity.args.as_slice()) {
        ("LINE", [_, origin, direction]) => Ok(Image::Line {
            origin: point(source, reference(origin)?, dim, budget)?,
            direction: raw_vector(source, reference(direction)?, dim, budget)?,
        }),
        ("CIRCLE", [_, axis, radius]) => {
            let f = frame(source, reference(axis)?, dim, budget)?;
            let r = scalar(source, id, radius, budget)?;
            Ok(Image::Periodic {
                center: f.origin,
                cosine: scale(f.axes[0], r),
                sine: scale(f.axes[1], r),
            })
        }
        ("ELLIPSE", [_, axis, a, b]) => {
            let f = frame(source, reference(axis)?, dim, budget)?;
            Ok(Image::Periodic {
                center: f.origin,
                cosine: scale(f.axes[0], scalar(source, id, a, budget)?),
                sine: scale(f.axes[1], scalar(source, id, b, budget)?),
            })
        }
        _ => Err(ReadError::Unsupported),
    }
}
fn zero_direction_component(
    source: &Source<'_>,
    line_id: u32,
    axis: usize,
    budget: &mut Budget<'_>,
) -> Result<bool> {
    let [_, _, vector_id] = source.args(line_id, "LINE")? else {
        return Err(ReadError::Invalid);
    };
    let [_, direction, _] = source.args(reference(vector_id)?, "VECTOR")? else {
        return Err(ReadError::Invalid);
    };
    let direction = reference(direction)?;
    let [_, values] = source.args(direction, "DIRECTION")? else {
        return Err(ReadError::Invalid);
    };
    let values = list(values)?;
    if values.len() != 2 {
        return Err(ReadError::Unsupported);
    }
    let zero = scalar(source, direction, &values[axis], budget)?;
    Ok(zero.lo == 0. && zero.hi == 0.)
}
/// Actual cardinal-plane affine image with ORIGINAL rational cancellation
/// before conversion. Large in-plane gauges must not lose their correlation
/// by independently rounding the plane and UV anchors into interval boxes.
pub(super) fn plane_point(
    source: &Source<'_>,
    plane_id: u32,
    uv_point_id: u32,
    budget: &mut Budget<'_>,
) -> Result<V> {
    plane_lift(source, plane_id, budget)?.point(source, uv_point_id, budget)
}
/// The raw REAL frame of one `PLANE`, read once and applied to many UV
/// points: exactly what [`plane_point`] reads and lifts for each point.
pub(super) struct PlaneLift {
    frame: RawFrame,
    /// The exact placement origin of a cardinal frame, whose points are
    /// lifted by an exact translation instead.
    cardinal_origin: Option<[BigRational; 3]>,
}
pub(super) fn plane_lift(
    source: &Source<'_>,
    plane_id: u32,
    budget: &mut Budget<'_>,
) -> Result<PlaneLift> {
    budget.charge(1024)?;
    let [_, axis] = source.args(plane_id, "PLANE")? else {
        return Err(ReadError::Invalid);
    };
    let axis = reference(axis)?;
    let f = frame(source, axis, 3, budget)?;
    let cardinal = f
        .axes
        .into_iter()
        .flatten()
        .all(|c| c.lo == c.hi && [-1., 0., 1.].contains(&c.lo));
    let cardinal_origin = if cardinal {
        let [_, origin, _, _] = source.args(axis, "AXIS2_PLACEMENT_3D")? else {
            return Err(ReadError::Invalid);
        };
        Some(translation_origin(source, reference(origin)?, budget)?)
    } else {
        None
    };
    Ok(PlaneLift {
        frame: f,
        cardinal_origin,
    })
}
impl PlaneLift {
    pub(super) fn point(
        &self,
        source: &Source<'_>,
        uv_point_id: u32,
        budget: &mut Budget<'_>,
    ) -> Result<V> {
        budget.charge(64)?;
        if let Some(origin) = &self.cardinal_origin {
            return translated_from(source, origin, uv_point_id, self.frame.axes, budget);
        }
        Ok(lift(self.frame, point(source, uv_point_id, 2, budget)?))
    }
}
fn curve_origin_id(source: &Source<'_>, curve_id: u32) -> Result<u32> {
    let entity = source.simple(curve_id)?;
    match (entity.name, entity.args.as_slice()) {
        ("LINE", [_, origin, _]) => reference(origin),
        ("CIRCLE", [_, axis, _]) | ("ELLIPSE", [_, axis, _, _]) => {
            let [_, origin, _] = source.args(reference(axis)?, "AXIS2_PLACEMENT_2D")? else {
                return Err(ReadError::Invalid);
            };
            reference(origin)
        }
        _ => Err(ReadError::Unsupported),
    }
}
pub(super) fn actual_image(
    source: &Source<'_>,
    definition: &uv::Definition,
    surface: &Surface,
    budget: &mut Budget<'_>,
) -> Result<Image> {
    budget.charge(2048)?;
    let id = definition.surface;
    if let Surface::Plane(_) = surface {
        let [_, axis] = source.args(id, "PLANE")? else {
            return Err(ReadError::Invalid);
        };
        let f = frame(source, reference(axis)?, 3, budget)?;
        let origin = plane_point(
            source,
            id,
            curve_origin_id(source, definition.curve_id)?,
            budget,
        )?;
        return Ok(match raw_curve(source, definition.curve_id, 2, budget)? {
            Image::Line { direction, .. } => Image::Line {
                origin,
                direction: lift_vector(f, direction),
            },
            Image::Periodic { cosine, sine, .. } => Image::Periodic {
                center: origin,
                cosine: lift_vector(f, cosine),
                sine: lift_vector(f, sine),
            },
        });
    }
    // The coordinator separately proves ORIGINAL SI radian ownership for
    // angular surface parameters. A cached unit conversion equal to 1 is not
    // that proof; this additionally refuses a known non-radian reference.
    if source.angle != 1. {
        return Err(ReadError::Unsupported);
    }
    let (moving, origin, rate, decoded_sign) =
        if let Curve::Line { direction, .. } = &definition.geometry {
            let moving = if direction[0] != 0. && direction[1] == 0. {
                0
            } else if direction[0] == 0. && direction[1] != 0. {
                1
            } else {
                return Err(ReadError::Unsupported);
            };
            if !zero_direction_component(source, definition.curve_id, 1 - moving, budget)? {
                return Err(ReadError::Unsupported);
            }
            let Image::Line {
                origin,
                direction: raw_direction,
            } = raw_curve(source, definition.curve_id, 2, budget)?
            else {
                return Err(ReadError::Unsupported);
            };
            (
                moving,
                origin,
                raw_direction[moving],
                direction[moving].signum(),
            )
        } else {
            let segment =
                super::axis_segment(source, definition, budget)?.ok_or(ReadError::Unsupported)?;
            let rate = segment.endpoints[1][segment.moving] - segment.endpoints[0][segment.moving];
            let Curve::Nurbs { curve, .. } = &definition.geometry else {
                return Err(ReadError::Unsupported);
            };
            let decoded = curve.points()[1][segment.moving] - curve.points()[0][segment.moving];
            (segment.moving, segment.endpoints[0], rate, decoded.signum())
        };
    let sign = if rate.valid() && rate.lo > 0. {
        1.
    } else if rate.valid() && rate.hi < 0. {
        -1.
    } else {
        return Err(ReadError::Unsupported);
    };
    if sign != decoded_sign {
        return Err(ReadError::Unsupported);
    }
    // The original fixed coordinate is exactly constant. The moving rate and
    // origin reparameterize its analytic continuation. For finite rational
    // segments ONLY a finite caller earns the actual endpoint window; this
    // coefficient builder never grants natural-domain coverage by itself.
    let entity = source.simple(id)?;
    let result = match (surface, entity.name, entity.args.as_slice()) {
        (Surface::Cylinder { .. }, "CYLINDRICAL_SURFACE", [_, axis, r]) => {
            let f = frame(source, reference(axis)?, 3, budget)?;
            let r = scalar(source, id, r, budget)?;
            if moving == 0 {
                Image::Periodic {
                    center: add(f.origin, scale(f.axes[2], origin[1])),
                    cosine: scale(f.axes[0], r),
                    sine: scale(f.axes[1], r),
                }
            } else {
                Image::Line {
                    origin: add(
                        add(f.origin, scale(radial(f, origin[0])?, r)),
                        scale(f.axes[2], origin[1]),
                    ),
                    direction: scale(f.axes[2], rate),
                }
            }
        }
        (Surface::Cone { .. }, "CONICAL_SURFACE", [_, axis, r, angle]) => {
            let f = frame(source, reference(axis)?, 3, budget)?;
            let r = scalar(source, id, r, budget)?;
            let angle = scalar(source, id, angle, budget)?;
            let cos = trig(angle, true)?;
            if cos.lo <= 0. && cos.hi >= 0. {
                return Err(ReadError::Unsupported);
            }
            let tangent = trig(angle, false)? / cos;
            if moving == 0 {
                let rho = r + origin[1] * tangent;
                Image::Periodic {
                    center: add(f.origin, scale(f.axes[2], origin[1])),
                    cosine: scale(f.axes[0], rho),
                    sine: scale(f.axes[1], rho),
                }
            } else {
                let radial = radial(f, origin[0])?;
                let generator = add(f.axes[2], scale(radial, tangent));
                Image::Line {
                    origin: add(add(f.origin, scale(radial, r)), scale(generator, origin[1])),
                    direction: scale(generator, rate),
                }
            }
        }
        (Surface::Sphere { .. }, "SPHERICAL_SURFACE", [_, axis, r]) => {
            let f = frame(source, reference(axis)?, 3, budget)?;
            let r = scalar(source, id, r, budget)?;
            if moving == 0 {
                if origin[1].lo < -std::f64::consts::FRAC_PI_2
                    || origin[1].hi > std::f64::consts::FRAC_PI_2
                {
                    return Err(ReadError::Unsupported);
                }
                let rho = r * trig(origin[1], true)?;
                Image::Periodic {
                    center: add(f.origin, scale(f.axes[2], r * trig(origin[1], false)?)),
                    cosine: scale(f.axes[0], rho),
                    sine: scale(f.axes[1], rho),
                }
            } else {
                // Coefficients of the analytic continuation ONLY. Sphere V
                // actually covers a half-meridian. A finite-domain caller
                // must separately earn every queried V inside [-pi/2,pi/2].
                Image::Periodic {
                    center: f.origin,
                    cosine: scale(radial(f, origin[0])?, r),
                    sine: scale(f.axes[2], r),
                }
            }
        }
        (Surface::Torus { .. }, "TOROIDAL_SURFACE", [_, axis, major, minor]) => {
            let f = frame(source, reference(axis)?, 3, budget)?;
            let major = scalar(source, id, major, budget)?;
            let minor = scalar(source, id, minor, budget)?;
            if moving == 0 {
                let rho = major + minor * trig(origin[1], true)?;
                Image::Periodic {
                    center: add(f.origin, scale(f.axes[2], minor * trig(origin[1], false)?)),
                    cosine: scale(f.axes[0], rho),
                    sine: scale(f.axes[1], rho),
                }
            } else {
                let radial = radial(f, origin[0])?;
                Image::Periodic {
                    center: add(f.origin, scale(radial, major)),
                    cosine: scale(radial, minor),
                    sine: scale(f.axes[2], minor),
                }
            }
        }
        _ => return Err(ReadError::Unsupported),
    };
    Ok(match result {
        Image::Periodic {
            center,
            cosine,
            sine,
        } => Image::Periodic {
            center,
            cosine,
            sine: scale(sine, I::point(sign)),
        },
        line => line,
    })
}
fn periodic_difference(raw: Image, retained: Image) -> Result<I> {
    let (
        Image::Periodic {
            center: a,
            cosine: u,
            sine: v,
        },
        Image::Periodic {
            center: b,
            cosine: x,
            sine: y,
        },
    ) = (raw, retained)
    else {
        return Err(ReadError::Unsupported);
    };
    if ![a, u, v, b, x, y].into_iter().all(finite) {
        return Err(ReadError::Unsupported);
    }
    let error = norm(difference(a, b)) + norm(difference(u, x)) + norm(difference(v, y));
    if error.valid() {
        Ok(error)
    } else {
        Err(ReadError::Unsupported)
    }
}
fn line_difference(raw: Image, retained: Image, parameter: I) -> Result<I> {
    let (
        Image::Line {
            origin: a,
            direction: u,
        },
        Image::Line {
            origin: b,
            direction: v,
        },
    ) = (raw, retained)
    else {
        return Err(ReadError::Unsupported);
    };
    if ![a, u, b, v].into_iter().all(finite) || !parameter.valid() {
        return Err(ReadError::Unsupported);
    }
    let error = norm(difference(a, b)) + I::point(parameter.magnitude()) * norm(difference(u, v));
    if error.valid() {
        Ok(error)
    } else {
        Err(ReadError::Unsupported)
    }
}
fn edge_enclosure(edge: &Edge, budget: &mut Budget<'_>) -> Result<V> {
    budget.charge(1024)?;
    let parameter = I::new(edge.range[0], edge.range[1]);
    if !parameter.valid() {
        return Err(ReadError::Unsupported);
    }
    let result = match &edge.curve {
        Curve::Line { origin, direction } => {
            add(cached(*origin), scale(cached(*direction), parameter))
        }
        Curve::Circle { frame, radius } => add(
            cached(frame.origin),
            add(
                scale(cached(frame.x), I::new(-1., 1.) * I::point(*radius)),
                scale(cached(frame.y), I::new(-1., 1.) * I::point(*radius)),
            ),
        ),
        Curve::Ellipse {
            frame,
            semiaxes: [a, b],
        } => add(
            cached(frame.origin),
            add(
                scale(cached(frame.x), I::new(-1., 1.) * I::point(*a)),
                scale(cached(frame.y), I::new(-1., 1.) * I::point(*b)),
            ),
        ),
        _ => return Err(ReadError::Unsupported),
    };
    if finite(result) {
        Ok(result)
    } else {
        Err(ReadError::Unsupported)
    }
}
/// A finite line reserve covering BOTH the successful correspondence map and
/// every raw closest-line parameter for ANY retained point on this edge.
/// This latter bound is essential when transferring a point-to-WHOLE-line
/// lower witness; there is no uniform infinite-line Hausdorff bound.
pub(super) fn line_image_reserve(
    raw: Image,
    reference: Image,
    edge: &Edge,
    budget: &mut Budget<'_>,
) -> Result<f64> {
    budget.charge(2048)?;
    let (
        Image::Line {
            origin: raw_origin,
            direction: raw_direction,
        },
        Image::Line { origin, direction },
        Curve::Line {
            origin: p,
            direction: d,
        },
    ) = (raw, reference, &edge.curve)
    else {
        return Err(ReadError::Unsupported);
    };
    let denominator = dot(direction, direction);
    if !denominator.valid() || denominator.lo <= 0. {
        return Err(ReadError::Unsupported);
    }
    let alpha = dot(cached(*d), direction) / denominator;
    let beta = dot(difference(cached(*p), origin), direction) / denominator;
    if !alpha.valid() || !beta.valid() {
        return Err(ReadError::Unsupported);
    }
    let parameter = alpha * I::new(edge.range[0], edge.range[1]) + beta;
    if !parameter.valid() {
        return Err(ReadError::Unsupported);
    }
    // Include the exact floating proposal used by the analytic certificate as
    // well as the interval projection used by the affine/sweep certificate.
    // A genuine lower witness does not require same sense. The positive
    // admission path separately earns IP2; a negative-slope reserve remains
    // valid for transferring an actual geometric discrepancy.
    let proposed = if let Some([a, b]) = retained::line_map(reference, edge, budget)? {
        let value = I::point(a) * I::new(edge.range[0], edge.range[1]) + I::point(b);
        if !value.valid() {
            return Err(ReadError::Unsupported);
        }
        value.magnitude()
    } else {
        0.
    };
    let points = edge_enclosure(edge, budget)?;
    let distance = norm(difference(points, raw_origin));
    let length = norm(raw_direction);
    if !distance.valid() || !length.valid() || length.lo <= 0. {
        return Err(ReadError::Unsupported);
    }
    // For each retained p, the raw closest parameter satisfies
    // |q*| <= ||p - raw_origin|| / ||raw_direction||.
    let closest = distance / length;
    if !closest.valid() {
        return Err(ReadError::Unsupported);
    }
    let bound = parameter.magnitude().max(proposed).max(closest.hi);
    let reserve = line_difference(raw, reference, I::new(-bound, bound))?;
    Ok(reserve.hi)
}
/// Raw original PCURVE/surface image to the SAME retained image. For a finite
/// LINE this uses the exact positive map already proposed by its certificate.
/// Periodic coefficients cover the complete natural image and preserve sense.
pub(super) fn image_reserve(
    source: &Source<'_>,
    definition: &uv::Definition,
    surface: &Surface,
    edge: &Edge,
    budget: &mut Budget<'_>,
) -> Result<f64> {
    budget.charge(2048)?;
    let cached_image =
        retained::image(surface, &definition.geometry, budget)?.ok_or(ReadError::Unsupported)?;
    let raw = actual_image(source, definition, surface, budget)?;
    let reserve = match cached_image {
        Image::Periodic { .. } => periodic_difference(raw, cached_image)?,
        Image::Line { .. } => return line_image_reserve(raw, cached_image, edge, budget),
    };
    Ok(reserve.hi)
}
/// Original analytic 3D curve to its retained reference over the actual finite
/// 3D edge parameter range. Periodic shapes use a global coefficient envelope.
pub(super) fn curve_reserve(
    source: &Source<'_>,
    id: u32,
    edge: &Edge,
    budget: &mut Budget<'_>,
) -> Result<f64> {
    Ok(curve_packet(source, id, edge, budget)?.1)
}
/// The immutable Original image and its numeric transfer are prepared once
/// for one finite association. The caller retains this packet for its sense
/// check instead of re-reading the same source operators.
pub(super) fn curve_packet(
    source: &Source<'_>,
    id: u32,
    edge: &Edge,
    budget: &mut Budget<'_>,
) -> Result<(Image, f64)> {
    budget.charge(2048)?;
    let id = unwrapped(source, id, budget)?;
    let raw = raw_curve(source, id, 3, budget)?;
    let reserve = match &edge.curve {
        Curve::Line { origin, direction } => line_difference(
            raw,
            Image::Line {
                origin: cached(*origin),
                direction: cached(*direction),
            },
            I::new(edge.range[0], edge.range[1]),
        )?,
        Curve::Circle { frame, radius } => periodic_difference(
            raw,
            Image::Periodic {
                center: cached(frame.origin),
                cosine: scale(cached(frame.x), I::point(*radius)),
                sine: scale(cached(frame.y), I::point(*radius)),
            },
        )?,
        Curve::Ellipse {
            frame,
            semiaxes: [a, b],
        } => periodic_difference(
            raw,
            Image::Periodic {
                center: cached(frame.origin),
                cosine: scale(cached(frame.x), I::point(*a)),
                sine: scale(cached(frame.y), I::point(*b)),
            },
        )?,
        _ => return Err(ReadError::Unsupported),
    };
    Ok((raw, reserve.hi))
}

/// Original REAL IP2 sense is checked independently of positional tolerance.
/// This is needed only before positive admission, not to invalidate an actual
/// point-set discrepancy witness whose orientation is immaterial.
pub(super) fn sense(
    source: &Source<'_>,
    curve3d_id: u32,
    definition: &uv::Definition,
    surface: &Surface,
    _edge: &Edge,
    budget: &mut Budget<'_>,
) -> Result<()> {
    let id = unwrapped(source, curve3d_id, budget)?;
    let target = raw_curve(source, id, 3, budget)?;
    let image = actual_image(source, definition, surface, budget)?;
    sense_images(target, image, budget)
}
/// Only the alignment arithmetic is shared. Preparing the two Original
/// packets and earning their finite correspondence remain caller obligations.
pub(super) fn sense_images(target: Image, image: Image, budget: &mut Budget<'_>) -> Result<()> {
    budget.charge(2048)?;
    let alignment = match (target, image) {
        (Image::Line { direction: a, .. }, Image::Line { direction: b, .. }) => dot(a, b),
        (
            Image::Periodic {
                cosine: a, sine: b, ..
            },
            Image::Periodic {
                cosine: c, sine: d, ..
            },
        ) => dot(cross(a, b), cross(c, d)),
        _ => return Err(ReadError::Unsupported),
    };
    if !alignment.valid() || alignment.lo <= 0. {
        return Err(ReadError::Unsupported);
    }
    Ok(())
}

#[cfg(test)]
mod tests;
