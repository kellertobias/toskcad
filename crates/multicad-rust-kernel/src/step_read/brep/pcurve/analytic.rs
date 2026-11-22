//! Natural-domain analytic trim images, with parameter-independent coverage.
//!
//! A LINE domain is all real numbers and conics here have their full natural
//! domain. Bounded NURBS and TRIMMED_CURVE are deliberately not admitted by
//! a support-only comparison. Every successful map below is surjective onto
//! the analytic support; the finite 3D edge selects its corresponding UV arc.
use super::{Budget, Proof, Result, difference, norm, vector};
use crate::brep::{Curve, Edge, Frame, Point, Surface};
use crate::tapered_thread::interval::I;

type V = [I; 3];
#[derive(Clone, Copy)]
pub(super) enum Image {
    Line { origin: V, direction: V },
    Periodic { center: V, cosine: V, sine: V },
}
fn add(a: V, b: V) -> V {
    std::array::from_fn(|i| a[i] + b[i])
}
fn scale(a: V, b: I) -> V {
    a.map(|v| v * b)
}
fn dot(a: V, b: V) -> I {
    (0..3).fold(I::point(0.), |s, i| s + a[i] * b[i])
}
fn cross(a: V, b: V) -> V {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn midpoint(a: V) -> Point {
    a.map(I::midpoint)
}
fn valid(a: V) -> bool {
    a.into_iter().all(I::valid)
}
fn absolute_lower(a: I) -> f64 {
    if a.lo > 0. {
        a.lo
    } else if a.hi < 0. {
        -a.hi
    } else {
        0.
    }
}
fn lift(frame: &Frame, p: Point) -> V {
    add(
        vector(frame.origin),
        add(
            scale(vector(frame.x), I::point(p[0])),
            scale(vector(frame.y), I::point(p[1])),
        ),
    )
}
fn lift_vector(frame: &Frame, p: Point) -> V {
    add(
        scale(vector(frame.x), I::point(p[0])),
        scale(vector(frame.y), I::point(p[1])),
    )
}
fn radial(frame: &Frame, u: f64) -> V {
    let u = I::point(u);
    add(
        scale(vector(frame.x), u.cos()),
        scale(vector(frame.y), u.sin()),
    )
}
fn physical_periodic(curve: &Curve) -> Option<Image> {
    match curve {
        Curve::Circle { frame, radius } => Some(Image::Periodic {
            center: vector(frame.origin),
            cosine: scale(vector(frame.x), I::point(*radius)),
            sine: scale(vector(frame.y), I::point(*radius)),
        }),
        Curve::Ellipse {
            frame,
            semiaxes: [a, b],
        } => Some(Image::Periodic {
            center: vector(frame.origin),
            cosine: scale(vector(frame.x), I::point(*a)),
            sine: scale(vector(frame.y), I::point(*b)),
        }),
        _ => None,
    }
}
pub(super) fn image(
    surface: &Surface,
    uv: &Curve,
    budget: &mut Budget<'_>,
) -> Result<Option<Image>> {
    // Fixed-size coefficient construction: at most two trigonometric pairs,
    // one division and 3-vector arithmetic; no heap or iterative evaluator.
    budget.charge(1024)?;
    if let Surface::Plane(frame) = surface {
        return Ok(match uv {
            Curve::Line { origin, direction } => Some(Image::Line {
                origin: lift(frame, *origin),
                direction: lift_vector(frame, *direction),
            }),
            Curve::Circle { frame: uv, radius } => Some(Image::Periodic {
                center: lift(frame, uv.origin),
                cosine: scale(lift_vector(frame, uv.x), I::point(*radius)),
                sine: scale(lift_vector(frame, uv.y), I::point(*radius)),
            }),
            Curve::Ellipse {
                frame: uv,
                semiaxes: [a, b],
            } => Some(Image::Periodic {
                center: lift(frame, uv.origin),
                cosine: scale(lift_vector(frame, uv.x), I::point(*a)),
                sine: scale(lift_vector(frame, uv.y), I::point(*b)),
            }),
            _ => None,
        });
    }
    let Curve::Line {
        origin: [u, v, _],
        direction: [du, dv, _],
    } = uv
    else {
        return Ok(None);
    };
    // Exact zero means an actual isocurve. Near-zero diagonal lines are not
    // rounded into isocurves, whatever the requested geometric tolerance.
    let along_u = *du != 0. && *dv == 0.;
    let along_v = *du == 0. && *dv != 0.;
    if !along_u && !along_v {
        return Ok(None);
    }
    let result = match surface {
        Surface::Cylinder { frame, radius } => {
            if along_u {
                Image::Periodic {
                    center: add(vector(frame.origin), scale(vector(frame.z), I::point(*v))),
                    cosine: scale(vector(frame.x), I::point(*radius)),
                    sine: scale(vector(frame.y), I::point(*radius)),
                }
            } else {
                Image::Line {
                    origin: add(
                        add(
                            vector(frame.origin),
                            scale(radial(frame, *u), I::point(*radius)),
                        ),
                        scale(vector(frame.z), I::point(*v)),
                    ),
                    direction: scale(vector(frame.z), I::point(*dv)),
                }
            }
        }
        Surface::Cone {
            frame,
            radius,
            half_angle,
        } => {
            let angle = I::point(*half_angle);
            let cosine = angle.cos();
            if cosine.lo <= 0. && cosine.hi >= 0. {
                return Ok(None);
            }
            let slope = angle.sin() / cosine;
            if along_u {
                let rho = I::point(*radius) + I::point(*v) * slope;
                Image::Periodic {
                    center: add(vector(frame.origin), scale(vector(frame.z), I::point(*v))),
                    cosine: scale(vector(frame.x), rho),
                    sine: scale(vector(frame.y), rho),
                }
            } else {
                let r = radial(frame, *u);
                let tangent = add(scale(r, slope), vector(frame.z));
                Image::Line {
                    origin: add(
                        add(vector(frame.origin), scale(r, I::point(*radius))),
                        scale(tangent, I::point(*v)),
                    ),
                    direction: scale(tangent, I::point(*dv)),
                }
            }
        }
        Surface::Sphere { frame, radius } => {
            // STEP spherical V is restricted to [-pi/2,pi/2]. A V-line
            // therefore covers only HALF a great circle; support equality
            // with a full circle would not earn finite-edge coverage.
            if !along_u || v.abs() > std::f64::consts::FRAC_PI_2 {
                return Ok(None);
            }
            let r = I::point(*radius);
            let v = I::point(*v);
            let rho = r * v.cos();
            Image::Periodic {
                center: add(vector(frame.origin), scale(vector(frame.z), r * v.sin())),
                cosine: scale(vector(frame.x), rho),
                sine: scale(vector(frame.y), rho),
            }
        }
        Surface::Torus {
            frame,
            major,
            minor,
        } => {
            let m = I::point(*minor);
            if along_u {
                let v = I::point(*v);
                let rho = I::point(*major) + m * v.cos();
                Image::Periodic {
                    center: add(vector(frame.origin), scale(vector(frame.z), m * v.sin())),
                    cosine: scale(vector(frame.x), rho),
                    sine: scale(vector(frame.y), rho),
                }
            } else {
                let r = radial(frame, *u);
                Image::Periodic {
                    center: add(vector(frame.origin), scale(r, I::point(*major))),
                    cosine: scale(r, m),
                    sine: scale(vector(frame.z), m),
                }
            }
        }
        _ => return Ok(None),
    };
    // lambda = |moving rate| * q + constant is increasing in the ACTUAL
    // UV-curve parameter q. Retain its orientation in the periodic basis.
    let result = match result {
        Image::Periodic {
            center,
            cosine,
            sine,
        } => {
            let sign = if along_u { du.signum() } else { dv.signum() };
            Image::Periodic {
                center,
                cosine,
                sine: scale(sine, I::point(sign)),
            }
        }
        line => line,
    };
    Ok(Some(result))
}
pub(super) fn line_map(
    image: Image,
    edge: &Edge,
    budget: &mut Budget<'_>,
) -> Result<Option<[f64; 2]>> {
    budget.charge(512)?;
    let (
        Image::Line {
            origin: q,
            direction: d,
        },
        Curve::Line {
            origin: p,
            direction: v,
        },
    ) = (image, &edge.curve)
    else {
        return Ok(None);
    };
    if !valid(q) || !valid(d) {
        return Ok(None);
    }
    let (qm, dm) = (midpoint(q), midpoint(d));
    let length2 = dm.into_iter().map(|x| x * x).sum::<f64>();
    let directed_alpha = dot(vector(*v), d) / dot(d, d);
    if !directed_alpha.valid() || directed_alpha.lo <= 0. {
        return Ok(None);
    }
    let alpha = (0..3).map(|i| v[i] * dm[i]).sum::<f64>() / length2;
    let beta = (0..3).map(|i| (p[i] - qm[i]) * dm[i]).sum::<f64>() / length2;
    if !alpha.is_finite() || alpha <= 0. || !beta.is_finite() || norm(d).lo <= 0. {
        return Ok(None);
    }
    Ok(Some([alpha, beta]))
}
fn line_upper(image: Image, edge: &Edge, budget: &mut Budget<'_>) -> Result<Option<f64>> {
    budget.charge(1024)?;
    let Some([alpha, beta]) = line_map(image, edge, budget)? else {
        return Ok(None);
    };
    let (
        Image::Line {
            origin: q,
            direction: d,
        },
        Curve::Line {
            origin: p,
            direction: v,
        },
    ) = (image, &edge.curve)
    else {
        return Ok(None);
    };
    // The PROPOSED real map t=alpha*s+beta is increasing and bijective.
    // The independently enclosed actual projected tangent slope is positive.
    // Its residual is affine, hence its norm on the finite actual edge is
    // bounded by the maximum of the two endpoint interval norms.
    let residual_origin = difference(difference(vector(*p), q), scale(d, I::point(beta)));
    let residual_direction = difference(vector(*v), scale(d, I::point(alpha)));
    let mut upper = 0f64;
    for s in edge.range {
        let error = norm(add(residual_origin, scale(residual_direction, I::point(s))));
        // max() must not mask NaN/overflow from a failed norm enclosure.
        if !error.valid() {
            return Ok(None);
        }
        upper = upper.max(error.hi);
    }
    Ok(Some(upper))
}
pub(super) struct PeriodicMap {
    pub phase: f64,
    pub error: f64,
    pub speed: f64,
}
pub(super) fn periodic_map(
    image: Image,
    edge: &Edge,
    phase_shift: f64,
    budget: &mut Budget<'_>,
) -> Result<Option<PeriodicMap>> {
    budget.charge(2048)?;
    let Image::Periodic {
        center: b,
        cosine: u,
        sine: v,
    } = image
    else {
        return Ok(None);
    };
    let Some(Image::Periodic {
        center: a,
        cosine: x,
        sine: y,
    }) = physical_periodic(&edge.curve)
    else {
        return Ok(None);
    };
    if !valid(b) || !valid(u) || !valid(v) || norm(cross(u, v)).lo <= 0. {
        return Ok(None);
    }
    // Sense is independent of tolerance: even a tiny opposite-sense loop
    // cannot become valid merely because its diameter fits the C0 limit.
    let orientation = dot(cross(x, y), cross(u, v));
    if !orientation.valid() || orientation.lo <= 0. {
        return Ok(None);
    }
    // Solve a small Gram system only to propose a phase. Its accuracy never
    // authorizes geometry: all coefficients are checked independently below.
    let (um, vm, xm) = (midpoint(u), midpoint(v), midpoint(x));
    let d = |a: Point, b: Point| (0..3).map(|i| a[i] * b[i]).sum::<f64>();
    let (uu, uv, vv) = (d(um, um), d(um, vm), d(vm, vm));
    let determinant = uu * vv - uv * uv;
    let (xu, xv) = (d(xm, um), d(xm, vm));
    let c = (xu * vv - xv * uv) / determinant;
    let s = (xv * uu - xu * uv) / determinant;
    let phi = s.atan2(c) + phase_shift;
    if !phi.is_finite() {
        return Ok(None);
    }
    let (cosine, sine) = (I::point(phi).cos(), I::point(phi).sin());
    let first = add(scale(u, cosine), scale(v, sine));
    let second = add(scale(u, -sine), scale(v, cosine));
    let center = norm(difference(a, b));
    let error = center + norm(difference(x, first)) + norm(difference(y, second));
    // theta_source = theta_target+phi is increasing and surjective on the
    // full periodic supports. A reflected map is NEVER an IP2 certificate.
    // The bound uses |sin theta|,|cos theta|<=1 over the WHOLE support.
    let speed = (norm(u) + norm(v)).hi;
    Ok((error.valid() && speed.is_finite()).then_some(PeriodicMap {
        phase: phi,
        error: error.hi,
        speed,
    }))
}
fn periodic_upper(image: Image, edge: &Edge, budget: &mut Budget<'_>) -> Result<Option<f64>> {
    Ok(periodic_map(image, edge, 0., budget)?.map(|map| map.error))
}

fn physical_point(edge: &Edge, t: f64) -> Option<V> {
    match &edge.curve {
        Curve::Line { origin, direction } => {
            Some(add(vector(*origin), scale(vector(*direction), I::point(t))))
        }
        curve => {
            let Image::Periodic {
                center,
                cosine,
                sine,
            } = physical_periodic(curve)?
            else {
                return None;
            };
            Some(add(
                center,
                add(
                    scale(cosine, I::point(t).cos()),
                    scale(sine, I::point(t).sin()),
                ),
            ))
        }
    }
}
fn support_lower(image: Image, point: V) -> f64 {
    match image {
        Image::Line { origin, direction } => {
            let denominator = norm(direction);
            if !denominator.valid() || denominator.lo <= 0. {
                return 0.;
            }
            let numerator = norm(cross(difference(point, origin), direction));
            (numerator.lo / denominator.hi).next_down().max(0.)
        }
        Image::Periodic {
            center,
            cosine: u,
            sine: v,
        } => {
            // Gershgorin encloses the Gram eigenvalues, hence the radius of
            // EVERY point in this actual (possibly slightly skew) ellipse.
            let (uu, uv, vv) = (dot(u, u), dot(u, v), dot(v, v));
            let upper = (I::point(uu.hi.max(vv.hi)) + I::point(uv.magnitude())).sqrt();
            let lower_squared = (uu.lo.min(vv.lo) - uv.magnitude()).next_down().max(0.);
            let lower = I::point(lower_squared).sqrt();
            let delta = difference(point, center);
            let distance = norm(delta);
            if !upper.valid() || !lower.valid() || !distance.valid() {
                return 0.;
            }
            let radial = (distance.lo - upper.hi)
                .next_down()
                .max((lower.lo - distance.hi).next_down())
                .max(0.);
            // This additional witness is distance to the ENTIRE ellipse's
            // containing plane, not a same-parameter curve residual.
            let normal = cross(u, v);
            let length = norm(normal);
            let plane = if length.valid() && length.lo > 0. {
                (absolute_lower(dot(delta, normal)) / length.hi)
                    .next_down()
                    .max(0.)
            } else {
                0.
            };
            radial.max(plane)
        }
    }
}
/// The caller must earn ownership and the complete natural image/domain of
/// these coefficients, including positive source parameter orientation. This
/// shared mathematical checker does not infer a STEP source from coefficients.
pub(super) fn certify_image(
    image: Image,
    edge: &Edge,
    tolerance: f64,
    budget: &mut Budget<'_>,
) -> Result<Proof> {
    if !tolerance.is_finite()
        || tolerance <= 0.
        || edge.range.into_iter().any(|x| !x.is_finite())
        || edge.range[0] >= edge.range[1]
    {
        return Ok(Proof::Unsupported);
    }
    let finite = match image {
        Image::Line { origin, direction } => valid(origin) && valid(direction),
        Image::Periodic {
            center,
            cosine,
            sine,
        } => valid(center) && valid(cosine) && valid(sine),
    };
    if !finite {
        return Ok(Proof::Unsupported);
    }
    let upper = match image {
        Image::Line { .. } => line_upper(image, edge, budget)?,
        Image::Periodic { .. } => periodic_upper(image, edge, budget)?,
    };
    if let Some(upper) = upper
        && upper.is_finite()
        && upper <= tolerance
    {
        return Ok(Proof::Certified(upper));
    }
    // A few points may prove a lower witness, but can NEVER prove agreement.
    // Each point is interval-evaluated from the actual 3D analytic operator.
    let middle = edge.range[0].midpoint(edge.range[1]);
    for parameter in [edge.range[0], middle, edge.range[1]] {
        budget.charge(1024)?;
        if let Some(point) = physical_point(edge, parameter)
            && valid(point)
        {
            let lower = support_lower(image, point);
            if lower.is_finite() && lower > tolerance {
                return Ok(Proof::Disagreement(lower));
            }
        }
    }
    Ok(Proof::Unsupported)
}
/// Compare a caller-earned full natural periodic image. The caller must derive
/// these SAME actual source coefficients, oriented for POSITIVE source UV
/// parameter, and prove full-period coverage; this
/// mathematical helper does not create source ownership or a trim domain.
pub(super) fn periodic(
    center: [I; 3],
    cosine: [I; 3],
    sine: [I; 3],
    edge: &Edge,
    tolerance: f64,
    budget: &mut Budget<'_>,
) -> Result<Proof> {
    budget.charge(64)?;
    certify_image(
        Image::Periodic {
            center,
            cosine,
            sine,
        },
        edge,
        tolerance,
        budget,
    )
}
pub(super) fn certify(
    surface: &Surface,
    uv: &Curve,
    edge: &Edge,
    tolerance: f64,
    budget: &mut Budget<'_>,
) -> Result<Proof> {
    budget.charge(64)?;
    let Some(image) = image(surface, uv, budget)? else {
        return Ok(Proof::Unsupported);
    };
    certify_image(image, edge, tolerance, budget)
}

#[cfg(test)]
mod tests;
