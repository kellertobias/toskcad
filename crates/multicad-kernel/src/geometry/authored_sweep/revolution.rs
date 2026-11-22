//! Exact revolutions of certified curved meridians about an axis in their
//! sketch plane, for signed partial and full turns.
//!
//! Each meridian curve sweeps its own exact support: lines become planes,
//! cylinders or cones, circular arcs become tori and splines become surfaces
//! of revolution of the source spline. Every meridian vertex traces one
//! circular rail. A partial turn keeps planar start and end caps; a full
//! turn closes each wall on its meridian seam, and each hole becomes a void.
use super::*;
use crate::RevolutionAxis;

const REVOLVE: &str = "Rust curved revolution";

/// Rotation about the unit axis `a` through `origin` by `angle`.
fn rotation(origin: P, a: P, angle: f64) -> Affine {
    let (s, c) = angle.portable_sin_cos();
    let m: [P; 3] = std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            let identity = f64::from(u8::from(i == j));
            let skew = match (i, j) {
                (0, 1) => -a[2],
                (0, 2) => a[1],
                (1, 0) => a[2],
                (1, 2) => -a[0],
                (2, 0) => -a[1],
                (2, 1) => a[0],
                _ => 0.,
            };
            c * identity + s * skew + (1. - c) * a[i] * a[j]
        })
    });
    let turned = m.map(|row| dot(row, origin));
    Affine {
        m,
        t: sub(origin, turned),
    }
}

/// The point at the middle of a path and the unit tangent in its authored
/// direction there.
fn middle(path: &Path, tolerance: f64) -> Result<(P, P), GeometryError> {
    let (point, tangent) = match path {
        Path::Line([a, b]) => (scale(add(*a, *b), 0.5), sub(*b, *a)),
        Path::Conic {
            center,
            axes,
            range,
        } => {
            let t = 0.5 * (range[0] + range[1]);
            let (sin, cos) = t.portable_sin_cos();
            (
                add(*center, add(scale(axes[0], cos), scale(axes[1], sin))),
                add(scale(axes[0], -sin), scale(axes[1], cos)),
            )
        }
        Path::Spline { .. } => {
            let Curve::Nurbs { curve, .. } = path.curve(tolerance)? else {
                unreachable!("spline path");
            };
            let [a, b] = curve.range();
            curve.evaluate(0.5 * (a + b))
        }
    };
    Ok((point, unit(tangent)?))
}

/// The meridian frame of a revolution: axis point, unit axis and the radial
/// direction on the material's side.
#[derive(Clone, Copy)]
struct Meridian {
    origin: P,
    axis: P,
    radial: P,
}

impl Meridian {
    fn coordinates(self, point: P) -> [f64; 2] {
        let d = sub(point, self.origin);
        [dot(d, self.radial), dot(d, self.axis)]
    }

    fn radial_at(self, angle: f64) -> P {
        let (s, c) = angle.portable_sin_cos();
        add(
            scale(self.radial, c),
            scale(cross(self.axis, self.radial), s),
        )
    }

    fn frame(self, height: f64) -> Frame {
        Frame {
            origin: add(self.origin, scale(self.axis, height)),
            x: self.radial,
            y: cross(self.axis, self.radial),
            z: self.axis,
        }
    }

    /// The exact support swept by one meridian path.
    fn surface(self, path: &Path, tolerance: f64) -> Result<Surface, GeometryError> {
        match path {
            Path::Line([a, b]) => {
                let ([r0, z0], [r1, z1]) = (self.coordinates(*a), self.coordinates(*b));
                let scale_mm = (r1 - r0).portable_hypot(z1 - z0);
                if (r1 - r0).abs() <= 1e-12 * scale_mm {
                    return Ok(Surface::Cylinder {
                        frame: self.frame(z0),
                        radius: 0.5 * (r0 + r1),
                    });
                }
                if (z1 - z0).abs() <= 1e-12 * scale_mm {
                    let mut frame = self.frame(0.5 * (z0 + z1));
                    frame.origin = *a;
                    return Ok(Surface::Plane(frame));
                }
                // The cone opens from its narrower circle.
                let (z, radius, height) = if r1 > r0 {
                    (z0, r0, z1 - z0)
                } else {
                    (z1, r1, z0 - z1)
                };
                let mut frame = self.frame(z);
                if height < 0. {
                    frame.z = frame.z.map(|x| -x);
                    frame.y = cross(frame.z, frame.x);
                }
                Ok(Surface::Cone {
                    frame,
                    radius,
                    half_angle: ((r1 - r0).abs() / height.abs()).portable_atan(),
                })
            }
            Path::Conic { center, axes, .. } => {
                let [rc, zc] = self.coordinates(*center);
                if rc.abs() <= tolerance {
                    return Ok(Surface::Sphere {
                        frame: self.frame(zc),
                        radius: norm(axes[0]),
                    });
                }
                Ok(Surface::Torus {
                    frame: self.frame(zc),
                    major: rc,
                    minor: norm(axes[0]),
                })
            }
            Path::Spline { .. } => Ok(Surface::Revolution {
                curve: path.curve(tolerance)?,
                origin: self.origin,
                axis: self.axis,
            }),
        }
    }
}

impl Solid {
    /// Revolve a certified section whose meridian holds circular arcs or
    /// splines. The section lies strictly on one side of the axis.
    pub(in crate::geometry) fn rust_authored_revolution(
        profile: &CurvedProfile,
        axis: RevolutionAxis,
        angle: f64,
        feature: EntityId,
        origins: &[Origin],
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        policy.validate()?;
        axis.validate(profile, angle, policy)?;
        certify(profile, feature, origins, policy)?;
        let tolerance = policy.linear_tolerance_mm;
        let frame = profile.frame;
        let a = unit(axis.direction)?;
        let mut radial = unit(cross(a, frame.normal()))?;
        // Radius is affine in sketch coordinates; bound it on the true curves.
        let reach = |radial: P| {
            let offset = dot(sub(frame.origin, axis.origin_mm), radial);
            linear_range(
                &profile.outer,
                [dot(radial, frame.u), dot(radial, frame.v)],
                policy,
            )
            .map(|[low, high]| [offset + low, offset + high])
        };
        let mut range = reach(radial)?;
        if range[1] <= 0. {
            radial = radial.map(|x| -x);
            range = reach(radial)?;
        }
        if range[0] < -tolerance {
            return Err(input("Revolution profile crosses the axis"));
        }
        let meridian = Meridian {
            origin: axis.origin_mm,
            axis: a,
            radial,
        };
        let full = angle.abs() == TAU;
        let turn = rotation(axis.origin_mm, a, angle);
        // Orientation in meridian coordinates (radius, height).
        let handed = dot(cross(radial, a), frame.normal()).signum();
        let mut first = 0;
        let rings: Vec<Ring> = std::iter::once(&profile.outer)
            .chain(&profile.holes)
            .enumerate()
            .map(|(index, wire)| {
                let area = wire
                    .iter()
                    .map(|edge| twice_signed_area(edge, policy))
                    .sum::<Result<f64, _>>()?;
                let ring = Ring {
                    first,
                    len: wire.len(),
                    forward: (area * handed > 0.) == (index == 0),
                    outer: index == 0,
                };
                first += wire.len();
                Ok(ring)
            })
            .collect::<Result<_, GeometryError>>()?;
        let start: Vec<_> = edges(profile).map(|e| Path::sketch(frame, e)).collect();
        let end: Vec<_> = start.iter().map(|path| path.mapped(turn)).collect();
        let count = start.len();
        let on_axis = |p: P| meridian.coordinates(p)[0] <= tolerance;
        let mut along_axis = vec![false; count];
        for (index, (path, edge)) in start.iter().zip(edges(profile)).enumerate() {
            let [low, _] = linear_range(
                std::slice::from_ref(edge),
                [dot(radial, frame.u), dot(radial, frame.v)],
                policy,
            )?;
            if low + dot(sub(frame.origin, axis.origin_mm), radial) > tolerance {
                continue;
            }
            let ends = path.ends(tolerance)?;
            let touching = ends.map(on_axis);
            along_axis[index] = matches!(path, Path::Line(_)) && touching == [true, true];
            // Lines reach the axis only at an end. A curve may end on it, but
            // a tangent or interior contact would pinch the solid.
            let interior = match path {
                Path::Line(_) => false,
                Path::Conic {
                    center,
                    axes,
                    range,
                } => {
                    let along = [dot(axes[0], radial), dot(axes[1], radial)];
                    let nearest = (-along[1]).portable_atan2(-along[0]);
                    let travel = (nearest - range[0]).rem_euclid(TAU);
                    let inside = travel > 1e-9 && travel < range[1] - range[0] - 1e-9;
                    let (sin, cos) = nearest.portable_sin_cos();
                    let reach = meridian.coordinates(*center)[0] + along[0] * cos + along[1] * sin;
                    inside && reach <= tolerance
                }
                Path::Spline { .. } => true,
            };
            if interior || !touching.iter().any(|&t| t) {
                return Err(error(
                    REVOLVE,
                    GeometryErrorKind::Unsupported,
                    "A meridian curve meets the axis away from its ends; split it there",
                ));
            }
        }
        // The axis may bound only the outer wire, along one connected run.
        for ring in &rings {
            let touching: Vec<bool> = (0..ring.len)
                .map(|k| {
                    let edge = ring.edge(k);
                    start[edge]
                        .ends(tolerance)
                        .map(|ends| on_axis(ends[usize::from(!ring.forward)]))
                })
                .collect::<Result<_, _>>()?;
            let runs = (0..ring.len)
                .filter(|&k| {
                    let previous = (k + ring.len - 1) % ring.len;
                    touching[k] && !(touching[previous] && along_axis[ring.edge(previous)])
                })
                .count();
            if (!ring.outer && runs > 0) || runs > 1 {
                return Err(input(
                    "Revolution profile touches the axis at separate places and would pinch",
                ));
            }
        }
        let surfaces = start
            .iter()
            .map(|path| meridian.surface(path, tolerance))
            .collect::<Result<Vec<_>, _>>()?;
        // A full turn of a radial line is a planar disc or annulus bounded
        // by its rails alone; it keeps no seam.
        let seamless: Vec<bool> = (0..count)
            .map(|edge| full && (along_axis[edge] || surfaces[edge].is_plane()))
            .collect();

        let mut boundary = Boundary::default();
        let mut ring_vertices = Vec::with_capacity(rings.len());
        for ring in &rings {
            let mut vertices = Vec::with_capacity(ring.len);
            for k in 0..ring.len {
                let edge = ring.edge(k);
                let point = start[edge].ends(tolerance)?[usize::from(!ring.forward)];
                boundary.vertices.push(point);
                let low = boundary.vertices.len() - 1;
                let high = if full || on_axis(point) {
                    low
                } else {
                    boundary.vertices.push(turn.point(point));
                    low + 1
                };
                vertices.push([low, high]);
            }
            ring_vertices.push(vertices);
        }
        // Meridian edges at both ends; a full turn shares one seam.
        let mut meridians = [vec![usize::MAX; count], vec![usize::MAX; count]];
        let mut rails = Vec::with_capacity(rings.len());
        for (ring, vertices) in rings.iter().zip(&ring_vertices) {
            for k in 0..ring.len {
                let edge = ring.edge(k);
                if seamless[edge] {
                    continue;
                }
                let (here, next) = (vertices[k], vertices[(k + 1) % ring.len]);
                // An edge on the axis is shared by both caps.
                let shared = full || along_axis[edge];
                for side in 0..if shared { 1 } else { 2 } {
                    let path = if side == 0 { &start[edge] } else { &end[edge] };
                    let [from, to] = if ring.forward {
                        [here[side], next[side]]
                    } else {
                        [next[side], here[side]]
                    };
                    let points = [from, to].map(|v| boundary.vertices[v]);
                    boundary.edges.push(
                        brep::Edge::new(path.curve(tolerance)?, [from, to], points, tolerance)
                            .map_err(failure)?,
                    );
                    meridians[side][edge] = boundary.edges.len() - 1;
                }
                if shared {
                    meridians[1][edge] = meridians[0][edge];
                }
            }
            // Rails run counter-clockwise about the axis; `toward` says
            // whether that is the sweep direction.
            let mut ring_rails = Vec::with_capacity(ring.len);
            for &[low, high] in vertices {
                let point = boundary.vertices[low];
                if on_axis(point) {
                    ring_rails.push(None);
                    continue;
                }
                let [radius, height] = meridian.coordinates(point);
                let center = add(axis.origin_mm, scale(a, height));
                let x = unit(sub(point, center))?;
                let circle = Curve::Circle {
                    frame: Frame {
                        origin: center,
                        x,
                        y: cross(a, x),
                        z: a,
                    },
                    radius,
                };
                let ends = if angle > 0. { [low, high] } else { [high, low] };
                let points = ends.map(|v| boundary.vertices[v]);
                boundary
                    .edges
                    .push(brep::Edge::new(circle, ends, points, tolerance).map_err(failure)?);
                ring_rails.push(Some(boundary.edges.len() - 1));
            }
            rails.push(ring_rails);
        }
        let toward = angle > 0.;

        let mut faces = Vec::with_capacity(count + 2);
        let mut names = Vec::with_capacity(count + 2);
        if !full {
            for (cap, end) in [false, true].into_iter().enumerate() {
                let anchor = if end {
                    PlaneFrame {
                        origin: turn.point(frame.origin),
                        u: turn.vector(frame.u),
                        v: turn.vector(frame.v),
                    }
                } else {
                    frame
                };
                // Material leaves the start cap along the sweep and reaches
                // the end cap along it.
                let sweep = scale(
                    cross(a, meridian.radial_at(if end { angle } else { 0. })),
                    angle.signum(),
                );
                let along = (dot(anchor.normal(), sweep) > 0.) == end;
                let support = ReferenceFrame {
                    origin_mm: anchor.origin,
                    u: anchor.u,
                    v: if along {
                        anchor.v
                    } else {
                        anchor.v.map(|x| -x)
                    },
                };
                let loops = rings
                    .iter()
                    .map(|ring| {
                        let uses = (0..ring.len)
                            .map(|k| (meridians[cap][ring.edge(k)], ring.forward == end));
                        Loop {
                            uses: if end {
                                uses.collect()
                            } else {
                                uses.rev().collect()
                            },
                            vertex: None,
                        }
                    })
                    .collect();
                faces.push((
                    0,
                    brep::Face {
                        surface: Surface::Plane(Frame {
                            origin: support.origin_mm,
                            x: support.u,
                            y: support.v,
                            z: support.normal(),
                        }),
                        same_sense: true,
                        loops,
                    },
                ));
                names.push(FaceName {
                    origins: vec![Origin {
                        source: feature,
                        role: if end { "cap.end" } else { "cap.start" }.into(),
                    }],
                    support: Some(support),
                    anchor: Some(ReferenceFrame {
                        origin_mm: anchor.origin,
                        u: anchor.u,
                        v: anchor.v,
                    }),
                });
            }
        }
        let mut sides = vec![None; count];
        for (ring_index, (ring, ring_rails)) in rings.iter().zip(&rails).enumerate() {
            for k in 0..ring.len {
                let edge = ring.edge(k);
                if along_axis[edge] {
                    continue;
                }
                let next = (k + 1) % ring.len;
                // A full turn closes each hole into a void shell, described
                // as the solid it encloses and then reversed as a whole.
                let void = full && !ring.outer;
                let [after, before] = [ring_rails[next], ring_rails[k]];
                let loops = if seamless[edge] {
                    [(after, toward), (before, !toward)]
                        .into_iter()
                        .filter_map(|(rail, forward)| {
                            rail.map(|rail| Loop {
                                uses: vec![(rail, forward)],
                                vertex: None,
                            })
                        })
                        .collect()
                } else {
                    let uses = [
                        Some((meridians[0][edge], ring.forward)),
                        after.map(|rail| (rail, toward)),
                        Some((meridians[1][edge], !ring.forward)),
                        before.map(|rail| (rail, !toward)),
                    ];
                    // With the meridian tangent t, the outward normal n and
                    // the sweep s, t x s = -n for a positive turn: reverse the
                    // loop there so that it runs counter-clockwise about n
                    // even where winding loops cannot show their direction.
                    let mut uses: Vec<_> = uses.into_iter().flatten().collect();
                    if toward != void {
                        uses.reverse();
                        uses.iter_mut().for_each(|used| used.1 = !used.1);
                    }
                    vec![Loop { uses, vertex: None }]
                };
                let mut surface = surfaces[edge].clone();
                // Outward from the material at the middle of the wall.
                let (point, tangent) = middle(&start[edge], tolerance)?;
                let tangent = if ring.forward {
                    tangent
                } else {
                    tangent.map(|x| -x)
                };
                let half = rotation(axis.origin_mm, a, 0.5 * angle);
                let [dr, dz] = [dot(tangent, radial), dot(tangent, a)];
                let outward = sub(scale(meridian.radial_at(0.5 * angle), dz), scale(a, dr));
                let outward = if void { outward.map(|x| -x) } else { outward };
                let sample = half.point(point);
                let uv = surface.invert(sample, None);
                let (on, [su, sv]) = surface.evaluate(uv);
                if norm(sub(on, sample)) > tolerance.max(1e-9 * norm(sample)) {
                    return Err(error(
                        REVOLVE,
                        GeometryErrorKind::InvalidResult,
                        "A revolved wall does not contain its own meridian",
                    ));
                }
                let mut same_sense = dot(cross(su, sv), outward) > 0.;
                if full && !same_sense && !surface.is_plane() {
                    // A closed wall facing away from its natural normal
                    // cannot show its sense through winding loops. Revolve
                    // the meridian in the direction whose own normal faces
                    // out of the material instead.
                    for path in [start[edge].clone(), start[edge].reversed()] {
                        let candidate = Surface::Revolution {
                            curve: path.curve(tolerance)?,
                            origin: axis.origin_mm,
                            axis: a,
                        };
                        let uv = candidate.invert(sample, None);
                        let (_, [su, sv]) = candidate.evaluate(uv);
                        if dot(cross(su, sv), outward) > 0. {
                            surface = candidate;
                            same_sense = true;
                            break;
                        }
                    }
                }
                let mut support = None;
                let mut anchor = None;
                if let Surface::Plane(plane) = &mut surface {
                    // Planar annuli face outward and keep authored frames.
                    if !same_sense {
                        plane.y = plane.y.map(|x| -x);
                        plane.z = plane.z.map(|x| -x);
                        same_sense = true;
                    }
                    support = Some(ReferenceFrame {
                        origin_mm: plane.origin,
                        u: plane.x,
                        v: plane.y,
                    });
                    let origin_mm = start[edge].ends(tolerance)?[0];
                    let mut u = sub(frame.u, scale(a, dot(frame.u, a)));
                    if norm(u) < 1e-6 {
                        u = sub(frame.v, scale(a, dot(frame.v, a)));
                    }
                    let u = unit(u)?;
                    anchor = Some(ReferenceFrame {
                        origin_mm,
                        u,
                        v: cross(a, u),
                    });
                }
                sides[edge] = Some((
                    ring_index,
                    brep::Face {
                        surface,
                        same_sense,
                        loops,
                    },
                    FaceName {
                        origins: vec![
                            Origin {
                                source: feature,
                                role: format!(
                                    "side/{}/{}",
                                    origins[edge].source, origins[edge].role
                                ),
                            },
                            origins[edge].clone(),
                        ],
                        support,
                        anchor,
                    },
                ));
            }
        }
        // Edges on the axis sweep no wall.
        for (ring, face, name) in sides.into_iter().flatten() {
            faces.push((ring, face));
            names.push(name);
        }
        // A partial turn is one shell; a full turn closes every hole into a
        // void of its own.
        let shell_count = if full { rings.len() } else { 1 };
        let mut shells = vec![Vec::new(); shell_count];
        for (index, (ring, face)) in faces.into_iter().enumerate() {
            shells[if full { ring } else { 0 }].push(index);
            boundary.faces.push(face);
        }
        boundary.shells = shells
            .into_iter()
            .enumerate()
            .map(|(index, faces)| Shell {
                faces,
                reversed: index > 0,
            })
            .collect();
        let solid = brep::Solid::new(boundary, tolerance, None).map_err(failure)?;
        Solid::from_authored_brep(solid, feature, AuthoredNaming { faces: names }, policy)
    }
}
