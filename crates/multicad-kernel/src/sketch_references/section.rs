//! Exact plane sections of finite Rust faces for sketch references.
//!
//! Planar faces are cut along the line shared with the sketch plane. The
//! cut is split at the certified contacts of the face's own boundary edges
//! and each piece is kept only when a transverse in-face ray proves its
//! midpoint inside the face. Cylinder side faces are cut analytically, and
//! only where the body itself is a right circular cylinder. Every other
//! face family is a typed refusal, never a sampled guess.
use super::{
    SketchReferenceCatalog, SketchReferenceGeometry, SketchReferenceMethod, SketchReferenceResult,
    error, rust_element_reference,
};
use crate::{
    BSpline2d, CurvedEdge2d, GeometryError, GeometryErrorKind, GeometryPolicy, PlaneFrame, Solid,
};
use multicad_core::references::{
    AnalyticGeometry, ReferenceFrame, TopologyElement, TopologyKey, TopologyKind, TopologySnapshot,
};
use multicad_precision::PortableMath as _;

type V3 = [f64; 3];

fn dot(a: V3, b: V3) -> f64 {
    a.into_iter().zip(b).map(|(a, b)| a * b).sum()
}
fn cross(a: V3, b: V3) -> V3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn sub(a: V3, b: V3) -> V3 {
    std::array::from_fn(|i| a[i] - b[i])
}
fn add(a: V3, b: V3) -> V3 {
    std::array::from_fn(|i| a[i] + b[i])
}
fn scale(a: V3, s: f64) -> V3 {
    a.map(|v| v * s)
}
fn norm(a: V3) -> f64 {
    dot(a, a).sqrt()
}
fn unit(a: V3) -> Option<V3> {
    let n = norm(a);
    (n.is_finite() && n > 1e-12).then(|| scale(a, 1. / n))
}
fn unsupported(message: &str) -> GeometryError {
    GeometryError {
        operation: "sketch reference",
        kind: GeometryErrorKind::Unsupported,
        message: message.into(),
    }
}
/// Angular tolerance for parallel and perpendicular unit directions.
const ANGLE: f64 = 1e-12;

struct Target {
    frame: PlaneFrame,
    normal: V3,
}
impl Target {
    fn new(frame: PlaneFrame) -> Result<Self, GeometryError> {
        let normal = unit(frame.normal()).ok_or_else(|| error("Invalid sketch plane"))?;
        Ok(Self { frame, normal })
    }
    fn local(&self, p: V3) -> [f64; 2] {
        let d = sub(p, self.frame.origin);
        [dot(d, self.frame.u), dot(d, self.frame.v)]
    }
    fn local_direction(&self, d: V3) -> [f64; 2] {
        [dot(d, self.frame.u), dot(d, self.frame.v)]
    }
    fn world(&self, p: [f64; 2]) -> V3 {
        add(
            self.frame.origin,
            add(scale(self.frame.u, p[0]), scale(self.frame.v, p[1])),
        )
    }
}

fn face_element(
    topology: &TopologySnapshot,
    index: u32,
) -> Result<&TopologyElement, GeometryError> {
    let key = TopologyKey {
        kind: TopologyKind::Face,
        ordinal: index,
    };
    topology
        .elements
        .iter()
        .find(|element| element.key == key)
        .ok_or_else(|| error("Selected Rust reference topology is missing"))
}

/// Exact geometric incidence: the analytic edges lying in a face plane.
/// Snapshot origins can be empty, so they cannot identify incidence. A
/// second coplanar face or an edge without analytic geometry would make
/// this incidence ambiguous and is refused.
fn planar_boundary<'a>(
    topology: &'a TopologySnapshot,
    face: &TopologyElement,
    origin: V3,
    normal: V3,
    tolerance: f64,
) -> Result<Vec<&'a TopologyElement>, GeometryError> {
    let on_plane = |p: V3| dot(sub(p, origin), normal).abs() <= tolerance;
    let in_plane = |frame: &ReferenceFrame| {
        unit(cross(frame.u, frame.v)).is_some_and(|n| norm(cross(n, normal)) <= ANGLE)
            && on_plane(frame.origin_mm)
    };
    let mut edges = vec![];
    for element in &topology.elements {
        match (element.key.kind, element.geometry) {
            (TopologyKind::Face, Some(AnalyticGeometry::Plane { frame }))
                if element.key != face.key && in_plane(&frame) =>
            {
                return Err(unsupported(
                    "Another coplanar Rust face makes the face boundary ambiguous",
                ));
            }
            (TopologyKind::Edge, None) => {
                return Err(unsupported(
                    "Rust face sections require analytic boundary edges",
                ));
            }
            (TopologyKind::Edge, Some(AnalyticGeometry::Line { a_mm, b_mm }))
                if on_plane(a_mm) && on_plane(b_mm) =>
            {
                edges.push(element);
            }
            (
                TopologyKind::Edge,
                Some(
                    AnalyticGeometry::CircularArc { frame, .. }
                    | AnalyticGeometry::EllipticArc { frame, .. },
                ),
            ) if in_plane(&frame) => edges.push(element),
            (TopologyKind::Edge, Some(AnalyticGeometry::InvoluteArc { .. })) => {
                return Err(unsupported(
                    "Rust face sections do not yet bound faces by involute edges",
                ));
            }
            _ => {}
        }
    }
    Ok(edges)
}

/// Finite end vertices of open boundary edges.
fn vertices(edges: &[&TopologyElement]) -> Vec<V3> {
    let mut result = vec![];
    for edge in edges {
        match edge.geometry {
            Some(AnalyticGeometry::Line { a_mm, b_mm }) => result.extend([a_mm, b_mm]),
            Some(
                AnalyticGeometry::CircularArc {
                    frame,
                    radius_mm,
                    start_rad,
                    sweep_rad,
                }
                | AnalyticGeometry::EllipticArc {
                    frame,
                    semiaxes_mm: [radius_mm, _],
                    start_rad,
                    sweep_rad,
                },
            ) if sweep_rad.abs() < std::f64::consts::TAU - 1e-9 => {
                let semiaxes = match edge.geometry {
                    Some(AnalyticGeometry::EllipticArc { semiaxes_mm, .. }) => semiaxes_mm,
                    _ => [radius_mm; 2],
                };
                for angle in [start_rad, start_rad + sweep_rad] {
                    let (sin, cos) = angle.portable_sin_cos();
                    result.push(add(
                        frame.origin_mm,
                        add(
                            scale(frame.u, semiaxes[0] * cos),
                            scale(frame.v, semiaxes[1] * sin),
                        ),
                    ));
                }
            }
            _ => {}
        }
    }
    result
}

struct Contact {
    point: [f64; 2],
    support: Vec<TopologyKey>,
}

/// Transverse contacts of all boundary edges with one plane. `None` means
/// the plane runs along a boundary edge or touches one tangentially.
fn contacts(
    edges: &[&TopologyElement],
    frame: PlaneFrame,
    policy: GeometryPolicy,
    allow_tangent: bool,
) -> Result<Option<Vec<Contact>>, GeometryError> {
    let mut result = vec![];
    for edge in edges {
        let catalog =
            rust_element_reference(edge, frame, SketchReferenceMethod::Intersect, policy)?;
        if catalog.coincident {
            return Ok(None);
        }
        for hit in catalog.results {
            if hit.tangent && !allow_tangent {
                return Ok(None);
            }
            if hit.tangent {
                continue;
            }
            let SketchReferenceGeometry::Point(point) = hit.geometry else {
                return Err(error("Rust edge plane contact is not a point"));
            };
            result.push(Contact {
                point,
                support: vec![edge.key],
            });
        }
    }
    Ok(Some(result))
}

/// In-face ray parity at `point`, using rays that avoid vertices and tangents.
fn inside(
    edges: &[&TopologyElement],
    corners: &[V3],
    point: V3,
    face: (V3, V3, V3),
    policy: GeometryPolicy,
) -> Result<bool, GeometryError> {
    let tolerance = policy.linear_tolerance_mm;
    let (u, v, normal) = face;
    for angle in [0.3711_f64, 1.2345, 2.0187, 2.9031, 3.7777, 4.6093, 5.4419] {
        let (sin, cos) = angle.portable_sin_cos();
        let Some(direction) = unit(add(scale(u, cos), scale(v, sin))) else {
            continue;
        };
        // The cutting plane contains the ray and the face normal.
        let frame = PlaneFrame {
            origin: point,
            u: direction,
            v: normal,
        };
        // A probe whose contacts cannot be certified proves nothing; the
        // next independent direction is tried instead.
        let Ok(Some(hits)) = contacts(edges, frame, policy, false) else {
            continue;
        };
        let ray = |p: [f64; 2]| p[0];
        if hits
            .iter()
            .any(|hit| ray(hit.point).abs() <= tolerance * 2.)
            || hits.iter().any(|hit| {
                let world = add(point, scale(direction, hit.point[0]));
                corners
                    .iter()
                    .any(|corner| norm(sub(*corner, world)) <= tolerance * 2.)
            })
        {
            continue;
        }
        let crossings = hits.iter().filter(|hit| ray(hit.point) > 0.).count();
        return Ok(crossings % 2 == 1);
    }
    Err(unsupported(
        "No transverse in-face ray could classify the Rust face section",
    ))
}

fn planar(
    topology: &TopologySnapshot,
    face: &TopologyElement,
    plane: ReferenceFrame,
    target: &Target,
    method: SketchReferenceMethod,
    policy: GeometryPolicy,
) -> Result<SketchReferenceCatalog, GeometryError> {
    let tolerance = policy.linear_tolerance_mm;
    let (Some(u), Some(v)) = (unit(plane.u), unit(plane.v)) else {
        return Err(error("Invalid Rust face plane"));
    };
    let normal = unit(cross(u, v)).ok_or_else(|| error("Invalid Rust face plane"))?;
    let edges = planar_boundary(topology, face, plane.origin_mm, normal, tolerance)?;
    if edges.is_empty() {
        return Err(error("Rust face has no boundary edges"));
    }
    let parallel = norm(cross(normal, target.normal)) <= ANGLE;
    let offset = dot(sub(plane.origin_mm, target.frame.origin), target.normal);
    let coplanar = parallel && offset.abs() <= tolerance;
    match method {
        SketchReferenceMethod::Project => Err(unsupported(
            "Rust face references support coplanar boundaries and plane intersections",
        )),
        SketchReferenceMethod::Coplanar => {
            if !coplanar {
                return Err(error(
                    "Selected finite Rust face is not coplanar with the sketch plane",
                ));
            }
            let mut results = vec![];
            for edge in &edges {
                for mut result in
                    rust_element_reference(edge, target.frame, method, policy)?.results
                {
                    result.support.push(face.key);
                    result.support.sort();
                    result.support.dedup();
                    results.push(result);
                }
            }
            Ok(SketchReferenceCatalog {
                results,
                coincident: false,
            })
        }
        SketchReferenceMethod::Intersect => {
            if parallel {
                return Ok(SketchReferenceCatalog {
                    results: vec![],
                    coincident: coplanar,
                });
            }
            if let [rim] = edges.as_slice()
                && let Some(AnalyticGeometry::CircularArc {
                    frame: circle,
                    radius_mm,
                    sweep_rad,
                    ..
                }) = rim.geometry
                && (sweep_rad.abs() - std::f64::consts::TAU).abs() <= 1e-9
            {
                return disk(
                    face.key,
                    rim.key,
                    circle.origin_mm,
                    radius_mm,
                    normal,
                    target,
                    policy,
                );
            }
            let Some(mut hits) = contacts(&edges, target.frame, policy, true)? else {
                return Err(unsupported(
                    "The sketch plane runs along a boundary edge of the selected Rust face",
                ));
            };
            let direction = target.local_direction(cross(target.normal, normal));
            let along = |p: [f64; 2]| p[0] * direction[0] + p[1] * direction[1];
            hits.sort_by(|a, b| along(a.point).total_cmp(&along(b.point)));
            let mut merged: Vec<Contact> = vec![];
            for hit in hits {
                if let Some(last) = merged.last_mut()
                    && (last.point[0] - hit.point[0]).portable_hypot(last.point[1] - hit.point[1])
                        <= tolerance * 2.
                {
                    last.support.extend(hit.support);
                    continue;
                }
                merged.push(hit);
            }
            let corners = vertices(&edges);
            let mut pieces: Vec<(usize, usize)> = vec![];
            for i in 1..merged.len() {
                let [a, b] = [merged[i - 1].point, merged[i].point];
                let middle = target.world([(a[0] + b[0]) * 0.5, (a[1] + b[1]) * 0.5]);
                if inside(&edges, &corners, middle, (u, v, normal), policy)? {
                    match pieces.last_mut() {
                        Some(piece) if piece.1 == i - 1 => piece.1 = i,
                        _ => pieces.push((i - 1, i)),
                    }
                }
            }
            let mut results = vec![];
            for (first, last) in pieces {
                let curve = CurvedEdge2d::Line {
                    a: merged[first].point,
                    b: merged[last].point,
                };
                curve.data(policy, 0)?;
                let mut support = vec![face.key];
                support.extend(&merged[first].support);
                support.extend(&merged[last].support);
                support.sort();
                support.dedup();
                results.push(SketchReferenceResult {
                    geometry: SketchReferenceGeometry::Curve(curve),
                    support,
                    source_parameter: None,
                    source_interval: None,
                    tangent: false,
                    crossing: None,
                });
            }
            Ok(SketchReferenceCatalog {
                results,
                coincident: false,
            })
        }
    }
}

/// A planar face bounded by one full circle is a disk; its cut by a
/// nonparallel plane is an exact chord, a tangent touch or nothing.
fn disk(
    face: TopologyKey,
    rim: TopologyKey,
    center: V3,
    radius: f64,
    normal: V3,
    target: &Target,
    policy: GeometryPolicy,
) -> Result<SketchReferenceCatalog, GeometryError> {
    let tolerance = policy.linear_tolerance_mm;
    let along = unit(cross(target.normal, normal)).ok_or_else(|| error("Invalid face plane"))?;
    // In-face direction from the center toward the cut line.
    let toward = unit(cross(normal, along)).ok_or_else(|| error("Invalid face plane"))?;
    let slope = dot(toward, target.normal);
    let distance = -dot(sub(center, target.frame.origin), target.normal) / slope;
    let mut results = vec![];
    if distance.abs() < radius - tolerance {
        let foot = add(center, scale(toward, distance));
        let half = ((radius - distance.abs()) * (radius + distance.abs())).sqrt();
        let curve = CurvedEdge2d::Line {
            a: target.local(sub(foot, scale(along, half))),
            b: target.local(add(foot, scale(along, half))),
        };
        let mut result = curve_result(face, curve, false, policy)?;
        result.support.push(rim);
        result.support.sort();
        results.push(result);
    }
    Ok(SketchReferenceCatalog {
        results,
        coincident: false,
    })
}

/// The two full rims of a right circular cylinder side face.
struct Cylinder {
    base: V3,
    axis: V3,
    height: f64,
    radius: f64,
    x: V3,
    y: V3,
}

/// An axial cylinder body has exactly one curved face, bounded by its two
/// full rims; the body type proves the surface, the rims give its placement.
fn cylinder(topology: &TopologySnapshot, tolerance: f64) -> Option<Cylinder> {
    let rims: Vec<_> = topology
        .elements
        .iter()
        .filter(|edge| edge.key.kind == TopologyKind::Edge)
        .filter_map(|edge| match edge.geometry {
            Some(AnalyticGeometry::CircularArc {
                frame,
                radius_mm,
                sweep_rad,
                ..
            }) if (sweep_rad.abs() - std::f64::consts::TAU).abs() <= 1e-9 => {
                Some((frame, radius_mm))
            }
            _ => None,
        })
        .collect();
    let [(bottom, radius), (top, other)] = rims.as_slice() else {
        return None;
    };
    let delta = sub(top.origin_mm, bottom.origin_mm);
    let height = norm(delta);
    let axis = unit(delta)?;
    let x = unit(bottom.u)?;
    let y = unit(cross(axis, x))?;
    let normal = |f: &ReferenceFrame| unit(cross(f.u, f.v));
    ((radius - other).abs() <= tolerance
        && height > tolerance
        && [bottom, top]
            .into_iter()
            .all(|f| normal(f).is_some_and(|n| norm(cross(n, axis)) <= ANGLE))
        && dot(x, axis).abs() <= ANGLE)
        .then_some(Cylinder {
            base: bottom.origin_mm,
            axis,
            height,
            radius: *radius,
            x,
            y,
        })
}

fn curve_result(
    face: TopologyKey,
    curve: CurvedEdge2d,
    tangent: bool,
    policy: GeometryPolicy,
) -> Result<SketchReferenceResult, GeometryError> {
    curve.data(policy, 0)?;
    Ok(SketchReferenceResult {
        geometry: SketchReferenceGeometry::Curve(curve),
        support: vec![face],
        source_parameter: None,
        source_interval: None,
        tangent,
        crossing: None,
    })
}

fn cylindrical(
    c: Cylinder,
    face: TopologyKey,
    target: &Target,
    method: SketchReferenceMethod,
    policy: GeometryPolicy,
) -> Result<SketchReferenceCatalog, GeometryError> {
    let tolerance = policy.linear_tolerance_mm;
    if method != SketchReferenceMethod::Intersect {
        return Err(unsupported(
            "Rust cylinder face references support plane intersections",
        ));
    }
    let n = target.normal;
    let slope = dot(n, c.axis);
    let offset = dot(sub(c.base, target.frame.origin), n);
    let mut results = vec![];
    if slope.abs() <= ANGLE {
        // Axis parallel to the plane: generator lines.
        let foot = sub(c.base, scale(n, offset));
        let distance = offset.abs();
        let line = |start: V3| CurvedEdge2d::Line {
            a: target.local(start),
            b: target.local(add(start, scale(c.axis, c.height))),
        };
        if (distance - c.radius).abs() <= tolerance {
            results.push(curve_result(face, line(foot), true, policy)?);
        } else if distance < c.radius {
            let half = ((c.radius - distance) * (c.radius + distance)).sqrt();
            let across = unit(cross(n, c.axis)).ok_or_else(|| error("Invalid cylinder axis"))?;
            for sign in [-1., 1.] {
                results.push(curve_result(
                    face,
                    line(add(foot, scale(across, sign * half))),
                    false,
                    policy,
                )?);
            }
        }
    } else if norm(cross(n, c.axis)) <= ANGLE {
        // Plane normal to the axis: one full circle inside the face span.
        let along = -offset / slope;
        if (-tolerance..=c.height + tolerance).contains(&along) {
            results.push(curve_result(
                face,
                CurvedEdge2d::Circle {
                    center: target.local(add(c.base, scale(c.axis, along))),
                    radius_mm: c.radius,
                },
                false,
                policy,
            )?);
        }
    } else {
        // Oblique plane: P(t) = base + r(cos t x + sin t y) + h(t) axis.
        let center_height = -offset / slope;
        let [hx, hy] = [c.x, c.y].map(|d| -c.radius * dot(d, n) / slope);
        let reach = hx.portable_hypot(hy);
        let [low, high] = [center_height - reach, center_height + reach];
        if high < -tolerance || low > c.height + tolerance {
            return Ok(SketchReferenceCatalog {
                results,
                coincident: false,
            });
        }
        if low < -tolerance || high > c.height + tolerance {
            return Err(unsupported(
                "Rust cylinder references cannot yet trim an oblique section at the face rims",
            ));
        }
        let center = target.local(add(c.base, scale(c.axis, center_height)));
        let [b, d] = [(c.x, hx), (c.y, hy)].map(|(radial, h)| {
            target.local_direction(add(scale(radial, c.radius), scale(c.axis, h)))
        });
        // Exact affine image of the four-span rational quadratic circle.
        let w = std::f64::consts::FRAC_1_SQRT_2;
        let unit_poles = [
            [1., 0.],
            [1., 1.],
            [0., 1.],
            [-1., 1.],
            [-1., 0.],
            [-1., -1.],
            [0., -1.],
            [1., -1.],
            [1., 0.],
        ];
        let poles = unit_poles
            .map(|[x, y]| std::array::from_fn(|i| center[i] + x * b[i] + y * d[i]))
            .to_vec();
        let ellipse = CurvedEdge2d::BSpline(BSpline2d {
            degree: 2,
            poles,
            knots: vec![0., 0.25, 0.5, 0.75, 1.],
            multiplicities: vec![3, 2, 2, 2, 3],
            weights: Some(vec![1., w, 1., w, 1., w, 1., w, 1.]),
        });
        results.push(curve_result(face, ellipse, false, policy)?);
    }
    Ok(SketchReferenceCatalog {
        results,
        coincident: false,
    })
}

pub(super) fn face(
    solid: &Solid,
    topology: &TopologySnapshot,
    index: u32,
    frame: PlaneFrame,
    method: SketchReferenceMethod,
    policy: GeometryPolicy,
) -> Result<SketchReferenceCatalog, GeometryError> {
    let target = Target::new(frame)?;
    let element = face_element(topology, index)?;
    match element.geometry {
        Some(AnalyticGeometry::Plane { frame: plane }) => {
            planar(topology, element, plane, &target, method, policy)
        }
        None if solid.axial_cylinder_frame().is_some() => {
            let c = cylinder(topology, policy.linear_tolerance_mm).ok_or_else(|| {
                unsupported("Rust cylinder face rims do not form one right circular cylinder")
            })?;
            cylindrical(c, element.key, &target, method, policy)
        }
        _ => Err(unsupported(
            "Rust face references support planar faces and right circular cylinder faces",
        )),
    }
}

pub(super) fn whole(
    solid: &Solid,
    topology: &TopologySnapshot,
    frame: PlaneFrame,
    method: SketchReferenceMethod,
    policy: GeometryPolicy,
) -> Result<SketchReferenceCatalog, GeometryError> {
    if method != SketchReferenceMethod::Intersect {
        return Err(error(
            "A whole-body sketch reference must intersect the sketch plane",
        ));
    }
    let mut catalog = SketchReferenceCatalog {
        results: vec![],
        coincident: false,
    };
    for element in &topology.elements {
        if element.key.kind != TopologyKind::Face {
            continue;
        }
        let section = face(solid, topology, element.key.ordinal, frame, method, policy)?;
        catalog.coincident |= section.coincident;
        catalog.results.extend(section.results);
        if catalog.results.len() > 256 {
            return Err(error("Rust reference result budget exceeded"));
        }
    }
    Ok(catalog)
}
