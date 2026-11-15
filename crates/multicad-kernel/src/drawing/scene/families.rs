//! Scene members of the Rust body families.
use super::bore;
use super::prism::{Band, Prism};
use super::profile::Rings;
use super::*;
use multicad_core::references::{AnalyticGeometry, TopologyKind};
use multicad_rust_kernel::{normal_sweep, spur_gear::Edge};
use std::f64::consts::TAU;

/// The exact scene member of one leaf solid in a view, or `None` when its
/// family has no member yet. Each family contributes one constructor here.
pub(in crate::drawing) fn member(
    solid: &Solid,
    view: &View,
) -> Option<Result<Box<dyn Member>, GeometryError>> {
    // HOOK twisted: helical and herringbone gears (`scene/twisted.rs`), e.g.
    //     .or_else(|| twisted::member(solid, view))
    // HOOK bevel: spherical bevel gears (`scene/bevel.rs`).
    // HOOK worm: ZA worms and worm wheels (`scene/worm.rs`).
    gear(solid, view)
        .or_else(|| extrusion(solid, view))
        .or_else(|| bore::Bore::new(solid, view).map(|bore| Ok(Box::new(bore) as Box<dyn Member>)))
        .or_else(|| super::bevel::member(solid, view))
        .or_else(|| super::helix::member(solid, view))
        .or_else(|| super::pipe::member(solid, view))
}

/// Draw whole solids through the scene engine when every leaf has a member.
pub(in crate::drawing) fn project_solids(
    solids: &[DrawingSolid<'_>],
    options: DrawingProjection,
) -> Option<Result<Vec<ProjectedVector>, GeometryError>> {
    let view = View::new(options.frame);
    let mut policy = GeometryPolicy::default();
    let mut members = Vec::new();
    for (source, solid) in projection_members(solids) {
        let built = member(solid, &view)?;
        let leaf = solid.geometry_policy();
        policy.linear_tolerance_mm = policy.linear_tolerance_mm.max(leaf.linear_tolerance_mm);
        policy.minimum_feature_mm = policy.minimum_feature_mm.min(leaf.minimum_feature_mm);
        members.push(built.map(|member| SceneMember {
            source,
            member,
            hatch: solids[source].hatch,
        }));
    }
    Some(
        members
            .into_iter()
            .collect::<Result<Vec<_>, _>>()
            .and_then(|members| project(&members, options, policy)),
    )
}

fn placed(
    view: &View,
    frame: PlaneFrame,
    bands: Vec<Band>,
    solid: &Solid,
) -> Result<Box<dyn Member>, GeometryError> {
    Ok(Box::new(Prism::new(
        view.point(frame.origin),
        [frame.u, frame.v, frame.normal()].map(|axis| view.vector(axis)),
        bands,
        solid.geometry_policy().linear_tolerance_mm,
    )?))
}

fn gear(solid: &Solid, view: &View) -> Option<Result<Box<dyn Member>, GeometryError>> {
    if let Some((body, frame)) = solid.spur_gear_body() {
        let profile = body.profile();
        let bore: Option<Vec<Edge>> = profile.bore().map(|bore| vec![Edge::Planar(bore)]);
        let mut rings: Rings = vec![profile.outline().iter().map(|edge| edge.edge).collect()];
        rings.extend(bore.clone());
        let mut bands = vec![Band {
            rings,
            axis: [0., profile.height_mm()],
        }];
        for extension in profile.hub_extensions().into_iter().flatten() {
            let mut rings: Rings = vec![vec![Edge::Planar(normal_sweep::Edge::Circle {
                center: [0.; 2],
                radius_mm: extension.radius_mm,
                clockwise: false,
            })]];
            rings.extend(bore.clone());
            let [a, b] = [extension.gear_z_mm, extension.end_z_mm];
            bands.push(Band {
                rings,
                axis: [a.min(b), a.max(b)],
            });
        }
        return Some(placed(view, frame, bands, solid));
    }
    let (body, frame) = solid.internal_ring_body()?;
    let ring = &body.ring;
    let bands = vec![Band {
        rings: vec![
            vec![Edge::Planar(ring.rim())],
            ring.cavity().iter().map(|edge| edge.edge).collect(),
        ],
        axis: [0., ring.height_mm()],
    }];
    Some(placed(view, frame, bands, solid))
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Right prisms recovered from the exact topology of the extruded families:
/// every edge lies in one of two parallel cap planes or joins them along the
/// extrusion, and both caps carry the same profile.
fn extrusion(solid: &Solid, view: &View) -> Option<Result<Box<dyn Member>, GeometryError>> {
    if solid.rust_prism_frame().is_none()
        && solid.axial_cylinder_frame().is_none()
        && solid.annulus_body().is_none()
        && solid.rounded_slot_body().is_none()
        && solid.orthogonal_body().is_none()
    {
        return None;
    }
    let tolerance = solid.geometry_policy().linear_tolerance_mm;
    let attempt = || -> Result<Option<Box<dyn Member>>, GeometryError> {
        let Some(topology) = solid.rust_topology() else {
            return Ok(None);
        };
        let edges: Vec<AnalyticGeometry> = match topology?
            .elements
            .into_iter()
            .filter(|element| element.key.kind == TopologyKind::Edge)
            .map(|element| element.geometry)
            .collect::<Option<Vec<_>>>()
        {
            Some(edges) => edges,
            None => return Ok(None),
        };
        let mut normals: Vec<[f64; 3]> = Vec::new();
        for face in solid.faces()? {
            if let Some((_, normal)) = face.plane
                && !normals
                    .iter()
                    .any(|known| dot(*known, normal).abs() > 1. - 1e-12)
            {
                normals.push(normal);
            }
        }
        for axis in normals {
            if let Some(rings_and_levels) = caps(&edges, axis, tolerance) {
                let (rings, levels, [x, y]) = rings_and_levels;
                return Ok(Some(Box::new(Prism::new(
                    view.point([0.; 3]),
                    [x, y, axis].map(|direction| view.vector(direction)),
                    vec![Band {
                        rings,
                        axis: levels,
                    }],
                    tolerance,
                )?)));
            }
        }
        Ok(None)
    };
    attempt().transpose()
}

type Caps = (Rings, [f64; 2], [[f64; 3]; 2]);

fn caps(edges: &[AnalyticGeometry], axis: [f64; 3], tolerance: f64) -> Option<Caps> {
    let seed = if axis[0].abs() < 0.9 {
        [1., 0., 0.]
    } else {
        [0., 1., 0.]
    };
    let cross = |a: [f64; 3], b: [f64; 3]| -> [f64; 3] {
        std::array::from_fn(|i| a[(i + 1) % 3] * b[(i + 2) % 3] - a[(i + 2) % 3] * b[(i + 1) % 3])
    };
    let unit = |v: [f64; 3]| {
        let length = dot(v, v).sqrt();
        v.map(|c| c / length)
    };
    let y = unit(cross(axis, seed));
    let x = cross(y, axis);
    let flat = |p: [f64; 3]| [dot(p, x), dot(p, y)];
    let mut levels: Vec<(f64, Vec<normal_sweep::Edge>)> = Vec::new();
    for edge in edges {
        let (level, planar) = match *edge {
            AnalyticGeometry::Line { a_mm, b_mm } => {
                let [la, lb] = [dot(a_mm, axis), dot(b_mm, axis)];
                let [pa, pb] = [flat(a_mm), flat(b_mm)];
                if (la - lb).abs() <= tolerance {
                    (la, normal_sweep::Edge::Line { a: pa, b: pb })
                } else if (pa[0] - pb[0]).hypot(pa[1] - pb[1]) <= tolerance {
                    continue;
                } else {
                    return None;
                }
            }
            AnalyticGeometry::CircularArc {
                frame,
                radius_mm,
                start_rad,
                sweep_rad,
            } => {
                let facing = dot(frame.normal(), axis);
                if facing.abs() < 1. - 1e-12 {
                    return None;
                }
                let at = |angle: f64| {
                    let (sin, cos) = angle.sin_cos();
                    flat(std::array::from_fn(|i| {
                        frame.origin_mm[i] + radius_mm * (cos * frame.u[i] + sin * frame.v[i])
                    }))
                };
                let center = flat(frame.origin_mm);
                let clockwise = facing * sweep_rad < 0.;
                let planar = if sweep_rad.abs() >= TAU - 1e-12 {
                    normal_sweep::Edge::Circle {
                        center,
                        radius_mm,
                        clockwise,
                    }
                } else {
                    normal_sweep::Edge::Arc {
                        center,
                        a: at(start_rad),
                        b: at(start_rad + sweep_rad),
                        clockwise,
                    }
                };
                (dot(frame.origin_mm, axis), planar)
            }
            _ => return None,
        };
        match levels
            .iter_mut()
            .find(|(known, _)| (*known - level).abs() <= tolerance)
        {
            Some((_, list)) => list.push(planar),
            None => levels.push((level, vec![planar])),
        }
    }
    levels.sort_by(|a, b| a.0.total_cmp(&b.0));
    let [(low, base), (high, lid)] = <[_; 2]>::try_from(levels).ok()?;
    // Both caps must carry the same curves for the walls to be a prism's.
    let near = |a: [f64; 2], b: [f64; 2]| (a[0] - b[0]).hypot(a[1] - b[1]) <= 4. * tolerance;
    let same = |a: &normal_sweep::Edge, b: &normal_sweep::Edge| {
        std::mem::discriminant(a) == std::mem::discriminant(b)
            && a.radius_mm()
                .zip(b.radius_mm())
                .is_none_or(|(r, s)| (r - s).abs() <= tolerance)
            && near(a.point(0.5), b.point(0.5))
            && ((near(a.start(), b.start()) && near(a.end(), b.end()))
                || (near(a.start(), b.end()) && near(a.end(), b.start())))
    };
    if base.len() != lid.len()
        || !base
            .iter()
            .all(|edge| lid.iter().any(|twin| same(edge, twin)))
    {
        return None;
    }
    // Chain the cap edges into closed rings.
    let mut pool = base;
    let mut rings: Rings = Vec::new();
    while let Some(first) = pool.pop() {
        let mut ring = vec![first];
        while !near(ring[0].start(), ring.last()?.end())
            || (ring.len() == 1 && first.start() != first.end() && false)
        {
            let end = ring.last()?.end();
            let next = pool
                .iter()
                .position(|edge| near(edge.start(), end) || near(edge.end(), end))?;
            let edge = pool.swap_remove(next);
            ring.push(if near(edge.start(), end) {
                edge
            } else {
                edge.reversed()
            });
        }
        rings.push(ring.into_iter().map(Edge::Planar).collect());
    }
    Some((rings, [low, high], [x, y]))
}
