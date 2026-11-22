//! A straight-walled pocket entering one cap of an analytic normal sweep.
//!
//! The pocket outline is a closed line ring lying strictly inside the
//! profile material apart from the holes it encloses. Above its floor the
//! section is the profile minus the pocket and the holes clear of it; below
//! the floor it is the whole profile. The enclosed holes therefore end on the
//! floor, which is the pocket outline with those holes removed. Every face is
//! the exact plane or cylinder of a profile edge, the pocket or a cap.
use super::*;

/// The walls of one analytic ring between two levels and its oriented
/// bounds on a face below (`bounds[0]`) and above (`bounds[1]`) it.
fn ring_walls(
    file: &mut StepFile,
    ring: &[Edge],
    levels: [f64; 2],
    placement: &Placement,
    material_inside: bool,
    outer_bound: [bool; 2],
) -> Result<(Vec<Entity>, [Entity; 2]), StepError> {
    let normal = placement.normal;
    if let [
        Edge::Circle {
            center, radius_mm, ..
        },
    ] = ring
    {
        let trim = CylindricalTrim::new(
            file,
            levels.map(|z| placement.point_f64(*center, z)),
            normal,
            placement.u,
            *radius_mm,
        )?;
        let bounds = [false, true].map(|top| {
            // The orientation of a circle follows its material side; the
            // bound kind only marks the face's outer loop.
            let edge = file.oriented_edge(trim.rims[usize::from(top)], top == material_inside);
            file.bound(&[edge], outer_bound[usize::from(top)])
        });
        return Ok((vec![trim.wall(file, material_inside)], bounds));
    }
    let ring: Vec<_> = ring.iter().copied().map(ProfileEdge::Analytic).collect();
    let low = Rim::new(file, &ring, placement, levels[0])?;
    let high = Rim::new(file, &ring, placement, levels[1])?;
    let bounds = [
        low.bound(file, false, outer_bound[0]),
        high.bound(file, true, outer_bound[1]),
    ];
    let rails = low
        .vertices
        .iter()
        .zip(&high.vertices)
        .map(|(&a, &b)| line_edge(file, a, b))
        .collect::<Result<Vec<_>, _>>()?;
    let mut walls = Vec::with_capacity(ring.len());
    for (i, edge) in ring.iter().enumerate() {
        let next = (i + 1) % ring.len();
        let boundary = [
            low.edges[i],
            (rails[next], true),
            (high.edges[i].0, !high.edges[i].1),
            (rails[i], false),
        ]
        .map(|(edge, forward)| file.oriented_edge(edge, forward));
        let bound = file.bound(&boundary, true);
        let face = match (edge, &low.supports[i]) {
            (
                ProfileEdge::Analytic(Edge::Arc {
                    center,
                    a,
                    clockwise,
                    ..
                }),
                WallSupport::Cylinder(axis),
            ) => {
                let [x, y] = std::array::from_fn(|j| a[j] - center[j]);
                let radius = x.portable_hypot(y);
                let axis = if *clockwise {
                    let center = file.point(placement.point_f64(*center, levels[0]));
                    file.axis(
                        center,
                        normal,
                        placement.vector([x / radius, y / radius, 0.]),
                    )
                } else {
                    *axis
                };
                let cylinder =
                    file.entity(format_args!("CYLINDRICAL_SURFACE('',{axis},{radius:.16E})"));
                file.face_on(&[bound], cylinder, !*clockwise)
            }
            (ProfileEdge::Analytic(Edge::Line { .. }), WallSupport::Plane) => {
                let tangent = direction(low.vertices[i].point, low.vertices[next].point)?;
                file.face(
                    &[bound],
                    low.vertices[i].cartesian,
                    cross(tangent, normal),
                    tangent,
                )
            }
            _ => return Err(StepError::InvalidGeometry),
        };
        walls.push(face);
    }
    Ok((walls, bounds))
}

/// Encode `sweep` with a pocket of outline `pocket` (a clockwise line ring,
/// oriented like a hole) entering the upper cap when `top`, else the lower
/// one, down to the placement level `floor`. `enclosed[i]` marks hole ring
/// `i + 1` as lying inside the pocket outline; the others lie clear of it.
#[allow(clippy::too_many_arguments)]
pub(in crate::step) fn encode(
    file: &mut StepFile,
    sweep: &Sweep,
    pocket: &[Edge],
    enclosed: &[bool],
    floor: f64,
    top: bool,
    placement: Placement,
    reverse: bool,
) -> Result<Entity, StepError> {
    let height = sweep.height_mm();
    let levels = if reverse { [-height, 0.] } else { [0., height] };
    let rings = sweep.rings();
    if !(levels[0] < floor && floor < levels[1])
        || enclosed.len() + 1 != rings.len()
        || pocket.len() < 3
        || pocket.iter().any(|edge| !matches!(edge, Edge::Line { .. }))
        || crate::normal_sweep::signed_area(pocket) >= 0.
    {
        return Err(StepError::InvalidGeometry);
    }
    let normal = placement.normal;
    // Faces facing down at `levels[0]` and up at `levels[1]`, and the floor.
    let mut caps = [Vec::new(), Vec::new()];
    let mut floor_bounds = Vec::new();
    let mut walls = Vec::new();
    let pocketed = if top {
        [floor, levels[1]]
    } else {
        [levels[0], floor]
    };
    let solid = if top {
        [levels[0], floor]
    } else {
        [floor, levels[1]]
    };
    // The pocket outline bounds the floor from outside.
    let floor_end = usize::from(top);
    let (faces, bounds) = ring_walls(
        file,
        pocket,
        pocketed,
        &placement,
        false,
        if top { [true, false] } else { [false, true] },
    )?;
    walls.extend(faces);
    let [low, high] = bounds;
    // The pocket mouth joins the cap after the profile rings.
    let mouth = if top {
        floor_bounds.push(low);
        high
    } else {
        floor_bounds.push(high);
        low
    };
    for (index, ring) in rings.iter().enumerate() {
        let outer = index == 0;
        let inside = !outer && enclosed[index - 1];
        let span = if inside { solid } else { levels };
        let (faces, [low, high]) = ring_walls(file, ring, span, &placement, outer, [outer, outer])?;
        walls.extend(faces);
        if !inside {
            caps[0].push(low);
            caps[1].push(high);
        } else if top {
            caps[0].push(low);
            floor_bounds.push(high);
        } else {
            floor_bounds.push(low);
            caps[1].push(high);
        }
    }
    caps[floor_end].push(mouth);
    let supports = levels.map(|z| file.point(placement.point_f64([0.; 2], z)));
    let floor_support = file.point(placement.point_f64([0.; 2], floor));
    let up = normal;
    let down = normal.map(|n| -n);
    let mut faces = vec![
        file.face(&caps[0], supports[0], down, placement.u),
        file.face(&caps[1], supports[1], up, placement.u),
        file.face(
            &floor_bounds,
            floor_support,
            if floor_end == 1 { up } else { down },
            placement.u,
        ),
    ];
    faces.extend(walls);
    Ok(file.solid(&faces))
}
