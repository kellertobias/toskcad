//! Tessellate consolidated stock faces and preserve both rims of each drill.

use super::stadium::StadiumRing;
use super::*;
use crate::{
    helical_thread::Thread,
    helical_thread_mesh::{open_blind_cavity, open_counterbored_cavity},
    orthogonal_boundary::{Boundary, BoundaryFace, point_key},
    orthogonal_pockets::{FacePlane, MAX_POCKETS, OrthogonalPocket, OrthogonalThroughBore},
    orthogonal_slots::OrthogonalSlot,
};
use std::collections::BTreeMap;

type Point = [f64; 3];
type SharedCorners = BTreeMap<(u32, [u64; 3]), usize>;

/// Tessellate the faces of an orthogonal solid with certified round pockets.
/// Face IDs preserve the original boundary order, then one wall/floor pair
/// per pocket. Only adjacent faces of the same stock component share corners.
pub fn pocketed_orthogonal(
    boundary: &Boundary,
    pockets: &[OrthogonalPocket],
    options: Tessellation,
) -> Result<DisplayMesh, MeshError> {
    slotted_orthogonal(boundary, pockets, &[], options)
}

/// Tessellate orthogonal stock with both circular and rounded-slot recesses.
/// Slots append four profile walls followed by the blind floor.
pub fn slotted_orthogonal(
    boundary: &Boundary,
    pockets: &[OrthogonalPocket],
    slots: &[OrthogonalSlot],
    options: Tessellation,
) -> Result<DisplayMesh, MeshError> {
    machined_orthogonal(boundary, pockets, slots, &[], None, options)
}

/// Tessellate orthogonal stock with blind pockets and paired through bores.
/// A bore appends one wall face after all blind pocket faces.
pub fn bored_orthogonal(
    boundary: &Boundary,
    pockets: &[OrthogonalPocket],
    through: &[OrthogonalThroughBore],
    options: Tessellation,
) -> Result<DisplayMesh, MeshError> {
    machined_orthogonal(boundary, pockets, &[], through, None, options)
}

/// One certified through bore whose middle interval is replaced by an exact
/// internal helical wall. `origin` and `axes` place the thread tool in stock
/// coordinates; its axis is `axes[0] x axes[1]` and must be the bore axis.
#[derive(Clone, Copy)]
pub struct ThreadedThroughBore<'a> {
    pub bore: usize,
    pub wall: BandWall<'a>,
    pub origin: Point,
    pub axes: [Point; 2],
}

/// The wall that replaces a bore's middle interval: the exact helical
/// thread, or a display-only smooth cylinder of the same tool length.
#[derive(Clone, Copy)]
pub enum BandWall<'a> {
    Thread(&'a Thread),
    Cylinder { radius_mm: f64, length_mm: f64 },
}
impl BandWall<'_> {
    fn length_mm(self) -> f64 {
        match self {
            Self::Thread(thread) => thread.spec().length_mm,
            Self::Cylinder { length_mm, .. } => length_mm,
        }
    }
    fn inner_radius_mm(self) -> f64 {
        match self {
            Self::Thread(thread) => thread.root_radius_mm(),
            Self::Cylinder { radius_mm, .. } => radius_mm,
        }
    }
    /// Outward wall in tool coordinates with counterclockwise rims at z=0
    /// and z=length; faces are numbered from 2 like the thread envelope.
    fn open(self, options: Tessellation) -> Result<(DisplayMesh, [Vec<u32>; 2]), MeshError> {
        let (radius, length) = match self {
            Self::Thread(thread) => return crate::helical_thread_mesh::open_wall(thread, options),
            Self::Cylinder {
                radius_mm,
                length_mm,
            } => (radius_mm, length_mm),
        };
        if !(radius.is_finite() && radius > 0. && length.is_finite() && length > 0.) {
            return Err(MeshError::UnresolvedBoundary);
        }
        let count = options.arc_segments(radius, TAU)?;
        let mut mesh = DisplayMesh {
            vertices: Vec::with_capacity(6 * count),
            indices: Vec::with_capacity(6 * count),
            face_indices: Vec::with_capacity(2 * count),
        };
        for z in [0., length] {
            for i in 0..count {
                let (sin, cos) = (TAU * i as f64 / count as f64).portable_sin_cos();
                mesh.vertices.extend([radius * cos, radius * sin, z]);
            }
        }
        for i in 0..count {
            let next = (i + 1) % count;
            mesh.triangle([i, next, count + i], 2);
            mesh.triangle([next, count + next, count + i], 2);
        }
        let ring = |base: usize| (0..count).map(|i| (base + i) as u32).collect();
        Ok((mesh, [ring(0), ring(count)]))
    }
}

/// Tessellate bored orthogonal stock where one bore carries a buried thread
/// band. Face IDs follow [`bored_orthogonal`]; the threaded bore's own wall ID
/// is its entry collar. Then follow its exit collar, the steps at thread
/// tool caps 0 and 1, and the helical bands in thread face order from 2.
pub fn threaded_bored_orthogonal(
    boundary: &Boundary,
    pockets: &[OrthogonalPocket],
    through: &[OrthogonalThroughBore],
    threaded: ThreadedThroughBore<'_>,
    options: Tessellation,
) -> Result<DisplayMesh, MeshError> {
    machined_orthogonal(boundary, pockets, &[], through, Some(threaded), options)
}

/// Stock ordinals in a blind tap shell. The mouth keeps face zero so the
/// helical band and blind floor ordinals match the original thread envelope.
pub fn blind_tap_stock_face(face: usize, mouth: usize, thread_faces: usize) -> u32 {
    if face == mouth {
        0
    } else {
        (thread_faces + face - usize::from(face > mouth)) as u32
    }
}

/// Close an open exact helical cavity against a certified orthogonal stock
/// face. Both sides of the mouth reuse the same indexed helical rim.
pub fn tapped_orthogonal(
    boundary: &Boundary,
    thread: &Thread,
    pocket: OrthogonalPocket,
    [u, v]: [Point; 2],
    options: Tessellation,
) -> Result<DisplayMesh, MeshError> {
    tapped_stock(boundary, thread, pocket, [u, v], options, None)
}
/// A counterbored blind tap shares all finite shoulder and stock-mouth vertices.
pub fn counterbored_tapped_orthogonal(
    boundary: &Boundary,
    seat: &crate::thread_counterbore::Counterbore,
    pocket: OrthogonalPocket,
    axes: [Point; 2],
    options: Tessellation,
) -> Result<DisplayMesh, MeshError> {
    tapped_stock(boundary, seat.thread(), pocket, axes, options, Some(seat))
}
fn tapped_stock(
    boundary: &Boundary,
    thread: &Thread,
    pocket: OrthogonalPocket,
    [u, v]: [Point; 2],
    options: Tessellation,
    seat: Option<&crate::thread_counterbore::Counterbore>,
) -> Result<DisplayMesh, MeshError> {
    let invalid = MeshError::UnresolvedBoundary;
    let mouth = boundary.faces.get(pocket.face).ok_or(invalid)?;
    let mouth_plane = FacePlane::from_face(mouth).ok_or(invalid)?;
    let spec = thread.spec();
    let dot = |a: Point, b: Point| a.into_iter().zip(b).map(|(a, b)| a * b).sum::<f64>();
    let inward = [
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    ];
    if spec.length_mm != pocket.depth_mm
        || pocket.center_mm[mouth_plane.axis] != mouth_plane.level
        || u[mouth_plane.axis] != 0.
        || v[mouth_plane.axis] != 0.
        || ![
            dot(u, u) - 1.,
            dot(v, v) - 1.,
            dot(u, v),
            1. + dot(inward, mouth.normal),
        ]
        .into_iter()
        .all(|v| v.abs() <= 1e-10)
        || !mouth_plane.contains_disk(
            mouth,
            pocket.center_mm,
            seat.map_or(spec.major_radius_mm, |s| s.radius_mm()),
        )
    {
        return Err(invalid);
    }
    let mut normal = mouth.normal;
    normal[mouth_plane.axis] *= -1.;
    let (mut mesh, rim) = match seat {
        Some(seat) => open_counterbored_cavity(seat, options)?,
        None => open_blind_cavity(thread, options)?,
    };
    for point in mesh.vertices.as_chunks_mut::<3>().0 {
        let [x, y, z] = *point;
        *point = std::array::from_fn(|i| pocket.center_mm[i] + u[i] * x + v[i] * y + normal[i] * z);
    }
    let mut corners = BTreeMap::new();
    let thread_faces = seat.map_or(2 + 4 * spec.starts as usize, |s| {
        s.cavity_face_count() as usize
    });
    for (i, face) in boundary.faces.iter().enumerate() {
        let plane = FacePlane::from_face(face).ok_or(invalid)?;
        let mut region = stock_region(&mut mesh, &mut corners, face, plane)?;
        let is_mouth = i == pocket.face;
        if is_mouth {
            region.ring(rim.iter().map(|&index| {
                let index = index as usize;
                let point = mesh.vertices.as_chunks::<3>().0[index];
                (plane.project(point), index)
            }));
        }
        let holes = face.rings.len() - 1 + usize::from(is_mouth);
        let triangulation = triangulate_region(region.points, &region.edges, holes)?;
        if mesh
            .face_indices
            .len()
            .saturating_add(triangulation.triangles.len())
            > options.triangle_limit
        {
            return Err(MeshError::BudgetExceeded);
        }
        for triangle in triangulation.triangles {
            mesh.triangle(
                triangle.map(|vertex| region.indices[vertex]),
                blind_tap_stock_face(i, pocket.face, thread_faces),
            );
        }
    }
    Ok(mesh)
}

fn stock_region(
    mesh: &mut DisplayMesh,
    corners: &mut SharedCorners,
    face: &BoundaryFace,
    plane: FacePlane,
) -> Result<FaceRegion, MeshError> {
    let mut region = FaceRegion::default();
    for ring in &face.rings {
        let mut projected = Vec::with_capacity(ring.len());
        for &point in ring {
            let uv = plane.project(point);
            if !point.iter().all(|x| x.is_finite()) || plane.lift(uv) != point {
                return Err(MeshError::UnresolvedBoundary);
            }
            let index = stock_corner(mesh, corners, face, point);
            projected.push((uv, index));
        }
        region.ring(projected.into_iter());
    }
    Ok(region)
}

fn machined_orthogonal(
    boundary: &Boundary,
    pockets: &[OrthogonalPocket],
    slots: &[OrthogonalSlot],
    through: &[OrthogonalThroughBore],
    threaded: Option<ThreadedThroughBore<'_>>,
    options: Tessellation,
) -> Result<DisplayMesh, MeshError> {
    options.check()?;
    if boundary.faces.is_empty() || threaded.is_some_and(|t| t.bore >= through.len()) {
        return Err(MeshError::UnresolvedBoundary);
    }
    if pockets
        .len()
        .saturating_add(slots.len())
        .saturating_add(through.len())
        > MAX_POCKETS
        || boundary.faces.len() > options.triangle_limit
    {
        return Err(MeshError::BudgetExceeded);
    }

    let mut openings = vec![FaceOpenings::default(); boundary.faces.len()];
    let mut counts = Vec::with_capacity(pockets.len());
    let mut pocket_triangles = 0;
    for (index, pocket) in pockets.iter().enumerate() {
        if !pocket.radius_mm.is_finite()
            || pocket.radius_mm <= 0.
            || !pocket.depth_mm.is_finite()
            || pocket.depth_mm <= 0.
            || pocket.center_mm.iter().any(|x| !x.is_finite())
        {
            return Err(MeshError::UnresolvedBoundary);
        }
        let count = options.arc_segments(pocket.radius_mm, TAU)?;
        let face = openings
            .get_mut(pocket.face)
            .ok_or(MeshError::UnresolvedBoundary)?;
        face.circles.push(index);
        counts.push(count);
        // The mouth contributes count + 2 triangles on its parent face.
        // The floor and wall contribute another 3 * count.
        pocket_triangles += 4 * count + 2;
        if pocket_triangles > options.triangle_limit {
            return Err(MeshError::BudgetExceeded);
        }
    }
    let mut slot_rings = Vec::with_capacity(slots.len());
    for (index, slot) in slots.iter().enumerate() {
        let face = boundary
            .faces
            .get(slot.face)
            .ok_or(MeshError::UnresolvedBoundary)?;
        let profile = slot
            .profile(face)
            .map_err(|_| MeshError::UnresolvedBoundary)?;
        let ring = StadiumRing::new(&profile, options)?;
        if !slot.depth_mm.is_finite() || slot.depth_mm <= 0. {
            return Err(MeshError::UnresolvedBoundary);
        }
        pocket_triangles += 4 * ring.points.len() + 2;
        if pocket_triangles > options.triangle_limit {
            return Err(MeshError::BudgetExceeded);
        }
        openings[slot.face].slots.push(index);
        slot_rings.push(ring);
    }
    let mut bore_counts = Vec::with_capacity(through.len());
    for (index, bore) in through.iter().enumerate() {
        check_bore(boundary, bore)?;
        let count = options.arc_segments(bore.radius_mm, TAU)?;
        pocket_triangles += 4 * count + 4;
        if pocket_triangles > options.triangle_limit {
            return Err(MeshError::BudgetExceeded);
        }
        bore_counts.push(count);
        for (side, face) in bore.faces.into_iter().enumerate() {
            openings[face].bores.push((index, side));
        }
    }
    let mut total = pocket_triangles;
    for face in &boundary.faces {
        if face.rings.is_empty() || face.rings.iter().any(|ring| ring.len() < 3) {
            return Err(MeshError::UnresolvedBoundary);
        }
        total += face.rings.iter().map(Vec::len).sum::<usize>() + 2 * face.rings.len() - 4;
        if total > options.triangle_limit {
            return Err(MeshError::BudgetExceeded);
        }
    }
    let mut mesh = DisplayMesh {
        vertices: Vec::new(),
        indices: Vec::with_capacity(3 * total),
        face_indices: Vec::with_capacity(total),
    };
    let mut corners = BTreeMap::new();
    let mut mouths = vec![0..0; pockets.len()];
    let mut slot_mouths = vec![0..0; slots.len()];
    let mut bore_mouths = vec![[0..0, 0..0]; through.len()];
    for (face_index, face) in boundary.faces.iter().enumerate() {
        let plane = FacePlane::from_face(face).ok_or(MeshError::UnresolvedBoundary)?;
        let mut region = stock_region(&mut mesh, &mut corners, face, plane)?;
        let opening = &openings[face_index];
        for &index in &opening.circles {
            let pocket = pockets[index];
            mouths[index] = round_mouth(
                &mut mesh,
                &mut region,
                plane,
                pocket.center_mm,
                pocket.radius_mm,
                counts[index],
                false,
            )?;
        }
        for &index in &opening.slots {
            let ring = &slot_rings[index];
            let start = mesh.vertices.len() / 3;
            region.ring(ring.points.iter().rev().enumerate().map(|(i, &uv)| {
                mesh.vertices.extend(plane.lift(uv));
                (uv, start + i)
            }));
            slot_mouths[index] = start..start + ring.points.len();
        }
        for &(index, side) in &opening.bores {
            let bore = through[index];
            bore_mouths[index][side] = round_mouth(
                &mut mesh,
                &mut region,
                plane,
                bore.centers_mm[side],
                bore.radius_mm,
                bore_counts[index],
                side == 1,
            )?;
        }
        let holes =
            face.rings.len() + opening.circles.len() + opening.slots.len() + opening.bores.len()
                - 1;
        for triangle in triangulate_region(region.points, &region.edges, holes)?.triangles {
            mesh.triangle(
                triangle.map(|index| region.indices[index]),
                face_index as u32,
            );
        }
    }
    for (index, (pocket, mouth)) in pockets.iter().zip(mouths).enumerate() {
        let normal = boundary.faces[pocket.face].normal;
        let wall_face = (boundary.faces.len() + 2 * index) as u32;
        add_recess_surfaces(
            &mut mesh,
            mouth,
            normal,
            pocket.depth_mm,
            pocket.center_mm,
            std::iter::repeat(wall_face),
            wall_face + 1,
        )?;
    }
    let first_slot_face = boundary.faces.len() + 2 * pockets.len();
    for (index, ((slot, ring), mouth)) in slots.iter().zip(slot_rings).zip(slot_mouths).enumerate()
    {
        let face = &boundary.faces[slot.face];
        let [a, b] = slot.centers_mm;
        let center = std::array::from_fn(|i| a[i].midpoint(b[i]));
        let first_wall = (first_slot_face + 5 * index) as u32;
        // The mouth ring winds clockwise around the stock face.
        let walls = ring
            .sides
            .iter()
            .rev()
            .cycle()
            .skip(1)
            .map(|&side| first_wall + side);
        add_recess_surfaces(
            &mut mesh,
            mouth,
            face.normal,
            slot.depth_mm,
            center,
            walls,
            first_wall + 4,
        )?;
    }
    let first_bore_face = first_slot_face + 5 * slots.len();
    let first_thread_face = (first_bore_face + through.len()) as u32;
    for (index, [mouth, exit]) in bore_mouths.into_iter().enumerate() {
        let count = mouth.len();
        let face = (first_bore_face + index) as u32;
        if let Some(threaded) = threaded.filter(|t| t.bore == index) {
            threaded_wall(
                &mut mesh,
                [mouth, exit],
                through[index],
                threaded,
                [face, first_thread_face],
                options,
            )?;
            continue;
        }
        for i in 0..count {
            let next = (i + 1) % count;
            let near = [i, next].map(|j| mouth.start + j);
            let far = [i, next].map(|j| exit.start + (count - j) % count);
            wall(&mut mesh, far, near, face);
        }
    }
    Ok(mesh)
}

/// Entry collar, planar step, inward helical wall, step and exit collar.
/// Collar rings repeat the mouth samples along the exact bore axis; steps
/// share the bore rings and the helical rims vertex for vertex.
fn threaded_wall(
    mesh: &mut DisplayMesh,
    [mouth, exit]: [Range<usize>; 2],
    bore: OrthogonalThroughBore,
    threaded: ThreadedThroughBore<'_>,
    [entry_face, first]: [u32; 2],
    options: Tessellation,
) -> Result<(), MeshError> {
    let invalid = MeshError::UnresolvedBoundary;
    let dot = |a: Point, b: Point| a.into_iter().zip(b).map(|(a, b)| a * b).sum::<f64>();
    let [c0, c1] = bore.centers_mm;
    let depth = bore.depth_mm();
    let d: Point = std::array::from_fn(|i| (c1[i] - c0[i]) / depth);
    let axis = d.iter().position(|x| x.abs() == 1.).ok_or(invalid)?;
    let [u, v] = threaded.axes;
    let n = [
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    ];
    let along = dot(n, d);
    let length = threaded.wall.length_mm();
    let offset: Point = std::array::from_fn(|i| threaded.origin[i] - c0[i]);
    let start = dot(offset, d);
    let radial = std::array::from_fn::<f64, 3, _>(|i| offset[i] - start * d[i]);
    let end = start + along * length;
    let [low, high] = [start.min(end), start.max(end)];
    if ![dot(u, u) - 1., dot(v, v) - 1., dot(u, v), along.abs() - 1.]
        .into_iter()
        .all(|x| x.abs() <= 1e-10)
        || dot(radial, radial).sqrt() > 1e-9 * (1. + depth)
        || low <= 0.
        || high >= depth
        || bore.radius_mm >= threaded.wall.inner_radius_mm()
    {
        return Err(invalid);
    }
    let count = mouth.len();
    let level = |s: f64| c0[axis] + d[axis] * s;
    // Ring vertex j shares the angular sample of mouth vertex j.
    let ring = |mesh: &mut DisplayMesh, s: f64| {
        let first = mesh.vertices.len() / 3;
        for j in 0..count {
            let mut point: Point =
                std::array::from_fn(|i| mesh.vertices[3 * (mouth.start + j) + i]);
            point[axis] = level(s);
            mesh.vertices.extend(point);
        }
        first
    };
    let lower = ring(mesh, low);
    let upper = ring(mesh, high);
    for i in 0..count {
        let next = (i + 1) % count;
        let pair = |base: usize| [base + i, base + next];
        wall(mesh, pair(lower), pair(mouth.start), entry_face);
        let far = [i, next].map(|j| exit.start + (count - j) % count);
        wall(mesh, far, pair(upper), first);
    }
    let (wall_mesh, rims) = threaded.wall.open(options)?;
    if mesh
        .face_indices
        .len()
        .saturating_add(wall_mesh.face_indices.len())
        > options.triangle_limit
    {
        return Err(MeshError::BudgetExceeded);
    }
    let base = mesh.vertices.len() / 3;
    for point in wall_mesh.vertices.as_chunks::<3>().0 {
        let [x, y, z] = *point;
        mesh.vertices.extend::<[f64; 3]>(std::array::from_fn(|i| {
            threaded.origin[i] + u[i] * x + v[i] * y + n[i] * z
        }));
    }
    for (triangle, &face) in wall_mesh
        .indices
        .as_chunks::<3>()
        .0
        .iter()
        .zip(&wall_mesh.face_indices)
    {
        let face = face.checked_sub(2).ok_or(invalid)?;
        // The cut wall faces the bore axis, opposite to the external wall.
        let [a, b, c] = triangle.map(|i| base + i as usize);
        mesh.triangle([a, c, b], first + 3 + face);
    }
    // Each step faces into the tool: +n at cap 0, -n at cap 1.
    for (cap, rim) in rims.iter().enumerate() {
        let s = if cap == 0 { start } else { end };
        let outward = if cap == 0 { along } else { -along };
        let plane = FacePlane {
            axis,
            sign: outward * d[axis],
            level: level(s),
        };
        let collar = if s == low { lower } else { upper };
        let mut region = FaceRegion::default();
        let signed = |points: &[([f64; 2], usize)]| {
            (0..points.len())
                .map(|i| {
                    let [a, b] = [points[i].0, points[(i + 1) % points.len()].0];
                    a[0] * b[1] - a[1] * b[0]
                })
                .sum::<f64>()
        };
        let mut outer: Vec<_> = rim
            .iter()
            .map(|&i| {
                let index = base + i as usize;
                let point = std::array::from_fn(|k| mesh.vertices[3 * index + k]);
                (plane.project(point), index)
            })
            .collect();
        let mut inner: Vec<_> = (0..count)
            .map(|j| {
                let index = collar + j;
                let point = std::array::from_fn(|k| mesh.vertices[3 * index + k]);
                (plane.project(point), index)
            })
            .collect();
        if signed(&outer) < 0. {
            outer.reverse();
        }
        if signed(&inner) > 0. {
            inner.reverse();
        }
        region.ring(outer.into_iter());
        region.ring(inner.into_iter());
        let triangles = triangulate_region(region.points, &region.edges, 1)?.triangles;
        if mesh.face_indices.len().saturating_add(triangles.len()) > options.triangle_limit {
            return Err(MeshError::BudgetExceeded);
        }
        for triangle in triangles {
            mesh.triangle(
                triangle.map(|vertex| region.indices[vertex]),
                first + 1 + cap as u32,
            );
        }
    }
    Ok(())
}

#[derive(Clone, Default)]
struct FaceOpenings {
    circles: Vec<usize>,
    slots: Vec<usize>,
    bores: Vec<(usize, usize)>,
}

#[derive(Default)]
struct FaceRegion {
    points: Vec<Point2<f64>>,
    edges: Vec<[usize; 2]>,
    indices: Vec<usize>,
}

impl FaceRegion {
    fn ring(&mut self, ring: impl ExactSizeIterator<Item = ([f64; 2], usize)>) {
        let first = self.points.len();
        let end = first + ring.len();
        for (i, ([u, v], vertex)) in (first..).zip(ring) {
            self.points.push(Point2::new(u, v));
            self.edges
                .push([i, if i + 1 == end { first } else { i + 1 }]);
            self.indices.push(vertex);
        }
    }
}

fn check_bore(boundary: &Boundary, bore: &OrthogonalThroughBore) -> Result<(), MeshError> {
    let [front, back] = bore.faces.map(|face| {
        boundary
            .faces
            .get(face)
            .ok_or(MeshError::UnresolvedBoundary)
    });
    let [front, back] = [front?, back?];
    let normal = front.normal;
    let [a, b] = bore.centers_mm;
    let movement: Point = std::array::from_fn(|i| b[i] - a[i]);
    let depth: f64 = movement
        .into_iter()
        .zip(normal)
        .map(|(distance, normal)| -distance * normal)
        .sum();
    if front.component != back.component
        || back.normal != normal.map(|n| -n)
        || bore.radius_mm <= 0.
        || !bore.radius_mm.is_finite()
        || !depth.is_finite()
        || depth <= 0.
        || movement != normal.map(|n| -depth * n)
    {
        return Err(MeshError::UnresolvedBoundary);
    }
    Ok(())
}

fn round_mouth(
    mesh: &mut DisplayMesh,
    region: &mut FaceRegion,
    plane: FacePlane,
    center: Point,
    radius: f64,
    count: usize,
    opposite: bool,
) -> Result<Range<usize>, MeshError> {
    let uv = plane.project(center);
    if plane.lift(uv) != center {
        return Err(MeshError::UnresolvedBoundary);
    }
    let start = mesh.vertices.len() / 3;
    region.ring((0..count).map(|i| {
        // Opposite normals reverse the planar V axis. Reverse the sample
        // index as well to use the identical world offset on both mouths.
        let sample = if opposite { (count - i) % count } else { i };
        let (sin, cos) = (-TAU * sample as f64 / count as f64).portable_sin_cos();
        let sign = if opposite { -1. } else { 1. };
        let point = [uv[0] + radius * cos, uv[1] + sign * radius * sin];
        mesh.vertices.extend(plane.lift(point));
        (point, start + i)
    }));
    Ok(start..start + count)
}

fn stock_corner(
    mesh: &mut DisplayMesh,
    corners: &mut SharedCorners,
    face: &BoundaryFace,
    point: Point,
) -> usize {
    *corners
        .entry((face.component, point_key(point)))
        .or_insert_with(|| {
            let index = mesh.vertices.len() / 3;
            mesh.vertices.extend(point);
            index
        })
}

fn add_recess_surfaces(
    mesh: &mut DisplayMesh,
    mouth: Range<usize>,
    normal: Point,
    depth: f64,
    center: Point,
    wall_faces: impl Iterator<Item = u32>,
    floor_face: u32,
) -> Result<(), MeshError> {
    let count = mouth.len();
    let floor = mesh.vertices.len() / 3;
    let floor_center = std::array::from_fn(|i| center[i] - normal[i] * depth);
    if !floor_center.iter().all(|x| x.is_finite()) || floor_center == center {
        return Err(MeshError::UnresolvedBoundary);
    }
    for vertex in mouth.clone() {
        let point: Point = std::array::from_fn(|i| mesh.vertices[3 * vertex + i]);
        mesh.vertices
            .extend((0..3).map(|i| point[i] - normal[i] * depth));
    }
    mesh.vertices.extend(floor_center);
    for (i, wall_face) in (0..count).zip(wall_faces) {
        let next = (i + 1) % count;
        let bottom @ [a, b] = [floor + i, floor + next];
        wall(
            mesh,
            bottom,
            [mouth.start + i, mouth.start + next],
            wall_face,
        );
        mesh.triangle([floor + count, b, a], floor_face);
    }
    Ok(())
}

#[cfg(test)]
mod tests;
