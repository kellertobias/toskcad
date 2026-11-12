use super::*;
use crate::normal_relief::Kind;
use crate::spline_sweep::Sweep;
/// Reconstructed retained source curves; no authored recipe or feature roles.
#[derive(Clone, Debug, PartialEq)]
pub struct ImportedSplineCapRelief {
    pub origin_mm: Point,
    pub u: Point,
    pub v: Point,
    pub relief: Relief,
}
pub fn read_spline_cap_relief(
    input: &str,
    tolerance: f64,
) -> Result<(Readback, ImportedSplineCapRelief)> {
    compound::read_members_by(input, tolerance, true, true, |shell, id, _| {
        read_member(shell, id)
    })?
    .pop()
    .ok_or(ReadError::Invalid)
}
pub fn read_spline_cap_relief_compound(
    input: &str,
    tolerance: f64,
) -> Result<Vec<(Readback, ImportedSplineCapRelief)>> {
    compound::read_members_by(input, tolerance, true, false, |shell, id, _| {
        read_member(shell, id)
    })
}
pub(super) fn error(error: crate::spline_sweep::ConstructionError) -> ReadError {
    match error {
        crate::spline_sweep::ConstructionError::Limit => ReadError::Limit,
        _ => ReadError::Invalid,
    }
}
fn opposite(face: &FaceWire, edge: u32) -> Result<u32> {
    require(face.holes.is_empty() && face.outer.len() == 4)?;
    let at = face
        .outer
        .iter()
        .position(|e| e.edge == edge)
        .ok_or(ReadError::Invalid)?;
    Ok(face.outer[(at + 2) % 4].edge)
}
pub(super) fn other(owners: &BTreeMap<u32, Vec<usize>>, edge: u32, face: usize) -> Result<usize> {
    let ids = owners.get(&edge).ok_or(ReadError::Invalid)?;
    let [a, b] = ids[..] else {
        return Err(ReadError::Invalid);
    };
    if a == face && b != face {
        Ok(b)
    } else if b == face && a != face {
        Ok(a)
    } else {
        Err(ReadError::Invalid)
    }
}
fn far_rim(
    faces: &[FaceWire],
    owners: &BTreeMap<u32, Vec<usize>>,
    source: usize,
    ring: &[Traversal],
    work: &mut usize,
) -> Result<Vec<u32>> {
    ring.iter()
        .map(|edge| {
            charge(work, 1)?;
            let wall = other(owners, edge.edge, source)?;
            opposite(&faces[wall], edge.edge)
        })
        .collect()
}
fn cap_for_rim(
    faces: &[FaceWire],
    owners: &BTreeMap<u32, Vec<usize>>,
    source: usize,
    ring: &[Traversal],
    far: &[u32],
) -> Result<usize> {
    let edge = ring.first().ok_or(ReadError::Invalid)?;
    let wall = other(owners, edge.edge, source)?;
    let cap = other(owners, far[0], wall)?;
    let mut actual = faces[cap].outer.iter().map(|e| e.edge).collect::<Vec<_>>();
    actual.sort_unstable();
    let mut expected = far.to_vec();
    expected.sort_unstable();
    require(actual == expected)?;
    Ok(cap)
}
fn section(
    shell: &Shell<'_, '_>,
    edges: &Edges,
    plane: &Plane,
    face: &FaceWire,
    reverse: bool,
    height: f64,
) -> Result<Sweep> {
    let outline =
        spline_sweep::reconstruction::section_edges(shell, edges, plane, &face.outer, reverse)?;
    let holes = face
        .holes
        .iter()
        .map(|ring| spline_sweep::reconstruction::section_edges(shell, edges, plane, ring, reverse))
        .collect::<Result<Vec<_>>>()?;
    Sweep::from_verified_rings(outline, holes, height, shell.tolerance)
        .map_err(|_| ReadError::Invalid)
}
pub(in crate::step_read) fn read_member(
    shell: &mut Shell<'_, '_>,
    solid: u32,
) -> Result<(Readback, ImportedSplineCapRelief)> {
    let ids = shell.member_faces(solid)?;
    if ids.len() > 2 + 6 * crate::spline_cap_relief::MAX_EDGES {
        return Err(ReadError::Limit);
    }
    let mut edges = Edges::new();
    let faces = ids
        .iter()
        .map(|&id| spline_sweep::read_face(shell, id, &mut edges))
        .collect::<Result<Vec<_>>>()?;
    require(edges.values().all(|edge| edge.uses == [1, 1]))?;
    let spline = faces
        .iter()
        .find(|face| {
            shell
                .simple(face.surface)
                .is_ok_and(|e| e.name == "SURFACE_OF_LINEAR_EXTRUSION")
        })
        .ok_or(ReadError::Unsupported)?;
    let (_, axis) = shell.spline_extrusion(spline.surface)?;
    let tol = shell.tolerance;
    let mut caps = Vec::new();
    let mut work = MAX_WORK;
    for (i, face) in faces.iter().enumerate() {
        charge(&mut work, 1)?;
        if shell.simple(face.surface)?.name == "PLANE" {
            let plane = spline_sweep::face_plane(shell, face)?;
            if length(cross(plane.normal, axis)) <= 1e-10 {
                caps.push((i, plane));
            }
        }
    }
    let bottom = caps
        .iter()
        .filter(|(_, p)| dot(p.normal, axis) < 0.)
        .collect::<Vec<_>>();
    let &[&(bottom, ref plane)] = bottom.as_slice() else {
        return Err(ReadError::Unsupported);
    };
    let frame = Frame::new(plane.origin, plane.u, plane.v.map(|n| -n))?;
    let [u, v, n] = frame.axes;
    let determinant = dot(u, cross(v, n));
    // Pull actual retained controls through the true inverse, not a transpose.
    let section_plane = Plane {
        origin: frame.origin,
        u: cross(v, n).map(|x| x / determinant),
        v: cross(n, u).map(|x| x / determinant),
        normal: n,
    };
    let mut owners: BTreeMap<u32, Vec<usize>> = BTreeMap::new();
    for (i, face) in faces.iter().enumerate() {
        for ring in std::iter::once(&face.outer).chain(&face.holes) {
            for edge in ring {
                charge(&mut work, 1)?;
                owners.entry(edge.edge).or_default().push(i);
            }
        }
    }
    let far = far_rim(&faces, &owners, bottom, &faces[bottom].outer, &mut work)?;
    let entry = cap_for_rim(&faces, &owners, bottom, &faces[bottom].outer, &far)?;
    let entry_plane = spline_sweep::face_plane(shell, &faces[entry])?;
    let h = frame.local(entry_plane.origin)[2];
    require(h > tol && aligned(entry_plane.normal, n, h, tol / 16.))?;
    let base = section(shell, &edges, &section_plane, &faces[bottom], true, h)?;
    let mut stock_top = BTreeSet::new();
    for ring in std::iter::once(&faces[bottom].outer).chain(&faces[bottom].holes) {
        stock_top.extend(far_rim(&faces, &owners, bottom, ring, &mut work)?);
    }
    let mut motifs = Vec::new();
    let mut kind = None;
    let mut used = BTreeSet::new();
    for ring in &faces[entry].holes {
        charge(&mut work, ring.len())?;
        if ring.iter().all(|e| stock_top.contains(&e.edge)) {
            continue;
        }
        let far = far_rim(&faces, &owners, entry, ring, &mut work)?;
        let cap = cap_for_rim(&faces, &owners, entry, ring, &far)?;
        require(used.insert(cap))?;
        let plane = spline_sweep::face_plane(shell, &faces[cap])?;
        let level = frame.local(plane.origin)[2];
        let delta = level - h;
        require(
            delta.is_finite()
                && delta.abs() > tol
                && aligned(plane.normal, n, h.max(level.abs()), tol / 16.),
        )?;
        let current = if delta > 0. {
            Kind::Raised
        } else {
            Kind::Recessed
        };
        if kind.is_some_and(|prior| prior != current) {
            return Err(ReadError::Unsupported);
        }
        kind = Some(current);
        motifs.push(section(
            shell,
            &edges,
            &section_plane,
            &faces[cap],
            false,
            delta.abs(),
        )?);
    }
    let kind = kind.ok_or(ReadError::Unsupported)?;
    let minimum = tol * 4.;
    let relief = Relief::new(base, motifs, kind, tol, minimum).map_err(error)?;
    // Reparse into a clean member shell and prove the complete reconstructed
    // physical graph. Candidate recognition alone admits no geometry.
    let mut clean = Shell::new(shell.entities, tol, true);
    verify_member(&mut clean, solid, &relief, [frame.origin, u, v])?;
    let seams = owners
        .iter()
        .filter_map(|(&id, faces)| matches!(faces.as_slice(),[a,b] if a==b).then_some(id))
        .collect();
    let readback = spline_sweep::metrics(shell, edges, &seams, ids.len(), relief.volume_mm3())?;
    Ok((
        readback,
        ImportedSplineCapRelief {
            origin_mm: frame.origin,
            u,
            v,
            relief,
        },
    ))
}
