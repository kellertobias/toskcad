//! Cylindrical stock joined to the original rims of a through internal
//! thread. The same retained rim ownership, planar cap triangulation and
//! counted connectivity proof as the hex nut; the six exact flats become an
//! inscribed chord polygon whose sagitta and normal turn meet the request.
use super::*;
use std::borrow::Cow;

/// Chord count of the inscribed stock polygon: the sagitta `R(1-cos(pi/n))`
/// is within the requested deflection and the normal of every wall facet is
/// within half the requested angle of the exact cylinder normal across it.
pub fn sleeve_segments(radius_mm: f64, tessellation: Tessellation) -> Result<usize, Error> {
    tessellation.check().map_err(|_| Error::Invalid)?;
    if !radius_mm.is_finite() || radius_mm <= 0. {
        return Err(Error::Invalid);
    }
    let ratio = (tessellation.deflection / radius_mm).min(1.);
    let chord = (1. - ratio).portable_acos();
    let by_sagitta = if chord > 0. {
        (std::f64::consts::PI / chord).ceil()
    } else {
        f64::INFINITY
    };
    let by_angle = (std::f64::consts::PI / (tessellation.angle / 2.)).ceil();
    let n = by_sagitta.max(by_angle).max(8.);
    if !n.is_finite() || n > 1_000_000. {
        return Err(Error::Limit);
    }
    // Directed verification: sagitta and half turn from the chosen count.
    let half = I::point(std::f64::consts::PI) / I::point(n);
    let half = I::new(half.lo.next_down(), half.hi.next_up());
    let sagitta = I::point(radius_mm) * (I::point(1.) - half.cos());
    if sagitta.hi > tessellation.deflection || half.hi > tessellation.angle / 2. {
        // One more chord removes floating directed rounding at the frontier.
        let n = n + 1.;
        let half = I::point(std::f64::consts::PI) / I::point(n);
        let half = I::new(half.lo.next_down(), half.hi.next_up());
        let sagitta = I::point(radius_mm) * (I::point(1.) - half.cos());
        if sagitta.hi > tessellation.deflection || half.hi > tessellation.angle / 2. {
            return Err(Error::Unresolved);
        }
        return Ok(n as usize);
    }
    Ok(n as usize)
}

/// Join a retained through-bore trace to cylindrical stock of `radius_mm`
/// over the source's own axial cap interval. The bore's cap fans are
/// replaced by annuli that keep the source cap face ids; the chord wall is
/// face `max source face + 1`. Every new facet normal is checked against
/// the exact support and the complete source must clear the inscribed circle.
pub fn join_sleeve<'a, 's>(
    source: &'a TracedMesh<'s>,
    radius_mm: f64,
    options: MeshOptions,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<JoinedMeshCandidate<'a, 's>, Error> {
    options.tessellation.check().map_err(|_| Error::Invalid)?;
    if options.maximum_work == 0
        || options.maximum_work > 512_000_000
        || options.maximum_bytes == 0
        || options.maximum_bytes > 64 * 1024 * 1024
        || !radius_mm.is_finite()
        || radius_mm <= 0.
    {
        return Err(Error::Invalid);
    }
    let source = ThreadBoundaryMesh::Retained(source);
    let held = source
        .bytes()
        .checked_add(options.additional_resident_bytes)
        .ok_or(Error::Limit)?;
    let mut work = Work {
        used: 0,
        held,
        peak: 0,
        options,
        interrupt,
    };
    work.storage(0)?;
    work.charge(1)?;
    let rings = rims(source, &mut work)?;
    work.held = work
        .held
        .checked_add(
            rings
                .vertices
                .iter()
                .map(|r| r.capacity() * 4)
                .sum::<usize>(),
        )
        .ok_or(Error::Limit)?;
    work.storage(0)?;
    if !rings.height[0].is_finite()
        || !rings.height[1].is_finite()
        || rings.height[0] >= rings.height[1]
    {
        return Err(Error::Unsupported);
    }
    let (mesh, connectivity, distance) = join_rims(
        Cow::Borrowed(source.mesh()),
        rings,
        radius_mm,
        None,
        &mut work,
    )?;
    Ok(JoinedMeshCandidate {
        source,
        mesh,
        connectivity,
        distance,
        work: work.used,
        peak: work.peak,
    })
}

/// Join a dedicated authored-law bore surface (no cap fans; rims at the
/// two cap planes, counterclockwise) to cylindrical stock. `faces` gives
/// the start cap, end cap and wall face ids; `source_distance_mm` is the
/// law mesh's proven facet distance. Returns the closed mesh, its
/// connectivity proof and the complete distance bound.
pub fn join_sleeve_law(
    law: crate::axial_pipe_thread::authored::LawMesh,
    heights: [f64; 2],
    faces: [u32; 3],
    radius_mm: f64,
    options: MeshOptions,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<(DisplayMesh, MeshProof, f64), Error> {
    options.tessellation.check().map_err(|_| Error::Invalid)?;
    if options.maximum_work == 0
        || options.maximum_work > 512_000_000
        || options.maximum_bytes == 0
        || options.maximum_bytes > 64 * 1024 * 1024
        || !radius_mm.is_finite()
        || radius_mm <= 0.
        || law.proof.is_some()
    {
        return Err(Error::Invalid);
    }
    // The law mesh is consumed: it becomes the joined mesh in place.
    let mut work = Work {
        used: 0,
        held: options.additional_resident_bytes,
        peak: 0,
        options,
        interrupt,
    };
    work.storage(0)?;
    work.charge(1)?;
    let rings = Rims {
        vertices: law.rims,
        caps: [faces[0], faces[1]],
        height: heights,
        unchanged_source_bound: law.deviation_mm,
        // Rim chords are edges of law facets: the same proven distance.
        distance: law.deviation_mm,
    };
    join_rims(
        Cow::Owned(law.mesh),
        rings,
        radius_mm,
        Some(faces[2]),
        &mut work,
    )
}

/// Shared stock join. `wall` is the explicit wall face for a source without
/// cap fans; otherwise the source's cap fans are replaced and the wall is
/// face `max source face + 1`.
fn join_rims(
    source: Cow<'_, DisplayMesh>,
    mut rings: Rims,
    radius_mm: f64,
    wall: Option<u32>,
    work: &mut Work<'_>,
) -> Result<(DisplayMesh, MeshProof, f64), Error> {
    let options = work.options;
    let fans = wall.is_none();
    // An owned fan-free source is extended in place, without a copy or a
    // compaction pass; a borrowed trace keeps both.
    let in_place = !fans && matches!(source, Cow::Owned(_));
    let source_mesh: &DisplayMesh = &source;
    if !rings.height[0].is_finite()
        || !rings.height[1].is_finite()
        || rings.height[0] >= rings.height[1]
    {
        return Err(Error::Unsupported);
    }
    let n = sleeve_segments(radius_mm, options.tessellation)?;
    // Each cap fan has one triangle per rim coedge; an annulus between a
    // rim of m and the n-gon has n+m triangles, the wall has 2n. A source
    // without fans gains the rim counts as well.
    let rim_triangles = if wall.is_some() {
        rings.vertices[0].len() + rings.vertices[1].len()
    } else {
        0
    };
    let expected_triangles = source_mesh
        .face_indices
        .len()
        .checked_add(n.checked_mul(4).ok_or(Error::Limit)?)
        .and_then(|t| t.checked_add(rim_triangles))
        .ok_or(Error::Limit)?;
    if expected_triangles > options.tessellation.triangle_limit {
        return Err(Error::Limit);
    }
    let mut cap_height_error = 0f64;
    for (end, ring) in rings.vertices.iter().enumerate() {
        for &id in ring {
            work.charge(1)?;
            cap_height_error = cap_height_error.max(
                (I::point(source_mesh.vertices[3 * id as usize + 2]) - I::point(rings.height[end]))
                    .magnitude(),
            );
        }
    }
    rings.distance = (I::point(rings.distance) + I::point(cap_height_error)).hi;
    let containment_distance = rings.unchanged_source_bound.max(rings.distance);
    let zero = I::point(0.);
    let step = I::point(std::f64::consts::TAU) / I::point(n as f64);
    let step = I::new(step.lo.next_down(), step.hi.next_up());
    let inscribed = I::point(radius_mm) * (step / I::point(2.)).cos();
    for p in source_mesh.vertices.as_chunks::<3>().0 {
        work.charge(16)?;
        // The literal triangles lie in the convex hull of their vertices and
        // the complete actual source lies within the inherited C0 tube. The
        // inscribed disc is inside the chord polygon, so this strict margin
        // keeps the complete actual bore inside the emitted stock.
        let r = (I::point(p[0]).square() + I::point(p[1]).square()).sqrt();
        if (r + I::point(containment_distance)).hi >= inscribed.lo {
            return Err(Error::Unsupported);
        }
    }
    let [lower, upper] = rings.height.map(I::point);
    let count = source_mesh.vertices.len() / 3 + 2 * n;
    let extra_triangles = expected_triangles + 8;
    work.storage(if in_place {
        count * 24 + extra_triangles * 16 + 52 * extra_triangles + 1024 * 1024
    } else {
        count * 24 + extra_triangles * 16 + count * 4 + 52 * extra_triangles + 1024 * 1024
    })?;
    if in_place {
        for &face in &source_mesh.face_indices {
            if face == rings.caps[0] || face == rings.caps[1] || Some(face) == wall {
                return Err(Error::Invalid);
            }
        }
    }
    let wall = wall.unwrap_or_else(|| {
        source_mesh
            .face_indices
            .iter()
            .copied()
            .max()
            .expect("qualified mesh has faces")
            + 1
    });
    let mut mesh = if in_place {
        let mut mesh = source.into_owned();
        mesh.vertices.reserve_exact(6 * n);
        let more = extra_triangles - mesh.face_indices.len();
        mesh.indices.reserve_exact(3 * more);
        mesh.face_indices.reserve_exact(more);
        mesh
    } else {
        let mut mesh = DisplayMesh {
            vertices: Vec::with_capacity(3 * count),
            indices: Vec::with_capacity(3 * extra_triangles),
            face_indices: Vec::with_capacity(extra_triangles),
        };
        mesh.vertices.extend_from_slice(&source_mesh.vertices);
        for (triangle, &face) in source_mesh
            .indices
            .as_chunks::<3>()
            .0
            .iter()
            .zip(&source_mesh.face_indices)
        {
            work.charge(16)?;
            if face == rings.caps[0] || face == rings.caps[1] {
                if fans {
                    continue;
                }
                return Err(Error::Invalid);
            }
            // A traced thread is oriented as material; the stock's outward
            // normal on its void is the reverse. A law bore already faces it.
            let mut ids = *triangle;
            if fans {
                ids.swap(1, 2);
            }
            mesh.indices.extend(ids);
            mesh.face_indices.push(face);
        }
        mesh
    };
    let base = mesh.vertices.len() / 3;
    let mut rounding = 0f64;
    let mut corners = Vec::with_capacity(n);
    for i in 0..n {
        work.charge(32)?;
        let angle = step * I::point(i as f64);
        corners.push([
            I::point(radius_mm) * angle.cos(),
            I::point(radius_mm) * angle.sin(),
        ]);
    }
    for z in [lower, upper] {
        for xy in &corners {
            work.charge(32)?;
            let p = [xy[0], xy[1], z];
            let literal = p.map(|v| v.lo + (v.hi - v.lo) * 0.5);
            if literal.iter().any(|v| !v.is_finite()) {
                return Err(Error::Unsupported);
            }
            let error = (0..3)
                .map(|k| (I::point(literal[k]) - p[k]).square())
                .fold(zero, |a, b| a + b)
                .sqrt()
                .hi;
            if !error.is_finite() {
                return Err(Error::Unsupported);
            }
            rounding = rounding.max(error);
            mesh.vertices.extend(literal);
        }
    }
    // A chord wall facet is within its sagitta of the exact cylinder. New cap
    // patches incur the actual rim/plane transfer; the union bound is the
    // maximum over the disjoint face families.
    let half = step / I::point(2.);
    let sagitta = (I::point(radius_mm) * (I::point(1.) - half.cos())).hi;
    let cap_distance = (I::point(rings.distance) + I::point(rounding)).hi;
    let wall_distance = (I::point(sagitta) + I::point(rounding)).hi;
    let distance = rings
        .unchanged_source_bound
        .max(cap_distance)
        .max(wall_distance);
    if distance > options.tessellation.deflection {
        return Err(Error::Unsupported);
    }
    // The wall facet normal across a chord is the radial direction at its
    // middle; the exact normal turns by at most half a step on either side,
    // which sleeve_segments kept within half the requested angle.
    for i in 0..n {
        let j = (i + 1) % n;
        let middle = step * (I::point(i as f64) + I::point(0.5));
        let normal = [middle.cos(), middle.sin(), zero];
        let a = (base + i) as u32;
        let b = (base + j) as u32;
        let c = b + n as u32;
        let d = a + n as u32;
        work.triangle(&mut mesh, [a, b, c], wall, normal)?;
        work.triangle(&mut mesh, [a, c, d], wall, normal)?;
    }
    let outer0 = (base as u32..(base + n) as u32).collect::<Vec<_>>();
    let outer1 = ((base + n) as u32..(base + 2 * n) as u32).collect::<Vec<_>>();
    plane(
        &mut mesh,
        &outer0,
        &rings.vertices[0],
        false,
        rings.caps[0],
        work,
    )?;
    plane(
        &mut mesh,
        &outer1,
        &rings.vertices[1],
        true,
        rings.caps[1],
        work,
    )?;
    if in_place {
        // No apex exists: every vertex is used. Certify the mesh as is.
        if mesh.face_indices.len() != expected_triangles {
            return Err(Error::Unresolved);
        }
        let resident = mesh_bytes(&mesh);
        let connectivity = manufacturing::certify_counted(&mesh, &mut |n, scratch| {
            if let Some(scratch) = scratch {
                work.storage(resident.checked_add(scratch).ok_or(Error::Limit)?)?;
            }
            work.charge(n)
        })
        .map_err(|error| match error {
            manufacturing::CountedError::Mesh(_) => Error::Unresolved,
            manufacturing::CountedError::Interrupted(error) => error,
        })?;
        return Ok((mesh, connectivity, distance));
    }
    // Remove former cap apices only by actual index use; no coordinate weld.
    work.storage(mesh_bytes(&mesh) * 2 + count * 4 + 52 * mesh.face_indices.len())?;
    let mut remap = vec![u32::MAX; count];
    let mut compact = DisplayMesh {
        vertices: Vec::with_capacity(mesh.vertices.len()),
        indices: Vec::with_capacity(mesh.indices.len()),
        face_indices: mesh.face_indices.clone(),
    };
    let old_mesh_bytes = mesh_bytes(&mesh);
    for id in mesh.indices {
        work.charge(1)?;
        let new = &mut remap[id as usize];
        if *new == u32::MAX {
            *new = (compact.vertices.len() / 3) as u32;
            compact
                .vertices
                .extend_from_slice(&mesh.vertices[3 * id as usize..3 * id as usize + 3]);
        }
        compact.indices.push(*new);
    }
    if compact.face_indices.len() != expected_triangles {
        return Err(Error::Unresolved);
    }
    let resident = old_mesh_bytes + mesh_bytes(&compact) + remap.capacity() * 4;
    let connectivity = manufacturing::certify_counted(&compact, &mut |n, scratch| {
        if let Some(scratch) = scratch {
            work.storage(resident.checked_add(scratch).ok_or(Error::Limit)?)?;
        }
        work.charge(n)
    })
    .map_err(|error| match error {
        manufacturing::CountedError::Mesh(_) => Error::Unresolved,
        manufacturing::CountedError::Interrupted(error) => error,
    })?;
    Ok((compact, connectivity, distance))
}
