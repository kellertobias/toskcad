//! Union of two retained boundaries that meet only across one shared planar
//! face from opposite sides, such as the two halves of a two-sided sweep.
//!
//! The shared face is removed from both boundaries and its edges are joined
//! once. Each operand must lie in the closed half-space behind its shared
//! face, so the operands' interiors are disjoint and the result keeps every
//! other face and name of both operands unchanged.
use super::*;
use multicad_core::references::{AnalyticGeometry, TopologyKind, TopologySnapshot};

const UNION: &str = "Rust plane-glued union";

fn refusal(message: &str) -> GeometryError {
    error(UNION, GeometryErrorKind::Unsupported, message)
}

fn distance(a: P, b: P) -> f64 {
    norm(sub(a, b))
}

/// A rotation taking `normal` to `+z`, as rows.
fn rotation(normal: P) -> Result<[P; 3], GeometryError> {
    let z = unit(normal)?;
    let seed = if z[0].abs() < 0.9 {
        [1., 0., 0.]
    } else {
        [0., 1., 0.]
    };
    let x = unit(cross(seed, z))?;
    Ok([x, cross(z, x), z])
}

/// The outward unit normal and a point of a planar face.
pub(super) fn plane(solid: &brep::Solid, face: usize) -> Option<(P, P)> {
    let face = &solid.boundary().faces[face];
    let Surface::Plane(frame) = &face.surface else {
        return None;
    };
    let sign = if face.same_sense { 1. } else { -1. };
    Some((frame.origin, scale(frame.z, sign)))
}

/// Points that identify an edge's geometry independently of direction.
fn probes(boundary: &Boundary, edge: usize) -> [P; 3] {
    let edge = &boundary.edges[edge];
    let [a, b] = edge.range;
    let at = |s: f64| edge.curve.evaluate(a + (b - a) * s).0;
    [at(0.), at(1. / 3.), at(2. / 3.)]
}

fn same_edge(a: &Boundary, ea: usize, b: &Boundary, eb: usize, tolerance: f64) -> Option<bool> {
    let [pa, qa, ra] = probes(a, ea);
    let [pb, qb, rb] = probes(b, eb);
    let close = |x: P, y: P| distance(x, y) <= tolerance;
    let ends = |e: &brep::Edge, boundary: &Boundary| e.vertices.map(|v| boundary.vertices[v]);
    let [a0, a1] = ends(&a.edges[ea], a);
    let [b0, b1] = ends(&b.edges[eb], b);
    if close(a0, b0) && close(a1, b1) && close(qa, qb) && close(ra, rb) && close(pa, pb) {
        return Some(true);
    }
    // Opposite direction: the probes at one third and two thirds swap.
    (close(a0, b1) && close(a1, b0) && close(qa, rb) && close(ra, qb)).then_some(false)
}

fn face_edges(boundary: &Boundary, face: usize) -> Vec<usize> {
    boundary.faces[face]
        .loops
        .iter()
        .flat_map(|trim| trim.uses.iter().map(|&(edge, _)| edge))
        .collect()
}

/// Match every edge of face `fa` of `a` to one edge of face `fb` of `b`:
/// the index in `b` and whether it runs the same way.
fn matching(
    a: &Boundary,
    fa: usize,
    b: &Boundary,
    fb: usize,
    tolerance: f64,
) -> Option<Vec<(usize, usize, bool)>> {
    let ea = face_edges(a, fa);
    let mut eb = face_edges(b, fb);
    if ea.len() != eb.len()
        || a.faces[fa].loops.len() != b.faces[fb].loops.len()
        || a.faces[fa].loops.iter().any(|trim| trim.vertex.is_some())
        || b.faces[fb].loops.iter().any(|trim| trim.vertex.is_some())
    {
        return None;
    }
    let mut pairs = Vec::with_capacity(ea.len());
    for edge in ea {
        let found = eb
            .iter()
            .position(|&other| same_edge(a, edge, b, other, tolerance).is_some())?;
        let other = eb.swap_remove(found);
        let forward = same_edge(a, edge, b, other, tolerance)?;
        pairs.push((edge, other, forward));
    }
    Some(pairs)
}

/// Whether `solid` lies in the closed half-space behind its face `face`, and
/// faces that do not border it stay clear of the plane by more than the
/// tessellation can hide.
pub(super) fn behind(
    solid: &brep::Solid,
    face: usize,
    tolerance: f64,
) -> Result<bool, GeometryError> {
    let (origin, normal) = plane(solid, face).ok_or_else(|| refusal("Glued face is not planar"))?;
    let rows = rotation(normal)?;
    let turn = |p: P| rows.map(|row| dot(row, p));
    let level = dot(rows[2], origin);
    let moved = solid.moved(&turn, &turn).map_err(failure)?;
    if moved.bounds_mm()[1][2] > level + tolerance {
        return Ok(false);
    }
    let boundary = solid.boundary();
    let glued = face_edges(boundary, face);
    let bordering: Vec<bool> = boundary
        .faces
        .iter()
        .enumerate()
        .map(|(index, other)| {
            index == face
                || other
                    .loops
                    .iter()
                    .flat_map(|trim| &trim.uses)
                    .any(|(edge, _)| glued.contains(edge))
        })
        .collect();
    let [low, high] = solid.bounds_mm();
    let deflection = (1e-4 * distance(low, high)).max(10. * tolerance);
    let mesh = brep::mesh::tessellate(
        solid,
        brep::mesh::Options {
            deflection_mm: deflection,
            angle_rad: 0.5,
            triangle_limit: 2_000_000,
        },
        None,
    )
    .map_err(failure)?;
    let clear = mesh
        .triangles
        .iter()
        .zip(&mesh.faces)
        .all(|(triangle, &f)| {
            bordering[f as usize]
                || triangle.iter().all(|&v| {
                    level - dot(rows[2], mesh.positions[v as usize]) > deflection + tolerance
                })
        });
    Ok(clear)
}

fn names(snapshot: &TopologySnapshot, face: usize) -> FaceName {
    let element = snapshot
        .elements
        .iter()
        .find(|e| e.key.kind == TopologyKind::Face && e.key.ordinal as usize == face);
    let Some(element) = element else {
        return FaceName {
            origins: vec![],
            support: None,
            anchor: None,
        };
    };
    let anchored = element.frames.first();
    // The anchored origin stays first so that the frame keeps its owner.
    let mut origins: Vec<_> = anchored.map(|f| f.origin.clone()).into_iter().collect();
    origins.extend(
        element
            .origins
            .iter()
            .filter(|o| Some(*o) != anchored.map(|f| &f.origin))
            .cloned(),
    );
    FaceName {
        origins,
        support: match element.geometry {
            Some(AnalyticGeometry::Plane { frame }) => Some(frame),
            _ => None,
        },
        anchor: anchored.map(|f| f.frame),
    }
}

impl Solid {
    /// A union of two retained boundaries joined across one shared planar
    /// face, or `None` when either operand keeps another representation.
    pub(in crate::geometry) fn glued_brep_union(
        &self,
        other: &Self,
        operation: BooleanOperation,
        policy: GeometryPolicy,
    ) -> Option<Result<Option<Self>, GeometryError>> {
        if operation != BooleanOperation::Union {
            return None;
        }
        let (a, b) = (self.imported_brep_body()?, other.imported_brep_body()?);
        Some(self.glue(other, a, b, policy).map(Some))
    }

    fn glue(
        &self,
        other: &Self,
        a: &ImportedBrepBody,
        b: &ImportedBrepBody,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        if policy != self.policy || policy != other.policy {
            return Err(refusal(
                "A glued union requires both operands' geometry policy",
            ));
        }
        let tolerance = policy.linear_tolerance_mm;
        let (sa, sb) = (&a.solid, &b.solid);
        let (ba, bb) = (sa.boundary(), sb.boundary());
        let mut candidates = Vec::new();
        for fa in 0..ba.faces.len() {
            let Some((oa, na)) = plane(sa, fa) else {
                continue;
            };
            for fb in 0..bb.faces.len() {
                let Some((ob, nb)) = plane(sb, fb) else {
                    continue;
                };
                if dot(na, nb) < -1. + 1e-12
                    && dot(sub(ob, oa), na).abs() <= tolerance
                    && let Some(pairs) = matching(ba, fa, bb, fb, tolerance)
                {
                    candidates.push((fa, fb, pairs));
                }
            }
        }
        let [(fa, fb, pairs)] = <[_; 1]>::try_from(candidates).map_err(|_| {
            refusal("A glued union requires exactly one shared planar face of equal boundary")
        })?;
        if !behind(sa, fa, tolerance)? || !behind(sb, fb, tolerance)? {
            return Err(refusal(
                "Glued operands must lie on opposite sides of their shared face",
            ));
        }

        // B's vertices and edges on the shared face become A's.
        let mut vertex_map = vec![usize::MAX; bb.vertices.len()];
        let mut edge_map = vec![None; bb.edges.len()];
        for &(ea, eb, forward) in &pairs {
            edge_map[eb] = Some((ea, forward));
            let [a0, a1] = ba.edges[ea].vertices;
            let [b0, b1] = bb.edges[eb].vertices;
            let [m0, m1] = if forward { [a0, a1] } else { [a1, a0] };
            for (vb, va) in [(b0, m0), (b1, m1)] {
                if vertex_map[vb] != usize::MAX && vertex_map[vb] != va {
                    return Err(refusal("Shared face vertices do not correspond"));
                }
                vertex_map[vb] = va;
            }
        }
        let mut merged = Boundary {
            vertices: ba.vertices.clone(),
            edges: ba.edges.clone(),
            faces: Vec::new(),
            shells: Vec::new(),
        };
        for (index, &point) in bb.vertices.iter().enumerate() {
            if vertex_map[index] == usize::MAX {
                merged.vertices.push(point);
                vertex_map[index] = merged.vertices.len() - 1;
            }
        }
        let mut edges = vec![(0, true); bb.edges.len()];
        for (index, edge) in bb.edges.iter().enumerate() {
            edges[index] = match edge_map[index] {
                Some(mapped) => mapped,
                None => {
                    merged.edges.push(brep::Edge {
                        vertices: edge.vertices.map(|v| vertex_map[v]),
                        ..edge.clone()
                    });
                    (merged.edges.len() - 1, true)
                }
            };
        }
        let mut face_map = [vec![None; ba.faces.len()], vec![None; bb.faces.len()]];
        let snapshots = [
            self.topology_snapshot(EntityId::new(), 0)?,
            other.topology_snapshot(EntityId::new(), 0)?,
        ];
        let mut naming = Vec::new();
        for (side, (boundary, skip)) in [(ba, fa), (bb, fb)].into_iter().enumerate() {
            for (index, face) in boundary.faces.iter().enumerate() {
                if index == skip {
                    continue;
                }
                let mut face = face.clone();
                if side == 1 {
                    for trim in &mut face.loops {
                        for used in &mut trim.uses {
                            let (edge, forward) = edges[used.0];
                            *used = (edge, used.1 == forward);
                        }
                        trim.vertex = trim.vertex.map(|v| vertex_map[v]);
                    }
                }
                merged.faces.push(face);
                face_map[side][index] = Some(merged.faces.len() - 1);
                naming.push(names(&snapshots[side], index));
            }
        }
        // The exterior shells join; voids keep their own shells.
        let mut exterior = Vec::new();
        let mut voids = Vec::new();
        for (side, boundary) in [ba, bb].into_iter().enumerate() {
            for (index, shell) in boundary.shells.iter().enumerate() {
                let faces: Vec<_> = shell
                    .faces
                    .iter()
                    .filter_map(|&face| face_map[side][face])
                    .collect();
                if index == 0 {
                    if shell.reversed {
                        return Err(refusal("Glued operands need outward exterior shells"));
                    }
                    exterior.extend(faces);
                } else {
                    voids.push(Shell {
                        faces,
                        reversed: shell.reversed,
                    });
                }
            }
        }
        merged.shells.push(Shell {
            faces: exterior,
            reversed: false,
        });
        merged.shells.extend(voids);
        let solid = brep::Solid::new(merged, tolerance, None).map_err(failure)?;
        let expected = sa.volume_mm3() + sb.volume_mm3();
        let slack = sa.volume_error_mm3() + sb.volume_error_mm3() + solid.volume_error_mm3();
        if (solid.volume_mm3() - expected).abs() > slack.max(1e-9 * expected) {
            return Err(error(
                UNION,
                GeometryErrorKind::InvalidResult,
                "Glued union volume differs from its operands' sum",
            ));
        }
        Solid::from_authored_brep(solid, a.feature, AuthoredNaming { faces: naming }, policy)
    }
}
