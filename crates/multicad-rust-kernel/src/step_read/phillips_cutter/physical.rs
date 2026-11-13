//! Literal finite cavity parsing only. This is not a geometric support-incidence,
//! mass, self-intersection, Phillips gauge or native-publication certificate.
use super::*;

mod branches;
mod endpoints;
pub use branches::{
    EdgeSection as IntersectionSection, EdgeTube as IntersectionEdgeTube,
    Options as IntersectionTubeOptions, Tubes as IntersectionTubes,
};
mod incidence;
pub use endpoints::{EndpointCouplingBounds, EndpointCouplingOptions};
mod literal;
mod supports;
pub use supports::{
    CubicSupportBounds, CubicSupportOptions, EdgeSupportBound, ImplicitSupportEnclosure,
};
#[cfg(test)]
mod tests;
pub use literal::{Axis, Curve, Surface};

/// The actual head face remains part of the material shell. Only this explicitly
/// identified inner bound is reversed into a virtual cavity closure.
#[derive(Clone, Debug)]
pub struct Mouth {
    pub material_face: u32,
    pub material_surface: u32,
    pub material_bound: u32,
    pub material_forward: bool,
    pub material_uses: Vec<(u32, bool)>,
    pub material_loop: MaterialLoop,
}
#[derive(Clone, Debug)]
pub struct MaterialUse {
    pub entity: u32,
    pub edge_entity: u32,
    pub forward: bool,
}
#[derive(Clone, Debug)]
pub struct MaterialLoop {
    pub bound_entity: u32,
    pub loop_entity: u32,
    pub outside: bool,
    pub forward: bool,
    pub uses: Vec<MaterialUse>,
}
#[derive(Clone, Debug)]
pub struct Vertex {
    pub entity: u32,
    pub point_entity: u32,
    pub point: Point,
}
#[derive(Clone, Debug)]
pub struct Edge {
    pub entity: u32,
    pub vertices: [usize; 2],
    pub geometry: u32,
    pub same_sense: bool,
    pub owners: [usize; 2],
    pub curve: Curve,
}
#[derive(Clone, Debug)]
pub struct Face {
    /// None denotes only the virtual Mouth. Every material side face retains
    /// its actual entity and original sense independently of the cavity view.
    pub material_entity: Option<u32>,
    pub material_forward: bool,
    pub material_loop: Option<MaterialLoop>,
    pub cavity_forward: bool,
    pub surface_entity: u32,
    pub support: Surface,
    pub uses: Vec<(usize, bool)>,
}
/// Immutable parsed physical data. Private fields prevent this structural
/// receipt from being substituted for a qualified geometric body.
#[derive(Clone, Debug)]
pub struct Boundary {
    mouth: Mouth,
    vertices: Vec<Vertex>,
    edges: Vec<Edge>,
    faces: Vec<Face>,
}
impl Boundary {
    pub fn mouth(&self) -> &Mouth {
        &self.mouth
    }
    pub fn vertices(&self) -> &[Vertex] {
        &self.vertices
    }
    pub fn edges(&self) -> &[Edge] {
        &self.edges
    }
    pub fn faces(&self) -> &[Face] {
        &self.faces
    }
}

/// Parse the literal connected finite cavity component of one material screw.
/// The returned data has structural provenance only. No source recipe, support
/// intersection, volume, watertight geometry or manufacturing claim is made.
pub fn read_phillips_screw_boundary(
    input: &str,
    tolerance_mm: f64,
    work_limit: usize,
) -> std::result::Result<Boundary, PhillipsReadError> {
    read_phillips_screw_boundary_with_interrupt(input, tolerance_mm, work_limit, None)
}
pub fn read_phillips_screw_boundary_with_interrupt(
    input: &str,
    tolerance_mm: f64,
    work_limit: usize,
    interrupt: Option<&dyn Fn() -> bool>,
) -> std::result::Result<Boundary, PhillipsReadError> {
    let mut work = Work {
        remaining: work_limit,
        interrupt,
        cancelled: Cell::new(false),
    };
    let result = compound::read_members_by(input, tolerance_mm, true, true, |shell, id, _| {
        work.step(1)?;
        let faces = shell.member_faces(id)?;
        read_member(shell, &faces, &mut work)
    });
    if work.cancelled.get() {
        return Err(PhillipsReadError::Cancelled);
    }
    let mut members = result.map_err(PhillipsReadError::from)?;
    members.pop().ok_or(PhillipsReadError::Invalid)
}

struct MaterialFace {
    entity: u32,
    surface: u32,
    forward: bool,
    bounds: Vec<u32>,
    uses: Vec<(u32, bool)>,
}
fn material_loop(shell: &Shell<'_, '_>, bound: u32) -> Result<MaterialLoop> {
    let (outside, wire, forward) = shell.bound(bound, false)?;
    let [Value::String(_), uses] = shell.args(wire, "EDGE_LOOP")? else {
        return Err(ReadError::Invalid);
    };
    let uses = refs(uses)?
        .into_iter()
        .map(|entity| {
            let [
                Value::String(_),
                Value::Derived,
                Value::Derived,
                edge,
                sense,
            ] = shell.args(entity, "ORIENTED_EDGE")?
            else {
                return Err(ReadError::Invalid);
            };
            Ok(MaterialUse {
                entity,
                edge_entity: reference(edge)?,
                forward: crate::step_read::sense(sense)?,
            })
        })
        .collect::<Result<Vec<_>>>()?;
    Ok(MaterialLoop {
        bound_entity: bound,
        loop_entity: wire,
        outside,
        forward,
        uses,
    })
}
fn material_faces(
    shell: &mut Shell<'_, '_>,
    ids: &[u32],
    work: &mut Work<'_>,
) -> Result<Vec<MaterialFace>> {
    // Check the original material graph BEFORE any virtual-cap decomposition.
    helical::profiled::closed_incidence(shell, ids)?;
    ids.iter()
        .map(|&entity| {
            work.step(1)?;
            let (surface, forward, bounds) = helical::face(shell, entity)?;
            let mut uses = Vec::new();
            for &bound in &bounds {
                let (outside, _, _) = shell.bound(bound, false)?;
                uses.extend(helical::ring(shell, bound, outside)?);
            }
            Ok(MaterialFace {
                entity,
                surface,
                forward,
                bounds,
                uses,
            })
        })
        .collect()
}
fn owners(faces: &[MaterialFace]) -> Result<BTreeMap<u32, [usize; 2]>> {
    let mut map = BTreeMap::<u32, Vec<usize>>::new();
    for (i, face) in faces.iter().enumerate() {
        for &(edge, _) in &face.uses {
            map.entry(edge).or_default().push(i);
        }
    }
    map.into_iter()
        .map(|(edge, owners)| {
            let [a, b] = owners.as_slice() else {
                return Err(ReadError::Invalid);
            };
            Ok((edge, [*a, *b]))
        })
        .collect()
}
fn component(
    faces: &[MaterialFace],
    owners: &BTreeMap<u32, [usize; 2]>,
    head: usize,
    hole: &[(u32, bool)],
    work: &mut Work<'_>,
) -> Result<Option<Vec<usize>>> {
    let cut: BTreeSet<_> = hole.iter().map(|u| u.0).collect();
    let mut pending = Vec::new();
    for &edge in &cut {
        let pair = owners[&edge];
        let Some(other) = pair.into_iter().find(|&i| i != head) else {
            return Err(ReadError::Invalid);
        };
        pending.push(other);
    }
    let mut found = BTreeSet::new();
    while let Some(i) = pending.pop() {
        work.step(1)?;
        if i == head {
            return Ok(None);
        }
        if !found.insert(i) {
            continue;
        }
        if found.len() > 28 {
            return Ok(None);
        }
        for &(edge, _) in &faces[i].uses {
            if !cut.contains(&edge) {
                pending.extend(owners[&edge].into_iter().filter(|&j| j != i));
            }
        }
    }
    Ok((found.len() == 28).then(|| found.into_iter().collect()))
}
pub(crate) fn read_member(
    shell: &mut Shell<'_, '_>,
    ids: &[u32],
    work: &mut Work<'_>,
) -> Result<Boundary> {
    let material = material_faces(shell, ids, work)?;
    let owners = owners(&material)?;
    let mut candidates = Vec::new();
    for (index, face) in material.iter().enumerate() {
        work.step(1)?;
        if face.bounds.len() != 2 || !face.forward || shell.simple(face.surface)?.name != "PLANE" {
            continue;
        }
        let hole = helical::ring(shell, face.bounds[1], false)?;
        if let Some(component) = component(&material, &owners, index, &hole, work)? {
            candidates.push((index, hole, component));
        }
    }
    let [(index, hole, component)] = candidates.as_slice() else {
        return Err(ReadError::Unsupported);
    };
    let head = &material[*index];
    let mouth = Mouth {
        material_face: head.entity,
        material_surface: head.surface,
        material_bound: head.bounds[1],
        material_forward: head.forward,
        material_uses: hole.clone(),
        material_loop: material_loop(shell, head.bounds[1])?,
    };
    let mut views = vec![SourceFaceView {
        surface: head.surface,
        forward: head.forward,
        uses: hole.iter().rev().map(|&(id, sense)| (id, !sense)).collect(),
    }];
    let mut entities = vec![None];
    for &i in component {
        let face = &material[i];
        require(face.bounds.len() == 1)?;
        let uses = helical::ring(shell, face.bounds[0], true)?;
        views.push(SourceFaceView {
            surface: face.surface,
            forward: !face.forward,
            uses: uses
                .into_iter()
                .rev()
                .map(|(id, sense)| (id, !sense))
                .collect(),
        });
        entities.push(Some((
            face.entity,
            face.forward,
            material_loop(shell, face.bounds[0])?,
        )));
    }
    incidence::read(shell, mouth, &views, &entities, work)
}
