//! Selected literal analytic boundary evidence. This stage proves finite schema
//! and manifold incidence, but grants no occupied material or Native admission.
use super::*;
use owned_insert::CompoundDocument;
mod chains;
#[cfg(test)]
mod tests;
pub use chains::{BoreChain, Chains};

#[derive(Clone, Copy)]
pub struct Options {
    pub maximum_work: usize,
    /// Additional query storage; the borrowed document has its own ledger.
    pub maximum_bytes: usize,
}
impl Default for Options {
    fn default() -> Self {
        Self {
            maximum_work: 16_000_000,
            maximum_bytes: 64 * 1024 * 1024,
        }
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Error {
    Invalid,
    Unsupported,
    Limit,
    Cancelled,
}
impl From<ReadError> for Error {
    fn from(value: ReadError) -> Self {
        match value {
            ReadError::Invalid => Self::Invalid,
            ReadError::Unsupported => Self::Unsupported,
            ReadError::Limit => Self::Limit,
        }
    }
}
#[derive(Clone, Copy)]
pub struct LiteralNumber<'input> {
    entity: u32,
    ordinal: usize,
    spelling: &'input str,
}
impl<'input> LiteralNumber<'input> {
    pub fn entity(&self) -> u32 {
        self.entity
    }
    pub fn ordinal(&self) -> usize {
        self.ordinal
    }
    pub fn spelling(&self) -> &'input str {
        self.spelling
    }
}
/// Original owners and scalar spellings stay live with the neutral document.
/// No source recipe, binary64 equality, or header supplies their identity.
pub struct Boundary<'document, 'input> {
    document: &'document CompoundDocument<'input>,
    item: u32,
    solid: u32,
    dependencies: Vec<u32>,
    face_entities: Vec<u32>,
    numbers: Vec<LiteralNumber<'input>>,
    shell: Shell<'document, 'input>,
    faces: Vec<Face>,
    work: usize,
    bytes: usize,
}
impl<'document, 'input> Boundary<'document, 'input> {
    pub fn document(&self) -> &'document CompoundDocument<'input> {
        self.document
    }
    pub fn item_entity(&self) -> u32 {
        self.item
    }
    pub fn solid_entity(&self) -> u32 {
        self.solid
    }
    pub fn dependencies(&self) -> &[u32] {
        &self.dependencies
    }
    pub fn face_entities(&self) -> &[u32] {
        &self.face_entities
    }
    pub fn numeric_leaves(&self) -> &[LiteralNumber<'input>] {
        &self.numbers
    }
    pub fn edge_count(&self) -> usize {
        self.shell.edges.len()
    }
    pub fn vertex_count(&self) -> usize {
        self.shell.points.len()
    }
    pub fn planar_face_count(&self) -> usize {
        self.faces
            .iter()
            .filter(|f| matches!(f, Face::Planar { .. }))
            .count()
    }
    pub fn conical_face_count(&self) -> usize {
        self.faces
            .iter()
            .filter(|f| matches!(f, Face::Conical(_)))
            .count()
    }
    pub fn cylindrical_face_count(&self) -> usize {
        self.faces
            .iter()
            .filter(|f| matches!(f, Face::Cylindrical { .. }))
            .count()
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    /// Conservative query allocation bound, including traversal and Shell
    /// scratch. This excludes the separately retained neutral document.
    pub fn query_budget_bytes(&self) -> usize {
        self.bytes
    }
}
struct Budget<'a> {
    work: usize,
    maximum: usize,
    interrupt: &'a dyn Fn() -> bool,
}
impl Budget<'_> {
    fn charge(&mut self, n: usize) -> std::result::Result<(), Error> {
        if (self.interrupt)() {
            return Err(Error::Cancelled);
        }
        self.work = self.work.checked_add(n).ok_or(Error::Limit)?;
        if self.work > self.maximum {
            return Err(Error::Limit);
        }
        Ok(())
    }
}
fn references(
    values: &[Value<'_>],
    pending: &mut Vec<u32>,
    budget: &mut Budget<'_>,
) -> std::result::Result<(), Error> {
    for value in values {
        budget.charge(1)?;
        match value {
            Value::Reference(id) => {
                if pending.len() == 16_384 {
                    return Err(Error::Limit);
                }
                pending.push(*id);
            }
            Value::List(values) => references(values, pending, budget)?,
            Value::Typed(_) => return Err(Error::Unsupported),
            _ => {}
        }
    }
    Ok(())
}
fn supported(name: &str) -> bool {
    matches!(
        name,
        "MANIFOLD_SOLID_BREP"
            | "CLOSED_SHELL"
            | "ADVANCED_FACE"
            | "FACE_OUTER_BOUND"
            | "FACE_BOUND"
            | "EDGE_LOOP"
            | "ORIENTED_EDGE"
            | "EDGE_CURVE"
            | "VERTEX_POINT"
            | "CARTESIAN_POINT"
            | "AXIS2_PLACEMENT_3D"
            | "DIRECTION"
            | "LINE"
            | "VECTOR"
            | "PLANE"
            | "CIRCLE"
            | "CYLINDRICAL_SURFACE"
            | "CONICAL_SURFACE"
    )
}
/// Decode the actual local selected member. Mapped occurrences and other
/// schemas remain separate; none is normalized into a supported ordinary item.
pub fn decode<'document, 'input>(
    document: &'document CompoundDocument<'input>,
    item: u32,
    options: Options,
    interrupt: &dyn Fn() -> bool,
) -> std::result::Result<Boundary<'document, 'input>, Error> {
    let mut budget = Budget {
        work: 0,
        maximum: options.maximum_work,
        interrupt,
    };
    budget.charge(0)?;
    if options.maximum_work == 0
        || options.maximum_work > 16_000_000
        || options.maximum_bytes == 0
        || options.maximum_bytes > 64 * 1024 * 1024
    {
        return Err(Error::Invalid);
    }
    // This small analytic schema has at most64 faces,256 edges,4096 original
    // dependency records and16384 numeric leaves. Eight MiB covers their live
    // pinned BTree node minima, Vec doubling capacities and Shell scratch;
    // document input/AST are borrowed, not copied into this query budget.
    const BYTES: usize = 8 * 1024 * 1024;
    if options.maximum_bytes < BYTES {
        return Err(Error::Limit);
    }
    budget.charge(document.represented_members().len())?;
    let member = document
        .represented_members()
        .iter()
        .find(|m| m.item_entity() == item)
        .ok_or(Error::Invalid)?;
    if member.mapping_link().is_some() {
        return Err(Error::Unsupported);
    }
    let solid = member.local_solid_entity();
    let records = document.records();
    let mut pending = vec![solid];
    let mut selected = BTreeSet::new();
    let mut numbers = Vec::new();
    while let Some(id) = pending.pop() {
        budget.charge(64)?;
        if !selected.insert(id) {
            continue;
        }
        if selected.len() > 4096 {
            return Err(Error::Limit);
        }
        let [entity] = records.get(&id).ok_or(Error::Invalid)?.as_slice() else {
            return Err(Error::Unsupported);
        };
        if !supported(entity.name) {
            return Err(Error::Unsupported);
        }
        let source = entity.source.ok_or(Error::Invalid)?;
        budget.charge(source.len())?;
        // Bound the tokenizer's temporary Vec capacity before it runs.
        if source.len() > 65_536 {
            return Err(Error::Limit);
        }
        let leaves = part21::numeric_spellings(entity).ok_or(Error::Invalid)?;
        if leaves.len() > 16_384 - numbers.len() {
            return Err(Error::Limit);
        }
        numbers.extend(
            leaves
                .into_iter()
                .enumerate()
                .map(|(ordinal, spelling)| LiteralNumber {
                    entity: id,
                    ordinal,
                    spelling,
                }),
        );
        references(&entity.args, &mut pending, &mut budget)?;
    }
    let mut shell = Shell::new(records, 1e-7, true);
    let mut context = None;
    for components in records.values() {
        budget.charge(components.len() + 1)?;
        if let [entity] = components.as_slice()
            && entity.name == "ADVANCED_BREP_SHAPE_REPRESENTATION"
        {
            let [_, _, value] = entity.args.as_slice() else {
                return Err(Error::Invalid);
            };
            if context.replace(reference(value)?).is_some() {
                return Err(Error::Unsupported);
            }
        }
    }
    // check_units scans component names once more to detect angular supports.
    budget.charge(records.values().map(|c| c.len() + 1).sum())?;
    let context = context.ok_or(Error::Invalid)?;
    for component in shell.components(context)? {
        if component.name == "GLOBAL_UNIT_ASSIGNED_CONTEXT" {
            let [units] = component.args.as_slice() else {
                return Err(Error::Invalid);
            };
            if list(units)?.len() > 16 {
                return Err(Error::Limit);
            }
        }
    }
    shell.check_units(context)?;
    let face_entities = shell.member_faces(solid)?;
    if face_entities.len() > 64 {
        return Err(Error::Limit);
    }
    let edge_count = selected
        .iter()
        .filter(|id| records[id][0].name == "EDGE_CURVE")
        .count();
    if edge_count > 256 {
        return Err(Error::Limit);
    }
    let mut bounds = BTreeSet::new();
    let mut loops = BTreeSet::new();
    let mut coedges = 0usize;
    for id in &face_entities {
        budget.charge(64)?;
        let [_, face_bounds, _, _] = shell.args(*id, "ADVANCED_FACE")? else {
            return Err(Error::Invalid);
        };
        let face_bounds = list(face_bounds)?;
        if face_bounds.len() > 512 {
            return Err(Error::Limit);
        }
        for value in face_bounds {
            let bound = reference(value)?;
            if !bounds.insert(bound) {
                return Err(Error::Invalid);
            }
            let entity = shell.simple(bound)?;
            if !matches!(entity.name, "FACE_BOUND" | "FACE_OUTER_BOUND") {
                return Err(Error::Invalid);
            }
            let [_, wire, _] = entity.args.as_slice() else {
                return Err(Error::Invalid);
            };
            let wire = reference(wire)?;
            if !loops.insert(wire) {
                return Err(Error::Invalid);
            }
            let [_, uses] = shell.args(wire, "EDGE_LOOP")? else {
                return Err(Error::Invalid);
            };
            coedges = coedges.checked_add(list(uses)?.len()).ok_or(Error::Limit)?;
            if coedges > 512 {
                return Err(Error::Limit);
            }
        }
    }
    // Buy the finite support/trim/incidence bound once on the same ledger.
    // Sum of all wire pair counts is <=coedges^2. The128 factor covers each
    // pair's scalar geometry and lookup predicates; linear terms cover curve
    // normalization, traversal, face fans and their pinned BTree operations.
    let support_work = 128 * coedges * coedges + 4096 * coedges + 8192 * face_entities.len();
    budget.charge(support_work)?;
    let faces = shell.validate_faces(face_entities.clone())?;
    budget.charge(0)?;
    if faces.iter().any(|f| {
        !matches!(
            f,
            Face::Planar { .. } | Face::Cylindrical { .. } | Face::Conical(_)
        )
    }) {
        return Err(Error::Unsupported);
    }
    Ok(Boundary {
        document,
        item,
        solid,
        dependencies: selected.into_iter().collect(),
        face_entities,
        numbers,
        shell,
        faces,
        work: budget.work,
        bytes: BYTES,
    })
}
