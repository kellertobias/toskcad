//! Product structure: which shape representations a product uses and where
//! each occurrence of a component sits.
//!
//! Representations joined by a relationship without a transformation
//! describe one shape and form one node. A node is placed in its parent by
//! a next-assembly-usage occurrence through its context-dependent shape
//! representation, or by a mapped item. A node used several times is read
//! once; every occurrence places that one definition.
//!
//! Only rigid, orientation-preserving placements are admitted. A mirroring
//! or scaling operator, a placement between items outside the related
//! representations, a cycle, an occurrence without a placement and any
//! solid that no occurrence reaches are refused rather than read otherwise.
use super::{Source, check, cross, normalized};
use crate::brep::Point;
use crate::step_read::part21::{Entity, Value};
use crate::step_read::{ReadError, Result, reference};
use std::collections::{BTreeMap, BTreeSet};

/// The deepest nesting of occurrences below a root.
pub const MAX_DEPTH: usize = 32;

/// A rigid, orientation-preserving placement: a point `p` of the placed
/// shape lies at `origin + axes[0] p[0] + axes[1] p[1] + axes[2] p[2]`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Placement {
    pub origin: Point,
    pub axes: [Point; 3],
}

impl Placement {
    pub const IDENTITY: Self = Self {
        origin: [0.; 3],
        axes: [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
    };

    pub fn vector(&self, v: Point) -> Point {
        std::array::from_fn(|i| {
            self.axes[0][i] * v[0] + self.axes[1][i] * v[1] + self.axes[2][i] * v[2]
        })
    }

    pub fn point(&self, p: Point) -> Point {
        let v = self.vector(p);
        std::array::from_fn(|i| v[i] + self.origin[i])
    }

    /// The placement, in this placement's parent, of a shape placed by
    /// `inner` in this placement's coordinates.
    pub fn then(&self, inner: &Self) -> Self {
        Self {
            origin: self.point(inner.origin),
            axes: inner.axes.map(|axis| self.vector(axis)),
        }
    }

    pub fn inverse(&self) -> Self {
        let axes = std::array::from_fn(|i| std::array::from_fn(|j| self.axes[j][i]));
        let turned = Self {
            origin: [0.; 3],
            axes,
        };
        Self {
            origin: turned.vector(self.origin).map(|x| -x),
            axes,
        }
    }

    /// Largest coordinate difference of the axes and of the origins.
    fn distance(&self, other: &Self) -> [f64; 2] {
        let largest = |a: &[f64], b: &[f64]| {
            a.iter()
                .zip(b)
                .map(|(x, y)| (x - y).abs())
                .fold(0., f64::max)
        };
        [
            largest(self.axes.as_flattened(), other.axes.as_flattened()),
            largest(&self.origin, &other.origin),
        ]
    }
}

/// One use of a node inside another.
#[derive(Clone, Debug, PartialEq)]
pub struct Instance {
    pub node: usize,
    /// The child's coordinates placed in the parent's.
    pub placement: Placement,
    /// The occurrence name, or its identifier when it has no name.
    pub name: String,
}

/// One shape definition: its own solids and the definitions it places.
#[derive(Clone, Debug, PartialEq)]
pub struct Node<S> {
    /// The product name, or empty for a shape that no product defines.
    pub name: String,
    pub solids: Vec<S>,
    pub children: Vec<Instance>,
}

/// One placed use of a node's solid.
#[derive(Clone, Debug, PartialEq)]
pub struct Occurrence {
    pub node: usize,
    pub solid: usize,
    pub placement: Placement,
    /// The index of the root and of the instance chosen at every level.
    pub path: Vec<usize>,
}

/// Definitions and their placed uses. Every node is reached from a root
/// and no node contains itself.
#[derive(Clone, Debug, PartialEq)]
pub struct Assembly<S> {
    pub nodes: Vec<Node<S>>,
    pub roots: Vec<usize>,
}

impl<S> Assembly<S> {
    /// Every solid occurrence, depth first: a node's own solids, then its
    /// children in order.
    pub fn occurrences(&self) -> Vec<Occurrence> {
        let mut found = Vec::new();
        for (index, &root) in self.roots.iter().enumerate() {
            self.visit(root, Placement::IDENTITY, &mut vec![index], &mut found);
        }
        found
    }

    fn visit(
        &self,
        node: usize,
        placement: Placement,
        path: &mut Vec<usize>,
        found: &mut Vec<Occurrence>,
    ) {
        let entry = &self.nodes[node];
        for solid in 0..entry.solids.len() {
            found.push(Occurrence {
                node,
                solid,
                placement,
                path: path.clone(),
            });
        }
        for (index, child) in entry.children.iter().enumerate() {
            path.push(index);
            self.visit(child.node, placement.then(&child.placement), path, found);
            path.pop();
        }
    }

    pub(super) fn try_map<T, E>(
        self,
        mut f: impl FnMut(S) -> std::result::Result<T, E>,
    ) -> std::result::Result<Assembly<T>, E> {
        let mut nodes = Vec::with_capacity(self.nodes.len());
        for node in self.nodes {
            nodes.push(Node {
                name: node.name,
                solids: node
                    .solids
                    .into_iter()
                    .map(&mut f)
                    .collect::<std::result::Result<_, E>>()?,
                children: node.children,
            });
        }
        Ok(Assembly {
            nodes,
            roots: self.roots,
        })
    }
}

/// The text of a string attribute, with doubled apostrophes undone.
fn text(value: &Value<'_>) -> Result<String> {
    match value {
        Value::String(text) => Ok(text.replace("''", "'")),
        _ => Err(ReadError::Invalid),
    }
}

/// Exactly `N` attributes.
fn arity<'a, 'i, const N: usize>(args: &'a [Value<'i>]) -> Result<&'a [Value<'i>; N]> {
    args.try_into().ok().ok_or(ReadError::Invalid)
}

/// The attributes of `id` when it is a simple `name` entity, or `None` for
/// any other entity. A missing entity is invalid.
fn kind<'s, 'i>(source: &'s Source<'i>, id: u32, name: &str) -> Result<Option<&'s [Value<'i>]>> {
    match source.args(id, name) {
        Ok(args) => Ok(Some(args)),
        Err(ReadError::Invalid) => Err(ReadError::Invalid),
        Err(_) => Ok(None),
    }
}

fn pair(args: &[Value<'_>]) -> Result<[u32; 2]> {
    let [_, _, first, second] = arity(args)?;
    Ok([reference(first)?, reference(second)?])
}

/// Union of representations related without a transformation; the group
/// of a representation is the smallest id it is joined with.
struct Groups(BTreeMap<u32, u32>);

impl Groups {
    fn of(&mut self, id: u32) -> u32 {
        let mut id = id;
        loop {
            let up = self.0[&id];
            if up == id {
                return id;
            }
            let above = self.0[&up];
            self.0.insert(id, above);
            id = up;
        }
    }

    fn join(&mut self, a: u32, b: u32) {
        let (a, b) = (self.of(a), self.of(b));
        self.0.insert(a.max(b), a.min(b));
    }
}

/// A representation: its items and its context.
struct Representation<'a, 'i> {
    items: &'a [Value<'i>],
    context: u32,
}

impl Source<'_> {
    /// The placement an item or operator states, as a map from its own
    /// coordinates to its representation's.
    fn placement(&self, id: u32) -> Result<Placement> {
        if let Some(placement) = self.complex_operator(id)? {
            return Ok(placement);
        }
        let entity = self.simple(id)?;
        match entity.name {
            "AXIS2_PLACEMENT_3D" => {
                let frame = self.frame(id)?;
                Ok(Placement {
                    origin: frame.origin,
                    axes: [frame.x, frame.y, frame.z],
                })
            }
            "CARTESIAN_TRANSFORMATION_OPERATOR_3D" => self.operator(&entity.args),
            _ => Err(ReadError::Unsupported),
        }
    }

    /// ISO 10303-42 `cartesian_transformation_operator_3d`: the axes are
    /// derived by `base_axis`, which may yield a left-handed set. Mirroring
    /// and scaling are refused.
    fn operator(&self, args: &[Value<'_>]) -> Result<Placement> {
        let [_, _, _, axis1, axis2, origin, scale, axis3] = arity(args)?;
        self.operator_attributes([axis1, axis2, origin, scale, axis3])
    }

    /// The same operator as an ISO 10303-21 complex instance, its attributes
    /// split over its supertypes in alphabetical order:
    /// `(CARTESIAN_TRANSFORMATION_OPERATOR(axis1, axis2, origin, scale)
    /// CARTESIAN_TRANSFORMATION_OPERATOR_3D(axis3)
    /// FUNCTIONALLY_DEFINED_TRANSFORMATION(name, description)
    /// GEOMETRIC_REPRESENTATION_ITEM() REPRESENTATION_ITEM(name))`.
    /// `None` when `id` is not exactly that instance.
    fn complex_operator(&self, id: u32) -> Result<Option<Placement>> {
        let components = self.components(id)?;
        let names = [
            "CARTESIAN_TRANSFORMATION_OPERATOR",
            "CARTESIAN_TRANSFORMATION_OPERATOR_3D",
            "FUNCTIONALLY_DEFINED_TRANSFORMATION",
            "GEOMETRIC_REPRESENTATION_ITEM",
            "REPRESENTATION_ITEM",
        ];
        if components.len() != names.len()
            || components
                .iter()
                .zip(names)
                .any(|(entity, name)| entity.name != name)
        {
            return Ok(None);
        }
        let [axis1, axis2, origin, scale] = arity(&components[0].args)?;
        let [axis3] = arity(&components[1].args)?;
        let [_, _] = arity(&components[2].args)?;
        let [] = arity(&components[3].args)?;
        let [_] = arity(&components[4].args)?;
        self.operator_attributes([axis1, axis2, origin, scale, axis3])
            .map(Some)
    }

    fn operator_attributes(
        &self,
        [axis1, axis2, origin, scale, axis3]: [&Value<'_>; 5],
    ) -> Result<Placement> {
        let optional = |value: &Value<'_>| match value {
            Value::Omitted => Ok(None),
            value => self.direction(reference(value)?).map(Some),
        };
        match scale {
            Value::Omitted => {}
            Value::Number(scale) => check(*scale == 1., ReadError::Unsupported)?,
            _ => return Err(ReadError::Invalid),
        }
        let z = optional(axis3)?.unwrap_or([0., 0., 1.]);
        let x = project(optional(axis1)?, z, None)?;
        let y = project(optional(axis2)?, z, Some(x))?;
        let handed = cross(x, y);
        let right = handed.iter().zip(z).map(|(a, b)| a * b).sum::<f64>();
        check(right > 0., ReadError::Unsupported)?;
        Ok(Placement {
            origin: self.point(reference(origin)?)?,
            axes: [x, y, z],
        })
    }
}

/// `first_proj_axis` (no `x`) and `second_proj_axis` of ISO 10303-42:
/// `arg` with its parts along `z` and `x` removed, or the default axis.
fn project(arg: Option<Point>, z: Point, x: Option<Point>) -> Result<Point> {
    let dot = |a: Point, b: Point| a.iter().zip(b).map(|(p, q)| p * q).sum::<f64>();
    let along = |v: Point, axis: Point| {
        let d = dot(v, axis);
        std::array::from_fn(|i| v[i] - d * axis[i])
    };
    let arg = match (arg, x) {
        (Some(arg), _) => arg,
        (None, Some(x)) => cross(z, x),
        (None, None) if z[0].abs() < 1. => [1., 0., 0.],
        (None, None) => [0., 1., 0.],
    };
    let v = along(arg, z);
    normalized(x.map_or(v, |x| along(v, x)))
}

/// Read the product structure of an assembly file. Solids are identified by
/// their entity ids.
pub(super) fn read(source: &Source<'_>) -> Result<Assembly<u32>> {
    let mut representations = BTreeMap::new();
    let mut identities = Vec::new();
    let mut transformed = BTreeMap::new();
    let mut shapes = Vec::new();
    let mut usages = BTreeMap::new();
    let mut contexts = BTreeMap::new();
    let mut mapped = BTreeSet::new();
    for (&id, components) in &source.entities {
        let find = |name: &str| components.iter().find(|c| c.name == name);
        for entity in components {
            if let (true, [_, Value::List(items), context]) = (
                entity.name.ends_with("SHAPE_REPRESENTATION"),
                entity.args.as_slice(),
            ) {
                representations.insert(
                    id,
                    Representation {
                        items,
                        context: reference(context)?,
                    },
                );
            }
        }
        match components.as_slice() {
            [Entity { name, args, .. }] => match *name {
                "SHAPE_REPRESENTATION_RELATIONSHIP" => identities.push(pair(args)?),
                "SHAPE_DEFINITION_REPRESENTATION" => {
                    let [definition, representation] = arity(args)?;
                    shapes.push([reference(definition)?, reference(representation)?]);
                }
                "NEXT_ASSEMBLY_USAGE_OCCURRENCE" => {
                    let [label, name, _, parent, child, _] = arity(args)?;
                    let name = Some(text(name)?).filter(|name| !name.is_empty());
                    usages.insert(
                        id,
                        (
                            [reference(parent)?, reference(child)?],
                            name.map_or_else(|| text(label), Ok)?,
                        ),
                    );
                }
                "CONTEXT_DEPENDENT_SHAPE_REPRESENTATION" => {
                    let [relation, shape] = arity(args)?;
                    let shape = source.args(reference(shape)?, "PRODUCT_DEFINITION_SHAPE")?;
                    let [_, _, usage] = arity(shape)?;
                    let previous = contexts.insert(reference(usage)?, reference(relation)?);
                    check(previous.is_none(), ReadError::Unsupported)?;
                }
                "MAPPED_ITEM" => {
                    mapped.insert(id);
                }
                _ => {}
            },
            _ => {
                if let Some(relation) = find("REPRESENTATION_RELATIONSHIP") {
                    let reps = pair(&relation.args)?;
                    match find("REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION") {
                        Some(Entity { args, .. }) => {
                            let [transformation] = arity(args)?;
                            transformed.insert(id, (reps, reference(transformation)?));
                        }
                        None => identities.push(reps),
                    }
                }
            }
        }
    }

    let mut groups = Groups(representations.keys().map(|&id| (id, id)).collect());
    for [a, b] in identities {
        if representations.contains_key(&a) && representations.contains_key(&b) {
            groups.join(a, b);
        }
    }
    let mut group = |id: u32| -> Result<u32> {
        check(representations.contains_key(&id), ReadError::Invalid)?;
        Ok(groups.of(id))
    };

    // The shape of every product definition and the name of its product.
    let mut defined = BTreeMap::new();
    let mut names = BTreeMap::new();
    for [shape, representation] in shapes {
        // A shape of anything but a product definition, such as of an
        // occurrence or a property, places nothing; a solid that only such a
        // shape holds is refused below.
        let Some(args) = kind(source, shape, "PRODUCT_DEFINITION_SHAPE")? else {
            continue;
        };
        let [_, _, definition] = arity(args)?;
        let definition = reference(definition)?;
        let Some(args) = kind(source, definition, "PRODUCT_DEFINITION")? else {
            continue;
        };
        let [_, _, formation, _] = arity(args)?;
        let shape = group(representation)?;
        check(
            *defined.entry(definition).or_insert(shape) == shape,
            ReadError::Unsupported,
        )?;
        let formation = source.simple(reference(formation)?)?;
        let product = match (formation.name, formation.args.as_slice()) {
            (
                "PRODUCT_DEFINITION_FORMATION"
                | "PRODUCT_DEFINITION_FORMATION_WITH_SPECIFIED_SOURCE",
                [_, _, product, ..],
            ) => reference(product)?,
            _ => return Err(ReadError::Invalid),
        };
        let [_, name, _, _] = arity(source.args(product, "PRODUCT")?)?;
        names.entry(shape).or_insert(text(name)?);
    }

    // Occurrences of products and their placements.
    let mut children: BTreeMap<u32, Vec<(u32, Placement, String)>> = BTreeMap::new();
    let mut components = BTreeSet::new();
    for (usage, ([parent, child], name)) in usages {
        let relation = contexts.remove(&usage).ok_or(ReadError::Unsupported)?;
        let (reps, transformation) = transformed.remove(&relation).ok_or(ReadError::Invalid)?;
        // A dangling definition is invalid; a product without a shape
        // cannot be placed.
        let shape = |definition| {
            source.args(definition, "PRODUCT_DEFINITION")?;
            defined
                .get(&definition)
                .copied()
                .ok_or(ReadError::Unsupported)
        };
        let [parent_shape, child_shape] = [shape(parent)?, shape(child)?];
        let ends = [group(reps[0])?, group(reps[1])?];
        let reversed = match ends {
            _ if ends == [child_shape, parent_shape] => false,
            _ if ends == [parent_shape, child_shape] => true,
            _ => return Err(ReadError::Invalid),
        };
        let placement = relate(source, &representations, reps, transformation)?;
        let placement = if reversed {
            placement.inverse()
        } else {
            placement
        };
        components.insert(child);
        children
            .entry(parent_shape)
            .or_default()
            .push((child_shape, placement, name));
    }
    // A placement must place an occurrence, and an occurrence context must
    // belong to one.
    check(
        contexts.is_empty() && transformed.is_empty(),
        ReadError::Unsupported,
    )?;

    // Mapped items place another representation inside their own, once.
    let mut held = BTreeSet::new();
    for (&id, representation) in &representations {
        for item in representation.items {
            let item = reference(item)?;
            if !mapped.contains(&item) {
                continue;
            }
            check(held.insert(item), ReadError::Unsupported)?;
            let [name, map, target] = arity(source.args(item, "MAPPED_ITEM")?)?;
            let [origin, shape] = arity(source.args(reference(map)?, "REPRESENTATION_MAP")?)?;
            let origin = source.placement(reference(origin)?)?;
            let target = source.placement(reference(target)?)?;
            let placement = agreed(origin, target, source.tolerance)?;
            let shape = group(reference(shape)?)?;
            children
                .entry(group(id)?)
                .or_default()
                .push((shape, placement, text(name)?));
        }
    }
    // A mapped item that no representation holds is never placed.
    check(held == mapped, ReadError::Unsupported)?;

    let mut roots = Vec::new();
    for (definition, shape) in &defined {
        if !components.contains(definition) && !roots.contains(shape) {
            roots.push(*shape);
        }
    }
    check(!roots.is_empty(), ReadError::Unsupported)?;

    // Number the nodes reached from the roots, refusing cycles.
    let mut index = BTreeMap::new();
    let mut order = Vec::new();
    let mut open = BTreeSet::new();
    for &root in &roots {
        number(root, 0, &children, &mut index, &mut order, &mut open)?;
    }
    let solids: BTreeSet<u32> = source.solids.iter().copied().collect();
    let mut nodes = Vec::with_capacity(order.len());
    let mut placed = BTreeSet::new();
    for &shape in &order {
        let mut own = Vec::new();
        for (&id, representation) in &representations {
            if group(id)? != shape {
                continue;
            }
            check(
                source
                    .components(representation.context)?
                    .iter()
                    .any(|part| part.name == "GLOBAL_UNIT_ASSIGNED_CONTEXT"),
                ReadError::Unsupported,
            )?;
            for item in representation.items {
                let item = reference(item)?;
                if solids.contains(&item) {
                    own.push(item);
                    placed.insert(item);
                }
            }
        }
        nodes.push(Node {
            name: names.get(&shape).cloned().unwrap_or_default(),
            solids: own,
            children: children
                .get(&shape)
                .into_iter()
                .flatten()
                .map(|(child, placement, name)| Instance {
                    node: index[child],
                    placement: *placement,
                    name: name.clone(),
                })
                .collect(),
        });
    }
    // Every solid has a place.
    check(placed == solids, ReadError::Unsupported)?;
    let assembly = Assembly {
        roots: roots.iter().map(|root| index[root]).collect(),
        nodes,
    };
    check(uses(&assembly) <= super::MAX_SOLIDS, ReadError::Limit)?;
    Ok(assembly)
}

/// The number of placed solids, saturating above the solid budget.
fn uses<S>(assembly: &Assembly<S>) -> usize {
    fn count<S>(assembly: &Assembly<S>, node: usize, known: &mut [Option<usize>]) -> usize {
        if let Some(count) = known[node] {
            return count;
        }
        let entry = &assembly.nodes[node];
        let mut total = entry.solids.len();
        for child in &entry.children {
            total = total
                .saturating_add(count(assembly, child.node, known))
                .min(super::MAX_SOLIDS + 1);
        }
        known[node] = Some(total);
        total
    }
    let mut known = vec![None; assembly.nodes.len()];
    let mut total = 0usize;
    for &root in &assembly.roots {
        total = total.saturating_add(count(assembly, root, &mut known));
    }
    total
}

/// Number `shape` and the shapes below it, depth first. A shape met again
/// while it is still open contains itself.
fn number(
    shape: u32,
    depth: usize,
    children: &BTreeMap<u32, Vec<(u32, Placement, String)>>,
    index: &mut BTreeMap<u32, usize>,
    order: &mut Vec<u32>,
    open: &mut BTreeSet<u32>,
) -> Result<()> {
    check(!open.contains(&shape), ReadError::Invalid)?;
    if index.contains_key(&shape) {
        return Ok(());
    }
    check(depth <= MAX_DEPTH, ReadError::Limit)?;
    index.insert(shape, order.len());
    order.push(shape);
    open.insert(shape);
    for (child, ..) in children.get(&shape).into_iter().flatten() {
        number(*child, depth + 1, children, index, order, open)?;
    }
    open.remove(&shape);
    Ok(())
}

/// The placement of `reps[0]`'s coordinates in `reps[1]`'s.
fn relate(
    source: &Source<'_>,
    representations: &BTreeMap<u32, Representation<'_, '_>>,
    reps: [u32; 2],
    transformation: u32,
) -> Result<Placement> {
    if let Some(placement) = source.complex_operator(transformation)? {
        return Ok(placement);
    }
    let entity = source.simple(transformation)?;
    match (entity.name, entity.args.as_slice()) {
        ("ITEM_DEFINED_TRANSFORMATION", [_, _, first, second]) => {
            let items = [reference(first)?, reference(second)?];
            for (item, rep) in items.iter().zip(reps) {
                let holds = representations[&rep]
                    .items
                    .iter()
                    .any(|value| reference(value) == Ok(*item));
                check(holds, ReadError::Invalid)?;
            }
            let [origin, target] = [source.placement(items[0])?, source.placement(items[1])?];
            agreed(origin, target, source.tolerance)
        }
        ("CARTESIAN_TRANSFORMATION_OPERATOR_3D", args) => source.operator(args),
        _ => Err(ReadError::Unsupported),
    }
}

/// The placement that carries `origin` onto `target`.
///
/// ISO 10303-43 maps the first item onto the second, `target · origin⁻¹`;
/// OpenCascade 7.9 applies `origin⁻¹ · target`. The two agree whenever the
/// origin item is the identity, as producers write it; a placement on which
/// they disagree has no single reading and is refused.
fn agreed(origin: Placement, target: Placement, tolerance: f64) -> Result<Placement> {
    let standard = target.then(&origin.inverse());
    let [axes, offset] = standard.distance(&origin.inverse().then(&target));
    check(axes <= 1e-12 && offset <= tolerance, ReadError::Unsupported)?;
    Ok(standard)
}

#[cfg(test)]
mod tests;
