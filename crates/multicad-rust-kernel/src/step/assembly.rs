//! Product structure written back: every definition once, in its own
//! coordinates, with the occurrences that place it.
use super::*;
use crate::brep::Solid;
use crate::step_read::brep::{Assembly, Placement};

/// A Part 21 string with its apostrophes doubled.
struct Text<'a>(&'a str);

impl fmt::Display for Text<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "'{}'", self.0.replace('\'', "''"))
    }
}

/// A written definition: its name, its own solids and its placed children.
type Written<'a> = (&'a str, &'a [Solid], Vec<(usize, Placement, &'a str)>);

fn axis(file: &mut StepFile, placement: &Placement) -> Entity {
    let origin = file.point(placement.origin);
    file.axis(origin, placement.axes[2], placement.axes[0])
}

/// Write `assembly` with its roots placed by `placement`.
///
/// Every node becomes one product whose representation holds its own
/// solids and the placements of its children. Every instance becomes a
/// next-assembly-usage occurrence whose item-defined transformation carries
/// the child's origin onto its placement in the parent. Roots that are not
/// at the identity are held by one unnamed product that places them.
pub fn write_assembly(
    assembly: &Assembly<Solid>,
    placement: Placement,
) -> Result<String, StepError> {
    let mut file = StepFile::new();
    let application = file.entity(format_args!(
        "APPLICATION_CONTEXT('core data for automotive mechanical design processes')"
    ));
    file.entity(format_args!("APPLICATION_PROTOCOL_DEFINITION('international standard','automotive_design',2000,{application})"));
    let product_context = file.entity(format_args!(
        "PRODUCT_CONTEXT('',{application},'mechanical')"
    ));
    let definition_context = file.entity(format_args!(
        "PRODUCT_DEFINITION_CONTEXT('part definition',{application},'design')"
    ));
    let length = file.entity(format_args!(
        "(LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.))"
    ));
    let angle = file.entity(format_args!(
        "(NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.))"
    ));
    let solid_angle = file.entity(format_args!(
        "(NAMED_UNIT(*) SI_UNIT($,.STERADIAN.) SOLID_ANGLE_UNIT())"
    ));
    let uncertainty = file.entity(format_args!(
        "UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-7),{length},'distance_accuracy_value','')"
    ));
    let context = file.entity(format_args!("(GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT(({uncertainty})) GLOBAL_UNIT_ASSIGNED_CONTEXT(({length},{angle},{solid_angle})) REPRESENTATION_CONTEXT('',''))"));
    let origin = axis(&mut file, &Placement::IDENTITY);

    // A wrapper places roots that are not at the identity.
    let wrapped = placement != Placement::IDENTITY;
    let mut nodes: Vec<Written<'_>> = assembly
        .nodes
        .iter()
        .map(|node| {
            let children = node
                .children
                .iter()
                .map(|child| (child.node, child.placement, child.name.as_str()))
                .collect();
            (node.name.as_str(), node.solids.as_slice(), children)
        })
        .collect();
    if wrapped {
        let roots = assembly
            .roots
            .iter()
            .map(|&root| (root, placement, ""))
            .collect();
        nodes.push(("", &[], roots));
    }

    let mut representations = Vec::with_capacity(nodes.len());
    let mut placements = Vec::with_capacity(nodes.len());
    for (_, solids, children) in &nodes {
        let mut items = vec![origin];
        for solid in *solids {
            items.push(super::brep::encode(&mut file, solid)?);
        }
        let mut placed = Vec::with_capacity(children.len());
        for (_, placement, _) in children {
            let item = axis(&mut file, placement);
            placed.push(item);
            items.push(item);
        }
        let kind = if solids.is_empty() {
            "SHAPE_REPRESENTATION"
        } else {
            "ADVANCED_BREP_SHAPE_REPRESENTATION"
        };
        representations
            .push(file.entity(format_args!("{kind}('',{},{context})", Entities(&items))));
        placements.push(placed);
    }
    let mut definitions = Vec::with_capacity(nodes.len());
    for ((name, ..), representation) in nodes.iter().zip(&representations) {
        let name = Text(name);
        let product = file.entity(format_args!(
            "PRODUCT({name},{name},'',({product_context}))"
        ));
        let formation = file.entity(format_args!(
            "PRODUCT_DEFINITION_FORMATION('','',{product})"
        ));
        let definition = file.entity(format_args!(
            "PRODUCT_DEFINITION('design','',{formation},{definition_context})"
        ));
        let shape = file.entity(format_args!("PRODUCT_DEFINITION_SHAPE('','',{definition})"));
        file.entity(format_args!(
            "SHAPE_DEFINITION_REPRESENTATION({shape},{representation})"
        ));
        file.entity(format_args!(
            "PRODUCT_RELATED_PRODUCT_CATEGORY('part',$,({product}))"
        ));
        definitions.push(definition);
    }
    let mut count = 0;
    for (parent, (_, _, children)) in nodes.iter().enumerate() {
        for ((child, _, name), item) in children.iter().zip(&placements[parent]) {
            count += 1;
            let usage = file.entity(format_args!(
                "NEXT_ASSEMBLY_USAGE_OCCURRENCE('{count}',{},'',{},{},$)",
                Text(name),
                definitions[parent],
                definitions[*child]
            ));
            let shape = file.entity(format_args!(
                "PRODUCT_DEFINITION_SHAPE('Placement','',{usage})"
            ));
            let transformation = file.entity(format_args!(
                "ITEM_DEFINED_TRANSFORMATION('','',{origin},{item})"
            ));
            let relation = file.entity(format_args!(
                "(REPRESENTATION_RELATIONSHIP('','',{},{}) REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION({transformation}) SHAPE_REPRESENTATION_RELATIONSHIP())",
                representations[*child],
                representations[parent]
            ));
            file.entity(format_args!(
                "CONTEXT_DEPENDENT_SHAPE_REPRESENTATION({relation},{shape})"
            ));
        }
    }
    file.data.push_str("ENDSEC;\nEND-ISO-10303-21;\n");
    if writer::exceeds_readback_budget(&file) {
        return Err(StepError::Limit);
    }
    Ok(file.data)
}
