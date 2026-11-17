//! Schema-helper contracts over reparsed physical bytes. Some malformed
//! metadata is rejected by the neutral reader before this private validator.
use super::*;

#[test]
fn genuine_physical_metadata_roles_and_arities_cannot_authorize_other_records() {
    let actual = super::super::super::tests::candidate(false);
    let arena = crate::step_read::part21::parse(&actual).unwrap();
    let id = |name: &str| {
        *arena
            .iter()
            .find(|(_, parts)| parts.iter().any(|e| e.name == name))
            .unwrap()
            .0
    };
    let app = id("APPLICATION_CONTEXT");
    let product = id("PRODUCT");
    let context = id("PRODUCT_CONTEXT");
    let formation = id("PRODUCT_DEFINITION_FORMATION");
    let definition = id("PRODUCT_DEFINITION");
    let definition_context = id("PRODUCT_DEFINITION_CONTEXT");
    let shape = id("PRODUCT_DEFINITION_SHAPE");
    let representation = id("ADVANCED_BREP_SHAPE_REPRESENTATION");
    let solid = id("MANIFOLD_SOLID_BREP");
    let schema = reference(&arena[&representation][0].args[2]).unwrap();
    let length_unit = id("LENGTH_UNIT");
    let angle_unit = id("PLANE_ANGLE_UNIT");
    let parametric_context = id("PARAMETRIC_REPRESENTATION_CONTEXT");
    let cases = [
        (
            "SHAPE_DEFINITION_REPRESENTATION",
            "SHAPE_DEFINITION_REPRESENTATION()".to_owned(),
        ),
        (
            "ADVANCED_BREP_SHAPE_REPRESENTATION",
            format!("ADVANCED_BREP_SHAPE_REPRESENTATION('',(#{solid},#{length_unit}),#{schema})"),
        ),
        (
            "ADVANCED_BREP_SHAPE_REPRESENTATION",
            format!(
                "ADVANCED_BREP_SHAPE_REPRESENTATION('',(#{solid},#{solid}),#{parametric_context})"
            ),
        ),
        (
            "DEFINITIONAL_REPRESENTATION",
            format!("DEFINITIONAL_REPRESENTATION('',(),#{parametric_context})"),
        ),
        ("APPLICATION_CONTEXT", "APPLICATION_CONTEXT(1.)".to_owned()),
        ("APPLICATION_CONTEXT", "UNSUPPORTED_METADATA('')".to_owned()),
        (
            "GEOMETRIC_REPRESENTATION_CONTEXT",
            "GEOMETRIC_REPRESENTATION_CONTEXT(4.)".to_owned(),
        ),
        (
            "GEOMETRIC_REPRESENTATION_CONTEXT",
            "GEOMETRIC_REPRESENTATION_CONTEXT(4)".to_owned(),
        ),
        (
            "GEOMETRIC_REPRESENTATION_CONTEXT",
            "GEOMETRIC_REPRESENTATION_CONTEXT()".to_owned(),
        ),
        (
            "REPRESENTATION_CONTEXT",
            "REPRESENTATION_CONTEXT(1.,'')".to_owned(),
        ),
        (
            "GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT",
            "GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT(())".to_owned(),
        ),
        (
            "GLOBAL_UNIT_ASSIGNED_CONTEXT",
            "GLOBAL_UNIT_ASSIGNED_CONTEXT(())".to_owned(),
        ),
        (
            "GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT",
            "GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT()".to_owned(),
        ),
        (
            "GLOBAL_UNIT_ASSIGNED_CONTEXT",
            "GLOBAL_UNIT_ASSIGNED_CONTEXT()".to_owned(),
        ),
        ("LENGTH_UNIT", "LENGTH_UNIT(1.)".to_owned()),
        ("PLANE_ANGLE_UNIT", "PLANE_ANGLE_UNIT(1.)".to_owned()),
        ("SOLID_ANGLE_UNIT", "SOLID_ANGLE_UNIT(1.)".to_owned()),
        ("NAMED_UNIT", "NAMED_UNIT()".to_owned()),
        ("SI_UNIT", "SI_UNIT(.KILO.,.METRE.)".to_owned()),
        (
            "UNCERTAINTY_MEASURE_WITH_UNIT",
            "UNCERTAINTY_MEASURE_WITH_UNIT()".to_owned(),
        ),
        (
            "UNCERTAINTY_MEASURE_WITH_UNIT",
            format!("UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(-1.),#{length_unit},'','')"),
        ),
        (
            "UNCERTAINTY_MEASURE_WITH_UNIT",
            format!("UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(0.001),#{angle_unit},'','')"),
        ),
        ("PRODUCT_CONTEXT", format!("PRODUCT_CONTEXT(1.,#{app},'')")),
        (
            "PRODUCT_CONTEXT",
            format!("PRODUCT_CONTEXT('',#{product},'')"),
        ),
        ("PRODUCT", format!("PRODUCT('',1.,'',(#{context}))")),
        ("PRODUCT", "PRODUCT('','','',())".to_owned()),
        ("PRODUCT", format!("PRODUCT('','','',(#{app}))")),
        (
            "PRODUCT_DEFINITION_FORMATION",
            format!("PRODUCT_DEFINITION_FORMATION(1.,'',#{product})"),
        ),
        (
            "PRODUCT_DEFINITION_FORMATION",
            format!("PRODUCT_DEFINITION_FORMATION('','',#{app})"),
        ),
        (
            "PRODUCT_DEFINITION",
            format!("PRODUCT_DEFINITION(1.,'',#{formation},#{definition_context})"),
        ),
        (
            "PRODUCT_DEFINITION",
            format!("PRODUCT_DEFINITION('','',#{product},#{definition_context})"),
        ),
        (
            "PRODUCT_DEFINITION_SHAPE",
            format!("PRODUCT_DEFINITION_SHAPE(1.,'',#{definition})"),
        ),
        (
            "SHAPE_DEFINITION_REPRESENTATION",
            format!("SHAPE_DEFINITION_REPRESENTATION(#{product},#{representation})"),
        ),
        (
            "SHAPE_DEFINITION_REPRESENTATION",
            format!("SHAPE_DEFINITION_REPRESENTATION(#{shape},#{product})"),
        ),
        (
            "ADVANCED_BREP_SHAPE_REPRESENTATION",
            format!("ADVANCED_BREP_SHAPE_REPRESENTATION(1.,(#{solid}),#{schema})"),
        ),
        (
            "ADVANCED_BREP_SHAPE_REPRESENTATION",
            format!("ADVANCED_BREP_SHAPE_REPRESENTATION('',(#{solid}),#{schema})"),
        ),
        (
            "ADVANCED_BREP_SHAPE_REPRESENTATION",
            format!("ADVANCED_BREP_SHAPE_REPRESENTATION('',(#{solid},#{solid}),#{schema})"),
        ),
        (
            "ADVANCED_BREP_SHAPE_REPRESENTATION",
            format!("ADVANCED_BREP_SHAPE_REPRESENTATION('',(#{product},#{solid}),#{schema})"),
        ),
        (
            "DEFINITIONAL_REPRESENTATION",
            format!("DEFINITIONAL_REPRESENTATION(1.,(),#{schema})"),
        ),
    ];
    for (kind, replacement) in cases {
        let owners = arena
            .iter()
            .filter_map(|(&owner, parts)| {
                parts
                    .iter()
                    .find(|e| e.name == kind)
                    .map(|component| (owner, component))
            })
            .collect::<Vec<_>>();
        assert!(!owners.is_empty());
        for (owner, component) in owners {
            let source = component.source.unwrap();
            let start = source.as_ptr() as usize - actual.as_ptr() as usize;
            let mut changed = actual.clone();
            changed.replace_range(start..start + source.len(), &replacement);
            let records = crate::step_read::part21::parse(&changed).unwrap();
            let mut budget =
                Budget::new(crate::step_read::npt_stock::Options::default(), None).unwrap();
            assert!(
                matches!(
                    check(&records, &records[&owner], &mut budget),
                    Err(Error::Invalid | Error::Unsupported)
                ),
                "{kind}: {replacement}"
            );
        }
    }
    let items = list(&arena[&representation][0].args[1]).unwrap();
    let axis = items
        .iter()
        .map(|v| reference(v).unwrap())
        .find(|id| arena[id].iter().any(|e| e.name == "AXIS2_PLACEMENT_3D"))
        .unwrap();
    let origin = reference(&arena[&axis][0].args[1]).unwrap();
    let normal = reference(&arena[&axis][0].args[2]).unwrap();
    let reference_direction = reference(&arena[&axis][0].args[3]).unwrap();
    for (owner, replacement) in [
        (axis, "AXIS2_PLACEMENT_3D('')".to_owned()),
        (
            axis,
            format!("AXIS2_PLACEMENT_3D('',#{product},#{product},#{product})"),
        ),
        (origin, "CARTESIAN_POINT('',(1.,0.,0.))".to_owned()),
        (normal, "DIRECTION('',(0.,0.,-1.))".to_owned()),
        (reference_direction, "DIRECTION('',(-1.,0.,0.))".to_owned()),
    ] {
        let source = arena[&owner][0].source.unwrap();
        let start = source.as_ptr() as usize - actual.as_ptr() as usize;
        let mut changed = actual.clone();
        changed.replace_range(start..start + source.len(), &replacement);
        let records = crate::step_read::part21::parse(&changed).unwrap();
        let mut budget =
            Budget::new(crate::step_read::npt_stock::Options::default(), None).unwrap();
        assert!(matches!(
            check(&records, &records[&representation], &mut budget),
            Err(Error::Invalid | Error::Unsupported)
        ));
    }
    // Empty component input is a separate helper contract; the physical
    // parser rejects an empty record before it could reach this validator.
    let mut budget = Budget::new(crate::step_read::npt_stock::Options::default(), None).unwrap();
    assert!(matches!(
        check(&arena, &[], &mut budget),
        Err(Error::Invalid)
    ));
    // Reuse the untouched original raw records after every rejected mutation.
    for parts in arena.values() {
        let mut budget =
            Budget::new(crate::step_read::npt_stock::Options::default(), None).unwrap();
        check(&arena, parts, &mut budget).unwrap();
    }
}

#[test]
fn actual_optional_protocol_metadata_keeps_edition_and_owner_roles() {
    let mut actual = super::super::super::tests::candidate(false);
    let arena = crate::step_read::part21::parse(&actual).unwrap();
    let owner = |name| {
        *arena
            .iter()
            .find(|(_, p)| p.iter().any(|e| e.name == name))
            .unwrap()
            .0
    };
    let app = owner("APPLICATION_CONTEXT");
    let product = owner("PRODUCT");
    let protocol = 999_999_997;
    let category = 999_999_998;
    assert!(!arena.contains_key(&protocol) && !arena.contains_key(&category));
    let protocol_body = format!(
        "APPLICATION_PROTOCOL_DEFINITION('international standard','automotive_design',2000,#{app})"
    );
    let category_body = format!("PRODUCT_RELATED_PRODUCT_CATEGORY('part',$,(#{product}))");
    let extra = format!("#{protocol}={protocol_body};\n#{category}={category_body};\n");
    let end = actual.rfind("ENDSEC;").unwrap();
    actual.insert_str(end, &extra);
    let records = crate::step_read::part21::parse(&actual).unwrap();
    for id in [protocol, category] {
        let mut b = Budget::new(crate::step_read::npt_stock::Options::default(), None).unwrap();
        check(&records, &records[&id], &mut b).unwrap();
    }
    for (id, replacement) in [
        (
            protocol,
            format!("APPLICATION_PROTOCOL_DEFINITION('unknown','automotive_design',2000,#{app})"),
        ),
        (
            protocol,
            format!(
                "APPLICATION_PROTOCOL_DEFINITION('international standard','automotive_design',1999,#{app})"
            ),
        ),
        (
            protocol,
            format!(
                "APPLICATION_PROTOCOL_DEFINITION('international standard','automotive_design',2000,#{product})"
            ),
        ),
        (
            category,
            format!("PRODUCT_RELATED_PRODUCT_CATEGORY('unknown',$,(#{product}))"),
        ),
        (
            category,
            "PRODUCT_RELATED_PRODUCT_CATEGORY('part',$,())".to_owned(),
        ),
        (
            category,
            format!("PRODUCT_RELATED_PRODUCT_CATEGORY('part',$,(#{app}))"),
        ),
    ] {
        let source = records[&id][0].source.unwrap();
        let start = source.as_ptr() as usize - actual.as_ptr() as usize;
        let mut changed = actual.clone();
        changed.replace_range(start..start + source.len(), &replacement);
        let parsed = crate::step_read::part21::parse(&changed).unwrap();
        let mut b = Budget::new(crate::step_read::npt_stock::Options::default(), None).unwrap();
        assert!(matches!(
            check(&parsed, &parsed[&id], &mut b),
            Err(Error::Invalid | Error::Unsupported)
        ));
    }
}
