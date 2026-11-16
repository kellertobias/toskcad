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
    let schema = id("GEOMETRIC_REPRESENTATION_CONTEXT");
    let cases = [
        ("APPLICATION_CONTEXT", "APPLICATION_CONTEXT(1.)".to_owned()),
        ("APPLICATION_CONTEXT", "UNSUPPORTED_METADATA('')".to_owned()),
        (
            "GEOMETRIC_REPRESENTATION_CONTEXT",
            "GEOMETRIC_REPRESENTATION_CONTEXT(4.)".to_owned(),
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
        ("LENGTH_UNIT", "LENGTH_UNIT(1.)".to_owned()),
        ("PLANE_ANGLE_UNIT", "PLANE_ANGLE_UNIT(1.)".to_owned()),
        ("SOLID_ANGLE_UNIT", "SOLID_ANGLE_UNIT(1.)".to_owned()),
        ("NAMED_UNIT", "NAMED_UNIT()".to_owned()),
        ("SI_UNIT", "SI_UNIT(.KILO.,.METRE.)".to_owned()),
        (
            "UNCERTAINTY_MEASURE_WITH_UNIT",
            "UNCERTAINTY_MEASURE_WITH_UNIT()".to_owned(),
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
        let owner = id(kind);
        let component = arena[&owner].iter().find(|e| e.name == kind).unwrap();
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
    // Reuse the untouched original raw records after every rejected mutation.
    for parts in arena.values() {
        let mut budget =
            Budget::new(crate::step_read::npt_stock::Options::default(), None).unwrap();
        check(&arena, parts, &mut budget).unwrap();
    }
}
