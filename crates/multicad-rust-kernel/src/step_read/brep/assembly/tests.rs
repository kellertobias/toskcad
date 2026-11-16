//! Product structure against OpenCascade's two-body assembly export and
//! tampered or extended copies of it.
use super::*;
use crate::brep::{Error, Solid};
use crate::step::write_assembly;
use crate::step_read::brep::{has_product_structure, read_assembly, read_boundaries, read_solids};
use std::f64::consts::PI;

const TOLERANCE: f64 = 1e-7;

/// OpenCascade 7.9.3's assembly of a 10 x 10 x 5 mm box (#31, shape #36)
/// and a radius 4 mm, 8 mm cylinder at (30, 0) (#380, shape #385) under a
/// root product (#5, shape #10).
const TWO: &str =
    include_str!("../../../../../multicad-app/tests/fixtures/exchange/corpus/occt-two-bodies.step");
const FLAT: &str =
    include_str!("../../../../../multicad-app/tests/fixtures/exchange/corpus/occt-cylinder.step");

const BOX: f64 = 500.;
const CYLINDER: f64 = 128. * PI;

fn swap(text: &str, from: &str, to: &str) -> String {
    assert!(text.contains(from), "{from}");
    text.replacen(from, to, 1)
}

/// Append entities to the data section.
fn with(text: &str, extra: &str) -> String {
    let end = text.rfind("ENDSEC;").unwrap();
    format!("{}{extra}\n{}", &text[..end], &text[end..])
}

/// List `items` after the root representation's own.
fn rooted(text: &str, items: &str) -> String {
    swap(
        text,
        "#10 = SHAPE_REPRESENTATION('',(#11,#15,#19),#23);",
        &format!("#10 = SHAPE_REPRESENTATION('',(#11,#15,#19,{items}),#23);"),
    )
}

fn axis(id: u32, origin: Point, z: Point, x: Point) -> String {
    let p = |v: Point| format!("({:?},{:?},{:?})", v[0], v[1], v[2]);
    format!(
        "#{id}=AXIS2_PLACEMENT_3D('',#{},#{},#{});#{}=CARTESIAN_POINT('',{});#{}=DIRECTION('',{});#{}=DIRECTION('',{});",
        id + 1,
        id + 2,
        id + 3,
        id + 1,
        p(origin),
        id + 2,
        p(z),
        id + 3,
        p(x)
    )
}

/// An occurrence of `child` in `parent` (product definitions) whose
/// relationship relates `reps` by the items `placed`, using ids from `id`.
fn occurrence(
    id: u32,
    name: &str,
    [parent, child]: [u32; 2],
    reps: [u32; 2],
    placed: [u32; 2],
) -> String {
    format!(
        "#{id}=NEXT_ASSEMBLY_USAGE_OCCURRENCE('{id}','{name}','',#{parent},#{child},$);\
         #{}=PRODUCT_DEFINITION_SHAPE('','',#{id});\
         #{}=ITEM_DEFINED_TRANSFORMATION('','',#{},#{});\
         #{}=(REPRESENTATION_RELATIONSHIP('','',#{},#{})REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION(#{})SHAPE_REPRESENTATION_RELATIONSHIP());\
         #{}=CONTEXT_DEPENDENT_SHAPE_REPRESENTATION(#{},#{});",
        id + 1,
        id + 2,
        placed[0],
        placed[1],
        id + 3,
        reps[0],
        reps[1],
        id + 2,
        id + 4,
        id + 3,
        id + 1
    )
}

/// A product with an empty shape at ids from `id`: (definition, shape).
fn product(id: u32, name: &str) -> (String, [u32; 2]) {
    (
        format!(
            "#{id}=PRODUCT('{name}','{name}','',(#8));#{}=PRODUCT_DEFINITION_FORMATION('','',#{id});\
             #{}=PRODUCT_DEFINITION('design','',#{},#9);#{}=PRODUCT_DEFINITION_SHAPE('','',#{});\
             #{}=SHAPE_REPRESENTATION('',(#11),#23);#{}=SHAPE_DEFINITION_REPRESENTATION(#{},#{});",
            id + 1,
            id + 2,
            id + 1,
            id + 3,
            id + 2,
            id + 4,
            id + 5,
            id + 3,
            id + 4
        ),
        [id + 2, id + 4],
    )
}

/// `levels` empty products below the root, each placed `fan` times in the
/// one above; the last places the box `fan` times.
fn tree(levels: u32, fan: u32) -> String {
    let mut extra = String::new();
    let mut above = [5, 10];
    let mut next = 50_000;
    for level in 0..=levels {
        let below = if level == levels {
            [31, 36]
        } else {
            let (text, below) = product(30_000 + 10 * level, &format!("level {level}"));
            extra += &text;
            below
        };
        for _ in 0..if level == 0 { 1 } else { fan } {
            extra += &occurrence(
                next,
                "",
                [above[0], below[0]],
                [below[1], above[1]],
                [11, 11],
            );
            next += 5;
        }
        above = below;
    }
    with(TWO, &extra)
}

fn read(text: &str) -> Assembly<Solid> {
    read_assembly(text, TOLERANCE, None).unwrap()
}

fn refused(text: &str) -> Error {
    read_assembly(text, TOLERANCE, None).unwrap_err()
}

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-9 * b.abs().max(1.)
}

fn same(a: &Placement, b: &Placement) -> bool {
    let [axes, origin] = a.distance(b);
    axes <= 1e-12 && origin <= 1e-9
}

fn volumes(assembly: &Assembly<Solid>) -> Vec<f64> {
    assembly
        .placed()
        .unwrap()
        .iter()
        .map(Solid::volume_mm3)
        .collect()
}

/// The repeated box placed by `#20000`: 90 degrees about z at (100, 0, 0).
fn repeated() -> String {
    let text = with(
        TWO,
        &(axis(20_000, [100., 0., 0.], [0., 0., 1.], [0., 1., 0.])
            + &occurrence(20_010, "again", [5, 31], [36, 10], [11, 20_000])),
    );
    rooted(&text, "#20000")
}

const TURNED: Placement = Placement {
    origin: [100., 0., 0.],
    axes: [[0., 1., 0.], [-1., 0., 0.], [0., 0., 1.]],
};

#[test]
fn every_product_is_placed_once_at_its_stated_occurrence() {
    let assembly = read(TWO);
    assert_eq!(assembly.roots, [0]);
    assert_eq!(assembly.nodes.len(), 3);
    let root = &assembly.nodes[0];
    assert_eq!(root.name, "Open CASCADE STEP translator 7.9 9");
    assert!(root.solids.is_empty());
    // Unnamed occurrences are known by their identifiers.
    let names: Vec<_> = root.children.iter().map(|c| c.name.as_str()).collect();
    assert_eq!(names, ["1", "2"]);
    assert_eq!(
        assembly.nodes[root.children[1].node].name,
        "Open CASCADE STEP translator 7.9 9.2"
    );
    let occurrences = assembly.occurrences();
    assert_eq!(occurrences.len(), 2);
    assert_eq!(occurrences[1].path, [0, 1]);
    assert!(
        occurrences
            .iter()
            .all(|o| o.placement == Placement::IDENTITY && o.solid == 0)
    );
    let found = volumes(&assembly);
    assert!(
        close(found[0], BOX) && close(found[1], CYLINDER),
        "{found:?}"
    );
    // The flat readers see the placed solids, or refuse the structure.
    assert_eq!(read_solids(TWO, TOLERANCE, None).unwrap().len(), 2);
    assert_eq!(
        read_boundaries(TWO, TOLERANCE).unwrap_err(),
        ReadError::Unsupported
    );
    assert!(has_product_structure(TWO));
    assert!(!has_product_structure(FLAT));
    assert!(!has_product_structure("not a step file"));
    // A flat file is one unnamed root holding its solids.
    let flat = read(FLAT);
    assert_eq!(flat.nodes.len(), 1);
    assert!(flat.nodes[0].name.is_empty() && flat.nodes[0].children.is_empty());
}

#[test]
fn a_repeated_component_shares_one_definition() {
    let assembly = read(&repeated());
    assert_eq!(assembly.nodes.len(), 3);
    let occurrences = assembly.occurrences();
    assert_eq!(occurrences.len(), 3);
    assert_eq!(occurrences[0].node, occurrences[2].node);
    assert!(same(&occurrences[2].placement, &TURNED));
    assert_eq!(assembly.nodes[0].children[2].name, "again");
    let placed = assembly.placed().unwrap();
    assert!(close(placed[2].volume_mm3(), BOX));
    let [low, high] = placed[2].bounds_mm();
    let expected = [[90., 0., 0.], [100., 10., 5.]];
    for (found, expected) in [low, high].iter().flatten().zip(expected.iter().flatten()) {
        assert!((found - expected).abs() < 1e-9, "{low:?} {high:?}");
    }
}

#[test]
fn nested_occurrences_compose_their_placements() {
    // The cylinder moves into the box, which the root offsets by (5, 6, 7);
    // the box turns the cylinder by 90 degrees about z.
    let text = swap(
        TWO,
        "NEXT_ASSEMBLY_USAGE_OCCURRENCE('2','','',#5,#380,$)",
        "NEXT_ASSEMBLY_USAGE_OCCURRENCE('2','','',#31,#380,$)",
    );
    let text = swap(&text, "('','',#385,#10)", "('','',#385,#36)");
    let text = swap(
        &text,
        "#36 = ADVANCED_BREP_SHAPE_REPRESENTATION('',(#11,#37),#367);",
        "#36 = ADVANCED_BREP_SHAPE_REPRESENTATION('',(#11,#37,#19),#367);",
    );
    let text = swap(
        &text,
        "#16 = CARTESIAN_POINT('',(0.,0.,0.));",
        "#16 = CARTESIAN_POINT('',(5.,6.,7.));",
    );
    let text = swap(
        &text,
        "#22 = DIRECTION('',(1.,0.,-0.));",
        "#22 = DIRECTION('',(0.,1.,0.));",
    );
    let assembly = read(&text);
    let occurrences = assembly.occurrences();
    assert_eq!(occurrences.len(), 2);
    assert_eq!(occurrences[1].path, [0, 0, 0]);
    let expected = Placement {
        origin: [5., 6., 7.],
        axes: TURNED.axes,
    };
    assert!(same(&occurrences[1].placement, &expected));
    let cylinder = &assembly.placed().unwrap()[1];
    let [low, high] = cylinder.bounds_mm();
    // The axis at (30, 0) turns to (0, 30) and moves by (5, 6, 7).
    assert!(
        (low[0] - 1.).abs() < 1e-9 && (high[1] - 40.).abs() < 1e-9,
        "{low:?} {high:?}"
    );
    assert!((low[2] - 7.).abs() < 1e-9);
}

#[test]
fn a_relationship_stated_from_the_parent_reads_the_same() {
    let text = repeated();
    let text = swap(
        &text,
        "ITEM_DEFINED_TRANSFORMATION('','',#11,#20000)",
        "ITEM_DEFINED_TRANSFORMATION('','',#20000,#11)",
    );
    let text = swap(
        &text,
        "('','',#36,#10)REPRESENTATION_RELATIONSHIP_WITH",
        "('','',#10,#36)REPRESENTATION_RELATIONSHIP_WITH",
    );
    let occurrences = read(&text).occurrences();
    assert!(same(&occurrences[2].placement, &TURNED));
}

#[test]
fn mapped_items_place_their_representation() {
    let mapping = axis(20_000, [100., 0., 0.], [0., 0., 1.], [0., 1., 0.])
        + "#20010=REPRESENTATION_MAP(#11,#385);#20011=MAPPED_ITEM('copy',#20010,#20000);";
    let assembly = read(&rooted(&with(TWO, &mapping), "#20011"));
    let occurrences = assembly.occurrences();
    assert_eq!(occurrences.len(), 3);
    assert_eq!(occurrences[1].node, occurrences[2].node);
    assert_eq!(assembly.nodes[0].children[2].name, "copy");
    assert!(same(&occurrences[2].placement, &TURNED));
    // An operator may be the target; an empty representation that no
    // product defines is an unnamed definition.
    let operator = "#20020=CARTESIAN_POINT('',(0.,0.,9.));\
        #20021=CARTESIAN_TRANSFORMATION_OPERATOR_3D('','','',$,$,#20020,$,$);\
        #20022=SHAPE_REPRESENTATION('',(#11),#23);#20023=REPRESENTATION_MAP(#11,#20022);\
        #20024=MAPPED_ITEM('',#20023,#20021);";
    let assembly = read(&rooted(&with(TWO, operator), "#20024"));
    let child = &assembly.nodes[0].children[2];
    assert!(assembly.nodes[child.node].name.is_empty());
    assert_eq!(child.placement.origin, [0., 0., 9.]);
    // An item listed twice would be placed twice.
    let twice = rooted(&with(TWO, &mapping), "#20011,#20011");
    assert_eq!(refused(&twice), Error::Unsupported);
    // A mapped origin that is not the identity reads differently under
    // ISO 10303-43 and OpenCascade.
    let skew = mapping.replace("REPRESENTATION_MAP(#11,", "REPRESENTATION_MAP(#20030,")
        + &axis(20_030, [1., 0., 0.], [0., 0., 1.], [0., 1., 0.]);
    assert_eq!(
        refused(&rooted(&with(TWO, &skew), "#20011")),
        Error::Unsupported
    );
}

#[test]
fn identity_relationships_join_representations_into_one_definition() {
    // The cylinder's product defines a chain of shapes related to its solid.
    let text = swap(
        TWO,
        "#378 = SHAPE_DEFINITION_REPRESENTATION(#379,#385);",
        "#378 = SHAPE_DEFINITION_REPRESENTATION(#379,#20020);",
    );
    let chain = "#20020=SHAPE_REPRESENTATION('',(#11),#484);\
        #20021=SHAPE_REPRESENTATION_RELATIONSHIP('','',#20020,#20022);\
        #20022=SHAPE_REPRESENTATION('',(#11),#484);\
        #20023=(REPRESENTATION_RELATIONSHIP('','',#20022,#385)SHAPE_REPRESENTATION_RELATIONSHIP());\
        #20024=SHAPE_REPRESENTATION_RELATIONSHIP('','',#20020,#2);";
    let assembly = read(&with(&text, chain));
    assert_eq!(assembly.nodes.len(), 3);
    let found = volumes(&assembly);
    assert!(close(found[1], CYLINDER));
}

#[test]
fn transformation_operators_are_rigid_and_refuse_mirrors_and_scales() {
    let placed = |operator: &str| {
        let text = swap(
            TWO,
            "REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION(#491)",
            "REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION(#20030)",
        );
        let directions = "#20031=DIRECTION('',(0.,1.,0.));#20032=DIRECTION('',(-1.,0.,0.));\
            #20033=DIRECTION('',(0.,0.,1.));#20034=CARTESIAN_POINT('',(100.,0.,0.));\
            #20035=DIRECTION('',(1.,0.,0.));";
        let text = with(&text, &format!("{directions}#20030={operator};"));
        read_assembly(&text, TOLERANCE, None).map(|a| a.occurrences()[1].placement)
    };
    let turned =
        placed("CARTESIAN_TRANSFORMATION_OPERATOR_3D('','','',#20031,#20032,#20034,1.,#20033)");
    assert!(same(&turned.unwrap(), &TURNED));
    // Omitted axes follow base_axis.
    let derived = placed("CARTESIAN_TRANSFORMATION_OPERATOR_3D('','','',#20031,$,#20034,$,$)");
    assert!(same(&derived.unwrap(), &TURNED));
    let default =
        placed("CARTESIAN_TRANSFORMATION_OPERATOR_3D('','','',$,$,#20034,$,#20035)").unwrap();
    assert_eq!(default.axes, [[0., 1., 0.], [0., 0., 1.], [1., 0., 0.]]);
    let refusals = [
        // Mirroring: the second axis opposes z cross x.
        ("#20031,#20035,#20034,1.,#20033", Error::Unsupported),
        ("#20031,#20032,#20034,2.,#20033", Error::Unsupported),
        ("#20031,#20032,#20034,'1',#20033", Error::Invalid),
        // The first axis along z leaves no x.
        ("#20033,$,#20034,$,#20033", Error::Invalid),
    ];
    for (args, kind) in refusals {
        let operator = format!("CARTESIAN_TRANSFORMATION_OPERATOR_3D('','','',{args})");
        assert_eq!(placed(&operator).unwrap_err(), kind, "{args}");
    }
    assert_eq!(
        placed("CARTESIAN_TRANSFORMATION_OPERATOR_3D('','',#20031,#20032,#20034,1.,#20033)")
            .unwrap_err(),
        Error::Invalid
    );
    // An operator may also be an item of an item-defined transformation.
    let text = swap(
        TWO,
        "ITEM_DEFINED_TRANSFORMATION('','',#11,#19)",
        "ITEM_DEFINED_TRANSFORMATION('','',#11,#20030)",
    );
    let text = with(
        &rooted(&text, "#20030"),
        "#20031=DIRECTION('',(0.,1.,0.));#20034=CARTESIAN_POINT('',(100.,0.,0.));\
         #20030=CARTESIAN_TRANSFORMATION_OPERATOR_3D('','','',#20031,$,#20034,$,$);",
    );
    assert!(same(&read(&text).occurrences()[1].placement, &TURNED));
}

#[test]
fn item_placements_must_agree_and_belong_to_their_representations() {
    let items = |first: u32, second: u32, extra: &str| {
        let text = swap(
            &repeated(),
            "ITEM_DEFINED_TRANSFORMATION('','',#11,#20000)",
            &format!("ITEM_DEFINED_TRANSFORMATION('','',#{first},#{second})"),
        );
        let text = swap(
            &text,
            "#36 = ADVANCED_BREP_SHAPE_REPRESENTATION('',(#11,#37),#367);",
            "#36 = ADVANCED_BREP_SHAPE_REPRESENTATION('',(#11,#37,#20040),#367);",
        );
        read_assembly(&with(&text, extra), TOLERANCE, None).map(|a| a.occurrences()[2].placement)
    };
    // An origin item shifted along the turning axis reads the same either
    // way; shifted across it, the two readings part.
    let shifted = axis(20_040, [0., 0., 1.], [0., 0., 1.], [1., 0., 0.]);
    let moved = items(20_040, 20_000, &shifted).unwrap();
    let expected = Placement {
        origin: [100., 0., -1.],
        axes: TURNED.axes,
    };
    assert!(same(&moved, &expected), "{moved:?}");
    for across in [
        axis(20_040, [1., 0., 0.], [0., 0., 1.], [1., 0., 0.]),
        axis(20_040, [1., 0., 0.], [0., 0., 1.], [0., 1., 0.]),
    ] {
        assert_eq!(
            items(20_040, 20_000, &across).unwrap_err(),
            Error::Unsupported
        );
    }
    // Each item must be an item of its representation; a solid places
    // nothing.
    assert_eq!(items(11, 15, &shifted).map(|_| ()), Ok(()));
    assert_eq!(items(20_000, 11, &shifted).unwrap_err(), Error::Invalid);
    assert_eq!(items(37, 20_000, &shifted).unwrap_err(), Error::Unsupported);
    // A transformation must be an item-defined one or an operator.
    let text = swap(
        TWO,
        "REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION(#491)",
        "REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION(#12)",
    );
    assert_eq!(refused(&text), Error::Unsupported);
}

#[test]
fn malformed_and_unplaced_structure_is_refused_with_its_kind() {
    let cases = [
        // A product placed inside itself.
        (
            with(TWO, &occurrence(20_000, "", [31, 31], [36, 36], [11, 11])),
            Error::Invalid,
        ),
        // Every product is a component, so nothing is a root.
        (
            with(TWO, &occurrence(20_000, "", [31, 5], [10, 36], [11, 11])),
            Error::Unsupported,
        ),
        // An occurrence without a placement.
        (
            with(
                TWO,
                "#20000=NEXT_ASSEMBLY_USAGE_OCCURRENCE('3','','',#5,#31,$);",
            ),
            Error::Unsupported,
        ),
        // A placement of a product definition, not of an occurrence.
        (
            with(
                TWO,
                "#20000=CONTEXT_DEPENDENT_SHAPE_REPRESENTATION(#373,#4);",
            ),
            Error::Unsupported,
        ),
        // Two placements of one occurrence.
        (
            with(
                TWO,
                "#20000=CONTEXT_DEPENDENT_SHAPE_REPRESENTATION(#490,#375);",
            ),
            Error::Unsupported,
        ),
        // A transformed relationship that places no occurrence.
        (
            with(
                TWO,
                "#20000=(REPRESENTATION_RELATIONSHIP('','',#36,#10)REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION(#374)SHAPE_REPRESENTATION_RELATIONSHIP());",
            ),
            Error::Unsupported,
        ),
        // An occurrence context whose relationship has no transformation.
        (
            swap(
                TWO,
                "REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION(#374) \n",
                "",
            ),
            Error::Invalid,
        ),
        // A dangling definition, and one that is not a product definition.
        (
            swap(TWO, "('1','','',#5,#31,$)", "('1','','',#5,#99999,$)"),
            Error::Invalid,
        ),
        (
            swap(TWO, "('1','','',#5,#31,$)", "('1','','',#5,#4,$)"),
            Error::Unsupported,
        ),
        // A component without a shape.
        (
            with(
                TWO,
                &("#20000=PRODUCT_DEFINITION('design','',#32,#35);".to_owned()
                    + &occurrence(20_010, "", [5, 20_000], [36, 10], [11, 11])),
            ),
            Error::Unsupported,
        ),
        // A relationship between other shapes than the occurrence's.
        (
            swap(TWO, "('','',#36,#10)", "('','',#36,#385)"),
            Error::Invalid,
        ),
        // A shape that is not a representation.
        (
            swap(
                TWO,
                "SHAPE_DEFINITION_REPRESENTATION(#4,#10)",
                "SHAPE_DEFINITION_REPRESENTATION(#4,#2)",
            ),
            Error::Invalid,
        ),
        // Two shapes of one product.
        (
            with(TWO, "#20000=SHAPE_DEFINITION_REPRESENTATION(#30,#385);"),
            Error::Unsupported,
        ),
        // A product definition without a formation of a product.
        (
            swap(
                TWO,
                "PRODUCT_DEFINITION('design','',#32,#35)",
                "PRODUCT_DEFINITION('design','',#33,#35)",
            ),
            Error::Invalid,
        ),
        // An occurrence name that is not a string, or too few attributes.
        (
            swap(TWO, "('1','','',#5,#31,$)", "('1',$,'',#5,#31,$)"),
            Error::Invalid,
        ),
        (
            swap(TWO, "('1','','',#5,#31,$)", "('1','','',#5,#31)"),
            Error::Invalid,
        ),
        // A representation whose context states no units.
        (
            with(
                &swap(TWO, "(#11,#15,#19),#23);", "(#11,#15,#19),#20000);"),
                "#20000=(GEOMETRIC_REPRESENTATION_CONTEXT(3)REPRESENTATION_CONTEXT('',''));",
            ),
            Error::Unsupported,
        ),
        // A component in metres beside millimetres.
        (
            swap(
                TWO,
                "#485 = ( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.) );",
                "#485 = ( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT($,.METRE.) );",
            ),
            Error::Unsupported,
        ),
        // Nesting deeper than the structure budget, or too many placed solids.
        (tree(MAX_DEPTH as u32 + 1, 1), Error::Limit),
        (tree(10, 2), Error::Limit),
    ];
    for (index, (text, kind)) in cases.iter().enumerate() {
        assert_eq!(refused(text), *kind, "case {index}");
    }
    // A shape stated for an occurrence or a property places nothing: the
    // same shape as the product's reads, another is left without a place,
    // and a dangling one is invalid.
    let occurrence_shape = with(TWO, "#20000=SHAPE_DEFINITION_REPRESENTATION(#375,#36);");
    assert_eq!(read(&occurrence_shape).occurrences().len(), 2);
    let property = "#20001=PROPERTY_DEFINITION('','',#5);\
        #20000=SHAPE_DEFINITION_REPRESENTATION(#20001,#20002);\
        #20002=SHAPE_REPRESENTATION('',(#11),#23);";
    assert_eq!(read(&with(TWO, property)).occurrences().len(), 2);
    let elsewhere = swap(
        TWO,
        "#378 = SHAPE_DEFINITION_REPRESENTATION(#379,#385);",
        "#378 = SHAPE_DEFINITION_REPRESENTATION(#492,#385);",
    );
    assert_eq!(refused(&elsewhere), Error::Unsupported);
    for dangling in [
        "(#99999,#385)",
        "(#20000,#385);#20000=PRODUCT_DEFINITION_SHAPE('','',#99999)",
    ] {
        let text = with(
            TWO,
            &format!("#20001=SHAPE_DEFINITION_REPRESENTATION{dangling};"),
        );
        assert_eq!(refused(&text), Error::Invalid, "{dangling}");
    }
    // The cylinder's shape, occurrence and placement removed: its solid is
    // left without a place.
    let mut text = TWO.to_owned();
    for entity in [
        "#378 = SHAPE_DEFINITION_REPRESENTATION(#379,#385);",
        "#489 = CONTEXT_DEPENDENT_SHAPE_REPRESENTATION(#490,#492);",
        "#493 = NEXT_ASSEMBLY_USAGE_OCCURRENCE('2','','',#5,#380,$);",
    ] {
        text = swap(&text, entity, "");
    }
    assert_eq!(refused(&text), Error::Unsupported);
    // A shared definition within the budgets is counted once per use.
    assert_eq!(read(&tree(3, 2)).occurrences().len(), 10);
}

#[test]
fn written_structure_reads_back_with_the_same_definitions_and_placements() {
    let source = read(&repeated());
    let text = write_assembly(&source, Placement::IDENTITY).unwrap();
    let back = read(&text);
    assert_eq!(back.nodes.len(), source.nodes.len());
    let names = |node: &Node<Solid>| {
        node.children
            .iter()
            .map(|c| (c.node, c.name.clone()))
            .collect::<Vec<_>>()
    };
    for (a, b) in back.nodes.iter().zip(&source.nodes) {
        assert_eq!(a.name, b.name);
        assert_eq!(a.solids.len(), b.solids.len());
        assert_eq!(names(a), names(b));
    }
    let (a, b) = (back.occurrences(), source.occurrences());
    assert_eq!(a.len(), b.len());
    for (a, b) in a.iter().zip(&b) {
        assert_eq!((a.node, a.solid, &a.path), (b.node, b.solid, &b.path));
        assert!(same(&a.placement, &b.placement));
    }
    // Moved roots are held by one unnamed product that places them.
    let moved = write_assembly(&source, TURNED).unwrap();
    let wrapped = read(&moved);
    assert_eq!(wrapped.roots.len(), 1);
    let root = &wrapped.nodes[wrapped.roots[0]];
    assert!(root.name.is_empty() && root.children.len() == 1);
    let occurrences = wrapped.occurrences();
    assert!(same(&occurrences[2].placement, &TURNED.then(&TURNED)));
    // Apostrophes survive.
    let quoted = swap(&repeated(), "'again'", "'it''s'");
    let text = write_assembly(&read(&quoted), Placement::IDENTITY).unwrap();
    assert_eq!(read(&text).nodes[0].children[2].name, "it's");
}

#[test]
fn a_structure_beyond_the_readback_budget_is_not_written() {
    let mut source = read(TWO);
    let child = source.nodes[0].children[0].clone();
    source.nodes[0].children = vec![child; 120_000];
    assert_eq!(
        write_assembly(&source, Placement::IDENTITY).unwrap_err(),
        crate::step::StepError::Limit
    );
}
