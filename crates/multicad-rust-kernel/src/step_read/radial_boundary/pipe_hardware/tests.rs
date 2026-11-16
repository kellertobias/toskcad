use super::*;
use crate::step::pipe_hardware::{self as writer, test_fixture::engineering};
use crate::step_read::literal_member_binding::{Token, token};
use std::ops::Range;

struct Damage {
    label: &'static str,
    range: Range<usize>,
    replacement: String,
    expected: GeometryError,
}

// Only called on original components already accepted by the real neutral
// parser. This cursor locates literal fields to alter; it grants no admission.
fn arguments(source: &str) -> Vec<Range<usize>> {
    let mut rest = source;
    let mut depth = 0;
    let mut start = 0;
    let mut fields = Vec::new();
    while let Some(next) = token(&mut rest) {
        if let Token::Literal(value) = next {
            let offset = value.as_ptr() as usize - source.as_ptr() as usize;
            match value {
                "(" => {
                    depth += 1;
                    if depth == 1 {
                        start = source.len() - rest.len();
                    }
                }
                ")" => {
                    depth -= 1;
                    if depth == 0 {
                        fields.push(start..offset);
                        break;
                    }
                }
                "," if depth == 1 => {
                    fields.push(start..offset);
                    start = source.len() - rest.len();
                }
                _ => {}
            }
        }
    }
    assert_eq!(depth, 0);
    fields
}
fn named<'a, 'i>(
    arena: &'a BTreeMap<u32, Vec<part21::Entity<'i>>>,
    name: &str,
) -> &'a part21::Entity<'i> {
    arena.values().flatten().find(|e| e.name == name).unwrap()
}
fn record_range(entity: &part21::Entity<'_>, input: &str) -> Range<usize> {
    let source = entity.source.unwrap();
    let start = source.as_ptr() as usize - input.as_ptr() as usize;
    assert_eq!(input.get(start..start + source.len()), Some(source));
    start..start + source.len()
}
fn fields(entity: &part21::Entity<'_>, input: &str) -> Vec<Range<usize>> {
    let body = record_range(entity, input);
    let args = arguments(&input[body.clone()]);
    assert_eq!(args.len(), entity.args.len());
    args.into_iter()
        .map(|r| body.start + r.start..body.start + r.end)
        .collect()
}

#[test]
fn physical_malformed_geometry_and_unsupported_curve_families_cannot_earn_material() {
    let source = engineering();
    let candidate = writer::propose(&source, writer::Options::default(), None).unwrap();
    let path =
        std::env::temp_dir().join(format!("codec-schema-corpus-{}.step", std::process::id()));
    std::fs::write(&path, candidate.as_str()).unwrap();
    drop(candidate);
    drop(source);
    let actual = std::fs::read_to_string(&path).unwrap();
    let cases = {
        let arena = part21::parse(&actual).unwrap();
        let spline = fields(named(&arena, "B_SPLINE_CURVE_WITH_KNOTS"), &actual);
        let tensor = fields(named(&arena, "B_SPLINE_SURFACE_WITH_KNOTS"), &actual);
        let line = named(&arena, "LINE");
        let vector = &arena[&reference(&line.args[2]).unwrap()][0];
        let direction = &arena[&reference(&vector.args[1]).unwrap()][0];
        let line = fields(line, &actual);
        let vector = fields(vector, &actual);
        let direction = fields(direction, &actual);
        let circle = named(&arena, "CIRCLE");
        let circle_fields = fields(circle, &actual);
        let circle_range = record_range(circle, &actual);
        let surface_rows = arguments(&actual[tensor[3].clone()]);
        assert!(surface_rows.len() > 2);
        let second = tensor[3].start + surface_rows[1].start..tensor[3].start + surface_rows[1].end;
        let shell = fields(named(&arena, "CLOSED_SHELL"), &actual);
        let shell_members = arguments(&actual[shell[1].clone()]);
        assert_eq!(shell_members.len(), 13);
        let fewer = shell_members[..12]
            .iter()
            .map(|r| &actual[shell[1].start + r.start..shell[1].start + r.end])
            .collect::<Vec<_>>()
            .join(",");
        let face_record = named(&arena, "ADVANCED_FACE");
        let first_bound = &arena[&refs(&face_record.args[1]).unwrap()[0]][0];
        let actual_loop = &arena[&reference(&first_bound.args[1]).unwrap()][0];
        let oriented = refs(&actual_loop.args[1]).unwrap()[0];
        let repeated = (0..81)
            .map(|_| format!("#{oriented}"))
            .collect::<Vec<_>>()
            .join(",");
        let loop_fields = fields(actual_loop, &actual);
        let edge_record = named(&arena, "EDGE_CURVE");
        let vertex_record = &arena[&reference(&edge_record.args[1]).unwrap()][0];
        assert_eq!(vertex_record.name, "VERTEX_POINT");
        let referenced_point = &arena[&reference(&vertex_record.args[1]).unwrap()][0];
        assert_eq!(referenced_point.name, "CARTESIAN_POINT");
        let point_fields = fields(referenced_point, &actual);
        let edge_fields = fields(edge_record, &actual);
        let inner = named(&arena, "FACE_BOUND");
        let inner_name = record_range(inner, &actual).start
            ..record_range(inner, &actual).start + inner.name.len();
        let face = fields(face_record, &actual);
        let mut out = Vec::new();
        let mut add = |label, range, replacement: &str, expected| {
            out.push(Damage {
                label,
                range,
                replacement: replacement.to_owned(),
                expected,
            })
        };
        add(
            "spline closed flag",
            spline[4].clone(),
            ".T.",
            GeometryError::Invalid,
        );
        add(
            "spline missing poles",
            spline[2].clone(),
            "$",
            GeometryError::Invalid,
        );
        add(
            "spline multiplicities scalar",
            spline[6].clone(),
            "$",
            GeometryError::Invalid,
        );
        add(
            "spline knots scalar",
            spline[7].clone(),
            "$",
            GeometryError::Invalid,
        );
        add(
            "spline multiplicity cardinality",
            spline[6].clone(),
            "(10,10)",
            GeometryError::Unsupported,
        );
        add(
            "spline knot cardinality",
            spline[7].clone(),
            "(0,1)",
            GeometryError::Unsupported,
        );
        add(
            "tensor closed flag",
            tensor[5].clone(),
            ".T.",
            GeometryError::Invalid,
        );
        add(
            "tensor missing rows",
            tensor[3].clone(),
            "$",
            GeometryError::Invalid,
        );
        add(
            "tensor second row scalar",
            second.clone(),
            "$",
            GeometryError::Invalid,
        );
        add(
            "tensor second row width",
            second,
            "(#1)",
            GeometryError::Invalid,
        );
        add(
            "line record schema",
            line[0].clone(),
            "$",
            GeometryError::Invalid,
        );
        add(
            "vector missing magnitude",
            vector[2].clone(),
            "$",
            GeometryError::Invalid,
        );
        add(
            "direction record schema",
            direction[0].clone(),
            "$",
            GeometryError::Invalid,
        );
        add(
            "line direction changes actual endpoints",
            direction[1].clone(),
            "(0.,0.,0.)",
            GeometryError::Unsupported,
        );
        add(
            "ellipse family unsupported",
            circle_range,
            &format!("ELLIPSE('',{},1.)", &actual[circle_fields[1].clone()]),
            GeometryError::Unsupported,
        );
        add(
            "missing actual shell face",
            shell[1].clone(),
            &format!("({fewer})"),
            GeometryError::Unsupported,
        );
        add(
            "empty actual face boundary",
            face[1].clone(),
            "()",
            GeometryError::Unsupported,
        );
        add(
            "actual referenced vertex exceeds coordinate domain",
            point_fields[1].clone(),
            "(1000000000001.,0.,0.)",
            GeometryError::Limit,
        );
        add(
            "edge curve record schema",
            edge_fields[0].clone(),
            "$",
            GeometryError::Invalid,
        );
        add(
            "duplicate shoulder outer bound",
            inner_name,
            "FACE_OUTER_BOUND",
            GeometryError::Invalid,
        );
        add(
            "actual edge loop exceeds finite incidence cap",
            loop_fields[1].clone(),
            &format!("({repeated})"),
            GeometryError::Limit,
        );
        let record_end = record_range(referenced_point, &actual).end;
        add(
            "physical syntax truncation crosses referenced point record",
            record_end - 1..record_end,
            "",
            GeometryError::Invalid,
        );
        let insert = actual.rfind("ENDSEC;").unwrap();
        add(
            "typed auxiliary subULP numeric alias",
            insert..insert,
            "#999999=CARTESIAN_POINT('',(LENGTH_MEASURE(1.0000000000000000000001),0.,0.));\n",
            GeometryError::Unsupported,
        );
        out
    }; // The genuine original neutral arena and all range scaffolding die here.
    let resident = actual.capacity()
        + cases.capacity() * size_of::<Damage>()
        + cases
            .iter()
            .map(|c| c.replacement.capacity())
            .sum::<usize>()
        + 4096;
    for case in &cases {
        let mut changed = actual.clone();
        changed.replace_range(case.range.clone(), &case.replacement);
        assert!(
            changed != actual,
            "physical mutation must alter {}",
            case.label
        );
        std::fs::write(&path, &changed).unwrap();
        drop(changed);
        let damaged = std::fs::read_to_string(&path).unwrap();
        let result = read(
            &damaged,
            1e-7,
            1e-6,
            Options {
                additional_resident_bytes: resident + damaged.capacity() - damaged.len(),
                ..Default::default()
            },
            None,
        );
        assert_eq!(
            result.err(),
            Some(case.expected),
            "physical malformed case: {}",
            case.label
        );
    }
    drop(cases);
    std::fs::write(&path, &actual).unwrap();
    let healthy = std::fs::read_to_string(&path).unwrap();
    std::fs::remove_file(&path).unwrap();
    let retry = read(
        &healthy,
        1e-7,
        1e-6,
        Options {
            additional_resident_bytes: actual.capacity() + healthy.capacity() - healthy.len()
                + 4096,
            ..Default::default()
        },
        None,
    )
    .unwrap();
    assert_eq!(retry.original_input(), actual);
    assert_eq!(retry.material().head_height_mm(), [6., 9.]);
    assert_eq!(
        retry
            .material()
            .capped_thread()
            .classify([0., 0., 3.], None)
            .unwrap(),
        crate::radial_tensor::Location::Inside
    );
}

pub(crate) struct InteriorChange {
    pub text: String,
    pub tensor_entity: u32,
    pub point_entity: u32,
    pub pole_index: usize,
    pub old_value: f64,
    pub new_value: f64,
}
// Literal preparation only. The changed document must independently earn actual
// occupied material in the caller; this structure is not a Source certificate.
pub(crate) fn exclusive_interior_nextup(input: &str) -> InteriorChange {
    fn count(value: &Value<'_>, counts: &mut BTreeMap<u32, usize>) {
        match value {
            Value::Reference(id) => *counts.entry(*id).or_default() += 1,
            Value::List(values) => {
                for v in values {
                    count(v, counts);
                }
            }
            Value::Typed(e) => {
                for v in &e.args {
                    count(v, counts);
                }
            }
            _ => {}
        }
    }
    let arena = part21::parse(input).unwrap();
    let mut counts = BTreeMap::new();
    for e in arena.values().flatten() {
        for v in &e.args {
            count(v, &mut counts);
        }
    }
    let mut chosen = None;
    for (&tensor_id, components) in &arena {
        for e in components
            .iter()
            .filter(|e| e.name == "B_SPLINE_SURFACE_WITH_KNOTS")
        {
            assert_eq!(e.args[1], Value::Number(9.));
            assert_eq!(e.args[2], Value::Number(9.));
            if let Value::List(rows) = &e.args[3] {
                for (i, row) in rows.iter().enumerate() {
                    // Exclude boundary and first derivative/seam controls in
                    // EVERY degree9 cell, not merely the global grid corners.
                    if !(2..=7).contains(&(i % 9)) {
                        continue;
                    }
                    if let Value::List(poles) = row {
                        for (j, pole) in poles.iter().enumerate() {
                            if !(2..=7).contains(&(j % 9)) {
                                continue;
                            }
                            let id = reference(pole).unwrap();
                            if counts[&id] == 1 {
                                chosen = Some((tensor_id, id, i * poles.len() + j, i, j));
                                break;
                            }
                        }
                    }
                    if chosen.is_some() {
                        break;
                    }
                }
            }
            if chosen.is_some() {
                break;
            }
        }
        if chosen.is_some() {
            break;
        }
    }
    let (tensor_entity, point_entity, pole_index, row, col) =
        chosen.expect("exclusive actual interior tensor pole");
    assert_eq!(
        counts[&point_entity], 1,
        "the actual ID occurs only in that tensor interior, never a boundary curve, vertex, seam, cap or another tensor"
    );
    let point = &arena[&point_entity][0];
    assert_eq!(point.name, "CARTESIAN_POINT");
    let coords = fields(point, input)[1].clone();
    let values = arguments(&input[coords.clone()]);
    assert_eq!(values.len(), 3);
    let scalar = coords.start + values[0].start..coords.start + values[0].end;
    let old_value = input[scalar.clone()].parse::<f64>().unwrap();
    let new_value = old_value.next_up();
    assert!(old_value.is_finite() && new_value.is_finite() && new_value > old_value);
    let literal = writer::numbers::literal(new_value).unwrap();
    let mut text = input.to_owned();
    text.replace_range(scalar, &literal);
    assert!(text != input);
    eprintln!(
        "EXCLUSIVE_INTERIOR tensor{tensor_entity} point{point_entity} row{row} col{col} pole{pole_index} oldbits{} newbits{} actual_references1",
        old_value.to_bits(),
        new_value.to_bits()
    );
    InteriorChange {
        text,
        tensor_entity,
        point_entity,
        pole_index,
        old_value,
        new_value,
    }
}
