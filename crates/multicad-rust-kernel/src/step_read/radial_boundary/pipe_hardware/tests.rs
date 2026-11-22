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
        // Genuine distinct face/edge records retain their original incidence;
        // substituting another actual support must not grant occupied material.
        let tensor_faces: Vec<_> = arena
            .values()
            .flatten()
            .filter(|e| {
                e.name == "ADVANCED_FACE"
                    && arena[&reference(&e.args[2]).unwrap()][0].name
                        == "B_SPLINE_SURFACE_WITH_KNOTS"
            })
            .collect();
        assert!(tensor_faces.len() > 1);
        let support0 = reference(&tensor_faces[0].args[2]).unwrap();
        let support1 = reference(&tensor_faces[1].args[2]).unwrap();
        assert_ne!(support0, support1);
        let aliased_face = fields(tensor_faces[1], &actual)[2].clone();
        let spline_edges: Vec<_> = arena
            .values()
            .flatten()
            .filter(|e| {
                e.name == "EDGE_CURVE"
                    && arena[&reference(&e.args[3]).unwrap()][0].name == "B_SPLINE_CURVE_WITH_KNOTS"
            })
            .collect();
        assert!(spline_edges.len() > 1);
        let curve0 = reference(&spline_edges[0].args[3]).unwrap();
        let other_edge = spline_edges
            .iter()
            .find(|e| reference(&e.args[3]).unwrap() != curve0)
            .unwrap();
        assert_ne!(
            [
                reference(&spline_edges[0].args[1]).unwrap(),
                reference(&spline_edges[0].args[2]).unwrap()
            ],
            [
                reference(&other_edge.args[1]).unwrap(),
                reference(&other_edge.args[2]).unwrap()
            ]
        );
        let aliased_edge = fields(other_edge, &actual)[3].clone();
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
        add(
            "distinct tensor face borrows another actual support ID",
            aliased_face,
            &format!("#{support0}"),
            GeometryError::Unresolved,
        );
        add(
            "distinct spline edge borrows another actual curve ID",
            aliased_edge,
            &format!("#{curve0}"),
            GeometryError::Invalid,
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
// Copy actual physical records under fresh IDs. The resulting graphs are
// deliberately over the finite decoder frontier, not earned material objects.
fn physical_incidence_frontier(actual: &str, unique_vertices: bool) -> String {
    use std::fmt::Write as _;
    let arena = part21::parse(actual).unwrap();
    let loops: Vec<_> = arena
        .values()
        .flatten()
        .filter(|e| e.name == "EDGE_LOOP" && refs(&e.args[1]).unwrap().len() == 4)
        .take(if unique_vertices { 1 } else { 2 })
        .collect();
    assert_eq!(loops.len(), if unique_vertices { 1 } else { 2 });
    let mut next = *arena.keys().next_back().unwrap() + 1;
    let mut records = String::new();
    let mut replacements = Vec::new();
    let mut copied_edges = 0;
    let mut copied_vertices = 0;
    for (index, loop_record) in loops.into_iter().enumerate() {
        let original = refs(&loop_record.args[1]).unwrap();
        let length = if unique_vertices {
            52
        } else if index == 0 {
            80
        } else {
            52
        };
        assert!(length <= 80 && length % original.len() == 0);
        let mut oriented_ids = Vec::new();
        for offset in 0..length {
            let oriented = &arena[&original[offset % original.len()]][0];
            assert_eq!(oriented.name, "ORIENTED_EDGE");
            let edge = &arena[&reference(&oriented.args[3]).unwrap()][0];
            assert_eq!(edge.name, "EDGE_CURVE");
            let mut body = edge.source.unwrap().to_owned();
            let args = arguments(&body);
            if unique_vertices {
                for k in [2, 1] {
                    let vertex = &arena[&reference(&edge.args[k]).unwrap()][0];
                    assert_eq!(vertex.name, "VERTEX_POINT");
                    let point = reference(&vertex.args[1]).unwrap();
                    assert_eq!(arena[&point][0].name, "CARTESIAN_POINT");
                    writeln!(records, "#{next}=VERTEX_POINT('',#{point});").unwrap();
                    body.replace_range(args[k].clone(), &format!("#{next}"));
                    next += 1;
                    copied_vertices += 1;
                }
            }
            let edge_id = next;
            next += 1;
            writeln!(records, "#{edge_id}={body};").unwrap();
            copied_edges += 1;
            let mut body = oriented.source.unwrap().to_owned();
            body.replace_range(arguments(&body)[3].clone(), &format!("#{edge_id}"));
            writeln!(records, "#{next}={body};").unwrap();
            oriented_ids.push(format!("#{next}"));
            next += 1;
        }
        replacements.push((
            fields(loop_record, actual)[1].clone(),
            format!("({})", oriented_ids.join(",")),
        ));
    }
    assert_eq!(copied_edges, if unique_vertices { 52 } else { 132 });
    assert_eq!(copied_vertices, if unique_vertices { 104 } else { 0 });
    let mut changed = actual.to_owned();
    changed.insert_str(actual.rfind("ENDSEC;").unwrap(), &records);
    replacements.sort_by_key(|(range, _)| std::cmp::Reverse(range.start));
    for (range, value) in replacements {
        changed.replace_range(range, &value);
    }
    // Independently parse the actual resulting file and prove the added IDs
    // really own the advertised distinct vertex/edge records.
    let parsed = part21::parse(&changed).unwrap();
    for (kind, added) in [
        ("EDGE_CURVE", copied_edges),
        ("VERTEX_POINT", copied_vertices),
    ] {
        let count = |map: &BTreeMap<u32, Vec<part21::Entity<'_>>>| {
            map.values().flatten().filter(|e| e.name == kind).count()
        };
        assert_eq!(count(&parsed), count(&arena) + added);
    }
    changed
}

#[test]
fn physical_graph_vertex_and_edge_frontiers_refuse_before_material_and_retry() {
    let source = engineering();
    let candidate = writer::propose(&source, writer::Options::default(), None).unwrap();
    let path = std::env::temp_dir().join(format!(
        "codec-incidence-frontier-{}.step",
        std::process::id()
    ));
    std::fs::write(&path, candidate.as_str()).unwrap();
    drop(candidate);
    drop(source);
    let actual = std::fs::read_to_string(&path).unwrap();
    for unique_vertices in [false, true] {
        let changed = physical_incidence_frontier(&actual, unique_vertices);
        assert_ne!(changed, actual);
        std::fs::write(&path, &changed).unwrap();
        drop(changed); // All construction arenas and buffers die before read.
        let damaged = std::fs::read_to_string(&path).unwrap();
        assert_eq!(
            read(
                &damaged,
                1e-7,
                1e-6,
                Options {
                    additional_resident_bytes: actual.capacity() + damaged.capacity()
                        - damaged.len()
                        + 4096,
                    ..Default::default()
                },
                None
            )
            .err(),
            Some(GeometryError::Limit),
            "unique vertices={unique_vertices}"
        );
    }
    std::fs::write(&path, &actual).unwrap();
    let restored = std::fs::read_to_string(&path).unwrap();
    std::fs::remove_file(&path).unwrap();
    let retry = read(
        &restored,
        1e-7,
        1e-6,
        Options {
            additional_resident_bytes: actual.capacity() + restored.capacity() - restored.len()
                + 4096,
            ..Default::default()
        },
        None,
    )
    .unwrap();
    assert_eq!(retry.original_input(), actual);
    assert_eq!(
        retry
            .material()
            .capped_thread()
            .classify([0., 0., 3.], None)
            .unwrap(),
        radial::Location::Inside
    );
}

// Physical preparation only: clone a real point under a fresh ID so changing
// its curve use cannot silently change interned tensor/vertex/seam owners.
pub(crate) fn exclusive_curve_interior_nextup(input: &str) -> String {
    let arena = part21::parse(input).unwrap();
    let (&curve_id, components) = arena
        .iter()
        .find(|(_, components)| {
            components.len() == 1 && components[0].name == "B_SPLINE_CURVE_WITH_KNOTS"
        })
        .unwrap();
    let curve = &components[0];
    assert_eq!(curve.args[1], Value::Number(9.));
    let poles = refs(&curve.args[2]).unwrap();
    let pole_index = 2;
    assert!(poles.len() >= 10 && (2..=7).contains(&(pole_index % 9)));
    let old_point_id = poles[pole_index];
    let point = &arena[&old_point_id][0];
    assert_eq!(point.name, "CARTESIAN_POINT");
    let coords = fields(point, input)[1].clone();
    let values = arguments(&input[coords.clone()]);
    let scalar = coords.start + values[0].start..coords.start + values[0].end;
    let old_value = input[scalar.clone()].parse::<f64>().unwrap();
    let new_value = old_value.next_up();
    assert!(old_value.is_finite() && new_value.is_finite() && new_value > old_value);
    let new_id = *arena.keys().next_back().unwrap() + 1;
    assert!(!arena.contains_key(&new_id));
    let body_range = record_range(point, input);
    let mut new_body = input[body_range.clone()].to_owned();
    new_body.replace_range(
        scalar.start - body_range.start..scalar.end - body_range.start,
        &writer::numbers::literal(new_value).unwrap(),
    );
    let pole_list = fields(curve, input)[2].clone();
    let pole_fields = arguments(&input[pole_list.clone()]);
    assert_eq!(pole_fields.len(), poles.len());
    let reference_range = pole_list.start + pole_fields[pole_index].start
        ..pole_list.start + pole_fields[pole_index].end;
    assert_eq!(input[reference_range.clone()], format!("#{old_point_id}"));
    let mut text = input.to_owned();
    text.insert_str(
        input.rfind("ENDSEC;").unwrap(),
        &format!("#{new_id}={new_body};\n"),
    );
    text.replace_range(reference_range, &format!("#{new_id}"));
    let changed = part21::parse(&text).unwrap();
    assert_eq!(changed.len(), arena.len() + 1);
    for (&id, original) in &arena {
        if id != curve_id {
            assert_eq!(
                original.iter().map(|e| e.source).collect::<Vec<_>>(),
                changed[&id].iter().map(|e| e.source).collect::<Vec<_>>()
            );
        }
    }
    let actual_curve = &changed[&curve_id][0];
    assert_eq!(curve.name, actual_curve.name);
    for (k, (a, b)) in curve.args.iter().zip(&actual_curve.args).enumerate() {
        if k != 2 {
            assert_eq!(a, b);
        }
    }
    let after = refs(&actual_curve.args[2]).unwrap();
    assert_eq!(after.len(), poles.len());
    for (i, (a, b)) in poles.iter().zip(&after).enumerate() {
        assert_eq!(*b, if i == pole_index { new_id } else { *a });
    }
    eprintln!(
        "EXCLUSIVE_CURVE_INTERIOR curve{curve_id} originalpoint{old_point_id} freshpoint{new_id} pole{pole_index} oldbits{} newbits{} actual_curve_ref_changes1 all_other_records_unchangedtrue",
        old_value.to_bits(),
        new_value.to_bits()
    );
    text
}

#[test]
fn physical_exclusive_curve_control_refusal_keeps_original_source_usable() {
    let source = engineering();
    let before = source.shell().input().clone();
    let snapshot_budget = source.receipt().retained_result_bytes_upper_bound;
    let candidate = writer::propose(
        &source,
        writer::Options {
            additional_resident_bytes: snapshot_budget,
            ..Default::default()
        },
        None,
    )
    .unwrap();
    let path =
        std::env::temp_dir().join(format!("codec-curve-control-{}.step", std::process::id()));
    std::fs::write(&path, candidate.as_str()).unwrap();
    drop(candidate);
    let original = std::fs::read_to_string(&path).unwrap();
    let changed = exclusive_curve_interior_nextup(&original);
    std::fs::write(&path, &changed).unwrap();
    drop(changed);
    let changed = std::fs::read_to_string(&path).unwrap();
    // This actual changed curve is refused before material admission. Never
    // pass it to the independent correspondence helper as an earned source.
    assert_eq!(
        read(
            &changed,
            1e-7,
            1e-6,
            Options {
                additional_resident_bytes: source.receipt().retained_result_bytes_upper_bound
                    + snapshot_budget
                    + original.capacity()
                    + changed.capacity()
                    - changed.len()
                    + 4096,
                ..Default::default()
            },
            None
        )
        .err(),
        Some(GeometryError::Unresolved)
    );
    assert_eq!(source.shell().input(), &before);
    drop(before);
    drop(changed);
    std::fs::write(&path, &original).unwrap();
    let restored = std::fs::read_to_string(&path).unwrap();
    std::fs::remove_file(&path).unwrap();
    let retry = writer::qualify_actual_file(
        &source,
        &restored,
        1e-7,
        1e-6,
        writer::QualificationOptions {
            additional_resident_bytes: original.capacity() + restored.capacity() - restored.len()
                + 4096,
            ..Default::default()
        },
        None,
    )
    .unwrap();
    assert_eq!(retry.as_str(), original);
    assert_eq!(
        retry
            .material()
            .material()
            .capped_thread()
            .classify([0., 0., 3.], None)
            .unwrap(),
        radial::Location::Inside
    );
}

pub(crate) fn actual_cone_first_face_order(input: &str) -> String {
    let arena = part21::parse(input).unwrap();
    let (&shell_id, shell) = arena
        .iter()
        .find(|(_, c)| c.len() == 1 && c[0].name == "CLOSED_SHELL")
        .unwrap();
    let shell = &shell[0];
    let original = refs(&shell.args[1]).unwrap();
    assert_eq!(original.len(), 13);
    let cone_indices: Vec<_> = original
        .iter()
        .enumerate()
        .filter_map(|(i, id)| {
            let face = &arena[id][0];
            assert_eq!(face.name, "ADVANCED_FACE");
            (arena[&reference(&face.args[2]).unwrap()][0].name == "CONICAL_SURFACE").then_some(i)
        })
        .collect();
    let [index] = cone_indices.as_slice() else {
        panic!("one actual cone face");
    };
    assert!(*index > 0);
    let mut ordered = original.clone();
    let cone = ordered.remove(*index);
    ordered.insert(0, cone);
    assert_eq!(
        ordered.iter().collect::<std::collections::BTreeSet<_>>(),
        original.iter().collect::<std::collections::BTreeSet<_>>()
    );
    let mut text = input.to_owned();
    text.replace_range(
        fields(shell, input)[1].clone(),
        &format!(
            "({})",
            ordered
                .iter()
                .map(|id| format!("#{id}"))
                .collect::<Vec<_>>()
                .join(",")
        ),
    );
    assert_ne!(text, input);
    let changed = part21::parse(&text).unwrap();
    assert_eq!(changed.len(), arena.len());
    for (&id, original) in &arena {
        if id != shell_id {
            assert_eq!(
                original.iter().map(|e| e.source).collect::<Vec<_>>(),
                changed[&id].iter().map(|e| e.source).collect::<Vec<_>>()
            );
        }
    }
    assert_eq!(refs(&changed[&shell_id][0].args[1]).unwrap(), ordered);
    eprintln!(
        "ACTUAL_CONE_FIRST shell{shell_id} cone_face{cone} moved_from{index} faces13 all_other_records_unchangedtrue"
    );
    text
}

pub(crate) struct CircleChange {
    pub text: String,
    pub entity: u32,
    pub old_radius: f64,
    pub new_radius: f64,
}
// Only actual physical bytes change. This value carries no geometry admission.
pub(crate) fn actual_circle_radius_nextup(input: &str) -> CircleChange {
    let arena = part21::parse(input).unwrap();
    let (&id, components) = arena
        .iter()
        .find(|(_, c)| c.len() == 1 && c[0].name == "CIRCLE")
        .unwrap();
    let circle = &components[0];
    let Value::Number(old_radius) = circle.args[2] else {
        panic!("actual radius literal");
    };
    let new_radius = old_radius.next_up();
    assert!(old_radius > 0. && new_radius.is_finite() && new_radius > old_radius);
    let mut text = input.to_owned();
    text.replace_range(
        fields(circle, input)[2].clone(),
        &writer::numbers::literal(new_radius).unwrap(),
    );
    assert_ne!(text, input);
    let changed = part21::parse(&text).unwrap();
    assert_eq!(arena.len(), changed.len());
    for (&other, before) in &arena {
        if other != id {
            assert_eq!(
                before.iter().map(|e| e.source).collect::<Vec<_>>(),
                changed[&other].iter().map(|e| e.source).collect::<Vec<_>>()
            );
        }
    }
    assert_eq!(circle.args[0..2], changed[&id][0].args[0..2]);
    assert_eq!(changed[&id][0].args[2], Value::Number(new_radius));
    eprintln!(
        "ACTUAL_CIRCLE_RADIUS entity{id} oldbits{} newbits{} all_other_records_unchangedtrue",
        old_radius.to_bits(),
        new_radius.to_bits()
    );
    CircleChange {
        text,
        entity: id,
        old_radius,
        new_radius,
    }
}
