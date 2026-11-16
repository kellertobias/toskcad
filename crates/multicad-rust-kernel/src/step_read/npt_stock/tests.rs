use super::*;
mod guards;
use crate::axial_pipe_thread::authored::{AuthoredSection, AuthoredSpec, AuthoredThread};
use crate::step::npt_stock::{CurveOptions, CurvePlan, SurfacePlan, write_source_candidate};
use crate::step_read::owned_insert::DocumentOptions;
pub(super) fn fixture(quarter: bool) -> String {
    let (tpi, e1, l1, l2, od, max) = if quarter {
        (18., 0.49163, 0.2278, 0.4018, 0.540, 0.088)
    } else {
        (27., 0.37360, 0.1615, 0.2639, 0.405, 0.096)
    };
    let pitch = 25.4 / tpi;
    let k = 1. / 32.;
    let m = 3f64.sqrt();
    let slopes = 1. / (m - k) + 1. / (m + k);
    let depth = pitch / slopes - (0.033 + max) * pitch;
    let source = AuthoredThread::new(
        AuthoredSpec {
            major_radius_mm: (e1 * 25.4 - 2. * k * l1 * 25.4 + depth) / 2.,
            pitch_mm: pitch,
            length_mm: l2 * 25.4,
            radial_taper: k,
            section: AuthoredSection::NptFlat {
                depth_mm: depth,
                flat_width_mm: (pitch - depth * slopes) / 2.,
            },
            crest_radius_limit_mm: Some(od * 25.4 / 2.),
        },
        1e-7,
        1e-6,
    )
    .unwrap()
    .stock_clip_plan(Default::default(), None)
    .unwrap();
    let options = CurveOptions::default();
    let curves = CurvePlan::new(&source, options, None).unwrap();
    let plan = SurfacePlan::new(curves, options.into(), None).unwrap();
    write_source_candidate(&plan, options.maximum_work, options.maximum_bytes, None)
        .unwrap()
        .into_string()
}

#[test]
fn public_literal_npt_import_refuses_changed_units_and_exact_noninteger_dimensions() {
    use crate::step_read::{ReadError, SourceImport, read_source_compound};
    let original = fixture(false);
    let before = original.clone();
    for (find, replace, expected) in [
        (
            "SI_UNIT(.MILLI.,.METRE.)",
            "SI_UNIT(.CENTI.,.METRE.)",
            ReadError::Unsupported,
        ),
        (
            "GEOMETRIC_REPRESENTATION_CONTEXT(3)",
            "GEOMETRIC_REPRESENTATION_CONTEXT(4)",
            ReadError::Unsupported,
        ),
        (
            "GEOMETRIC_REPRESENTATION_CONTEXT(3)",
            "GEOMETRIC_REPRESENTATION_CONTEXT(3.000000000000000000001)",
            ReadError::Invalid,
        ),
    ] {
        let changed = original.replacen(find, replace, 1);
        assert_ne!(
            changed, original,
            "the actual unit/context record must change"
        );
        assert!(matches!(read_source_compound(&changed, 1e-7), Err(error) if error == expected));
    }
    assert_eq!(original, before);
    assert_eq!("3.000000000000000000001".parse::<f64>().unwrap(), 3.);
    // The unchanged millimetre input must still import, retaining actual
    // original records and independently bounded material after its input dies.
    let material = {
        let mut imported = read_source_compound(&original, 1e-7).unwrap();
        assert_eq!(imported.len(), 1);
        let SourceImport::LiteralNpt(material) = imported.pop().unwrap() else {
            panic!("the original literal NPT file lost its owned source")
        };
        material
    };
    drop(original);
    assert!(material.point_count() > 0 && material.source_entity_count() > material.point_count());
    assert!(material.local_volume_mm3()[0] > 470. && material.local_volume_mm3()[1] < 471.);
    let point = material.point(0).unwrap();
    assert!(material.original_component(point.entity(), 0).is_some());
}

#[test]
fn actual_packet_parses_numeric_leaves_coupled_boundary_controls_and_closed_owning_graph() {
    for quarter in [false, true] {
        let text = fixture(quarter);
        let document = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
        let solid = document.represented_members()[0].item_entity();
        let parsed = Parsed::parse(&document, solid, Default::default(), None).unwrap();
        eprintln!(
            "NPT DECODE quarter={quarter} work={} querybytes={} documentbytes={}",
            parsed.work_used(),
            parsed.query_budget_bytes(),
            document.retained_budget_bytes()
        );
        assert!(std::ptr::eq(parsed.document(), &document));
        assert_eq!(parsed.solid_entity(), solid);
        assert_eq!(parsed.vertices().len(), 12);
        assert_eq!(parsed.curves().len(), 20);
        assert_eq!(parsed.tensors().len(), 8);
        assert_eq!(parsed.edges().len(), 20);
        assert_eq!(parsed.faces().len(), 10);
        for point in parsed.points() {
            assert!(point.entity() > 0);
            for scalar in point.coordinates() {
                assert!(!scalar.spelling().is_empty());
                let b = scalar.bounds();
                assert!(b[0] <= scalar.value() && scalar.value() <= b[1]);
            }
        }
        for vertex in parsed.vertices() {
            assert!(vertex.entity() > 0);
            assert!(vertex.point_index() < parsed.points().len());
        }
        for curve in parsed.curves() {
            assert!(curve.entity() > 0);
            assert_eq!(curve.point_indices().len(), 11 * curve.spans() + 1);
        }
        for tensor in parsed.tensors() {
            assert!(tensor.entity() > 0);
            assert_eq!(
                tensor.point_indices().len(),
                tensor.spans().map(|n| 11 * n + 1).into_iter().product()
            );
        }
        for edge in parsed.edges() {
            assert!(edge.entity() > 0);
            assert!(edge.curve_index() < parsed.curves().len());
            assert!(
                edge.vertex_indices()
                    .iter()
                    .all(|v| *v < parsed.vertices().len())
            );
        }
        let mut caps = 0;
        for face in parsed.faces() {
            assert!(face.entity() > 0);
            for coedge in face.coedges() {
                assert!(coedge.edge_index() < 20);
                let _ = coedge.forward();
            }
            match face.support() {
                Support::Tensor(index) => {
                    assert!(index < 8);
                    assert!(!face.forward());
                    let mut sides = face.tensor_sides().unwrap();
                    sides.sort_unstable();
                    assert_eq!(sides, [0, 1, 2, 3]);
                }
                Support::Plane { height } => {
                    caps += 1;
                    assert!(height < parsed.points().len());
                    assert_eq!(face.tensor_sides(), None);
                }
            }
        }
        assert_eq!(caps, 2);
    }
}
pub(super) fn record_id(text: &str, fragment: &str, occurrence: usize) -> u32 {
    let record = text
        .lines()
        .filter(|line| line.starts_with('#') && line.contains(fragment))
        .nth(occurrence)
        .unwrap();
    record[1..record.find('=').unwrap()].parse().unwrap()
}
pub(super) fn change_record(text: &str, id: u32, change: impl FnOnce(&str) -> String) -> String {
    let prefix = format!("#{id}=");
    let line = text.lines().find(|line| line.starts_with(&prefix)).unwrap();
    let changed = change(line);
    assert_ne!(line, changed);
    text.replacen(line, &changed, 1)
}
fn first_reference(text: &str) -> u32 {
    let start = text.find('#').unwrap() + 1;
    let end = start + text[start..].bytes().take_while(u8::is_ascii_digit).count();
    text[start..end].parse().unwrap()
}
fn rejected(text: &str) -> Error {
    let document = CompoundDocument::parse(text, DocumentOptions::default(), None).unwrap();
    let solid = document.represented_members()[0].item_entity();
    match Parsed::parse(&document, solid, Default::default(), None) {
        Ok(_) => panic!("malformed actual geometry was accepted"),
        Err(error) => error,
    }
}
#[test]
fn actual_header_labels_are_irrelevant_and_equivalent_decimal_formats_remain_equal() {
    let text = fixture(false);
    let renamed = text
        .replace("MultiCAD", "Different CAD 8172")
        .replace("2D SPACE", "Other parameter space");
    let point = text
        .lines()
        .find(|line| line.contains("=CARTESIAN_POINT(") && line.matches(',').count() == 2)
        .unwrap();
    let id = point[1..point.find('=').unwrap()].parse().unwrap();
    let changed = change_record(&renamed, id, |line| {
        line.replacen("0.0000000000000000E0", "+0.000E00017", 1)
    });
    let document = CompoundDocument::parse(&changed, DocumentOptions::default(), None).unwrap();
    let parsed = Parsed::parse(
        &document,
        document.represented_members()[0].item_entity(),
        Default::default(),
        None,
    )
    .unwrap();
    assert_eq!(parsed.faces().len(), 10);
    let one = text
        .lines()
        .find(|line| {
            line.contains("=CARTESIAN_POINT(")
                && line.matches(',').count() == 2
                && line.contains("1.0000000000000000E0")
        })
        .unwrap();
    let id = one[1..one.find('=').unwrap()].parse().unwrap();
    let changed = change_record(&text, id, |line| {
        line.replacen("1.0000000000000000E0", "1000.000E-3", 1)
    });
    let document = CompoundDocument::parse(&changed, DocumentOptions::default(), None).unwrap();
    assert!(
        Parsed::parse(
            &document,
            document.represented_members()[0].item_entity(),
            Default::default(),
            None
        )
        .is_ok()
    );
}
#[test]
fn actual_literal_mutations_cannot_use_midpoint_equality_or_unsupported_geometry_as_binding() {
    let text = fixture(false);
    let normal = record_id(&text, "=DIRECTION('',(0.,0.,1.))", 0);
    let changed = change_record(&text, normal, |line| {
        line.replacen("(0.,0.,1.)", "(1E-400,0.,1.)", 1)
    });
    assert_eq!(rejected(&changed), Error::Unsupported);
    let tensor = record_id(&text, "=B_SPLINE_SURFACE_WITH_KNOTS(", 0);
    assert_eq!(
        rejected(&change_record(&text, tensor, |line| line.replacen(
            "('',11,11,",
            "('',10,11,",
            1
        ))),
        Error::Unsupported
    );
    let curve = record_id(&text, "=B_SPLINE_CURVE_WITH_KNOTS('',11,", 0);
    assert_eq!(
        rejected(&change_record(&text, curve, |line| line.replacen("('',11,", "('',10,", 1))),
        Error::Unsupported
    );
    assert_eq!(
        rejected(&change_record(&text, curve, |line| line.replacen("(12,", "(13,", 1))),
        Error::Unsupported
    );
    assert_eq!(
        rejected(&change_record(&text, curve, |line| line.replacen(
            "0.0000000000000000E0",
            "1E-400",
            1
        ))),
        Error::Unsupported
    );
    let edge = record_id(&text, "=EDGE_CURVE(", 0);
    assert_eq!(
        rejected(&change_record(&text, edge, |line| line.replace(",.T.)", ",.F.)"))),
        Error::Unsupported
    );
    let bound = record_id(&text, "=FACE_OUTER_BOUND(", 0);
    assert_eq!(
        rejected(&change_record(&text, bound, |line| line.replace(",.T.)", ",.F.)"))),
        Error::Unsupported
    );
    let oriented = record_id(&text, "=ORIENTED_EDGE(", 0);
    assert_eq!(
        rejected(&change_record(&text, oriented, |line| line.replace(",.T.)", ",.F.)"))),
        Error::Invalid
    );
    let point = record_id(&text, "=CARTESIAN_POINT(", 0);
    assert_eq!(
        rejected(&change_record(&text, point, |line| {
            let _ = line;
            format!("#{point}=CARTESIAN_POINT('',(1E8,0.,0.));")
        })),
        Error::Limit
    );
    let axis = record_id(&text, "=AXIS2_PLACEMENT_3D(", 0);
    let line = text
        .lines()
        .find(|line| line.starts_with(&format!("#{axis}=")))
        .unwrap();
    let origin = first_reference(line.split_once('=').unwrap().1);
    assert_eq!(
        rejected(&change_record(&text, origin, |line| {
            let _ = line;
            format!("#{origin}=CARTESIAN_POINT('',(1E-400,0.,0.));")
        })),
        Error::Unsupported
    );
}
#[test]
fn independent_graph_queries_preserve_source_and_real_resource_and_cancel_boundaries() {
    let text = fixture(false);
    let document = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
    let solid = document.represented_members()[0].item_entity();
    let options = Options::default();
    let result = Parsed::parse(&document, solid, options, None).unwrap();
    for bad in [
        Options {
            maximum_work: result.work_used() - 1,
            ..options
        },
        Options {
            maximum_bytes: result.query_budget_bytes() - 1,
            ..options
        },
    ] {
        assert!(matches!(
            Parsed::parse(&document, solid, bad, None),
            Err(Error::Limit)
        ));
    }
    assert_eq!(
        Parsed::parse(
            &document,
            solid,
            Options {
                maximum_work: result.work_used(),
                maximum_bytes: result.query_budget_bytes()
            },
            None
        )
        .unwrap()
        .work_used(),
        result.work_used()
    );
    for bad in [
        Options {
            maximum_work: 0,
            ..options
        },
        Options {
            maximum_work: options.maximum_work + 1,
            ..options
        },
        Options {
            maximum_bytes: 0,
            ..options
        },
        Options {
            maximum_bytes: options.maximum_bytes + 1,
            ..options
        },
    ] {
        assert!(matches!(
            Parsed::parse(&document, solid, bad, None),
            Err(Error::Invalid)
        ));
    }
    assert!(matches!(
        Parsed::parse(
            &document,
            solid,
            Options {
                maximum_bytes: 1,
                ..options
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        Parsed::parse(&document, solid, options, Some(&|| true)),
        Err(Error::Cancelled)
    ));
    let polls = std::cell::Cell::new(0usize);
    let count = || {
        polls.set(polls.get() + 1);
        false
    };
    Parsed::parse(&document, solid, options, Some(&count)).unwrap();
    let seen = std::cell::Cell::new(0usize);
    let cancel = || {
        seen.set(seen.get() + 1);
        seen.get() == polls.get()
    };
    assert!(matches!(
        Parsed::parse(&document, solid, options, Some(&cancel)),
        Err(Error::Cancelled)
    ));
    assert_eq!(document.literal_input(), text);
}
pub(super) fn fields(record: &str) -> Vec<&str> {
    let start = record.find('(').unwrap();
    let end = record.rfind(')').unwrap();
    let body = &record[start + 1..end];
    let mut result = Vec::new();
    let mut depth = 0usize;
    let mut quoted = false;
    let mut begin = 0;
    let bytes = body.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        match bytes[i] {
            b'\'' => {
                if quoted && bytes.get(i + 1) == Some(&b'\'') {
                    i += 1;
                } else {
                    quoted = !quoted;
                }
            }
            b'(' if !quoted => depth += 1,
            b')' if !quoted => depth -= 1,
            b',' if !quoted && depth == 0 => {
                result.push(&body[begin..i]);
                begin = i + 1;
            }
            _ => {}
        }
        i += 1;
    }
    result.push(&body[begin..]);
    result
}
pub(super) fn with_fields(record: &str, fields: &[&str]) -> String {
    let start = record.find('(').unwrap();
    let end = record.rfind(')').unwrap();
    format!(
        "{}{}{}",
        &record[..start + 1],
        fields.join(","),
        &record[end..]
    )
}
fn refuses_pipeline(text: &str, solid: u32) -> bool {
    let Ok(document) = CompoundDocument::parse(text, DocumentOptions::default(), None) else {
        return true;
    };
    Parsed::parse(&document, solid, Default::default(), None).is_err()
}
#[test]
fn actual_step_schema_mutation_corpus_refuses_missing_wrong_and_ragged_geometry_fields() {
    let text = fixture(false);
    let baseline = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
    let solid = baseline.represented_members()[0].item_entity();
    let mut representatives = Vec::new();
    for fragment in [
        "=CARTESIAN_POINT(",
        "=VERTEX_POINT(",
        "=B_SPLINE_CURVE_WITH_KNOTS('',11,",
        "=B_SPLINE_CURVE_WITH_KNOTS('',1,",
        "=B_SPLINE_SURFACE_WITH_KNOTS(",
        "=AXIS2_PLACEMENT_3D(",
        "=DIRECTION(",
        "=PLANE(",
        "=DEFINITIONAL_REPRESENTATION(",
        "=PCURVE(",
        "=SURFACE_CURVE(",
        "=EDGE_CURVE(",
        "=ORIENTED_EDGE(",
        "=EDGE_LOOP(",
        "=FACE_OUTER_BOUND(",
        "=ADVANCED_FACE(",
        "=CLOSED_SHELL(",
        "=MANIFOLD_SOLID_BREP(",
    ] {
        representatives.push(record_id(&text, fragment, 0));
    }
    let point2 = text
        .lines()
        .find(|line| line.contains("=CARTESIAN_POINT(") && line.matches(',').count() == 2)
        .unwrap();
    representatives.push(point2[1..point2.find('=').unwrap()].parse().unwrap());
    for id in representatives {
        let prefix = format!("#{id}=");
        let record = text.lines().find(|line| line.starts_with(&prefix)).unwrap();
        let arguments = fields(record);
        for index in 0..arguments.len() {
            let mut changed = arguments.clone();
            changed[index] = "$";
            let output = change_record(&text, id, |line| with_fields(line, &changed));
            assert!(
                refuses_pipeline(&output, solid),
                "wrong field {index} of{record}"
            );
            let mut changed = arguments.clone();
            changed.remove(index);
            let output = change_record(&text, id, |line| with_fields(line, &changed));
            assert!(
                refuses_pipeline(&output, solid),
                "missing field {index} of{record}"
            );
            if arguments[index].starts_with('(') && arguments[index].ends_with(')') {
                let list = fields(arguments[index]);
                if !list.is_empty() {
                    let mut shortened = list.clone();
                    shortened.pop();
                    let value = format!("({})", shortened.join(","));
                    let mut changed = arguments.clone();
                    changed[index] = &value;
                    let output = change_record(&text, id, |line| with_fields(line, &changed));
                    assert!(
                        refuses_pipeline(&output, solid),
                        "short list {index} of{record}"
                    );
                }
            }
        }
    }
    assert!(matches!(
        Parsed::parse(&baseline, solid + 1, Default::default(), None),
        Err(Error::Unsupported)
    ));
}
