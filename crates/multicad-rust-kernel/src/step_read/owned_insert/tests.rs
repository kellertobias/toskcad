use super::*;
fn make_file(left: bool) -> String {
    make_file_case(left, 4.)
}
pub(crate) fn graph(left: bool, length: f64) -> crate::wood_insert::Boundary {
    use crate::{
        authored_flat_thread::{AuthoredThread, SourceOptions},
        helical_thread::Spec,
        wood_insert::{BoundaryOptions, Insert},
    };
    let spec = Spec {
        major_radius_mm: 5.,
        pitch_mm: 2.,
        length_mm: length,
        depth_mm: 1.,
        crest_width_mm: 0.4,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts: 1,
        left_hand: left,
        root_radius_mm: None,
    };
    let options = SourceOptions {
        position_limit_mm: 1e-7,
        normal_angle_limit_rad: 1e-7,
    };
    let outer = AuthoredThread::new(spec, options).unwrap();
    let inner = AuthoredThread::new(
        Spec {
            major_radius_mm: 1.5,
            pitch_mm: 0.5,
            depth_mm: (3. - (2.459 + 2.599) / 2.) / 2.,
            crest_width_mm: 0.25 - (3. - (2.675 + 2.775) / 2.) / 3_f64.sqrt(),
            ..spec
        },
        options,
    )
    .unwrap();
    Insert::new(outer, inner, 0.5, None)
        .unwrap()
        .boundary_with_interrupt(BoundaryOptions::default(), None)
        .unwrap()
}
fn make_file_case(left: bool, length: f64) -> String {
    use crate::step::insert::{
        CurveOptions, DomainOptions, DomainPlan, SurfacePlan, TopologyPlan, write_local_candidate,
    };
    let body = graph(left, length);
    let d = DomainPlan::prepare(&body, DomainOptions::default(), None).unwrap();
    let s = SurfacePlan::prepare(&d, CurveOptions::default(), None).unwrap();
    let t = TopologyPlan::prepare(&s, 16_000_000, 64 * 1024 * 1024, None).unwrap();
    write_local_candidate(&t, 16_000_000, 64 * 1024 * 1024, None).unwrap()
}
pub(super) fn file(left: bool) -> String {
    static FILES: std::sync::OnceLock<[String; 2]> = std::sync::OnceLock::new();
    FILES.get_or_init(|| [make_file(false), make_file(true)])[usize::from(left)].clone()
}
pub(super) fn record_change(text: &str, id: u32, change: impl FnOnce(&str) -> String) -> String {
    let prefix = format!("#{id}=");
    let line = text.lines().find(|line| line.starts_with(&prefix)).unwrap();
    let next = change(line);
    assert_ne!(line, next);
    text.replacen(line, &next, 1)
}
#[test]
fn actual_file_retains_original_scalar_carriers_and_complete_represented_graph() {
    for left in [false, true] {
        let text = file(left);
        let p = parse(&text, Options::default(), None).unwrap();
        assert_eq!(p.literal_input(), text);
        assert_eq!(p.planes().len(), 2);
        assert_eq!(p.faces().len(), p.tensors().len() + 2);
        assert!(p.work_used() <= 16_000_000);
        assert!(p.retained_budget_bytes() <= 64 * 1024 * 1024);
        assert!(p.literal_entity_count() > p.points().len());
        for point in p.points() {
            assert_eq!(p.point(point.entity()).unwrap().entity(), point.entity());
            for v in point.coordinates() {
                assert_eq!(v.spelling().parse::<f64>().unwrap(), v.value());
                let [a, b] = v.bounds();
                assert!(a <= v.value() && b >= v.value());
            }
        }
        eprintln!(
            "ACTUAL_PARSE left{left} faces{} edges{} points{} work{} bytes{}",
            p.faces().len(),
            p.edges().len(),
            p.points().len(),
            p.work_used(),
            p.retained_budget_bytes()
        );
    }
}

#[test]
fn actual_local_file_has_literal_dual_support_and_full_uv_domain_binding() {
    for left in [false, true] {
        let text = file(left);
        let parsed = parse(&text, Options::default(), None).unwrap();
        let proof = parsed.support_bindings(Options::default(), None).unwrap();
        assert!(std::ptr::eq(proof.source(), &parsed));
        assert_eq!(proof.tensors().len(), parsed.tensors().len());
        assert_eq!(proof.edge_face_indices().len(), parsed.edges().len());
        eprintln!(
            "ACTUAL_BIND left{left} work{} bytes{}",
            proof.work_used(),
            proof.retained_budget_bytes()
        );
    }
}

#[test]
fn actual_raw_coefficients_earn_closed_genus_one_material_without_source_headers() {
    for left in [false, true] {
        let text = file(left);
        let parsed = parse(&text, Options::default(), None).unwrap();
        let bindings = parsed.support_bindings(Options::default(), None).unwrap();
        let proof = bindings.material_bounds(Options::default(), None).unwrap();
        assert!(std::ptr::eq(proof.binding(), &bindings));
        assert!(proof.minimum_signed_jacobian() > 0.);
        assert!(proof.minimum_signed_cap_derivative() > 0.);
        assert!(proof.component_radius_mm()[0][0] > proof.component_radius_mm()[1][1]);
        eprintln!(
            "ACTUAL_MATERIAL left{left} work{} bytes{} radii{:?}",
            proof.work_used(),
            proof.retained_budget_bytes(),
            proof.component_radius_mm()
        );
    }
}

#[test]
fn original_schema_and_sub_ulp_surface_numeric_leaves_cannot_be_replaced_by_float_bits() {
    let text = file(false);
    let changed = text.replacen(
        "B_SPLINE_CURVE_WITH_KNOTS('',11,",
        "B_SPLINE_CURVE_WITH_KNOTS('',11.0000000000000001,",
        1,
    );
    assert!(matches!(
        parse(&changed, Options::default(), None),
        Err(Error::Invalid)
    ));
    let p = parse(&text, Options::default(), None).unwrap();
    let proof = p.support_bindings(Options::default(), None).unwrap();
    let cap = &p.faces()[proof.cap_face_indices()[0]];
    let edge = p.edge(cap.loops()[0].coedges()[0].edge).unwrap();
    let pc = edge
        .pcurves()
        .iter()
        .find(|pc| pc.support_entity == cap.surface_entity())
        .unwrap();
    let curve = p.curve(pc.curve_entity).unwrap();
    let (point, index) = curve
        .pole_entities()
        .iter()
        .find_map(|&id| {
            let point = p.point(id).unwrap();
            point
                .coordinates()
                .iter()
                .position(|v| {
                    v.value().abs() >= 2.
                        && v.value().fract() == 0.
                        && v.spelling().contains(".0000000000000000E")
                })
                .map(|i| (point, i))
        })
        .unwrap();
    let old = point.coordinates()[index];
    let replacement = old
        .spelling()
        .replacen(".0000000000000000E", ".0000000000000001E", 1);
    assert_eq!(replacement.parse::<f64>().unwrap(), old.value());
    let changed = record_change(&text, point.entity(), |line| {
        line.replacen(old.spelling(), &replacement, 1)
    });
    let changed = parse(&changed, Options::default(), None).unwrap();
    assert!(matches!(
        changed.support_bindings(Options::default(), None),
        Err(Error::Unsupported)
    ));
}
#[test]
fn raw_tensor_controls_uv_underflow_and_tilted_plane_axes_are_not_healed() {
    let text = file(false);
    let p = parse(&text, Options::default(), None).unwrap();
    let tensor = &p.tensors()[0];
    let rows = tensor.rows();
    let a = rows[0][0];
    let b = rows[0][1];
    let changed = record_change(&text, tensor.entity(), |line| {
        line.replacen(&format!("#{a},#{b},"), &format!("#{a},#{a},"), 1)
    });
    let p1 = parse(&changed, Options::default(), None).unwrap();
    assert!(matches!(
        p1.support_bindings(Options::default(), None),
        Err(Error::Unsupported)
    ));
    let uv = p
        .points()
        .iter()
        .find(|point| point.coordinates().len() == 2 && point.coordinates()[0].is_zero())
        .unwrap();
    let changed = record_change(&text, uv.entity(), |line| {
        line.replacen(uv.coordinates()[0].spelling(), "1.E-1000", 1)
    });
    let p1 = parse(&changed, Options::default(), None).unwrap();
    assert_eq!(p1.point(uv.entity()).unwrap().coordinates()[0].value(), 0.);
    assert!(matches!(
        p1.support_bindings(Options::default(), None),
        Err(Error::Unsupported)
    ));
    let normal = p.planes()[0].direction_entities()[0];
    let changed = record_change(&text, normal, |line| {
        line.replacen("0.0000000000000000E0", "1.E-1000", 1)
    });
    let p1 = parse(&changed, Options::default(), None).unwrap();
    assert!(matches!(
        p1.support_bindings(Options::default(), None),
        Err(Error::Unsupported)
    ));
}
fn interior_control(p: &Parsed<'_>, b: &SupportBindings<'_, '_>) -> u32 {
    p.edges()
        .iter()
        .enumerate()
        .find_map(|(i, e)| {
            if b.edge_face_indices()[i]
                .iter()
                .all(|&f| matches!(p.faces()[f].surface(), SurfaceIndex::Tensor(_)))
            {
                let c = p.curve(e.curve_entity()).unwrap();
                if c.degree() == 11 {
                    let id = c.pole_entities()[5];
                    let point = p.point(id).unwrap();
                    if point.coordinates()[0].value().abs() > point.coordinates()[1].value().abs() {
                        return Some(id);
                    }
                }
            }
            None
        })
        .unwrap()
}
#[test]
fn actual_changed_shared_carriers_are_evaluated_and_folded_geometry_is_refused() {
    let text = file(false);
    let p = parse(&text, Options::default(), None).unwrap();
    let b = p.support_bindings(Options::default(), None).unwrap();
    let id = interior_control(&p, &b);
    let x = p.point(id).unwrap().coordinates()[0];
    for folded in [false, true] {
        let value = if folded {
            -1000. * x.value().signum()
        } else {
            x.value() + 0.0001
        };
        let changed = record_change(&text, id, |line| {
            line.replacen(x.spelling(), &format!("{value:.16E}"), 1)
        });
        let parsed = parse(&changed, Options::default(), None).unwrap();
        let binding = parsed.support_bindings(Options::default(), None).unwrap();
        let actual = binding.material_bounds(Options::default(), None);
        if folded {
            assert!(matches!(
                actual,
                Err(Error::Unsupported) | Err(Error::Limit)
            ));
        } else {
            let actual = actual.unwrap();
            assert!(actual.minimum_signed_jacobian() > 0.);
            assert_eq!(parsed.point(id).unwrap().coordinates()[0].value(), value);
        }
    }
}
#[test]
fn actual_face_orientation_and_cap_hole_ownership_are_geometry_not_labels() {
    let text = file(false);
    let p = parse(&text, Options::default(), None).unwrap();
    let b = p.support_bindings(Options::default(), None).unwrap();
    let face = &p.faces()[0];
    let changed = record_change(&text, face.entity(), |line| {
        if line.ends_with(",.F.);") {
            line.strip_suffix(",.F.);").unwrap().to_owned() + ",.T.);"
        } else {
            line.strip_suffix(",.T.);").unwrap().to_owned() + ",.F.);"
        }
    });
    let p1 = parse(&changed, Options::default(), None).unwrap();
    assert!(matches!(
        p1.support_bindings(Options::default(), None),
        Err(Error::Unsupported)
    ));
    let cap = &p.faces()[b.cap_face_indices()[0]];
    let [outer, inner] = cap.loops() else {
        panic!("actual annular cap")
    };
    let a = outer.entity();
    let c = inner.entity();
    let changed = record_change(&text, outer.bound_entity(), |line| {
        line.replacen(&format!("#{a},"), &format!("#{c},"), 1)
    });
    let changed = record_change(&changed, inner.bound_entity(), |line| {
        line.replacen(&format!("#{c},"), &format!("#{a},"), 1)
    });
    let p1 = parse(&changed, Options::default(), None).unwrap();
    let b1 = p1.support_bindings(Options::default(), None).unwrap();
    assert!(matches!(
        b1.material_bounds(Options::default(), None),
        Err(Error::Unsupported)
    ));
}

#[test]
fn actual_metric_nonbinary_edit_lengths_keep_physical_readback_within_fixed_geometry_budget() {
    for (left, length) in [(false, 4.3), (true, 4.7)] {
        let text = make_file_case(left, length);
        let p = parse(&text, Options::default(), None).unwrap();
        let b = p.support_bindings(Options::default(), None).unwrap();
        eprintln!(
            "NONBINARY_SOURCE {left}/{length} faces{} bindwork{} bytes{}",
            p.faces().len(),
            b.work_used(),
            b.retained_budget_bytes()
        );
        let proof = b.material_bounds(Options::default(), None).unwrap();
        eprintln!(
            "NONBINARY_MATERIAL {left}/{length} work{}",
            proof.work_used()
        );
    }
}
#[test]
fn original_graph_accessors_preserve_all_literal_and_semantic_identities() {
    let text = file(false);
    let p = parse(&text, Options::default(), None).unwrap();
    for curve in p.curves() {
        assert_eq!(p.curve(curve.entity()).unwrap().entity(), curve.entity());
        assert_eq!(curve.knots().len(), curve.multiplicities().len());
    }
    for plane in p.planes() {
        assert!(
            p.literal_input()
                .contains(&format!("#{}=AXIS2_PLACEMENT_3D", plane.axis_entity()))
        );
        assert!(p.point(plane.origin_entity()).is_some());
        assert!(
            p.literal_input()
                .contains(&format!("#{}=PLANE", plane.entity()))
        );
    }
    for edge in p.edges() {
        assert_eq!(p.edge(edge.entity()).unwrap().entity(), edge.entity());
        assert!(
            p.literal_input()
                .contains(&format!("#{}=SURFACE_CURVE", edge.surface_curve_entity()))
        );
        for v in edge.vertex_entities() {
            assert_eq!(p.vertex(v).unwrap().entity, v);
        }
    }
    for tensor in p.tensors() {
        assert!(
            p.literal_input()
                .contains(&format!("#{}=B_SPLINE_SURFACE_WITH_KNOTS", tensor.entity()))
        );
        assert_eq!(tensor.knots_v().len(), tensor.multiplicities_v().len());
    }
    assert!(p.point(u32::MAX).is_none());
    assert!(p.curve(u32::MAX).is_none());
    assert!(p.edge(u32::MAX).is_none());
    assert!(p.vertex(u32::MAX).is_none());
    let b = p.support_bindings(Options::default(), None).unwrap();
    for t in b.tensors() {
        for v in t.corner_vertex_entities() {
            assert!(p.vertex(v).is_some());
        }
    }
    let m = b.material_bounds(Options::default(), None).unwrap();
    assert!(m.outer_component() < 2);
    for (i, face) in p.faces().iter().enumerate() {
        assert_eq!(
            m.tensor_component(i).is_some(),
            matches!(face.surface(), SurfaceIndex::Tensor(_))
        );
    }
    assert!(m.tensor_component(usize::MAX).is_none());
    let bounds = m.bounds_mm();
    assert_eq!(bounds[2], [0., 4.]);
    for tensor in p.tensors() {
        for row in tensor.rows() {
            for &id in row {
                for (i, v) in p.point(id).unwrap().coordinates().iter().enumerate() {
                    assert!(bounds[i][0] <= v.value() && v.value() <= bounds[i][1]);
                }
            }
        }
    }
}
#[test]
fn actual_raw_reader_resource_and_cancellation_failures_publish_no_proof() {
    let text = file(false);
    for options in [
        Options {
            maximum_work: 0,
            ..Options::default()
        },
        Options {
            maximum_work: 16_000_001,
            ..Options::default()
        },
        Options {
            retained_budget_bytes: 0,
            ..Options::default()
        },
        Options {
            retained_budget_bytes: 64 * 1024 * 1024 + 1,
            ..Options::default()
        },
    ] {
        assert!(matches!(parse(&text, options, None), Err(Error::Invalid)));
    }
    for options in [
        Options {
            maximum_work: 1,
            ..Options::default()
        },
        Options {
            retained_budget_bytes: 1,
            ..Options::default()
        },
    ] {
        assert!(matches!(parse(&text, options, None), Err(Error::Limit)));
    }
    assert!(matches!(
        parse(&text, Options::default(), Some(&|| true)),
        Err(Error::Cancelled)
    ));
    let polls = std::cell::Cell::new(0);
    let stop = || {
        polls.set(polls.get() + 1);
        polls.get() > 100
    };
    assert!(matches!(
        parse(&text, Options::default(), Some(&stop)),
        Err(Error::Cancelled)
    ));
    let p = parse(&text, Options::default(), None).unwrap();
    for options in [
        Options {
            maximum_work: 0,
            ..Options::default()
        },
        Options {
            retained_budget_bytes: 0,
            ..Options::default()
        },
    ] {
        assert!(matches!(
            p.support_bindings(options, None),
            Err(Error::Invalid)
        ));
    }
    assert!(matches!(
        p.support_bindings(
            Options {
                maximum_work: p.work_used(),
                ..Options::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        p.support_bindings(
            Options {
                retained_budget_bytes: p.retained_budget_bytes(),
                ..Options::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        p.support_bindings(Options::default(), Some(&|| true)),
        Err(Error::Cancelled)
    ));
    let b = p.support_bindings(Options::default(), None).unwrap();
    assert!(matches!(
        b.material_bounds(
            Options {
                maximum_work: 0,
                ..Options::default()
            },
            None
        ),
        Err(Error::Invalid)
    ));
    assert!(matches!(
        b.material_bounds(
            Options {
                retained_budget_bytes: 0,
                ..Options::default()
            },
            None
        ),
        Err(Error::Invalid)
    ));
    assert!(matches!(
        b.material_bounds(
            Options {
                maximum_work: b.work_used(),
                ..Options::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        b.material_bounds(
            Options {
                retained_budget_bytes: b.retained_budget_bytes(),
                ..Options::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        b.material_bounds(Options::default(), Some(&|| true)),
        Err(Error::Cancelled)
    ));
    for (before, after) in [
        (super::super::ReadError::Invalid, Error::Invalid),
        (super::super::ReadError::Unsupported, Error::Unsupported),
        (super::super::ReadError::Limit, Error::Limit),
    ] {
        assert_eq!(Error::from(before), after);
    }
    for (before, after) in [
        (crate::step::insert::Error::Invalid, Error::Invalid),
        (crate::step::insert::Error::Unresolved, Error::Unsupported),
        (crate::step::insert::Error::Limit, Error::Limit),
        (crate::step::insert::Error::Cancelled, Error::Cancelled),
    ] {
        assert_eq!(Error::from(before), after);
    }
}
pub(super) fn argument_change(text: &str, id: u32, index: usize, value: &str) -> String {
    record_change(text, id, |line| {
        let open = line.find('(').unwrap();
        let inside = &line[open + 1..line.len() - 2];
        let mut fields = Vec::new();
        let mut depth = 0usize;
        let mut quoted = false;
        let mut start = 0;
        for (i, b) in inside.bytes().enumerate() {
            if b == b'\'' {
                quoted = !quoted;
            }
            if !quoted {
                if b == b'(' {
                    depth += 1;
                }
                if b == b')' {
                    depth -= 1;
                }
                if b == b',' && depth == 0 {
                    fields.push(&inside[start..i]);
                    start = i + 1;
                }
            }
        }
        fields.push(&inside[start..]);
        fields[index] = value;
        format!("{}({});", &line[..open], fields.join(","))
    })
}
#[test]
fn actual_malformed_schema_records_and_missing_graph_references_are_atomic_refusals() {
    let text = file(false);
    let p = parse(&text, Options::default(), None).unwrap();
    let c = &p.curves()[0];
    let t = &p.tensors()[0];
    let plane = &p.planes()[0];
    let face = &p.faces()[0];
    let edge = &p.edges()[0];
    let pc = &edge.pcurves()[0];
    let vertex = &p.vertices()[0];
    let point = p
        .points()
        .iter()
        .find(|p| p.coordinates().len() == 3)
        .unwrap();
    let bound = &face.loops()[0];
    let cases = [
        (point.entity(), 0, ".T."),
        (point.entity(), 1, "(0.,0.,0.,0.)"),
        (c.entity(), 1, "3"),
        (c.entity(), 2, "()"),
        (c.entity(), 3, ".POLYLINE_FORM."),
        (c.entity(), 6, "(1,1)"),
        (c.entity(), 7, "(0.,0.)"),
        (c.entity(), 7, "(0.,0.5,1.)"),
        (t.entity(), 1, "2"),
        (t.entity(), 3, "()"),
        (t.entity(), 8, "(1,1)"),
        (t.entity(), 10, "(0.,0.5,1.)"),
        (t.entity(), 11, "(0.,0.)"),
        (plane.entity(), 1, "#4294967295"),
        (plane.axis_entity(), 2, "$"),
        (face.entity(), 1, "()"),
        (face.entity(), 3, ".U."),
        (bound.bound_entity(), 2, ".F."),
        (edge.entity(), 4, ".F."),
        (edge.surface_curve_entity(), 3, ".PCURVE_S1."),
        (pc.entity, 2, "#4294967295"),
        (pc.representation, 1, "()"),
        (vertex.entity, 1, "#4294967295"),
    ];
    for (id, index, value) in cases {
        let changed = argument_change(&text, id, index, value);
        assert!(
            parse(&changed, Options::default(), None).is_err(),
            "actualrecord#{id} arg{index}={value}"
        );
    }
    let changed = argument_change(
        &text,
        edge.entity(),
        1,
        &format!("#{}", edge.vertex_entities()[1]),
    );
    assert!(matches!(
        parse(&changed, Options::default(), None),
        Err(Error::Unsupported)
    ));
}
#[test]
fn coherently_translated_actual_z_planes_and_carriers_keep_their_own_datum() {
    let text = file(false);
    let p = parse(&text, Options::default(), None).unwrap();
    let mut changed = text.clone();
    for point in p.points() {
        if point.coordinates().len() != 3 {
            continue;
        }
        let c = point.coordinates();
        let z = c[2].value() - 10.;
        let value = format!("({},{},{z:.16E})", c[0].spelling(), c[1].spelling());
        changed = argument_change(&changed, point.entity(), 1, &value);
    }
    let p = parse(&changed, Options::default(), None).unwrap();
    let b = p.support_bindings(Options::default(), None).unwrap();
    assert_eq!(b.cap_heights().map(|h| h.value()), [-10., -6.]);
    let m = b.material_bounds(Options::default(), None).unwrap();
    assert_eq!(m.bounds_mm()[2], [-10., -6.]);
    assert!(m.minimum_signed_jacobian() > 0.);
}
#[test]
fn actual_face_and_plane_entity_orders_do_not_define_material_or_datum() {
    let text = file(false);
    let p = parse(&text, Options::default(), None).unwrap();
    let shell = p
        .entities
        .iter()
        .find_map(|(&id, e)| matches!(e.as_slice(),[e] if e.name=="CLOSED_SHELL").then_some(id))
        .unwrap();
    let reversed = p
        .faces()
        .iter()
        .rev()
        .map(|f| format!("#{}", f.entity()))
        .collect::<Vec<_>>()
        .join(",");
    let changed = argument_change(&text, shell, 1, &format!("({reversed})"));
    let p1 = parse(&changed, Options::default(), None).unwrap();
    let b1 = p1.support_bindings(Options::default(), None).unwrap();
    let m1 = b1.material_bounds(Options::default(), None).unwrap();
    assert_eq!(m1.outer_component(), 1);
    assert_eq!(m1.bounds_mm()[2], [0., 4.]);
    let [a, c] = [p.planes()[0].entity(), p.planes()[1].entity()];
    let mut changed = String::with_capacity(text.len());
    let bytes = text.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] == b'#' {
            let start = i + 1;
            let mut end = start;
            while end < bytes.len() && bytes[end].is_ascii_digit() {
                end += 1;
            }
            let id = text[start..end].parse::<u32>().unwrap();
            let id = if id == a {
                c
            } else if id == c {
                a
            } else {
                id
            };
            changed.push_str(&format!("#{id}"));
            i = end;
        } else {
            changed.push(bytes[i] as char);
            i += 1;
        }
    }
    let p1 = parse(&changed, Options::default(), None).unwrap();
    let b1 = p1.support_bindings(Options::default(), None).unwrap();
    assert_eq!(b1.cap_heights().map(|h| h.value()), [0., 4.]);
    assert!(b1.material_bounds(Options::default(), None).is_ok());
    let limited = Options {
        retained_budget_bytes: p.retained_budget_bytes() - 1,
        ..Options::default()
    };
    assert!(matches!(parse(&text, limited, None), Err(Error::Limit)));
}
