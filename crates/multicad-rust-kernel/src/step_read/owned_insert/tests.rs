use super::*;
fn make_file(left: bool) -> String {
    use crate::{
        authored_flat_thread::{AuthoredThread, SourceOptions},
        helical_thread::Spec,
        step::insert::{
            CurveOptions, DomainOptions, DomainPlan, SurfacePlan, TopologyPlan,
            write_local_candidate,
        },
        wood_insert::{BoundaryOptions, Insert},
    };
    let spec = Spec {
        major_radius_mm: 5.,
        pitch_mm: 2.,
        length_mm: 4.,
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
    let body = Insert::new(outer, inner, 0.5, None)
        .unwrap()
        .boundary_with_interrupt(BoundaryOptions::default(), None)
        .unwrap();
    let d = DomainPlan::prepare(&body, DomainOptions::default(), None).unwrap();
    let s = SurfacePlan::prepare(&d, CurveOptions::default(), None).unwrap();
    let t = TopologyPlan::prepare(&s, 16_000_000, 64 * 1024 * 1024, None).unwrap();
    write_local_candidate(&t, 16_000_000, 64 * 1024 * 1024, None).unwrap()
}
fn file(left: bool) -> String {
    static FILES: std::sync::OnceLock<[String; 2]> = std::sync::OnceLock::new();
    FILES.get_or_init(|| [make_file(false), make_file(true)])[usize::from(left)].clone()
}
fn record_change(text: &str, id: u32, change: impl FnOnce(&str) -> String) -> String {
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
