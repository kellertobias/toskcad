use super::*;
use crate::step::insert::{
    CurveOptions, DomainOptions, DomainPlan, SurfacePlan, TopologyPlan, write_local_candidate,
};
#[test]
fn actual_source_rows_incidence_and_decoded_mass_keep_separate_authority() {
    for (left, length) in [(false, 4.), (true, 4.), (false, 4.3)] {
        let body = super::super::tests::graph(left, length);
        let d = DomainPlan::prepare(&body, DomainOptions::default(), None).unwrap();
        let s = SurfacePlan::prepare(&d, CurveOptions::default(), None).unwrap();
        let t = TopologyPlan::prepare(&s, 16_000_000, 64 * 1024 * 1024, None).unwrap();
        let text = write_local_candidate(&t, 16_000_000, 64 * 1024 * 1024, None).unwrap();
        let p = parse(&text, Options::default(), None).unwrap();
        let b = p.support_bindings(Options::default(), None).unwrap();
        let c = b
            .source_correspondence(&t, CorrespondenceOptions::default(), None)
            .unwrap();
        assert!(std::ptr::eq(c.actual(), &b));
        assert!(std::ptr::eq(c.expected(), &t));
        assert!(c.position_error_mm() > 0. && c.position_error_mm() < 1e-7);
        assert!(c.row_derivative_error_mm() > 0. && c.row_derivative_error_mm() < 1e-6);
        assert!(c.ruling_derivative_error_mm() >= 2. * c.position_error_mm());
        assert!(c.work_used() <= 64_000_000 && c.retained_budget_bytes() <= 64 * 1024 * 1024);
        let material = b.material_bounds(Options::default(), None).unwrap();
        let volume = material
            .occupied_volume(FluxOptions::default(), None)
            .unwrap();
        assert!(volume.volume_mm3()[0] > 0.);
        eprintln!(
            "SOURCE_CORRESP {left}/{length} pos{} rowd{} rulingd{} work{} bytes{}",
            c.position_error_mm(),
            c.row_derivative_error_mm(),
            c.ruling_derivative_error_mm(),
            c.work_used(),
            c.retained_budget_bytes()
        );
    }
}
#[test]
fn actual_geometry_correspondence_limits_do_not_forge_a_mass_or_normal_receipt() {
    let body = super::super::tests::graph(false, 4.);
    let d = DomainPlan::prepare(&body, DomainOptions::default(), None).unwrap();
    let s = SurfacePlan::prepare(&d, CurveOptions::default(), None).unwrap();
    let t = TopologyPlan::prepare(&s, 16_000_000, 64 * 1024 * 1024, None).unwrap();
    let text = write_local_candidate(&t, 16_000_000, 64 * 1024 * 1024, None).unwrap();
    let p = parse(&text, Options::default(), None).unwrap();
    let b = p.support_bindings(Options::default(), None).unwrap();
    let c = b
        .source_correspondence(&t, CorrespondenceOptions::default(), None)
        .unwrap();
    for options in [
        CorrespondenceOptions {
            position_limit_mm: f64::NAN,
            ..CorrespondenceOptions::default()
        },
        CorrespondenceOptions {
            maximum_work: 0,
            ..CorrespondenceOptions::default()
        },
        CorrespondenceOptions {
            retained_budget_bytes: 0,
            ..CorrespondenceOptions::default()
        },
    ] {
        assert!(matches!(
            b.source_correspondence(&t, options, None),
            Err(Error::Invalid)
        ));
    }
    assert!(matches!(
        b.source_correspondence(
            &t,
            CorrespondenceOptions {
                maximum_work: c.work_used() - 1,
                ..CorrespondenceOptions::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        b.source_correspondence(
            &t,
            CorrespondenceOptions {
                retained_budget_bytes: 1,
                ..CorrespondenceOptions::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        b.source_correspondence(
            &t,
            CorrespondenceOptions {
                position_limit_mm: c.position_error_mm() / 2.,
                ..CorrespondenceOptions::default()
            },
            None
        ),
        Err(Error::Unsupported)
    ));
    assert!(matches!(
        b.source_correspondence(&t, CorrespondenceOptions::default(), Some(&|| true)),
        Err(Error::Cancelled)
    ));
}
#[test]
fn actual_changed_coefficients_are_measured_and_geometry_equivalent_reorder_is_not_guessed() {
    let body = super::super::tests::graph(false, 4.);
    let d = DomainPlan::prepare(&body, DomainOptions::default(), None).unwrap();
    let s = SurfacePlan::prepare(&d, CurveOptions::default(), None).unwrap();
    let t = TopologyPlan::prepare(&s, 16_000_000, 64 * 1024 * 1024, None).unwrap();
    let text = write_local_candidate(&t, 16_000_000, 64 * 1024 * 1024, None).unwrap();
    let p = parse(&text, Options::default(), None).unwrap();
    let pole = p
        .tensors()
        .iter()
        .flat_map(|t| t.rows())
        .flat_map(|r| r[1..r.len() - 1].iter())
        .map(|id| p.point(*id).unwrap())
        .find(|p| {
            let z = p.coordinates()[2].value();
            z > 0. && z < 4.
        })
        .unwrap();
    for shift in [1e-9, 1e-4] {
        let v = pole.coordinates();
        let xyz = format!(
            "({:.16E},{},{})",
            v[0].value() + shift,
            v[1].spelling(),
            v[2].spelling()
        );
        let changed = super::super::tests::argument_change(&text, pole.entity(), 1, &xyz);
        let pp = parse(&changed, Options::default(), None).unwrap();
        let b = pp.support_bindings(Options::default(), None).unwrap();
        assert!(b.material_bounds(Options::default(), None).is_ok());
        let receipt = b.source_correspondence(&t, CorrespondenceOptions::default(), None);
        if shift < 1e-7 {
            let c = receipt.unwrap();
            assert!(c.position_error_mm() >= shift * 0.99);
            assert!(matches!(
                b.source_correspondence(
                    &t,
                    CorrespondenceOptions {
                        position_limit_mm: 1e-10,
                        ..CorrespondenceOptions::default()
                    },
                    None
                ),
                Err(Error::Unsupported)
            ));
        } else {
            assert!(matches!(receipt, Err(Error::Unsupported)));
        }
    }
    // Cyclic loop starts are physically equivalent; canonical source-known
    // correspondence deliberately refuses rather than infer a new owner map.
    let wire = &p.faces()[0].loops()[0];
    let mut ids = refs(&p.entities[&wire.entity()][0].args[1]).unwrap();
    ids.rotate_left(1);
    let changed = super::super::tests::argument_change(
        &text,
        wire.entity(),
        1,
        &format!(
            "({})",
            ids.iter()
                .map(|id| format!("#{id}"))
                .collect::<Vec<_>>()
                .join(",")
        ),
    );
    let pp = parse(&changed, Options::default(), None).unwrap();
    let b = pp.support_bindings(Options::default(), None).unwrap();
    assert!(b.material_bounds(Options::default(), None).is_ok());
    assert!(matches!(
        b.source_correspondence(&t, CorrespondenceOptions::default(), None),
        Err(Error::Unsupported)
    ));
    let shell = p
        .entities
        .iter()
        .find_map(|(&id, r)| matches!(r.as_slice(),[e]if e.name=="CLOSED_SHELL").then_some(id))
        .unwrap();
    let reversed = p
        .faces()
        .iter()
        .rev()
        .map(|f| format!("#{}", f.entity()))
        .collect::<Vec<_>>()
        .join(",");
    let changed = super::super::tests::argument_change(&text, shell, 1, &format!("({reversed})"));
    let pp = parse(&changed, Options::default(), None).unwrap();
    let b = pp.support_bindings(Options::default(), None).unwrap();
    assert!(b.material_bounds(Options::default(), None).is_ok());
    assert!(matches!(
        b.source_correspondence(&t, CorrespondenceOptions::default(), None),
        Err(Error::Unsupported)
    ));
}
#[test]
fn actual_count_and_numeric_failures_keep_the_source_query_bounded() {
    let body = super::super::tests::graph(false, 4.);
    let d = DomainPlan::prepare(&body, DomainOptions::default(), None).unwrap();
    let s = SurfacePlan::prepare(&d, CurveOptions::default(), None).unwrap();
    let t = TopologyPlan::prepare(&s, 16_000_000, 64 * 1024 * 1024, None).unwrap();
    let text = write_local_candidate(&t, 16_000_000, 64 * 1024 * 1024, None).unwrap();
    let p = parse(&text, Options::default(), None).unwrap();
    let b = p.support_bindings(Options::default(), None).unwrap();
    let c = b
        .source_correspondence(&t, CorrespondenceOptions::default(), None)
        .unwrap();
    assert!(matches!(
        b.source_correspondence(
            &t,
            CorrespondenceOptions {
                derivative_limit_mm: c.row_derivative_error_mm() / 2.,
                ..CorrespondenceOptions::default()
            },
            None
        ),
        Err(Error::Unsupported)
    ));
    assert!(matches!(
        b.source_correspondence(
            &t,
            CorrespondenceOptions {
                retained_budget_bytes: c.retained_budget_bytes() - 1,
                ..CorrespondenceOptions::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
    let other = super::super::tests::graph(false, 4.3);
    let od = DomainPlan::prepare(&other, DomainOptions::default(), None).unwrap();
    let os = SurfacePlan::prepare(&od, CurveOptions::default(), None).unwrap();
    let ot = TopologyPlan::prepare(&os, 16_000_000, 64 * 1024 * 1024, None).unwrap();
    assert!(matches!(
        b.source_correspondence(&ot, CorrespondenceOptions::default(), None),
        Err(Error::Unsupported)
    ));
    let pole = p
        .tensors()
        .iter()
        .flat_map(|t| t.rows())
        .flat_map(|r| r[1..r.len() - 1].iter())
        .map(|id| p.point(*id).unwrap())
        .find(|p| {
            let z = p.coordinates()[2].value();
            z > 0. && z < 4.
        })
        .unwrap();
    let xyz = pole.coordinates();
    let changed = super::super::tests::argument_change(
        &text,
        pole.entity(),
        1,
        &format!(
            "(1.7976931348623157E308,{},{})",
            xyz[1].spelling(),
            xyz[2].spelling()
        ),
    );
    let pp = parse(&changed, Options::default(), None).unwrap();
    let bb = pp.support_bindings(Options::default(), None).unwrap();
    assert!(matches!(
        bb.source_correspondence(&t, CorrespondenceOptions::default(), None),
        Err(Error::Unsupported)
    ));
}
#[test]
fn actual_numeric_knot_leaves_do_not_mint_exact_source_parameter_identity() {
    let original = super::super::tests::file(false);
    let header = &original[..original.find("DATA;").unwrap() + 5];
    for (literal, wanted, accepted) in [
        ("0.", 0., true),
        ("1.", 1., true),
        ("5.0000000000000000E-1", 0.5, true),
        ("0.50000000000000001", 0.5, false),
        ("0.5", 0.5, false),
    ] {
        let input = format!("{header}#1=ACTUAL_REAL({literal});ENDSEC;END-ISO-10303-21;");
        let entities = part21::parse(&input).unwrap();
        let e = simple(&entities, 1).unwrap();
        let mut w = Ledger {
            used: 0,
            remaining: 16_000_000,
            interrupt: None,
        };
        let mut n = Numeric::new(e, &mut w).unwrap();
        let actual = n.scalar(&e.args[0]).unwrap();
        n.finish().unwrap();
        assert_eq!(knot(actual, wanted), accepted);
    }
    assert!(matches!(sum(f64::MAX, f64::MAX), Err(Error::Unsupported)));
}
fn rename(text: &str, a: u32, b: u32) -> String {
    let mut out = String::with_capacity(text.len());
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
            out.push_str(&format!(
                "#{}",
                if id == a {
                    b
                } else if id == b {
                    a
                } else {
                    id
                }
            ));
            i = end;
        } else {
            out.push(bytes[i] as char);
            i += 1;
        }
    }
    out
}
#[test]
fn actual_identity_renaming_and_wire_arities_never_infer_source_owners() {
    let body = super::super::tests::graph(false, 4.);
    let d = DomainPlan::prepare(&body, DomainOptions::default(), None).unwrap();
    let s = SurfacePlan::prepare(&d, CurveOptions::default(), None).unwrap();
    let t = TopologyPlan::prepare(&s, 16_000_000, 64 * 1024 * 1024, None).unwrap();
    let text = write_local_candidate(&t, 16_000_000, 64 * 1024 * 1024, None).unwrap();
    let p = parse(&text, Options::default(), None).unwrap();
    let row = p
        .edges()
        .iter()
        .position(|e| p.curve(e.curve_entity()).unwrap().degree() == 11)
        .unwrap();
    let ruling = p
        .edges()
        .iter()
        .position(|e| p.curve(e.curve_entity()).unwrap().degree() == 1)
        .unwrap();
    let pair = (0..s.rows().len())
        .flat_map(|a| (a + 1..s.rows().len()).map(move |b| (a, b)))
        .find(|&(a, b)| {
            t.edges()[a]
                .ends
                .iter()
                .any(|v| t.edges()[b].ends.contains(v))
        })
        .unwrap();
    for (a, b) in [(row, ruling), pair] {
        let changed = rename(&text, p.edges()[a].entity(), p.edges()[b].entity());
        let pp = parse(&changed, Options::default(), None).unwrap();
        let bb = pp.support_bindings(Options::default(), None).unwrap();
        assert!(bb.material_bounds(Options::default(), None).is_ok());
        assert!(matches!(
            bb.source_correspondence(&t, CorrespondenceOptions::default(), None),
            Err(Error::Unsupported)
        ));
    }
    let faces = p.faces();
    let (a, b) = (0..faces.len())
        .flat_map(|a| (a + 1..faces.len()).map(move |b| (a, b)))
        .find(|&(a, b)| {
            matches!(faces[a].surface(), SurfaceIndex::Tensor(_))
                && matches!(faces[b].surface(), SurfaceIndex::Tensor(_))
                && faces[a].same_sense() == faces[b].same_sense()
                && faces[a].loops()[0].coedges().len() != faces[b].loops()[0].coedges().len()
        })
        .unwrap();
    let shell = p
        .entities
        .iter()
        .find_map(|(&id, r)| matches!(r.as_slice(),[e]if e.name=="CLOSED_SHELL").then_some(id))
        .unwrap();
    let mut ids = faces.iter().map(Face::entity).collect::<Vec<_>>();
    ids.swap(a, b);
    let changed = super::super::tests::argument_change(
        &text,
        shell,
        1,
        &format!(
            "({})",
            ids.iter()
                .map(|id| format!("#{id}"))
                .collect::<Vec<_>>()
                .join(",")
        ),
    );
    let pp = parse(&changed, Options::default(), None).unwrap();
    let bb = pp.support_bindings(Options::default(), None).unwrap();
    assert!(bb.material_bounds(Options::default(), None).is_ok());
    assert!(matches!(
        bb.source_correspondence(&t, CorrespondenceOptions::default(), None),
        Err(Error::Unsupported)
    ));
    let b = p.support_bindings(Options::default(), None).unwrap();
    let polls = std::cell::Cell::new(0);
    let stop = || {
        polls.set(polls.get() + 1);
        polls.get() == t.edges().len() + 3
    };
    assert!(matches!(
        b.source_correspondence(&t, CorrespondenceOptions::default(), Some(&stop)),
        Err(Error::Cancelled)
    ));
}
#[test]
fn actual_refined_carrier_grid_and_changed_literal_knots_remain_explicit() {
    let body = super::super::tests::graph(false, 4.);
    let d = DomainPlan::prepare(&body, DomainOptions::default(), None).unwrap();
    let coarse = SurfacePlan::prepare(&d, CurveOptions::default(), None).unwrap();
    let c = TopologyPlan::prepare(&coarse, 16_000_000, 64 * 1024 * 1024, None).unwrap();
    let fine = SurfacePlan::prepare(
        &d,
        CurveOptions {
            position_error_mm: 1e-11,
            derivative_error_mm: 1e-10,
            ..CurveOptions::default()
        },
        None,
    )
    .unwrap();
    assert!(fine.rows().iter().any(|r| r.knots().len() > 2));
    let t = TopologyPlan::prepare(&fine, 16_000_000, 64 * 1024 * 1024, None).unwrap();
    let text = write_local_candidate(&t, 16_000_000, 64 * 1024 * 1024, None).unwrap();
    let p = parse(&text, Options::default(), None).unwrap();
    let b = p.support_bindings(Options::default(), None).unwrap();
    assert!(
        b.source_correspondence(&t, CorrespondenceOptions::default(), None)
            .is_ok()
    );
    assert!(matches!(
        b.source_correspondence(&c, CorrespondenceOptions::default(), None),
        Err(Error::Unsupported)
    ));
    let changed = text.replace("5.0000000000000000E-1", "0.50000000000000001");
    assert_ne!(changed, text);
    let pp = parse(&changed, Options::default(), None).unwrap();
    let bb = pp.support_bindings(Options::default(), None).unwrap();
    assert!(bb.material_bounds(Options::default(), None).is_ok());
    assert!(matches!(
        bb.source_correspondence(&t, CorrespondenceOptions::default(), None),
        Err(Error::Unsupported)
    ));
}
