use super::*;
const BASE: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../multicad-app/tests/fixtures/exchange/corpus/occt-cone-frustum.step"
));
fn input(end: bool) -> String {
    let z = if end { "-1." } else { "0." };
    let z1 = if end { "0." } else { "1." };
    let u = if end { "1." } else { "0." };
    let extra = format!(
        "\n#900001=B_SPLINE_CURVE_WITH_KNOTS('',2,(#900010,#900011,#900012),.UNSPECIFIED.,.F.,.F.,(3,3),(0.,6.323173064558),.UNSPECIFIED.);\n#900010=CARTESIAN_POINT('',(0.,0.,0.));\n#900011=CARTESIAN_POINT('',(0.5,0.,0.));\n#900012=CARTESIAN_POINT('',(1.,1.,0.));\n#900020=CARTESIAN_POINT('',(0.,0.,0.));\n#900021=CARTESIAN_POINT('',(0.5,0.25,0.));\n#900030=VERTEX_POINT('',#900020);\n#900031=VERTEX_POINT('',#900021);\n#900032=EDGE_CURVE('',#900030,#900031,#900001,.T.);\n#900040=B_SPLINE_SURFACE_WITH_KNOTS('',1,2,((#900050,#900051,#900052),(#900053,#900054,#900055)),.UNSPECIFIED.,.F.,.F.,.F.,(2,2),(3,3),(0.,1.),(0.,6.323173064558),.UNSPECIFIED.);\n#900050=CARTESIAN_POINT('',(0.,0.,{z}));\n#900051=CARTESIAN_POINT('',(0.5,0.,{z}));\n#900052=CARTESIAN_POINT('',(1.,1.,{z}));\n#900053=CARTESIAN_POINT('',(0.,0.,{z1}));\n#900054=CARTESIAN_POINT('',(0.5,0.,{z1}));\n#900055=CARTESIAN_POINT('',(1.,1.,{z1}));\n#900100=PCURVE('',#900040,#900101);\n#900101=DEFINITIONAL_REPRESENTATION('',(#900102),#900103);\n#900103=(GEOMETRIC_REPRESENTATION_CONTEXT(2)PARAMETRIC_REPRESENTATION_CONTEXT()REPRESENTATION_CONTEXT('',''));\n#900102=B_SPLINE_CURVE_WITH_KNOTS('',1,(#900110,#900111),.UNSPECIFIED.,.F.,.F.,(2,2),(0.,3.161586532279),.UNSPECIFIED.);\n#900110=CARTESIAN_POINT('',({u},0.));\n#900111=CARTESIAN_POINT('',({u},3.161586532279));\n"
    );
    let at = BASE.rfind("ENDSEC;").unwrap();
    format!("{}{}{}", &BASE[..at], extra, &BASE[at..])
}
fn pair(
    source: &Source<'_>,
    work: &mut Budget<'_>,
) -> (super::super::super::uv::Definition, Surface, Edge) {
    let d = super::super::super::uv::definition(source, 900100, work).unwrap();
    let s = source.surface(d.surface).unwrap();
    let e = super::super::super::edge(
        source,
        super::super::super::Binding {
            edge: 900032,
            curve: 900001,
            pcurve: 900100,
        },
        work,
    )
    .unwrap();
    (d, s, e)
}
#[test]
fn original_nonbinary_domain_partial_boundary_both_sides_and_real_refusals() {
    for end in [false, true] {
        let text = input(end);
        let source = Source::parse(&text, 1e-7).unwrap();
        let mut work = Budget::new(512_000_000, None);
        let (d, s, e) = pair(&source, &mut work);
        let proof = certify(&source, 900001, &d, &s, &e, 1e-7, &mut work).unwrap();
        assert!(
            matches!(proof,Some(super::super::super::Proof::Certified(v))if v<1e-9),
            "unqualified portable finite boundary"
        );
        // Genuine changed ORIGINAL row, same defining schema/knots/endpoints.
        let id = if end { 900054 } else { 900051 };
        let old = format!("#{id}=CARTESIAN_POINT('',(0.5,0.,0.))");
        let bad = text.replace(&old, &format!("#{id}=CARTESIAN_POINT('',(0.5,0.,0.001))"));
        assert_ne!(bad, text);
        let source = Source::parse(&bad, 1e-7).unwrap();
        let (d, s, e) = pair(&source, &mut work);
        assert!(
            certify(&source, 900001, &d, &s, &e, 1e-7, &mut work)
                .unwrap()
                .is_none()
        );
    }
    let text = input(false);
    for bad in [
        text.replace(
            "#900111=CARTESIAN_POINT('',(0.,3.161586532279))",
            "#900111=CARTESIAN_POINT('',(0.00000000000000000001,3.161586532279))",
        ),
        text.replace(
            "#900110=CARTESIAN_POINT('',(0.,0.))",
            "#900110=CARTESIAN_POINT('',(0.,1E-1000))",
        ),
    ] {
        assert_ne!(bad, text);
        let source = Source::parse(&bad, 1e-7).unwrap();
        let mut b = Budget::new(512_000_000, None);
        let (d, s, e) = pair(&source, &mut b);
        assert!(
            certify(&source, 900001, &d, &s, &e, 1e-7, &mut b)
                .unwrap()
                .is_none()
        );
    }
}
#[test]
fn measured_work_frontier_observed_last_cancellation_and_original_retry() {
    let text = input(false);
    let source = Source::parse(&text, 1e-7).unwrap();
    let (d, s, e) = pair(&source, &mut Budget::new(512_000_000, None));
    let run = |b: &mut Budget<'_>| certify(&source, 900001, &d, &s, &e, 1e-7, b);
    let mut b = Budget::new(512_000_000, None);
    let expected = run(&mut b).unwrap();
    assert!(expected.is_some());
    let used = b.used();
    assert!(matches!(
        run(&mut Budget::new(used - 1, None)),
        Err(ReadError::Limit)
    ));
    assert!(run(&mut Budget::new(used, None)).unwrap().is_some());
    let calls = std::cell::Cell::new(0);
    let observe = || {
        calls.set(calls.get() + 1);
        false
    };
    run(&mut Budget::new(512_000_000, Some(&observe))).unwrap();
    let n = calls.get();
    calls.set(0);
    let cancel = || {
        calls.set(calls.get() + 1);
        calls.get() == n
    };
    let mut b = Budget::new(512_000_000, Some(&cancel));
    assert!(matches!(run(&mut b), Err(ReadError::Limit)));
    assert!(b.cancelled.get());
    let mut b = Budget::new(512_000_000, None);
    assert!(run(&mut b).unwrap().is_some());
    assert_eq!(b.used(), used);
}
#[test]
fn actual_sma_partial2453_both_original_surface_bindings() {
    let path = "/Applications/KiCad/KiCad.app/Contents/SharedSupport/3dmodels/Connector_Coaxial.3dshapes/SMA_Samtec_SMA-J-P-H-ST-EM1_EdgeMount.step";
    let Ok(text) = std::fs::read_to_string(path) else {
        eprintln!("unavailable: external KiCad SMA fixture");
        return;
    };
    let source = Source::parse(&text, 1e-7).unwrap();
    let mut b = Budget::new(512_000_000, None);
    for pc in [2489, 2495] {
        let d = super::super::super::uv::definition(&source, pc, &mut b).unwrap();
        let s = source.surface(d.surface).unwrap();
        let e = super::super::super::edge(
            &source,
            super::super::super::Binding {
                edge: 2453,
                curve: 2457,
                pcurve: pc,
            },
            &mut b,
        )
        .unwrap();
        let proof = certify(&source, 2457, &d, &s, &e, 1e-7, &mut b).unwrap();
        assert!(
            matches!(proof,Some(super::super::super::Proof::Certified(v))if v<1e-7),
            "unqualified actual PC{pc}"
        );
        eprintln!(
            "actual SMA edge2453 PC{pc} upper={} aggregatework={}",
            match proof {
                Some(super::super::super::Proof::Certified(v)) => v,
                _ => panic!("unqualified PC{pc}"),
            },
            b.used()
        );
    }
}

#[test]
fn actual_raw_boundary_and_natural_domain_shadows_cannot_be_cache_authority() {
    let text = input(true);
    let bad = text
        .replace(
            "#900110=CARTESIAN_POINT('',(1.,0.))",
            "#900110=CARTESIAN_POINT('',(1.0000000000000000001,0.))",
        )
        .replace(
            "#900111=CARTESIAN_POINT('',(1.,3.161586532279))",
            "#900111=CARTESIAN_POINT('',(1.0000000000000000001,3.161586532279))",
        );
    assert_ne!(text, bad);
    let source = Source::parse(&bad, 1e-7).unwrap();
    let mut b = Budget::new(512_000_000, None);
    let (d, s, e) = pair(&source, &mut b);
    assert!(
        certify(&source, 900001, &d, &s, &e, 1e-7, &mut b)
            .unwrap()
            .is_none()
    );
    let text = input(false);
    let bad = text
        .replace(
            "#900110=CARTESIAN_POINT('',(0.,0.))",
            "#900110=CARTESIAN_POINT('',(0.,-1E-20))",
        )
        .replace("(2,2),(0.,3.161586532279)", "(2,2),(-1E-20,3.161586532279)");
    assert_ne!(text, bad);
    let source = Source::parse(&bad, 1e-7).unwrap();
    let mut b = Budget::new(512_000_000, None);
    let (d, s, e) = pair(&source, &mut b);
    assert!(
        certify(&source, 900001, &d, &s, &e, 1e-7, &mut b)
            .unwrap()
            .is_none()
    );
}

#[test]
fn opposite_original_uv_sense_is_not_a_finite_boundary_permission() {
    let text = input(false);
    let bad = text.replace("1,(#900110,#900111)", "1,(#900111,#900110)");
    assert_ne!(text, bad);
    let source = Source::parse(&bad, 1e-7).unwrap();
    let mut b = Budget::new(512_000_000, None);
    let (d, s, e) = pair(&source, &mut b);
    assert!(
        certify(&source, 900001, &d, &s, &e, 1e-7, &mut b)
            .unwrap()
            .is_none()
    );
}
