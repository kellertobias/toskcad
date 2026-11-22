//! Parsed Original STEP operators, not fabricated cached curves or ranges.
use super::*;
use std::cell::Cell;
const BASE: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../multicad-app/tests/fixtures/exchange/corpus/occt-cone-frustum.step"
));
fn input(segment: Option<(f64, f64, f64)>) -> String {
    let mut records = String::new();
    let mut rows = Vec::new();
    let expanded: Vec<f64> = std::iter::repeat_n(0., 10)
        .chain(std::iter::repeat_n(0.5, 7))
        .chain(std::iter::repeat_n(1., 10))
        .collect();
    for u in 0..2 {
        let mut ids = Vec::new();
        for v in 0..17 {
            let knots = &expanded[v + 1..v + 10];
            let sum = knots.iter().sum::<f64>();
            let y = sum / 9.;
            let y2 = (sum * sum - knots.iter().map(|x| x * x).sum::<f64>()) / 72.;
            let id = 990010 + u * 17 + v;
            ids.push(format!("#{id}"));
            records += &format!(
                "#{id}=CARTESIAN_POINT('',({u}.,{y:.17e},{:.17e}));\n",
                u as f64 * y + 0.125 * y2
            );
        }
        rows.push(format!("({})", ids.join(",")));
    }
    records += &format!(
        "#990001=B_SPLINE_SURFACE_WITH_KNOTS('',1,9,({}),.UNSPECIFIED.,.F.,.F.,.F.,(2,2),(10,7,10),(0.,1.),(0.,0.5,1.),.UNSPECIFIED.);\n",
        rows.join(",")
    );
    let (uv, target) = if let Some((a, b, offset)) = segment {
        (
            vec![[0., 0.5], [1., 0.5]],
            vec![
                [a, 0.5, 0.5 * a + 0.03125 + offset],
                [b, 0.5, 0.5 * b + 0.03125 + offset],
            ],
        )
    } else {
        (
            vec![[0.125, 0.], [0.5, 0.], [0.875, 1.]],
            vec![
                [0.125, 0., 0.],
                [0.3125, 0., 0.],
                [0.5, 1. / 6., 1. / 48.],
                [0.6875, 0.5, 0.25],
                [0.875, 1., 1.],
            ],
        )
    };
    for (i, p) in uv.iter().enumerate() {
        records += &format!(
            "#{}=CARTESIAN_POINT('',({:.17e},{:.17e}));\n",
            990050 + i,
            p[0],
            p[1]
        );
    }
    for (i, p) in target.iter().enumerate() {
        records += &format!(
            "#{}=CARTESIAN_POINT('',({:.17e},{:.17e},{:.17e}));\n",
            990060 + i,
            p[0],
            p[1],
            p[2]
        );
    }
    let ids = |first, n| {
        (0..n)
            .map(|i| format!("#{}", first + i))
            .collect::<Vec<_>>()
            .join(",")
    };
    let p = uv.len() - 1;
    let q = target.len() - 1;
    records += &format!(
        "#990002=B_SPLINE_CURVE_WITH_KNOTS('',{p},({}),.UNSPECIFIED.,.F.,.F.,({},{}) ,(0.,1.),.UNSPECIFIED.);\n",
        ids(990050, uv.len()),
        p + 1,
        p + 1
    );
    records += &format!(
        "#990003=B_SPLINE_CURVE_WITH_KNOTS('',{q},({}),.UNSPECIFIED.,.F.,.F.,({},{}) ,(0.,1.),.UNSPECIFIED.);\n",
        ids(990060, target.len()),
        q + 1,
        q + 1
    );
    records += "#990080=PCURVE('',#990001,#990081);\n#990081=DEFINITIONAL_REPRESENTATION('',(#990002),#990082);\n#990082=(GEOMETRIC_REPRESENTATION_CONTEXT(2) PARAMETRIC_REPRESENTATION_CONTEXT() REPRESENTATION_CONTEXT('',''));\n#990083=SURFACE_CURVE('',#990003,(#990080),.PCURVE_S1.);\n#990084=VERTEX_POINT('',#990060);\n";
    records += &format!(
        "#990085=VERTEX_POINT('',#{});\n#990090=EDGE_CURVE('',#990084,#990085,#990083,.T.);\n",
        990060 + q
    );
    let at = BASE.rfind("ENDSEC;").unwrap();
    format!(
        "{}{}{}",
        &BASE[..at],
        records.to_ascii_uppercase(),
        &BASE[at..]
    )
}
fn actual(
    source: &Source<'_>,
    work: &mut Budget<'_>,
) -> (super::super::super::uv::Definition, Surface, Edge) {
    let binding = super::super::super::bindings(source, work)
        .unwrap()
        .into_iter()
        .find(|b| b.edge == 990090)
        .unwrap();
    let def = super::super::super::uv::definition(source, binding.pcurve, work).unwrap();
    let surface = source.surface(def.surface).unwrap();
    let edge = super::super::super::edge(source, binding, work).unwrap();
    (def, surface, edge)
}
#[test]
fn actual_original_curved_trace_crosses_a_nonrational_surface_join() {
    let text = input(None);
    let source = Source::parse(&text, 1e-7).unwrap();
    let (mut setup, mut work) = (
        Budget::new(512_000_000, None),
        Budget::new(512_000_000, None),
    );
    let (def, surface, edge) = actual(&source, &mut setup);
    assert!(
        matches!(certify(&source,990003,&def,&surface,&edge,1e-7,&mut work),Ok(Some(Proof::Certified(e))) if e<=1e-7)
    );
    // The exact UV V(t)=t² crosses the interior V=.5 knot at irrational t.
    // Alter an actual interior 3D coefficient while retaining both endpoints.
    let changed = text.replace("2.08333333333333322E-2", "1.02083333333333326E0");
    assert_ne!(changed, text);
    let source = Source::parse(&changed, 1e-7).unwrap();
    let (def, surface, edge) = actual(&source, &mut Budget::new(512_000_000, None));
    assert!(!matches!(
        certify(
            &source,
            990003,
            &def,
            &surface,
            &edge,
            1e-7,
            &mut Budget::new(512_000_000, None)
        ),
        Ok(Some(Proof::Certified(_)))
    ));
    let outside = text.replace(
        "8.75000000000000000E-1,1.00000000000000000E0",
        "1.25000000000000000E0,1.00000000000000000E0",
    );
    assert_ne!(outside, text);
    let source = Source::parse(&outside, 1e-7).unwrap();
    let (def, surface, edge) = actual(&source, &mut Budget::new(512_000_000, None));
    assert!(!matches!(
        certify(
            &source,
            990003,
            &def,
            &surface,
            &edge,
            1e-7,
            &mut Budget::new(512_000_000, None)
        ),
        Ok(Some(Proof::Certified(_)))
    ));
}
#[test]
fn finite_used_segment_requires_an_increasing_subsegment_inside_the_uv_image() {
    for (a, b, offset, positive) in [
        (0.2, 0.8, 0., true),
        (0.2, 1.2, 0., false),
        (0.2, 0.8, 1., false),
    ] {
        let text = input(Some((a, b, offset)));
        let source = Source::parse(&text, 1e-7).unwrap();
        let (def, surface, edge) = actual(&source, &mut Budget::new(512_000_000, None));
        let result = certify(
            &source,
            990003,
            &def,
            &surface,
            &edge,
            1e-7,
            &mut Budget::new(512_000_000, None),
        );
        if positive {
            assert!(matches!(result,Ok(Some(Proof::Certified(e))) if e<=1e-7));
        } else {
            assert!(matches!(result,Ok(Some(Proof::Disagreement(e))) if e>1e-7));
        }
    }
}
#[test]
fn actual_source_proof_work_frontier_observed_late_cancel_and_same_input_retry() {
    let text = input(None);
    let source = Source::parse(&text, 1e-7).unwrap();
    let (def, surface, edge) = actual(&source, &mut Budget::new(512_000_000, None));
    let run = |work: &mut Budget<'_>| certify(&source, 990003, &def, &surface, &edge, 1e-7, work);
    let mut full = Budget::new(512_000_000, None);
    let result = run(&mut full).unwrap();
    let used = full.used();
    assert!(matches!(result, Some(Proof::Certified(_))));
    assert!(matches!(
        run(&mut Budget::new(used - 1, None)),
        Err(ReadError::Limit)
    ));
    let retry = run(&mut Budget::new(used, None)).unwrap();
    let bits = |proof: &Option<Proof>| match proof {
        Some(Proof::Certified(e)) => e.to_bits(),
        _ => panic!("no certificate"),
    };
    assert_eq!(bits(&result), bits(&retry));
    let calls = Cell::new(0);
    let observe = || {
        calls.set(calls.get() + 1);
        false
    };
    run(&mut Budget::new(used, Some(&observe))).unwrap();
    let last = calls.get();
    calls.set(0);
    let cancel = || {
        calls.set(calls.get() + 1);
        calls.get() == last
    };
    let mut stopped = Budget::new(used, Some(&cancel));
    assert!(matches!(run(&mut stopped), Err(ReadError::Limit)));
    assert!(stopped.cancelled.get());
    assert_eq!(
        bits(&run(&mut Budget::new(used, None)).unwrap()),
        bits(&retry)
    );
}
