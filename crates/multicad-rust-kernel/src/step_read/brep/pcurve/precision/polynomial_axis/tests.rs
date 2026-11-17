use super::*;
const BASE: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../multicad-app/tests/fixtures/exchange/corpus/occt-cone-frustum.step"
));
fn input(reverse: bool, constant: bool) -> String {
    let mut extra = String::new();
    let mut ids = Vec::new();
    for i in 0..6 {
        let id = 900010 + i;
        let u = if reverse { 5 - i } else { i };
        let v = if !constant && i == 3 {
            "1.000000000000000000000001"
        } else {
            "1."
        };
        extra += &format!("#{id}=CARTESIAN_POINT('',({u}.,{v}));\n");
        ids.push(format!("#{id}"));
    }
    extra += &format!(
        "#900001=B_SPLINE_CURVE_WITH_KNOTS('',3,({}),.UNSPECIFIED.,.F.,.F.,(4,2,4),(0.,0.25,1.),.UNSPECIFIED.);\n#900100=PCURVE('',#31,#900101);\n#900101=DEFINITIONAL_REPRESENTATION('',(#900001),#900102);\n#900102=(GEOMETRIC_REPRESENTATION_CONTEXT(2)PARAMETRIC_REPRESENTATION_CONTEXT()REPRESENTATION_CONTEXT('',''));\n",
        ids.join(",")
    );
    let at = BASE.rfind("ENDSEC;").unwrap();
    format!("{}{}{}", &BASE[..at], extra, &BASE[at..])
}
fn run(source: &Source<'_>, budget: &mut Budget<'_>) -> Result<Option<AxisSegment>> {
    let d = super::super::super::uv::definition(source, 900100, budget)?;
    segment(source, &d, budget)
}
#[test]
fn original_multispan_polynomial_axis_earns_both_orientations_and_rejects_alias() {
    for reverse in [false, true] {
        let text = input(reverse, true);
        let source = Source::parse(&text, 1e-7).unwrap();
        let mut budget = Budget::new(512_000_000, None);
        let axis = run(&source, &mut budget).unwrap().unwrap();
        assert_eq!(axis.moving, 0);
        assert_eq!(
            axis.endpoints[0][0].midpoint(),
            if reverse { 5. } else { 0. }
        );
        assert_eq!(
            axis.endpoints[1][0].midpoint(),
            if reverse { 0. } else { 5. }
        );
        let used = budget.used();
        assert!(matches!(
            run(&source, &mut Budget::new(used - 1, None)),
            Err(ReadError::Limit)
        ));
        assert!(
            run(&source, &mut Budget::new(used, None))
                .unwrap()
                .is_some()
        );
        let callbacks = std::cell::Cell::new(0);
        let observe = || {
            callbacks.set(callbacks.get() + 1);
            false
        };
        run(&source, &mut Budget::new(512_000_000, Some(&observe))).unwrap();
        let last = callbacks.get();
        callbacks.set(0);
        let cancel = || {
            callbacks.set(callbacks.get() + 1);
            callbacks.get() == last
        };
        let mut b = Budget::new(512_000_000, Some(&cancel));
        assert!(matches!(run(&source, &mut b), Err(ReadError::Limit)));
        assert!(b.cancelled.get());
        assert!(
            run(&source, &mut Budget::new(512_000_000, None))
                .unwrap()
                .is_some()
        );
    }
    let bad = input(false, false);
    let source = Source::parse(&bad, 1e-7).unwrap();
    assert!(
        run(&source, &mut Budget::new(512_000_000, None))
            .unwrap()
            .is_none()
    );
    let text = input(false, true).replace(
        "#900012=CARTESIAN_POINT('',(2.,1.))",
        "#900012=CARTESIAN_POINT('',(-100.,1.))",
    );
    let source = Source::parse(&text, 1e-7).unwrap();
    assert!(
        run(&source, &mut Budget::new(512_000_000, None))
            .unwrap()
            .is_none()
    );
}
#[test]
fn installed_dip_original_four_circle_trims_are_measured_separately() {
    let path = "/Applications/KiCad/KiCad.app/Contents/SharedSupport/3dmodels/Package_DIP.3dshapes/DIP-32_W7.62mm.step";
    let Ok(text) = std::fs::read_to_string(path) else {
        eprintln!("unavailable: external KiCad DIP fixture");
        return;
    };
    let source = Source::parse(&text, 1e-7).unwrap();
    let mut budget = Budget::new(512_000_000, None);
    let mut checked = 0;
    for binding in super::super::super::bindings(&source, &mut budget).unwrap() {
        if ![37456, 37500, 39276, 39320].contains(&binding.pcurve) {
            continue;
        }
        let e = super::super::super::edge(&source, binding, &mut budget).unwrap();
        let d = super::super::super::uv::definition(&source, binding.pcurve, &mut budget).unwrap();
        let surface = source.surface(d.surface).unwrap();
        let proof = super::super::super::finite_axis::certify(
            &source,
            binding.curve,
            &d,
            &surface,
            &e,
            1e-7,
            &mut budget,
        )
        .unwrap();
        assert!(
            matches!(proof,Some(super::super::super::Proof::Certified(error))if error<=1e-7),
            "unqualified circle edge{}pc{}",
            binding.edge,
            binding.pcurve
        );
        checked += 1;
    }
    assert_eq!(checked, 4);
    eprintln!("actual DIP circle trims={checked} work={}", budget.used());
}
