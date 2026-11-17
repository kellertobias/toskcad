use super::*;
use crate::step_read::brep::pcurve::uv;
const TOL: f64 = 1e-7;
const BASE: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../multicad-app/tests/fixtures/exchange/corpus/occt-cone-frustum.step"
));
// Geometry records are byte-exact from immutable original M6 SHA256
// 96aa15409ae6b5a2ef87e5ae7aabfb0636f6518110cfe63d6116a4aba2357907.
// Only the ordinary 2D context is remapped and the unrelated association
// omitted. This trim island is not claimed to be the complete original solid.
const RECORDS: &str = include_str!("fixtures/m6-original-trim.stepfrag");
fn input() -> String {
    let at = BASE.rfind("ENDSEC;").unwrap();
    format!("{}\n{}{}", &BASE[..at], RECORDS, &BASE[at..])
}
fn run(
    text: &str,
    resource: usize,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<(Option<Proof>, usize)> {
    let source = Source::parse(text, TOL)?;
    let mut setup = Budget::new(512_000_000, None);
    let d = uv::definition(&source, 24271, &mut setup)?;
    let surface = source.surface(d.surface)?;
    let curve = source.curve(24367, 0)?;
    let points = [source.point(24068)?, source.point(24070)?];
    let edge = Edge::new(curve, [24069, 24071], points, TOL).map_err(|_| ReadError::Unsupported)?;
    let mut work = Budget::new(resource, interrupt);
    let result = certify(&source, 24367, &d, &surface, &edge, TOL, &mut work)?;
    Ok((result, work.used()))
}
fn lower(text: &str) -> f64 {
    match run(text, 512_000_000, None).unwrap().0 {
        Some(Proof::Disagreement(v)) => v,
        _ => panic!("entire original slabs should separate"),
    }
}
#[test]
fn actual_original_m6_trim_images_have_a_whole_plane_gap() {
    let text = input();
    let value = lower(&text);
    let expected = num_rational::BigRational::new(
        16_568_542_494_923_788_i64.into(),
        100_000_000_000_000_000_i64.into(),
    );
    assert!(value > 0.16568 && value < 0.16569);
    assert!(num_rational::BigRational::from_float(value).unwrap() <= expected);
    let (_, used) = run(&text, 512_000_000, None).unwrap();
    assert!(matches!(run(&text, used - 1, None), Err(ReadError::Limit)));
    assert!(matches!(
        run(&text, used, None).unwrap().0,
        Some(Proof::Disagreement(_))
    ));
    let calls = std::cell::Cell::new(0);
    let observe = || {
        calls.set(calls.get() + 1);
        false
    };
    run(&text, 512_000_000, Some(&observe)).unwrap();
    let final_call = calls.get();
    calls.set(0);
    let cancel = || {
        calls.set(calls.get() + 1);
        calls.get() == final_call
    };
    assert!(matches!(
        run(&text, 512_000_000, Some(&cancel)),
        Err(ReadError::Limit)
    ));
    assert_eq!(lower(&text).to_bits(), value.to_bits());
}
fn change_v(text: &str, all: bool, value: &str) -> String {
    let source = Source::parse(text, TOL).unwrap();
    let mut work = Budget::new(512_000_000, None);
    let fields = fields(&source, 24269, 1, &mut work).unwrap();
    let mut output = text.to_owned();
    for (i, control) in list(fields.controls).unwrap().iter().enumerate() {
        if !all && i > 0 {
            break;
        }
        let id = reference(control).unwrap();
        let [_, values] = source.args(id, "CARTESIAN_POINT").unwrap() else {
            panic!()
        };
        let values = list(values).unwrap();
        let u = spelling(&source, id, &values[0], &mut work).unwrap();
        let old = format!("#{id}={};", source.simple(id).unwrap().source.unwrap());
        let new = format!("#{id}=CARTESIAN_POINT('',({u},{value}));");
        assert!(output.contains(&old));
        output = output.replacen(&old, &new, 1);
    }
    output
}
#[test]
fn overlapping_or_varying_planes_cannot_borrow_a_negative_or_positive_certificate() {
    let text = input();
    let overlap = change_v(&text, true, "4.0000000000000002E-1");
    assert!(run(&overlap, 512_000_000, None).unwrap().0.is_none());
    let varying = change_v(&text, false, "0.");
    assert!(run(&varying, 512_000_000, None).unwrap().0.is_none());
    assert!(matches!(
        run(&text, 512_000_000, None).unwrap().0,
        Some(Proof::Disagreement(_))
    ));
}
#[test]
fn original_slabs_survive_an_independent_noncardinal_literal_pose() {
    let text = input();
    let source = Source::parse(&text, TOL).unwrap();
    let mut posed = text.clone();
    let columns = [[0.8, 0.6, 0.], [-0.36, 0.48, 0.8], [0.48, -0.64, 0.6]];
    let mut ids = std::collections::BTreeSet::new();
    for record in RECORDS.lines() {
        if let Some((id, _)) = record.strip_prefix('#').and_then(|r| r.split_once('=')) {
            ids.insert(id.parse::<u32>().unwrap());
        }
    }
    for id in ids {
        let [entity] = source.components(id).unwrap() else {
            continue;
        };
        if !matches!(entity.name, "CARTESIAN_POINT" | "DIRECTION") {
            continue;
        }
        let [_, values] = entity.args.as_slice() else {
            panic!()
        };
        let values = list(values).unwrap();
        if values.len() != 3 {
            continue;
        }
        let p: [f64; 3] = std::array::from_fn(|i| {
            let Value::Number(v) = values[i] else {
                panic!()
            };
            v
        });
        let offset = if entity.name == "CARTESIAN_POINT" {
            [17., -8., 13.]
        } else {
            [0.; 3]
        };
        let q: [f64; 3] =
            std::array::from_fn(|i| offset[i] + (0..3).map(|j| columns[j][i] * p[j]).sum::<f64>());
        let old = format!("#{id}={};", entity.source.unwrap());
        let new = format!(
            "#{id}={}('',({:.17E},{:.17E},{:.17E}));",
            entity.name, q[0], q[1], q[2]
        );
        assert!(posed.contains(&old));
        posed = posed.replacen(&old, &new, 1);
    }
    assert!(lower(&posed) > 0.16568);
}
