use super::test_fixture::{engineering, engineering_height};
use super::*;
use crate::{
    axial_pipe_thread::authored::hardware::shell::{self, occupied},
    radial_tensor::{Curve, Surface},
    step_read::pipe_hardware as reader,
};
#[test]
fn actual_engineering_head_file_retains_non_power_of_two_domains_after_arena_drop() {
    let source = engineering();
    let candidate = propose(&source, Options::default(), None).unwrap();
    let proposal = candidate.receipt();
    let path = std::env::var_os("MULTICAD_CODEC_ARTIFACT_DIR")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(std::env::temp_dir)
        .join(format!("engineering-npt-head-{}.step", std::process::id()));
    std::fs::write(&path, candidate.as_str()).unwrap();
    drop(candidate);
    let actual = std::fs::read_to_string(&path).unwrap();
    let resident =
        source.receipt().retained_result_bytes_upper_bound + actual.capacity() - actual.len();
    println!(
        "ENGINEERING_PROPOSAL bytes{} resident{resident} receipt{proposal:?}",
        actual.len()
    );
    let qualified = qualify_actual_file(
        &source,
        &actual,
        1e-7,
        1e-6,
        QualificationOptions {
            additional_resident_bytes: actual.capacity() - actual.len(),
            ..Default::default()
        },
        None,
    )
    .unwrap();
    let qualification = qualified.receipt();
    assert_eq!(qualification.raw_operator_position_delta_mm, 0.);
    assert_eq!(qualification.mapped_operator_derivative_delta_mm, 0.);
    assert!(qualification.matched_original_records > 30_000);
    assert!(qualification.aggregate_work_reserved <= 272_000_000);
    assert!(qualification.correspondence_peak_bytes_reserved <= 64 * 1024 * 1024);
    println!("QUALIFIED_HEAD {qualification:?}");
    let decoded = qualified.into_material();
    println!(
        "ENGINEERING_HEAD bytes{} proposal{proposal:?} reader{:?} retained{}",
        actual.len(),
        decoded.receipt(),
        decoded.retained_bytes_upper_bound()
    );
    assert_eq!(decoded.original_input(), actual);
    assert_eq!(decoded.material().head_height_mm(), [6., 9.]);
    assert!(decoded.parameter_maps().iter().any(|map|matches!(map.mapping(),reader::ParameterMapping::Spline{spans}if !spans.is_power_of_two())));
    assert!(decoded.parameter_maps().iter().any(
        |map| matches!(map.mapping(),reader::ParameterMapping::Tensor{spans}if spans.contains(&47))
    ));
    assert!(decoded.receipt().aggregate_work_reserved <= 224_000_000);
    assert!(decoded.receipt().parser_peak_bytes_reserved <= 64 * 1024 * 1024);
    assert!(decoded.receipt().decode_peak_bytes_reserved <= 64 * 1024 * 1024);
    assert!(decoded.receipt().occupied_peak_bytes_reserved <= 64 * 1024 * 1024);
    // Independent cone slant/axial chain witness. Forgetting sec(alpha)
    // changes the axial tangent even though all scalar tables would match.
    let cone = decoded
        .parameter_maps()
        .iter()
        .find(|m| matches!(m.mapping(), reader::ParameterMapping::Cone { .. }))
        .unwrap();
    let reader::Target::Surface(index) = cone.target() else {
        panic!("cone surface")
    };
    let Surface::Cone {
        semi_angle_rad: angle,
        ..
    } = decoded.material().shell().input().surfaces[index]
    else {
        panic!("actual cone")
    };
    let theta = 0.37f64;
    let raw = [
        theta.cos() * angle.sin(),
        theta.sin() * angle.sin(),
        angle.cos(),
    ];
    let pull = cone
        .reference_derivative_bounds(raw.map(|x| [x.next_down(), x.next_up()]), 1)
        .unwrap();
    let expected = [theta.cos() * angle.tan(), theta.sin() * angle.tan(), 1.];
    for k in 0..3 {
        assert!(pull[k][0] <= expected[k] && expected[k] <= pull[k][1]);
    }
    assert!((raw[2] - expected[2]).abs() > 1e-4);
    assert!(
        decoded
            .material()
            .capped_thread()
            .input()
            .curves
            .iter()
            .any(|c| matches!(c, Curve::Spline { spans: 47, .. }))
    );
    independent_spline_chain(&source, &decoded);
    independent_mapped_geometry_consumers(&decoded);
    drop(source);
    drop(actual);
    // Actual record bytes and owned original operators survive every parser
    // owner and original engineering-source owner, not only cached metadata.
    assert!(decoded.original_input().starts_with("ISO-10303-21;"));
    assert_eq!(
        decoded
            .material()
            .capped_thread()
            .classify([0., 0., 3.], None)
            .unwrap(),
        crate::radial_tensor::Location::Inside
    );
}
#[test]
fn actual_file_changes_cannot_borrow_engineering_source_or_scalar_cache_admission() {
    use crate::{
        axial_pipe_thread::Error,
        step_read::literal_member_binding::{RetainedRecords, Token, token},
    };
    let source = engineering();
    let c = propose(&source, Options::default(), None).unwrap();
    let components = (0..c.record_count()).map(|i| c.component(i, 0).unwrap());
    let spline = components
        .clone()
        .find(|s| s.starts_with("B_SPLINE_CURVE_WITH_KNOTS") && s.ends_with(",47),.UNSPECIFIED.)"))
        .unwrap()
        .to_owned();
    let tensor = components
        .clone()
        .find(|s| s.starts_with("B_SPLINE_SURFACE_WITH_KNOTS"))
        .unwrap();
    let tensor_prefix = "B_SPLINE_SURFACE_WITH_KNOTS('',9,9,(";
    assert!(tensor.starts_with(tensor_prefix));
    let first_row_start = tensor_prefix.len();
    let first_row_end = first_row_start + tensor[first_row_start..].find(')').unwrap() + 1;
    let first_row = &tensor[first_row_start..first_row_end];
    let row_needle = format!("{tensor_prefix}{first_row}");
    let vector = components
        .clone()
        .find(|s| s.starts_with("VECTOR("))
        .unwrap()
        .to_owned();
    let inner = components
        .clone()
        .find(|s| s.starts_with("FACE_BOUND("))
        .unwrap()
        .to_owned();
    let circle = components
        .clone()
        .find(|s| s.starts_with("CIRCLE("))
        .unwrap();
    let refs = |mut text: &str| {
        let mut ids = Vec::new();
        while let Some(t) = token(&mut text) {
            if let Token::Reference(id) = t {
                ids.push(id);
            }
        }
        ids
    };
    let axis = refs(circle)[0];
    let axis = c.component(c.record_index(axis).unwrap(), 0).unwrap();
    let direction = refs(axis)[2];
    let direction = c
        .component(c.record_index(direction).unwrap(), 0)
        .unwrap()
        .to_owned();
    let mut point_leaf = None;
    for component in components.filter(|s| s.starts_with("CARTESIAN_POINT(")) {
        let mut rest = component;
        while let Some(t) = token(&mut rest) {
            if let Token::Literal(t) = t
                && t.contains(".E-")
                && t.parse::<f64>().is_ok_and(|x| x != 0.)
            {
                point_leaf = Some(t.to_owned());
                break;
            }
        }
        if point_leaf.is_some() {
            break;
        }
    }
    let leaf = point_leaf.unwrap();
    let (m, e) = leaf.split_once(".E-").unwrap();
    let alias = format!("{m}1.E-{}", e.parse::<usize>().unwrap() + 1);
    assert_ne!(leaf, alias);
    assert_eq!(leaf.parse::<f64>().unwrap(), alias.parse::<f64>().unwrap());
    let cases = vec![
        (leaf, alias, Error::Unsupported),
        (
            "B_SPLINE_CURVE_WITH_KNOTS('',9,".into(),
            "B_SPLINE_CURVE_WITH_KNOTS('',8,".into(),
            Error::Unsupported,
        ),
        (
            "B_SPLINE_CURVE_WITH_KNOTS('',9,".into(),
            "B_SPLINE_CURVE_WITH_KNOTS('',9.E0,".into(),
            Error::Invalid,
        ),
        (
            "B_SPLINE_SURFACE_WITH_KNOTS('',9,9,".into(),
            "B_SPLINE_SURFACE_WITH_KNOTS('',8,9,".into(),
            Error::Unsupported,
        ),
        (
            row_needle.clone(),
            format!("{tensor_prefix}$"),
            Error::Invalid,
        ),
        (
            row_needle,
            format!("{tensor_prefix}(#1)"),
            Error::Unsupported,
        ),
        (
            vector.clone(),
            vector.replacen(",1.)", ",2.)", 1),
            Error::Unsupported,
        ),
        (
            spline.clone(),
            spline.replacen(",47),.UNSPECIFIED.)", ",48),.UNSPECIFIED.)", 1),
            Error::Unsupported,
        ),
        (
            inner.clone(),
            inner.replacen(".T.)", ".F.)", 1),
            Error::Invalid,
        ),
        (
            direction,
            "DIRECTION('',(0.,1.E0,0.))".into(),
            Error::Unsupported,
        ),
        (
            "GEOMETRIC_REPRESENTATION_CONTEXT(3)".into(),
            "GEOMETRIC_REPRESENTATION_CONTEXT(3.000000000000000000001)".into(),
            Error::Invalid,
        ),
        (
            "SI_UNIT(.MILLI.,.METRE.)".into(),
            "SI_UNIT(.CENTI.,.METRE.)".into(),
            Error::Unsupported,
        ),
    ];
    let case_bytes = cases
        .iter()
        .map(|(a, b, _)| a.capacity() + b.capacity())
        .sum::<usize>()
        + cases.capacity() * std::mem::size_of::<(String, String, Error)>();
    let actual = c.into_string();
    let outside =
        source.receipt().retained_result_bytes_upper_bound + actual.capacity() + case_bytes;
    for (old, new, error) in &cases {
        let changed = actual.replacen(old, new, 1);
        assert!(
            changed != actual,
            "mutation must change the actual STEP record: {old}"
        );
        assert!(
            matches!(reader::read(&changed,1e-7,1e-6,reader::Options{additional_resident_bytes:outside,..Default::default()},None),Err(e)if e==*error)
        );
    }
    let at = actual.rfind("ENDSEC;").unwrap();
    let orphan = format!(
        "{}#999999=VERTEX_POINT('',#1);\n{}",
        &actual[..at],
        &actual[at..]
    );
    assert!(matches!(
        reader::read(
            &orphan,
            1e-7,
            1e-6,
            reader::Options {
                additional_resident_bytes: outside,
                ..Default::default()
            },
            None
        ),
        Err(Error::Invalid)
    ));
    drop(orphan);
    drop(actual);
    let outside = source.receipt().retained_result_bytes_upper_bound + case_bytes;
    let different = engineering_height(4.);
    let mut changed = propose(
        &different,
        Options {
            additional_resident_bytes: outside,
            ..Default::default()
        },
        None,
    )
    .unwrap()
    .into_string();
    changed.shrink_to_fit();
    drop(different);
    let admitted = reader::read(
        &changed,
        1e-7,
        1e-6,
        reader::Options {
            additional_resident_bytes: outside,
            ..Default::default()
        },
        None,
    )
    .unwrap();
    assert_eq!(admitted.material().head_height_mm(), [6., 10.]);
    drop(admitted);
    assert!(matches!(
        qualify_actual_file(
            &source,
            &changed,
            1e-7,
            1e-6,
            QualificationOptions {
                additional_resident_bytes: changed.capacity() - changed.len() + case_bytes,
                ..Default::default()
            },
            None
        ),
        Err(Error::Unsupported)
    ));
}

fn independent_spline_chain(source: &occupied::OccupiedHead, decoded: &reader::ReadHead) {
    use crate::step_read::literal_member_binding::{RetainedRecords, Token, token};
    let mapping = decoded
        .parameter_maps()
        .iter()
        .find(|m| matches!(m.mapping(), reader::ParameterMapping::Spline { spans: 47 }))
        .unwrap();
    let reader::Target::Curve(index) = mapping.target() else {
        panic!("spline curve")
    };
    let shell::Curve::Original(Curve::Spline { spans, poles }) =
        &decoded.material().shell().input().curves[index]
    else {
        panic!("actual spline")
    };
    let component = decoded
        .records
        .component(
            decoded
                .records
                .record_index(mapping.actual_entity())
                .unwrap(),
            0,
        )
        .unwrap();
    let mut rest = component;
    let mut numeric = Vec::new();
    while let Some(t) = token(&mut rest) {
        if let Token::Literal(t) = t
            && t.as_bytes()
                .first()
                .is_some_and(|b| matches!(b, b'+' | b'-' | b'0'..=b'9'))
        {
            numeric.push(t.parse::<f64>().unwrap());
        }
    }
    assert_eq!(numeric.len(), 1 + 2 * (spans + 1));
    let mult = &numeric[1..spans + 2];
    let knots = &numeric[spans + 2..];
    let expanded = knots
        .iter()
        .zip(mult)
        .flat_map(|(&k, &m)| std::iter::repeat_n(k, m as usize))
        .collect::<Vec<_>>();
    let u = 0.375;
    let actual_parameter = u * (*spans as f64);
    let derivative = poles
        .windows(2)
        .enumerate()
        .map(|(i, p)| {
            let width = expanded[i + 10] - expanded[i + 1];
            std::array::from_fn(|k| 9. * (p[1][k] - p[0][k]) / width)
        })
        .collect::<Vec<[f64; 3]>>();
    let step = de_boor(
        &derivative,
        &expanded[1..expanded.len() - 1],
        8,
        actual_parameter,
    );
    let original = source
        .shell()
        .input()
        .curves
        .iter()
        .find_map(|c| match c {
            shell::Curve::Original(Curve::Spline { spans: n, poles: p })
                if n == spans && p == poles =>
            {
                Some(p)
            }
            _ => None,
        })
        .unwrap();
    let cell = actual_parameter.floor() as usize;
    let t = actual_parameter - cell as f64;
    let controls = &original[cell * 9..=cell * 9 + 9];
    let expected = std::array::from_fn::<_, 3, _>(|k| {
        (0..9)
            .map(|i| {
                let binomial = (0..i).fold(1., |v, j| v * (8 - j) as f64 / (j + 1) as f64);
                binomial
                    * t.powi(i as i32)
                    * (1. - t).powi((8 - i) as i32)
                    * 9.
                    * (controls[i + 1][k] - controls[i][k])
                    * (*spans as f64)
            })
            .sum::<f64>()
    });
    let pulled = mapping
        .reference_derivative_bounds(step.map(|x| [x.next_down(), x.next_up()]), 0)
        .unwrap();
    for k in 0..3 {
        assert!((0.5 * (pulled[k][0] + pulled[k][1]) - expected[k]).abs() < 1e-9);
    }
    assert!((0..3).any(|k| (step[k] - expected[k]).abs() > 0.1));
}
fn de_boor(points: &[[f64; 3]], knots: &[f64], degree: usize, u: f64) -> [f64; 3] {
    let span = (degree..points.len())
        .find(|&i| knots[i] <= u && u < knots[i + 1])
        .unwrap();
    let mut d = (0..=degree)
        .map(|j| points[span - degree + j])
        .collect::<Vec<_>>();
    for r in 1..=degree {
        for j in (r..=degree).rev() {
            let i = span - degree + j;
            let a = (u - knots[i]) / (knots[i + degree - r + 1] - knots[i]);
            d[j] = std::array::from_fn(|k| (1. - a) * d[j - 1][k] + a * d[j][k]);
        }
    }
    d[degree]
}

#[test]
fn actual_file_qualification_exact_frontiers_cancel_and_same_source_retry() {
    use crate::axial_pipe_thread::Error;
    let source = engineering();
    let mut actual = propose(&source, Options::default(), None)
        .unwrap()
        .into_string();
    actual.shrink_to_fit();
    public_option_contracts(&source, &actual);
    let calls = std::cell::Cell::new(0usize);
    let count = || {
        calls.set(calls.get() + 1);
        false
    };
    let baseline = qualify_actual_file(
        &source,
        &actual,
        1e-7,
        1e-6,
        QualificationOptions::default(),
        Some(&count),
    )
    .unwrap();
    let receipt = baseline.receipt();
    let original_records = receipt.matched_original_records;
    drop(baseline);
    assert!(matches!(
        propose(
            &source,
            Options {
                maximum_work: receipt.proposal_work_used - 1,
                additional_resident_bytes: actual.capacity(),
                ..Default::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
    let exact_proposal = propose(
        &source,
        Options {
            maximum_work: receipt.proposal_work_used,
            additional_resident_bytes: actual.capacity(),
            ..Default::default()
        },
        None,
    )
    .unwrap();
    assert_eq!(
        exact_proposal.receipt().work_used,
        receipt.proposal_work_used
    );
    assert!(
        exact_proposal.as_str() == actual,
        "late writer refusal must leave a deterministic same-source retry"
    );
    drop(exact_proposal);
    // Independently count the two genuine public phases to target cancellation
    // inside the subsequent dense record binder, without an injected error or
    // test-only bypass in production. Residency changes do not change visits.
    let phase_calls = std::cell::Cell::new(0usize);
    let phase_count = || {
        phase_calls.set(phase_calls.get() + 1);
        false
    };
    drop(
        propose(
            &source,
            Options {
                additional_resident_bytes: actual.capacity(),
                ..Default::default()
            },
            Some(&phase_count),
        )
        .unwrap(),
    );
    let proposal_calls = phase_calls.get();
    assert!(proposal_calls > 2);
    let visited = std::cell::Cell::new(0usize);
    let late_face_cancel = || {
        visited.set(visited.get() + 1);
        visited.get() == proposal_calls - 2
    };
    // The final physical face precedes CLOSED_SHELL and MANIFOLD_SOLID_BREP.
    // This requests real public cancellation, not an injected serializer error.
    assert!(matches!(
        propose(
            &source,
            Options {
                maximum_work: receipt.proposal_work_used,
                additional_resident_bytes: actual.capacity(),
                ..Default::default()
            },
            Some(&late_face_cancel)
        ),
        Err(Error::Cancelled)
    ));
    assert_eq!(visited.get(), proposal_calls - 2);
    let after_cancel = propose(
        &source,
        Options {
            maximum_work: receipt.proposal_work_used,
            additional_resident_bytes: actual.capacity(),
            ..Default::default()
        },
        None,
    )
    .unwrap();
    assert!(after_cancel.as_str() == actual);
    drop(after_cancel);
    phase_calls.set(0);
    drop(
        reader::read(
            &actual,
            1e-7,
            1e-6,
            reader::Options {
                additional_resident_bytes: source.receipt().retained_result_bytes_upper_bound
                    + actual.capacity()
                    - actual.len(),
                ..Default::default()
            },
            Some(&phase_count),
        )
        .unwrap(),
    );
    let reader_calls = phase_calls.get();
    let binding_start = proposal_calls + reader_calls + 1;
    assert!(binding_start < calls.get());
    for remaining in [0, 1] {
        assert!(matches!(
            qualify_actual_file(
                &source,
                &actual,
                1e-7,
                1e-6,
                QualificationOptions {
                    maximum_work: receipt.proposal_work_used
                        + receipt.independent_read_work_reserved
                        + remaining,
                    ..Default::default()
                },
                None
            ),
            Err(Error::Limit)
        ));
    }
    let exact = QualificationOptions {
        maximum_work: receipt.aggregate_work_reserved,
        maximum_bytes: receipt
            .proposal_peak_bytes_reserved
            .max(receipt.read_peak_bytes_reserved)
            .max(receipt.correspondence_peak_bytes_reserved),
        additional_resident_bytes: actual.capacity() - actual.len(),
    };
    println!(
        "EXACT_QUALIFICATION_FRONTIER {exact:?} interrupt_checks{}",
        calls.get()
    );
    for options in [
        QualificationOptions {
            maximum_work: exact.maximum_work - 1,
            ..exact
        },
        QualificationOptions {
            maximum_bytes: exact.maximum_bytes - 1,
            ..exact
        },
    ] {
        assert!(matches!(
            qualify_actual_file(&source, &actual, 1e-7, 1e-6, options, None),
            Err(Error::Limit)
        ));
    }
    // These checkpoints cross proposal, cold parsing, independent material and
    // correspondence phases. No partial QualifiedStep is returned on refusal.
    for stop in [1, calls.get() / 2, binding_start, calls.get() - 1] {
        let visited = std::cell::Cell::new(0usize);
        let cancel = || {
            let n = visited.get() + 1;
            visited.set(n);
            n >= stop
        };
        assert!(matches!(
            qualify_actual_file(&source, &actual, 1e-7, 1e-6, exact, Some(&cancel)),
            Err(Error::Cancelled)
        ));
        assert_eq!(visited.get(), stop);
    }
    let retry = qualify_actual_file(&source, &actual, 1e-7, 1e-6, exact, None).unwrap();
    assert_eq!(
        retry.receipt().aggregate_work_reserved,
        receipt.aggregate_work_reserved
    );
    assert_eq!(retry.receipt().matched_original_records, original_records);
    assert_eq!(retry.as_str(), actual);
    assert_eq!(retry.material().material().head_height_mm(), [6., 9.]);
    assert_eq!(source.head_height_mm(), [6., 9.]);
}

fn independent_mapped_geometry_consumers(decoded: &reader::ReadHead) {
    use crate::tapered_thread::interval::I;
    use num_rational::BigRational;
    let q = |v| BigRational::from_float(v).unwrap();
    let contains = |bounds: [f64; 2], exact: BigRational| {
        assert!(q(bounds[0]) <= exact && exact <= q(bounds[1]));
    };
    let input = decoded.material().shell().input();
    let mut families = [false; 4];
    for map in decoded.parameter_maps() {
        assert!(matches!(
            map.reference_derivative_bounds([[0., 0.]; 3], 2),
            Err(crate::axial_pipe_thread::Error::Invalid)
        ));
        assert!(
            map.reference_derivative_bounds([[f64::NAN, f64::NAN]; 3], 0)
                .is_err()
        );
        assert!(map.reference_derivative_bounds([[1., -1.]; 3], 0).is_err());
        match map.mapping() {
            reader::ParameterMapping::Tensor { spans } => {
                families[0] = true;
                let reader::Target::Surface(index) = map.target() else {
                    unreachable!()
                };
                let Surface::Tensor(t) = &input.surfaces[index] else {
                    unreachable!()
                };
                let first = t.poles()[0];
                for axis in 0..2 {
                    // Integer-knot de Boor endpoint derivative has local span
                    // width1. Reference uniform knots have width1/n instead.
                    let next = t.poles()[if axis == 0 { 9 * spans[1] + 1 } else { 1 }];
                    let raw = std::array::from_fn(|k| {
                        (I::point(next[k]) - I::point(first[k])) * I::point(9.)
                    });
                    let mapped = map
                        .reference_derivative_bounds(raw.map(|v| [v.lo, v.hi]), axis)
                        .unwrap();
                    for k in 0..3 {
                        contains(
                            mapped[k],
                            (q(next[k]) - q(first[k]))
                                * BigRational::from_integer((9 * spans[axis]).into()),
                        );
                    }
                }
            }
            reader::ParameterMapping::Line { .. } => {
                families[1] = true;
                let reader::Target::Curve(index) = map.target() else {
                    unreachable!()
                };
                let shell::Curve::Line(ends) = input.curves[index] else {
                    unreachable!()
                };
                let a = input.vertices[ends[0]];
                let z = input.vertices[ends[1]];
                let delta = std::array::from_fn::<_, 3, _>(|k| I::point(z[k]) - I::point(a[k]));
                let length = delta
                    .into_iter()
                    .fold(I::point(0.), |s, d| s + d.square())
                    .sqrt();
                let unit = delta.map(|d| d / length);
                let mapped = map
                    .reference_derivative_bounds(unit.map(|v| [v.lo, v.hi]), 0)
                    .unwrap();
                for k in 0..3 {
                    contains(mapped[k], q(z[k]) - q(a[k]));
                }
                assert!(map.reference_derivative_bounds([[0., 0.]; 3], 1).is_err());
            }
            reader::ParameterMapping::AngularIdentity => {
                families[2] = true;
                let reader::Target::Curve(index) = map.target() else {
                    unreachable!()
                };
                let shell::Curve::Original(Curve::Circle { axis, radius_mm }) =
                    &input.curves[index]
                else {
                    unreachable!()
                };
                let tangent = [0., axis.normal[2] * radius_mm, 0.];
                let mapped = map
                    .reference_derivative_bounds(tangent.map(|v| [v, v]), 0)
                    .unwrap();
                for k in 0..3 {
                    contains(mapped[k], q(tangent[k]));
                }
                assert!(map.reference_derivative_bounds([[0., 0.]; 3], 1).is_err());
            }
            reader::ParameterMapping::Cone { .. } => {
                families[3] = true;
                let reader::Target::Surface(index) = map.target() else {
                    unreachable!()
                };
                let Surface::Cone {
                    axis, radius_mm, ..
                } = &input.surfaces[index]
                else {
                    unreachable!()
                };
                let tangent = [0., axis.normal[2] * radius_mm, 0.];
                let mapped = map
                    .reference_derivative_bounds(tangent.map(|v| [v, v]), 0)
                    .unwrap();
                for k in 0..3 {
                    contains(mapped[k], q(tangent[k]));
                }
            }
            reader::ParameterMapping::Spline { .. } => {}
        }
    }
    assert_eq!(families, [true; 4]);
}

fn public_option_contracts(source: &occupied::OccupiedHead, actual: &str) {
    use crate::axial_pipe_thread::Error;
    for options in [
        Options {
            maximum_work: 0,
            ..Default::default()
        },
        Options {
            maximum_work: 32_000_001,
            ..Default::default()
        },
        Options {
            maximum_bytes: 0,
            ..Default::default()
        },
        Options {
            maximum_bytes: 64 * 1024 * 1024 + 1,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            propose(source, options, None),
            Err(Error::Invalid)
        ));
    }
    for options in [
        QualificationOptions {
            maximum_work: 0,
            ..Default::default()
        },
        QualificationOptions {
            maximum_work: 272_000_001,
            ..Default::default()
        },
        QualificationOptions {
            maximum_bytes: 0,
            ..Default::default()
        },
        QualificationOptions {
            maximum_bytes: 64 * 1024 * 1024 + 1,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            qualify_actual_file(source, actual, 1e-7, 1e-6, options, None),
            Err(Error::Invalid)
        ));
    }
    for options in [
        reader::Options {
            maximum_work: 0,
            ..Default::default()
        },
        reader::Options {
            maximum_work: 224_000_001,
            ..Default::default()
        },
        reader::Options {
            maximum_bytes: 0,
            ..Default::default()
        },
        reader::Options {
            maximum_bytes: 64 * 1024 * 1024 + 1,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            reader::read(actual, 1e-7, 1e-6, options, None),
            Err(Error::Invalid)
        ));
    }
    assert!(matches!(
        reader::read(
            actual,
            1e-7,
            1e-6,
            reader::Options::default(),
            Some(&|| true)
        ),
        Err(Error::Cancelled)
    ));
    assert!(matches!(
        propose(
            source,
            Options {
                maximum_bytes: 1,
                ..Default::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        propose(
            source,
            Options {
                additional_resident_bytes: usize::MAX,
                ..Default::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        reader::read(
            actual,
            1e-7,
            1e-6,
            reader::Options {
                additional_resident_bytes: usize::MAX,
                ..Default::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        qualify_actual_file(
            source,
            actual,
            1e-7,
            1e-6,
            QualificationOptions {
                additional_resident_bytes: usize::MAX,
                ..Default::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
}

#[test]
fn earned_finite_line_refuses_rounding_an_unrepresentable_endpoint_difference() {
    use super::test_fixture::engineering_dimensions;
    use num_rational::BigRational;
    let source = engineering_dimensions(6.1, 100.);
    let [lo, hi] = source.head_height_mm();
    assert_eq!([lo, hi], [6.1, 106.1]);
    let q = |v| BigRational::from_float(v).unwrap();
    assert_ne!(q(hi) - q(lo), q(hi - lo));
    assert!(matches!(
        propose(&source, Options::default(), None),
        Err(crate::axial_pipe_thread::Error::Unsupported)
    ));
    assert_eq!(source.head_height_mm(), [lo, hi]);
    assert_eq!(
        source.capped_thread().classify([0., 0., 3.], None).unwrap(),
        crate::radial_tensor::Location::Inside
    );
    drop(source);
    let healthy = engineering();
    let retry = propose(&healthy, Options::default(), None).unwrap();
    assert!(retry.as_str().contains("MANIFOLD_SOLID_BREP"));
}
