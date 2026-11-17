use super::*;
use crate::{
    Rational,
    rounded_thread::{Spec, Thread},
    step::{Writer, write_rounded_thread_blind_bore},
};
fn source(left_hand: bool) -> BlindBore {
    let p = 1.27;
    let a = 27.5_f64.to_radians();
    let r = p * a.cos() / (12. * (1. - a.sin()));
    BlindBore::new(
        Thread::new(Spec {
            major_radius_mm: 3.175,
            pitch_mm: p,
            length_mm: 3. * p,
            depth_mm: p / (3. * a.tan()),
            crest_radius_mm: r,
            root_radius_mm: r,
            flank_angle_deg: 27.5,
            starts: 1,
            left_hand,
        })
        .unwrap(),
        Rational::new(127, 20).unwrap(),
        Rational::new(127, 25).unwrap(),
    )
    .unwrap()
}
fn value(v: &Value<'_>) -> String {
    match v {
        Value::String(s) => format!("'{s}'"),
        Value::Enum(s) => format!(".{s}."),
        Value::Reference(i) => format!("#{i}"),
        Value::Number(n) => format!("{n:.16E}"),
        Value::List(a) => format!("({})", a.iter().map(value).collect::<Vec<_>>().join(",")),
        Value::Omitted => "$".into(),
        Value::Derived => "*".into(),
        Value::Typed(e) => entity(e),
    }
}
fn entity(e: &Entity<'_>) -> String {
    format!(
        "{}({})",
        e.name,
        e.args.iter().map(value).collect::<Vec<_>>().join(",")
    )
}
fn rewrite<'a>(file: &'a str, edit: impl FnOnce(&mut BTreeMap<u32, Vec<Entity<'a>>>)) -> String {
    let mut entities = part21::parse(file).unwrap();
    edit(&mut entities);
    let data = file.find("DATA;").unwrap() + 5;
    let end = file.rfind("ENDSEC;").unwrap();
    let body = entities
        .iter()
        .map(|(id, parts)| {
            let s = if parts.len() == 1 {
                entity(&parts[0])
            } else {
                format!(
                    "({})",
                    parts.iter().map(entity).collect::<Vec<_>>().join(" ")
                )
            };
            format!("#{id}={s};")
        })
        .collect::<Vec<_>>()
        .join("\n");
    format!("{}\n{body}\n{}", &file[..data], &file[end..])
}
fn nth<'m, 's>(
    map: &'m mut BTreeMap<u32, Vec<Entity<'s>>>,
    name: &str,
    n: usize,
) -> &'m mut Entity<'s> {
    map.values_mut()
        .filter_map(|p| p.first_mut())
        .filter(|e| e.name == name)
        .nth(n)
        .unwrap()
}

#[test]
fn actual_nine_face_blind_files_keep_original_rounded_walls_and_separate_stock_top() {
    for hand in [false, true] {
        let b = source(hand);
        for f in [
            [[0.; 3], [1., 0., 0.], [0., 1., 0.]],
            [[17., -8., 13.], [0.8, 0.6, 0.], [-0.36, 0.48, 0.8]],
        ] {
            let file = write_rounded_thread_blind_bore(&b, f[0], f[1], f[2], 1e-7).unwrap();
            let volume =
                verify_rounded_thread_blind_bore(&file, &b, f[0], f[1], f[2], 1e-7).unwrap();
            let expected = std::f64::consts::PI
                * b.stock_radius_mm().as_f64().powi(2)
                * b.stock_height_mm().as_f64()
                - b.cutter().volume_mm3();
            assert!((volume - expected).abs() < 1e-6, "{volume} vs{expected}");
            let ast = part21::parse(&file).unwrap();
            assert_eq!(
                ast.values()
                    .flatten()
                    .filter(|e| e.name == "ADVANCED_FACE")
                    .count(),
                9
            );
            assert_eq!(
                ast.values().flatten().filter(|e| e.name == "PLANE").count(),
                3
            );
            assert_eq!(
                ast.values()
                    .flatten()
                    .filter(|e| e.name == "CYLINDRICAL_SURFACE")
                    .count(),
                1
            );
            assert_eq!(
                ast.values()
                    .flatten()
                    .filter(|e| e.name == "B_SPLINE_SURFACE_WITH_KNOTS")
                    .count(),
                5
            );
        }
    }
}

#[test]
fn physical_floor_top_wall_and_original_tensor_mutations_refuse_with_original_retry() {
    let b = source(false);
    let f = [[0.; 3], [1., 0., 0.], [0., 1., 0.]];
    let good = write_rounded_thread_blind_bore(&b, f[0], f[1], f[2], 1e-7).unwrap();
    for case in 0..6 {
        let bad = rewrite(&good, |map| match case {
            0 => nth(map, "ADVANCED_FACE", 6).args[3] = Value::Enum("F"),
            1 => {
                let top = refs(&nth(map, "ADVANCED_FACE", 7).args[1]).unwrap();
                nth(map, "ADVANCED_FACE", 6).args[1] =
                    Value::List(top.into_iter().map(Value::Reference).collect());
            }
            2 => nth(map, "CYLINDRICAL_SURFACE", 0).args[2] = Value::Number(6.5),
            3 => {
                let Value::List(v) = &mut nth(map, "CLOSED_SHELL", 0).args[1] else {
                    panic!()
                };
                let Value::Reference(first) = v[0] else {
                    panic!()
                };
                v.push(Value::Reference(first));
            }
            4 => {
                let Value::Reference(surface) = nth(map, "ADVANCED_FACE", 7).args[2] else {
                    panic!()
                };
                let Value::Reference(axis) = map[&surface][0].args[1] else {
                    panic!()
                };
                let Value::Reference(point) = map[&axis][0].args[1] else {
                    panic!()
                };
                let Value::List(coords) = &mut map.get_mut(&point).unwrap()[0].args[1] else {
                    panic!()
                };
                coords[2] = Value::Number(5.5);
            }
            5 => {
                let Value::List(rows) = &nth(map, "B_SPLINE_SURFACE_WITH_KNOTS", 0).args[3] else {
                    panic!()
                };
                let Value::List(row) = &rows[rows.len() / 2] else {
                    panic!()
                };
                let Value::Reference(point) = row[row.len() / 2] else {
                    panic!()
                };
                let Value::List(coords) = &mut map.get_mut(&point).unwrap()[0].args[1] else {
                    panic!()
                };
                let Value::Number(x) = coords[0] else {
                    panic!()
                };
                coords[0] = Value::Number(x + 0.001);
            }
            _ => unreachable!(),
        });
        assert_ne!(bad, good);
        // The independently reconstructed repaired interpretation must reject
        // the same actual broken floor/top/stock/tensor graphs.
        assert!(
            read_repaired_rounded_blind_bore(&bad, 1e-7).is_err(),
            "repaired case {case}"
        );
        assert!(read_repaired_rounded_blind_bore(&good, 1e-7).is_ok());

        assert!(
            verify_rounded_thread_blind_bore(&bad, &b, f[0], f[1], f[2], 1e-7).is_err(),
            "case{case}"
        );
    }
    assert!(verify_rounded_thread_blind_bore(&good, &b, f[0], f[1], f[2], 1e-7).is_ok());
}

#[test]
fn failed_blind_append_preserves_existing_physical_member() {
    let b = source(true);
    let f = [[0.; 3], [1., 0., 0.], [0., 1., 0.]];
    let mut writer = Writer::new();
    writer
        .append_rounded_thread_blind_bore(&b, f[0], f[1], f[2], 1e-7)
        .unwrap();
    assert!(
        writer
            .append_rounded_thread_blind_bore(&b, f[0], f[1], f[2], 0.)
            .is_err()
    );
    assert!(
        writer
            .append_rounded_thread_blind_bore(&b, f[0], [2., 0., 0.], f[2], 1e-7)
            .is_err()
    );
    let file = writer.finish().unwrap();
    assert!(verify_rounded_thread_blind_bore(&file, &b, f[0], f[1], f[2], 1e-7).is_ok());
    assert_eq!(b.stock_height_mm(), Rational::new(127, 25).unwrap());
}

#[test]
fn large_stock_decimal_aliases_cannot_earn_sub_ulp_correspondence() {
    let small = source(false);
    let b = BlindBore::new(
        *small.cutter(),
        Rational::new(8_000_000, 1).unwrap(),
        small.stock_height_mm(),
    )
    .unwrap();
    let f = [[0.; 3], [1., 0., 0.], [0., 1., 0.]];
    let tolerance = 1e-10;
    // Genuine public emission, followed by a coherent actual-file change of
    // only the outer stock radius and its cardinal seam coordinates.
    let ordinary = write_rounded_thread_blind_bore(&small, f[0], f[1], f[2], 1e-7).unwrap();
    let stock = format!("{:.16E}", small.stock_radius_mm().as_f64());
    assert!(ordinary.matches(&stock).count() >= 3);
    let needle = format!("{:.16E}", 8_000_000_f64);
    let good = ordinary.replace(&stock, &needle);
    let alias = "8000000.0000000002";
    assert_eq!(alias.parse::<f64>().unwrap(), 8_000_000.);
    assert!(good.matches(&needle).count() >= 3);
    let bad = good.replace(&needle, alias);
    assert_ne!(good, bad);
    // Actual stock circles/cylinder and their cardinal seam vertices change
    // coherently by2e-10, while every neutral-parser cached value is identical.
    assert_eq!(part21::parse(&good).unwrap(), part21::parse(&bad).unwrap());
    if let Ok(path) = std::env::var("MULTICAD_BLIND_ALIAS_DIAGNOSTIC") {
        std::fs::write(format!("{path}-good.step"), &good).unwrap();
        std::fs::write(format!("{path}-alias.step"), &bad).unwrap();
    }
    assert_eq!(
        verify_rounded_thread_blind_bore(&bad, &b, f[0], f[1], f[2], tolerance),
        Err(ReadError::Limit)
    );
    assert_eq!(
        verify_rounded_thread_blind_bore(&good, &b, f[0], f[1], f[2], tolerance),
        Err(ReadError::Limit)
    );
    assert_eq!(
        write_rounded_thread_blind_bore(&b, f[0], f[1], f[2], tolerance),
        Err(crate::step::StepError::Limit)
    );
    assert!(verify_rounded_thread_blind_bore(&ordinary, &small, f[0], f[1], f[2], 1e-7).is_ok());
}

#[test]
fn raw_axis_conditioning_and_actual_near_unit_support_require_earned_transfer() {
    let b = source(false);
    let f = [[0.; 3], [1., 0., 0.], [0., 1., 0.]];
    let good = write_rounded_thread_blind_bore(&b, f[0], f[1], f[2], 1e-7).unwrap();
    let tiny = rewrite(&good, |map| {
        nth(map, "DIRECTION", 0).args[1] = Value::List(vec![
            Value::Number(1e-310),
            Value::Number(0.),
            Value::Number(0.),
        ]);
    });
    assert_eq!(
        verify_rounded_thread_blind_bore(&tiny, &b, f[0], f[1], f[2], 1e-7),
        Err(ReadError::Unsupported)
    );
    let projected = rewrite(&good, |map| {
        let Value::Reference(id) = nth(map, "AXIS2_PLACEMENT_3D", 0).args[3] else {
            panic!()
        };
        map.get_mut(&id).unwrap()[0].args[1] = Value::List(vec![
            Value::Number(1e-9),
            Value::Number(0.),
            Value::Number(1.),
        ]);
    });
    assert_eq!(
        verify_rounded_thread_blind_bore(&projected, &b, f[0], f[1], f[2], 1e-7),
        Err(ReadError::Unsupported)
    );
    // Gram admission alone accepts this source frame, but its finite stock
    // circle transfer exceeds the requested tolerance. Existing complete
    // decoded outer-frame comparison must refuse rather than normalize it.
    assert!(write_rounded_thread_blind_bore(&b, f[0], [1. + 2e-11, 0., 0.], f[2], 1e-10).is_err());
    assert!(verify_rounded_thread_blind_bore(&good, &b, f[0], f[1], f[2], 1e-7).is_ok());
    let unit_with_tiny_component = [1., f64::from_bits(1), 0.];
    let supported =
        write_rounded_thread_blind_bore(&b, f[0], unit_with_tiny_component, f[2], 1e-7).unwrap();
    assert!(
        verify_rounded_thread_blind_bore(
            &supported,
            &b,
            f[0],
            unit_with_tiny_component,
            f[2],
            1e-7,
        )
        .is_ok()
    );
}
