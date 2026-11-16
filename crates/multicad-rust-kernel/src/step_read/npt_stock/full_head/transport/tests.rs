use super::*;
use num_rational::BigRational as Q;
use num_traits::One;
fn q(x: f64) -> Q {
    Q::from_float(x).unwrap()
}
fn contains(b: [f64; 2], q: &Q) {
    assert!(
        self::q(b[0]) <= *q && *q <= self::q(b[1]),
        "{b:?} excludes {q}"
    );
}
fn candidate(quarter: bool) -> String {
    let input = super::super::super::tests::fixture(quarter);
    let doc = CompoundDocument::parse(&input, Default::default(), None).unwrap();
    let parsed = Parsed::parse(
        &doc,
        doc.represented_members()[0].item_entity(),
        Default::default(),
        None,
    )
    .unwrap();
    let cap = parsed
        .faces()
        .iter()
        .filter_map(|f| match f.support() {
            Support::Plane { height } => {
                Some((parsed.points()[height].coordinates()[2].value(), f.entity()))
            }
            _ => None,
        })
        .max_by(|a, b| a.0.total_cmp(&b.0))
        .unwrap()
        .1;
    let head = HeadInput::regular_hex(&parsed, cap, 18., 3.).unwrap();
    let proof = head::certify(&parsed, head, 1e-7, 1e-6, Default::default(), None).unwrap();
    crate::step::npt_stock::full_head::propose(&proof, Default::default(), None)
        .unwrap()
        .as_str()
        .to_owned()
}
fn change(input: &str, id: u32, body: &str) -> String {
    let start = input.find(&format!("#{id}=")).unwrap();
    let end = start + input[start..].find(';').unwrap() + 1;
    format!("{}#{id}={body};{}", &input[..start], &input[end..])
}
pub(super) fn mapped(quarter: bool, pose: usize) -> String {
    let text = candidate(quarter);
    let doc = CompoundDocument::parse(&text, Default::default(), None).unwrap();
    let (&top, parts) = doc
        .records()
        .iter()
        .find(|(_, p)| matches!(p.as_slice(),[e] if e.name=="ADVANCED_BREP_SHAPE_REPRESENTATION"))
        .unwrap();
    let solid = doc.represented_members()[0].item_entity();
    let [_, Value::List(items), Value::Reference(context)] = parts[0].args.as_slice() else {
        panic!()
    };
    let axis = items
        .iter()
        .find_map(|v| match v {
            Value::Reference(id) if *id != solid => Some(*id),
            _ => None,
        })
        .unwrap();
    let first = *doc.records().keys().next_back().unwrap() + 1;
    let [rep, map, item, op, origin, x, y, z] =
        std::array::from_fn::<_, 8, _>(|i| first + i as u32);
    let body = format!("ADVANCED_BREP_SHAPE_REPRESENTATION('',(#{item},#{axis}),#{context})");
    let mut out = change(&text, top, &body);
    let directions = match pose {
        0 => [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        1 => [[1, 0, 0], [0, -1, 0], [0, 0, -1]],
        _ => [[3, 4, 0], [4, -3, 0], [0, 0, -5]],
    };
    let (origin_text, scale) = if pose == 0 {
        ("0.,0.,0.", 1)
    } else {
        ("2.,-3.,7.", 2)
    };
    let mut extra = format!(
        "#{rep}=SHAPE_REPRESENTATION('',(#{solid},#{axis}),#{context});\n#{map}=REPRESENTATION_MAP(#{axis},#{rep});\n#{item}=MAPPED_ITEM('',#{map},#{op});\n#{op}=(CARTESIAN_TRANSFORMATION_OPERATOR(#{x},#{y},#{origin},{scale}.)CARTESIAN_TRANSFORMATION_OPERATOR_3D(#{z})FUNCTIONALLY_DEFINED_TRANSFORMATION('','')GEOMETRIC_REPRESENTATION_ITEM()REPRESENTATION_ITEM(''));\n#{origin}=CARTESIAN_POINT('',({origin_text}));\n"
    );
    for (id, p) in [x, y, z].into_iter().zip(directions) {
        writeln!(
            extra,
            "#{id}=DIRECTION('',({}.,{}.,{}.));",
            p[0], p[1], p[2]
        )
        .unwrap();
    }
    let end = out.rfind("ENDSEC;").unwrap();
    out.insert_str(end, &extra);
    out
}
#[test]
fn actual_mapped_full_head_retains_physical_geometry_and_earns_file_operator_after_all_arenas_drop()
{
    for quarter in [false, true] {
        for pose in 0..3 {
            let owned = {
                let input: Arc<str> = Arc::from(mapped(quarter, pose));
                MappedActualHead::read(input, Default::default(), None).unwrap()
            };
            let original =
                CompoundDocument::parse(owned.actual_file(), Default::default(), None).unwrap();
            let local = CompoundDocument::parse(
                owned.local_geometry().actual_file(),
                Default::default(),
                None,
            )
            .unwrap();
            // Every actual physical component, including scalar spellings and
            // UV ownership, survives with the SAME actual entity ID.
            for (&id, parts) in local.records() {
                if parts.iter().any(|p| {
                    matches!(
                        p.name,
                        "CARTESIAN_POINT"
                            | "DIRECTION"
                            | "B_SPLINE_SURFACE_WITH_KNOTS"
                            | "B_SPLINE_CURVE_WITH_KNOTS"
                            | "PCURVE"
                            | "ADVANCED_FACE"
                            | "EDGE_CURVE"
                            | "PLANE"
                            | "LINE"
                            | "VECTOR"
                    )
                }) {
                    assert_eq!(parts.len(), original.records()[&id].len());
                    for (a, z) in parts.iter().zip(&original.records()[&id]) {
                        assert_eq!(a.source, z.source);
                    }
                }
            }
            let member = original.represented_members()[0];
            assert_eq!(owned.member_entity(), member.item_entity());
            assert_eq!(
                owned.operator_entity(),
                member.mapping_link().unwrap().operator_entity()
            );
            drop(original);
            drop(local);
            let frame = owned.actual_frame();
            let origin = if pose == 0 {
                [0., 0., 0.]
            } else {
                [2., -3., 7.]
            };
            let scale = if pose == 0 { 1. } else { 2. };
            for (b, v) in frame.origin_mm().into_iter().zip(origin) {
                contains(b, &q(v));
            }
            contains(frame.scale_bounds(), &q(scale));
            let fifth = Q::one() / Q::from_integer(5.into());
            let columns = match pose {
                0 => [
                    [q(1.), q(0.), q(0.)],
                    [q(0.), q(1.), q(0.)],
                    [q(0.), q(0.), q(1.)],
                ],
                1 => [
                    [q(1.), q(0.), q(0.)],
                    [q(0.), q(-1.), q(0.)],
                    [q(0.), q(0.), q(-1.)],
                ],
                _ => [
                    [&fifth * q(3.), &fifth * q(4.), q(0.)],
                    [&fifth * q(4.), &fifth * q(-3.), q(0.)],
                    [q(0.), q(0.), q(-1.)],
                ],
            };
            for (bounds, col) in frame.axes().into_iter().zip(&columns) {
                for (b, v) in bounds.into_iter().zip(col) {
                    contains(b, v);
                }
            }
            let bounds = owned.world_sewn_enclosure_mm();
            for vertex in &owned.local_geometry().head().vertices {
                for i in 0..3 {
                    let value = (0..3).fold(q(origin[i]), |s, j| {
                        s + q(scale) * &columns[j][i] * q(vertex[j])
                    });
                    contains(bounds[i], &value);
                }
            }
            // Independent exact rational extrema of the COMPLETE retained
            // enclosure, not only the finite head vertices.
            let local_bounds = owned.local_geometry().local_sewn_enclosure_mm();
            for mask in 0..8 {
                let point = std::array::from_fn::<_, 3, _>(|j| q(local_bounds[j][(mask >> j) & 1]));
                for i in 0..3 {
                    let value = (0..3).fold(q(origin[i]), |sum, j| {
                        sum + q(scale) * &columns[j][i] * &point[j]
                    });
                    contains(bounds[i], &value);
                }
            }
            let volume = owned.local_geometry().local_sewn_volume_mm3();
            for v in volume {
                contains(owned.world_sewn_volume_mm3(), &(q(v) * q(scale).pow(3)));
            }
            eprintln!(
                "MAPPED_HEAD quarter={quarter} pose={pose} {:?}",
                owned.receipt()
            );
        }
    }
}
#[test]
fn actual_mapped_head_refuses_orphan_geometry_reflection_and_unknown_operator_components() {
    let text = mapped(false, 2);
    let doc = CompoundDocument::parse(&text, Default::default(), None).unwrap();
    let link = doc.represented_members()[0].mapping_link().unwrap();
    let op = link.operator_entity();
    let parts = &doc.records()[&op];
    let combined = parts.iter().map(|p| p.source.unwrap()).collect::<String>();
    let unknown = change(&text, op, &format!("({combined}EXTRA_OWNER(#{op}))"));
    assert!(matches!(
        MappedActualHead::read(Arc::from(unknown), Default::default(), None),
        Err(Error::Unsupported)
    ));
    let y = parts
        .iter()
        .find(|e| e.name == "CARTESIAN_TRANSFORMATION_OPERATOR")
        .unwrap();
    let Value::Reference(y) = y.args[1] else {
        panic!()
    };
    let reflected = change(&text, y, "DIRECTION('',(-4.,3.,0.))");
    assert!(matches!(
        MappedActualHead::read(Arc::from(reflected), Default::default(), None),
        Err(Error::Unsupported)
    ));
    let face = doc
        .records()
        .values()
        .find(|p| matches!(p.as_slice(),[e] if e.name=="ADVANCED_FACE"))
        .unwrap()[0]
        .source
        .unwrap();
    let fresh = *doc.records().keys().next_back().unwrap() + 1;
    let mut orphan = text.clone();
    let end = orphan.rfind("ENDSEC;").unwrap();
    orphan.insert_str(end, &format!("#{fresh}={face};\n"));
    assert!(matches!(
        MappedActualHead::read(Arc::from(orphan), Default::default(), None),
        Err(Error::Unsupported)
    ));
    let huge_operator = combined.replacen(
        &format!(
            "#{},2.)",
            match parts
                .iter()
                .find(|e| e.name == "CARTESIAN_TRANSFORMATION_OPERATOR")
                .unwrap()
                .args[2]
            {
                Value::Reference(v) => v,
                _ => panic!(),
            }
        ),
        &format!(
            "#{},1.E300)",
            match parts
                .iter()
                .find(|e| e.name == "CARTESIAN_TRANSFORMATION_OPERATOR")
                .unwrap()
                .args[2]
            {
                Value::Reference(v) => v,
                _ => panic!(),
            }
        ),
        1,
    );
    assert_ne!(huge_operator, combined);
    let overflow = change(&text, op, &format!("({huge_operator})"));
    assert!(matches!(
        MappedActualHead::read(Arc::from(overflow), Default::default(), None),
        Err(Error::Unsupported)
    ));
    let good = MappedActualHead::read(Arc::from(text), Default::default(), None).unwrap();
    assert!(good.world_sewn_volume_mm3()[0] > 0.);
}
#[test]
fn mapped_head_work_bytes_and_cancel_frontiers_are_atomic() {
    use std::cell::Cell;
    let input: Arc<str> = Arc::from(mapped(false, 2));
    let good = MappedActualHead::read(Arc::clone(&input), Default::default(), None).unwrap();
    let r = good.receipt();
    for (work, bytes, ok) in [
        (r.aggregate_work - 1, 64 * 1024 * 1024, false),
        (r.aggregate_work, 64 * 1024 * 1024, true),
        (272_000_000, r.accounted_peak_bytes - 1, false),
        (272_000_000, r.accounted_peak_bytes, true),
    ] {
        let value = MappedActualHead::read(
            Arc::clone(&input),
            OwnedHeadOptions {
                resources: Options {
                    maximum_work: work,
                    maximum_bytes: bytes,
                    ..Default::default()
                },
                ..Default::default()
            },
            None,
        );
        assert_eq!(value.is_ok(), ok);
        if !ok {
            assert!(matches!(value, Err(Error::Limit)));
        }
    }
    let count = Cell::new(0);
    let observe = || {
        count.set(count.get() + 1);
        false
    };
    MappedActualHead::read(Arc::clone(&input), Default::default(), Some(&observe)).unwrap();
    for at in [1, count.get() / 2, count.get()] {
        let seen = Cell::new(0);
        let cancel = || {
            seen.set(seen.get() + 1);
            seen.get() >= at
        };
        assert!(matches!(
            MappedActualHead::read(Arc::clone(&input), Default::default(), Some(&cancel)),
            Err(Error::Cancelled)
        ));
    }
    let retry = MappedActualHead::read(input, Default::default(), None).unwrap();
    assert_eq!(retry.receipt().aggregate_work, r.aggregate_work);
    assert_eq!(retry.world_sewn_volume_mm3(), good.world_sewn_volume_mm3());
}
