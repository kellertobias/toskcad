use super::*;
use crate::step_read::owned_insert::DocumentOptions;
#[test]
fn actual_decimal_tensor_charts_and_cap_rims_earn_local_projection_without_source_metadata() {
    for quarter in [false, true] {
        let text = super::super::tests::fixture(quarter);
        let document = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
        let parsed = Parsed::parse(
            &document,
            document.represented_members()[0].item_entity(),
            Default::default(),
            None,
        )
        .unwrap();
        let projection = parsed.projection(Default::default(), None).unwrap();
        eprintln!(
            "NPT LITERAL PROJECTION quarter={quarter} work={} nodes={} leaves={} querybytes={}",
            projection.work_used(),
            projection.nodes(),
            projection.leaves(),
            projection.query_budget_bytes()
        );
        assert!(std::ptr::eq(projection.parsed(), &parsed));
        assert!(
            projection
                .cap_face_indices()
                .iter()
                .all(|i| *i < parsed.faces().len())
        );
        assert!(projection.nodes() >= projection.leaves() && projection.leaves() >= 216);
        assert!(projection.work_used() > parsed.work_used());
    }
}
#[test]
fn actual_literal_occupied_volume_and_cap_integrals_use_original_decimal_geometry() {
    for quarter in [false, true] {
        let text = super::super::tests::fixture(quarter);
        let document = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
        let parsed = Parsed::parse(
            &document,
            document.represented_members()[0].item_entity(),
            Default::default(),
            None,
        )
        .unwrap();
        let occupied = parsed
            .projection(Default::default(), None)
            .unwrap()
            .occupied(Default::default(), None)
            .unwrap();
        eprintln!(
            "NPT LITERAL OCCUPIED quarter={quarter} volume={:?} areas={:?} work={} querybytes={}",
            occupied.volume_mm3(),
            occupied.cap_areas_mm2(),
            occupied.work_used(),
            occupied.query_budget_bytes()
        );
        assert!(
            occupied.volume_mm3()[0] > 0.
                && occupied.volume_mm3()[1] - occupied.volume_mm3()[0] < 1e-5
        );
        assert!(
            occupied
                .cap_areas_mm2()
                .iter()
                .all(|area| area[0] > 0. && area[1] - area[0] < 1e-7)
        );
        assert!(std::ptr::eq(occupied.projection().parsed(), &parsed));
        let direct = direct_flux(&parsed, 17);
        let check = direct_flux(&parsed, 20);
        assert!((direct - check).abs() < 1e-10, "{direct} {check}");
        assert!(
            occupied.volume_mm3()[0] <= direct && direct <= occupied.volume_mm3()[1],
            "{direct} {:?}",
            occupied.volume_mm3()
        );
    }
}
fn changed_record(text: &str, id: u32, change: impl FnOnce(&str) -> String) -> String {
    let prefix = format!("#{id}=");
    let record = text.lines().find(|line| line.starts_with(&prefix)).unwrap();
    let changed = change(record);
    assert_ne!(record, changed);
    text.replacen(record, &changed, 1)
}
#[test]
fn actual_face_senses_must_match_the_proved_source_orientation() {
    let text = super::super::tests::fixture(false);
    let document = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
    let parsed = Parsed::parse(
        &document,
        document.represented_members()[0].item_entity(),
        Default::default(),
        None,
    )
    .unwrap();
    for face in parsed.faces() {
        let wrong = changed_record(&text, face.entity(), |record| {
            if face.forward() {
                record.replacen(",.T.)", ",.F.)", 1)
            } else {
                record.replacen(",.F.)", ",.T.)", 1)
            }
        });
        let document = CompoundDocument::parse(&wrong, DocumentOptions::default(), None).unwrap();
        let parsed = Parsed::parse(
            &document,
            document.represented_members()[0].item_entity(),
            Default::default(),
            None,
        );
        match face.support() {
            Support::Tensor(_) => assert!(matches!(parsed, Err(Error::Unsupported))),
            Support::Plane { .. } => assert!(matches!(
                parsed.unwrap().projection(Default::default(), None),
                Err(Error::Unsupported)
            )),
        }
    }
}
#[test]
fn actual_interior_z_coefficients_must_earn_the_global_cap_range() {
    let text = super::super::tests::fixture(false);
    let document = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
    let parsed = Parsed::parse(
        &document,
        document.represented_members()[0].item_entity(),
        Default::default(),
        None,
    )
    .unwrap();
    let tensor = &parsed.tensors()[0];
    let point = parsed.points()[tensor.point_indices()[11 * tensor.spans()[1] + 2]];
    let coordinates = point.coordinates();
    let wrong = changed_record(&text, point.entity(), |_| {
        format!(
            "#{}=CARTESIAN_POINT('',({},{},20.));",
            point.entity(),
            coordinates[0].spelling(),
            coordinates[1].spelling()
        )
    });
    let document = CompoundDocument::parse(&wrong, DocumentOptions::default(), None).unwrap();
    let parsed = Parsed::parse(
        &document,
        document.represented_members()[0].item_entity(),
        Default::default(),
        None,
    )
    .unwrap();
    assert!(matches!(
        parsed.projection(Default::default(), None),
        Err(Error::Unsupported)
    ));
}
#[test]
fn actual_mixed_z_derivative_is_retained_in_the_flux_remainder() {
    let poles = std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            [
                I::point(11.),
                I::point(-(i as f64 + j as f64)),
                I::point((i + 2 * j) as f64),
            ]
        })
    });
    let mut budget = Budget::new(Default::default(), None).unwrap();
    let matrix = super::super::moments::Moments::new(&mut budget).unwrap();
    let flux = cell_flux(&poles, 0., 22., &matrix, &mut budget).unwrap();
    // x=11, y=-11(u+v), z=11u+22v: J=-1331 exactly.
    // Dropping z_u would incorrectly return 1331 instead of 1331/2.
    assert!(flux.lo <= 665.5 && flux.hi >= 665.5, "{flux:?}");
    assert!(flux.lo > 665.49 && flux.hi < 1996.51);
}
fn basis(n: usize, x: f64) -> Vec<f64> {
    let mut values = vec![0.; n + 1];
    values[0] = 1.;
    for degree in 1..=n {
        for i in (1..=degree).rev() {
            values[i] = values[i] * (1. - x) + values[i - 1] * x;
        }
        values[0] *= 1. - x;
    }
    values
}
fn quadrature(n: usize) -> Vec<(f64, f64)> {
    (0..n)
        .map(|i| {
            let mut x = (std::f64::consts::PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
            let mut derivative = 0.;
            for _ in 0..64 {
                let (mut a, mut b) = (1., x);
                for k in 2..=n {
                    (a, b) = (
                        b,
                        ((2 * k - 1) as f64 * x * b - (k - 1) as f64 * a) / k as f64,
                    );
                }
                derivative = n as f64 * (x * b - a) / (x * x - 1.);
                let next = x - b / derivative;
                if next == x {
                    break;
                }
                x = next;
            }
            (
                (1. + x) * 0.5,
                1. / ((1. - x * x) * derivative * derivative),
            )
        })
        .collect()
}
// Test-only independent direct evaluation of the complete literal polynomial J.
// A degree11 tensor gives J degree <=32 in each variable, so 17-point Gauss
// quadrature is exact in real arithmetic. Compare a separate 20-point rule too.
fn direct_flux(parsed: &Parsed<'_, '_>, order: usize) -> f64 {
    let samples: Vec<_> = quadrature(order)
        .into_iter()
        .map(|(x, w)| (basis(11, x), basis(10, x), w))
        .collect();
    let mut sum = 0.;
    let mut compensation = 0.;
    for tensor in &parsed.tensors {
        for u in 0..tensor.spans[0] {
            for v in 0..tensor.spans[1] {
                let height = 11 * tensor.spans[1] + 1;
                let p: [[[f64; 3]; 12]; 12] = std::array::from_fn(|i| {
                    std::array::from_fn(|j| {
                        parsed.points[tensor.poles[(11 * u + i) * height + 11 * v + j]]
                            .coordinates
                            .map(|s| s.value())
                    })
                });
                for (bu, du, wu) in &samples {
                    for (bv, dv, wv) in &samples {
                        let (mut point, mut a, mut b) = ([0.; 3], [0.; 3], [0.; 3]);
                        for i in 0..12 {
                            for j in 0..12 {
                                for axis in 0..3 {
                                    point[axis] += bu[i] * bv[j] * p[i][j][axis];
                                    if i < 11 {
                                        a[axis] += 11.
                                            * du[i]
                                            * bv[j]
                                            * (p[i + 1][j][axis] - p[i][j][axis]);
                                    }
                                    if j < 11 {
                                        b[axis] += 11.
                                            * bu[i]
                                            * dv[j]
                                            * (p[i][j + 1][axis] - p[i][j][axis]);
                                    }
                                }
                            }
                        }
                        let jacobian = (point[0] * a[1] - point[1] * a[0]) * b[2]
                            - (point[0] * b[1] - point[1] * b[0]) * a[2];
                        let term = -0.5 * wu * wv * jacobian - compensation;
                        let next = sum + term;
                        compensation = (next - sum) - term;
                        sum = next;
                    }
                }
            }
        }
    }
    sum
}
#[test]
fn literal_material_queries_enforce_resources_and_cancellation() {
    let text = super::super::tests::fixture(false);
    let document = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
    let parsed = Parsed::parse(
        &document,
        document.represented_members()[0].item_entity(),
        Default::default(),
        None,
    )
    .unwrap();
    let standard = MaterialOptions::default();
    for options in [
        MaterialOptions {
            maximum_work: 0,
            ..standard
        },
        MaterialOptions {
            maximum_work: 64_000_001,
            ..standard
        },
        MaterialOptions {
            maximum_bytes: 0,
            ..standard
        },
        MaterialOptions {
            maximum_bytes: 64 * 1024 * 1024 + 1,
            ..standard
        },
        MaterialOptions {
            maximum_nodes: 65_537,
            ..standard
        },
        MaterialOptions {
            maximum_depth: 33,
            ..standard
        },
    ] {
        assert!(matches!(
            parsed.projection(options, None),
            Err(Error::Invalid)
        ));
        let projection = parsed.projection(standard, None).unwrap();
        assert!(matches!(
            projection.occupied(options, None),
            Err(Error::Invalid)
        ));
    }
    for options in [
        MaterialOptions {
            maximum_work: 1,
            ..standard
        },
        MaterialOptions {
            maximum_bytes: 1,
            ..standard
        },
        MaterialOptions {
            maximum_nodes: 0,
            ..standard
        },
    ] {
        assert!(matches!(
            parsed.projection(options, None),
            Err(Error::Limit)
        ));
    }
    assert!(matches!(
        parsed.projection(
            MaterialOptions {
                maximum_depth: 0,
                ..standard
            },
            None
        ),
        Err(Error::Unsupported)
    ));
    let projection = parsed.projection(standard, None).unwrap();
    let work = projection.work_used();
    assert!(matches!(
        projection.occupied(
            MaterialOptions {
                maximum_work: work,
                ..standard
            },
            None
        ),
        Err(Error::Limit)
    ));
    let projection = parsed.projection(standard, None).unwrap();
    let bytes = projection.query_budget_bytes();
    assert!(matches!(
        projection.occupied(
            MaterialOptions {
                maximum_bytes: bytes - 1,
                ..standard
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        parsed.projection(standard, Some(&|| true)),
        Err(Error::Cancelled)
    ));
    let projection = parsed.projection(standard, None).unwrap();
    assert!(matches!(
        projection.occupied(standard, Some(&|| true)),
        Err(Error::Cancelled)
    ));
    // Real cancellation boundaries in cap-loop and tensor subdivision, without
    // injecting a fabricated prevalidated graph or numeric failure.
    let calls = std::cell::Cell::new(0usize);
    parsed
        .projection(
            standard,
            Some(&|| {
                calls.set(calls.get() + 1);
                false
            }),
        )
        .unwrap();
    for checkpoint in [2, calls.get() / 4, calls.get() / 2, calls.get() - 1] {
        let count = std::cell::Cell::new(0usize);
        let stop = || {
            count.set(count.get() + 1);
            count.get() >= checkpoint
        };
        assert!(matches!(
            parsed.projection(standard, Some(&stop)),
            Err(Error::Cancelled)
        ));
    }
}
#[test]
fn actual_literal_interior_warp_is_integrated_without_reconstructing_source_height() {
    let text = super::super::tests::fixture(false);
    let document = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
    let parsed = Parsed::parse(
        &document,
        document.represented_members()[0].item_entity(),
        Default::default(),
        None,
    )
    .unwrap();
    let tensor = &parsed.tensors()[0];
    let point = parsed.points()[tensor.point_indices()[11 * tensor.spans()[1] + 2]];
    let coordinates = point.coordinates();
    let warped = changed_record(&text, point.entity(), |_| {
        format!(
            "#{}=CARTESIAN_POINT('',({},{},{}));",
            point.entity(),
            coordinates[0].spelling(),
            coordinates[1].spelling(),
            coordinates[2].value() + 0.00002
        )
    });
    let document = CompoundDocument::parse(&warped, DocumentOptions::default(), None).unwrap();
    let parsed = Parsed::parse(
        &document,
        document.represented_members()[0].item_entity(),
        Default::default(),
        None,
    )
    .unwrap();
    let occupied = parsed
        .projection(Default::default(), None)
        .unwrap()
        .occupied(Default::default(), None)
        .unwrap();
    let direct = direct_flux(&parsed, 17);
    assert!(
        occupied.volume_mm3()[0] <= direct && direct <= occupied.volume_mm3()[1],
        "{direct} {:?}",
        occupied.volume_mm3()
    );
    assert!((direct - direct_flux(&parsed, 20)).abs() < 1e-10);
}
#[test]
fn decoded_actual_cell_frontiers_and_decimal_range_are_checked() {
    let text = super::super::tests::fixture(false);
    let document = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
    let parsed = Parsed::parse(
        &document,
        document.represented_members()[0].item_entity(),
        Default::default(),
        None,
    )
    .unwrap();
    let poles = controls(&parsed, &parsed.tensors[0], 0, 0);
    for (maximum_nodes, maximum_depth, expected) in
        [(0, 32, Error::Limit), (65_536, 0, Error::Unsupported)]
    {
        let mut state = State {
            budget: Budget::new(Default::default(), None).unwrap(),
            nodes: 0,
            leaves: 0,
            maximum_nodes,
            maximum_depth,
        };
        let cell = calculation_cell11(std::array::from_fn(|i| poles[i / 12][i % 12])).unwrap();
        assert_eq!(state.cell(cell), Err(expected));
    }
    let tensor = &parsed.tensors()[0];
    let point = parsed.points()[tensor.point_indices()[11 * tensor.spans()[1] + 2]];
    let coordinates = point.coordinates();
    let wrong = changed_record(&text, point.entity(), |_| {
        format!(
            "#{}=CARTESIAN_POINT('',({},{},-20.));",
            point.entity(),
            coordinates[0].spelling(),
            coordinates[1].spelling()
        )
    });
    let document = CompoundDocument::parse(&wrong, DocumentOptions::default(), None).unwrap();
    let decoded = Parsed::parse(
        &document,
        document.represented_members()[0].item_entity(),
        Default::default(),
        None,
    )
    .unwrap();
    assert!(matches!(
        decoded.projection(Default::default(), None),
        Err(Error::Unsupported)
    ));
    let out_of_range = changed_record(&text, point.entity(), |_| {
        format!(
            "#{}=CARTESIAN_POINT('',(9999999.9999999999,{},{}));",
            point.entity(),
            coordinates[1].spelling(),
            coordinates[2].spelling()
        )
    });
    let document =
        CompoundDocument::parse(&out_of_range, DocumentOptions::default(), None).unwrap();
    assert!(matches!(
        Parsed::parse(
            &document,
            document.represented_members()[0].item_entity(),
            Default::default(),
            None
        ),
        Err(Error::Limit)
    ));
    let record = text
        .lines()
        .find(|line| line.contains("=CLOSED_SHELL("))
        .unwrap();
    let shell: u32 = record[1..record.find('=').unwrap()].parse().unwrap();
    let reordered = changed_record(&text, shell, |_| {
        format!(
            "#{shell}=CLOSED_SHELL('',({}));",
            parsed
                .faces
                .iter()
                .rev()
                .map(|face| format!("#{}", face.id))
                .collect::<Vec<_>>()
                .join(",")
        )
    });
    let document = CompoundDocument::parse(&reordered, DocumentOptions::default(), None).unwrap();
    let parsed = Parsed::parse(
        &document,
        document.represented_members()[0].item_entity(),
        Default::default(),
        None,
    )
    .unwrap();
    assert!(
        parsed
            .projection(Default::default(), None)
            .unwrap()
            .occupied(Default::default(), None)
            .is_ok()
    );
}
fn warp_height(text: &str, parsed: &Parsed<'_, '_>, power: i32) -> String {
    let height = parsed
        .points
        .iter()
        .map(|p| p.coordinates[2].value())
        .fold(0_f64, f64::max);
    let records: std::collections::BTreeMap<_, _> = parsed
        .points
        .iter()
        .map(|point| {
            let c = point.coordinates();
            (
                point.entity(),
                format!(
                    "#{}=CARTESIAN_POINT('',({},{},{}));",
                    point.entity(),
                    c[0].spelling(),
                    c[1].spelling(),
                    height
                        * (0.01 * (c[2].value() / height)
                            + 0.99 * (c[2].value() / height).powi(power))
                ),
            )
        })
        .collect();
    text.lines()
        .map(|line| {
            line.strip_prefix('#')
                .and_then(|s| s.split_once('='))
                .and_then(|(id, _)| id.parse::<u32>().ok())
                .and_then(|id| records.get(&id))
                .map_or(line, String::as_str)
        })
        .collect::<Vec<_>>()
        .join("\n")
}
#[test]
fn actual_nonlinear_height_requires_a_positive_resolved_volume_interval() {
    let text = super::super::tests::fixture(false);
    let document = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
    let parsed = Parsed::parse(
        &document,
        document.represented_members()[0].item_entity(),
        Default::default(),
        None,
    )
    .unwrap();
    // A monotone height warp preserves cap endpoints and shared owners. Its
    // positive linear component keeps the actual closed-cell z_v resolved.
    let warped = warp_height(&text, &parsed, 32);
    let document = CompoundDocument::parse(&warped, DocumentOptions::default(), None).unwrap();
    let parsed = Parsed::parse(
        &document,
        document.represented_members()[0].item_entity(),
        Default::default(),
        None,
    )
    .unwrap();
    let projection = parsed.projection(Default::default(), None).unwrap();
    assert!(direct_flux(&parsed, 17) > 0.);
    assert!(matches!(
        projection.occupied(Default::default(), None),
        Err(Error::Unsupported)
    ));
}

pub(in crate::step_read::npt_stock) fn scale_xy(
    text: &str,
    document: &CompoundDocument<'_>,
    factor: i32,
) -> String {
    let mut cap_points = std::collections::BTreeSet::new();
    for entries in document.records().values() {
        if let [entity] = entries.as_slice()
            && entity.name == "B_SPLINE_CURVE_WITH_KNOTS"
            && matches!(entity.args.get(1),Some(Value::Number(n)) if *n==11.)
            && let Value::List(poles) = &entity.args[2]
        {
            for pole in poles {
                if let Value::Reference(id) = pole {
                    cap_points.insert(*id);
                }
            }
        }
    }
    // Integer decimal multiplication is exact, including original trailing
    // digits that disappear from the binary64 cache. Scale by factor*10^-162.
    let scaled = |text: &str| {
        let (mantissa, exponent) = text.split_once('E').unwrap_or((text, "0"));
        let fraction = mantissa
            .split_once('.')
            .map_or(0, |(_, fraction)| fraction.len());
        let digits = mantissa
            .replace('.', "")
            .parse::<num_bigint::BigInt>()
            .unwrap()
            * factor;
        format!(
            "{digits}E{}",
            exponent.parse::<i32>().unwrap() - 162 - fraction as i32
        )
    };
    let mut replacements = std::collections::BTreeMap::new();
    for (id, entries) in document.records() {
        if let [entity] = entries.as_slice()
            && entity.name == "CARTESIAN_POINT"
            && let Some(Value::List(coords)) = entity.args.get(1)
            && (coords.len() == 3 || (coords.len() == 2 && cap_points.contains(id)))
        {
            let numbers = part21::numeric_spellings(entity).unwrap();
            let mut values = vec![scaled(numbers[0]), scaled(numbers[1])];
            if coords.len() == 3 {
                let (mantissa, exponent) = numbers[2].split_once('E').unwrap_or((numbers[2], "0"));
                values.push(format!(
                    "{mantissa}E{}",
                    exponent.parse::<i32>().unwrap() + 3
                ));
            }
            replacements.insert(
                *id,
                format!("#{id}=CARTESIAN_POINT('',({}));", values.join(",")),
            );
        }
    }
    text.lines()
        .map(|line| {
            line.strip_prefix('#')
                .and_then(|s| s.split_once('='))
                .and_then(|(id, _)| id.parse::<u32>().ok())
                .and_then(|id| replacements.get(&id))
                .map_or(line, String::as_str)
        })
        .collect::<Vec<_>>()
        .join("\n")
}
#[test]
fn actual_subnormal_cap_areas_require_positive_resolved_integrals() {
    let text = super::super::tests::fixture(false);
    let document = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
    for (factor, accepted) in [(7, false), (9, true)] {
        // Exact decimal XY scaling preserves radial geometry; Z*1000 keeps
        // the actual local Jacobian above underflow while cap area is tiny.
        let scaled = scale_xy(&text, &document, factor);
        let document = CompoundDocument::parse(&scaled, DocumentOptions::default(), None).unwrap();
        let parsed = Parsed::parse(
            &document,
            document.represented_members()[0].item_entity(),
            Default::default(),
            None,
        )
        .unwrap();
        let projection = parsed.projection(Default::default(), None).unwrap();
        let occupied = projection.occupied(Default::default(), None);
        if accepted {
            let occupied = occupied.unwrap();
            assert!(
                occupied
                    .cap_areas_mm2()
                    .iter()
                    .all(|area| area[0] > 0. && area[1] < 1e-318)
            );
        } else {
            assert!(matches!(occupied, Err(Error::Unsupported)));
        }
    }
}
