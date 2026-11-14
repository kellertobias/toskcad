use super::*;
fn packet_text() -> String {
    let screw = super::super::super::tests::fixture(0, false, false, 2);
    crate::step::write_phillips_screw(&screw, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7).unwrap()
}
#[test]
fn actual_unqualified_later_tensor_metadata_has_independent_schema_refusals() {
    let text = packet_text();
    let body = read_phillips_body(&text, 1e-7, 1_000_000).unwrap();
    let shell = Shell::new(&body.entities, 1e-7, true);
    let (surface, _, _) = helical::face(&shell, body.tensor_face_entities()[4]).unwrap();
    let args = shell.args(surface, "B_SPLINE_SURFACE_WITH_KNOTS").unwrap();
    let rows = list(&args[3]).unwrap();
    let first = refs(&rows[0]).unwrap()[0];
    let short_rows = format!(
        "({})",
        std::iter::repeat_n(format!("(#{first})"), 12)
            .collect::<Vec<_>>()
            .join(",")
    );
    for replacement in [
        "0.".to_owned(),
        "(0.)".to_owned(),
        "()".to_owned(),
        format!("((#{first}))"),
        short_rows,
    ] {
        let changed = super::super::tests::argument(&text, surface, 3, &replacement);
        // Body deliberately retains later unqualified carriers. This exercises
        // metadata validation, never a manufactured Solid or a Common receipt.
        let body = read_phillips_body(&changed, 1e-7, 1_000_000).unwrap();
        let mut shell = Shell::new(&body.entities, 1e-7, true);
        let mut work = Work {
            remaining: WORK,
            interrupt: None,
            cancelled: Cell::new(false),
        };
        assert!(
            grids(&body, &mut shell, &mut work, BYTES).is_err(),
            "{replacement}"
        );
    }
}
#[test]
fn actual_partial_packet_ambiguous_knot_entry_refuses_without_using_interior_poles() {
    let text = packet_text();
    let body = read_phillips_body(&text, 1e-7, 1_000_000).unwrap();
    let shell = Shell::new(&body.entities, 1e-7, true);
    let (surface, _, _) = helical::face(&shell, body.tensor_face_entities()[4]).unwrap();
    let args = shell.args(surface, "B_SPLINE_SURFACE_WITH_KNOTS").unwrap();
    let mut ids = list(&args[3])
        .unwrap()
        .iter()
        .map(refs)
        .collect::<Result<Vec<_>>>()
        .unwrap();
    let endpoint = shell.point(ids[0][11]).unwrap();
    for (index, ambiguous) in [(12, false), (22, true)] {
        ids[0][index] = ids[0][11];
        let grid = format!(
            "({})",
            ids.iter()
                .map(|row| format!(
                    "({})",
                    row.iter()
                        .map(|id| format!("#{id}"))
                        .collect::<Vec<_>>()
                        .join(",")
                ))
                .collect::<Vec<_>>()
                .join(",")
        );
        let changed = super::super::tests::argument(&text, surface, 3, &grid);
        let body = read_phillips_body(&changed, 1e-7, 1_000_000).unwrap();
        let mut shell = Shell::new(&body.entities, 1e-7, true);
        let mut work = Work {
            remaining: WORK,
            interrupt: None,
            cancelled: Cell::new(false),
        };
        let grids = grids(&body, &mut shell, &mut work, BYTES).unwrap();
        let result = entry_index(&shell, &grids[4].ids[0], endpoint, &mut work);
        if ambiguous {
            assert!(matches!(result, Err(ReadError::Unsupported)))
        } else {
            assert_eq!(result.unwrap(), 11);
        }
    }
}
#[test]
fn real_qualified_chart_inverse_encloses_exact_independent_point_and_refuses_overflow() {
    use num_rational::BigRational as Q;
    let text = super::super::super::branches::representative::scaled_text(1.137, true);
    super::super::common_mouth::tests::derived_thread_probe::with_head(&text, &mut |head| {
        let chart = chart::Chart::new(head);
        let raw = head.body().literal_tensor_control_points(0).unwrap()[1][11];
        let [o, u, v, n] = head
            .chart()
            .columns
            .map(|p| p.map(|x| Q::from_float(x).unwrap()));
        let cross = |a: &[Q; 3], b: &[Q; 3]| {
            std::array::from_fn::<_, 3, _>(|i| {
                &a[(i + 1) % 3] * &b[(i + 2) % 3] - &a[(i + 2) % 3] * &b[(i + 1) % 3]
            })
        };
        let dot = |a: &[Q; 3], b: &[Q; 3]| {
            (0..3).fold(Q::from_integer(0.into()), |s, i| s + &a[i] * &b[i])
        };
        let d = std::array::from_fn(|i| Q::from_float(raw[i]).unwrap() - &o[i]);
        let det = dot(&u, &cross(&v, &n));
        let mut exact = [
            dot(&d, &cross(&v, &n)) / &det,
            dot(&d, &cross(&n, &u)) / &det,
            dot(&d, &cross(&u, &v)) / &det,
        ];
        exact[2] += Q::from_float(head.chart().outer_top_z).unwrap();
        let local = chart.local(raw).unwrap();
        for (x, q) in local.into_iter().zip(exact) {
            assert!(Q::from_float(x.lo).unwrap() <= q && q <= Q::from_float(x.hi).unwrap());
        }
        for p in [[f64::MAX; 3], [1e8, 0., 0.], [f64::NAN, 0., 0.]] {
            assert!(chart.local(p).is_err());
        }
    });
}
#[test]
fn actual_repaired_polynomial_scaled_beyond_retained_domain_is_refused_before_publication() {
    let text = super::super::super::branches::representative::scaled_text(1., false);
    super::super::common_mouth::tests::derived_thread_probe::with_head(&text, &mut |head| {
        let source = head
            .derived_thread_carrier(Default::default(), None)
            .unwrap();
        let tensor = source
            .source()
            .surfaces
            .iter()
            .find_map(|s| match s {
                radial::Surface::Tensor(t) => Some(t),
                _ => None,
            })
            .unwrap();
        let height = 9 * tensor.spans()[1] + 1;
        // This is an actual source end cell, with positive z near2mm. Both
        // transformations are finite scalings; the second overflows results.
        let q = std::array::from_fn::<_, 10, _>(|i| {
            std::array::from_fn::<_, 10, _>(|j| tensor.poles()[i * height + height - 10 + j])
        });
        proposal_domain(&q).unwrap();
        for scale in [1e8, f64::MAX] {
            let scaled = q.map(|row| row.map(|p| p.map(|x| x * scale)));
            assert!(matches!(proposal_domain(&scaled), Err(ReadError::Invalid)));
            assert!(matches!(
                RetainedTensor::new([1, 1], scaled.into_iter().flatten().collect()),
                Err(crate::tapered_thread::Error::Invalid)
            ));
        }
    });
}
