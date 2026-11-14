use super::*;
fn fixture() -> String {
    let source = super::super::super::tests::fixture(1, false, false, 1);
    crate::step::write_phillips_screw(
        &source,
        [105., -20., 60.],
        [0.8, 0.6, 0.],
        [0.36, -0.48, -0.8],
        1e-7,
    )
    .unwrap()
}
#[test]
fn original_points_and_curve_scalar_leaves_retain_body_identity_exact_tokens_and_numeric_enclosures()
 {
    let text = fixture();
    let body = read_phillips_body(&text, 1e-7, 1_000_000).unwrap();
    assert!(body.retained_budget_bytes() > text.len() + std::mem::size_of_val(&body));
    assert!(body.retained_budget_bytes() < 64 * 1024 * 1024);
    let shell = Shell::new(&body.entities, 1e-7, true);
    for vertex in body.cavity().vertices() {
        let leaf = body
            .exact_point_leaf(vertex.point_entity, 1000, None)
            .unwrap();
        assert!(std::ptr::eq(leaf.body(), &body));
        assert_eq!(leaf.entity(), vertex.point_entity);
        assert!(leaf.work_used() > 0 && leaf.work_used() < 1000);
        for ((token, bounds), value) in leaf
            .coordinate_spellings()
            .into_iter()
            .zip(leaf.coordinate_bounds_mm())
            .zip(vertex.point)
        {
            assert_eq!(token.parse::<f64>().unwrap(), value);
            assert!(bounds[0] < value && bounds[1] > value);
        }
    }
    let curve = body
        .cavity()
        .edges()
        .iter()
        .find_map(|e| {
            if let Curve::CubicIntersection { spline_entity, .. } = e.curve {
                Some(spline_entity)
            } else {
                None
            }
        })
        .unwrap();
    let numbers = body.original_numeric_leaves(curve, 100_000, None).unwrap();
    assert!(std::ptr::eq(numbers.body(), &body));
    assert_eq!(numbers.entity(), curve);
    assert!(numbers.retained_capacity_bytes() > std::mem::size_of_val(&numbers));
    assert!(numbers.work_used() > 0);
    assert_eq!(numbers.leaves()[0].spelling(), "3");
    assert!(numbers.leaves()[0].bounds()[0] < 3. && numbers.leaves()[0].bounds()[1] > 3.);
    // Real file perturbation below one binary64 ULP remains distinct symbolic
    // source, while the same outward floating evaluation enclosure stays valid.
    let point = body.cavity().vertices()[0].point_entity;
    let before = body.exact_point_leaf(point, 1000, None).unwrap();
    let token = before.coordinate_spellings()[0];
    let (m, e) = token.split_once('E').unwrap();
    let changed = format!("{m}0000000000000000000000000000000001E{e}");
    assert_eq!(
        token.parse::<f64>().unwrap(),
        changed.parse::<f64>().unwrap()
    );
    let raw = shell.simple(point).unwrap().source.unwrap();
    let altered = text.replace(raw, &raw.replacen(token, &changed, 1));
    let after = read_phillips_body(&altered, 1e-7, 1_000_000).unwrap();
    let after_leaf = after.exact_point_leaf(point, 1000, None).unwrap();
    assert_ne!(
        before.coordinate_spellings(),
        after_leaf.coordinate_spellings()
    );
    assert_eq!(
        before.coordinate_bounds_mm(),
        after_leaf.coordinate_bounds_mm()
    );
    // The source arena also retains typed unit-measure numeric leaves.
    let typed = body
        .entities
        .iter()
        .find_map(|(&id, e)| {
            (e.len() == 1 && e[0].args.iter().any(|v| matches!(v, Value::Typed(_)))).then_some(id)
        })
        .unwrap();
    assert!(
        !body
            .original_numeric_leaves(typed, 1000, None)
            .unwrap()
            .leaves()
            .is_empty()
    );
}
fn append(text: &str, record: &str) -> String {
    let mut result = text.to_owned();
    result.insert_str(text.rfind("ENDSEC;").unwrap(), record);
    result
}
#[test]
fn non_point_malformed_numeric_and_resource_queries_never_publish_leaf_authority() {
    let text = fixture();
    let body = read_phillips_body(&text, 1e-7, 1_000_000).unwrap();
    let id = body.cavity().vertices()[0].point_entity;
    assert!(body.original_numeric_leaves(u32::MAX, 1000, None).is_err());
    assert!(
        body.exact_point_leaf(body.material_faces[0], 1000, None)
            .is_err()
    );
    let complex = body
        .entities
        .iter()
        .find_map(|(&id, e)| (e.len() > 1).then_some(id))
        .unwrap();
    assert!(body.original_numeric_leaves(complex, 1000, None).is_err());
    for record in [
        "#10000000=CARTESIAN_POINT('');\n",
        "#10000000=CARTESIAN_POINT('',(1.,2.));\n",
        "#10000000=CARTESIAN_POINT('',(1.7976931348623157E308,0.,0.));\n",
    ] {
        let altered = append(&text, record);
        let parsed = read_phillips_body(&altered, 1e-7, 1_000_000).unwrap();
        assert!(parsed.exact_point_leaf(10000000, 1000, None).is_err());
    }
    assert!(body.exact_point_leaf(id, 0, None).is_err());
    assert!(body.exact_point_leaf(id, 1, None).is_err());
    assert!(matches!(
        body.exact_point_leaf(id, 1000, Some(&|| true)),
        Err(PhillipsReadError::Cancelled)
    ));
    let calls = Cell::new(0);
    let stop = || {
        calls.set(calls.get() + 1);
        calls.get() > 2
    };
    assert!(matches!(
        body.exact_point_leaf(id, 1000, Some(&stop)),
        Err(PhillipsReadError::Cancelled)
    ));
}
#[test]
fn whole_body_budget_is_charged_before_publication_and_unsupported_head_partition_refuses() {
    let text = fixture();
    let body = read_phillips_body(&text, 1e-7, 1_000_000).unwrap();
    let mut work = Work {
        remaining: 1_000_000,
        interrupt: None,
        cancelled: Cell::new(false),
    };
    assert_eq!(
        budget(&body, &mut work).unwrap(),
        body.retained_budget_bytes()
    );
    let cost = 1_000_000 - work.remaining;
    assert!(cost > 1000);
    assert!(read_phillips_body(&text, 1e-7, cost).is_err());
    work.remaining = 0;
    assert!(budget(&body, &mut work).is_err());
    let count = Cell::new(0);
    read_phillips_body_with_interrupt(
        &text,
        1e-7,
        1_000_000,
        Some(&|| {
            count.set(count.get() + 1);
            false
        }),
    )
    .unwrap();
    let full = count.get();
    count.set(0);
    let stop = || {
        count.set(count.get() + 1);
        count.get() + 10 >= full
    };
    assert!(matches!(
        read_phillips_body_with_interrupt(&text, 1e-7, 1_000_000, Some(&stop)),
        Err(PhillipsReadError::Cancelled)
    ));
    // Actual raw shell order no longer supplies the supported outer/cavity
    // decomposition. The geometry remains a structurally valid cavity, so this
    // is an unsupported body-family layout, not a geometric-invalidity claim.
    let raw = body
        .entities
        .iter()
        .find_map(|(&id, e)| (e.len() == 1 && e[0].name == "CLOSED_SHELL").then_some(id))
        .unwrap();
    let mut faces = body.material_faces.clone();
    let a = faces
        .iter()
        .position(|&f| f == body.cap.material_face)
        .unwrap();
    let b = faces
        .iter()
        .position(|&f| Some(f) == body.cavity().faces()[1].material_entity)
        .unwrap();
    faces.swap(a, b);
    let list = format!(
        "({})",
        faces
            .iter()
            .map(|f| format!("#{f}"))
            .collect::<Vec<_>>()
            .join(",")
    );
    let changed = super::super::tests::argument(&text, raw, 1, &list);
    read_phillips_screw_boundary(&changed, 1e-7, 1_000_000).unwrap();
    assert!(read_phillips_body(&changed, 1e-7, 1_000_000).is_err());
}

#[test]
fn exact_point_publication_has_its_own_charged_final_cancel_checkpoint() {
    let text = fixture();
    let body = read_phillips_body(&text, 1e-7, 1_000_000).unwrap();
    let id = body.cavity().vertices()[0].point_entity;
    let count = Cell::new(0);
    let numeric = body
        .original_numeric_leaves(
            id,
            1000,
            Some(&|| {
                count.set(count.get() + 1);
                false
            }),
        )
        .unwrap();
    let generic_calls = count.get();
    count.set(0);
    let interrupt = || {
        count.set(count.get() + 1);
        count.get() > generic_calls
    };
    assert!(matches!(
        body.exact_point_leaf(id, 1000, Some(&interrupt)),
        Err(PhillipsReadError::Cancelled)
    ));
    assert!(
        body.exact_point_leaf(id, numeric.work_used(), None)
            .is_err()
    );
    let leaf = body
        .exact_point_leaf(id, numeric.work_used() + 1, None)
        .unwrap();
    assert_eq!(leaf.work_used(), numeric.work_used() + 1);
}
