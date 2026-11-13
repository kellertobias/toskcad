use super::*;
#[test]
fn real_metric_body_keeps_its_actual_raw_tensors_and_explicit_virtual_caps() {
    for index in 0..6 {
        for left in [false, true] {
            let source = super::super::tests::fixture(index, left, false, 1);
            let (o, u, v) = if left {
                ([105., -20., 60.], [0.8, 0.6, 0.], [0.36, -0.48, -0.8])
            } else {
                ([0.; 3], [1., 0., 0.], [0., 1., 0.])
            };
            let text = crate::step::write_phillips_screw(&source, o, u, v, 1e-7).unwrap();
            let body = read_phillips_body(&text, 1e-7, 1_000_000)
                .unwrap_or_else(|e| panic!("index{index} left{left}: {e:?}"));
            assert_eq!(body.literal_input(), text);
            assert!(body.literal_entity_count() > 1000);
            assert_eq!(body.material_face_entities().len(), source.faces().len());
            assert_eq!(body.tensor_face_entities().len(), 4);
            assert_eq!(body.inferred_thread_spec().left_hand, left);
            assert!(
                (body.inferred_thread_spec().major_radius_mm
                    - source.thread().spec().major_radius_mm)
                    .abs()
                    < 1e-12
            );
            assert!(
                (body.inferred_thread_spec().length_mm - source.thread().spec().length_mm).abs()
                    < 1e-12
            );
            assert_eq!(
                body.profile_face_entities().len(),
                source.profile().segments().count()
            );
            assert_eq!(
                body.retained_meridian_header().segments().count(),
                source.profile().segments().count()
            );
            let cap = body.virtual_outer_cap();
            assert_eq!(cap.material_face, body.cavity().mouth().material_face);
            assert_eq!(cap.surface, body.cavity().mouth().material_surface);
            assert_eq!(
                cap.removed_inner_bound,
                body.cavity().mouth().material_bound
            );
            assert_ne!(cap.retained_outer_bound, cap.removed_inner_bound);
            assert_eq!(cap.removed_cavity_faces.len(), 28);
            assert!(body.literal_tensor_control_points(0).unwrap().len() > 11);
            assert!(body.literal_tensor_control_points(99).is_err());
            assert!(body.inferred_visible_interval_mm()[0] > 0.);
            assert!(
                body.inferred_thread_frame()
                    .into_iter()
                    .flatten()
                    .all(f64::is_finite)
            );
        }
    }
}
#[test]
fn circular_root_and_multistart_headers_keep_the_original_support_datum_and_law() {
    for left in [false, true] {
        let base = super::super::tests::fixture(1, left, false, 1);
        let mut spec = base.thread().spec();
        spec.major_radius_mm = 1.4635;
        spec.pitch_mm = 0.5;
        spec.length_mm = 9.;
        spec.depth_mm = 0.29908756729740676;
        spec.crest_width_mm = 0.07131009168581062;
        spec.root_radius_mm = Some(0.07216878364870323);
        let thread = crate::helical_thread::Thread::new(spec).unwrap();
        let source = crate::phillips_drive::PhillipsScrew::new(
            thread,
            *base.profile(),
            base.spec(),
            base.cavity().clone(),
            100_000,
            &|| false,
        )
        .unwrap();
        let text = crate::step::write_phillips_screw(
            &source,
            [105., -20., 60.],
            [0.8, 0.6, 0.],
            [0.36, -0.48, -0.8],
            1e-7,
        )
        .unwrap();
        let b = read_phillips_body(&text, 1e-7, 1_000_000).unwrap();
        let decoded = b.inferred_thread_spec();
        assert_eq!(decoded.left_hand, left);
        assert!((decoded.length_mm - 9.).abs() < 1e-12);
        assert!((decoded.root_radius_mm.unwrap() - spec.root_radius_mm.unwrap()).abs() < 1e-12);
        assert!((decoded.depth_mm - spec.depth_mm).abs() < 1e-12);
        assert!(b.inferred_visible_interval_mm()[0] > 0.);
        let source = super::super::tests::fixture(1, left, true, 2);
        let text =
            crate::step::write_phillips_screw(&source, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7)
                .unwrap();
        let b = read_phillips_body(&text, 1e-7, 1_000_000).unwrap();
        assert_eq!(b.tensor_face_entities().len(), 8);
        assert_eq!(b.inferred_thread_spec().starts, 2);
        assert_eq!(b.material_face_entities().len(), source.faces().len() + 1);
    }
}
#[test]
fn a_changed_literal_tensor_is_retained_and_never_replaced_by_the_inferred_header() {
    let source = super::super::tests::fixture(1, false, false, 1);
    let text =
        crate::step::write_phillips_screw(&source, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7)
            .unwrap();
    let b = read_phillips_body(&text, 1e-7, 1_000_000).unwrap();
    let shell = Shell::new(&b.entities, 1e-7, true);
    let (surface, _, _) = helical::face(&shell, b.tensor_face_entities()[0]).unwrap();
    let Value::List(rows) = &shell.args(surface, "B_SPLINE_SURFACE_WITH_KNOTS").unwrap()[3] else {
        panic!("literal rows")
    };
    let id = refs(&rows[1]).unwrap()[1];
    let p = shell.point(id).unwrap();
    let old = text
        .lines()
        .find(|l| l.starts_with(&format!("#{id}=")))
        .unwrap();
    let changed = text.replacen(
        old,
        &format!(
            "#{id}=CARTESIAN_POINT('',({:.16E},{:.16E},{:.16E}));",
            p[0] + 0.001,
            p[1],
            p[2]
        ),
        1,
    );
    let changed = read_phillips_body(&changed, 1e-7, 1_000_000).unwrap();
    assert_eq!(changed.inferred_thread_spec(), b.inferred_thread_spec());
    assert_eq!(
        changed.literal_tensor_control_points(0).unwrap()[1][1],
        [p[0] + 0.001, p[1], p[2]]
    );
    // This packet deliberately remains scaffolding and retains the changed net.
    // It is not a tensor-match, manufacturing or physical mass certificate.
    assert_ne!(
        changed.literal_tensor_control_points(0).unwrap(),
        b.literal_tensor_control_points(0).unwrap()
    );
    assert!(read_phillips_body(&text, 1e-7, 0).is_err());
    assert_eq!(
        read_phillips_body_with_interrupt(&text, 1e-7, 1_000_000, Some(&|| true)).unwrap_err(),
        PhillipsReadError::Cancelled
    );
}
fn argument(text: &str, id: u32, index: usize, replacement: &str) -> String {
    let line = text
        .lines()
        .find(|l| l.starts_with(&format!("#{id}=")))
        .unwrap();
    let mut depth = 0;
    let mut quoted = false;
    let mut field = 0;
    let mut start = line.find('(').unwrap() + 1;
    let mut range = None;
    let begin = start;
    for (i, c) in line.char_indices().skip_while(|(i, _)| *i < begin) {
        if c == '\'' {
            quoted = !quoted;
        }
        if quoted {
            continue;
        }
        if c == '(' {
            depth += 1;
        }
        if c == ')' {
            if depth == 0 {
                if field == index {
                    range = Some((start, i));
                }
                break;
            }
            depth -= 1;
        }
        if c == ',' && depth == 0 {
            if field == index {
                range = Some((start, i));
                break;
            }
            field += 1;
            start = i + 1;
        }
    }
    let (a, b) = range.unwrap();
    text.replacen(
        line,
        &format!("{}{}{}", &line[..a], replacement, &line[b..]),
        1,
    )
}
#[test]
fn actual_malformed_outer_metadata_and_orphan_topology_are_refused_without_rewriting_the_file() {
    let source = super::super::tests::fixture(1, false, false, 1);
    let text =
        crate::step::write_phillips_screw(&source, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7)
            .unwrap();
    let body = read_phillips_body(&text, 1e-7, 1_000_000).unwrap();
    let shell = Shell::new(&body.entities, 1e-7, true);
    let count = body.tensor_face_count;
    let (surface, _, _) = helical::face(&shell, body.material_faces[count]).unwrap();
    let line = text
        .lines()
        .find(|l| l.starts_with(&format!("#{surface}=")))
        .unwrap();
    let bad = text.replacen(line, &format!("#{surface}=PLANE('');"), 1);
    assert!(read_phillips_body(&bad, 1e-7, 1_000_000).is_err());
    let (tip, _, _) = helical::face(&shell, body.material_faces[count + 1]).unwrap();
    let plane_args = shell.args(tip, "PLANE").unwrap();
    let axis = reference(&plane_args[1]).unwrap();
    let origin = reference(&shell.args(axis, "AXIS2_PLACEMENT_3D").unwrap()[1]).unwrap();
    let p = shell.point(origin).unwrap();
    let old = text
        .lines()
        .find(|l| l.starts_with(&format!("#{origin}=")))
        .unwrap();
    let changed = text.replacen(
        old,
        &format!(
            "#{origin}=CARTESIAN_POINT('',({:.16E},{:.16E},{:.16E}));",
            p[0] + 0.001,
            p[1],
            p[2]
        ),
        1,
    );
    assert!(read_phillips_body(&changed, 1e-7, 1_000_000).is_err());
    let profile_face = body.profile_faces[0];
    let (profile_surface, _, _) = helical::face(&shell, profile_face).unwrap();
    let first = shell.simple(profile_surface).unwrap();
    let axis = reference(&first.args[1]).unwrap();
    let old = text
        .lines()
        .find(|l| l.starts_with(&format!("#{profile_surface}=")))
        .unwrap();
    let changed = text.replacen(old, &format!("#{profile_surface}=PLANE('',#{axis});"), 1);
    assert!(read_phillips_body(&changed, 1e-7, 1_000_000).is_err());
    let mut orphan = text.clone();
    orphan.insert_str(
        text.rfind("ENDSEC;").unwrap(),
        &format!("#10000000=VERTEX_POINT('',#{origin});\n"),
    );
    assert!(read_phillips_body(&orphan, 1e-7, 1_000_000).is_err());
    // Primitive prefix guard: actual parsed faces, no complete four-band period.
    assert!(inferred_thread(&shell, &body.material_faces[..3]).is_err());
    let excessive = vec![body.material_faces[0]; 36];
    assert!(inferred_thread(&shell, &excessive).is_err());
}
#[test]
fn an_unverified_later_period_net_remains_retrievable_as_literal_data_or_a_schema_error() {
    let source = super::super::tests::fixture(1, false, true, 2);
    let text =
        crate::step::write_phillips_screw(&source, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7)
            .unwrap();
    let body = read_phillips_body(&text, 1e-7, 1_000_000).unwrap();
    let shell = Shell::new(&body.entities, 1e-7, true);
    let (surface, _, _) = helical::face(&shell, body.tensor_face_entities()[4]).unwrap();
    let changed = argument(&text, surface, 3, "0.");
    let packet = read_phillips_body(&changed, 1e-7, 1_000_000).unwrap();
    assert!(packet.literal_tensor_control_points(4).is_err());
    // No primitive/header path claims this malformed later net was qualified.
    assert!(
        packet
            .literal_input()
            .contains("B_SPLINE_SURFACE_WITH_KNOTS")
    );
}
