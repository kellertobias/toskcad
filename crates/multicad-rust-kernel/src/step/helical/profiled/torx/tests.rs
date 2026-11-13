use super::*;
use crate::hexalobular_drive::{Cone, Drive, FiniteBinding, TangentSpec};
use std::f64::consts::PI;

fn payload(source: Drive) -> ConicalHexalobular {
    let crate::screw_recess::Recess::ConicalHexalobular(payload) =
        crate::screw_recess::Recess::conical_hexalobular(source, 1e-7, 100_000).unwrap()
    else {
        unreachable!()
    };
    payload
}
fn source() -> Drive {
    Drive::from_authored_tangent(
        TangentSpec {
            major_diameter_mm: 2.8,
            minor_diameter_mm: 2.05,
            outer_lobe_radius_mm: 0.28,
        },
        Cone {
            tip_depth_mm: 1.1 + 1.4 / 3_f64.sqrt(),
            depth_per_radius: 1. / 3_f64.sqrt(),
        },
    )
    .unwrap()
}

fn retained_payload(shift: f64) -> ConicalHexalobular {
    let source = source();
    let mut curves: [Vec<_>; 12] =
        std::array::from_fn(|i| source.intersection_carrier(i, 1e-9, 10_000).unwrap());
    let vertices = FiniteVertices {
        mouth_mm: source.arcs().map(|a| [a.start_mm[0], a.start_mm[1], 0.]),
        floor_mm: std::array::from_fn(|i| curves[i][0].poles_mm[0]),
        apex_mm: [0., 0., source.cone().tip_depth_mm],
    };
    curves[0][0].poles_mm[1][2] += shift;
    let posts = std::array::from_fn(|i| {
        let d: Point = std::array::from_fn(|j| vertices.floor_mm[i][j] - vertices.mouth_mm[i][j]);
        FinitePost {
            origin_mm: std::array::from_fn(|j| vertices.mouth_mm[i][j] + 7. * d[j]),
            direction_mm_per_parameter: d.map(|x| -2. * x),
            parameters: [3.5, 3.],
        }
    });
    let binding = FiniteBinding::new_with_posts(source, vertices, posts, curves, 0.01)
        .unwrap()
        .with_literal_cone_semi_angle((PI / 3.).next_up())
        .unwrap();
    let crate::screw_recess::Recess::ConicalHexalobular(payload) =
        crate::screw_recess::Recess::conical_hexalobular_from_binding(binding, 1., 100_000)
            .unwrap()
    else {
        unreachable!()
    };
    payload
}

#[test]
fn retained_encoder_keeps_actual_controls_post_supports_and_original_cone_angle() {
    let recess = retained_payload(2e-10);
    let binding = recess.finite_binding().unwrap();
    let p = Placement::new([10., -20., 30.], [0.8, 0.6, 0.], [0., 0., 1.]).unwrap();
    let top = -2.4;
    let mut file = StepFile::new();
    encode(&mut file, &p, &recess, top, 1e-7).unwrap();
    let curve = file
        .data
        .lines()
        .find(|line| line.contains("B_SPLINE_CURVE_WITH_KNOTS("))
        .unwrap();
    let refs = references(curve);
    let points = std::iter::once(binding.carriers()[0][0].poles_mm[0]).chain(
        binding.carriers()[0]
            .iter()
            .flat_map(|s| s.poles_mm[1..].iter().copied()),
    );
    let expected = points
        .map(|[x, y, z]| {
            let point = p
                .point_f64([x, y], top + z)
                .map(|x| if x == 0. { 0. } else { x });
            file.points[&point.map(f64::to_bits)].0.to_string()
        })
        .collect::<Vec<_>>();
    assert_eq!(refs[1..], expected);
    for post in binding.post_lines() {
        let [x, y, z] = post.origin_mm;
        assert!(
            file.points
                .contains_key(&p.point_f64([x, y], top + z).map(f64::to_bits))
        );
        let value = p.vector(post.direction_mm_per_parameter);
        assert!(
            file.data
                .contains(&format!("DIRECTION('',{})", Coordinates(value)))
        );
    }
    let angle = binding.literal_cone_semi_angle_rad().unwrap();
    assert_ne!(
        angle.to_bits(),
        recess
            .source()
            .cone()
            .depth_per_radius
            .recip()
            .atan()
            .to_bits()
    );
    assert!(file.data.contains(&format!(",0.,{angle:.16E})")));
    assert_eq!(recess.finite_binding().unwrap(), binding);
}

#[test]
fn retained_encoder_refuses_insufficient_precision_or_collapsed_world_posts_without_replacing_sources()
 {
    let too_wide = retained_payload(1e-4);
    let p = Placement::new([0.; 3], [1., 0., 0.], [0., 1., 0.]).unwrap();
    let mut file = StepFile::new();
    let before = file.data.clone();
    assert!(matches!(
        encode(&mut file, &p, &too_wide, -2.4, 1e-7),
        Err(StepError::Limit)
    ));
    assert_eq!(file.data, before);
    let recess = retained_payload(2e-10);
    let collapsed = Placement::new([1e300; 3], [1., 0., 0.], [0., 1., 0.]).unwrap();
    assert!(matches!(
        encode(&mut StepFile::new(), &collapsed, &recess, -2.4, 1e-7),
        Err(StepError::InvalidGeometry)
    ));
    assert!(encode(&mut StepFile::new(), &p, &recess, -2.4, 1e-7).is_ok());
}

#[test]
fn retained_encoder_keeps_positive_post_parameterization_without_a_literal_angle_override() {
    let previous = retained_payload(2e-10);
    let stored = previous.finite_binding().unwrap();
    let binding = FiniteBinding::new(
        stored.source().clone(),
        *stored.vertices(),
        stored.carriers().clone(),
        1e-7,
    )
    .unwrap();
    assert!(
        binding
            .post_lines()
            .iter()
            .all(|p| p.parameters == [0., 1.])
    );
    assert_eq!(binding.literal_cone_semi_angle_rad(), None);
    let crate::screw_recess::Recess::ConicalHexalobular(recess) =
        crate::screw_recess::Recess::conical_hexalobular_from_binding(binding, 1e-5, 100_000)
            .unwrap()
    else {
        unreachable!()
    };
    let p = Placement::new([0.; 3], [1., 0., 0.], [0., 1., 0.]).unwrap();
    assert!(encode(&mut StepFile::new(), &p, &recess, -2.4, 1e-7).is_ok());
}
#[test]
fn torx_encoder_retains_real_surfaces_contacts_and_apex() {
    let source = source();
    let recess = payload(source.clone());
    let p = Placement::new([10., -20., 30.], [0.8, 0.6, 0.], [0., 0., 1.]).unwrap();
    let mut file = StepFile::new();
    let (_mouth, faces) = encode(&mut file, &p, &recess, -2.4, 1e-7).unwrap();
    assert_eq!(faces.len(), 13);
    assert_eq!(file.data.matches("CYLINDRICAL_SURFACE(").count(), 12);
    assert_eq!(file.data.matches("CONICAL_SURFACE(").count(), 1);
    assert_eq!(file.data.matches("CIRCLE(").count(), 12);
    assert_eq!(file.data.matches("INTERSECTION_CURVE(").count(), 12);
    assert_eq!(file.data.matches("B_SPLINE_CURVE_WITH_KNOTS(").count(), 12);
    assert_eq!(file.data.matches("VERTEX_POINT(").count(), 25);
    assert_eq!(file.data.matches("EDGE_CURVE(").count(), 37);
    assert_eq!(file.data.matches("EVALUATED_DEGENERATE_PCURVE(").count(), 1);
    assert_eq!(file.data.matches("SEAM_CURVE(").count(), 1);
    assert_eq!(file.data.matches("SEAM_EDGE(").count(), 2);
    for arc in source.arcs() {
        let xy = arc.start_mm;
        for z in [
            -2.4,
            -2.4 + source.cone().tip_depth_mm - source.cone().depth_per_radius * xy[0].hypot(xy[1]),
        ] {
            let point = p.point_f64(xy, z).map(|v| if v == 0. { 0. } else { v });
            assert!(file.points.contains_key(&point.map(f64::to_bits)));
        }
    }
}

fn head_fixture() -> (String, Screw) {
    let thread = Thread::new(crate::helical_thread::Spec {
        major_radius_mm: 1.5,
        pitch_mm: 0.5,
        length_mm: 3.,
        depth_mm: 0.2,
        crest_width_mm: 0.0625,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts: 1,
        left_hand: false,
        root_radius_mm: None,
    })
    .unwrap();
    let recess = crate::screw_recess::Recess::conical_hexalobular(source(), 1e-7, 100_000).unwrap();
    let screw = Screw::new_with_recess(
        thread,
        1.5,
        crate::head_profile::Head::Crowned {
            radius_mm: 2.8,
            height_mm: 2.4,
            crown_radius_mm: 5.,
            top_land_radius_mm: 2.05,
            neck_radius_mm: 0.1,
        },
        1.,
        0.1,
        recess,
    )
    .unwrap();
    let p = Placement::new([0.; 3], [1., 0., 0.], [0., 1., 0.]).unwrap();
    let mut file = StepFile::new();
    let rim = CircularRim::new(&mut file, [0.; 3], p.normal, p.u, 1.5).unwrap();
    let closing = (rim.edge, rim.center);
    let join = Junction {
        crests: vec![],
        complements: vec![],
        vertex: rim.vertex,
        reference: p.u,
        positive: vec![(rim.edge, true)],
        outer: Some(rim),
    };
    let mut faces = head(&mut file, &p, &screw, join, 1e-7).unwrap();
    assert_eq!(faces.len(), screw.profile().segments().count() + 14);
    assert_eq!(file.data.matches("INTERSECTION_CURVE(").count(), 12);
    let end = circular_bound(&mut file, closing.0, true, true);
    faces.push(file.face(&[end], closing.1, p.normal, p.u));
    let solid = file.solid(&faces);
    file.finish(vec![solid]);
    (file.data, screw)
}

#[test]
fn public_pan_head_export_has_independently_decoded_curved_cavity() {
    let (text, screw) = head_fixture();
    crate::step_read::verify_finite_graph_fixture(&text).unwrap();
    assert!(crate::step::supports_profiled_screw(&screw, [0.; 3], 1e-7));
    assert!(!crate::step::supports_profiled_screw_in_frame(
        &screw,
        [0.; 3],
        [1., 0., 0.],
        [1., 0., 0.],
        1e-7
    ));
    let text = crate::step::write_profiled_screw(&screw, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7)
        .unwrap();
    let (_, imported) = crate::step_read::read_import_solid(&text, 1e-7).unwrap();
    assert!(crate::step_read::profiled_screw_import_matches(
        &imported,
        &screw,
        [[0.; 3], [1., 0., 0.], [0., 1., 0.]],
        1e-7
    ));
}

fn full_fixture() -> (String, Screw) {
    full_fixture_with(false, 1., [[0.; 3], [1., 0., 0.], [0., 1., 0.]])
}
fn full_fixture_with(left: bool, plain: f64, pose: [Point; 3]) -> (String, Screw) {
    let (_, original) = head_fixture();
    let mut spec = original.thread().spec();
    spec.left_hand = left;
    spec.starts = if left { 2 } else { 1 };
    if left {
        spec.major_radius_mm = 1.4635;
        spec.depth_mm = 0.29908756729740676;
        spec.crest_width_mm = 0.07131009168581062;
        spec.root_radius_mm = Some(0.07216878364870323);
    }
    let thread = crate::helical_thread::Thread::new(spec).unwrap();
    let mut head = original.spec().head.unwrap();
    if let crate::head_profile::Head::Crowned {
        ref mut neck_radius_mm,
        ..
    } = head
    {
        *neck_radius_mm = 0.125;
    }
    let screw =
        Screw::new_with_recess(thread, 1.5, head, plain, 0.1, original.recess().clone()).unwrap();
    full_fixture_source(screw, pose)
}
fn full_fixture_source(screw: Screw, pose: [Point; 3]) -> (String, Screw) {
    let mut source_spec = screw.spec();
    source_spec.head = None;
    let hex = Screw::from_source_with_union(
        *screw.thread(),
        *screw.profile(),
        source_spec,
        crate::screw_recess::Recess::canonical_flat_hex(2., 1.).unwrap(),
    )
    .unwrap();
    let mut text =
        crate::step::write_profiled_screw(&hex, pose[0], pose[1], pose[2], 1e-7).unwrap();
    let shell = text
        .lines()
        .find(|line| line.contains("=CLOSED_SHELL("))
        .unwrap()
        .to_owned();
    let mut faces = references(&shell)[1..].to_vec();
    let head_index = faces.len() - hex.profile().segments().count() - 8;
    let top = entity(&text, &faces[head_index]).to_owned();
    let old_mouth = references(&top)[2].clone();
    let mut file = StepFile::new();
    file.data.clear();
    file.next = text
        .lines()
        .filter_map(|line| {
            line.strip_prefix('#')?
                .split_once('=')?
                .0
                .parse::<usize>()
                .ok()
        })
        .max()
        .unwrap()
        + 1;
    let crate::screw_recess::Recess::ConicalHexalobular(recess) = screw.recess() else {
        panic!("actual Torx cavity");
    };
    let (mouth, new_faces) = encode(
        &mut file,
        &Placement::new(pose[0], pose[1], pose[2]).unwrap(),
        recess,
        screw.head_top_z(),
        1e-7,
    )
    .unwrap();
    text = text.replacen(
        &top,
        &top.replace(&format!(",#{old_mouth})"), &format!(",{mouth})")),
        1,
    );
    faces.truncate(faces.len() - 7);
    let mut all_faces = faces.iter().map(|id| format!("#{id}")).collect::<Vec<_>>();
    all_faces.extend(new_faces.iter().map(ToString::to_string));
    let replacement = format!(
        "{}=CLOSED_SHELL('',({}));",
        shell.split_once('=').unwrap().0,
        all_faces.join(",")
    );
    text = text.replacen(&shell, &replacement, 1);
    let end = text.rfind("ENDSEC;").unwrap();
    text.insert_str(end, &file.data);
    // The fixture replaces an existing cavity. Remove its now-unreferenced
    // topological entities; production ownership checks must still see one
    // complete shell with no orphaned legacy cavity.
    let solid = text
        .lines()
        .find(|line| line.contains("=MANIFOLD_SOLID_BREP("))
        .unwrap();
    let mut pending = vec![references(solid)[0].clone()];
    let mut reachable = std::collections::BTreeSet::new();
    let entities: std::collections::HashMap<_, _> = text
        .lines()
        .filter_map(|line| {
            let id = line.strip_prefix('#')?.split_once('=')?.0;
            Some((id, line))
        })
        .collect();
    while let Some(id) = pending.pop() {
        if reachable.insert(id.clone()) {
            pending.extend(references(entities[id.as_str()]).into_iter().skip(1));
        }
    }
    text = text
        .lines()
        .filter(|line| {
            let topology = [
                "ADVANCED_FACE(",
                "FACE_OUTER_BOUND(",
                "FACE_BOUND(",
                "EDGE_LOOP(",
                "ORIENTED_EDGE(",
                "EDGE_CURVE(",
                "VERTEX_POINT(",
            ]
            .into_iter()
            .any(|kind| line.contains(&format!("={kind}")));
            !topology || reachable.contains(&references(line)[0])
        })
        .collect::<Vec<_>>()
        .join("\n");
    (text, screw)
}

#[test]
fn independent_torx_cavity_reader_retains_complete_actual_supports_and_curves() {
    let (text, source) = full_fixture();
    let decoded = crate::step_read::read_torx_cavity_fixture(
        &text,
        source.thread(),
        [[0.; 3], [1., 0., 0.], [0., 1., 0.]],
    )
    .unwrap();
    let crate::screw_recess::Recess::ConicalHexalobular(recess) = decoded.recess() else {
        panic!("physical Torx cavity");
    };
    let binding = recess.finite_binding().unwrap();
    assert_eq!(binding.carriers().len(), 12);
    assert!(binding.literal_cone_semi_angle_rad().is_some());
    assert!(binding.carriers().iter().all(|c| !c.is_empty()));
    assert!(decoded.spec().head.is_none());
    assert!(decoded.recess().volume_bounds_mm3()[0] > 0.);
    let (readback, imported) = crate::step_read::read_import_solid(&text, 1e-7).unwrap();
    assert!(readback.volume_mm3 > 0.);
    assert!(crate::step_read::profiled_screw_import_matches(
        &imported,
        &source,
        [[0.; 3], [1., 0., 0.], [0., 1., 0.]],
        1e-7
    ));
    for (tip_shift, transport) in [(0.01, 0.), (0., 2e-8)] {
        let source_drive = source.recess().hexalobular_geometry().unwrap();
        let mut cone = source_drive.cone();
        cone.tip_depth_mm += tip_shift;
        let drive = Drive::from_retained(
            *source_drive.arcs(),
            cone,
            crate::head_profile::TransportError {
                numeric_conversion_mm: 0.,
                arithmetic_mm: transport,
            },
        )
        .unwrap();
        let recess = crate::screw_recess::Recess::conical_hexalobular(drive, 0.1, 100_000).unwrap();
        let changed = Screw::new_with_recess(
            *source.thread(),
            source.profile().end()[0],
            source.spec().head.unwrap(),
            source.spec().shank_length_mm,
            0.1,
            recess,
        )
        .unwrap();
        assert!(!crate::step_read::profiled_screw_import_matches(
            &imported,
            &changed,
            [[0.; 3], [1., 0., 0.], [0., 1., 0.]],
            1e-7
        ));
    }
}

#[test]
fn independent_torx_readback_rejects_changed_supports_carriers_and_finite_uses() {
    let (text, source) = full_fixture();
    let intersection = text
        .lines()
        .find(|line| line.contains("=INTERSECTION_CURVE("))
        .unwrap();
    let ir = references(intersection);
    let spline = entity(&text, &ir[1]);
    let poles = references(spline);
    let interior = entity(&text, &poles[2]);
    let cylinder = entity(&text, &ir[2]);
    let cr = references(cylinder);
    let cylinder_axis = entity(&text, &cr[1]);
    let axis_refs = references(cylinder_axis);
    let cylinder_origin = entity(&text, &axis_refs[1]);
    let cylinder_normal = entity(&text, &axis_refs[2]);
    let cone = entity(&text, &ir[3]);
    let cone_axis = entity(&text, &references(cone)[1]);
    let cone_origin = entity(&text, &references(cone_axis)[1]);
    let seam_use = text
        .lines()
        .find(|line| line.contains("=SEAM_EDGE("))
        .unwrap();
    let seam_refs = references(seam_use);
    let seam = entity(&text, &references(entity(&text, &seam_refs[1]))[3]);
    let selected = entity(&text, &seam_refs[2]);
    let definition = entity(&text, &references(selected)[2]);
    let polyline = entity(&text, &references(definition)[1]);
    let uv = entity(&text, &references(polyline)[1]);
    let replace_entity =
        |line: &str, body: &str| format!("{}={body};", line.split_once('=').unwrap().0);
    let point_shift = |line: &str, index: usize, delta: f64| {
        let coordinates = line
            .split_once("('',(")
            .unwrap()
            .1
            .split_once(')')
            .unwrap()
            .0;
        let mut values = coordinates
            .split(',')
            .map(|v| v.parse::<f64>().unwrap())
            .collect::<Vec<_>>();
        values[index] += delta;
        replace_entity(
            line,
            &format!(
                "CARTESIAN_POINT('',({}))",
                values
                    .iter()
                    .map(|v| format!("{v:.16E}"))
                    .collect::<Vec<_>>()
                    .join(",")
            ),
        )
    };
    let mut cases = vec![
        edited(&text, interior, point_shift(interior, 2, 1e-4)),
        edited(
            &text,
            cylinder_origin,
            point_shift(cylinder_origin, 0, 1e-8),
        ),
        edited(&text, cone_origin, point_shift(cone_origin, 0, 1e-8)),
        edited(
            &text,
            cylinder_normal,
            replace_entity(cylinder_normal, "DIRECTION('',(1.E-8,0.,1.))"),
        ),
        edited(
            &text,
            intersection,
            intersection.replace(
                &format!("(#{},#{})", ir[2], ir[3]),
                &format!("(#{0},#{0})", ir[3]),
            ),
        ),
        edited(&text, spline, spline.replacen("('',3,", "('',2,", 1)),
        edited(&text, spline, spline.replace(".F.,.F.", ".T.,.F.")),
        edited(&text, spline, spline.replace(",3,3,", ",2,3,")),
        edited(&text, seam_use, seam_use.replace(",.T.,", ",.F.,")),
        edited(&text, seam, seam.replace(".PCURVE_S1.", ".PCURVE_S2.")),
        edited(&text, uv, point_shift(uv, 1, 1e-8)),
    ];
    let mut wall = None;
    for line in text.lines().filter(|line| line.contains("=ADVANCED_FACE(")) {
        if references(line).last() == Some(&ir[2]) {
            wall = Some(line);
            break;
        }
    }
    let wall = wall.unwrap();
    let wall_bound = entity(&text, &references(wall)[1]);
    let wall_loop = entity(&text, &references(wall_bound)[1]);
    let post_use = entity(&text, &references(wall_loop)[1]);
    let post_edge = entity(&text, &references(post_use)[1]);
    let post_line = entity(&text, &references(post_edge)[3]);
    let post_vector = entity(&text, &references(post_line)[2]);
    let post_direction = entity(&text, &references(post_vector)[1]);
    let lower_use = entity(&text, &references(wall_loop)[2]);
    let lower_edge = entity(&text, &references(lower_use)[1]);
    let lower_refs = references(lower_edge);
    let lower_spline_first = references(spline)[1].clone();
    let seam_line = entity(&text, &references(seam)[1]);
    let shell = text
        .lines()
        .find(|line| line.contains("=CLOSED_SHELL("))
        .unwrap();
    let faces = references(shell);
    let top_face = entity(
        &text,
        &faces[faces.len() - source.profile().segments().count() - 14],
    );
    // Actual support/parser failures must remain distinct from a broad file
    // tolerance; these are the literal geometric entities in a closed shell.
    cases.extend([
        edited(
            &text,
            spline,
            spline.replacen(
                &format!("(#{lower_spline_first},"),
                &format!("(#{},", poles[2]),
                1,
            ),
        ),
        edited(&text, seam_line, seam_line.replace("('',", "(0.,")),
        edited(
            &text,
            lower_edge,
            lower_edge.replacen(
                &format!(",#{},", lower_refs[1]),
                &format!(",#{},", lower_refs[2]),
                1,
            ),
        ),
        edited(&text, cylinder, cylinder.replace("('',", "(0.,")),
        edited(&text, post_line, post_line.replace("('',", "(0.,")),
        edited(&text, post_vector, post_vector.replace("('',", "(0.,")),
        edited(
            &text,
            post_direction,
            post_direction.replace("('',", "(0.,"),
        ),
        edited(
            &text,
            post_vector,
            replace_entity(
                post_vector,
                &format!("VECTOR('',#{},1.E13)", references(post_vector)[1]),
            ),
        ),
        edited(&text, selected, selected.replace("('',", "(0.,")),
        edited(
            &text,
            intersection,
            intersection.replace(".CURVE_3D.", ".PCURVE_S1."),
        ),
        edited(&text, spline, spline.replace("),(0.,", "),('',")),
        edited(&text, spline, spline.replace("),(0.,", "),(0.5,")),
        edited(&text, cone, cone.replace(",0.,", ",1.,")),
        edited(&text, post_edge, post_edge.replace(",.T.);", ",.F.);")),
    ]);
    let top_surface = entity(&text, references(top_face).last().unwrap());
    cases.push(edited(
        &text,
        top_surface,
        top_surface.replace("('',", "(0.,"),
    ));
    let flipped = if wall.ends_with(",.F.);") {
        wall.replace(",.F.);", ",.T.);")
    } else {
        wall.replace(",.T.);", ",.F.);")
    };
    let mouth_use = entity(&text, &references(wall_loop)[4]);
    let mouth_edge = entity(&text, &references(mouth_use)[1]);
    let mut invalid_contour = edited(&text, mouth_edge, mouth_edge.replace(",.T.);", ",.F.);"));
    invalid_contour = invalid_contour.replacen(wall, &flipped, 1);
    cases.push(invalid_contour);
    let evaluated = text
        .lines()
        .find(|line| line.contains("=EVALUATED_DEGENERATE_PCURVE("))
        .unwrap();
    let evaluated_refs = references(evaluated);
    let apex_vertex = text
        .lines()
        .find(|line| {
            line.contains("=VERTEX_POINT(") && references(line).get(1) == Some(&evaluated_refs[0])
        })
        .unwrap();
    let mut nonzero_reference_cone = edited(&text, cone, cone.replace(",0.,", ",1.,"));
    nonzero_reference_cone = nonzero_reference_cone.replacen(
        apex_vertex,
        &apex_vertex.replace(
            &format!("#{}", evaluated_refs[0]),
            &format!("#{}", evaluated_refs[3]),
        ),
        1,
    );
    cases.push(nonzero_reference_cone);
    cases.push(edited(
        &text,
        wall,
        if wall.ends_with(",.F.);") {
            wall.replace(",.F.);", ",.T.);")
        } else {
            wall.replace(",.T.);", ",.F.);")
        },
    ));
    let mouth_circle = text
        .lines()
        .filter(|line| line.contains("=CIRCLE("))
        .find(|line| {
            let axis = entity(&text, &references(line)[1]);
            let point = entity(&text, &references(axis)[1]);
            point == cylinder_origin
        })
        .unwrap();
    let mouth_axis = entity(&text, &references(mouth_circle)[1]);
    let mouth_edge_id = references(mouth_edge)[0].clone();
    let mouth_top_use = text
        .lines()
        .find(|line| {
            line.contains("=ORIENTED_EDGE(")
                && references(line).get(1) == Some(&mouth_edge_id)
                && line.ends_with(",.T.);")
        })
        .unwrap();
    let top_loop = text
        .lines()
        .find(|line| {
            line.contains("=EDGE_LOOP(") && references(line).contains(&references(mouth_top_use)[0])
        })
        .unwrap();
    let loop_refs = references(top_loop);
    let mut reordered = loop_refs[1..].to_vec();
    reordered.swap(0, 1);
    cases.push(edited(
        &text,
        top_loop,
        replace_entity(
            top_loop,
            &format!(
                "EDGE_LOOP('',({}))",
                reordered
                    .iter()
                    .map(|id| format!("#{id}"))
                    .collect::<Vec<_>>()
                    .join(",")
            ),
        ),
    ));
    let floor_use = text
        .lines()
        .find(|line| {
            line.contains("=ORIENTED_EDGE(")
                && references(line).get(1) == Some(&lower_refs[0])
                && line.ends_with(",.F.);")
        })
        .unwrap();
    cases.push(edited(
        &text,
        floor_use,
        floor_use.replace(",.F.);", ",.T.);"),
    ));
    // Duplicate a circle origin so the shift cannot also move the cylinder.
    let new_point = "#1900000=CARTESIAN_POINT('',(0.,0.,0.));";
    let shifted = point_shift(cylinder_origin, 2, 1e-8).replacen(
        cylinder_origin.split_once('=').unwrap().0,
        "#1900000",
        1,
    );
    let modified_axis = mouth_axis.replace(&format!(",#{},", axis_refs[1]), ",#1900000,");
    let mut circle_offset = edited(&text, mouth_axis, modified_axis);
    circle_offset.insert_str(
        circle_offset.rfind("ENDSEC;").unwrap(),
        &format!("{new_point}\n"),
    );
    circle_offset = circle_offset.replace(new_point, &shifted);
    cases.push(circle_offset);
    for (i, case) in cases.iter().enumerate() {
        assert!(case != &text, "fixture {i} must change actual data");
        assert!(
            crate::step_read::read_import_solid(case, 1e-7).is_err(),
            "changed support/carrier/use {i}"
        );
        assert!(
            crate::step_read::read_torx_cavity_fixture(
                case,
                source.thread(),
                [[0.; 3], [1., 0., 0.], [0., 1., 0.]]
            )
            .is_err(),
            "local cavity must reject malformed data {i}"
        );
    }
    assert!(crate::step_read::read_import_solid(&text, 1e-7).is_ok());
    let mouth_refs = references(mouth_edge);
    let changed_edge = replace_entity(
        mouth_edge,
        &format!(
            "EDGE_CURVE('',#{},#{},#{},.F.)",
            mouth_refs[2], mouth_refs[1], mouth_refs[3]
        ),
    );
    let mut equivalent = edited(&text, mouth_edge, changed_edge);
    equivalent = equivalent.replacen(mouth_top_use, &mouth_top_use.replace(",.T.);", ",.F.);"), 1);
    equivalent = equivalent.replacen(mouth_use, &mouth_use.replace(",.F.);", ",.T.);"), 1);
    let (_, imported) = crate::step_read::read_import_solid(&equivalent, 1e-7).unwrap();
    assert!(crate::step_read::profiled_screw_import_matches(
        &imported,
        &source,
        [[0.; 3], [1., 0., 0.], [0., 1., 0.]],
        1e-7
    ));
}

#[test]
fn independently_decoded_torx_heads_keep_original_hands_datums_and_proper_world_frames() {
    for left in [false, true] {
        for plain in [0., 1.] {
            let pose = [
                [105.4080000120003, -20.544000016000386, 59.26000001500036],
                [
                    0.8000001131000105,
                    0.5999998491999037,
                    -2.5133344449321634e-7,
                ],
                [
                    0.3599997486665377,
                    -0.48000018850005965,
                    -0.7999999999999605,
                ],
            ];
            let (_, source) = full_fixture_with(left, plain, pose);
            let text = crate::step::write_profiled_screw(&source, pose[0], pose[1], pose[2], 1e-7)
                .unwrap();
            let (_, decoded) = crate::step_read::read_import_solid(&text, 1e-7).unwrap();
            assert!(
                crate::step_read::profiled_screw_import_matches(&decoded, &source, pose, 1e-7),
                "left={left} plain={plain}"
            );
            let (physical, error) = decoded.profiled_screw_geometry().unwrap();
            assert!(error > 0.);
            assert_eq!(physical.thread().spec().left_hand, left);
            assert_eq!(physical.has_head_thread_overlap(), plain == 0.);
            assert!(physical.recess().hexalobular_geometry().is_some());
            assert!(crate::step::supports_profiled_screw_in_frame(
                physical,
                decoded.origin_mm,
                decoded.u,
                decoded.v,
                1e-7
            ));
            let roundtrip = crate::step::write_profiled_screw(
                physical,
                decoded.origin_mm,
                decoded.u,
                decoded.v,
                1e-7,
            )
            .unwrap();
            assert!(crate::step_read::read_import_solid(&roundtrip, 1e-7).is_ok());
        }
    }
}

#[test]
fn far_placed_metric_torx_import_preserves_its_literal_finite_reexport_data() {
    let pose = [[100000., -100000., 100000.], [0.6, 0.8, 0.], [0., 0., 1.]];
    let (_, source) = full_fixture_with(true, 0., pose);
    let text = crate::step::write_profiled_screw(&source, pose[0], pose[1], pose[2], 1e-7).unwrap();
    let (_, decoded) = crate::step_read::read_import_solid(&text, 1e-7).unwrap();
    assert!(crate::step_read::profiled_screw_import_matches(
        &decoded, &source, pose, 1e-7
    ));
    let (physical, _) = decoded.profiled_screw_geometry().unwrap();
    let crate::screw_recess::Recess::ConicalHexalobular(recess) = physical.recess() else {
        panic!("actual finite Torx");
    };
    let binding = recess.finite_binding().unwrap();
    assert!(binding.source().source_error_mm() > source.recess().source_error_mm());
    let mut file = StepFile::new();
    encode(
        &mut file,
        &Placement::new(decoded.origin_mm, decoded.u, decoded.v).unwrap(),
        recess,
        physical.head_top_z(),
        1e-7,
    )
    .unwrap();
    assert_eq!(file.data.matches("INTERSECTION_CURVE(").count(), 12);
    let reexport =
        crate::step::write_profiled_screw(physical, decoded.origin_mm, decoded.u, decoded.v, 1e-7)
            .unwrap();
    let (_, reimport) = crate::step_read::read_import_solid(&reexport, 1e-7).unwrap();
    assert!(crate::step_read::profiled_screw_import_matches(
        &reimport,
        physical,
        [decoded.origin_mm, decoded.u, decoded.v],
        1e-7
    ));
}

#[test]
fn public_torx_export_refuses_noncanonical_finite_binding_even_with_a_small_trim_tube() {
    let (_, source) = full_fixture();
    let text =
        crate::step::write_profiled_screw(&source, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7)
            .unwrap();
    let (_, imported) = crate::step_read::read_import_solid(&text, 1e-7).unwrap();
    let (physical, _) = imported.profiled_screw_geometry().unwrap();
    let crate::screw_recess::Recess::ConicalHexalobular(recess) = physical.recess() else {
        panic!("finite Torx");
    };
    let binding = recess.finite_binding().unwrap();
    let mut invalid = Vec::new();
    for kind in 0..3 {
        let mut vertices = *binding.vertices();
        let mut posts = *binding.post_lines();
        match kind {
            0 => posts[0].origin_mm[0] += 2e-10,
            1 => vertices.apex_mm[0] += 2e-10,
            _ => vertices.mouth_mm[0][2] += 2e-10,
        }
        let binding = FiniteBinding::new_with_posts(
            binding.source().clone(),
            vertices,
            posts,
            binding.carriers().clone(),
            1e-7,
        )
        .unwrap()
        .with_literal_cone_semi_angle(binding.literal_cone_semi_angle_rad().unwrap())
        .unwrap();
        let recess =
            crate::screw_recess::Recess::conical_hexalobular_from_binding(binding, 0.1, 100_000)
                .unwrap();
        assert!(recess.source_error_mm() < 1e-7 / 8.);
        let screw = Screw::from_source_with_union(
            *physical.thread(),
            *physical.profile(),
            physical.spec(),
            recess,
        )
        .unwrap();
        assert!(!crate::step::supports_profiled_screw(&screw, [0.; 3], 1e-7));
        assert_eq!(
            crate::step::write_profiled_screw(&screw, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7),
            Err(StepError::Unsupported)
        );
        invalid.push(screw);
    }
    let mut writer = crate::step::Writer::new();
    for screw in invalid {
        assert!(
            writer
                .append_profiled_screw(&screw, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7)
                .is_err()
        );
    }
    writer
        .append_profiled_screw(&source, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7)
        .unwrap();
    let written = writer.finish().unwrap();
    assert_eq!(
        crate::step_read::read_import_compound(&written, 1e-7)
            .unwrap()
            .len(),
        1
    );
}

#[test]
fn public_torx_gate_preserves_retained_negative_post_gauges_and_refuses_unqualified_families() {
    let (_, base) = head_fixture();
    let thin_flat = Screw::new(
        *base.thread(),
        crate::profiled_screw::Spec {
            head: base.spec().head.unwrap(),
            shank_length_mm: 1.,
            socket_across_flats_mm: 2.,
            socket_depth_mm: 5e-7,
            minimum_wall_mm: 0.1,
        },
    )
    .unwrap();
    assert_eq!(
        crate::step::write_profiled_screw(&thin_flat, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7),
        Err(StepError::Limit)
    );
    let make = |recess| {
        Screw::new_with_recess(
            *base.thread(),
            1.5,
            base.spec().head.unwrap(),
            1.,
            0.1,
            recess,
        )
        .unwrap()
    };
    let retained = make(crate::screw_recess::Recess::ConicalHexalobular(
        retained_payload(2e-10),
    ));
    assert!(crate::step::supports_profiled_screw(
        &retained, [0.; 3], 1e-7
    ));
    let text =
        crate::step::write_profiled_screw(&retained, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7)
            .unwrap();
    let (_, imported) = crate::step_read::read_import_solid(&text, 1e-7).unwrap();
    assert!(crate::step_read::profiled_screw_import_matches(
        &imported,
        &retained,
        [[0.; 3], [1., 0., 0.], [0., 1., 0.]],
        1e-7
    ));
    let source = source();
    let arbitrary = Drive::from_retained(
        *source.arcs(),
        source.cone(),
        crate::head_profile::TransportError::ZERO,
    )
    .unwrap();
    let arbitrary =
        make(crate::screw_recess::Recess::conical_hexalobular(arbitrary, 1e-6, 100_000).unwrap());
    assert_eq!(
        crate::step::write_profiled_screw(&arbitrary, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7),
        Err(StepError::Unsupported)
    );
    let shallow = Drive::from_authored_tangent(
        TangentSpec {
            major_diameter_mm: 2.8,
            minor_diameter_mm: 2.05,
            outer_lobe_radius_mm: 0.28,
        },
        Cone {
            tip_depth_mm: 1.,
            depth_per_radius: f64::MIN_POSITIVE,
        },
    )
    .unwrap();
    let shallow =
        make(crate::screw_recess::Recess::conical_hexalobular(shallow, 1e-6, 100_000).unwrap());
    assert_eq!(
        crate::step::write_profiled_screw(&shallow, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7),
        Err(StepError::Unsupported)
    );
    let wide = make(crate::screw_recess::Recess::ConicalHexalobular(
        retained_payload(1e-4),
    ));
    assert_eq!(
        crate::step::write_profiled_screw(&wide, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7),
        Err(StepError::Limit)
    );
    let thin = Drive::from_authored_tangent(
        TangentSpec {
            major_diameter_mm: 2.8,
            minor_diameter_mm: 2.05,
            outer_lobe_radius_mm: 0.28,
        },
        Cone {
            tip_depth_mm: 1.4 * 5e-5 + 1e-9,
            depth_per_radius: 5e-5,
        },
    )
    .unwrap();
    let thin = make(crate::screw_recess::Recess::conical_hexalobular(thin, 1e-6, 100_000).unwrap());
    // The raw-support phase has no tolerance argument: when placement rounds
    // both ends to indistinguishable coordinates, it cannot certify an edge
    // orientation. The public writer also enforces its earlier global budget.
    let unresolved = Placement::new([1e12; 3], [1., 0., 0.], [0., 1., 0.]).unwrap();
    assert_eq!(
        preflight::placement(&thin, &unresolved),
        Err(StepError::Unsupported)
    );
    assert_eq!(
        crate::step::write_profiled_screw(&thin, [1e12; 3], [1., 0., 0.], [0., 1., 0.], 1e-7),
        Err(StepError::Limit)
    );
}

fn entity<'a>(text: &'a str, id: &str) -> &'a str {
    text.lines()
        .find(|line| line.starts_with(&format!("#{id}=")))
        .unwrap()
}
fn references(line: &str) -> Vec<String> {
    line.split('#')
        .skip(1)
        .map(|s| s.chars().take_while(char::is_ascii_digit).collect())
        .collect()
}
fn edited(text: &str, line: &str, new: String) -> String {
    text.replacen(line, &new, 1)
}

#[test]
fn finite_apex_graph_rejects_wrong_parameter_points_seam_refs_and_orphaned_topology() {
    let (text, _) = head_fixture();
    let apex = text
        .lines()
        .find(|l| l.contains("EVALUATED_DEGENERATE_PCURVE("))
        .unwrap();
    let a = references(apex);
    let cone = entity(&text, &a[1]);
    let definition = entity(&text, &a[2]);
    let d = references(definition);
    let curve = entity(&text, &d[1]);
    let context = entity(&text, &d[2]);
    let c = references(curve);
    let point = entity(&text, &c[1]);
    let seam_edge = text.lines().find(|l| l.contains("SEAM_EDGE(")).unwrap();
    let e = references(seam_edge);
    let edge = entity(&text, &e[1]);
    let edge_refs = references(edge);
    let seam = entity(&text, &edge_refs[3]);
    let selected = entity(&text, &e[2]);
    let mut cases = vec![
        edited(&text, apex, apex.replace("('',", "(0.,")),
        edited(&text, cone, cone.replace(",0.,", ",1.,")),
        edited(
            &text,
            cone,
            format!("{},0.);", cone.rsplit_once(',').unwrap().0),
        ),
        edited(&text, definition, definition.replace("('',", "(0.,")),
        edited(
            &text,
            definition,
            definition.replace(&format!("(#{})", d[1]), &format!("(#{0},#{0})", d[1])),
        ),
        edited(&text, context, context.replace("CONTEXT(2)", "CONTEXT(3)")),
        edited(&text, curve, curve.replace("('',", "(0.,")),
        edited(
            &text,
            curve,
            curve.replace(&format!("(#{},#{})", c[1], c[2]), &format!("(#{})", c[1])),
        ),
        edited(&text, point, point.replace("('',", "(0.,")),
        edited(
            &text,
            point,
            format!(
                "{}=CARTESIAN_POINT('',(0.,0.,0.));",
                point.split_once('=').unwrap().0
            ),
        ),
        edited(
            &text,
            point,
            format!(
                "{}=CARTESIAN_POINT('',(1.E100,0.));",
                point.split_once('=').unwrap().0
            ),
        ),
        edited(
            &text,
            point,
            format!(
                "{}=CARTESIAN_POINT('',(0.,0.01));",
                point.split_once('=').unwrap().0
            ),
        ),
        edited(&text, seam_edge, seam_edge.replace("('',", "(0.,")),
        edited(
            &text,
            seam_edge,
            seam_edge.replace(&format!(",#{});", e[2]), &format!(",#{});", a[2])),
        ),
        edited(&text, selected, selected.replace("=PCURVE(", "=PLANE(")),
        edited(&text, seam, seam.replace("('',", "(0.,")),
        edited(&text, seam, seam.replace(".PCURVE_S1.", ".CURVE_3D.")),
    ];
    let vertex = text.lines().find(|l| l.contains("VERTEX_POINT(")).unwrap();
    cases.push(edited(&text, vertex, vertex.replace("('',", "(0.,")));
    let loop_ = text
        .lines()
        .find(|line| line.contains("EDGE_LOOP("))
        .unwrap();
    cases.push(edited(&text, loop_, loop_.replace("('',", "(0.,")));
    let bound = text
        .lines()
        .find(|line| line.contains("FACE_OUTER_BOUND("))
        .unwrap();
    cases.push(edited(&text, bound, bound.replace(".T.", ".F.")));
    let first_point = text
        .lines()
        .find(|l| l.contains("=CARTESIAN_POINT("))
        .unwrap();
    let point_id = references(first_point)[0].clone();
    cases.push(edited(
        &text,
        apex,
        apex.replace(&format!(",#{});", a[3]), &format!(",#{point_id});")),
    ));
    cases.push(edited(&text, seam_edge, seam_edge.replace(".T.", ".F.")));
    let floor_face = text
        .lines()
        .find(|line| line.contains("ADVANCED_FACE(") && line.contains(&format!(",#{},", a[1])))
        .unwrap();
    let wrong_loop = edited(&text, seam_edge, seam_edge.replace(".T.", ".F."));
    cases.push(edited(
        &wrong_loop,
        floor_face,
        floor_face.replace(".F.", ".T."),
    ));
    let orphan = format!("#999999=SEAM_EDGE('',*,*,#{},.T.,#{});\n", e[1], e[2]);
    cases.push(text.replace("ENDSEC;\nEND-ISO", &format!("{orphan}ENDSEC;\nEND-ISO")));
    for (i, bad) in cases.iter().enumerate() {
        assert!(
            crate::step_read::verify_finite_graph_fixture(bad).is_err(),
            "malformed finite graph case {i}"
        );
    }
    assert!(crate::step_read::verify_finite_graph_fixture(&text).is_ok());
}

#[test]
fn torx_encoder_refuses_unresolvable_precision_and_collapsed_world_geometry() {
    let source = source();
    let recess = payload(source.clone());
    let snapshot = source.clone();
    let p = Placement::new([0.; 3], [1., 0., 0.], [0., 1., 0.]).unwrap();
    for tolerance in [0., -1., f64::NAN, f64::INFINITY] {
        assert!(matches!(
            encode(&mut StepFile::new(), &p, &recess, -2.4, tolerance),
            Err(StepError::InvalidGeometry)
        ));
    }
    assert!(matches!(
        encode(&mut StepFile::new(), &p, &recess, f64::NAN, 1e-7),
        Err(StepError::InvalidGeometry)
    ));
    assert!(matches!(
        encode(&mut StepFile::new(), &p, &recess, -2.4, 1e-20),
        Err(StepError::Unsupported)
    ));
    let shallow = Drive::from_authored_tangent(
        TangentSpec {
            major_diameter_mm: 2.8,
            minor_diameter_mm: 2.05,
            outer_lobe_radius_mm: 0.28,
        },
        Cone {
            tip_depth_mm: 1.,
            depth_per_radius: f64::MIN_POSITIVE,
        },
    )
    .unwrap();
    assert!(matches!(
        encode(&mut StepFile::new(), &p, &payload(shallow), -2.4, 1e-7),
        Err(StepError::Unsupported)
    ));
    let collapsed = Placement::new([1e300; 3], [1., 0., 0.], [0., 1., 0.]).unwrap();
    assert!(matches!(
        encode(&mut StepFile::new(), &collapsed, &recess, -2.4, 1e-7),
        Err(StepError::InvalidGeometry)
    ));
    let skewed = Placement {
        origin: [0.; 3],
        u: [0.5, 0., 0.],
        v: [0., 1., 0.],
        normal: [0., 0., 1.],
    };
    assert!(matches!(
        encode(&mut StepFile::new(), &skewed, &recess, -2.4, 1e-7),
        Err(StepError::InvalidGeometry)
    ));
    for error in [
        Error::Invalid,
        Error::Unresolved,
        Error::Limit,
        Error::Cancelled,
    ] {
        let actual = source_error(error);
        assert_eq!(
            actual,
            match error {
                Error::Invalid => StepError::InvalidGeometry,
                Error::Unresolved => StepError::Unsupported,
                Error::Limit | Error::Cancelled => StepError::Limit,
            }
        );
    }
    assert_eq!(source, snapshot);
    assert!(encode(&mut StepFile::new(), &p, &recess, -2.4, 1e-7).is_ok());
}
