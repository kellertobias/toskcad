use super::*;
use crate::{
    head_profile::{Head, Profile},
    helical_thread,
    phillips_drive::{Cavity, Cutter},
    profiled_screw::SourceSpec,
};
fn fixture(index: usize, left: bool, nominal: bool, starts: u32) -> PhillipsScrew {
    let (s, d) = crate::phillips_drive::test_fixture(index);
    let cavity = Cavity::new(Cutter::new(s, d).unwrap(), 1e-5, 100_000).unwrap();
    let (rs, pitch, r, h, crown, land, neck) = [
        (1., 0.4, 2., 1.6, 3.2, 1.25, 0.1),
        (1.5, 0.5, 2.8, 2.4, 3., 2.1, 0.1),
        (1.5, 0.5, 2.75, 1.65, 0., 0., 0.1),
        (2., 0.7, 4., 3.1, 6.5, 3., 0.2),
        (3., 1., 6., 4.6, 10., 4.2, 0.25),
        (4., 1.25, 8., 6., 13., 5.4, 0.4),
    ][index];
    let head = if index == 2 {
        Head::Countersunk {
            radius_mm: r,
            height_mm: h,
            included_angle_deg: 90.,
            edge_chamfer_mm: 0.05,
            neck_radius_mm: neck,
        }
    } else {
        Head::Crowned {
            radius_mm: r,
            height_mm: h,
            crown_radius_mm: crown,
            top_land_radius_mm: land,
            neck_radius_mm: neck,
        }
    };
    let thread = Thread::new(helical_thread::Spec {
        major_radius_mm: if nominal { rs } else { rs - 0.03 },
        pitch_mm: pitch,
        length_mm: 2.,
        depth_mm: 0.6 * pitch,
        crest_width_mm: 0.125 * pitch,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts,
        left_hand: left,
        root_radius_mm: None,
    })
    .unwrap();
    PhillipsScrew::new(
        thread,
        Profile::new_for_thread_union(head, rs, 0.).unwrap(),
        SourceSpec {
            head: Some(head),
            shank_length_mm: 0.,
            minimum_wall_mm: 1e-5,
        },
        cavity,
        100_000,
        &|| false,
    )
    .unwrap()
}
#[test]
fn actual_metric_phi_heads_keep_original_thread_and_closed_raw_opposite_uses() {
    for index in 0..6 {
        for left in [false, true] {
            let source = fixture(index, left, false, 1);
            let original = *source.thread();
            let (origin, u, v) = if left {
                ([105., -20., 60.], [0.8, 0.6, 0.], [0.36, -0.48, -0.8])
            } else {
                ([0.; 3], [1., 0., 0.], [0., 1., 0.])
            };
            let text = write_phillips_screw(&source, origin, u, v, 1e-7).unwrap();
            crate::step_read::verify_finite_graph_fixture(&text).unwrap();
            let proof = crate::step_read::verify_phillips_screw_source_known(
                &text, &source, origin, u, v, 1e-7,
            )
            .unwrap();
            assert_eq!(proof.body_face_count(), source.faces().len());
            assert_eq!(proof.material_cavity_face_count(), 28);
            assert_eq!(proof.reference_cavity().vertex_count(), 49);
            assert_eq!(text.matches("MANIFOLD_SOLID_BREP(").count(), 1);
            assert_eq!(text.matches("ADVANCED_FACE(").count(), source.faces().len());
            assert_eq!(text.matches("INTERSECTION_CURVE(").count(), 32);
            assert!(text.contains("TOROIDAL_SURFACE(") && text.contains("CONICAL_SURFACE("));
            if index == 0 {
                assert!(text.contains("ELLIPSE("));
            }
            assert_eq!(source.thread(), &original);
        }
    }
    for left in [false, true] {
        let source = fixture(1, left, true, 2);
        let text =
            write_phillips_screw(&source, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7).unwrap();
        crate::step_read::verify_finite_graph_fixture(&text).unwrap();
        let proof = crate::step_read::verify_phillips_screw_source_known(
            &text,
            &source,
            [0.; 3],
            [1., 0., 0.],
            [0., 1., 0.],
            1e-7,
        )
        .unwrap();
        assert_eq!(proof.body_face_count(), source.faces().len() + 1);
        assert_eq!(
            text.matches("ADVANCED_FACE(").count(),
            source.faces().len() + 1
        );
    }
}
#[test]
fn callback_cancellation_rolls_back_real_thread_and_cavity_entities_before_retry() {
    use std::cell::Cell;
    let source = fixture(1, false, false, 1);
    let [u, v] = [[1., 0., 0.], [0., 1., 0.]];
    let origin = [10., 20., 30.];
    let count = Cell::new(0);
    write_phillips_screw_with_interrupt(
        &source,
        origin,
        u,
        v,
        1e-7,
        Some(&|| {
            count.set(count.get() + 1);
            false
        }),
    )
    .unwrap();
    let total = count.get();
    assert!(total > 100);
    let mut writer = Writer::new();
    writer
        .append_phillips_screw(&source, [0.; 3], u, v, 1e-7)
        .unwrap();
    writer
        .append_phillips_screw(&source, origin, u, v, 1e-7)
        .unwrap();
    let expected = writer.finish().unwrap();
    for at in [0, 1, total / 4, total / 2, total * 3 / 4, total - 1] {
        let mut writer = Writer::new();
        writer
            .append_phillips_screw(&source, [0.; 3], u, v, 1e-7)
            .unwrap();
        let count = Cell::new(0);
        assert_eq!(
            writer.append_phillips_screw_with_interrupt(
                &source,
                origin,
                u,
                v,
                1e-7,
                Some(&|| {
                    let n = count.get();
                    count.set(n + 1);
                    n >= at
                })
            ),
            Err(Error::Cancelled)
        );
        writer
            .append_phillips_screw(&source, origin, u, v, 1e-7)
            .unwrap();
        assert_eq!(writer.finish().unwrap(), expected);
    }
}
fn changed_point(text: &str, id: &str) -> String {
    let line = text
        .lines()
        .find(|l| l.starts_with(&format!("{id}=CARTESIAN_POINT(")))
        .unwrap();
    let coords = line
        .split("('',(")
        .nth(1)
        .unwrap()
        .split("))")
        .next()
        .unwrap();
    let mut values = coords
        .split(',')
        .map(|v| v.parse::<f64>().unwrap())
        .collect::<Vec<_>>();
    values[0] += 1e-4;
    text.replacen(
        line,
        &format!(
            "{id}=CARTESIAN_POINT('',({:.16E},{:.16E},{:.16E}));",
            values[0], values[1], values[2]
        ),
        1,
    )
}
#[test]
fn real_cavity_supports_every_cubic_and_material_face_senses_bind_to_the_source() {
    let source = fixture(0, false, false, 1);
    let origin = [0.; 3];
    let u = [1., 0., 0.];
    let v = [0., 1., 0.];
    let text = write_phillips_screw(&source, origin, u, v, 1e-7).unwrap();
    let check = |s: &str| {
        crate::step_read::verify_phillips_screw_source_known(s, &source, origin, u, v, 1e-7)
    };
    check(&text).unwrap();
    let cubics: Vec<_> = text
        .lines()
        .filter(|l| l.contains("=B_SPLINE_CURVE_WITH_KNOTS('',3,"))
        .collect();
    assert_eq!(cubics.len(), 32);
    for cubic in cubics {
        let list = cubic
            .split("('',3,(")
            .nth(1)
            .unwrap()
            .split(')')
            .next()
            .unwrap();
        let id = list.split(',').nth(1).unwrap();
        assert!(
            check(&changed_point(&text, id)).is_err(),
            "actual cubic {cubic}"
        );
    }
    let torus = text
        .lines()
        .rfind(|l| l.contains("=TOROIDAL_SURFACE("))
        .unwrap();
    let fields: Vec<_> = torus.split(',').collect();
    let radius = fields[2].parse::<f64>().unwrap();
    let replacement = format!(
        "{},{},{:.16E},{}",
        fields[0],
        fields[1],
        radius + 1e-4,
        fields[3]
    );
    assert!(check(&text.replacen(torus, &replacement, 1)).is_err());
    let face = text
        .lines()
        .rfind(|l| l.contains("=ADVANCED_FACE("))
        .unwrap();
    let flipped = if face.ends_with(".T.);") {
        face.replace(".T.);", ".F.);")
    } else {
        face.replace(".F.);", ".T.);")
    };
    let changed = text.replacen(face, &flipped, 1);
    crate::step_read::verify_finite_graph_fixture(&changed).unwrap();
    assert!(
        check(&changed).is_err(),
        "face sense cannot repair raw incidence"
    );
    assert!(
        crate::step_read::verify_phillips_screw_source_known(
            &text,
            &source,
            [1., 0., 0.],
            u,
            v,
            1e-7
        )
        .is_err()
    );
    for tolerance in [0., f64::NAN, 1e-20] {
        assert!(
            crate::step_read::verify_phillips_screw_source_known(
                &text, &source, origin, u, v, tolerance
            )
            .is_err()
        );
    }
}
#[test]
fn a_different_valid_actual_cavity_cannot_match_through_its_shared_outer_head() {
    let source = fixture(1, false, false, 1);
    let text = write_phillips_screw(&source, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7).unwrap();
    let mut spec = source.cavity().source().spec();
    spec.mouth_diameter_mm += 0.1;
    let cutter = Cutter::new(spec, source.cavity().source().domain()).unwrap();
    let changed = PhillipsScrew::new(
        *source.thread(),
        *source.profile(),
        source.spec(),
        Cavity::new(cutter, 1e-5, 100_000).unwrap(),
        100_000,
        &|| false,
    )
    .unwrap();
    assert!(
        crate::step_read::verify_phillips_screw_source_known(
            &text,
            &changed,
            [0.; 3],
            [1., 0., 0.],
            [0., 1., 0.],
            1e-7
        )
        .is_err()
    );
}

#[test]
fn public_writer_refuses_invalid_precision_and_preserves_error_taxonomy() {
    use crate::phillips_drive::Error as SourceError;
    for (input, expected) in [
        (SourceError::Invalid, StepError::InvalidGeometry),
        (SourceError::Unresolved, StepError::Limit),
        (SourceError::Limit, StepError::Limit),
        (SourceError::Cancelled, StepError::Limit),
    ] {
        assert_eq!(source_error(input), expected);
    }
    for (input, expected) in [
        (Error::InvalidGeometry, StepError::InvalidGeometry),
        (Error::InvalidPlacement, StepError::InvalidPlacement),
        (Error::Unsupported, StepError::Unsupported),
        (Error::Limit, StepError::Limit),
        (Error::Cancelled, StepError::Limit),
    ] {
        assert_eq!(plan_error(input), expected);
    }
    let source = fixture(1, false, false, 1);
    for tolerance in [0., f64::NAN, 1e-30] {
        assert!(
            write_phillips_screw(&source, [0.; 3], [1., 0., 0.], [0., 1., 0.], tolerance).is_err()
        );
    }
}
#[test]
fn cavity_primitive_refuses_an_extra_actual_bound_before_support_evaluation() {
    let source = fixture(1, false, false, 1);
    let text = write_phillips_screw(&source, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7).unwrap();
    crate::step_read::verify_phillips_cavity_fixture(&text, &source).unwrap();
    let face = text
        .lines()
        .rev()
        .filter(|l| l.contains("=ADVANCED_FACE("))
        .nth(1)
        .unwrap();
    let start = face.find("('',(").unwrap() + 5;
    let end = start + face[start..].find(')').unwrap();
    let bound = &face[start..end];
    let changed_face = format!("{},{bound}{}", &face[..end], &face[end..]);
    // Keep the parser and actual entity references: this guard is tested at the
    // cavity primitive, while the whole-shell public gate rejects repeated uses.
    let changed = text.replacen(face, &changed_face, 1);
    assert_eq!(
        crate::step_read::verify_phillips_cavity_fixture(&changed, &source),
        Err(crate::step_read::ReadError::Unsupported)
    );
    assert!(
        crate::step_read::verify_phillips_screw_source_known(
            &changed,
            &source,
            [0.; 3],
            [1., 0., 0.],
            [0., 1., 0.],
            1e-7
        )
        .is_err()
    );
}
