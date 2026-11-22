use super::*;
use crate::{
    head_profile::{Curve, Head, Profile},
    helical_thread::{Spec as ThreadSpec, Thread},
    phillips_drive::{Cavity, Cutter, PhillipsScrew},
    profiled_screw::SourceSpec,
};
pub(in crate::step_read::phillips_cutter::physical::body) fn additional_text(head: Head) -> String {
    let (s, d) = crate::phillips_drive::test_fixture(0);
    let cavity = Cavity::new(Cutter::new(s, d).unwrap(), 1e-5, 100_000).unwrap();
    let thread = Thread::new(ThreadSpec {
        major_radius_mm: 1.47,
        pitch_mm: 0.4,
        length_mm: 2.,
        depth_mm: 0.24,
        crest_width_mm: 0.05,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts: 1,
        left_hand: false,
        root_radius_mm: None,
    })
    .unwrap();
    let screw = PhillipsScrew::new(
        thread,
        Profile::new_for_thread_union(head, 1.5, 0.).unwrap(),
        SourceSpec {
            head: Some(head),
            shank_length_mm: 0.,
            minimum_wall_mm: 1e-5,
        },
        cavity,
        100_000,
        &|| false,
    )
    .unwrap();
    crate::step::write_phillips_screw(&screw, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7).unwrap()
}
fn text() -> String {
    super::super::super::branches::representative::scaled_text(1., false)
}
fn chart(profile: &Profile) -> CommonMouthChart {
    CommonMouthChart {
        columns: [
            [0., 0., profile.top()[1]],
            [1., 0., 0.],
            [0., 1., 0.],
            [0., 0., 1.],
        ],
        outer_top_z: profile.top()[1],
    }
}
#[test]
fn actual_raw_support_parser_and_cross_association_refusals_use_real_file_entities() {
    let text = text();
    let body = read_phillips_body(&text, 1e-7, 1_000_000).unwrap();
    let shell = Shell::new(&body.entities, body.tolerance, true);
    let profile = &body.profile;
    let source = profile
        .closed_repaired_meridian(Default::default(), None)
        .unwrap();
    let ids = body.profile_face_entities();
    let mut sphere = None;
    let mut cylinder = None;
    for (i, (&face, s)) in ids.iter().zip(profile.segments()).enumerate() {
        let (surface, _, _) = helical::face(&shell, face).unwrap();
        let raw = supports::read(&shell, surface).unwrap();
        let (k, bound) = supports::distance(
            &raw,
            chart(profile),
            s,
            source.segment_certificate(i).unwrap().support_transfer_mm,
        )
        .unwrap();
        assert!(bound > 0. && bound < 1e-7);
        if k == RawHeadSupportKind::Sphere {
            sphere = Some((surface, s));
        }
        if k == RawHeadSupportKind::Cylinder {
            cylinder = Some((surface, s));
        }
    }
    // A real polynomial thread surface is outside this ANALYTIC head decoder.
    let (surface, _, _) = helical::face(&shell, body.tensor_face_entities()[0]).unwrap();
    assert!(matches!(
        supports::read(&shell, surface),
        Err(ReadError::Unsupported)
    ));
    let (sphere_id, sphere_segment) = sphere.unwrap();
    let (cylinder_id, cylinder_segment) = cylinder.unwrap();
    let sphere_raw = supports::read(&shell, sphere_id).unwrap();
    assert!(matches!(
        supports::distance(&sphere_raw, chart(profile), cylinder_segment, 0.),
        Err(ReadError::Unsupported)
    ));
    assert!(
        supports::top_cap_distance(
            &sphere_raw,
            chart(profile),
            profile.top()[0],
            profile.top()[1]
        )
        .is_err()
    );
    let cylinder_raw = supports::read(&shell, cylinder_id).unwrap();
    assert!(matches!(
        supports::distance(&cylinder_raw, chart(profile), sphere_segment, 0.),
        Err(ReadError::Unsupported)
    ));
    // Properly rotate a real proposed cylinder relative to the unchanged raw
    // target: a valid parameter-basis fallback must not create a close proof.
    let mut turned = chart(profile);
    turned.columns[1] = [0., 0., 1.];
    turned.columns[2] = [0., 1., 0.];
    turned.columns[3] = [-1., 0., 0.];
    let (_, bound) = supports::distance(&cylinder_raw, turned, cylinder_segment, 0.).unwrap();
    assert!(bound > 0.1);
    // A genuine affine scale is outside the near-rigid Gram interval method;
    // it must refuse instead of taking sqrt of a nonpositive lower eigenbound.
    let mut scaled = chart(profile);
    scaled.columns[1] = [2., 0., 0.];
    assert!(supports::distance(&sphere_raw, scaled, sphere_segment, 0.).is_err());
    assert_eq!(body.literal_input(), text);
}
#[test]
fn actual_radius_and_axis_mutations_cannot_be_hidden_by_support_comparison() {
    let text = text();
    let body = read_phillips_body(&text, 1e-7, 1_000_000).unwrap();
    let shell = Shell::new(&body.entities, body.tolerance, true);
    let (index, s) = body
        .profile
        .segments()
        .enumerate()
        .find(|(_, s)| matches!(s.curve, Curve::Arc { .. }))
        .unwrap();
    let (surface, _, _) = helical::face(&shell, body.profile_faces[index]).unwrap();
    let entity = shell.simple(surface).unwrap();
    let words = part21::numeric_spellings(entity).unwrap();
    let record = text
        .lines()
        .find(|l| l.starts_with(&format!("#{surface}=")))
        .unwrap();
    for new_radius in [-1., s.start[0] + 10.] {
        let changed_record = record.replacen(words[0], &format!("{new_radius:.17E}"), 1);
        assert_ne!(record, changed_record);
        let changed = text.replacen(record, &changed_record, 1);
        let entities = part21::parse(&changed).unwrap();
        let input = Shell::new(&entities, 1e-7, true);
        if new_radius < 0. {
            assert!(supports::read(&input, surface).is_err());
        } else {
            let raw = supports::read(&input, surface).unwrap();
            let (_, error) = supports::distance(&raw, chart(&body.profile), s, 0.).unwrap();
            assert!(error > 1.);
        }
    }
    let axis = reference(&entity.args[1]).unwrap();
    let normal = reference(&body.entities[&axis][0].args[2]).unwrap();
    let record = text
        .lines()
        .find(|l| l.starts_with(&format!("#{normal}=")))
        .unwrap();
    let changed = text.replacen(record, &format!("#{normal}=DIRECTION('',(0.,0.,0.));"), 1);
    let entities = part21::parse(&changed).unwrap();
    let input = Shell::new(&entities, 1e-7, true);
    assert!(supports::read(&input, surface).is_err());
    assert_eq!(body.literal_input(), text);
}
#[test]
fn actual_cone_sheet_and_angular_domain_refusals_preserve_original_data() {
    let text = additional_text(Head::Countersunk {
        radius_mm: 3.,
        height_mm: 1.65,
        included_angle_deg: 90.,
        edge_chamfer_mm: 0.1,
        neck_radius_mm: 0.1,
    });
    let body = read_phillips_body(&text, 1e-7, 1_000_000).unwrap();
    let shell = Shell::new(&body.entities, body.tolerance, true);
    let (mut found, mut shifted) = (false, false);
    for (&face, s) in body.profile_faces.iter().zip(body.profile.segments()) {
        let (surface, _, _) = helical::face(&shell, face).unwrap();
        let raw = supports::read(&shell, surface).unwrap();
        if !matches!(raw.kind, supports::Kind::Cone(..)) {
            continue;
        }
        found = true;
        let mut chart = chart(&body.profile);
        chart.columns[0][2] += 100.;
        if supports::distance(&raw, chart, s, 0.).is_err() {
            shifted = true;
        }
        let record = text
            .lines()
            .find(|l| l.starts_with(&format!("#{surface}=")))
            .unwrap();
        let words = part21::numeric_spellings(shell.simple(surface).unwrap()).unwrap();
        let changed_record = record.replacen(words[1], "2.00000000000000000E+00", 1);
        assert_ne!(record, changed_record);
        let changed = text.replacen(record, &changed_record, 1);
        let entities = part21::parse(&changed).unwrap();
        let input = Shell::new(&entities, 1e-7, true);
        assert!(supports::read(&input, surface).is_err());
    }
    assert!(found && shifted);
    assert_eq!(body.literal_input(), text);
}
