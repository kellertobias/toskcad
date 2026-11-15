use super::*;
use crate::radial_tensor::{BoundaryOwnership, RadialBoundary};

fn material(owned: BoundaryOwnership) -> RadialBoundary {
    owned
        .bind_surfaces(1e-7, None)
        .unwrap()
        .certify_axial(1e-7, 1e-6, None)
        .unwrap()
        .certify_caps(None)
        .unwrap()
        .certify_ordering(None)
        .unwrap()
        .certify_material(1e-7, None)
        .unwrap()
}

#[test]
fn mixed_reader_preserves_analytic_receipts_and_actual_finite_source_ownership() {
    let body = material(
        BoundaryOwnership::new(crate::radial_tensor::tests::cam::wedge(&[1, 2]), 1e-7, None)
            .unwrap(),
    );
    let posed = body
        .transformed_with_interrupt(
            [[25., 10., 5.], [0.8, 0.6, 0.], [-0.36, 0.48, 0.8]],
            1e-7,
            1e-6,
            None,
        )
        .unwrap();
    let mut writer = crate::step::Writer::new();
    writer.append_radial_boundary(&posed, 1e-7).unwrap();
    writer
        .append_box(
            [2, 3, 4].map(|v| crate::Rational::new(v, 1).unwrap()),
            [-20., 0., 0.],
            [1., 0., 0.],
            [0., 1., 0.],
        )
        .unwrap();
    let text = writer.finish().unwrap();
    let members = read_source_compound(&text, 1e-7).unwrap();
    assert_eq!(members.len(), 2);
    let mut retained = 0;
    let mut analytic = 0;
    for member in members {
        match member {
            SourceImport::LiteralNpt(_) => {
                panic!("this actual fixture contains only radial and analytic members")
            }
            SourceImport::Radial(owned) => {
                retained += 1;
                let restored = material(*owned);
                assert_eq!(
                    crate::step::write_radial_boundary(&restored, 1e-7).unwrap(),
                    crate::step::write_radial_boundary(&posed, 1e-7).unwrap()
                );
                assert_eq!(restored.input().faces.len(), posed.input().faces.len());
            }
            SourceImport::Analytic(value) => {
                analytic += 1;
                assert_eq!(value.0.volume_mm3, 24.);
                assert_eq!(value.0.solid_count, 1);
            }
        }
    }
    assert_eq!((retained, analytic), (1, 1));
    assert!(read_source_compound("not STEP", 1e-7).is_err());
    assert!(read_source_compound(&text, 0.).is_err());
    let malformed = text.replacen(
        "B_SPLINE_SURFACE_WITH_KNOTS('',9,9",
        "B_SPLINE_SURFACE_WITH_KNOTS('',9,8",
        1,
    );
    assert_ne!(malformed, text);
    assert!(read_source_compound(&malformed, 1e-7).is_err());
}

#[test]
fn mixed_reader_keeps_existing_void_shell_import() {
    use crate::orthogonal::{Box3, Operation, boolean};
    let solid = boolean(
        &[Box3 {
            min: [0.; 3],
            max: [15., 10., 8.],
        }],
        &[Box3 {
            min: [2., 3., 1.],
            max: [7., 7., 4.],
        }],
        Operation::Subtract,
        1e-7,
    )
    .unwrap()
    .unwrap();
    let boundary = crate::orthogonal_boundary::build(solid.quads()).unwrap();
    let text = crate::step::write_orthogonal(&boundary, |p| p).unwrap();
    let mut members = read_source_compound(&text, 1e-7).unwrap();
    let SourceImport::Analytic(member) = members.pop().unwrap() else {
        panic!()
    };
    assert_eq!(member.0.volume_mm3, 1140.);
    assert_eq!(member.0.face_count, 12);
}

#[test]
fn retained_member_verification_integrates_the_decoded_source_and_refuses_changed_values() {
    let body = material(
        BoundaryOwnership::new(crate::radial_tensor::tests::cam::wedge(&[1, 2]), 1e-7, None)
            .unwrap(),
    );
    let input = crate::step::write_radial_boundary(&body, 1e-7).unwrap();
    let expected = [ExpectedStepMember::RetainedRadial {
        boundary: &body,
        minimum_feature_mm: 1e-6,
    }];
    assert!(format!("{:?}", expected[0]).contains("RadialBoundary"));
    let mut verified = verify_authored_compound(&input, &expected, 1e-7).unwrap();
    let VerifiedStepMember::RetainedRadial {
        properties,
        bounds_mm,
    } = verified.pop().unwrap()
    else {
        panic!()
    };
    assert!(properties.volume_mm3[0] > 0. && properties.volume_mm3[1] >= properties.volume_mm3[0]);
    assert_eq!(properties.faces.len(), body.input().faces.len());
    assert!(bounds_mm[0][0] < bounds_mm[1][0]);
    let moved = body
        .transformed_with_interrupt([[1., 0., 0.], [1., 0., 0.], [0., 1., 0.]], 1e-7, 1e-6, None)
        .unwrap();
    assert!(
        verify_authored_compound(
            &input,
            &[ExpectedStepMember::RetainedRadial {
                boundary: &moved,
                minimum_feature_mm: 1e-6
            }],
            1e-7
        )
        .is_err()
    );
    assert!(
        verify_authored_compound(
            &input,
            &[ExpectedStepMember::RetainedRadial {
                boundary: &body,
                minimum_feature_mm: 0.
            }],
            1e-7
        )
        .is_err()
    );
    let tensor = body
        .input()
        .surfaces
        .iter()
        .find_map(|s| {
            if let crate::radial_tensor::Surface::Tensor(t) = s {
                Some(t)
            } else {
                None
            }
        })
        .unwrap();
    let p = tensor.poles()[tensor.dimensions()[1] + 1];
    let coordinates = |p: [f64; 3]| format!("({:.16E},{:.16E},{:.16E})", p[0], p[1], p[2]);
    let mut q = p;
    q[2] += 1000.;
    let changed = input.replacen(&coordinates(p), &coordinates(q), 1);
    assert_ne!(changed, input);
    assert!(verify_authored_compound(&changed, &expected, 1e-7).is_err());
}

#[test]
fn retained_volume_only_verification_keeps_geometry_identity_and_interval_authority() {
    let body = material(
        BoundaryOwnership::new(crate::radial_tensor::tests::cam::wedge(&[1, 2]), 1e-7, None)
            .unwrap(),
    );
    let input = crate::step::write_radial_boundary(&body, 1e-7).unwrap();
    let expected = [ExpectedStepMember::RetainedRadialVolume {
        boundary: &body,
        minimum_feature_mm: 1e-6,
    }];
    assert!(format!("{:?}", expected[0]).contains("RadialBoundary"));
    let mut verified = verify_authored_compound(&input, &expected, 1e-7).unwrap();
    let VerifiedStepMember::RetainedRadialVolume { volume, bounds_mm } = verified.pop().unwrap()
    else {
        panic!()
    };
    assert!(volume.volume_mm3[0] > 0. && volume.volume_mm3[1] >= volume.volume_mm3[0]);
    assert!(bounds_mm[0][0] < bounds_mm[1][0]);
    let moved = body
        .transformed_with_interrupt([[1., 0., 0.], [1., 0., 0.], [0., 1., 0.]], 1e-7, 1e-6, None)
        .unwrap();
    assert!(
        verify_authored_compound(
            &input,
            &[ExpectedStepMember::RetainedRadialVolume {
                boundary: &moved,
                minimum_feature_mm: 1e-6
            }],
            1e-7
        )
        .is_err()
    );
    assert!(
        verify_authored_compound(
            &input,
            &[ExpectedStepMember::RetainedRadialVolume {
                boundary: &body,
                minimum_feature_mm: 0.
            }],
            1e-7
        )
        .is_err()
    );
    let changed = input.replace("ADVANCED_FACE(", "INVALID_FACE(");
    assert!(verify_authored_compound(&changed, &expected, 1e-7).is_err());
}
