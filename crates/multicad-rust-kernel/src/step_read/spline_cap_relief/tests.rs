use super::*;
use crate::{
    normal_sweep::Edge as A,
    spline_cap_relief::{Kind, tests::fixture},
    spline_sweep::{Cubic, Sweep},
};
const FRAME: [Point; 3] = [[0.; 3], [1., 0., 0.], [0., 1., 0.]];
fn write(relief: &Relief) -> String {
    crate::step::write_spline_cap_relief(relief, FRAME[0], FRAME[1], FRAME[2]).unwrap()
}
fn changed(relief: &Relief) -> Relief {
    let m = &relief.motifs()[0];
    let mut holes = m.holes().to_vec();
    let Edge::Cubic(c) = &holes[0][0] else {
        panic!()
    };
    let mut poles = c.poles().to_vec();
    poles[1][1] += 0.01;
    holes[0][0] = Edge::Cubic(Cubic::new(poles, c.knots().to_vec(), c.reversed()).unwrap());
    let m = Sweep::new_rings(m.outline().to_vec(), holes, m.height_mm(), 1e-7).unwrap();
    Relief::new(relief.base().clone(), vec![m], relief.kind(), 1e-7, 1e-4).unwrap()
}
#[test]
fn generic_import_and_full_matcher_preserve_curves_under_record_loop_and_member_order_changes() {
    for kind in [Kind::Raised, Kind::Recessed] {
        let relief = fixture(kind);
        let source = write(&relief);
        let mut reverse = source.clone();
        let (before, data) = source.split_once("DATA;").unwrap();
        let (data, after) = data.split_once("ENDSEC;").unwrap();
        reverse.clear();
        reverse.push_str(before);
        reverse.push_str("DATA;\n");
        reverse.push_str(&data.lines().rev().collect::<Vec<_>>().join("\n"));
        reverse.push_str("\nENDSEC;");
        reverse.push_str(after);
        for text in [&source, &reverse] {
            let (proof, imported) = read_import_solid(text, 1e-7).unwrap();
            assert_eq!(proof.solid_count, 1);
            assert!((proof.volume_mm3 - relief.volume_mm3()).abs() < 1e-8);
            assert!(spline_cap_relief_import_matches(
                &imported, &relief, FRAME, 1e-7
            ));
            assert!(!spline_cap_relief_import_matches(
                &imported,
                &changed(&relief),
                FRAME,
                1e-7
            ));
            for tol in [0., f64::NAN] {
                assert!(!spline_cap_relief_import_matches(
                    &imported, &relief, FRAME, tol
                ));
            }
            let mut frame = FRAME;
            frame[0][0] = 1.;
            assert!(!spline_cap_relief_import_matches(
                &imported, &relief, frame, 1e-7
            ));
            frame = FRAME;
            frame[1][0] = f64::NAN;
            assert!(!spline_cap_relief_import_matches(
                &imported, &relief, frame, 1e-7
            ));
            let mut altered = imported.clone();
            altered.u[0] = f64::NAN;
            assert!(!spline_cap_relief_import_matches(
                &altered, &relief, FRAME, 1e-7
            ));
            altered = imported.clone();
            altered.u[0] += 4e-11;
            assert!(!spline_cap_relief_import_matches(
                &altered, &relief, FRAME, 1e-7
            ));
            altered = imported.clone();
            altered.profile = ImportProfile::Cylinder { radius_mm: 2. };
            assert!(!spline_cap_relief_import_matches(
                &altered, &relief, FRAME, 1e-7
            ));
            assert!(read_analytic_solid(text, 1e-7).is_ok());
        }
        let (mut outline, mut holes) = (
            relief.motifs()[0].outline().to_vec(),
            relief.motifs()[0].holes().to_vec(),
        );
        outline.rotate_left(1);
        holes[0].rotate_left(2);
        holes.reverse();
        let rotated = Sweep::new_rings(outline, holes, 1., 1e-7).unwrap();
        let rotated = Relief::new(relief.base().clone(), vec![rotated], kind, 1e-7, 1e-4).unwrap();
        let (_, imported) = read_import_solid(&write(&rotated), 1e-7).unwrap();
        assert!(spline_cap_relief_import_matches(
            &imported, &relief, FRAME, 1e-7
        ));
    }
}
#[test]
fn arcs_and_multiple_profiles_keep_one_solid_and_independent_member_ownership() {
    let source = fixture(Kind::Raised);
    let original = &source.motifs()[0];
    let outer = vec![
        Edge::Analytic(A::Line {
            a: [2., 2.],
            b: [7., 2.],
        }),
        Edge::Analytic(A::Arc {
            a: [7., 2.],
            b: [9.5, 4.5],
            center: [7., 4.5],
            clockwise: false,
        }),
        Edge::Analytic(A::Arc {
            a: [9.5, 4.5],
            b: [7., 7.],
            center: [7., 4.5],
            clockwise: false,
        }),
        Edge::Analytic(A::Line {
            a: [7., 7.],
            b: [2., 7.],
        }),
        Edge::Analytic(A::Line {
            a: [2., 7.],
            b: [2., 2.],
        }),
    ];
    let motif = Sweep::from_verified_rings(outer, original.holes().to_vec(), 1., 1e-7).unwrap();
    let relief = Relief::new(
        source.base().clone(),
        vec![motif],
        Kind::Recessed,
        1e-7,
        1e-4,
    )
    .unwrap();
    let text = write(&relief);
    let (proof, decoded) = read_import_solid(&text, 1e-7).unwrap();
    assert_eq!(proof.arcs.len(), 4);
    assert!(spline_cap_relief_import_matches(
        &decoded, &relief, FRAME, 1e-7
    ));
    let mut writer = crate::step::Writer::new();
    writer
        .append_spline_cap_relief(&relief, FRAME[0], FRAME[1], FRAME[2])
        .unwrap();
    writer
        .append_spline_cap_relief(&source, [30., 0., 0.], FRAME[1], FRAME[2])
        .unwrap();
    let text = writer.finish().unwrap();
    let members = read_spline_cap_relief_compound(&text, 1e-7).unwrap();
    assert_eq!(members.len(), 2);
    assert!(read_spline_cap_relief(&text, 1e-7).is_err());
    assert_eq!(read_import_compound(&text, 1e-7).unwrap().len(), 2);
}

fn with_shell<T>(text: &str, f: impl FnOnce(&mut Shell<'_, '_>, u32) -> T) -> T {
    let entities = parse_entities(text, 1e-7).unwrap();
    let mut shell = Shell::new(&entities, 1e-7, true);
    let solid = shell.represented_solids(true).unwrap()[0];
    f(&mut shell, solid)
}
fn plane_origin(text: &str, face_index: usize, z: f64) -> String {
    let (axis, normal, u, point, id) = with_shell(text, |shell, solid| {
        let ids = shell.member_faces(solid).unwrap();
        let (surface, _, _) = helical::face(shell, ids[face_index]).unwrap();
        let [Value::String(_), axis] = shell.args(surface, "PLANE").unwrap() else {
            panic!()
        };
        let axis = reference(axis).unwrap();
        let [Value::String(_), point, normal, u] = shell.args(axis, "AXIS2_PLACEMENT_3D").unwrap()
        else {
            panic!()
        };
        (
            axis,
            reference(normal).unwrap(),
            reference(u).unwrap(),
            shell.point(reference(point).unwrap()).unwrap(),
            shell.entities.keys().max().unwrap() + 1,
        )
    });
    let row = text
        .lines()
        .find(|line| line.starts_with(&format!("#{axis}=")))
        .unwrap();
    let mut output = text.replacen(
        row,
        &format!("#{axis}=AXIS2_PLACEMENT_3D('',#{id},#{normal},#{u});"),
        1,
    );
    let at = output.find("ENDSEC;\nEND-ISO").unwrap();
    output.insert_str(
        at,
        &format!(
            "#{id}=CARTESIAN_POINT('',({:.16E},{:.16E},{z:.16E}));\n",
            point[0], point[1]
        ),
    );
    output
}
#[test]
fn every_known_curve_class_and_untrusted_ownership_refuses_incompatible_finite_geometry() {
    let relief = fixture(Kind::Raised);
    let text = write(&relief);
    let frame = Frame::new(FRAME[0], FRAME[1], FRAME[2]).unwrap();
    with_shell(&text, |shell, solid| {
        let mut actual = Edges::new();
        for id in shell.member_faces(solid).unwrap() {
            spline_sweep::read_face(shell, id, &mut actual).unwrap();
        }
        let cubic = actual
            .values()
            .find(|e| matches!(e.geometry, Geometry::Cubic { .. }))
            .unwrap();
        let circle = actual
            .values()
            .find(|e| matches!(e.geometry, Geometry::Analytic(Curve::Circle(_))))
            .unwrap();
        let line = actual
            .values()
            .find(|e| matches!(e.geometry, Geometry::Analytic(Curve::Line)))
            .unwrap();
        let expected = relief
            .boundary()
            .edges
            .iter()
            .find(|e| {
                matches!(
                    e.curve,
                    ExpectedCurve::Profile {
                        edge: Edge::Cubic(_),
                        ..
                    }
                )
            })
            .unwrap();
        assert!(curve_matches(
            &cubic.geometry,
            &expected.curve,
            true,
            frame,
            1e-7
        ));
        assert!(!curve_matches(
            &cubic.geometry,
            &expected.curve,
            false,
            frame,
            1e-7
        ));
        let ExpectedCurve::Profile {
            edge: Edge::Cubic(c),
            z_mm,
        } = &expected.curve
        else {
            panic!()
        };
        let altered = ExpectedCurve::Profile {
            edge: Edge::Cubic(Cubic::new(c.poles().to_vec(), vec![0., 2.], c.reversed()).unwrap()),
            z_mm: *z_mm,
        };
        assert!(!curve_matches(&cubic.geometry, &altered, true, frame, 1e-7));
        let circle_source = relief
            .boundary()
            .edges
            .iter()
            .find(|e| {
                matches!(
                    e.curve,
                    ExpectedCurve::Profile {
                        edge: Edge::Analytic(A::Circle { .. }),
                        ..
                    }
                )
            })
            .unwrap();
        assert!(!curve_matches(
            &line.geometry,
            &circle_source.curve,
            true,
            frame,
            1e-7
        ));
        let expected_line = ExpectedCurve::Profile {
            edge: Edge::Analytic(A::Line {
                a: [0., 0.],
                b: [1., 0.],
            }),
            z_mm: 0.,
        };
        assert!(!curve_matches(
            &circle.geometry,
            &expected_line,
            true,
            frame,
            1e-7
        ));
        assert!(!curve_matches(
            &circle.geometry,
            &expected.curve,
            true,
            frame,
            1e-7
        ));
        assert!(!curve_matches(
            &circle.geometry,
            &ExpectedCurve::Line {
                a_mm: [0.; 3],
                b_mm: [1.; 3]
            },
            true,
            frame,
            1e-7
        ));
    });
    let mut shifted = FRAME;
    shifted[0][0] = 1.;
    assert_eq!(
        verify_spline_cap_relief(&text, &relief, shifted[0], shifted[1], shifted[2], 1e-7).err(),
        Some(ReadError::Invalid)
    );
    assert_eq!(
        verify_spline_cap_relief(&text, &changed(&relief), FRAME[0], FRAME[1], FRAME[2], 1e-7)
            .err(),
        Some(ReadError::Invalid)
    );
    let mut work = 1;
    charge(&mut work, 1).unwrap();
    assert_eq!(charge(&mut work, 1), Err(ReadError::Limit));
    assert_eq!(work, 0);
}
#[test]
fn blind_candidate_support_depth_resource_and_shell_ownership_failures_are_typed() {
    let relief = fixture(Kind::Raised);
    let source = write(&relief);
    let low = source
        .lines()
        .find(|line| line.contains("=ADVANCED_FACE("))
        .unwrap();
    let wrong = source.replacen(low, &low.replace(",.T.);", ",.F.);"), 1);
    assert_eq!(
        read_spline_cap_relief(&wrong, 1e-7).err(),
        Some(ReadError::Unsupported)
    );
    let cap = relief
        .boundary()
        .faces
        .iter()
        .position(|f| {
            matches!(
                f.role,
                crate::spline_cap_relief::FaceRole::MotifCap { motif: 0 }
            )
        })
        .unwrap();
    assert!(read_spline_cap_relief(&plane_origin(&source, cap, 5.), 1e-7).is_err());
    let (vector, direction_id) = with_shell(&source, |shell, solid| {
        let faces = shell.member_faces(solid).unwrap();
        let face = faces
            .into_iter()
            .find(|&id| spline_sweep::has_extrusion_support(shell, id))
            .unwrap();
        let (surface, _, _) = helical::face(shell, face).unwrap();
        let [Value::String(_), _, vector] =
            shell.args(surface, "SURFACE_OF_LINEAR_EXTRUSION").unwrap()
        else {
            panic!()
        };
        let vector = reference(vector).unwrap();
        let [Value::String(_), direction, _] = shell.args(vector, "VECTOR").unwrap() else {
            panic!()
        };
        (vector, reference(direction).unwrap())
    });
    let _ = vector;
    let row = source
        .lines()
        .find(|l| l.starts_with(&format!("#{direction_id}=")))
        .unwrap();
    let altered = source.replacen(
        row,
        &format!("#{direction_id}=DIRECTION('',(1.,0.,0.));"),
        1,
    );
    assert_eq!(
        verify_spline_cap_relief(&altered, &relief, FRAME[0], FRAME[1], FRAME[2], 1e-7).err(),
        Some(ReadError::Invalid)
    );
    let row = source
        .lines()
        .find(|l| l.contains("=CLOSED_SHELL("))
        .unwrap();
    let (id, _) = row.split_once('=').unwrap();
    let face = with_shell(&source, |s, id| s.member_faces(id).unwrap()[0]);
    let ids = std::iter::repeat_n(
        format!("#{face}"),
        3 + 6 * crate::spline_cap_relief::MAX_EDGES,
    )
    .collect::<Vec<_>>()
    .join(",");
    // A spline-supported shell below the legacy reader's minimum face
    // count is outside that family, not a resource refusal (CAD-359); the
    // new provider must not turn it into an import.
    let short_ids = with_shell(&source, |s, solid| {
        let faces = s.member_faces(solid).unwrap();
        let spline = *faces
            .iter()
            .find(|&&f| spline_sweep::has_extrusion_support(s, f))
            .unwrap();
        std::iter::once(spline)
            .chain(faces.into_iter().filter(|&f| f != spline).take(3))
            .map(|f| format!("#{f}"))
            .collect::<Vec<_>>()
            .join(",")
    });
    let short = source.replacen(row, &format!("{id}=CLOSED_SHELL('',({short_ids}));"), 1);
    assert_eq!(
        read_import_solid(&short, 1e-7).err(),
        Some(ReadError::Unsupported)
    );
    assert_eq!(
        read_analytic_solid(&short, 1e-7).err(),
        Some(ReadError::Unsupported)
    );
    let huge = source.replacen(row, &format!("{id}=CLOSED_SHELL('',({ids}));"), 1);
    let cause = read_spline_cap_relief(&huge, 1e-7).err().unwrap();
    assert_eq!(cause, ReadError::Limit);
    assert_eq!(
        compound::spline_error(cause, ReadError::Invalid),
        ReadError::Limit
    );
    assert_eq!(
        compound::spline_error(ReadError::Unsupported, ReadError::Invalid),
        ReadError::Invalid
    );
    assert_eq!(
        compound::spline_error(ReadError::Invalid, ReadError::Unsupported),
        ReadError::Unsupported
    );
    let mut owners = BTreeMap::from([(1, vec![0])]);
    assert_eq!(recognize::other(&owners, 1, 0), Err(ReadError::Invalid));
    owners.insert(1, vec![1, 2]);
    assert_eq!(recognize::other(&owners, 1, 0), Err(ReadError::Invalid));
    owners.insert(1, vec![0, 1]);
    assert_eq!(recognize::other(&owners, 1, 0), Ok(1));
    assert_eq!(
        recognize::error(crate::spline_sweep::ConstructionError::Limit),
        ReadError::Limit
    );
    assert_eq!(
        recognize::error(crate::spline_sweep::ConstructionError::Unresolved),
        ReadError::Invalid
    );
}

#[test]
fn retained_stock_through_holes_remain_distinct_and_mixed_depth_proposals_refuse() {
    let source = fixture(Kind::Raised);
    let base = Sweep::new_rings(
        source.base().outline().to_vec(),
        vec![vec![Edge::Analytic(A::Circle {
            center: [14., 14.],
            radius_mm: 1.,
            clockwise: true,
        })]],
        5.,
        1e-7,
    )
    .unwrap();
    let p = [[12., 2.], [15., 2.], [15., 5.], [12., 5.]];
    let second = Sweep::new(
        (0..4)
            .map(|i| {
                Edge::Analytic(A::Line {
                    a: p[i],
                    b: p[(i + 1) % 4],
                })
            })
            .collect(),
        None,
        1.5,
        1e-7,
    )
    .unwrap();
    let relief = Relief::new(
        base,
        vec![source.motifs()[0].clone(), second],
        Kind::Raised,
        1e-7,
        1e-4,
    )
    .unwrap();
    let text = write(&relief);
    let (_, decoded) = read_spline_cap_relief(&text, 1e-7).unwrap();
    assert_eq!(decoded.relief.base().holes().len(), 1);
    assert_eq!(decoded.relief.motifs().len(), 2);
    let cap = relief
        .boundary()
        .faces
        .iter()
        .position(|f| {
            matches!(
                f.role,
                crate::spline_cap_relief::FaceRole::MotifCap { motif: 1 }
            )
        })
        .unwrap();
    assert_eq!(
        read_spline_cap_relief(&plane_origin(&text, cap, 3.5), 1e-7).err(),
        Some(ReadError::Unsupported)
    );
    let text = write(&source);
    let to_flip = with_shell(&text, |shell, solid| {
        let ids = shell.member_faces(solid).unwrap();
        let mut bounds = Vec::new();
        for (i, face) in source.boundary().faces.iter().enumerate() {
            let (_, _, current) = helical::face(shell, ids[i]).unwrap();
            if matches!(
                face.role,
                crate::spline_cap_relief::FaceRole::MotifCap { .. }
                    | crate::spline_cap_relief::FaceRole::CounterCap { .. }
                    | crate::spline_cap_relief::FaceRole::MotifWall { .. }
            ) {
                bounds.extend(current);
            } else if matches!(
                face.role,
                crate::spline_cap_relief::FaceRole::BaseCap { top: true }
            ) {
                bounds.push(current[1]);
            }
        }
        bounds
    });
    let mut changed = text.clone();
    for id in to_flip {
        let row = text
            .lines()
            .find(|l| l.starts_with(&format!("#{id}=")))
            .unwrap();
        changed = changed.replacen(row, &row.replace(",.T.);", ",.F.);"), 1);
    }
    assert_eq!(
        read_spline_cap_relief(&changed, 1e-7).err(),
        Some(ReadError::Invalid)
    );
}
