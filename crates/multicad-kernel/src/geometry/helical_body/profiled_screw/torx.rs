//! Nominal finite curved-face evaluations. Source-family uncertainty stays
//! separate; these properties alone do not admit a solid or a STEP boundary.
pub(in crate::geometry::helical_body) mod display;

use super::*;
use multicad_rust_kernel::hexalobular_drive::NominalSurfaceEvaluation;

pub(super) fn qualify_source(
    screw: &Screw,
    frame: PlaneFrame,
    policy: GeometryPolicy,
) -> Result<(), GeometryError> {
    let source = screw
        .recess()
        .hexalobular_geometry()
        .ok_or_else(|| input("A curved recess requires its physical contour and cone"))?;
    let reserve = screw.recess().source_error_mm();
    // Spend the same complete source reserve as profiled_options. A retained
    // carrier tube is separate from the arithmetic-only head-meridian bound.
    if !reserve.is_finite() || (4. * reserve).next_up() > policy.linear_tolerance_mm {
        return Err(input("Curved recess source precision is unresolved"));
    }
    policy.length(source.cone().tip_depth_mm)?;
    let floor_drop = (source.cone().depth_per_radius * source.maximum_radius_bound_mm()).next_up();
    let shortest_wall =
        ((source.cone().tip_depth_mm - floor_drop).next_down() - reserve).next_down();
    policy.length(shortest_wall)?;
    for arc in source.arcs() {
        policy.length((arc.radius_mm - 2. * reserve).next_down())?;
        // Include both angular endpoint transfers and radius uncertainty;
        // policy feature conditioning keeps radius above forty source errors.
        policy.length(
            ((arc.radius_mm * arc.sweep_rad.abs()).next_down() - 16. * reserve).next_down(),
        )?;
    }
    // This complete local envelope includes every retained trim and finite
    // cone point. Checking its corners bounds all placed coordinates.
    let radius = (source.maximum_radius_bound_mm() + reserve).next_up();
    let top = screw.head_top_z();
    for x in [-radius, radius] {
        for y in [-radius, radius] {
            for z in [
                (top - reserve).next_down(),
                ((top + source.cone().tip_depth_mm).next_up() + reserve).next_up(),
            ] {
                policy.point(position(frame, [x, y, z]))?;
            }
        }
    }
    Ok(())
}

fn evaluated(
    evaluation: NominalSurfaceEvaluation,
    top: f64,
    policy: GeometryPolicy,
    scale: f64,
) -> Result<(Point, f64), GeometryError> {
    let [lo, hi] = evaluation.area_mm2;
    if !lo.is_finite()
        || !hi.is_finite()
        || lo <= 0.
        || hi < lo
        || hi - lo > policy.linear_tolerance_mm * scale / 16.
        || evaluation.center_mm.into_iter().any(|[a, b]| {
            !a.is_finite() || !b.is_finite() || b < a || b - a > policy.linear_tolerance_mm / 16.
        })
    {
        return Err(input("Unresolved nominal curved recess face evaluation"));
    }
    let mut center = evaluation.center_mm.map(|[a, b]| a.midpoint(b));
    center[2] += top;
    Ok((center, lo.midpoint(hi)))
}

pub(super) fn extend_face_properties(
    screw: &Screw,
    faces: &mut Vec<(Point, f64)>,
    policy: GeometryPolicy,
) -> Result<(), GeometryError> {
    let source = screw
        .recess()
        .hexalobular_geometry()
        .expect("sealed physical recess variants");
    let scale = source
        .maximum_radius_bound_mm()
        .max(source.cone().tip_depth_mm);
    let properties = source
        .nominal_surface_evaluations(policy.linear_tolerance_mm * scale / 256., 1_000_000)
        .map_err(|cause| {
            error(
                OPERATION,
                GeometryErrorKind::InvalidResult,
                format!("Cannot evaluate actual nominal curved recess surfaces: {cause:?}"),
            )
        })?;
    let top = screw.head_top_z();
    let (mouth_center, mouth_area) = evaluated(properties.mouth, top, policy, scale)?;
    let floor = evaluated(properties.floor, top, policy, scale)?;
    let walls = properties
        .walls
        .into_iter()
        .map(|wall| evaluated(wall, top, policy, scale))
        .collect::<Result<Vec<_>, _>>()?;
    let (removed_center, removed_area) = faces[0];
    let radius = screw.profile().end()[0];
    let shoulder = PI * radius * radius - removed_area;
    faces[0] = (
        [
            -removed_center[0] * removed_area / shoulder,
            -removed_center[1] * removed_area / shoulder,
            removed_center[2],
        ],
        shoulder,
    );
    let segments: Vec<_> = screw.profile().segments().collect();
    for role in screw.faces().into_iter().skip(faces.len()) {
        let property = match role {
            ScrewFace::HeadTop => {
                let area = PI * screw.profile().top()[0].powi(2) - mouth_area;
                (
                    [
                        -mouth_center[0] * mouth_area / area,
                        -mouth_center[1] * mouth_area / area,
                        top,
                    ],
                    area,
                )
            }
            ScrewFace::ProfileWall(i) => {
                let segment = segments[i as usize];
                let mut center = segment.center_mm();
                center[2] -= screw.spec().shank_length_mm;
                (center, segment.area_mm2())
            }
            ScrewFace::SocketFloor => floor,
            ScrewFace::SocketWall(i) => walls[i as usize],
            _ => unreachable!("Original helical face properties already exist"),
        };
        if property.1 <= 0.
            || !property.1.is_finite()
            || !property.0.into_iter().all(f64::is_finite)
        {
            return Err(input("Invalid nominal profiled screw face property"));
        }
        faces.push(property);
    }
    Ok(())
}

/// The owning screw retains the actual circle/cone intersection source. Its
/// nonplanar floor edges intentionally have no planar analytic substitute.
pub(super) fn extend_boundary(screw: &Screw, add: &mut Boundary<'_, ScrewFace>) {
    let source = screw
        .recess()
        .hexalobular_geometry()
        .expect("sealed recess");
    let top = screw.head_top_z();
    for (i, arc) in source.arcs().iter().enumerate() {
        let wall = ScrewFace::SocketWall(i as u32);
        let previous = ScrewFace::SocketWall(((i + 11) % 12) as u32);
        let [x, y] = arc.start_mm;
        let binding = match screw.recess() {
            multicad_rust_kernel::screw_recess::Recess::ConicalHexalobular(recess) => {
                recess.finite_binding()
            }
            _ => unreachable!("curved recess source"),
        };
        let nominal = [
            [x, y, 0.],
            [
                x,
                y,
                source.cone().tip_depth_mm - x.hypot(y) * source.cone().depth_per_radius,
            ],
        ];
        let [a, b] = binding
            .map_or(nominal, |bound| {
                [bound.vertices().mouth_mm[i], bound.vertices().floor_mm[i]]
            })
            .map(|mut p| {
                p[2] += top;
                p
            });
        let [line_a, line_b] = binding.map_or([a, b], |bound| {
            let post = bound.post_lines()[i];
            post.parameters.map(|t| {
                let mut p = std::array::from_fn(|k| {
                    post.origin_mm[k] + t * post.direction_mm_per_parameter[k]
                });
                p[2] += top;
                p
            })
        });
        add.line(
            format!("screw.drive/post/{i}"),
            line_a,
            line_b,
            &[previous, wall],
        );
        add.circle_at(
            format!("screw.drive/rim/0/{i}"),
            [arc.center_mm[0], arc.center_mm[1], top],
            arc.radius_mm,
            [arc.start_angle_rad, arc.sweep_rad],
            &[ScrewFace::HeadTop, wall],
        );
        add.surface_edge(
            format!("screw.drive/rim/1/{i}"),
            &[ScrewFace::SocketFloor, wall],
        );
        add.vertex(
            format!("screw.drive/corner/0/{i}"),
            a,
            &[ScrewFace::HeadTop, previous, wall],
        );
        add.vertex(
            format!("screw.drive/corner/1/{i}"),
            b,
            &[ScrewFace::SocketFloor, previous, wall],
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use multicad_rust_kernel::{
        head_profile::Head,
        helical_thread::Spec,
        hexalobular_drive::{Cone, Drive, TangentSpec},
        screw_recess::Recess,
    };
    pub(super) fn source(left: bool, plain: f64) -> Screw {
        let drive = Drive::from_authored_tangent(
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
        .unwrap();
        let recess = Recess::conical_hexalobular(drive, 1e-7, 100_000).unwrap();
        let thread = Thread::new(Spec {
            major_radius_mm: 1.5,
            pitch_mm: 0.5,
            length_mm: 3.,
            depth_mm: 0.2,
            crest_width_mm: 0.0625,
            lower_flank_deg: 30.,
            upper_flank_deg: 30.,
            starts: if left { 2 } else { 1 },
            left_hand: left,
            root_radius_mm: None,
        })
        .unwrap();
        Screw::new_with_recess(
            thread,
            1.5,
            Head::Crowned {
                radius_mm: 2.8,
                height_mm: 2.4,
                crown_radius_mm: 5.,
                top_land_radius_mm: 2.05,
                neck_radius_mm: 0.125,
            },
            plain,
            0.1,
            recess,
        )
        .unwrap()
    }
    #[test]
    fn retained_curves_vertices_and_post_supports_keep_literal_native_identity() {
        use multicad_rust_kernel::hexalobular_drive::{FiniteBinding, FinitePost, FiniteVertices};
        let original = source(false, 0.);
        let drive = original.recess().hexalobular_geometry().unwrap().clone();
        let mouth = drive.arcs().map(|a| [a.start_mm[0], a.start_mm[1], 0.]);
        let floor = mouth.map(|[x, y, _]| {
            [
                x,
                y,
                drive.cone().tip_depth_mm - drive.cone().depth_per_radius * x.hypot(y) + 1e-9,
            ]
        });
        let mut carriers =
            std::array::from_fn(|i| drive.intersection_carrier(i, 1e-10, 10_000).unwrap());
        for spans in &mut carriers {
            for span in spans {
                for p in &mut span.poles_mm {
                    p[2] += 1e-9;
                }
            }
        }
        let mut posts = std::array::from_fn(|i| FinitePost {
            origin_mm: mouth[i],
            direction_mm_per_parameter: std::array::from_fn(|k| floor[i][k] - mouth[i][k]),
            parameters: [0., 1.],
        });
        posts[0].origin_mm[0] += 2e-10;
        let actual_line = posts[0];
        let binding = FiniteBinding::new_with_posts(
            drive.clone(),
            FiniteVertices {
                mouth_mm: mouth,
                floor_mm: floor,
                apex_mm: [0., 0., drive.cone().tip_depth_mm + 1e-9],
            },
            posts,
            carriers,
            1e-8,
        )
        .unwrap();
        let recess = Recess::conical_hexalobular_from_binding(binding, 1e-5, 100_000).unwrap();
        let screw = Screw::new_with_recess(
            *original.thread(),
            original.profile().end()[0],
            original.spec().head.unwrap(),
            0.,
            0.1,
            recess,
        )
        .unwrap();
        let body = HelicalBody::new(
            *screw.thread(),
            EntityId::new(),
            HelicalForm::ProfiledScrew(Box::new(screw.clone())),
            None,
            GeometryPolicy::default(),
        )
        .unwrap();
        let snapshot = body
            .topology(PlaneFrame::default(), body.feature, 0)
            .unwrap();
        for (i, expected_floor) in floor.iter().enumerate() {
            let vertex = snapshot
                .elements
                .iter()
                .find(|e| {
                    e.key.kind == TopologyKind::Vertex
                        && e.origins
                            .contains(&body.origin(format!("screw.drive/corner/1/{i}")))
                })
                .unwrap();
            let Some(AnalyticGeometry::Point { position_mm: p }) = vertex.geometry else {
                panic!("literal floor vertex")
            };
            assert_eq!(
                p,
                [
                    expected_floor[0],
                    expected_floor[1],
                    expected_floor[2] + screw.head_top_z()
                ]
            );
        }
        let post = snapshot
            .elements
            .iter()
            .find(|e| {
                e.key.kind == TopologyKind::Edge
                    && e.origins
                        .contains(&body.origin("screw.drive/post/0".into()))
            })
            .unwrap();
        let Some(AnalyticGeometry::Line { a_mm, b_mm }) = post.geometry else {
            panic!("retained post support")
        };
        let expected = actual_line.parameters.map(|t| {
            let mut p = std::array::from_fn(|k| {
                actual_line.origin_mm[k] + t * actual_line.direction_mm_per_parameter[k]
            });
            p[2] += screw.head_top_z();
            p
        });
        assert_eq!([a_mm, b_mm], expected);
        assert_ne!(a_mm[0], mouth[0][0]);
        let mut remaining = MAX_DISPLAY_EDGE_POINTS;
        let edges =
            display::edges(&body, &screw, PlaneFrame::default(), 0.01, &mut remaining).unwrap();
        for (i, edge) in edges.iter().enumerate() {
            assert_eq!(
                edge.points_mm[0],
                [floor[i][0], floor[i][1], floor[i][2] + screw.head_top_z()]
            );
        }
        assert!(
            Solid::rust_profiled_screw_from_source(
                screw,
                EntityId::new(),
                PlaneFrame::default(),
                GeometryPolicy::default()
            )
            .is_err()
        );
    }
    #[test]
    fn nominal_curved_faces_preserve_actual_cone_and_circles_without_planar_reclassification() {
        let policy = GeometryPolicy::default();
        for left in [false, true] {
            for plain in [0., 3.] {
                let source = source(left, plain);
                let original = source.clone();
                let body = HelicalBody::new(
                    *source.thread(),
                    EntityId::new(),
                    HelicalForm::ProfiledScrew(Box::new(source.clone())),
                    None,
                    policy,
                )
                .unwrap();
                assert_eq!(body.local_faces.len(), source.faces().len());
                let snapshot = body
                    .topology(PlaneFrame::default(), body.feature, 0)
                    .unwrap();
                let floor = snapshot
                    .elements
                    .iter()
                    .find(|e| {
                        e.key.kind == TopologyKind::Face
                            && e.origins.contains(&body.origin("screw.drive/end".into()))
                    })
                    .unwrap();
                assert_eq!(floor.boundary.len(), 12);
                assert!(floor.geometry.is_none() && floor.frames.is_empty());
                for i in 0..12 {
                    let origin = body.origin(format!("screw.drive/rim/1/{i}"));
                    let edge = snapshot
                        .elements
                        .iter()
                        .find(|e| e.key.kind == TopologyKind::Edge && e.origins.contains(&origin))
                        .unwrap();
                    assert!(edge.geometry.is_none());
                    assert_eq!(
                        snapshot
                            .elements
                            .iter()
                            .filter(|e| e.key.kind == TopologyKind::Face
                                && e.boundary.contains(&origin))
                            .count(),
                        2
                    );
                    let contact =
                        source.recess().hexalobular_geometry().unwrap().arcs()[i].start_mm;
                    let vertex = snapshot
                        .elements
                        .iter()
                        .find(|e| {
                            e.key.kind == TopologyKind::Vertex
                                && e.origins
                                    .contains(&body.origin(format!("screw.drive/corner/1/{i}")))
                        })
                        .unwrap();
                    let Some(AnalyticGeometry::Point {
                        position_mm: actual,
                    }) = vertex.geometry
                    else {
                        panic!("actual cone contact")
                    };
                    assert_eq!(actual[0], contact[0]);
                    assert_eq!(actual[1], contact[1]);
                    assert!(actual[2] > source.head_top_z());
                }

                let properties = source
                    .recess()
                    .hexalobular_geometry()
                    .unwrap()
                    .nominal_surface_evaluations(1e-8, 1_000_000)
                    .unwrap();
                for (index, role) in source.faces().into_iter().enumerate() {
                    let (center, area) = body.local_faces[index];
                    assert!(area.is_finite() && area > 0.);
                    let expected = match role {
                        ScrewFace::SocketFloor => Some(properties.floor),
                        ScrewFace::SocketWall(i) => Some(properties.walls[i as usize]),
                        _ => None,
                    };
                    if let Some(expected) = expected {
                        assert!(
                            expected.area_mm2[0] - 1e-8 <= area
                                && area <= expected.area_mm2[1] + 1e-8
                        );
                        for (axis, [lo, hi]) in expected.center_mm.into_iter().enumerate() {
                            let actual =
                                center[axis] - if axis == 2 { source.head_top_z() } else { 0. };
                            assert!(lo - 1e-8 <= actual && actual <= hi + 1e-8);
                        }
                        assert!(plane(&source, role, PlaneFrame::default()).is_none());
                    }
                }
                let floor_index = source
                    .faces()
                    .iter()
                    .position(|role| *role == ScrewFace::SocketFloor)
                    .unwrap();
                let top_index = source
                    .faces()
                    .iter()
                    .position(|role| *role == ScrewFace::HeadTop)
                    .unwrap();
                let mouth = properties.mouth.area_mm2[0].midpoint(properties.mouth.area_mm2[1]);
                let (floor_center, floor_area) = body.local_faces[floor_index];
                let slope = source
                    .recess()
                    .hexalobular_geometry()
                    .unwrap()
                    .cone()
                    .depth_per_radius;
                assert!((floor_area / mouth - (1. + slope * slope).sqrt()).abs() < 1e-8);
                assert!(
                    (body.local_faces[top_index].1 + mouth
                        - PI * source.profile().top()[0].powi(2))
                    .abs()
                        < 1e-8
                );
                assert!(floor_center[2] > source.head_top_z());
                assert!(floor_center[2] < source.head_top_z() + source.recess().maximum_depth_mm());
                assert!(floor_center[0].abs() < 1e-8 && floor_center[1].abs() < 1e-8);
                assert!(plane(&source, ScrewFace::HeadTop, PlaneFrame::default()).is_some());
                assert_eq!(source, original);
                let published = Solid::rust_profiled_screw_from_source(
                    source.clone(),
                    EntityId::new(),
                    PlaneFrame::default(),
                    policy,
                )
                .unwrap();
                let mesh = published
                    .manufacturing_mesh_with_receipt(crate::MeshOptions::default())
                    .unwrap();
                assert_eq!(mesh.edge_pairs as usize * 2, mesh.mesh.indices.len());
                let directory = tempfile::tempdir().unwrap();
                let path = directory.path().join("physical-torx.step");
                published.write_step(&path).unwrap();
                let decoded = Solid::import_rust_step(
                    &std::fs::read_to_string(path).unwrap(),
                    EntityId::new(),
                    policy,
                )
                .unwrap();
                let (physical, _) = decoded.rust_profiled_screw_geometry().unwrap();
                assert!(physical.spec().head.is_none());
                assert!(physical.recess().flat_hex_geometry().is_none());
                let multicad_rust_kernel::screw_recess::Recess::ConicalHexalobular(recess) =
                    physical.recess()
                else {
                    panic!("physical curved cavity");
                };
                assert!(recess.finite_binding().is_some());
                assert!(
                    (decoded.volume_mm3().unwrap() - published.volume_mm3().unwrap()).abs() < 1e-5
                );
                assert_eq!(source, original);
            }
        }
    }
}
