//! Nominal finite curved-face evaluations. Source-family uncertainty stays
//! separate; these properties alone do not admit a solid or a STEP boundary.
pub(in crate::geometry::helical_body) mod display;

use super::*;
use multicad_rust_kernel::hexalobular_drive::NominalSurfaceEvaluation;

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
        let a = [x, y, top];
        let b = [
            x,
            y,
            top + source.cone().tip_depth_mm - x.hypot(y) * source.cone().depth_per_radius,
        ];
        add.line(format!("screw.drive/post/{i}"), a, b, &[previous, wall]);
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
                // Source evaluations do not admit unresolved native incidence/exchange.
                assert!(
                    Solid::rust_profiled_screw_from_source(
                        source,
                        EntityId::new(),
                        PlaneFrame::default(),
                        policy
                    )
                    .is_err()
                );
            }
        }
    }
}
