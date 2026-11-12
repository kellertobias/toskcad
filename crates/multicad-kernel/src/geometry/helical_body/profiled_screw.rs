//! Policy admission and finite face identities for actual profiled heads.
use super::{
    boundary::{Boundary, edge_midpoint, hex_wall_plane, placed_frame},
    *,
};
use multicad_rust_kernel::{
    head_profile::{Curve, Role},
    profiled_screw::{Face as ScrewFace, Screw, Spec as ScrewSpec},
};

impl Solid {
    /// Current retained head meridian and true helical thread, in its placed
    /// thread frame. Source meridian Z is shifted by minus the plain length.
    pub fn rust_profiled_screw_geometry(&self) -> Option<(&Screw, PlaneFrame)> {
        let SolidBody::Helical(body, frame) = &self.body else {
            return None;
        };
        let HelicalForm::ProfiledScrew(screw) = &body.form else {
            return None;
        };
        Some((screw, *frame))
    }

    pub fn rust_profiled_machine_screw(
        thread: HelicalThread,
        root_radius_mm: Option<f64>,
        screw: ScrewSpec,
        feature: EntityId,
        underhead: PlaneFrame,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        policy.validate()?;
        underhead.checked(policy)?;
        policy.length(screw.socket_across_flats_mm)?;
        policy.length(screw.socket_depth_mm)?;
        policy.length(screw.minimum_wall_mm)?;
        if screw.shank_length_mm != 0. {
            policy.length(screw.shank_length_mm)?;
        }
        let frame = PlaneFrame {
            origin: position(underhead, [0., 0., screw.shank_length_mm]),
            ..underhead
        };
        Self::rust_helical_body(thread, root_radius_mm, feature, frame, policy, |thread| {
            let screw = Screw::new(thread, screw).map_err(|cause| match cause {
                ThreadError::Invalid => {
                    input("Invalid profiled screw head, neck, shank or hex drive")
                }
                ThreadError::Limit => error(
                    OPERATION,
                    GeometryErrorKind::Limit,
                    "Profiled screw exceeds the geometry envelope",
                ),
            })?;
            let reserve = screw.profile().rounding_error_mm();
            if reserve > policy.linear_tolerance_mm / 8. {
                return Err(input(
                    "Head source arithmetic exceeds the requested geometry precision",
                ));
            }
            policy.length(screw.profile().top()[0])?;
            policy.length(-screw.profile().top()[1])?;
            for segment in screw.profile().segments() {
                if segment.area_mm2() <= 0.
                    || !segment.area_mm2().is_finite()
                    || !segment.center_mm().into_iter().all(f64::is_finite)
                {
                    return Err(input("Unrepresentable profiled head surface properties"));
                }
                for p in [segment.start, segment.end] {
                    policy.point(position(underhead, [p[0], 0., p[1]]))?;
                }
                if let Curve::Arc {
                    center,
                    radius_mm: r,
                    sweep,
                } = segment.curve
                {
                    policy.length(r)?;
                    let angle = (segment.start[1] - center[1]).atan2(segment.start[0] - center[0]);
                    for (fraction, endpoint) in [(0., segment.start), (1., segment.end)] {
                        let t = angle + fraction * sweep;
                        let actual = [center[0] + r * t.cos(), center[1] + r * t.sin()];
                        let residual = (actual[0] - endpoint[0]).hypot(actual[1] - endpoint[1]);
                        if residual + reserve > policy.linear_tolerance_mm / 4. {
                            return Err(input(
                                "Head arc endpoints exceed the requested trim precision",
                            ));
                        }
                    }
                } else {
                    policy.length(
                        (segment.end[0] - segment.start[0])
                            .hypot(segment.end[1] - segment.start[1]),
                    )?;
                }
            }
            Ok(HelicalForm::ProfiledScrew(Box::new(screw)))
        })
    }
}

pub(super) fn extend_face_properties(screw: &Screw, faces: &mut Vec<(Point, f64)>) {
    let spec = screw.spec();
    let radius = screw.thread().spec().major_radius_mm;
    let (removed_center, removed_area) = faces[0];
    let shoulder = PI * radius.powi(2) - removed_area;
    faces[0] = (
        removed_center.map(|x| -x * removed_area / shoulder),
        shoulder,
    );
    let top = screw.head_top_z();
    let socket_area = 3_f64.sqrt() * spec.socket_across_flats_mm.powi(2) / 2.;
    let segments: Vec<_> = screw.profile().segments().collect();
    for face in screw.faces().into_iter().skip(faces.len()) {
        faces.push(match face {
            ScrewFace::HeadTop => (
                [0., 0., top],
                PI * screw.profile().top()[0].powi(2) - socket_area,
            ),
            ScrewFace::ProfileWall(i) => {
                let segment = segments[i as usize];
                let mut center = segment.center_mm();
                center[2] -= spec.shank_length_mm;
                (center, segment.area_mm2())
            }
            ScrewFace::SocketFloor => ([0., 0., top + spec.socket_depth_mm], socket_area),
            ScrewFace::SocketWall(i) => (
                edge_midpoint(
                    spec.socket_across_flats_mm,
                    i,
                    top + spec.socket_depth_mm / 2.,
                ),
                spec.socket_across_flats_mm * spec.socket_depth_mm / 3_f64.sqrt(),
            ),
            _ => unreachable!("Original helical face properties already exist"),
        });
    }
}

pub(super) fn plane(screw: &Screw, role: ScrewFace, frame: PlaneFrame) -> Option<ReferenceFrame> {
    let axial = |z, sign| placed_frame(frame, [0., 0., z], [1., 0., 0.], [0., sign, 0.]);
    let top = screw.head_top_z();
    match role {
        ScrewFace::HeadTop => Some(axial(top, -1.)),
        ScrewFace::ProfileWall(i) => {
            let segment = screw
                .profile()
                .segments()
                .nth(i as usize)
                .expect("source face");
            (segment.role == Role::Bearing)
                .then(|| axial(segment.start[1] - screw.spec().shank_length_mm, 1.))
        }
        ScrewFace::SocketFloor => Some(axial(top + screw.spec().socket_depth_mm, -1.)),
        ScrewFace::SocketWall(i) => Some(hex_wall_plane(
            frame,
            screw.spec().socket_across_flats_mm,
            i,
            top,
            -1.,
        )),
        _ => None,
    }
}

fn face_name(screw: &Screw, role: ScrewFace) -> String {
    match role {
        ScrewFace::HeadTop => "screw.head/top.land".into(),
        ScrewFace::ProfileWall(i) => format!(
            "screw.head/{}",
            screw
                .profile()
                .segments()
                .nth(i as usize)
                .expect("source face")
                .role
                .source_role()
        ),
        ScrewFace::SocketFloor => "screw.drive/end".into(),
        ScrewFace::SocketWall(i) => format!("screw.drive/edge/{i}"),
        _ => unreachable!("Original helical face identities already exist"),
    }
}

pub(super) fn extend_topology(
    body: &HelicalBody,
    screw: &Screw,
    frame: PlaneFrame,
    elements: &mut Vec<TopologyElement>,
) {
    let roles = screw.faces();
    let count = body.rails.len();
    let segments: Vec<_> = screw.profile().segments().collect();
    for (ordinal, &role) in roles.iter().enumerate().skip(2 + count) {
        let origin = body.origin(face_name(screw, role));
        let plane = plane(screw, role, frame);
        elements.push(TopologyElement {
            key: TopologyKey {
                kind: TopologyKind::Face,
                ordinal: ordinal as u32,
            },
            origins: BTreeSet::from([origin.clone()]),
            boundary: BTreeSet::new(),
            frames: plane
                .into_iter()
                .map(|frame| AnchoredFrame {
                    origin: origin.clone(),
                    frame,
                })
                .collect(),
            geometry: plane.map(|frame| AnalyticGeometry::Plane { frame }),
        });
    }
    let mut add = Boundary::new(body, frame, &roles, elements, false);
    let last = ScrewFace::ProfileWall(segments.len() as u32 - 1);
    let shoulder = ScrewFace::ThreadShoulder;
    let rate = screw.thread().helical_rate();
    for i in 0..count {
        if (body.band(i) - 2) % 4 == 1 {
            let edge = body.rim_origin(0, i);
            add.face(shoulder).boundary.remove(&edge);
            add.face(last).boundary.insert(edge);
        } else {
            let phase = body.rails[i];
            add.arc(
                format!("screw.head/thread.junction/{}", body.rail_name(i)),
                0.,
                screw.thread().spec().major_radius_mm,
                [-phase / rate, (phase - body.next_phase(i)) / rate],
                &[shoulder, last],
            );
        }
    }
    for (i, segment) in segments.iter().enumerate() {
        let before = if i == 0 {
            ScrewFace::HeadTop
        } else {
            ScrewFace::ProfileWall(i as u32 - 1)
        };
        let current = ScrewFace::ProfileWall(i as u32);
        add.arc(
            format!("screw.head/{}/start", segment.role.source_role()),
            segment.start[1] - screw.spec().shank_length_mm,
            segment.start[0],
            [0., TAU],
            &[before, current],
        );
    }
    add.hex_prism(
        "screw.drive",
        screw.spec().socket_across_flats_mm,
        [
            screw.head_top_z(),
            screw.head_top_z() + screw.spec().socket_depth_mm,
        ],
        [ScrewFace::HeadTop, ScrewFace::SocketFloor],
        ScrewFace::SocketWall,
    );
}
