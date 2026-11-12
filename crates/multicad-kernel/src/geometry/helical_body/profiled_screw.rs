//! Policy admission and finite face identities for actual profiled heads.
pub(super) mod torx;

use super::{
    boundary::{Boundary, placed_frame},
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
        Self::rust_profiled_machine_screw_impl(
            thread,
            root_radius_mm,
            screw,
            feature,
            underhead,
            policy,
            false,
            None,
        )
    }

    /// Retain the original thread phase while the true head covers its prefix.
    /// Exchange remains gated until finite trimmed-surface readback is available.
    pub fn rust_profiled_machine_screw_with_head_thread_union(
        thread: HelicalThread,
        root_radius_mm: Option<f64>,
        screw: ScrewSpec,
        feature: EntityId,
        underhead: PlaneFrame,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        Self::rust_profiled_machine_screw_impl(
            thread,
            root_radius_mm,
            screw,
            feature,
            underhead,
            policy,
            true,
            None,
        )
    }

    /// Preserve the nominal shank and the fitted thread, including head overlap.
    pub fn rust_profiled_machine_screw_with_shank(
        thread: HelicalThread,
        root_radius_mm: Option<f64>,
        shank_radius_mm: f64,
        screw: ScrewSpec,
        feature: EntityId,
        underhead: PlaneFrame,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        policy.length(shank_radius_mm)?;
        Self::rust_profiled_machine_screw_impl(
            thread,
            root_radius_mm,
            screw,
            feature,
            underhead,
            policy,
            true,
            Some(shank_radius_mm),
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn rust_profiled_machine_screw_impl(
        thread: HelicalThread,
        root_radius_mm: Option<f64>,
        screw: ScrewSpec,
        feature: EntityId,
        underhead: PlaneFrame,
        policy: GeometryPolicy,
        union: bool,
        shank_radius_mm: Option<f64>,
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
            let shank = shank_radius_mm.unwrap_or(thread.spec().major_radius_mm);
            let screw = if union {
                Screw::new_with_shank_thread_union(thread, shank, screw)
            } else {
                Screw::new_with_shank_radius(thread, shank, screw)
            }
            .map_err(|cause| match cause {
                ThreadError::Invalid => {
                    input("Invalid profiled screw head, neck, shank or hex drive")
                }
                ThreadError::Limit => error(
                    OPERATION,
                    GeometryErrorKind::Limit,
                    "Profiled screw exceeds the geometry envelope",
                ),
            })?;
            qualify_source(&screw, frame, policy)?;
            Ok(HelicalForm::ProfiledScrew(Box::new(screw)))
        })
    }

    /// Admit an already validated physical source in its original thread frame.
    /// Retained meridians, source uncertainty and absent authoring metadata are
    /// preserved; no head recipe or placement is inferred.
    pub fn rust_profiled_screw_from_source(
        source: Screw,
        feature: EntityId,
        thread_frame: PlaneFrame,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        policy.validate()?;
        thread_frame.checked(policy)?;
        let geometry = *source.thread();
        let spec = geometry.spec();
        checked_thread(
            HelicalThread {
                major_radius_mm: spec.major_radius_mm,
                pitch_mm: spec.pitch_mm,
                length_mm: spec.length_mm,
                depth_mm: spec.depth_mm,
                crest_width_mm: spec.crest_width_mm,
                lower_flank_deg: spec.lower_flank_deg,
                upper_flank_deg: spec.upper_flank_deg,
                radial_taper: 0.,
                starts: spec.starts,
                left_hand: spec.left_hand,
            },
            spec.root_radius_mm,
            policy,
        )?;
        Self::publish_helical_geometry(geometry, feature, thread_frame, policy, |_| {
            qualify_source(&source, thread_frame, policy)?;
            Ok(HelicalForm::ProfiledScrew(Box::new(source)))
        })
    }
}

fn qualify_source(
    screw: &Screw,
    frame: PlaneFrame,
    policy: GeometryPolicy,
) -> Result<(), GeometryError> {
    let (drive, socket_depth) = screw.recess().flat_hex_geometry().ok_or_else(|| {
        error(
            OPERATION,
            GeometryErrorKind::Unsupported,
            "Conical curved recesses need the qualified native mesh and exchange pipeline",
        )
    })?;
    let spec = screw.spec();
    policy.length(drive.across_flats_mm())?;
    policy.length(socket_depth)?;
    policy.length(spec.minimum_wall_mm)?;
    for (i, [x, y]) in drive.corners_mm().into_iter().enumerate() {
        policy.point(position(frame, [x, y, screw.head_top_z()]))?;
        policy.length(
            (drive.side_length_mm(i) - 2. * screw.profile().rounding_error_mm()).next_down(),
        )?;
    }
    if spec.shank_length_mm != 0. {
        policy.length(spec.shank_length_mm)?;
    }
    let reserve = screw.profile().rounding_error_mm();
    let shank = screw.profile().end()[0];
    let major = screw.thread().spec().major_radius_mm;
    if shank > major {
        policy.length(((shank - major).next_down() - 2. * reserve).next_down())?;
    }
    if screw.has_head_thread_overlap() {
        let [start, end] = screw.visible_thread_interval_mm();
        policy.length(((end - start).next_down() - 2. * reserve).next_down())?;
    }
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
            policy.point(position(
                frame,
                [p[0], 0., p[1] - screw.spec().shank_length_mm],
            ))?;
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
                (segment.end[0] - segment.start[0]).hypot(segment.end[1] - segment.start[1]),
            )?;
        }
    }
    Ok(())
}

pub(super) fn extend_face_properties(
    screw: &Screw,
    faces: &mut Vec<(Point, f64)>,
    policy: GeometryPolicy,
) -> Result<(), GeometryError> {
    let Some((drive, socket_depth)) = screw.recess().flat_hex_geometry() else {
        return torx::extend_face_properties(screw, faces, policy);
    };
    let spec = screw.spec();
    let radius = screw.profile().end()[0];
    let (removed_center, removed_area) = faces[0];
    let shoulder = PI * radius.powi(2) - removed_area;
    faces[0] = (
        [
            -removed_center[0] * removed_area / shoulder,
            -removed_center[1] * removed_area / shoulder,
            removed_center[2],
        ],
        shoulder,
    );
    let top = screw.head_top_z();
    let socket_area = drive.area_mm2();
    let [cx, cy] = drive.center_mm();
    let segments: Vec<_> = screw.profile().segments().collect();
    for face in screw.faces().into_iter().skip(faces.len()) {
        faces.push(match face {
            ScrewFace::HeadTop => {
                let area = PI * screw.profile().top()[0].powi(2) - socket_area;
                (
                    [-cx * socket_area / area, -cy * socket_area / area, top],
                    area,
                )
            }
            ScrewFace::ProfileWall(i) => {
                let segment = segments[i as usize];
                let mut center = segment.center_mm();
                center[2] -= spec.shank_length_mm;
                (center, segment.area_mm2())
            }
            ScrewFace::SocketFloor => ([cx, cy, top + socket_depth], socket_area),
            ScrewFace::SocketWall(i) => {
                let [a, b] = drive.edge_mm(i as usize);
                (
                    [
                        a[0].midpoint(b[0]),
                        a[1].midpoint(b[1]),
                        top + socket_depth / 2.,
                    ],
                    drive.side_length_mm(i as usize) * socket_depth,
                )
            }
            _ => unreachable!("Original helical face properties already exist"),
        });
    }
    Ok(())
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
        ScrewFace::SocketFloor => {
            let (_, socket_depth) = screw.recess().flat_hex_geometry()?;
            Some(axial(top + socket_depth, -1.))
        }
        ScrewFace::SocketWall(i) => {
            let (drive, _) = screw.recess().flat_hex_geometry()?;
            let [a, b] = drive.edge_mm(i as usize);
            let d = [b[0] - a[0], b[1] - a[1]];
            let length = d[0].hypot(d[1]);
            Some(placed_frame(
                frame,
                [a[0], a[1], top],
                [-d[0] / length, -d[1] / length, 0.],
                [0., 0., 1.],
            ))
        }
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
    let start = screw.visible_thread_interval_mm()[0];
    if screw.profile().end()[0] > screw.thread().spec().major_radius_mm {
        add.arc(
            "screw.head/thread.junction/outer".into(),
            start,
            screw.profile().end()[0],
            [0., TAU],
            &[shoulder, last],
        );
    } else {
        for i in 0..count {
            if (body.band(i) - 2) % 4 == 1 {
                let edge = body.rim_origin(0, i);
                add.face(shoulder).boundary.remove(&edge);
                add.face(last).boundary.insert(edge);
            } else {
                let phase = body.rails[i];
                add.arc(
                    format!("screw.head/thread.junction/{}", body.rail_name(i)),
                    start,
                    screw.thread().spec().major_radius_mm,
                    [(start - phase) / rate, (phase - body.next_phase(i)) / rate],
                    &[shoulder, last],
                );
            }
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
    if let Some((drive, socket_depth)) = screw.recess().flat_hex_geometry() {
        add.polygon_prism(
            "screw.drive",
            drive.corners_mm(),
            [screw.head_top_z(), screw.head_top_z() + socket_depth],
            [ScrewFace::HeadTop, ScrewFace::SocketFloor],
            ScrewFace::SocketWall,
        );
    } else {
        torx::extend_boundary(screw, &mut add);
    }
}
