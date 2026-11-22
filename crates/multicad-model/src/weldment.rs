//! Weldment members: the profile section swept along the member axis
//! as an exact prism (hollow sections keep their wall), and member cut
//! records measured from each member Part's final solid in its member frame.
use crate::{GeometryBackend, Model, ModelError, error};
use multicad_core::{
    document::{Document, DocumentError, EntityData, EntityId},
    materials::{self, ProfileShape},
    model::{
        self, ConnectionKind, EndCut, EndTreatment, MemberDefinition, MemberEnd, MemberFrame,
        MemberPath, PathSegment, ResolvedConnection, ResolvedMember, SectionEdge, SectionLoop,
    },
    parameters::EvaluatedParameters,
    references::Origin,
    sketch::{Frame, Geometry},
};
use multicad_kernel::{
    BooleanOperation, CopedTube, CopedTubeNames, CurvedEdge2d, CurvedProfile, GeometryErrorKind,
    GeometryPolicy, PlaneFrame, RigidTransform, Solid,
};
use multicad_precision::PortableMath as _;
use multicad_solver::SketchEvaluationContext;
use std::collections::BTreeMap;

const TOLERANCE_MM: f64 = 1e-3;

/// Solved sketch-line paths, evaluated once per sketch per compile.
pub(crate) struct Paths<'a> {
    context: &'a SketchEvaluationContext<'a>,
    sketches: BTreeMap<EntityId, (Frame, Vec<Geometry>)>,
}
impl<'a> Paths<'a> {
    pub(crate) fn new(context: &'a SketchEvaluationContext<'a>) -> Self {
        Self {
            context,
            sketches: BTreeMap::new(),
        }
    }
    /// The solved world segment of a member path.
    pub(crate) fn segment(
        &mut self,
        member: &MemberDefinition,
        owner: EntityId,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<PathSegment, ModelError> {
        let (sketch, line) = match &member.path {
            MemberPath::SketchLine { sketch, line } => (*sketch, *line),
            path => {
                return model::point_path(path, self.context.parameters())
                    .map_err(|e| error(Some(owner), e))?
                    .ok_or_else(|| error(Some(owner), "Member path is not a point pair"));
            }
        };
        if let std::collections::btree_map::Entry::Vacant(entry) = self.sketches.entry(sketch) {
            let evaluated = self
                .context
                .evaluate(sketch, None, interrupt)
                .map_err(|e| error(Some(owner), format!("Member path sketch: {e}")))?;
            let geometry = evaluated
                .geometry()
                .map_err(|e| error(Some(owner), format!("Member path sketch: {e}")))?;
            entry.insert((evaluated.frame, geometry));
        }
        let (frame, geometry) = &self.sketches[&sketch];
        let found = geometry.iter().find_map(|g| match g {
            Geometry::Line { id, a, b } if *id == line => Some((*a, *b)),
            _ => None,
        });
        let (a, b) = found.ok_or_else(|| {
            error(
                Some(owner),
                "The member's sketch line no longer exists or is not a line; choose another path",
            )
        })?;
        Ok(PathSegment {
            start: frame.point(a),
            end: frame.point(b),
            reference: Some(frame.n),
        })
    }
}

fn edges(section_loop: &SectionLoop) -> Vec<CurvedEdge2d> {
    match section_loop {
        SectionLoop::Polygon(points) => (0..points.len())
            .map(|i| CurvedEdge2d::Line {
                a: points[i],
                b: points[(i + 1) % points.len()],
            })
            .collect(),
        SectionLoop::Curved(edges) => edges
            .iter()
            .map(|edge| match *edge {
                SectionEdge::Line { a, b } => CurvedEdge2d::Line { a, b },
                SectionEdge::Arc { center, a, b } => CurvedEdge2d::Arc {
                    center,
                    a,
                    b,
                    clockwise: false,
                },
            })
            .collect(),
        SectionLoop::Circle { center, radius_mm } => vec![CurvedEdge2d::Circle {
            center: *center,
            radius_mm: *radius_mm,
        }],
    }
}
/// The member section placed at axial distance `s` from the path start.
/// With `envelope`, the solid outside grown by that clearance, without the
/// bore; otherwise the real section with its holes.
pub(crate) fn section_profile(
    member: &ResolvedMember,
    s: f64,
    envelope: Option<f64>,
) -> CurvedProfile {
    let holes = match envelope {
        None => member.section.holes().iter().map(edges).collect(),
        Some(_) => vec![],
    };
    CurvedProfile {
        frame: PlaneFrame {
            origin: member.section_origin(s),
            u: member.frame.x,
            v: member.frame.y,
        },
        outer: edges(&member.section.outer(envelope.unwrap_or(0.))),
        holes,
    }
}
fn origins(
    owner: EntityId,
    member: &ResolvedMember,
    profile: &CurvedProfile,
    prefix: &str,
) -> Vec<Origin> {
    let outer = (0..profile.outer.len()).map(|i| {
        member.section.slotted_outer_role(i).map_or_else(
            || format!("{prefix}.outer.{i}"),
            |role| format!("{prefix}.{role}"),
        )
    });
    let inner = profile
        .holes
        .iter()
        .flatten()
        .enumerate()
        .map(|(i, _)| format!("{prefix}.inner.{i}"));
    outer
        .chain(inner)
        .map(|role| Origin {
            source: owner,
            role,
        })
        .collect()
}
/// An exact prism of the member section from axial distance `from` over
/// `length` along the axis.
pub(crate) fn prism(
    member: &ResolvedMember,
    from: f64,
    length: f64,
    envelope: Option<f64>,
    owner: EntityId,
    prefix: &str,
    policy: GeometryPolicy,
) -> Result<Solid, ModelError> {
    let profile = section_profile(member, from, envelope);
    Solid::extrude_curved_profile_vector_named(
        &profile,
        member.frame.axis.map(|a| a * length),
        owner,
        &origins(owner, member, &profile, prefix),
        policy,
    )
    .map_err(|e| error(Some(owner), e))
}

/// Resolve a member feature with its solved path.
pub(crate) fn resolve(
    document: &Document,
    parameters: &EvaluatedParameters,
    feature: EntityId,
    paths: &mut Paths,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<ResolvedMember, ModelError> {
    let e = |e: DocumentError| error(Some(feature), e);
    let source = model::read_member(document, feature).map_err(e)?;
    let segment = paths.segment(&source, feature, interrupt)?;
    model::resolve_member(document, parameters, &source, segment).map_err(e)
}

/// One connection treating an end of a member.
#[derive(Clone, Debug)]
pub(crate) struct Joint {
    pub connection: EntityId,
    pub partner: EntityId,
    pub partner_member: ResolvedMember,
    pub resolved: ResolvedConnection,
    pub treatment: EndTreatment,
}
/// Every connection that cuts an end of the member Part `part`. Two
/// connections cutting the same end are refused.
pub(crate) fn joints(
    document: &Document,
    parameters: &EvaluatedParameters,
    part: EntityId,
    member: &ResolvedMember,
    other: &mut dyn FnMut(EntityId) -> Result<ResolvedMember, ModelError>,
) -> Result<Vec<Joint>, ModelError> {
    let mut out: Vec<Joint> = vec![];
    for (id, connection) in
        model::member_connections(document, part).map_err(|e| error(Some(part), e))?
    {
        let e = |e: DocumentError| error(Some(id), e);
        let index = usize::from(connection.members[1] == part);
        let partner = connection.members[1 - index];
        let feature = model::member_feature(document, partner).ok_or_else(|| {
            error(
                Some(id),
                "The connected Part is no longer a weldment member",
            )
        })?;
        let partner_member = other(feature)?;
        if member.section.is_item_line_8_40x40() || partner_member.section.is_item_line_8_40x40() {
            return Err(error(
                Some(id),
                "Slotted extrusion members use Mechanical profile connectors; welded end cuts are unavailable",
            ));
        }
        let gap = connection.gap_mm(parameters).map_err(e)?;
        let (a, b) = if index == 0 {
            (member, &partner_member)
        } else {
            (&partner_member, member)
        };
        let resolved = model::resolve_connection(connection.kind, a, b, gap).map_err(e)?;
        let Some(treatment) = resolved.treatments[index] else {
            continue;
        };
        if let Some(previous) = out.iter().find(|j| j.treatment.end == treatment.end) {
            let name = |id: EntityId| document.get(id).map_or("".into(), |e| e.name.clone());
            return Err(error(
                Some(id),
                format!(
                    "{} and {} both cut the {:?} end of {}; remove one connection",
                    name(previous.connection),
                    name(id),
                    treatment.end,
                    name(part)
                )
                .replace("Start end", "start")
                .replace("End end", "end"),
            ));
        }
        out.push(Joint {
            connection: id,
            partner,
            partner_member,
            resolved,
            treatment,
        });
    }
    Ok(out)
}
fn diagonal(m: &ResolvedMember) -> f64 {
    m.section.width_mm.portable_hypot(m.section.height_mm)
}
fn plane_frame(point: [f64; 3], normal: [f64; 3]) -> PlaneFrame {
    let seed = if normal[0].abs() < 0.9 {
        [1., 0., 0.]
    } else {
        [0., 1., 0.]
    };
    let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    let d = dot(seed, normal);
    let u = [
        seed[0] - normal[0] * d,
        seed[1] - normal[1] * d,
        seed[2] - normal[2] * d,
    ];
    let len = dot(u, u).sqrt();
    let u = u.map(|c| c / len);
    let v = [
        normal[1] * u[2] - normal[2] * u[1],
        normal[2] * u[0] - normal[0] * u[2],
        normal[0] * u[1] - normal[1] * u[0],
    ];
    PlaneFrame {
        origin: point,
        u,
        v,
    }
}

/// Remove the half-space behind a plane: material on the side `normal`
/// points to is kept.
fn cut_plane(
    solid: Solid,
    point: [f64; 3],
    normal: [f64; 3],
    size: f64,
    owner: EntityId,
    role: &str,
    policy: GeometryPolicy,
) -> Result<Option<Solid>, ModelError> {
    let corners = [[-size, -size], [size, -size], [size, size], [-size, size]];
    let profile = CurvedProfile {
        frame: plane_frame(point, normal),
        outer: (0..4)
            .map(|i| CurvedEdge2d::Line {
                a: corners[i],
                b: corners[(i + 1) % 4],
            })
            .collect(),
        holes: vec![],
    };
    let origins: Vec<_> = (0..4)
        .map(|i| Origin {
            source: owner,
            role: format!("{role}.side.{i}"),
        })
        .collect();
    let displacement = normal.map(|n| -n * size);
    let tool = if solid.uses_rust_geometry() {
        Solid::rust_extrude_linear_named(&profile, displacement, 0., owner, &origins, policy)
    } else {
        Solid::extrude_curved_profile_vector_named(&profile, displacement, owner, &origins, policy)
    }
    .map_err(|e| error(Some(owner), e))?;
    solid
        .boolean(&tool, BooleanOperation::Subtract, policy)
        .map_err(|e| error(Some(owner), e))
}

fn rust_member(
    member: &ResolvedMember,
    joints: &[Joint],
    id: EntityId,
    policy: GeometryPolicy,
) -> Result<Solid, ModelError> {
    let from = member.start_offset_mm();
    let mut stations = [from, from + member.length_mm()];
    let mut miters = [None; 2];
    for joint in joints {
        let EndCut::Plane { point, normal } = joint.treatment.cut else {
            return rust_clipped_member(member, joints, id, policy);
        };
        let end = usize::from(joint.treatment.end == MemberEnd::End);
        if joint.treatment.kind == ConnectionKind::Miter {
            miters[end] = Some((joint.connection, plane_frame(point, normal)));
        } else {
            stations[end] = member.axial(point);
        }
    }
    let profile = section_profile(member, stations[0], None);
    let origins = origins(id, member, &profile, "member");
    let straight = std::iter::once(&profile.outer)
        .chain(&profile.holes)
        .flatten()
        .all(|edge| matches!(edge, CurvedEdge2d::Line { .. }));
    let Some((connection, _)) = miters.iter().flatten().next().copied() else {
        let height = stations[1] - stations[0];
        if height <= 0. {
            return Err(error(
                Some(joints.first().map_or(id, |joint| joint.connection)),
                "Connection cuts leave no member between the end planes",
            ));
        }
        if straight {
            return Solid::rust_extrude_linear_named(
                &profile,
                member.frame.axis.map(|component| component * height),
                0.,
                id,
                &origins,
                policy,
            )
            .map_err(|cause| error(Some(id), cause));
        }
        return crate::rust_backend::stock_profile(id, &profile, &origins, height, policy);
    };
    if !straight {
        return Err(error(
            Some(connection),
            "Rust weldment miters require a straight-sided section",
        ));
    }
    let planes = std::array::from_fn(|end| {
        miters[end].map_or_else(
            || PlaneFrame {
                origin: member.section_origin(stations[end]),
                ..profile.frame
            },
            |(_, plane)| plane,
        )
    });
    let trimmed = Solid::rust_extrude_between_planes_named(
        &profile,
        member.frame.axis,
        planes,
        id,
        &origins,
        policy,
    );
    match trimmed {
        Ok(solid) => Ok(solid),
        Err(cause) if cause.kind == GeometryErrorKind::Input => {
            // A short member can lose part of its section after the first
            // miter. A generator-preserving sweep cannot represent that
            // intermediate shape; clip the actual extended stock instead.
            rust_clipped_member(member, joints, id, policy)
        }
        Err(cause) => Err(error(Some(connection), cause)),
    }
}

fn cut_stock_extent(member: &ResolvedMember, joints: &[Joint]) -> (f64, f64, f64) {
    let mut from = member.start_offset_mm();
    let mut to = from + member.length_mm();
    let mut size = member.length_mm() + diagonal(member);
    for joint in joints {
        let sin = (joint.resolved.angle_deg.to_radians().portable_sin()).max(0.25);
        let margin = 2. * (diagonal(member) + diagonal(&joint.partner_member)) / sin
            + joint.resolved.gap_mm
            + 1.;
        match joint.treatment.end {
            MemberEnd::Start => from = from.min(-margin),
            MemberEnd::End => to = to.max(member.path_length_mm + margin),
        }
        size = size.max(joint.partner_member.length_mm() + 4. * margin);
    }
    (from, to, 2. * size + 100.)
}

/// The center and radius of a profile loop that is one full circle.
fn circle(profile: &CurvedProfile, edges: &[CurvedEdge2d]) -> Option<([f64; 3], f64)> {
    let [CurvedEdge2d::Circle { center, radius_mm }] = edges else {
        return None;
    };
    let f = &profile.frame;
    let world = std::array::from_fn(|i| f.origin[i] + center[0] * f.u[i] + center[1] * f.v[i]);
    Some((world, *radius_mm))
}

/// A round hollow branch with one cope onto a round receiving member and a
/// square far end: the exact fish-mouth boundary. `None` leaves every other
/// case to the clipped stock.
fn rust_round_cope(
    member: &ResolvedMember,
    joints: &[Joint],
    id: EntityId,
    policy: GeometryPolicy,
) -> Option<Result<Solid, ModelError>> {
    let [joint] = joints else {
        return None;
    };
    let EndCut::Cope { grow_mm, .. } = joint.treatment.cut else {
        return None;
    };
    let partner = &joint.partner_member;
    let from = member.start_offset_mm();
    let to = from + member.length_mm();
    let profile = section_profile(member, from, None);
    let (center, outer) = circle(&profile, &profile.outer)?;
    let [hole] = profile.holes.as_slice() else {
        return None;
    };
    let (hole_center, inner) = circle(&profile, hole)?;
    let envelope = section_profile(partner, partner.start_offset_mm(), Some(grow_mm));
    let (point, radius) = circle(&envelope, &envelope.outer)?;
    let magnitude = center
        .into_iter()
        .chain(point)
        .map(f64::abs)
        .fold(outer.max(radius), f64::max);
    let rounding = 64. * f64::EPSILON * magnitude;
    if distance(center, hole_center) > rounding {
        return None;
    }
    // Where the branch axis meets the receiving axis.
    let (a, b) = (member.frame.axis, partner.frame.axis);
    let w = sub(point, center);
    let (ab, aw, bw) = (dot(a, b), dot(a, w), dot(b, w));
    let denominator = 1. - ab * ab;
    if denominator <= 1e-6 {
        return None;
    }
    let s = (aw - ab * bw) / denominator;
    let u = (ab * aw - bw) / denominator;
    let on_branch = add(center, scale(a, s));
    if distance(on_branch, add(point, scale(b, u))) > rounding {
        return None;
    }
    // The finite envelope must span all the material the cope removes:
    // between the plane through the receiving axis and the saddle, that
    // material projects onto the receiving axis within
    // (R + |cos| Rc) / sin of the axes' meeting point.
    let sin = denominator.sqrt();
    let reach = (outer + ab.abs() * radius) / sin + policy.linear_tolerance_mm;
    let along = partner.axial(on_branch);
    let span = partner.start_offset_mm();
    if along - reach < span || along + reach > span + partner.length_mm() {
        return None;
    }
    let (axis, length) = match joint.treatment.end {
        MemberEnd::Start => (a, to - from - s),
        MemberEnd::End => (scale(a, -1.), s - (to - from)),
    };
    let tube = CopedTube {
        origin_mm: on_branch,
        axis,
        outer_radius_mm: outer,
        inner_radius_mm: inner,
        length_mm: length,
        envelope_axis: b,
        envelope_radius_mm: radius,
    };
    let roles = origins(id, member, &profile, "member");
    let [outer_role, inner_role] = [&roles[0], &roles[1]];
    let role = |role: String| Origin { source: id, role };
    let names = CopedTubeNames {
        outer: outer_role.clone(),
        inner: inner_role.clone(),
        end: role(match joint.treatment.end {
            MemberEnd::Start => "member.end".into(),
            MemberEnd::End => "member.start".into(),
        }),
        saddle: role(format!("connection.{}.cope", joint.connection)),
    };
    Some(
        Solid::rust_coped_round_tube(&tube, &names, id, policy)
            .map_err(|cause| error(Some(joint.connection), cause)),
    )
}

fn rust_clipped_member(
    member: &ResolvedMember,
    joints: &[Joint],
    id: EntityId,
    policy: GeometryPolicy,
) -> Result<Solid, ModelError> {
    if let Some(coped) = rust_round_cope(member, joints, id, policy) {
        return coped;
    }
    let (from, to, size) = cut_stock_extent(member, joints);
    let profile = section_profile(member, from, None);
    let straight = |profile: &CurvedProfile| {
        profile
            .outer
            .iter()
            .chain(profile.holes.iter().flatten())
            .all(|edge| matches!(edge, CurvedEdge2d::Line { .. }))
    };
    if !straight(&profile) {
        return Err(error(
            Some(
                joints
                    .iter()
                    .find(|joint| matches!(joint.treatment.cut, EndCut::Cope { .. }))
                    .map_or(id, |joint| joint.connection),
            ),
            "Rust coped ends require straight-sided member sections",
        ));
    }
    let mut solid = Solid::rust_extrude_linear_named(
        &profile,
        member.frame.axis.map(|a| a * (to - from)),
        0.,
        id,
        &origins(id, member, &profile, "member"),
        policy,
    )
    .map_err(|cause| error(Some(id), cause))?;
    for joint in joints {
        let (point, normal) = match joint.treatment.cut {
            EndCut::Plane { point, normal } => (point, normal),
            EndCut::Cope {
                grow_mm,
                point,
                normal,
            } => {
                let partner = &joint.partner_member;
                let profile = section_profile(partner, partner.start_offset_mm(), Some(grow_mm));
                if !straight(&profile) {
                    return Err(error(
                        Some(joint.connection),
                        "Rust coped ends require a straight-sided receiving section",
                    ));
                }
                let envelope = Solid::rust_extrude_linear_named(
                    &profile,
                    partner.frame.axis.map(|a| a * partner.length_mm()),
                    0.,
                    id,
                    &origins(
                        id,
                        partner,
                        &profile,
                        &format!("connection.{}.cope", joint.connection),
                    ),
                    policy,
                )
                .map_err(|cause| error(Some(joint.connection), cause))?;
                solid = solid
                    .boolean(&envelope, BooleanOperation::Subtract, policy)
                    .map_err(|cause| error(Some(joint.connection), cause))?
                    .ok_or_else(|| {
                        error(
                            Some(joint.connection),
                            "The connection cut removes the whole member",
                        )
                    })?;
                (point, normal)
            }
        };
        solid = cut_plane(
            solid,
            point,
            normal,
            size,
            id,
            &format!("connection.{}", joint.connection),
            policy,
        )?
        .ok_or_else(|| {
            error(
                Some(joint.connection),
                "The connection cut removes the whole member",
            )
        })?;
    }
    if solid
        .solid_count()
        .map_err(|cause| error(Some(id), cause))?
        != 1
    {
        return Err(error(
            Some(id),
            "Connection cuts split the member into separate pieces; check the joint geometry",
        ));
    }
    Ok(solid)
}

/// Build a member solid with its connection end cuts; returns it with its
/// resolved member.
pub(crate) fn build(
    document: &Document,
    parameters: &EvaluatedParameters,
    id: EntityId,
    policy: GeometryPolicy,
    backend: GeometryBackend,
    paths: &mut Paths,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<(Solid, ResolvedMember), ModelError> {
    let member = resolve(document, parameters, id, paths, interrupt)?;
    let source = model::read_member(document, id).map_err(|e| error(Some(id), e))?;
    let joints = joints(
        document,
        parameters,
        source.output.body,
        &member,
        &mut |feature| resolve(document, parameters, feature, paths, interrupt),
    )?;
    if backend == GeometryBackend::Rust {
        return Ok((rust_member(&member, &joints, id, policy)?, member));
    }
    // Connected ends start long enough to be cut back to their joint.
    let (from, to, size) = cut_stock_extent(&member, &joints);
    let mut solid = prism(&member, from, to - from, None, id, "member", policy)?;
    for joint in &joints {
        let role = format!("connection.{}", joint.connection);
        let removed = || {
            error(
                Some(joint.connection),
                "The connection cut removes the whole member",
            )
        };
        solid = match joint.treatment.cut {
            EndCut::Plane { point, normal } => {
                cut_plane(solid, point, normal, size, id, &role, policy)?.ok_or_else(removed)?
            }
            EndCut::Cope {
                grow_mm,
                point,
                normal,
            } => {
                let partner = &joint.partner_member;
                let envelope = prism(
                    partner,
                    partner.start_offset_mm(),
                    partner.length_mm(),
                    Some(grow_mm),
                    id,
                    &format!("{role}.cope"),
                    policy,
                )?;
                let coped = solid
                    .boolean(&envelope, BooleanOperation::Subtract, policy)
                    .map_err(|e| error(Some(joint.connection), e))?
                    .ok_or_else(removed)?;
                cut_plane(coped, point, normal, size, id, &role, policy)?.ok_or_else(removed)?
            }
        };
    }
    if solid.solid_count().map_err(|e| error(Some(id), e))? != 1 {
        return Err(error(
            Some(id),
            "Connection cuts split the member into separate pieces; check the joint geometry",
        ));
    }
    Ok((solid, member))
}

/// The world-to-member transform: x, y and the axis become X, Y and Z, with
/// the path start at the origin.
pub(crate) fn into_member(frame: &MemberFrame, start: [f64; 3]) -> RigidTransform {
    let d = |a: [f64; 3]| -(a[0] * start[0] + a[1] * start[1] + a[2] * start[2]);
    RigidTransform {
        rotation: [frame.x, frame.y, frame.axis],
        translation_mm: [d(frame.x), d(frame.y), d(frame.axis)],
    }
}

/// How one end of a member is finished.
#[derive(Clone, Debug, PartialEq)]
pub struct EndRecord {
    /// None for a square end without a connection.
    pub kind: Option<model::ConnectionKind>,
    pub connection: Option<EntityId>,
    /// The member Part at the other side of the connection.
    pub partner: Option<EntityId>,
    /// Plane cuts: the cut face angle from square (0 = square, 45 = a
    /// right-angle miter). Copes: the branch angle to the receiving axis.
    pub angle_deg: f64,
    /// False for a continuing member cut square to cover its partner.
    pub terminates: bool,
    /// Plane cuts more than 0.01° from square: where the cut end reaches
    /// farthest along the axis. None for square ends and copes.
    pub long_point: Option<LongPoint>,
}
impl EndRecord {
    fn square() -> Self {
        Self {
            kind: None,
            connection: None,
            partner: None,
            angle_deg: 0.,
            terminates: true,
            long_point: None,
        }
    }
}
/// The long point of an angled plane cut: the side of the section where the
/// cut end reaches farthest along the member axis.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LongPoint {
    /// World unit direction from the axis toward the long point, square to
    /// the axis.
    pub direction: [f64; 3],
    /// The section face it lies on, in section coordinates (width x, height
    /// y): "+x", "-x", "+y" or "-y", or a corner such as "+x-y" when the cut
    /// runs diagonally across the section.
    pub face: &'static str,
}
impl LongPoint {
    /// The long point of a plane cut that keeps the side `normal` points to.
    /// For a section offset q the cut face lies at s(q) = -(n_q · q) / n_a
    /// along the axis, so at either end the material reaches farthest on the
    /// side of the normal's component square to the axis.
    pub(crate) fn of_plane(frame: &MemberFrame, normal: [f64; 3]) -> Option<Self> {
        let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        let along = dot(normal, frame.axis);
        let q = [0, 1, 2].map(|i| normal[i] - along * frame.axis[i]);
        let length = dot(q, q).sqrt();
        if length < 0.01_f64.to_radians().portable_sin() {
            return None;
        }
        let direction = q.map(|c| c / length);
        let (x, y) = (dot(direction, frame.x), dot(direction, frame.y));
        let face = if (x.abs() - y.abs()).abs() < 1e-6 {
            match (x > 0., y > 0.) {
                (true, true) => "+x+y",
                (true, false) => "+x-y",
                (false, true) => "-x+y",
                (false, false) => "-x-y",
            }
        } else if x.abs() > y.abs() {
            if x > 0. { "+x" } else { "-x" }
        } else if y > 0. {
            "+y"
        } else {
            "-y"
        };
        Some(Self { direction, face })
    }
}

/// One member's cut record, measured from its current Part solid.
#[derive(Clone, Debug, PartialEq)]
pub struct MemberRecord {
    pub part: EntityId,
    pub feature: EntityId,
    pub name: String,
    pub material: Option<EntityId>,
    pub stock: EntityId,
    pub stock_name: String,
    pub designation: String,
    pub shape: ProfileShape,
    /// Section bounding box: width (x) and height (y).
    pub section_mm: [f64; 2],
    pub wall_mm: Option<f64>,
    /// Sketch line or point pair the member follows.
    pub path: MemberPath,
    /// Path length plus both end extensions.
    pub nominal_length_mm: f64,
    /// Overall length along the axis, measured from the solid: the cut
    /// length from the farthest point of one end to the farthest point of
    /// the other.
    pub length_mm: f64,
    /// Measured axial range of the solid from the path start.
    pub span_mm: [f64; 2],
    /// The start and end treatments.
    pub ends: [EndRecord; 2],
    /// Measured section size in the member frame.
    pub measured_section_mm: [f64; 2],
    /// Section rotation about the axis, and the anchor on the path.
    pub rotation_deg: f64,
    pub anchor: [model::AxisAnchor; 2],
    /// World directions of the member axis and the section width and height.
    pub frame: MemberFrame,
    pub path_start: [f64; 3],
    pub volume_mm3: f64,
    /// Whether the overall length fits one stock bar.
    pub fits_stock: bool,
    /// Review findings: section or stock assignment that disagree.
    pub issues: Vec<String>,
}

impl Model {
    /// Cut records for every evaluated member Part, in document order. The
    /// length convention is the overall length along the member axis.
    pub fn member_cut_list(&self, document: &Document) -> Result<Vec<MemberRecord>, ModelError> {
        let parameters =
            multicad_core::parameters::evaluate(document).map_err(|e| error(None, e))?;
        let mut out = vec![];
        for entity in model::members(document) {
            let id = entity.id;
            let e = |e: DocumentError| error(Some(id), e);
            let source = model::read_member(document, id).map_err(e)?;
            let part = source.output.body;
            let Some(segment) = self.snapshot().member_paths.get(&id) else {
                continue;
            };
            let Some(solid) = self
                .snapshot()
                .bodies
                .get(&part)
                .and_then(|head| self.solid(*head))
            else {
                continue;
            };
            let member =
                model::resolve_member(document, &parameters, &source, *segment).map_err(e)?;
            let b = solid
                .transformed(into_member(&member.frame, member.path.start))
                .map_err(|e| error(Some(id), e))?
                .bounds()
                .map_err(|e| error(Some(id), e))?;
            let length = b.max[2] - b.min[2];
            let snapshot = self.snapshot();
            let joints = joints(document, &parameters, part, &member, &mut |feature| {
                let e = |e: DocumentError| error(Some(feature), e);
                let other = model::read_member(document, feature).map_err(e)?;
                let segment = snapshot
                    .member_paths
                    .get(&feature)
                    .ok_or_else(|| error(Some(feature), "Connected member was not evaluated"))?;
                model::resolve_member(document, &parameters, &other, *segment).map_err(e)
            })?;
            let mut ends = [EndRecord::square(), EndRecord::square()];
            for j in &joints {
                ends[usize::from(j.treatment.end == model::MemberEnd::End)] = EndRecord {
                    kind: Some(j.treatment.kind),
                    connection: Some(j.connection),
                    partner: Some(j.partner),
                    angle_deg: j.treatment.angle_deg,
                    terminates: j.treatment.terminates,
                    long_point: match j.treatment.cut {
                        EndCut::Plane { normal, .. } if j.treatment.angle_deg > 0.01 => {
                            LongPoint::of_plane(&member.frame, normal)
                        }
                        _ => None,
                    },
                };
            }
            let measured = [b.max[0] - b.min[0], b.max[1] - b.min[1]];
            let section = &member.section;
            let mut issues = vec![];
            if (measured[0] - section.width_mm).abs() > TOLERANCE_MM
                || (measured[1] - section.height_mm).abs() > TOLERANCE_MM
            {
                issues.push(format!(
                    "Measured section {:.3} × {:.3} mm differs from the {} × {} mm stock",
                    measured[0], measured[1], section.width_mm, section.height_mm
                ));
            }
            if let Some((_, properties)) = materials::part_properties(document, part).map_err(e)?
                && properties.stock.as_ref().map(|s| s.id) != Some(member.stock)
            {
                issues.push("The Part's stock assignment differs from its member stock".into());
            }
            let Some(EntityData::Part { material, .. }) = document.get(part).map(|e| &e.data)
            else {
                continue;
            };
            out.push(MemberRecord {
                part,
                feature: id,
                name: document
                    .get(part)
                    .map(|e| e.name.clone())
                    .unwrap_or_default(),
                material: *material,
                stock: member.stock,
                stock_name: document
                    .get(member.stock)
                    .map(|e| e.name.clone())
                    .unwrap_or_default(),
                designation: section.designation.clone(),
                shape: section.shape,
                section_mm: [section.width_mm, section.height_mm],
                wall_mm: section.wall_mm,
                path: source.path.clone(),
                nominal_length_mm: member.length_mm(),
                length_mm: length,
                span_mm: [b.min[2], b.max[2]],
                ends,
                measured_section_mm: measured,
                rotation_deg: member.rotation_rad.to_degrees(),
                anchor: source.anchor,
                frame: member.frame,
                path_start: member.path.start,
                volume_mm3: solid.volume_mm3().map_err(|e| error(Some(id), e))?,
                fits_stock: length <= section.stock_length_mm + TOLERANCE_MM,
                issues,
            });
        }
        Ok(out)
    }
}

fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|i| a[i] + b[i])
}
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|i| a[i] - b[i])
}
fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    a.map(|x| x * s)
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    (0..3).map(|i| a[i] * b[i]).sum()
}
fn distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    let d = sub(a, b);
    dot(d, d).sqrt()
}
