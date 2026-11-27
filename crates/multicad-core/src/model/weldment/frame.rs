//! Frames drawn as lines (CAD-439, design `docs/architecture/frame-drawing.md`).
//!
//! A `weldment.frame` extension holds intent only: a base plane, nodes with
//! position definitions, members between nodes, and per-node joint
//! overrides. Every frame member is an ordinary `weldment.member` feature
//! with its own Part whose path is `MemberPath::FrameMember`, so cut lists,
//! welds and drawings work unchanged. Member feature and Part ids are created
//! once and kept through every edit of the frame.
//!
//! Joints are not stored: each node derives its joint from the current
//! geometry and sections (§5.4) and feeds the existing miter, butt and cope
//! builders. An override is stored only when the user changes a node; an
//! override that no longer applies is reported at the node, never replaced.
use super::*;
use crate::{model::common::number, parameters::EvaluatedParameters};
use multicad_precision::PortableMath as _;
use resolve::vec3::*;
use std::collections::BTreeMap;

pub const FRAME_KIND: &str = "weldment.frame";
pub const FRAME_CAPABILITY: &str = "weldment.frame.1";
pub const MAX_FRAME_NODES: usize = 1024;
/// Directions closer than this are parallel (collinear chains, legs).
const PARALLEL: f64 = 1e-9;

/// The frame's base plane. Node coordinates `(u, v, w)` are measured along
/// `u`, `v = normal × u` and `normal` from `origin`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FramePlane {
    pub origin: [Dimension; 3],
    pub u: [f64; 3],
    pub normal: [f64; 3],
}
impl FramePlane {
    /// The world XY plane at the origin.
    pub fn xy() -> Self {
        let zero = || super::super::length(0.);
        Self {
            origin: [zero(), zero(), zero()],
            u: [1., 0., 0.],
            normal: [0., 0., 1.],
        }
    }
    pub fn v(&self) -> [f64; 3] {
        cross(self.normal, self.u)
    }
}

/// Where a node is. Positions are evaluated in dependency order; a cycle is
/// refused.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum NodePosition {
    /// Plane coordinates; `w` along the plane normal, 0 when omitted.
    Absolute {
        u: Dimension,
        v: Dimension,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        w: Option<Dimension>,
    },
    /// From another node: `length` along the in-plane `angle` (from +u
    /// toward +v), then `rise` along the plane normal.
    Polar {
        from: EntityId,
        length: Dimension,
        angle: Dimension,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        rise: Option<Dimension>,
    },
    /// A station on a member: `distance` from its start, or from its end.
    Along {
        member: EntityId,
        distance: Dimension,
        #[serde(default, skip_serializing_if = "std::ops::Not::not")]
        from_end: bool,
    },
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FrameNode {
    pub id: EntityId,
    pub position: NodePosition,
}
/// Which side of its line a member's section lies on (§3.3). `Center` uses
/// the member's own grid anchor; `Inner` and `Outer` look toward or away
/// from the frame centroid; `Left` and `Right` look along the member with
/// the plane normal up.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum FrameSide {
    #[default]
    Center,
    Inner,
    Outer,
    Left,
    Right,
}
impl FrameSide {
    fn is_center(&self) -> bool {
        *self == Self::Center
    }
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FrameMember {
    pub id: EntityId,
    pub start: EntityId,
    pub end: EntityId,
    /// None follows the frame default side.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub side: Option<FrameSide>,
}
/// A stored change to one node's derived joint. Every field left empty
/// stays automatic.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct JointOverride {
    pub node: EntityId,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub kind: Option<ConnectionKind>,
    /// The member that runs through the node.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub through: Option<EntityId>,
    /// Weld gap at this node; the frame default otherwise.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gap: Option<Dimension>,
}
impl JointOverride {
    pub fn is_automatic(&self) -> bool {
        self.kind.is_none() && self.through.is_none() && self.gap.is_none()
    }
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FrameDefinition {
    pub plane: FramePlane,
    pub nodes: Vec<FrameNode>,
    pub members: Vec<FrameMember>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub joints: Vec<JointOverride>,
    /// Default side of every member that does not set its own.
    #[serde(default, skip_serializing_if = "FrameSide::is_center")]
    pub side: FrameSide,
    /// Frame-wide weld gap; 0 mm when omitted (Q9).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gap: Option<Dimension>,
}

/// One segment of a drawn chain: in-plane length and angle, then rise
/// along the plane normal.
#[derive(Clone, Debug, PartialEq)]
pub struct ChainSegment {
    pub length: Dimension,
    pub angle: Dimension,
    pub rise: Option<Dimension>,
}

impl FrameDefinition {
    pub fn new(plane: FramePlane) -> Self {
        Self {
            plane,
            nodes: vec![],
            members: vec![],
            joints: vec![],
            side: FrameSide::Center,
            gap: None,
        }
    }
    pub fn node(&self, id: EntityId) -> Option<&FrameNode> {
        self.nodes.iter().find(|n| n.id == id)
    }
    pub fn member(&self, id: EntityId) -> Option<&FrameMember> {
        self.members.iter().find(|m| m.id == id)
    }
    pub fn joint(&self, node: EntityId) -> Option<&JointOverride> {
        self.joints.iter().find(|j| j.node == node)
    }
    fn dimensions(&self) -> Vec<(&Dimension, Unit)> {
        let mut out: Vec<(&Dimension, Unit)> = self
            .plane
            .origin
            .iter()
            .map(|d| (d, Unit::Millimeter))
            .collect();
        for node in &self.nodes {
            match &node.position {
                NodePosition::Absolute { u, v, w } => {
                    out.extend([u, v].into_iter().chain(w).map(|d| (d, Unit::Millimeter)));
                }
                NodePosition::Polar {
                    length,
                    angle,
                    rise,
                    ..
                } => {
                    out.push((length, Unit::Millimeter));
                    out.push((angle, Unit::Degree));
                    out.extend(rise.iter().map(|d| (d, Unit::Millimeter)));
                }
                NodePosition::Along { distance, .. } => out.push((distance, Unit::Millimeter)),
            }
        }
        out.extend(
            self.joints
                .iter()
                .filter_map(|j| j.gap.as_ref())
                .chain(&self.gap)
                .map(|d| (d, Unit::Millimeter)),
        );
        out
    }
    /// Structural checks that need no evaluation. Joint overrides naming a
    /// node or member that is gone stay valid source: they are reported at
    /// evaluation instead of being dropped.
    pub fn validate(&self) -> Result<(), DocumentError> {
        let unit_ok = |v: [f64; 3]| v.iter().all(|c| c.is_finite()) && (norm(v) - 1.).abs() < 1e-9;
        if !unit_ok(self.plane.u)
            || !unit_ok(self.plane.normal)
            || dot(self.plane.u, self.plane.normal).abs() > 1e-9
        {
            return Err(invalid(
                "A frame plane needs unit u and normal directions square to each other",
            ));
        }
        if self.nodes.len() > MAX_FRAME_NODES || self.members.len() > MAX_FRAME_MEMBERS {
            return Err(invalid(format!(
                "A frame holds at most {MAX_FRAME_NODES} nodes and {MAX_FRAME_MEMBERS} members"
            )));
        }
        let mut ids = BTreeSet::new();
        for id in self
            .nodes
            .iter()
            .map(|n| n.id)
            .chain(self.members.iter().map(|m| m.id))
        {
            if id.is_nil() || !ids.insert(id) {
                return Err(invalid(
                    "Frame node and member ids must be unique and non-nil",
                ));
            }
        }
        for member in &self.members {
            if member.start == member.end {
                return Err(invalid("A frame member needs two different nodes"));
            }
            if self.node(member.start).is_none() || self.node(member.end).is_none() {
                return Err(invalid("A frame member must join two nodes of its frame"));
            }
        }
        for node in &self.nodes {
            match &node.position {
                NodePosition::Polar { from, .. } if self.node(*from).is_none() => {
                    return Err(invalid("A polar node must start from a node of its frame"));
                }
                NodePosition::Along { member, .. } if self.member(*member).is_none() => {
                    return Err(invalid(
                        "A station node must lie on a member of its frame; delete the station first",
                    ));
                }
                _ => {}
            }
        }
        let mut joints = BTreeSet::new();
        for joint in &self.joints {
            if !joints.insert(joint.node) {
                return Err(invalid("A frame node has at most one joint override"));
            }
        }
        for (d, unit) in self.dimensions() {
            check_dimension(d, unit)?;
        }
        Ok(())
    }
    pub fn inputs(&self) -> Vec<Link> {
        parameter_links(self.dimensions().into_iter().map(|(d, _)| d).collect())
    }
    pub(crate) fn entity_data(&self) -> Result<EntityData, DocumentError> {
        self.validate()?;
        Ok(EntityData::Extension {
            definition: definition(FRAME_KIND, self)?,
            inputs: self.inputs(),
        })
    }

    /// World positions of every node.
    pub fn evaluate(
        &self,
        parameters: &EvaluatedParameters,
    ) -> Result<EvaluatedFrame, DocumentError> {
        self.validate()?;
        let mm = |d: &Dimension| number(d, Unit::Millimeter, parameters);
        let origin = [
            mm(&self.plane.origin[0])?,
            mm(&self.plane.origin[1])?,
            mm(&self.plane.origin[2])?,
        ];
        let (u, n) = (self.plane.u, self.plane.normal);
        let v = self.plane.v();
        let mut done: BTreeMap<EntityId, [f64; 3]> = BTreeMap::new();
        let mut visiting = BTreeSet::new();
        fn place(
            frame: &FrameDefinition,
            id: EntityId,
            axes: ([f64; 3], [f64; 3], [f64; 3], [f64; 3]),
            parameters: &EvaluatedParameters,
            done: &mut BTreeMap<EntityId, [f64; 3]>,
            visiting: &mut BTreeSet<EntityId>,
        ) -> Result<[f64; 3], DocumentError> {
            if let Some(p) = done.get(&id) {
                return Ok(*p);
            }
            if !visiting.insert(id) {
                return Err(invalid(
                    "Frame node positions depend on each other in a cycle; place one of them absolutely",
                ));
            }
            let (origin, u, v, n) = axes;
            let mm = |d: &Dimension| number(d, Unit::Millimeter, parameters);
            let node = frame
                .node(id)
                .ok_or_else(|| invalid("A frame position names a node that does not exist"))?;
            let p = match &node.position {
                NodePosition::Absolute { u: a, v: b, w } => add(
                    origin,
                    add(
                        add(scale(u, mm(a)?), scale(v, mm(b)?)),
                        scale(n, w.as_ref().map_or(Ok(0.), mm)?),
                    ),
                ),
                NodePosition::Polar {
                    from,
                    length,
                    angle,
                    rise,
                } => {
                    let base = place(frame, *from, axes, parameters, done, visiting)?;
                    let (l, a) = (mm(length)?, number(angle, Unit::Radian, parameters)?);
                    let r = rise.as_ref().map_or(Ok(0.), mm)?;
                    if l < 0. {
                        return Err(invalid(
                            "A polar node needs a length of 0 mm or more; turn the angle instead",
                        ));
                    }
                    let (s, c) = a.portable_sin_cos();
                    add(
                        base,
                        add(add(scale(u, l * c), scale(v, l * s)), scale(n, r)),
                    )
                }
                NodePosition::Along {
                    member,
                    distance,
                    from_end,
                } => {
                    let m = frame.member(*member).ok_or_else(|| {
                        invalid("A station node's member was deleted; delete the station first")
                    })?;
                    let a = place(frame, m.start, axes, parameters, done, visiting)?;
                    let b = place(frame, m.end, axes, parameters, done, visiting)?;
                    let length = distance3(a, b);
                    let d = mm(distance)?;
                    if !(0. ..=length).contains(&d) || length < PATH_TOLERANCE_MM {
                        return Err(invalid(format!(
                            "A station {d} mm along a {length:.3} mm member lies outside it"
                        )));
                    }
                    let t = if *from_end { length - d } else { d } / length;
                    add(a, scale(sub(b, a), t))
                }
            };
            visiting.remove(&id);
            done.insert(id, p);
            Ok(p)
        }
        for node in &self.nodes {
            place(
                self,
                node.id,
                (origin, u, v, n),
                parameters,
                &mut done,
                &mut visiting,
            )?;
        }
        // Members join only at shared nodes: two separate nodes at one
        // place would leave members overlapping without a joint.
        let placed: Vec<_> = done.iter().collect();
        for (i, (_, a)) in placed.iter().enumerate() {
            if placed[i + 1..]
                .iter()
                .any(|(_, b)| distance3(**a, **b) < PATH_TOLERANCE_MM)
            {
                return Err(invalid(
                    "Two frame nodes are at the same place, so their members do not meet at a joint; join them at one node",
                ));
            }
        }
        for m in &self.members {
            if distance3(done[&m.start], done[&m.end]) < PATH_TOLERANCE_MM {
                return Err(invalid(
                    "A frame member has zero length; its two nodes coincide",
                ));
            }
        }
        // The in-plane centroid of all nodes: the interior for side rules.
        let count = done.len().max(1) as f64;
        let sum = done.values().fold([0.; 3], |s, p| add(s, *p));
        let mean = scale(sum, 1. / count);
        let centroid = sub(mean, scale(n, dot(sub(mean, origin), n)));
        Ok(EvaluatedFrame {
            origin,
            u,
            v,
            n,
            nodes: done,
            centroid,
        })
    }

    // --- Typed edits shared by the UI and MCP ------------------------------

    pub fn add_node(&mut self, position: NodePosition) -> EntityId {
        let id = EntityId::new();
        self.nodes.push(FrameNode { id, position });
        id
    }
    pub fn add_member(&mut self, start: EntityId, end: EntityId) -> EntityId {
        let id = EntityId::new();
        self.members.push(FrameMember {
            id,
            start,
            end,
            side: None,
        });
        id
    }
    /// Draw a chain of members from node `from` (§4.1): one polar node per
    /// segment, then a member to it. A segment that continues the previous
    /// member straight on extends that member instead (Q5), so a collinear
    /// chain stays one member. `close` joins the last node back to `from`.
    /// Returns the members created, in order.
    pub fn add_chain(
        &mut self,
        parameters: &EvaluatedParameters,
        from: EntityId,
        segments: &[ChainSegment],
        close: bool,
    ) -> Result<Vec<EntityId>, DocumentError> {
        if self.node(from).is_none() {
            return Err(invalid("A chain starts at a node of its frame"));
        }
        if segments.is_empty() {
            return Err(invalid("A chain needs at least one segment"));
        }
        let mut created = vec![];
        let mut at = from;
        let mut previous: Option<EntityId> = None;
        for segment in segments {
            let node = self.add_node(NodePosition::Polar {
                from: at,
                length: segment.length.clone(),
                angle: segment.angle.clone(),
                rise: segment.rise.clone(),
            });
            let continued = previous.filter(|m| self.continues_straight(parameters, *m, node));
            match continued {
                Some(member) => {
                    let m = self.members.iter_mut().find(|x| x.id == member).unwrap();
                    m.end = node;
                }
                None => {
                    let member = self.add_member(at, node);
                    created.push(member);
                    previous = Some(member);
                }
            }
            at = node;
        }
        if close {
            if at == from {
                return Err(invalid("The chain already ends where it starts"));
            }
            created.push(self.add_member(at, from));
        }
        self.evaluate(parameters)?;
        Ok(created)
    }
    /// Whether ending `member` at `node` instead keeps it straight: the new
    /// node lies on the member's line beyond its end, and nothing else
    /// meets the member's current end.
    fn continues_straight(
        &self,
        parameters: &EvaluatedParameters,
        member: EntityId,
        node: EntityId,
    ) -> bool {
        let Ok(frame) = self.evaluate(parameters).map(|f| f.nodes) else {
            return false;
        };
        let Some(m) = self.member(member) else {
            return false;
        };
        let (a, b, c) = (frame[&m.start], frame[&m.end], frame[&node]);
        let (Some(d1), Some(d2)) = (unit(sub(b, a)), unit(sub(c, b))) else {
            return false;
        };
        let others = self
            .members
            .iter()
            .any(|x| x.id != member && (x.start == m.end || x.end == m.end));
        let station = self
            .nodes
            .iter()
            .any(|n| matches!(n.position, NodePosition::Along { .. }) && n.id == m.end);
        !others && !station && norm(cross(d1, d2)) < PARALLEL && dot(d1, d2) > 0.
    }
    /// A station on `member`, to start a T branch or cross an X.
    pub fn add_station(
        &mut self,
        member: EntityId,
        distance: Dimension,
        from_end: bool,
    ) -> Result<EntityId, DocumentError> {
        if self.member(member).is_none() {
            return Err(invalid("A station lies on a member of its frame"));
        }
        Ok(self.add_node(NodePosition::Along {
            member,
            distance,
            from_end,
        }))
    }
    /// Rewrite one node's own definition; nodes that depend on it follow.
    pub fn set_node(
        &mut self,
        node: EntityId,
        position: NodePosition,
    ) -> Result<(), DocumentError> {
        let n = self
            .nodes
            .iter_mut()
            .find(|n| n.id == node)
            .ok_or_else(|| invalid("Choose a node of this frame"))?;
        n.position = position;
        Ok(())
    }
    /// Store, change or (when automatic) clear a node's joint override.
    pub fn set_joint(&mut self, joint: JointOverride) {
        self.joints.retain(|j| j.node != joint.node);
        if !joint.is_automatic() {
            self.joints.push(joint);
        }
    }
    /// Remove a member. Joint overrides naming it stay and are reported as
    /// stale; a station on it must be deleted first.
    pub fn remove_member(&mut self, member: EntityId) -> Result<(), DocumentError> {
        if self.member(member).is_none() {
            return Err(invalid("Choose a member of this frame"));
        }
        if self
            .nodes
            .iter()
            .any(|n| matches!(n.position, NodePosition::Along { member: m, .. } if m == member))
        {
            return Err(invalid(
                "A station node lies on this member; delete the branch at the station first",
            ));
        }
        self.members.retain(|m| m.id != member);
        Ok(())
    }
}
fn distance3(a: [f64; 3], b: [f64; 3]) -> f64 {
    distance(a, b)
}

/// Evaluated node positions of a frame, in world millimeters.
#[derive(Clone, Debug, PartialEq)]
pub struct EvaluatedFrame {
    pub origin: [f64; 3],
    pub u: [f64; 3],
    pub v: [f64; 3],
    pub n: [f64; 3],
    pub nodes: BTreeMap<EntityId, [f64; 3]>,
    pub centroid: [f64; 3],
}

pub fn frames(document: &Document) -> impl Iterator<Item = &Entity> {
    document.entities().filter(
        |e| matches!(&e.data, EntityData::Extension { definition, .. } if definition.type_id == FRAME_KIND),
    )
}
pub fn read_frame(document: &Document, id: EntityId) -> Result<FrameDefinition, DocumentError> {
    let Some(EntityData::Extension { definition, inputs }) = document.get(id).map(|e| &e.data)
    else {
        return Err(match document.get(id) {
            Some(_) => invalid(format!("{id} is not a {FRAME_KIND}")),
            None => DocumentError::Missing(id),
        });
    };
    let value: FrameDefinition = parse(definition, FRAME_KIND)?;
    value.validate()?;
    let expected = value.inputs();
    if inputs.len() != expected.len() || expected.iter().any(|l| !inputs.contains(l)) {
        return Err(invalid("Frame source and canonical input links disagree"));
    }
    Ok(value)
}
/// The member feature carrying frame member `member`, if any.
pub fn frame_member_feature(
    document: &Document,
    frame: EntityId,
    member: EntityId,
) -> Option<EntityId> {
    members(document).map(|e| e.id).find(|id| {
        read_member(document, *id)
            .is_ok_and(|m| m.path == MemberPath::FrameMember { frame, member })
    })
}
/// Every member feature of a frame with its frame member id.
pub fn frame_member_features(document: &Document, frame: EntityId) -> Vec<(EntityId, EntityId)> {
    members(document)
        .filter_map(|e| match read_member(document, e.id).ok()?.path {
            MemberPath::FrameMember { frame: f, member } if f == frame => Some((member, e.id)),
            _ => None,
        })
        .collect()
}

/// The world segment of a frame member. In-plane and inclined members take
/// the plane normal as their section height reference; members along the
/// normal (legs) take the plane `u` direction.
pub fn frame_member_path(
    document: &Document,
    parameters: &EvaluatedParameters,
    frame: EntityId,
    member: EntityId,
) -> Result<PathSegment, DocumentError> {
    let definition = read_frame(document, frame)?;
    let evaluated = definition.evaluate(parameters)?;
    member_segment(&definition, &evaluated, member)
}
fn member_segment(
    definition: &FrameDefinition,
    frame: &EvaluatedFrame,
    member: EntityId,
) -> Result<PathSegment, DocumentError> {
    let m = definition
        .member(member)
        .ok_or_else(|| invalid("This member was removed from its frame"))?;
    let (start, end) = (frame.nodes[&m.start], frame.nodes[&m.end]);
    let axis = unit(sub(end, start)).ok_or_else(|| invalid("Zero-length frame member"))?;
    let reference = if norm(cross(axis, frame.n)) < 1e-6 {
        frame.u
    } else {
        frame.n
    };
    Ok(PathSegment {
        start,
        end,
        reference: Some(reference),
    })
}
/// The section anchor of a frame member with a side rule: for each section
/// axis that leans toward the chosen side, the section's face on the line
/// is the one away from that side, so the section lies on it. Axes square
/// to the side keep `anchor_mm`.
pub(super) fn frame_anchor(
    document: &Document,
    parameters: &EvaluatedParameters,
    frame_id: EntityId,
    member: EntityId,
    member_frame: &MemberFrame,
    section: &ProfileSection,
    anchor_mm: [f64; 2],
) -> Result<[f64; 2], DocumentError> {
    let definition = read_frame(document, frame_id)?;
    let m = definition
        .member(member)
        .ok_or_else(|| invalid("This member was removed from its frame"))?;
    let side = m.side.unwrap_or(definition.side);
    if side == FrameSide::Center {
        return Ok(anchor_mm);
    }
    let frame = definition.evaluate(parameters)?;
    let (a, b) = (frame.nodes[&m.start], frame.nodes[&m.end]);
    let axis = member_frame.axis;
    let in_plane = |x: [f64; 3]| sub(x, scale(frame.n, dot(x, frame.n)));
    let toward = match side {
        FrameSide::Inner | FrameSide::Outer => {
            let mid = scale(add(a, b), 0.5);
            let r = in_plane(sub(frame.centroid, mid));
            let r = sub(r, scale(axis, dot(r, axis)));
            let r = unit(r);
            if side == FrameSide::Outer {
                r.map(|r| scale(r, -1.))
            } else {
                r
            }
        }
        FrameSide::Left => unit(cross(frame.n, axis)),
        FrameSide::Right => unit(cross(axis, frame.n)),
        FrameSide::Center => None,
    };
    let Some(d) = toward else {
        return Ok(anchor_mm);
    };
    let size = [section.width_mm, section.height_mm];
    let mut out = anchor_mm;
    for (k, axis) in [member_frame.x, member_frame.y].into_iter().enumerate() {
        let c = dot(axis, d);
        if c > 1e-6 {
            out[k] = 0.;
        } else if c < -1e-6 {
            out[k] = size[k];
        }
    }
    Ok(out)
}

/// One derived joint between two frame members at a node, in the order
/// `resolve_connection` takes them: first the mitre partner, terminating
/// member or branch; second the other mitre partner, the through member or
/// the receiving member.
#[derive(Clone, Debug, PartialEq)]
pub struct FramePair {
    pub node: EntityId,
    pub kind: ConnectionKind,
    pub members: [EntityId; 2],
    pub gap_mm: f64,
}
fn bigger(a: &ResolvedMember, b: &ResolvedMember) -> bool {
    let size = |m: &ResolvedMember| {
        (
            m.section.width_mm.max(m.section.height_mm),
            m.section.area_mm2(),
        )
    };
    size(a) > size(b)
}
fn is_round(m: &ResolvedMember) -> bool {
    m.section.shape == crate::materials::ProfileShape::Round
}
/// The derived joints at every node of a frame that touches `focus`
/// (every node when `focus` is none), with the frame's overrides applied.
/// A node whose joint cannot be made is an error naming the reason; a
/// stale override is an error too, never silently dropped.
pub fn frame_pairs(
    document: &Document,
    definition: &FrameDefinition,
    parameters: &EvaluatedParameters,
    members: &BTreeMap<EntityId, ResolvedMember>,
    focus: Option<EntityId>,
) -> Result<Vec<FramePair>, DocumentError> {
    let name = |member: EntityId| {
        frame_member_part_name(document, member).unwrap_or_else(|| "a member".into())
    };
    let gap = |joint: Option<&JointOverride>| -> Result<f64, DocumentError> {
        let d = joint
            .and_then(|j| j.gap.as_ref())
            .or(definition.gap.as_ref());
        let g = d.map_or(Ok(0.), |d| number(d, Unit::Millimeter, parameters))?;
        if g < 0. {
            return Err(invalid("A weld gap must be 0 mm or more"));
        }
        Ok(g)
    };
    let mut out = vec![];
    for joint in &definition.joints {
        if definition.node(joint.node).is_none() {
            return Err(invalid(
                "A joint override names a node that was deleted; remove the override",
            ));
        }
    }
    for node in &definition.nodes {
        let ends: Vec<EntityId> = definition
            .members
            .iter()
            .filter(|m| m.start == node.id || m.end == node.id)
            .map(|m| m.id)
            .collect();
        let station = match node.position {
            NodePosition::Along { member, .. } => Some(member),
            _ => None,
        };
        let incident: Vec<EntityId> = ends.iter().copied().chain(station).collect();
        if focus.is_some_and(|f| !incident.contains(&f)) {
            continue;
        }
        let joint = definition.joint(node.id);
        if incident.len() < 2 {
            if joint.is_some_and(|j| j.through.is_some() || j.kind.is_some()) {
                return Err(invalid(
                    "A joint override is stored at a node where fewer than two members meet now; remove the override",
                ));
            }
            continue;
        }
        let resolved = |id: EntityId| {
            members
                .get(&id)
                .ok_or_else(|| invalid("A frame member has no member feature"))
        };
        let gap_mm = gap(joint)?;
        // The through member: stored, a station's member, or derived.
        let through = match joint.and_then(|j| j.through) {
            Some(t) if incident.contains(&t) => Some(t),
            Some(t) => {
                return Err(invalid(format!(
                    "The joint override at this node makes {} run through, but it no longer meets the node; change or remove the override",
                    if definition.member(t).is_some() {
                        name(t)
                    } else {
                        "a deleted member".into()
                    }
                )));
            }
            None => station,
        };
        let kind_override = joint.and_then(|j| j.kind);
        if kind_override == Some(ConnectionKind::Miter) && through.is_some() {
            return Err(invalid(
                "A mitre has no through member; choose butt or cope, or clear the through member",
            ));
        }
        let through = match through {
            Some(t) => Some(t),
            None if ends.len() >= 3 => {
                // Rule 8: the member along the plane normal (a leg) runs
                // through; otherwise the largest section.
                let frame_n = definition.plane.normal;
                let leg = ends.iter().copied().find(|id| {
                    resolved(*id).is_ok_and(|m| norm(cross(m.frame.axis, frame_n)) < 1e-6)
                });
                let largest = ends.iter().copied().try_fold(
                    None::<EntityId>,
                    |best, id| -> Result<_, DocumentError> {
                        Ok(match best {
                            Some(b) if !bigger(resolved(id)?, resolved(b)?) => Some(b),
                            _ => Some(id),
                        })
                    },
                )?;
                leg.or(largest)
            }
            None => None,
        };
        match through {
            Some(t) => {
                for branch in ends.iter().copied().filter(|b| *b != t) {
                    let (bm, tm) = (resolved(branch)?, resolved(t)?);
                    let kind = kind_override.unwrap_or(
                        if is_round(bm)
                            || tm.section.shape != crate::materials::ProfileShape::Rectangular
                        {
                            ConnectionKind::Cope
                        } else {
                            ConnectionKind::Butt
                        },
                    );
                    out.push(FramePair {
                        node: node.id,
                        kind,
                        members: [branch, t],
                        gap_mm,
                    });
                }
            }
            None => {
                let [a, b] = [ends[0], ends[1]];
                let (am, bm) = (resolved(a)?, resolved(b)?);
                // Rule 2: a straight run with nothing else needs no joint.
                if norm(cross(am.frame.axis, bm.frame.axis)) < PARALLEL && kind_override.is_none() {
                    continue;
                }
                let same = am.stock == bm.stock && am.section == bm.section;
                let (first, second) = if bigger(am, bm) { (b, a) } else { (a, b) };
                let kind = match kind_override {
                    Some(kind) => kind,
                    // Rule 4: round tubes cope, larger through; equal
                    // rounds at a corner mitre.
                    None if (is_round(am) || is_round(bm)) && !(same && is_round(am)) => {
                        ConnectionKind::Cope
                    }
                    // Rule 5: mirroring sections mitre (Q4).
                    None if same => ConnectionKind::Miter,
                    // Rule 6: different sections butt, larger through.
                    None => ConnectionKind::Butt,
                };
                let members = match kind {
                    ConnectionKind::Miter => [a, b],
                    // A round branch copes onto the other member.
                    ConnectionKind::Cope if is_round(am) && !is_round(bm) => [a, b],
                    ConnectionKind::Cope if is_round(bm) && !is_round(am) => [b, a],
                    _ => [first, second],
                };
                out.push(FramePair {
                    node: node.id,
                    kind,
                    members,
                    gap_mm,
                });
            }
        }
    }
    Ok(out)
}
/// The Part name of a frame member, for messages.
fn frame_member_part_name(document: &Document, member: EntityId) -> Option<String> {
    let feature = members(document).find(|e| {
        read_member(document, e.id).is_ok_and(
            |m| matches!(m.path, MemberPath::FrameMember { member: x, .. } if x == member),
        )
    })?;
    let body = read_member(document, feature.id).ok()?.output.body;
    document.get(body).map(|p| p.name.clone())
}

pub(super) fn validate_frames(document: &Document) -> Result<(), DocumentError> {
    let mut carried = BTreeSet::new();
    for e in members(document) {
        let member = read_member(document, e.id)?;
        if let MemberPath::FrameMember { frame, member: m } = member.path {
            let definition = read_frame(document, frame)?;
            if definition.member(m).is_none() {
                return Err(invalid(
                    "A frame member feature names a member its frame no longer has",
                ));
            }
            if !carried.insert((frame, m)) {
                return Err(invalid(
                    "A frame member carries more than one member feature",
                ));
            }
        }
    }
    for e in frames(document) {
        let definition = read_frame(document, e.id)?;
        for m in &definition.members {
            if !carried.contains(&(e.id, m.id)) {
                return Err(invalid(format!(
                    "Frame {} has a member without its member feature and Part",
                    e.name
                )));
            }
        }
    }
    Ok(())
}

/// One node's derived joint as the Frame tool shows it: the joint kind and
/// whether its members are reversed from automatic, or the reason it cannot
/// be built.
#[derive(Clone, Debug, PartialEq)]
pub struct NodeJoint {
    pub node: EntityId,
    pub pairs: Vec<FramePair>,
    pub overridden: bool,
    pub problem: Option<String>,
}
/// The derived joint at every node where two or more members meet, checked
/// against the same cut builders the evaluator uses, on a document that
/// already holds the frame and its member features.
pub fn frame_joint_report(
    document: &Document,
    parameters: &EvaluatedParameters,
    frame: EntityId,
) -> Result<Vec<NodeJoint>, DocumentError> {
    let definition = read_frame(document, frame)?;
    let mut resolved = BTreeMap::new();
    for (member, feature) in frame_member_features(document, frame) {
        let source = read_member(document, feature)?;
        let segment = frame_member_path(document, parameters, frame, member)?;
        resolved.insert(
            member,
            resolve_member(document, parameters, &source, segment)?,
        );
    }
    let mut out = vec![];
    for node in &definition.nodes {
        let mut single = definition.clone();
        // Report each node on its own so one refusal does not hide others.
        single.joints.retain(|j| j.node == node.id);
        let incident = definition
            .members
            .iter()
            .filter(|m| m.start == node.id || m.end == node.id)
            .count()
            + usize::from(matches!(node.position, NodePosition::Along { .. }));
        if incident < 2 && definition.joint(node.id).is_none() {
            continue;
        }
        single.nodes.retain(|n| n.id == node.id);
        let pairs = node_pairs(
            document,
            &definition,
            &single,
            parameters,
            &resolved,
            node.id,
        );
        let overridden = definition.joint(node.id).is_some();
        match pairs {
            Ok(pairs) => {
                let problem = pairs.iter().find_map(|p| {
                    resolve_connection(
                        p.kind,
                        &resolved[&p.members[0]],
                        &resolved[&p.members[1]],
                        p.gap_mm,
                    )
                    .err()
                    .map(|e| match e {
                        DocumentError::Invalid(m) => m,
                        other => other.to_string(),
                    })
                });
                if pairs.is_empty() && problem.is_none() && !overridden {
                    continue;
                }
                out.push(NodeJoint {
                    node: node.id,
                    pairs,
                    overridden,
                    problem,
                });
            }
            Err(e) => out.push(NodeJoint {
                node: node.id,
                pairs: vec![],
                overridden,
                problem: Some(match e {
                    DocumentError::Invalid(m) => m,
                    other => other.to_string(),
                }),
            }),
        }
    }
    Ok(out)
}
/// The pairs at one node: `single` holds only that node and its override,
/// with the full definition's members.
fn node_pairs(
    document: &Document,
    full: &FrameDefinition,
    single: &FrameDefinition,
    parameters: &EvaluatedParameters,
    resolved: &BTreeMap<EntityId, ResolvedMember>,
    node: EntityId,
) -> Result<Vec<FramePair>, DocumentError> {
    let mut probe = full.clone();
    probe.joints = single.joints.clone();
    Ok(frame_pairs(document, &probe, parameters, resolved, None)?
        .into_iter()
        .filter(|p| p.node == node)
        .collect())
}
