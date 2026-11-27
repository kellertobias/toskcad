//! Frame tools (CAD-439): create and edit `weldment.frame` sources through
//! the same typed frame edits and `prepare_frame` transaction as the UI.
//! Every member is built natively, with its derived node joints, before the
//! one transaction commits; the cut list is the weldment cut list.
use super::{Quantity, Request, ToolResult, commit, compile, id, revision_guard, snapshot};
use crate::{
    error::ToolError,
    host::{Host, with},
};
use multicad_core::{
    constraint_sketch::Dimension,
    document::{Document, EntityId},
    model::{
        self, ChainSegment, ConnectionKind, FrameDefinition, FrameMemberRequest, FramePlane,
        FrameRequest, FrameSide, JointOverride, NodePosition, angle, length,
    },
    parameters::EvaluatedParameters,
};
use schemars::JsonSchema;
use serde::Deserialize;
use serde_json::{Value, json};

fn q(q: &Quantity, d: &Document, e: &EvaluatedParameters) -> Result<Dimension, ToolError> {
    Ok(q.dimension(d, e)?.0)
}
fn opt(
    v: &Option<Quantity>,
    d: &Document,
    e: &EvaluatedParameters,
) -> Result<Option<Dimension>, ToolError> {
    v.as_ref().map(|v| q(v, d, e)).transpose()
}

#[derive(Clone, Copy, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum FrameSideArg {
    Center,
    Inner,
    Outer,
    Left,
    Right,
}
impl From<FrameSideArg> for FrameSide {
    fn from(s: FrameSideArg) -> Self {
        match s {
            FrameSideArg::Center => Self::Center,
            FrameSideArg::Inner => Self::Inner,
            FrameSideArg::Outer => Self::Outer,
            FrameSideArg::Left => Self::Left,
            FrameSideArg::Right => Self::Right,
        }
    }
}
#[derive(Clone, Copy, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum JointKindArg {
    Miter,
    Butt,
    Cope,
}
impl From<JointKindArg> for ConnectionKind {
    fn from(k: JointKindArg) -> Self {
        match k {
            JointKindArg::Miter => Self::Miter,
            JointKindArg::Butt => Self::Butt,
            JointKindArg::Cope => Self::Cope,
        }
    }
}
/// One drawn segment: in-plane length along `angle_deg` (from the plane u
/// axis toward v), then `rise` along the plane normal.
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct FrameSegmentArg {
    pub length: Quantity,
    pub angle_deg: f64,
    #[serde(default)]
    pub rise: Option<Quantity>,
}
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct FramePlaneArg {
    pub origin_mm: [f64; 3],
    /// Unit in-plane u direction.
    pub u: [f64; 3],
    /// Unit plane normal, square to u.
    pub normal: [f64; 3],
}
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct CreateFrameArgs {
    pub expected_revision: u64,
    pub name: String,
    /// Project profile stock for every member drawn.
    pub stock: String,
    /// Base plane; the world XY plane when omitted.
    #[serde(default)]
    pub plane: Option<FramePlaneArg>,
    /// Default side of the section on each line: center (default), inner or
    /// outer of the frame, or left or right looking along each member.
    #[serde(default)]
    pub side: Option<FrameSideArg>,
    /// Frame-wide weld gap (default 0 mm). A gap shortens the members.
    #[serde(default)]
    pub gap: Option<Quantity>,
    /// The first node, in plane coordinates.
    pub start: StartArg,
    /// The chain drawn from the first node.
    pub segments: Vec<FrameSegmentArg>,
    /// Join the last node back to the first.
    #[serde(default)]
    pub close: bool,
    /// Member names are "<prefix> N" (default "Member").
    #[serde(default)]
    pub name_prefix: Option<String>,
}
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct StartArg {
    pub u: Quantity,
    pub v: Quantity,
    #[serde(default)]
    pub w: Option<Quantity>,
}
#[allow(clippy::large_enum_variant)]
#[derive(Deserialize, JsonSchema)]
#[serde(tag = "op", rename_all = "snake_case", deny_unknown_fields)]
pub enum FrameOperationArg {
    /// Draw a chain from an existing node; a straight continuation extends
    /// the previous member.
    AddChain {
        from_node: String,
        segments: Vec<FrameSegmentArg>,
        #[serde(default)]
        close: bool,
    },
    /// A station on a member, to start a T branch or cross an X.
    AddStation {
        member: String,
        distance: Quantity,
        #[serde(default)]
        from_end: bool,
    },
    /// Rewrite one node's own position: absolute plane coordinates, or polar
    /// from another node. Nodes that depend on it follow.
    SetNode {
        node: String,
        #[serde(default)]
        from_node: Option<String>,
        #[serde(default)]
        u: Option<Quantity>,
        #[serde(default)]
        v: Option<Quantity>,
        #[serde(default)]
        w: Option<Quantity>,
        #[serde(default)]
        length: Option<Quantity>,
        #[serde(default)]
        angle_deg: Option<f64>,
        #[serde(default)]
        rise: Option<Quantity>,
    },
    /// Override a node's joint; omit every field to make it automatic again.
    SetJoint {
        node: String,
        #[serde(default)]
        kind: Option<JointKindArg>,
        #[serde(default)]
        through: Option<String>,
        #[serde(default)]
        gap: Option<Quantity>,
    },
    /// Change one member's side, or the frame default without `member`.
    SetSide {
        #[serde(default)]
        member: Option<String>,
        side: FrameSideArg,
    },
    /// Delete a member with its Part. Overrides naming it are reported.
    RemoveMember { member: String },
}
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct EditFrameArgs {
    pub expected_revision: u64,
    pub frame: String,
    /// Stock for members the operations create; required when they do.
    #[serde(default)]
    pub stock: Option<String>,
    pub operations: Vec<FrameOperationArg>,
    #[serde(default)]
    pub name_prefix: Option<String>,
}

fn segments(
    args: &[FrameSegmentArg],
    d: &Document,
    e: &EvaluatedParameters,
) -> Result<Vec<ChainSegment>, ToolError> {
    args.iter()
        .map(|s| {
            Ok(ChainSegment {
                length: q(&s.length, d, e)?,
                angle: angle(s.angle_deg),
                rise: opt(&s.rise, d, e)?,
            })
        })
        .collect()
}
fn frame_json(frame: EntityId, definition: &FrameDefinition, document: &Document) -> Value {
    let carried: std::collections::BTreeMap<_, _> = model::frame_member_features(document, frame)
        .into_iter()
        .collect();
    json!({
        "frame": frame,
        "nodes": definition.nodes.iter().map(|n| n.id).collect::<Vec<_>>(),
        "frame_members": definition.members.iter().map(|m| {
            let feature = carried.get(&m.id);
            let part = feature.and_then(|f| model::read_member(document, *f).ok()).map(|x| x.output.body);
            json!({ "member": m.id, "start": m.start, "end": m.end, "feature": feature, "part": part })
        }).collect::<Vec<_>>(),
    })
}
/// Prepare, build natively with every derived joint, then commit.
fn apply(
    host: &dyn Host,
    call: &Request,
    revision: u64,
    request: FrameRequest,
    created: Vec<EntityId>,
) -> ToolResult {
    let call_owned = call.clone();
    let request_owned = request.clone();
    let prepared = with(host, move |target| {
        revision_guard(target, revision)?;
        Ok::<_, ToolError>(
            model::prepare_frame(
                target.session(),
                revision,
                call_owned.actor.clone(),
                &request_owned,
            )?
            .prepared,
        )
    })??;
    let model = compile(prepared.document(), revision, call)?;
    let records = model.member_cut_list(prepared.document())?;
    let mut out = frame_json(request.frame, &request.definition, prepared.document());
    out["created_members"] = json!(created);
    let call = call.clone();
    let committed = with(host, move |target| commit(target, &call, prepared))??;
    out["revision"] = json!(committed);
    out["members"] = json!(records.len());
    Ok(out)
}
fn member_requests(
    created: &[EntityId],
    stock: Option<EntityId>,
    prefix: &str,
    first: usize,
) -> Result<Vec<FrameMemberRequest>, ToolError> {
    if created.is_empty() {
        return Ok(vec![]);
    }
    let stock = stock.ok_or_else(|| ToolError::invalid("Give stock for the new members"))?;
    Ok(created
        .iter()
        .enumerate()
        .map(|(i, m)| FrameMemberRequest::new(*m, format!("{prefix} {}", first + i + 1), stock))
        .collect())
}

pub(super) fn create_frame(host: &dyn Host, call: &Request, args: CreateFrameArgs) -> ToolResult {
    let (document, revision, evaluated) = snapshot(host)?;
    if args.expected_revision != revision {
        return Err(ToolError::new(
            "stale_revision",
            format!(
                "Request is based on revision {}, but the document is at {revision}",
                args.expected_revision
            ),
        )
        .at_revision(revision));
    }
    let (d, e) = (&document, &evaluated);
    let plane = match &args.plane {
        Some(p) => FramePlane {
            origin: p.origin_mm.map(length),
            u: p.u,
            normal: p.normal,
        },
        None => FramePlane::xy(),
    };
    let mut definition = FrameDefinition::new(plane);
    definition.side = args.side.map_or(FrameSide::Center, Into::into);
    definition.gap = opt(&args.gap, d, e)?;
    let start = definition.add_node(NodePosition::Absolute {
        u: q(&args.start.u, d, e)?,
        v: q(&args.start.v, d, e)?,
        w: opt(&args.start.w, d, e)?,
    });
    let created = definition.add_chain(e, start, &segments(&args.segments, d, e)?, args.close)?;
    let prefix = args.name_prefix.as_deref().unwrap_or("Member");
    let members = member_requests(&created, Some(id(&args.stock)?), prefix, 0)?;
    let request = FrameRequest {
        frame: EntityId::new(),
        name: args.name.clone(),
        definition,
        members,
    };
    apply(host, call, revision, request, created)
}

pub(super) fn edit_frame(host: &dyn Host, call: &Request, args: EditFrameArgs) -> ToolResult {
    let (document, revision, evaluated) = snapshot(host)?;
    if args.expected_revision != revision {
        return Err(ToolError::new(
            "stale_revision",
            format!(
                "Request is based on revision {}, but the document is at {revision}",
                args.expected_revision
            ),
        )
        .at_revision(revision));
    }
    let (d, e) = (&document, &evaluated);
    let frame = id(&args.frame)?;
    let mut definition = model::read_frame(d, frame)?;
    let name = d.get(frame).map(|x| x.name.clone()).unwrap_or_default();
    let mut created = vec![];
    let mut stations = vec![];
    for op in &args.operations {
        match op {
            FrameOperationArg::AddChain {
                from_node,
                segments: s,
                close,
            } => {
                created.extend(definition.add_chain(
                    e,
                    id(from_node)?,
                    &segments(s, d, e)?,
                    *close,
                )?);
            }
            FrameOperationArg::AddStation {
                member,
                distance,
                from_end,
            } => {
                stations.push(definition.add_station(id(member)?, q(distance, d, e)?, *from_end)?)
            }
            FrameOperationArg::SetNode {
                node,
                from_node,
                u,
                v,
                w,
                length: l,
                angle_deg,
                rise,
            } => {
                let position = match (from_node, u, v) {
                    (Some(from), None, None) => NodePosition::Polar {
                        from: id(from)?,
                        length: l
                            .as_ref()
                            .map(|x| q(x, d, e))
                            .transpose()?
                            .unwrap_or_else(|| length(0.)),
                        angle: angle(angle_deg.unwrap_or(0.)),
                        rise: opt(rise, d, e)?,
                    },
                    (None, Some(u), Some(v)) => NodePosition::Absolute {
                        u: q(u, d, e)?,
                        v: q(v, d, e)?,
                        w: opt(w, d, e)?,
                    },
                    _ => {
                        return Err(ToolError::invalid(
                            "Give from_node (with length, angle_deg, rise) or u and v (with w)",
                        ));
                    }
                };
                definition.set_node(id(node)?, position)?;
            }
            FrameOperationArg::SetJoint {
                node,
                kind,
                through,
                gap,
            } => definition.set_joint(JointOverride {
                node: id(node)?,
                kind: kind.map(Into::into),
                through: through.as_deref().map(id).transpose()?,
                gap: opt(gap, d, e)?,
            }),
            FrameOperationArg::SetSide { member, side } => match member {
                Some(m) => {
                    let m = id(m)?;
                    let member = definition
                        .members
                        .iter_mut()
                        .find(|x| x.id == m)
                        .ok_or_else(|| ToolError::invalid("Choose a member of this frame"))?;
                    member.side = Some((*side).into());
                }
                None => definition.side = (*side).into(),
            },
            FrameOperationArg::RemoveMember { member } => definition.remove_member(id(member)?)?,
        }
    }
    let prefix = args.name_prefix.as_deref().unwrap_or("Member");
    let existing = model::frame_member_features(d, frame).len();
    let stock = args.stock.as_deref().map(id).transpose()?;
    let members = member_requests(&created, stock, prefix, existing)?;
    let request = FrameRequest {
        frame,
        name,
        definition,
        members,
    };
    let mut out = apply(host, call, revision, request, created)?;
    out["created_stations"] = json!(stations);
    Ok(out)
}
