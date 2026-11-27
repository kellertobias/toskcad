//! Model → Metal → Frame (CAD-440): draw a frame of profile stock as lines.
//!
//! The form holds a working copy of the `weldment.frame` source. Clicks on
//! the drawing plane add nodes and members (snapping to nodes and 15°
//! steps); typed segments add polar nodes whose length keeps any
//! `=$parameter`, at an angle relative to the previous member. Joints are
//! derived by the core at every node and drawn as badges; clicking a node
//! opens its joint editor. Nothing is written until Apply, which commits
//! one `prepare_frame` transaction; Discard keeps the source.
use super::{
    FormStamp, OutlineKey, PreviewState, TerminationViewport,
    dimension::DimensionField,
    form::{self, FormAction, ModelEvent, ModelRequest},
    value_entry::format_number,
};
use crate::{controls, egui};
use multicad_core::{
    commands::{Actor, Session},
    document::{Document, EntityId},
    model::{
        self, ChainSegment, ConnectionKind, FrameDefinition, FrameMemberRequest, FramePlane,
        FrameRequest, FrameSide, JointOverride, ModelPrepared, NodeJoint, NodePosition, angle,
        length,
    },
    parameters::{EvaluatedParameters, Unit},
    sketch::Frame,
};
use multicad_precision::PortableMath as _;
use std::collections::BTreeMap;

const PICK_PX: f32 = 9.;
const ANGLE_STEP_DEG: f64 = 15.;
const ANGLE_SNAP_DEG: f64 = 4.;

/// The typed frame intent sent to native preview and Apply.
#[derive(Clone, Debug)]
pub struct FrameCommand {
    pub request: FrameRequest,
    pub editing: bool,
}
impl FrameCommand {
    pub(super) fn prepare(&self, session: &Session, actor: Actor) -> Result<ModelPrepared, String> {
        model::prepare_frame(session, session.revision(), actor, &self.request)
            .map_err(|e| e.to_string())
    }
    /// Whether `feature` is a member feature of the edited frame.
    pub(super) fn edits(&self, document: &Document, feature: EntityId) -> bool {
        self.editing
            && model::frame_member_features(document, self.request.frame)
                .iter()
                .any(|(_, f)| *f == feature)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum DrawPlane {
    Top,
    Front,
    Side,
}
impl DrawPlane {
    const ALL: [Self; 3] = [Self::Top, Self::Front, Self::Side];
    fn label(self) -> &'static str {
        match self {
            Self::Top => "Top (XY)",
            Self::Front => "Front (XZ)",
            Self::Side => "Side (YZ)",
        }
    }
    fn plane(self) -> FramePlane {
        let mut plane = FramePlane::xy();
        (plane.u, plane.normal) = match self {
            Self::Top => ([1., 0., 0.], [0., 0., 1.]),
            Self::Front => ([1., 0., 0.], [0., -1., 0.]),
            Self::Side => ([0., 1., 0.], [1., 0., 0.]),
        };
        plane
    }
}
/// Planes a frame can be drawn on: the three datum planes, each level at
/// its elevation (kept as its own dimension, so a parameter-driven level
/// moves the frame), and every work or sketch plane of the document.
pub(crate) fn plane_choices(document: &Document) -> Vec<(String, FramePlane)> {
    let mut out: Vec<(String, FramePlane)> = DrawPlane::ALL
        .into_iter()
        .map(|p| (p.label().to_owned(), p.plane()))
        .collect();
    for level in model::levels(document) {
        if let Ok(definition) = model::read_level(document, level.id) {
            let mut plane = DrawPlane::Top.plane();
            plane.origin[2] = definition.elevation.clone();
            out.push((format!("Level {}", level.name), plane));
        }
    }
    for entity in document.entities() {
        if let multicad_core::document::EntityData::Plane {
            origin_mm, normal, ..
        } = &entity.data
        {
            let n = *normal;
            let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
            if len < 1e-12 {
                continue;
            }
            let n = n.map(|c| c / len);
            // u: world X square to the normal, or world Y for planes facing X.
            let seed = if n[0].abs() > 0.9 {
                [0., 1., 0.]
            } else {
                [1., 0., 0.]
            };
            let d = seed[0] * n[0] + seed[1] * n[1] + seed[2] * n[2];
            let u = [seed[0] - d * n[0], seed[1] - d * n[1], seed[2] - d * n[2]];
            let ul = (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt();
            out.push((
                format!("Plane {}", entity.name),
                FramePlane {
                    origin: origin_mm.map(length),
                    u: u.map(|c| c / ul),
                    normal: n,
                },
            ));
        }
    }
    out
}
const SIDES: [(FrameSide, &str); 6] = [
    (FrameSide::Inner, "Inner edge on the line"),
    (FrameSide::Outer, "Outer edge on the line"),
    (FrameSide::Center, "Centered on the line"),
    (FrameSide::Left, "Left of the line"),
    (FrameSide::Right, "Right of the line"),
    (FrameSide::Centroid, "Centroid on the line"),
];
fn side_label(side: FrameSide) -> &'static str {
    SIDES
        .iter()
        .find(|(s, _)| *s == side)
        .map_or("Centered on the line", |(_, l)| l)
}

/// What a joint report was computed for: the form state, the document
/// revision, and a node being moved with its rounded position.
type ReportKey = (u64, u64, Option<(EntityId, [i64; 2])>);

pub struct FramePanel {
    stamp: FormStamp,
    frame: EntityId,
    name: String,
    editing: bool,
    definition: FrameDefinition,
    stock: Option<EntityId>,
    /// Feature and Part ids of members new in this edit, created once so
    /// Preview and Apply commit the same identities.
    new_members: BTreeMap<EntityId, FrameMemberRequest>,
    /// Member names of applied members, for labels.
    names: BTreeMap<EntityId, String>,
    drawing: bool,
    /// The next viewport click picks a planar face as the drawing plane.
    picking_face: bool,
    face_label: Option<String>,
    /// The nodes of the chain being drawn, with the member that reached each.
    chain: Vec<(EntityId, Option<EntityId>)>,
    hover: Option<[f64; 2]>,
    moving: Option<EntityId>,
    selected_node: Option<EntityId>,
    selected_member: Option<EntityId>,
    start: [DimensionField; 2],
    segment_length: DimensionField,
    segment_angle: DimensionField,
    gap_text: String,
    report: Option<(ReportKey, Result<Vec<NodeJoint>, String>)>,
    /// Placed section outlines at both ends of every member, from the
    /// same candidate the badges come from.
    sections: Vec<[Vec<[f64; 3]>; 2]>,
    dragging: bool,
    drag_from: Option<[f64; 2]>,
    error: Option<String>,
    counter: usize,
}

impl FramePanel {
    pub fn new(document: &Document, name: String) -> Result<Self, String> {
        let mut definition = FrameDefinition::new(DrawPlane::Top.plane());
        definition.side = FrameSide::Inner;
        let stock = super::weldment::profile_stock_specs(document, false)
            .first()
            .map(|s| s.0);
        Self::with(name, EntityId::new(), definition, stock, false)
    }
    fn with(
        name: String,
        frame: EntityId,
        definition: FrameDefinition,
        stock: Option<EntityId>,
        editing: bool,
    ) -> Result<Self, String> {
        let field = |v: f64, unit: &'static str| {
            DimensionField::from_dimension(
                &multicad_core::constraint_sketch::Dimension::Literal {
                    value: v,
                    unit: unit.into(),
                },
                unit,
            )
        };
        Ok(Self {
            stamp: FormStamp::new(),
            frame,
            name,
            editing,
            definition,
            stock,
            new_members: BTreeMap::new(),
            names: BTreeMap::new(),
            drawing: false,
            chain: vec![],
            hover: None,
            moving: None,
            selected_node: None,
            selected_member: None,
            start: [field(0., "mm")?, field(0., "mm")?],
            segment_length: field(1000., "mm")?,
            segment_angle: field(0., "deg")?,
            gap_text: String::new(),
            report: None,
            sections: vec![],
            dragging: false,
            drag_from: None,
            error: None,
            counter: 0,
            picking_face: false,
            face_label: None,
        })
    }
    /// Reopen the frame that member `feature` belongs to.
    pub fn existing(document: &Document, feature: EntityId) -> Result<Self, String> {
        let member = model::read_member(document, feature).map_err(|e| e.to_string())?;
        let model::MemberPath::FrameMember { frame, .. } = member.path else {
            return Err("This member is not part of a frame".into());
        };
        let definition = model::read_frame(document, frame).map_err(|e| e.to_string())?;
        let name = document
            .get(frame)
            .map(|e| e.name.clone())
            .unwrap_or_default();
        let mut panel = Self::with(name, frame, definition, Some(member.stock), true)?;
        for (m, f) in model::frame_member_features(document, frame) {
            if let Ok(source) = model::read_member(document, f) {
                let part = document
                    .get(source.output.body)
                    .map(|e| e.name.clone())
                    .unwrap_or_default();
                panel.names.insert(m, part);
            }
        }
        panel.counter = panel.names.len();
        Ok(panel)
    }
    pub fn stamp(&self) -> FormStamp {
        self.stamp
    }
    pub fn label(&self) -> &'static str {
        if self.editing { "Edit frame" } else { "Frame" }
    }
    pub(super) fn select_stock(&mut self, stock: EntityId) {
        self.stock = Some(stock);
        self.changed();
    }
    fn changed(&mut self) {
        self.stamp.changed();
        self.error = None;
    }
    fn member_name(&self, member: EntityId) -> String {
        self.new_members
            .get(&member)
            .map(|r| r.name.clone())
            .or_else(|| self.names.get(&member).cloned())
            .unwrap_or_else(|| "Member".into())
    }

    pub fn request(&self, document: &Document, enabled: bool) -> Result<FrameCommand, String> {
        self.request_for(&self.definition, document, enabled)
    }
    fn request_for(
        &self,
        definition: &FrameDefinition,
        document: &Document,
        enabled: bool,
    ) -> Result<FrameCommand, String> {
        if !enabled {
            return Err("Native modeling is unavailable".into());
        }
        if self.name.trim().is_empty() {
            return Err("Name the frame".into());
        }
        if definition.members.is_empty() {
            return Err("Draw at least one member: click two points or type a segment".into());
        }
        for m in self.new_members.values() {
            model::read_profile_stock(document, m.stock).map_err(|e| e.to_string())?;
        }
        Ok(FrameCommand {
            request: FrameRequest {
                frame: self.frame,
                name: self.name.trim().into(),
                definition: definition.clone(),
                members: definition
                    .members
                    .iter()
                    .filter_map(|m| self.new_members.get(&m.id).cloned())
                    .collect(),
            },
            editing: self.editing,
        })
    }
    pub fn model_request(
        &self,
        document: &Document,
        enabled: bool,
    ) -> Result<ModelRequest, String> {
        Ok(ModelRequest::Frame(Box::new(
            self.request(document, enabled)?,
        )))
    }

    // --- Geometry ------------------------------------------------------------

    fn plane_frame(&self, parameters: &EvaluatedParameters) -> Frame {
        let p = &self.definition.plane;
        let origin = std::array::from_fn(|i| {
            p.origin[i]
                .value_in(Unit::Millimeter, Some(parameters))
                .unwrap_or(0.)
        });
        Frame {
            origin,
            u: p.u,
            v: p.v(),
            n: p.normal,
        }
    }
    fn positions(
        &self,
        definition: &FrameDefinition,
        parameters: &EvaluatedParameters,
    ) -> Option<BTreeMap<EntityId, [f64; 3]>> {
        definition.evaluate(parameters).ok().map(|f| f.nodes)
    }
    fn to_plane(frame: &Frame, p: [f64; 3]) -> [f64; 2] {
        let d = [
            p[0] - frame.origin[0],
            p[1] - frame.origin[1],
            p[2] - frame.origin[2],
        ];
        let dot = |a: [f64; 3]| d[0] * a[0] + d[1] * a[1] + d[2] * a[2];
        [dot(frame.u), dot(frame.v)]
    }
    /// The in-plane direction of the last chain member, in degrees.
    fn last_angle(&self, parameters: &EvaluatedParameters) -> Option<f64> {
        let (_, member) = self.chain.last()?;
        let m = self.definition.member((*member)?)?;
        let frame = self.plane_frame(parameters);
        let nodes = self.positions(&self.definition, parameters)?;
        let (a, b) = (
            Self::to_plane(&frame, nodes[&m.start]),
            Self::to_plane(&frame, nodes[&m.end]),
        );
        Some((b[1] - a[1]).portable_atan2(b[0] - a[0]).to_degrees())
    }
    fn snap(
        &self,
        raw: [f64; 2],
        nodes: &[(EntityId, [f64; 2])],
        last: Option<[f64; 2]>,
        directions: &[f64],
        tolerance: f64,
    ) -> ([f64; 2], Option<EntityId>) {
        if let Some((id, p)) = nodes
            .iter()
            .map(|(id, p)| (id, p, (p[0] - raw[0]).portable_hypot(p[1] - raw[1])))
            .filter(|(.., d)| *d <= tolerance)
            .min_by(|a, b| a.2.total_cmp(&b.2))
            .map(|(id, p, _)| (*id, *p))
        {
            return (p, Some(id));
        }
        if let Some(l) = last {
            let (dx, dy) = (raw[0] - l[0], raw[1] - l[1]);
            let r = dx.portable_hypot(dy);
            let a = dy.portable_atan2(dx).to_degrees();
            // 15° steps from the plane axes and from every member drawn
            // (which includes parallel and perpendicular to each), and the
            // exact directions of existing members.
            let wrap = |d: f64| (d + 540.).rem_euclid(360.) - 180.;
            let step = std::iter::once(0.)
                .chain(directions.iter().copied())
                .map(|base| base + (wrap(a - base) / ANGLE_STEP_DEG).round() * ANGLE_STEP_DEG)
                .min_by(|x, y| wrap(a - x).abs().total_cmp(&wrap(a - y).abs()))
                .unwrap_or(a);
            if r > 0. && wrap(a - step).abs() <= ANGLE_SNAP_DEG {
                let r = r.round();
                let (s, c) = step.to_radians().portable_sin_cos();
                return ([l[0] + r * c, l[1] + r * s], None);
            }
        }
        ([raw[0].round(), raw[1].round()], None)
    }

    // --- Edits ---------------------------------------------------------------

    fn register(&mut self, member: EntityId) -> Result<(), String> {
        let stock = self
            .stock
            .ok_or("Choose profile stock before drawing members")?;
        self.counter += 1;
        let name = format!("{} member {}", self.name.trim(), self.counter);
        self.new_members
            .insert(member, FrameMemberRequest::new(member, name, stock));
        Ok(())
    }
    fn absolute(&mut self, p: [f64; 2]) -> EntityId {
        self.definition.add_node(NodePosition::Absolute {
            u: length(p[0]),
            v: length(p[1]),
            w: None,
        })
    }
    /// A click on the plane while drawing: start, extend or close the chain.
    fn draw_click(&mut self, point: [f64; 2], hit: Option<EntityId>) -> Result<bool, String> {
        let node = match hit {
            Some(node) => node,
            None => self.absolute(point),
        };
        let Some(&(last, _)) = self.chain.last() else {
            self.chain.push((node, None));
            return Ok(true);
        };
        if last == node {
            return Ok(false);
        }
        let closes = self.chain.len() >= 3 && self.chain[0].0 == node;
        let member = self.definition.add_member(last, node);
        self.register(member)?;
        self.chain.push((node, Some(member)));
        if closes {
            self.drawing = false;
            self.chain.clear();
        }
        Ok(true)
    }
    /// Add the typed segment: its length (kept as typed, with any
    /// `=$parameter`) at the typed angle relative to the previous member.
    fn add_segment(
        &mut self,
        document: &Document,
        parameters: &EvaluatedParameters,
    ) -> Result<(), String> {
        let length_dim = self.segment_length.resolve(document)?;
        let relative = self
            .segment_angle
            .resolve(document)?
            .value_in(Unit::Degree, Some(parameters))
            .map_err(|e| e.to_string())?;
        if self.chain.is_empty() {
            let u = self.start[0].resolve(document)?;
            let v = self.start[1].resolve(document)?;
            let node = self
                .definition
                .add_node(NodePosition::Absolute { u, v, w: None });
            self.chain.push((node, None));
        }
        let absolute = self.last_angle(parameters).unwrap_or(0.) + relative;
        let from = self.chain.last().expect("chain start").0;
        let before: Vec<EntityId> = self.definition.members.iter().map(|m| m.id).collect();
        let created = self
            .definition
            .add_chain(
                parameters,
                from,
                &[ChainSegment {
                    length: length_dim,
                    angle: angle(absolute),
                    rise: None,
                }],
                false,
            )
            .map_err(|e| e.to_string())?;
        let node = self.definition.nodes.last().expect("new node").id;
        for m in &created {
            self.register(*m)?;
        }
        // A straight continuation extended the previous member instead.
        let member = created.first().copied().or_else(|| {
            self.definition
                .members
                .iter()
                .find(|m| before.contains(&m.id) && m.end == node)
                .map(|m| m.id)
        });
        self.chain.push((node, member));
        Ok(())
    }
    fn close_loop(&mut self) -> Result<(), String> {
        if self.chain.len() < 3 {
            return Err("Draw at least two members before closing the loop".into());
        }
        let (first, last) = (self.chain[0].0, self.chain.last().unwrap().0);
        let member = self.definition.add_member(last, first);
        self.register(member)?;
        self.chain.clear();
        self.drawing = false;
        Ok(())
    }
    /// Remove the last node of the chain with the member that reached it.
    fn remove_last(&mut self) -> bool {
        let Some((node, member)) = self.chain.pop() else {
            return false;
        };
        if let Some(member) = member {
            let _ = self.definition.remove_member(member);
            self.new_members.remove(&member);
        }
        let used = self
            .definition
            .members
            .iter()
            .any(|m| m.start == node || m.end == node)
            || self.definition.nodes.iter().any(|n| match &n.position {
                NodePosition::Polar { from, .. } => *from == node,
                _ => false,
            });
        if !used {
            self.definition.nodes.retain(|n| n.id != node);
        }
        true
    }
    /// Enter, Escape or Backspace while drawing; `None` when the key is the
    /// form's (Apply or Cancel).
    pub(super) fn outline_key(&mut self, key: OutlineKey) -> Option<bool> {
        if self.moving.take().is_some() {
            return Some(false);
        }
        if !self.drawing && self.chain.is_empty() {
            return None;
        }
        let changed = match key {
            OutlineKey::RemoveLast => self.remove_last(),
            OutlineKey::Finish | OutlineKey::Cancel => {
                // Enter ends the chain; Esc cancels the pending segment.
                self.drawing = false;
                self.chain.clear();
                self.hover = None;
                false
            }
        };
        if changed {
            self.changed();
        }
        Some(changed)
    }
    fn reverse_member(&mut self, member: EntityId) {
        if let Some(m) = self.definition.members.iter_mut().find(|m| m.id == member) {
            std::mem::swap(&mut m.start, &mut m.end);
            // Left and right look along the member: reversing keeps the
            // section on the same physical side.
            m.side = m.side.map(|s| match s {
                FrameSide::Left => FrameSide::Right,
                FrameSide::Right => FrameSide::Left,
                other => other,
            });
        }
    }

    // --- Joint badges ----------------------------------------------------------

    /// Recompute joint badges and section outlines when the source, the
    /// document or a moving node changed. A moving node is evaluated at its
    /// live position, so badges and sections follow a drag.
    fn refresh_report(&mut self, session: &Session, probe: Option<(EntityId, [f64; 2])>) {
        let key: ReportKey = (
            self.stamp.revision,
            session.revision(),
            probe.map(|(n, p)| (n, [p[0].round() as i64, p[1].round() as i64])),
        );
        if self.report.as_ref().is_some_and(|(k, _)| *k == key) {
            return;
        }
        let mut definition = self.definition.clone();
        if let Some((node, p)) = probe {
            let _ = definition.set_node(
                node,
                NodePosition::Absolute {
                    u: length(p[0]),
                    v: length(p[1]),
                    w: None,
                },
            );
        }
        let mut sections = vec![];
        let result = (|| -> Result<Vec<NodeJoint>, String> {
            let command = self.request_for(&definition, session.document(), true)?;
            let actor = Actor::local(multicad_core::commands::Channel::Ui);
            let prepared = command.prepare(session, actor)?;
            let document = prepared.prepared.document();
            let parameters = session.parameters();
            for (member, feature) in model::frame_member_features(document, self.frame) {
                let source = model::read_member(document, feature).map_err(|e| e.to_string())?;
                let segment = model::frame_member_path(document, parameters, self.frame, member)
                    .map_err(|e| e.to_string())?;
                let m = model::resolve_member(document, parameters, &source, segment)
                    .map_err(|e| e.to_string())?;
                let start = m.start_offset_mm();
                sections.push([
                    m.outline_world(start),
                    m.outline_world(start + m.length_mm()),
                ]);
            }
            model::frame_joint_report(document, parameters, self.frame).map_err(|e| e.to_string())
        })();
        self.sections = sections;
        self.report = Some((key, result));
    }
    fn badge(joint: &NodeJoint) -> (String, bool) {
        if joint.problem.is_some() {
            return ("!".into(), true);
        }
        let mut text = match joint.pairs.first().map(|p| p.kind) {
            Some(ConnectionKind::Miter) => "M".to_owned(),
            Some(ConnectionKind::Butt) => "B→".to_owned(),
            Some(ConnectionKind::Cope) => "C".to_owned(),
            None => "–".to_owned(),
        };
        if joint.overridden {
            text.push('*');
        }
        (text, false)
    }
    /// The joint shown for a node, if any.
    pub fn joint(&self, node: EntityId) -> Option<&NodeJoint> {
        match &self.report {
            Some((_, Ok(joints))) => joints.iter().find(|j| j.node == node),
            _ => None,
        }
    }
    /// Draw, pick and move on the plane; returns whether the source changed.
    pub fn viewport(&mut self, session: &Session, view: TerminationViewport<'_>) -> bool {
        if self.picking_face {
            if view.clicked && view.pointer.is_some_and(|p| view.rect.contains(p)) {
                self.picking_face = false;
                match Self::face_plane(session.document(), &view) {
                    Ok((plane, label)) => {
                        self.definition.plane = plane;
                        self.face_label = Some(label);
                        self.changed();
                        return true;
                    }
                    Err(e) => self.error = Some(e),
                }
            }
            return false;
        }
        let parameters = session.parameters();
        let frame = self.plane_frame(parameters);
        let tolerance = (PICK_PX / view.camera.scale.max(1e-6)) as f64;
        let nodes_world = self
            .positions(&self.definition, parameters)
            .unwrap_or_default();
        let nodes: Vec<(EntityId, [f64; 2])> = nodes_world
            .iter()
            .map(|(id, p)| (*id, Self::to_plane(&frame, *p)))
            .collect();
        let inside = view.pointer.filter(|p| view.rect.contains(*p));
        let raw = inside.and_then(|p| view.camera.on_plane(p, view.rect, frame));
        let last = self
            .chain
            .last()
            .and_then(|(n, _)| nodes.iter().find(|(id, _)| id == n))
            .map(|(_, p)| *p);
        let at = |id: EntityId| nodes.iter().find(|(n, _)| *n == id).map(|(_, p)| *p);
        let directions: Vec<f64> = self
            .definition
            .members
            .iter()
            .filter_map(|m| {
                let (a, b) = (at(m.start)?, at(m.end)?);
                Some((b[1] - a[1]).portable_atan2(b[0] - a[0]).to_degrees())
            })
            .collect();
        // A node being moved never snaps onto itself.
        let targets: Vec<(EntityId, [f64; 2])> = nodes
            .iter()
            .filter(|(id, _)| Some(*id) != self.moving)
            .copied()
            .collect();
        let anchor = if self.drawing { last } else { None };
        let snapped = raw.map(|r| self.snap(r, &targets, anchor, &directions, tolerance));
        self.hover = snapped.map(|s| s.0);
        let mut changed = false;
        // Dragging a node: press on it, move, release to place it. A press
        // that does not move stays a click that selects the node.
        let (pressed, released) = view
            .painter
            .ctx()
            .input(|i| (i.pointer.primary_pressed(), i.pointer.primary_released()));
        if pressed
            && !self.drawing
            && self.moving.is_none()
            && inside.is_some()
            && let Some((point, Some(node))) = snapped
        {
            self.moving = Some(node);
            self.dragging = true;
            self.drag_from = Some(point);
        }
        let mut consumed = false;
        if self.dragging && released {
            self.dragging = false;
            let from = self.drag_from.take();
            let node = self.moving.take();
            if let (Some(node), Some(from)) = (node, from) {
                match self.hover {
                    Some(to) if (to[0] - from[0]).portable_hypot(to[1] - from[1]) >= 0.5 => {
                        changed |= self.place(node, to, parameters);
                    }
                    // A press that did not move is a click on the node.
                    _ => self.select_node(node, parameters),
                }
                consumed = true;
            }
        }
        if view.clicked && inside.is_some() && !consumed {
            match (raw, snapped) {
                (None, _) | (_, None) => {
                    self.error = Some(
                        "That click does not hit the drawing plane; turn the view toward it or choose another plane"
                            .into(),
                    );
                }
                (Some(_), Some((point, hit))) => {
                    if let Some(node) = self.moving.take() {
                        changed |= self.place(node, point, parameters);
                    } else if self.drawing {
                        match self.draw_click(point, hit) {
                            Ok(c) => changed |= c,
                            Err(e) => self.error = Some(e),
                        }
                    } else if let Some(node) = hit {
                        self.select_node(node, parameters);
                    } else {
                        self.selected_member = self.member_at(point, &nodes, tolerance);
                        self.selected_node = None;
                    }
                }
            }
        }
        if changed {
            self.changed();
        }
        let probe = self.moving.zip(self.hover);
        self.refresh_report(session, probe);
        self.paint(&view, &frame, &nodes);
        changed
    }
    fn select_node(&mut self, node: EntityId, parameters: &EvaluatedParameters) {
        self.selected_node = Some(node);
        self.selected_member = None;
        self.gap_text = self
            .definition
            .joint(node)
            .and_then(|j| j.gap.as_ref())
            .map(|g| {
                format!(
                    "{} mm",
                    g.value_in(Unit::Millimeter, Some(parameters)).unwrap_or(0.)
                )
            })
            .unwrap_or_default();
    }
    /// The plane of the planar model face under the pointer: its semantic
    /// face frame, with the frame's u axis.
    fn face_plane(
        document: &Document,
        view: &TerminationViewport<'_>,
    ) -> Result<(FramePlane, String), String> {
        let scene = view
            .scene
            .ok_or("Wait for the current model before picking a face")?;
        let pointer = view.pointer.ok_or("Point at a face")?;
        let ray = view
            .camera
            .ray(pointer, view.rect)
            .ok_or("Point at a face")?;
        let hit = scene
            .pick_face_section(ray, view.section)
            .ok_or("That click hits no model face; click a flat face of a body")?;
        let occurrence = scene
            .instances()
            .iter()
            .find(|i| i.id == hit.instance && i.geometry == hit.geometry)
            .ok_or("That face is no longer shown")?;
        if occurrence.transform != Default::default() {
            return Err("Faces of placed occurrences cannot be drawn on yet; pick a face of an unplaced body".into());
        }
        let snapshot = scene
            .geometry_snapshot(hit.geometry)
            .ok_or("The face has no current topology; wait for the model to rebuild")?;
        let frame = snapshot
            .elements
            .iter()
            .find(|e| e.key == hit.key)
            .and_then(|e| e.frames.first())
            .ok_or("Choose a flat face; curved faces cannot hold a frame plane")?
            .frame;
        frame.validate()?;
        let n = frame.normal();
        let name = document
            .entities()
            .find_map(|e| match &e.data {
                multicad_core::document::EntityData::Part { features, .. }
                    if features.contains(&hit.geometry) =>
                {
                    Some(e.name.clone())
                }
                _ => None,
            })
            .unwrap_or_else(|| "a body".into());
        Ok((
            FramePlane {
                origin: frame.origin_mm.map(length),
                u: frame.u,
                normal: n,
            },
            format!("Face of {name}"),
        ))
    }
    /// Rewrite a node's own position to an absolute plane point; nodes that
    /// depend on it follow.
    fn place(&mut self, node: EntityId, point: [f64; 2], parameters: &EvaluatedParameters) -> bool {
        let position = NodePosition::Absolute {
            u: length(point[0]),
            v: length(point[1]),
            w: None,
        };
        let mut probe = self.definition.clone();
        match probe
            .set_node(node, position)
            .and_then(|_| probe.evaluate(parameters).map(|_| ()))
        {
            Ok(()) => {
                self.definition = probe;
                true
            }
            Err(e) => {
                self.error = Some(e.to_string());
                false
            }
        }
    }
    fn member_at(
        &self,
        point: [f64; 2],
        nodes: &[(EntityId, [f64; 2])],
        tolerance: f64,
    ) -> Option<EntityId> {
        let at = |id: EntityId| nodes.iter().find(|(n, _)| *n == id).map(|(_, p)| *p);
        self.definition
            .members
            .iter()
            .filter_map(|m| {
                let (a, b) = (at(m.start)?, at(m.end)?);
                let (dx, dy) = (b[0] - a[0], b[1] - a[1]);
                let l2 = dx * dx + dy * dy;
                if l2 <= 0. {
                    return None;
                }
                let t = (((point[0] - a[0]) * dx + (point[1] - a[1]) * dy) / l2).clamp(0., 1.);
                let d = (point[0] - a[0] - t * dx).portable_hypot(point[1] - a[1] - t * dy);
                (d <= tolerance).then_some((m.id, d))
            })
            .min_by(|a, b| a.1.total_cmp(&b.1))
            .map(|(id, _)| id)
    }
    fn paint(&self, view: &TerminationViewport<'_>, frame: &Frame, nodes: &[(EntityId, [f64; 2])]) {
        let world = |p: [f64; 2]| {
            std::array::from_fn::<f32, 3, _>(|i| {
                (frame.origin[i] + p[0] * frame.u[i] + p[1] * frame.v[i]) as f32
            })
        };
        let screen = |p: [f64; 2]| {
            let s = view.camera.project(world(p), view.rect);
            egui::pos2(s[0], s[1])
        };
        let at = |id: EntityId| {
            if self.moving == Some(id) {
                return self.hover;
            }
            nodes.iter().find(|(n, _)| *n == id).map(|(_, p)| *p)
        };
        let line = egui::Color32::from_rgb(64, 196, 255);
        let accent = egui::Color32::from_rgb(255, 176, 64);
        let bad = egui::Color32::from_rgb(255, 80, 80);
        let font = egui::FontId::proportional(11.);
        // The real placed section at both ends of every member, joined
        // along its edges: the section preview while drawing.
        let section = egui::Stroke::new(
            1.,
            egui::Color32::from_rgba_unmultiplied(160, 200, 230, 150),
        );
        let to_screen = |p: &[f64; 3]| {
            let s = view.camera.project(p.map(|v| v as f32), view.rect);
            egui::pos2(s[0], s[1])
        };
        for [a, b] in &self.sections {
            for outline in [a, b] {
                let points: Vec<egui::Pos2> = outline.iter().map(to_screen).collect();
                view.painter.add(egui::Shape::closed_line(points, section));
            }
            let step = (a.len() / 8).max(1);
            for (p, q) in a.iter().zip(b).step_by(step) {
                view.painter
                    .line_segment([to_screen(p), to_screen(q)], section);
            }
        }
        for m in &self.definition.members {
            let (Some(a), Some(b)) = (at(m.start), at(m.end)) else {
                continue;
            };
            let selected = self.selected_member == Some(m.id);
            let color = if selected { accent } else { line };
            let (sa, sb) = (screen(a), screen(b));
            view.painter.line_segment(
                [sa, sb],
                egui::Stroke::new(if selected { 3. } else { 2. }, color),
            );
            let l = (b[0] - a[0]).portable_hypot(b[1] - a[1]);
            let dir = (sb - sa).normalized();
            view.painter.text(
                sa + (sb - sa) / 2. + egui::vec2(-dir.y, dir.x) * 12.,
                egui::Align2::CENTER_CENTER,
                format!("{} mm", format_number(l.round())),
                font.clone(),
                color,
            );
        }
        for (id, p) in nodes {
            let p = if self.moving == Some(*id) {
                self.hover.unwrap_or(*p)
            } else {
                *p
            };
            let s = screen(p);
            let selected = self.selected_node == Some(*id);
            view.painter.circle_filled(
                s,
                if selected { 5. } else { 3.5 },
                if selected { accent } else { line },
            );
            if let Some(joint) = self.joint(*id) {
                let (text, refused) = Self::badge(joint);
                let pos = s + egui::vec2(10., -14.);
                let galley = view
                    .painter
                    .layout_no_wrap(text, font.clone(), egui::Color32::WHITE);
                let rect = egui::Rect::from_min_size(pos, galley.size()).expand(3.);
                view.painter.rect_filled(
                    rect,
                    3.,
                    if refused {
                        bad
                    } else {
                        egui::Color32::from_rgb(40, 60, 90)
                    },
                );
                view.painter.galley(pos, galley, egui::Color32::WHITE);
            }
        }
        // The pending segment from the last chain node.
        if self.drawing
            && let (Some(hover), Some((last, _))) = (self.hover, self.chain.last())
            && let Some(l) = at(*last)
        {
            view.painter
                .line_segment([screen(l), screen(hover)], egui::Stroke::new(1.5, accent));
            let length = (hover[0] - l[0]).portable_hypot(hover[1] - l[1]);
            view.painter.text(
                screen(hover) + egui::vec2(12., 12.),
                egui::Align2::LEFT_TOP,
                format!("{} mm", format_number(length.round())),
                font,
                accent,
            );
        }
    }

    // --- Inspector ---------------------------------------------------------------

    /// The stock a frame member is made of: requested, or its feature's.
    fn member_stock(&self, document: &Document, member: EntityId) -> Option<EntityId> {
        if let Some(r) = self.new_members.get(&member) {
            return Some(r.stock);
        }
        let feature = model::frame_member_feature(document, self.frame, member)?;
        model::read_member(document, feature).ok().map(|m| m.stock)
    }
    /// The selected member's section with its offset and roll as callouts
    /// linked to their fields (CAD-395).
    fn placement_diagram(
        &mut self,
        ui: &mut egui::Ui,
        document: &Document,
        parameters: &EvaluatedParameters,
        member: EntityId,
    ) -> bool {
        use super::diagram::{DiagramLink, Notation};
        use crate::materials::profile_picker as picker;
        const ACROSS: &str = "frame.offset.across";
        const UP: &str = "frame.offset.up";
        const ROLL: &str = "frame.roll";
        let Some(spec) = self
            .member_stock(document, member)
            .and_then(|id| document.get(id))
            .and_then(|e| multicad_core::materials::read(e).ok().flatten())
            .and_then(|d| match d {
                multicad_core::materials::LibraryDefinition::Stock(s) => Some(s.specification),
                _ => None,
            })
        else {
            return false;
        };
        let Some(section) = picker::section(&spec) else {
            return false;
        };
        let Some(m) = self.definition.member(member).cloned() else {
            return false;
        };
        let mm = |d: &multicad_core::constraint_sketch::Dimension, unit| {
            d.value_in(unit, Some(parameters)).unwrap_or(0.)
        };
        let mut across = m
            .offset
            .as_ref()
            .map_or(0., |o| mm(&o[0], Unit::Millimeter));
        let mut up = m
            .offset
            .as_ref()
            .map_or(0., |o| mm(&o[1], Unit::Millimeter));
        let mut roll = m.roll.as_ref().map_or(0., |r| mm(r, Unit::Degree));
        let (w, h) = (section.width_mm as f32, section.height_mm as f32);
        let size = w.max(h);
        let diagram = picker::art(&section)
            .leader(
                ACROSS,
                [w, h / 2.],
                [w + 0.35 * size, h * 0.2],
                Notation::Length,
                "Offset across",
            )
            .leader(
                UP,
                [w / 2., h],
                [-0.45 * size, h + 0.25 * size],
                Notation::Length,
                "Offset up",
            )
            .leader(
                ROLL,
                [w * 0.85, h * 0.85],
                [w + 0.3 * size, h + 0.3 * size],
                Notation::Angle,
                "Roll",
            );
        let id = ui.make_persistent_id("frame-placement-link");
        let mut link: DiagramLink = ui.data(|d| d.get_temp(id)).unwrap_or_default();
        let value = |key| match key {
            ACROSS => format!("{} mm", format_number(across)),
            UP => format!("{} mm", format_number(up)),
            _ => format!("{}°", format_number(roll)),
        };
        link.diagram(ui, "frame-placement-diagram", &diagram, 160., &value);
        let mut changed = false;
        let signed = -1e6..=1e6;
        changed |= picker::linked_value(
            ui,
            &mut link,
            ACROSS,
            "Offset across",
            &mut across,
            "Shift the section along its width after placing it",
            signed.clone(),
            " mm",
        );
        changed |= picker::linked_value(
            ui,
            &mut link,
            UP,
            "Offset up",
            &mut up,
            "Shift the section along its height after placing it",
            signed,
            " mm",
        );
        changed |= picker::linked_value(
            ui,
            &mut link,
            ROLL,
            "Roll",
            &mut roll,
            "Turn the section about the member axis",
            -360.0..=360.0,
            "°",
        );
        ui.data_mut(|d| d.insert_temp(id, link));
        if changed && let Some(target) = self.definition.members.iter_mut().find(|x| x.id == member)
        {
            target.offset = (across != 0. || up != 0.).then(|| [length(across), length(up)]);
            target.roll = (roll != 0.).then(|| angle(roll));
        }
        changed
    }

    pub fn show(
        &mut self,
        ui: &mut egui::Ui,
        document: &Document,
        enabled: bool,
        preview: PreviewState<'_>,
    ) -> Option<ModelEvent> {
        let parameters = &match multicad_core::parameters::evaluate(document) {
            Ok(p) => p,
            Err(e) => {
                ui.colored_label(ui.visuals().error_fg_color, e.to_string());
                return None;
            }
        };
        controls::inspector_layout(ui);
        ui.strong(if self.editing { "Edit frame" } else { "Frame" });
        ui.small("Draw profile members as lines. Joints are made at every node: equal sections mitre, different sections butt with the larger through, round tubes cope.");
        ui.label("Frame name");
        let mut changed = ui
            .add_enabled(
                !self.editing,
                egui::TextEdit::singleline(&mut self.name).desired_width(ui.available_width()),
            )
            .changed();
        // Plane: chosen before the first node.
        ui.label("Drawing plane");
        let choices = plane_choices(document);
        let current = choices
            .iter()
            .find(|(_, p)| *p == self.definition.plane)
            .map(|(label, _)| label.clone());
        let empty = self.definition.nodes.is_empty();
        ui.add_enabled_ui(empty, |ui| {
            controls::inspector_combo(ui, "frame-plane")
                .selected_text(
                    current
                        .clone()
                        .or_else(|| self.face_label.clone())
                        .unwrap_or_else(|| "Custom plane".into()),
                )
                .show_ui(ui, |ui| {
                    for (label, plane) in &choices {
                        if ui
                            .selectable_label(current.as_ref() == Some(label), label)
                            .clicked()
                        {
                            self.definition.plane = plane.clone();
                            self.face_label = None;
                            changed = true;
                        }
                    }
                });
            let label = if self.picking_face {
                "Click a planar face…"
            } else {
                "Pick a planar face"
            };
            if ui
                .button(label)
                .on_hover_text("Then click a flat face of a model body to draw on it")
                .clicked()
            {
                self.picking_face = !self.picking_face;
                self.drawing = false;
            }
        });
        // Profile stock with pictures and plain names.
        ui.label("Profile stock");
        let stocks = super::weldment::profile_stock_specs(document, false);
        let selected = self
            .stock
            .and_then(|id| stocks.iter().find(|(s, ..)| *s == id))
            .map_or("Choose profile stock".to_owned(), |(_, spec, name)| {
                crate::materials::profile_picker::description(spec).unwrap_or_else(|| name.clone())
            });
        let create = controls::inspector_combo(ui, "frame-stock")
            .selected_text(selected)
            .show_ui(ui, |ui| {
                ui.set_min_width(260.);
                for (id, spec, name) in &stocks {
                    let is = self.stock == Some(*id);
                    if crate::materials::profile_picker::stock_row(ui, is, spec, name).clicked()
                        && !is
                    {
                        self.stock = Some(*id);
                        changed = true;
                    }
                }
                ui.button("New profile stock…").clicked()
            })
            .inner
            .unwrap_or(false);
        if create {
            return Some(ModelEvent::CreateStock(
                crate::materials::ToolStock::Profile,
            ));
        }
        ui.label("Section placement");
        controls::inspector_combo(ui, "frame-side")
            .selected_text(side_label(self.definition.side))
            .show_ui(ui, |ui| {
                for (side, label) in SIDES {
                    if ui
                        .selectable_label(self.definition.side == side, label)
                        .clicked()
                    {
                        self.definition.side = side;
                        changed = true;
                    }
                }
            });
        ui.separator();
        // Selected member.
        if let Some(member) = self.selected_member
            && self.definition.member(member).is_some()
        {
            ui.separator();
            ui.strong(self.member_name(member));
            let own = self.definition.member(member).and_then(|m| m.side);
            controls::inspector_combo(ui, "frame-member-side")
                .selected_text(own.map_or("Frame placement", side_label))
                .show_ui(ui, |ui| {
                    if ui
                        .selectable_label(own.is_none(), "Frame placement")
                        .clicked()
                    {
                        if let Some(m) = self.definition.members.iter_mut().find(|m| m.id == member)
                        {
                            m.side = None;
                        }
                        changed = true;
                    }
                    for (side, label) in SIDES {
                        if ui.selectable_label(own == Some(side), label).clicked() {
                            if let Some(m) =
                                self.definition.members.iter_mut().find(|m| m.id == member)
                            {
                                m.side = Some(side);
                            }
                            changed = true;
                        }
                    }
                });
            ui.horizontal(|ui| {
                if ui.button("Reverse").clicked() {
                    self.reverse_member(member);
                    changed = true;
                }
                if ui.button("Delete member").clicked() {
                    match self.definition.remove_member(member) {
                        Ok(()) => {
                            self.new_members.remove(&member);
                            self.selected_member = None;
                            changed = true;
                        }
                        Err(e) => self.error = Some(e.to_string()),
                    }
                }
            });
            changed |= self.placement_diagram(ui, document, parameters, member);
        }
        // Joint editor for the selected node.
        if let Some(node) = self.selected_node
            && self.definition.node(node).is_some()
        {
            ui.separator();
            ui.strong("Joint");
            if ui
                .button("Move node")
                .on_hover_text("Then click the node's new place on the plane")
                .clicked()
            {
                self.moving = Some(node);
            }
            if let Some(problem) = self.joint(node).and_then(|j| j.problem.clone()) {
                ui.colored_label(ui.visuals().error_fg_color, problem);
            }
            let mut joint = self
                .definition
                .joint(node)
                .cloned()
                .unwrap_or(JointOverride {
                    node,
                    ..Default::default()
                });
            let before = joint.clone();
            ui.horizontal_wrapped(|ui| {
                for (kind, label) in [
                    (None, "Automatic"),
                    (Some(ConnectionKind::Miter), "Mitre"),
                    (Some(ConnectionKind::Butt), "Butt"),
                    (Some(ConnectionKind::Cope), "Cope"),
                ] {
                    if ui.selectable_label(joint.kind == kind, label).clicked() {
                        joint.kind = kind;
                        if kind == Some(ConnectionKind::Miter) {
                            joint.through = None;
                        }
                    }
                }
            });
            let incident: Vec<EntityId> = self
                .definition
                .members
                .iter()
                .filter(|m| m.start == node || m.end == node)
                .map(|m| m.id)
                .collect();
            ui.label("Through member");
            controls::inspector_combo(ui, "frame-joint-through")
                .selected_text(joint.through.map_or("Automatic".to_owned(), |m| {
                    format!("{} through", self.member_name(m))
                }))
                .show_ui(ui, |ui| {
                    if ui
                        .selectable_label(joint.through.is_none(), "Automatic")
                        .clicked()
                    {
                        joint.through = None;
                    }
                    for m in &incident {
                        let label = format!("{} through", self.member_name(*m));
                        if ui
                            .selectable_label(joint.through == Some(*m), label)
                            .clicked()
                        {
                            joint.through = Some(*m);
                        }
                    }
                });
            ui.label("Weld gap (mm, empty = frame default 0 mm)");
            let gap = ui.add(
                egui::TextEdit::singleline(&mut self.gap_text)
                    .id_salt("frame-joint-gap")
                    .hint_text("0"),
            );
            if gap.changed() {
                joint.gap = if self.gap_text.trim().is_empty() {
                    None
                } else {
                    let mut field = self.segment_length.clone();
                    field.text = self.gap_text.clone();
                    field.resolve(document).ok()
                };
            }
            if joint != before {
                self.definition.set_joint(joint);
                changed = true;
            }
        }
        if self.moving.is_some() {
            ui.small("Click the new position on the plane; Esc keeps the node where it is.");
        }
        ui.separator();
        // Drawing.
        let draw_label = if self.drawing {
            "Stop drawing"
        } else {
            "Draw in viewport"
        };
        if ui.button(draw_label).clicked() {
            self.drawing = !self.drawing;
            self.chain.clear();
        }
        if self.drawing {
            ui.small("Click nodes on the plane; click the first node to close the loop. Backspace removes the last node, Esc cancels the pending segment, Enter ends the chain.");
        }
        if self.chain.is_empty() {
            ui.horizontal(|ui| {
                changed |= self.start[0].show_compact(
                    ui,
                    "frame-start-u",
                    "Start u",
                    "First node, along the plane's u axis",
                );
                changed |= self.start[1].show_compact(
                    ui,
                    "frame-start-v",
                    "Start v",
                    "First node, along the plane's v axis",
                );
            });
        }
        ui.label(format!("Segment length ({})", self.segment_length.unit()));
        let length_response = ui.add(
            egui::TextEdit::singleline(&mut self.segment_length.text)
                .id_salt("frame-segment-length")
                .hint_text("1000, 1 m, 39 in or =$width"),
        );
        ui.label("Turn from the previous member (deg, counterclockwise)");
        let angle_response = ui.add(
            egui::TextEdit::singleline(&mut self.segment_angle.text)
                .id_salt("frame-segment-angle")
                .hint_text("90"),
        );
        // Enter in the length field adds the segment; the turn is set first.
        let entered = length_response.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter));
        let _ = angle_response;
        ui.horizontal(|ui| {
            if ui.button("Add segment").clicked() || entered {
                match self.add_segment(document, parameters) {
                    Ok(()) => changed = true,
                    Err(e) => self.error = Some(e),
                }
            }
            if ui
                .add_enabled(self.chain.len() >= 3, egui::Button::new("Close loop"))
                .clicked()
            {
                match self.close_loop() {
                    Ok(()) => changed = true,
                    Err(e) => self.error = Some(e),
                }
            }
        });
        ui.small(format!(
            "{} members · {} nodes",
            self.definition.members.len(),
            self.definition.nodes.len()
        ));
        ui.separator();
        form::preview_status(ui, self.stamp, preview);
        let validation = self.request(document, enabled);
        if let Err(error) = &validation {
            ui.colored_label(ui.visuals().error_fg_color, error);
        }
        if let Some((_, Err(error))) = &self.report
            && validation.is_ok()
        {
            ui.colored_label(ui.visuals().error_fg_color, error);
        }
        let action = form::actions(ui, validation.is_ok(), changed, self.stamp, preview);
        if let Some(error) = &self.error {
            ui.colored_label(ui.visuals().error_fg_color, error);
        }
        if action == Some(FormAction::Cancel) {
            return Some(ModelEvent::Cancel);
        }
        if changed {
            self.changed();
            if action != Some(FormAction::Preview) {
                return Some(ModelEvent::Changed);
            }
        }
        match action {
            Some(FormAction::Preview) => match self.model_request(document, enabled) {
                Ok(request) => Some(ModelEvent::Preview {
                    stamp: self.stamp,
                    request,
                }),
                Err(error) => {
                    self.error = Some(error);
                    None
                }
            },
            Some(FormAction::Apply)
                if preview.ready(self.stamp) && self.request(document, enabled).is_ok() =>
            {
                Some(ModelEvent::Apply)
            }
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::modeling::diagram::harness::{self, WIDTHS};
    use multicad_core::{
        commands::Channel,
        document::{Entity, EntityData},
        materials::{LibraryDefinition, ProfileShape, Stock, StockKind, prepare_library},
    };

    fn session() -> Session {
        let mut s = Session::new(Document::new("Frame panel").unwrap()).unwrap();
        let steel = Entity::new(
            "Steel",
            EntityData::Material {
                density_kg_m3: Some(7850.),
                color_rgba: [120, 124, 130, 255],
                section_hatch: None,
            },
        );
        let stock = LibraryDefinition::Stock(Stock {
            material: steel.id,
            specification: StockKind::Profile {
                designation: "L 40×40×4".into(),
                length_mm: 6000.,
                width_mm: 40.,
                height_mm: 40.,
                wall_mm: Some(4.),
                shape: ProfileShape::Angle,
                section: None,
                corner_radius_mm: None,
                web_mm: None,
            },
        })
        .entity(EntityId::new(), "L 40×40×4")
        .unwrap();
        s.commit(prepare_library(&s, 0, Actor::local(Channel::Ui), &[steel, stock]).unwrap())
            .unwrap();
        s
    }

    #[test]
    fn the_placement_diagram_links_offset_and_roll_callouts_to_their_fields() {
        let session = session();
        let document = session.document().clone();
        for width in WIDTHS {
            let mut panel = FramePanel::new(&document, "Frame".into()).unwrap();
            panel.add_segment(&document, session.parameters()).unwrap();
            panel.selected_member = Some(panel.definition.members[0].id);
            let ctx = egui::Context::default();
            let mut show = |ui: &mut egui::Ui| {
                panel.show(ui, &document, true, PreviewState::Idle);
            };
            harness::edit_through_callout(&ctx, width, "Offset across 0 mm", "-15", &mut show);
            harness::edit_through_callout(&ctx, width, "∠ Roll 0°", "90", &mut show);
            let m = &panel.definition.members[0];
            assert_eq!(
                m.offset.as_ref().map(|o| o[0].clone()),
                Some(length(-15.)),
                "{width}"
            );
            assert_eq!(m.roll, Some(angle(90.)), "{width}");
        }
    }
}
