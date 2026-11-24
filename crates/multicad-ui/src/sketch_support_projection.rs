//! Snap-time projection of the supporting face (CAD-383). While a draw tool
//! is active in a face-attached sketch, the corners and edges of that face are
//! offered as snap targets. When the gesture commits, each corner or edge it
//! placed a point on is captured as a canonical projected reference (the
//! CAD-128 mechanism): a corner's output point is reused by identity and a
//! point on an edge lies exactly on it. The captures are prepared separately
//! against the same base, merged into one candidate, extended with the drawn
//! geometry and entered sizes, checked by the asynchronous native candidate
//! validation and committed once. Nothing changes on failure or cancellation.
//!
//! The face and the sketch are both in the producer's local coordinates, so
//! the capture is exact whatever pose the part or occurrence has.
use super::*;
use inline_input::Field;
use multicad_core::{
    constraint_sketch::Dimension,
    references::{
        self, AnalyticGeometry, GeometryIntent, ReferenceReport, Resolution, TopologyKey,
        TopologyKind,
    },
    sketch_draw::{DrawInput, DrawRequest, PointInput},
    sketch_external::{self, Behavior, CaptureInput, CaptureRequest, Method},
};
use snap::SupportEdge;

/// Draw tools whose point placements may land on a face corner or edge.
/// Fillet, Chamfer and Text place no points.
pub(super) fn projects(tool: DrawTool) -> bool {
    !matches!(tool, DrawTool::Fillet | DrawTool::Chamfer | DrawTool::Text)
}

/// A polyline's first vertex (position, point) and its committed segments.
type PolylineChain = (Option<([f64; 2], EntityId)>, usize);

/// A drawing onto face references, waiting for its native check.
pub(super) struct Pending {
    prepared: Prepared,
    requested: bool,
    validated: Option<Result<(), String>>,
    created: Vec<EntityId>,
    /// The solved candidate sketch, to continue a polyline after commit.
    solved: Option<multicad_core::constraint_sketch::Sketch>,
    /// A polyline's first vertex and segment count, to continue the chain.
    polyline: Option<PolylineChain>,
}

/// The point inputs of a draw request a placement may sit on.
fn point_inputs(input: &mut DrawInput) -> Vec<&mut PointInput> {
    match input {
        DrawInput::BSpline { poles, .. } => poles.iter_mut().collect(),
        DrawInput::Polyline { vertices, .. } => vertices.iter_mut().collect(),
        DrawInput::Line { start, end } => vec![start, end],
        DrawInput::Rectangle { a, b } => vec![a, b],
        DrawInput::Circle { center, .. } => vec![center],
        DrawInput::Arc {
            center, start, end, ..
        } => vec![center, start, end],
        DrawInput::CenterRectangle { corner, .. }
        | DrawInput::CenterPointsRectangle { corner, .. } => vec![corner],
        DrawInput::ThreePointRectangle { a, b, .. } => vec![a, b],
        DrawInput::TwoPointCircle { .. } | DrawInput::ThreePointCircle { .. } => vec![],
    }
}

/// The explicit refusal for a captured face reference that does not resolve
/// to exactly one current native element.
pub(super) fn resolution_error(resolution: &Resolution) -> Option<String> {
    match resolution {
        Resolution::Resolved { .. } => None,
        Resolution::Missing | Resolution::MissingSource => Some(
            "Missing face reference: the snapped corner or edge has no current native geometry"
                .into(),
        ),
        Resolution::Ambiguous { candidates } => Some(format!(
            "Ambiguous face reference ({} projections); choose one with Project reference",
            candidates.len()
        )),
        other => Some(format!("The face reference does not resolve: {other:?}")),
    }
}

impl SketchEditor {
    /// The face producer and its current native topology for a sketch on a
    /// face, when the display scene is current.
    fn support_source(
        &self,
        document: &Document,
        revision: u64,
        sketch: EntityId,
    ) -> Option<(EntityId, references::TopologySnapshot)> {
        let Some(EntityData::Sketch { plane, .. }) = document.get(sketch).map(|e| &e.data) else {
            return None;
        };
        let Some(EntityData::Plane {
            anchor: Some(link), ..
        }) = document.get(*plane).map(|e| &e.data)
        else {
            return None;
        };
        link.geometry.as_ref()?;
        let producer = link.target;
        if !matches!(
            document.get(producer).map(|e| &e.data),
            Some(EntityData::Feature {
                suppressed: false,
                ..
            })
        ) {
            return None;
        }
        let scene = self.display_scene.as_deref().filter(|scene| {
            scene.stamp().document == document.id() && scene.stamp().revision == revision
        })?;
        Some((producer, scene.geometry_snapshot(producer)?))
    }
    /// Vertices of the face producer that lie in the sketch plane, in sketch
    /// coordinates: the corners a Coplanar projection captures exactly.
    pub(super) fn support_corners(
        &self,
        session: &Session,
        sketch: EntityId,
    ) -> Vec<(TopologyKey, [f64; 2])> {
        if self.tool.is_none_or(|tool| !projects(tool)) {
            return vec![];
        }
        let revision = session.revision();
        let Some((_, snapshot)) = self.support_source(session.document(), revision, sketch) else {
            return vec![];
        };
        let Ok((_, frame, _)) = self.read_source(session.document(), sketch, revision) else {
            return vec![];
        };
        snapshot
            .elements
            .iter()
            .filter_map(|element| match element.geometry {
                Some(AnalyticGeometry::Point { position_mm }) => {
                    let off: f64 = (0..3)
                        .map(|i| (position_mm[i] - frame.origin[i]) * frame.n[i])
                        .sum();
                    (off.abs() <= 1e-6 && GeometryIntent::capture(&snapshot, element.key).is_ok())
                        .then(|| (element.key, frame.coordinates(position_mm)))
                }
                _ => None,
            })
            .collect()
    }
    /// Straight and circular edges of the face producer that lie in the
    /// sketch plane, in sketch coordinates: the edges a Coplanar projection
    /// captures exactly.
    pub(super) fn support_edges(
        &self,
        session: &Session,
        sketch: EntityId,
    ) -> Vec<(TopologyKey, SupportEdge)> {
        if self.tool.is_none_or(|tool| !projects(tool)) {
            return vec![];
        }
        let revision = session.revision();
        let Some((_, snapshot)) = self.support_source(session.document(), revision, sketch) else {
            return vec![];
        };
        let Ok((_, frame, _)) = self.read_source(session.document(), sketch, revision) else {
            return vec![];
        };
        let dot = |a: [f64; 3], b: [f64; 3]| (0..3).map(|i| a[i] * b[i]).sum::<f64>();
        let on = |p: [f64; 3]| {
            dot(std::array::from_fn(|i| p[i] - frame.origin[i]), frame.n).abs() <= 1e-6
        };
        snapshot
            .elements
            .iter()
            .filter(|element| element.key.kind == TopologyKind::Edge)
            .filter_map(|element| {
                let edge = match *element.geometry.as_ref()? {
                    AnalyticGeometry::Line { a_mm, b_mm } => (on(a_mm) && on(b_mm)).then(|| {
                        SupportEdge::Line(frame.coordinates(a_mm), frame.coordinates(b_mm))
                    })?,
                    AnalyticGeometry::CircularArc {
                        frame: arc,
                        radius_mm,
                        start_rad,
                        sweep_rad,
                    } => {
                        let normal = [
                            arc.u[1] * arc.v[2] - arc.u[2] * arc.v[1],
                            arc.u[2] * arc.v[0] - arc.u[0] * arc.v[2],
                            arc.u[0] * arc.v[1] - arc.u[1] * arc.v[0],
                        ];
                        let along = dot(normal, frame.n);
                        if !on(arc.origin_mm) || along.abs() < 1. - 1e-9 || radius_mm <= 1e-9 {
                            return None;
                        }
                        let center = frame.coordinates(arc.origin_mm);
                        let at = |angle: f64| {
                            let (s, c) = (angle.portable_sin(), angle.portable_cos());
                            frame.coordinates(std::array::from_fn(|i| {
                                arc.origin_mm[i] + radius_mm * (c * arc.u[i] + s * arc.v[i])
                            }))
                        };
                        let angle =
                            |p: [f64; 2]| (p[1] - center[1]).portable_atan2(p[0] - center[0]);
                        let sweep = sweep_rad.abs().min(std::f64::consts::TAU);
                        // Counter-clockwise in the sketch, whatever the arc frame.
                        let start = if sweep >= std::f64::consts::TAU - 1e-9 {
                            0.
                        } else if along * sweep_rad > 0. {
                            angle(at(start_rad))
                        } else {
                            angle(at(start_rad + sweep_rad))
                        };
                        SupportEdge::Arc {
                            center,
                            radius: radius_mm,
                            start,
                            sweep,
                        }
                    }
                    _ => return None,
                };
                GeometryIntent::capture(&snapshot, element.key)
                    .is_ok()
                    .then_some((element.key, edge))
            })
            .collect()
    }
    /// One candidate holding a projected reference per face corner or edge.
    /// Each capture is prepared against the current document and the
    /// captures are merged into one transaction, so none can change another
    /// or existing source. Returns the candidate, its reference report and,
    /// per key in order, the projected point of a corner.
    fn capture_support(
        &self,
        session: &Session,
        sketch: EntityId,
        keys: &[TopologyKey],
        channel: Channel,
    ) -> Result<(Prepared, ReferenceReport, Vec<Option<EntityId>>), String> {
        let revision = session.revision();
        let (producer, snapshot) = self
            .support_source(session.document(), revision, sketch)
            .ok_or("The supporting face's current native geometry is unavailable")?;
        let snapshots = self
            .display_scene
            .as_ref()
            .map(|scene| scene.reference_snapshots())
            .unwrap_or_default();
        let document = session.document();
        let mut previews = Vec::new();
        let mut captured: Vec<(TopologyKey, EntityId, Option<EntityId>)> = Vec::new();
        for key in keys {
            if captured.iter().any(|(k, ..)| k == key) {
                continue;
            }
            if key.kind == TopologyKind::Face {
                return Err("Only a corner or an edge of the supporting face is projected".into());
            }
            if !snapshot.elements.iter().any(|element| element.key == *key) {
                return Err(resolution_error(&Resolution::Missing).expect("missing"));
            }
            let intent = GeometryIntent::capture(&snapshot, *key)
                .map_err(|e| format!("Unsupported face reference: {e}"))?;
            let request = CaptureRequest::new(
                document,
                revision,
                CaptureInput {
                    sketch,
                    source: producer,
                    intent,
                    method: Method::Coplanar,
                    behavior: Behavior::ReferenceOnly,
                },
                &snapshots,
            )
            .map_err(|e| e.to_string())?;
            let point = if key.kind == TopologyKind::Vertex {
                Some(
                    *request
                        .output
                        .points()
                        .first()
                        .ok_or("The face corner has no projected point")?,
                )
            } else {
                None
            };
            let step = sketch_external::capture_edits(document, revision, &request, &snapshots)
                .map_err(|e| e.to_string())?;
            previews.push(document.edited(&step).map_err(|e| e.to_string())?);
            captured.push((*key, request.reference, point));
        }
        let previews: Vec<&Document> = previews.iter().collect();
        let edits = sketch_external::merge_capture_previews(document, sketch, &previews)
            .map_err(|e| e.to_string())?;
        let prepared = session
            .prepare(Request {
                expected_revision: revision,
                actor: Actor::local(channel),
                label: "Project face reference".into(),
                edits,
            })
            .map_err(|e| e.to_string())?;
        let report = references::inspect(prepared.document(), revision, &snapshots);
        for (_, reference, _) in &captured {
            if let Some(error) = report
                .entries
                .iter()
                .filter(|entry| entry.address.owner == *reference)
                .find_map(|entry| resolution_error(&entry.resolution))
            {
                return Err(error);
            }
        }
        let points = keys
            .iter()
            .map(|key| {
                captured
                    .iter()
                    .find(|(k, ..)| k == key)
                    .and_then(|(.., point)| *point)
            })
            .collect();
        Ok((prepared, report, points))
    }
    /// Whether a draw request places a point on a face corner or edge that
    /// the current gesture snapped to.
    pub(super) fn support_requested(&self, request: &DrawRequest) -> bool {
        let mut input = request.input.clone();
        !self.support_used.is_empty()
            && point_inputs(&mut input).into_iter().any(|point| {
                matches!(point, PointInput::Position(p)
                    if self.support_used.iter().any(|(_, at)| at == p))
            })
    }
    /// Prepares a drawing whose point placements `used` sit on face corners
    /// or edges: they become projected references, a corner's point is
    /// reused, and the references, the drawing and the entered `sizes` form
    /// one candidate awaiting native validation.
    pub(super) fn draw_on_support(
        &mut self,
        session: &Session,
        request: &DrawRequest,
        used: &[(TopologyKey, [f64; 2])],
        sizes: Vec<(Field, Dimension)>,
        tool: Option<DrawTool>,
        channel: Channel,
    ) -> Result<(), String> {
        if !self.editing {
            return Err("Create or open a sketch first".into());
        }
        if self.has_pending_edits() {
            return Err("Apply or discard pending edits before drawing".into());
        }
        let sketch = self.selected.ok_or("Select a sketch")?;
        let mut request = request.clone();
        let mut keys: Vec<TopologyKey> = Vec::new();
        for point in point_inputs(&mut request.input) {
            if let PointInput::Position(p) = point
                && let Some((key, _)) = used.iter().find(|(_, at)| at == p)
                && !keys.contains(key)
            {
                keys.push(*key);
            }
        }
        let (base, report, points) = self.capture_support(session, sketch, &keys, channel)?;
        for point in point_inputs(&mut request.input) {
            if let PointInput::Position(p) = point
                && let Some((key, _)) = used.iter().find(|(_, at)| at == p)
                && let Some(Some(id)) = keys.iter().position(|k| k == key).map(|i| points[i])
            {
                *point = PointInput::Existing(id);
            }
        }
        let drawn = multicad_solver::draw::prepare_draw_on_with(
            session,
            session.revision(),
            Actor::local(channel),
            sketch,
            &base,
            Some(&report),
            &request,
            |source, created| inline_input::size_constraints(source, created, &sizes).map(|_| ()),
        )
        .map_err(|e| e.to_string())?;
        let prepared = drawn.edit.preview.ok_or_else(|| {
            format!(
                "Sketch rejected: {:?}; nothing was placed",
                drawn.edit.evaluation.diagnostics.state
            )
        })?;
        let polyline = (tool == Some(DrawTool::Polyline))
            .then_some((self.polyline_start, self.polyline_segments));
        self.projection = Some(Pending {
            prepared,
            requested: false,
            validated: None,
            created: drawn.created_curves,
            solved: drawn.edit.evaluation.candidate,
            polyline,
        });
        self.last_draw_line = None;
        self.clear_drawing_progress();
        self.message = "Projecting the face reference · checking it against the model…".into();
        Ok(())
    }
    /// Prepares the drawing whose `start`/`end` placements are on face
    /// corners: the corners become projected references and the drawing
    /// reuses their points, in one candidate awaiting native validation.
    #[cfg(test)]
    pub(super) fn draw_on_corners(
        &mut self,
        session: &Session,
        tool: DrawTool,
        start: (PointInput, Option<TopologyKey>),
        end: (PointInput, Option<TopologyKey>),
        radius_mm: f64,
        channel: Channel,
    ) -> Result<(), String> {
        let used: Vec<(TopologyKey, [f64; 2])> = [&start, &end]
            .into_iter()
            .filter_map(|(input, key)| match (input, key) {
                (PointInput::Position(p), Some(key)) => Some((*key, *p)),
                _ => None,
            })
            .collect();
        let (a, b) = (start.0, end.0);
        let input = match tool {
            DrawTool::Line => DrawInput::Line { start: a, end: b },
            DrawTool::Rectangle => DrawInput::Rectangle { a, b },
            DrawTool::Circle => DrawInput::Circle {
                center: a,
                radius_mm,
            },
            _ => return Err("This tool does not place points on face corners".into()),
        };
        let request = DrawRequest {
            input,
            construction: self.construction,
        };
        self.draw_on_support(session, &request, &used, vec![], Some(tool), channel)
    }
    pub fn projection_pending(&self) -> bool {
        self.projection.is_some()
    }
    /// Escape withdraws a queued or running projection; nothing is created.
    pub(super) fn cancel_projection(&mut self) -> bool {
        if self.projection.take().is_none() {
            return false;
        }
        self.external.cancel();
        self.message = "Face reference projection cancelled; nothing was placed".into();
        true
    }
    /// Requests and collects the native candidate validation.
    pub(super) fn poll_projection(
        &mut self,
        session: &Session,
        resources: &multicad_core::project_file::ProjectResources,
    ) {
        let Some(pending) = &mut self.projection else {
            return;
        };
        if session.validate_preview(&pending.prepared).is_err() {
            self.projection = None;
            self.external.cancel();
            self.message = "The document changed while projecting; nothing was placed".into();
            return;
        }
        if !pending.requested {
            pending.requested = true;
            if let Err(error) =
                self.external
                    .validate_candidate(session, &pending.prepared, resources)
            {
                pending.validated = Some(Err(error));
            }
            return;
        }
        if pending.validated.is_none()
            && let Some(result) = self.external.candidate_validation()
        {
            pending.validated = Some(result.clone());
        }
    }
    /// Commits a validated projection, or reports why nothing changed.
    /// Returns whether the projection finished (either way).
    pub fn finish_projection(&mut self, session: &mut Session) -> bool {
        let Some(result) = self.projection.as_ref().and_then(|p| p.validated.clone()) else {
            return false;
        };
        let pending = self.projection.take().expect("checked projection");
        let Pending {
            prepared,
            created,
            solved,
            polyline,
            ..
        } = pending;
        match result.and_then(|()| {
            session
                .commit(prepared)
                .map(|_| ())
                .map_err(|e| e.to_string())
        }) {
            Ok(()) => {
                self.element = created.first().copied();
                self.revision = session.revision();
                self.inline.texts.iter_mut().for_each(String::clear);
                self.message =
                    "Face reference projected and geometry committed · one Undo step".into();
                // A polyline continues from the end of the committed segment.
                let line = created
                    .last()
                    .zip(solved.as_ref())
                    .and_then(|(id, sketch)| {
                        let curve = sketch.curves.iter().find(|c| c.id == *id)?;
                        let multicad_core::constraint_sketch::CurveKind::Line { start, end } =
                            curve.kind
                        else {
                            return None;
                        };
                        let at = |p: EntityId| {
                            sketch
                                .points
                                .iter()
                                .find(|q| q.id == p)
                                .map(|q| q.position_mm)
                        };
                        Some([(start, at(start)?), (end, at(end)?)])
                    });
                if let Some((first, segments)) = polyline
                    && self.tool == Some(DrawTool::Polyline)
                    && let Some([start, end]) = line
                    && first.is_none_or(|(_, id)| id != end.0)
                {
                    self.polyline_start = first.or(Some((start.1, start.0)));
                    self.polyline_segments = segments + 1;
                    self.first = Some(end.1);
                    self.first_anchor = Some(end.0);
                }
            }
            Err(error) => {
                self.message = format!("{error}; nothing was placed");
            }
        }
        true
    }
}
