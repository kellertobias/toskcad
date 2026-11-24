//! Snap-time projection of the supporting face (CAD-383). While a two-point
//! draw tool is active in a face-attached sketch, the corners of that face are
//! offered as snap targets. Placing a point on one first captures the corner as
//! a canonical projected reference (the CAD-128 mechanism) and reuses its
//! output point; the capture and the drawn geometry are prepared as ONE
//! transaction, checked by the asynchronous native candidate validation and
//! committed only when it passes. Nothing changes on failure.
//!
//! The corner and the sketch are both in the producer's local coordinates, so
//! the capture is exact whatever pose the part or occurrence has.
use super::*;
use multicad_core::{
    references::{self, AnalyticGeometry, GeometryIntent, ReferenceReport, TopologyKey},
    sketch_draw::DrawRequest,
    sketch_external::{self, Behavior, CaptureInput, CaptureRequest, Method},
};

/// Draw tools whose point placements may land on a face corner.
pub(super) fn projects(tool: DrawTool) -> bool {
    matches!(
        tool,
        DrawTool::Line | DrawTool::Rectangle | DrawTool::Circle
    )
}

pub(super) struct Pending {
    prepared: Prepared,
    requested: bool,
    validated: Option<Result<(), String>>,
    created: Option<EntityId>,
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
    /// One candidate that captures each corner as a projected reference.
    /// Returns it with the reference report for the candidate document and
    /// the captured output point of each corner, in order.
    fn capture_corners(
        &self,
        session: &Session,
        sketch: EntityId,
        corners: &[TopologyKey],
        channel: Channel,
    ) -> Result<(Prepared, ReferenceReport, Vec<EntityId>), String> {
        let revision = session.revision();
        let (producer, snapshot) = self
            .support_source(session.document(), revision, sketch)
            .ok_or("The supporting face's current native geometry is unavailable")?;
        let snapshots = self
            .display_scene
            .as_ref()
            .map(|scene| scene.reference_snapshots())
            .unwrap_or_default();
        let mut document = session.document().clone();
        let mut edits = Vec::new();
        let mut captured: Vec<(TopologyKey, EntityId)> = Vec::new();
        for key in corners {
            if captured.iter().any(|(k, _)| k == key) {
                continue;
            }
            let intent = GeometryIntent::capture(&snapshot, *key)?;
            let request = CaptureRequest::new(
                &document,
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
            let point = *request
                .output
                .points()
                .first()
                .ok_or("The face corner has no projected point")?;
            let step = sketch_external::capture_edits(&document, revision, &request, &snapshots)
                .map_err(|e| e.to_string())?;
            document = document.edited(&step).map_err(|e| e.to_string())?;
            edits.extend(step);
            captured.push((*key, point));
        }
        let prepared = session
            .prepare(Request {
                expected_revision: revision,
                actor: Actor::local(channel),
                label: "Project face corner".into(),
                edits,
            })
            .map_err(|e| e.to_string())?;
        let report = references::inspect(prepared.document(), revision, &snapshots);
        let points = corners
            .iter()
            .map(|key| {
                captured
                    .iter()
                    .find(|(k, _)| k == key)
                    .map(|(_, id)| *id)
                    .expect("captured")
            })
            .collect();
        Ok((prepared, report, points))
    }
    /// Prepares the drawing whose `start`/`end` placements are on face
    /// corners: the corners become projected references and the drawing
    /// reuses their points, in one candidate awaiting native validation.
    pub(super) fn draw_on_corners(
        &mut self,
        session: &Session,
        tool: DrawTool,
        start: (multicad_core::sketch_draw::PointInput, Option<TopologyKey>),
        end: (multicad_core::sketch_draw::PointInput, Option<TopologyKey>),
        radius_mm: f64,
        channel: Channel,
    ) -> Result<(), String> {
        use multicad_core::sketch_draw::{DrawInput, PointInput};
        if self.has_pending_edits() {
            return Err("Apply or discard pending edits before drawing".into());
        }
        let sketch = self.selected.ok_or("Select a sketch")?;
        let corners: Vec<TopologyKey> = [start.1, end.1].into_iter().flatten().collect();
        let (base, report, points) = self.capture_corners(session, sketch, &corners, channel)?;
        let mut points = points.into_iter();
        let mut point = |input: PointInput, corner: Option<TopologyKey>| {
            corner
                .and_then(|_| points.next())
                .map(PointInput::Existing)
                .unwrap_or(input)
        };
        let (a, b) = (point(start.0, start.1), point(end.0, end.1));
        let input = match tool {
            DrawTool::Line => DrawInput::Line { start: a, end: b },
            DrawTool::Rectangle => DrawInput::Rectangle { a, b },
            DrawTool::Circle => DrawInput::Circle {
                center: a,
                radius_mm,
            },
            _ => return Err("This tool does not place points on face corners".into()),
        };
        let drawn = multicad_solver::draw::prepare_draw_on(
            session,
            session.revision(),
            Actor::local(channel),
            sketch,
            &base,
            Some(&report),
            &DrawRequest {
                input,
                construction: self.construction,
            },
        )
        .map_err(|e| e.to_string())?;
        let prepared = drawn.edit.preview.ok_or_else(|| {
            format!(
                "Sketch rejected: {:?}; nothing was placed",
                drawn.edit.evaluation.diagnostics.state
            )
        })?;
        self.projection = Some(Pending {
            prepared,
            requested: false,
            validated: None,
            created: drawn.created_curves.first().copied(),
        });
        self.clear_drawing_progress();
        self.message = "Projecting face corner · checking the native reference…".into();
        Ok(())
    }
    pub fn projection_pending(&self) -> bool {
        self.projection.is_some()
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
        match result.and_then(|()| {
            session
                .commit(pending.prepared)
                .map(|_| ())
                .map_err(|e| e.to_string())
        }) {
            Ok(()) => {
                self.element = pending.created;
                self.revision = session.revision();
                self.message =
                    "Face corner projected and geometry committed · one Undo step".into();
            }
            Err(error) => {
                self.message = format!("{error}; nothing was placed");
            }
        }
        true
    }
}
