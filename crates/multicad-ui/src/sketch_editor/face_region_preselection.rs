//! Selection-first Model interaction for a finished sketch: on a finite face,
//! or on a datum plane (CAD-400).
use super::{Document, EntityData, EntityId, Pos2, Rect, Session, SketchEditor};
use crate::viewport::Scene;
use multicad_core::sketch_regions::supporting_face;

impl SketchEditor {
    pub(super) fn profile_pointer_allowed(&self, document: &Document, sketch: EntityId) -> bool {
        let developed = super::developed_sheet::chart(document, sketch).is_some();
        (self.editing
            || self.acquiring_profiles()
            || (supporting_face(document, sketch).is_some()
                || !developed
                || self.developed_viewport.is_some() && developed)
                && self.overlay_visible(sketch))
            && self.tool.is_none()
            && !self.moving
            && !self.picking
            && !self.measurement.active()
            && (!self.has_pending_edits() || self.acquiring_profiles())
    }

    pub(super) fn profile_fill_point(
        &self,
        session: &Session,
        rect: Rect,
        pointer: Pos2,
    ) -> Option<(EntityId, super::Frame, [f64; 2])> {
        let sketch = self.selected?;
        let presentation = self.sketch_evaluation.current(sketch)?;
        if !presentation
            .evaluated
            .is_current(session.document(), session.revision())
        {
            return None;
        }
        let frame = self.display_frame(sketch, presentation.evaluated.frame);
        let point = self.camera.on_plane(pointer, rect, frame)?;
        // Cached display fills only decide which input lane owns this click.
        // The sketch worker classifies the exact region and its boundaries.
        presentation
            .fills
            .iter()
            .any(|fill| fill.contains(point))
            .then_some((sketch, frame, point))
    }

    pub(super) fn face_region_point(
        &self,
        session: &Session,
        scene: &Scene,
        rect: Rect,
        pointer: Pos2,
    ) -> Option<(EntityId, [f64; 2])> {
        if scene.stamp().document != session.document().id()
            || scene.stamp().revision != session.revision()
        {
            return None;
        }
        let (sketch, frame, point) = self.profile_fill_point(session, rect, pointer)?;
        // Do not select a profile through foreground geometry or through an
        // opening in the supporting face.
        let ray = self.camera.ray(pointer, rect)?;
        let hit = scene.pick_face_section(ray, self.modeling_section());
        let point_world = frame.point(point);
        let distance: f64 = point_world
            .into_iter()
            .zip(ray.origin())
            .zip(ray.direction())
            .map(|((point, origin), axis)| (point - origin) * axis)
            .sum();
        if supporting_face(session.document(), sketch).is_none() {
            // A datum-plane sketch has no supporting face; only geometry in
            // front of its plane hides the area.
            return hit
                .is_none_or(|hit| {
                    hit.distance_mm >= distance - super::face_profiles::COPLANAR_DISTANCE_MM
                })
                .then_some((sketch, point));
        }
        ((hit?.distance_mm - distance).abs() <= super::face_profiles::COPLANAR_DISTANCE_MM)
            .then_some((sketch, point))
    }
    /// Tool-first Extrude or Revolve with no sketch chosen yet (CAD-400): the
    /// nearest visible finished sketch whose closed area lies under the
    /// pointer and is not hidden behind model geometry.
    pub(super) fn sketch_area_under_pointer(
        &self,
        session: &Session,
        scene: Option<&Scene>,
        rect: Rect,
        pointer: Pos2,
    ) -> Option<(EntityId, [f64; 2])> {
        let document = session.document();
        let ray = self.camera.ray(pointer, rect)?;
        let hit = scene
            .filter(|scene| {
                scene.stamp().document == document.id()
                    && scene.stamp().revision == session.revision()
            })
            .and_then(|scene| scene.pick_face_section(ray, self.modeling_section()));
        let mut best: Option<(f64, EntityId, [f64; 2])> = None;
        for entity in document.entities() {
            let id = entity.id;
            if !matches!(entity.data, EntityData::Sketch { .. })
                || super::developed_sheet::chart(document, id).is_some()
                || !self.overlay_visible(id)
                || !self.tree.show_sources() && self.solid_source_sketches.contains(&id)
            {
                continue;
            }
            let Some(presentation) = self
                .sketch_evaluation
                .current(id)
                .filter(|p| p.evaluated.is_current(document, session.revision()))
            else {
                continue;
            };
            let frame = self.display_frame(id, presentation.evaluated.frame);
            let Some(point) = self.camera.on_plane(pointer, rect, frame) else {
                continue;
            };
            if !presentation.fills.iter().any(|fill| fill.contains(point)) {
                continue;
            }
            let distance: f64 = frame
                .point(point)
                .into_iter()
                .zip(ray.origin())
                .zip(ray.direction())
                .map(|((point, origin), axis)| (point - origin) * axis)
                .sum();
            let tolerance = super::face_profiles::COPLANAR_DISTANCE_MM;
            let visible = if supporting_face(document, id).is_some() {
                hit.as_ref()
                    .is_some_and(|hit| (hit.distance_mm - distance).abs() <= tolerance)
            } else {
                hit.as_ref()
                    .is_none_or(|hit| hit.distance_mm >= distance - tolerance)
            };
            if visible && best.is_none_or(|(nearest, _, _)| distance < nearest) {
                best = Some((distance, id, point));
            }
        }
        best.map(|(_, id, point)| (id, point))
    }
    /// Why a sketch offers no closed area (CAD-400): its open, crossing or
    /// self-intersecting boundaries, as the region diagnostics state them.
    pub(super) fn closed_area_refusal(
        &self,
        session: &Session,
        sketch: EntityId,
    ) -> Option<String> {
        let presentation = self.sketch_evaluation.current(sketch).filter(|p| {
            p.evaluated
                .is_current(session.document(), session.revision())
        })?;
        if !presentation.regions.regions.is_empty() {
            return None;
        }
        let mut reasons: Vec<&str> = vec![];
        for diagnostic in &presentation.regions.diagnostics {
            if !reasons.contains(&diagnostic.message.as_str()) {
                reasons.push(&diagnostic.message);
            }
        }
        let name = session
            .document()
            .get(sketch)
            .map_or("The sketch", |entity| entity.name.as_str());
        Some(if reasons.is_empty() {
            format!("{name} has no closed area to use as a profile")
        } else {
            format!(
                "{name} has no closed area to use as a profile: {}",
                reasons.join("; ")
            )
        })
    }
    /// Refusal for the visible finished sketch nearest the pointer's plane
    /// point when a tool-first profile click finds no closed area.
    pub(super) fn refusal_under_pointer(
        &self,
        session: &Session,
        rect: Rect,
        pointer: Pos2,
    ) -> Option<String> {
        let document = session.document();
        document
            .entities()
            .filter(|entity| {
                matches!(entity.data, EntityData::Sketch { .. })
                    && super::developed_sheet::chart(document, entity.id).is_none()
                    && self.overlay_visible(entity.id)
                    && (self.tree.show_sources()
                        || !self.solid_source_sketches.contains(&entity.id))
            })
            .find_map(|entity| {
                let presentation = self.sketch_evaluation.current(entity.id)?;
                let frame = self.display_frame(entity.id, presentation.evaluated.frame);
                self.camera.on_plane(pointer, rect, frame)?;
                self.closed_area_refusal(session, entity.id)
            })
    }
}
