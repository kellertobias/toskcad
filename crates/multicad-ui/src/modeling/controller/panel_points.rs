//! Viewport picking for point-to-point placement of a new or reopened
//! rectangular panel draft (CAD-375). The source is a draft snap point; the
//! target is a supported snap point on other canonical geometry, including
//! translated or rotated Instances, or the world origin.
use super::*;
use crate::modeling::{cabinet::PanelPointPick, placement::points};
use std::collections::BTreeSet;

impl ModelToolController {
    /// The panel form is waiting for a point-to-point source or target click.
    pub fn panel_point_picking(&self) -> bool {
        matches!(&self.form, Some(ModelForm::Cabinet(form)) if form.point_pick().is_some())
    }

    /// Show snap candidates and markers; a click stages the pick. Returns a
    /// status message after a successful click. Failures keep the last
    /// valid draft and pick state.
    pub fn panel_point_viewport(
        &mut self,
        session: &Session,
        resources: &ProjectResources,
        view: super::super::TerminationViewport<'_>,
    ) -> Result<Option<String>, String> {
        let Some(ModelForm::Cabinet(form)) = &self.form else {
            return Ok(None);
        };
        let Some(pick) = form.point_pick() else {
            return Ok(None);
        };
        let scene = view.scene.filter(|scene| {
            scene.is_full_model()
                && scene.stamp().document == session.document().id()
                && scene.stamp().revision == session.revision()
        });
        let (candidates, source) = match pick {
            PanelPointPick::Source => (form.source_points(session), None),
            PanelPointPick::Target { source } => {
                let Some(scene) = scene else {
                    return if view.clicked {
                        Err("Wait for the current model before choosing a target point".into())
                    } else {
                        Ok(None)
                    };
                };
                // A reopened panel never targets its own canonical geometry.
                let own = points::moving_instances(
                    session.document(),
                    scene,
                    &BTreeSet::from([form.body()]),
                );
                let mut candidates = points::snap_points(scene, |id| !own.contains(&id), |p| p);
                candidates.push(points::SnapPoint {
                    owner: None,
                    canonical: [0.; 3],
                    shown: [0.; 3],
                    label: "World origin",
                });
                let source = form
                    .source_point(session, source)
                    .map(|point| (point.shown, point.label));
                (candidates, source)
            }
        };
        let hovered = view
            .pointer
            .and_then(|pointer| points::nearest(view.camera, view.rect, pointer, &candidates))
            .filter(|(point, _)| {
                view.section
                    .is_none_or(|section| section.contains(point.shown))
            });
        let hover = hovered.map(|(point, _)| (point.shown, point.label));
        match pick {
            PanelPointPick::Source => {
                points::paint_markers(view.painter, view.camera, view.rect, hover, None)
            }
            PanelPointPick::Target { .. } => {
                points::paint_markers(view.painter, view.camera, view.rect, source, hover)
            }
        }
        if !view.clicked {
            return Ok(None);
        }
        let Some((point, _)) = hovered else {
            return Err(match pick {
                PanelPointPick::Source => "Click a highlighted corner or edge midpoint of the new panel",
                PanelPointPick::Target { .. } => "Click a highlighted vertex or edge midpoint on another part, or the world origin",
            }
            .into());
        };
        if let (PanelPointPick::Target { .. }, Some(scene)) = (pick, scene) {
            self.verify_context(session, resources, scene)?;
        }
        let Some(ModelForm::Cabinet(form)) = &mut self.form else {
            return Ok(None);
        };
        match pick {
            PanelPointPick::Source => {
                let label = form.pick_source(session, point.canonical)?;
                Ok(Some(format!(
                    "Panel source: {} · now click a target point",
                    label.to_lowercase()
                )))
            }
            PanelPointPick::Target { .. } => {
                form.pick_target(session, point.canonical, point.label)?;
                self.invalidate();
                // The form stays open if the preview provider reports an error.
                let _ = self.preview(session, resources, Channel::Ui);
                Ok(Some(format!(
                    "Panel placed on {} · Enter applies, Escape cancels",
                    point.label.to_lowercase()
                )))
            }
        }
    }
}
