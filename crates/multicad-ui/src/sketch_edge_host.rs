//! Edge acquisition uses the current native upstream producer, including saved edits.
use super::*;
use crate::viewport::{PickFilter, PickOptions, Scene, SceneFeed};
use multicad_core::{
    model::BodyInput,
    project_file::{ContentHash, ProjectResources},
    references::{self, Resolution, TopologyKind},
};
use std::{collections::BTreeSet, sync::Arc};

impl SketchEditor {
    pub fn begin_edge_treatment(
        &mut self,
        session: &Session,
        resources: &ProjectResources,
        feed: &SceneFeed,
        channel: Channel,
    ) -> Result<(), String> {
        if !self.primitive_ready() || self.editing || self.tree_selection.len() != 1 {
            return Err("Select one current native Part definition and finish pending edits before treating its edges".into());
        }
        let body = *self.tree_selection.first().expect("one selected Part");
        if !matches!(
            session.document().get(body).map(|e| &e.data),
            Some(EntityData::Part { .. })
        ) {
            return Err("Fillet / Chamfer edits a Part definition and its linked occurrences. Select the Part rather than a placed Instance".into());
        }
        let scene = feed
            .visible()
            .filter(|_| feed.is_current(session.document().id(), session.revision()))
            .ok_or("Wait for the current native model")?;
        // Only an edge hit seeds Fillet / Chamfer; a face set (CAD-152) is
        // deliberately ignored here.
        let edges = self
            .solid
            .selected_hit()
            .filter(|hit| hit.instance == body && hit.key.kind == TopologyKind::Edge)
            .and_then(|hit| hit.intent.clone())
            .into_iter()
            .collect();
        self.finish()?;
        self.model
            .begin_edge_treatment(session, resources, scene, body, edges, channel)?;
        self.message="Fillet / Chamfer · select native Part edges; Shift adds or removes. Linked occurrences follow this definition".into();
        Ok(())
    }
    pub(super) fn edge_source_scene(
        &self,
        session: &Session,
        feed: &SceneFeed,
    ) -> Result<Arc<Scene>, String> {
        let input = self
            .model
            .edge_treatment_input()
            .ok_or("Start Fillet / Chamfer first")?;
        let scene = feed
            .visible()
            .filter(|_| feed.is_current(session.document().id(), session.revision()))
            .ok_or("Wait for current native geometry before selecting edges")?;
        if !scene.is_full_model()
            || scene.model_evaluation().is_none_or(|info| {
                ContentHash::document(session.document()).ok() != Some(info.source_hash)
            })
        {
            return Err("Edge selection requires the complete source-bound native model".into());
        }
        // CAD-409: every displayed body stays drawn and pickable; only the
        // edited body is shown from its frozen upstream producer.
        scene.source_in_context(input.body, input.producer)
    }
    pub(super) fn edge_viewport(
        &mut self,
        ui: &egui::Ui,
        response: &egui::Response,
        session: &Session,
        resources: &ProjectResources,
        feed: &SceneFeed,
        rect: Rect,
    ) {
        if !self.model.edge_treatment_picking() {
            return;
        }
        let Ok(scene) = self.edge_source_scene(session, feed) else {
            return;
        };
        let Some(input) = self.model.edge_treatment_input().cloned() else {
            return;
        };
        let snapshot = scene.geometry_snapshot(input.producer);
        let kind = self
            .model
            .edge_treatment_pick_kind()
            .unwrap_or(TopologyKind::Edge);
        let hovered = response.hover_pos().and_then(|pointer| {
            scene.pick_screen_with(
                self.camera,
                rect,
                pointer,
                PickOptions {
                    filter: if kind == TopologyKind::Face {
                        PickFilter::Faces
                    } else {
                        PickFilter::Edges
                    },
                    include_hidden: false,
                    section: self.solid.section.plane(),
                },
            )
        });
        // CAD-409: before the first edge, a pick on another body (edge or
        // face) switches the edited body; afterwards those picks are ignored.
        let other_body = response.clicked().then(|| {
            hovered.clone().or_else(|| {
                let pointer = response.interact_pointer_pos()?;
                scene.pick_screen_with(
                    self.camera,
                    rect,
                    pointer,
                    PickOptions {
                        filter: PickFilter::Faces,
                        include_hidden: false,
                        section: self.solid.section.plane(),
                    },
                )
            })
        });
        if let Some(Some(hit)) = other_body.filter(|hit| {
            hit.as_ref()
                .is_some_and(|hit| hit.instance != input.producer)
        }) {
            self.retarget_edge_treatment(session, resources, feed, &hit, kind);
            return;
        }
        if response.clicked() {
            if let Some(hit) = &hovered {
                if let Some(capture) = &producer_capture(&scene, &input, hit) {
                    let current = self
                        .model
                        .edge_treatment_pick_captures()
                        .unwrap_or_default();
                    let same = current.iter().position(|intent| {
                        matches!(snapshot.as_ref().map(|snapshot| references::resolve(intent, snapshot, session.revision())), Some(Resolution::Resolved {key, ..}) if key == hit.key)
                    });
                    let modifiers = crate::selection_input::click_modifiers(ui, response);
                    let edges = if kind == TopologyKind::Edge
                        && (modifiers.shift || modifiers.command)
                    {
                        let mut edges = current.to_vec();
                        if let Some(index) = same {
                            edges.remove(index);
                        } else {
                            edges.push(capture.clone());
                        }
                        edges
                    } else {
                        vec![same.map_or_else(|| capture.clone(), |index| current[index].clone())]
                    };
                    if let Err(error) = self.model.set_edge_treatment_pick_captures(
                        session,
                        resources,
                        feed.visible().expect("current canonical scene"),
                        edges,
                    ) {
                        self.message = error;
                    }
                } else {
                    self.message="This native element has no persistent geometry name; select another or repair its source".into();
                }
            } else {
                self.message = if kind == TopologyKind::Face {
                    "Choose the adjacent side face for this chamfer edge"
                } else {
                    "Pick an edge of a body"
                }
                .into();
            }
        }
        if kind == TopologyKind::Face {
            let highlight = hovered.as_ref().map(|hit| hit.key).or_else(|| {
                self.model
                    .edge_treatment_pick_captures()?
                    .last()
                    .and_then(|capture| {
                        match references::resolve(capture, snapshot.as_ref()?, session.revision()) {
                            Resolution::Resolved { key, .. } => Some(key),
                            _ => None,
                        }
                    })
            });
            self.solid
                .highlight_reference(highlight.map(|key| crate::viewport::Selection::Face {
                    instance: input.producer,
                    ordinal: key.ordinal,
                }));
            return;
        }
        let selected = self
            .model
            .edge_treatment_edges()
            .unwrap_or_default()
            .iter()
            .filter_map(|capture| {
                match references::resolve(capture, snapshot.as_ref()?, session.revision()) {
                    Resolution::Resolved { key, .. } => Some(key.ordinal),
                    _ => None,
                }
            })
            .collect::<BTreeSet<_>>();
        let hover = hovered.as_ref().map(|hit| hit.key.ordinal);
        EdgeView {
            camera: self.camera,
            rect,
            section: self.solid.section.plane(),
            painter: &ui.painter().with_clip_rect(rect),
        }
        .paint(&scene, &input, &selected, hover);
    }
}
impl SketchEditor {
    /// An edge or face picked on another displayed body (CAD-409): before
    /// the first edge a new treatment moves there; a side-face pick must
    /// stay on the edited body.
    fn retarget_edge_treatment(
        &mut self,
        session: &Session,
        resources: &ProjectResources,
        feed: &SceneFeed,
        hit: &crate::viewport::Hit,
        kind: TopologyKind,
    ) {
        if kind == TopologyKind::Face {
            self.message = "Choose a side face of the body whose edge you picked".into();
            return;
        }
        let Some(scene) = feed.visible() else {
            return;
        };
        let edges = if hit.key.kind == TopologyKind::Edge {
            let Some(capture) = scene
                .geometry_snapshot(hit.geometry)
                .and_then(|snapshot| references::GeometryIntent::capture(&snapshot, hit.key).ok())
            else {
                self.message =
                    "This edge has no persistent geometry name; pick another edge".into();
                return;
            };
            vec![capture]
        } else {
            vec![]
        };
        self.message =
            match self
                .model
                .retarget_edge_treatment(session, resources, scene, hit.instance, edges)
            {
                Ok(name) => format!("Fillet / Chamfer · {name}"),
                Err(error) => error,
            };
    }
}
/// A pick on the edited body, named against its frozen upstream producer.
fn producer_capture(
    scene: &Scene,
    input: &BodyInput,
    hit: &crate::viewport::Hit,
) -> Option<references::GeometryIntent> {
    let mut snapshot = scene.geometry_snapshot(hit.geometry)?;
    snapshot.owner = input.producer;
    references::GeometryIntent::capture(&snapshot, hit.key).ok()
}
struct EdgeView<'a> {
    camera: Camera,
    rect: Rect,
    section: Option<crate::viewport::SectionPlane>,
    painter: &'a egui::Painter,
}
impl EdgeView<'_> {
    fn paint(
        &self,
        scene: &Scene,
        input: &BodyInput,
        selected: &BTreeSet<u32>,
        hover: Option<u32>,
    ) {
        let Some(geometry) = scene.geometries().get(&input.producer) else {
            return;
        };
        let transform = scene
            .instances()
            .iter()
            .find(|instance| instance.id == input.producer)
            .map(|instance| instance.transform)
            .unwrap_or_default();
        for edge in geometry
            .edges()
            .iter()
            .filter(|edge| selected.contains(&edge.index) || hover == Some(edge.index))
        {
            let color = if selected.contains(&edge.index) {
                Color32::from_rgb(80, 220, 160)
            } else {
                Color32::YELLOW
            };
            for pair in edge.points_mm.windows(2) {
                let pair = [transform.point(pair[0]), transform.point(pair[1])];
                let segment = if let Some(plane) = self.section {
                    plane.clip_segment(pair[0], pair[1])
                } else {
                    Some([pair[0], pair[1]])
                };
                if let Some([a, b]) = segment
                    && let Some(segment) = self.camera.segment(a, b, self.rect)
                {
                    self.painter.line_segment(segment, Stroke::new(2.5, color));
                }
            }
        }
    }
}
