//! Point-to-point placement of a rectangular panel draft (CAD-375): a source
//! snap point on the draft envelope and a target snap point on other
//! geometry stage one translation of the panel origin. Plane, length, width
//! and stock thickness are never changed by the pick.
use super::{CabinetPanel, Kind};
use crate::modeling::{
    dimension::DimensionField,
    placement::points::{self, SnapPoint},
    primitive::gizmo,
};
use multicad_core::{commands::Session, model::ResolvedPrimitive};

type Point = [f64; 3];

/// Which draft snap point is the source. Indices address
/// [`CabinetPanel::source_points`], whose order is fixed, so the pick follows
/// later edits of the draft size.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(in crate::modeling) enum PanelPointPick {
    Source,
    Target { source: usize },
}

/// The last staged pair, kept while the draft still realises it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(in crate::modeling) struct PanelAlignment {
    pub source: usize,
    pub target: Point,
    pub target_label: &'static str,
}

/// Envelope edges as corner-index pairs (bits 0, 1, 2 = high U, V, N).
const EDGES: [(usize, usize); 12] = [
    (0, 1),
    (2, 3),
    (4, 5),
    (6, 7),
    (0, 2),
    (1, 3),
    (4, 6),
    (5, 7),
    (0, 4),
    (1, 5),
    (2, 6),
    (3, 7),
];

fn source_points(shape: &ResolvedPrimitive) -> Vec<SnapPoint> {
    let corners = gizmo::corners(shape);
    let point = |canonical: Point, label| SnapPoint {
        owner: None,
        canonical,
        shown: canonical,
        label,
    };
    corners
        .iter()
        .map(|corner| point(*corner, "Panel corner"))
        .chain(EDGES.iter().map(|(a, b)| {
            point(
                std::array::from_fn(|i| (corners[*a][i] + corners[*b][i]) * 0.5),
                "Panel edge midpoint",
            )
        }))
        .collect()
}

impl CabinetPanel {
    pub(in crate::modeling) fn point_pick(&self) -> Option<PanelPointPick> {
        self.point_pick
    }

    /// Supported source snap points on the current draft envelope.
    pub(in crate::modeling) fn source_points(&self, session: &Session) -> Vec<SnapPoint> {
        self.panel_shape(session)
            .map_or_else(Vec::new, |shape| source_points(&shape))
    }

    /// The source point of a pick or alignment on the current draft.
    pub(in crate::modeling) fn source_point(
        &self,
        session: &Session,
        index: usize,
    ) -> Option<SnapPoint> {
        self.source_points(session).get(index).copied()
    }

    #[cfg(test)]
    pub(in crate::modeling) fn alignment(&self) -> Option<PanelAlignment> {
        self.alignment
    }

    pub(in crate::modeling) fn toggle_point_pick(&mut self) {
        self.point_pick = match self.point_pick {
            None if self.kind == Kind::Panel => Some(PanelPointPick::Source),
            _ => None,
        };
    }

    pub(in crate::modeling) fn pick_source(
        &mut self,
        session: &Session,
        point: Point,
    ) -> Result<&'static str, String> {
        if self.point_pick != Some(PanelPointPick::Source) {
            return Err("Choose Points in the Inspector first".into());
        }
        let points = self.source_points(session);
        let (index, picked) = points
            .iter()
            .enumerate()
            .find(|(_, candidate)| candidate.canonical == point)
            .ok_or("The chosen source is no longer on the panel draft")?;
        self.point_pick = Some(PanelPointPick::Target { source: index });
        Ok(picked.label)
    }

    /// Translate the draft so its source point coincides with `target`.
    /// A failure leaves the previous staged origin and pick state intact.
    pub(in crate::modeling) fn pick_target(
        &mut self,
        session: &Session,
        target: Point,
        target_label: &'static str,
    ) -> Result<(), String> {
        let Some(PanelPointPick::Target { source }) = self.point_pick else {
            return Err("Choose a source point on the new panel first".into());
        };
        if target.iter().any(|value| !value.is_finite()) {
            return Err("The chosen target point is not finite".into());
        }
        let shape = self
            .panel_shape(session)
            .ok_or("Fix the panel inputs before placing it by points")?;
        let ResolvedPrimitive::Box { frame, .. } = &shape else {
            unreachable!("panel drafts resolve as boxes")
        };
        let point = source_points(&shape)
            .get(source)
            .ok_or("The chosen source is no longer on the panel draft")?
            .canonical;
        let previous = self.origin.clone();
        for (axis, field) in self.origin.iter_mut().enumerate() {
            let next = frame.origin[axis] + target[axis] - point[axis];
            set_exact(field, frame.origin[axis], next);
        }
        let placed = self.source_point(session, source).map(|p| p.canonical);
        if placed.is_none_or(|placed| (0..3).any(|i| (placed[i] - target[i]).abs() > 1e-6)) {
            self.origin = previous;
            return Err("The point-to-point placement could not be represented exactly".into());
        }
        self.point_pick = None;
        self.alignment = Some(PanelAlignment {
            source,
            target,
            target_label,
        });
        self.stamp.changed();
        self.error = None;
        self.guide = None;
        Ok(())
    }

    /// Paint the staged pair while the draft still places its source on the
    /// target; any later edit that moves it away clears the stale pair.
    pub(in crate::modeling) fn paint_alignment(
        &mut self,
        session: &Session,
        painter: &crate::egui::Painter,
        camera: crate::camera::Camera,
        rect: crate::egui::Rect,
    ) {
        let Some(alignment) = self.alignment else {
            return;
        };
        let Some(source) = self.source_point(session, alignment.source) else {
            return;
        };
        if (0..3).any(|i| (source.canonical[i] - alignment.target[i]).abs() > 1e-6) {
            self.alignment = None;
            return;
        }
        points::paint_markers(
            painter,
            camera,
            rect,
            Some((source.canonical, source.label)),
            Some((alignment.target, alignment.target_label)),
        );
    }
}

#[cfg(test)]
impl CabinetPanel {
    pub(in crate::modeling) fn origin_texts(&self) -> Vec<String> {
        self.origin.iter().map(|field| field.text.clone()).collect()
    }
}

fn set_exact(field: &mut DimensionField, before: f64, after: f64) {
    if (after - before).abs() > 1e-12 {
        field.set_exact_literal(after, "mm");
    }
}

#[cfg(test)]
mod tests {
    use super::super::tests::session;
    use super::*;

    #[test]
    fn source_points_cover_corners_and_edge_midpoints_in_a_stable_order() {
        let s = session();
        let form = CabinetPanel::panel(s.document(), "Panel 1".into());
        let points = form.source_points(&s);
        assert_eq!(points.len(), 20);
        assert_eq!(points[0].canonical, [0.; 3]);
        assert_eq!(points[7].canonical, [800., 400., 18.]);
        assert_eq!(points[8].canonical, [400., 0., 0.]);
        assert_eq!(points[8].label, "Panel edge midpoint");
        assert!(points.iter().all(|p| p.owner.is_none()));
    }

    #[test]
    fn picking_source_and_target_moves_only_the_origin() {
        let s = session();
        let mut form = CabinetPanel::panel(s.document(), "Panel 1".into());
        let before = (
            form.length.text.clone(),
            form.panel_width.text.clone(),
            form.plane,
        );
        assert!(
            form.pick_source(&s, [800., 400., 18.]).is_err(),
            "not picking yet"
        );
        form.toggle_point_pick();
        assert!(form.pick_target(&s, [1., 2., 3.], "Vertex").is_err());
        assert!(
            form.pick_source(&s, [1., 2., 3.]).is_err(),
            "not on the draft"
        );
        assert_eq!(
            form.pick_source(&s, [800., 0., 18.]).unwrap(),
            "Panel corner"
        );
        let stamp = form.stamp();
        form.pick_target(&s, [1000., 250., 720.], "Vertex").unwrap();
        assert_ne!(form.stamp(), stamp);
        assert_eq!(form.point_pick(), None);
        let source = form.source_point(&s, 5).unwrap().canonical;
        assert_eq!(source, [1000., 250., 720.]);
        assert_eq!(form.origin[0].text, "200");
        assert_eq!(form.origin[1].text, "250");
        assert_eq!(form.origin[2].text, "702");
        assert_eq!(
            (
                form.length.text.clone(),
                form.panel_width.text.clone(),
                form.plane
            ),
            before
        );
        assert!(form.alignment().is_some());
    }

    #[test]
    fn a_failed_target_keeps_the_staged_origin_and_pick() {
        let s = session();
        let mut form = CabinetPanel::panel(s.document(), "Panel 1".into());
        form.toggle_point_pick();
        form.pick_source(&s, [0.; 3]).unwrap();
        let origin: Vec<_> = form.origin.iter().map(|f| f.text.clone()).collect();
        assert!(form.pick_target(&s, [f64::NAN, 0., 0.], "Vertex").is_err());
        // An invalid draft cannot be placed either.
        form.length.text = "=$missing".into();
        assert!(form.pick_target(&s, [5., 5., 5.], "Vertex").is_err());
        assert_eq!(
            form.origin
                .iter()
                .map(|f| f.text.clone())
                .collect::<Vec<_>>(),
            origin
        );
        assert_eq!(
            form.point_pick(),
            Some(PanelPointPick::Target { source: 0 })
        );
        form.toggle_point_pick();
        assert_eq!(form.point_pick(), None);
    }
}
