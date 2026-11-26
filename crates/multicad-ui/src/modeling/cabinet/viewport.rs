//! Live viewport draft of a rectangular panel (CAD-374): the finished
//! envelope with corner resize handles and dimension inputs that edit the
//! same staged Inspector fields. Thickness belongs to the panel stock, so the
//! handles never change it and the inputs say where it comes from.
use super::{CabinetCommand, CabinetPanel, Kind, panel_stocks};
use crate::{
    egui,
    modeling::{
        dimension::DimensionField,
        primitive::{
            PrimitiveViewportResult,
            corner_inputs::{CornerEdit, CornerField},
            gizmo,
        },
    },
};
use multicad_core::{
    commands::Session,
    model::{self, PanelShape, PanelSource, ResolvedPrimitive},
    sketch::Frame,
};

/// The panel thickness axis (N) is stock-controlled.
pub(super) const LOCKED: [bool; 3] = [false, false, true];

impl CabinetPanel {
    /// The resolved finished envelope of a rectangular panel draft, as a box
    /// along the panel's (u, v, thickness normal) frame.
    pub(in crate::modeling) fn panel_shape(&self, session: &Session) -> Option<ResolvedPrimitive> {
        if self.kind != Kind::Panel {
            return None;
        }
        let Ok(CabinetCommand::Panel { panel, .. }) = self.request(session.document(), true) else {
            return None;
        };
        if !matches!(panel.source, PanelSource::Rectangle { .. }) {
            return None;
        }
        let resolved =
            model::resolve_panel(session.document(), session.parameters(), &panel).ok()?;
        let PanelShape::Box(shape) = resolved.shape else {
            return None;
        };
        (shape.thickness_mm > 0.).then(|| ResolvedPrimitive::Box {
            size_mm: [shape.length_mm, shape.width_mm, shape.thickness_mm],
            frame: Frame {
                origin: shape.origin,
                u: shape.u,
                v: shape.v,
                n: shape.normal(),
            },
        })
    }

    pub(in crate::modeling) fn envelope(&self, session: &Session) -> Option<[[f64; 3]; 8]> {
        Some(gizmo::corners(&self.panel_shape(session)?))
    }

    pub(in crate::modeling) fn dragging(&self) -> bool {
        self.resize.active()
    }

    #[cfg(test)]
    pub(in crate::modeling) fn corner_input_active(&self) -> bool {
        self.kind == Kind::Panel && self.corner_inputs.active()
    }

    pub(in crate::modeling) fn corner_input_focused(&self) -> bool {
        self.kind == Kind::Panel && self.corner_inputs.focused()
    }

    pub(in crate::modeling) fn open_dimension_inputs(&mut self) {
        if self.kind == Kind::Panel {
            self.corner_inputs.open(7);
        }
    }

    /// A read-only explanation of the stock-controlled thickness.
    fn thickness_note(&self, session: &Session) -> Option<String> {
        let stock = self.stock?;
        panel_stocks(session.document())
            .into_iter()
            .find(|(id, ..)| *id == stock)
            .map(|(_, name, thickness)| {
                format!(
                    "Thickness {thickness} mm comes from the panel stock \u{201c}{name}\u{201d}. \
                     Choose other stock in the Inspector or edit it in Project \u{2192} Materials."
                )
            })
    }

    /// Snap panel corner drags to the displayed grid spacing (CAD-402).
    pub(in crate::modeling) fn set_drag_snap(&mut self, snap: super::super::drag_snap::DragSnap) {
        self.resize.set_snap(snap);
    }

    pub(in crate::modeling) fn viewport(
        &mut self,
        ui: &mut egui::Ui,
        session: &Session,
        camera: crate::camera::Camera,
        rect: egui::Rect,
        start_allowed: bool,
        guide: bool,
    ) -> PrimitiveViewportResult {
        if self.kind != Kind::Panel {
            return PrimitiveViewportResult::default();
        }
        let current = self.panel_shape(session);
        let valid = current.is_some();
        if valid {
            self.last_shape = current.clone();
        } else {
            let snap = self.resize.snap();
            self.resize = gizmo::ResizeGizmo::locked(LOCKED);
            self.resize.set_snap(snap);
        }
        let Some(shape) = current.or_else(|| self.last_shape.clone()) else {
            return PrimitiveViewportResult::default();
        };
        let over_inputs = ui.input(|input| {
            input
                .pointer
                .interact_pos()
                .is_some_and(|position| self.corner_inputs.contains_pointer(position))
        });
        let gesture = self.resize.show(
            ui,
            &shape,
            camera,
            rect,
            valid && start_allowed && !over_inputs,
            guide || !valid,
        );
        if let Some(corner) = gesture.clicked_corner {
            self.corner_inputs.open(corner);
        }
        let mut result = PrimitiveViewportResult {
            input_owned: gesture.input_owned || over_inputs,
            changed: gesture
                .changed
                .is_some_and(|next| self.resize_draft(&shape, &next)),
            released: gesture.released,
            confirm: false,
        };
        if self.corner_inputs.active() {
            let anchor = if result.changed {
                self.panel_shape(session).unwrap_or(shape)
            } else {
                shape
            };
            let notes: Vec<String> = self.thickness_note(session).into_iter().collect();
            let mut fields = [
                CornerField {
                    label: "Length (u)",
                    value: &mut self.length,
                },
                CornerField {
                    label: "Width (v)",
                    value: &mut self.panel_width,
                },
            ];
            let CornerEdit { changed, confirm } = self.corner_inputs.show(
                ui,
                self.body,
                rect,
                camera,
                gizmo::corners(&anchor),
                session,
                &mut fields,
                &notes,
            );
            result.changed |= changed;
            result.confirm = confirm && self.panel_shape(session).is_some();
        }
        if result.changed {
            self.stamp.changed();
            self.error = None;
            self.guide = None;
        }
        self.paint_alignment(session, &ui.painter_at(rect), camera, rect);
        result
    }

    /// Write a dragged envelope into the staged origin, length and width.
    fn resize_draft(&mut self, before: &ResolvedPrimitive, after: &ResolvedPrimitive) -> bool {
        let (
            ResolvedPrimitive::Box {
                size_mm: old_size,
                frame: old_frame,
            },
            ResolvedPrimitive::Box {
                size_mm: new_size,
                frame: new_frame,
            },
        ) = (before, after)
        else {
            return false;
        };
        let mut changed = false;
        for (index, field) in self.origin.iter_mut().enumerate() {
            changed |= set_length(field, old_frame.origin[index], new_frame.origin[index]);
        }
        changed |= set_length(&mut self.length, old_size[0], new_size[0]);
        changed |= set_length(&mut self.panel_width, old_size[1], new_size[1]);
        changed
    }
}

fn set_length(field: &mut DimensionField, before: f64, after: f64) -> bool {
    if (after - before).abs() <= 1e-8 {
        return false;
    }
    field.set_exact_literal(after, "mm");
    true
}

#[cfg(test)]
mod tests {
    use super::super::tests::session;
    use super::*;
    use crate::modeling::PreviewState;
    use multicad_core::{
        document::{Entity, EntityData},
        model::PanelPlane,
    };

    const RECT: egui::Rect = egui::Rect {
        min: egui::Pos2::ZERO,
        max: egui::pos2(1024., 700.),
    };

    fn top_view() -> crate::camera::Camera {
        crate::camera::Camera {
            yaw: 0.,
            pitch: std::f32::consts::FRAC_PI_2,
            scale: 0.8,
            target: [400., 200., 0.],
            perspective: false,
        }
    }

    fn project(camera: crate::camera::Camera, point: [f64; 3]) -> egui::Pos2 {
        let [x, y, _] = camera.project(point.map(|value| value as f32), RECT);
        egui::pos2(x, y)
    }

    fn button(pos: egui::Pos2, pressed: bool) -> egui::Event {
        egui::Event::PointerButton {
            pos,
            pressed,
            button: egui::PointerButton::Primary,
            modifiers: egui::Modifiers::NONE,
        }
    }

    fn frame(
        form: &mut CabinetPanel,
        session: &Session,
        context: &egui::Context,
        camera: crate::camera::Camera,
        events: Vec<egui::Event>,
        start_allowed: bool,
    ) -> PrimitiveViewportResult {
        let mut result = None;
        context
            .run_ui(
                egui::RawInput {
                    screen_rect: Some(RECT),
                    events,
                    focused: true,
                    ..Default::default()
                },
                |ui| {
                    result = Some(form.viewport(ui, session, camera, RECT, start_allowed, true));
                },
            )
            .drop_without_applying_deltas();
        result.unwrap()
    }

    fn size(form: &CabinetPanel, session: &Session) -> ([f64; 3], [f64; 3]) {
        let Some(ResolvedPrimitive::Box { size_mm, frame }) = form.panel_shape(session) else {
            panic!("a valid panel draft resolves to a box envelope");
        };
        (size_mm, frame.origin)
    }

    fn near(actual: [f64; 3], expected: [f64; 3]) {
        for (a, e) in actual.into_iter().zip(expected) {
            assert!((a - e).abs() < 1e-4, "{actual:?} != {expected:?}");
        }
    }

    #[test]
    fn panel_draft_envelope_follows_plane_and_stock_thickness_without_preview() {
        let s = session();
        let mut form = CabinetPanel::panel(s.document(), "Panel 1".into());
        let (size_mm, origin) = size(&form, &s);
        near(size_mm, [800., 400., 18.]);
        near(origin, [0.; 3]);
        // XZ panels grow their thickness toward −Y, exactly as the kernel does.
        form.plane = PanelPlane::Xz;
        let corners = form.envelope(&s).unwrap();
        near(corners[7], [800., -18., 400.]);
        // Other stock changes only the thickness.
        let thin = panel_stocks(s.document())
            .into_iter()
            .find(|(_, name, _)| name == "Ply 6")
            .unwrap()
            .0;
        form.select_stock(thin, s.document());
        near(size(&form, &s).0, [800., 400., 6.]);
        assert!(form.thickness_note(&s).unwrap().contains("6 mm"));
        assert!(form.thickness_note(&s).unwrap().contains("Ply 6"));
        // Carcass drafts and drafts without stock show no panel envelope.
        form.stock = None;
        assert!(form.panel_shape(&s).is_none());
        let carcass = CabinetPanel::new(s.document(), "Carcass".into());
        assert!(carcass.envelope(&s).is_none());
    }

    #[test]
    fn panel_corner_drag_edits_length_and_width_but_never_stock_thickness() {
        let s = session();
        let before = s.document().clone();
        let mut form = CabinetPanel::panel(s.document(), "Panel 1".into());
        let context = egui::Context::default();
        let camera = top_view();
        let corner = project(camera, [800., 400., 18.]);
        let target = project(camera, [900., 450., 18.]);
        let request = form.request(s.document(), true).unwrap();
        let stamp = form.stamp();
        for allowed in [false, true] {
            frame(&mut form, &s, &context, camera, vec![], allowed);
            let press = frame(
                &mut form,
                &s,
                &context,
                camera,
                vec![egui::Event::PointerMoved(corner), button(corner, true)],
                allowed,
            );
            assert_eq!(press.input_owned, allowed);
            let moved = frame(
                &mut form,
                &s,
                &context,
                camera,
                vec![egui::Event::PointerMoved(target)],
                allowed,
            );
            assert_eq!(moved.changed, allowed);
            let release = frame(
                &mut form,
                &s,
                &context,
                camera,
                vec![button(target, false)],
                allowed,
            );
            assert_eq!(release.released, allowed);
            if !allowed {
                assert_eq!(
                    format!("{:?}", form.request(s.document(), true).unwrap()),
                    format!("{request:?}")
                );
                assert_eq!(form.stamp(), stamp);
            }
        }
        assert_ne!(form.stamp(), stamp, "a drag stages a new draft revision");
        let (size_mm, origin) = size(&form, &s);
        near(size_mm, [900., 450., 18.]);
        near(origin, [0.; 3]);
        // The Inspector fields are the staged values the drag changed.
        assert_eq!(form.length.text, "900");
        assert_eq!(form.panel_width.text, "450");

        // The low corner moves the origin in the panel plane only.
        let low = project(camera, [0., 0., 18.]);
        let shifted = project(camera, [50., 30., 18.]);
        frame(
            &mut form,
            &s,
            &context,
            camera,
            vec![egui::Event::PointerMoved(low), button(low, true)],
            true,
        );
        frame(
            &mut form,
            &s,
            &context,
            camera,
            vec![egui::Event::PointerMoved(shifted)],
            true,
        );
        frame(
            &mut form,
            &s,
            &context,
            camera,
            vec![button(shifted, false)],
            true,
        );
        let (size_mm, origin) = size(&form, &s);
        near(size_mm, [850., 420., 18.]);
        near(origin, [50., 30., 0.]);

        // A side view along V would make thickness the drag axis; it stays.
        let side = crate::camera::Camera {
            yaw: 0.,
            pitch: 0.,
            ..camera
        };
        let edge = project(side, [900., 450., 18.]);
        frame(&mut form, &s, &context, side, vec![], true);
        frame(
            &mut form,
            &s,
            &context,
            side,
            vec![egui::Event::PointerMoved(edge), button(edge, true)],
            true,
        );
        frame(
            &mut form,
            &s,
            &context,
            side,
            vec![egui::Event::PointerMoved(edge + egui::vec2(40., -40.))],
            true,
        );
        frame(
            &mut form,
            &s,
            &context,
            side,
            vec![button(edge + egui::vec2(40., -40.), false)],
            true,
        );
        assert_eq!(size(&form, &s).0[2], 18.);
        assert_eq!(s.document(), &before, "drafting never edits the document");
    }

    #[test]
    fn panel_viewport_inputs_tab_between_length_and_width_and_confirm_once() {
        let mut s = session();
        let mut document = s.document().clone();
        document
            .insert(Entity::new(
                "shelf_width",
                EntityData::Parameter {
                    expression: "380".into(),
                    unit: "mm".into(),
                    exposed: true,
                    inputs: vec![],
                },
            ))
            .unwrap();
        s = Session::new(document.clone()).unwrap();
        let mut form = CabinetPanel::panel(s.document(), "Panel 1".into());
        form.open_dimension_inputs();
        assert!(form.corner_input_active());
        let context = egui::Context::default();
        context.enable_accesskit();
        let camera = top_view();
        let focus = || context.memory(|memory| memory.focused());
        let key = |form: &mut CabinetPanel, key, modifiers| {
            let mut result = None;
            for pressed in [true, false] {
                let current = frame(
                    form,
                    &s,
                    &context,
                    camera,
                    vec![
                        egui::Event::ModifiersChanged(if pressed {
                            modifiers
                        } else {
                            egui::Modifiers::NONE
                        }),
                        egui::Event::Key {
                            key,
                            physical_key: None,
                            pressed,
                            repeat: false,
                            modifiers,
                        },
                    ],
                    true,
                );
                if pressed {
                    result = Some(current);
                }
            }
            result.unwrap()
        };
        frame(&mut form, &s, &context, camera, vec![], true);
        frame(&mut form, &s, &context, camera, vec![], true);
        let length_focus = focus().expect("new panel inputs focus Length");
        frame(
            &mut form,
            &s,
            &context,
            camera,
            vec![egui::Event::Text("1200 mm".into())],
            true,
        );
        assert_eq!(form.length.text, "1200 mm");
        key(&mut form, egui::Key::Tab, egui::Modifiers::NONE);
        let width_focus = focus().expect("Tab focuses Width");
        assert_ne!(length_focus, width_focus);
        // An unresolved formula never confirms and leaves Length intact.
        frame(
            &mut form,
            &s,
            &context,
            camera,
            vec![egui::Event::Text("=$missing".into())],
            true,
        );
        assert!(!key(&mut form, egui::Key::Enter, egui::Modifiers::NONE).confirm);
        assert_eq!(form.length.text, "1200 mm");
        key(&mut form, egui::Key::Tab, egui::Modifiers::SHIFT);
        assert_eq!(focus(), Some(length_focus));
        key(&mut form, egui::Key::Tab, egui::Modifiers::NONE);
        assert_eq!(focus(), Some(width_focus));
        frame(
            &mut form,
            &s,
            &context,
            camera,
            vec![egui::Event::Text("=$shelf_width + 2".into())],
            true,
        );
        assert_eq!(form.panel_width.text, "=$shelf_width + 2");
        near(size(&form, &s).0, [1200., 382., 18.]);
        // Thickness is explained, not editable, in the viewport inputs.
        let output = context.run_ui(
            egui::RawInput {
                screen_rect: Some(RECT),
                ..Default::default()
            },
            |ui| {
                form.viewport(ui, &s, camera, RECT, true, true);
            },
        );
        let nodes = &output
            .platform_output
            .accesskit_update
            .as_ref()
            .unwrap()
            .nodes;
        assert!(nodes.iter().any(|(_, node)| {
            node.value()
                .or(node.label())
                .is_some_and(|text| text.contains("Thickness 18 mm comes from the panel stock"))
        }));
        output.drop_without_applying_deltas();
        assert!(key(&mut form, egui::Key::Enter, egui::Modifiers::NONE).confirm);
        assert_eq!(s.document(), &document);

        // The Inspector shows the same staged values.
        let inspector = egui::Context::default();
        inspector.enable_accesskit();
        let output = inspector.run_ui(
            egui::RawInput {
                screen_rect: Some(RECT),
                ..Default::default()
            },
            |ui| {
                form.show(ui, s.document(), true, PreviewState::Idle);
            },
        );
        let nodes = &output
            .platform_output
            .accesskit_update
            .as_ref()
            .unwrap()
            .nodes;
        for text in ["1200 mm", "=$shelf_width + 2"] {
            assert!(
                nodes.iter().any(|(_, node)| node.value() == Some(text)),
                "{text}"
            );
        }
        output.drop_without_applying_deltas();
    }
}
