//! CAD-430: position (origin corner X/Y/Z) and orientation (face direction
//! preset and turn about it) of an existing rectangle panel, shown where the
//! Part is now. Apply sends one canonical placement edit through the shared
//! core API; the viewport draft handles keep these values current until
//! the user types or picks here.
//!
//! The fields use the inspector's `DimensionField` widget; switching to the
//! shared inspector field widget (CAD-416) only changes `show`.
use super::{DimensionField, field};
use crate::{controls, egui};
use multicad_core::{
    document::{Document, EntityId},
    model::{self, PanelFacing, PanelPlacementRequest, PanelPlane, RigidPose},
    parameters,
};

pub(super) struct PlacementFields {
    pub position: [DimensionField; 3],
    pub facing: Option<PanelFacing>,
    pub rotation: DimensionField,
    /// The Part's Move / rotate after its panel feature, when opened.
    pose: RigidPose,
    /// What the draft showed last, so unedited fields follow it.
    shown: Option<([f64; 3], PanelPlane)>,
    /// The user typed or picked here: Apply places the panel.
    pub edited: bool,
}
pub(super) const POSITION: [&str; 3] = ["Position X", "Position Y", "Position Z"];
pub(super) const FACING: &str = "Face direction";
pub(super) const ROTATION: &str = "Rotation about face";

fn facing_label(facing: Option<PanelFacing>) -> &'static str {
    facing.map_or("Turned off the principal directions", PanelFacing::label)
}
impl PlacementFields {
    /// Where the panel Part is now, or None when it cannot be placed here.
    pub fn open(document: &Document, part: EntityId, feature: EntityId) -> Option<Self> {
        let p = parameters::evaluate(document).ok()?;
        let placed = model::panel_placement(document, &p, part).ok()?;
        let pose = model::panel_pose_after(document, &p, part, feature).ok()?;
        Some(Self {
            position: placed.origin_mm.map(|v| field(v, "mm")),
            facing: placed.facing,
            rotation: field(placed.rotation_deg, "deg"),
            pose,
            shown: None,
            edited: false,
        })
    }
    /// Show the draft's panel-feature origin and plane where the Part puts
    /// them, unless the user has typed or picked a placement here.
    pub fn follow(&mut self, origin_mm: Option<[f64; 3]>, plane: PanelPlane) {
        let Some(origin_mm) = origin_mm else { return };
        if self.edited || self.shown == Some((origin_mm, plane)) {
            return;
        }
        let first = self.shown.is_none();
        self.shown = Some((origin_mm, plane));
        if first {
            // Opened from the document: the fields already show it exactly.
            return;
        }
        let (u, v) = plane.axes();
        let placed = model::placement_of(
            self.pose.point(origin_mm),
            self.pose.vector(u),
            self.pose.vector(v),
        );
        self.position = placed.origin_mm.map(|v| field(v, "mm"));
        self.facing = placed.facing;
        self.rotation = field(placed.rotation_deg, "deg");
    }
    pub fn show(&mut self, ui: &mut egui::Ui) -> bool {
        let mut changed = false;
        ui.label("Placement");
        for (f, label) in self.position.iter_mut().zip(POSITION) {
            changed |= f.show(
                ui,
                label,
                "World position of the panel's origin corner, where the panel is now",
            );
        }
        ui.label(FACING);
        controls::inspector_combo(ui, "panel-facing")
            .selected_text(facing_label(self.facing))
            .show_ui(ui, |ui| {
                for facing in PanelFacing::ALL {
                    changed |= ui
                        .selectable_value(&mut self.facing, Some(facing), facing.label())
                        .changed();
                }
            });
        changed |= self.rotation.show(
            ui,
            ROTATION,
            "Turn of the panel about its face direction, counter-clockwise looking at the face",
        );
        if self.facing.is_none() {
            ui.small("Choose a face direction to place this panel from here.");
        }
        self.edited |= changed;
        changed
    }
    /// The placement edit, or None while nothing was typed or picked here.
    pub fn request(&self, document: &Document) -> Result<Option<PanelPlacementRequest>, String> {
        if !self.edited {
            return Ok(None);
        }
        let facing = self
            .facing
            .ok_or("Choose the direction the panel face looks along")?;
        Ok(Some(PanelPlacementRequest {
            position: [
                self.position[0].resolve(document)?,
                self.position[1].resolve(document)?,
                self.position[2].resolve(document)?,
            ],
            facing,
            rotation: self.rotation.resolve(document)?,
        }))
    }
}
