//! Tool-first Extrude / Revolve (CAD-400): the form opens before any profile
//! is picked, and every viewer pick replaces its profiles live while the
//! distance, angle and operation fields stay as the user set them. Nothing
//! is committed until the form's Apply or Enter.
use super::*;
use multicad_core::model::FaceInput;

/// The profiles one Extrude / Revolve form reads: closed areas of one
/// sketch, or coplanar planar model faces.
#[derive(Clone, Debug, PartialEq)]
pub enum ProfileTarget {
    Sketch {
        sketch: EntityId,
        regions: Vec<RegionIntent>,
    },
    Faces(Vec<FaceInput>),
}

impl ModelToolController {
    /// A new Extrude / Revolve form still waiting for its first profile:
    /// it neither previews nor applies.
    pub(super) fn profiles_missing(&self) -> bool {
        match &self.form {
            Some(ModelForm::Extrude(panel)) => panel.profiles_missing(),
            Some(ModelForm::Revolve(panel)) => panel.profiles_missing(),
            _ => false,
        }
    }
    /// Open the Extrude or Revolve form with no profile yet. Preview and
    /// Apply refuse with `missing` until [`Self::retarget_profiles`].
    pub fn begin_pending_profiles(
        &mut self,
        session: &Session,
        resources: &ProjectResources,
        tool: super::super::ModelingTool,
        missing: String,
        channel: Channel,
    ) -> Result<(), String> {
        match tool {
            super::super::ModelingTool::Extrude => {
                let mut panel = ExtrudePanel::default();
                panel.begin_pending(missing)?;
                self.begin_form(
                    session,
                    resources,
                    channel,
                    ModelForm::Extrude(Box::new(panel)),
                )?;
                self.operations = Some(OperationControls::new(session.document(), None));
                self.identities = BooleanIdentities::new(session.document());
            }
            super::super::ModelingTool::Revolve => self.begin_form(
                session,
                resources,
                channel,
                ModelForm::Revolve(Box::new(RevolvePanel::pending(missing)?)),
            )?,
            _ => return Err("Choose a profile modeling tool".into()),
        }
        Ok(())
    }

    /// Replace the profiles of the open new Extrude / Revolve form. `None`
    /// empties them; Preview and Apply then refuse with `missing`. The
    /// display guide is cleared and must be set again by the host.
    pub fn retarget_profiles(
        &mut self,
        session: &Session,
        target: Option<ProfileTarget>,
        missing: String,
    ) -> Result<(), String> {
        let document = session.document();
        let mut names = GeneratedNames::new(document);
        match self.form.as_mut() {
            Some(ModelForm::Extrude(panel)) => {
                let request = target.map(|target| {
                    let mut request = match target {
                        ProfileTarget::Sketch { sketch, regions } => {
                            ExtrudeRequest::new(sketch, regions, Default::default())
                        }
                        ProfileTarget::Faces(faces) => {
                            ExtrudeRequest::from_faces(faces, Default::default())
                        }
                    };
                    for output in &mut request.outputs {
                        output.name = names.numbered("Extruded part", " extrusion");
                    }
                    request
                });
                let before = panel.profile_sketch();
                let after = request.as_ref().map(|request| request.sketch);
                panel.retarget(request, missing)?;
                // Operation targets and the through-sheet mode depend on the
                // profile sketch; keep the user's choices while it is the same.
                if let Some(sketch) = after
                    && before != Some(sketch)
                {
                    let operations = OperationControls::new(
                        document,
                        (sketch != EntityId::nil()).then_some(sketch),
                    );
                    panel.set_through_sheet(operations.through_sheet());
                    self.operations = Some(operations);
                }
            }
            Some(ModelForm::Revolve(panel)) => {
                let request = target.map(|target| {
                    let mut request = match target {
                        ProfileTarget::Sketch { sketch, regions } => {
                            RevolveRequest::new(sketch, regions, Default::default())
                        }
                        ProfileTarget::Faces(faces) => {
                            RevolveRequest::from_faces(faces, Default::default())
                        }
                    };
                    for output in &mut request.outputs {
                        output.name = names.numbered("Revolved part", " feature");
                    }
                    request
                });
                panel.retarget(request, missing);
            }
            _ => return Err("Start Extrude or Revolve first".into()),
        }
        self.operation.guide = None;
        self.operation.drag = None;
        self.operation.confirm = None;
        self.keyboard_confirm = None;
        self.error = None;
        Ok(())
    }
}
