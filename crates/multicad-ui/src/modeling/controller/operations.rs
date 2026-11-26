//! The cut discovery reply cannot authorize Apply. It only selects targets for
//! a second, complete transaction evaluated on the original source Session.
use super::*;
use crate::modeling::extrude_operation::DiscoveryMode;
use multicad_core::model::{AutomaticOperation, Operation};

impl ModelToolController {
    pub(super) fn operation_changed(&mut self) {
        let through_sheet = self
            .operations
            .as_ref()
            .is_some_and(OperationControls::through_sheet);
        if through_sheet {
            self.operation.entry = None;
            self.operation.drag = None;
            self.operation.confirm = None;
            self.operation.touched = false;
        }
        if let Some(ModelForm::Extrude(panel)) = &mut self.form {
            panel.set_through_sheet(through_sheet);
            panel.changed();
        }
        self.invalidate();
    }
    pub fn set_extrusion_context(
        &mut self,
        session: &Session,
        resources: &ProjectResources,
        scene: &Scene,
    ) -> Result<(), String> {
        let result = (|| {
            self.verify_context(session, resources, scene)?;
            self.operations
                .as_mut()
                .ok_or("Start a new extrusion before choosing operation targets")?
                .set_context(session.document(), session.revision(), scene)
        })();
        self.context_result(scene.stamp(), &result);
        result
    }
    pub(super) fn submit_extrusion(
        &mut self,
        session: &Session,
        resources: &ProjectResources,
        channel: Channel,
        form: FormStamp,
        extrusion: &ExtrudeRequest,
    ) -> Result<(), String> {
        self.current(session, resources, true)?;
        if self.stamp() != Some(form) {
            return Err("Extrusion inputs changed; Preview the current values".into());
        }
        let controls = self
            .operations
            .as_ref()
            .ok_or("Start a new extrusion first")?;
        if controls.automatic() {
            if !controls.needs_discovery() {
                // Nothing to touch: New body without a native round trip.
                self.operations
                    .as_mut()
                    .expect("checked operation controls")
                    .decide(
                        AutomaticOperation::NewBody,
                        model::CutTargetSelection {
                            defaults: vec![],
                            targets: vec![],
                            not_intersecting: vec![],
                        },
                    );
                let request = self.new_body_request(extrusion, true);
                return self.submit_final(session, resources, channel, form, &request);
            }
            // Each preview decides again from the current distance and
            // direction until the author picks an operation (CAD-404).
            return self.start_discovery(
                session,
                resources,
                channel,
                form,
                extrusion,
                DiscoveryMode::Automatic,
            );
        }
        match controls.operation {
            Operation::NewBody => {
                let request = self.new_body_request(extrusion, false);
                self.submit_final(session, resources, channel, form, &request)
            }
            Operation::Cut if controls.through_sheet() => {
                let target = controls.sheet_target.as_ref().map_err(Clone::clone)?;
                let [output] = extrusion.outputs.as_slice() else {
                    return Err(
                        "Choose one simple closed sketch region for a through-sheet cut".into(),
                    );
                };
                let request = ModelRequest::SheetProfileCut {
                    feature: output.feature,
                    definition: Box::new(model::SheetProfileCutDefinition {
                        input: target.input.clone(),
                        output: model::BodyOutput {
                            body: target.input.body,
                        },
                        profile: model::ProfileSource {
                            sketch: extrusion.sketch,
                            region: output.region.clone(),
                        },
                        face: target.face,
                    }),
                };
                self.submit_final(session, resources, channel, form, &request)
            }
            Operation::Cut => self.start_discovery(
                session,
                resources,
                channel,
                form,
                extrusion,
                DiscoveryMode::Cut,
            ),
            Operation::Join | Operation::Intersect => {
                if controls.automatic_join() {
                    return self.start_discovery(
                        session,
                        resources,
                        channel,
                        form,
                        extrusion,
                        DiscoveryMode::Join,
                    );
                }
                let targets = controls.explicit_targets(session.document())?;
                let request = self.identities.compose(
                    extrusion.clone(),
                    controls.operation,
                    targets,
                    controls.keep_tools,
                )?;
                self.submit_final(
                    session,
                    resources,
                    channel,
                    form,
                    &ModelRequest::ExtrudeBoolean(Box::new(request)),
                )
            }
            Operation::Split => {
                Err("Split extrusion outputs are not available in this form".into())
            }
        }
    }
    fn new_body_request(&mut self, extrusion: &ExtrudeRequest, automatic: bool) -> ModelRequest {
        let mut group = self.identities.group_draft();
        group.automatic = automatic;
        ModelRequest::GroupedExtrude(Box::new(multicad_core::model::GroupedExtrudeRequest {
            extrusion: extrusion.clone(),
            group,
        }))
    }
    /// Ask the native model which bodies the extrusion tools affect. Join and
    /// the automatic default use the active scope only; Cut adds its explicit
    /// overrides.
    fn start_discovery(
        &mut self,
        session: &Session,
        resources: &ProjectResources,
        channel: Channel,
        form: FormStamp,
        extrusion: &ExtrudeRequest,
        mode: DiscoveryMode,
    ) -> Result<(), String> {
        let controls = self
            .operations
            .as_ref()
            .ok_or("Start a new extrusion first")?;
        let (include, exclude) = if mode != DiscoveryMode::Cut {
            (vec![], vec![])
        } else {
            (
                controls.include.iter().copied().collect(),
                controls.exclude.iter().copied().collect(),
            )
        };
        let cut = CutDiscovery::prepare(
            session,
            Actor::local(channel),
            extrusion.clone(),
            controls.component,
            include,
            exclude,
            mode,
        )?;
        self.generation = self
            .generation
            .checked_add(1)
            .ok_or("Model preview generation limit reached")?;
        let stamp = SceneStamp {
            document: session.document().id(),
            revision: session.revision(),
            generation: self.generation,
        };
        self.backend
            .as_mut()
            .ok_or("The native modeling preview provider is unavailable")?
            .request_evaluation(cut.prepared.document(), resources, stamp, &cut.options)?;
        self.discovery = Some(Discovery {
            cut,
            form,
            stamp,
            channel,
        });
        self.candidate = None;
        self.error = None;
        Ok(())
    }
    pub(super) fn receive_discovery(
        &mut self,
        session: &Session,
        resources: &ProjectResources,
        reply: crate::viewport::SceneReply,
    ) -> Result<(), String> {
        let discovery = self
            .discovery
            .take()
            .ok_or("The cut query is no longer active")?;
        self.current(session, resources, true)?;
        if self.stamp() != Some(discovery.form)
            || ContentHash::document(discovery.cut.prepared.document())? != discovery.cut.content
        {
            return Err("Extrusion controls changed; Preview the current cut".into());
        }
        let scene = reply.result?;
        let selections: Vec<_> = discovery
            .cut
            .selections(session, &scene, discovery.stamp)?
            .into_iter()
            .cloned()
            .collect();
        if discovery.cut.mode() == DiscoveryMode::Automatic {
            return self.receive_automatic(session, resources, discovery, selections);
        }
        let [selection] = <[_; 1]>::try_from(selections)
            .map_err(|_| "Cut target results do not match the requested probe")?;
        let controls = self
            .operations
            .as_mut()
            .ok_or("The extrusion operation was cancelled")?;
        let operation = if discovery.cut.mode() == DiscoveryMode::Join {
            Operation::Join
        } else {
            Operation::Cut
        };
        if controls.operation != operation
            || (operation == Operation::Join && !controls.automatic_join())
        {
            return Err("The extrusion operation changed; Preview again".into());
        }
        controls.receive(selection.clone());
        let targets = if operation == Operation::Join {
            // Never pick one of several plausible bodies silently.
            match selection.single_target() {
                Ok(target) => vec![target.clone()],
                Err(model::AutomaticTargetError::None) => vec![],
                Err(model::AutomaticTargetError::Several(several)) => {
                    return Err(format!(
                        "Several bodies touch or overlap the extrusion: {}; choose explicit targets",
                        several
                            .iter()
                            .map(|target| controls.label(target.body))
                            .collect::<Vec<_>>()
                            .join(", ")
                    ));
                }
            }
        } else {
            selection.targets
        };
        let mut request = self.identities.compose(
            discovery.cut.extrusion().clone(),
            operation,
            targets,
            controls.keep_tools,
        )?;
        if operation == Operation::Cut {
            request
                .group
                .as_mut()
                .expect("pinned operation record")
                .captured_scope = Some(multicad_core::model::ExtrudeCapturedScope {
                component: controls.component,
                include: controls.include.iter().copied().collect(),
                exclude: controls.exclude.iter().copied().collect(),
            });
        }
        self.submit_final(
            session,
            resources,
            discovery.channel,
            discovery.form,
            &ModelRequest::ExtrudeBoolean(Box::new(request)),
        )
    }
    /// Turn the two automatic probes into the default operation and submit
    /// it, unless the author has picked an operation meanwhile.
    fn receive_automatic(
        &mut self,
        session: &Session,
        resources: &ProjectResources,
        discovery: Discovery,
        selections: Vec<model::CutTargetSelection>,
    ) -> Result<(), String> {
        let [contact, volume] = <[_; 2]>::try_from(selections)
            .map_err(|_| "Automatic operation results do not match the requested probes")?;
        let decision = AutomaticOperation::decide(&contact, &volume);
        let controls = self
            .operations
            .as_mut()
            .ok_or("The extrusion operation was cancelled")?;
        if !controls.decide(decision.clone(), contact) {
            return Err("The extrusion operation changed; Preview again".into());
        }
        let keep_tools = controls.keep_tools;
        let component = controls.component;
        let extrusion = discovery.cut.extrusion().clone();
        let request = match (decision.operation(), decision.target()) {
            (operation @ (Operation::Join | Operation::Cut), Some(target)) => {
                let mut request = self.identities.compose(
                    extrusion.clone(),
                    operation,
                    vec![target.clone()],
                    keep_tools,
                )?;
                let group = request.group.as_mut().expect("pinned operation record");
                group.automatic = true;
                if operation == Operation::Cut {
                    group.captured_scope = Some(multicad_core::model::ExtrudeCapturedScope {
                        component,
                        include: vec![],
                        exclude: vec![],
                    });
                }
                ModelRequest::ExtrudeBoolean(Box::new(request))
            }
            _ => self.new_body_request(&extrusion, true),
        };
        let cut = decision.operation() == Operation::Cut;
        self.submit_final(
            session,
            resources,
            discovery.channel,
            discovery.form,
            &request,
        )?;
        if cut {
            self.cut_fallback = Some(super::CutFallback {
                form: discovery.form,
                channel: discovery.channel,
                extrusion,
            });
        }
        Ok(())
    }
    /// An automatic Cut whose native result removes all material of its
    /// target falls back to New body with a note. An explicit Cut keeps the
    /// error, because the author asked for exactly that operation.
    pub(super) fn automatic_cut_fallback(
        &mut self,
        session: &Session,
        resources: &ProjectResources,
        form: FormStamp,
    ) {
        if !self
            .error
            .as_deref()
            .is_some_and(|error| error.contains("removed all material"))
        {
            return;
        }
        let Some(fallback) = self.cut_fallback.take().filter(|f| f.form == form) else {
            return;
        };
        if !self
            .operations
            .as_mut()
            .is_some_and(OperationControls::cut_removes_all)
        {
            return;
        }
        self.error = None;
        let request = self.new_body_request(&fallback.extrusion, true);
        if let Err(error) = self.submit_final(
            session,
            resources,
            fallback.channel,
            fallback.form,
            &request,
        ) {
            self.error = Some(error);
        }
    }
}
