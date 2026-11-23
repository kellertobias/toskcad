//! Native preflight and write receipts for one project's print plates.
use super::*;
use multicad_core::commands::SessionStamp;

/// Native preflight results kept for reselecting a plate. Each is keyed by
/// the complete session stamp, plate and arrangement it answered.
const MAX_RECEIPTS: usize = 16;

#[derive(Clone, PartialEq)]
struct Source {
    session: SessionStamp,
    document: EntityId,
    revision: u64,
    plate: EntityId,
    definition: PrintPlateDefinition,
}

struct Job {
    id: EntityId,
    source: Source,
    actor: Actor,
    expected: Option<PrintKey>,
}

#[derive(Default)]
pub(crate) struct PlateJobs {
    backend: Option<Box<dyn PrintBackend>>,
    source: Option<Source>,
    job: Option<Job>,
    due: Option<f64>,
    preview: Option<PrintPreview>,
    error: Option<String>,
    notice: Option<(String, bool)>,
    /// Completed previews of this session revision, most recent last. Only
    /// an identical source (session, revision, plate and definition) reuses
    /// one; writing still rebuilds and checks the current source natively.
    receipts: Vec<(Source, PrintPreview)>,
}

impl PlateJobs {
    pub(crate) fn set_backend(&mut self, backend: Option<Box<dyn PrintBackend>>) {
        self.backend = backend;
        self.source = None;
        self.preview = None;
        self.error = None;
        self.receipts.clear();
    }

    pub(crate) fn has_backend(&self) -> bool {
        self.backend.is_some()
    }

    pub(crate) fn busy(&self) -> bool {
        self.job.is_some()
    }

    pub(crate) fn scheduled(&self) -> bool {
        self.due.is_some()
    }

    pub(crate) fn writing(&self) -> bool {
        self.job.as_ref().is_some_and(|job| job.expected.is_some())
    }

    pub(crate) fn preview(&self) -> Option<&PrintPreview> {
        self.preview.as_ref()
    }

    pub(crate) fn error(&self) -> Option<&str> {
        self.error.as_deref()
    }

    pub(crate) fn take_notice(&mut self) -> Option<(String, bool)> {
        self.notice.take()
    }

    pub(crate) fn cancel(&mut self) {
        self.due = None;
        if let Some(job) = &self.job
            && let Some(backend) = &mut self.backend
        {
            backend.cancel(job.id);
        }
    }

    pub(crate) fn retry(&mut self, time: f64) {
        self.error = None;
        self.due = Some(time);
    }

    pub(crate) fn update(
        &mut self,
        s: &Session,
        plate: EntityId,
        definition: &PrintPlateDefinition,
        time: f64,
        allow_io: bool,
    ) {
        if self.source.as_ref().is_none_or(|source| {
            source.session != s.stamp() || source.plate != plate || source.definition != *definition
        }) {
            if !self.writing() {
                self.cancel();
            }
            let source = Source {
                session: s.stamp(),
                document: s.document().id(),
                revision: s.revision(),
                plate,
                definition: definition.clone(),
            };
            // A receipt from another revision can never answer this one.
            self.receipts
                .retain(|(receipt, _)| receipt.session == source.session);
            self.preview = self
                .receipts
                .iter()
                .find(|(receipt, _)| *receipt == source)
                .map(|(_, preview)| preview.clone());
            self.error = None;
            self.due = (self.preview.is_none() && !definition.placements.is_empty())
                .then_some(time + 0.15);
            self.source = Some(source);
        }
        if allow_io && !self.busy() && self.due.is_some_and(|due| due <= time) {
            self.due = None;
            if let Err(error) = self.request(s) {
                self.error = Some(error);
            }
        }
    }

    fn request(&mut self, s: &Session) -> Result<(), String> {
        let source = self.source.as_ref().ok_or("Choose a print plate")?;
        let actor = Actor::local(Channel::Ui);
        let id = self
            .backend
            .as_mut()
            .ok_or("Native print plate service is unavailable")?
            .preview_plate(
                s.document(),
                s.revision(),
                &source.definition,
                actor.clone(),
            )
            .map_err(|error| error.to_string())?;
        self.job = Some(Job {
            id,
            source: source.clone(),
            actor,
            expected: None,
        });
        Ok(())
    }

    pub(crate) fn can_write(
        &self,
        s: &Session,
        plate: EntityId,
        definition: &PrintPlateDefinition,
    ) -> bool {
        !self.busy()
            && self.source.as_ref().is_some_and(|source| {
                source.session == s.stamp()
                    && source.plate == plate
                    && source.definition == *definition
            })
            && self.preview.as_ref().is_some_and(PrintPreview::can_write)
    }

    pub(crate) fn write(
        &mut self,
        s: &Session,
        plate: EntityId,
        definition: &PrintPlateDefinition,
        path: &Path,
    ) -> Result<(), String> {
        if !self.can_write(s, plate, definition) {
            return Err("Wait for the current plate to pass native preflight".into());
        }
        let path = super::three_mf_path(path)?;
        let expected = self
            .preview
            .as_ref()
            .expect("ready plate preview")
            .key
            .clone();
        expected
            .source
            .check(s.document(), s.revision())
            .map_err(|error| error.to_string())?;
        let actor = Actor::local(Channel::Ui);
        let id = self
            .backend
            .as_mut()
            .ok_or("Native print plate service is unavailable")?
            .write_plate(
                s.document(),
                s.revision(),
                definition,
                &expected,
                &path,
                actor.clone(),
            )
            .map_err(|error| error.to_string())?;
        self.job = Some(Job {
            id,
            source: self.source.clone().expect("ready plate source"),
            actor,
            expected: Some(expected),
        });
        self.error = None;
        Ok(())
    }

    pub(crate) fn poll(&mut self, s: &Session) {
        while let Some(reply) = self.backend.as_mut().and_then(|backend| backend.poll()) {
            if self.job.as_ref().is_none_or(|job| job.id != reply.id) {
                continue;
            }
            let job = self.job.take().expect("matching job");
            let current =
                self.source.as_ref() == Some(&job.source) && job.source.session == s.stamp();
            let write = job.expected.is_some();
            if !current && !write {
                continue;
            }
            if reply.document != job.source.document
                || reply.revision != job.source.revision
                || reply.actor != job.actor
            {
                self.failed(
                    "Native print plate reply did not match its request".into(),
                    write,
                );
                continue;
            }
            match (reply.result, job.expected) {
                (PrintResult::Preview(Ok(preview)), None)
                    if preview.requested_by == job.actor
                        && preview.request
                            == job.source.definition.summary_request(s.document())
                        && preview.selection.stamp == preview.key.source
                        && preview.key.source.check(s.document(), s.revision()).is_ok() =>
                {
                    self.receipts.retain(|(receipt, _)| {
                        receipt.session == job.source.session && *receipt != job.source
                    });
                    if self.receipts.len() == MAX_RECEIPTS {
                        self.receipts.remove(0);
                    }
                    self.receipts.push((job.source, (*preview).clone()));
                    self.preview = Some(*preview);
                    self.error = None;
                }
                (PrintResult::Preview(Err(error)), None)
                | (PrintResult::Write(Err(error)), Some(_)) => {
                    if error.kind == ErrorKind::Cancelled && write {
                        self.notice = Some(("Print plate write canceled".into(), false));
                    } else {
                        self.failed(error.to_string(), write);
                    }
                }
                (PrintResult::Write(Ok(receipt)), Some(expected))
                    if receipt.actor == job.actor && receipt.key == expected =>
                {
                    let path = receipt
                        .files
                        .first()
                        .map_or("the selected folder", |file| file.path.as_str());
                    self.notice = Some((
                        format!("Wrote {} print plate file(s) · {path}", receipt.files.len()),
                        false,
                    ));
                }
                _ => self.failed(
                    "Native print plate receipt did not match its request".into(),
                    write,
                ),
            }
        }
    }

    fn failed(&mut self, error: String, write: bool) {
        if write {
            self.notice = Some((error, true));
        } else {
            self.error = Some(error);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use multicad_core::{
        commands::Request,
        document::{Document, Entity, EntityData, Mutation},
        project_file::ContentHash,
    };
    use std::{cell::RefCell, rc::Rc};

    type Previews = Vec<(EntityId, u64, PrintPlateDefinition, Actor)>;
    #[derive(Default)]
    struct State {
        previews: Previews,
        writes: Vec<(u64, PrintKey)>,
        replies: VecDeque<PrintReply>,
    }
    struct Backend(Rc<RefCell<State>>);
    impl PrintBackend for Backend {
        fn preview(
            &mut self,
            _: &Document,
            _: u64,
            _: &PrintRequest,
            _: Actor,
        ) -> Result<EntityId, ExportError> {
            unreachable!("plates request plate previews")
        }
        fn write(
            &mut self,
            _: &Document,
            _: u64,
            _: &PrintRequest,
            _: &PrintKey,
            _: &Path,
            _: Actor,
        ) -> Result<EntityId, ExportError> {
            unreachable!("plates request plate writes")
        }
        fn preview_plate(
            &mut self,
            _: &Document,
            revision: u64,
            plate: &PrintPlateDefinition,
            actor: Actor,
        ) -> Result<EntityId, ExportError> {
            let id = EntityId::new();
            let mut state = self.0.borrow_mut();
            state.previews.push((id, revision, plate.clone(), actor));
            Ok(id)
        }
        fn write_plate(
            &mut self,
            _: &Document,
            revision: u64,
            _: &PrintPlateDefinition,
            expected: &PrintKey,
            _: &Path,
            _: Actor,
        ) -> Result<EntityId, ExportError> {
            self.0
                .borrow_mut()
                .writes
                .push((revision, expected.clone()));
            Ok(EntityId::new())
        }
        fn cancel(&mut self, _: EntityId) -> bool {
            true
        }
        fn poll(&mut self) -> Option<PrintReply> {
            self.0.borrow_mut().replies.pop_front()
        }
        fn pending(&self) -> bool {
            false
        }
    }

    fn session() -> (Session, EntityId) {
        let mut s = Session::new(Document::new("Rack").unwrap()).unwrap();
        let part = Entity::new(
            "Rack body",
            EntityData::Part {
                features: vec![],
                material: None,
                tags: Default::default(),
            },
        );
        let id = part.id;
        s.apply(Request {
            expected_revision: 0,
            actor: Actor::local(Channel::Ui),
            label: "Part".into(),
            edits: vec![Mutation::Insert(part)],
        })
        .unwrap();
        (s, id)
    }

    fn plate(part: EntityId, x: f64) -> PrintPlateDefinition {
        PrintPlateDefinition {
            bed_mm: [350., 350.],
            placements: vec![PrintPlatePlacement {
                source: PartRef {
                    part,
                    occurrence_path: vec![],
                    root_assembly: None,
                },
                position_mm: [x, 4.],
                rotation_deg: [-90., 0., 0.],
            }],
            ..Default::default()
        }
    }

    /// Answer the latest native request as a writable preview.
    fn answer(s: &Session, state: &Rc<RefCell<State>>) {
        let (id, revision, plate, actor) = state.borrow().previews.last().unwrap().clone();
        let stamp = ExportStamp::capture(s.document(), revision).unwrap();
        let preview = PrintPreview {
            requested_by: actor.clone(),
            key: PrintKey {
                source: stamp.clone(),
                plan_hash: ContentHash::of(&serde_json::to_vec(&plate).unwrap()),
            },
            provider: "test".into(),
            request: plate.summary_request(s.document()),
            selection: PrintSelection {
                stamp,
                parts: vec![],
                skipped: vec![],
                processes: vec![],
            },
            objects: vec![PrintObject {
                part: plate.placements[0].source.part,
                occurrence_path: vec![],
                name: "Rack body".into(),
                role: PrintRole::Finished,
                bounds_mm: [[0.; 3], [1.; 3]],
                plate_pose: None,
                volume_mm3: 1.,
                triangles: 12,
                mesh_volume_mm3: 1.,
            }],
            files: vec![],
            processes: vec![],
            findings: vec![],
            blocked: None,
        };
        state.borrow_mut().replies.push_back(PrintReply {
            id,
            document: s.document().id(),
            revision,
            actor,
            result: PrintResult::Preview(Ok(Box::new(preview))),
        });
    }

    /// Select `plate` and let any due native check run and finish.
    fn select(
        jobs: &mut PlateJobs,
        s: &Session,
        state: &Rc<RefCell<State>>,
        id: EntityId,
        plate: &PrintPlateDefinition,
        time: &mut f64,
    ) -> bool {
        let before = state.borrow().previews.len();
        jobs.update(s, id, plate, *time, true);
        *time += 1.;
        jobs.update(s, id, plate, *time, true);
        let requested = state.borrow().previews.len() > before;
        if requested {
            answer(s, state);
            jobs.poll(s);
        }
        requested
    }

    #[test]
    fn reselecting_a_plate_reuses_its_receipt_until_source_or_arrangement_changes() {
        let (mut s, part) = session();
        let state = Rc::new(RefCell::new(State::default()));
        let mut jobs = PlateJobs::default();
        jobs.set_backend(Some(Box::new(Backend(state.clone()))));
        let (first, second) = (EntityId::new(), EntityId::new());
        let (rack, other) = (plate(part, 4.), plate(part, 40.));
        let mut time = 0.;
        // The first selection runs native preflight.
        assert!(select(&mut jobs, &s, &state, first, &rack, &mut time));
        assert!(jobs.can_write(&s, first, &rack));
        let receipt = jobs.preview().unwrap().clone();
        assert!(select(&mut jobs, &s, &state, second, &other, &mut time));
        assert!(!jobs.can_write(&s, first, &rack));
        // Reselecting the same plate at the same revision is ready at once,
        // with the identical receipt and no native request.
        jobs.update(&s, first, &rack, time, true);
        assert!(!jobs.scheduled() && !jobs.busy());
        assert!(jobs.can_write(&s, first, &rack));
        assert_eq!(jobs.preview(), Some(&receipt));
        assert_eq!(state.borrow().previews.len(), 2);
        // Moving or turning the placement is another arrangement.
        let mut moved = rack.clone();
        moved.placements[0].position_mm[0] += 1.;
        assert!(select(&mut jobs, &s, &state, first, &moved, &mut time));
        let mut turned = rack.clone();
        turned.placements[0].rotation_deg[2] = 90.;
        assert!(select(&mut jobs, &s, &state, first, &turned, &mut time));
        // The same arrangement under another plate id is checked again.
        assert!(select(&mut jobs, &s, &state, second, &rack, &mut time));
        // Back to a receipted arrangement: reused.
        assert!(!select(&mut jobs, &s, &state, first, &rack, &mut time));
        assert!(jobs.can_write(&s, first, &rack));
        // Any source change (a new revision) discards every receipt.
        s.apply(Request {
            expected_revision: s.revision(),
            actor: Actor::local(Channel::Mcp),
            label: "Edit".into(),
            edits: vec![Mutation::RenameDocument("Rack v2".into())],
        })
        .unwrap();
        assert!(!jobs.can_write(&s, first, &rack));
        jobs.update(&s, first, &rack, time, true);
        assert!(jobs.scheduled() && jobs.preview().is_none());
        assert!(jobs.receipts.is_empty());
        time += 1.;
        jobs.update(&s, first, &rack, time, true);
        assert_eq!(state.borrow().previews.last().unwrap().1, s.revision());
        answer(&s, &state);
        jobs.poll(&s);
        assert!(jobs.can_write(&s, first, &rack));
        // Writing from a reused receipt still asks the native service to
        // rebuild and check the current source against its key.
        assert!(!select(&mut jobs, &s, &state, first, &rack, &mut time));
        let path = std::env::temp_dir().join("plate.3mf");
        jobs.write(&s, first, &rack, &path).unwrap();
        let writes = state.borrow().writes.clone();
        assert_eq!(writes.len(), 1);
        assert_eq!(
            writes[0],
            (s.revision(), jobs.preview().unwrap().key.clone())
        );
    }

    #[test]
    fn receipts_are_bounded_and_a_new_backend_forgets_them() {
        let (s, part) = session();
        let state = Rc::new(RefCell::new(State::default()));
        let mut jobs = PlateJobs::default();
        jobs.set_backend(Some(Box::new(Backend(state.clone()))));
        let id = EntityId::new();
        let mut time = 0.;
        for n in 0..=MAX_RECEIPTS {
            assert!(select(
                &mut jobs,
                &s,
                &state,
                id,
                &plate(part, n as f64),
                &mut time
            ));
        }
        assert_eq!(jobs.receipts.len(), MAX_RECEIPTS);
        // The oldest arrangement was evicted; the newest are kept.
        assert!(select(
            &mut jobs,
            &s,
            &state,
            id,
            &plate(part, 0.),
            &mut time
        ));
        assert!(!select(
            &mut jobs,
            &s,
            &state,
            id,
            &plate(part, MAX_RECEIPTS as f64),
            &mut time
        ));
        jobs.set_backend(Some(Box::new(Backend(state.clone()))));
        assert!(jobs.receipts.is_empty());
        assert!(select(
            &mut jobs,
            &s,
            &state,
            id,
            &plate(part, 1.),
            &mut time
        ));
    }
}
