//! CAD-347 headless response harness. It runs the production backend
//! composition (`install_backends`) and real `Application` frames through an
//! egui context, writing the response trace that
//! `tools/measure-response-trace.mjs` reads. There is no window, GPU upload,
//! presentation or operating-system input: actions are invoked as their
//! buttons invoke them and a result counts as shown when the frame that
//! accepted it finishes building. It is not a native interactive measurement.
//!
//! ```sh
//! HOME=$ISOLATED MULTICAD_RESPONSE_TRACE=$TRACE CAD347_PROJECT=examples/toolcase.tcd \
//!   cargo test -p multicad-app --release --bin multicad response_harness -- --ignored --nocapture
//! ```
//! Optional: `CAD347_PARAMETER=name` and `CAD347_VALUES=a,b,...` edit one
//! parameter per value; `CAD347_TOOLS=thread,gear` preview, cancel and preview
//! again with each modeling tool's defaults.
use super::{files, install_backends};
use multicad_core::{commands::Channel, response_trace as trace};
use multicad_model::GeometryBackend;
use multicad_ui::{
    application::Application,
    egui,
    modeling::ModelingTool,
    registry::{Command, Mode},
};
use std::{
    path::PathBuf,
    time::{Duration, Instant},
};

const TIMEOUT: Duration = Duration::from_secs(300);

struct Harness {
    app: Application,
    ctx: egui::Context,
    clock: Instant,
}

impl Harness {
    fn frame(&mut self) {
        let input = egui::RawInput {
            screen_rect: Some(egui::Rect::from_min_size(
                egui::Pos2::ZERO,
                egui::vec2(1728., 1080.),
            )),
            time: Some(self.clock.elapsed().as_secs_f64()),
            ..Default::default()
        };
        let app = &mut self.app;
        let _ = self.ctx.run_ui(input, |ui| app.show_frame(ui));
        // About one 120 Hz display interval between frames.
        std::thread::sleep(Duration::from_millis(8));
    }
    /// Frames until `done`, returning wall time from the call.
    fn until(&mut self, what: &str, mut done: impl FnMut(&mut Application) -> bool) -> Duration {
        let start = Instant::now();
        loop {
            self.frame();
            if done(&mut self.app) {
                return start.elapsed();
            }
            assert!(
                start.elapsed() < TIMEOUT,
                "{what} timed out; status {:?}",
                self.app.status()
            );
        }
    }
    fn action(&mut self, label: &str) -> Instant {
        trace::mark("harness-action", |record| {
            record.insert("label".into(), label.into());
        });
        Instant::now()
    }
    fn invoke(&mut self, command: Command) -> Result<(), String> {
        self.app.invoke(command, &[], Channel::Ui)
    }
    fn scene_current(app: &mut Application) -> bool {
        let p = app.current_mut();
        p.scene.error().is_none()
            && p.scene
                .is_current(p.session.document().id(), p.session.revision())
    }
    fn report(&self, label: &str, elapsed: Duration) {
        let (status, error) = self.app.status();
        eprintln!(
            "{label}: {:.1} ms; mode {:?}; status {status:?}{}",
            elapsed.as_secs_f64() * 1000.,
            self.app.mode(),
            if error { " (error)" } else { "" }
        );
    }
}

fn env(name: &str) -> Option<String> {
    std::env::var(name).ok().filter(|value| !value.is_empty())
}

#[test]
#[ignore = "timing harness; run explicitly with MULTICAD_RESPONSE_TRACE and CAD347_PROJECT"]
fn response_harness() {
    let path = trace::install_from_env()
        .unwrap()
        .expect("set MULTICAD_RESPONSE_TRACE");
    eprintln!("trace {}", path.display());
    let project = PathBuf::from(env("CAD347_PROJECT").expect("set CAD347_PROJECT"));
    let mut app = Application::welcome(Box::new(files::ProjectStore::default())).unwrap();
    install_backends(&mut app, GeometryBackend::Rust).unwrap();
    let mut h = Harness {
        app,
        ctx: egui::Context::default(),
        clock: Instant::now(),
    };

    let start = h.action("open project");
    h.app.request_open_documents([project.clone()]);
    h.until("project open", |app| {
        app.projects
            .iter()
            .any(|p| p.path.as_deref() == Some(&project))
    });
    let opened = start.elapsed();
    h.report("open project (document)", opened);
    h.invoke(Command::Mode(Mode::Model)).unwrap();
    let shown = h.until("first model", Harness::scene_current);
    h.report("first model after open", opened + shown);

    if let Some(name) = env("CAD347_PARAMETER") {
        h.invoke(Command::Mode(Mode::Parameters)).unwrap();
        h.frame();
        for value in env("CAD347_VALUES").unwrap_or_default().split(',') {
            let before = h.app.current().session.revision();
            let p = h.app.current_mut();
            let names: Vec<String> = p.editor.rows.iter().map(|row| row.name.clone()).collect();
            let row = p
                .editor
                .rows
                .iter_mut()
                .find(|row| row.name == name)
                .unwrap_or_else(|| panic!("no parameter {name}; parameters {names:?}"));
            row.expression = value.into();
            let start = h.action("parameter apply");
            h.invoke(Command::Apply).unwrap();
            assert_eq!(h.app.current().session.revision(), before + 1);
            let elapsed = start.elapsed() + h.until("parameter rebuild", Harness::scene_current);
            h.report(&format!("parameter {name} = {value}"), elapsed);
        }
    }

    let has_pages = multicad_core::drawings::books(h.app.current().session.document())
        .next()
        .is_some();
    if has_pages {
        for (label, leave) in [("drawings cold", false), ("drawings warm", true)] {
            if leave {
                h.invoke(Command::Mode(Mode::Model)).unwrap();
                h.until("model", Harness::scene_current);
            }
            let start = h.action(label);
            h.invoke(Command::Mode(Mode::Drawings)).unwrap();
            let elapsed = start.elapsed()
                + h.until(label, |app| {
                    let p = app.current_mut();
                    p.drawings.preview_outcome(&p.session).is_some()
                });
            let p = h.app.current();
            eprintln!(
                "  {label} outcome {:?}",
                p.drawings.preview_outcome(&p.session)
            );
            h.report(label, elapsed);
        }
        h.invoke(Command::Mode(Mode::Model)).unwrap();
        h.until("model", Harness::scene_current);
    }

    for tool in env("CAD347_TOOLS").unwrap_or_default().split(',') {
        let tool = match tool {
            "thread" => ModelingTool::Thread,
            "gear" => ModelingTool::Gear,
            "" => continue,
            other => panic!("unknown tool {other}"),
        };
        let revision = h.app.current().session.revision();
        h.action("open tool");
        h.invoke(Command::Modeling(tool)).unwrap();
        h.frame();
        // Preview, cancel at once, then preview to completion: the previous
        // model must stay current and the document unchanged throughout.
        let start = h.action("preview then cancel");
        let preview = h.invoke(Command::Preview);
        h.frame();
        h.invoke(Command::Discard).unwrap();
        let cancelled =
            start.elapsed() + h.until("cancel", |app| !app.current().sketch.model_active());
        h.report(
            &format!("{tool:?} preview {preview:?} then cancel"),
            cancelled,
        );
        assert_eq!(h.app.current().session.revision(), revision);
        assert!(
            Harness::scene_current(&mut h.app),
            "previous model was lost"
        );

        h.action("open tool");
        h.invoke(Command::Modeling(tool)).unwrap();
        h.frame();
        let start = h.action("preview");
        let preview = h.invoke(Command::Preview);
        let elapsed = start.elapsed()
            + h.until("preview", |app| {
                let p = app.current();
                preview.is_err()
                    || p.sketch.model_error().is_some()
                    || p.sketch.model_can_apply(&p.session, &p.resources)
            });
        let p = h.app.current();
        eprintln!(
            "  {tool:?} preview request {preview:?}; can apply {}; error {:?}",
            p.sketch.model_can_apply(&p.session, &p.resources),
            p.sketch.model_error()
        );
        h.report(&format!("{tool:?} preview"), elapsed);
        h.invoke(Command::Discard).unwrap();
        h.until("discard", |app| !app.current().sketch.model_active());
        assert_eq!(h.app.current().session.revision(), revision);
    }
    // Let the trace's last frame and worker records land.
    for _ in 0..4 {
        h.frame();
    }
}
