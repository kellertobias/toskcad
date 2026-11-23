//! Opt-in operator response trace (CAD-347).
//!
//! When `MULTICAD_RESPONSE_TRACE` names a file at application start, the UI,
//! worker and paint paths append one JSON object per line: observed inputs,
//! frame boundaries, timed stages (preparation, compilation, queries,
//! publication, GPU upload) and the moment each frame's draw commands are
//! encoded immediately before presentation. All times are seconds on one
//! monotonic clock that starts with the trace, so records from different
//! threads compare directly. Without the variable every call is a single
//! atomic load and nothing is allocated or written.
//!
//! The trace observes the application from inside the process. Input times
//! are when a frame first sees an event, not the operating-system event time,
//! and the last stamp is the encoded frame, not scan-out on the display.
use serde_json::{Map, Value};
use std::{
    fs::File,
    io::{LineWriter, Write},
    path::Path,
    sync::{Mutex, OnceLock},
    time::Instant,
};

/// Environment variable naming the trace file.
pub const ENV: &str = "MULTICAD_RESPONSE_TRACE";
/// Record format version; bump when a field changes meaning.
pub const VERSION: u64 = 1;

/// Pipeline stage of a timed span.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Stage {
    /// Main-thread work that builds a request (snapshots, validation, forms).
    Preparation,
    /// Exact source evaluation on a worker.
    Compilation,
    /// Typed questions answered from an evaluated model.
    Query,
    /// Turning results into display data and accepting them on the UI thread.
    Publication,
    /// Uploading display buffers to the GPU while a frame is painted.
    Upload,
}

impl Stage {
    pub const fn name(self) -> &'static str {
        match self {
            Self::Preparation => "preparation",
            Self::Compilation => "compilation",
            Self::Query => "query",
            Self::Publication => "publication",
            Self::Upload => "upload",
        }
    }
}

struct Tracer {
    epoch: Instant,
    out: Mutex<Box<dyn Write + Send>>,
}

impl Tracer {
    fn now(&self) -> f64 {
        self.epoch.elapsed().as_secs_f64()
    }
    fn write(&self, mut record: Map<String, Value>) {
        if let Some(name) = std::thread::current().name() {
            record.insert("thread".into(), name.into());
        }
        let Ok(mut line) = serde_json::to_vec(&Value::Object(record)) else {
            return;
        };
        line.push(b'\n');
        if let Ok(mut out) = self.out.lock() {
            // A failed diagnostic write must never disturb the application.
            let _ = out.write_all(&line);
        }
    }
}

static TRACER: OnceLock<Tracer> = OnceLock::new();

/// Start tracing to the file named by [`ENV`], if any. Returns the path when
/// tracing is active. Later calls keep the first destination.
pub fn install_from_env() -> Result<Option<std::path::PathBuf>, String> {
    let Some(path) = std::env::var_os(ENV).filter(|path| !path.is_empty()) else {
        return Ok(None);
    };
    let path = std::path::PathBuf::from(path);
    install_file(&path)?;
    Ok(Some(path))
}

/// Start tracing to a new file. An existing file is replaced.
pub fn install_file(path: &Path) -> Result<(), String> {
    let file = File::create(path)
        .map_err(|error| format!("Cannot create response trace {}: {error}", path.display()))?;
    install(Box::new(LineWriter::new(file)));
    Ok(())
}

fn install(out: Box<dyn Write + Send>) {
    let tracer = Tracer {
        epoch: Instant::now(),
        out: Mutex::new(out),
    };
    if TRACER.set(tracer).is_ok() {
        let build = crate::BuildInfo::current();
        mark("start", |record| {
            record.insert("version".into(), VERSION.into());
            record.insert("pid".into(), std::process::id().into());
            record.insert("product_version".into(), build.version.into());
            record.insert("commit".into(), build.commit.into());
            record.insert(
                "unix_time".into(),
                std::time::SystemTime::now()
                    .duration_since(std::time::UNIX_EPOCH)
                    .map_or(0., |d| d.as_secs_f64())
                    .into(),
            );
        });
    }
}

/// Whether a trace destination is active.
#[inline]
pub fn enabled() -> bool {
    TRACER.get().is_some()
}

/// Seconds on the trace clock, or `None` without a trace.
pub fn now() -> Option<f64> {
    TRACER.get().map(Tracer::now)
}

/// Append an instantaneous record. `fields` runs only while tracing.
pub fn mark(phase: &str, fields: impl FnOnce(&mut Map<String, Value>)) {
    let Some(tracer) = TRACER.get() else {
        return;
    };
    let mut record = Map::new();
    record.insert("phase".into(), phase.into());
    record.insert("t".into(), tracer.now().into());
    fields(&mut record);
    tracer.write(record);
}

/// A timed stage written when dropped, so early returns and errors are kept.
#[must_use = "a span measures until it is dropped"]
pub struct Span {
    stage: Stage,
    label: &'static str,
    start: f64,
    fields: Map<String, Value>,
}

impl Span {
    /// Attach a field, for example a revision, generation or outcome.
    pub fn set(&mut self, key: &str, value: impl Into<Value>) {
        self.fields.insert(key.into(), value.into());
    }
}

impl Drop for Span {
    fn drop(&mut self) {
        let Some(tracer) = TRACER.get() else {
            return;
        };
        let mut record = std::mem::take(&mut self.fields);
        record.insert("phase".into(), "span".into());
        record.insert("stage".into(), self.stage.name().into());
        record.insert("label".into(), self.label.into());
        record.insert("start".into(), self.start.into());
        record.insert("end".into(), tracer.now().into());
        tracer.write(record);
    }
}

/// Begin a timed stage; `None` without a trace.
pub fn span(stage: Stage, label: &'static str) -> Option<Span> {
    TRACER.get().map(|tracer| Span {
        stage,
        label,
        start: tracer.now(),
        fields: Map::new(),
    })
}

/// Attach a field to an optional span.
pub fn set(span: &mut Option<Span>, key: &str, value: impl Into<Value>) {
    if let Some(span) = span {
        span.set(key, value);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;

    #[derive(Clone, Default)]
    struct Shared(Arc<Mutex<Vec<u8>>>);
    impl Write for Shared {
        fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
            self.0.lock().unwrap().extend_from_slice(bytes);
            Ok(bytes.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }

    // One process-wide tracer: this is the only test that installs it.
    #[test]
    fn records_are_one_json_object_per_line_on_one_monotonic_clock() {
        assert!(!enabled());
        assert!(span(Stage::Compilation, "ignored").is_none());
        mark("ignored", |_| {
            panic!("fields must not be built without a trace")
        });
        let out = Shared::default();
        install(Box::new(out.clone()));
        install(Box::new(Shared::default()));
        assert!(enabled());
        {
            let mut outer = span(Stage::Compilation, "compile");
            set(&mut outer, "revision", 7);
            let inner = span(Stage::Query, "cut targets");
            mark("scene-request", |record| {
                record.insert("generation".into(), 3.into());
            });
            drop(inner);
        }
        let worker = std::thread::Builder::new()
            .name("worker".into())
            .spawn(|| drop(span(Stage::Publication, "scene")))
            .unwrap();
        worker.join().unwrap();
        let text = String::from_utf8(out.0.lock().unwrap().clone()).unwrap();
        let rows: Vec<Value> = text
            .lines()
            .map(|line| serde_json::from_str(line).unwrap())
            .collect();
        assert_eq!(rows.len(), 5, "{text}");
        assert_eq!(rows[0]["phase"], "start");
        assert_eq!(rows[0]["version"], VERSION);
        assert_eq!(rows[1]["phase"], "scene-request");
        assert_eq!(rows[1]["generation"], 3);
        let (query, compile, publish) = (&rows[2], &rows[3], &rows[4]);
        assert_eq!(query["stage"], "query");
        assert_eq!(compile["stage"], "compilation");
        assert_eq!(compile["label"], "compile");
        assert_eq!(compile["revision"], 7);
        assert_eq!(publish["thread"], "worker");
        let time = |row: &Value, key: &str| row[key].as_f64().unwrap();
        assert!(time(compile, "start") <= time(query, "start"));
        assert!(time(query, "end") <= time(compile, "end"));
        assert!(time(&rows[1], "t") <= time(query, "end"));
        assert!(time(compile, "end") <= time(publish, "start"));
        assert!(now().unwrap() >= time(publish, "end"));
    }
}
