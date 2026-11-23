//! CAD-330 timing harness: the saved Werkzeugkoffer "Rack on its back" plate
//! through the print worker the native app uses, and the same build split by
//! stage. Headless; it excludes window, frame scheduling and model display
//! work, so it is not a native interactive measurement.
//!
//! cargo test -p multicad-app --release --test rack_plate_timing -- --ignored --nocapture
use multicad_core::{
    commands::{Actor, Channel},
    document::Document,
    exchange::read_print_plate,
    project_file,
};
use multicad_model::{
    GeometryBackend,
    export::{Exporter, PrintWorker},
};
use std::{
    path::Path,
    time::{Duration, Instant},
};

fn toolcase() -> Document {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/toolcase.tcd");
    project_file::read(
        &mut std::fs::File::open(path).unwrap(),
        Default::default(),
        &project_file::supported_capabilities(),
    )
    .unwrap()
    .project
    .document
}

#[test]
#[ignore = "timing harness; run explicitly in release"]
fn rack_plate_first_check_timings() {
    let document = toolcase();
    let id = document
        .entities()
        .find(|entity| entity.name == "Rack on its back")
        .unwrap()
        .id;
    let plate = read_print_plate(&document, id).unwrap();
    let actor = || Actor::local(Channel::Ui);
    // A fresh worker per sample, as for the first selection in a new session.
    for sample in 0..5 {
        let worker =
            PrintWorker::with_geometry_backend(Exporter::default(), GeometryBackend::Rust).unwrap();
        let mut rounds = vec![];
        for _ in 0..2 {
            let started = Instant::now();
            let job = worker.preview_plate(&document, 0, &plate, actor()).unwrap();
            let reply = loop {
                if let Some(reply) = worker.poll() {
                    break reply;
                }
                std::thread::sleep(Duration::from_micros(200));
            };
            assert_eq!(reply.id, job);
            let multicad_core::exchange::PrintResult::Preview(Ok(preview)) = reply.result else {
                panic!("plate preview failed");
            };
            assert!(preview.can_write());
            rounds.push(started.elapsed());
        }
        let (preview, stages) = Exporter::default()
            .print_plate_preview_timed(&document, 0, &plate, &actor(), None)
            .unwrap();
        assert!(preview.can_write());
        eprintln!(
            "sample {sample}: worker first {:?} repeat {:?}; stages {stages:?}",
            rounds[0], rounds[1]
        );
    }
}
