use super::*;
use crate::step::degree11_phillips::{self as exchange, Error, Options};
use crate::step_read::phillips_cutter::physical::body::{
    derived_thread::actual_m3, read_degree11_phillips_material,
};
fn read(text: &str) -> OwnedDegree11Material {
    read_degree11_phillips_material(text, Default::default(), None)
        .unwrap()
        .into_material()
}
fn limited<T>(value: std::result::Result<T, Error>) -> bool {
    matches!(
        value,
        Err(Error::Limit)
            | Err(Error::Reader(crate::step_read::PhillipsDegree11ReadError {
                cause: PhillipsReadError::Limit,
                ..
            }))
    )
}
#[test]
fn actual_retained_d11_singleton_proposal_earns_distinct_decoded_material_both_hands() {
    for left in [false, true] {
        let original = actual_m3(true, left);
        let source = read(&original);
        let qualified = exchange::export_retained(&source, Options::default(), None).unwrap();
        assert_eq!(qualified.as_str(), original);
        let decoded = qualified.qualification().material();
        assert_ne!(
            source.thread_faces()[0].source().poles().as_ptr(),
            decoded.thread_faces()[0].source().poles().as_ptr()
        );
        let p = qualified.qualification().verification();
        eprintln!(
            "D11_RETAINED_STEP left={left} bytes={} aggregate={} source={} reader={} decoded={} fresh={} encoding={} known_bytes={}",
            qualified.as_str().len(),
            p.aggregate_reserved_work,
            p.original_source_prefix_work,
            p.reader_reserved_work,
            p.actual_decoded_material_proof_work,
            p.fresh_correspondence_work,
            p.proposal_work,
            p.known_resident_peak_bytes
        );
        assert_eq!(p.original_source_prefix_work, source.work_used());
        assert_eq!(p.actual_decoded_material_proof_work, decoded.work_used());
        assert_eq!(
            p.aggregate_reserved_work,
            p.original_source_prefix_work
                + p.reader_reserved_work
                + p.fresh_correspondence_work
                + p.proposal_work
        );
        assert!(p.maximum_world_distance_mm <= 1e-7 && p.maximum_world_chart_derivative_mm <= 1e-7);
        assert!(p.independent_repaired_volume_discrepancy_mm3 <= 1e-5);
        let actual =
            exchange::qualify_actual_file(&source, qualified.as_str(), Options::default(), None)
                .unwrap();
        assert_eq!(
            actual.material().world_volume_mm3(),
            decoded.world_volume_mm3()
        );
        assert_eq!(actual.verification().proposal_work, 0);
        let detached = actual.into_material();
        drop(qualified);
        drop(source);
        drop(original);
        assert!(detached.local_volume_mm3()[0] > 0.);
        assert!(
            detached
                .thread_faces()
                .iter()
                .all(|f| !f.source().poles().is_empty())
        );
        assert!(detached.original_input().starts_with("ISO-10303-21;"));
    }
}
#[test]
fn retained_d11_proposal_limits_actual_stale_file_and_late_cancel_never_publish() {
    let text = actual_m3(true, false);
    let source = read(&text);
    let options = Options::default();
    let calls = Cell::new(0usize);
    let count = || {
        calls.set(calls.get() + 1);
        false
    };
    let candidate = exchange::export_retained(&source, options, Some(&count)).unwrap();
    let earned = candidate.qualification().verification();
    for bad in [
        Options {
            maximum_output_bytes: text.len() - 1,
            ..options
        },
        Options {
            maximum_proposal_work: earned.proposal_work - 1,
            ..options
        },
        Options {
            maximum_aggregate_work: earned.aggregate_reserved_work - 1,
            ..options
        },
        Options {
            maximum_known_resident_bytes: earned.known_resident_peak_bytes - 1,
            ..options
        },
    ] {
        assert!(limited(exchange::export_retained(&source, bad, None)));
    }
    let exact = exchange::export_retained(
        &source,
        Options {
            maximum_output_bytes: text.len(),
            maximum_proposal_work: earned.proposal_work,
            maximum_aggregate_work: earned.aggregate_reserved_work,
            maximum_known_resident_bytes: earned.known_resident_peak_bytes,
            ..options
        },
        None,
    )
    .unwrap();
    assert_eq!(exact.as_str(), text);
    assert_eq!(
        exact.qualification().verification().aggregate_reserved_work,
        earned.aggregate_reserved_work
    );
    let mut bad = options;
    bad.correspondence.maximum_work = 512_000_001;
    assert!(matches!(
        exchange::export_retained(&source, bad, None),
        Err(Error::Invalid)
    ));
    let stale = actual_m3(false, false);
    assert!(matches!(
        exchange::qualify_actual_file(&source, &stale, options, None),
        Err(Error::Unsupported)
    ));
    assert!(exchange::qualify_actual_file(&source, "ISO-10303-21;BROKEN", options, None).is_err());
    assert!(matches!(
        exchange::export_retained(&source, options, Some(&|| true)),
        Err(Error::Cancelled)
    ));
    let late = Cell::new(0usize);
    let cancel = || {
        late.set(late.get() + 1);
        late.get() >= calls.get()
    };
    assert!(matches!(
        exchange::export_retained(&source, options, Some(&cancel)),
        Err(Error::Cancelled)
    ));
    let retry = exchange::export_retained(&source, options, None).unwrap();
    assert_eq!(retry.as_str(), candidate.as_str());
    assert_eq!(source.original_input(), text);
}
