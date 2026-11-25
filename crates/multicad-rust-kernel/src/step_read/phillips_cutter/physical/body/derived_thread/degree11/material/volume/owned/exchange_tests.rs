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
fn authored_ok(p: exchange::AuthoredVerification) {
    assert!(p.maximum_world_distance_mm <= 1e-7);
    assert!(p.maximum_world_chart_derivative_mm <= 1e-7);
    assert!(
        p.independent_volume_difference_mm3
            .into_iter()
            .all(|v| v.is_finite() && v.abs() <= 1e-5)
    );
    assert!(p.decoded_repaired_volume_bounds_mm3[0] > 0.);
    assert!(p.decoded_repaired_volume_bounds_mm3[1] >= p.decoded_repaired_volume_bounds_mm3[0]);
    assert!(
        p.aggregate_reserved_work <= exchange::AuthoredOptions::default().maximum_aggregate_work
    );
    assert!(p.fresh_correspondence_work > 0);
}
#[test]
fn authored_fitted_export_qualifies_its_own_bytes_against_the_author_both_hands() {
    use crate::step_read::phillips_cutter::physical::body::derived_thread::actual_m3_parts;
    for left in [false, true] {
        let (author, frame) = actual_m3_parts(true, left);
        let candidate =
            exchange::export_authored(&author, frame, 1e-7, Default::default(), None).unwrap();
        // The proposal is exactly the deterministic authored writer output.
        let [origin, u, v] = frame;
        assert_eq!(
            candidate.as_str(),
            crate::step::write_phillips_screw(&author, origin, u, v, 1e-7).unwrap()
        );
        let p = candidate.qualification().verification();
        authored_ok(p);
        assert!(p.proposal_work >= 2 * candidate.as_str().len());
        // Author volume and decoded repaired volume agree within the bound.
        let author_volume = author.source_volume_bounds().unwrap().volume_mm3();
        let decoded = p.decoded_repaired_volume_bounds_mm3;
        assert!(author_volume[0] - decoded[1] <= 1e-5 && decoded[0] - author_volume[1] <= 1e-5);
        // Fresh qualification of the same actual bytes repeats the proof.
        let again = exchange::qualify_authored_file(
            &author,
            frame,
            candidate.as_str(),
            Default::default(),
            None,
        )
        .unwrap();
        authored_ok(again.verification());
        assert_eq!(again.verification().proposal_work, 0);
    }
}
#[test]
fn authored_qualification_refuses_other_hands_moved_frames_and_bad_options() {
    use crate::step_read::phillips_cutter::physical::body::derived_thread::actual_m3_parts;
    let (author, frame) = actual_m3_parts(true, false);
    let text = actual_m3(true, false);
    // A valid fitted file of the other hand is different geometry.
    let (other, _) = actual_m3_parts(true, true);
    assert!(
        exchange::qualify_authored_file(&other, frame, &text, Default::default(), None).is_err()
    );
    // The nominal (flat-root) author is not this fitted file either.
    let (flat, _) = actual_m3_parts(false, false);
    assert!(
        exchange::qualify_authored_file(&flat, frame, &text, Default::default(), None).is_err()
    );
    // One micrometre of occurrence shift exceeds the 1e-7 mm correspondence.
    let mut moved = frame;
    moved[0][0] += 1e-3;
    assert!(
        exchange::qualify_authored_file(&author, moved, &text, Default::default(), None).is_err()
    );
    // Non-orthonormal frames, invalid tolerances and oversized budgets refuse
    // before any read.
    let mut skew = frame;
    skew[1] = [1., 1e-3, 0.];
    assert!(matches!(
        exchange::qualify_authored_file(&author, skew, &text, Default::default(), None),
        Err(Error::Invalid)
    ));
    assert!(matches!(
        exchange::export_authored(&author, frame, 1e-6, Default::default(), None),
        Err(Error::Invalid)
    ));
    let over = exchange::AuthoredOptions {
        maximum_aggregate_work: exchange::AuthoredOptions::default().maximum_aggregate_work + 1,
        ..Default::default()
    };
    assert!(matches!(
        exchange::qualify_authored_file(&author, frame, &text, over, None),
        Err(Error::Invalid)
    ));
    let tight = exchange::AuthoredOptions {
        maximum_aggregate_work: 1_000_000,
        ..Default::default()
    };
    assert!(limited(exchange::qualify_authored_file(
        &author, frame, &text, tight, None
    )));
    let tiny = exchange::AuthoredOptions {
        maximum_output_bytes: text.len() - 1,
        ..Default::default()
    };
    assert!(limited(exchange::export_authored(
        &author, frame, 1e-7, tiny, None
    )));
    assert!(matches!(
        exchange::qualify_authored_file(&author, frame, &text, Default::default(), Some(&|| true)),
        Err(Error::Cancelled)
    ));
}
#[test]
fn authored_candidate_replays_into_a_writer_whose_actual_bytes_requalify() {
    use crate::step_read::phillips_cutter::physical::body::derived_thread::actual_m3_parts;
    let (author, frame) = actual_m3_parts(true, true);
    let candidate =
        exchange::export_authored(&author, frame, 1e-7, Default::default(), None).unwrap();
    let mut writer = crate::step::Writer::new();
    let receipt = writer
        .append_qualified_authored_degree11_phillips(&candidate, Default::default(), None)
        .unwrap();
    assert!(receipt.matched_entities > 0);
    assert_eq!(receipt.candidate, candidate.qualification().verification());
    let text = writer.finish().unwrap();
    let replayed =
        exchange::qualify_authored_file(&author, frame, &text, Default::default(), None).unwrap();
    authored_ok(replayed.verification());
    // A cancelled append leaves the writer unchanged.
    let mut cancelled = crate::step::Writer::new();
    assert!(matches!(
        cancelled.append_qualified_authored_degree11_phillips(
            &candidate,
            Default::default(),
            Some(&|| true)
        ),
        Err(Error::Cancelled)
    ));
    // Same writer, appended afterwards: byte-identical to the clean append.
    cancelled
        .append_qualified_authored_degree11_phillips(&candidate, Default::default(), None)
        .unwrap();
    assert_eq!(cancelled.finish().unwrap(), text);
}
#[test]
fn authored_routing_admits_contained_countersunk_recesses_and_their_export_qualifies() {
    use crate::step_read::phillips_cutter::physical::body::derived_thread::m3_parts;
    let (author, frame) = m3_parts(true, false, 1.5, 1., true);
    assert!(matches!(
        author.spec().head,
        Some(crate::head_profile::Head::Countersunk { .. })
    ));
    // The routing predicate is not a certificate: when it admits, the actual
    // bytes must still qualify against the author.
    if exchange::supports_author(&author) {
        let candidate =
            exchange::export_authored(&author, frame, 1e-7, Default::default(), None).unwrap();
        authored_ok(candidate.qualification().verification());
    }
    // A Pan author always routes.
    let (pan, _) = m3_parts(true, false, 1.5, 1., false);
    assert!(exchange::supports_author(&pan));
}
#[test]
fn retained_d11_moved_occurrence_replays_literal_local_records_through_one_mapped_item() {
    let original = actual_m3(true, false);
    let source = read(&original);
    // A proper rigid occurrence with irrational axes.
    let (s, c) = (0.6_f64, 0.8_f64);
    let occurrence = [[10., -20., 30.], [c, s, 0.], [-s, c, 0.], [0., 0., 1.]];
    let placed =
        exchange::export_retained_placed(&source, occurrence, Options::default(), None).unwrap();
    assert_eq!(placed.as_str().matches("MAPPED_ITEM(").count(), 1);
    let p = placed.qualification().verification();
    assert!(p.maximum_world_distance_mm <= 1e-7 && p.maximum_world_chart_derivative_mm <= 1e-7);
    assert!(p.independent_repaired_volume_discrepancy_mm3 <= 1e-5);
    let fresh = exchange::qualify_placed_actual_file(
        &source,
        placed.as_str(),
        occurrence,
        Options::default(),
        None,
    )
    .unwrap();
    assert_eq!(
        fresh.material().world_volume_mm3(),
        placed.qualification().material().world_volume_mm3()
    );
    // The identity export is the original flat file; the mapped and flat
    // contracts do not accept each other's bytes.
    assert!(
        exchange::qualify_actual_file(&source, placed.as_str(), Options::default(), None).is_err()
    );
    let identity = [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
    assert_eq!(
        exchange::qualify_placed_actual_file(
            &source,
            &original,
            identity,
            Options::default(),
            None
        )
        .err(),
        Some(Error::Unsupported)
    );
    // A one-micrometre different occurrence, a mirrored or skewed occurrence,
    // another source and cancellation are refused.
    let mut moved = occurrence;
    moved[0][0] += 1e-3;
    assert!(
        exchange::qualify_placed_actual_file(
            &source,
            placed.as_str(),
            moved,
            Options::default(),
            None
        )
        .is_err()
    );
    let mirrored = [occurrence[0], occurrence[1], occurrence[2], [0., 0., -1.]];
    let skewed = [occurrence[0], [1., 0., 0.], [0.5, 1., 0.], [0., 0., 1.]];
    for bad in [mirrored, skewed] {
        assert_eq!(
            exchange::export_retained_placed(&source, bad, Options::default(), None).err(),
            Some(Error::Invalid)
        );
    }
    let other = read(&actual_m3(true, true));
    assert!(
        exchange::qualify_placed_actual_file(
            &other,
            placed.as_str(),
            occurrence,
            Options::default(),
            None
        )
        .is_err()
    );
    assert_eq!(
        exchange::export_retained_placed(&source, occurrence, Options::default(), Some(&|| true))
            .err(),
        Some(Error::Cancelled)
    );
}
