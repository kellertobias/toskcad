use super::*;
use crate::step_read::owned_insert::DocumentOptions;
fn upper(parsed: &Parsed<'_, '_>) -> u32 {
    parsed
        .faces
        .iter()
        .filter_map(|f| match f.support {
            Support::Plane { height } => Some((parsed.points[height].coordinates[2].value(), f.id)),
            _ => None,
        })
        .max_by(|a, b| a.0.total_cmp(&b.0))
        .unwrap()
        .1
}
fn signature(parsed: &Parsed<'_, '_>) -> u64 {
    use std::hash::{Hash, Hasher};
    let mut h = std::collections::hash_map::DefaultHasher::new();
    for p in &parsed.points {
        p.id.hash(&mut h);
        for x in p.coordinates {
            x.spelling().hash(&mut h);
            x.value().to_bits().hash(&mut h);
            x.bounds().map(f64::to_bits).hash(&mut h);
        }
    }
    for t in &parsed.tensors {
        t.id.hash(&mut h);
        t.spans.hash(&mut h);
        t.poles.hash(&mut h);
    }
    for e in &parsed.edges {
        e.id.hash(&mut h);
        e.curve.hash(&mut h);
        e.vertices.hash(&mut h);
    }
    h.finish()
}
#[test]
fn canonical_clipped_d11_head_requalifies_original_thread_caps_and_ids() {
    for quarter in [false, true] {
        let text = super::super::tests::fixture(quarter);
        let document = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
        let parsed = Parsed::parse(
            &document,
            document.represented_members()[0].item_entity(),
            Default::default(),
            None,
        )
        .unwrap();
        let before = signature(&parsed);
        let cap = upper(&parsed);
        let input = HeadInput::regular_hex(&parsed, cap, 18., 3.).unwrap();
        let proof = certify(&parsed, input, 1e-7, 1e-6, HeadOptions::default(), None).unwrap();
        assert_eq!(proof.head().auxiliary_cap_entity, cap);
        assert!(std::ptr::eq(
            proof.capped_thread().projection().parsed(),
            &parsed
        ));
        assert_eq!(
            proof.capped_thread().projection().parsed().tensors().len(),
            8
        );
        assert!(proof.capped_thread().volume_mm3()[0] > 0.);
        assert!(proof.minimum_wall_mm() > 1e-6);
        assert!(proof.join_sewing_mm() <= 1e-7 / 8.);
        assert!(proof.head_stations_mm()[1] > proof.head_stations_mm()[0]);
        assert_eq!(signature(&parsed), before);
        let receipt = proof.receipt();
        assert_eq!(
            receipt.aggregate_work,
            receipt.inherited_document_work
                + receipt.material_work_with_parsed_prefix
                + receipt.new_head_work
        );
        assert!(receipt.aggregate_work <= 96_000_000);
        assert!(receipt.accounted_peak_bytes <= 64 * 1024 * 1024);
        eprintln!(
            "CLIPPED_D11_HEAD quarter={quarter} sourcepoints={} work={} peak={} collar={} wall={}",
            parsed.points.len(),
            receipt.aggregate_work,
            receipt.accounted_peak_bytes,
            proof.join_sewing_mm(),
            proof.minimum_wall_mm()
        );
    }
}

#[test]
fn actual_clipped_head_refuses_wrong_rim_star_intersection_gap_and_bad_supports() {
    let text = super::super::tests::fixture(false);
    let document = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
    let parsed = Parsed::parse(
        &document,
        document.represented_members()[0].item_entity(),
        Default::default(),
        None,
    )
    .unwrap();
    let cap = upper(&parsed);
    let before = signature(&parsed);
    let input = HeadInput::regular_hex(&parsed, cap, 18., 3.).unwrap();
    let mut cases = Vec::new();
    let mut bad = input.clone();
    bad.faces[0].inner[0].1 = !bad.faces[0].inner[0].1;
    cases.push(("rim sense", bad));
    let mut bad = input.clone();
    bad.faces[0].inner[0].0 = u32::MAX;
    cases.push(("foreign rim", bad));
    let mut bad = input.clone();
    bad.edges[0] = bad.edges[1];
    cases.push(("broken fan", bad));
    let mut bad = input.clone();
    bad.faces[2].plane.normal = [0.; 3];
    cases.push(("zero normal", bad));
    let mut bad = input.clone();
    bad.faces[2].plane.origin[0] = f64::NAN;
    cases.push(("NaN plane", bad));
    let mut bad = input.clone();
    bad.faces[2].plane.reference = bad.faces[2].plane.normal;
    cases.push(("parallel frame", bad));
    let mut bad = input.clone();
    bad.vertices[6][0] += 0.001;
    cases.push(("nonvertical head", bad));
    let mut bad = input.clone();
    for p in &mut bad.vertices {
        p[2] += 1e-5;
    }
    for face in &mut bad.faces {
        face.plane.origin[2] += 1e-5;
    }
    cases.push(("detached contact", bad));
    let mut bad = input.clone();
    let angles = [0_f64, 130., 260., 30., 160., 290.].map(f64::to_radians);
    for (i, p) in bad.vertices.iter_mut().enumerate() {
        p[0] = 20. * angles[i % 6].cos();
        p[1] = 20. * angles[i % 6].sin();
    }
    for i in 0..6 {
        let j = (i + 1) % 6;
        let a = bad.vertices[i];
        let b = bad.vertices[j];
        let dx = b[0] - a[0];
        let dy = b[1] - a[1];
        bad.faces[i + 2].plane = Plane {
            origin: a,
            normal: [dy, -dx, 0.],
            reference: [dx, dy, 0.],
        };
    }
    // All local turns are positive, and the central thread fits inside each
    // local wedge. Only the complete all-edge convexity check rejects this star.
    for i in 0..6 {
        let a = bad.vertices[i];
        let b = bad.vertices[(i + 1) % 6];
        let c = bad.vertices[(i + 2) % 6];
        assert!((b[0] - a[0]) * (c[1] - b[1]) - (b[1] - a[1]) * (c[0] - b[0]) > 0.);
    }
    cases.push(("self-intersecting star", bad));
    cases.push((
        "insufficient full hull clearance",
        HeadInput::regular_hex(&parsed, cap, 9., 3.).unwrap(),
    ));
    cases.push((
        "thin head",
        HeadInput::regular_hex(&parsed, cap, 18., 1e-8).unwrap(),
    ));
    let lower = parsed
        .faces
        .iter()
        .find(|f| matches!(f.support, Support::Plane { .. }) && f.id != cap)
        .unwrap()
        .id;
    cases.push((
        "wrong auxiliary cap",
        HeadInput::regular_hex(&parsed, lower, 18., 3.).unwrap(),
    ));
    for (label, bad) in cases {
        assert!(
            certify(&parsed, bad, 1e-7, 1e-6, HeadOptions::default(), None).is_err(),
            "{label}"
        );
        assert_eq!(signature(&parsed), before);
    }
    assert!(certify(&parsed, input, 1e-7, 1e-6, HeadOptions::default(), None).is_ok());
}

#[test]
fn actual_raw_cap_sense_is_requalified_and_cannot_inherit_a_parsed_material_claim() {
    let text = super::super::tests::fixture(false);
    let document = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
    let parsed = Parsed::parse(
        &document,
        document.represented_members()[0].item_entity(),
        Default::default(),
        None,
    )
    .unwrap();
    let cap = upper(&parsed);
    let record = document.records().get(&cap).unwrap()[0].source.unwrap();
    let wrong = record.replacen(".T.", ".F.", 1);
    assert_ne!(wrong, record);
    let changed = text.replacen(record, &wrong, 1);
    assert_ne!(changed, text);
    let invalid_document =
        CompoundDocument::parse(&changed, DocumentOptions::default(), None).unwrap();
    let actual = Parsed::parse(
        &invalid_document,
        invalid_document.represented_members()[0].item_entity(),
        Default::default(),
        None,
    )
    .unwrap();
    let head = HeadInput::regular_hex(&actual, cap, 18., 3.).unwrap();
    assert!(matches!(
        certify(&actual, head, 1e-7, 1e-6, HeadOptions::default(), None),
        Err(Error::Unsupported)
    ));
    assert!(
        certify(
            &parsed,
            HeadInput::regular_hex(&parsed, cap, 18., 3.).unwrap(),
            1e-7,
            1e-6,
            HeadOptions::default(),
            None
        )
        .is_ok()
    );
}

#[test]
fn clipped_head_counts_real_source_and_head_memory_work_cancel_and_retry() {
    let text = super::super::tests::fixture(false);
    let document = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
    let parsed = Parsed::parse(
        &document,
        document.represented_members()[0].item_entity(),
        Default::default(),
        None,
    )
    .unwrap();
    let cap = upper(&parsed);
    let before = signature(&parsed);
    let input = || HeadInput::regular_hex(&parsed, cap, 18., 3.).unwrap();
    let calls = std::cell::Cell::new(0usize);
    let observe = || {
        calls.set(calls.get() + 1);
        false
    };
    let positive = certify(
        &parsed,
        input(),
        1e-7,
        1e-6,
        HeadOptions::default(),
        Some(&observe),
    )
    .unwrap();
    let receipt = positive.receipt();
    let total = calls.get();
    drop(positive);
    assert!(receipt.borrowed_document_bytes >= text.len());
    assert!(receipt.head_input_bytes > 0);
    for options in [
        HeadOptions {
            maximum_work: receipt.aggregate_work - 1,
            ..Default::default()
        },
        HeadOptions {
            maximum_bytes: receipt.accounted_peak_bytes - 1,
            ..Default::default()
        },
        HeadOptions {
            additional_resident_bytes: usize::MAX,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            certify(&parsed, input(), 1e-7, 1e-6, options, None),
            Err(Error::Limit)
        ));
    }
    let exact = HeadOptions {
        maximum_work: receipt.aggregate_work,
        maximum_bytes: receipt.accounted_peak_bytes,
        ..Default::default()
    };
    assert!(certify(&parsed, input(), 1e-7, 1e-6, exact, None).is_ok());
    for stop in [1, total / 2, total] {
        calls.set(0);
        let cancel = || {
            calls.set(calls.get() + 1);
            calls.get() >= stop
        };
        assert!(matches!(
            certify(
                &parsed,
                input(),
                1e-7,
                1e-6,
                HeadOptions::default(),
                Some(&cancel)
            ),
            Err(Error::Cancelled)
        ));
        assert_eq!(signature(&parsed), before);
    }
    let mut capacity = input();
    capacity.vertices.reserve_exact(5000);
    let accounted = storage(&capacity).unwrap();
    assert!(accounted > receipt.head_input_bytes);
    assert!(matches!(
        certify(
            &parsed,
            capacity,
            1e-7,
            1e-6,
            HeadOptions {
                maximum_bytes: accounted - 1,
                ..Default::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
    let head = input();
    let bytes = storage(&head).unwrap();
    let resident = 4096;
    let proof = certify(
        &parsed,
        head,
        1e-7,
        1e-6,
        HeadOptions {
            additional_resident_bytes: resident,
            ..Default::default()
        },
        None,
    )
    .unwrap();
    assert_eq!(proof.receipt().head_input_bytes, bytes);
    assert_eq!(
        proof.receipt().accounted_peak_bytes,
        receipt.accounted_peak_bytes + resident
    );
    drop(proof);
    for (tolerance, minimum) in [
        (0., 1e-6),
        (f64::NAN, 1e-6),
        (1e-7, 0.),
        (1e-7, f64::INFINITY),
    ] {
        assert!(matches!(
            certify(
                &parsed,
                input(),
                tolerance,
                minimum,
                HeadOptions::default(),
                None
            ),
            Err(Error::Invalid)
        ));
    }
    assert!(certify(&parsed, input(), 1e-7, 1e-6, HeadOptions::default(), None).is_ok());
    assert_eq!(signature(&parsed), before);
}

#[test]
fn clipped_head_public_options_and_actual_clearance_thresholds_are_enforced() {
    let text = super::super::tests::fixture(true);
    let document = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
    let parsed = Parsed::parse(
        &document,
        document.represented_members()[0].item_entity(),
        Default::default(),
        None,
    )
    .unwrap();
    let cap = upper(&parsed);
    let input = || HeadInput::regular_hex(&parsed, cap, 18., 3.).unwrap();
    for options in [
        HeadOptions {
            maximum_work: 0,
            ..Default::default()
        },
        HeadOptions {
            maximum_work: MAX_WORK + 1,
            ..Default::default()
        },
        HeadOptions {
            maximum_bytes: 0,
            ..Default::default()
        },
        HeadOptions {
            maximum_bytes: MAX_BYTES + 1,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            certify(&parsed, input(), 1e-7, 1e-6, options, None),
            Err(Error::Invalid)
        ));
    }
    for (across, height) in [(0., 3.), (f64::NAN, 3.), (18., 0.), (18., f64::INFINITY)] {
        assert!(matches!(
            HeadInput::regular_hex(&parsed, cap, across, height),
            Err(Error::Invalid)
        ));
    }
    assert!(HeadInput::regular_hex(&parsed, u32::MAX, 18., 3.).is_err());
    let tensor = parsed
        .faces
        .iter()
        .find(|f| matches!(f.support, Support::Tensor(_)))
        .unwrap()
        .id;
    assert!(matches!(
        HeadInput::regular_hex(&parsed, tensor, 18., 3.),
        Err(Error::Unsupported)
    ));
    // This actual quarter-inch retained hull clears the eighteen-mm hex by
    // about1.8347mm: an1.8mm requirement passes and1.9mm must refuse.
    assert!(certify(&parsed, input(), 1e-7, 1.8, HeadOptions::default(), None).is_ok());
    assert!(matches!(
        certify(&parsed, input(), 1e-7, 1.9, HeadOptions::default(), None),
        Err(Error::Unsupported)
    ));
}
