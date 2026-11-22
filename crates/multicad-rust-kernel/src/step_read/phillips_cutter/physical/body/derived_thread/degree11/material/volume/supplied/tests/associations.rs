use super::*;
use crate::step_read::{
    literal_member_binding::{Token, token},
    part21,
};

// Renumber the actual complete STEP graph, including declarations, topology,
// PCURVEs and metadata. The validated token cursor leaves literals untouched.
fn independent_ids(input: &str) -> (String, u32) {
    let original = part21::parse(input).unwrap();
    let offset = original.keys().next_back().unwrap().checked_add(1).unwrap();
    let mut rest = input;
    let mut changes = Vec::new();
    while let Some(value) = token(&mut rest) {
        if let Token::Reference(id) = value {
            assert!(original.contains_key(&id));
            let end = input.len() - rest.len();
            let spelling = format!("#{id}");
            let start = end - spelling.len();
            assert_eq!(&input[start..end], spelling);
            changes.push((start..end, format!("#{}", id.checked_add(offset).unwrap())));
        }
    }
    let mut text = input.to_owned();
    for (range, replacement) in changes.into_iter().rev() {
        text.replace_range(range, &replacement);
    }
    let rebased = part21::parse(&text).unwrap();
    assert_eq!(original.len(), rebased.len());
    assert!(original.keys().all(|id| !rebased.contains_key(id)));
    let (mut a, mut b) = (input, text.as_str());
    loop {
        match (token(&mut a), token(&mut b)) {
            (None, None) => break,
            (Some(Token::Reference(a)), Some(Token::Reference(b))) => {
                assert_eq!(a.checked_add(offset), Some(b))
            }
            (Some(Token::Literal(a)), Some(Token::Literal(b))) => assert_eq!(a, b),
            _ => panic!("only complete actual entity IDs may change"),
        }
    }
    (text, offset)
}
fn extent(material: &Degree11ClosedMaterial<'_, '_, '_>) -> I {
    let mut radius = I::point(0.);
    let mut z = I::point(0.);
    for i in 0..material.head().segment_count() {
        let b = material.head().certified_interval_bounds(i, [0., 1.]);
        radius = maximum([radius, abs(interval(b.position_mm[0]))]);
        z = maximum([z, abs(interval(b.position_mm[1]))]);
    }
    maximum([
        chart::norm([radius, radius, z]),
        chart::norm(
            material
                .thread()
                .local_enclosure_bounds_mm()
                .map(|x| I::point(x[0].abs().max(x[1].abs()))),
        ),
    ])
}
fn fresh_work() -> Work<'static> {
    Work {
        remaining: WORK_LIMIT,
        interrupt: None,
        cancelled: std::cell::Cell::new(false),
    }
}
fn signature(value: &ThreadBounds) -> Vec<[f64; 2]> {
    std::iter::once(value.delta)
        .chain(value.derivative)
        .chain(value.speed)
        .chain(
            value
                .rims
                .iter()
                .flat_map(|r| [r.delta, r.derivative, r.speed]),
        )
        .map(|x| [x.lo, x.hi])
        .collect()
}
#[test]
fn actual_independently_earned_material_and_receipt_owners_reject_mixed_face_ids() {
    let (author, frame) = actual_m3_parts(true, false);
    let original =
        crate::step::write_phillips_screw(&author, frame[0], frame[1], frame[2], 1e-7).unwrap();
    let (rebased, offset) = independent_ids(&original);
    let path_a = std::env::temp_dir().join(format!("association-a-{}.step", std::process::id()));
    let path_b = std::env::temp_dir().join(format!("association-b-{}.step", std::process::id()));
    std::fs::write(&path_a, &original).unwrap();
    std::fs::write(&path_b, &rebased).unwrap();
    drop(original);
    drop(rebased);
    let a = std::fs::read_to_string(&path_a).unwrap();
    let b = std::fs::read_to_string(&path_b).unwrap();
    std::fs::remove_file(&path_a).unwrap();
    std::fs::remove_file(&path_b).unwrap();
    with_head(&a, &mut |head_a| {
        let occupied_a = head_a
            .degree11_thread_carrier(Default::default(), None)
            .unwrap()
            .certify_occupied(Default::default(), None)
            .unwrap();
        let material_a = head_a
            .closed_degree11_material(&occupied_a, Default::default(), None)
            .unwrap();
        let volume_a = material_a.volume_bounds(Default::default(), None).unwrap();
        let matched_a = volume_a
            .match_supplied_author(&author, frame, Default::default(), None)
            .unwrap();
        let world_a = chart::compare(frame, material_a.chart()).unwrap();
        let extent_a = extent(&material_a);
        let baseline = thread_bounds(
            &material_a,
            matched_a.supplied_thread(),
            &world_a,
            extent_a,
            &mut fresh_work(),
        )
        .unwrap();
        with_head(&b, &mut |head_b| {
            let occupied_b = head_b
                .degree11_thread_carrier(Default::default(), None)
                .unwrap()
                .certify_occupied(Default::default(), None)
                .unwrap();
            let material_b = head_b
                .closed_degree11_material(&occupied_b, Default::default(), None)
                .unwrap();
            let volume_b = material_b.volume_bounds(Default::default(), None).unwrap();
            let matched_b = volume_b
                .match_supplied_author(&author, frame, Default::default(), None)
                .unwrap();
            let world_b = chart::compare(frame, material_b.chart()).unwrap();
            let extent_b = extent(&material_b);
            assert_eq!(material_a.thread().carrier().faces().len(), 4);
            assert_eq!(matched_b.supplied_thread().transfers().len(), 4);
            for (left, right) in matched_a
                .supplied_thread()
                .transfers()
                .iter()
                .zip(matched_b.supplied_thread().transfers())
            {
                assert_eq!(
                    left.face_entity.checked_add(offset),
                    Some(right.face_entity)
                );
                assert_eq!(
                    left.surface_entity.checked_add(offset),
                    Some(right.surface_entity)
                );
                assert_eq!(
                    left.original_rail_entities.map(|id| id + offset),
                    right.original_rail_entities
                );
            }
            assert!(matches!(
                thread_bounds(
                    &material_a,
                    matched_b.supplied_thread(),
                    &world_a,
                    extent_a,
                    &mut fresh_work()
                ),
                Err(PhillipsReadError::Unsupported)
            ));
            let own_a = thread_bounds(
                &material_a,
                matched_a.supplied_thread(),
                &world_a,
                extent_a,
                &mut fresh_work(),
            )
            .unwrap();
            assert_eq!(signature(&own_a), signature(&baseline));
            let own_b = thread_bounds(
                &material_b,
                matched_b.supplied_thread(),
                &world_b,
                extent_b,
                &mut fresh_work(),
            )
            .unwrap();
            assert_eq!(signature(&own_b), signature(&baseline));
            assert_eq!(volume_a.local_volume_mm3(), volume_b.local_volume_mm3());
            assert_eq!(matched_a.supplied_thread().body().literal_input(), a);
            assert_eq!(matched_b.supplied_thread().body().literal_input(), b);
            eprintln!(
                "ACTUAL_ASSOCIATION two_complete_material_matchestrue offset{offset} mixed_faces_refusedtrue own_retriestrue"
            );
        });
    });
}

#[test]
fn actual_multistart_material_and_single_start_receipt_refuse_count_mismatch() {
    let (single, frame) = actual_m3_parts(true, false);
    let mut spec = single.thread().spec();
    spec.starts = 2;
    let multi = crate::phillips_drive::PhillipsScrew::new(
        crate::helical_thread::Thread::new(spec).unwrap(),
        *single.profile(),
        single.spec(),
        single.cavity().clone(),
        100_000,
        &|| false,
    )
    .unwrap();
    let a = crate::step::write_phillips_screw(&single, frame[0], frame[1], frame[2], 1e-7).unwrap();
    let b = crate::step::write_phillips_screw(&multi, frame[0], frame[1], frame[2], 1e-7).unwrap();
    let path = std::env::temp_dir().join(format!(
        "association-multistart-{}.step",
        std::process::id()
    ));
    std::fs::write(&path, &b).unwrap();
    drop(b);
    let b = std::fs::read_to_string(&path).unwrap();
    std::fs::remove_file(&path).unwrap();
    with_head(&a, &mut |single_head| {
        let occupied = single_head
            .degree11_thread_carrier(Default::default(), None)
            .unwrap()
            .certify_occupied(Default::default(), None)
            .unwrap();
        let material = single_head
            .closed_degree11_material(&occupied, Default::default(), None)
            .unwrap();
        let volume = material.volume_bounds(Default::default(), None).unwrap();
        let supplied = volume
            .match_supplied_author(&single, frame, Default::default(), None)
            .unwrap();
        let world = chart::compare(frame, material.chart()).unwrap();
        let baseline = thread_bounds(
            &material,
            supplied.supplied_thread(),
            &world,
            extent(&material),
            &mut fresh_work(),
        )
        .unwrap();
        with_head(&b, &mut |multi_head| {
            let occupied = multi_head
                .degree11_thread_carrier(Default::default(), None)
                .unwrap()
                .certify_occupied(Default::default(), None)
                .unwrap();
            let multi_material = multi_head
                .closed_degree11_material(&occupied, Default::default(), None)
                .unwrap();
            let multi_volume = multi_material
                .volume_bounds(Default::default(), None)
                .unwrap();
            assert_eq!(multi_material.thread().carrier().faces().len(), 8);
            assert_eq!(supplied.supplied_thread().transfers().len(), 4);
            let multi_world = chart::compare(frame, multi_material.chart()).unwrap();
            assert!(matches!(
                thread_bounds(
                    &multi_material,
                    supplied.supplied_thread(),
                    &multi_world,
                    extent(&multi_material),
                    &mut fresh_work()
                ),
                Err(PhillipsReadError::Unsupported)
            ));
            let retry = thread_bounds(
                &material,
                supplied.supplied_thread(),
                &world,
                extent(&material),
                &mut fresh_work(),
            )
            .unwrap();
            assert_eq!(signature(&retry), signature(&baseline));
            let retry_volume = multi_material
                .volume_bounds(Default::default(), None)
                .unwrap();
            assert_eq!(
                multi_volume.local_volume_mm3(),
                retry_volume.local_volume_mm3()
            );
            assert_eq!(multi_material.body().literal_input(), b);
            eprintln!(
                "ACTUAL_ASSOCIATION multistart_material_faces8 single_receipt_faces4 count_refusedtrue real_material_volume_retrytrue"
            );
        });
    });
}
