use super::*;
use crate::step::insert::{CurveOptions, coupled::CoupledPlan};
fn change_record(text: &str, id: u32, change: impl FnOnce(&str) -> String) -> String {
    let prefix = format!("#{id}=");
    let line = text.lines().find(|s| s.starts_with(&prefix)).unwrap();
    let changed = change(line);
    assert_ne!(line, changed);
    text.replacen(line, &changed, 1)
}
#[test]
fn actual_shared_parallel_row_fold_is_not_hidden_by_the_strip_shortcut() {
    use crate::step_read::owned_insert::{CompoundDocument, DocumentOptions, Options};
    use num_traits::Zero;
    let source = crate::step_read::owned_insert::tests::graph(false, 4.);
    let plan = CoupledPlan::prepare(&source, CurveOptions::default(), None).unwrap();
    let top = CoupledTopology::prepare(&plan, 16_000_000, 64 * 1024 * 1024, None).unwrap();
    let text = top
        .write_mapped_candidate(
            [0.; 3],
            [1., 0., 0.],
            [0., 1., 0.],
            16_000_000,
            64 * 1024 * 1024,
            None,
        )
        .unwrap();
    let document = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
    let item = document.represented_members()[0].item_entity();
    let packet = document
        .member_packet(item, Options::default(), None)
        .unwrap();
    let p = packet.parsed();
    let original = p.support_bindings(Options::default(), None).unwrap();
    original.material_bounds(Options::default(), None).unwrap();
    let integer = |s: &str| BigInt::parse_bytes(s.trim_end_matches('.').as_bytes(), 10).unwrap();
    let weights = [1_i64, 11, 55, 165, 330, 462, 462, 330, 165, 55, 11, 1];
    let mut change = None;
    for t in p.tensors() {
        if !t.rows()[0].iter().zip(t.rows()[1]).all(|(a, b)| {
            (0..2).all(|axis| {
                p.point(*a).unwrap().coordinates()[axis].spelling()
                    == p.point(*b).unwrap().coordinates()[axis].spelling()
            })
        }) {
            continue;
        }
        let mut before = BigInt::zero();
        for (i, weight) in weights.into_iter().enumerate() {
            let a = p.point(t.rows()[0][i]).unwrap();
            let b = p.point(t.rows()[1][i]).unwrap();
            before += (integer(b.coordinates()[2].spelling())
                - integer(a.coordinates()[2].spelling()))
                * weight;
        }
        let point = p.point(t.rows()[1][5]).unwrap();
        let after = &before - integer(point.coordinates()[2].spelling()) * weights[5];
        // Exact midpoint Bernstein sums from actual encoded INTEGER leaves,
        // independent of interval-Jacobian implementation or source recipe.
        if before > BigInt::zero() && after < BigInt::zero() {
            change = Some((
                point.entity(),
                point.coordinates()[0].spelling(),
                point.coordinates()[1].spelling(),
            ));
            break;
        }
    }
    let (id, x, y) = change
        .expect("actual finite constant-radius strip with an independently demonstrated fold");
    let changed = change_record(&text, id, |_| {
        format!("#{id}=CARTESIAN_POINT('',({x},{y},0.));")
    });
    let d = CompoundDocument::parse(&changed, DocumentOptions::default(), None).unwrap();
    let packet = d.member_packet(item, Options::default(), None).unwrap();
    let b = packet
        .parsed()
        .support_bindings(Options::default(), None)
        .unwrap();
    assert!(b.material_bounds(Options::default(), None).is_err());
}
#[test]
fn real_mapped_pose_and_original_scale_are_not_discarded_or_guessed() {
    use crate::step_read::owned_insert::{CompoundDocument, DocumentOptions, Options};
    let source = crate::step_read::owned_insert::tests::graph(false, 4_f64.next_up());
    let plan = CoupledPlan::prepare(&source, CurveOptions::default(), None).unwrap();
    let top = CoupledTopology::prepare(&plan, 16_000_000, 64 * 1024 * 1024, None).unwrap();
    let origin = [40., -18., 7.];
    let u = [0.8, 0.6, 0.];
    let v = [-0.36, 0.48, 0.8];
    let text = top
        .write_mapped_candidate(origin, u, v, 16_000_000, 64 * 1024 * 1024, None)
        .unwrap();
    let document = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
    let member = document.represented_members()[0];
    let op = member.mapping_link().unwrap().operator_entity();
    let packet = document
        .member_packet(member.item_entity(), Options::default(), None)
        .unwrap();
    let mapping = packet
        .parsed()
        .physical_mapping(op, Options::default(), None)
        .unwrap();
    assert_eq!(mapping.handedness().unwrap(), 1);
    for (i, x) in origin.into_iter().enumerate() {
        let a = mapping.origin_mm()[i];
        assert!(a[0] <= x && a[1] >= x);
    }
    let bindings = packet
        .parsed()
        .support_bindings(Options::default(), None)
        .unwrap();
    let material = bindings.material_bounds(Options::default(), None).unwrap();
    let world = mapping.map_box_mm(material.bounds_mm()).unwrap();
    assert!(
        world
            .iter()
            .all(|a| a[0].is_finite() && a[1].is_finite() && a[0] < a[1])
    );
    let prefix = format!("#{op}=");
    let record = text.lines().find(|l| l.starts_with(&prefix)).unwrap();
    let attributes = record
        .split("CARTESIAN_TRANSFORMATION_OPERATOR(")
        .nth(1)
        .unwrap()
        .split(')')
        .next()
        .unwrap();
    let scale = attributes.rsplit(',').next().unwrap();
    for invalid in ["0.", "-1."] {
        let changed = change_record(&text, op, |s| {
            s.replacen(&format!(",{scale})"), &format!(",{invalid})"), 1)
        });
        let changed = CompoundDocument::parse(&changed, DocumentOptions::default(), None).unwrap();
        let p = changed
            .member_packet(member.item_entity(), Options::default(), None)
            .unwrap();
        assert!(
            p.parsed()
                .physical_mapping(op, Options::default(), None)
                .is_err()
        );
    }
    let missing = change_record(&text, op, |_| String::new());
    assert!(CompoundDocument::parse(&missing, DocumentOptions::default(), None).is_err());
    // A different positive scale is valid DIFFERENT geometry. The physical
    // reader must retain that fact; a future source matcher must not silently
    // restore the intended grid unit or authorize Native publication from it.
    let changed = change_record(&text, op, |s| s.replacen(&format!(",{scale})"), ",1.)", 1));
    let changed = CompoundDocument::parse(&changed, DocumentOptions::default(), None).unwrap();
    let p = changed
        .member_packet(member.item_entity(), Options::default(), None)
        .unwrap();
    let different = p
        .parsed()
        .physical_mapping(op, Options::default(), None)
        .unwrap();
    assert!(different.scale_bounds()[0] > plan.coefficient_unit_mm() * 1e20);
}
#[test]
fn real_coupled_append_budget_pose_and_late_cancel_leave_existing_members_unchanged() {
    use std::cell::Cell;
    let source = crate::step_read::owned_insert::tests::graph(false, 4.);
    let plan = CoupledPlan::prepare(&source, CurveOptions::default(), None).unwrap();
    let top = CoupledTopology::prepare(&plan, 16_000_000, 64 * 1024 * 1024, None).unwrap();
    let seed = || {
        let mut w = Writer::new();
        w.append_sphere(
            crate::sphere::Sphere::new(1.).unwrap(),
            [20., 0., 0.],
            [1., 0., 0.],
            [0., 1., 0.],
        )
        .unwrap();
        w
    };
    let before = seed().finish().unwrap();
    for (origin, u, v, work, bytes, error) in [
        (
            [f64::NAN, 0., 0.],
            [1., 0., 0.],
            [0., 1., 0.],
            16_000_000,
            64 * 1024 * 1024,
            Error::Invalid,
        ),
        (
            [0.; 3],
            [1., 0., 0.],
            [1., 0., 0.],
            16_000_000,
            64 * 1024 * 1024,
            Error::Unresolved,
        ),
        (
            [0.; 3],
            [1., 0., 0.],
            [0., 1., 0.],
            64_000_000,
            64 * 1024 * 1024,
            Error::Invalid,
        ),
        (
            [0.; 3],
            [1., 0., 0.],
            [0., 1., 0.],
            top.work_used() + 1,
            64 * 1024 * 1024,
            Error::Limit,
        ),
        (
            [0.; 3],
            [1., 0., 0.],
            [0., 1., 0.],
            16_000_000,
            1,
            Error::Limit,
        ),
    ] {
        let mut w = seed();
        assert_eq!(
            w.append_coupled_insert_candidate(&top, origin, u, v, work, bytes, None),
            Err(error)
        );
        assert_eq!(w.finish().unwrap(), before);
    }
    let polls = Cell::new(0usize);
    let count = || {
        polls.set(polls.get() + 1);
        false
    };
    let mut good = seed();
    good.append_coupled_insert_candidate(
        &top,
        [0.; 3],
        [1., 0., 0.],
        [0., 1., 0.],
        16_000_000,
        64 * 1024 * 1024,
        Some(&count),
    )
    .unwrap();
    let good = good.finish().unwrap();
    assert!(polls.get() > 1000);
    let stop = polls.get() - 1;
    polls.set(0);
    let cancel = || {
        polls.set(polls.get() + 1);
        polls.get() >= stop
    };
    let mut retry = seed();
    assert_eq!(
        retry.append_coupled_insert_candidate(
            &top,
            [0.; 3],
            [1., 0., 0.],
            [0., 1., 0.],
            16_000_000,
            64 * 1024 * 1024,
            Some(&cancel)
        ),
        Err(Error::Cancelled)
    );
    retry
        .append_coupled_insert_candidate(
            &top,
            [0.; 3],
            [1., 0., 0.],
            [0., 1., 0.],
            16_000_000,
            64 * 1024 * 1024,
            None,
        )
        .unwrap();
    assert_eq!(retry.finish().unwrap(), good);
}
#[test]
fn actual_mapped_numerator_candidates_keep_literal_material_inside_original_query_budget() {
    for (left, length) in [
        (false, 4.),
        (false, 4_f64.next_up()),
        (true, 4_f64.next_up()),
        (false, 4_f64.next_down()),
        (true, 8.),
    ] {
        let source = crate::step_read::owned_insert::tests::graph(left, length);
        let plan = CoupledPlan::prepare(&source, CurveOptions::default(), None).unwrap();
        let top = CoupledTopology::prepare(&plan, 16_000_000, 64 * 1024 * 1024, None).unwrap();
        let text = top
            .write_mapped_candidate(
                [0.; 3],
                [1., 0., 0.],
                [0., 1., 0.],
                16_000_000,
                64 * 1024 * 1024,
                None,
            )
            .unwrap();
        eprintln!("NUMERATOR_FILE {left}/{length:?}: {}bytes", text.len());
        let document = crate::step_read::owned_insert::CompoundDocument::parse(
            &text,
            Default::default(),
            None,
        )
        .unwrap();
        let member = document.represented_members()[0];
        let packet = document
            .member_packet(member.item_entity(), Default::default(), None)
            .unwrap();
        let parsed = packet.parsed();
        let mapping = parsed
            .physical_mapping(
                member.mapping_link().unwrap().operator_entity(),
                Default::default(),
                None,
            )
            .unwrap();
        let expected = plan.coefficient_unit_mm();
        let actual = mapping.scale_bounds();
        assert!(actual[0] <= expected && actual[1] >= expected);
        let bindings = parsed.support_bindings(Default::default(), None).unwrap();
        eprintln!(
            "NUMERATOR_BINDINGS {left}/{length:?}: {}work",
            bindings.work_used()
        );
        let material = bindings.material_bounds(Default::default(), None).unwrap();
        eprintln!(
            "NUMERATOR_MATERIAL {left}/{length:?}: {}work",
            material.work_used()
        );
        let scale = mapping.volume_scale_bounds().unwrap();
        let raw_width = (I::point(5e-7) / I::point(scale[1])).lo;
        let volume = material
            .occupied_volume(
                crate::step_read::owned_insert::FluxOptions {
                    maximum_width_mm3: raw_width,
                    ..Default::default()
                },
                None,
            )
            .unwrap();
        let [lo, hi] = volume.volume_mm3();
        let physical = I::new(lo, hi) * I::new(scale[0], scale[1]);
        eprintln!("NUMERATOR_PHYSICAL_VOLUME {left}/{length:?}: {physical:?}");
        assert!(physical.valid() && physical.lo > 0. && physical.hi < 1000.);
        assert!(physical.width() <= 1e-6);
    }
}
