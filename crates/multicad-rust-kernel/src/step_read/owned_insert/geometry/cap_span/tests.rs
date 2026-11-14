use super::*;
use num_bigint::BigInt;
fn cache(p: &Parsed<'_>) -> Vec<Option<[I; 3]>> {
    p.points()
        .iter()
        .map(|p| {
            (p.coordinates().len() == 3).then(|| {
                std::array::from_fn(|i| {
                    let [a, b] = p.coordinates()[i].bounds();
                    I::new(a, b)
                })
            })
        })
        .collect()
}
fn decimal(n: BigInt) -> String {
    let digits = (n * BigInt::from(5).pow(100)).to_string();
    format!(
        "{}.{}",
        &digits[..digits.len() - 100],
        &digits[digits.len() - 100..]
    )
}
#[test]
fn immutable_parsed_near_parallel_span_requires_true_exact_angular_fallback() {
    let text = crate::step_read::owned_insert::tests::file(false);
    let original = parse(&text, Options::default(), None).unwrap();
    let binding = original.support_bindings(Options::default(), None).unwrap();
    let co = &original.faces()[binding.cap_face_indices()[0]].loops()[0].coedges()[0];
    let id = original.edge(co.edge).unwrap().curve_entity();
    let curve = original.curve(id).unwrap();
    let mut changed = text.clone();
    for (i, &point) in curve.pole_entities()[..12].iter().enumerate() {
        let x = BigInt::from(8 + i) << 97_usize;
        let y = &x + BigInt::from(i);
        let z = original.point(point).unwrap().coordinates()[2].spelling();
        changed = crate::step_read::owned_insert::tests::record_change(&changed, point, |_| {
            format!(
                "#{point}=CARTESIAN_POINT('',({},{},{z}));",
                decimal(x),
                decimal(y)
            )
        });
    }
    // This is the actual numeric-stage contract only: the modified surrounding
    // supports/wires have NOT earned material, incidence or decoded mass.
    let actual = parse(&changed, Options::default(), None).unwrap();
    let cache = cache(&actual);
    let mut w = Work::new(16_000_000, None).unwrap();
    let finite = ParsedCapSpan::load(&actual, id, 0, &cache, &mut w).unwrap();
    assert!(finite.points()[0][0].lo <= 1. && finite.points()[0][0].hi >= 1.);
    let leaf = actual
        .point(actual.curve(id).unwrap().pole_entities()[1])
        .unwrap();
    assert_eq!(leaf.coordinates()[0].value(), leaf.coordinates()[1].value());
    assert_ne!(
        leaf.coordinates()[0].spelling(),
        leaf.coordinates()[1].spelling()
    );
    let mut memory = Memory::new(0, 64 * 1024 * 1024).unwrap();
    let mut stack = Vec::new();
    let lower = finite
        .angular_lower(1., &mut stack, &mut memory, &mut w)
        .unwrap();
    assert!(lower > 0. && lower <= 11. * 2_f64.powi(-100));
    // Interrupt real finite-curve queries throughout the numeric fallback,
    // including original-leaf lookup. Each request owns fresh scratch/work.
    let calls = std::cell::Cell::new(0_usize);
    let count = || {
        calls.set(calls.get() + 1);
        false
    };
    let mut counted = Work::new(16_000_000, Some(&count)).unwrap();
    let mut memory_counted = Memory::new(0, 64 * 1024 * 1024).unwrap();
    stack.clear();
    assert_eq!(
        finite
            .angular_lower(1., &mut stack, &mut memory_counted, &mut counted)
            .unwrap(),
        lower
    );
    let total = calls.get();
    for at in (0..128)
        .map(|i| total * i / 128)
        .chain(0..total.min(128))
        .chain(total.saturating_sub(128)..total)
    {
        let n = std::cell::Cell::new(0);
        let stop = || {
            let old = n.get();
            n.set(old + 1);
            old >= at
        };
        let mut limited = Work::new(16_000_000, Some(&stop)).unwrap();
        let mut memory = Memory::new(0, 64 * 1024 * 1024).unwrap();
        stack.clear();
        assert!(matches!(
            finite.angular_lower(1., &mut stack, &mut memory, &mut limited),
            Err(Error::Cancelled)
        ));
    }
    stack.clear();

    assert!(matches!(
        finite.angular_lower(-1., &mut stack, &mut memory, &mut w),
        Err(Error::Unsupported)
    ));
    assert!(matches!(
        finite.angular_lower(0., &mut stack, &mut memory, &mut w),
        Err(Error::Invalid)
    ));
    assert!(matches!(
        ParsedCapSpan::load(&actual, u32::MAX, 0, &cache, &mut w),
        Err(Error::Invalid)
    ));
    assert!(matches!(
        ParsedCapSpan::load(&actual, id, usize::MAX, &cache, &mut w),
        Err(Error::Unsupported)
    ));
    let line = actual.curves().iter().find(|c| c.degree() == 1).unwrap();
    assert!(matches!(
        ParsedCapSpan::load(&actual, line.entity(), 0, &cache, &mut w),
        Err(Error::Unsupported)
    ));
}

#[test]
fn parsed_span_queries_preserve_finite_domain_and_caller_resource_scopes() {
    use std::cell::Cell;
    let text = crate::step_read::owned_insert::tests::file(false);
    let p = parse(&text, Options::default(), None).unwrap();
    let id = p
        .curves()
        .iter()
        .find(|c| {
            c.degree() == 11 && p.point(c.pole_entities()[0]).unwrap().coordinates().len() == 3
        })
        .unwrap()
        .entity();
    let cache = cache(&p);
    let mut w = Work::new(16_000_000, None).unwrap();
    let finite = ParsedCapSpan::load(&p, id, 0, &cache, &mut w).unwrap();
    let empty = vec![None; p.points().len()];
    assert!(matches!(
        ParsedCapSpan::load(&p, id, 0, &empty, &mut w),
        Err(Error::Invalid)
    ));
    let mut short = Work::new(1, None).unwrap();
    assert!(matches!(
        ParsedCapSpan::load(&p, id, 0, &cache, &mut short),
        Err(Error::Limit)
    ));
    let stop = || true;
    let mut cancelled = Work::new(16_000_000, Some(&stop)).unwrap();
    assert!(matches!(
        ParsedCapSpan::load(&p, id, 0, &cache, &mut cancelled),
        Err(Error::Cancelled)
    ));
    let mut initial_work = Work::new(16_000_000, None).unwrap();
    let mut initial_memory = Memory::new(0, 64 * 1024 * 1024).unwrap();
    let mut initial_stack = Vec::new();
    let sign = if finite
        .angular_lower(
            1.,
            &mut initial_stack,
            &mut initial_memory,
            &mut initial_work,
        )
        .is_ok()
    {
        1.
    } else {
        -1.
    };
    let calls = Cell::new(0_usize);
    let count = || {
        calls.set(calls.get() + 1);
        false
    };
    let mut complete = Work::new(16_000_000, Some(&count)).unwrap();
    let mut memory = Memory::new(0, 64 * 1024 * 1024).unwrap();
    let mut stack = Vec::new();
    let baseline = finite
        .angular_lower(sign, &mut stack, &mut memory, &mut complete)
        .unwrap();
    for at in (0..128)
        .map(|i| calls.get() * i / 128)
        .chain([calls.get() - 1])
    {
        let n = Cell::new(0);
        let stop = || {
            let old = n.get();
            n.set(old + 1);
            old >= at
        };
        let mut w = Work::new(16_000_000, Some(&stop)).unwrap();
        let mut memory = Memory::new(0, 64 * 1024 * 1024).unwrap();
        stack.clear();
        assert!(matches!(
            finite.angular_lower(sign, &mut stack, &mut memory, &mut w),
            Err(Error::Cancelled)
        ));
    }
    stack.clear();
    let mut enough = Work::new(complete.used, None).unwrap();
    let mut memory = Memory::new(0, 64 * 1024 * 1024).unwrap();
    assert_eq!(
        finite
            .angular_lower(sign, &mut stack, &mut memory, &mut enough)
            .unwrap(),
        baseline
    );
}
