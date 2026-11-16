use super::*;
#[test]
fn cached_constants_preserve_exact_original_mm_and_account_retention() {
    let mut m = Memory::new(19, 64 * 1024 * 1024).unwrap();
    let mut w = Work::new(16_000_000, None).unwrap();
    let mut cache = DecimalCache::new(&mut m, &mut w).unwrap();
    for twos in [0, 1, 7, 128, 1074] {
        for fives in [0, 1, 9, 64] {
            for numerator in [-37, 0, 1, 1234567] {
                let value = N(BigRational::new(
                    numerator.into(),
                    (BigInt::one() << twos) * BigInt::from(5_u32).pow(fives),
                ));
                let before = m.used;
                let old = decimal(&value, &mut m, &mut w).unwrap();
                let expected = old.as_str().to_owned();
                drop(old);
                m.used = before;
                for _ in 0..2 {
                    let before = m.used;
                    let text = cache.format(&value, &mut m, &mut w).unwrap();
                    assert_eq!(text.as_str(), expected);
                    let output = text.retained_bytes();
                    drop(text);
                    m.used -= output;
                    assert!(m.used >= before);
                }
            }
        }
    }
    cache.release(&mut m);
    assert_eq!(m.used, 19);
}
#[test]
fn refusal_restores_only_memory_and_keeps_cache_atomic() {
    for budget in [1, 10, 100, 1000] {
        let mut m = Memory::new(0, 1_000_000).unwrap();
        let mut w = Work::new(1_000_000, None).unwrap();
        let mut cache = DecimalCache::new(&mut m, &mut w).unwrap();
        let before = m.used;
        let mut limited = Work::new(budget, None).unwrap();
        let value = N(BigRational::new(7.into(), BigInt::one() << 1074));
        assert!(matches!(
            cache.format(&value, &mut m, &mut limited),
            Err(Error::Limit)
        ));
        assert_eq!(m.used, before);
        assert!(cache.powers.is_empty());
        cache.release(&mut m);
        assert_eq!(m.used, 0);
    }
    let mut m = Memory::new(0, 1_000_000).unwrap();
    let mut w = Work::new(1_000_000, None).unwrap();
    let mut cache = DecimalCache::new(&mut m, &mut w).unwrap();
    let before = m.used;
    assert!(matches!(
        cache.format(&N(BigRational::new(1.into(), 3.into())), &mut m, &mut w),
        Err(Error::Unresolved)
    ));
    assert_eq!(m.used, before);
    cache.release(&mut m);
}
#[test]
fn real_allocation_bit_and_cancellation_frontiers_are_atomic() {
    use std::cell::Cell;
    for bytes in [0, 1, size_of::<DecimalCache>()] {
        let mut m = Memory::new(0, bytes).unwrap();
        let mut w = Work::new(1_000_000, None).unwrap();
        assert!(matches!(
            DecimalCache::new(&mut m, &mut w),
            Err(Error::Limit)
        ));
        assert_eq!(m.used, 0);
    }
    let q = N(BigRational::new(7.into(), BigInt::one() << 128));
    let execute = |stop: usize| {
        let calls = Cell::new(0);
        let cancel = || {
            calls.set(calls.get() + 1);
            calls.get() == stop
        };
        let mut m = Memory::new(0, 1_000_000).unwrap();
        let mut w = Work::new(1_000_000, Some(&cancel)).unwrap();
        let result = (|| {
            let mut c = DecimalCache::new(&mut m, &mut w)?;
            let a = c.format(&q, &mut m, &mut w)?;
            let b = c.format(&q, &mut m, &mut w)?;
            drop((a, b));
            c.release(&mut m);
            Ok::<_, Error>(())
        })();
        (result, calls.get())
    };
    let (result, calls) = execute(usize::MAX);
    assert_eq!(result, Ok(()));
    for stop in 1..=calls {
        assert_eq!(execute(stop).0, Err(Error::Cancelled));
    }
    let mut m = Memory::new(0, 1_000_000).unwrap();
    let mut w = Work::new(1_000_000, None).unwrap();
    let mut c = DecimalCache::new(&mut m, &mut w).unwrap();
    let before = m.used;
    let huge = N(BigRational::from_integer(BigInt::one() << 16385));
    assert!(matches!(c.format(&huge, &mut m, &mut w), Err(Error::Limit)));
    assert_eq!(m.used, before);
    c.release(&mut m);
}
#[test]
fn successful_number_does_not_require_optional_cache_retention() {
    let q = N(BigRational::new(7.into(), BigInt::one() << 128));
    let run = |budget| {
        let mut m = Memory::new(0, 1_000_000).unwrap();
        let mut initial = Work::new(1_000_000, None).unwrap();
        let mut c = DecimalCache::new(&mut m, &mut initial).unwrap();
        let mut w = Work::new(budget, None).unwrap();
        let result = c.format(&q, &mut m, &mut w).map(|s| s.as_str().to_owned());
        (result, c.powers.len())
    };
    let (mut lo, mut hi) = (1_usize, 100_000_usize);
    while lo + 1 < hi {
        let mid = lo.midpoint(hi);
        if run(mid).0.is_ok() {
            hi = mid
        } else {
            lo = mid
        }
    }
    let (result, retained) = run(hi);
    assert!(result.is_ok());
    assert_eq!(retained, 0);
    assert!(matches!(run(lo).0, Err(Error::Limit)));
}
