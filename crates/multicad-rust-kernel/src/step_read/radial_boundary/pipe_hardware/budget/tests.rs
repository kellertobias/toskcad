use super::*;
use crate::{Rational, step::Writer, step_read::part21};
fn actual_box() -> String {
    let mut writer = Writer::new();
    writer
        .append_box(
            [2, 3, 4].map(|n| Rational::new(n, 1).unwrap()),
            [0.; 3],
            [1., 0., 0.],
            [0., 1., 0.],
        )
        .unwrap();
    writer.finish().unwrap()
}
fn vectors(v: &Vec<Value<'_>>) -> (usize, usize) {
    let mut total = v.capacity() * size_of::<Value<'_>>();
    let mut largest = total;
    for value in v {
        let part = match value {
            Value::List(v) => Some(vectors(v)),
            Value::Typed(e) => Some(vectors(&e.args)),
            _ => None,
        };
        if let Some((n, m)) = part {
            total += n;
            largest = largest.max(m);
        }
    }
    (total, largest)
}
#[test]
fn actual_cold_file_accounts_vectors_and_literal_comments_before_parser_allocation() {
    let ordinary = actual_box();
    let basic = reserve(&ordinary, 64_000_000, 64 * 1024 * 1024, None).unwrap();
    let comment = format!("/* {} #999=FAKE((9,8)) */", "λ".repeat(5000));
    let actual = ordinary
        .replacen(
            "FILE_NAME('MultiCAD'",
            "FILE_NAME('μ ''#999=FAKE((1))'' /*literal*/'",
            1,
        )
        .replacen("DATA;", &format!("{comment}DATA;"), 1);
    let plan = reserve(&actual, 64_000_000, 64 * 1024 * 1024, None).unwrap();
    assert_eq!(
        plan.represented_record_count,
        basic.represented_record_count
    );
    assert_eq!(plan.references, basic.references);
    let arena = part21::parse(&actual).unwrap();
    let mut bytes = actual.len();
    let mut largest = 0;
    for records in arena.values() {
        let n = records.capacity() * size_of::<Entity<'_>>();
        bytes += n;
        largest = largest.max(n);
        for e in records {
            let (n, m) = vectors(&e.args);
            bytes += n;
            largest = largest.max(m);
        }
    }
    let node = 20 * size_of::<usize>() + 11 * (size_of::<u32>() + size_of::<Vec<Entity<'_>>>());
    bytes += (1 + arena.len() / 5) * node
        + size_of::<std::collections::BTreeMap<u32, Vec<Entity<'_>>>>();
    assert!(plan.retained_arena_bytes_reserved >= bytes);
    assert!(plan.parser_peak_bytes_reserved >= bytes + largest + 64 * node + 64 * 1024);
    assert_eq!(
        plan.aggregate_work_reserved,
        plan.scan_work_used + plan.delegated_parser_work_reserved
    );
    let exact = reserve(
        &actual,
        plan.aggregate_work_reserved,
        plan.parser_peak_bytes_reserved,
        None,
    )
    .unwrap();
    assert_eq!(exact.aggregate_work_reserved, plan.aggregate_work_reserved);
    assert!(matches!(
        reserve(
            &actual,
            plan.aggregate_work_reserved - 1,
            plan.parser_peak_bytes_reserved,
            None
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        reserve(
            &actual,
            plan.aggregate_work_reserved,
            plan.parser_peak_bytes_reserved - 1,
            None
        ),
        Err(Error::Limit)
    ));
}
#[test]
fn actual_cold_file_rejects_malformed_depth_tokens_and_cancellation_without_parser_side_effects() {
    let healthy = actual_box();
    let insert = |record: &str| {
        let at = healthy.rfind("ENDSEC;").unwrap();
        format!("{}{record}\n{}", &healthy[..at], &healthy[at..])
    };
    for (text, expected) in [
        (insert("#999999=(A()A());"), Error::Invalid),
        (
            insert(&format!("#999999=X('{}');", "x".repeat(65_537))),
            Error::Limit,
        ),
        (
            insert(&format!("#999999={}();", "X".repeat(65_537))),
            Error::Limit,
        ),
        (
            insert(&format!(
                "#999999=X({}0{});",
                "(".repeat(25),
                ")".repeat(25)
            )),
            Error::Limit,
        ),
        (
            healthy.replacen("DATA;", "DATA('unsupported');", 1),
            Error::Unsupported,
        ),
        (format!("{healthy} trailing"), Error::Invalid),
        (
            healthy.replacen("HEADER;", "HEADER;/* unclosed", 1),
            Error::Invalid,
        ),
    ] {
        assert_eq!(
            reserve(&text, 64_000_000, 64 * 1024 * 1024, None).unwrap_err(),
            expected
        );
    }
    let calls = std::cell::Cell::new(0usize);
    let count = || {
        calls.set(calls.get() + 1);
        false
    };
    let plan = reserve(&healthy, 64_000_000, 64 * 1024 * 1024, Some(&count)).unwrap();
    for stop in [1, calls.get() / 2, calls.get() - 1] {
        let visited = std::cell::Cell::new(0usize);
        let cancel = || {
            let n = visited.get() + 1;
            visited.set(n);
            n >= stop
        };
        assert!(matches!(
            reserve(&healthy, 64_000_000, 64 * 1024 * 1024, Some(&cancel)),
            Err(Error::Cancelled)
        ));
        assert_eq!(visited.get(), stop);
    }
    assert_eq!(
        reserve(&healthy, 64_000_000, 64 * 1024 * 1024, None)
            .unwrap()
            .aggregate_work_reserved,
        plan.aggregate_work_reserved
    );
}

#[test]
fn long_whitespace_is_charged_in_bounded_cancellable_chunks() {
    let ordinary = actual_box();
    let padded = format!("{}{}", " ".repeat(20_000), ordinary);
    let calls = std::cell::Cell::new(0usize);
    let cancel = || {
        calls.set(calls.get() + 1);
        calls.get() == 3
    };
    assert!(matches!(
        reserve(&padded, 64_000_000, 64 * 1024 * 1024, Some(&cancel)),
        Err(Error::Cancelled)
    ));
    assert_eq!(calls.get(), 3);
    let plan = reserve(&padded, 64_000_000, 64 * 1024 * 1024, None).unwrap();
    assert_eq!(
        plan.represented_record_count,
        part21::parse(&ordinary).unwrap().len()
    );
    assert!(plan.scan_work_used >= padded.len());
}
