use super::*;
use crate::{Rational, step::Writer};

#[test]
fn genuine_record_arena_rejects_truncated_byte_owner_then_retains_full_original() {
    let mut writer = Writer::new();
    writer
        .append_box(
            [2, 3, 4].map(|v| Rational::new(v, 1).unwrap()),
            [0.; 3],
            [1., 0., 0.],
            [0., 1., 0.],
        )
        .unwrap();
    let file = writer.finish().unwrap();
    let path = std::env::temp_dir().join(format!("record-owner-box-{}.step", std::process::id()));
    std::fs::write(&path, &file).unwrap();
    drop(file);
    let actual = std::fs::read_to_string(&path).unwrap();
    std::fs::remove_file(&path).unwrap();
    let arena = part21::parse(&actual).unwrap();
    let member = *arena
        .iter()
        .find(|(_, entities)| {
            entities
                .iter()
                .any(|e| matches!(e.name, "FACETED_BREP" | "MANIFOLD_SOLID_BREP"))
        })
        .unwrap()
        .0;
    let first_source = arena.values().next().unwrap()[0].source.unwrap();
    let end = first_source.as_ptr() as usize - actual.as_ptr() as usize + first_source.len();
    assert!(end < actual.len());
    let truncated = &actual[..end - 1];
    // Independent helper provenance contract only: all Entity.source leaves
    // are genuine parser output. Public read always supplies the full owner;
    // this check must not be reported as a reachable public decoder failure.
    assert!(matches!(
        Records::from_actual(truncated, &arena, member, None),
        Err(GeometryError::Invalid)
    ));
    let expected = arena
        .iter()
        .flat_map(|(&id, entities)| {
            entities
                .iter()
                .enumerate()
                .map(move |(index, e)| (id, index, e.source.unwrap().to_owned()))
        })
        .collect::<Vec<_>>();
    let retained = Records::from_actual(&actual, &arena, member, None).unwrap();
    assert_eq!(retained.record_count(), arena.len());
    assert_eq!(retained.member, member);
    assert_ne!(retained.input().as_ptr(), actual.as_ptr());
    assert_eq!(retained.input(), actual);
    drop(arena);
    drop(actual);
    for (id, index, original) in expected {
        assert_eq!(
            retained.component(retained.record_index(id).unwrap(), index),
            Some(original.as_str())
        );
    }
    assert!(retained.input().starts_with("ISO-10303-21;"));
    assert!(retained.input().ends_with("END-ISO-10303-21;\n"));
    assert!(retained.bytes().unwrap() > retained.input().len());
}
