//! Private schema contracts over actual reparsed physical record bytes.
//! Neutral parsing may reject these before the full-file decoder reaches here.
use super::*;

#[test]
fn actual_graph_record_helpers_reject_wrong_shapes_and_preserve_valid_retry() {
    let actual = super::super::tests::candidate(false);
    let original = crate::step_read::part21::parse(&actual).unwrap();
    for (name, replacements) in [
        (
            "FACE_OUTER_BOUND",
            vec!["FACE_OUTER_BOUND('',#1,.F.)", "FACE_OUTER_BOUND()"],
        ),
        ("EDGE_LOOP", vec!["EDGE_LOOP()", "EDGE_LOOP('',())"]),
        ("ORIENTED_EDGE", vec!["ORIENTED_EDGE('')"]),
        ("VERTEX_POINT", vec!["VERTEX_POINT()"]),
        ("PLANE", vec!["PLANE()"]),
        ("AXIS2_PLACEMENT_3D", vec!["AXIS2_PLACEMENT_3D('')"]),
        ("EDGE_CURVE", vec!["EDGE_CURVE()"]),
    ] {
        let owner = if name == "AXIS2_PLACEMENT_3D" {
            let parts = original
                .values()
                .find(|parts| parts.iter().any(|e| e.name == "PLANE"))
                .unwrap();
            reference(&parts[0].args[1]).unwrap()
        } else if name == "EDGE_LOOP" || name == "ORIENTED_EDGE" {
            let bound = original
                .values()
                .find(|parts| parts[0].name == "FACE_OUTER_BOUND")
                .unwrap();
            let loop_id = reference(&bound[0].args[1]).unwrap();
            if name == "EDGE_LOOP" {
                loop_id
            } else {
                reference(&list(&original[&loop_id][0].args[1]).unwrap()[0]).unwrap()
            }
        } else {
            *original
                .iter()
                .find(|(_, parts)| parts.iter().any(|e| e.name == name))
                .unwrap()
                .0
        };
        let components = &original[&owner];
        let component = components.iter().find(|e| e.name == name).unwrap();
        let source = component.source.unwrap();
        let start = source.as_ptr() as usize - actual.as_ptr() as usize;
        for replacement in replacements {
            let mut changed = actual.clone();
            changed.replace_range(start..start + source.len(), replacement);
            let records = crate::step_read::part21::parse(&changed).unwrap();
            let mut budget =
                Budget::new(crate::step_read::npt_stock::Options::default(), None).unwrap();
            let result = match name {
                "VERTEX_POINT" => vertex(&records, owner, true, &mut budget).map(|_| ()),
                "PLANE" => plane(&records, owner, &mut budget).map(|_| ()),
                "AXIS2_PLACEMENT_3D" => {
                    let plane_id = original
                        .iter()
                        .find(|(_, parts)| {
                            parts.iter().any(|e| {
                                e.name == "PLANE" && reference(&e.args[1]).ok() == Some(owner)
                            })
                        })
                        .unwrap()
                        .0;
                    plane(&records, *plane_id, &mut budget).map(|_| ())
                }
                "EDGE_CURVE" => edge(&records, owner, &mut budget).map(|_| ()),
                _ => {
                    let bound = if name == "FACE_OUTER_BOUND" {
                        owner
                    } else {
                        let loop_id = if name == "EDGE_LOOP" {
                            owner
                        } else {
                            *original
                                .iter()
                                .find(|(_, parts)| {
                                    parts.iter().any(|e| {
                                        e.name == "EDGE_LOOP"
                                            && list(&e.args[1])
                                                .unwrap()
                                                .iter()
                                                .any(|v| reference(v).ok() == Some(owner))
                                    })
                                })
                                .unwrap()
                                .0
                        };
                        *original
                            .iter()
                            .find(|(_, parts)| {
                                parts.iter().any(|e| {
                                    e.name == "FACE_OUTER_BOUND"
                                        && reference(&e.args[1]).ok() == Some(loop_id)
                                })
                            })
                            .unwrap()
                            .0
                    };
                    wire(&records, bound, true, &mut budget).map(|_| ())
                }
            };
            assert!(
                matches!(result, Err(Error::Invalid | Error::Unsupported)),
                "{name}: {replacement}"
            );
        }
    }
    let mut budget = Budget::new(crate::step_read::npt_stock::Options::default(), None).unwrap();
    let mut exact_vertices = 0;
    let mut inexact_vertices = 0;
    for (&id, parts) in &original {
        match parts[0].name {
            "FACE_OUTER_BOUND" => {
                wire(&original, id, true, &mut budget).unwrap();
            }
            "VERTEX_POINT" => {
                vertex(&original, id, false, &mut budget).unwrap();
                match vertex(&original, id, true, &mut budget) {
                    Ok(_) => exact_vertices += 1,
                    Err(Error::Unsupported) => inexact_vertices += 1,
                    Err(e) => panic!("unexpected actual vertex error: {e:?}"),
                }
            }
            "PLANE" => {
                plane(&original, id, &mut budget).unwrap();
            }
            "EDGE_CURVE" => {
                edge(&original, id, &mut budget).unwrap();
            }
            _ => {}
        }
    }
    assert!(exact_vertices > 0);
    assert!(inexact_vertices > 0);
}

#[test]
fn actual_graph_consumption_reserves_its_tree_and_preserves_cancel_retry() {
    let actual = super::super::tests::candidate(false);
    let records = crate::step_read::part21::parse(&actual).unwrap();
    let tree = (1 + records.len() / 5) * (20 * size_of::<usize>() + 11 * size_of::<u32>());
    let options = crate::step_read::npt_stock::Options {
        maximum_bytes: crate::step_read::npt_stock::SCRATCH + tree - 1,
        ..Default::default()
    };
    let mut bytes = Budget::new(options, None).unwrap();
    assert_eq!(consume(&records, &mut bytes), Err(Error::Limit));
    let stop = || true;
    let mut cancelled = Budget::new(Default::default(), Some(&stop)).unwrap();
    assert_eq!(consume(&records, &mut cancelled), Err(Error::Cancelled));
    // The single-member parser budget is smaller than this full-head graph.
    let mut limited = Budget::new(Default::default(), None).unwrap();
    assert_eq!(consume(&records, &mut limited), Err(Error::Limit));
    let retry = DecodedFile::parse(&actual, Default::default(), None).unwrap();
    assert_eq!(retry.actual_file(), actual);
}
