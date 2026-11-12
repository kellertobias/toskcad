//! Unknown closed material is not silently interpreted as a known head recipe.
use super::*;

fn references(value: &Value<'_>, pending: &mut Vec<u32>) {
    match value {
        Value::Reference(id) => pending.push(*id),
        Value::List(values) => {
            for value in values {
                references(value, pending);
            }
        }
        _ => {}
    }
}
fn topology(name: &str) -> bool {
    matches!(
        name,
        "MANIFOLD_SOLID_BREP"
            | "CLOSED_SHELL"
            | "ADVANCED_FACE"
            | "FACE_OUTER_BOUND"
            | "FACE_BOUND"
            | "EDGE_LOOP"
            | "ORIENTED_EDGE"
            | "EDGE_CURVE"
            | "VERTEX_POINT"
    )
}
fn rewrite(text: &str, id: u32, record: String) -> String {
    let key = format!("#{id}=");
    let line = text.lines().find(|l| l.starts_with(&key)).unwrap();
    text.replacen(line, &format!("{key}{record}"), 1)
}

#[test]
fn opaque_closed_head_is_not_misidentified_from_its_thread_or_tail_count() {
    let thread = Thread::new(Spec {
        major_radius_mm: 1.5,
        pitch_mm: 0.5,
        length_mm: 0.5,
        depth_mm: 0.2,
        crest_width_mm: 0.0625,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts: 1,
        left_hand: false,
        root_radius_mm: None,
    })
    .unwrap();
    let screw = crate::profiled_screw::Screw::new(
        thread,
        crate::profiled_screw::Spec {
            head: crate::head_profile::Head::Countersunk {
                radius_mm: 3.,
                height_mm: 1.65,
                included_angle_deg: 90.,
                edge_chamfer_mm: 0.,
                neck_radius_mm: 0.1,
            },
            shank_length_mm: 1.,
            socket_across_flats_mm: 2.,
            socket_depth_mm: 0.6,
            minimum_wall_mm: 0.1,
        },
    )
    .unwrap();
    assert_eq!(screw.profile().segments().count(), 4);
    let text = crate::step::write_profiled_screw(&screw, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7)
        .unwrap();
    let data = parse_entities(&text, 1e-7).unwrap();
    let shell = Shell::new(&data, 1e-7, true);
    let solid = shell.represented_solids(true).unwrap()[0];
    let faces = shell.member_faces(solid).unwrap();
    let (surface, forward, bounds) = face(&shell, faces[6]).unwrap();
    assert!(forward && bounds.len() == 2);
    // Filling the socket leaves the same real meridian and thread, with a
    // complete circular top cap. The resulting opaque head is outside the
    // supported driven-screw family, though its material boundary is closed.
    let text = rewrite(
        &text,
        faces[6],
        format!("ADVANCED_FACE('',(#{}),#{surface},.T.);", bounds[0]),
    );
    let closed = reference(&shell.args(solid, "MANIFOLD_SOLID_BREP").unwrap()[1]).unwrap();
    let kept = faces[..faces.len() - 7]
        .iter()
        .map(|id| format!("#{id}"))
        .collect::<Vec<_>>()
        .join(",");
    let text = rewrite(&text, closed, format!("CLOSED_SHELL('',({kept}));"));
    let data = parse_entities(&text, 1e-7).unwrap();
    let shell = Shell::new(&data, 1e-7, true);
    let mut pending = vec![solid];
    let mut live = BTreeSet::new();
    while let Some(id) = pending.pop() {
        if !live.insert(id) {
            continue;
        }
        for entity in shell.components(id).unwrap() {
            for value in &entity.args {
                references(value, &mut pending);
            }
        }
    }
    let dead = data
        .iter()
        .filter(|(id, entities)| !live.contains(id) && entities.iter().any(|e| topology(e.name)))
        .map(|(id, _)| format!("#{id}="))
        .collect::<Vec<_>>();
    let text = text
        .lines()
        .filter(|line| !dead.iter().any(|prefix| line.starts_with(prefix)))
        .map(|line| format!("{line}\n"))
        .collect::<String>();
    compound::read_members_by(&text, 1e-7, true, true, |shell, solid, _| {
        let faces = shell.member_faces(solid)?;
        helical::profiled::closed_incidence(shell, &faces)?;
        require(faces.len() == 11)?;
        require(read_member(shell, solid, &faces)?.is_none())
    })
    .unwrap();
    assert!(crate::step_read::read_import_solid(&text, 1e-7).is_err());
}
