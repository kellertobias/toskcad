use super::*;
use crate::{
    head_profile::{Head, Profile},
    helical_thread::{self, Thread},
    phillips_drive::{Cavity, Cutter, PhillipsScrew},
    profiled_screw::SourceSpec,
};
pub(super) fn fixture(index: usize, left: bool, nominal: bool, starts: u32) -> PhillipsScrew {
    let (s, d) = crate::phillips_drive::test_fixture(index);
    let cavity = Cavity::new(Cutter::new(s, d).unwrap(), 1e-5, 100_000).unwrap();
    let (rs, pitch, r, h, crown, land, neck) = [
        (1., 0.4, 2., 1.6, 3.2, 1.25, 0.1),
        (1.5, 0.5, 2.8, 2.4, 3., 2.1, 0.1),
        (1.5, 0.5, 2.75, 1.65, 0., 0., 0.1),
        (2., 0.7, 4., 3.1, 6.5, 3., 0.2),
        (3., 1., 6., 4.6, 10., 4.2, 0.25),
        (4., 1.25, 8., 6., 13., 5.4, 0.4),
    ][index];
    let head = if index == 2 {
        Head::Countersunk {
            radius_mm: r,
            height_mm: h,
            included_angle_deg: 90.,
            edge_chamfer_mm: 0.05,
            neck_radius_mm: neck,
        }
    } else {
        Head::Crowned {
            radius_mm: r,
            height_mm: h,
            crown_radius_mm: crown,
            top_land_radius_mm: land,
            neck_radius_mm: neck,
        }
    };
    let thread = Thread::new(helical_thread::Spec {
        major_radius_mm: if nominal { rs } else { rs - 0.03 },
        pitch_mm: pitch,
        length_mm: 2.,
        depth_mm: 0.6 * pitch,
        crest_width_mm: 0.125 * pitch,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts,
        left_hand: left,
        root_radius_mm: None,
    })
    .unwrap();
    PhillipsScrew::new(
        thread,
        Profile::new_for_thread_union(head, rs, 0.).unwrap(),
        SourceSpec {
            head: Some(head),
            shank_length_mm: 0.,
            minimum_wall_mm: 1e-5,
        },
        cavity,
        100_000,
        &|| false,
    )
    .unwrap()
}

#[test]
fn actual_metric_material_boundaries_parse_without_any_authored_source_argument() {
    for index in 0..6 {
        for left in [false, true] {
            let source = fixture(index, left, false, 1);
            let (origin, u, v) = if left {
                ([105., -20., 60.], [0.8, 0.6, 0.], [0.36, -0.48, -0.8])
            } else {
                ([0.; 3], [1., 0., 0.], [0., 1., 0.])
            };
            let text = crate::step::write_phillips_screw(&source, origin, u, v, 1e-7).unwrap();
            let b = read_phillips_screw_boundary(&text, 1e-7, 1_000_000).unwrap();
            assert_eq!(
                (b.vertices().len(), b.edges().len(), b.faces().len()),
                (49, 76, 29)
            );
            assert!(b.faces()[0].material_entity.is_none());
            assert_eq!(b.faces()[0].surface_entity, b.mouth().material_surface);
            assert!(
                b.mouth().material_face > 0
                    && b.mouth().material_bound > 0
                    && b.mouth().material_forward
            );
            assert_eq!(
                b.edges()
                    .iter()
                    .filter(|e| matches!(e.curve, Curve::CubicIntersection { .. }))
                    .count(),
                32
            );
            assert_eq!(
                b.edges()
                    .iter()
                    .any(|e| matches!(e.curve, Curve::Ellipse { .. })),
                index == 0
            );
            assert!(
                b.faces()
                    .iter()
                    .all(|f| f.material_entity.is_none() || f.material_forward != f.cavity_forward)
            );
        }
    }
}
fn text() -> String {
    let source = fixture(0, false, false, 1);
    crate::step::write_phillips_screw(&source, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7).unwrap()
}
#[test]
fn literal_nonincident_data_stays_literal_and_is_not_mislabeled_a_geometry_certificate() {
    let text = text();
    let boundary = read_phillips_screw_boundary(&text, 1e-7, 1_000_000).unwrap();
    let (pole, old) = boundary
        .edges()
        .iter()
        .find_map(|e| match &e.curve {
            Curve::CubicIntersection {
                pole_entities,
                spans,
                ..
            } => Some((pole_entities[1], spans[0].poles_mm[1])),
            _ => None,
        })
        .unwrap();
    let line = text
        .lines()
        .find(|l| l.starts_with(&format!("#{pole}=")))
        .unwrap();
    let replacement = format!(
        "#{pole}=CARTESIAN_POINT('',({:.16E},{:.16E},{:.16E}));",
        old[0] + 0.001,
        old[1],
        old[2]
    );
    let changed = text.replacen(line, &replacement, 1);
    let actual = read_phillips_screw_boundary(&changed, 1e-7, 1_000_000).unwrap();
    let point = actual
        .edges()
        .iter()
        .find_map(|e| match &e.curve {
            Curve::CubicIntersection {
                pole_entities,
                spans,
                ..
            } if pole_entities[1] == pole => Some(spans[0].poles_mm[1]),
            _ => None,
        })
        .unwrap();
    assert_eq!(point, [old[0] + 0.001, old[1], old[2]]);
    // This is intentionally only a parsed graph. The changed literal cubic is
    // retained; no source-known geometric fit or nominal replacement is invoked.
    assert_ne!(point, old);
    let line_edge = boundary
        .edges()
        .iter()
        .find(|e| matches!(e.curve, Curve::Line { .. }))
        .unwrap();
    let Curve::Line { vector_entity, .. } = line_edge.curve else {
        panic!("line selected")
    };
    let vector = text
        .lines()
        .find(|l| l.starts_with(&format!("#{vector_entity}=")))
        .unwrap();
    let prefix = vector.rsplit_once(',').unwrap().0;
    let changed = text.replacen(vector, &format!("{prefix},3.);"), 1);
    let actual = read_phillips_screw_boundary(&changed, 1e-7, 1_000_000).unwrap();
    assert!(
        actual
            .edges()
            .iter()
            .any(|e| matches!(e.curve, Curve::Line { magnitude: 3., .. }))
    );
}
#[test]
fn original_raw_loops_and_resource_and_cancellation_limits_are_checked_first() {
    let text = text();
    for tolerance in [0., f64::NAN] {
        assert!(read_phillips_screw_boundary(&text, tolerance, 1_000_000).is_err());
    }
    assert_eq!(
        read_phillips_screw_boundary(&text, 1e-7, 0).unwrap_err(),
        PhillipsReadError::Limit
    );
    assert_eq!(
        read_phillips_screw_boundary_with_interrupt(&text, 1e-7, 1_000_000, Some(&|| true))
            .unwrap_err(),
        PhillipsReadError::Cancelled
    );
    let edge = text
        .lines()
        .find(|l| l.contains("=ORIENTED_EDGE("))
        .unwrap();
    let replacement = if edge.ends_with(".T.);") {
        edge.replace(".T.);", ".F.);")
    } else {
        edge.replace(".F.);", ".T.);")
    };
    assert!(
        read_phillips_screw_boundary(&text.replacen(edge, &replacement, 1), 1e-7, 1_000_000)
            .is_err()
    );
    let top = read_phillips_screw_boundary(&text, 1e-7, 1_000_000)
        .unwrap()
        .mouth()
        .material_face;
    let face = text
        .lines()
        .find(|l| l.starts_with(&format!("#{top}=")))
        .unwrap();
    assert_eq!(
        read_phillips_screw_boundary(
            &text.replacen(face, &face.replace(".T.);", ".F.);"), 1),
            1e-7,
            1_000_000
        )
        .unwrap_err(),
        PhillipsReadError::Unsupported
    );
    let calls = Cell::new(0);
    assert_eq!(
        read_phillips_screw_boundary_with_interrupt(
            &text,
            1e-7,
            1_000_000,
            Some(&|| {
                let n = calls.get();
                calls.set(n + 1);
                n >= 60
            })
        )
        .unwrap_err(),
        PhillipsReadError::Cancelled
    );
    assert!(calls.get() > 60);
}
fn with_shell<R>(text: &str, f: impl FnOnce(&mut Shell<'_, '_>) -> R) -> R {
    let entities = parse_entities(text, 1e-7).unwrap();
    f(&mut Shell::new(&entities, 1e-7, true))
}
fn changed_entity(text: &str, id: u32, body: &str) -> String {
    let old = text
        .lines()
        .find(|l| l.starts_with(&format!("#{id}=")))
        .unwrap();
    text.replacen(old, &format!("#{id}={body}"), 1)
}
#[test]
fn reversed_literal_edge_parameterization_preserves_the_same_finite_component() {
    let text = text();
    let b = read_phillips_screw_boundary(&text, 1e-7, 1_000_000).unwrap();
    let edge = b
        .edges()
        .iter()
        .find(|e| matches!(e.curve, Curve::CubicIntersection { .. }))
        .unwrap();
    let [a, c] = edge.vertices.map(|i| b.vertices()[i].entity);
    let mut changed = changed_entity(
        &text,
        edge.entity,
        &format!("EDGE_CURVE('',#{c},#{a},#{},.F.);", edge.geometry),
    );
    for line in text
        .lines()
        .filter(|l| l.contains("=ORIENTED_EDGE(") && l.contains(&format!(",#{},", edge.entity)))
    {
        let flipped = if line.ends_with(".T.);") {
            line.replace(".T.);", ".F.);")
        } else {
            line.replace(".F.);", ".T.);")
        };
        changed = changed.replacen(line, &flipped, 1);
    }
    let actual = read_phillips_screw_boundary(&changed, 1e-7, 1_000_000).unwrap();
    let actual = actual
        .edges()
        .iter()
        .find(|e| e.entity == edge.entity)
        .unwrap();
    assert!(!actual.same_sense);
    assert_eq!(actual.vertices, [edge.vertices[1], edge.vertices[0]]);
}
#[test]
fn primitive_readers_reject_malformed_actual_wire_and_vertex_records() {
    let text = text();
    let b = read_phillips_screw_boundary(&text, 1e-7, 1_000_000).unwrap();
    let mouth = b.mouth();
    let loop_bad = changed_entity(&text, mouth.material_loop.loop_entity, "EDGE_LOOP('');");
    assert!(with_shell(&loop_bad, |s| material_loop(s, mouth.material_bound)).is_err());
    let oriented_bad = changed_entity(
        &text,
        mouth.material_loop.uses[0].entity,
        "ORIENTED_EDGE('');",
    );
    assert!(with_shell(&oriented_bad, |s| material_loop(s, mouth.material_bound)).is_err());
    let views: Vec<_> = b
        .faces()
        .iter()
        .map(|f| SourceFaceView {
            surface: f.surface_entity,
            forward: f.cavity_forward,
            uses: f
                .uses
                .iter()
                .map(|&(e, s)| (b.edges()[e].entity, s))
                .collect(),
        })
        .collect();
    let entities: Vec<_> = b
        .faces()
        .iter()
        .map(|f| {
            f.material_entity
                .map(|id| (id, f.material_forward, f.material_loop.clone().unwrap()))
        })
        .collect();
    let run = |text: &str, views: &[SourceFaceView]| {
        with_shell(text, |s| {
            incidence::read(
                s,
                mouth.clone(),
                views,
                &entities,
                &mut Work {
                    remaining: 1_000_000,
                    interrupt: None,
                    cancelled: Cell::new(false),
                },
            )
        })
    };
    let edge_bad = changed_entity(&text, b.edges()[0].entity, "EDGE_CURVE('');");
    assert!(run(&edge_bad, &views).is_err());
    let vertex_bad = changed_entity(&text, b.vertices()[0].entity, "VERTEX_POINT('');");
    assert!(run(&vertex_bad, &views).is_err());
    let mut missing = views;
    missing[0].uses.pop();
    assert!(run(&text, &missing).is_err());
}
#[test]
fn malformed_literal_line_and_ellipse_are_not_hidden_by_a_valid_cell_complex() {
    let text = text();
    let b = read_phillips_screw_boundary(&text, 1e-7, 1_000_000).unwrap();
    let direction = b
        .edges()
        .iter()
        .find_map(|e| match e.curve {
            Curve::Line {
                direction_entity, ..
            } => Some(direction_entity),
            _ => None,
        })
        .unwrap();
    let changed = changed_entity(&text, direction, "DIRECTION('',(0.,0.,0.));");
    assert!(read_phillips_screw_boundary(&changed, 1e-7, 1_000_000).is_err());
    let ellipse = b
        .edges()
        .iter()
        .find(|e| matches!(e.curve, Curve::Ellipse { .. }))
        .unwrap();
    let line = text
        .lines()
        .find(|l| l.starts_with(&format!("#{}=ELLIPSE", ellipse.geometry)))
        .unwrap();
    let fields: Vec<_> = line.split(',').collect();
    let changed = changed_entity(
        &text,
        ellipse.geometry,
        &format!("ELLIPSE('',{},-1.,{}", fields[1], fields[3]),
    );
    assert!(read_phillips_screw_boundary(&changed, 1e-7, 1_000_000).is_err());
}
#[test]
fn source_free_cut_search_refuses_unpaired_holes_partial_cuts_and_oversized_components() {
    let face = |i, uses| MaterialFace {
        entity: i,
        surface: 0,
        forward: true,
        bounds: vec![],
        uses,
    };
    let seam = vec![face(0, vec![(7, true), (7, false)])];
    let pairs = owners(&seam).unwrap();
    assert_eq!(pairs[&7], [0, 0]);
    let mut work = Work {
        remaining: 10_000,
        interrupt: None,
        cancelled: Cell::new(false),
    };
    assert!(component(&seam, &pairs, 0, &[(7, true)], &mut work).is_err());
    assert!(owners(&[face(0, vec![(7, true)])]).is_err());
    let triangle = vec![
        face(0, vec![(0, false), (2, true)]),
        face(1, vec![(0, true), (1, true)]),
        face(2, vec![(1, false), (2, false)]),
    ];
    assert!(
        component(
            &triangle,
            &owners(&triangle).unwrap(),
            0,
            &[(0, false)],
            &mut work
        )
        .unwrap()
        .is_none()
    );
    let mut chain = vec![face(0, vec![(0, true)])];
    for i in 1..30 {
        let mut uses = vec![(i - 1, false)];
        if i < 29 {
            uses.push((i, true));
        }
        chain.push(face(i, uses));
    }
    assert!(
        component(&chain, &owners(&chain).unwrap(), 0, &[(0, true)], &mut work)
            .unwrap()
            .is_none()
    );
}
