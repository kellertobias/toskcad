use super::*;
use crate::{
    Circle, Rational, chamfered_bores::BoreChamfer, polygon::PolygonPrism,
    step::write_chamfered_polygon_prism,
};

const Z: Point = [0., 0., 1.];
const ORIGIN: Point = [0.; 3];
const TOLERANCE: f64 = 1e-7;

fn fixture() -> String {
    let n = |value| Rational::new(value, 1).unwrap();
    let ring = [[0, 0], [18, 0], [18, 12], [0, 12]]
        .map(|point| point.map(n))
        .to_vec();
    let through = Circle::new([n(9), n(6)], n(2)).unwrap();
    let prism = PolygonPrism::with_circular_holes(vec![ring], vec![through], n(6)).unwrap();
    write_chamfered_polygon_prism(
        &prism,
        &[false, true].map(|end| BoreChamfer {
            hole: 0,
            end,
            size_mm: 1.,
        }),
        ORIGIN,
        [1., 0., 0.],
        [0., 1., 0.],
        false,
    )
    .unwrap()
}

#[test]
fn countersink_import_requires_all_walls_in_their_geometric_sequence() {
    let source = fixture();
    let entities = parse_entities(&source, TOLERANCE).unwrap();
    let mut shell = Shell::new(&entities, TOLERANCE, true);
    let mut faces = shell.validated_faces().unwrap();
    let mut edges = shell.indexed_edges(faces.len());
    let (rim, mouth) = edges[0]
        .iter()
        .find_map(|&edge| match shell.edges[&edge].curve {
            Curve::Circle(circle) => Some((edge, circle)),
            _ => None,
        })
        .unwrap();
    let cone = shell.adjacent_face(rim, 0).unwrap();
    let read_bore = |axis, caps, mut covered| {
        shell.conical_seated_bore(&faces, &edges, caps, (rim, mouth), axis, &mut covered)
    };
    assert!(read_bore(Z, [0, 1], BTreeSet::new()).is_some());
    assert!(read_bore(Z, [0, 1], BTreeSet::from([cone])).is_none());
    assert!(read_bore([1., 0., 0.], [0, 1], BTreeSet::new()).is_none());
    // A cone after the cylinder must terminate at the upper stock cap.
    assert!(read_bore(Z, [0, usize::MAX], BTreeSet::new()).is_none());

    // Rims may appear later than the seam in another writer's edge order.
    edges[cone].sort_by_key(|edge| !matches!(shell.edges[edge].curve, Curve::Line));
    assert!(
        shell
            .match_countersunk_polygon(&faces, &edges, 0, ORIGIN, Z)
            .is_some()
    );
    let far = *edges[cone]
        .iter()
        .find(|&&edge| edge != rim && matches!(shell.edges[&edge].curve, Curve::Circle(_)))
        .unwrap();
    let wall_edges = edges[cone].clone();
    edges[cone].retain(|&edge| edge != far);
    assert!(
        shell
            .match_countersunk_polygon(&faces, &edges, 0, ORIGIN, Z)
            .is_none()
    );
    edges[cone] = wall_edges;
    for edge in [rim, far] {
        let uses = shell.edges.get_mut(&edge).unwrap();
        let saved = uses.uses.clone();
        uses.uses
            .retain(|use_| use_.face == if edge == rim { 0 } else { cone });
        assert!(
            shell
                .match_countersunk_polygon(&faces, &edges, 0, ORIGIN, Z)
                .is_none()
        );
        shell.edges.get_mut(&edge).unwrap().uses = saved;
    }

    // A disjoint face can never be silently ignored by recognition.
    faces.push(Face::Planar {
        origin: ORIGIN,
        normal: Z,
        area: 1.,
    });
    assert!(
        shell
            .match_countersunk_polygon(&faces, &edges, 0, ORIGIN, Z)
            .is_none()
    );
    faces.pop();

    // Changing a serialized rim without changing its cone support is invalid
    // before recognition; non-45-degree *valid* cones are admitted separately.
    let needle = format!(",{:.16E})", mouth.radius_mm);
    let line = source
        .lines()
        .find(|l| l.contains("=CIRCLE(") && l.contains(&needle))
        .unwrap();
    let changed = line.replacen(&needle, &format!(",{:.16E})", mouth.radius_mm + 0.25), 1);
    let corrupt = source.replacen(line, &changed, 1);
    assert!(read_import_solid(&corrupt, TOLERANCE).is_err());

    // Reversing the surface turns the inward wall into a positive feature.
    let Face::Conical(wall) = &mut faces[cone] else {
        unreachable!()
    };
    wall.outward = true;
    assert!(shell.countersunk_polygon(&faces).is_none());
}

#[test]
fn raw_candidates_with_a_missing_far_circle_or_an_unclassified_cap_curve_are_refused() {
    let source = fixture();
    let entities = parse_entities(&source, TOLERANCE).unwrap();
    let mut shell = Shell::new(&entities, TOLERANCE, true);
    let faces = shell.validated_faces().unwrap();
    let mut edges = shell.indexed_edges(faces.len());
    let (rim, mouth) = edges[0]
        .iter()
        .find_map(|&id| {
            if let Curve::Circle(c) = shell.edges[&id].curve {
                Some((id, c))
            } else {
                None
            }
        })
        .unwrap();
    let cone = shell.adjacent_face(rim, 0).unwrap();
    let seam = *edges[cone]
        .iter()
        .find(|&&id| matches!(shell.edges[&id].curve, Curve::Line))
        .unwrap();
    let other = *shell
        .edges
        .iter()
        .find(|(id, e)| **id != seam && matches!(e.curve, Curve::Line))
        .unwrap()
        .0;
    let original = edges[cone].clone();
    // These are untrusted reader candidates, not forged accepted solid proofs.
    // Arity alone cannot substitute for an actual second circular rim.
    edges[cone] = vec![rim, seam, other];
    assert!(
        shell
            .conical_seated_bore(
                &faces,
                &edges,
                [0, 1],
                (rim, mouth),
                Z,
                &mut BTreeSet::new()
            )
            .is_none()
    );
    edges[cone] = original;
    let id = shell.edges.keys().max().unwrap() + 1;
    shell.edges.insert(
        id,
        Edge {
            vertices: shell.edges[&rim].vertices,
            curve: Curve::Semicircle(SemicircleEdge {
                center_mm: mouth.center_mm,
                normal: Z,
                radius_mm: mouth.radius_mm,
                radial_mid: [1., 0., 0.],
            }),
            uses: Vec::new(),
        },
    );
    // A partially recognized cap must not discard an unsupported finite trim.
    edges[0].push(id);
    assert!(
        shell
            .match_countersunk_polygon(&faces, &edges, 0, ORIGIN, Z)
            .is_none()
    );
}
