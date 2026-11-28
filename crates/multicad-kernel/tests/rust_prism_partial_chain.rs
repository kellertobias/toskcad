//! CAD-421: fillets and chamfers on part of a tangent edge chain of an
//! extruded profile.
//!
//! The reported bodies are an L whose inner corner is rounded in the sketch
//! and a plate whose rectangular pocket has rounded corners. On a cap, the
//! straight edges and the rounded corner of the inner side form one smooth
//! chain. Selecting only part of the chain must still build.
//!
//! Behaviour: the blend stops where the selected edge ends. A flat end face,
//! perpendicular to the edge there, closes it, and the unselected tangent
//! neighbour keeps its sharp edge.
//!
//! Every volume is checked against an independent slice oracle written in
//! this file. A horizontal section of a cap-treated prism at blend inset
//! `delta` is each profile ring with every treated edge offset inward by
//! `delta`. At a sharp vertex, consecutive offset curves meet at their
//! intersection nearest the vertex. At a tangent vertex each edge ends on
//! the common normal at its own offset, and a straight chord joins the two
//! ends. The oracle takes the exact area of that section by Green's theorem
//! on lines and arcs, then integrates it over the blend height with
//! Gauss-Legendre quadrature.
use multicad_core::{document::EntityId, references::*};
use multicad_kernel::*;
use multicad_precision::PortableMath as _;
use std::f64::consts::{FRAC_PI_2, PI};

#[derive(Clone, Copy, Debug)]
enum Seg {
    Line {
        a: [f64; 2],
        b: [f64; 2],
    },
    Arc {
        c: [f64; 2],
        a: [f64; 2],
        b: [f64; 2],
        cw: bool,
    },
}

impl Seg {
    fn start(self) -> [f64; 2] {
        match self {
            Seg::Line { a, .. } | Seg::Arc { a, .. } => a,
        }
    }
    fn end(self) -> [f64; 2] {
        match self {
            Seg::Line { b, .. } | Seg::Arc { b, .. } => b,
        }
    }
    fn curved(self) -> CurvedEdge2d {
        match self {
            Seg::Line { a, b } => CurvedEdge2d::Line { a, b },
            Seg::Arc { c, a, b, cw } => CurvedEdge2d::Arc {
                center: c,
                a,
                b,
                clockwise: cw,
            },
        }
    }
    /// Unit tangent along the ring at the point `p` of the edge.
    fn tangent(self, p: [f64; 2]) -> [f64; 2] {
        match self {
            Seg::Line { a, b } => unit(sub(b, a)),
            Seg::Arc { c, cw, .. } => {
                let r = unit(sub(p, c));
                if cw { [r[1], -r[0]] } else { [-r[1], r[0]] }
            }
        }
    }
    /// Unit normal into the material (left of the ring direction).
    fn inward(self, p: [f64; 2]) -> [f64; 2] {
        let t = self.tangent(p);
        [-t[1], t[0]]
    }
}

fn sub(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] - b[0], a[1] - b[1]]
}
fn len(a: [f64; 2]) -> f64 {
    a[0].portable_hypot(a[1])
}
fn unit(a: [f64; 2]) -> [f64; 2] {
    [a[0] / len(a), a[1] / len(a)]
}
fn cross(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}
fn shifted(p: [f64; 2], d: [f64; 2], s: f64) -> [f64; 2] {
    [p[0] + s * d[0], p[1] + s * d[1]]
}

fn line(a: [f64; 2], b: [f64; 2]) -> Seg {
    Seg::Line { a, b }
}

/// An L of 20 by 20 with arms 8 wide; the inner corner is rounded with
/// radius 3 in the sketch. Edges 2 (straight), 3 (rounded corner) and 4
/// (straight) form the inner tangent chain.
fn rounded_l() -> Vec<Seg> {
    vec![
        line([0., 0.], [20., 0.]),
        line([20., 0.], [20., 8.]),
        line([20., 8.], [11., 8.]),
        Seg::Arc {
            c: [11., 11.],
            a: [11., 8.],
            b: [8., 11.],
            cw: true,
        },
        line([8., 11.], [8., 20.]),
        line([8., 20.], [0., 20.]),
        line([0., 20.], [0., 0.]),
    ]
}

/// A 30 by 20 plate and its pocket from (10, 6) to (20, 14), whose corners
/// are rounded with radius 2. The pocket ring runs clockwise: edge 4 is its
/// bottom side, edge 5 the rounded corner after it.
fn pocket_plate() -> Vec<Vec<Seg>> {
    let outer = vec![
        line([0., 0.], [30., 0.]),
        line([30., 0.], [30., 20.]),
        line([30., 20.], [0., 20.]),
        line([0., 20.], [0., 0.]),
    ];
    let r = 2.;
    let (x0, x1, y0, y1) = (10., 20., 6., 14.);
    let arc = |c: [f64; 2], a: [f64; 2], b: [f64; 2]| Seg::Arc { c, a, b, cw: true };
    let pocket = vec![
        line([x1 - r, y0], [x0 + r, y0]),
        arc([x0 + r, y0 + r], [x0 + r, y0], [x0, y0 + r]),
        line([x0, y0 + r], [x0, y1 - r]),
        arc([x0 + r, y1 - r], [x0, y1 - r], [x0 + r, y1]),
        line([x0 + r, y1], [x1 - r, y1]),
        arc([x1 - r, y1 - r], [x1 - r, y1], [x1, y1 - r]),
        line([x1, y1 - r], [x1, y0 + r]),
        arc([x1 - r, y0 + r], [x1, y0 + r], [x1 - r, y0]),
    ];
    vec![outer, pocket]
}

/// A 24 by 16 plate whose four corners are rounded with radius 3 in the
/// sketch: a smooth outer chain of four straight edges and four arcs.
fn rounded_plate() -> Vec<Seg> {
    let r = 3.;
    let (w, d) = (24., 16.);
    let arc = |c: [f64; 2], a: [f64; 2], b: [f64; 2]| Seg::Arc { c, a, b, cw: false };
    vec![
        line([r, 0.], [w - r, 0.]),
        arc([w - r, r], [w - r, 0.], [w, r]),
        line([w, r], [w, d - r]),
        arc([w - r, d - r], [w, d - r], [w - r, d]),
        line([w - r, d], [r, d]),
        arc([r, d - r], [r, d], [0., d - r]),
        line([0., d - r], [0., r]),
        arc([r, r], [0., r], [r, 0.]),
    ]
}

/// Inward offset curve of one profile edge: a line (point, direction) or a
/// circle (center, radius).
#[derive(Clone, Copy)]
enum Offset {
    Line([f64; 2], [f64; 2]),
    Circle([f64; 2], f64),
}

fn offset(seg: Seg, delta: f64) -> Offset {
    match seg {
        Seg::Line { a, b } => Offset::Line(shifted(a, seg.inward(a), delta), sub(b, a)),
        Seg::Arc { c, a, cw, .. } => {
            let r = len(sub(a, c));
            Offset::Circle(c, if cw { r + delta } else { r - delta })
        }
    }
}

/// Intersections of two offset curves.
fn meet(p: Offset, q: Offset) -> Vec<[f64; 2]> {
    match (p, q) {
        (Offset::Line(a, d), Offset::Line(b, e)) => {
            let det = cross(d, e);
            if det.abs() < 1e-14 {
                return vec![];
            }
            vec![shifted(a, d, cross(sub(b, a), e) / det)]
        }
        (Offset::Line(a, d), Offset::Circle(c, r)) | (Offset::Circle(c, r), Offset::Line(a, d)) => {
            let u = unit(d);
            let w = sub(a, c);
            let along = w[0] * u[0] + w[1] * u[1];
            let disc = along * along - (w[0] * w[0] + w[1] * w[1] - r * r);
            if disc < 0. {
                return vec![];
            }
            [-1., 1.]
                .map(|s| shifted(a, u, -along + s * disc.sqrt()))
                .to_vec()
        }
        (Offset::Circle(c0, r0), Offset::Circle(c1, r1)) => {
            let d = sub(c1, c0);
            let l = len(d);
            let a = (r0 * r0 - r1 * r1 + l * l) / (2. * l);
            let h2 = r0 * r0 - a * a;
            if h2 < 0. {
                return vec![];
            }
            let m = shifted(c0, d, a / l);
            [-1., 1.]
                .map(|s| {
                    let h = s * h2.sqrt() / l;
                    [m[0] - h * d[1], m[1] + h * d[0]]
                })
                .to_vec()
        }
    }
}

/// Where the offset of `seg` by `own` ends at the vertex it shares with
/// `other` (offset by `theirs`).
fn end_point(seg: Seg, other: Seg, vertex: [f64; 2], own: f64, theirs: f64) -> [f64; 2] {
    let tangent = cross(seg.tangent(vertex), other.tangent(vertex)).abs() < 1e-9;
    if tangent || (own == 0. && theirs == 0.) {
        return shifted(vertex, seg.inward(vertex), own);
    }
    *meet(offset(seg, own), offset(other, theirs))
        .iter()
        .min_by(|a, b| len(sub(**a, vertex)).total_cmp(&len(sub(**b, vertex))))
        .expect("offsets meet")
}

/// Area of the profile rings with each edge offset inward by its delta.
fn section_area(rings: &[Vec<Seg>], deltas: &[Vec<f64>]) -> f64 {
    let mut twice = 0.;
    for (ring, deltas) in rings.iter().zip(deltas) {
        let n = ring.len();
        let ends: Vec<[[f64; 2]; 2]> = (0..n)
            .map(|k| {
                let (previous, next) = ((k + n - 1) % n, (k + 1) % n);
                [
                    end_point(
                        ring[k],
                        ring[previous],
                        ring[k].start(),
                        deltas[k],
                        deltas[previous],
                    ),
                    end_point(ring[k], ring[next], ring[k].end(), deltas[k], deltas[next]),
                ]
            })
            .collect();
        for k in 0..n {
            let [p, q] = ends[k];
            twice += match (ring[k], offset(ring[k], deltas[k])) {
                (Seg::Arc { cw, .. }, Offset::Circle(c, r)) => {
                    let angle = |x: [f64; 2]| (x[1] - c[1]).portable_atan2(x[0] - c[0]);
                    let mut sweep = angle(q) - angle(p);
                    sweep = if cw {
                        -(-sweep).rem_euclid(2. * PI)
                    } else {
                        sweep.rem_euclid(2. * PI)
                    };
                    r * r * sweep + c[0] * (q[1] - p[1]) - c[1] * (q[0] - p[0])
                }
                _ => cross(p, q),
            };
            // Chord to the start of the next edge: empty unless one side of
            // a tangent vertex is treated.
            twice += cross(q, ends[(k + 1) % n][0]);
        }
    }
    twice / 2.
}

fn gauss_legendre(n: usize) -> Vec<(f64, f64)> {
    (1..=n)
        .map(|i| {
            let mut x = (PI * (i as f64 - 0.25) / (n as f64 + 0.5)).portable_cos();
            loop {
                let (mut p0, mut p1) = (1., x);
                for k in 2..=n {
                    let k = k as f64;
                    (p0, p1) = (p1, ((2. * k - 1.) * x * p1 - (k - 1.) * p0) / k);
                }
                let dp = n as f64 * (x * p1 - p0) / (x * x - 1.);
                let step = p1 / dp;
                x -= step;
                if step.abs() < 1e-16 {
                    return (x, 2. / ((1. - x * x) * dp * dp));
                }
            }
        })
        .collect()
}

#[derive(Clone, Copy, PartialEq)]
enum Kind {
    Fillet,
    Chamfer,
}

/// Oracle volume. `top` and `bottom` list the treated `(ring, edge)` pairs.
fn oracle(
    rings: &[Vec<Seg>],
    height: f64,
    size: f64,
    kind: Kind,
    top: &[(usize, usize)],
    bottom: &[(usize, usize)],
) -> f64 {
    let deltas = |marked: &[(usize, usize)], delta: f64| -> Vec<Vec<f64>> {
        rings
            .iter()
            .enumerate()
            .map(|(r, ring)| {
                (0..ring.len())
                    .map(|k| if marked.contains(&(r, k)) { delta } else { 0. })
                    .collect()
            })
            .collect()
    };
    let base = section_area(rings, &deltas(&[], 0.));
    let band = |marked: &[(usize, usize)]| {
        if marked.is_empty() {
            return base * size;
        }
        gauss_legendre(64)
            .iter()
            .map(|&(x, w)| {
                let (delta, weight) = match kind {
                    Kind::Fillet => {
                        let theta = FRAC_PI_2 * (x + 1.) / 2.;
                        (
                            size * (1. - theta.portable_cos()),
                            size * theta.portable_cos() * FRAC_PI_2 / 2.,
                        )
                    }
                    Kind::Chamfer => (size * (x + 1.) / 2., size / 2.),
                };
                w * weight * section_area(rings, &deltas(marked, delta))
            })
            .sum()
    };
    base * (height - 2. * size) + band(top) + band(bottom)
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-7 * expected.abs().max(1.),
        "{actual} != {expected} (difference {})",
        actual - expected
    );
}

fn prism(rings: &[Vec<Seg>], height: f64) -> (Solid, EntityId) {
    let id = EntityId::new();
    let count = rings.iter().map(Vec::len).sum::<usize>();
    let origins: Vec<_> = (0..count)
        .map(|i| Origin {
            source: id,
            role: format!("edge/{i}"),
        })
        .collect();
    let curved = |ring: &Vec<Seg>| ring.iter().map(|s| s.curved()).collect::<Vec<_>>();
    let solid = Solid::rust_extrude_normal_profile_named(
        &CurvedProfile {
            frame: PlaneFrame::default(),
            outer: curved(&rings[0]),
            holes: rings[1..].iter().map(curved).collect(),
        },
        height,
        id,
        &origins,
        Default::default(),
    )
    .unwrap();
    (solid, id)
}

fn on_segment(seg: Seg, p: [f64; 3]) -> bool {
    let q = [p[0], p[1]];
    match seg {
        Seg::Line { a, b } => {
            let d = sub(b, a);
            let t = ((q[0] - a[0]) * d[0] + (q[1] - a[1]) * d[1]) / (d[0] * d[0] + d[1] * d[1]);
            (0. ..=1.).contains(&t) && cross(d, sub(q, a)).abs() / len(d) < 1e-6
        }
        Seg::Arc { c, a, b, cw } => {
            // On the circle and between the ends in the arc direction.
            let (u, v, w) = (sub(a, c), sub(q, c), sub(b, c));
            let turn = |x: [f64; 2], y: [f64; 2]| {
                let angle = cross(x, y).portable_atan2(x[0] * y[0] + x[1] * y[1]);
                if cw { -angle } else { angle }.rem_euclid(2. * PI)
            };
            (len(v) - len(u)).abs() < 1e-6 && turn(u, v) <= turn(u, w)
        }
    }
}

fn midpoint(geometry: AnalyticGeometry) -> Option<[f64; 3]> {
    match geometry {
        AnalyticGeometry::Line { a_mm, b_mm } => {
            Some(std::array::from_fn(|i| (a_mm[i] + b_mm[i]) / 2.))
        }
        AnalyticGeometry::CircularArc {
            frame,
            radius_mm,
            start_rad,
            sweep_rad,
        } => {
            let t = start_rad + sweep_rad / 2.;
            Some(std::array::from_fn(|i| {
                frame.origin_mm[i]
                    + radius_mm * (t.portable_cos() * frame.u[i] + t.portable_sin() * frame.v[i])
            }))
        }
        _ => None,
    }
}

/// The cap edges at height `z` on profile edges `(ring, edge)`.
fn select(
    solid: &Solid,
    owner: EntityId,
    rings: &[Vec<Seg>],
    picks: &[(f64, usize, usize)],
    feature: EntityId,
) -> Vec<TreatmentEdge> {
    let snapshot = solid.topology_snapshot(owner, 0).unwrap();
    picks
        .iter()
        .enumerate()
        .map(|(n, &(z, r, i))| {
            let found = snapshot.elements.iter().find(|e| {
                e.key.kind == TopologyKind::Edge
                    && e.geometry.is_some_and(|g| {
                        let flat = match g {
                            AnalyticGeometry::Line { a_mm, b_mm } => {
                                (a_mm[2] - z).abs() < 1e-9 && (b_mm[2] - z).abs() < 1e-9
                            }
                            _ => true,
                        };
                        flat && midpoint(g)
                            .is_some_and(|m| (m[2] - z).abs() < 1e-9 && on_segment(rings[r][i], m))
                    })
            });
            TreatmentEdge {
                key: found.expect("selected edge").key,
                origin: Origin {
                    source: feature,
                    role: format!("treated/{n}"),
                },
            }
        })
        .collect()
}

fn treat(
    rings: &[Vec<Seg>],
    height: f64,
    picks: &[(f64, usize, usize)],
    treatment: EdgeTreatment,
) -> Result<Solid, GeometryError> {
    let (solid, owner) = prism(rings, height);
    let feature = EntityId::new();
    let edges = select(&solid, owner, rings, picks, feature);
    solid.treat_edges_named(&edges, treatment, feature, Default::default())
}

fn step_round_trip(solid: &Solid) -> String {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("treated.step");
    solid.write_step(&path).unwrap();
    let text = std::fs::read_to_string(&path).unwrap();
    let imported = Solid::import_rust_step(&text, EntityId::new(), Default::default()).unwrap();
    assert_eq!(imported.solid_count().unwrap(), 1);
    near(imported.volume_mm3().unwrap(), solid.volume_mm3().unwrap());
    text
}

fn closed_mesh(solid: &Solid) {
    let receipt = solid
        .manufacturing_mesh_with_receipt(MeshOptions {
            linear_deflection_mm: 0.01,
            ..Default::default()
        })
        .unwrap();
    assert_eq!(receipt.edge_pairs as usize * 2, receipt.mesh.indices.len());
}

/// Build, check one closed solid against the oracle, and round-trip STEP.
fn check(
    rings: &[Vec<Seg>],
    height: f64,
    kind: Kind,
    size: f64,
    top: &[(usize, usize)],
    bottom: &[(usize, usize)],
) -> String {
    let picks: Vec<_> = top
        .iter()
        .map(|&(r, i)| (height, r, i))
        .chain(bottom.iter().map(|&(r, i)| (0., r, i)))
        .collect();
    let treatment = match kind {
        Kind::Fillet => EdgeTreatment::Fillet { radius_mm: size },
        Kind::Chamfer => EdgeTreatment::Chamfer { distance_mm: size },
    };
    let result = treat(rings, height, &picks, treatment).unwrap();
    assert_eq!(result.solid_count().unwrap(), 1);
    near(
        result.volume_mm3().unwrap(),
        oracle(rings, height, size, kind, top, bottom),
    );
    closed_mesh(&result);
    step_round_trip(&result)
}

#[test]
fn reported_l_single_inner_straight_edge_fillet_ends_in_flat_faces() {
    let l = vec![rounded_l()];
    let text = check(&l, 6., Kind::Fillet, 1., &[(0, 2)], &[]);
    // One cylinder blend, one rounded-corner side; the blend stops at the
    // rounded corner with a flat end face instead of running into it.
    assert_eq!(text.matches("CYLINDRICAL_SURFACE").count(), 2);
    assert!(!text.contains("TOROIDAL_SURFACE"));
    // The other straight edge of the chain, alone.
    check(&l, 6., Kind::Fillet, 1., &[(0, 4)], &[]);
}

#[test]
fn reported_l_partial_chains_fillet_and_chamfer_on_both_caps() {
    let l = vec![rounded_l()];
    for kind in [Kind::Fillet, Kind::Chamfer] {
        // One edge, two edges (straight and rounded), the rounded corner
        // alone, and the whole chain.
        for top in [
            &[(0, 2)][..],
            &[(0, 2), (0, 3)],
            &[(0, 3)],
            &[(0, 2), (0, 3), (0, 4)],
        ] {
            check(&l, 6., kind, 1., top, &[]);
            check(&l, 6., kind, 1., &[], top);
        }
        // Different partial chains on the two caps.
        check(&l, 6., kind, 1., &[(0, 2)], &[(0, 3), (0, 4)]);
    }
}

#[test]
fn reported_pocket_single_inner_side_fillet_and_chamfer() {
    let plate = pocket_plate();
    check(&plate, 5., Kind::Fillet, 1., &[(1, 0)], &[]);
    check(&plate, 5., Kind::Chamfer, 1., &[(1, 0)], &[]);
    check(&plate, 5., Kind::Fillet, 1., &[], &[(1, 2)]);
    check(&plate, 5., Kind::Chamfer, 1., &[], &[(1, 2)]);
}

#[test]
fn pocket_partial_and_whole_chains_on_both_caps() {
    let plate = pocket_plate();
    let all: Vec<_> = (0..8).map(|i| (1, i)).collect();
    for kind in [Kind::Fillet, Kind::Chamfer] {
        for top in [&[(1, 0), (1, 1)][..], &[(1, 1)], &[(1, 0), (1, 4)], &all] {
            check(&plate, 5., kind, 1., top, &[]);
            check(&plate, 5., kind, 1., &[], top);
        }
        // Outer and pocket edges together on both caps.
        check(&plate, 5., kind, 0.8, &[(0, 0), (1, 6)], &[(0, 2), (1, 3)]);
    }
}

#[test]
fn convex_rounded_corners_partial_chain_on_both_caps() {
    let plate = vec![rounded_plate()];
    for kind in [Kind::Fillet, Kind::Chamfer] {
        check(&plate, 5., kind, 1., &[(0, 0)], &[]);
        check(&plate, 5., kind, 1., &[(0, 0), (0, 1)], &[(0, 5)]);
        check(&plate, 5., kind, 1., &[(0, 1), (0, 2), (0, 3)], &[(0, 4)]);
    }
}

#[test]
fn unsupported_partial_chains_refuse_in_plain_language() {
    // The rounded corner has the fillet radius: the blend would collapse to
    // a point exactly where it has to end.
    let plate = vec![rounded_plate()];
    let refused = treat(
        &plate,
        5.,
        &[(5., 0, 1)],
        EdgeTreatment::Fillet { radius_mm: 3. },
    )
    .err()
    .expect("refused");
    assert_eq!(refused.kind, GeometryErrorKind::Unsupported);
    assert!(
        refused.message.contains("select the whole tangent chain"),
        "{}",
        refused.message
    );
    assert!(!refused.message.contains("Rust"), "{}", refused.message);
    // Too large for the body height.
    let l = vec![rounded_l()];
    let refused = treat(
        &l,
        1.5,
        &[(1.5, 0, 2)],
        EdgeTreatment::Fillet { radius_mm: 2. },
    )
    .err()
    .expect("refused");
    assert!(
        refused
            .message
            .contains("too large for the height of the body"),
        "{}",
        refused.message
    );
}
