//! CAD-411: fillets and chamfers on selected edges of extruded profiles.
//!
//! The reported body is a circle around a box face corner, clipped by the
//! face to two straight edges and one arc, then extruded as a new body.
//!
//! Every volume is checked against an independent slice oracle written in
//! this file. A horizontal section of a cap-treated prism at depth `delta`
//! below the start of the blend is the profile with each treated edge offset
//! inward by `delta`. Consecutive offset curves meet at their intersection
//! nearest the original vertex. The oracle takes the exact area of that
//! section by Green's theorem on lines and arcs and integrates it over the
//! blend height with Gauss-Legendre quadrature. The kernel instead integrates
//! its trimmed boundary surfaces, so the two share no construction code.
use multicad_core::{document::EntityId, references::*};
use multicad_kernel::{rust_kernel::Rational, *};
use multicad_precision::PortableMath as _;
use std::f64::consts::{FRAC_PI_2, PI};

#[derive(Clone, Copy, Debug)]
enum Seg {
    Line { a: [f64; 2], b: [f64; 2] },
    Arc { c: [f64; 2], a: [f64; 2], b: [f64; 2], cw: bool },
}

impl Seg {
    fn start(self) -> [f64; 2] {
        match self {
            Seg::Line { a, .. } | Seg::Arc { a, .. } => a,
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
}

fn sub(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] - b[0], a[1] - b[1]]
}
fn len(a: [f64; 2]) -> f64 {
    a[0].portable_hypot(a[1])
}
fn cross(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}

/// The test body: the part of the circle about `(2, 3)` with radius 8 that
/// lies inside the face quadrant `x, y >= 0` of a box whose corner is the
/// origin.
fn corner_segment() -> Vec<Seg> {
    let (c, r): ([f64; 2], f64) = ([2., 3.], 8.);
    let x = c[0] + (r * r - c[1] * c[1]).sqrt();
    let y = c[1] + (r * r - c[0] * c[0]).sqrt();
    vec![
        Seg::Line {
            a: [0., 0.],
            b: [x, 0.],
        },
        Seg::Arc {
            c,
            a: [x, 0.],
            b: [0., y],
            cw: false,
        },
        Seg::Line {
            a: [0., y],
            b: [0., 0.],
        },
    ]
}

fn rectangle(w: f64, d: f64) -> Vec<Seg> {
    let p = [[0., 0.], [w, 0.], [w, d], [0., d]];
    (0..4)
        .map(|i| Seg::Line {
            a: p[i],
            b: p[(i + 1) % 4],
        })
        .collect()
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
        Seg::Line { a, b } => {
            let d = sub(b, a);
            let n = [-d[1] / len(d), d[0] / len(d)];
            Offset::Line([a[0] + delta * n[0], a[1] + delta * n[1]], d)
        }
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
            let t = cross(sub(b, a), e) / det;
            vec![[a[0] + t * d[0], a[1] + t * d[1]]]
        }
        (Offset::Line(a, d), Offset::Circle(c, r)) | (Offset::Circle(c, r), Offset::Line(a, d)) => {
            let u = [d[0] / len(d), d[1] / len(d)];
            let w = sub(a, c);
            let along = w[0] * u[0] + w[1] * u[1];
            let disc = along * along - (w[0] * w[0] + w[1] * w[1] - r * r);
            if disc < 0. {
                return vec![];
            }
            [-1., 1.]
                .map(|s| {
                    let t = -along + s * disc.sqrt();
                    [a[0] + t * u[0], a[1] + t * u[1]]
                })
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
            let m = [c0[0] + a * d[0] / l, c0[1] + a * d[1] / l];
            [-1., 1.]
                .map(|s| {
                    let h = s * h2.sqrt();
                    [m[0] - h * d[1] / l, m[1] + h * d[0] / l]
                })
                .to_vec()
        }
    }
}

/// Area of the profile with each edge `i` offset inward by `deltas[i]`.
fn section_area(profile: &[Seg], deltas: &[f64]) -> f64 {
    let n = profile.len();
    let corners: Vec<[f64; 2]> = (0..n)
        .map(|k| {
            let previous = (k + n - 1) % n;
            let vertex = profile[k].start();
            let candidates = meet(
                offset(profile[previous], deltas[previous]),
                offset(profile[k], deltas[k]),
            );
            if candidates.is_empty() {
                // Tangent continuation: both offsets touch at the shifted
                // vertex.
                let step = deltas[k].max(deltas[previous]);
                let normal = match offset(profile[k], 0.) {
                    Offset::Line(_, d) => [-d[1] / len(d), d[0] / len(d)],
                    Offset::Circle(c, _) => {
                        let s = if let Seg::Arc { cw: true, .. } = profile[k] {
                            -1.
                        } else {
                            1.
                        };
                        let r = sub(c, vertex);
                        [s * r[0] / len(r), s * r[1] / len(r)]
                    }
                };
                return [vertex[0] + step * normal[0], vertex[1] + step * normal[1]];
            }
            *candidates
                .iter()
                .min_by(|a, b| len(sub(**a, vertex)).total_cmp(&len(sub(**b, vertex))))
                .unwrap()
        })
        .collect();
    let mut twice = 0.;
    for k in 0..n {
        let [p, q] = [corners[k], corners[(k + 1) % n]];
        twice += match (profile[k], offset(profile[k], deltas[k])) {
            (Seg::Arc { cw, .. }, Offset::Circle(c, r)) => {
                let angle = |x: [f64; 2]| (x[1] - c[1]).portable_atan2(x[0] - c[0]);
                let mut sweep = angle(q) - angle(p);
                if cw {
                    sweep = -(-sweep).rem_euclid(2. * PI);
                } else {
                    sweep = sweep.rem_euclid(2. * PI);
                }
                if r < 1e-12 {
                    sweep = 0.;
                }
                r * r * sweep + c[0] * (q[1] - p[1]) - c[1] * (q[0] - p[0])
            }
            _ => p[0] * q[1] - q[0] * p[1],
        };
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
                    let w = 2. / ((1. - x * x) * dp * dp);
                    return (x, w);
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

/// Oracle volume: `top[i]` / `bottom[i]` mark treated cap edges.
fn oracle(profile: &[Seg], height: f64, size: f64, kind: Kind, top: &[bool], bottom: &[bool]) -> f64 {
    let base = section_area(profile, &vec![0.; profile.len()]);
    let band = |mask: &[bool]| {
        if !mask.contains(&true) {
            return base * size;
        }
        let nodes = gauss_legendre(64);
        nodes
            .iter()
            .map(|&(x, w)| {
                let (delta, weight) = match kind {
                    Kind::Fillet => {
                        let theta = FRAC_PI_2 * (x + 1.) / 2.;
                        (size * (1. - theta.portable_cos()), size * theta.portable_cos() * FRAC_PI_2 / 2.)
                    }
                    Kind::Chamfer => (size * (x + 1.) / 2., size / 2.),
                };
                let deltas: Vec<_> = mask.iter().map(|&m| if m { delta } else { 0. }).collect();
                w * weight * section_area(profile, &deltas)
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

fn prism(profile: &[Seg], height: f64) -> (Solid, EntityId) {
    let id = EntityId::new();
    let origins: Vec<_> = (0..profile.len())
        .map(|i| Origin {
            source: id,
            role: format!("edge/{i}"),
        })
        .collect();
    let solid = Solid::rust_extrude_normal_profile_named(
        &CurvedProfile {
            frame: PlaneFrame::default(),
            outer: profile.iter().map(|s| s.curved()).collect(),
            holes: vec![],
        },
        height,
        id,
        &origins,
        Default::default(),
    )
    .unwrap();
    (solid, id)
}

/// Cap edges at height `z` whose midpoint lies on profile edge `i`, and the
/// vertical edge through profile vertex `k`.
enum Pick {
    Cap(f64, usize),
    Side(usize),
}

fn on_segment(seg: Seg, p: [f64; 3]) -> bool {
    let q = [p[0], p[1]];
    match seg {
        Seg::Line { a, b } => {
            let d = sub(b, a);
            let t = ((q[0] - a[0]) * d[0] + (q[1] - a[1]) * d[1]) / (d[0] * d[0] + d[1] * d[1]);
            (0. ..=1.).contains(&t) && cross(d, sub(q, a)).abs() / len(d) < 1e-6
        }
        Seg::Arc { c, a, .. } => (len(sub(q, c)) - len(sub(a, c))).abs() < 1e-6,
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
                frame.origin_mm[i] + radius_mm * (t.portable_cos() * frame.u[i] + t.portable_sin() * frame.v[i])
            }))
        }
        _ => None,
    }
}

fn select(
    solid: &Solid,
    owner: EntityId,
    profile: &[Seg],
    picks: &[Pick],
) -> (Vec<TreatmentEdge>, EntityId) {
    let snapshot = solid.topology_snapshot(owner, 0).unwrap();
    let feature = EntityId::new();
    let edges = picks
        .iter()
        .enumerate()
        .map(|(n, pick)| {
            let found = snapshot.elements.iter().find(|e| {
                e.key.kind == TopologyKind::Edge
                    && e.geometry.is_some_and(|g| match (pick, g) {
                        (Pick::Cap(z, i), g) => midpoint(g).is_some_and(|m| {
                            (m[2] - z).abs() < 1e-9
                                && on_segment(profile[*i], m)
                                && match g {
                                    AnalyticGeometry::Line { a_mm, b_mm } => {
                                        (a_mm[2] - z).abs() < 1e-9 && (b_mm[2] - z).abs() < 1e-9
                                    }
                                    _ => true,
                                }
                        }),
                        (Pick::Side(k), AnalyticGeometry::Line { a_mm, b_mm }) => {
                            let v = profile[*k].start();
                            [a_mm, b_mm].iter().all(|p| len(sub([p[0], p[1]], v)) < 1e-9)
                                && (a_mm[2] - b_mm[2]).abs() > 1e-6
                        }
                        _ => false,
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
        .collect();
    (edges, feature)
}

fn treat(
    solid: &Solid,
    owner: EntityId,
    profile: &[Seg],
    picks: &[Pick],
    treatment: EdgeTreatment,
) -> Result<Solid, GeometryError> {
    let (edges, feature) = select(solid, owner, profile, picks);
    solid.treat_edges_named(&edges, treatment, feature, Default::default())
}

fn step_text(solid: &Solid) -> (String, Solid) {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("treated.step");
    solid.write_step(&path).unwrap();
    let text = std::fs::read_to_string(&path).unwrap();
    let imported = Solid::import_rust_step(&text, EntityId::new(), Default::default()).unwrap();
    (text, imported)
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

fn fillet(radius_mm: f64) -> EdgeTreatment {
    EdgeTreatment::Fillet { radius_mm }
}

fn chamfer(distance_mm: f64) -> EdgeTreatment {
    EdgeTreatment::Chamfer { distance_mm }
}

#[test]
fn reported_corner_segment_top_edges_fillet_is_exact_and_round_trips_step() {
    let profile = corner_segment();
    let (solid, owner) = prism(&profile, 5.);
    let all = [Pick::Cap(5., 0), Pick::Cap(5., 1), Pick::Cap(5., 2)];
    let result = treat(&solid, owner, &profile, &all, fillet(1.)).unwrap();
    assert_eq!(result.solid_count().unwrap(), 1);
    let expected = oracle(&profile, 5., 1., Kind::Fillet, &[true; 3], &[false; 3]);
    near(result.volume_mm3().unwrap(), expected);
    closed_mesh(&result);
    let (text, imported) = step_text(&result);
    assert!(text.contains("TOROIDAL_SURFACE"), "arc edge blend is a torus");
    assert_eq!(text.matches("CYLINDRICAL_SURFACE").count(), 3, "two line blends and the arc side");
    assert_eq!(imported.solid_count().unwrap(), 1);
    near(imported.volume_mm3().unwrap(), expected);
}

#[test]
fn reported_corner_segment_top_edges_chamfer_is_exact_and_round_trips_step() {
    let profile = corner_segment();
    let (solid, owner) = prism(&profile, 5.);
    let all = [Pick::Cap(5., 0), Pick::Cap(5., 1), Pick::Cap(5., 2)];
    let result = treat(&solid, owner, &profile, &all, chamfer(1.)).unwrap();
    assert_eq!(result.solid_count().unwrap(), 1);
    let expected = oracle(&profile, 5., 1., Kind::Chamfer, &[true; 3], &[false; 3]);
    near(result.volume_mm3().unwrap(), expected);
    closed_mesh(&result);
    let (text, imported) = step_text(&result);
    assert!(text.contains("CONICAL_SURFACE"), "arc edge chamfer is a cone");
    near(imported.volume_mm3().unwrap(), expected);
}

#[test]
fn corner_segment_single_and_partial_cap_edges() {
    let profile = corner_segment();
    let (solid, owner) = prism(&profile, 5.);
    let cases: [(&[usize], f64); 4] = [(&[1], 1.), (&[0], 1.5), (&[2, 1], 0.75), (&[0, 2], 1.)];
    for (edges, size) in cases {
        let picks: Vec<_> = edges.iter().map(|&i| Pick::Cap(5., i)).collect();
        let mut mask = [false; 3];
        edges.iter().for_each(|&i| mask[i] = true);
        for kind in [Kind::Fillet, Kind::Chamfer] {
            let treatment = match kind {
                Kind::Fillet => fillet(size),
                Kind::Chamfer => chamfer(size),
            };
            let result = treat(&solid, owner, &profile, &picks, treatment)
                .unwrap_or_else(|e| panic!("{edges:?}: {e}"));
            let expected = oracle(&profile, 5., size, kind, &mask, &[false; 3]);
            near(result.volume_mm3().unwrap(), expected);
            let (_, imported) = step_text(&result);
            near(imported.volume_mm3().unwrap(), expected);
        }
    }
}

#[test]
fn corner_segment_both_caps_and_bottom_only() {
    let profile = corner_segment();
    let (solid, owner) = prism(&profile, 5.);
    let both: Vec<_> = (0..3)
        .flat_map(|i| [Pick::Cap(5., i), Pick::Cap(0., i)])
        .collect();
    let result = treat(&solid, owner, &profile, &both, fillet(1.)).unwrap();
    near(
        result.volume_mm3().unwrap(),
        oracle(&profile, 5., 1., Kind::Fillet, &[true; 3], &[true; 3]),
    );
    closed_mesh(&result);
    let bottom = treat(&solid, owner, &profile, &[Pick::Cap(0., 1)], chamfer(0.5)).unwrap();
    near(
        bottom.volume_mm3().unwrap(),
        oracle(&profile, 5., 0.5, Kind::Chamfer, &[false; 3], &[false, true, false]),
    );
}

#[test]
fn vertical_segment_edge_rounds_the_profile_exactly() {
    let profile = corner_segment();
    let (solid, owner) = prism(&profile, 5.);
    // The 90 degree corner at the box corner.
    let result = treat(&solid, owner, &profile, &[Pick::Side(0)], fillet(1.)).unwrap();
    let base = section_area(&profile, &[0.; 3]) * 5.;
    near(result.volume_mm3().unwrap(), base - (1. - PI / 4.) * 5.);
    let cut = treat(&solid, owner, &profile, &[Pick::Side(0)], chamfer(1.)).unwrap();
    near(cut.volume_mm3().unwrap(), base - 0.5 * 5.);
}

fn box_solid(size: [i64; 3]) -> (Solid, EntityId) {
    let id = EntityId::new();
    let solid = Solid::rust_box(
        id.uuid_bits(),
        size.map(|s| Rational::new(s, 1).unwrap()),
        PlaneFrame::default(),
        Default::default(),
    )
    .unwrap();
    (solid, id)
}

#[test]
fn box_single_partial_and_corner_edge_sets() {
    let (solid, owner) = box_solid([30, 20, 10]);
    let rect = rectangle(30., 20.);
    let cases: [&[usize]; 4] = [&[0], &[0, 1], &[0, 1, 2], &[0, 1, 2, 3]];
    for edges in cases {
        let picks: Vec<_> = edges.iter().map(|&i| Pick::Cap(10., i)).collect();
        let mut mask = [false; 4];
        edges.iter().for_each(|&i| mask[i] = true);
        for (kind, treatment) in [(Kind::Fillet, fillet(2.)), (Kind::Chamfer, chamfer(2.))] {
            let result = treat(&solid, owner, &rect, &picks, treatment)
                .unwrap_or_else(|e| panic!("{edges:?}: {e}"));
            assert_eq!(result.solid_count().unwrap(), 1);
            let expected = oracle(&rect, 10., 2., kind, &mask, &[false; 4]);
            near(result.volume_mm3().unwrap(), expected);
            closed_mesh(&result);
            let (_, imported) = step_text(&result);
            near(imported.volume_mm3().unwrap(), expected);
        }
    }
}

#[test]
fn box_corner_with_vertical_and_both_top_edges_meets_in_a_sphere() {
    let (solid, owner) = box_solid([30, 20, 10]);
    let rect = rectangle(30., 20.);
    // Vertex 1 is (30, 0): its vertical edge and the two top edges there.
    let picks = [Pick::Side(1), Pick::Cap(10., 0), Pick::Cap(10., 1)];
    let result = treat(&solid, owner, &rect, &picks, fillet(2.)).unwrap();
    let rounded = vec![
        Seg::Line {
            a: [0., 0.],
            b: [28., 0.],
        },
        Seg::Arc {
            c: [28., 2.],
            a: [28., 0.],
            b: [30., 2.],
            cw: false,
        },
        Seg::Line {
            a: [30., 2.],
            b: [30., 20.],
        },
        Seg::Line {
            a: [30., 20.],
            b: [0., 20.],
        },
        Seg::Line {
            a: [0., 20.],
            b: [0., 0.],
        },
    ];
    let expected = oracle(
        &rounded,
        10.,
        2.,
        Kind::Fillet,
        &[true, true, true, false, false],
        &[false; 5],
    );
    near(result.volume_mm3().unwrap(), expected);
    let (text, imported) = step_text(&result);
    assert!(text.contains("SPHERICAL_SURFACE"));
    near(imported.volume_mm3().unwrap(), expected);
}

#[test]
fn unsupported_selections_refuse_in_plain_language() {
    let profile = corner_segment();
    let (solid, owner) = prism(&profile, 5.);
    let all = [Pick::Cap(5., 0), Pick::Cap(5., 1), Pick::Cap(5., 2)];
    // Larger than the body is tall.
    let error = treat(&solid, owner, &profile, &all, fillet(6.)).err().unwrap();
    assert!(error.to_string().contains("too large"), "{error}");
    // An L-shaped profile has a concave corner between the two top edges.
    let l = {
        let p = [[0., 0.], [20., 0.], [20., 8.], [8., 8.], [8., 20.], [0., 20.]];
        (0..6)
            .map(|i| Seg::Line {
                a: p[i],
                b: p[(i + 1) % 6],
            })
            .collect::<Vec<_>>()
    };
    let (solid, owner) = prism(&l, 5.);
    let error = treat(&solid, owner, &l, &[Pick::Cap(5., 2), Pick::Cap(5., 3)], fillet(1.))
        .err()
        .unwrap();
    assert_eq!(error.kind, GeometryErrorKind::Unsupported);
    assert!(error.to_string().contains("concave"), "{error}");
    assert!(!error.to_string().contains("Rust"), "{error}");
}
