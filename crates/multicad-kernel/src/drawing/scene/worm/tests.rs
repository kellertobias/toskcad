//! The worm members proven against the provider's own membership predicate
//! and surface laws, independently of the scene engine.
use super::*;
use crate::{GeneratedZaWormWheel, RigidTransform, ZaWorm};
use multicad_core::document::EntityId;
use multicad_rust_kernel::worm_gear::Containment;

fn worm(m: f64, length: f64, bore: f64) -> Solid {
    Solid::rust_za_worm_named(
        ZaWorm {
            axial_module_mm: m,
            threaded_length_mm: length * m,
            bore_radius_mm: bore,
            axial_thinning_mm: 0.02 * m,
        },
        EntityId::new(),
        GeometryPolicy::default(),
    )
    .unwrap()
}

fn wheel(m: f64, teeth: u32, bore: f64) -> Solid {
    Solid::rust_generated_za_worm_wheel_named(
        GeneratedZaWormWheel {
            teeth,
            axial_module_mm: m,
            face_width_mm: 1.6 * m,
            bore_radius_mm: bore,
            pitch_thinning_mm: 0.03 * m,
        },
        EntityId::new(),
        GeometryPolicy::default(),
    )
    .unwrap()
}

fn posed(solid: Solid) -> Solid {
    let (s, c) = 0.3_f64.sin_cos();
    solid
        .transformed(RigidTransform {
            rotation: [[c, -s, 0.], [s, c, 0.], [0., 0., 1.]],
            translation_mm: [1.5, -2., 0.75],
        })
        .unwrap()
}

fn solids() -> Vec<(&'static str, Solid)> {
    vec![
        ("worm", worm(1., 20., 0.)),
        ("bored worm", posed(worm(1.5, 16., 2.5))),
        ("wheel", wheel(1., 32, 0.)),
        ("bored wheel", posed(wheel(1.25, 40, 6.))),
    ]
}

fn views() -> Vec<(&'static str, PlaneFrame)> {
    let origin = [0.3, -0.2, 0.1];
    let s = 0.5_f64.sqrt();
    vec![
        (
            "axial",
            PlaneFrame {
                origin,
                u: [1., 0., 0.],
                v: [0., 1., 0.],
            },
        ),
        (
            "reversed",
            PlaneFrame {
                origin,
                u: [-1., 0., 0.],
                v: [0., 1., 0.],
            },
        ),
        (
            "side",
            PlaneFrame {
                origin,
                u: [0., 1., 0.],
                v: [0., 0., 1.],
            },
        ),
        (
            "top",
            PlaneFrame {
                origin,
                u: [1., 0., 0.],
                v: [0., 0., 1.],
            },
        ),
        (
            "oblique",
            PlaneFrame {
                origin,
                u: [s, -s, 0.],
                v: [s / 3_f64.sqrt(), s / 3_f64.sqrt(), -2. * s / 3_f64.sqrt()],
            },
        ),
        (
            "tilted",
            PlaneFrame {
                origin,
                u: [1., 0., 0.],
                v: [0., 0.96, 0.28],
            },
        ),
    ]
}

fn built(solid: &Solid, frame: PlaneFrame) -> Body {
    body(solid, &View::new(frame)).unwrap().unwrap()
}

fn classify(solid: &Solid, local: V3, tolerance: f64) -> Containment {
    match solid.rust_worm() {
        Some((worm, _)) => worm.classify(local, tolerance),
        None => solid
            .rust_worm_wheel()
            .unwrap()
            .0
            .classify(local, tolerance),
    }
}

struct Random(u64);

impl Random {
    fn next(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 11) as f64 / (1_u64 << 53) as f64
    }
}

#[test]
fn sight_intervals_are_the_membership_predicate_with_roots_on_the_surfaces() {
    let mut random = Random(7);
    for (name, solid) in solids() {
        for (label, frame) in views() {
            let body = built(&solid, frame);
            let inner = &body.0;
            let [lo, hi] = body.bounds();
            let mut inside = 0;
            let mut ends = 0;
            for _ in 0..400 {
                let paper = [0, 1].map(|k| lo[k] + (hi[k] - lo[k]) * random.next());
                let spans = body.sight(paper).unwrap();
                for pair in spans.windows(2) {
                    assert!(pair[0][1] < pair[1][0], "{name} {label} {spans:?}");
                }
                for span in &spans {
                    assert!(span[1] > span[0]);
                    for depth in span {
                        let local = inner.local([paper[0], paper[1], *depth]);
                        assert_eq!(
                            classify(&solid, local, 1e-9),
                            Containment::Boundary,
                            "{name} {label} {paper:?} {depth}"
                        );
                        ends += 1;
                    }
                }
                for _ in 0..24 {
                    let depth = lo[2] - 0.5 + (hi[2] - lo[2] + 1.) * random.next();
                    let local = inner.local([paper[0], paper[1], depth]);
                    let found = spans.iter().any(|[a, b]| *a < depth && depth < *b);
                    let expected = match classify(&solid, local, 1e-6) {
                        Containment::Inside => true,
                        Containment::Outside => false,
                        Containment::Boundary => continue,
                    };
                    assert_eq!(
                        found, expected,
                        "{name} {label} {paper:?} {depth} {spans:?}"
                    );
                    inside += usize::from(found);
                }
            }
            assert!(inside > 200 && ends > 100, "{name} {label} {inside} {ends}");
        }
    }
}

/// Distance from a point to a polyline.
fn distance(point: P3, line: &[P3]) -> f64 {
    line.windows(2)
        .map(|pair| {
            let d = sub(pair[1], pair[0]);
            let s = (dot3(sub(point, pair[0]), d) / dot3(d, d).max(1e-300)).clamp(0., 1.);
            norm(sub(point, add(pair[0], scaled(d, s))))
        })
        .fold(f64::INFINITY, f64::min)
}

#[test]
fn candidates_lie_on_the_boundary_and_outlines_are_complete() {
    for (name, solid) in solids() {
        for (label, frame) in views() {
            let body = built(&solid, frame);
            let inner = &body.0;
            let candidates = body.candidates(0.002).unwrap();
            let mut outlines: Vec<Vec<P3>> = Vec::new();
            for candidate in &candidates {
                let stops = candidate.path.flatten(0.002);
                assert!(stops[0] == 0. && *stops.last().unwrap() == 1.);
                let points: Vec<P3> = stops.iter().map(|&t| candidate.path.point(t)).collect();
                for (i, pair) in stops.windows(2).enumerate() {
                    assert!(pair[1] > pair[0]);
                    for s in [0., 0.37, 1.] {
                        let p = candidate.path.point(pair[0] + (pair[1] - pair[0]) * s);
                        assert_eq!(
                            classify(&solid, inner.local(p), 1e-8),
                            Containment::Boundary,
                            "{name} {label} {:?} {i}",
                            candidate.kind
                        );
                    }
                    // Chords stay within the tolerance of the path.
                    let mid = candidate.path.point(0.5 * (pair[0] + pair[1]));
                    let chord = scaled(add(points[i], points[i + 1]), 0.5);
                    assert!(norm(sub(mid, chord)) <= 0.002, "{name} {label}");
                }
                // Every part is injective in the paper along one paper axis.
                let monotone = |k: usize| {
                    points.windows(2).all(|p| p[1][k] >= p[0][k] - 1e-9)
                        || points.windows(2).all(|p| p[1][k] <= p[0][k] + 1e-9)
                };
                assert!(
                    monotone(0) || monotone(1),
                    "{name} {label} {:?}",
                    candidate.kind
                );
                if candidate.kind == Kind::Outline {
                    outlines.push(points);
                }
            }
            // The exact normal is perpendicular to the sight line along
            // every traced outline.
            let shell = &inner.shell;
            let half = shell.half();
            let mut traced = 0;
            let mut crossings = 0;
            for patch in 0..shell.patches() {
                for path in body.traced(patch, Level::Outline).unwrap() {
                    for i in 0..=20 {
                        let jet = path.jet(f64::from(i) / 20.);
                        let facing = dot3(unit(jet.normal()), inner.toward);
                        assert!(facing.abs() < 1e-9, "{name} {label} {patch} {facing}");
                    }
                    traced += 1;
                }
                // Dense sign changes of `n · d` on the patch all lie on a
                // candidate outline.
                let [nu, nz] = shell.grid(patch).map(|n| 5 * n + 3);
                let facing =
                    |u: f64, z: f64| dot3(unit(shell.eval(patch, u, z).normal()), inner.toward);
                let node = |i: usize, j: usize| {
                    (
                        i as f64 / nu as f64,
                        -half + 2. * half * j as f64 / nz as f64,
                    )
                };
                for i in 0..nu {
                    for j in 0..=nz {
                        let [(ua, z), (ub, _)] = [node(i, j), node(i + 1, j)];
                        let [fa, fb] = [facing(ua, z), facing(ub, z)];
                        if fa.abs() < 1e-7 || fb.abs() < 1e-7 || (fa < 0.) == (fb < 0.) {
                            continue;
                        }
                        let u = root(|u| facing(u, z), ua, ub);
                        let point = inner.view(shell.eval(patch, u, z).p);
                        let gap = outlines
                            .iter()
                            .map(|line| distance(point, line))
                            .fold(f64::INFINITY, f64::min);
                        assert!(gap < 0.004, "{name} {label} {patch} {u} {z} {gap}");
                        crossings += 1;
                    }
                }
            }
            let axial = inner.toward[0].hypot(inner.toward[1]) < 1e-9;
            assert!(axial || (traced > 0 && crossings > 0), "{name} {label}");
        }
    }
}

/// Signed area of oriented pieces by Green's theorem on fine chords.
fn area(pieces: &[Arc<dyn Path>]) -> f64 {
    let mut total = 0.;
    for piece in pieces {
        let points: Vec<P3> = piece
            .flatten(2e-5)
            .into_iter()
            .map(|t| piece.point(t))
            .collect();
        for pair in points.windows(2) {
            total += 0.5 * (pair[0][0] * pair[1][1] - pair[1][0] * pair[0][1]);
        }
    }
    total
}

/// Area of the membership predicate over the plane at a depth: exact
/// transitions along rows, midpoint rule across them. The rows are turned
/// in the paper so that no straight boundary runs along a row.
fn quadrature(solid: &Solid, body: &Body, depth: f64) -> f64 {
    let inner = &body.0;
    let [lo, hi] = body.bounds();
    // The wheel predicate bisects its section at every call.
    let (rows, columns) = if solid.rust_worm().is_some() {
        (1500, 3000)
    } else {
        (400, 800)
    };
    let centre = [0, 1].map(|k| 0.5 * (lo[k] + hi[k]));
    let reach = 0.5 * (hi[0] - lo[0]).hypot(hi[1] - lo[1]) + 1e-3;
    let (sin, cos) = 0.37_f64.sin_cos();
    let inside = |x: f64, y: f64| {
        let paper = [
            centre[0] + cos * x - sin * y,
            centre[1] + sin * x + cos * y,
            depth,
        ];
        classify(solid, inner.local(paper), 0.) == Containment::Inside
    };
    let mut total = 0.;
    let height = 2. * reach / rows as f64;
    let width = 2. * reach / columns as f64;
    for row in 0..rows {
        let y = -reach + (row as f64 + 0.5) * height;
        let mut previous = (-reach - width, false);
        let mut entered = 0.;
        for column in 0..=columns + 1 {
            let x = -reach + column as f64 * width;
            let state = inside(x, y);
            if state != previous.1 {
                let (mut a, mut b) = (previous.0, x);
                for _ in 0..50 {
                    let mid = 0.5 * (a + b);
                    if inside(mid, y) == previous.1 {
                        a = mid;
                    } else {
                        b = mid;
                    }
                }
                if state {
                    entered = a;
                } else {
                    total += (a - entered) * height;
                }
            }
            previous = (x, state);
        }
    }
    total
}

#[test]
fn sections_are_closed_oriented_boundaries_of_the_membership_predicate() {
    for (name, solid) in solids() {
        for (label, frame) in views() {
            let body = built(&solid, frame);
            let inner = &body.0;
            let [lo, hi] = body.bounds();
            for fraction in [0.31, 0.5, 0.83] {
                let depth = lo[2] + (hi[2] - lo[2]) * fraction;
                let pieces = body.cut(depth).unwrap();
                assert!(!pieces.is_empty(), "{name} {label} {fraction}");
                let ends: Vec<[P3; 2]> =
                    pieces.iter().map(|p| [p.point(0.), p.point(1.)]).collect();
                for (i, piece) in pieces.iter().enumerate() {
                    for step in 0..=64 {
                        let t = f64::from(step) / 64.;
                        let p = piece.point(t);
                        assert!((p[2] - depth).abs() < 1e-9, "{name} {label}");
                        assert_eq!(
                            classify(&solid, inner.local(p), 1e-8),
                            Containment::Boundary,
                            "{name} {label} {fraction} {i} {t}"
                        );
                    }
                    // Every end continues in exactly one other start.
                    let next = ends
                        .iter()
                        .filter(|other| norm(sub(other[0], ends[i][1])) < 1e-7)
                        .count();
                    assert_eq!(next, 1, "{name} {label} {fraction} {i} {:?}", ends[i]);
                    // Material lies on the left of the travel direction.
                    let [a, b] = [piece.point(0.5 - 1e-4), piece.point(0.5 + 1e-4)];
                    let travel = unit(sub(b, a));
                    let mid = piece.point(0.5);
                    let side = |sign: f64| {
                        let p = [
                            mid[0] - sign * 2e-4 * travel[1],
                            mid[1] + sign * 2e-4 * travel[0],
                            depth,
                        ];
                        classify(&solid, inner.local(p), 0.)
                    };
                    assert_eq!(
                        side(1.),
                        Containment::Inside,
                        "{name} {label} {fraction} {i}"
                    );
                    assert_eq!(
                        side(-1.),
                        Containment::Outside,
                        "{name} {label} {fraction} {i}"
                    );
                }
                let exact = area(&pieces);
                let sampled = quadrature(&solid, &body, depth);
                assert!(
                    (exact - sampled).abs() < 1e-3 * sampled,
                    "{name} {label} {fraction} {exact} {sampled}"
                );
            }
        }
    }
}

#[test]
fn the_axial_plane_section_of_a_worm_is_its_straight_sided_profile() {
    let solid = worm(1., 20., 2.);
    let (source, _) = solid.rust_worm().unwrap();
    let thread = source.thread();
    let frame = PlaneFrame {
        origin: [0.; 3],
        u: [0., 0., 1.],
        v: [1., 0., 0.],
    };
    let body = built(&solid, frame);
    let pieces = body.cut(0.).unwrap();
    let slopes = [1, 3].map(|band| thread.slope_at_phase(0., band));
    let mut flanks = 0;
    for piece in &pieces {
        // Paper x is the axis and paper y the signed radius.
        let points: Vec<P3> = (0..=8).map(|i| piece.point(f64::from(i) / 8.)).collect();
        let [a, b] = [points[0], points[8]];
        for (i, p) in points.iter().enumerate() {
            let s = i as f64 / 8.;
            let chord = add(scaled(a, 1. - s), scaled(b, s));
            let off = (p[0] - chord[0]) * (b[1] - a[1]) - (p[1] - chord[1]) * (b[0] - a[0]);
            assert!(off.abs() < 1e-12 * 36., "{off}");
            assert!(p[2].abs() < 5e-12);
        }
        let radial = (b[1].abs() - a[1].abs()).abs();
        if radial > 1. && (b[0] - a[0]).abs() > 0.1 {
            // A flank: the radius changes by the profile slope per axial mm.
            let slope = (b[1].abs() - a[1].abs()) / (b[0] - a[0]);
            assert!(
                slopes.iter().any(|s| (s - slope).abs() < 1e-11 * s.abs()),
                "{slope} {slopes:?}"
            );
            flanks += 1;
        }
    }
    // Two flanks per turn on both sides of the axis, less the cut turns.
    assert!(flanks >= 4 * 5, "{flanks}");
}
